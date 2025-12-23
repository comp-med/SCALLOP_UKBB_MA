#!/usr/bin/env Rscript

## script to run MR for SCALLOP pQTLs and FinnGen stats
## Maik Pietzner 23/11/2023
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("<path to dir>")

## --> packages required <-- ##

require(data.table)
require(doMC)
require(MendelianRandomization)

## --> import parameters <-- ##

## get the ID of the protein
olink.id <- args[1]

cat("Run MR for", olink.id, "\n")

#-----------------------------------------#
##--       load fine-mapped pQTLs      --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Reading summary statistics in \n")

## import regional statistics with credible set assignment 
res.fine  <- fread("<path to dir>/SCALLOP.pQTLs.MR.QCed.updated.20240125.txt")
## keep only those needed
res.fine  <- res.fine[ MA_prot_id == olink.id]
## drop secondary signals in trans
res.fine  <- res.fine[ cis_trans == "Cis" | signal.group == 1]

## import proxy snps (SNP_A is the respective pQTL and SNP_B are all proxies; Allele1 is the effect allele)
proxy.snp <- fread(paste0("<path to dir>/", olink.id, ".txt"))
## subset to set of SNPs QC'ed
proxy.snp <- proxy.snp[ MarkerName_SNP_A %in% res.fine$MarkerName ]

## create signal column
proxy.snp <- as.data.table(proxy.snp)
proxy.snp <- proxy.snp[order(SNP_A, -R2)]
proxy.snp[, r2.order := 1:.N, by="SNP_A"]

## repair poor coding at the Cam end
proxy.snp[, rsid := ifelse(rsid == "SNP_B", SNP_B, rsid)]

## write SNPs to file to be queried
write.table(c("SNP", proxy.snp$rsid), paste0("tmpdir/", olink.id, ".snplist"), col.names = F, row.names = F, quote = F)

cat("------------------------------\n")

#-------------------------------------------#
##-- loop through phenotypes of interest --##
#-------------------------------------------#

## import phenotypes to be tested
finn.dis <- fread("input/FinnGen.MR.priority.20231129.txt")
## import broader definition
finn.tmp <- fread("<path to dir>/Summary.FinnGen.Release8.Phenotypes.20231129.txt")

## do in parallel
registerDoMC(10)

## run through
res.mr   <- mclapply(finn.dis$processed_data_name, function(x){
  
  print(x)
  
  ## get the relevant associations (A2 is the effect allele)
  finn.stat <- fread(cmd=paste0("zgrep -wF -f tmpdir/", olink.id,".snplist <path to dir>", x, ".tsv.gz"))
  
  ## only if enough SNPs found
  if(nrow(finn.stat) > 0){
    
    ## merge with proxy SNPs
    finn.stat <- merge(proxy.snp, finn.stat, by.x="rsid", by.y="SNP", suffixes = c(".scallop", ".finn"))
    finn.stat <- as.data.table(finn.stat)
    
    ## keep only best proxy for each lead SNP (SNP_A)
    finn.stat <- finn.stat[ order(SNP_A, -R2)]
    finn.stat[, ind := 1:.N, by="SNP_A"]
    finn.stat <- finn.stat[ ind == 1]
    ## add cis/trans information
    finn.stat <- merge(finn.stat, res.fine[, c("MarkerName", "cis_trans", "count", "interaction.trans", "num.efo.parent", "cis.pav.flag")], by.x="MarkerName_SNP_A", by.y="MarkerName")
    
    ## align effect estimates (this currently ignores INDELs --> check!)
    finn.stat[, Effect.finn := ifelse(toupper(Allele1) == A2, BETA, -BETA)]
    
    ## code protein increasing (to have valid MR-Egger estimates)
    finn.stat[, Effect.finn := sign(Effect)*Effect.finn]
    finn.stat[, Effect.pqtl := abs(Effect)]
    
    ## --> check for reverse confounding <-- ##
    
    ## pQTLs from here https://rdrr.io/github/MRCfinn/TwoSampleMR/src/R/add_rsq.r
    finn.stat[, r2.pqtl := apply(finn.stat[, c("Effect.pqtl", "StdErr", "TotalSampleSize"), with=F], 1, function(k){
      ## F-value
      Fval <- (k[1]/k[2])^2
      ## estimated variance
      R2   <- Fval / (k[3]-2+Fval)
      return(sqrt(R2) * sign(k[1]))
    })]
    
    
    ## import function to do so (adopted from TwoSampleMR package)
    source("scripts/compute_rsq_binary.R")
    ## new column for allele frequencies according to effect allele
    finn.stat[, eaf.pqtl := ifelse(Effect > 0, Freq1, 1 - Freq1)]
    ## continuous
    finn.stat[, r2.finn := apply(finn.stat[, c("Effect.finn", "eaf.pqtl"), with=F], 1, function(k){
      ## case number
      ncase      <- finn.tmp$n_cases[ which(finn.tmp$processed_data_name == x)]
      ## control number
      ncontrol   <- finn.tmp$n_control[ which(finn.tmp$processed_data_name == x)]
      ## prevalence
      prevalence <- ncase/ncontrol
      ## compute
      rsq        <- get_r_from_lor(k[1], k[2], ncase, ncontrol, prevalence)
      return(rsq)
    })]
    ## account for very low explained variance
    finn.stat[, r2.finn := ifelse(r2.finn < 0, 1e-4, r2.finn)]
    ## add effective case/control sample size
    finn.stat[, n_effective.finn := 2 / (1/finn.tmp$n_cases[ which(finn.tmp$processed_data_name == x)] + 1/finn.tmp$n_control[ which(finn.tmp$processed_data_name == x)])]
    ## compute index
    finn.stat[ , steiger.pval := apply(finn.stat[, c("r2.pqtl", "r2.finn", "TotalSampleSize", "n_effective.finn"), with=F], 1, function(k){
      ## perform test
      st <- psych::r.test( n = k[3], n2 = k[4], r12 = sqrt(k[1]), r34 = sqrt(k[2]))
      ## return p-value
      return(stats::pnorm(-abs(st$z)) * 2)
    })]
    
    ## drop possible SNPs (do not drop cis SNPs)
    finn.stat <- finn.stat[ !(r2.finn > r2.pqtl & steiger.pval < 1e-3 & cis_trans == "Trans")]
    
    ## proceed only if anything is left
    if(nrow(finn.stat) > 0){
      
      ## loop through different sets of instruments
      res      <- lapply(c("cis", "cis_pav", "trans", "trans_spec", "trans_inter", "cis_trans", "cis_trans_spec", "cis_trans_inter", "cis_trans_pleio", "trans_pleio"), function(k){
        
        ## create MR input
        if(k == "cis"){
          tmp <- finn.stat[ cis_trans == "Cis" ]
        }else if(k == "cis_pav"){
          tmp <- finn.stat[ cis_trans == "Cis" & cis.pav.flag == F ]
        }else if(k == "trans"){
          tmp <- finn.stat[ cis_trans == "Trans" ]
        }else if(k == "trans_spec"){
          tmp <- finn.stat[ cis_trans == "Trans" & count <= 5 ]
        }else if(k == "cis_trans"){
          tmp <- finn.stat
        }else if(k == "cis_trans_spec"){
          tmp <- finn.stat[ cis_trans == "Cis" | (cis_trans == "Trans" & count <= 5)]
        }else if(k == "cis_trans_pleio"){
          tmp <- finn.stat[ cis_trans == "Cis" | (cis_trans == "Trans" & count <= 5 & num.efo.parent <= 5)]
        }else if(k == "trans_inter"){
          tmp <- finn.stat[ interaction.trans == T]
        }else if(k == "trans_pleio"){
          tmp <- finn.stat[ cis_trans == "Trans" & count <= 5 & num.efo.parent <= 5]
        }else{
          tmp <- finn.stat[ interaction.trans == T | cis_trans == "Cis"]
        }
        
        print(k)
        
        ## --> preform LOO to identify influential variants (>1 variants) <-- ##
        
        ## do only if >1 instrument
        if(nrow(tmp) > 2){
          ## run LOO using simple IVW
          res.loo              <- lapply(1:nrow(tmp), function(c){
            
            ## perpare input
            mr.input             <- mr_input(bx = tmp[-c]$Effect.pqtl, bxse = tmp[-c]$StdErr,
                                             by = tmp[-c]$Effect.finn, byse = tmp[-c]$SE)
            ## run ivw
            tmp.ivw              <- mr_ivw(mr.input)
            ## return relevant effect stats
            return(data.table(loo.MarkerName = tmp$MarkerName_SNP_A[c], 
                              beta=tmp.ivw@Estimate, se=tmp.ivw@StdError, pval=tmp.ivw@Pvalue, 
                              q2=tmp.ivw@Heter.Stat[1], pval.q2=tmp.ivw@Heter.Stat[2]))
          })
          ## combine
          res.loo              <- rbindlist(res.loo, fill=T)
          
          ## evidence for excess heterogeneity
          res.loo[, het.excess := q2 < (median(q2)-3*sd(q2))]
          
          ## remove SNP(s) causing excess heterogeneity
          if(sum(res.loo$het.excess) > 0){
            tmp <- tmp[-which(MarkerName_SNP_A %in% res.loo[ het.excess == T]$loo.MarkerName)]
          }
          
        }else{
          ## dummy results to ease storing the results later on
          res.loo <- tmp
          res.loo[, het.excess := F]
        }
        
        ## proceed only if at least one SNP
        if(nrow(tmp) %in% c(1,2)){
          
          ## create MR input
          mr.input             <- mr_input(bx = tmp$Effect.pqtl, bxse = tmp$StdErr,
                                           by = tmp$Effect.finn, byse = tmp$SE)
          
          ## add Q2 estimate
          tmp.mr                <- mr_ivw(mr.input)
          
          ## create results data frame
          mr.results            <- data.frame(beta.IVW = tmp.mr@Estimate, se.IVW = tmp.mr@StdError, cil.IVW = tmp.mr@CILower, ciu.IVW = tmp.mr@CIUpper, pval.IVW = tmp.mr@Pvalue, 
                                              nsnps.loo=nrow(tmp), nsnps.all=nrow(res.loo),
                                              outlying.snps = paste(res.loo[het.excess == T]$loo.MarkerName, collapse = "|"))
          
          ## add to results
          mr.results$q2         <- tmp.mr@Heter.Stat[1]
          mr.results$pval.q2    <- tmp.mr@Heter.Stat[2]
          
          ## add type of MR
          mr.results$type       <- k
          ## add protein and outcome
          mr.results$MA_prot_id <- olink.id
          mr.results$finn.id     <- x
          ## plotting indicator
          mr.results$plot.ind   <- mr.results$pval.IVW < .05
          
          ## return results; only if independent
          if(!is.na(mr.results$se.IVW)){
            return(mr.results)
          }
          
        }else if(nrow(tmp) > 2){
          
          ## create MR input
          mr.input                 <- mr_input(bx = tmp$Effect.pqtl, bxse = tmp$StdErr,
                                               by = tmp$Effect.finn, byse = tmp$SE)
          
          ## perform MR
          mr.results               <- mr_allmethods(mr.input, method = "main") 
          
          ## convert into data frame 
          mr.results               <- mr.results@Values
          ## rename 
          names(mr.results)        <- c("Method", "beta", "se", "cil", "ciu", "pval")
          mr.results$ind           <- 1
          ## indicator for plotting
          plot.ind                 <- min(mr.results$pval, na.rm=T) < 0.05
          ## reshape
          mr.results               <- reshape(mr.results, timevar = "Method", direction = "wide", idvar = "ind")
          ## omit brackets from the output
          names(mr.results)        <- gsub("\\(|\\)", "", names(mr.results))
          names(mr.results)        <- gsub(" |-", ".", names(mr.results))
          
          ## add Q2 estimate
          tmp.mr                   <- mr_ivw(mr.input)
          ## add to results
          mr.results$q2            <- tmp.mr@Heter.Stat[1]
          mr.results$pval.q2       <- tmp.mr@Heter.Stat[2]
          
          ## add type of MR
          mr.results$type          <- k
          ## add protein and outcome
          mr.results$MA_prot_id    <- olink.id
          mr.results$finn.id       <- x
          ## how many snps (account for possibly removed SNPs)
          mr.results$nsnps.all     <- nrow(res.loo)
          mr.results$nsnps.loo     <- nrow(tmp)
          ## plotting indicator
          mr.results$plot.ind      <- plot.ind
          
          ## flag potentially outlying SNPs
          mr.results$outlying.snps <- paste(res.loo[het.excess == T]$loo.MarkerName, collapse = "|")
          
          ## return results
          return(mr.results)
          
        }
      })
      ## combine everything
      res      <- do.call(plyr::rbind.fill, res)
      
      #----------------------------#
      ##--     possible plot    --##
      #----------------------------#
      
      if(nrow(res) > 0){
        if(sum(res$plot.ind) > 0){
          ## open device c("cis", "trans", "trans_spec", "cis_trans", "cis_trans_spec")
          png(paste0("graphics/", olink.id, ".", x, ".png"), width = 16, height = 11.6, units = "cm", res = 300)
          ## parameters
          par(mar=c(1.5,1.5,.5,.5), cex.lab=.5, cex.axis=.5, mgp=c(.6,0,0), tck=-.01, mfrow=c(2,3))
          
          ## loop through
          for(k in c("cis", "trans", "trans_spec", "trans_pleio", "cis_trans_spec", "cis_trans_pleio")){
            
            ## empty plot
            plot(c(0, max(finn.stat$Effect.pqtl + 1.96*finn.stat$StdErr, na.rm=T)),
                 c(min(finn.stat$Effect.finn - 1.96*finn.stat$SE, na.rm=T), max(finn.stat$Effect.finn + 1.96*finn.stat$SE, na.rm=T)),
                 xlab="Effect protein", ylab="Effect outcome", xaxt="n", yaxt="n", type="n")
            ## axis
            axis(1, lwd=.5); axis(2, lwd=.5); abline(v=0, lwd=.5); abline(h=0, lwd=.5)
            ## title
            mtext(k, cex=.5)
            
            ## get data to be plotted
            if(k == "cis"){
              tmp <- finn.stat[ cis_trans == "Cis" ]
            }else if(k == "trans"){
              tmp <- finn.stat[ cis_trans == "Trans" ]
            }else if(k == "trans_spec"){
              tmp <- finn.stat[ cis_trans == "Trans" & count <= 5 ]
            }else if(k == "cis_trans"){
              tmp <- finn.stat
            }else if(k == "cis_trans_spec"){
              tmp <- finn.stat[ cis_trans == "Cis" | (cis_trans == "Trans" & count <= 5)]
            }else if(k == "cis_trans_pleio"){
              tmp <- finn.stat[ cis_trans == "Cis" | (cis_trans == "Trans" & count <= 5 & num.efo.parent <= 5)]
            }else if(k == "trans_inter"){
              tmp <- finn.stat[ interaction.trans == T]
            }else if(k == "trans_pleio"){
              tmp <- finn.stat[ cis_trans == "Trans" & count <= 5 & num.efo.parent <= 5]
            }else{
              tmp <- finn.stat[ interaction.trans == T | cis_trans == "Cis"]
            }
            
            ## plot if any
            if(nrow(tmp) > 0){
              ## confidence intervals
              arrows(tmp$Effect.pqtl - 1.96*tmp$StdErr, tmp$Effect.finn, tmp$Effect.pqtl + 1.96*tmp$StdErr, tmp$Effect.finn,
                     lwd=.3, length = 0)
              arrows(tmp$Effect.pqtl, tmp$Effect.finn - 1.96*tmp$SE, tmp$Effect.pqtl, tmp$Effect.finn + 1.96*tmp$SE,
                     lwd=.3, length = 0)
              ## point estimates
              points(tmp$Effect.pqtl, tmp$Effect.finn, cex=.7, bg=ifelse(tmp$cis_trans == "Cis", "#F95700", "#00A4CC"), lwd=.3, pch=21)
              
              ## add MR estimates
              tmp <- subset(res, type == k)
              if(nrow(tmp) > 0){
                for(l in c("IVW", "Simple.median", "Weighted.median")){
                  ## test if included
                  if(length(grep(l, names(tmp))) > 0){
                    if(!is.na(tmp[, paste0("beta.", l)])) abline(a=0, b=tmp[, paste0("beta.", l)], lwd=.5, lty=2)
                  }
                }
                ## MR-Egger
                if(length(grep("intercept", names(tmp))) > 0){
                  if(!is.na(tmp$beta.intercept)) abline(a=tmp$beta.intercept, b=tmp$beta.MR.Egger, lty=3, lwd=.5)
                } 
              }
            }
          }
          dev.off()
        }
      }
      
      ## return
      return(res)
      
    }
  }
}, mc.cores = 10)
# })
## combine everything
res.mr      <- rbindlist(res.mr, fill = T)

## store results
write.table(res.mr, paste0("output/MR.results.", olink.id, ".txt"), sep="\t", row.names = F)
