#!/usr/bin/env Rscript

## script to run coloc for Olink regions
## Maik Pietzner 15/12/2023
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("<path to dir>")

## --> packages required <-- ##

require(data.table)
require(susieR)
require(coloc)
require(doMC)
require(Rfast)

## --> import parameters <-- ##

olink <- args[1]
chr.s <- as.numeric(args[2])
pos.s <- as.numeric(args[3])
pos.e <- as.numeric(args[4])

cat("Run coloc with", olink, chr.s, pos.s, pos.e, "\n")

#-----------------------------------------#
##--     load regional assoc stats     --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Reading summary statistics in \n")

## import regional statistics with credible set assignment 
res.olink <- fread(paste("<path to dir>",
                         olink, chr.s, pos.s, pos.e, "txt", sep="."))
## recode X-chromosome
res.olink[, chr := ifelse(chr == "X", 23, chr)]

cat("Found", nrow(res.olink), "entries \n")
cat("------------------------------\n")

#-----------------------------------------#
##--            LD-matrix              --##
#-----------------------------------------#

## write list of SNPs to be queried to file 
write.table(res.olink$rsid, paste("tmpdir/snplist", olink, chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)

## snp data to be queried (positions are in build37, might differ from sentinel variant file with positions in build38)
tmp.z        <- as.data.frame(res.olink[, c("rsid", "chr", "pos", "NEA", "EA")])
names(tmp.z) <- c("rsid", "chromosome", "position", "allele1", "allele2")

## adopt chromosome if needed
if(chr.s < 10){
  print("tada")
  tmp.z$chromosome <- paste0("0", tmp.z$chromosome)
}else if(chr.s == 23){
  tmp.z$chromosome <- "X"
}

## check for input
print(head(tmp.z))

## write to file
write.table(tmp.z, paste("tmpdir/snpz", olink, chr.s, pos.s, pos.e, "z", sep="."), row.names = F, quote = F)

## --> create master file for LDstore2 <-- ##

## assign entries
m.file      <- data.frame(z=paste("tmpdir/snpz", olink, chr.s, pos.s, pos.e, "z", sep="."),
                          bgen=paste("tmpdir/filtered", olink, chr.s, pos.s, pos.e, "bgen", sep="."),
                          bgi=paste("tmpdir/filtered", olink, chr.s, pos.s, pos.e, "bgen.bgi", sep="."),
                          ld=paste("tmpdir/ld", olink, chr.s, pos.s, pos.e, "ld", sep="."),
                          incl="/sc-projects/sc-proj-computational-medicine/data/UK_biobank/genotypes/sample_inclusion/qctool_pass_EUR_panukbb_all_unrelated.incl",
                          n_samples=ifelse(chr.s != 23, 487409, 486757))

## write to file
write.table(m.file, paste("tmpdir/master", olink, chr.s, pos.s, pos.e, "z", sep="."), sep=";", row.names = F, quote = F)

cat("--------------------------------------\n")
cat("computing LD matrix\n")

## submit the job
system(paste("./scripts/get_LD_matrix_ldstore.sh", chr.s, pos.s, pos.e, olink))

## read in matrix
ld         <- fread(paste("tmpdir/ld", olink, chr.s, pos.s, pos.e, "ld", sep="."), data.table = F)
## print for checking: N.B.: This assumes alleles are aligned to UKB ones! (this was done during fine-mapping on the Cambridge server)
print(ld[1:5,1:5])

## create identifier column in results to keep the mapping to the LD matrix
res.olink[, snp.id := 1:nrow(res.olink)]

cat("Done\n")
cat("--------------------------------------\n")

#------------------------------------------#
##--     obtain top SNPs and proxies    --##
#------------------------------------------#

## align order
res.olink  <- as.data.table(res.olink)
## order, keep in mind that cs -1 means not included in the credible set
res.olink  <- res.olink[ order(-cs, -pip) ]
## create indicator to select
res.olink[, ind := 1:.N, by="cs"]
## get top SNPs
top.snp    <- res.olink[ ind == 1 & cs > 0]
## add LD columns: be careful 'cs' does not mean that numbers match, but is a legacy from fine-map
for(j in top.snp$cs){
  res.olink[, paste0("R2.", j)] <- ld[ res.olink$snp.id, top.snp$snp.id[which(top.snp$cs == j)]]^2
}

## get all proxy SNPs
proxy.snps <- paste0("R2.", top.snp$cs)
proxy.snps <- apply(res.olink[, ..proxy.snps], 1, function(x){
  ifelse(sum(x >= .8) > 0, T, F)
})
## may include SNPs not in the credible set
proxy.snps <- res.olink[proxy.snps]

## store the LD pattern across top SNPs (convert to data frame to ease downstream operations)
ld.top.snps <- as.data.frame(res.olink)
ld.top.snps <- lapply(top.snp$cs, function(x){
  ## get all SNPs and corresponding LD
  tmp        <- paste0("R2.", x)
  tmp        <- ld.top.snps[which(ld.top.snps[, tmp] >= .8), c("MarkerName", "id", "rsid", tmp)]
  ## edit names
  names(tmp) <- c("MarkerName.proxy", "id.proxy", "rsid.proxy", "R2")
  print(tmp)
  ## add top SNP
  tmp        <- merge(as.data.frame(top.snp[x, c("MarkerName", "id", "rsid", "cs")]), tmp, suffix=c(".lead", ".proxy"))
  ## do some renaming to ease downstream coding
  names(tmp) <- c("MarkerName.lead", "id.lead", "rsid.lead", "cs", "MarkerName.proxy", "id.proxy", "rsid.proxy", "R2.proxy")
  ## return
  return(tmp)
})
## combine everything
ld.top.snps <- do.call(rbind, ld.top.snps)

## write SNPs to file to be queried
write.table(unique(c("SNP", ld.top.snps$rsid.proxy)), paste0("tmpdir/", olink, ".snplist"), col.names = F, row.names = F, quote = F)

#-------------------------------------------#
##-- loop through phenotypes of interest --##
#-------------------------------------------#

## import phenotypes to be tested
finn.dis <- fread("input/FinnGen.MR.priority.20231129.txt")

## do in parallel
registerDoMC(5)

## run through: finn.dis$processed_data_name
res.mr   <- mclapply(finn.dis$processed_data_name, function(x){
  
  print(x)
  
  ## get the relevant associations (A2 is the effect allele)
  finn.stat <- fread(cmd=paste0("zgrep -wF -f tmpdir/", olink,".snplist <path to dir>/", x, ".tsv.gz"))
  
  ## proceed only if at least suggestive evidence
  if(nrow(finn.stat) > 0 & min(finn.stat$P, na.rm = T)  < 1e-4 ){
    
    ## import complete stats
    finn.stats <- paste0("zcat <path to dir>/", x, ".tsv.gz | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                         " '{if(($1 == chr && $9 >= low && $9 <= upp) || NR == 1) print $0}'")
    finn.stats <- data.table::fread(cmd = finn.stats, sep = "\t")
    ## make unique for some reason
    finn.stats <- unique(finn.stats)
    
    ## create MarkerName to enable mapping to the LD file
    finn.stats[, MarkerName := paste0(CHR, ":", BP, ":", pmin(A1, A2), "_", pmax(A1, A2))]
    
    #-----------------------------#
    ##-- combined set of stats --##
    #-----------------------------#
    
    ## combine
    res.all    <- merge(res.olink, finn.stats[, c("MarkerName", "A1", "A2", "BETA", "SE", "P")], by="MarkerName", suffixes = c(".pQTL", ".finngen"))
    
    ## align effect estimates (this currently ignores INDELs --> check!)
    res.all[, Effect.finn := ifelse(toupper(EA) == A2, BETA, -BETA)]
    res.all[, StdErr.finn := SE]
    
    #-----------------------------#
    ##--       naive coloc     --##
    #-----------------------------#
    
    ## import function to do so
    source("scripts/naive_coloc.R")
    
    ## run coloc
    res.naive  <- naive.coloc(res.all, ld, x)
    ## add type of coloc
    res.naive[, type := "naive"]
    
    #-----------------------------#
    ##--       susie coloc     --##
    #-----------------------------#
    
    ## import function to do so
    source("scripts/susie_coloc.R")
    
    ## run coloc
    res.susie  <- susie.coloc(res.all, ld, top.snp, x)
    ## add type of coloc
    if(!is.null(res.susie) > 0){
      res.susie[, type := "susie"]
    }  
    
    #-----------------------------#
    ##--    create return      --##
    #-----------------------------#
    
    ## combine
    res.coloc  <- plyr::rbind.fill(res.naive, res.susie)
    
    ## return results
    return(res.coloc)
  }
}, mc.cores = 5)
# })
print(res.mr)
## combine
res.mr <- rbindlist(res.mr, fill=T)

## store results
write.table(res.mr, paste0("output_coloc//coloc.results.", olink, ".txt"), sep="\t", row.names = F)

## do some cleaning
system(paste0("rm tmpdir/*.", olink, ".", chr.s, ".", pos.s, ".", pos.e,".*"))
