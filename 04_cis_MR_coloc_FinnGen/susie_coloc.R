####################################################
## function to compute colocalisation assuming more
## than one causal variant with susie coloc

susie.coloc <- function(res.all, ld, top.snp, x){
  
  ## 'res.all' -- data set containing merged and aligned statistics
  ## 'ld'      -- corresponding LD matrix
  ## 'top.snp' -- top snps for pQTL credible sets

  #-----------------------------------------#
  ##-- 	        prepare coloc            --##
  #-----------------------------------------#
  
  ## order by position
  res.all      <- as.data.table(res.all)
  res.all      <- res.all[order(pos)]
  
  ## edit names
  rownames(ld) <- colnames(ld)
  
  ## prepare input
  D1           <- list(beta=res.all$Effect, varbeta=res.all$StdErr^2, 
                       type="quant", 
                       sdY=1, 
                       N=max(res.all$TotalSampleSize, na.rm=T),
                       MAF=res.all$MAF,
                       snp=paste0("V", res.all$snp.id),
                       position=1:nrow(res.all),
                       ## subset the LD matrix to what is needed
                       LD=as.matrix(ld[res.all$snp.id, res.all$snp.id]))
  
  ## binary outcome
  D2          <- list(beta=res.all$Effect.finn, varbeta=res.all$StdErr.finn^2, 
                      type="cc",
                      # s=tr.info$ncase/(tr.info$ncontrol+tr.info$ncase), 
                      # N=tr.info$sample_size,
                      MAF=res.all$MAF,
                      snp=paste0("V", res.all$snp.id),
                      position=1:nrow(res.all),
                      ## subset the LD matrix to what is needed
                      LD=as.matrix(ld[res.all$snp.id, res.all$snp.id]))
  
  ## Olink protein
  set.seed(42)
  susie.olink <- tryCatch(
    {
      runsusie(D1, 
               # r2.prune = .25, 
               max_iter = 10000,
               L=ifelse(nrow(top.snp) == 1, 2, nrow(top.snp)))
    }, error=function(e){
      return(NA)
    })
  
  ## trait
  susie.trait <- tryCatch(
    {
      runsusie(D2, 
               # r2.prune = .25, 
               max_iter = 10000,
               L=5)
    }, error=function(e){
      return(NA)
    })
  
  ## check whether both have outcome
  if(length(susie.olink) > 1 & length(susie.trait) > 1){
    
    ## additional check, whether both traits have at least one credible set to be tested
    if(length(summary(susie.olink)$cs) > 0 & length(summary(susie.trait)$cs) > 0){
      
      #-----------------------------------------#
      ##-- 	            run coloc            --##
      #-----------------------------------------#
      
      ## run coloc with susie input
      res.coloc        <- coloc.susie(susie.olink, susie.trait, p12 = 5e-6)
      
      #-----------------------------------------#
      ##--   cross-check with fine-mapping   --##
      #-----------------------------------------#
      
      ## add LD with fine-mapped variants 
      res.coloc        <- as.data.table(res.coloc$summary)
      
      print(res.coloc)
      
      ## some processing to ease downstream mapping
      res.coloc[, hit1 := as.numeric(gsub("V", "", hit1))]
      res.coloc[, hit2 := as.numeric(gsub("V", "", hit2))]
      
      ## compute ld between selected lead hits
      res.coloc$ld.top <- apply(res.coloc[, c("hit1", "hit2"), drop=F], 1, function(k) ld[k[1], k[2]]^2)
      
      #-----------------------------------------#
      ##--     add additional information    --##
      #-----------------------------------------#
      
      ## add effect estimates (restrict to selected pQTLs, which means no effect estimates for the regional lead for the trait are taken forward)
      res.coloc                   <- merge(res.coloc, res.all[, c("snp.id", "MarkerName", "rsid", "EA", "NEA", "MAF", "Effect", "StdErr", "Pvalue", "pip", "cs", "Effect.finn",
                                                                  "StdErr.finn", "P")], 
                                           by.x = "hit1", by.y = "snp.id")
      
      ## add phecode
      res.coloc$olink             <- olink
      ## add trait
      res.coloc$finn.id           <- x
      
      ## add rsids and LD with other variants
      res.coloc                   <- merge(res.coloc, res.olink[, c("snp.id", "rsid")], by.x="hit2", by.y = "snp.id", suffixes = c(".olink", ".finn"))
      
      #-----------------------------------------#
      ##-- 	        draw selected            --##
      #-----------------------------------------#
      
      if(sum(res.coloc$PP.H4.abf > .7) > 0 | sum(res.coloc$ld.top > .8) > 0){
        
        source("scripts/plot_locus_compare.R")
        png(paste("graphics/susie", olink, x, chr.s, pos.s, pos.e, "png", sep="."), width=16, height=8, units="cm", res=200)
        par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
        ## more complex layout for gene assignment
        layout(matrix(c(1,1,1,2,3,4),3,2), heights = c(.43,.37,.2))
        plot.locus.compare(res.all, res.coloc, ld, a.vars=unique(c(res.coloc$rsid.olink, res.coloc$rsid.finn)))
        dev.off()
      }
      
      ## do some renaming
      names(res.coloc)            <- c("snp.id.protein", "snp.id.finn", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "idx1", "idx2",
                                       "ld.top", "MarkerName", "rsid.protein", "EA", "NEA", "MAF", "Effect.protein", "StdErr.protein", "Pvalue.protein", "pip", "cs", 
                                       "Effect.finn", "StdErr.finn", "Pvalue.finn",  "olink", "finn.id", "rsid.finn")
      
      ## add R2 with top credible set
      res.coloc                   <- merge(res.coloc, res.all[, c("MarkerName", grep("R2", names(res.all), value=T)), with=F], by="MarkerName")
      
      ## write results to file
      return(res.coloc)
      
    }else{
      return(NULL)
    }
    
  }else{
    return(NULL)
  }
  
}
  