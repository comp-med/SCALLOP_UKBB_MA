#!/usr/bin/env Rscript
## get the arguments from the command line that are the protein to start with and which to finsih
args <- commandArgs(trailingOnly = T)
## Check we have given arguments
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  

## give path to the libraries
.libPaths(c( "/R/path/4.1", 
             "/opt/R/4.1.0/lib/R/library"))

## loadin the packages
require(data.table)
require(susieR)
require(doMC)
require(Rfast)
library(tidyverse)
library(unixtools)

set.tempdir("/Temp_LD_Dir/")

## Load in the BIM file
bim_file <- fread("/LD_REFERENCE/ukb.20k.LDref.bim")
## set column names 
colnames(bim_file) <- c("Chr", "rsid", "centimorgan", "bp", "Allele1", "Allele2")
## Getthe regions we want
file_list_processed <- list.files("/path_to_dir/Post_MA_Files/Processed/", recursive = T)
file_list_processed_reg <- file_list_processed[grep("Regions", file_list_processed)]
MA_reg <- data.frame()
for (i in file_list_processed_reg) {
  temp_sum <- fread(paste0("/Processed/", i))
  MA_reg <- rbind(MA_reg, temp_sum)
}
MA_reg <- MA_reg %>% filter(str_detect(protein_model, pattern = "FE1")) 
MA_reg <- MA_reg %>% mutate(region_comb = paste0(Chr, "-", Start, "-", End), 
                            remove =case_when(Chr == 6 & (Start > 23000000 & End <35000000)~ "Yes", T ~ "No")) %>%
  filter(remove == "No")


## tidy function to get LD matrix
ld_mat_tidy <- function (variants, bfile, plink_bin, with_alleles = TRUE) 
{
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                  "sh")
  fn <- tempfile()
  write.table(data.frame(variants), file = fn, row.names = F, 
              col.names = F, quote = F)
  fun1 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --extract ", shQuote(fn, 
                                                                      type = shell), " --make-just-bim ", " --keep-allele-order ", 
                 " --out ", shQuote(fn, type = shell))
  system(fun1)
  bim <- read.table(paste0(fn, ".bim"), stringsAsFactors = FALSE)
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --extract ", shQuote(fn, 
                                                                      type = shell), " --r square ", " --keep-allele-order ", 
                 " --out ", shQuote(fn, type = shell))
  system(fun2)
  res <- read.table(paste0(fn, ".ld"), header = FALSE) %>% 
    as.matrix
  if (with_alleles) {
    rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, 
                                            bim$V6, sep = "_")
  }
  else {
    rownames(res) <- colnames(res) <- bim$V2
  }
  ## tidy up after ourselves
  file.remove(list.files(tempdir(), pattern = ".ld", recursive = T, full.names = T))
  file.remove(list.files(tempdir(), pattern = ".bim", recursive = T, full.names = T))
  file.remove(list.files(tempdir(), pattern = ".log", recursive = T, full.names = T))
  file.remove(list.files(tempdir(), pattern = ".nosex", recursive = T, full.names = T))
  return(res)
}

region_prot <- unique(MA_reg$protein_model)
cat("Running ", as.numeric(args[2]) -as.numeric(args[1]), "Regions. \n This will begin with ", region_prot[as.numeric(args[1])], "and finish with ", region_prot[as.numeric(args[2])], "\n")

for (prot in as.numeric(args[1]): as.numeric(args[2])) {
  
  cat("Running....", region_prot[prot], "\n")
  ### create a directory to store our results
  out_dir = paste0(str_replace(region_prot[prot], pattern = "_FE1", replacement = ""), "/")
  mainDir_out = "/hpath_to_dir/Post_MA_Files/Susie_Res/"
  dir.create(file.path(mainDir_out, out_dir), showWarnings = FALSE)
  
  
  ## create a file path to load in the UKBB summary stats
  temp_path <- paste0("/path_to_dir/UKBB_SCALLOP_MA/", 
                      str_replace(region_prot[prot], pattern = "_FE1", replacement = ""), "/")
  ## get list of files for that protein
  file.list <- list.files(temp_path)
  ## Load in the UKBB summary stats
  UKBB_temp <- fread(paste0(temp_path,file.list[grep("^(?=.*txt.gz)(?!.*QC)", file.list, perl=TRUE)]))
  ## Extract the region we want to begin with 
  temp_region_dat <- MA_reg %>% filter(protein_model == region_prot[prot])
  
  for (i in 1:dim(temp_region_dat)[1]) {
    
    out_file_sub_dir <- paste0(temp_region_dat$Chr[i], "_", 
                               temp_region_dat$Start[i], "_", 
                               temp_region_dat$End[i])
    cat("Running...", out_file_sub_dir, "for ", out_dir)
    
    ## extract the SNPs within the region we want
    UKBB_temp_reg <- UKBB_temp %>% filter(CHROM == temp_region_dat$Chr[i] & 
                                            between(pos, temp_region_dat$Start[i], 
                                                    temp_region_dat$End[i])) %>%
      left_join(., bim_file, by = "rsid") %>% mutate(BETA_Susie = case_when(ALLELE1 == Allele1 ~ BETA, 
                                                                            T ~ -BETA))
    
    
    ## Get the LD matrix
    ld_temp_mat <- ld_mat_tidy(bfile = "/LD_REFERENCE/ukb.20k.LDref", 
                               plink_bin = "/plink",UKBB_temp_reg$rsid)
    
    ## Rename the columns/rows to match the dataset
    row.names(ld_temp_mat) <- data.frame(str_split_fixed(row.names(ld_temp_mat), "_", 2))[,1]
    colnames(ld_temp_mat) <- data.frame(str_split_fixed(colnames(ld_temp_mat), "_", 2))[,1]
    
    
    ## drop those with all missing
    ## format the LD matrix
    ld_temp_mat <- ld_temp_mat[,!colnames(ld_temp_mat) %in% names(which(colSums(is.na(as.matrix(ld_temp_mat))) > dim(ld_temp_mat)-1)) ]
    ld_temp_mat <- ld_temp_mat[complete.cases(ld_temp_mat),]
    
    ## Keep only the rows where we have SNPs and LD
    UKBB_temp_reg <- subset(UKBB_temp_reg, rsid %in% row.names(ld_temp_mat))
    
    ## align order
    UKBB_temp_reg <- as.data.table(UKBB_temp_reg)
    UKBB_temp_reg <- UKBB_temp_reg[ order(pos) ]
    ## order ld matrix accordingly
    ld_temp_mat        <- ld_temp_mat[ UKBB_temp_reg$rsid, UKBB_temp_reg$rsid]
    
    
    ## do in parallel
    registerDoMC(6)
    
    ## run fine mapping to obtain 95%-credible sets (throws an error due to possible rounding errors in LD matrix)
    set.seed(42)
    ## run through different values for L and catch possible errors
    res.fine <- mclapply(2:10, function(x){
      
      ## run SuSiE
      tmp <- tryCatch(
        {
          susie_rss(UKBB_temp_reg$BETA_Susie/UKBB_temp_reg$SE, as.matrix(ld_temp_mat), L = x, coverage = .95, min_abs_corr=.1, max_iter = 10000)
        }, error=function(e){
          return(list(pip=rep(NA, nrow(UKBB_temp_reg)),
                      sets=list(cs=NA),
                      converged=F))
        })
      ## add L
      tmp$L <- x
      ## return
      return(tmp)
      
    }, mc.cores=6)
    ## reduce to entries with at least one credible set
    jj       <- unlist(lapply(res.fine, function(x) length(na.omit(x$sets$cs))))
    ## delete
    res.fine <- res.fine[which(jj > 0)]
    
    
    ## proceed only if any
    if(length(res.fine) > 0){
      
      #------------------------------------------#
      ##--        prune sets if needed        --##
      #------------------------------------------#
      
      ## do LD assessment
      ld.top <- mclapply(res.fine, function(x){
        
        ## add PIPs to ease selection of SNPs
        tmp.olink    <- summary(x)$vars
        ## subset to those in credible sets
        tmp.olink    <- as.data.table(subset(tmp.olink, cs > 0))
        ## only top SNPs
        tmp.olink    <- tmp.olink[order(cs, -variable_prob)]
        ## get only the top
        tmp.olink[, ind := 1:.N, by="cs"]
        tmp.olink    <- tmp.olink[ ind == 1]
        ## get the names
        tmp.olink[, id := sapply(tmp.olink$variable, function(x) rownames(ld_temp_mat)[x])]
        ## generate LD
        top.ld       <- ld_temp_mat[tmp.olink$id, tmp.olink$id, drop=F]^2
        ## identify possible sets in LD (r2>0.25); set diagonal to zero to ease downstream analysis
        diag(top.ld) <- 0
        top.ld       <- reshape2::melt(top.ld)
        ## subset to possible problematic candidates
        top.ld       <- subset(top.ld, value >= .25)
        ## return
        return(top.ld)
        
      }, mc.cores=6)
      
      ## find the maximum set of unrelated variants
      jj       <- unlist(lapply(ld.top, nrow))  
      ## subset accordingly
      if(sum( jj > 0) > 0){
        res.fine <- res.fine[-which(jj > 0)]
      }
      
      ## add to the results, if any
      if(length(res.fine) > 0){
        
        dir.create(file.path(paste0(mainDir_out,"/", out_dir), out_file_sub_dir), showWarnings = FALSE)
        ## take only the last one
        res.fine  <- res.fine[[length(res.fine)]]
        
        ## add the information on pip and cs to the summary statistics
        tmp        <- summary(res.fine)$vars
        tmp$id     <- sapply(tmp$variable, function(x) rownames(ld_temp_mat)[x])
        names(tmp) <- c("variable", "pip", "cs", "rsid")
        UKBB_temp_reg  <- left_join(UKBB_temp_reg, tmp[, c("rsid", "pip", "cs")], by="rsid")
        
        ## write all results to file
        write.table(UKBB_temp_reg, paste0(file.path(paste0(mainDir_out,"/", out_dir), out_file_sub_dir),"/Results_with_PIP.txt"), sep="\t", row.names = F)
        
        ## write top credible sets to separate file (N.B. variants not in any credible set are indicated by -1)
        UKBB_temp_reg  <- as.data.table(subset(UKBB_temp_reg, cs > 0))
        ## create indicator of strongest
        UKBB_temp_reg  <- UKBB_temp_reg[ order(cs, -pip)]
        UKBB_temp_reg[, ind := 1:.N, by="cs"]
        
        ## write all results to file
        write.table(UKBB_temp_reg[ind == 1], paste0(file.path(paste0(mainDir_out,"/", out_dir), out_file_sub_dir),"/Top_Hits_with_PIP.txt"), sep="\t", row.names = F)
        
      }else{
        
        cat("\n-----------------------\n")  
        cat("Found no credible sets for",  out_file_sub_dir, "for ", out_dir, "\n")
        
        ## create empty file to be able to check for successful run jobs later on
        write.table(data.frame(id=NA, cs=NA), paste0(file.path(paste0(mainDir_out,"/", "Non_Converge_Check/")),gsub(x = paste0(out_dir,out_file_sub_dir ), pattern = "/", replacement = "_"),".txt"), sep="\t", row.names = F)
        
      }
      
    }else{
      
      cat("\n-----------------------\n")  
      cat("Found no credible sets for",  out_file_sub_dir, "for ", out_dir, "\n")
      
      ## create empty file to be able to check for successful run jobs later on
      write.table(data.frame(id=NA, cs=NA), paste0(file.path(paste0(mainDir_out,"/", "Non_Converge_Check/")),gsub(x = paste0(out_dir,out_file_sub_dir ), pattern = "/", replacement = "_"),".txt"), sep="\t", row.names = F)
      
    }
    
  }
}
