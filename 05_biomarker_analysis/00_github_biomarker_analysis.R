########################################
#### protein - disease associations ####
#### Maik Pietzner       25/05/2023 ####
########################################

rm(list=ls())
setwd("<path to dir>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(arrow)
require(doMC)
require(igraph)
require(arrow)

####################################
####      import covariates     ####
####################################

## An example list of columns:
cl.select       <- c("f.eid", "f.21022.0.0", "f.31.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", "f.20116.0.0", "f.1558.0.0", paste0("f.22009.0.", 1:10))
## names to assign
cl.names        <- c("f.eid", "age", "sex", "baseline_date", "centre", "bmi", "smoking", "alcohol", paste0("pc", 1:10))
## import data from the main release
ukb.dat         <- read_parquet("<path to dir>", col_select = cl.select)
## change names
names(ukb.dat)  <- cl.names

## define factor label for smoking
ukb.dat$smoking <- as.character(ukb.dat$smoking)
ukb.dat$smoking <- ifelse(ukb.dat$smoking == "Prefer not to answer", "Previous", ukb.dat$smoking)
ukb.dat$smoking <- factor(ukb.dat$smoking, levels = c("Never", "Previous", "Current"))

## define factor level for alcohol
ukb.dat$alcohol <- as.character(ukb.dat$alcohol)
ukb.dat$alcohol <- ifelse(ukb.dat$alcohol == "Prefer not to answer", "Three or four times a week", ukb.dat$alcohol)
ukb.dat$alcohol <- factor(ukb.dat$alcohol, levels = c("Never", "Special occasions only", "One to three times a month", "Once or twice a week", "Three or four times a week",
                                                      "Daily or almost daily"))
## convert to data table
ukb.dat         <- as.data.table(ukb.dat)

## --> get cystatin C <-- ##

## import QC'ed biomarker levels
ukb.bio         <- read_parquet("/<path to dir>", col_select=c("f.eid", "log_cysc_cleaned"))
## add to the data
ukb.dat         <- merge(ukb.dat, ukb.bio)

####################################
####     import protein data    ####
####################################

## import protein data
ukb.prot                <- fread("<path to dir>")
## import label
lab.prot                <- fread("<path to dir>")

## apply normalization
ukb.prot                <- as.data.frame(ukb.prot)
ukb.prot[, lab.prot$id] <- apply(ukb.prot[, lab.prot$id], 2, function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
})

## subset ukb data to same number of participants
ukb.dat                 <- ukb.dat[ f.eid %in% ukb.prot$eid]
gc(reset=T)

####################################
####     import phecode data    ####
####################################

## import
ukb.phe <- fread("<path to dir>", sep="\t", header = T)
## get a good right censoring date
quantile(ukb.phe$date, c(.05,.5,.975,.99), type=1) ## 2020-08-31 as end date

## subset to people with protein data
ukb.phe <- ukb.phe[ f.eid %in% ukb.prot$eid ]
## convert to wide format
ukb.phe <- dcast(ukb.phe, f.eid ~ phecode, value.var = c("date", "resource"))
## keep only what is of immediate interest
ukb.phe <- ukb.phe[, c("f.eid", grep("date", names(ukb.phe), value=T)), with=F]
## import death dates for censoring
ukb.dod <- fread("<path to dir>/death.txt")
## add date of death (keep only first instance)
ukb.phe <- merge(ukb.phe, ukb.dod[ ins_index == 0 , c("eid", "date_of_death")], all.x=T, by.x="f.eid", by.y="eid")

## import label
lab.phe <- fread("<path to dir>", sep="\t", header=T)
## add identifier to match with data set
lab.phe[, id := paste0("date_", phecode)]
## subset to what is in the data
lab.phe <- lab.phe[ id %in% names(ukb.phe) ]
## add sex-specific coding
lab.tmp <- fread("<path to dir>")
## create variable to map to UKBB data set
lab.tmp[, id := paste0("date_", phecode)]
## add sex-specific information to the data
lab.phe <- merge(lab.phe, lab.tmp[, c("id", "sex")], by="id")
## replace missing ones
lab.phe[, sex := ifelse(sex == "", "Both", sex)]

## convert to data frame to ease coding
ukb.phe <- as.data.frame(ukb.phe)
## add non-cases and baseline data
ukb.phe <- merge(ukb.phe, ukb.dat[, c("f.eid", "baseline_date")], all=T)
## convert death date
ukb.phe$date_of_death <- as.IDate(ukb.phe$date_of_death, format = "%d/%m/%Y")
ukb.phe$baseline_date <- as.IDate(ukb.phe$baseline_date)

## create case and date variable: choose
for(j in lab.phe$id){
  ## define case status (binary)
  ukb.phe[, gsub("date", "bin", j)]  <- ifelse(is.na(ukb.phe[, j]), 0, 
                                               ifelse(ukb.phe[, j] > ukb.phe[, "baseline_date"], 1, NA))
  ## define new date
  ukb.phe[, gsub("date", "surv", j)] <- as.IDate(ifelse(is.na(ukb.phe[, j]) & is.na(ukb.phe[, "date_of_death"]), as.IDate("2020-08-31"), 
                                                        ifelse(is.na(ukb.phe[, j]) & !is.na(ukb.phe[, "date_of_death"]), ukb.phe[, "date_of_death"], ukb.phe[, j])))
}
## add sex to count cases
ukb.phe <- merge(ukb.phe, ukb.dat[, c("f.eid", "sex")])

## compute case numbers
lab.phe[, inc_cases := apply(lab.phe[, c("id", "sex")], 1, function(x){
  if(x[2] == "Both"){
    sum(ukb.phe[, gsub("date", "bin", x[1])] == 1, na.rm=T)
  }else{
    sum(ukb.phe[ukb.phe$sex == x[2], gsub("date", "bin", x[1])] == 1, na.rm=T)
  }
})]

## at least 50 incident cases
lab.phe <- lab.phe[ inc_cases >= 50]

## keep only what is really needed
ukb.phe <- ukb.phe[, c("f.eid", "baseline_date", paste0("bin_", lab.phe$phecode), paste0("surv_", lab.phe$phecode))]
gc(reset=T)

####################################
####   export for association   ####
####################################

## covariates
fwrite(ukb.dat, "UKB.covariates.20240624.txt", sep="\t", row.names=F, na = NA)
## proteins
fwrite(ukb.prot, "UKB.proteins.20240624.txt", sep="\t", row.names=F, na = NA)
## phecodes
fwrite(ukb.phe, "UKB.phecodes.20240624.txt", sep="\t", row.names=F, na = NA)
## label for phecodes
fwrite(lab.phe, "label.phecodes.20240624.txt", sep="\t", row.names=F, na=NA, quote = F)
fwrite(lab.prot, "label.Olink.proteins.20240624.txt", sep="\t", row.names=F, na=NA)

####################################
####       import results       ####
####################################

#----------------------#
##--  first version --##
#----------------------#

## check for all results files
ii <- dir("../output/")
ii <- grep("surv.crude", ii, value = T)
## seems to be all there

## import in parallel
registerDoMC(12)
## run the import; careful, some have not run through!; also, drop outcome w/o effect estimates for early years
res.surv <- mclapply(ii, function(x){
  ## import
  tmp <- fread(paste0("../output/", x))
  ## report only if non-missing
  ii  <- grep("Error", tmp)
  if(length(ii) == 0){
    return(tmp)
  }
}, mc.cores = 12)
## combine
res.surv <- rbindlist(res.surv, fill = T)
## add phecode and protein information
res.surv <- merge(lab.phe[, c("phecode", "category", "cl", "phenotype")], res.surv)
## drop disease examlpes not fully covered
ii       <- table(res.surv$phecode)
ii       <- names(ii[ ii == 1463])
res.surv <- subset(res.surv, phecode %in% ii)
## drop one additional associations with inconsistent time outcomes
res.surv <- as.data.table(res.surv)
res.surv <- res.surv[ se.2yr != 0 & se.5yr != 0 & se.10yr != 0 & se.18yr != 0]

#######################################################################################################
#######################################################################################################
####                                  REVISION - CELL 02/07/2025                                   ####
#######################################################################################################
#######################################################################################################

######################################################
####           import genetic results             ####
######################################################

## --> import MVP results from scratch <-- ##

## poorly mapped phecodes to FinnGen
bad.mvp.finngen <- fread("FinnGen_R7_BAD_MAPPING.txt")
## updated MVP: Mendelian randomizaiton (be careful; this omits associations not meeting significance!!)
res.mvp.mr.upd  <- fread("SCALLOP_MR_rerun_significant_hits_16MAY2024.txt")
## updated MVP: Colocalisation
res.mvp.cc.upd  <- fread("scallop_rerun_coloc_29MAY2024.txt")
## import original results
res.mvp.mr      <- fread("aggregated_SCALLOP_MR_28JAN2024.txt")
## drop bad phecodes
res.mvp.mr      <- res.mvp.mr[ !(PHENOTYPE %in% bad.mvp.finngen$Phenotype)]
## add updated results
res.mvp.mr      <- rbind(res.mvp.mr, res.mvp.mr.upd)
## import original phecode coloc results
res.mvp.cc      <- fread("phecode_scallop_coloc.txt")
## combine and filter
res.mvp.cc      <- res.mvp.cc[  !(Phenotype %in% bad.mvp.finngen$Phenotype)]
res.mvp.cc      <- rbind(res.mvp.cc, res.mvp.cc.upd, fill = T)
res.mvp.cc[, GENE_NAME_SeqID := NULL]
## make unique
res.mvp.cc      <- unique(res.mvp.cc)

## combine
res.mvp.mr      <- res.mvp.mr[ order(GENE_NAME, PHENOTYPE, -FDR)]
res.mvp.mr[, ind := 1:.N, by = c("GENE_NAME", "PHENOTYPE")]
res.gen         <- merge(res.mvp.mr[ ind == 1, .(PHENOTYPE, GENE_NAME, BETA, SE, PVAL, FDR)], res.mvp.cc, 
                         by.x = c("PHENOTYPE", "GENE_NAME"), by.y = c("Phenotype", "gene_name"), all.x = T)
## drop non-phecode results
res.gen         <- res.gen[ substr(PHENOTYPE, 1, 3) == "phe"]

## make comparable
res.gen[, id := tolower(GENE_NAME)]
res.gen[, phecode := as.numeric(gsub("_", ".", gsub("phe_", "", PHENOTYPE)))]

## make unique (some duplications due to coloc)
res.gen         <- res.gen[ order(GENE_NAME, PHENOTYPE, -FDR)]
res.gen[, ind := 1:.N, by = c("GENE_NAME", "PHENOTYPE")]
## drop
res.gen         <- res.gen[ ind == 1]

## import protein label from SCALLOP
prot.scallop    <- fread("<path to dir>/MA_protein_list_and_locations.txt")

######################################################
####      import results from extended models     ####
######################################################

#----------------------#
##--   adj version  --##
#----------------------#

## check for all results files
ii       <- dir("../output/")
ii       <- grep("crude", ii, value = T, invert = T)
## seems to be all there

## import in parallel
registerDoMC(12)
## run the import; careful, some have not run through!; also, drop outcome w/o effect estimates for early years
res.adj  <- mclapply(ii, function(x){
  ## import
  tmp <- fread(paste0("../output/", x))
  ## report only if non-missing
  ii  <- grep("Error", tmp)
  if(length(ii) == 0){
    return(tmp)
  }
}, mc.cores = 12)
## combine
res.adj  <- rbindlist(res.adj, fill = T)
## add phecode and protein information
res.adj  <- merge(lab.phe[, c("phecode", "category", "cl", "phenotype")], res.adj)
## drop disease examples not fully covered
ii       <- table(res.adj$phecode)
ii       <- names(ii[ ii == 1463])
## n = 682 out of 844
res.adj  <- subset(res.adj, phecode %in% ii)
## drop one additional associations with inconsistent time outcomes
res.adj  <- as.data.table(res.adj)
res.adj  <- res.adj[ se.2yr != 0 & se.5yr != 0 & se.10yr != 0 & se.18yr != 0]
## n = 984,223

## check again for completeness
ii       <- table(res.adj$phecode)
ii       <- names(ii[ ii == 1463])
## subset accordingly
res.adj  <- res.adj[ phecode %in% ii]

######################################################
####            look at prevalent diseases        ####
######################################################

#--------------------------------------#
##--    prepare data accordingly    --##
#--------------------------------------#

## create new variables
ukb.phe <- as.data.frame(ukb.phe)

## loop through all endpoints accordingly
for(j in lab.phe$id){
  ## create variable for people with the disease at baseline
  ukb.phe[, paste0("prev_", gsub("date_", "",j))] <- ifelse(is.na(ukb.phe[, paste0("bin_", gsub("date_", "",j))]), 1, 0)
}

## convert back to data table
ukb.phe <- as.data.table(ukb.phe)

## export for association analysis
fwrite(ukb.phe, "UKB.phecodes.prevalence.20250702.txt", sep="\t", row.names=F, na = NA)

#--------------------------------------#
##--          import results        --##
#--------------------------------------#

## check for all results files
ii <- dir("../output/")
ii <- grep("prev", ii, value = T)

## import in parallel
registerDoMC(12)
## run the import; careful, some have not run through!; also, drop outcome w/o effect estimates for early years
res.prev <- mclapply(ii, function(x){
  ## import
  tmp <- fread(paste0("../output/", x))
  ## report only if non-missing
  ii  <- grep("Error", tmp)
  if(length(ii) == 0){
    return(tmp)
  }
}, mc.cores = 12)
## combine
res.prev <- rbindlist(res.prev, fill = T)
## restrict to outcomes also in incident diseases
res.prev <- res.prev[ phecode %in% res.surv$phecode]
## add phecode and protein information
res.prev <- merge(lab.phe[, c("phecode", "category", "cl", "phenotype")], res.prev)

## adopt names
names(res.prev) <- c("phecode", "category", "cl", "phenotype", "id", "beta.prev", "se.prev", "pval.prev", "nevent.prev", "nall.prev")

######################################################
####              combine everything              ####
######################################################

#--------------------------------------#
##--   merge with survival results  --##
#--------------------------------------#

## define intersection of phecodes across all data sets to be included
phe.shared   <- Reduce(intersect, list(res.surv$phecode, res.gen$phecode, res.adj$phecode, res.prev$phecode))
## n = 517 phecodes passing qc in all four data sets

## add those; restrict to outcomes in both; account for the fact that some proteins had no cis-pQTL
res.comb     <- merge(res.surv[ phecode %in% phe.shared & id %in% tolower(prot.scallop$Assay)], 
                      res.gen, by.x = c("phecode", "id"),
                      by.y = c("phecode", "id"),
                      all.x = T)

## indicator of significance: 9.371298e-08
res.comb[, sig.surv.crude := pval < .05/(nrow(res.comb))]

## same for genetics: 1.442148e-07
res.comb[, sig.mr := PVAL < .05/nrow(res.comb[ !is.na(PVAL)])]

## indicate high-conf genetic support
res.comb[, high.conf.gen := PP.H4.abf >= .8 & sig.mr == T]

## add survival results adjusted for other confounders
res.comb     <- merge(res.comb[, .(phecode, id, category, cl, phenotype, beta, se, pval, nevent, nall, p.resid.prot, sig.surv.crude, 
                                   sig.mr, high.conf.gen, BETA, SE, PVAL, PP.H0.abf, PP.H1.abf, PP.H3.abf, PP.H4.abf)],
                      res.adj[, .(phecode, id, beta, se, pval, nevent, nall, p.resid.prot)], 
                      by = c("phecode", "id"), 
                      suffixes = c(".surv.crude", ".surv.adj"))

## add results from prevalent phecodes
res.comb     <- merge(res.comb, res.prev, by = c("phecode", "id", "category", "phenotype", "cl"))

######################################################
####     pQTL enrichment of biomarker signature   ####
######################################################

#---------------------------#
##-- import pQTL results --##
#---------------------------#

## scallop results
res.scallop <- fread("<path to dir>/Results.SCALLOP.GWAS.catalog.overlap.pathway.genes.20250711.txt")
## pQTLs matched to the GWAS catalog
res.gwas    <- fread("<path to dir>/Results.R2.groupGWAS.catalog.pruned.overlap.pathway.genes.20250711.txt")

#---------------------------#
##--  perform enrichm.   --##
#---------------------------#

## --> w/o additional adjustment <-- ##

## pQTLs with ≥10 proteins
r2.test     <- res.gwas[ R2.count >= 10 ]$R2.group

## run in parallel
registerDoMC(10)

## run the actual enrichment
res.enrich.crude <- lapply(unique(res.comb$phecode), function(x){
  
  ## get all associated proteins
  prot.foreground <- res.comb[ phecode == x & sig.surv.crude == T]$id
  
  ## proceed only if enough findings
  if(length(prot.foreground) >= 5){
    
    ## iterate through pQTLs
    res <- mclapply(r2.test, function(k){
      
      ## get pQTL associated biomarker
      prot.pqtl <- tolower(res.scallop[ R2.group == k]$Protein)
      
      ## biomarker and pQTL associated
      d1    <- length(intersect(prot.pqtl, prot.foreground))
      ## not biomarker and pQTL associated
      d2    <- length(prot.pqtl[!(prot.pqtl %in% prot.foreground)])
      ## biomarker and not pQTL associated
      d3    <- length(prot.foreground[!(prot.foreground %in% prot.pqtl)])
      ## not biomarker and not pQTL
      d4    <- length(unique(prot.scallop[ !(tolower(Assay) %in% c(prot.pqtl, prot.foreground))]$Assay))
      
      ## test for enrichment
      enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
      
      ## return information needed
      return(data.table(phecode=x, R2.group = k, or = enr$estimate, pval = enr$p.value, 
                        intersection=paste(intersect(prot.pqtl, prot.foreground), collapse = "|"), 
                        d1=d1, d2=d2, d3=d3, d4=d4))
      
    }, mc.cores = 10)
    return(rbindlist(res))
  }
})
## combine all results
res.enrich.crude <- rbindlist(res.enrich.crude)

## add FDR correction
res.enrich.crude[, fdr := p.adjust(pval, method = "BH")]

## add labels
res.enrich.crude <- merge(res.enrich.crude, unique(res.comb[, .(phecode, phenotype, category, cl)]))
res.enrich.crude <- merge(res.enrich.crude, res.gwas[, .(R2.group, R2.count, protein.profile, cis.trans, candidate.gene.classifier, candidate.gene.score.top, mapped_trait, pqtl.category, drug.indicator, drug.gene, drug.id, drug.name, drug.gwas.evidence)],
                          by = "R2.group")

## --> w/ additional adjustment <-- ##

## run in parallel
registerDoMC(10)

## run the actual enrichment
res.enrich.adj   <- lapply(unique(res.comb$phecode), function(x){
  
  ## get all associated proteins
  prot.foreground <- res.comb[ phecode == x & sig.surv.adj == T]$id
  
  ## proceed only if enough findings
  if(length(prot.foreground) >= 5){
    
    ## iterate through pQTLs
    res <- mclapply(r2.test, function(k){
      
      ## get pQTL associated biomarker
      prot.pqtl <- tolower(res.scallop[ R2.group == k]$Protein)
      
      ## biomarker and pQTL associated
      d1    <- length(intersect(prot.pqtl, prot.foreground))
      ## not biomarker and pQTL associated
      d2    <- length(prot.pqtl[!(prot.pqtl %in% prot.foreground)])
      ## biomarker and not pQTL associated
      d3    <- length(prot.foreground[!(prot.foreground %in% prot.pqtl)])
      ## not biomarker and not pQTL
      d4    <- length(unique(prot.scallop[ !(tolower(Assay) %in% c(prot.pqtl, prot.foreground))]$Assay))
      
      ## test for enrichment
      enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
      
      ## return information needed
      return(data.table(phecode=x, R2.group = k, or = enr$estimate, pval = enr$p.value, 
                        intersection=paste(intersect(prot.pqtl, prot.foreground), collapse = "|"), 
                        d1=d1, d2=d2, d3=d3, d4=d4))
      
    }, mc.cores = 10)
    return(rbindlist(res))
  }
})
## combine all results
res.enrich.adj   <- rbindlist(res.enrich.adj)

## add FDR correction
res.enrich.adj[, fdr := p.adjust(pval, method = "BH")]

## add labels
res.enrich.adj <- merge(res.enrich.adj, unique(res.comb[, .(phecode, phenotype, category, cl)]))
res.enrich.adj <- merge(res.enrich.adj, res.gwas[, .(R2.group, R2.count, protein.profile, cis.trans, candidate.gene.classifier, candidate.gene.score.top, mapped_trait, pqtl.category, drug.indicator, drug.gene, drug.id, drug.name, drug.gwas.evidence)],
                        by = "R2.group")

#---------------------------#
##--     link to GWAS    --##
#---------------------------#

## import all entries and generate one large table
ot.diseases <- dir("<path to dir>")
ot.diseases <- grep("part", ot.diseases, value=T)
ot.diseases <- lapply(ot.diseases, function(x) read_parquet(paste0("<path to dir>", x)))
ot.diseases <- lapply(ot.diseases, as.data.table)
ot.diseases <- rbindlist(ot.diseases)
## transform some columns to make life easier
ot.diseases[, dbXRefs := sapply(dbXRefs, function(x) paste(x, collapse = "|"))]
ot.diseases[, directLocationIds := sapply(directLocationIds, function(x) paste(x, collapse = "|"))]
ot.diseases[, obsoleteTerms := sapply(obsoleteTerms, function(x) paste(x, collapse = "|"))]
ot.diseases[, parents := sapply(parents, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasBroadSynonym := sapply(synonyms.hasBroadSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasExactSynonym := sapply(synonyms.hasExactSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasNarrowSynonym := sapply(synonyms.hasNarrowSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, synonyms.hasRelatedSynonym := sapply(synonyms.hasRelatedSynonym, function(x) paste(x, collapse = "|"))]
ot.diseases[, descendants := sapply(descendants, function(x) paste(x, collapse = "|"))]
ot.diseases[, children := sapply(children, function(x) paste(x, collapse = "|"))]
ot.diseases[, therapeuticAreas := sapply(therapeuticAreas, function(x) paste(x, collapse = "|"))]
ot.diseases[, indirectLocationIds := sapply(indirectLocationIds, function(x) paste(x, collapse = "|"))]
ot.diseases[, ancestors := sapply(ancestors, function(x) paste(x, collapse = "|"))]

## harmonize codings
ot.diseases[, dbXRefs := gsub("\\.", "", dbXRefs)]

#-----------------------------#
##-- map OT ID to phecodes --##
#-----------------------------#

## import phecode mapping
lab.phecodes <- fread("<path to dir>")

## create an unpacked set of ot.diseases (inflate)
tmp.ot       <- rbindlist(mclapply(1:nrow(ot.diseases), function(x) return(data.table(id=ot.diseases$id[x], db.ref = strsplit(ot.diseases$dbXRefs[x], "\\|")[[1]]))))

## add to phecode label
lab.phe[, ot.disease.mapping := sapply(phecode, function(x){
  
  ## get relevant code
  jj      <- which(lab.phecodes$id == paste0("date_", x))
  
  ## do only if any mapping
  if(length(jj) > 0){
    ## get all relevant codes
    icd10   <- sapply(strsplit(lab.phecodes$icd10[jj], "\\|")[[1]], function(c) paste0("ICD10:", c))
    snomed1 <- sapply(strsplit(lab.phecodes$snomed[jj], "\\|")[[1]], function(c) paste0("SNOMEDCT:", c))
    snomed2 <- sapply(strsplit(lab.phecodes$snomed[jj], "\\|")[[1]], function(c) paste0("SCTID:", c))
    
    ## return relevant mappings
    return(paste(unique(tmp.ot[ db.ref %in% c(icd10, snomed1, snomed2)]$id), collapse = "|"))
  }else{
    return("")
  }
  
})]

## import phecode to EFO assignment
phecode.efo <- fread("<path to dir>")

## map EFO terms more broadly to phecodes
icd10.efo   <- read.delim("<path to dir>", sep= "\t", fill = T)
icd10.efo   <- as.data.table(icd10.efo)

## add to label
lab.phe[, efo.mapping := sapply(phecode, function(x){
  
  ## test if looking at a phecode
  if(paste0("date_", x) %in% lab.phecodes$id){
    ## get the relevant icd-10 codes
    icd10 <- strsplit(lab.phecodes[ which(id == paste0("date_", x))]$icd10, "\\|")
    ## shrink to three characters
    icd10 <- sapply(icd10, function(c) substr(c, 1, 3))
    print(icd10)
    ## return all matching EFO terms
    return(paste(icd10.efo[ICD10_CODE.SELF_REPORTED_TRAIT_FIELD_CODE %in% icd10]$MAPPED_TERM_URI, collapse = "|"))
  }else{
    return("")
  }
  
})]

## --> add indicator, whether the pQTL has been linked to disease <-- ##

## overall enrichment
res.enrich.crude[, gwas.evidence.pQTL := apply(res.enrich.crude[, .(R2.group, phecode)], 1, function(x){
  
  ## not done for MHC region!
  if(x[1] != 0){
    ## identify relevant EFO terms (if any)
    efo.gwas <- strsplit(res.gwas[ R2.group == x[1]]$efo.term, "\\|\\||\\|")[[1]]
    efo.gwas <- sapply(efo.gwas, function(k) gsub(":", "_", k))
    ## get efo for phecode
    efo.phe  <- c(strsplit(lab.phe[ phecode == x[2]]$ot.disease.mapping, "\\||, ")[[1]],
                  strsplit(lab.phe[ phecode == x[2]]$efo.mapping, "\\||, ")[[1]])
    ## return results
    return(paste( intersect(efo.gwas, efo.phe), collapse = "|"))
  }else{
    return("")
  }
  
})]

## adjusted enrichment
res.enrich.adj[, gwas.evidence.pQTL := apply(res.enrich.adj[, .(R2.group, phecode)], 1, function(x){
  
  ## not done for MHC region!
  if(x[1] != 0){
    ## identify relevant EFO terms (if any)
    efo.gwas <- strsplit(res.gwas[ R2.group == x[1]]$efo.term, "\\|\\||\\|")[[1]]
    efo.gwas <- sapply(efo.gwas, function(k) gsub(":", "_", k))
    ## get efo for phecode
    efo.phe  <- c(strsplit(lab.phe[ phecode == x[2]]$ot.disease.mapping, "\\||, ")[[1]],
                  strsplit(lab.phe[ phecode == x[2]]$efo.mapping, "\\||, ")[[1]])
    ## return results
    return(paste( intersect(efo.gwas, efo.phe), collapse = "|"))
  }else{
    return("")
  }
  
})]

######################################################
####              new summary figure              ####
######################################################

## --> create variables to ease plotting <-- ##

## indicator support: crude survival models --> genetics
res.comb[, indicator.surv.crude.gen := ifelse(sig.surv.crude == T & sign(beta.surv.crude) == sign(BETA) & PVAL < .05 & PP.H4.abf > .8, "concordant",
                                              ifelse(sig.surv.crude == T & sign(beta.surv.crude) != sign(BETA) & PVAL < .05 & PP.H4.abf > .8, "discordant", "no support"))]
## account for possible missing values
res.comb[, indicator.surv.crude.gen := ifelse(is.na(indicator.surv.crude.gen), "no support", indicator.surv.crude.gen)]
## indicator support: genetics --> crude survival models
res.comb[, indicator.gen.surv.crude := ifelse(high.conf.gen == T & sign(beta.surv.crude) == sign(BETA) & pval.surv.crude < .05, "concordant",
                                              ifelse(high.conf.gen == T & sign(beta.surv.crude) != sign(BETA) & pval.surv.crude < .05, "discordant", "no support"))]
## account for possible missing values
res.comb[, indicator.gen.surv.crude := ifelse(is.na(indicator.gen.surv.crude), "no support", indicator.gen.surv.crude)]


## indicator support: adjust survival models --> genetics
res.comb[, indicator.surv.adj.gen := ifelse(sig.surv.adj == T & sign(beta.surv.adj) == sign(BETA) & PVAL < .05 & PP.H4.abf > .8, "concordant",
                                            ifelse(sig.surv.adj == T & sign(beta.surv.adj) != sign(BETA) & PVAL < .05 & PP.H4.abf > .8, "discordant", "no support"))]
## account for possible missing values
res.comb[, indicator.surv.adj.gen := ifelse(is.na(indicator.surv.adj.gen), "no support", indicator.surv.adj.gen)]
## indicator support: genetics --> adjusted survival models
res.comb[, indicator.gen.surv.adj := ifelse(high.conf.gen == T & sign(beta.surv.adj) == sign(BETA) & pval.surv.adj < .05, "concordant",
                                            ifelse(high.conf.gen == T & sign(beta.surv.adj) != sign(BETA) & pval.surv.adj< .05, "discordant", "no support"))]
## account for possible missing values
res.comb[, indicator.gen.surv.adj := ifelse(is.na(indicator.gen.surv.adj), "no support", indicator.gen.surv.adj)]


## indicator support: prevalent models --> genetics
res.comb[, indicator.prev.gen := ifelse(sig.prev == T & sign(beta.prev) == sign(BETA) & PVAL < .05 & PP.H4.abf > .8, "concordant",
                                        ifelse(sig.prev == T & sign(beta.prev) != sign(BETA) & PVAL < .05 & PP.H4.abf > .8, "discordant", "no support"))]
## account for possible missing values
res.comb[, indicator.prev.gen := ifelse(is.na(indicator.prev.gen), "no support", indicator.prev.gen)]
## indicator support: genetics --> adjusted survival models
res.comb[, indicator.gen.prev := ifelse(high.conf.gen == T & sign(beta.prev) == sign(BETA) & pval.prev < .05, "concordant",
                                        ifelse(high.conf.gen == T & sign(beta.prev) != sign(BETA) & pval.prev < .05, "discordant", "no support"))]
## account for possible missing values
res.comb[, indicator.gen.prev := ifelse(is.na(indicator.gen.prev), "no support", indicator.gen.prev)]


## --> define any convergence (survival) <-- ##

## apply
res.comb[, convergent.example := indicator.surv.crude.gen == "concordant" | indicator.gen.surv.crude == "concordant"]

#---------------------------------------#
##--  summary figure - cis genetics  --##
#---------------------------------------#

## open device
pdf("../graphics/Summary.Convergence.MR.obs.biomarker.20251204.pdf", width = 6.3, height = 6.3*(2/3))

## graphical parameters
par(mar=c(1.5,1.5,.5,.5), cex.axis = .5, cex.lab = .5, tck = -.01, mgp = c(.6,0,0), lwd = .5, bty = "n", xaxs = "i", yaxs = "i")

## define layout
layout(matrix(1:6,2,3, byrow = T), heights = c(.2,.8), widths = c(.46,.27,.27))

#---------------------------------#
##--    fraction of support    --##
#---------------------------------#

## space for names
par(mar=c(1.5,10,.5,.5))

## --> crude survival models <-- ##

## empty plot (two rectangles are needed)
plot(c(0, 1), c(0,2), type = "n", xlab="Percentage (Incident models)", ylab="", xaxt="n", yaxt="n")
## add axis
axis(1, at=seq(0,1,.25), labels = seq(0,1,.25)*100, lwd=.5)
## get plotting coordinates
pm <- par("usr")

## add first rectangle: genetics --> survial models (Living Coral (#FC766AFF), Storm Gray (#B0B8B4FF) and Forest Biome (#184A45FF))
tmp <- table(res.comb[ high.conf.gen == T]$indicator.gen.surv.crude)
rect(c(0, cumsum(tmp)[1:2])/sum(tmp), 1.1, cumsum(tmp)/sum(tmp), 1.6, lwd=.3, border = "white", col = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))
## add legend
# text(pm[1], 1.9, labels=paste(paste0(names(tmp), " (n=",tmp, ")"), collapse = " - "), cex=.4, pos=4, xpd = NA)
legend(pm[1], 2.1, bty = "n", cex=.4, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 2,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))

## add second rectangle: survial models --> genetics (Living Coral (#FC766AFF), Storm Gray (#B0B8B4FF) and Forest Biome (#184A45FF))
tmp <- table(res.comb[ sig.surv.crude == T]$indicator.surv.crude.gen)
rect(c(0, cumsum(tmp)[1:2])/sum(tmp), .1, cumsum(tmp)/sum(tmp), .6, lwd=.3, border = "white", col = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))
## add legend
legend(pm[1], 1.1, bty = "n", cex=.4, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 2,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))

## add text
text(pm[1], c(.45, 1.45), cex=.6, xpd = NA, pos=2, offset = .1,
     labels = c("Observational -> Genetics", "Genetics -> Observational"))

## --> adjusted survival models <-- ##

## reduce space for names
par(mar=c(1.5,.5,.5,.5))

## empty plot (two rectangles are needed)
plot(c(0, 1), c(0,2), type = "n", xlab="Percentage (Extended incident models)", ylab="", xaxt="n", yaxt="n")
## add axis
axis(1, at=seq(0,1,.25), labels = seq(0,1,.25)*100, lwd=.5)
## get plotting coordinates
pm <- par("usr")

## add first rectangle: genetics --> survial models (Living Coral (#FC766AFF), Storm Gray (#B0B8B4FF) and Forest Biome (#184A45FF))
tmp <- table(res.comb[ high.conf.gen == T]$indicator.gen.surv.adj)
rect(c(0, cumsum(tmp)[1:2])/sum(tmp), 1.1, cumsum(tmp)/sum(tmp), 1.6, lwd=.3, border = "white", col = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))
## add legend
legend(pm[1], 2.1, bty = "n", cex=.4, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 2,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))

## add second rectangle: survial models --> genetics (Living Coral (#FC766AFF), Storm Gray (#B0B8B4FF) and Forest Biome (#184A45FF))
tmp <- table(res.comb[ sig.surv.adj == T]$indicator.surv.adj.gen)
rect(c(0, cumsum(tmp)[1:2])/sum(tmp), .1, cumsum(tmp)/sum(tmp), .6, lwd=.3, border = "white", col = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))
legend(pm[1], 1.1, bty = "n", cex=.4, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 2,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))

## --> prevalent models <-- ##

## empty plot (two rectangles are needed)
plot(c(0, 1), c(0,2), type = "n", xlab="Percentage (Prevalent models)", ylab="", xaxt="n", yaxt="n")
## add axis
axis(1, at=seq(0,1,.25), labels = seq(0,1,.25)*100, lwd=.5)
## get plotting coordinates
pm <- par("usr")

## add first rectangle: genetics --> survial models (Living Coral (#FC766AFF), Storm Gray (#B0B8B4FF) and Forest Biome (#184A45FF))
tmp <- table(res.comb[ high.conf.gen == T]$indicator.gen.prev)
rect(c(0, cumsum(tmp)[1:2])/sum(tmp), 1.1, cumsum(tmp)/sum(tmp), 1.6, lwd=.3, border = "white", col = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))
## add legend
legend(pm[1], 2.1, bty = "n", cex=.4, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 2,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))

## add second rectangle: survial models --> genetics (Living Coral (#FC766AFF), Storm Gray (#B0B8B4FF) and Forest Biome (#184A45FF))
tmp <- table(res.comb[ sig.prev == T]$indicator.prev.gen)
rect(c(0, cumsum(tmp)[1:2])/sum(tmp), .1, cumsum(tmp)/sum(tmp), .6, lwd=.3, border = "white", col = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))
## add legend
legend(pm[1], 1.1, bty = "n", cex=.4, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 2,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#FC766AFF", adjustcolor("#FC766AFF", .5), "#B0B8B4FF"))

#---------------------------------#
##--    converging examples    --##
#---------------------------------#

## select examples with survival and prevalent support 
tmp <- res.comb[(indicator.surv.crude.gen == "concordant" & indicator.prev.gen == "concordant") | 
                  (indicator.gen.surv.crude == "concordant" & indicator.gen.prev == "concordant")]

## order accrodingly for plotting
tmp <- tmp[ order(id, phecode)]
tmp[, srt := 1:.N]

## --> genetic prediction <-- ##

## space for names
par(mar=c(1.5,10,1.5,.5))

## empty plot
plot(c(-1.1,.8), c(.5, nrow(tmp)+.5), ylim = rev(c(.5, nrow(tmp)+.5)), type = "n", xlab = "Effect per 1 s.d. increase in genetically predicted protein levels",
     ylab = "", xaxt = "n", yaxt = "n")
## add x-axis
axis(1, lwd=.5)
## get borders
pm <- par("usr")
## rectangles to divide (by protein)
rect(pm[1], tmp[, .(left = min(srt)), by = "id"]$left - .5, 
     pm[2], tmp[, .(right = max(srt)), by = "id"]$right + .5 , col = c("white", "grey90"),
     border = NA)
## zero crossing
abline(v=0, lwd=.5)

## add confidence intervals
arrows(tmp$BETA - 1.96 * tmp$SE, 1:nrow(tmp), tmp$BETA + 1.96 * tmp$SE, 1:nrow(tmp), length = 0,
       lwd=.3)
## add point estimates
points(tmp$BETA, 1:nrow(tmp), pch = 22, lwd=.3, col = "grey20", cex = .5, bg = "#FC766AFF")
## add names
text(pm[1], 1:nrow(tmp), cex=.45, xpd=NA, pos = 2, offset = .1, labels = paste0(toupper(tmp$id), " - ", tmp$phenotype))
## add title
mtext("Mendelian randomization", cex=.5)


## --> Survial models <-- ##

## reduce space for names
par(mar=c(1.5,.5,1.5,.5))

## empty plot
plot(c(-.7,1.1), c(.5, nrow(tmp)+.5), ylim = rev(c(.5, nrow(tmp)+.5)), type = "n", xlab = "Effect per 1 s.d. increase in measured protein levels",
     ylab = "", xaxt = "n", yaxt = "n")
## add x-axis
axis(1, lwd=.5)
## get borders
pm <- par("usr")
## rectangles to divide (by protein)
rect(pm[1], tmp[, .(left = min(srt)), by = "id"]$left - .5, 
     pm[2], tmp[, .(right = max(srt)), by = "id"]$right + .5 , col = c("white", "grey90"),
     border = NA)
## zero crossing
abline(v=0, lwd=.5)

## add confidence intervals: minimal adjustment
arrows(tmp$beta.surv.crude - 1.96 * tmp$se.surv.crude, 1:nrow(tmp)-.25, tmp$beta.surv.crude + 1.96 * tmp$se.surv.crude, 1:nrow(tmp)-.25, length = 0,
       lwd=.3)
## add confidence intervals: extended adjustment
arrows(tmp$beta.surv.adj - 1.96 * tmp$se.surv.adj, 1:nrow(tmp)+.25, tmp$beta.surv.adj + 1.96 * tmp$se.surv.adj, 1:nrow(tmp)+.25, length = 0,
       lwd=.3)
## add point estimates: minimal adjustment
points(tmp$beta.surv.crude, 1:nrow(tmp)-.25, pch = 22, lwd=.3, col = "grey20", cex = .5,
       bg = ifelse(tmp$indicator.gen.surv.crude == "concordant" | tmp$indicator.surv.crude.gen == "concordant", "#FC766AFF",
                   ifelse(tmp$indicator.gen.surv.crude == "disconcordant" & tmp$indicator.surv.crude.gen == "disconcordant", adjustcolor("#FC766AFF", .5), "#B0B8B4FF")))
## add point estimates: maximal adjustment
points(tmp$beta.surv.adj, 1:nrow(tmp)+.25, pch = 21, lwd=.3, col = "grey20", cex = .5,
       bg = ifelse(tmp$indicator.gen.surv.adj == "concordant" | tmp$indicator.surv.adj.gen == "concordant", "#FC766AFF",
                   ifelse(tmp$indicator.gen.surv.adj == "disconcordant" & tmp$indicator.surv.adj.gen == "disconcordant", adjustcolor("#FC766AFF", .5), "#B0B8B4FF")))
## add title
mtext("Incident disease risk", cex=.5)
## specific legend
legend("topright", cex=.5, bty="n", lty=0,
       pch=c(22,21), pt.lwd=.3, pt.cex = .8,
       legend = c("age + sex", ".. + confounder"))


## --> genetic prediction <-- ##

## empty plot
plot(c(-1.1,1.3), c(.5, nrow(tmp)+.5), ylim = rev(c(.5, nrow(tmp)+.5)), type = "n", xlab = "Effect per 1 s.d. increase in protein levels",
     ylab = "", xaxt = "n", yaxt = "n")
## add x-axis
axis(1, lwd=.5)
## get borders
pm <- par("usr")
## rectangles to divide (by protein)
rect(pm[1], tmp[, .(left = min(srt)), by = "id"]$left - .5, 
     pm[2], tmp[, .(right = max(srt)), by = "id"]$right + .5 , col = c("white", "grey90"),
     border = NA)
## zero crossing
abline(v=0, lwd=.5)

## add confidence intervals
arrows(tmp$beta.prev - 1.96 * tmp$se.prev, 1:nrow(tmp), tmp$beta.prev + 1.96 * tmp$se.prev, 1:nrow(tmp), length = 0,
       lwd=.3)
## add point estimates
points(tmp$beta.prev, 1:nrow(tmp), pch = 22, lwd=.3, col = "grey20", cex = .5,
       bg = ifelse(tmp$indicator.gen.prev == "concordant" | tmp$indicator.prev.gen == "concordant", "#FC766AFF",
                   ifelse(tmp$indicator.gen.prev == "disconcordant" & tmp$indicator.prev.gen == "disconcordant", adjustcolor("#FC766AFF", .5), "#B0B8B4FF")))
## add title
mtext("Prevalent disease risk", cex=.5)

## close device
dev.off()


#---------------------------------------#
##--  summary figure - trans enrich. --##
#---------------------------------------#

## open device
pdf("../graphics/pQTL.enrichment.obs.biomarker.dot.plot.20250715.pdf", width = 6.3, height = 6.3*(1/3))

## graphical parameters
par(mar=c(5,1.5,1.5,.5), cex.axis = .5, cex.lab = .5, tck = -.01, mgp = c(.6,0,0), lwd = .5, bty = "n", xaxs = "i", yaxs = "i")

## get what should be drawn
tmp <- res.enrich.crude[or > 1 & fdr < .05 & gwas.evidence.pQTL != ""]
## cap enrichment
tmp[, or.plot := ifelse(is.finite(or) & or < 100, or, 100)]
## order by enrichment
tmp <- tmp[ order(candidate.gene.classifier, fdr)]
## add sorting
tmp[, srt := 1:.N]

## create color vector for plotting?

## empty plot
plot(c(.5, nrow(tmp)+.5), c(0, 17), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = expression(-log[10]("FDR")))
## add axis
axis(2, lwd=.5)

## devide by protein target
foo <- tmp[, .(left = min(srt), right = max(srt), mid = mean(srt)), by = c("candidate.gene.classifier", "drug.id")]
rect(foo$left-.5, 0, foo$right+.5, 17, border = NA, col = c("white", "grey90"))

## add estimates
points(1:nrow(tmp), -log10(tmp$fdr), pch = 21, cex=log10(tmp$or.plot), lwd=.3,
       bg = adjustcolor(ifelse(tmp$cis.trans == "cis", "#F95700", "#00A4CC"), .6), xpd = NA)

## add legend
legend("topright", bty = "n", cex = .5, lty = 0, pch = 22, 
       pt.lwd = .3, pt.bg = c("#F95700", "#00A4CC"),
       legend = c("cis-pQTL", "trans-pQTL"))

## get plotting coordinates
pm <- par("usr")

## add some text
text(1:nrow(tmp), pm[3] - (pm[4]-pm[3])*.05, pos = 2, xpd = NA, srt = 60, offset = 0, 
     labels = tmp$phenotype, cex=.35)

## add label
text(foo$mid, pm[4] + (pm[4]-pm[3])*c(0,.08,.16,.24), pos = 3, xpd = NA, offset = 0, 
     labels = foo$candidate.gene.classifier,
     cex=.35, font = ifelse(is.na(foo$drug.id), 3, 4))

## legend for Odds ratio
legend(pm[1]+(pm[2]-pm[1])*.05, pm[4], bty = "n", lty=0, cex = .5, pch = 21,
       pt.cex = log10(c(5,10,50,100)),
       title = "Fold enrichment", legend = c(5,10,50,100),
       ncol=4, pt.bg = rgb(0,0,0,.2))


## close device
dev.off()

######################################################
####          look into discrepant examples       ####
######################################################

## create relevant columns
prot.scallop[, phecode.analysis := tolower(Assay) %in% res.comb$id]

## how many cis/trans pQTL
prot.scallop <- merge(prot.scallop, res.scallop[, .(num.cis.pQTL = sum(cis_trans == "Cis"),
                                                    num.trans.pQTL = sum(cis_trans != "Cis")), by = "MA_prot_id"],
                      by.x = "OID_MA", by.y = "MA_prot_id")
## define categories
prot.scallop[, pQTL.category := ifelse(num.cis.pQTL > 0 & num.trans.pQTL == 0, "cis",
                                       ifelse(num.cis.pQTL > 0 & num.trans.pQTL > 0, "cis+trans", "trans"))]

#------------------------------------------------------#
##-- check for genetic connection of top rank assoc --##
#------------------------------------------------------#

## create tmp file to ease mapping
tmp                <- rbindlist(lapply(1:nrow(lab.phe), function(x) return(data.table(phecode = lab.phe$phecode[x],
                                                                                      efo.term = unique(c( strsplit(lab.phe$ot.disease.mapping[x], "\\||, ")[[1]],
                                                                                                           strsplit(lab.phe$efo.mapping[x], "\\||, ")[[1]]))))))

## create matching phecode column to GWAS cat results
res.gwas[, phecode.efo := sapply(efo.term, function(x){
  
  ## get all phecodes
  efo.gwas <- strsplit(x, "\\|\\||\\|")[[1]]
  efo.gwas <- sapply(efo.gwas, function(k) gsub(":", "_", k))
  
  ## return matching phecodes
  return(paste(sort(unique(tmp[ efo.term %in% efo.gwas]$phecode)), collapse = "|"))
  
})]

## add indicator, whether the protein could be linked genetically to the disease (based on GWAS cat.)
tmp                <- rbindlist(lapply(1:nrow(res.gwas), function(x){
  ## extract entries
  prot <- tolower(strsplit(res.gwas$protein.profile[x], "\\|")[[1]])
  phe  <- strsplit(res.gwas$phecode.efo[x], "\\|")[[1]]
  ## return (if any)
  if(length(phe) > 0){
    return(data.table(R2.group = res.gwas$R2.group[x], expand.grid(id = prot, phecode = phe, stringsAsFactors = F)))
  }
}))

## drop MHC region
tmp                <- tmp[ R2.group != 0]

## compress R2 groups to ease merging
tmp                <- tmp[, .(R2.group = paste(unique(R2.group), collapse = "|"),
                              R2.num = length(R2.group)),
                          by = c("id", "phecode")]
## add cis/trans information
tmp[, cis.trans := apply(tmp[, .(id, R2.group)], 1, function(x){
  
  ## split R2 groups
  r2 <- strsplit(x[2], "\\|")[[1]]
  ## return info
  return(paste(sort(unique(tolower(res.scallop[ R2.group %in% r2 & Protein == toupper(x[1])]$cis_trans))), collapse = "|"))
  
})]
## change data type to enable merging
tmp[, phecode := as.numeric(phecode)]

## add to results
res.comb           <- merge(res.comb, tmp, by = c("id", "phecode"), all.x = T)

## --> more likely to persist adjustment? <-- ##

## table
tmp <- res.comb[ sig.surv.crude == T]
table(tmp$sig.surv.adj, ifelse(tmp$R2.num > 0 & !is.na(tmp$R2.num), 1, 0))
## more likely 
fisher.test(table(tmp$sig.surv.adj, ifelse(tmp$R2.num > 0 & !is.na(tmp$R2.num), 1, 0)))$p.value
## odds ratio: 1.677; p-value: 1.721998e-60
  
