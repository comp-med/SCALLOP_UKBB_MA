#################################################
#### Assignment of trans-pQTLs               ####
#### Maik Pietzner                20/05/2024 ####
#################################################

rm(list=ls())
setwd("<path to dir>")
options(stringsAsFactors=F)
load(".RData")

## --> packages required <-- ##

require(data.table)
require(doMC)
require(readxl)
require(caret)
require(igraph)
require(pROC)
require(tidyverse)
require(gprofiler2)
require(ndjson)
require(circlize)
require(arrow)
require(pscl)

##############################################
####      import protein statistics       ####
##############################################

## import all results (no cleaning as for MR)
res.scallop <- fread("<path to dir>/Finemapped_results_SCALLOP_UKBB_R2_results.txt")

## import LD-proxies
r2.proxies  <- fread("<path to dir>/SCALLOP_UKBB_MA_all_finemapped_r0.6.formatted_withMarkerName.ld")
r2.proxies  <- unique(r2.proxies)

###############################################
####          map to build 38              ####
###############################################

#-------------------------------------#
##--      lift over to build 38    --##
#-------------------------------------#

## R package to do so
require(liftOver)

## edit chromosome
r2.proxies[, CHR_B.hg39 := ifelse( CHR_B == 23, "X", CHR_B)]

## import chain
chain              <- import.chain("<path to dir>/hg19ToHg38.over.chain")

## create GRanges object
grObject           <- GRanges(seqnames = paste0("chr", r2.proxies$CHR_B.hg39), ranges=IRanges(start = r2.proxies$BP_B, end = r2.proxies$BP_B, names=r2.proxies$SNP_B))

## now do the mapping
tmp                <- as.data.table(liftOver(grObject, chain))
## rename
names(tmp)         <- c("group", "SNP_B", "seqnames", "BP_B.hg38", "end", "width", "strand")
## add to the data
r2.proxies         <- merge(r2.proxies, unique(tmp[, c("SNP_B", "BP_B.hg38")]), all.x=T, by="SNP_B")


##############################################
####           import GWAS catalog        ####
##############################################

## import download from 20/05/2023
gwas.catalogue           <- fread("gwas_catalog_v1.0.2-associations_e111_r2024-05-05.tsv", sep="\t", header=T)

## rename some columns
names(gwas.catalogue)[8] <- "TRAIT"

## prune GWAS catalogue data
gwas.catalogue           <- gwas.catalogue[ !is.na(`OR or BETA`) & is.finite(`OR or BETA`) ]
## generate risk allele and drop everything w/o this information
gwas.catalogue[, riskA := sapply(`STRONGEST SNP-RISK ALLELE`, function(x) strsplit(x,"-")[[1]][2])] 
gwas.catalogue[, riskA := trimws(riskA, which = "b")] 
## drop interaction entries
ii                       <- grep("[0-9]", gwas.catalogue$riskA)
gwas.catalogue           <- gwas.catalogue[-ii,]
## only genome-wide significant ones
gwas.catalogue           <- gwas.catalogue[ PVALUE_MLOG > 7.3 ]

## create another entry to possible merge on (careful, genome build 38 mapping)
gwas.catalogue[, snp.id := paste0(ifelse(CHR_ID == "X", 23, CHR_ID), ":", CHR_POS)]

###############################################
####  assign each variant to GWAS catalog  ####
###############################################

#-------------------------------------#
##--        create mapping         --##
#-------------------------------------#

## get all unique variants from regional results
res.var     <- unique(res.scallop[, .(MarkerName, ID, CHROM, position_hg19, Allele1_SCALLOP_MA, Allele2_UKBB_SCALLOP_MA)])

## do in parallel
require(doMC)
registerDoMC(12)

## go through each variant and 1) identify all proxies, 2) map to GWAS catalog findings, and 3) reduce redundancy
res.gwas    <- mclapply(1:nrow(res.var), function(x){
  
  ## get all possible proxies
  snp              <- r2.proxies[ MarkerName_SNP_A == res.var$MarkerName[x]]
  
  ## create snp id to optimize merging multiple mappings
  snp$snp.id       <- paste0(snp$CHR_A, ":", snp$BP_B.hg38)
  
  ## create two versions of mapping
  snp.rsid         <- merge(snp, gwas.catalogue, by.x="SNP_B", by.y="SNPS")
  snp.pos          <- merge(snp, gwas.catalogue, by = "snp.id")
  ## edit
  snp.rsid$snp.id   <- snp.rsid$snp.id.x
  snp.rsid$snp.id.x <- snp.rsid$snp.id.y <- NULL
  snp.pos$SNPS     <- snp.pos$SNP_B
  snp.rsid$SNPS    <- snp.rsid$SNP_B
  ## combine
  snp              <- unique(rbind(snp.rsid, snp.pos))
  
  ## prepare return
  if(nrow(snp) > 0){
    ## sort 
    snp              <- snp[order(TRAIT, SNP_A, -R2)]
    ## create indicator
    snp[, ind := 1:.N, by=c("TRAIT", "SNP_A")]
    ## keep only one finding per trait
    snp              <- snp[ind == 1]
    
    ## report summary back
    snp              <- data.table(rsid.gwas=paste(sort(unique(snp$SNP_B)), collapse = "||"),
                                   trait_reported=paste(sort(unique(snp$TRAIT)), collapse = "||"),
                                   mapped_trait=paste(sort(unique(snp$MAPPED_TRAIT)), collapse = "||"),
                                   study_id=paste(sort(unique(snp$`STUDY ACCESSION`)), collapse = "||"),
                                   source_gwas=paste(sort(unique(snp$PUBMEDID)), collapse = "||"),
                                   num_reported=nrow(snp))
    
  }else{
    
    ## report summary back
    snp              <- data.table(rsid.gwas="",
                                   trait_reported="",
                                   mapped_trait="",
                                   study_id="",
                                   source_gwas="",
                                   num_reported=0)
  }
  
  ## return data set
  return(data.table(res.var[x,], snp))
}, mc.cores = 10)
## combine 
res.gwas    <- rbindlist(res.gwas)
save.image()

## add to fine-mapping results
res.scallop <- merge(res.scallop, res.gwas)

###############################################
####    import candidate gene assignment   ####
###############################################

# ## import trans-pQTL candidate gene assignment
# trans.gene  <- data.table(read_excel("<path to dir>/trans_pqtl_causal_gene_assignment_scores_with_MHC_20240510.xlsx")) 
# ## reduce to unique entries
# tmp         <- unique(trans.gene[, .(pqtl_MarkerName, Protein_name, pheno)])
# trans.gene  <- lapply(1:nrow(tmp), function(x){
#   ## all entries needed
#   foo <- trans.gene[ pqtl_MarkerName == tmp$pqtl_MarkerName[x] & Protein_name == tmp$Protein_name[x] & pheno == tmp$pheno[x]]
#   print(foo)
#   ## create a single line return
#   return(data.table(tmp[x,], candidate_gene_name=paste(foo$candidate_gene_name, collapse = "|"),
#                     candidate_gene_ensembl_id=paste(foo$candidate_gene_ensembl_id, collapse = "|"),
#                     test=foo$test[1], score=foo$score[1]))
# })
# trans.gene  <- rbindlist(trans.gene)
# ## add to SCALLOP
# res.scallop <- merge(res.scallop, trans.gene, by.x = c("MarkerName", "Protein", "MA_prot_id"), by.y = c("pqtl_MarkerName", "Protein_name", "pheno"), all.x=T)
# ## import cis-pQTL assignment
# res.scallop[, candidate_gene_name := ifelse(is.na(candidate_gene_name), Protein, candidate_gene_name)]

###############################################
####        import pathway members         ####
###############################################

## compile all Olink genes with at least one genetic signal
olink.pathway <- unique(res.scallop[, .(Protein, Uniprot_ID)])
## n = 995

## split odd names (e.g., complexes)
olink.pathway <- lapply(1:nrow(olink.pathway), function(x){
  return(data.table(Protein=strsplit(olink.pathway$Protein[x], "_")[[1]],
                    Uniprot_ID=strsplit(olink.pathway$Uniprot_ID[x], "_")[[1]]))
})
## combine again
olink.pathway <- rbindlist(olink.pathway)

#---------------------------------#
##--     protein complexes     --##
#---------------------------------#

## import complexes
protein.complexes <- import_omnipath_complexes(resources = c("CORUM", "hu.MAP", "hu.MAP2"))
## spread out to ease mapping with olink data
protein.complexes <- lapply(1:nrow(protein.complexes), function(x){
  ## use only gene symbols to map
  return(data.table(protein.complexes[x,], gene.symbol=strsplit(protein.complexes$components_genesymbols[x], "_")[[1]]))
})
## combine back into a list
protein.complexes <- rbindlist(protein.complexes)

## add to SCALLOP
olink.pathway[, protein.complex.genes := sapply(Protein, function(x){
  ## get all unique complexes in which the gene is involved in
  tmp <- unique(protein.complexes[ gene.symbol == x]$components_genesymbols)
  ## return all unique genes (if any)
  if(length(tmp) > 0){
    return(paste(unique(unlist(lapply(tmp, function(k) strsplit(k, "_")[[1]]))), collapse = "|"))
  }else{
    return(NA)
  }
})]

#---------------------------------#
##--  receptor - ligand pairs  --##
#---------------------------------#

## import ligand/receptor interactions (with reference)
ligand.receptor <- OmnipathR::import_ligrecextra_interactions()
ligand.receptor <- as.data.table(ligand.receptor)
## subset to those with references
ligand.receptor <- ligand.receptor[ n_references > 0]

## add to olink
olink.pathway[, ligand.receptor.pair := sapply(Protein, function(x){
  ## get all unique complexes in which the gene is involved in
  tmp <- subset(ligand.receptor, source_genesymbol == x | target_genesymbol == x,)
  ## return all unique genes (if any)
  if(nrow(tmp) > 0){
    tmp <- unique(c(tmp$source_genesymbol, tmp$target_genesymbol))
    tmp <- unlist(lapply(tmp, function(k) strsplit(k, "_")[[1]]))
    return(paste(tmp, collapse = "|"))
  }else{
    return(NA)
  }
})]

#---------------------------------#
##--        KEGG pathways      --##
#---------------------------------#

## We get entrez ids and their pathways.
kegg.pathway        <- getGeneKEGGLinks(species="hsa")
## This is to get the gene symbols using entrez ids
kegg.pathway$Symbol <- mapIds(org.Hs.eg.db, kegg.pathway$GeneID, column="SYMBOL", keytype="ENTREZID")
# pathway names
pathway_names       <- getKEGGPathwayNames(species="hsa")
## obtain all human KEGG pathways
kegg.pathway        <-  merge(kegg.pathway, pathway_names, by="PathwayID")
kegg.pathway        <-  as.data.table(kegg.pathway)
## drop some very large, but unspecific pathways
kegg.pathway        <- kegg.pathway[ !(Description %in% c("Metabolic pathways - Homo sapiens (human)", "Pathways in cancer - Homo sapiens (human)",
                                                          "Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)", "MicroRNAs in cancer - Homo sapiens (human)",
                                                          "Coronavirus disease - COVID-19 - Homo sapiens (human)", "Chemical carcinogenesis - reactive oxygen species - Homo sapiens (human)",
                                                          "Proteoglycans in cancer - Homo sapiens (human)", "Chemical carcinogenesis - receptor activation - Homo sapiens (human)",
                                                          "Transcriptional misregulation in cancer - Homo sapiens (human)"))]

## add to olink
olink.pathway[, kegg.pathway.genes := sapply(Protein, function(x){
  ## get pathways the gene is involved in
  tmp <- subset(kegg.pathway, Symbol == x)$PathwayID
  print(tmp)
  ## return all unique genes (if any)
  if(length(tmp) > 0){
    return(paste(unique(kegg.pathway[ PathwayID %in% tmp]$Symbol), collapse = "|"))
  }else{
    return(NA)
  }
})]

#---------------------------------#
##--     Reactome pathways     --##
#---------------------------------#

## get pathways
reactome.pathway <- getSchemaClass(class = "Pathway", species = "human", all = TRUE)
## drop disease pathways
reactome.pathway <- as.data.table(reactome.pathway)
reactome.pathway <- reactome.pathway[ isInDisease == F]
## drop top level pathways
reactome.pathway <- reactome.pathway[ schemaClass == "Pathway"]
## get associated genes
reactome.pathway <- mclapply(1:nrow(reactome.pathway), function(x){
  ## obtain relevant genes
  tmp <- event2Ids(reactome.pathway$stId[x])
  ## return information
  return(data.table(reactome.pathway[x,], gene.symbol=tmp$geneSymbol))
}, mc.cores = 10)
## combine
reactome.pathway <- rbindlist(reactome.pathway, fill = T)
reactome.pathway[, name := NULL ]
## make unique
reactome.pathway <- unique(reactome.pathway)
## exclude some very large, possibly unspecific pathways
tail(sort(table(reactome.pathway$displayName)))

## add to olink
olink.pathway[, reactome.pathway.genes := sapply(Protein, function(x){
  ## get pathways the gene is involved in
  tmp <- subset(reactome.pathway, gene.symbol== x)$stId
  print(tmp)
  ## return all unique genes (if any)
  if(length(tmp) > 0){
    return(paste(unique(reactome.pathway[ stId %in% tmp]$gene.symbol), collapse = "|"))
  }else{
    return(NA)
  }
})]

###############################################
####   add information to SCALLOP results  ####
###############################################

## get closets genes for each of the variants
gene.pos    <- fread("<path to dir>/Genes.GRCh37.complete.txt")
## add to unique pQTLs
res.var[, closets.genes := apply(res.var[, .(CHROM, position_hg19)], 1, function(x){
  ## get all genes mapping in a 2Mb window either side
  tmp <- gene.pos[ chromosome_name == x[1] & start_position >= as.numeric(x[2])-1e6 & end_position <= as.numeric(x[2])+1e6 ]
  ## return genes (if any)
  return(paste0(tmp$external_gene_name, collapse = "|"))
})]

## add to SCALLOP
res.scallop <- merge(res.scallop, res.var)

#---------------------------------#
##--     protein complexes     --##
#---------------------------------#

## add a column indicating, whether there is a protein-complex gene among the closest for the detected pQTLs
res.scallop[, protein.complex.gene := apply(res.scallop[, .(Protein, closets.genes)], 1, function(x){
  ## get all possible protein complexes for the protein target
  pr   <- strsplit(x[1], "_")[[1]]
  tmp  <- strsplit(olink.pathway[ Protein %in% pr]$protein.complex.genes, "\\|")[[1]]
  ## included in closest genes?
  cl   <- strsplit(x[2], "\\|")[[1]]
  print(intersect(tmp, cl))
  ## return possible match
  return(paste(intersect(tmp, cl), collapse = "|"))
})]

#---------------------------------#
##--      ligand receptor      --##
#---------------------------------#

## add a column indicating, whether there is a protein-complex gene among the closest for the detected pQTLs
res.scallop[, ligand.receptor.gene := apply(res.scallop[, .(Protein, closets.genes)], 1, function(x){
  ## get all possible protein complexes for the protein target
  pr   <- strsplit(x[1], "_")[[1]]
  tmp  <- strsplit(olink.pathway[ Protein %in% pr]$ligand.receptor.pair, "\\|")[[1]]
  ## included in closest genes?
  cl   <- strsplit(x[2], "\\|")[[1]]
  print(intersect(tmp, cl))
  ## return possible match
  return(paste(intersect(tmp, cl), collapse = "|"))
})]

#---------------------------------#
##--     KEGG pathway gene     --##
#---------------------------------#

## add a column indicating, whether there is a protein-complex gene among the closest for the detected pQTLs
res.scallop[, kegg.pathway.gene := apply(res.scallop[, .(Protein, closets.genes)], 1, function(x){
  ## get all possible protein complexes for the protein target
  pr   <- strsplit(x[1], "_")[[1]]
  tmp  <- strsplit(olink.pathway[ Protein %in% pr]$kegg.pathway.genes, "\\|")[[1]]
  ## included in closest genes?
  cl   <- strsplit(x[2], "\\|")[[1]]
  print(intersect(tmp, cl))
  ## return possible match
  return(paste(intersect(tmp, cl), collapse = "|"))
})]

#---------------------------------#
##--    Reactome pathway gene  --##
#---------------------------------#

## add a column indicating, whether there is a protein-complex gene among the closest for the detected pQTLs
res.scallop[, reactome.pathway.gene := apply(res.scallop[, .(Protein, closets.genes)], 1, function(x){
  ## get all possible protein complexes for the protein target
  pr   <- strsplit(x[1], "_")[[1]]
  tmp  <- strsplit(olink.pathway[ Protein %in% pr]$reactome.pathway.genes, "\\|")[[1]]
  ## included in closest genes?
  cl   <- strsplit(x[2], "\\|")[[1]]
  print(intersect(tmp, cl))
  ## return possible match
  return(paste(intersect(tmp, cl), collapse = "|"))
})]

###############################################
####     refine GWAS catalog assignment    ####
###############################################

## --> drop plasma protein related measures <-- ##

## get all EFO terms to be excluded - 'protein measurement' http://www.ebi.ac.uk/efo/EFO_0004747
efo.exclude           <- gwasrapidd::get_child_efo("EFO_0004747")
efo.exclude           <- paste0("http://www.ebi.ac.uk/efo/", efo.exclude$EFO_0004747)
## drop proteomic and related studies
gwas.catalogue.pruned <- gwas.catalogue[ !(MAPPED_TRAIT_URI %in% c(efo.exclude, "http://www.ebi.ac.uk/efo/EFO_0004747"))]
## some additional manual cleaning (e.g,. )
gwas.catalogue.pruned <- gwas.catalogue.pruned[ !(MAPPED_TRAIT %in% c("insulin like growth factor measurement", "protein measurement", "sex hormone-binding globulin measurement"))]
gwas.catalogue.pruned <- gwas.catalogue.pruned[ !(PUBMEDID %in% c("36349687", "35870639", "30134952"))]
## prune those with missing position
gwas.catalogue.pruned <- gwas.catalogue.pruned[ CHR_POS != ""]
## drop those without mapped trait
gwas.catalogue.pruned <- gwas.catalogue.pruned[ MAPPED_TRAIT != ""]

## --> map to Olink entries <-- ##

## do in parallel
require(doMC)
registerDoMC(12)

## go through each variant and 1) identify all proxies, 2) map to GWAS catalog findings, and 3) reduce redundancy
res.gwas.pruned       <- mclapply(1:nrow(res.var), function(x){
  
  ## get all possible proxies
  snp              <- r2.proxies[ MarkerName_SNP_A == res.var$MarkerName[x]]
  
  ## create snp id to optimize merging multiple mappings
  snp$snp.id       <- paste0(snp$CHR_A, ":", snp$BP_B.hg38)
  
  ## create two versions of mapping
  snp.rsid         <- merge(snp, gwas.catalogue.pruned, by.x="SNP_B", by.y="SNPS")
  snp.pos          <- merge(snp, gwas.catalogue.pruned, by = "snp.id")
  ## edit
  snp.rsid$snp.id   <- snp.rsid$snp.id.x
  snp.rsid$snp.id.x <- snp.rsid$snp.id.y <- NULL
  snp.pos$SNPS     <- snp.pos$SNP_B
  snp.rsid$SNPS    <- snp.rsid$SNP_B
  ## combine
  snp              <- unique(rbind(snp.rsid, snp.pos))
  
  ## prepare return
  if(nrow(snp) > 0){
    ## sort 
    snp              <- snp[order(TRAIT, SNP_A, -R2)]
    ## create indicator
    snp[, ind := 1:.N, by=c("TRAIT", "SNP_A")]
    ## keep only one finding per trait
    snp              <- snp[ind == 1]
    
    ## report summary back
    snp              <- data.table(rsid.gwas=paste(sort(unique(snp$SNP_B)), collapse = "||"),
                                   mapped_trait=paste(sort(unique(snp$MAPPED_TRAIT)), collapse = "||"),
                                   num_reported=length(unique(snp$MAPPED_TRAIT)))
    
  }else{
    
    ## report summary back
    snp              <- data.table(rsid.gwas="",
                                   mapped_trait="",
                                   num_reported=0)
  }
  
  ## return data set
  return(data.table(res.var[x,], snp))
}, mc.cores = 10)
## combine 
res.gwas.pruned       <- rbindlist(res.gwas.pruned)
save.image()

#----------------------------------#
##--       add annotations      --##
#----------------------------------#

## whether r2.group tags cis or trans-pQTLs
res.gwas.pruned[, cis.trans := sapply(R2.group, function(x){
  ifelse("Cis" %in% res.scallop[ R2.group == x]$cis_trans, "cis", "trans")
})]

## --> collate causal gene assignment and whether genes have a biological link
## to the associated protein target <-- 

## create column for gene assignment
registerDoMC(10)
tmp.gene <- mclapply(unique(res.scallop$R2.group), function(x){
  
  ## get all relevant entries
  tmp <- res.scallop[ R2.group == x]
  
  print(x)
  
  ## collate evidence across locus
  if("Cis" %in% tmp$cis_trans){
    ## find cis-target
    jj <- which(tmp$cis_trans == "Cis")[1]
    return(data.table(R2.group=x, candidate.gene=tmp$Protein[jj], candidate.gene.score=2, ligand.receptor.gene=tmp$Protein[jj], protein.complex.gene=tmp$Protein[jj],
                      pathway.gene=tmp$Protein[jj], refined.biological.gene=tmp$Protein[jj]))
  }else{
    ## get candidate gene with highest score
    if(sum(is.na(tmp$score)) == nrow(tmp)){
      can.gene <- unique(unlist(lapply(tmp$candidate_gene_name, function(k) strsplit(k, "\\|")[[1]])))
    }else{
      can.gene <- unique(unlist(lapply(tmp[ score == max(tmp$score, na.rm=T) ]$candidate_gene_name, function(k) strsplit(k, "\\|")[[1]])))
    }
    ## check whether any of those coincide with ligand receptor pairs
    lr.gene  <- unique(unlist(lapply(tmp$ligend.receptor.gene, function(k) strsplit(k, "\\|")[[1]])))
    ## check whether any of those coincide with protein complex genes
    pc.gene  <- unique(unlist(lapply(tmp$protein.complex.gene, function(k) strsplit(k, "\\|")[[1]])))
    ## check for pathway genes
    pw.gene  <- unique(unlist(lapply(c(tmp$kegg.pathway.gene, tmp$reactome.pathway.gene), function(k) strsplit(k, "\\|")[[1]])))
    ## compute refined biological gene
    rb.gene  <- table(c(lr.gene, pc.gene, can.gene))
    ## add evidence from each resource
    rb.gene  <- names(rb.gene[ rb.gene == max(rb.gene)])
    
    ## return
    return(data.table(R2.group=x, 
                      candidate.gene=paste(can.gene, collapse = "|"), 
                      candidate.gene.score=ifelse(sum(is.na(tmp$score)) != nrow(tmp), max(tmp$score, na.rm = T), .25),
                      ligand.receptor.gene=paste(lr.gene, collapse = "|"), 
                      protein.complex.gene=paste(pc.gene, collapse = "|"),
                      pathway.gene=paste(pw.gene, collapse = "|"),
                      refined.biological.gene=paste(rb.gene, collapse = "|")))
  }
  
  
}, mc.cores=10)
## combine
tmp.gene <- rbindlist(tmp.gene)

## --> create set to plot <-- ##

## count R2 groups
tmp                   <- as.data.table(table(res.scallop$R2.group))
names(tmp)            <- c("R2.group", "R2.count")
tmp$R2.group          <- as.numeric(tmp$R2.group)
res.scallop           <- merge(res.scallop, tmp, by="R2.group")

## add to GWAS catalog results
res.gwas.pruned       <- merge(res.gwas.pruned, unique(res.scallop[, .(MarkerName, R2.group, R2.count)]))
## keep only one proxy per R2 group (order to keep best pleiotropy mapping)
res.gwas.pruned       <- res.gwas.pruned[ order(R2.group, -num_reported)]
res.gwas.pruned[, ind := 1:.N, by="R2.group"]
res.gwas.pruned       <- res.gwas.pruned[ ind == 1]
## n = 10461

## import genomic coordinates from GWAS file
gwas.coord            <- fread("<path to dir>")
## order
gwas.coord            <- gwas.coord[ order(chromosome, position)]
## get end positions for different chromosomes
chr.dat               <- gwas.coord[, .(chr.start=1, chr.end=max(position)), by="chromosome"]
## add plotting index
chr.dat[, tmp := cumsum(as.numeric(chr.end))]
chr.dat$plt.start     <- c(1, chr.dat$tmp[-nrow(chr.dat)] + 1)
chr.dat$plt.end       <- chr.dat$plt.start + chr.dat$chr.end - 1
## add position of mids chromosome
chr.dat[, plt.mid := plt.start + (plt.end - plt.start)/2]

## order and create plotting coordinates
res.gwas.pruned[, chromosome := ifelse(CHROM == "X", 23, as.numeric(CHROM))]
## add plotting coordinates
res.gwas.pruned       <- merge(res.gwas.pruned, chr.dat[, .(chromosome, plt.start)], by="chromosome")
## add plotting coordinates
res.gwas.pruned[, plt.srt := plt.start + position_hg19 - 1]
## remove what is no longer needed
rm(gwas.coord); gc(reset=T)

## add gene assignment
res.gwas.pruned       <- merge(res.gwas.pruned, tmp.gene, by="R2.group")
## create indicator about gene
res.gwas.pruned[, gene.indicator := ifelse(cis.trans == "cis", "cis", 
                                           ifelse(ligand.receptor.gene != "" | protein.complex.gene != "", "trans - interaction",
                                                  ifelse(pathway.gene != "", "trans - pathway", "trans - systemic")))]

##############################################
####     prep data for gene assignment    ####
##############################################

## import data matrix provided by Mine
trans.pqtl.class                           <- fread("SCALLOP_UKBB_MA_trans_pqtls_annotations_for_all_genes_within_1Mb_20240523.txt")

## replace NAs with zeros
trans.pqtl.class[ is.na(trans.pqtl.class)] <- 0

## recode ExWAS results
trans.pqtl.class[, ExWAS_burden_lowest_significant_pvalue := ifelse(ExWAS_burden_lowest_significant_pvalue == 0, 0, -log10(ExWAS_burden_lowest_significant_pvalue))]
## recode VEP worst consequence (1 - most severe)
trans.pqtl.class[, VEP_pQTL_worst_rank := ifelse(VEP_pQTL_worst_rank == 0, 28, VEP_pQTL_worst_rank)]
trans.pqtl.class[, VEP_proxy_worst_rank := ifelse(VEP_proxy_worst_rank == 0, 28, VEP_proxy_worst_rank)]

## add R2.groups to create fair splits
trans.pqtl.class                           <- merge(trans.pqtl.class, unique(res.scallop[, .(MarkerName, R2.group)]), by.x="pqtl_MarkerName", by.y="MarkerName")
## add gene density for each
tmp                                        <- trans.pqtl.class[, .(gene.density=length(hgnc_symbol)), by="pqtl_MarkerName"]
trans.pqtl.class                           <- merge(trans.pqtl.class, tmp)

#-------------------------------------#
##--    Ligand/Receptor pairs      --##
#-------------------------------------#

## transform matrix to ease merging
tmp              <- rbindlist(apply(olink.pathway[, .(Protein, ligand.receptor.pair)], 1, function(x) return(data.table(Protein = x[1], ligand.receptor.pair = strsplit(x[2], "\\|")[[1]]))))
## create indicator
tmp[, ligand.receptor.gene := !is.na(ligand.receptor.pair)]
## merge to the trans-pQTL assignment
trans.pqtl.class <- merge(trans.pqtl.class, unique(tmp), by.x=c("Protein", "hgnc_symbol"), by.y=c("Protein", "ligand.receptor.pair"), all.x=T)
## replace missing values
trans.pqtl.class[, ligand.receptor.gene := ifelse(is.na(ligand.receptor.gene), F, ligand.receptor.gene)]

#-------------------------------------#
##--        Protein complex        --##
#-------------------------------------#

## transform matrix to ease merging
tmp              <- rbindlist(apply(olink.pathway[, .(Protein, protein.complex.genes)], 1, function(x) return(data.table(Protein = x[1], protein.complex.genes = strsplit(x[2], "\\|")[[1]]))))
## create indicator
tmp[, protein.complex.gene := !is.na(protein.complex.genes)]
## merge to the trans-pQTL assignment
trans.pqtl.class <- merge(trans.pqtl.class, unique(tmp), by.x=c("Protein", "hgnc_symbol"), by.y=c("Protein", "protein.complex.genes"), all.x=T)
## replace missing values
trans.pqtl.class[, protein.complex.gene := ifelse(is.na(protein.complex.gene), F, protein.complex.gene)]

#-------------------------------------#
##--          KEGG pathway         --##
#-------------------------------------#

## add pathway size to KEGG summary
tmp              <- as.data.table(table(kegg.pathway$PathwayID))
kegg.pathway     <- merge(kegg.pathway, tmp, by.x="PathwayID", by.y="V1")

## re-design to also take pathway size into account
tmp              <- rbindlist(lapply(unique(trans.pqtl.class$Protein), function(x){
  ## get all pathway the protein is involved in
  ii <- kegg.pathway[ Symbol == x ]$PathwayID
  ## only if anything
  if(length(ii) > 0){
    ## get all the pathways
    foo <- kegg.pathway[ PathwayID %in% ii ]
    ## get the size of each involved pathway
    return(data.table(Protein = x, kegg.pathway.pair = foo$Symbol, N = foo$N))
  }else{
    ## set pathway size to the entire genome
    return(data.table(Protein = x, kegg.pathway.pair = NA, N = 2e5))
  }
})) 
## drop missing values
tmp              <- tmp[ !is.na(kegg.pathway.pair)]
## clean to keep only the smallest pathway assignment (if any)
tmp              <- tmp[ order(Protein, kegg.pathway.pair, N)]
tmp              <- tmp[, ind := 1:.N, by=c("Protein", "kegg.pathway.pair")]
tmp              <- tmp[ ind == 1]
## create indicator
tmp[, kegg.pathway.gene := T]
## merge to the trans-pQTL assignment
trans.pqtl.class <- merge(trans.pqtl.class, unique(tmp[, .(Protein, kegg.pathway.pair, N, kegg.pathway.gene)]), by.x=c("Protein", "hgnc_symbol"), by.y=c("Protein", "kegg.pathway.pair"), all.x=T)
## replace missing values
trans.pqtl.class[, kegg.pathway.gene := ifelse(is.na(kegg.pathway.gene), F, kegg.pathway.gene)]
## adopt missing pathway size
trans.pqtl.class[, N.kegg.pathway := ifelse(is.na(N), 2e4, N)]
## delete what is no longer needed
trans.pqtl.class[, N := NULL ]

#-------------------------------------#
##--        REACTOME pathway       --##
#-------------------------------------#

## add pathway size to KEGG summary
tmp              <- as.data.table(table(reactome.pathway$stId))
reactome.pathway <- merge(reactome.pathway, tmp, by.x="stId", by.y="V1")

## redesign to also take pathway size into account
tmp              <- rbindlist(lapply(unique(trans.pqtl.class$Protein), function(x){
  ## get all pathway the protein is involved in
  ii <- reactome.pathway[ gene.symbol == x ]$stId
  ## only if anything
  if(length(ii) > 0){
    ## get all the pathways
    foo <- reactome.pathway[ stId %in% ii ]
    ## get the size of each involved pathway
    return(data.table(Protein = x, reactome.pathway.pair = foo$gene.symbol, N = foo$N))
  }else{
    ## set pathway size to the entire genome
    return(data.table(Protein = x, reactome.pathway.pair = NA, N = 2e4))
  }
})) 
## drop missing values
tmp              <- tmp[ !is.na(reactome.pathway.pair)]
## clean to keep only the smallest pathway assignment (if any)
tmp              <- tmp[ order(Protein, reactome.pathway.pair, N)]
tmp              <- tmp[, ind := 1:.N, by=c("Protein", "reactome.pathway.pair")]
tmp              <- tmp[ ind == 1]
## create indicator
tmp[, reactome.pathway.gene := T]
## merge to the trans-pQTL assignment
trans.pqtl.class <- merge(trans.pqtl.class, unique(tmp[, .(Protein, reactome.pathway.pair, N, reactome.pathway.gene)]), by.x=c("Protein", "hgnc_symbol"), by.y=c("Protein", "reactome.pathway.pair"), all.x=T)
## replace missing values
trans.pqtl.class[, reactome.pathway.gene := ifelse(is.na(reactome.pathway.gene), F, reactome.pathway.gene)]
## adopt missing pathway size
trans.pqtl.class[, N.reactome.pathway := ifelse(is.na(N), 2e4, N)]
## delete what is no longer needed
trans.pqtl.class[, N := NULL ]

#--------------------------------------#
##-- define a set of positive genes --##
#--------------------------------------#

## define pseudo 'true-positive' set
trans.pqtl.class[, pseudo.true.positive := ifelse(VEP_pQTL_worst_rank <= 12 | ligand.receptor.gene == T | protein.complex.gene == T, 1, 0)]

##############################################
####            build classifier          ####
##############################################

## adjust coding of variables
trans.pqtl.class[, kegg.pathway.gene := ifelse(kegg.pathway.gene == T, 1, 0)]
trans.pqtl.class[, reactome.pathway.gene := ifelse(reactome.pathway.gene == T, 1, 0)]
trans.pqtl.class[, chromosome_name := as.numeric(ifelse(chromosome_name == "X", 23, chromosome_name))]
trans.pqtl.class[, protein.complex.gene := ifelse(protein.complex.gene == T, 1, 0)]
trans.pqtl.class[, ligand.receptor.gene := ifelse(ligand.receptor.gene == T, 1, 0)]


## define features to be used
pred.feat        <- c("chromosome_name", "pqtl_pos", "distance_to_start_of_the_gene", "distance_to_end_of_the_gene", "shortest_distance_from_the_gene", "shortest_distance_from_the_gene_v2",
                      "VEP_pQTL_worst_rank", "VEP_proxy_worst_rank", "VEP_proxy_R2", grep("GTEX", names(trans.pqtl.train), value=T), "kegg.pathway.gene", "N.kegg.pathway",
                      "reactome.pathway.gene", "N.reactome.pathway", "gene.density", "ExWAS_burden_lowest_significant_pvalue", "protein.complex.gene", "ligand.receptor.gene")

## import function to run the classifier
source("../scripts/gene_annotator.R")

#-------------------------------------------#
##-- ligand/receptor and protein complex --##
#-------------------------------------------#

## define pseudo true positive set of genes
trans.pqtl.class[, pseudo.true.positive.lr.pc := ifelse(ligand.receptor.gene == 1 | protein.complex.gene == 1, 1, 0)]

## build classifier (omit features used to set true positive set)
pqtl.gene.lr.pc  <- gen.annotator.rf(trans.pqtl.class, pred.feat[-which(pred.feat %in% c("protein.complex.gene", "ligand.receptor.gene"))], "pseudo.true.positive.lr.pc") 

#-------------------------------------------#
##--        functional variants          --##
#-------------------------------------------#

## define pseudo true positive set of genes
trans.pqtl.class[, pseudo.true.positive.vep := ifelse(VEP_pQTL_worst_rank <= 12, 1, 0)]

## build classifier (omit features used to set true positive set)
pqtl.gene.vep    <- gen.annotator.rf(trans.pqtl.class, pred.feat[-which(pred.feat %in% c("VEP_pQTL_worst_rank", "shortest_distance_from_the_gene", "shortest_distance_from_the_gene_v2"))], "pseudo.true.positive.vep") 

#-------------------------------------------#
##--          ExWAS annotation           --##
#-------------------------------------------#

## define pseudo true positive set of genes (increase threshold to meet significance)
trans.pqtl.class[, pseudo.true.positive.exwas := ifelse(ExWAS_burden_lowest_significant_pvalue >= 8.7, 1, 0)]

## build classifier (omit features used to set true positive set)
pqtl.gene.exwas  <- gen.annotator.rf(trans.pqtl.class, pred.feat[-which(pred.feat %in% c("ExWAS_burden_lowest_significant_pvalue"))], "pseudo.true.positive.exwas") 

##############################################
####         examine classifier           ####
##############################################

#-----------------------------#
##--  variable importance  --##
#-----------------------------#

## get variable importance
pqtl.class.imp <- lapply(c("lr.pc", "vep", "exwas"), function(x){
  
  ## get the importance
  tmp <- lapply(get(paste0("pqtl.gene.", x)), function(k){
    ## get importance across all ten runs
    tmp <- varImp(k$pqtl.rf)$importance
    ## return data frame
    return(data.table(predictor=rownames(tmp), var.imp=tmp))
  })
  
  ## combine into single data table
  pqtl.class.imp        <- tmp %>% reduce(inner_join, by="predictor")
  ## edit names
  names(pqtl.class.imp) <- c("predictor", paste("var.imp", 1:10, sep="."))
  ## add median (relative) importance
  pqtl.class.imp[, var.imp.med := apply(pqtl.class.imp[, paste("var.imp", 1:10, sep="."), with=F], 1, median)]
  ## order accordingly
  pqtl.class.imp        <- pqtl.class.imp[ order(-var.imp.med)]
  ## add type of model
  pqtl.class.imp[, model := x]
  ## return
  return(pqtl.class.imp)
  
})
## combine
pqtl.class.imp <- rbindlist(pqtl.class.imp)

#-----------------------------#
##--        plotting       --##
#-----------------------------#

## open device
png("../graphics/pQTL.multiple.gene.annotator.RF.VEP.LigRec.ProComp.20240528.png", width = 16, height = 16, res=300, units = "cm")

## graphical parameters
par(mar=c(1.5,1.5,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, mgp=c(.6,0,0), xaxs="i", yaxs="i")
## define layout for the plot
layout(matrix(c(1,2,2,3,4,4,5,6,6), 3, 3, byrow = T))

## loop through all different instances of predictions
for(j in c("lr.pc", "vep", "exwas")){
  
  ## --> ROCs across validation sets <-- ##
  
  ## adopt figure margins
  par(mar=c(1.5,1.5,.5,.5))
  
  ## empty figure
  plot(c(0,1), c(0,1), type="n", xaxt="n", yaxt="n", xlab="Specificity", ylab="Sensitivity",
       xlim=c(1,0))
  ## add axis and line of random guessing
  axis(1, lwd=.5); axis(2, lwd=.5); abline(a=1, b=-1, lwd=.5)
  ## add reach ROC
  lapply(get(paste0("pqtl.gene.", j)), function(x) points(x$pqtl.roc$specificities, x$pqtl.roc$sensitivities, type="l", lwd=.5, 
                                                          col=colorspace::adjust_transparency("red2", .3)))
  ## add model
  legend("bottomright", cex=.5, bty="n", legend = j)
  
  ## --> variable importance <-- ##
  
  ## get relevant data to ease coding
  tmp <- pqtl.class.imp[ model == j]
  tmp <- tmp[order(-var.imp.med)]
  
  ## increase margin to add names
  par(mar=c(7,1.5,.5,.5))
  
  ## empty figure
  plot(c(.5, nrow(tmp)+.5), c(0,105), type="n", xaxt="n", yaxt="n", xlab="", ylab="Relative importance")
  ## add axis and line of random guessing
  axis(2, lwd=.5); axis(1, lwd=.5, at=1:nrow(tmp), labels = NA)
  ## add median
  arrows(1:nrow(tmp)-.4, tmp[ model == j]$var.imp.med, 1:nrow(tmp)+.4, tmp$var.imp.med, length = 0,
         lwd=1, col="red2")
  ## add individual estimates
  for(j in 1:nrow(tmp)) points(jitter(rep(j, 10), .2), t(tmp[ j, paste0("var.imp.", 1:10), with=T]), cex=.3, 
                               col=colorspace::adjust_transparency("red2", .3), pch=20)
  ## add predictor labels
  pm <- par("usr")
  text(1:nrow(tmp), pm[3]-(pm[4]-pm[3])*.05, cex=.4, pos=2, offset=0, 
       labels=gsub("GTEX_coloc_", "", tmp$predictor), srt=90, 
       xpd=NA)
  
}

## close device
dev.off()


##############################################
####       apply annotation at scale      ####
##############################################

## do some renaming
names(trans.pqtl.class) <- gsub("-", "_", names(trans.pqtl.class))

#---------------------------#
##--    obtain scores    --##
#---------------------------#

## do in parallel
registerDoMC(10)

## --> LR/PC prediction <-- ##
pred.lr.pc <- mclapply(1:10, function(x){
  ## do the prediction
  tmp <- predict(pqtl.gene.lr.pc[[x]]$pqtl.rf, newdata = as.data.frame(trans.pqtl.class[, ..pred.feat]), type = "prob")
  ## return prob for 1
  return(tmp[, 2])
}, mc.cores=10)
## combine
pred.lr.pc <- do.call(cbind, pred.lr.pc)
## convert to data table
pred.lr.pc <- data.table(trans.pqtl.class[, .(pqtl_MarkerName, pheno, R2.group, Protein, hgnc_symbol, ensembl_gene_id, pseudo.true.positive.lr.pc)], pred.lr.pc)
## create median values
pred.lr.pc[, median.score.lr.pc := apply(pred.lr.pc[, paste0("V", 1:10), with=F], 1, median)]

## --> VEP prediction <-- ##
pred.vep   <- mclapply(1:10, function(x){
  ## do the prediction
  tmp <- predict(pqtl.gene.vep[[x]]$pqtl.rf, newdata = as.data.frame(trans.pqtl.class[, ..pred.feat]), type = "prob")
  ## return prob for 1
  return(tmp[, 2])
}, mc.cores=10)
## combine
pred.vep <- do.call(cbind, pred.vep)
## convert to data table
pred.vep <- data.table(trans.pqtl.class[, .(pqtl_MarkerName, pheno, R2.group, Protein, hgnc_symbol, ensembl_gene_id, pseudo.true.positive.vep)], pred.vep)
## create median values
pred.vep[, median.score.vep := apply(pred.vep[, paste0("V", 1:10), with=F], 1, median)]

## --> ExWAS prediction <-- ##
pred.exwas   <- mclapply(1:10, function(x){
  ## do the prediction
  tmp <- predict(pqtl.gene.exwas[[x]]$pqtl.rf, newdata = as.data.frame(trans.pqtl.class[, ..pred.feat]), type = "prob")
  ## return prob for 1
  return(tmp[, 2])
}, mc.cores=10)
## combine
pred.exwas <- do.call(cbind, pred.exwas)
## convert to data table
pred.exwas <- data.table(trans.pqtl.class[, .(pqtl_MarkerName, pheno, R2.group, Protein, hgnc_symbol, ensembl_gene_id, pseudo.true.positive.exwas)], pred.exwas)
## create median values
pred.exwas[, median.score.exwas := apply(pred.exwas[, paste0("V", 1:10), with=F], 1, median)]

#---------------------------#
##--   add to the data   --##
#---------------------------#

## add to pQTL class data
trans.pqtl.class <- merge(trans.pqtl.class, pred.lr.pc[, .(pqtl_MarkerName, pheno, hgnc_symbol, median.score.lr.pc)], by=c("pqtl_MarkerName", "pheno", "hgnc_symbol"))
trans.pqtl.class <- merge(trans.pqtl.class, pred.vep[, .(pqtl_MarkerName, pheno, hgnc_symbol, median.score.vep)], by=c("pqtl_MarkerName", "pheno", "hgnc_symbol"))
trans.pqtl.class <- merge(trans.pqtl.class, pred.exwas[, .(pqtl_MarkerName, pheno, hgnc_symbol, median.score.exwas)], by=c("pqtl_MarkerName", "pheno", "hgnc_symbol"))

## compute sum across all three categories
trans.pqtl.class[, candidate.gene.score.sum := median.score.exwas + median.score.vep + median.score.lr.pc]

##############################################
####          map to SCALLOP data         ####
##############################################

## add to SCALLOP results
tmp <- mclapply(1:nrow(res.scallop), function(x){
  
  print(x)
  
  ## get the relevant subset of findings
  tmp <- trans.pqtl.class[ pqtl_MarkerName == res.scallop$MarkerName[x] & pheno == res.scallop$MA_prot_id[x] ]
  
  ## proceed only if any
  if(nrow(tmp) > 0){
    ## order by aggregated scores
    tmp     <- tmp[order(-candidate.gene.score.sum)]
    ## define possible gaps in assignments
    top.sel <- tmp$candidate.gene.score.sum[-nrow(tmp)] - tmp$candidate.gene.score.sum[-1]  
    ## get the largest difference
    ii       <- which.max(top.sel)
    ## get only those proteins
    if(any(ii)){
      tmp      <- tmp[1:ii]
    }
    ## return results
    return(data.table(res.scallop[x,], 
                      candidate.gene.classifier=paste(tmp$hgnc_symbol, collapse = "|"), 
                      candidate.gene.score.top=tmp$candidate.gene.score[1],
                      candidate.gene.score.all=paste(tmp$candidate.gene.score.sum, collapse = "|"),
                      candidate.gene.number=nrow(tmp)))
  }else{
    return(res.scallop[x,])
  }
}, mc.cores = 10)
tmp <- rbindlist(tmp, fill = T)
## careful, does not capture findings within the MHC region!

## sort
tmp         <- tmp[ order(-candidate.gene.score.top)]

## take care of cis-loci and MHC region
tmp[, candidate.gene.classifier := ifelse(cis_trans == "Cis", Protein, candidate.gene.classifier)]
tmp[, candidate.gene.classifier := ifelse(is.na(candidate.gene.classifier), candidate_gene_name, candidate.gene.classifier)]

## add to SCALLOP results
res.scallop <- tmp
res.scallop[, candidate.gene.classifier := ifelse(cis_trans == "Trans" & is.na(score) , "", candidate.gene.classifier)]

#----------------------------------#
##--      pathway enrichment    --##
#----------------------------------#

## import package
require(gprofiler2)

## --> perform pathway enrichment per protein <-- ##

## run for each ID
pathway.protein <- lapply(unique(res.scallop$pheno), function(x){
  
  print(x)
  
  ## get all relevant candidate genes, exclude cis gene and those that a very pleiotropic
  tmp <- unique(res.scallop[ pheno == x & cis_trans != "Cis" , candidate.gene.classifier])
  ## resolve ambiguous loci
  tmp <- unique(unlist(lapply(tmp, function(k) strsplit(k, "\\|")[[1]][1])))
  
  ## proceed only if more than 3 genes
  if(length(tmp) >= 3){
    
    ## do the enrichment
    en.res <- gost(query = tmp, 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                   measure_underrepresentation = FALSE, evcodes = TRUE, 
                   user_threshold = 0.05, correction_method = "fdr", 
                   domain_scope = "annotated",
                   numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
    ## return results
    return(data.table(pheno = x, en.res$result))
  }  
  
})
## combine everything into one large data frame
pathway.protein <- do.call(plyr::rbind.fill, pathway.protein)
## convert to data table
pathway.protein <- as.data.table(pathway.protein)
## drop findings with no enriched pathways
pathway.protein <- pathway.protein[ !(is.na(p_value))]
## minimum interaction size
pathway.protein <- pathway.protein[ intersection_size > 2]
## drop unspecific terms
pathway.protein <- pathway.protein[ term_size > 2 & term_size < 400]
## add more information on proteins
pathway.protein <- merge(unique(res.scallop[, .(Protein, pheno)]), pathway.protein)
## add fold change
pathway.protein[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]
## delete some colums
pathway.protein[, parents := NULL]

## create indicator whether the protein participates in the relevant pathway 
pathway.protein[, cis.included := apply(pathway.protein[, .(Protein, term_id, source)], 1, function(x){
  
  ## check whether the term occurs in the Olink mapping
  if(x[3] == "REAC"){
    nrow(reactome.pathway[ gene.symbol == x[1] & stId == gsub("REAC:", "", x[2], fixed=T)])
  }else{
    nrow(kegg.pathway[ Symbol == x[1] & gsub("hsa", "", PathwayID) == gsub("KEGG:", "", x[2], fixed=T)])
  }
  
})]

## plot some (reduce redundancy)
source("../scripts/plot_enrichment.R")

## create compressed version
pathway.protein.compressed <- mclapply(unique(pathway.protein$pheno), function(x) plot.enrich(pathway.protein,  s.set = x, compress = T), mc.cores = 10)
pathway.protein.compressed <- rbindlist(pathway.protein.compressed)

## --> pathway enrichment for pleiotropic pathways <-- ##

## create column for gene assignment
tmp.gene  <- mclapply(unique(res.scallop$R2.group), function(x){
  
  ## get all relevant entries
  tmp <- res.scallop[ R2.group == x]
  
  print(x)
  
  ## collate evidence across locus
  if("Cis" %in% tmp$cis_trans){
    ## define which is cis
    jj <- which(tmp$cis_trans == "Cis")[1]
    return(data.table(R2.group=x, candidate.gene.classifier=tmp$Protein[jj]))
  }else{
    ## get candidate gene with highest score
    if(sum(is.na(tmp$score)) == nrow(tmp)){
      can.gene <- unique(unlist(lapply(tmp$candidate_gene_name, function(k) strsplit(k, "\\|")[[1]])))
      ## return
      return(data.table(R2.group=x, candidate.gene.classifier=paste(can.gene, collapse = "|")))
    }else{
      ## split to average for each candidate gene (use all gene scores instead of top)
      can.gene <- trans.pqtl.class[ R2.group == x , .(can.score = mean(candidate.gene.score.sum)), by="hgnc_symbol"]
      ## order
      can.gene <- can.gene[order(-can.score)]
      ## define possible gaps in assignments
      top.sel  <- can.gene$can.score[-nrow(can.gene)] - can.gene$can.score[-1]  
      ## get the largest difference
      ii       <- which.max(top.sel)
      ## get only those proteins
      if(any(ii)){
        can.gene <- can.gene[1:ii]
      }
      ## return results
      return(data.table(R2.group=x, 
                        candidate.gene.classifier=paste(can.gene$hgnc_symbol, collapse = "|"), 
                        candidate.gene.score.top=can.gene$can.score[1],
                        candidate.gene.score.all=paste(can.gene$can.score, collapse = "|"),
                        candidate.gene.number=nrow(can.gene)))
    }
    
  }
  
  
}, mc.cores = 10)
## combine
tmp.gene  <- rbindlist(tmp.gene, fill = T)
## add R2 count
tmp.gene  <- merge(tmp.gene, unique(res.scallop[, .(R2.group, R2.count)])) 

#-----------------------------------#
##-- cell-type/tissue enrichment --##
#-----------------------------------#

## import annotations
hpa.tissue <- fread("<path to dir>/HPA.all.genes.tissue.enhanced.20240501.txt")
hpa.cells  <- fread("<path to dir>/HPA.all.genes.cell.enhanced.20240501.txt")
## import including anatomical ordering
tissues    <- fread("<path to dir>/Olink.tissue.annotation.summary.20240423.txt")
cell.type  <- fread("<path to dir>/Olink.cell.type.annotation.summary.20240423.txt")

## adjust 'olink.included' to reflect SCALLOP subset
hpa.tissue[, olink.included := ifelse(Gene %in% res.scallop$Protein, 1, 0)]
hpa.cells[, olink.included := ifelse(Gene %in% res.scallop$Protein, 1, 0)]

## --> enrichment of tissues among trans genes <-- ##

## do in parallel
registerDoMC(10)

## test whether proteins explained by certain factors are enriched for 
## expression in certain tissues
trans.gene.tissue <- mclapply(unique(res.scallop$pheno), function(x){
  
  print(x)
  
  ## get the proteins of interest (get only top gene at each locus!)
  prot.diff <- unique(unlist(lapply(res.scallop[ pheno == x & cis_trans == "Trans" ]$candidate.gene.classifier, function(x) strsplit(x, "\\|")[[1]][1])))
  
  ## only if at least three different genes
  if(length(prot.diff) > 3){
    ## test for the enrichment across different tissues
    enr       <- lapply(tissues$tissue, function(k){
      
      ## selected and tissue specific
      d1    <- nrow(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
      ## not selected and tissue specific
      d2    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
      ## selected and not tissue specific
      d3    <- nrow(hpa.tissue[ Gene %in% prot.diff & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
      ## not selected and not tissue specific
      d4    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
      
      ## test for enrichment
      enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
      
      ## return information needed
      return(data.table(tissue=k, pheno=x, or=enr$estimate, pval=enr$p.value, 
                        intersection=paste(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                        d1=d1, d2=d2, d3=d3, d4=d4))
    })
    enr       <- rbindlist(enr, fill=T)
    return(enr)
  }
  
}, mc.cores=10)
## combine results
trans.gene.tissue <- rbindlist(trans.gene.tissue, fill = T)
## add more information on proteins
trans.gene.tissue <- merge(unique(res.scallop[, .(Protein, pheno)]), trans.gene.tissue)
## compute FDR
trans.gene.tissue[, fdr := p.adjust(pval, method = "BH")]
## prune for single proteins (drop pheno_ID targeting the same protein)
trans.gene.tissue <- trans.gene.tissue[ order(Protein, tissue, pval)]
trans.gene.tissue[, ind := 1:.N, by=c("Protein", "tissue")]
## keep only those
trans.gene.tissue <- trans.gene.tissue[ ind == 1]

## add colour and ordering
trans.gene.tissue <- merge(trans.gene.tissue, tissues[, .(tissue, group, colour, label, order)])

## --> enrichment of cell-types among trans genes <-- ##

## test whether proteins explained by certain factors are enriched for 
## expression in certain tissues
trans.gene.cells  <- mclapply(unique(res.scallop$pheno), function(x){
  
  print(x)
  
  ## get the proteins of interest (take only top ranking gene at each locus)
  prot.diff <- unique(unlist(lapply(res.scallop[ pheno == x & cis_trans == "Trans" ]$candidate.gene.classifier, function(x) strsplit(x, "\\|")[[1]][1])))
  
  ## do only if at least three genes
  if(length(prot.diff) >= 3){
    
    ## test for the enrichment across different tissues
    enr       <- lapply(cell.type$cell.type, function(k){
      
      ## selected and tissue specific
      d1    <- nrow(hpa.cells[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
      ## not selected and tissue specific
      d2    <- nrow(hpa.cells[ !(Gene %in% prot.diff) & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
      ## selected and not tissue specific
      d3    <- nrow(hpa.cells[ Gene %in% prot.diff & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
      ## not selected and not tissue specific
      d4    <- nrow(hpa.cells[ !(Gene %in% prot.diff) & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
      
      ## test for enrichment
      enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
      
      ## return information needed
      return(data.table(cell.type=k, pheno=x, or=enr$estimate, pval=enr$p.value, 
                        intersection=paste(hpa.cells[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                        d1=d1, d2=d2, d3=d3, d4=d4))
    })
    enr       <- rbindlist(enr, fill=T)
    return(enr)
  }
  
}, mc.cores=10)
## combine results
trans.gene.cells  <- rbindlist(trans.gene.cells, fill = T)
## add more information on proteins
trans.gene.cells  <- merge(unique(res.scallop[, .(Protein, pheno)]), trans.gene.cells)
## compute FDR
trans.gene.cells[, fdr := p.adjust(pval, method = "BH")]
## prune for single proteins (drop pheno_ID targeting the same protein)
trans.gene.cells <- trans.gene.cells[ order(Protein, cell.type, pval)]
trans.gene.cells[, ind := 1:.N, by=c("Protein", "cell.type")]
## keep only those
trans.gene.cells <- trans.gene.cells[ ind == 1]

## add colour and ordering
trans.gene.cells <- merge(trans.gene.cells, cell.type[, .(cell.type, category, colour, label, order)])

#----------------------------------#
##--      druggable proteins    --##
#----------------------------------#

## import mapping
drug.genes  <- fread("../../target_list_grch37_final.txt")

## add to SCALLOP whether cis-gene is druggable (include druggable state)
res.scallop[, cis.druggable := apply(res.scallop[, .(Protein, Uniprot_ID)], 1, function(x){
  ## search for possible entries in the druggable genome
  symbol  <- strsplit(x[1], "_")[[1]]
  uniprot <- strsplit(x[2], "_")[[1]]
  ## anything in the druggable table
  tmp     <- drug.genes[ molecue_summary != "" & (gene_display_label %in% symbol | accession %in% uniprot)]
  ## return
  return(ifelse(nrow(tmp) > 0, max(tmp$max_molecule_phase_for_target, na.rm=T), -1))
})]

## add to SCALLOP whether trans-gene is druggable (include druggable state)
res.scallop[, trans.druggable := sapply(res.scallop$candidate.gene.classifier, function(x){
  ## search for possible entries in the druggable genome
  symbol  <- strsplit(x, "_")[[1]]
  ## anything in the druggable table
  tmp     <- drug.genes[ molecue_summary != "" & gene_display_label %in% symbol ]
  ## return
  return(ifelse(nrow(tmp) > 0, max(tmp$max_molecule_phase_for_target, na.rm=T), -1))
})]

## add information to pathway results
pathway.protein.compressed <- merge(pathway.protein.compressed, unique(res.scallop[, .(pheno, cis.druggable)]))
pathway.protein            <- merge(pathway.protein, unique(res.scallop[, .(pheno, cis.druggable)]))

#------------------------------------#
##-- summary figure for the draft --##
#------------------------------------#

## combine
png("../graphics/Trans.pqtl.annotation.pathway.enrichment.png", width = 16, height = 16/3, res=300, units = "cm")
## graphical parameters
par(mar=c(1.5,1.5,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, mgp=c(.6,0,0), xaxs="i", yaxs="i")
## define layout for the plot
layout(matrix(c(1,2,2), 1, 3, byrow = T))

## --> distribution gene scores <-- ##

hist(res.scallop$candidate.gene.score.top, breaks=50, lwd=.5, col=colorspace::adjust_transparency("orange", .5), xlab="Candidate gene score", ylab="Frequency",
     main="")

## --> enriched pathways <-- ##

## create relevant data
tmp <- pathway.protein.compressed[, .(count = length(pheno), count.drug = sum(cis.druggable >= 0)), by="term_name"]
## order
tmp <- tmp[order(-count)]
## chhose only 3 or more
tmp <- tmp[ count > 2]

## empty plot
plot(c(.5, nrow(tmp)+.5), c(0,80), xlab="Pathway", ylab="Number of proteins enriched for trans genes", type="n", xaxt="n", yaxt="n")
## add axis
axis(2, lwd=.5)
## pathways - all "#00A4CC"
rect(1:nrow(tmp)-.4, 0, 1:nrow(tmp)+.4, tmp$count, col=colorspace::adjust_transparency("#F95700", .5), lwd=.1, border="grey30")
## pathways - druggable"#00A4CC"
rect(1:nrow(tmp)-.4, 0, 1:nrow(tmp)+.4, tmp$count.drug, col=colorspace::adjust_transparency("#00A4CC", .5), lwd=.1, border="grey30")
## add pathway names
text(1:nrow(tmp), 1, srt=90, pos=4, labels = sapply(tmp$term_name, function(x) paste(strwrap(stringr::str_to_sentence(x), 70), collapse="\n")), cex=.4,
     font=1, offset=0, xpd=NA, col="grey30")
## add legend
legend("topright", bty="n", cex=.5, pch=22, pt.cex=1, pt.lwd=.1, col="grey30", pt.bg=c("#F95700", "#00A4CC"),
       legend = c("protein targets", "druggable protein targets"))
## close
dev.off()

## --> tissue/cell-types <-- ##

## start device
pdf("../graphics/Summary.trans.genes.tissue.cell.type.enrichment.202405230.pdf", width = 6.3, height = 3.15)
## define graphical parameters
par(mfrow=c(1,2), mar=rep(0,4))

#----------------------------#
##--   tissue enrichment  --##
#----------------------------#

circos.clear()

## create input
tmp               <- trans.gene.tissue[ fdr < .05 & or > 1 & d1 > 1, .(label, Protein, colour, order)]
## dummy
tmp[, value := 1]
## new order
tmp                <- tmp[, .(label, Protein, value, order, colour)]
tmp                <- tmp[ order(order)]
## define vector to customize the plot
chord.order        <- rev(c(unique(tmp$label), rev(unique(tmp$Protein))))
## colour
chord.color        <- rev(c(tissues[ label %in% tmp$label]$colour, rep("grey80", length(unique(tmp$Protein)))) )
## assign names
names(chord.color) <- chord.order

## define gap parameters
circos.par(gap.after = c(rep(.5, length(unique(tmp[[2]]))-1), 3, 
                         rep(1.5, length(unique(tmp[[1]]))-1), 3),
           start.degree = 160)

## create diagram
chordDiagram(tmp, 
             ## define order based on grouping of labels
             order = chord.order,
             ## define colour vector
             grid.col = chord.color,
             annotationTrack = "grid", 
             annotationTrackHeight = .01,
             ## larger layer for labels
             preAllocateTracks = list(track.height = .35),
             ## colour links according to characteristics
             col = tmp$colour)

## change orientation of labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = .25)
  # circos.axis(h = "top", labels.cex = .3, major.tick.length = .01, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()

#----------------------------#
##--    cell enrichment   --##
#----------------------------#

## create input
tmp               <- trans.gene.cells[ fdr < .05 & or > 1 & d1 > 1, .(label, Protein, colour, order)]
## dummy
tmp[, value := 1]
## new order
tmp                <- tmp[, .(label, Protein, value, order, colour)]
tmp                <- tmp[ order(order)]
## define vector to customize the plot
chord.order        <- rev(c(unique(tmp$label), rev(unique(tmp$Protein))))
## colour
chord.color        <- rev(c(cell.type[ label %in% tmp$label]$colour, rep("grey80", length(unique(tmp$Protein)))) )
## assign names
names(chord.color) <- chord.order

## define gap parameters
circos.par(gap.after = c(rep(.5, length(unique(tmp[[2]]))-1), 3, 
                         rep(1.5, length(unique(tmp[[1]]))-1), 3),
           start.degree = 160)

## create diagram
chordDiagram(tmp, 
             ## define order based on grouping of labels
             order = chord.order,
             ## define colour vector
             grid.col = chord.color,
             annotationTrack = "grid", 
             annotationTrackHeight = .01,
             ## larger layer for labels
             preAllocateTracks = list(track.height = .35),
             ## colour links according to characteristics
             col = tmp$colour)

## change orientation of labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = .2)
  # circos.axis(h = "top", labels.cex = .3, major.tick.length = .01, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()

## close
dev.off()

##############################################
####   add gene assignments to GWAS cat.  ####
##############################################

## add to GWAS catalog results
res.gwas.pruned <- merge(res.gwas.pruned, tmp.gene[, .(R2.group, candidate.gene.classifier, candidate.gene.score.top, candidate.gene.score.all, candidate.gene.number)])

## create second count based on refined set of GWAS catalog entries

#--------------------------------------------#
##-- add parent EFO terms to GWAS catalog --##
#--------------------------------------------#

require(ontologyIndex)

## import EFO ontology
efo.ontology <- get_ontology("efo.obo")

## add to pruned GWAS catalog results
gwas.catalogue.pruned[, efo.term := sapply(MAPPED_TRAIT_URI, function(x){
  ## split
  x <- strsplit(x, ",\\s*")[[1]]
  ## extract what is needed
  x <- sapply(x, function(k) sub(".*\\/([^\\/]+)$", "\\1", k))
  ## return results
  return(paste(x, collapse = "|"))
})]
## make consistent
gwas.catalogue.pruned[, efo.term := gsub("_", ":", efo.term)]

## call parent terms (at least for those where it is available)
gwas.catalogue.pruned[, efo.parents := sapply(efo.term, function(x){
  ## get all possible entries
  efo <- strsplit(x, "\\|")[[1]]
  ## get parents for each
  efo <- unique(unlist(lapply(efo, function(k) efo.ontology$parents[[k]])))
  ## return results
  if(length(efo) > 0){
    return(paste(efo, collapse = "|"))
  }else{
    return(paste(x, collapse = "|"))
  }
})]

## second step of refinement, get all EFO parents and rank occurrence
tmp.efo <- as.data.table(table(unlist(mclapply(gwas.catalogue.pruned$efo.parents, function(x) strsplit(x, "\\|")[[1]], mc.cores=10))))
## add label
tmp.efo[, label := efo.ontology$name[V1]]

## go through and choose the efo term most frequently found as parent (take original term otherwise)
gwas.catalogue.pruned[, efo.parent.pruned := sapply(efo.parents, function(x){
  ## get all relevant terms
  efo <- strsplit(x, "\\|")[[1]]
  ## get the subset
  efo <- tmp.efo[ V1 %in% efo]
  ## return the most frequent one
  return(efo$V1[which.max(efo$N)])
})]

#------------------------------------------------------#
##-- redo r2 mapping to find number of associations --##
#------------------------------------------------------#

## add to pQTL mapping (need to sort out mapping...)
## do in parallel
require(doMC)
registerDoMC(12)

## go through each variant and 1) identify all proxies, 2) map to GWAS catalog findings, and 3) reduce redundancy
res.gwas.efo    <- mclapply(1:nrow(res.var), function(x){
  
  ## get all possible proxies
  snp              <- r2.proxies[ MarkerName_SNP_A == res.var$MarkerName[x]]
  
  ## create snp id to optimize merging multiple mappings
  snp$snp.id       <- paste0(snp$CHR_A, ":", snp$BP_B.hg38)
  
  ## create two versions of mapping
  snp.rsid         <- merge(snp, gwas.catalogue.pruned, by.x="SNP_B", by.y="SNPS")
  snp.pos          <- merge(snp, gwas.catalogue.pruned, by = "snp.id")
  ## edit
  snp.rsid$snp.id   <- snp.rsid$snp.id.x
  snp.rsid$snp.id.x <- snp.rsid$snp.id.y <- NULL
  snp.pos$SNPS     <- snp.pos$SNP_B
  snp.rsid$SNPS    <- snp.rsid$SNP_B
  ## combine
  snp              <- unique(rbind(snp.rsid, snp.pos))
  
  ## prepare return
  if(nrow(snp) > 0){
    ## sort 
    snp              <- snp[order(TRAIT, SNP_A, -R2)]
    ## create indicator
    snp[, ind := 1:.N, by=c("TRAIT", "SNP_A")]
    ## keep only one finding per trait
    snp              <- snp[ind == 1]
    
    ## report summary back
    snp              <- data.table(rsid.gwas=paste(sort(unique(snp$SNP_B)), collapse = "||"),
                                   efo.term=paste(sort(unique(snp$efo.term)), collapse = "||"),
                                   efo.parent=paste(sort(unique(snp$efo.parent.pruned)), collapse = "||"),
                                   num.efo.parent=length(unique(snp$efo.parent.pruned)))
    
  }else{
    
    ## report summary back
    snp              <- data.table(rsid.gwas="",
                                   efo.term="",
                                   efo.parent="",
                                   num.efo.parent=0)
  }
  
  ## return data set
  return(data.table(res.var[x,], snp))
}, mc.cores = 10)
## combine 
res.gwas.efo    <- rbindlist(res.gwas.efo)

## add to the previous iteration focused on mapped traits
res.gwas.pruned <- merge(res.gwas.pruned, unique(res.gwas.efo[, .(MarkerName, efo.term, efo.parent, num.efo.parent)]), by="MarkerName")

#------------------------------------------------------#
##--             introduce pQTL categories          --##
#------------------------------------------------------#

## create arbritary grouping for now
res.gwas.pruned[, pqtl.category := ifelse(R2.count <= 5 & num.efo.parent <= 5, "specific",
                                          ifelse(R2.count <= 5 & num.efo.parent > 5, "systemic.pleiotropy",
                                                 ifelse(R2.count > 5 & num.efo.parent <= 5, "protein.pleiotropy", "unspecific.pleiotropy")))]

#-------------------------------------------------------#
##-- export results for MR analysis and more general --##
#-------------------------------------------------------#

## export including GWAS catalog annotation
write.table(res.gwas.pruned, "SCALLOP.R2.groups.GWAS.catalog.pleiotropy.20240627.txt", sep="\t", row.names = F)

#######################################################################################################
#######################################################################################################
####                                  REVISION - CELL 02/07/2025                                   ####
#######################################################################################################
#######################################################################################################

################################################
####    phenotypic relevance trans-pQTLs    ####
################################################

## export to combine with biomarker analysis
write.table(res.scallop, "Results.SCALLOP.GWAS.catalog.overlap.pathway.genes.20250711.txt", sep="\t", row.names = F)
write.table(res.gwas.pruned, "Results.R2.groupGWAS.catalog.pruned.overlap.pathway.genes.20250711.txt", sep="\t", row.names = F)

#----------------------------#
##-- enrichment cis/trans --##
#----------------------------#

## import genomic position of Olink targets (careful: coordinates are in build 37)
prot.label <- fread("MA_protein_list_and_locations.txt")
## expand by gene
prot.long  <- prot.label[, 1:17]
# Add row identifier
prot.long[, row_id := .I]
## edit one protein ensembl id
prot.long[, ensembl_gene_id := gsub("ENSG00000205595", "ENSG00000109321", ensembl_gene_id)]

# Split rows based on pipe separator - corrected approach
prot.long <- prot.long[, {
  # Find columns with pipe separators
  pipe_cols <- which(sapply(.SD, function(x) any(grepl("\\|", x))))
  
  if(length(pipe_cols) > 0) {
    # Get the first pipe column to determine split length
    first_pipe_col <- names(.SD)[pipe_cols[1]]
    split_length <- length(strsplit(.SD[[first_pipe_col]], "\\|")[[1]])
    
    # Create result list
    result_list <- list()
    
    # Process each new row
    for(i in 1:split_length) {
      new_row <- copy(.SD)
      
      # Update pipe columns
      for(col_idx in pipe_cols) {
        col_name <- names(.SD)[col_idx]
        split_vals <- strsplit(.SD[[col_name]], "\\|")[[1]]
        if(i <= length(split_vals)) {
          new_row[[col_name]] <- split_vals[i]
        } else {
          new_row[[col_name]] <- split_vals[length(split_vals)]
        }
      }
      
      result_list[[i]] <- new_row
    }
    
    rbindlist(result_list)
  } else {
    .SD
  }
}, by = row_id]

## add position in build 38; do two times due to missing identifiers
prot.long <- merge(prot.long, unique(human.genes[, .(gene_id, start, end)]), by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)
## now based on ENSEMBL
prot.long <- merge(prot.long, unique(human.genes[, .(gene_name, start, end)]), by.x = "hgnc_symbol", by.y = "gene_name", all.x = T, suffixes = c(".ensembl", ".hgnc"))
## fuse
prot.long[, start.hg38 := ifelse(!is.na(start.ensembl), start.ensembl, start.hgnc)]
prot.long[, end.hg38 := ifelse(!is.na(end.ensembl), end.ensembl, end.hgnc)]
## drop columns
prot.long[, start.ensembl := NULL]
prot.long[, end.ensembl := NULL]
prot.long[, start.hgnc := NULL]
prot.long[, end.hgnc := NULL]
## clean redundancies
prot.long <- unique(prot.long)

## create indicator for the GWAS catalog, whether the variant falls into a potential cis-pQTL region
gwas.catalogue.pruned[, protein.region := apply(gwas.catalogue.pruned[, .(CHR_ID, CHR_POS)], 1, function(x){
  ## get all proteins falling into that region
  tmp <- prot.long[ chromosome_name == x[1] & start.hg38 - 1e6 <= as.numeric(x[2])  & end.hg38 + 1e6 >= as.numeric(x[2])]
  ## return mapping proteins
  return(paste(sort(unique(tmp$Assay)), collapse = "|"))
})]
## add indicator whether the variant was linked to a cis-pQTL
jj <- unique(unlist(lapply(res.gwas.pruned[ cis.trans == "cis"]$rsid.gwas, function(x) strsplit(x, "\\|\\|")[[1]])))
gwas.catalogue.pruned[, cis.pQTL.linked := SNPS %in% jj]
## look out for oddities
gwas.catalogue.pruned[ protein.region == "" & cis.pQTL.linked == T]

#-------------------------------#
##-- probability of GWAS hit --##
#-------------------------------#

## enrichment of trans-pQTLs?
summary(glm(ifelse(num_reported > 0, 1, 0) ~ cis.trans, res.gwas.pruned[ R2.group != 0], family = binomial(link = "logit")))$coefficients

## number reported in the GWAS catalog; model as negbinomial: omit MHC region
summary(zeroinfl(num_reported ~ cis.trans , dist = "negbin", data = res.gwas.pruned[ R2.group != 0]))
## clear evidence that trans pQTLs are associated with more reported traits, but not if any
summary(zeroinfl(num.efo.parent ~ cis.trans , dist = "negbin", data = res.gwas.pruned[ R2.group != 0]))
## same when collapsing into parent terms

## numbers for the paper
summary(zeroinfl(num_reported ~ cis.trans , dist = "negbin", data = res.gwas.pruned[ R2.group != 0]))$coefficients

## look whether pleiotropy in the proteome matters: omit MHC region
summary(zeroinfl(num_reported ~ cis.trans * R2.count , dist = "negbin", data = res.gwas.pruned[ R2.group != 0]))
## evidence that stratification is needed
summary(zeroinfl(num_reported ~ R2.count , dist = "negbin", data = res.gwas.pruned[ cis.trans == "cis" & R2.group != 0]))
## R2 count is associated with higher number of reported traits in the GWAS catalog
summary(zeroinfl(num_reported ~ R2.count , dist = "negbin", data = res.gwas.pruned[ cis.trans == "trans" & R2.group != 0]))

##########################################################
####     additional analysis pleiotropy section       ####
##########################################################

#----------------------------------------#
##--  tissue/cell enrichment testing  --##
#----------------------------------------#

## loop thorugh all categories
enrich.pleio.tissue <- lapply(unique(res.gwas.pruned$pqtl.category), function(x){
  
  ## get all relevant candidate genes, exclude cis gene and those that a very pleiotropic
  tmp <- unique(res.gwas.pruned[ pqtl.category == x , candidate.gene.classifier])
  ## resolve ambiguous loci
  tmp <- unique(unlist(lapply(tmp, function(k) strsplit(k, "\\|")[[1]][1])))
  
  ## define all genes prioritised as background
  bg  <- unique(res.gwas.pruned[, candidate.gene.classifier])
  ## resolve ambiguous loci
  bg  <- unique(unlist(lapply(bg, function(k) strsplit(k, "\\|")[[1]][1])))
  
  ## run enrichment across tissues
  res <- lapply(tissues$tissue, function(k){
    
    ## selected and tissue specific
    d1    <- nrow(hpa.tissue[ Gene %in% tmp & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## not selected and tissue specific
    d2    <- nrow(hpa.tissue[ Gene %in% bg & !(Gene %in% tmp) & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## selected and not tissue specific
    d3    <- nrow(hpa.tissue[ Gene %in% tmp & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    ## not selected and not tissue specific
    d4    <- nrow(hpa.tissue[ Gene %in% bg & !(Gene %in% tmp) & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))])
    
    ## test for enrichment
    enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
    
    ## return information needed
    return(data.table(tpqtl.category = x, tissue=k, or=enr$estimate, pval=enr$p.value, 
                      intersection=paste(hpa.tissue[ Gene %in% tmp & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                      d1=d1, d2=d2, d3=d3, d4=d4))
  })
  res <- rbindlist(res, fill=T)
  return(res)
})
## combine and check
enrich.pleio.tissue <- rbindlist(enrich.pleio.tissue)
## add FDR
enrich.pleio.tissue[, fdr := p.adjust(pval, method =  "BH")]

## loop through all categories
enrich.pleio.cells  <- lapply(unique(res.gwas.pruned$pqtl.category), function(x){
  
  ## get all relevant candidate genes, exclude cis gene and those that a very pleiotropic
  tmp <- unique(res.gwas.pruned[ pqtl.category == x , candidate.gene.classifier])
  ## resolve ambiguous loci
  tmp <- unique(unlist(lapply(tmp, function(k) strsplit(k, "\\|")[[1]][1])))
  
  ## define all genes prioritised as background
  bg  <- unique(res.gwas.pruned[, candidate.gene.classifier])
  ## resolve ambiguous loci
  bg  <- unique(unlist(lapply(bg, function(k) strsplit(k, "\\|")[[1]][1])))
  
  ## run enrichment across tissues
  res <- lapply(cell.type$cell.type, function(k){
    
    ## selected and tissue specific
    d1    <- nrow(hpa.cells[ Gene %in% tmp & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched")& !is.na(eval(as.name(k))))])
    ## not selected and tissue specific
    d2    <- nrow(hpa.cells[ Gene %in% bg & !(Gene %in% tmp) & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## selected and not tissue specific
    d3    <- nrow(hpa.cells[ Gene %in% tmp & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    ## not selected and not tissue specific
    d4    <- nrow(hpa.cells[ Gene %in% bg & !(Gene %in% tmp) & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))])
    
    ## test for enrichment
    enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
    
    ## return information needed
    return(data.table(tpqtl.category = x, cell.type=k, or=enr$estimate, pval=enr$p.value, 
                      intersection=paste(hpa.cells[ Gene %in% tmp & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                      d1=d1, d2=d2, d3=d3, d4=d4))
  })
  res <- rbindlist(res, fill=T)
  return(res)
})
## combine and check
enrich.pleio.cells  <- rbindlist(enrich.pleio.cells)
## add FDR
enrich.pleio.cells[, fdr := p.adjust(pval, method =  "BH")]

#----------------------------------------#
##--    Protein-Protein Interaction   --##
#----------------------------------------#

## --> effector gene assignments <-- ##

## define reference
res.gwas.pruned[, pqtl.category.factor := factor(pqtl.category, levels = c("unspecific.pleiotropy", "specific", "protein.pleiotropy", "systemic.pleiotropy"))]
res.gwas.pruned[, pqtl.category.factor := factor(pqtl.category, levels = c( "specific", "unspecific.pleiotropy", "protein.pleiotropy", "systemic.pleiotropy"))]
table(res.gwas.pruned$pqtl.category.factor)

## any difference in trans-pQTLs for pathways
summary(glm(ifelse(pathway.gene != "", 1, 0) ~ pqtl.category.factor + R2.count, res.gwas.pruned[ cis.trans == "trans"], family = binomial(link = "logit")))$coefficients

## depletion of complex partners
summary(glm(ifelse(protein.complex.gene != "", 1, 0) ~ pqtl.category.factor + R2.count, res.gwas.pruned[ cis.trans == "trans"], family = binomial(link = "logit")))$coefficients

## enrichment of ligand receptor pairs
summary(glm(ifelse(ligand.receptor.gene != "", 1, 0) ~ pqtl.category.factor + R2.count, res.gwas.pruned[ cis.trans == "trans"], family = binomial(link = "logit")))$coefficients

#----------------------------------------#
##--     protein pathway enrichment   --##
#----------------------------------------#

## test pathway enrichment by R2 group
r2.group.pathway   <- lapply(1:nrow(res.gwas.pruned), function(x){
  
  ## proceed only if 3 or more proteins
  if(res.gwas.pruned$R2.count[x] >= 3){
    
    ## get protein profile
    prot.fore <- strsplit(res.gwas.pruned$protein.profile[x], "\\|")[[1]]
    
    ## do the enrichment
    en.res   <- gost(query = prot.fore, 
                     organism = "hsapiens", ordered_query = FALSE, 
                     multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                     measure_underrepresentation = FALSE, evcodes = TRUE, 
                     user_threshold = 0.05, correction_method = "fdr", 
                     custom_bg = unique(prot.meta$Protein),
                     domain_scope = "annotated",
                     numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
    ## return results
    return(data.table(res.gwas.pruned[x, .(R2.group, R2.count, candidate.gene.classifier, candidate.gene.score.all, cis.trans, protein.profile, pqtl.category)], en.res$result))
  }
})
## collate results
r2.group.pathway   <- rbindlist(r2.group.pathway, fill=T)
## add fold change
r2.group.pathway[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]
## delete some colums
r2.group.pathway[, parents := NULL]
## delete entries with no results
r2.group.pathway   <- r2.group.pathway[ !is.na(p_value)]
save.image()

## add GWAS and other results
r2.group.pathway   <- merge(r2.group.pathway, res.gwas.pruned)
## filter
r2.group.pathway   <- r2.group.pathway[ R2.group != 0 & intersection_size > 2]

## some numbers for the manuscript
length(unique(r2.group.pathway[ mapped_trait != ""]$R2.group))

## write to file
write.table(r2.group.pathway[ mapped_trait != "", .(R2.group, candidate.gene.classifier, cis.trans, protein.profile, pqtl.category, fc, p_value, term_id, term_name, intersection_size, mapped_trait)], 
            "Results.pQTL.protein.pathway.enrichment.GWAS.catalog.overlap.20250807.txt", sep = "\t", row.names = F)

##########################################################
####         trans-pQTLs and secreted proteins        ####
##########################################################

#-------------------------------------#
##--         HPA annotation        --##
#-------------------------------------#

## add to protein annotation table
prot.meta[, hpa.tissue.enhanced := sapply(Protein, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(max( hpa.tissue[ Gene %in% x]$tissue.spec))
})]

## add to protein annotation table
prot.meta[, hpa.cell.enhanced := sapply(Protein, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(max( hpa.cells[ Gene %in% x]$cell.spec))
})]

## import some more information
hpa.mapping    <- fread("/sc-projects/sc-proj-computational-medicine/data/07_public/01_Human_Protein_Atlas/data/hpa_all_2024-03-20.csv")

## add specifics - Secretome
prot.meta[, hpa.secretome := sapply(Protein, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(paste( hpa.mapping[ Gene %in% x]$`Secretome location`, collapse = " | "))
})]
## generate simplified variable for analysis
prot.meta[, hpa.secretome.binary := ifelse(hpa.secretome == "", 0, 1)]

## add specifics - `Subcellular location`
prot.meta[, hpa.subsellular.location := sapply(Protein, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(paste( hpa.mapping[ Gene %in% x]$`Subcellular location`, collapse = " | "))
})]

## add specifics - `Subcellular main location`
prot.meta[, hpa.subsellular.main.location := sapply(Protein, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(paste( hpa.mapping[ Gene %in% x]$`Subcellular main location`, collapse = " | "))
})]

## add specifics - `Blood concentration - Conc. blood IM [pg/L]`
prot.meta[, hpa.blood.conc.immuno := as.numeric(sapply(Protein, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(as.character(mean( hpa.mapping[ Gene %in% x]$`Blood concentration - Conc. blood IM [pg/L]`, na.rm=T)))
}))]

## add specifics - `Blood concentration - Conc. blood MS [pg/L]`
prot.meta[, hpa.blood.conc.mass.spec := as.numeric(sapply(Protein, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector: avoid issues with reporting large numbers
  return(as.character(mean( hpa.mapping[ Gene %in% x]$`Blood concentration - Conc. blood MS [pg/L]`, na.rm=T)))
}))]

#-------------------------------------#
##--         run analysis          --##
#-------------------------------------#

## parse some
prot.meta[, hpa.tissue.enhanced := ifelse(!is.finite(hpa.tissue.enhanced), NA, hpa.tissue.enhanced)]
prot.meta[, hpa.cell.enhanced := ifelse(!is.finite(hpa.cell.enhanced), NA, hpa.tissue.enhanced)]

## add to scallop results
res.scallop       <- merge(res.scallop, prot.meta[, .(MA_prot_id, hpa.secretome, hpa.secretome.binary, hpa.tissue.enhanced, hpa.cell.enhanced)],
                           by = "MA_prot_id")
## helper variables
res.scallop[, ind := 1]
res.scallop[, ind.R2 := paste0("R2.", R2.group)]
tmp               <- dcast(res.scallop, Protein ~ ind.R2, value.var = "ind")
## combine with protein information
tmp               <- merge(prot.meta, tmp, all.x = T)
## careful, this implies, that proteins with no genetic signal are set to missing

## run analysis by R2-group
res.secreted.pQTL <- lapply(unique(res.scallop$R2.group), function(x){
  
  ## get all entries
  jj <- res.scallop[ R2.group == x]
  
  ## progress
  print(jj)
  
  ## proceed only if at least five
  if(nrow(jj) > 4){
    
    ## run analysis
    ff <- summary(glm(paste0("hpa.secretome.binary ~ R2.", x), tmp, family = binomial(link = "logit")))$coefficients
    ## return results
    return(data.table(R2.group = x, beta=ff[2,1], se=ff[2,2], pval=ff[2,4]))
  }
})
## combine
res.secreted.pQTL <- rbindlist(res.secreted.pQTL)

##########################################################
####     heritability and protein characteristics     ####
##########################################################

#----------------------------#
##--    import results    --##
#----------------------------#

## import results generated by Mine
res.heritability <- fread("/<path to dir>/All_protein_targets_heritability.txt")

## simple plot
plot(sqrt(h2_polygenic_background) ~  I(sqrt(VE_cis_and_trans)), res.heritability)
## looks like there is a natural breakpoint, at which depencies changes
lw1 <- loess(sqrt(h2_polygenic_background) ~  I(sqrt(VE_cis_and_trans)), res.heritability,  span = 1)
j   <- order(res.heritability$VE_cis_and_trans)
lines(sqrt(res.heritability$VE_cis_and_trans)[j],lw1$fitted[j],col="red",lwd=3)
## use max of LOESS fit
lw1$x[which.max(lw1$fitted)]
## 0.34

## create new variables for changepoint detection
res.heritability[, sqrt.h2 := sqrt(h2_polygenic_background)]
res.heritability[, sqrt.var := sqrt(VE_cis_and_trans)]
## order for plotting
res.heritability <- res.heritability[order(sqrt.var)]

## create indicator: 0.34^2 = 0.1156
res.heritability[, poly.break := ifelse(sqrt.var >= 0.34, 1, 0)]

#----------------------------#
##-- protein mapping tab. --##
#----------------------------#

## import summary of proteins
tmp              <- fread("MA_protein_list_and_locations.txt")
## combine with heritability estimates
res.heritability <- merge(tmp, res.heritability, by.x = "OID_MA", by.y = "pheno")

#-------------------------------------#
##--         HPA annotation        --##
#-------------------------------------#

## add to protein annotation table
res.heritability[, hpa.tissue.enhanced := sapply(Assay, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(max( hpa.tissue[ Gene %in% x]$tissue.spec))
})]
## fix some
res.heritability[, hpa.tissue.enhanced := ifelse(!is.finite(hpa.tissue.enhanced), 0, 1)]

## add to protein annotation table
res.heritability[, hpa.cell.enhanced := sapply(Assay, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(max( hpa.cells[ Gene %in% x]$cell.spec))
})]
## fix some
res.heritability[, hpa.cell.enhanced := ifelse(!is.finite(hpa.cell.enhanced), 0, 1)]

## add specifics - Secretome
res.heritability[, hpa.secretome := sapply(Assay, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(paste( hpa.mapping[ Gene %in% x]$`Secretome location`, collapse = " | "))
})]
## generate simplified variable for analysis
res.heritability[, hpa.secretome.binary := ifelse(hpa.secretome == "", 0, 1)]

## add specifics - `Blood concentration - Conc. blood IM [pg/L]`
res.heritability[, hpa.blood.conc.immuno := as.numeric(sapply(Assay, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector
  return(as.character(mean( hpa.mapping[ Gene %in% x]$`Blood concentration - Conc. blood IM [pg/L]`, na.rm=T)))
}))]

## add specifics - `Blood concentration - Conc. blood MS [pg/L]`
res.heritability[, hpa.blood.conc.mass.spec := as.numeric(sapply(Assay, function(x){
  ## get all possible genes
  x <- strsplit(x, "_")[[1]]
  ## create a vector: avoid issues with reporting large numbers
  return(as.character(mean( hpa.mapping[ Gene %in% x]$`Blood concentration - Conc. blood MS [pg/L]`, na.rm=T)))
}))]

## simple analysis
hpa.heritability <- expand.grid(hpa.variable = c("hpa.secretome.binary", "hpa.tissue.enhanced", "hpa.cell.enhanced"), 
                                h2.variable = c("h2_polygenic_background", "VE_cis", "VE_trans", "VE_cis_and_trans", "h2_and_VE_all"), 
                                stringsAsFactors = F)
## simple testing
hpa.heritability <- lapply(1:nrow(hpa.heritability), function(x){
  
  ## get the relevant varibales
  hpa <- hpa.heritability$hpa.variable[x]
  h2  <- hpa.heritability$h2.variable[x]
  
  ## run wilcoxon rank sum test
  wt  <- wilcox.test(as.formula(paste(h2, hpa, sep = "~")), res.heritability)
  
  ## return values
  return(data.table(hpa.heritability[x,], p.value = wt$p.value, 
                    med.0 = median(unlist(res.heritability[ eval(as.name(hpa)) == 0, ..h2]), na.rm = T),
                    med.1 = median(unlist(res.heritability[ eval(as.name(hpa)) == 1, ..h2]), na.rm = T)))
  
  
})
## combine
hpa.heritability <- rbindlist(hpa.heritability)

#---------------------------------------#
##--  analysis break point proteins  --##
#---------------------------------------#

## --> pathway analysis <-- ##

## run pathway-based analysis
path.heritability <- gost(query = unique(unlist(lapply(res.heritability[ poly.break == 1]$Assay, function(x) strsplit(x, "_")[[1]]))), 
                          organism = "hsapiens", ordered_query = FALSE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                          measure_underrepresentation = FALSE, evcodes = TRUE, 
                          custom_bg = unique(unlist(lapply(res.heritability$Assay, function(x) strsplit(x, "_")[[1]]))),
                          user_threshold = 0.05, correction_method = "fdr", 
                          domain_scope = "annotated",
                          numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
## only one term
path.heritability <- as.data.table(path.heritability$result)
path.heritability[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]

## --> tissue analysis <-- ##

## run pathway-based analysis
tissue.heritability <- lapply(tissues$tissue, function(k){
  
  ## define differential proteins
  prot.diff <- unique(unlist(lapply(res.heritability[ poly.break == 1]$Assay, function(x) strsplit(x, "_")[[1]])))
  
  ## selected and tissue specific
  d1    <- nrow(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  ## not selected and tissue specific
  d2    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  ## selected and not tissue specific
  d3    <- nrow(hpa.tissue[ Gene %in% prot.diff & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  ## not selected and not tissue specific
  d4    <- nrow(hpa.tissue[ !(Gene %in% prot.diff) & !(`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  
  ## test for enrichment
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  
  ## return information needed
  return(data.table(tissue=k, or=enr$estimate, pval=enr$p.value, 
                    intersection=paste(hpa.tissue[ Gene %in% prot.diff & (`RNA tissue specificity` %in% c("Tissue enhanced", "Tissue enriched") & !is.na(eval(as.name(k)))) & olink.included == 1]$Gene, collapse = "|"), 
                    d1=d1, d2=d2, d3=d3, d4=d4))
})
tissue.heritability <- rbindlist(tissue.heritability, fill=T)

## --> cell-type enrichment <-- ##

cell.heritability   <- lapply(cell.type$cell.type, function(k){
  
  ## define differential proteins
  prot.diff <- unique(unlist(lapply(res.heritability[ poly.break == 1]$Assay, function(x) strsplit(x, "_")[[1]])))
  
  ## selected and tissue specific
  d1    <- nrow(hpa.cells[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  ## not selected and tissue specific
  d2    <- nrow(hpa.cells[ !(Gene %in% prot.diff) & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  ## selected and not tissue specific
  d3    <- nrow(hpa.cells[ Gene %in% prot.diff & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  ## not selected and not tissue specific
  d4    <- nrow(hpa.cells[ !(Gene %in% prot.diff) & !(`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k)))) & olink.included == 1])
  
  ## test for enrichment
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  
  ## return information needed
  return(data.table(cell.type=k, or=enr$estimate, pval=enr$p.value, 
                    intersection=paste(hpa.cells[ Gene %in% prot.diff & (`RNA single cell type specificity` %in% c("Cell type enhanced", "Cell type enriched") & !is.na(eval(as.name(k))))]$Gene, collapse = "|"), 
                    d1=d1, d2=d2, d3=d3, d4=d4))
})
cell.heritability   <- rbindlist(cell.heritability, fill=T)

##########################################################
####                   updated figure                 ####
##########################################################

## create column to do so
res.gwas.pruned[, effector.gene.category := ifelse(ligand.receptor.gene != "", "LR",
                                                   ifelse(protein.complex.gene != "", "PC",
                                                          ifelse(pathway.gene != "", "PW", "none")))]
## define levels
res.gwas.pruned[, effector.gene.category := factor(effector.gene.category, levels = c("LR", "PC", "PW", "none"))]

pdf("../graphics/Revised.FigXX.summary.trans.pQTL.pleiotropy.20250916.pdf", width = 6.3, height = 6.3)
## graphical parameters
par(mar=c(1.5,1.5,3.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, mgp=c(.6,0,0), xaxs="i", bty = "o", yaxs = "i")
## define layout to place opposing findings
layout(matrix(c(1,1,1,2:4,5,5,5), 3, 3, byrow = T), heights = c(1/3,1/3,1/3), widths = c(1/3,1/3,1/3))

#----------------------------------#
##--      Miami-like plot       --##
#----------------------------------#

## empty plot
plot(c(.5, chr.dat$plt.end[nrow(chr.dat)]+.5), c(-sqrt(600), sqrt(600)), type="n", xlab="Chromosomal position", ylab="Association count", xaxt="n", yaxt="n")
## add rectangles to divide chromosomes
pm <- par("usr")
rect(chr.dat$plt.start-.5, pm[3], chr.dat$plt.end+.5, pm[4], border = NA, col=c("white", "grey90"))
## add null
abline(h=0, lwd=.5)
## add axis
axis(2, at=c(sqrt(c(500,250,100,50))*-1, sqrt(c(0, 50, 100, 250, 500))), labels = c(500,250,100,50, 0, 50, 100, 250, 500), lwd=.5)
axis(1, at=chr.dat$plt.mid, labels = chr.dat$chromosome, lwd=.5)
## add estimates for proteins
points(res.gwas.pruned$plt.srt, sqrt(res.gwas.pruned$R2.count), type = "h", lwd=.3, col=ifelse(res.gwas.pruned$cis.trans == "cis", "#F95700", "#00A4CC"))
## add estiamtes for GWAS catalog traits
points(res.gwas.pruned$plt.srt, -sqrt(res.gwas.pruned$num_reported), type = "h", lwd=.3, col=ifelse(res.gwas.pruned$cis.trans == "cis", "#F95700", "#00A4CC"))
## add labels
legend("topleft", bty="n", lty=0, pch=NA, cex=.5, legend = "Protein targets")
legend("bottomleft", bty="n", lty=0, pch=NA, cex=.5, legend = "GWAS catalog")
legend("topright", bty="n", lty=1, lwd=.5, cex=.5, col=c("#F95700", "#00A4CC"), legend = c("cis", "trans"))

## add cis-protein genes with more than 20 associated EFO parent terms
foo          <- res.gwas.pruned[ cis.trans == "cis" & num.efo.parent >= 10]
## make unique by gene name
foo          <- foo[ order(candidate.gene.classifier, -num.efo.parent)]
foo[, ind := 1:.N, by="candidate.gene.classifier"]
foo          <- foo[ind == 1]
## order
foo          <- foo[order(plt.srt)]
## define plotting coordinates
pm  <- par("usr")

## define the coordinates of plotting the boxes
p.ii <- foo$plt.srt
## now get the distance between each point
d.ii <- p.ii[-1] - p.ii[-length(p.ii)]
## width of one box
w    <- strwidth("B")*1.2
## get the gaps
g.ii <- c(0,which(d.ii > w), length(p.ii))
for(j in 2:length(g.ii)){
  ## how large is the gap
  k  <- g.ii[j] - g.ii[j-1]
  if(k == 1){
    ## now set the labels
    text(p.ii[g.ii[j]], pm[4]+strheight("B")*.5, labels = foo$candidate.gene.classifier[g.ii[j]], cex=.4, srt=90, xpd=NA,
         pos=4, offset = 0, font=3, col = "#F95700")
    ## add arrow
    arrows(p.ii[g.ii[j]], pm[4]+strheight("B")*.5, p.ii[g.ii[j]],
           pm[4], lwd=.3, length=0, xpd=NA)
    # arrows(p.ii[g.ii[j]], pm[4], p.ii[g.ii[j]],
    #        foo$id[g.ii[j]]+2, lwd=.3, length=0, xpd=NA, lty=2)
  }else{
    ## define mean of the points in the gap
    m  <- mean(p.ii[(g.ii[j-1]+1):g.ii[j]])
    ## length of the whole block
    bw <- k*strwidth("B")*.5
    
    print(foo$candidate.gene[(g.ii[j-1]+1):g.ii[j]])
    print(p.ii[(g.ii[j-1]+1):g.ii[j]])
    print(m)
    print(seq(m-(bw/2), m+(bw/2), length.out=k))
    
    ## now set the labels
    text(seq(m-(bw/2), m+(bw/2), length.out=k), pm[4]+strheight("B")*.5, 
         labels = foo$candidate.gene.classifier[(g.ii[j-1]+1):g.ii[j]], cex=.4, srt=90, xpd=NA,
         pos=4, offset = 0, font=3, col = "#F95700")
    ## add arrows
    arrows(seq(m-(bw/2), m+(bw/2), length.out=k), pm[4]+strheight("B")*.5, p.ii[(g.ii[j-1]+1):g.ii[j]],
           pm[4], lwd=.3, length=0, xpd=NA)
    # arrows(p.ii[(g.ii[j-1]+1):g.ii[j]], pm[4], p.ii[(g.ii[j-1]+1):g.ii[j]],
    #        foo$id[(g.ii[j-1]+1):g.ii[j]]+2, lwd=.3, length=0, xpd=NA, lty=2)
  }
}

#----------------------------------#
##--   Opposing counts - MAP    --##
#----------------------------------#

## graphical parameters
par(mar=c(1.5,1.5,.5,.5))

## empty plot
plot(sqrt(c(0,600)), sqrt(c(0, 600)), type="n", ylab="Count GWAS catalog - mapped trait", xlab="Count protein targets", xaxt="n", yaxt="n", 
     lwd=.5, cex=.4, pch=20)
## add axis
axis(1, at=sqrt(c(0, 10, 50, 100, 250, 500)), labels = c(0, 10, 50, 100, 250, 500), lwd=.5)
axis(2, at=sqrt(c(0, 10, 50, 100, 250, 500)), labels = c(0, 10, 50, 100, 250, 500), lwd=.5)
## add points
points(sqrt(res.gwas.pruned[ R2.group != 0]$R2.count), sqrt(res.gwas.pruned[ R2.group != 0]$num_reported), cex=.3, pch=20, 
       col=ifelse(res.gwas.pruned$cis.trans == "cis", "#F95700", "#00A4CC"))
## annotate some extremes
ii <- which((res.gwas.pruned$R2.count > 150 | res.gwas.pruned$num_reported > 200) & res.gwas.pruned$R2.group != 0)
text(sqrt(res.gwas.pruned$R2.count[ii]), sqrt(res.gwas.pruned$num_reported[ii]), cex=.4, font=3, labels = res.gwas.pruned$candidate.gene.classifier[ii],
     xpd=NA)

#----------------------------------#
##--       pQTL categories      --##
#----------------------------------#

## adopt some parameters
par(mar=c(6,1.5,.5,.5))

## prepare the data for plotting
tmp <- table(res.gwas.pruned[pqtl.category != "specific"]$pqtl.category, 
             res.gwas.pruned[pqtl.category != "specific"]$cis.trans)
## reorder
tmp <- tmp[ c("unspecific.pleiotropy", "systemic.pleiotropy", "protein.pleiotropy"),]

## empty plot
plot(c(.5,3.5), c(0,1200), type="n", xlab="", ylab="Number of loci", xaxt="n", yaxt="n")
## add axis
axis(1, lwd=.5, at=1:3, labels=rownames(tmp), las=2, hadj = 1.2); axis(2, lwd=.5)
## add plots
rect(1:3-.4, 0, 1:3+.4, rowSums(tmp), border="white", lwd=.4, col="#00A4CC")
rect(1:3-.4, 0, 1:3+.4, tmp[1,], border="white", lwd=.4, col="#F95700")

#----------------------------------#
##--    causal gene fraction    --##
#----------------------------------#

## graphical parameters
par(mar=c(1.5,.5,.5,.5), bty = "n")

## empty plot
plot(c(0,1), c(.4,4.5), type = "n", xlab=" Fraction [%]", ylab = "", xaxt = "n", yaxt = "n")
## add axis
axis(1, at =c(0,.25,.5,.75,1), labels = c(0,25,50,75,100), lwd=.5)
## plotting coordinates
pm <- par("usr")

## add fractions: specific QTLs; White (#F1F3FFFF), Light Pink (#F7CED7FF), Pink (#F99FC9FF) and Dark Pink (#EF6079FF); Deep Blue (#2460A7FF), Northern Sky (#85B3D1FF), Baby Blue (#B3C7D6FF) and Coffee (#D9B48FFF)
tmp <- table(res.gwas.pruned[ cis.trans == "trans" & pqtl.category == "specific"]$effector.gene.category)
rect(c(0, cumsum(tmp)[1:3])/sum(tmp), .5, cumsum(tmp)/sum(tmp), 1, lwd=.3, border = "white", col = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))
## add legend
legend(pm[1], 1.4, bty = "n", cex=.6, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 4,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))

## add fractions: protein.pleiotropy QTLs; White (#F1F3FFFF), Light Pink (#F7CED7FF), Pink (#F99FC9FF) and Dark Pink (#EF6079FF)
tmp <- table(res.gwas.pruned[ cis.trans == "trans" & pqtl.category == "protein.pleiotropy"]$effector.gene.category)
rect(c(0, cumsum(tmp)[1:3])/sum(tmp), 1.5, cumsum(tmp)/sum(tmp), 2, lwd=.3, border = "white", col = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))
## add legend
legend(pm[1], 2.4, bty = "n", cex=.6, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 4,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))

## add fractions: protein.pleiotropy QTLs; White (#F1F3FFFF), Light Pink (#F7CED7FF), Pink (#F99FC9FF) and Dark Pink (#EF6079FF)
tmp <- table(res.gwas.pruned[ cis.trans == "trans" & pqtl.category == "systemic.pleiotropy"]$effector.gene.category)
rect(c(0, cumsum(tmp)[1:3])/sum(tmp), 2.5, cumsum(tmp)/sum(tmp), 3, lwd=.3, border = "white", col = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))
## add legend
legend(pm[1], 3.4, bty = "n", cex=.6, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 4,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))

## add fractions: protein.pleiotropy QTLs; White (#F1F3FFFF), Light Pink (#F7CED7FF), Pink (#F99FC9FF) and Dark Pink (#EF6079FF)
tmp <- table(res.gwas.pruned[ cis.trans == "trans" & pqtl.category == "unspecific.pleiotropy"]$effector.gene.category)
rect(c(0, cumsum(tmp)[1:3])/sum(tmp), 3.5, cumsum(tmp)/sum(tmp), 4, lwd=.3, border = "white", col = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))
## add legend
legend(pm[1], 4.4, bty = "n", cex=.6, pch = 22, pt.cex = .8, pt.lwd = .3, xpd = NA, ncol = 4,
       legend = paste0(names(tmp), " (n=",tmp, ")"),
       pt.bg = c("#2460A7FF", "#85B3D1FF", "#B3C7D6FF", "#D9B48FFF"))


#----------------------------------#
##--       Enrichment bubble    --##
#----------------------------------#

par(mar=c(1.5,1.5,.5,.5), bty = "o")

## empty plot
plot(c(.5, chr.dat$plt.end[nrow(chr.dat)]+.5), c(0,25), type="n", xlab="Chromosomal position", ylab=expression(-log[10]("FDR")), xaxt="n", yaxt="n")
## add rectangles to divide chromosomes
pm <- par("usr")
rect(chr.dat$plt.start-.5, pm[3], chr.dat$plt.end+.5, pm[4], border = NA, col=c("white", "grey90"))
## add axis
axis(2, lwd=.5)
axis(1, at=chr.dat$plt.mid, labels = chr.dat$chromosome, lwd=.5)
## add estimates for proteins
points(r2.group.pathway$plt.srt, -log10(r2.group.pathway$p_value),  lwd=.3, bg=adjustcolor(ifelse(r2.group.pathway$cis.trans == "cis", "#F95700", "#00A4CC"), .4),
       cex = log10(r2.group.pathway$fc), pch = 21)

## add legend
legend("topright", bty = "n", cex = .7, lty = 0, pch = 22, 
       pt.lwd = .3, pt.bg = c("#F95700", "#00A4CC"),
       legend = c("cis-pQTL", "trans-pQTL"))

## legend for Odds ratio
legend(pm[1]+(pm[2]-pm[1])*.2, pm[4], bty = "n", lty=0, cex = .7, pch = 21,
       pt.cex = log10(c(5,10,50,100)),
       title = "Fold enrichment", legend = c(5,10,50,100),
       ncol=4, pt.bg = rgb(0,0,0,.2))

## close device
dev.off()

#-------------------------------#
##--   results heritability  --##
#-------------------------------#

## create PNG for the Supplement
png("../graphics/Supplement.Heritability.Break.point.20250915.png", res = 300 , width = 8, height = 6, units = "cm")
## define plotting margins
par(mar=c(1.5,1.5,.5,.5), tck = -.01, mgp = c(.6,0,0), cex.axis = .5, cex.lab = .5, lwd=.5)

## oppose polygenic heritability and major loci
plot(sqrt.h2 ~ sqrt.var, res.heritability, cex = .5, pch = 21, bg = "grey80", lwd=.3,
     xlab = "Explained variance pQTL [%]",
     ylab = "Polygenic heritability (h2) excluding pQTL",
     xaxt = "n", yaxt = "n")
## axis
axis(1, at=sqrt(c(0,1,5,10,20,50)/100), labels = c(0,1,5,10,20,50), lwd = .5)
axis(2, at=sqrt(c(0,1,5,10,20,50)/100), labels = c(0,1,5,10,20,50)/100, lwd = .5)
## add LOESS fit
lw1 <- loess(sqrt.h2 ~ sqrt.var, res.heritability,  span = 1)
j   <- order(res.heritability$sqrt.var)
lines(res.heritability$sqrt.var[j],lw1$fitted[j],col="red",lwd=3)
## add breakpoint
abline(v=0.34, lwd=1, lty = 2)

## close device
dev.off()
