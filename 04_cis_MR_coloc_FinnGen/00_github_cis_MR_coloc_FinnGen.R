###################################################
#### perform FinnGen-wide MR for SCALLOP pQTL  ####
#### Maik Pietzner                  28/11/2023 ####
###################################################

rm(list=ls())
setwd("<path to dir>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(doMC)
require(readxl)
require(metafor)
require(colorspace)
require(biomaRt)
require(MendelianRandomization)

#############################################
####        import pQTLs to be run       ####
#############################################

## import summary of pQTLs
prot.meta       <- fread("<path to dir>/Protein.meta.SCALLOP.cis.trans.MR.pipeline.txt")

## import FinnGen data to be incorporated
def.finn        <- fread("<path to dir>/Summary.FinnGen.Release8.Phenotypes.20231129.txt")
def.finn        <- def.finn[ number.regional.sentinels > 0]

## import all fine-mapping results
res.fine        <- fread("<path to dir>/SCALLOP.pQTLs.MR.QCed.20231123.txt")

## add data on GWAS catalog associations to prune MR
res.gwas.pruned <- fread("<path to dir>/SCALLOP.R2.groups.GWAS.catalog.pleiotropy.20240613.txt")
## add the relevant part only
res.fine        <- merge(res.fine, res.gwas.pruned[, .(R2.group, num.efo.parent)])
## add information to protein meta data
prot.meta[, MarkerName.Trans.non.pleio := sapply(MA_prot_id, function(x) nrow(res.fine[ MA_prot_id == x & num.efo.parent <= 5]))]
## create another category containing only SNP with absence of pleiotropy
prot.meta[, MarkerName.Trans.non.pleio.both := sapply(MA_prot_id, function(x) nrow(res.fine[ MA_prot_id == x & num.efo.parent <= 5 & count <= 5]))]

## include flag for PAVs
scallop.pav     <- fread("SCALLOP_UKBB_MA_PAV_mapping_20240215.txt")
res.fine        <- merge(res.fine, scallop.pav[, .(pheno, MarkerName, PAV.lead.variant, PAV.in.LD)], all.x = T, by=c("pheno", "MarkerName"))
res.fine[, cis.pav.flag := ifelse(cis_trans == "Cis" & (PAV.lead.variant == 1 | PAV.in.LD == 1), T, F)]

## drop low freq variants
res.fine[, low.freq.variant := ifelse(Freq1_UKBB_SCALLOP_MA < .01 | Freq1_UKBB_SCALLOP_MA > .99, T, F)]
res.fine        <- res.fine[ low.freq.variant != T ]

## write updated list to file
write.table(res.fine, "<path to dir>/SCALLOP.pQTLs.MR.QCed.updated.20240125.txt", sep="\t", row.names = F)

## define number of elgible SNPs for each subset (consider only lead signals)
prot.meta[, MarkerName.Cis := sapply(MA_prot_id, function(x) nrow(res.fine[ MA_prot_id == x & cis_trans == "Cis"]))]
prot.meta[, MarkerName.Cis.no.pav := sapply(MA_prot_id, function(x) nrow(res.fine[ MA_prot_id == x & cis_trans == "Cis" & cis.pav.flag == F]))]
prot.meta[, MarkerName.Trans := sapply(MA_prot_id, function(x) nrow(res.fine[ MA_prot_id == x & cis_trans == "Trans" & signal.group == 1]))]
prot.meta[, MarkerName.Trans.specific := sapply(MA_prot_id, function(x) nrow(res.fine[ MA_prot_id == x & cis_trans == "Trans" & signal.group == 1 & count <= 5]))]
prot.meta[, MarkerName.Trans.non.pleio.both := sapply(MA_prot_id, function(x) nrow(res.fine[ MA_prot_id == x & num.efo.parent <= 5 & count <= 5 & signal.group == 1 & cis_trans == "Trans"]))]

#############################################
####          perform cis-coloc          ####
#############################################

## import Olink targets with at least one credible set
olink.targets <- fread("<path to dir>/Olink.targets.with.cis.credible.set.20240129.txt")

## store input file for the pipeline
write.table(olink.targets[, c("OID_MA", "chr", "region_start", "region_end")], "Olink.coloc.targets.txt", sep="\t", col.names = F, row.names = F, quote = F)

#----------------------------#
##--    import results    --##
#----------------------------#

## get all the files produced
ii        <- dir("../output_coloc/")

## import results
res.coloc <- lapply(ii, function(x){ fread(paste0("../output_coloc/", x))})
## combine
res.coloc <- rbindlist(res.coloc, fill = T)
## drop missing values
res.coloc <- res.coloc[ !is.na(PP.H0.abf)]
## add chromosome name
res.coloc[, chromosome_name := sub("^([0-9]+):.*", "\\1", MarkerName) ]
## add labels for proteins (duplicates some findings; allow for multiple cis mappings)
res.coloc <- merge(res.coloc, unique(olink.targets[, c("OID_MA", "Assay", "UniProt", "chromosome_name", "Panel.x")]), by.x=c("olink", "chromosome_name"), by.y=c("OID_MA", "chromosome_name"), all.x=T)
## add FinnGen explanation
res.coloc <- merge(res.coloc, def.finn[, c("processed_data_name", "LONGNAME", "NAME")], by.x = "finn.id", by.y = "processed_data_name")
## split by type of coloc and fuse with MR results
res.naive <- res.coloc[ type == "naive" ]
res.susie <- res.coloc[ type == "susie"]

## do some cleaning (top protein signal should be conserved)
res.naive <- res.naive[ R2.1 >= .8 ]
## same for the conditional coloc
res.susie[, R2.max := pmax(R2.1, R2.2, R2.3, R2.4, R2.4, R2.5, R2.6, R2.7, R2.7, R2.8, R2.9, R2.10, na.rm=T)]
res.susie <- res.susie[ R2.max >= .8 ]
## second line filter for potentially artificial results
res.susie <- res.susie[ as.numeric(Pvalue.protein) < 1e-6 ]

## indicate whether significant finding
res.susie[ PP.H4.abf >= .8 & ld.top < .5]

################################################################################
####                        coloc-first approach                            ####
################################################################################

## how many robust cis.coloc findings
res.naive.robust <- res.naive[ ld.top >= .8 & PP.H4.abf >= .8]

## same for credible set coloc
res.susie.robust <- res.susie[ ld.top >= .8 & PP.H4.abf >= .8]
## n = 343 (includes some duplications; take only strongest variant forward)
res.susie.robust <- res.susie.robust[ order(olink, finn.id, snp.id.protein, -abs(Effect.finn/StdErr.finn))]
res.susie.robust[, ind := 1:.N, by=c("olink", "finn.id", "snp.id.protein")]
res.susie.robust <- res.susie.robust[ ind == 1]

## add cis-mr estimates based on single variant
res.naive.robust <- lapply(1:nrow(res.naive.robust), function(x){
  ## prepare
  res <- mr_input(bx = res.naive.robust$Effect.protein[x], bxse = res.naive.robust$StdErr.protein[x],
                  by = res.naive.robust$Effect.finn[x], byse = res.naive.robust$StdErr.finn[x])
  ## run
  res <- mr_ivw(res)
  ## return results
  return(data.table(res.naive.robust[x,], beta.IVW = res@Estimate, se.IVW = res@StdError, pval.IVW = res@Pvalue, cil.IVW = res@CILower, ciu.IVW = res@CIUpper))
})
## combine again
res.naive.robust <- rbindlist(res.naive.robust)

## add cis-mr estimates based on single variant
res.susie.robust <- lapply(1:nrow(res.susie.robust), function(x){
  ## prepare
  res <- mr_input(bx = res.susie.robust$Effect.protein[x], bxse = res.susie.robust$StdErr.protein[x],
                  by = res.susie.robust$Effect.finn[x], byse = res.susie.robust$StdErr.finn[x])
  ## run
  res <- mr_ivw(res)
  ## return results
  return(data.table(res.susie.robust[x,], beta.IVW = res@Estimate, se.IVW = res@StdError, pval.IVW = res@Pvalue, cil.IVW = res@CILower, ciu.IVW = res@CIUpper))
})
## combine again
res.susie.robust <- rbindlist(res.susie.robust)

################################################################################
####                compare cis-coloc with cis-MR discovery                 ####
################################################################################

## get all cis results into one
jj            <- dir("../output/")

## subset to cis-findings only
registerDoMC(10)
res.cis <- mclapply(jj, function(x){ tmp <- fread(paste0("../output/", x))}, mc.cores=10) 
## combine everything
res.cis <- rbindlist(res.cis, fill=T)
res.cis <- as.data.table(res.cis)
## subset
res.cis <- res.cis[ type %in% c("cis", "cis_pav")]
## restrict to findings with fine-mapping support and targets outside of the extended MHC region 
res.cis <- res.cis[ MA_prot_id %in% olink.targets[ !(chr == 6 & region_start >= 25500000 & region_end <= 34000000)]$OID_MA ]
## add study description
res.cis <- merge(res.cis, def.finn[, c("processed_data_name", "LONGNAME", "NAME")], by.x = "finn.id", by.y = "processed_data_name")
## impose QC flags for MR
res.cis[, qc.hetero := is.na(pval.q2) | pval.q2 > 1e-4]
res.cis[, qc.pleio  := is.na(pval.intercept) | pval.intercept > 1e-3]
## directional concordant
res.cis[, directional.concordant := sign(beta.IVW) + sign(beta.Simple.median) + sign(beta.MR.Egger) + sign(beta.Weighted.median)]
res.cis[, directional.concordant := is.na(directional.concordant) | directional.concordant %in% c(-4,4)]
## significance IVW: 6.599131e-08
res.cis[, sig.ivw := pval.IVW < .05/(nrow(res.cis)/2)]
## reshape by type
res.cis <- tidyr::pivot_wider(res.cis, id_cols = c("MA_prot_id", "finn.id"), names_from = "type", 
                              values_from = c(grep("IVW", names(res.cis), value=T), "qc.pleio", "qc.hetero", "directional.concordant", "sig.ivw", "nsnps.all", "nsnps.loo", "outlying.snps"), names_sep = ".")
res.cis <- as.data.table(res.cis)
## merge the number of SNPs
res.cis <- merge(res.cis, prot.meta)
## add ac flag for retained SNPs
res.cis[, snp.conserved.cis := ifelse(nsnps.all.cis/MarkerName.Cis >= 2/3, T, F)]
res.cis[, snp.conserved.cis_pav := ifelse(nsnps.all.cis_pav/MarkerName.Cis.no.pav >= 2/3, T, F)]
## indicator high-confidence results
res.cis[, high.confidence.cis := (qc.hetero.cis + qc.pleio.cis + directional.concordant.cis + sig.ivw.cis + snp.conserved.cis) == 5]
res.cis[, high.confidence.cis_pav := (qc.hetero.cis_pav + qc.pleio.cis_pav + directional.concordant.cis_pav + sig.ivw.cis_pav + snp.conserved.cis_pav) == 5]
## add study description
res.cis <- merge(res.cis, def.finn[, c("processed_data_name", "LONGNAME", "NAME")], by.x = "finn.id", by.y = "processed_data_name")

## --> combine with cis-coloc findings <-- ##

## create intermediate SuSiE file d
tmp.susie        <- res.susie[ order( olink, finn.id, -PP.H4.abf)]
tmp.susie[, ind := 1:.N, by=c("olink", "finn.id")]
tmp.susie        <- tmp.susie[ ind == 1]
tmp.susie        <- tmp.susie[ as.numeric(Pvalue.finn) < 1e-5 ]
## combine
res.cis.combined <- merge(res.cis, tmp.susie, by.x=c("finn.id", "MA_prot_id", "LONGNAME", "NAME"), by.y=c("finn.id", "olink", "LONGNAME", "NAME"), 
                          all = T, suffixes = c(".MR", ".susie"))
## some duplications..
res.cis.combined <- res.cis.combined[ !is.na(pval.IVW.cis)]
## careful some multiplications
res.cis.combined[, id.key := paste(finn.id, MA_prot_id, sep = "$")]

## create indicator for cis-coloc
res.cis.combined[, cis.coloc := ifelse(is.na(PP.H4.abf) | PP.H4.abf < .8 | ld.top < .8, F, T)]

## plot cis concordance
tmp.mr           <- res.cis.combined[ high.confidence.cis == T | cis.coloc == T ]

## indicator what type of example
res.cis.combined[, type.example := ifelse(cis.coloc == T & high.confidence.cis == T, "coloc+MR", 
                                          ifelse(cis.coloc == T, "coloc",
                                                 ifelse(high.confidence.cis == T & nsnps.all.cis >= 2 & is.na(PP.H4.abf) == T, "MR", "none")))]

#-----------------------------------------#
##--           summary figure          --##
#-----------------------------------------#

## create colour vector
cat.col          <- data.table(Category = c("I Certain infectious and parasitic diseases (AB1_)",
                                            "II Neoplasms from hospital discharges (CD2_)",
                                            "II Neoplasms, from cancer register (ICD-O-3)",
                                            "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)",
                                            "IV Endocrine, nutritional and metabolic diseases (E4_)",
                                            "V Mental and behavioural disorders (F5_)",
                                            "VI Diseases of the nervous system (G6_)",
                                            "VII Diseases of the eye and adnexa (H7_)",
                                            "VIII Diseases of the ear and mastoid process (H8_)",
                                            "IX Diseases of the circulatory system (I9_)",
                                            "X Diseases of the respiratory system (J10_)",
                                            "XI Diseases of the digestive system (K11_)",
                                            "XII Diseases of the skin and subcutaneous tissue (L12_)",
                                            "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
                                            "XIV Diseases of the genitourinary system (N14_)",
                                            "XV Pregnancy, childbirth and the puerperium (O15_)",
                                            "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)",
                                            "XVIII Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R18_)",
                                            "XIX Injury, poisoning and certain other consequences of external causes (ST19_)", 
                                            "Composite endpoints"),
                               cl = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(20), cat.srt=1:20)

## add to results data frame
def.finn[, Category := ifelse(is.na(Category), "Composite endpoints", Category)]
def.finn         <- merge(def.finn, cat.col)
## create order
def.finn         <- def.finn[ order(cat.srt, NAME)]
def.finn[, finn.srt := 1:nrow(def.finn)]
res.cis.combined <- merge(res.cis.combined, def.finn[, .(processed_data_name, Category, cl, cat.srt, finn.srt)], 
                          by.x = "finn.id", by.y = "processed_data_name")

## --> add single point MR estimates wherever possible <-- ##

## do in parallel
registerDoMC(10)

## go through
tmp              <- mclapply(1:nrow(res.cis.combined), function(x){
  
  ## do only if coloc was performed
  if(!(is.na(res.cis.combined$Pvalue.protein[x]))){
    ## create MR input
    mr.input <- mr_input(bx = res.cis.combined$Effect.protein[x], bxse = res.cis.combined$StdErr.protein[x],
                         by = res.cis.combined$Effect.finn[x], byse = res.cis.combined$StdErr.finn[x])
    ## run the MR
    mr.res   <- mr_ivw(mr.input)
    ## return results
    return(data.table(res.cis.combined[x,], 
                      beta.mr.single=mr.res@Estimate, se.mr.single=mr.res@StdError, 
                      pval.mr.single=mr.res@Pvalue))
  }else{
    return(res.cis.combined[x,])
  }
  
}, mc.cores = 10)
## combine
tmp              <- rbindlist(tmp, fill = T)
## re-assign
res.cis.combined <- tmp

## create combined p-value for plotting
res.cis.combined[, log10p.plot := apply(res.cis.combined[, .(beta.IVW.cis, se.IVW.cis, beta.mr.single, se.mr.single)], 1, function(x){
  ## compute minimum z-score
  zscore <- min(c(x[1]/x[2], x[3]/x[4]), na.rm = T)
  ## return approximate log10p values
  return(-pchisq((zscore)^2, df=1, lower.tail=F, log.p=T)/log(10))
})]

## --> protein coordinates across the genome <-- ##

## order protein targets
olink.targets    <- olink.targets[ order(chr, start_position)]
## add order
olink.targets[, oid.srt := 1:nrow(olink.targets)]
## subset and make unique
tmp              <- unique(olink.targets[, .(OID_MA, oid.srt)])
tmp[, ind := 1:.N, by="OID_MA"]
tmp              <- tmp[ ind == 1]
## add to the summary results
res.cis.combined <- merge(res.cis.combined, tmp[, .(OID_MA, oid.srt)], by.x="MA_prot_id", by.y="OID_MA")

## --> prep the plot <-- ##

## create temp data for pseudo chromosomes
chr.dat        <- do.call(data.frame, aggregate(oid.srt ~ chr, olink.targets, function(x) c(min(x), mean(x), max(x)))) 
names(chr.dat) <- c("chr", "start", "mid", "end")

## same for finn IDs
tmp            <- do.call(data.frame, aggregate(finn.srt ~ cat.srt, def.finn, function(x) c(min(x), mean(x), max(x)))) 
names(tmp)     <- c("cat.srt", "start", "mid", "end")
## add to the colouring data set
cat.col        <- merge(cat.col, tmp)

## get only those of interest
tmp       <- res.cis.combined[ cis.coloc == T | (high.confidence.cis == T & nsnps.all.cis >= 2 & is.na(PP.H4.abf))]
## simple count
tmp.count <- tmp[, .(count = length(finn.id), oid.srt=unique(oid.srt), Protein=unique(Protein), Category=length(unique(Category))), by="MA_prot_id"]

## start graphical device
pdf("../graphics_transfer//Summary.cis.coloc.MR.FinnGen.20240628.pdf", width = 6.3, height = 4)
par(mar=c(.1,1.5,2.5,.5), mgp=c(.4,0,0), cex.axis=.4, cex.lab=.4, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
## define layout
layout(matrix(1:2, 2, 1), heights = c(.3,.7))

#---------------------------------#
##-- upper part for pleiotropy --##
#---------------------------------#

## empty plot
plot(c(chr.dat$start[1], chr.dat$end[23]), c(0,30), xaxt="n", yaxt="n",
     xlab="", ylab="Number\nassociated diseases", type="n")
## add axis
axis(2, lwd=.5)
## rectanlges to divide chromosomes
rect(chr.dat$start, 0, chr.dat$end, 50, col=c("white", "grey90"), border = NA)

## add regional estimates
arrows(tmp.count$oid.srt, 0, tmp.count$oid.srt, tmp.count$count,  length = 0, lwd=.5, col="grey50")

## add legend
pm <- par("usr")
legend(pm[1]+(pm[2]-pm[1])*.5, pm[4], bty="n", lty=0, pch=21, pt.bg="white", cex=.4, pt.lwd=.3,
       pt.cex=rev(1:4)/10, legend = rev(1:4), title = "ICD-10 chapters")

## add genes with at least 5 annotations
foo          <- tmp.count[ count > 1]
## make unique by gene name
foo          <- foo[ order(Protein, -count)]
foo[, ind := 1:.N, by="Protein"]
foo          <- foo[ind == 1]
## order
foo          <- foo[order(foo$oid.srt)]
## define plotting coordinates
pm  <- par("usr")

## define the coordinates of plotting the boxes
p.ii <- foo$oid.srt
## now get the distance between each point
d.ii <- p.ii[-1] - p.ii[-length(p.ii)]
## width of one box
w    <- strwidth("B")*1.7
## get the gaps
g.ii <- c(0,which(d.ii > w), length(p.ii))
for(j in 2:length(g.ii)){
  ## how large is the gap
  k  <- g.ii[j] - g.ii[j-1]
  if(k == 1){
    ## now set the labels
    text(p.ii[g.ii[j]], pm[4]+strheight("B")*.5, labels = foo$Protein[g.ii[j]], cex=.3, srt=90, xpd=NA,
         pos=4, offset = 0, font=3)
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
         labels = foo$Protein[(g.ii[j-1]+1):g.ii[j]], cex=.3, srt=90, xpd=NA,
         pos=4, offset = 0, font=3)
    ## add arrows
    arrows(seq(m-(bw/2), m+(bw/2), length.out=k), pm[4]+strheight("B")*.5, p.ii[(g.ii[j-1]+1):g.ii[j]],
           pm[4], lwd=.3, length=0, xpd=NA)
    # arrows(p.ii[(g.ii[j-1]+1):g.ii[j]], pm[4], p.ii[(g.ii[j-1]+1):g.ii[j]],
    #        foo$id[(g.ii[j-1]+1):g.ii[j]]+2, lwd=.3, length=0, xpd=NA, lty=2)
  }
}

## add categories on top
points(tmp.count$oid.srt, tmp.count$count, pch=21, cex=tmp.count$Category/10,
       bg="white", lwd=.3, xpd=NA)


#----------------------------------#
##-- lower part for single vars --##
#----------------------------------#

## adapt plotting parameters
par(mar=c(3.8,1.5,.1,.5))

## empty plot
plot(c(chr.dat$start[1], chr.dat$end[23]), c(0,max(cat.col$end)), xaxt="n", yaxt="n",
     xlab="Chromosomal position", ylab="", type="n", ylim=rev(c(0,max(cat.col$end))))
## add axis
axis(1, lwd=.5, at=chr.dat$mid, labels=1:23)
## get plotting coordinates
pm <- par("usr")
## add phenotype categories
rect(pm[1], cat.col$start, pm[2], cat.col$end, col=colorspace::lighten(cat.col$cl, .6), border=cat.col$cl, lwd=.1)
abline(v=chr.dat$start, lwd=.3, lty=2, col="white")
## add points
points(tmp$oid.srt, tmp$finn.srt, pch=20, col=ifelse(tmp$type.example == "coloc+MR", "black",
                                                     ifelse(tmp$type.example == "coloc", "grey30", "grey70")), 
       cex=.3, lwd=.1, xpd=NA)

## add legend
legend(pm[1]-(pm[2]-pm[1])*.02, pm[3]-(pm[4]-pm[3])*.1, lty=0, pch=22, pt.lwd=.1, bty="n", xpd=NA,
       pt.bg=c(cat.col$cl, "black", "grey30", "grey70"), legend = c(cat.col$Category, "coloc+MR", "coloc", "MR"),
       cex=.3, ncol=3, pt.cex=.5)
## close device
dev.off()

################################################################################
####                        cis+trans-MR discovery                          ####
################################################################################

## get all cis results into one
jj            <- dir("../output/")

## subset to cis-findings only
registerDoMC(10)
res.cis.trans <- mclapply(jj, function(x){ tmp <- fread(paste0("../output/", x))}, mc.cores=10) 
## combine everything
res.cis.trans <- rbindlist(res.cis.trans, fill=T)
res.cis.trans <- as.data.table(res.cis.trans)
## subset
res.cis.trans <- res.cis.trans[ type %in% c("cis_trans_spec", "cis_trans_pleio")]
## retain only proteins with at least two specific trans_pQTLs
res.cis.trans <- res.cis.trans[ MA_prot_id %in% prot.meta[MarkerName.Trans.specific >= 3 & MarkerName.Cis > 0 ]$MA_prot_id]
## fine-mapping results and outside of MHC region
res.cis.trans <- res.cis.trans[ MA_prot_id %in% olink.targets[ !(chr == 6 & region_start >= 25500000 & region_end <= 34000000)]$OID_MA]
## create QC measures
res.cis.trans[, qc.hetero := is.na(pval.q2) | pval.q2 > 1e-4]
res.cis.trans[, qc.pleio  := is.na(pval.intercept) | pval.intercept > 1e-3]
## directional concordant
res.cis.trans[, directional.concordant := sign(beta.IVW) + sign(beta.Simple.median) + sign(beta.MR.Egger) + sign(beta.Weighted.median)]
res.cis.trans[, directional.concordant := is.na(directional.concordant) | directional.concordant %in% c(-4,4)]
## significance IVW: 9.811463e-08
res.cis.trans[, sig.ivw := pval.IVW < .05/(nrow(res.cis.trans)/2)]
## reshape by type
res.cis.trans <- tidyr::pivot_wider(res.cis.trans, id_cols = c("MA_prot_id", "finn.id"), names_from = "type", 
                                    values_from = c(grep("IVW", names(res.cis.trans), value=T), "qc.pleio", "qc.hetero", "directional.concordant", "sig.ivw", "nsnps.all", "nsnps.loo", "outlying.snps"), names_sep = ".")
res.cis.trans <- as.data.table(res.cis.trans)
## merge the number of SNPs
res.cis.trans <- merge(res.cis.trans, prot.meta)
## add ac flag for retained SNPs
res.cis.trans[, snp.conserved.cis_trans_spec := ifelse(nsnps.all.cis_trans_spec/(MarkerName.Trans.specific + MarkerName.Cis) >= 2/3, T, F)]
res.cis.trans[, snp.conserved.cis_trans_pleio := ifelse(nsnps.all.cis_trans_pleio/(MarkerName.Trans.non.pleio.both + MarkerName.Cis) >= 2/3, T, F)]

## add study description
res.cis.trans <- merge(res.cis.trans, def.finn[, c("processed_data_name", "LONGNAME", "NAME")], by.x = "finn.id", by.y = "processed_data_name")
res.cis.trans <- as.data.table(res.cis.trans)

## declare high confidence findings
res.cis.trans[, high.confidence.cis_trans_spec := (sig.ivw.cis_trans_spec + qc.hetero.cis_trans_spec + qc.pleio.cis_trans_spec + directional.concordant.cis_trans_spec + snp.conserved.cis_trans_spec) == 5]
res.cis.trans[, high.confidence.cis_trans_pleio := (sig.ivw.cis_trans_pleio + qc.hetero.cis_trans_pleio + qc.pleio.cis_trans_pleio+ directional.concordant.cis_trans_pleio + snp.conserved.cis_trans_pleio) == 5]

## check whether there is support in cis
res.cis.trans <- merge(res.cis.trans, res.cis.combined[, c("MA_prot_id", "finn.id", grep("\\.cis", names(res.cis.combined), value=T), "cis.coloc"), with=F], all.x=T, by=c("MA_prot_id", "finn.id"))

################################################################################
####                          trans-MR discovery                            ####
################################################################################

## get all cis results into one
jj      <- dir("../output/")

## subset to cis-findings only
res.trans <- mclapply(jj, function(x){ tmp <- fread(paste0("../output/", x))}, mc.cores = 10) 
## combine everything
res.trans <- rbindlist(res.trans, fill=T)
res.trans <- as.data.table(res.trans)
## subset
res.trans <- res.trans[ type %in% c("trans_spec", "trans_pleio")]
## retain only proteins with at least two specific trans_pQTLs (n=549 OIDs)
res.trans <- res.trans[ MA_prot_id %in% prot.meta[MarkerName.Trans.specific >= 3]$MA_prot_id]
## create QC measures
res.trans[, qc.hetero := is.na(pval.q2) | pval.q2 > 1e-4]
res.trans[, qc.pleio  := is.na(pval.intercept) | pval.intercept > 1e-3]
## directional concordant
res.trans[, directional.concordant := sign(beta.IVW) + sign(beta.Simple.median) + sign(beta.MR.Egger) + sign(beta.Weighted.median)]
res.trans[, directional.concordant := is.na(directional.concordant) | directional.concordant %in% c(-4,4)]
## significance IVW: 1.436235e-07
res.trans[, sig.ivw := pval.IVW < .05/(nrow(res.trans)/2)]
## reshape by type
res.trans <- tidyr::pivot_wider(res.trans, id_cols = c("MA_prot_id", "finn.id"), names_from = "type", 
                                values_from = c(grep("IVW", names(res.trans), value=T), "qc.pleio", "qc.hetero", "directional.concordant", "sig.ivw", "nsnps.all", "nsnps.loo", "outlying.snps"), names_sep = ".")
res.trans <- as.data.table(res.trans)
## merge the number of SNPs
res.trans <- merge(res.trans, prot.meta)
## add ac flag for retained SNPs
res.trans[, snp.conserved.trans_spec := ifelse(nsnps.all.trans_spec/MarkerName.Trans.specific >= 2/3, T, F)]
res.trans[, snp.conserved.trans_pleio := ifelse(nsnps.all.trans_pleio/MarkerName.Trans.non.pleio.both >= 2/3, T, F)]

## add study description
res.trans <- merge(res.trans, def.finn[, c("processed_data_name", "LONGNAME", "NAME")], by.x = "finn.id", by.y = "processed_data_name")
res.trans <- as.data.table(res.trans)

## declare high confidence findings
res.trans[, high.confidence.trans_spec := (sig.ivw.trans_spec + qc.hetero.trans_spec + qc.pleio.trans_spec + directional.concordant.trans_spec + snp.conserved.trans_spec) == 5]
res.trans[, high.confidence.trans_pleio := (sig.ivw.trans_pleio + qc.hetero.trans_pleio + qc.pleio.trans_pleio + directional.concordant.trans_pleio + snp.conserved.trans_pleio) == 5 & nsnps.all.trans_pleio > 2]

## check whether there is support in cis
res.trans <- merge(res.trans, res.cis.combined[, c("MA_prot_id", "finn.id", grep("\\.cis", names(res.cis.combined), value=T), "cis.coloc"), with=F], all.x=T, by=c("MA_prot_id", "finn.id"))

################################################################################
####                      convergence cis/trans MR                          ####
################################################################################

## create combined data
res.combined <- merge(res.cis.combined, res.trans[, c("MA_prot_id", "finn.id", grep("trans_spec|trans_pleio", names(res.trans), value=T)), with=F], all.x=T, by=c("MA_prot_id", "finn.id"))

## subset to variants of interest
tmp.mr       <- res.combined[ (cis.coloc == T & pval.IVW.cis < .05 & !is.na(pval.IVW.trans_pleio)) | (high.confidence.cis == T & nsnps.all.cis >= 2 & is.na(PP.H4.abf) & !is.na(pval.IVW.trans_pleio)) ]

## define MR estimate to be used based on type of example
tmp.mr[, beta.IVW.cis_fused := ifelse( type.example == "coloc", beta.mr.single, beta.IVW.cis)]
tmp.mr[, se.IVW.cis_fused := ifelse( type.example == "coloc", se.mr.single, se.IVW.cis)]

## add heterogeneity estimates
tmp.mr       <- lapply(1:nrow(tmp.mr), function(x){
  
  ## compute MA across cis and specific trans
  ma.spec  <- rma(yi=c(tmp.mr$beta.IVW.cis_fused[x], tmp.mr$beta.IVW.trans_spec[x]),
                  sei=c(tmp.mr$se.IVW.cis_fused[x], tmp.mr$se.IVW.trans_spec[x]),
                  method="FE")
  
  ## compute MA across cis and specific trans
  ma.pleio  <- rma(yi=c(tmp.mr$beta.IVW.cis_fused[x], tmp.mr$beta.IVW.trans_pleio[x]),
                   sei=c(tmp.mr$se.IVW.cis_fused[x], tmp.mr$se.IVW.trans_pleio[x]),
                   method="FE")
  
  ## return results
  return(data.table(tmp.mr[x, ], 
                    heterogeneity_trans_spec=ma.spec$I2, pval.hetero_trans_spec=ma.spec$QEp,
                    heterogeneity_trans_pleio=ma.pleio$I2, pval.hetero_trans_pleio=ma.pleio$QEp))
  
})
## combine again
tmp.mr  <- rbindlist(tmp.mr)

##########################################################################################################
##########################################################################################################
####                                    REVISION CELL - 07/07/2025                                    ####
##########################################################################################################
##########################################################################################################

################################################
####       cis-trans MR findings - naive    ####
################################################

#-------------------------------------#
##-- findings cis+trans w/o filter --##
#-------------------------------------#

## get all cis results into one
jj            <- dir("../output/")

## subset to cis-findings only
registerDoMC(10)
res.cis.trans.naive <- mclapply(jj, function(x){ tmp <- fread(paste0("../output/", x))}, mc.cores=10) 
## combine everything
res.cis.trans.naive <- rbindlist(res.cis.trans.naive, fill=T)
res.cis.trans.naive <- as.data.table(res.cis.trans.naive)
## subset
res.cis.trans.naive <- res.cis.trans.naive[ type %in% c("cis", "cis_trans", "trans")]
## retain only proteins with at least five trans_pQTLs (n= OIDs)
res.cis.trans.naive <- res.cis.trans.naive[ MA_prot_id %in% prot.meta[MarkerName.Trans >= 5]$MA_prot_id]
## create QC measures
res.cis.trans.naive[, qc.hetero := is.na(pval.q2) | pval.q2 > 1e-4]
res.cis.trans.naive[, qc.pleio  := is.na(pval.intercept) | pval.intercept > 1e-3]
## directional concordant
res.cis.trans.naive[, directional.concordant := sign(beta.IVW) + sign(beta.Simple.median) + sign(beta.MR.Egger) + sign(beta.Weighted.median)]
res.cis.trans.naive[, directional.concordant := is.na(directional.concordant) | directional.concordant %in% c(-4,4)]
## significance IVW:  7.324852e-08
res.cis.trans.naive[, sig.ivw := pval.IVW < .05/(nrow(res.cis.trans.naive)/3)]

## reshape by type
res.cis.trans.naive <- tidyr::pivot_wider(res.cis.trans.naive, id_cols = c("MA_prot_id", "finn.id"), names_from = "type", 
                                          values_from = c(grep("IVW", names(res.cis.trans.naive), value=T), "qc.pleio", "qc.hetero", "directional.concordant", "sig.ivw", "nsnps.all", "nsnps.loo", "outlying.snps"), names_sep = ".")
res.cis.trans.naive <- as.data.table(res.cis.trans.naive)
## merge the number of SNPs
res.cis.trans.naive <- merge(res.cis.trans.naive, prot.meta)

## add study description
res.cis.trans.naive <- merge(res.cis.trans.naive, def.finn[, c("processed_data_name", "LONGNAME", "NAME")], by.x = "finn.id", by.y = "processed_data_name")
res.cis.trans.naive <- as.data.table(res.cis.trans.naive)

#--------------------------------------#
##-- benchmark against cis-examples --##
#--------------------------------------#

## create combined data
res.review          <- merge(res.cis.combined, res.cis.trans.naive[, c("MA_prot_id", "finn.id", grep("cis_trans|trans", names(res.cis.trans.naive), value=T)), with=F], all.x=T, by=c("MA_prot_id", "finn.id"))

## take forward all cis-examples and cis/trans
res.review          <- res.review[ (cis.coloc == T & pval.IVW.cis < .05 & !is.na(pval.IVW.trans)) | (high.confidence.cis == T & nsnps.all.cis >= 2 & is.na(PP.H4.abf) & !is.na(pval.IVW.trans)) | sig.ivw.cis_trans == T | sig.ivw.trans == T]

## define MR estimate to be used based on type of example
res.review[, beta.IVW.cis_fused := ifelse( type.example == "coloc", beta.mr.single, beta.IVW.cis)]
res.review[, se.IVW.cis_fused := ifelse( type.example == "coloc", se.mr.single, se.IVW.cis)]

## add heterogeneity estimates
res.review          <- lapply(1:nrow(res.review), function(x){
  
  ## compute MA across cis and cis/trans
  ma.cis.trans  <- rma(yi=c(res.review$beta.IVW.cis_fused[x], res.review$beta.IVW.cis_trans[x]),
                       sei=c(res.review$se.IVW.cis_fused[x], res.review$se.IVW.cis_trans[x]),
                       method="FE")
  
  ## compute MA across trans
  ma.trans     <- rma(yi=c(res.review$beta.IVW.cis_fused[x], res.review$beta.IVW.trans[x]),
                      sei=c(res.review$se.IVW.cis_fused[x], res.review$se.IVW.trans[x]),
                      method="FE")
  
  ## return results
  return(data.table(res.review[x, ], 
                    heterogeneity_cis_trans=ma.cis.trans$I2, pval.hetero_cis_trans=ma.cis.trans$QEp,
                    heterogeneity_trans=ma.trans$I2, pval.hetero_trans=ma.trans$QEp))
  
})
## combine again
res.review          <- rbindlist(res.review)
