####################################
#### PRE-PROCESSING OF DATASETS ####
####################################

## LOAD LIBRARIES

library(SummarizedExperiment)
library(GenomicRanges)
library(rexposome)
library(minfi)
library(stringr)
library(SNPassoc)
library(tidyverse)
library(arsenal) 


## LOAD FILES

#Methylome
load("/home/isglobal.lan/ncarreras/data/WS_HELIX/HELIX_preproc/methylation/Final_data/methylome_subcohort_residuals2beta_v4.RData")

#Transcriptome
load("/home/isglobal.lan/ncarreras/data/WS_HELIX/HELIX_preproc/gene_expression/Final_data/transcriptome_subcohort_f1_residuals3_v3.Rdata")

#Exposome
load("/home/isglobal.lan/ncarreras/data/WS_HELIX/HELIX_preproc/exposome/FinalDataset/imppostnatal_v3.Rdata")
load("/home/isglobal.lan/ncarreras/data/WS_HELIX/HELIX_preproc/exposome/FinalDataset/imppregnancy_v3.Rdata")

#Inversions (scClassDF)
load("/home/isglobal.lan/ncarreras/data/WS_HELIX/HELIX_preproc/gwas/inversions/scoreInvHapClassDFHELIX.Rdata")


## SAMPLEID --> HELIXID

#Methylome IDs
methy_HelixID <- pData(methylome_subcohort_residuals2beta)[,1]
methylome_Helix <- methylome_subcohort_residuals2beta
colnames(methylome_Helix) <- methy_HelixID

#Transcriptome IDs
trans_HelixID <- pData(transcriptome_subcohort_f1_residuals3)[,2]
transcriptome_Helix <- transcriptome_subcohort_f1_residuals3
colnames(transcriptome_Helix) <- trans_HelixID

#Exposome data
#First, get only one imputation
imppreg <- toES(imppregnancy, rid = 1)
imppost <- toES(imppostnatal, rid = 1)

pregn_HelixID <- pData(imppreg)[,3]
sampleNames(imppreg) <- pregn_HelixID
post_HelixID <- pData(imppost)[,3]
sampleNames(imppost) <- post_HelixID


## ADD THE INVERSIONS IN THE pData AND DISCARD THE NON CAUCASIAN SAMPLES

#Convert scClassDF into a table
inversions <- as.matrix(scClassDF)

#Convert inversions to SNP data
inversions.s <- setupSNP(inversions, 1:ncol(inversions),sep="/")
rownames(inversions.s) <- rownames(inversions)

#Sort inversion genotypes as N/N, N/I, and I/I
for (inv in colnames(inversions.s)){
  inversions.s[[inv]] <-  factor(inversions.s[[inv]], levels=c("N/N", "N/I", "I/I"))
}

add_inv_cauc <- function(dataset){
  inv_ID <- intersect(colnames(dataset),rownames(inversions.s))
  final <- dataset[,inv_ID]
  pData(final) <- cbind(pData(final),inversions.s[inv_ID, ])
  if (class(dataset)!="ExposomeSet"){
    final <- final[, which(final$h_ethnicity_cauc=="yes") ] #Discard non caucasian samples
  }
  if (class(dataset)=="ExposomeSet"){
    final <- final[, which(final$h_ethnicity_c_None=="Caucasian") ] #Discard non caucasian samples
  }
  return(final)
}

methy_final <- add_inv_cauc(methylome_Helix)
trans_final <- add_inv_cauc(transcriptome_Helix)
imppreg_final <- add_inv_cauc(imppreg)
imppost_final <- add_inv_cauc(imppost)


## SEX AND COHORT COVARIATES RENAME

#e3_sex/e3_sex_None --> sex
colnames(pData(methy_final))[7] <- "sex"
colnames(pData(trans_final))[6] <- "sex"
colnames(pData(imppreg_final))[53] <- "sex"
colnames(pData(imppost_final))[53] <- "sex"

#h_cohort --> cohort
colnames(pData(imppreg_final))[56] <- "cohort"
colnames(pData(imppost_final))[56] <- "cohort"


## REMOVE CpG SITES WITH SNPs

#Load dataset with masked CpG sites (https://zwdzwd.github.io/InfiniumAnnotation)
snps_cpgs <- read.table("/home/isglobal.lan/ncarreras/homews/TFM/Data_Tables/450Khg19SNP.tsv", header=TRUE, sep = "\t")

#Select CpG sites withuot SNPs
cpgs_unmask <- snps_cpgs[which(!snps_cpgs$MASK_general),]$probeID
methy_final <- methy_final[rownames(methy_final) %in% cpgs_unmask,]


## REMOVE TRANSCRIPTS WITHOUT GENE SYMBOL
dd <- featureData(trans_final)
genes <- rownames(dd[which(dd$GeneSymbol_Affy!=""),])
trans_final <- trans_final[genes,]


## SAVE FINAL DATASETS

save(methy_final,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Methy_final.Rdata")
save(trans_final,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Trans_final.Rdata")
save(imppreg_final,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppreg_final.Rdata")
save(imppost_final,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppost_final.Rdata")


## CREATE DATASETS FOR EACH INVERSION REGION +/-1 Mb

#Load inversion information
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")

region_inv <- function(inv){
  start <- as.numeric(start(inversionGR[inv,]))-1000000
  end <- as.numeric(end(inversionGR[inv,]))+1000000
  region <- paste(as.character(seqnames(inversionGR[inv,])),":",start,"-",end,sep="")
  return(region)
}

#Define the inversion region coordinates +/- 1Mb
region8 <- GRanges(region_inv("inv8_001"))
region16 <- GRanges(region_inv("inv16_009"))
region17 <- GRanges(region_inv("inv17_007"))

#Subset the features that overlap with the inversion region
methy8 <- subsetByOverlaps(methy_final, region8)
save(methy8,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
methy16 <- subsetByOverlaps(methy_final, region16)
save(methy16,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy16.Rdata")
methy17 <- subsetByOverlaps(methy_final, region17)
save(methy17,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy17.Rdata")

#For transcriptome, first transform ExpressionSet to RangedSummarizedExperiment
trans_ranged <- makeSummarizedExperimentFromExpressionSet(trans_final)

trans8 <- subsetByOverlaps(trans_ranged, region8)
save(trans8,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans8.Rdata")
trans16 <- subsetByOverlaps(trans_ranged, region16)
save(trans16,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans16.Rdata")
trans17 <- subsetByOverlaps(trans_ranged, region17)
save(trans17,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans17.Rdata")


## DESCRIBE THE DATA

#Methylation data
intersect_samples<-intersect(colnames(methy8),colnames(imppreg_final))
imppreg_methy <- imppreg_final[,intersect_samples]

phenos <- data.frame(sex=imppreg_methy$sex,
                     age=(imppreg_methy$hs_child_age_days_None)/365,
                     cohort=imppreg_methy$cohort,
                     trim_conception=imppreg_methy$h_trimcon_None,
                     parity=imppreg_methy$h_parity_None,
                     maternal_education=imppreg_methy$h_edumc_None,
                     maternal_bmi=imppreg_methy$h_mbmi_None,
                     maternal_age=imppreg_methy$h_age_None,
                     maternal_smoke=expos(imppreg_methy)$e3_asmokyn_p_None,
                     inv8p23.1=imppreg_methy$inv8_001,
                     inv16p11.2=imppreg_methy$inv16_009,
                     inv17q21.31=imppreg_methy$inv17_007)
rownames(phenos) <- colnames(imppreg_methy)

table_four <- tableby(~sex + age + cohort + trim_conception + parity +
                        maternal_education + maternal_bmi + maternal_age + 
                        maternal_smoke + inv8p23.1 + inv16p11.2 + inv17q21.31,
                      data = phenos) 

write2(table_four, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Data_description/Methylation_data.doc", output_format="word")


#Transcription data
intersect_samples<-intersect(colnames(trans8),colnames(imppost_final))
imppreg_trans <- imppreg_final[,intersect_samples]

phenos <- data.frame(sex=imppreg_trans$sex,
                     age=(imppreg_trans$hs_child_age_days_None)/365,
                     cohort=imppreg_trans$cohort,
                     trim_conception=imppreg_trans$h_trimcon_None,
                     parity=imppreg_trans$h_parity_None,
                     maternal_education=imppreg_trans$h_edumc_None,
                     maternal_bmi=imppreg_trans$h_mbmi_None,
                     maternal_age=imppreg_trans$h_age_None,
                     maternal_smoke=expos(imppreg_trans)$e3_asmokyn_p_None,
                     inv8p23.1=imppreg_trans$inv8_001,
                     inv16p11.2=imppreg_trans$inv16_009,
                     inv17q21.31=imppreg_trans$inv17_007)
rownames(phenos) <- colnames(imppreg_trans)

table_four <- tableby(~sex + age + cohort + trim_conception + parity +
                        maternal_education + maternal_bmi + maternal_age + 
                        maternal_smoke + inv8p23.1 + inv16p11.2 + inv17q21.31,
                      data = phenos)

write2(table_four, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Data_description/Transcription_data.doc", output_format="word")
