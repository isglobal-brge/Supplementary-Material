####################################
#### PRE-PROCESSING OF DATASETS ####
####################################

## LOAD LIBRARIES

library(SummarizedExperiment)
library(GenomicRanges)
library(rexposome)
library(minfi)
library(stringr)
library(plyr)
library(tidyverse)
library(arsenal) 
library(readxl)


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

#Principal components from GWAS
PCs <- as.data.frame(read_excel("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA_GWAS.xlsx",sheet = "metadataEUR"))
rownames(PCs) <- PCs$HelixID

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


## ADD THE INVERSIONS AND PCS FROM GWAS IN THE pData AND DISCARD THE NON-EUROPEAN SAMPLES

inversions <- select(scClassDF, -c("inv3_003","inv11_004","inv12_006","inv14_005")) #Discard the ones without 3 genotypes

#Sort inversion genotypes as NN, NI, and II
for (inv in colnames(inversions)){
  inversions[[inv]] <-  factor(inversions[[inv]], levels=c("NN", "NI", "II"))
  inversions[[inv]] <- revalue(inversions[[inv]], c("NN"="N/N", "NI"="N/I", "II"="I/I"))
}

add_inv_cauc <- function(dataset){
  EUR_ID <- intersect(colnames(dataset),rownames(PCs))
  final <- dataset[,EUR_ID]
  pData(final) <- cbind(pData(final),inversions[EUR_ID, ], PCs[EUR_ID,17:26])
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

#age_sample_years --> age
colnames(pData(methy_final))[10] <- "age"
colnames(pData(trans_final))[9] <- "age"

## REMOVE CpG SITES WITH SNPs

#Load dataset with masked CpG sites (https://zwdzwd.github.io/InfiniumAnnotation)
snps_cpgs <- read.table("/home/isglobal.lan/ncarreras/homews/TFM/Data_Tables/450Khg19SNP.tsv", header=TRUE, sep = "\t")

#Select CpG sites withuot SNPs
cpgs_unmask <- snps_cpgs[which(!snps_cpgs$MASK_general),]$probeID
methy_final <- methy_final[rownames(methy_final) %in% cpgs_unmask,]


## REMOVE TRANSCRIPTS WITHOUT GENE SYMBOL AND WITH CALLRATE <20%

dd <- featureData(trans_final)
genes <- rownames(dd[which(dd$GeneSymbol_Affy!="" & dd$CallRate>=20),])
trans_final <- trans_final[genes,]

## REPLACE THE GENE SYMBOL FOR THE ONES APPROVED BY HGNC

#Load txt file with HGNC Gene Symbols
hgnc <- read.table(file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/HGNC_gene_symbols.txt", 
                   header=T, sep="\t", quote="", fill=T)

#Create a list with all the previous and alias symbols for a gene
other_symbols <- function(num){
  app.symbol <- hgnc$Approved.symbol[num]
  prev.symbol <- strsplit(hgnc$Previous.symbols[num],", ")[[1]]
  ali.symbol <- strsplit(hgnc$Alias.symbols[num],", ")[[1]]
  other.symbol <- c(prev.symbol,ali.symbol)
  return(other.symbol)
}

hgnc_list <- sapply(1:nrow(hgnc),other_symbols)
names(hgnc_list) <- hgnc$Approved.symbol

#Create a function for changing the gene symbols in our datasets
change_symbol <- function(i, data){
  print(i)
  if (data=="trans"){
    gene_affy <- fData(trans_final)$GeneSymbol_Affy[i]
  }
  if (data=="methy"){
    gene_affy <- rowData(methy_final)$UCSC_RefGene_Name[i]
  }
  #Skip features without Gene Symbol annotated
  if (gene_affy==""){
    return("")
  }
  #Separate each gene symbol in a different character
  symbols <- strsplit(gene_affy,";")[[1]]
  #Through each gene, search the approved gene symbol in HGNC
  for (j in 1:length(symbols)){
    if(!symbols[j]%in%names(hgnc_list)){
      for (n in 1:length(hgnc_list)){
        if(symbols[j]%in%hgnc_list[[n]]){
          symbols[j] <- names(hgnc_list)[n]
          break
        }
      }
    }
  }
  #Join all the gene symbols
  if (length(symbols)>1){
    symbols <- paste0(symbols,collapse=";")
  }
  return(symbols)
}

##Transcription
symbols_list <- parallel::mclapply(1:nrow(trans_final),change_symbol, 
                                   mc.cores=10, data= "trans")

#Replace the gene symbols in the trans_final dataset
fData(trans_final)$GeneSymbol_Affy <- do.call(c, symbols_list)

#Methylation
symbols_list_methy <- parallel::mclapply(1:nrow(methy_final),change_symbol, 
                                         mc.cores=10, data= "methy")

#Replace the gene symbols in the methy_final dataset
rowData(methy_final)$UCSC_RefGene_Name <- do.call(c, symbols_list_methy)


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


## CREATE DESCRIPTIVE TABLES

#Methylation data
phenos_methy <- data.frame(sex=methy8$sex,
                           age=methy8$age,
                           cohort=methy8$cohort,
                           Inversion.8p23.1=methy8$inv8_001,
                           Inversion.16p11.2=methy8$inv16_009,
                           Inversion.17q21.31=methy8$inv17_007)

phenos_methy <- cbind(phenos_methy, pData(methy8)[,40:45])

table_methy <- tableby(formula(paste("~", paste(colnames(phenos_methy),collapse=" + "))),
                      data = phenos_methy) 

write2(table_methy, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Data_description/Methylation_data.doc", output_format="word")


#Transcription data
phenos_trans <- data.frame(sex=trans8$sex,
                     age=trans8$age,
                     cohort=trans8$cohort,
                     Inversion.8p23.1=trans8$inv8_001,
                     Inversion.16p11.2=trans8$inv16_009,
                     Inversion.17q21.31=trans8$inv17_007)

phenos_trans <- cbind(phenos_trans, colData(trans8)[,18:23])

#Select samples with cell type information
phenos_trans <- phenos_trans[which(!is.na(phenos_trans$NK_6)),]
table_trans <- tableby(formula(paste("~", paste(colnames(phenos_trans),collapse=" + "))),
                       data = phenos_trans) 

write2(table_trans, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Data_description/Transcription_data.doc", output_format="word")
