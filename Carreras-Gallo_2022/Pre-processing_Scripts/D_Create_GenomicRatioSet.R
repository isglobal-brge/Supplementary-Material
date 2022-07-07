############################
## CREATE GENOMICRATIOSET ##
############################

#Set "clean" or "whole" dataset
dataset <- "clean"

#Load libraries
library(minfi)
library(meffil)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#Set working directory
setwd("/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/Functional_Normalization")

#Load beta values
load(paste0("norm_beta_pc10_",dataset,".Robj"))

#Load metadata
load("/PROJECTES/GENOMICS/TruDiagnostic/metadata/metadata.Rdata")

#Load cellular composition
cellcounts <- read.table(paste0("/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/QC/cellcounts_",dataset,".txt"),
                         header=T)

#Load samplesheet to obtain the Slide variable
load("/PROJECTES/GENOMICS/TruDiagnostic/idat_files/idat_use/SampleSheet.Rdata")

#Select the overlapping IDs
IDs <- intersect(metadata$PatientID,colnames(norm.beta))
metadata <- metadata[which(metadata$PatientID%in%IDs),]
norm.beta <- norm.beta[,IDs]
samplesheet <- samplesheet[which(samplesheet$Sample_Name%in%IDs),]

gc(reset=TRUE)

dim(norm.beta)

#Add cellcounts to metadata
pheno <- merge(metadata, cellcounts, by="PatientID")

#Add Slide (as.numeric) to metadata
pheno$Slide <- as.numeric(as.factor(samplesheet$Slide))
rownames(pheno) <- pheno$PatientID

gc(reset=TRUE)

#Create the GenomicRatioSet with the beta values, the pheno data, and the default annotation by Illumina
GRset <- makeGenomicRatioSetFromMatrix(mat = norm.beta,
                                       rownames = rownames(norm.beta),
                                       pData = pheno,
                                       array = "IlluminaHumanMethylationEPIC", 
                                       annotation = "ilm10b4.hg19",
                                       mergeManifest = FALSE,
                                       what = "Beta")

#Garbage collector
gc(reset=TRUE)

### REMOVE CpG SITES WITH SNPS

#Load dataset with masked CpG sites (https://zwdzwd.github.io/InfiniumAnnotation)
snps_cpgs <- read.table("/PROJECTES/GENOMICS/TruDiagnostic/Processing_Scripts/Meffil/EPIC.hg19.manifest.tsv", header=TRUE, sep = "\t")

#Select CpG sites withuot SNPs
cpgs_unmask <- snps_cpgs[which(!snps_cpgs$MASK_general),]$probeID
dim(GRset)
GRset <- GRset[rownames(GRset) %in% cpgs_unmask,]
dim(GRset)

#Include the most important information about the CpG sites in the rowData()
rowData(GRset) <- getAnnotation(GRset)[,c("chr","pos","strand","UCSC_RefGene_Group","UCSC_RefGene_Name")]


### CHECK THE ANNOTATION OF THE CpG SITES (GENE SYMBOL AND GROUP)
#First, we remove the redundancy by selecting only the unique gene-group pairs
#Second, we rename the genes using the last version of HGNC

#Load txt file with HGNC Gene Symbols
hgnc <- read.table(file="/PROJECTES/GENOMICS/TruDiagnostic/Processing_Scripts/Meffil/HGNC_gene_symbols.txt", 
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

#Garbage collector
gc(reset=TRUE)

#Create a function to remove the redundancy and to get the HGNC gene symbols
redundancy_HGNC <- function(i){
  
  print(i)
  
  #Get the groups
  groups <- strsplit(rowData(GRset)$UCSC_RefGene_Group[i],";")[[1]]
  
  #Get the genes
  genes <- strsplit(rowData(GRset)$UCSC_RefGene_Name[i],";")[[1]]
  
  #Make "group;gene" pairs
  pairs <- paste(groups,genes,sep=";")
  
  #Select the non-repetitive pairs
  pairs <- unique(pairs)
  
  #Separate the genes and the groups and save them
  df <- do.call(rbind,strsplit(pairs,";"))
  groups <- paste(df[,1],collapse=";")
  genes <- paste(df[,2],collapse=";")
  
  #Check the HGNC last version for the genes
  #Separate each gene symbol in a different character
  symbols <- strsplit(genes,";")[[1]]
  
  #Stop fot the vectors without genes annotated
  if (length(symbols)==0){
    return("")
  }
  
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
  df <- data.frame(UCSC_RefGene_Group=groups,
             UCSC_RefGene_Name=genes,
             HGNC_GeneSymbol=symbols)
  return(df)
}

list.df <- parallel::mclapply(1:nrow(GRset),redundancy_HGNC,mc.cores=3)


#Substitute the previous columns of the rowData for the new ones
rowData(GRset) <- cbind(rowData(GRset)[,1:3],do.call(rbind,list.df))

#Save GRset
save(GRset,file=paste0("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_",dataset,".Rdata"))

GRset

# # CHANGE PHENOTYPE DATA
# dataset <- "clean_imp"
# load(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_",dataset,".Rdata"))
# load("/PROJECTES/GENOMICS/TruDiagnostic/metadata/metadata.Rdata")
# metadata <- metadata[which(metadata$PatientID%in%colnames(GRset)),]
# cellcounts <- read.table(paste0("/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/QC/cellcounts_clean.txt"),
#                          header=T)
# 
# pheno <- merge(metadata, cellcounts, by="PatientID")
# rownames(pheno) <- pheno$PatientID
# 
# pData(GRset) <- as(pheno,"DataFrame")
# gc(reset=TRUE)
# save(GRset,file=paste0("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_",dataset,".Rdata"))
