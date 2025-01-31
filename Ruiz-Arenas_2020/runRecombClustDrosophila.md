# Summary

Run recombClust in inversions in Drosophila. Pre-filter SNPs with MAF > 0.05 and call rate > 0.95 and individuals with call rate > 0.95.

## Load libraries

```{r}
library(recombClust)
library(GenomicRanges)
library(VariantAnnotation)
library(parallel)
library(snpStats)
library(BiocParallel)

fold <- "results/models/"
```


## Load inversion ranges

```{r}
invReg <- read.csv2("data/ranges.csv", header = TRUE) ## csv with inversion ranges
invReg$end <- as.numeric(gsub(",", "", invReg$End.Phys))
invReg$start <- as.numeric(gsub(",", "", invReg$Start.Phys))
invgr <- makeGRangesFromDataFrame(invReg, keep.extra.columns = TRUE)

names(invgr) <- invgr$Inversion 
```

## Define functions

```{r}
getRCmods <- function(range){
  message(names(range))
  message("Loading VCF")
  vcf <- readVcf("data/dgrp2.vcf.gz", "dm6", ScanVcfParam(which = range))
  
  snpsVCF <- genotypeToSnpMatrix(vcf)
  snpsVCF$genotypes <- snpsVCF$genotypes[, !snpsVCF$map$ignore]
  
  ## Remove individuals with low Call rate
  sumsr <- row.summary(snpsVCF$genotypes)
  inds <- rownames(snpsVCF$genotypes)[sumsr$Call.rate > 0.95]
  
  ## Remove missings in data
  sums <- col.summary(snpsVCF$genotypes[inds, ])
  snps <- colnames(snpsVCF$genotypes)[sums$MAF > 0.05 & sums$Call.rate == 1]
  
  vcf <- vcf[snps, inds]
  
  if (ncol(vcf) < 10){
    message("The region does not contain enough polymorphic SNPs")
    recomb <- NULL
    save(range, recomb, file = paste0(fold, names(range), "recombRes.Rdata")) 
    return(NULL)
  }
  
  ## Get phased genotypes
  geno <- geno(vcf)$GT
  phase <- lapply(1:ncol(geno), function(x){
    chr1 <- as.numeric(sapply(geno[, x], substring, 1, 1))
    chr2 <- as.numeric(sapply(geno[, x], substring, 3, 3))
    matrix(c(chr1, chr2), nrow = 2, byrow = TRUE)
  })
  phase <- Reduce(function(...) rbind(...), phase)
  rownames(phase) <- paste(rep(colnames(geno), each = 2), 1:2, sep = "_")
  colnames(phase) <- rownames(geno)


  message("Running recombClust")
  recomb <- runRecombClust(phase, annot = rowRanges(vcf), BPPARAM = MulticoreParam(20))
  
  save(range, recomb, file = paste0(fold, names(range), "recombRes.Rdata")) 
  NULL
}
```

Run recombClust. It stores the results in .Rdata files.

```{r}
lapply(names(invgr), function(x) getRCmods(invgr[x]))
```

