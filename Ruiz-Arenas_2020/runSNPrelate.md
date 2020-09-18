# Summary

Run SNPrelate in the regions described in the paper.

# Load libraries

```{r}
library(SNPRelate)
```

# Load data (data gathered with getFilesSNPrelate.sh)

```{r}
load("data/Samples_Pop1GK.Rdata")
load("data/LCTmodels.Rdata")
recombLCT <- recomb
load("data/TAR_reg_models.Rdata")
load("data/HumanInvsPCs.Rdata")
```


# Load Drosophila inversions

```{r}
invReg <- read.csv2("data/ranges.csv", header = TRUE)
invReg$end <- as.numeric(gsub(",", "", invReg$End.Phys))
invReg$start <- as.numeric(gsub(",", "", invReg$Start.Phys))

load("data/Dros2Lmodels.Rdata")
recomb2L <- recomb
load("data/Dros2Rmodels.Rdata")
recomb2R <- recomb
load("data/Dros3Rmodels.Rdata")
recomb3R <- recomb
```

# Convert VCF files to gds

```{r}
chrs <- c(paste0("data/chr", c(1, 2, 8, 17), ".vcf.gz"), "data/dgrp2.vcf.gz")
names(chrs) <- c(paste0("chr", c(1, 2, 8, 17)), "dros")

preproc <- "results/preproc/"

lapply(seq_len(length(chrs)), function(x)
  snpgdsVCF2GDS(chrs[x], paste0(preproc, names(chrs)[x], ".gds", method="biallelic.only"))
  )
```

# Run SNPrelate in humans

Run SNPrelate in European individuals of 1000 Genomes.

```{r}
EUR <- rownames(samp_pop)[samp_pop$superpop == "EUR"]

regionRanges <- data.frame(reg = c("chr1", "chr2", "chr8", "chr17"),
                            start = c(145550000, 135700000, 8055789, 43661775),
                            end = c(145750000, 136900000, 11980649, 44372665))

gdsFiles <- paste0(preproc, names(chrs)[1:4], ".gdsbiallelic.only")
names(gdsFiles) <- names(chrs)[1:4]
```

Run SNPrelate.

```{r}
getPCAs <- function(gdsFile, start, end, samples = NULL, chromosome = NULL){

  gds <- snpgdsOpen(gdsFile)
  pos <- read.gdsn(index.gdsn(gds, "snp.position"))
  
  if (!is.null(chromosome)){
    chr <-  read.gdsn(index.gdsn(gds, "snp.chromosome"))
    idx <- which(chr == chromosome & pos > start & pos < end)
  } else {
    idx <- which(pos > start & pos < end)
  }
  
  if (!is.null(samples)){
    selSamps <- intersect(samples, read.gdsn(index.gdsn(gds, "sample.id")))
    pca <- snpgdsPCA(gds, sample.id = selSamps, snp.id = idx, num.thread = 2, maf = 0.05)  
  } else {
    pca <- snpgdsPCA(gds, snp.id = idx, num.thread = 2, maf = 0.05)  
  }
  closefn.gds(gds)
  pca
}

genoPCAs <- lapply(names(gdsFiles), function(x) 
  getPCAs(gdsFiles[x], start = regionRanges[regionRanges$reg == x, "start"], 
          end = regionRanges[regionRanges$reg == x, "end"], samples = EUR))
```

# Run SNPrelate in Drosophila

```{r}
drosReg <- subset(invReg, Inversion %in% c("In(2L)t", "In(2R)NS", "In(3R)Mo"))

drosPCAs <- lapply(seq_len(nrow(drosReg)), function(x) 
  getPCAs("results/preproc/dros.gdsbiallelic.only", start = drosReg[x, "start"], 
          end = drosReg[x, "end"], chromosome = drosReg[x, "Chromosome"]))
```


# Get individuals classification

Convert recombClust at chromosome level to individuals classification.

```{r}
getDrosophilaClass <- function(class){
  
  names(class) <- substring(names(class), 1, nchar(names(class))-2)
  class <- class[!duplicated(names(class))]
  class
}

getGenoClass <- function(class){
  genoclass <- sapply(seq(1, length(class), 2), function(x) { 
    a <- paste(sort(class[x:(x+1)]), collapse = "-")
  })
  sampNames <- substring(names(class), 1, nchar(names(class))-2)
  names(genoclass) <- sampNames[!duplicated(sampNames)]
  genoclass
}


class8 <- kmeans(pc8_04$x[, 1], centers = 2, nstart = 1000)$cluster
names(class8) <-  paste(ind_ids, 1:2, sep = "_")
class17 <- kmeans(pc17_04$x[, 1], centers = 2, nstart = 1000)$cluster
names(class17) <- paste(ind_ids17, 1:2, sep = "_")



humanClass <- list(inv8 = getGenoClass(class8), 
                   inv17 = getGenoClass(class17),
                   lct = getGenoClass(recombLCT$class),
                   tar = getGenoClass(thousand_res$class))
drosClass <- list(In2L = getDrosophilaClass(recomb2L$class), 
                  In2R = getDrosophilaClass(recomb2R$class),
                  In3R = getDrosophilaClass(recomb3R$class))

save(drosPCAs, genoPCAs, humanClass, drosClass, file = "results/plotFiles/SNPrelate_results.Rdata")
```
