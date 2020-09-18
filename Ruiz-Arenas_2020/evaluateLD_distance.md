# Summary

This documents show how to generate Sup Fig S2, the association between linkage disequilibrium R2 difference and distance for inv8p23.1.

# Load libraries and data

```{r}
library(VariantAnnotation)
library(GenomicRanges)
library(snpStats)
library(ggplot2)
library(dplyr)
```

We need three files to do the computation:

1. inv8p23.1 genotypes
2. 1000 Genomes chromosome 8 VCF file
3. 1000 Genomes ancestry (to select European individuals)

```{r}
load("/SYNCRW10125/DATASETS/STUDY/1000GENOME/Samples_Pop1GK.Rdata") ### Load data.frame with 1000G ancestry
load("/home/cruiz/InversionSequencing/SandersROIs/invClust/invClustEURNoSeqDups.Rdata") ## Contains invsGeno 
```

# Compute LD (linkage disequilibrium) differences 

## Load genotype data

First, we loaded genotype data from inv8p23.1 regions in EUR individuals of 1000 Genomes. In this step, we used  `getVCFMatrix`, an in-house function based on `VariantAnnotation`, which loads the VCF in the SnpStats format. `ROIsGR["ROIno.8.3"]` represent the coordinates of inv8p23.1.

```{r}
EUR <- rownames(samp_pop)[samp_pop$pop %in% c("GBR", "FIN", "IBS", "TSI", "CEU")] ## Select 1000G from EUR ancestry

IIVCF <- NNVCF <- VCF <- getVCFmatrix(ROIsGR["ROIno.8.3"], EUR) ## Load genotypes from inversion region from EUR individuals of 1000Genomes
```

In this step, we make a two vectors with homozygous inverted and standard samples.

```{r}
invs8 <- invsGeno$ROIno.8.3 ### Inversion classification of inv8p23.1

## Select samples that are homozygous for the inversion or standard
II <- names(invs8[invs8 == "I/I"])
NN <-  names(invs8[invs8 == "NI/NI"])
```

We created two objects, one with II individuals (homozygous inverted) and another for NN (homozygous standard)

```{r}
IIVCF$genotypes <- IIVCF$genotypes[II, ]
NNVCF$genotypes <- NNVCF$genotypes[NN, ]

rownames(IIVCF$map) <- IIVCF$map$snp.names
rownames(NNVCF$map) <- NNVCF$map$snp.names
```

## Compute Linkage Disequilibrium with SNPstats

```{r}
ldII <- list(res = ld(IIVCF$genotypes, stats=c("R.squared", "D.prime"), 
                      depth = ncol(IIVCF$genotypes)), annot = IIVCF$map)

ldNN <- list(res = ld(NNVCF$genotypes, stats=c("R.squared", "D.prime"), 
                      depth = ncol(NNVCF$genotypes)), annot = NNVCF$map)
```

We computed the difference in R2 between II and NN individuals

```{r}
ldDiff <- ldNN 
ldDiff$res$R.squared <- ldII$res$R.squared - ldNN$res$R.squared
ldDiff$res$D.prime <- ldII$res$D.prime - ldNN$res$D.prime
```

## Prepare data

We added SNP names pairs to the LD data. 

```{r}
ldDiffMat <- ldDiff$res$R.squared

ldVecDiff <- ldDiffMat[upper.tri(ldDiffMat)]
rows <- ncol(ldDiffMat)
z <- sequence(rows)
names(ldVecDiff) <- paste(colnames(ldDiffMat)[unlist(lapply(1:(rows-1), function(x) 1:x), use.names = FALSE)], 
                      rownames(ldDiffMat)[rep(z[-1], times = z[-length(z)])], sep = "-")


pairs <- strsplit(names(ldVecDiff), "-")
mat <- matrix(unlist(pairs), ncol = 2, byrow = TRUE)
```

We added the distance between SNPs in the pairs and binned in four groups.

```{r}
distances <- VCF$map[mat[, 2], "position"] - VCF$map[mat[, 1], "position"]
names(distances) <- names(ldVecDiff)


df2 <- data.frame(DistanceBins = cut(distances, c(0, 1e3, 1e4, 1e5, 1e7)), 
                 DiffBins = cut(abs(ldVecDiff), seq(0, 1, 0.2)))
levels(df2$DistanceBins) <- c("0-1Kb", "1-10Kb", "10-100Kb", ">100Kb")
tabdf <- data.frame(prop.table(table(df2$DistanceBins, df2$DiffBins), margin = 2))
```

## Create plot

```{r}
ggplot(tabdf, aes(y = Freq*100, x = Var2, fill = Var1)) + geom_bar(stat = "identity") + 
  scale_y_continuous(name = "") +
  scale_fill_discrete(name = "SNPs Distance") +
  scale_x_discrete(labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")) +
  xlab(bquote(R^2~Difference)) + 
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.text.x  = element_text(angle=45, vjust=0.5, size=15),
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20))
```



