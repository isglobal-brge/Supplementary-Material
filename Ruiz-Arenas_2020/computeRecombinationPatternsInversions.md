# Summary

This document shows how to compute the recombination patterns with FastEPRR. We used this code to compute the recombination patterns for each recombClust group of inv8p23.1 in European samples of 1000 Genomes. It requires the R package `FastEPRR`, and the python script FASTEPPR_parser.py to process the output.

# Preprocess data

Select variants in the inversion regions from the 1000 Genomes VCF file.

```{bash}
KG=~/PublicData/STUDY/1000GENOME/VCF/ ## Directory to 1000 Genomes VCF files

## Shorten original 1000 Genomes VCFs to speed up analysis
vcftools --gzvcf $KG/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr 8 --from-bp 8000000 --to-bp 12000000 --out chr8 --recode 
```

Create required folders to run the method. 

```{bash}
mkdir NN
mkdir NN/step1
mkdir NN/step1/step1
mkdir NN/step2

mkdir II/step1/step1
mkdir II
mkdir II/step1/
mkdir II/step1/step1
mkdir II/step2

cp ~/data/CarlosRuiz/Inversions/LDclassifier/EuropeanSamples/PCsinvs8_17_EUR_pruned0.4.Rdata . ## Object with recombClust classification for inv8p23.1
```

# Run FastEPRR

Run FastEPRR on pre-filtered files. 

## Load libraries and files

```{r}
library(FastEPRR)
library(parallel)

load("PCsinvs8_17_EUR_pruned0.4.Rdata") ## Object with recombClust classification for inv8p23.1
```

## Recover cluster classification

This step gets recombClust classification for inversion 8. In the current version of recombClust, this step is already performed by `runRecombClust`.

```{r}
class8 <- kmeans(pc8_04$x[, 1], centers = 2, nstart = 1000)$cluster
names(class8) <- ind_ids
```

## Divide chromosomes by recombClust classification

Recombination patterns will be computed independently for each recombClust group. We considered that chromosomes classified as 1 are inverted while chromosomes classified as 2 are standard. In this step, we will define a vector to select the chromosomes in each of the groups.

### List of inverted chromosomes

```{r}
II <- sapply(seq(1, length(class8), 2), function(x) { 
  a <- paste(ifelse(class8[x:(x+1)] == 1, "1", "0"), collapse = ":")
  paste0(names(class8)[x], "[", a, "]")
})
### Remove samples with none of the chromosomes included
II <- II[!grepl("0:0", II, fixed = TRUE)]
II <- paste(II, collapse = ";")
```

### List of standard chromosomes
```{r}
NN <- sapply(seq(1, length(class8), 2), function(x) { 
  a <- paste(ifelse(class8[x:(x+1)] == 2, "1", "0"), collapse = ":")
  paste0(names(class8)[x], "[", a, "]")
})
NN <- NN[!grepl("0:0", NN, fixed = TRUE)]
NN <- paste(NN, collapse = ";")
```

## Apply algorithm

FastEPRR consists of 3 steps. They are run independently for each groups of chromosomes (inverted and standard). 

### Run step 1
```{r}
## NN
FastEPRR_VCF_step1(vcfFilePath = "./chr8.recode.vcf", 
                   winLength="80", stepLength = "50", 
                   idvlConsidered= NN,
                   srcOutputFilePath= "NN/step1/step1")

## II
FastEPRR_VCF_step1(vcfFilePath="./chr8.recode.vcf", 
                   winLength="80", stepLength = "50", 
                   idvlConsidered= II,
                   srcOutputFilePath= "./II/step1/step1")
```

### Run step 2

```{r}
## NN
mclapply(1:25, FastEPRR_VCF_step2, 
         srcFolderPath="./NN/step1/",
         jobNumber = 25, mc.cores = 20,
         DXOutputFolderPath = "./NN/step2")


## II
mclapply(1:25, FastEPRR_VCF_step2, 
         srcFolderPath="./II/step1",
         jobNumber = 25, mc.cores = 20,
         DXOutputFolderPath = "./II/step2")

```

### Run step3

```{r}
## NN
FastEPRR_VCF_step3(srcFolderPath = "./NN/step1",
                   DXFolderPath="./NN/step2", 
                   finalOutputFolderPath="./NN")

## II
FastEPRR_VCF_step3(srcFolderPath = "./II/step1",
                   DXFolderPath="./II/step2", 
                   finalOutputFolderPath="./II")

```

## Correct output

This python script converts FastEPRR output to a table format.

```{bash}
python ~/data/CarlosRuiz/Inversions/InversionSequencing/FASTEPPR_parser.py -i II/chr_8 -o II.txt
python ~/data/CarlosRuiz/Inversions/InversionSequencing/FASTEPPR_parser.py -i NN/chr_8 -o NN.txt
```
