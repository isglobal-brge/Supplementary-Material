# Summary

This document contains the script to generate most of the paper plots.

# Figure 3: Recombination patterns inversion 8


## Load libraries

```{r}
library(Gviz)
library(GenomicRanges)
```

## Load FastEPRR data

Files created with computeRecombinationPatternsInversions.Rmd

```{r}
II <- read.table("II.txt", header = TRUE)
NN <- read.table("NN.txt", header = TRUE)
```

## Load recombClust data

```{r}
load("recombClustPatterns.Rdata") ## Recombination patterns for inv8p23.1 estimated from recombClust
width <- 5e4
starts <- seq(8000000, 12000000, width)
```

## Generate annotation tracks

These tracks contain the ideogram and the transcripts in the region.

```{r}
basetracks <- list(IdeogramTrack(genome = "hg19", chromosome = "chr8"), 
                                 Gviz::GenomeAxisTrack())
range <- GRanges("chr8:8055789-11980649")
  
## Load transcripts data
data(dmrcatedata, envir = environment(), package = "DMRcatedata")

## Select only protein coding genes
txs <- subsetByOverlaps(tx.hg19, range)
txs <-  subset(txs, gene_type == "protein_coding")

genes <- Gviz::GeneRegionTrack(txs, name = "Transcripts", 
                      symbol = txs$gene_name, 
                      fill = "lightblue", 
                      gene = txs$gene_name,
                      showId = TRUE, geneSymbol = TRUE, cex.title = 0.7,
                      shape = "arrow", transcriptAnnotation = "symbol",
                      collapseTranscripts = TRUE, rotation.title = 0)
```

## Alves hotspots

Add hotspots defined by Alves and colleagues.

```{r}
Adf <- data.frame(chr = "chr8", start = c(9300000, 10800000), 
                  end = c(9350000, 10850000), 
                  group = c("Specific Standard Recombination Peak", "Specific Inverted Recombination Peak"))
AlvesGR <- makeGRangesFromDataFrame(Adf, keep.extra.columns = TRUE)

AlvesTrack <- Gviz::AnnotationTrack(AlvesGR,
                                    name = "Alves et al peaks", shape = "box",
                                    feature = c("std", "inv"),
                                    rotation.title = 0,
                                    cex.title = 0.7)
```

## FastEPPR track

```{r}

rownames(II) <- II$Start
rownames(NN) <- NN$Start
allNames <- union(rownames(II), rownames(NN))
II <- II[allNames, ]
rownames(II) <- allNames
II$NNRho <- NN[rownames(II), "Rho"]
II$Start <- II$Start * 1000
II$End <- II$End * 1000
II$Chr <- "chr8"

FGR <- makeGRangesFromDataFrame(II[, -(4:5)], keep.extra.columns = TRUE)
## Solve NAs manually
FGR$NNRho[20:22] <- (FGR$NNRho[23]-FGR$NNRho[19])/4*1:3 + FGR$NNRho[19]
FastEPPRtrack <- DataTrack(FGR, type ="a", name = "FastEPPR\nrecombination rate\n(cM/Mb)",
                           groups = c("Inverted", "Standard"),
                           col = c("#334491", "#349642"),
                           cex.title = 0.7)

FGR$Diff <- FGR$NNRho - FGR$Rho
FEPPRDifftrack <- DataTrack(FGR, type =c("a", "mountain"), baseline = 0, data =  FGR$Diff, 
                            name = "Recombination\ndifference (STD-INV)", 
                            col = "black",
                            fill.mountain = c("#536ee9", "#52e768"),
                            cex.title = 0.7)
```
                           
## RecombClust track

Track with the recombination patterns obtained with recombClust.

```{r}
df <- data.frame(chr = "chr8", start = starts[2:80], end = starts[2:80] + 5e4, 
                 RhoI = clusProbs[[2]], RhoN = clusProbs[[1]])
recombGR <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
recombTrack <- DataTrack(recombGR, type ="a", name = "Recomb proportion",
                         groups = c("Inverted", "Standard"),
                         col = c("#334491", "#349642"),
                         cex.title = 0.7)


recombGR$Diff <- recombGR$RhoN - recombGR$RhoI
recombDifftrack <- DataTrack(recombGR, type =c("a", "mountain"),  baseline = 0,
                             data =  recombGR$Diff, 
                             name = "Difference in\nrecomb proportion\n(STD-INV)", 
                             col = "black", 
                             fill.mountain = c("#536ee9", "#52e768"),
                             cex.title = 0.7)
```

## Make plot

```{r}
png("inv8RecombGviz.png", width = 30, height = 25, units = 'cm', res = 300)
gviz <- Gviz::plotTracks(c(basetracks, genes, AlvesTrack, FastEPPRtrack, FEPPRDifftrack, 
                   recombTrack, recombDifftrack), 
                 sizes = c(1, 1, 4, 1, 4, 4, 4, 4), 
                 from = start(range), to = end(range), 
                 groupAnnotation = "group",
                 std = "#349642", inv = "#334491", 
                 fontcolor.title = "black",
                 col.axis = "black",
                 background.title = "grey90")
dev.off()
```

# Figure 2: recombClust PCA of human and Drosophila chromosomal inversions

## Load libraries 

```{r}
library(ggplot2)
library(tidyr)
library(cowplot)
library(readxl)
library(dplyr)
library(GenomicRanges)
```

## Human inversions

These objects were obtained with an old version of recombClust. With the release version, both objects can be obtained with the `runRecombClust` function.

```{r}
load("Inversions/LDclassifier/EuropeanSamples/PCsinvs8_17_EUR.Rdata") ## Results from recombClust for inversions inv8p23.1 and inv17q21.31

pcList <- list(inv8 = pc8_04, inv17 = pc17_04) ## recombClust PCAs 
idsList <- list(inv8 = ind_ids, inv17 = ind_ids17) ## chromosome IDs
```

### inv8p23.1

Create data.frame with PCA values and inversion genotypes. These inversion genotypes were obtained with experimental methods.

```{r}
pc <- pcList[["inv8"]]$x
idsGenos <- idsList[["inv8"]]

df <- data.frame(pc[, 1:2])
df$idsGenos <- idsGenos
df$oriGenos <- genos[["inv8"]][df$idsGenos]
```

Create plot

```{r}
vars8 <- pcList$inv8$sdev^2/sum(pcList$inv8$sdev^2)
inv8 <- ggplot(df, aes(x = PC1, y = PC2, color = oriGenos)) +
  # Main points plotting
  geom_point(color = "grey", size = 3, alpha = 0.2) + ggtitle("inv8p23.1") +
  # Plot genotyped
  geom_point(data = df[df$oriGenos %in% c("STD/STD", "INV/INV"), ], aes(color = oriGenos), size = 3) + 
  # Change colors
  scale_color_manual(values = c("#334491", "#349642"), 
                     labels = c("Std", "Inv"), 
                     breaks = c("STD/STD", "INV/INV"), na.value = "black", 
                     name="Validated\nGenotypes") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     legend.position="none")  +
  scale_x_continuous(name = paste(sprintf("PC1 (%.2f", vars8[1]*100), "%)")) +
  scale_y_continuous(name = paste(sprintf("PC2 (%.2f", vars8[2]*100), "%)"))
```

### inv17q21.31

Create data.frame with PCA values and inversion genotypes. These inversion genotypes were obtained with experimental methods.

```{r}
pc <- pcList[["inv17"]]
idsGenos <- idsList[["inv17"]]

df <- data.frame(pc$x[, 1:2])
df$idsGenos <- idsGenos
df$oriGenos <- genos[["inv17"]][df$idsGenos]
```

Create plot

```{r}
vars17 <- pcList$inv17$sdev^2/sum(pcList$inv17$sdev^2)
inv17 <- ggplot(df, aes(x = PC1, y = PC2, color = oriGenos)) +
  # Main points plotting
    geom_point(color = "grey", size = 3, alpha = 0.2) + ggtitle("inv17q21.31") +
  # Plot genotyped
  geom_point(data = df[df$oriGenos %in% c("STD/STD", "INV/INV"), ], aes(color = oriGenos), size = 3) + 
  # Change colors
  scale_color_manual(values = c("#334491", "#349642"), 
                     labels = c("Std", "Inv"), 
                     breaks = c("STD/STD", "INV/INV"), na.value = "black", 
                     name="Validated\nAlleles") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) +
  scale_x_continuous(name = paste(sprintf("PC1 (%.2f", vars17[1]*100), "%)")) +
  scale_y_continuous(name = paste(sprintf("PC2 (%.2f", vars17[2]*100), "%)"))
```

## Drosophila inversions 

Load inversion genotypes.

```{r}
dfInvs <- read_excel("data/inv_genos.xlsx") ## Excel with Drosophila inversion genotypes 
dfInvs <- mutate(dfInvs, `DGRP Line` = gsub("DGRP", "line", `DGRP Line`))
dfInvs <- as.data.frame(dfInvs)

modFold  <- "results/models/"
```

### 2Lt 

Load recombClust results

```{r}
invName <- "In(2L)t"
load(paste0(modFold, invName, "recombRes.Rdata")) ## recombClust results (obtained with runRecombClustDrosophila.R)
```

Create data.frame with PCA values and experimental inversion genotypes.

```{r}
invClass <- dfInvs[, invName]
names(invClass) <- dfInvs[, 1]
  
class2 <- rep(invClass, each = 2)
names(class2) <- paste(rep(names(invClass), each = 2), 1:2, sep = "_")

  
df <- data.frame(recomb$pc$x)
df$class <- class2[rownames(df)]

df <- subset(df, !is.na(class) & class != "INV/ST")
```

Create plot

```{r}
vars <- recomb$pc$sdev^2/sum(recomb$pc$sdev^2)
inv2L <- ggplot(df, aes(x = PC1, y = PC2, color = class)) +
  ggtitle(invName) +
  # Plot genotyped
  geom_point(size = 3) + 
  # Change colors
  scale_color_manual(values = c("#334491", "#349642"), 
                     labels = c("Std", "Inv"), 
                     breaks = c("ST", "INV"), na.value = "black", 
                     guide = FALSE) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) +
  scale_x_continuous(name = paste(sprintf("PC1 (%.2f", vars[1]*100), "%)")) +
  scale_y_continuous(name = paste(sprintf("PC2 (%.2f", vars[2]*100), "%)"))
```

### 2RNS 

Load recombClust results

```{r}
invName <- "In(2R)NS"
load(paste0(modFold, invName, "recombRes.Rdata")) ## recombClust results (obtained with runRecombClustDrosophila.R)
```

Create data.frame with PCA values and experimental inversion genotypes.

```{r}
invClass <- dfInvs[, invName]
names(invClass) <- dfInvs[, 1]

class2 <- rep(invClass, each = 2)
names(class2) <- paste(rep(names(invClass), each = 2), 1:2, sep = "_")

df <- data.frame(recomb$pc$x)
df$class <- class2[rownames(df)]

df <- subset(df, !is.na(class) & class != "INV/ST")
```

Create plot

```{r}
vars <- recomb$pc$sdev^2/sum(recomb$pc$sdev^2)
inv2R <- ggplot(df, aes(x = PC1, y = PC2, color = class)) +
  ggtitle(invName) +
  # Plot genotyped
  geom_point(size = 3) + 
  # Change colors
  scale_color_manual(values = c("#334491", "#349642"), 
                     labels = c("Std", "Inv"), 
                     breaks = c("ST", "INV"), na.value = "black", 
                     guide = FALSE) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) +
  scale_x_continuous(name = paste(sprintf("PC1 (%.2f", vars[1]*100), "%)")) +
  scale_y_continuous(name = paste(sprintf("PC2 (%.2f", vars[2]*100), "%)"))
```

### 3RMo 

Load recombClust results

```{r}
invName <- "In(3R)Mo"
load(paste0(modFold, invName, "recombRes.Rdata")) ## recombClust results (obtained with runRecombClustDrosophila.R)
```

Create data.frame with PCA values and experimental inversion genotypes.

```{r}
invClass <- dfInvs[, invName]
names(invClass) <- dfInvs[, 1]

class2 <- rep(invClass, each = 2)
names(class2) <- paste(rep(names(invClass), each = 2), 1:2, sep = "_")

df <- data.frame(recomb$pc$x)
df$class <- class2[rownames(df)]

df <- subset(df, !is.na(class) & class != "INV/ST")
```

Create plot

```{r}
vars <- recomb$pc$sdev^2/sum(recomb$pc$sdev^2)

inv3R <- ggplot(df, aes(x = PC1, y = PC2, color = class)) +
  ggtitle(invName) +
  # Plot genotyped
  geom_point(size = 3) + 
  # Change colors
  scale_color_manual(values = c("#334491", "#349642"), 
                     labels = c("Std", "Inv"), 
                     breaks = c("ST", "INV"), na.value = "black", 
                     name="Validated\nAlleles") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())  +
  scale_x_continuous(name = paste(sprintf("PC1 (%.2f", vars[1]*100), "%)")) +
  scale_y_continuous(name = paste(sprintf("PC2 (%.2f", vars[2]*100), "%)"))
```

## Combine plots

Combine Drosophila and human plots

```{r}
dros <- plot_grid(inv2L, inv2R, inv3R , ncol = 3, labels = LETTERS[1:3], 
                  rel_widths = c(2, 2, 3))
hum <- plot_grid(inv8, inv17, ncol = 2, labels = c("D", "E"), 
                 rel_widths = c(3, 4))

png("RecombClustInversion.png", width = 20, height = 15, units = "cm", res = 300)
plot_grid(dros, hum, labels = "", ncol = 1)
dev.off()
```

# Sup Fig S8: recombClust All samples  

This plot contains recombClust analysis in inv8p23.1 and inv17q21.31 including all 1000 Genomes individuals.

## Load files

```{r}
load("PCApopulationsLDclassifier.Rdata") ## recombClust results for inv8p23.1 and inv17q21.31 including all 1000Genomes individuals
load("invFestgenos.Rdata") ## Experimental inversion classification 
load("Samples_Pop1GK.Rdata")  ## 1000G individuals' ancestry
genos <- genos[-5] ## Experimental inversion genotypes
names(genos) <- paste0("inv", c(7, "X", 8, 17))
```

Create lists with recombClust PCAs. These objects were obtained with an old version of recombClust. With the release version, both objects can be obtained with the `runRecombClust` function.

```{r}
pcList <- list(inv8 = pc8_mini, inv17 = pc17_mini) ## RecombClust PCAs
idsList <- list(inv8 = ind_ids, inv17 = ind_ids) ## Individuals IDs
```

## Load libraries

```{r}
library(ggplot2)
library(tidyr)
library(cowplot)
```

## Plot inversion status

Plot recombClust results with experimental inversion status

### inv8p23.1

```{r}
pc <- pcList[["inv8"]]
idsGenos <- idsList[["inv8"]]

df <- data.frame(pc[, 1:2])
df$idsGenos <- substring(idsGenos, 1, nchar(idsGenos) - 2)
df$oriGenos <- genos[["inv8"]][df$idsGenos]


recomb8 <- ggplot(df, aes(x = PC1, y = PC2, color = oriGenos)) +
  # Main points plotting
  geom_point(color = "grey") + 
  # Plot genotyped
  geom_point(data = df[df$oriGenos %in% c("STD/STD", "INV/INV"), ], aes(color = oriGenos), size = 3) + 
  # Change colors
  scale_color_manual(values = c("#009E73", "#0072B2"), 
                     labels = c("Std", "Inv"), 
                     breaks = c("STD/STD", "INV/INV"), na.value = "black", 
                     name="Inversion\nStatus") +
  ggtitle("inv8p23.1") +
  theme_bw() +   theme(plot.title = element_text(hjust = 0.5),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank())
```

### inv17q21.31

```{r}
pc <- pcList[["inv17"]]
idsGenos <- idsList[["inv17"]]

df <- data.frame(pc[, 1:2])
df$idsGenos <- substring(idsGenos, 1, nchar(idsGenos) - 2)
df$oriGenos <- genos[["inv17"]][df$idsGenos]

recomb17 <- ggplot(df, aes(x = PC1, y = PC2, color = oriGenos)) +
  # Main points plotting
  geom_point(color = "grey") + 
  # Plot genotyped
  geom_point(data = df[df$oriGenos %in% c("STD/STD", "INV/INV"), ], aes(color = oriGenos), size = 3) + 
  # Change colors
  scale_color_manual(values = c("#009E73", "#0072B2"), 
                     labels = c("Std", "Inv"), 
                     breaks = c("STD/STD", "INV/INV"), na.value = "black", 
                     name="Inversion\nStatus") +
  theme_bw() + ggtitle("inv17q21.31") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
```

## Plot ancestry

Plot ancestry in recombClust results

#### inv8p23.1

```{r}
pc <- pcList[["inv8"]]
idsGenos <- idsList[["inv8"]]
idsGenos <- substring(idsGenos, 1, nchar(idsGenos) - 2)

df <- data.frame(pc[, 1:2])
colnames(df) <- paste0("PC", 1:2)
df$ancestry <- samp_pop[idsGenos, "superpop"]

anc8 <- ggplot(df, aes(x = PC1, y = PC2)) +
  # Main points plotting
  geom_point(color = "grey") + 
  # Plot genotyped
  theme_bw() +
  geom_point(data = df[df$ancestry %in% c("AFR", "EUR", "EAS"),], aes(color = ancestry)) + 
  scale_color_manual(values = c("#E69F00", "#44AA99", "#56B4E9"), 
                     labels = c("African", "East Asian", "European"), 
                     breaks = c("AFR", "EAS", "EUR"),
                     name="Superpopulation") +
  ggtitle("inv8p23.1") +
  theme(plot.title = element_text(hjust = 0.5))
```

### inv17q21.31

```{r}
pc <- pcList[["inv17"]]
idsGenos <- idsList[["inv17"]]
idsGenos <- substring(idsGenos, 1, nchar(idsGenos) - 2)

df <- data.frame(pc[, 1:2])
colnames(df) <- paste0("PC", 1:2)
df$ancestry <- samp_pop[idsGenos, "superpop"]


anc17 <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(color = "grey") + theme_bw() +
  geom_point(data = df[df$ancestry != "SAS",], aes(color = ancestry)) + 
  scale_color_manual(values = c("#E69F00", "#F0E442", "#44AA99", "#56B4E9"), 
                     labels = c("African", "American", "East Asian", "European"), 
                     breaks = c("AFR", "AMR", "EAS", "EUR"),
                     name="Superpopulation")+
  ggtitle("inv17q21.31") +
  theme(plot.title = element_text(hjust = 0.5))
```

## Combined

Combine ancestry and inversion status plots

```{r}
res <- plot_grid(recomb8, anc8,  recomb17, anc17, labels = c("A", "", "B", ""))
png("RecombClustAncestries.jpg",  width = 3000, height = 2500, res = 300)
res
dev.off()
```

# Sup Fig S5: recombClust accuracy 

## Recombination vs mutation

Compare recombClust with invClust results.

### Load data

```{r}
load("simulatedPopulations.Rdata") ## Generated with generatedSyntheticDatasets.R
```


### Process recombClust results

From LD-mixture models, run PCA on recomb probabilities matrix.

```{r}
indsmat <- do.call(cbind, lapply(modelRes, `[[`, "r1"))
rownames(indsmat) <- 1:2000 
pc <- prcomp(indsmat)$x
```

Plot PCA with subpopulation colors.

```{r}
df <- data.frame(pc)
df$Population <- rep(c("A", "B"), each = 1000)

recomb <- ggplot(df, aes(x = PC1, y = PC2, color = Population)) + geom_point() +
  theme_bw() + ggtitle("Recombination") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cbbPalette[c(1, 3)], name = "Chromosomal\nSubpopulation")
```

### invClust results

Plot subpopulations on invClust PCA

```{r}
df <- data.frame(invClustGeno$datin$y)
colnames(df)[1:2] <- c("PC1", "PC2")
df$Population <- rep(c("A", "B"), each = 500)

mut <- ggplot(df, aes(x = PC1, y = PC2, color = Population)) + geom_point() +
  theme_bw() + ggtitle("Mutation") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cbbPalette[c(1, 3)], name = "Chromosomal\nSubpopulation")
```

### Combine

Combine recombClust and invClust plots.

```{r}
pcs <- plot_grid(recomb, mut, labels = "", ncol = 2)
```


## Simulated inversions

### Load data

```{r}
load("simulationRes.Rdata") ## Summary of simulated inversions results (from SimulateInversionsInvertFREGENE.R) 
```

### Summarize data

Plot accuracy per inversion length

```{r}
resdf$lengths <- factor(resdf$lengths, levels = c("50000", "100000", "250000", "500000", "1000000"))

dfbar <- aggregate(resdf$acc, by = list(length = resdf$lengths), mean, na.rm =TRUE)
dfbar$se <- aggregate(resdf$acc, by = list(length = resdf$lengths), 
                      function(x) (sd(x, na.rm =TRUE)/sqrt(length(x))))[, 2]
dfbar$xmin <- dfbar$x - 1.96*dfbar$se
dfbar$xmax <- dfbar$x + 1.96*dfbar$se

len <- ggplot(dfbar, aes(x = length, y = x, ymin = xmin, ymax = xmax)) + 
  geom_bar(stat = "identity", fill = "grey70") +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(limits = c(0, 1), name = "Accuracy") +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_discrete(name = "Inversion Length", 
                   labels = c("50Kb", "100Kb", "250Kb", "500Kb", "1Mb"))

```

## Combine

Combine synthetic population and simulated inversion plots

```{r}
res <- plot_grid(pcs, len, labels = c("A", "B"), ncol = 1)
png("Panel_RecombClust_accuracy.png", res = 300, width = 2485, height = 2665)
res
dev.off()
```

# Sup Fig S6: Inversion frequency

Plot of recombClust accuracy per inversion frequency.

```{r}
load("simulationRes.Rdata") ## Summary of simulated inversions results (from SimulateInversionsInvertFREGENE.R) 

resdf$freqs <- as.factor(resdf$freqs)

png("RecombClust_accuracy_invsFreq_BW.png", res = 300, width = 2000, height = 1100)
ggplot(resdf, aes(x = freqs, y = acc)) + geom_boxplot() +
  scale_y_continuous(limits = c(0.5, 1), name = "Accuracy") +
  theme(legend.position = "none") +
  scale_x_discrete(name = "Inversion Frequency") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
dev.off()
```

# Sup Fig S7: Inversion age

Plot of recombClust accuracy per inversion age.

```{r}
load("simulationRes.Rdata") ## Summary of simulated inversions results (from SimulateInversionsInvertFREGENE.R) 

png("RecombClust_accuracy_invsAge_BW.png", res = 300, width = 2000, height = 1600)
ggplot(resdf, aes(x = age, y = acc)) + geom_point() +
  scale_y_continuous(limits = c(0.5, 1), name = "Accuracy") +
  scale_x_continuous(name = "Inversion Age (in generations)") +
  theme_bw()  + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
dev.off()
```

# Sup Fig S1: Model Evaluation Figures  

Plots with recombClust accuracy in synthetic datasets.

## Load libraries and data

```{r}
library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

load("testProb0.Rdata") ## Synthetic datasets generated with generateSyntheticDatasets.R
```

## Mixture parameter initialization (Sup Fig S1B)

Compute accuracy for initial mixture parameter = 0.95

```{r}
allresHigh <- allresHigh[13:39]
dfHigh <- data.frame(prob0 = unlist(lapply(allresHigh, function(x) 
  sapply(x, function(y) y["prob"]))), 
  scenario = rep(names(allresHigh), lengths(allresHigh)),
  TP = unlist(lapply(allresHigh, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1:1000] > 0.5)))),
  FP = unlist(lapply(allresHigh, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1001:2000] > 0.5)))),
  FN = unlist(lapply(allresHigh, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1:1000] < 0.5)))),
  TN = unlist(lapply(allresHigh, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1001:2000] < 0.5)))))

```

Compute accuracy for initial mixture parameter = 0.05

```{r}
allresLow <- allresLow[13:39]
dfLow <- data.frame(prob0 = unlist(lapply(allresLow, function(x) 
  sapply(x, function(y) y["prob"]))), 
  scenario = rep(names(allresLow), lengths(allresLow)),
  TP = unlist(lapply(allresLow, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1:1000] > 0.5)))),
  FP = unlist(lapply(allresLow, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1001:2000] > 0.5)))),
  FN = unlist(lapply(allresLow, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1:1000] < 0.5)))),
  TN = unlist(lapply(allresLow, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1001:2000] < 0.5)))))
```


Create plot

```{r}
dfAll <- rbind(dfLow, dfHigh)
dfAll$p0 <- rep(c(0.05, 0.95), each = nrow(dfLow))
dfAll$SNlin <- dfAll$TN/(dfAll$TN + dfAll$FP)
dfAll$SNrec <- dfAll$TP/(dfAll$TP + dfAll$FN)
dfAll$Acc <- (dfAll$TP +  dfAll$TN)/(dfAll$TN + dfAll$FN + dfAll$TP + dfAll$FP)
dfAlltide <- gather(dfAll, "param", "value", 8:10)

prob0 <- dfAlltide %>% 
  group_by(param, p0) %>%
  summarize(m = mean(value), sd = sd(value), n = n()) %>%
  mutate(se = sd/sqrt(n)) %>%
  ggplot(aes(x = p0, y = m, color = param)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.1) +
  geom_line() +
  geom_point() + 
  scale_y_continuous(name = "Accuracy", limits = c(0, 1)) +
  scale_x_continuous(name = expression(pi[0])) +
  scale_color_discrete(name = "", 
                       labels = c("Global", "linkage", "recomb"),
                       breaks = c("Acc", "SNlin", "SNrec")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16), 
        panel.grid = element_blank(), 
        strip.background = element_rect(fill="white"))
```

  
## Intrablock Linakge (Sup Fig S1D)

Compute accuracy per number of SNP-blocks in high linkage.

```{r}
dfAlltide$combs <- gsub("NoLD_|R2_", "", dfAll$scenario)
Rg <- gregexpr("R", dfAlltide$combs)
Rgr <- regexpr("R", dfAlltide$combs)
dfAlltide$R <- lengths(Rg)
dfAlltide$R[dfAlltide$R == 1 & Rgr == -1] <- 0
```

Make plot

```{r}
LD <- dfAlltide %>% 
  group_by(param, R) %>%
  summarize(m = mean(value), sd = sd(value), n = n()) %>%
  mutate(se = sd/sqrt(n)) %>%
  ggplot(aes(x = R, y = m, color = param)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.1) +
  geom_line() +
  geom_point() + 
  scale_y_continuous(name = "Accuracy", limits = c(0, 1)) +
  scale_x_continuous(name = "Full Linkage Blocks") +
  scale_color_discrete(name = "", 
                       labels = c("Global", "linkage", "recomb"),
                       breaks = c("Acc", "SNlin", "SNrec")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16), 
        panel.grid = element_blank(), 
        strip.background = element_rect(fill="white"))
```



## Genetic Distance (Sup Fig S1C)

Compute accuracy for genetic distance (number of SNPs with different major allele between the subpopulations).

```{r}
load("testAllele.Rdata") ## Synthetic datasets generated with generatedSyntheticDatasets.R
 
allresAllele <- allresAllele[18:33]

dfAllele <- data.frame( 
  TP = unlist(lapply(allresAllele, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1:1000] > 0.5)))),
  FP = unlist(lapply(allresAllele, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1001:2000] > 0.5)))),
  FN = unlist(lapply(allresAllele, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1:1000] < 0.5)))),
  TN = unlist(lapply(allresAllele, function(x) 
    sapply(x, function(y) sum(y[["r1"]][1001:2000] < 0.5)))),
  scenario = rep(names(allresAllele), lengths(allresAllele)))


dfAllele$SNlin <- dfAllele$TN/(dfAllele$TN + dfAllele$FP)
dfAllele$SNrec <- dfAllele$TP/(dfAllele$TP + dfAllele$FN)
dfAllele$Acc <- (dfAllele$TP +  dfAllele$TN)/(dfAllele$TN + dfAllele$FN + dfAllele$TP + dfAllele$FP)

lowe <- gregexpr("[a-d]", dfAllele$scenario)
Rgr <- regexpr("[a-d]", dfAllele$scenario)
dfAllele$lower <- lengths(lowe)
dfAllele$lower[dfAllele$lower == 1 & Rgr == -1] <- 0


dfAlleletide <- gather(dfAllele, "param", "value", 6:8)

lowe <- gregexpr("[a-d]", dfAlleletide$scenario)
Rgr <- regexpr("[a-d]", dfAlleletide$scenario)
dfAlleletide$lower <- lengths(lowe)
dfAlleletide$lower[dfAlleletide$lower == 1 & Rgr == -1] <- 0

```

Create plot

```{r}
freq <- dfAlleletide %>% 
  group_by(param, lower) %>%
  summarize(m = mean(value), sd = sd(value), n = n()) %>%
  mutate(se = sd/sqrt(n)) %>%
  ggplot(aes(x = lower, y = m, color = param)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.1) +
  geom_line() +
  geom_point() + 
  scale_y_continuous(name = "Accuracy", limits = c(0, 1)) +
  scale_x_continuous(name = "Different Major Alleles") +
  scale_color_discrete(name = "", 
                       labels = c("Global", "linkage", "recomb"),
                       breaks = c("Acc", "SNlin", "SNrec")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16), 
        panel.grid = element_blank(), 
        strip.background = element_rect(fill="white"))
```

## Mixture proportions (Sup Fig S1A)

Compute accuracy per mixture proportion

```{r}
load("testCombs.Rdata")  ## Synthetic datasets generated with generateSyntheticDatasets.R
props <- seq(0.1, 0.9, 0.1)

getInd <- function(prop){
  if (prop > 0.5){
    prop <- 1 - prop
    ind <- prop*1000/(1 - prop)
  } else {
    ind <- 1000
  }
}

dfProps <- data.frame(
  TP =  c(unlist(lapply(seq_len(length(props)), function(i) 
    sapply(allresCombsProps[[2]][[i]], 
           function(y) sum(y[["r1"]][1:getInd(props[i])] > 0.5)
    ))), unlist(lapply(seq_len(length(props)), function(i) 
      sapply(allresCombsProps[[3]][[i]], 
             function(y) sum(y[["r1"]][-c(1:getInd(props[i]))] > 0.5))))), 
  FP =  c(unlist(lapply(seq_len(length(props)), function(i) 
    sapply(allresCombsProps[[2]][[i]], 
           function(y) sum(y[["r1"]][-c(1:getInd(props[i]))] > 0.5)
    ))), unlist(lapply(seq_len(length(props)), function(i) 
      sapply(allresCombsProps[[3]][[i]], 
             function(y) sum(y[["r1"]][1:getInd(props[i])] > 0.5))))), 
  TN =  c(unlist(lapply(seq_len(length(props)), function(i) 
    sapply(allresCombsProps[[2]][[i]], 
           function(y) sum(y[["r1"]][-c(1:getInd(props[i]))] < 0.5)
    ))), unlist(lapply(seq_len(length(props)), function(i) 
      sapply(allresCombsProps[[3]][[i]], 
             function(y) sum(y[["r1"]][1:getInd(props[i])] < 0.5))))), 
  FN =  c(unlist(lapply(seq_len(length(props)), function(i) 
    sapply(allresCombsProps[[2]][[i]], 
           function(y) sum(y[["r1"]][1:getInd(props[i])] < 0.5)
    ))), unlist(lapply(seq_len(length(props)), function(i) 
      sapply(allresCombsProps[[3]][[i]], 
             function(y) sum(y[["r1"]][-c(1:getInd(props[i]))] < 0.5))))),
  scenario = rep(names(allresCombsProps[2:3]), each = 200*9),
  prop = rep(c(props, 1 - props), each = 200))
dfProps$Proportion  <- as.factor(dfProps$prop)

dfProps$SNlin <- dfProps$TN/(dfProps$TN + dfProps$FP)
dfProps$SNrec <- dfProps$TP/(dfProps$TP + dfProps$FN)
dfProps$Acc <- (dfProps$TP +  dfProps$TN)/(dfProps$TN + dfProps$FN + dfProps$TP + dfProps$FP)
dfPropstide <- gather(dfProps, "param", "value", 8:10)
```

Create plot

```{r}
mixt <- dfPropstide %>% 
  group_by(param, Proportion) %>%
  summarize(m = mean(value), sd = sd(value), n = n()) %>%
  mutate(se = sd/sqrt(n)) %>%
  ggplot(aes(x = 1 - as.numeric(as.character(Proportion)), y = m, color = param)) + 
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.1) +
  geom_line() +
  geom_point() + 
  scale_y_continuous(name = "Accuracy", limits = c(0, 1)) +
  scale_x_continuous(name = expression(pi), 
                     breaks = round(seq(0, 1, 0.2), 1)) +
  scale_color_discrete(name = "", 
                       labels = c("Global", "linkage", "recomb"),
                       breaks = c("Acc", "SNlin", "SNrec")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16), 
        panel.grid = element_blank(), 
        strip.background = element_rect(fill="white"))
```

## Combine plots

```{r}
jpeg("model_evaluation.jpeg",  width = 2100, height = 1700, res = 300)
plot_grid(mixt, prob0, freq, LD, ncol = 2, labels = LETTERS[1:4])
dev.off()
```

# Sup Figs S3 and S4: Evaluate recombClust in simulated populations 

## Load libraries and data

```{r}
library(cluster)
library(dplyr)
library(parallel)
library(ggplot2)
library(cowplot)
library(tidyr)

load("simulations/SimBlockRecombResBlocks.Rdata") ## Simulated populations generated with simulatePopulations.Rmd
load("simulations/SimBlockRecombRes.Rdata") ## Simulated populations generated with simulatePopulations.Rmd
```


## Compute accuracy 

Compute accuracy per number of blocks or number of individuals in each subpopulation

```{r}
nInds <- c(10, 15, 20, 25, 30)
nBlocks <- c(10, 15, 20, 40, 70, 100) ## SNP block pairs per populations
nsim <- 1000
```

Define functions

```{r}
getClust <- function(recObj){
  indsmat <- do.call(cbind, lapply(recObj, `[[`, "r1"))
  pc <- prcomp(indsmat, rank. = 2)$x
  class <- kmeans(pc, centers = 2, nstart = 100)$cluster
  list(pc = pc, class = class)
}
computeSilhouette <- function(l){
  summary(silhouette(l$class, dist(l$pc)))$si.summary[4]
}
```

### Chromosome

Silhouette plots for different number of chromosomes in each subpopulation.

```{r}
mixPCs <- lapply(mixRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})
singlePCs <- lapply(singleRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})

mixSilhouete <- lapply(mixPCs, function(x) sapply(x, computeSilhouette))
singleSilhouete <- lapply(singlePCs, function(x) sapply(x, computeSilhouette))
indsDF <- data.frame(sil = c(unlist(mixSilhouete), unlist(singleSilhouete)),
                     inds = rep(nInds*2, each = nsim),
                     type = rep(c("Mixture", "Single"), each = nsim*length(nInds)))
```

Create plots

```{r}
indsPlot <- indsDF %>%
  group_by(inds, type) %>%
  summarize(mean = mean(sil),
            sd = sd(sil),
            se = sd/sqrt(n())) %>%
  ggplot(aes(x = inds, color = type, y = mean)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  geom_hline(yintercept = c(0.5, 0.7), linetype = "dashed") +
  scale_color_discrete(name = "Dataset") +
  scale_y_continuous(name = "Average silhouette", limits = c(0.25, 1)) + 
  scale_x_continuous(name = "Chromosomes")

indsEstsplot <- indsDF %>%
  mutate(sel = sil > 0.7) %>%
  group_by(inds, type) %>%
  summarize(mean = mean(sel)) %>%
  mutate(est = ifelse(type == "Mixture", "Power", "FPR")) %>%
  ggplot(aes(x = inds, color = est, y = mean)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  scale_color_discrete(name = "") +
  scale_x_continuous(name = "Chromosomes") +
  scale_y_continuous(name = "")
```

### SNP-blocks 

Silhouette plots for different number of SNP-blocks 

```{r}
mixBlockPCs <- lapply(mixBlockRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})
singleBlockPCs <- lapply(singleBlockRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})

mixBlockSilhouete <- lapply(mixBlockPCs, function(x) sapply(x, computeSilhouette))
singleBlockSilhouete <- lapply(singleBlockPCs, function(x) sapply(x, computeSilhouette))

blockDF <- data.frame(sil = c(unlist(mixBlockSilhouete), unlist(singleBlockSilhouete)),
                     blocks = rep(nBlocks, each = nsim),
                     type = rep(c("Mixture", "Single"), each = nsim*length(nBlocks)))
```

Make plots

```{r}
blockPlot <- blockDF %>%
  group_by(blocks, type) %>%
  summarize(mean = mean(sil),
            sd = sd(sil),
            se = sd/sqrt(n())) %>%
  ggplot(aes(x = blocks, color = type, y = mean)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean + se)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  geom_hline(yintercept = c(0.5, 0.7), linetype = "dashed") +
  scale_color_discrete(name = "Dataset") +
  scale_y_continuous(name = "Average silhouette", limits = c(0.25, 1)) + 
  scale_x_continuous(name = "Recombination points")

blockEstsplot <- blockDF %>%
  mutate(sel = sil > 0.7) %>%
  group_by(blocks, type) %>%
  summarize(mean = mean(sel)) %>%
  mutate(est = ifelse(type == "Mixture", "Power", "FPR")) %>%
  ggplot(aes(x = blocks, color = est, y = mean)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  scale_color_discrete(name = "") +
  scale_x_continuous(name = "Recombination points") +
  scale_y_continuous(name = "")
```

## Combine plots

```{r}
png("simulations/populationSimulations.png",  width = 2100, height = 1700, res = 300)
plot_grid(indsPlot, blockPlot, ncol = 1, labels = c("A", "B"))
dev.off()

png("simulations/populationSimulationsEstimates.png",  width = 2100, height = 1700, res = 300)
plot_grid(indsEstsplot, blockEstsplot, ncol = 1, labels = c("A", "B"))
dev.off()
```

# Figure 6: Compare recombClust with SNPRelate 

## Load library and data

```{r}
load("SNPRelate/results/plotFiles/SNPrelate_results.Rdata") ## SNPrelate results generated with runSNPrelate.Rmd

library(SNPRelate)
```

Create plots with SNPRelate results and displaying recombClust classification.

```{r}
levs <- c("1", "2", "1-1", "1-2", "2-2")

png("SNPRelate/results/plotFiles/SNPrelate.png",  width = 2100, height = 1700, res = 300)
layout(matrix(c(rep(1:3, each = 2), rep(4:7, each = 3)), ncol=6, byrow=TRUE))
plot(drosPCAs[[1]], 
     col = factor(drosClass$In2L[drosPCAs[[1]]$sample.id], levels = levs),
     main = "In(2L)t", 
     pch = 16)
plot(drosPCAs[[2]], 
     col = factor(drosClass$In2R[drosPCAs[[2]]$sample.id], levels = levs), 
     main = "In(2R)NS", 
     pch = 16)
plot(drosPCAs[[3]], 
     col = factor(drosClass$In3R[drosPCAs[[3]]$sample.id], levels = levs), 
     main = "In(3R)Mo", 
     pch = 16)
plot(genoPCAs[[3]],
     col = factor(humanClass$inv8[genoPCAs[[3]]$sample.id], levels = levs),
     main = "inv8p23.1", 
     pch = 16)
plot(genoPCAs[[4]],
     col = factor(humanClass$inv17[genoPCAs[[4]]$sample.id], levels = levs),
     main = "inv17q21.31", 
     pch = 16)
plot(genoPCAs[[2]],
     col = factor(humanClass$lct[genoPCAs[[2]]$sample.id], levels = levs), 
     main = "LCT region", 
     pch = 16)
plot(genoPCAs[[1]], 
     col = factor(humanClass$tar[genoPCAs[[1]]$sample.id], levels = levs),
     main = "TAR region", 
     pch = 16)
dev.off()
```
