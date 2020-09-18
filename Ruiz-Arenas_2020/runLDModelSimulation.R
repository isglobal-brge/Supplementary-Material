## Code to run LDclassifier to simulated inversions

args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
dist <- as.numeric(args[2])

# Load libraries
library(LDclassifier)
library(GenomicRanges)

## Load data
haplosraw <- data.matrix(read.table(paste0(root, ".haplotypes_0.dat")))
haplos <- haplosraw[-1, ]
rownames(haplos) <- 1:nrow(haplos)
annot <- data.frame(chromosome = 1, position = haplosraw[1, ], 
                    name = haplosraw[1, ], 
                    stringsAsFactors = FALSE)
rownames(annot) <- annot$name
colnames(haplos) <- rownames(annot)

## Read inv ids 
#### Ids are coded in a 0-index base!!!!
inv_ids <- read.table(paste0(root, ".log_ids.txt"))
ids <- rep(0, 2000)
ids[inv_ids$V1 + 1] <- 1

## Filter SNPs with MAF < 0.1
haplosfilt <- haplos[, colMeans(haplos) > 0.1 & colMeans(haplos) < 0.9]

## Select SNPs inside the inversion region 
haplosfilt <- haplosfilt[, as.numeric(colnames(haplosfilt)) > 750000 & as.numeric(colnames(haplosfilt)) < 750000 + dist]
annotfilt <- annot[as.character(colnames(haplosfilt)), ]
GRfilt <- makeGRangesFromDataFrame(annotfilt, start.field = "position", 
                                     end.field = "position")
message("Running models")
models <- runLDclassifier(haplosfilt, GRfilt, 
                            BPPARAM = BiocParallel::MulticoreParam(20))
indsmat <- do.call(cbind, lapply(models, `[[`, "r1"))
pc <- prcomp(indsmat)
class <- kmeans(pc$x[, 1], centers = 2, nstart = 1000)$cluster
save(pc, class, ids, file = paste0(root, ".modelres.Rdata"))
