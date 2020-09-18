#'#################################################################################
#'#################################################################################
#' Helper code to convert invertFREGENE data to plink
#'#################################################################################
#'#################################################################################
args = commandArgs(trailingOnly=TRUE)

root <- args[1]

library(snpStats)
genosraw <- data.matrix(read.table(paste0(root, ".genotypes_0.dat")))
genos <- genosraw[-1, ]
rownames(genos) <- 1:nrow(genos)
annot <- data.frame(chromosome = 1, position = genosraw[1, ], 
                    name = genosraw[1, ], 
                    stringsAsFactors = FALSE)
rownames(annot) <- annot$name
colnames(genos) <- rownames(annot)

snp <- new("SnpMatrix", genos + 1)
write.plink(paste0(root, ".genotypesplink"), snps = snp, chromosome = 1, position = annot$position)
