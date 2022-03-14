load("hpc/inputs/epimutations_found_GSE51032.rda")
subject  <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(found_GSE51032), seqnames.field = "chr", start.field = "start", end.field= "end")
epimutations <- found_GSE51032[, c("chr", "start", "end")]
epimutations$epi_id <- "grag"
colnames(epimutations) <- c("chromosome", "start", "end", "epi_id")

load("hpc/outputs/GSE51032/GSE51032beta.rda")
print("beta")
nrow(unique(res[, c("chromosome", "start", "end")]))
samples_beta <- res$sample
res <- unique(res[, c("chromosome", "start", "end")])
query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(res), seqnames.field = "chromosome", start.field = "start", end.field= "end")
hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")
query[S4Vectors::queryHits(hits)]
epimutations <- rbind(epimutations, res[, c("chromosome", "start", "end", "epi_id")])

load("hpc/outputs/GSE51032/GSE51032mahdistmcd.rda")
print("mahdistmcd")
samples_maha <- res$sample
nrow(unique(res[, c("chromosome", "start", "end")]))
res <- unique(res[, c("chromosome", "start", "end")])
query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(res), seqnames.field = "chromosome", start.field = "start", end.field= "end")
hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")
query[S4Vectors::queryHits(hits)]
epimutations <- rbind(epimutations, res[, c("chromosome", "start", "end", "epi_id")])

load("hpc/outputs/GSE51032/GSE51032manova.rda")
print("manova")
samples_manova <- res$sample

nrow(unique(res[, c("chromosome", "start", "end")]))
res <- unique(res[, c("chromosome", "start", "end")])
query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(res), seqnames.field = "chromosome", start.field = "start", end.field= "end")
hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")
query[S4Vectors::queryHits(hits)]
nrow(unique(res[, c("chromosome", "start", "end")]))
epimutations <- rbind(epimutations, res[, c("chromosome", "start", "end", "epi_id")])

load("hpc/outputs/GSE51032/GSE51032mlm.rda")
print("mlm")
samples_mlm <- res$sample
nrow(unique(res[, c("chromosome", "start", "end")]))
res <- unique(res[, c("chromosome", "start", "end")])
query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(res), seqnames.field = "chromosome", start.field = "start", end.field= "end")
hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")
query[S4Vectors::queryHits(hits)]
nrow(unique(res[, c("chromosome", "start", "end")]))
epimutations <- rbind(epimutations, res[, c("chromosome", "start", "end", "epi_id")])

load("hpc/outputs/GSE51032/GSE51032quantile.rda")
print("quantile")
samples_quantile <- res$sample

nrow(unique(res[, c("chromosome", "start", "end")]))
res <- unique(res[, c("chromosome", "start", "end")])
subject  <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(found_GSE51032), seqnames.field = "chr", start.field = "start", end.field= "end")
query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(res), seqnames.field = "chromosome", start.field = "start", end.field= "end")
hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")
query[S4Vectors::queryHits(hits)]
epimutations <- rbind(epimutations, res[, c("chromosome", "start", "end", "epi_id")])

load("hpc/outputs/GSE51032/GSE51032isoforest.rda")
print("isoforest")
samples_isoforest <- res$sample
nrow(unique(res[, c("chromosome", "start", "end")]))
res <- unique(res[, c("chromosome", "start", "end")])
subject  <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(found_GSE51032), seqnames.field = "chr", start.field = "start", end.field= "end")
query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(res), seqnames.field = "chromosome", start.field = "start", end.field= "end")
hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")
query[S4Vectors::queryHits(hits)]
epimutations <- rbind(epimutations, res[, c("chromosome", "start", "end", "epi_id")])


