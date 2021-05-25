library(epimutacions)

#1. Load the data
n <- Sys.getenv("N")
i <- Sys.getenv("I")

print(paste0("n=", n, ", i=", i)) 

input.filename <- paste0("inputs/GSE51032",".rda")
variable <- load(file=input.filename)
object   <- eval(parse(text=variable))
print(object)

load("inputs/epi_validated.rda")
methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "beta")

#2. Filtering inputs
control_samples  <- which(GSE51032$status == "control")
case_samples <- GSE51032[,"GSM1235784"]

#3. Simulation
rst <- do.call(rbind, lapply(seq_len(length(methods)), function(j) {
  rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii){
    samples_ncol <- sample(control_samples, size = n, replace = FALSE)
    control_panel <- GSE51032[,samples_ncol]
    epimutacions::epimutations(case_samples, control_panel,
                               method = methods[j], 
                               chr = "chr17",
                               start = 46000000,
                               end = 46020000)
    
  }))
}))


if(is.null(rst)){
  results <- NA
}else{
  query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(rst), seqnames.field = "chromosome", start.field = "start", end.field= "end")
  subject  <- epi_validated
  hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")

  if(length(hits) != 0){
    overlaps <- GenomicRanges::pintersect(query[S4Vectors::queryHits(hits)], subject[S4Vectors::subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(subject[S4Vectors::subjectHits(hits)])
  
  
  
    keep_rst <- S4Vectors::queryHits(hits)
    keep_validated <- S4Vectors::subjectHits(hits)
    rst <- as.data.frame(rst[keep_rst,])
    validated <- as.data.frame(epi_validated[keep_validated,])
    results <- cbind(rst,validated[,1:4])
    results$percent25 <- ifelse(percentOverlap >= 0.25, 1, 0)
    results$percent50 <- ifelse(percentOverlap >= 0.50, 1, 0)
    results$percent75 <- ifelse(percentOverlap >= 0.75, 1, 0)
    results$percent100 <- ifelse(percentOverlap >= 1, 1, 0)
    results$percent <- percentOverlap
  
    results <- results[,c(1:6, 16:24, 7:15)]
    colnames(results)[c(7:15)] <- c("chromosome_validated","start_validated","end_validated", "sz_validated", "25%", "50%", "75%", "100%", "%")
  }else{
    results <- NA
  }
}

#4. Save the results
output.filename <- paste0("outputs/outputfile-","TPR-", variable,"-", n, "-", i, ".rda")
save(results, file=output.filename)
print(paste0("Written to file '", output.filename, "'"))


