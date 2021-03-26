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

#2. Filtering inputs
case_samples <- object[,"GSM1235784"]
control_samples <- object[,GSE51032$cancer_type == "free"]

total_control_samples <- ncol(control_samples)
samples_ncol <- sample.int(total_control_samples, size = n, replace = FALSE)
control_samples <- control_samples[,samples_ncol] 
methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn")

#3. Simulation
rst <- do.call(rbind, lapply(seq_len(length(methods)), function(j) {
  rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii) {
    epimutations(case_samples[,ii], 
                 control_panel = control_samples, 
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
  
  results <- results[,c(1:6, 13:16, 7:8, 17:21, 9:12)]
  colnames(results)[c(7:10,13:17)]  <- c("chromosome_validated","start_validated","end_validated", "sz_validated", "25%", "50%", "75%", "100%", "%")
  }else{
    results <- NA
  }
}

#4. Save the results
output.filename <- paste0("outputs/outputfile-", variable,"-", n, "-", i, ".rda")
save(results, file=output.filename)
print(paste0("Written to file '", output.filename, "'"))


