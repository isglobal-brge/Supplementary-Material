library(epimutacions)

#1. Load the data
n <- Sys.getenv("N")
i <- Sys.getenv("I")

print(paste0("n=", n, ", i=", i)) 

input.filename <- paste0("inputs/GSE111629",".rda")
variable <- load(file=input.filename)
object   <- eval(parse(text=variable))
print(object)

load("inputs/epi_validated.rda")

#2. Filtering inputs
case_names <- c("GSM3035933", "GSM3035791", "GSM3035807","GSM3035685")
control_samples <- which(GSE111629$status == "control")
methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "beta")

#3. Simulation
rst <- do.call(rbind, lapply(seq_len(length(methods)), function(j){
    samples_ncol <- sample(control_samples, size = n, replace = FALSE)
    control_panel <- GSE111629[,samples_ncol]
    case_samples <- GSE111629[,case_names]
    epimutacions::epimutations(case_samples, control_panel, 
                               method = methods[j], 
                               chr = c("chr5", "chr5", "chr17","chr19"),
                               start = c(10240000, 67580000, 46000000, 11190000), 
                               end = c(10260000, 67600000, 46020000, 11210000))
}))

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
    colnames(results)[c(7:15)]  <- c("chromosome_validated","start_validated","end_validated", "sz_validated", "25%", "50%", "75%", "100%", "%")
  }else{
    results <- NA
  }


#4. Save the results
output.filename <- paste0("outputs/outputfile-", variable,"-", n, "-", i, ".rda")
results.name <- paste0("outputfile-","TPR-", variable,"-", n, "-", i)
save(results, file=output.filename)
print(paste0("Written to file '", output.filename, "'"))


