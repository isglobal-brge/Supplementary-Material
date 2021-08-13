n <- Sys.getenv("N")
i <- Sys.getenv("I")

print(paste0("n=", n, ", i=", i))

# 1. Generate the function
TPR <- function(file_name, case_sample, methods,  chr, start, end){
  #1.1 Load the results
  input.filename <- paste0("inputs/",file_name,".rda")
  variable <- load(file=input.filename)
  object   <- eval(parse(text=variable))
  #1.2 Filtering inputs
  control_samples  <- which(object$status == "control")
  case_samples <- object[,case_sample]
  #1.3 Epimutations detection
  results <- do.call(rbind, lapply(seq_len(length(methods)), function(j) {
    rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii){
    samples_ncol <- sample(control_samples, size = n, replace = FALSE)
    control_panel <- GSE51032[,samples_ncol]
    time <- system.time(res <- epimutacions::epimutations(case_samples, control_panel,
                                                          method = methods[j], 
                                                          chr = chr[ii],
                                                          start = start[ii],
                                                          end = end[ii]))
    res$time <- time[[2]]
    res
    
    }))
  }))
  
  #1.4. Save the results
  output.filename <- paste0("outputs/outputfile-","TPR-", variable,"-", n, "-", i, ".rda")
  save(results, file=output.filename)
  print(paste0("Written to file '", output.filename, "'"))
}

# 2. Calculate TPR for each cohort

## 2.1 GSE51032

TPR(file_name = "GSE51032", 
    case_sample = "GSM1235784", 
    methods = c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta"),  
    chr = "chr17", start = 46000000, end = 46020000)

## 2.2. GSE111629

TPR(file_name = "GSE111629", 
    case_sample = c("GSM3035933", "GSM3035791", "GSM3035807","GSM3035685"), 
    methods = c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta"),  
    chr = c("chr5", "chr5", "chr17","chr19"), 
    start = c(10240000, 67580000, 46000000, 11190000), 
    end = c(10260000, 67600000, 46020000, 11210000))



