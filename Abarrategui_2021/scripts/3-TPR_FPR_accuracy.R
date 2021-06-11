# 1. Function to create TPR, FPR and accuracy list
TPR_FPR_accuracy <- function(FPR_file_name, TPR_file_name, final_file_name){
# 1.1. Load the data files containing TPR, FPR and accuracy list 
  load(paste0("result_files/2-TPR_accuracy/", TPR_file_name, ".rda"))
  load(paste0("result_files/2-FPR/", FPR_file_name, ".rda"))
# 1.2. Create the table
  df <- do.call(list, lapply(seq_len(length(TPR)), function(i) {
    TPR_n <- as.data.frame(data.table::rbindlist(TPR[[i]]))
    FPR <- as.data.frame(data.table::rbindlist(FPR))
    out <- merge(TPR_n, FPR, by = c("method", "n"), sort = FALSE)
    out[,3:5] <- apply(out[,3:5], 2, function(x) as.numeric(x))
    out[,3:5] <- out[,3:5] * 100
    out
  }))
  names(df) <- names(TPR)
# 1.3. Save results  
  save(df, file = paste0("result_files/3-TPR_FPR_accuracy/",final_file_name, ".rda"))
}

# 2. Create tables
TPR_FPR_accuracy("GSE111629-FPR", "GSE111629-TPR", "TPR_FPR_accuracy_GSE111629")
TPR_FPR_accuracy("GSE51032-FPR", "GSE51032-TPR", "TPR_FPR_accuracy_GSE51032")
