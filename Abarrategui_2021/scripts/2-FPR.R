#1. Function to calculate the FPR 
FPR <- function(file_name){
#1.1. Load the data file containing the results list
  load(paste0("result_files/1-Result_list/", file_name, ".rda"))
#1.2. Set the parameters  
  control_n <- seq(20, 100, 10)
  methods <- c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "beta")
#1.3. Calculate the FPR 
  FPR <- do.call(list, lapply(seq(length(result_table)), function(n){
    rst <- do.call(rbind, lapply(seq_len(length(methods)), function(i){
      results <- grepl(methods[i],unlist(result_table[[n]]["epi_id"]))
      results <- result_table[[n]][results,]
      if(!all(is.na(results$pvalue)) & methods[i] == "manova" | methods[i] == "mlm"){
        results <- results[results$pvalue < 0.05/40408,]
      }else if(!all(is.na(results$pvalue)) & methods[i] == "isoforest"){
        results <- results[results$outlier_score > 0.7,]
      }
      if(all(is.na(results$cpg_n))){
        FPR <- 0
      }else{
        keep <- which(!is.na(results$cpg_n))
        results <- results[keep,]
        FPR <- nrow(unique(results))/400
      }
      
      table <- c(methods[i], as.numeric(control_n[n]), as.numeric(FPR))
      table
      
    }))
    rst <- as.data.frame(rst)
    colnames(rst) <- c("method", "n", "FPR")
    rst
  }))
  names(FPR) <-  paste0("n", seq(from = 20, to = 100, by = 10))
  
#1.4. Save results
save(FPR, file = paste0("C:result_files/2-FPR/", file_name, ".rda"))
}

# 2. Create FPR tables
FPR("GSE51032-FPR")
FPR("GSE111629-FPR")




