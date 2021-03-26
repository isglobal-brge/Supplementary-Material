load.dir<- "result_files/hpc_results_list"
file_name <- "GSE111629_out_chr_start_end"
load(paste0(load.dir, "/", file_name, ".rda"))


control_n <- seq(20, 100, 10)
simulation <- result_table
methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn")
FPR <- do.call(list, lapply(seq(length(simulation)), function(n){
  rst <- do.call(rbind, lapply(seq_len(length(methods)), function(i){
    results <- grepl(methods[i],unlist(simulation[[n]]["epi_id"]))
    results <- simulation[[n]][results,]
    if(methods[i] == "manova" | methods[i] == "mlm"){
      results <- results[results$pvalue < 0.05/40408,]
    }else if(methods[i] == "isoforest"){
      results <- results[results$outlier_score > 0.7,]
    }
    result <- nrow(unique(results[,c("sample","i")]))/400
    table <- c(methods[i], as.numeric(control_n[n]), as.numeric(result))
    table
    
    }))
  rst <- as.data.frame(rst)
  colnames(rst) <- c("method", "n", "FPR")
  rst
}))
names(FPR) <-  paste0("n", seq(from = 20, to = 100, by = 10))

save.dir <- "C:result_files/FPR/"
save(FPR, file = paste0(save.dir, "FPR_", file_name, ".rda"))
rm(list = ls())
