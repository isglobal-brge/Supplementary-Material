load.dir<- "result_files/hpc_results_list"
file_name <- "GSE51032_in_chr"
load(paste0(load.dir, "/", file_name, ".rda"))


simulation <- result_table
control_n <- seq(20, 100, 10)
methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn")
result <- do.call(list, lapply(seq(length(simulation)), function(n){
  rst <-  do.call(rbind, lapply(seq_len(length(methods)), function(i){
    results <- grepl(methods[i],unlist(simulation[[n]]["epi_id"]))
    results <- simulation[[n]][results,]
    if(methods[i] == "manova" | methods[i] == "mlm"){
      results <- results[results$pvalue < 0.05/40408,]
    }else if(methods[i] == "isoforest"){
      results <- results[results$outlier_score > 0.7,]
    }
    TPR <- nrow(results)/100
    percent <- mean(results[,17])
    table <- c(methods[i], as.numeric(control_n[n]), as.numeric(TPR), round(as.numeric(percent), digits = 3))
  }))
  rst <- as.data.frame(rst)
  colnames(rst) <- c("method", "n", "TPR", "accuracy")
  rst
}))
TPR <- result
names(TPR) <- paste0("n", seq(from = 20, to = 100, by = 10))
save.dir <- "result_files/TPR/"
save(TPR, file = paste0(save.dir, "TPR_", file_name, ".rda"))
rm(list =ls())

#############################################################################################
load.dir<- "result_files/hpc_results_list"
file_name <- "GSE111629_in_chr_start_end"
load(paste0(load.dir, "/", file_name, ".rda"))
load("hpc/inputs/epi_validated.rda")

epi_validated <- as.data.frame(epi_validated)
number <- c(100, 100, 100, 200)
control_n <- seq(20, 100, 10)

methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn")
regions <-  do.call(list, lapply(seq(nrow(epi_validated)), function(i){
  TPR <- do.call(list, lapply(seq(length(result_table)), function(n){
    rst <-  do.call(rbind, lapply(seq_len(length(methods)), function(j){
      results <- result_table[[n]]
      regions <- which(results$chromosome_validated == epi_validated$seqnames[i] & results$start_validated == epi_validated$start[i])
      regions <- results[regions,]
      keep_method <- grepl(methods[j],unlist(regions[,"epi_id"]))
      regions <- regions[keep_method,]
      if(methods[j] == "manova" | methods[j] == "mlm"){
        regions <- regions[regions$pvalue < 0.05/40408,]
      }else if(methods[j] == "isoforest"){
        regions <- regions[regions$outlier_score > 0.7,]
      }
      TPR <- nrow(regions)/number[i]
      percent <- mean(regions[,17])
      table <- c(methods[j], as.numeric(control_n[n]), as.numeric(TPR), round(as.numeric(percent), digits = 3))
    }))
    rst <- as.data.frame(rst)
    rownames(rst)
    colnames(rst) <- c("method", "n", "TPR", "accuracy")
    rst
  }))
  names(TPR) <- paste0("n", seq(from = 20, to = 100, by = 10))
  TPR
}))
names(regions) <- c("chr17_46018653_46019185", "chr19_11199850_11200147", "chr5_10249760_10251253", "chr5_67583971_67584381")
TPR <- regions

save.dir <- "result_files/TPR/"
save(TPR, file = paste0(save.dir, "TPR_", file_name, ".rda"))

