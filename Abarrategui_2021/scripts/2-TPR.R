#1. Function to calculate the TPR and accuracy
TPR <- function(file_name, epi_validated){
#1.1. Load the data file containing the results list  
  load(paste0("result_files/1-Result_list/", file_name, ".rda"))
#1.2. Set the parameters
  control_n <- seq(20, 100, 10)
  methods <- c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "beta")
#1.3. Calculate TPR and accuracy
regions <-  do.call(list, lapply(seq(nrow(epi_validated)), function(i){
  TPR <- do.call(list, lapply(seq(length(result_table)), function(n){
    rst <-  do.call(rbind, lapply(seq_len(length(methods)), function(j){
      df_n <- result_table[[n]]
      keep_region <- which(df_n$chromosome_validated == epi_validated$seqnames[i] & df_n$start_validated == epi_validated$start[i])
      region <- df_n[keep_region,]
      keep_method <- grepl(methods[j],unlist(region[,"epi_id"]))
      out <- region[keep_method,]
      if(methods[j] == "manova" | methods[j] == "mlm"){
        out <- out[out$pvalue < 0.05/40408,]
      }else if(methods[j] == "isoforest"){
        out <- out[out$outlier_score > 0.7,]
      }
      if(i == 4){
        TPR <- nrow(out)/200
      }else{
      TPR <- nrow(out)/100
      }
      percent <- mean(out[,"%"])
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
names(regions) <- paste0(epi_validated$seqnames, "_", epi_validated$start, "_", epi_validated$end)
TPR <- regions
#1.4. Save results
save(TPR, file = paste0("result_files/2-TPR_accuracy/", file_name, ".rda"))
}

# 2. Create TPR and accuracy tables 
## 2.1. load the validated epimutations
load("hpc/inputs/epi_validated.rda")
epi_validated <- as.data.frame(epi_validated)

## 2.2. Table for 'GSE111629' sample
GSE111629_regions <- epi_validated
TPR("GSE111629-TPR", GSE111629_regions)
## 2.3. Table for 'GSE51032' sample
GSE51032_region <- epi_validated[1,]
TPR("GSE51032-TPR", GSE51032_region)




