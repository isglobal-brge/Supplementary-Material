load("result_files/TPR/TPR_GSE51032_in_chr.rda")
load("result_files/FPR/FPR_GSE51032_out_chr.rda")
final.dataset.name <- "GSE51032_chr"
TPR <- as.data.frame(data.table::rbindlist(TPR))
FPR <- df <- as.data.frame(data.table::rbindlist(FPR))
df<- merge(TPR, FPR, by = c("method", "n"), sort = FALSE)
save.dir <-"result_files/TPR_FPR_accuracy/" 
save(df, file = paste0(save.dir, final.dataset.name, ".rda"))

####################################################
load("result_files/TPR/TPR_GSE111629_in_chr_start_end.rda")
load("result_files/FPR/FPR_GSE111629_out_chr_start_end.rda")
final.dataset.name <- "GSE111629_chr_start_end"

df <- do.call(list, lapply(seq_len(length(TPR)), function(i) {
  TPR_n <- as.data.frame(data.table::rbindlist(TPR[[i]]))
  FPR <- as.data.frame(data.table::rbindlist(FPR))
  merge(TPR_n, FPR, by = c("method", "n"), sort = FALSE)
}))

names(df) <- names(TPR)

save.dir <-"result_files/TPR_FPR_accuracy/" 
save(df, file = paste0(save.dir, final.dataset.name, ".rda"))
rm(list = ls())
