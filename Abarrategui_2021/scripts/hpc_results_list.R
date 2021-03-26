#1. Create the result list containing the validated epimutations
load.dir<- "hpc/outputs"
folder_name <- "GSE51032"
variable <- "GSE51032"
save.file.name <- "GSE51032_in_chr"


RES.n <- list()
for(n in seq(20, 100, 10)){
  RES.i <- list()
  for(i in 1:100){
    load(paste0(load.dir,"/", folder_name,"/", "outputfile", "-", variable, "-", n, "-", i,".rda"))
      RES.i[[length(RES.i)+1]]<- results 
   
}
  res <- as.data.frame(data.table::rbindlist(RES.i))
  RES.n[[length(RES.n)+1]] <- res
}

names(RES.n) <- paste0("n", seq(20, 100, 10))

result_table <- RES.n

save.dir <- "result_files/hpc_results_list"
save(result_table, file = paste0(save.dir, "/", save.file.name,".rda"))
rm(list = ls())


##########################################################################
# 2. Create the result list containing the not validated epimutations 
load.dir<- "hpc/outputs"
folder_name <- "GSE51032-out"
variable <- "GSE51032-out"
save.file.name <- "GSE51032-out_chr"


RES.n <- list()
for(n in seq(20, 100, 10)){
  RES.i <- list()
  for(i in 1:100){
    load(paste0(folder_name,"/", "outputfile", "-", variable, "-", n , "-", i,".rda"))
    try(if(!is.null(nrow(out))){
      out$i <- i
      RES.i[[length(RES.i)+1]]<- out 
    }, silent = TRUE)
 
  }
  res <- as.data.frame(data.table::rbindlist(RES.i))
  RES.n[[length(RES.n)+1]] <- res
}

names(RES.n) <- paste0("n", seq(20, 100, 10))

result_table <- RES.n

save.dir <- "result_files/hpc_results_list"
save(result_table, file = paste0(save.dir, "/", save.file.name,".rda"))
rm(list = ls())



