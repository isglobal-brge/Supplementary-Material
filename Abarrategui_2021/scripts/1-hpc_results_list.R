#1. Create the result list containing the validated epimutations

create_result_list<- function(load.dir, folder_name,rate,  variable, save_name, out_name){
    RES.n <- list()
    for(n in seq(20, 100, 10)){
          RES.i <- list()
          for(i in 1:100){
            load(paste0(load.dir,"/", folder_name,"/outputfile-", rate, "-", variable, "-", n, "-", i,".rda"))
            RES.i[[length(RES.i)+1]] <- eval(parse(text=out_name))}
          res <- as.data.frame(data.table::rbindlist(RES.i))
          RES.n[[length(RES.n)+1]] <- res
}

names(RES.n) <- paste0("n", seq(20, 100, 10))

result_table <- RES.n

save.dir <- "result_files/1-Result_list"
save(result_table, file = paste0(save.dir, "/", save_name,".rda"))
}

#2. Creating results
baseDir <- "hpc/outputs"
##2.1. TPR
create_result_list(load.dir = baseDir,
                   folder_name = "GSE111629-TPR", rate = "TPR", variable = "GSE111629", save_name = "GSE111629-TPR", out_name = "results")
create_result_list(load.dir = baseDir,
                   folder_name = "GSE51032-TPR", rate = "TPR", variable = "GSE51032", save_name = "GSE51032-TPR", out_name = "results")
##2.2. FPR
create_result_list(load.dir = baseDir,
                   folder_name = "GSE111629-FPR", rate = "FPR", variable = "GSE111629", save_name = "GSE111629-FPR", out_name = "out")
create_result_list(load.dir = baseDir,
                   folder_name = "GSE51032-FPR", rate = "FPR", variable = "GSE51032", save_name = "GSE51032-FPR", out_name = "out")
