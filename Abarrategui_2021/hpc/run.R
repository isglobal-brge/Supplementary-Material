
i <- Sys.getenv("I")

run <- function(file_name){
  
  load(file = paste0("inputs/",file_name,".rda"))
  load(paste0("inputs/names_",file_name,".rda"))
  case_samples <- object[, names]
  control_samples <- object[,!colnames(object) %in% names]
  
  time <- system.time(res <- epimutacions::epimutations(case_samples, control_samples, method = i))
  res$time <- time[[2]]
    
  output.filename <- paste0("outputs/", file_name, "-",  i, ".rda")
  save(res, file = output.filename)
  print(paste0("Written to file '", output.filename, "'"))
}

#a <- run("GSE51032")
b <- run("GSE111629")

