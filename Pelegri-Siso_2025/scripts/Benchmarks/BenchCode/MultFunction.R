## ##################################################################
## Benchmark Basic Functions
## 
## General Script to test matrix functions:
##      - Multiplication:
##          R implementation
##          BigDataStatMeth - memory version
##          BigDataStatMeth - hdf5 version
##          
##  Number of cores:
##      1, 2, 4, 6, 8, 12, 16, 20, 24
##      
##  
##  Author: Dolors Pelegrí
##  date: 02/2025
## 
## ##################################################################

library(BigDataStatMeth)
library(rhdf5)
library(microbenchmark)
library(ggplot2)
library(ggthemes)


dir.create(wrkdir, recursive = T, showWarnings = FALSE)
setwd(wrkdir)

dir.create("csv_results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("db", showWarnings = FALSE)
dir.create("hdf5_files", showWarnings = FALSE)

source(srcCrateMatrix)

datasets <- createHDF5FilesforBenchmarking(filename = "hdf5_files/MultFunctions.hdf5",
                                           groupname = "InputDatasets", datasetRows = size, datasetCols = size)

save(datasets, file = "db/datasetObject_Mult.RData")

finalBenchmark <- sapply( seq_len(nrow(datasets[["datasets"]])), function(l, ntimes) {
    benchCores <- sapply(ncores, function(nc, ntimes) {
        benchBlockSize <- sapply( blockSizes, function(bs, ntimes) {
            benchChecks <- sapply(check, function(f, ntimes) {
                
                file <- datasets[["file"]]
                group <- datasets[["group"]]
                dataset <- datasets[["datasets"]][l,"dsname"]
                nRows <- datasets[["datasets"]][l,"rows"]
                nCols <- datasets[["datasets"]][l,"cols"]
                
                print(paste0("Testing: ", nRows, " x ", nCols," with blocks of: ", bs, " ==> ", f))
                
                # if(f == "R" && nRows<=25000) {
                if(f == "R" && !bRcrash) {
                    print("R")
                    mat <-  h5read(file, paste(group, dataset, sep = "/"))
                    if(any(class(try(times <- microbenchmark::microbenchmark(MultR = mat %*% mat, times = ntimes, unit = "s"))) %in% "try-error")){
                        bRcrash = TRUE
                    }
                    
                } else {
                    if(f == "BDSM-mem") {
                        mat <-  h5read(file, paste(group, dataset, sep = "/"))
                        
                        times <- microbenchmark::microbenchmark(Mult_mem =  bdblockMult(mat, mat, paral = FALSE), times = ntimes, unit = "s")
                        
                    } else if( f == "BDSM-mem-Block") {
                        mat <-  h5read(file, paste(group, dataset, sep = "/"))
                        
                        times <- microbenchmark::microbenchmark(Mult_memBlock =  bdblockMult(mat, mat, paral = TRUE, threads = nc, block_size = bs, byBlocks = TRUE), times = ntimes, unit = "s")
                        
                    } else if(f == "BDSM-hdf5"){
                        times <- microbenchmark::microbenchmark(MultHDF5 =  
                                                                    bdblockmult_hdf5(filename = file, group = group, 
                                                                                          A = dataset, B = dataset, paral = TRUE,
                                                                                          outgroup = "results", outdataset = "res", 
                                                                                          block_size = bs, threads = nc, force = TRUE ),
                                                                times = ntimes, unit = "s")
                    } 
                }
                
                
                if(is.null(bs)) {
                    bs = "automatic";
                }
                
                if( !exists("times")) {
                    results <- cbind.data.frame( system = system, nCores = nc, funct = f, rows = nRows, cols = nCols, blockSize = bs,
                                                 expr = f, min = NA, lq = NA, mean = NA, median = NA, uq = NA, max = NA, 
                                                 repeats = ntimes, date = Sys.Date(), time = Sys.time())
                } else {
                    results <- cbind.data.frame( system = system, nCores = nc, funct = f, rows = nRows, cols = nCols, blockSize = bs, summary(times)[, c(1:7)], repeats = ntimes, date = Sys.Date(), time = Sys.time())    
                }
                
                if(!file.exists(paste0("csv_results/",st,"_Matirx-",funct,".csv"))) {
                    write.table(results, paste0("csv_results/",st,"_Matirx-",funct,".csv"),
                                append = FALSE,quote = FALSE,sep = "\t", row.names = FALSE)
                } else {
                    write.table(results, paste0("csv_results/",st,"_Matirx-",funct,".csv"),
                                append = TRUE, quote = FALSE,sep = "\t", row.names = FALSE, col.names = FALSE)
                }
                
                rm(mat); gc()
                
                return(results)
                
            }, simplify = FALSE, ntimes = repeats)
            return(do.call(rbind,benchChecks))
        }, simplify = FALSE, ntimes = repeats )
        return(do.call(rbind,benchBlockSize))
    }, simplify = FALSE, ntimes = repeats)
    return(do.call(rbind,benchCores))
}, simplify = FALSE, ntimes = repeats)

finalBenchmark <- do.call(rbind, finalBenchmark)

write.table(finalBenchmark, paste0("csv_results/",system, "_", st,"_Matirx-",funct,"_FINAL.csv"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)

# Load results from file
#..# finalBenchmark <- read.table(paste0("csv_results/",Sys.Date(),"Matirx-Sum_FINAL.csv"), sep = "\t", header = TRUE)
