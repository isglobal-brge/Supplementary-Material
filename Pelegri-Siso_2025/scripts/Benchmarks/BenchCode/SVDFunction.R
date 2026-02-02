## ##################################################################
## Benchmark Basic Functions
## 
## General Script to test matrix functions:
##      - SVD decomposition:
##          R implementation
##          BigDataStatMeth - hdf5 version
##          
##  Number of cores:
##      1, 2, 4, 6
##      
##  
##  Author: Dolors Pelegrí
##  date: 10/2025
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


# filename = "hdf5_files/SVDFunctions.hdf5"
# groupname = "InputDatasets" 
# datasetRows = size
# datasetCols = size

datasets <- createHDF5FilesforBenchmarking(filename = "hdf5_files/SVDFunctions.hdf5",
                                           groupname = "InputDatasets", datasetRows = row_sizes, datasetCols = col_sizes)

save(datasets, file = "db/datasetObject_SVD.RData")

finalBenchmark <- sapply( seq_len(nrow(datasets[["datasets"]])), function(l, ntimes) {
    benchCores <- sapply(ncores, function(nc, ntimes) {
            benchChecks <- sapply(check, function(f, ntimes) {
                
                file <- datasets[["file"]]
                group <- datasets[["group"]]
                dataset <- datasets[["datasets"]][l,"dsname"]
                nRows <- datasets[["datasets"]][l,"rows"]
                nCols <- datasets[["datasets"]][l,"cols"]
                
                print(paste0("Testing: ", nRows, " x ", nCols, " ==> ", f))
                
                if(f == "R" && !bRcrash && nc ==1) {
                    print("R")
                    mat <-  h5read(file, paste(group, dataset, sep = "/"))
                    if(any(class(try(times <- microbenchmark::microbenchmark(SVD_R = svd(mat) , times = ntimes, unit = "s"))) %in% "try-error")){
                        bRcrash = TRUE
                    }
                    rm(mat); gc()
                    
                } else {
                    
                    if(f == "BDSM-hdf5-full" && nc ==1) {
                        times <- microbenchmark::microbenchmark(SVDHDF5 =  
                                                                    bdSVD_hdf5(file, group = group, 
                                                                               dataset = dataset, method = "full",
                                                                               bcenter = TRUE, bscale = TRUE, 
                                                                               rankthreshold = 0.0, 
                                                                               overwrite  = TRUE, threads = nc ),
                                                                times = ntimes, unit = "s")
                        
                    } else if(f == "BDSM-hdf5-auto" && (nRows<10000 || (nRows >= 10000 && nc>1) )  ){
                        times <- microbenchmark::microbenchmark(SVDHDF5B_q1_k4 =  
                                                                    bdSVD_hdf5(file, group = group, 
                                                                               dataset = dataset, method = "auto",
                                                                               q = 1, k = 4,  bcenter = TRUE, bscale = TRUE, 
                                                                               rankthreshold = 0.0, 
                                                                               overwrite  = TRUE, threads = nc ),
                                                                SVDHDF5B_q1_k6 =  
                                                                    bdSVD_hdf5(file, group = group, 
                                                                               dataset = dataset, method = "auto",
                                                                               q = 1, k = 6,  bcenter = TRUE, bscale = TRUE, 
                                                                               rankthreshold = 0.0, 
                                                                               overwrite  = TRUE, threads = nc ),
                                                                SVDHDF5B_q2_k4 =  
                                                                    bdSVD_hdf5(file, group = group, 
                                                                               dataset = dataset, method = "auto",
                                                                               q = 2, k = 4,  bcenter = TRUE, bscale = TRUE, 
                                                                               rankthreshold = 0.0, 
                                                                               overwrite  = TRUE, threads = nc ),
                                                                SVDHDF5B_q2_k6 =  
                                                                    bdSVD_hdf5(file, group = group, 
                                                                               dataset = dataset, method = "auto",
                                                                               q = 2, k = 6,  bcenter = TRUE, bscale = TRUE, 
                                                                               rankthreshold = 0.0, 
                                                                               overwrite  = TRUE, threads = nc ),
                                                                SVDHDF5B_q3_k4 =  
                                                                    bdSVD_hdf5(file, group = group, 
                                                                               dataset = dataset, method = "auto",
                                                                               q = 3, k = 4,  bcenter = TRUE, bscale = TRUE, 
                                                                               rankthreshold = 0.0, 
                                                                               overwrite  = TRUE, threads = nc ),
                                                                times = ntimes, unit = "s")
                    } 
                    
                    
                }
                
                bs <- "automatic";
                
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
                
                return(results)
                
            }, simplify = FALSE, ntimes = repeats)
            return(do.call(rbind,benchChecks))
    }, simplify = FALSE, ntimes = repeats)
    
    return(do.call(rbind,benchCores))
}, simplify = FALSE, ntimes = repeats)

finalBenchmark <- do.call(rbind, finalBenchmark)

write.table(finalBenchmark, paste0("csv_results/",system, "_", st,"_Matirx-",funct,"_FINAL.csv"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)

# Load results from file
#..# finalBenchmark <- read.table(paste0("csv_results/",Sys.Date(),"Matirx-Sum_FINAL.csv"), sep = "\t", header = TRUE)

