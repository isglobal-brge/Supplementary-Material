## ##################################################################
## Create files with data to run benchmarks
## 
## Script to create files to be used with benchmarks, every
## file contains a set of matrices to be used to work with hdf5 functions
## 
## Parameters:
##  - filename
##  - datasets names
##  - datasets sizes
##      
##      
##  Author: Dolors Pelegri
##  date: 02/2025
## 
## ##################################################################


createHDF5FilesforBenchmarking <- function( filename, groupname, datasetnames = NULL, datasetRows = NULL, datasetCols = NULL) {
    
    if(is.null(datasetRows))
        datasetRows <- c(100, seq(from = 500, to = 2500, by = 500), seq(from = 5000, to = 20000, by = 5000) )
    
    if(is.null(datasetCols))
        datasetCols <- c(100, seq(from = 500, to = 2500, by = 500), seq(from = 5000, to = 20000, by = 5000) )
    
    if(length(datasetRows) != length(datasetCols))
        stop("Row and Col have different lengths")
    
    dims <- cbind.data.frame (rows = datasetRows, cols = datasetCols)
    
    if(is.null(datasetnames))
        datasetnames <- apply( dims, 1, function(d){ paste0("mat_",d["rows"],"x",d["cols"]) } )
    
    dims <- cbind.data.frame (rows = datasetRows, cols = datasetCols, dsname = datasetnames)
    
    if(file.exists(filename))
        file.remove( filename )
        
    sapply(seq_len(nrow(dims)), function(d, dsdims) {
        print(dsdims[d, "dsname"])
        set.seed(555)
        mat <- matrix( rnorm( dsdims[d, "rows"]*dsdims[d, "cols"], mean=0, sd=10), dsdims[d, "rows"], dsdims[d, "cols"]) 
        
        bdCreate_hdf5_matrix(filename = filename, 
                             object = mat, group = groupname, 
                             dataset = dsdims[d, "dsname"],
                             transp = FALSE,
                             overwriteFile = FALSE, 
                             overwriteDataset = TRUE, 
                             unlimited = FALSE)    
    }, dsdims = dims)
    
    message("File and dataset created")
    
    return(list(file = filename, 
                group = groupname,
                datasets = dims))
}
