## ##################################################################
## Benchmark Basic Functions
## 
## General Script to test matrix functions:
## 
##      - SVD:
##          R implementation
##          BigDataStatMeth - hdf5 version
##          BigDataStatMeth - hdf5 by Blocks version
##          
##  Number of cores:
##      1, 2, 4, 6
##      
##  k and q ( when working by blocks)
##      - k = 1 ; q = 4
##      - k = 2 ; q = 2
##  
##  Author: Dolors Pelegrí
##  date: 10/2025
## 
## ##################################################################
# devtools::reload(pkgload::inst("BigDataStatMeth"))

setwd("~/PhD/TREBALLANT/BigDataStatMeth_Benhmarks")
# source("code/BenchCode/createMatrixFile.R")

wrkdir <- "~/PhD/TREBALLANT/BigDataStatMeth_Benhmarks"
srcCrateMatrix <- "code/BenchCode/createMatrixFile.R"
srcBenchmarkSVD <-  "code/BenchCode/SVDFunction.R"

check <- c( "R", "BDSM-hdf5-full", "BDSM-hdf5-auto")
system <- "macOS"


repeats <- 2
ncores <- c(2, 4, 6)
bRcrash <-  FALSE

# SVD with Horizontal Matrices 
# ------------------------------

row_sizes <- c(rep(c(100, 500, 1000, 2500), each = 4)) # 35000 maàxim pot R amb RAM actual macbook (no suporta mida superior...)
col_sizes <- c( rep(seq(from = 10000, to = 100000, by = 30000), 4) ) # 35000 maàxim pot R amb RAM actual macbook (no suporta mida superior...)

st <- format(Sys.time(), "%Y-%m-%d__%H-%M")

funct <- c( "SVD_H")
source(srcBenchmarkSVD)

# SVD with Vertical Matrices 
# ------------------------------

row_sizes <- c( rep(seq(from = 10000, to = 100000, by = 30000), 4) ) # 35000 maàxim pot R amb RAM actual macbook (no suporta mida superior...)
col_sizes <- c(rep(c(100, 500, 1000, 2500), each = 4)) # 35000 maàxim pot R amb RAM actual macbook (no suporta mida superior...)

st <- format(Sys.time(), "%Y-%m-%d__%H-%M")

funct <- c( "SVD_V")
source(srcBenchmarkSVD)
