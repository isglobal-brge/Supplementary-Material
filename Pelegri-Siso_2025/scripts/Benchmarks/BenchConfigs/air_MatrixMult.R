## ##################################################################
## Benchmark Basic Functions
## 
## General Script to test matrix functions:
## 
##      - Sum:
##          R implementation
##          BigDataStatMeth - memory version
##          BigDataStatMeth - hdf5 version
##          
##      - Substract
##          R implementation
##          BigDataStatMeth - memory version
##          BigDataStatMeth - hdf5 version
##          
##      - Multiplication:
##          R implementation
##          BigDataStatMeth - memory version
##          BigDataStatMeth - hdf5 version
##          
##  Number of cores:
##      1, 2, 4, 6
##      
##  
##  Author: Dolors Pelegrí
##  e-mail: dolors.pelegri@isglobal.org
##  date: 10/2025
## 
## ##################################################################


setwd("~/PhD/TREBALLANT/BigDataStatMeth_Benhmarks")
source("code/BenchCode/createMatrixFile.R")

wrkdir <- "~/PhD/TREBALLANT/BigDataStatMeth_Benhmarks"
srcCrateMatrix <- "code/BenchCode/createMatrixFile.R"
srcBenchmarkSum <-  "code/BenchCode/SumFunction.R"
srcBenchmarkSubs <-  "code/BenchCode/SubstractFunction.R"
srcBenchmarkMult <-  "code/BenchCode/MultFunction.R"

check <- c( "R", "BDSM-mem", "BDSM-mem-Block", "BDSM-hdf5")
system <- "macOS"
size <- c(100, 1000, 2500, 5000, seq(from = 10000, to = 30000, by = 10000) ) # 35000 maàxim a R amb RAM actual macbook (no suporta mida superior...)
blockSizes <- c(500, 1000, 2500, 5000, 10000, NULL)
repeats <- 3
ncores <- c(1, 2, 4, 6)
bRcrash <-  FALSE
st <- format(Sys.time(), "%Y-%m-%d__%H-%M")

funct <- c( "Multiplication")
source(srcBenchmarkMult)
