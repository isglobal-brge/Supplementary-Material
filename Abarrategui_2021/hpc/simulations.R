# Simulations

## 1. Set parameters
n <- Sys.getenv("N")
i <- Sys.getenv("I")

print(paste0("n=", n, ", i=", i)) 

source("hpc/simulation_functions.R")

## 2. TPR

## 2.1 GSE51032

TPR(file_name = "GSE51032", 
    case_sample = "GSM1235784", 
    n = n, 
    methods = c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta"),  
    chr = "chr17", start = 46000000, end = 46020000)

## 2.2. GSE111629

TPR(file_name = "GSE111629", 
    case_sample = c("GSM3035933", "GSM3035791", "GSM3035807","GSM3035685"),
    n = n, 
    methods = c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta"),  
    chr = c("chr5", "chr5", "chr17","chr19"), 
    start = c(10240000, 67580000, 46000000, 11190000), 
    end = c(10260000, 67600000, 46020000, 11210000))

## 3. FPR
### 3.1 GSE51032

FPR(file_name = "GSE51032", 
    epi_found = "names_GSE51032",
    methods = c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta"),
    n = n)

### 3.2. GSE111629

FPR(file_name = "GSE111629", 
    epi_found = "names_GSE111629",
    methods = c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta"),
    n = n)