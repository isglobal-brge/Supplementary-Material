#0. load necessary data and functions
source("scripts/simulations_results_functions.R")
load("hpc/inputs/epi_validated.rda")

#1. Result list 
##1.1. TPR
create_result_list(folder_name = "GSE111629-TPR", rate = "TPR", variable = "GSE111629", save_name = "GSE111629-TPR", epi_validated = epi_validated)
create_result_list(folder_name = "GSE51032-TPR", rate = "TPR", variable = "GSE51032", save_name = "GSE51032-TPR", epi_validated = epi_validated)

create_result_list(folder_name = "ramr-GSE111629-TPR", 
                   rate = "TPR", 
                   variable = "GSE111629",
                   case_sample = c("GSM3035933", "GSM3035791", "GSM3035807","GSM3035685"),
                   save_name = "ramr-GSE111629-TPR", epi_validated = epi_validated)
create_result_list(folder_name = "ramr-GSE51032-TPR", 
                   rate = "TPR", 
                   variable = "GSE51032", 
                   case_sample = "GSM1235784",
                   save_name = "ramr-GSE51032-TPR", epi_validated = epi_validated)

##1.2. FPR
create_result_list(folder_name = "GSE111629-FPR", rate = "FPR", variable = "GSE111629", save_name = "GSE111629-FPR", epi_validated = NULL)
create_result_list(folder_name = "GSE51032-FPR", rate = "FPR", variable = "GSE51032", save_name = "GSE51032-FPR", epi_validated = NULL)

create_result_list(folder_name = "ramr-GSE111629-FPR", rate = "FPR", variable = "GSE111629", save_name = "ramr-GSE111629-FPR", epi_validated = NULL)
create_result_list(folder_name = "ramr-GSE51032-FPR", rate = "FPR", variable = "GSE51032", save_name = "ramr-GSE51032-FPR", epi_validated = NULL)

# 2. TPR, FPR and Accuracy

## 2.1. Table for 'GSE111629' sample
GSE111629_regions <- epi_validated
TPR("GSE111629-TPR", as.data.frame(GSE111629_regions))
FPR("GSE111629-FPR")
TPR_FPR_accuracy("GSE111629-FPR", "GSE111629-TPR", "TPR_FPR_accuracy_GSE111629")

TPR("ramr-GSE111629-TPR", as.data.frame(GSE111629_regions))
FPR("ramr-GSE111629-FPR")
TPR_FPR_accuracy("ramr-GSE111629-FPR", "ramr-GSE111629-TPR", "ramr-TPR_FPR_accuracy_GSE111629")

## 2.2. Table for 'GSE51032' sample
GSE51032_region <- epi_validated[1,]
TPR("GSE51032-TPR", as.data.frame(GSE51032_region))
FPR("GSE51032-FPR")
TPR_FPR_accuracy("GSE51032-FPR", "GSE51032-TPR", "TPR_FPR_accuracy_GSE51032")

TPR("ramr-GSE51032-TPR", as.data.frame(GSE51032_region))
FPR("ramr-GSE51032-FPR")
TPR_FPR_accuracy("ramr-GSE51032-FPR", "ramr-GSE51032-TPR", "ramr-TPR_FPR_accuracy_GSE51032")

#3. Create the graphics

graph("ramr-TPR_FPR_accuracy_GSE111629")
graph("ramr-TPR_FPR_accuracy_GSE51032")
