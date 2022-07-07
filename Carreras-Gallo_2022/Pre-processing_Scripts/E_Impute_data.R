#################
## IMPUTE DATA ##
#################

start <- Sys.time()

#Set "clean" or "whole" dataset
dataset <- "clean"

#Load libraries
library(minfi)
library(tidyverse)

#Load beta values
load(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_",dataset,".Rdata"))

#Impute data by median (it lasts ~ 7 minuts)
assay(GRset) <- assay(GRset) %>% 
  as.data.frame() %>%
  mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = T), x)) %>%
  as.matrix()

gc(reset=TRUE)

#Impute data using a regression-based imputation (it lasts for days)
#assay(GRset) <- t(methyLImp::methyLImp(t(assay(GRset)),min=0,max=0))

#Save new GRset with imputed data
save(GRset,file=paste0("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_",dataset,"_imp.Rdata"))

end <- Sys.time()
total <- end-start
total
