#################
## SampleSheet ##
#################

# Load package and set number of cores to 4

library(meffil)
library(readr)
library(stringr)
options(mc.cores=4)

#Set working directory

setwd("/PROJECTES/GENOMICS/TruDiagnostic/idat_files/idat_use")

#Create Sample Sheet

samplesheet <- meffil.create.samplesheet("./", recursive=TRUE)
dim(samplesheet)

#Load metadata
load("../../metadata/metadata.Rdata")

#Include sex and PID from metadata in the SampleSheet
samplesheet$PID <- NA
add_sex_PID <- function(row){
  id <- samplesheet$Sample_Name[row]
  
  #Prove that the id is included in the MetaData dataset
  if(id%in%metadata$PatientID){
    sex <- metadata[which(metadata$PatientID==id),]$Biological_Sex
    
    #Select only the first letter (M/F) instead of Male/Female -> needed for the QC
    samplesheet$Sex[row] <<- substr(sex,1,1)
    samplesheet$PID[row] <<- metadata[which(metadata$PatientID==id),]$PID
  }
}

sapply(1:nrow(samplesheet),add_sex_PID)

#Select only the idat files that have metadata
samplesheet <- samplesheet[!is.na(samplesheet$Sex),]

# Remove troubling samples
samplesheet <- samplesheet[!(samplesheet$Sample_Name=="205772280052_R06C01" | samplesheet$Sample_Name=="205772280052_R07C01" | samplesheet$Sample_Name=="205772280052_R08C01"),]

#For the patients who took repeated measures, select only one row
samplesheet <- samplesheet[!duplicated(samplesheet$PID),]

#Save Sample Sheet
save(samplesheet, file="SampleSheet.Rdata")

