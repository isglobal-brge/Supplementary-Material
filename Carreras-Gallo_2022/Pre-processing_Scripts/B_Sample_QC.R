###############
## Sample QC ##
###############

# Load package and set number of cores to 2

library(meffil)
library(readr)
library(stringr)
options(mc.cores=8)

#Set working directory
setwd("/PROJECTES/GENOMICS/TruDiagnostic/idat_files/idat_use")

#Load SampleSheet
load("SampleSheet.Rdata")

#Generate QC objects
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", 
                        verbose=TRUE,
                        detection.threshold	= 0.01,
                        bead.threshold = 3)

length(qc.objects)

setwd("/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/QC")

save(qc.objects,file="qc_objects_whole.Robj")


#Estimate cellular composition from the whole dataset
cc<-t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc<-data.frame(PatientID=row.names(cc),cc)

write.table(cc,paste0("cellcounts_whole.txt"),sep="\t",row.names=F,col.names=T,quote=F)


#Generating a QC summary and a QC report

qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.05,
  detectionp.samples.threshold          = 0.05,
  detectionp.cpgs.threshold             = 0.05, 
  beadnum.cpgs.threshold                = 0.05,
  sex.outlier.sd                        = 3
)

#QC summary
qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters,
  genotypes=NULL,
  verbose=TRUE
)

save(qc.summary, file="qcsummary_whole.Robj")

#QC report
meffil.qc.report(qc.summary, output.file="qc-report_whole.html")

#Removing bad samples
outlier <- qc.summary$bad.samples
table(outlier$issue)

#Select all the ouliers except the outliers based on control probes as there are 
#often a few probes for each control type only
index <- outlier$issue %in% c("Control probe (dye.bias)", 
                              "Methylated vs Unmethylated",
                              "X-Y ratio outlier",
                              "Low bead numbers",
                              "Detection p-value",
                              "Sex mismatch",
                              "Control probe (bisulfite1)",
                              "Control probe (bisulfite2)")

outlier <- outlier[index,]

#Remove bad samples
length(qc.objects)
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects)
save(qc.objects,file="qc_objects_clean.Robj")

#Rerun QC summary on clean dataset
qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters=qc.parameters,
  genotypes=NULL,
  verbose=TRUE
)

save(qc.summary, file="qcsummary_clean.Robj")

#Count number of bad cpgs
table(qc.summary$bad.cpgs$issue)

#Rerun QC report on clean dataset
meffil.qc.report(qc.summary, output.file="qc-report_clean.html")

#Estimate cellular composition from the clean dataset
cc<-t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc<-data.frame(PatientID=row.names(cc),cc)

write.table(cc,paste0("cellcounts_clean.txt"),sep="\t",row.names=F,col.names=T,quote=F)