##############################
## Functional normalization ##
##############################

#Set "clean" or "whole" dataset
dataset <- "clean"

#Load libraries
library(meffil)
library(naniar)
library(ggplot2)

#Set working directory
setwd("/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/QC")

#Load in the QC objects
load(paste0("qc_objects_",dataset,".Robj"))
length(qc.objects)
load(paste0("qcsummary_",dataset,".Robj"))

#Estimate the number of principal components to use
# y <- meffil.plot.pc.fit(qc.objects)
# ggsave(y$plot,filename=paste0("../Functional_Normalization/pc_fit_",dataset,".png"),height=6,width=6)

#Set the number of PCs to use going forward
pc <- 10

#Perform functional normalisation and remove bad CpGs and poor signal values
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects,file=paste0("../Functional_Normalization/norm_obj_pc",pc,"_",dataset,".Robj"))

setwd("/PROJECTES/GENOMICS/TruDiagnostic/idat_files/idat_use")

gds.filename <- paste0("../../prepro_files/Functional_Normalization/norm_beta_",dataset,".gds")

meffil.normalize.samples(
  norm.objects,
  just.beta=T,
  remove.poor.signal = T,
  cpglist.remove=qc.summary$bad.cpgs$name,
  gds.filename=gds.filename,
  verbose=T,
  mc.cores=1)

norm.beta <- meffil.gds.methylation(gds.filename)

#Garbage collector
gc(reset=TRUE)

##############
## MISSING ##
##############

#Percentage of missing
na_tot <- sum(is.na(norm.beta))
na_per <- na_tot/(nrow(norm.beta)*as.numeric(ncol(norm.beta)))*100
na_per

gc(reset=TRUE)

#Missing per CpG
miss_cpgs <- miss_var_summary(data.frame(t(norm.beta)), order = TRUE)

gc(reset=TRUE)

png(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Missing_values/Hist_CpGs_methylome_",dataset,".png"),
    res = 1200,
    width = 10, height = 3.25, units = "in", pointsize = 3)
par(mar=c(5,5,2,2))

ggplot(data=miss_cpgs, aes(x=1:nrow(miss_cpgs), y=pct_miss)) +
  geom_line(size=1, colour="red") +
  labs(title="Percentage of missing per CpG", 
       x="CpG sites", y = "Percentage of missing (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())

dev.off()

#Missing per ID
miss_IDs <- miss_case_summary(data.frame(t(norm.beta)), order = TRUE)
miss_IDs$case <- colnames(norm.beta)

gc(reset=TRUE)

png(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Missing_values/Hist_IDs_methylome_",dataset,".png"),
    res = 1200,
    width = 5, height = 3.25, units = "in", pointsize = 3)
par(mar=c(5,5,2,2))

ggplot(data=miss_IDs, aes(x=case, y=pct_miss)) +
  geom_bar(stat="identity") +
  labs(title="Percentage of missing per ID", 
       x="IDs", y = "Percentage of missing (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())

dev.off()

#Remove the probes and samples with more than 5% of NAs
bad_cpgs <- miss_cpgs[which(miss_cpgs$pct_miss>5),]$variable
length(bad_cpgs)
bad_ids <- miss_IDs[which(miss_IDs$pct_miss>5),]$case
bad_ids

norm.beta <- norm.beta[!rownames(norm.beta) %in% bad_cpgs, !colnames(norm.beta) %in% bad_ids]
norm.objects <- norm.objects[!names(norm.objects)%in% bad_ids]

gc(reset=TRUE)

#Percentage of missing after removal
na_tot <- sum(is.na(norm.beta))
na_per <- na_tot/(nrow(norm.beta)*as.numeric(ncol(norm.beta)))*100
na_per

gc(reset=TRUE)

#Save normalized betas
setwd("/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/Functional_Normalization")
save(norm.beta,file=paste0("norm_beta_pc",pc,"_",dataset,".Robj"))

#Generate normalization report
#Code batch variables as factors
for (i in 1:length(norm.objects)){
  norm.objects[[i]]$samplesheet$Slide<-as.factor(norm.objects[[i]]$samplesheet$Slide)
  norm.objects[[i]]$samplesheet$Sex<-as.factor(norm.objects[[i]]$samplesheet$Sex)
  norm.objects[[i]]$samplesheet$sentrix_row<-as.factor(norm.objects[[i]]$samplesheet$sentrix_row)
  norm.objects[[i]]$samplesheet$sentrix_col<-as.factor(norm.objects[[i]]$samplesheet$sentrix_col)
}

batch_var<-c("Slide","Sex","sentrix_row","sentrix_col")

norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables=batch_var,
  control.pcs=1:5,
  batch.pcs=1:5,
  batch.threshold=0.01
)

#Calculate PCs of normalized beta based on the 50000 most variable CpGs (by default) 
pcs <- meffil.methylation.pcs(norm.beta,verbose=T, full.obj = T)
save(pcs,file=paste0("pcs_norm_beta_",dataset,".Robj"))
norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs$x, parameters=norm.parameters)
meffil.normalization.report(norm.summary, output.file=paste0("normalization-report_",dataset,".html"))
