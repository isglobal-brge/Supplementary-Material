#################################
## SURROGATE VARIABLE ANALYSIS ##
#################################

library(SmartSVA)
library(minfi)
library(meffil)

#Load GenomicRatioSet
load("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_clean_imp.Rdata")

#Select only europeans
# GRset <- GRset[,which(GRset$Ethnicity=="Euro")]

#### ESTIMATE NUMBER OF SVs
# #Select the 50K CpGs most variables to estimate SVs
# var.sites <- meffil.most.variable.cpgs(assay(GRset), n=50000, autosomal=T)
# 
# #Obtain only the beta values with features at rows and samples at columns
# mat <- assay(GRset)[var.sites,]
# 
# #create a data.frame with the phenotype data
# pheno <- pData(GRset)
# 
# rm(GRset,var.sites)
# gc(reset = TRUE)
# 
# #Calculate the residual values
# mat.res <- t(resid(lm(t(mat) ~ as.factor(Marijuana)+Biological_Sex+Ethnicity+age_estimated+
#                         Neuropsychological_any+Cardiovascular_any+Respiratory_Disease_any+
#                         Endocrine_Disease_any+Tobacco_Use+Amphetamines+Benzodiazepines+
#                         Hallucinogens+MDMA+Drug_Alcohol_mother+
#                         Drug_Alcohol_father+Alcohol_Use_per_week,
#                    data=pheno)))
# 
# gc(reset=TRUE)
# 
# #Determine the number of SVs to be calculated
# n.sv <- isva::EstDimRMT(mat.res, FALSE)$dim + 1
# n.sv

#### RUN SVA
#Create the model
#Without Ethnicity when only Europeans selected
# mod <- model.matrix( ~ as.factor(Marijuana)+Biological_Sex+age+
#                        Neuropsychological_any+Cardiovascular_any+Respiratory_Disease_any+
#                        Endocrine_Disease_any+Tobacco_Use+Amphetamines+Benzodiazepines+
#                        Hallucinogens+MDMA+Drug_Alcohol_mother+
#                        Drug_Alcohol_father+Alcohol_Use_per_week, pData(GRset))

#With Ethnicity when the whole dataset selected
mod <- model.matrix( ~ as.factor(Marijuana)+Biological_Sex+Ethnicity+age+
                       Neuropsychological_any+Cardiovascular_any+Respiratory_Disease_any+
                       Endocrine_Disease_any+Tobacco_Use+Amphetamines+Benzodiazepines+
                       Hallucinogens+MDMA+Drug_Alcohol_mother+
                       Drug_Alcohol_father+Alcohol_Use_per_week, pData(GRset))

mod0 <- model.matrix(~1, data=pData(GRset))

#Run SVA
sv.obj <- smartsva.cpp(assay(GRset), mod, mod0, n.sv=60)

#Save sva object
save(sv.obj,file="/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/Batch/sv.obj.Rdata")

#load("/PROJECTES/GENOMICS/TruDiagnostic/prepro_files/Batch/sv.obj.Rdata")

SVs <- sv.obj$sv
colnames(SVs) <- paste0("SV",1:ncol(SVs))
mat <- t(getBeta(GRset))

#Estimate residuals for SVs
mat_res <- apply(mat, 2, function(x) {
  ex <- data.frame(cpg=x)
  ex <- cbind(ex,SVs)
  residuals(lm(formula(paste0("cpg ~ ",paste(colnames(SVs),collapse="+"))),
               data=ex,na.action=na.exclude))
}
)

gc(reset=TRUE)

GRset@assays@data$Beta <- as.matrix(t(mat_res))

#Save the new GRset
save(GRset, file="/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_SVA.Rdata")