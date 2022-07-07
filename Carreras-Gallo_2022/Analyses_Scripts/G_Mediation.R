#############################################################
## Mediation between Alcohol -> CpG -> High Blood Pressure ##
#############################################################

library(mediation)
library(MultiMed)
library(minfi)

setwd("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Uniform_analysis/Phenotypes/Mediation")
#Load GRset
dataset <- "SVA"
load(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_",dataset,".Rdata"))

#Load top cpgs from Alcohol EWAS
load("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Uniform_analysis/Alcohol/5levels/topcpgs_alcohol_5levels.Rdata")

substance <- "Alcohol_Use_per_week"
phenotype <- "High_Blood_Pressure"

#Select individuals who never drink and who drink regularly
GRset <- GRset[,which(GRset[[substance]]=="Never" | GRset[[substance]]=="Regularly")]
GRset[[substance]] <- droplevels(GRset[[substance]])

#### Only one CpG site as mediator
univariate <- function(mediator){
  data <- cbind(data,data.frame(assay(GRset)[mediator,]))
  colnames(data)[3] <- mediator
  data <- cbind(data,pData(GRset)[,covs])
  
  #Substance -> Phenotype
  design1 <- formula(paste0(phenotype," ~ ",substance,"+",paste0(covs,collapse="+")))
  mod1 <- glm(design1, data=data)

  #Substance -> Methylation
  design.M <- formula(paste0(mediator," ~ ",substance,"+",paste0(covs,collapse="+")))
  model.M <- glm(design.M, data=data)

  #Mediation
  design.Y <- formula(paste0(phenotype," ~ ",substance," + ",mediator,"+",paste0(covs,collapse="+")))
  model.Y <- glm(design.Y, data=data)

  #Test mediation
  res <- mediate(model.M, model.Y, treat=substance, mediator=mediator)

  df <- data.frame(IV_DV_Est = summary(mod1)$coefficients[2,1],
                   IV_DV_p = summary(mod1)$coefficients[2,4],
                   IV_med_Est = summary(model.M)$coefficients[2,1],
                   IV_med_p = summary(model.M)$coefficients[2,4],
                   IV_med_DV_Est = summary(model.Y)$coefficients[2,1],
                   IV_med_DV_p = summary(model.Y)$coefficients[2,4],
                   ACME_Est = res$d0,
                   ACME_p = res$d0.p,
                   ADE_Est = res$z0,
                   ADE_p = res$z0.p,
                   Total_Effect_Est = res$tau.coef,
                   Total_Effect_p = res$tau.p,
                   Prop_Mediated_Est = res$n0,
                   Prop_Mediated_p = res$n0.p)
}

data <- pData(GRset)[,c(substance,phenotype)]
covs <- c("Biological_Sex","age", "Ethnicity",
          "BMI","Level_of_Education","Slide",
          "Bcell","CD4T","CD8T","Eos","Mono","Neu","NK",
          "Tobacco_Use")

list <- lapply(topcpgs$CpG[1:6],univariate)
df <- do.call(rbind,list)
df$CpG <- topcpgs$CpG[1:6]
writexl::write_xlsx(df,"Univariate_mediation.xlsx")


#### Multivariant mediation

E <- as.numeric(GRset[[substance]])
mediators <- topcpgs$CpG[1:200]
M <- t(assay(GRset)[mediators,])
Y <- as.numeric(GRset[[phenotype]])

multi <- data.frame(medTest(E, M, Y, nperm = 500))
rownames(multi) <- 1:nrow(multi)
multi$CpG <- mediators

multi[which(multi$p<0.05),]

writexl::write_xlsx(multi,"Multivariate_mediation.xlsx")
