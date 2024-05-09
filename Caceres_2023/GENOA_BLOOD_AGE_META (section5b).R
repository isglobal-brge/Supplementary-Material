## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##X-Ra in blood 
#from GENOA study (GSE210255)
#############################

#### libraries
library(XRa)
library(SummarizedExperiment)
library(GEOquery)
library(lmerTest)
library(meta)

#load phenotypes
gsm <- getGEO("GSE210255",destdir ="./Data", AnnotGPL =TRUE)
gsm <- gsm[[1]]

phenobb <- pData(phenoData(gsm))

###Load immunce cell count in genoa study
load("./Data/cellcountsEPIC.RData")

phenobb <- cbind(phenobb, cellcountsEPIC)

###Methylation IDAT data was downloaded from 
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210256
###CpGs from chromosome X selected and saved in metXEPIC.RData
###############################################################

load("./Data/metXEPIC.RData")
met <- metX

met <- matrix(as.numeric(metX[,-1]), nrow=nrow(metX))
rownames(met) <- metX[,1]
colnames(met) <- rownames(phenobb)
cpgidtrue <- rownames(met)

#obtain beta values and phenotypes for females 
metfemale <- t(met)[phenobb$"gender:ch1"=="F", ]
phenofemale <- phenobb[phenobb$"gender:ch1"=="F", ]

#compute XRa and arrange relevant phenotype data
phenofemale$xra <- XRa(metfemale, cpgidtrue)
phenofemale$age <- as.numeric(phenofemale$"age(yrs):ch1") 
phenofemale$plate <- phenofemale$"plate:ch1"
phenofemale$age.dic <-  phenofemale$age > 65
phenofemale$pedigree <- phenofemale$"pedigreeid:ch1"

##Association between XRa and age
#################################
mod <- lmerTest::lmer(xra ~ age+ plate + Bcell + CD4T + CD8T + Mono + Neu + NK + (1|pedigree), 
                      data=phenofemale)

summary(mod)

mod2 <- lmerTest::lmer(xra ~ age.dic+ plate + Bcell + CD4T + CD8T + Mono + Neu + NK + (1|pedigree), 
                      data=phenofemale)

summary(mod2)

save(mod2, file="./Data/modGENOA.RData")

######meta analysis of XRA and age across TRUDIAGNOSTIC, GENOA y MESA (Figure 5A)
######################################################################

load("./Data/mod1Tru.RData")
study1 <- summary(mod1Tru)$coeff["age.dicTRUE", 1:2]

load("./Data/modGENOA.RData")
study2 <- summary(mod2)$coeff["age.dicTRUE", 1:2]

load("./Data/modMonocites.RData")
study3 <- summary(mod3)$coeff["age.dicTRUE", 1:2]

formeta <- rbind(study1, study2, study3)*10^3

matdat <- meta::metagen(formeta[,1], formeta[,2],
                        studlab= c("TruDiagnostic", "GENOA (GSE210255)",  "MESA (GSE56046)"),
                        level.ci = 0.95)


png("./Figures/Figure5A.png", width = 8, height = 8, units = 'in', res = 300) 
meta::forest(matdat, layout="JAMA",  
             leftlabs=c("Study", "beta 95%CI"), xlab="Effect of >65yrs on X-Ra (10e-3)",
             title="", color="grey")
dev.off()

