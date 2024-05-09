## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##X-Ra in blood 
#from TruDIagnostic biobank
########################### 


#### libraries
library(chrXRa)
library(SummarizedExperiment)
library(arm)
library(ggplot2)
library(smplot2)
library(ggpubr)
library(methylclock)
library(GEOquery)
library(visreg)
library(lmerTest)


#### TRUDIAGNOSTIC

### methylation data for X chromosome 
###(Data request varun@trudiagnostic.com)
#########################################

load("./Data/GRsetchrX.RData")

met <- assay(GRsetchrX)
cpgidtrue <- rownames(met)

phenobb <- pData(GRsetchrX)
subidspheno <- phenobb$PatientID
rownames(phenobb) <-  subidspheno

phenobb$sex <-  phenobb$Biological_Sex

metfemale <- t(met)[phenobb$sex=="Female", ]
phenofemale <- phenobb[phenobb$sex=="Female", ]


###########################################
#compute XRa and obtain relevant phenotypes
###########################################

phenofemale$xra <- XRa(metfemale, cpgidtrue)
phenofemale$menopause <- as.numeric(phenofemale$Menopause=="Yes")
phenofemale$age.dic <-  phenofemale$age > 65

#cancer status
phenofemale$cannum <- as.numeric(factor(phenofemale$Cancer_Diagnosis_any))-1


#########################
####Associations with age
#########################

mod <- lm(xra ~ age + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, 
           data = phenofemale, subset=which(phenofemale$cannum==0))

summary(mod)


####Associations with age>65 save result to validate in meta-analysis 
#performed in GENOA_BLOOD_AGE_META (Figure 2A)
##############################################

mod1Tru <- lm(xra ~ age.dic  
           + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, 
           data = phenofemale, subset=which(phenofemale$cannum==0))

summary(mod1Tru)

save(mod1Tru, file="./Data/mod1Tru.RData")


####XRa in age intervals (Figure 5B1)
#####################################

phenofemale$ageint <- cut(phenofemale$age, breaks = c(0,seq(60, 80, 10), 100))


dat <- data.frame(XRa=phenofemale$xra, ageint=phenofemale$ageint)
dat <- dat[]


plot1 <- ggplot(data = dat, mapping = aes(x = ageint, y = XRa, fill = ageint)) +
  sm_bar(errorbar_type = 'ci', point.params = list(size = 5, shape = 21,
                                                               color = 'transparent', alpha = 0.6,
                                                               position = sdamr::position_jitternudge(nudge.x = 0, 
                                                                                                      jitter.width = 0.5))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(9, 'Blues')[c(4,6,7,8,9)])+
  theme(axis.title.x = element_blank()) +
  sm_minimal()+
  theme(
        axis.text = element_text(size = 35),
        axis.title = element_text(size = 35),
        plot.title = element_text(size = 36),
        legend.text = element_text(size = 35))+
  labs(title = "TruDiagnostic", x = "Age", y = "X-Ra")


ggsave("./Figures/Figure5B1.png", plot1, width = 12, height = 10, dpi = 300)


##interaction with menopause
############################
mod <- glm(xra  ~ age*Menopause 
            + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, 
            data=phenofemale,subset = which(phenofemale$cannum==0))

summary(mod)

##stratification by menopause
mod <- glm(xra  ~ age 
           + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, 
           data=phenofemale, 
           subset = which(phenofemale$cannum==0 & phenofemale$Menopause=="Yes"))

summary(mod)

###################################
## Association with Telomere length
## Figure 5C
###################################

sel <-  rowSums(is.na(phenofemale[,c("age", "xra", "TelomereValues", "Bcell","CD4T", "CD8T", "Eos","Mono", "Neu" , "NK", "Slide")]))==0

mod <- lm(xra ~ TelomereValues + cannum + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, data=phenofemale[sel,])
summary(mod)

##plot residualized XRa
mod <- lm(xra ~ cannum + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, data=phenofemale[sel,])
summary(mod)
phenofemale$XRares <- NA
phenofemale$XRares[sel] <- mod$residuals

mod <- lm(xra ~ TelomereValues, data=phenofemale[sel,])

dat <- data.frame(XRa=phenofemale$XRares[sel], Telomere=phenofemale$TelomereValues[sel])

plot1 <- ggplot(dat, aes(x=Telomere, y=XRa)) + 
  geom_point(color="black", size=3, alpha=0.5)+
  geom_smooth(method=lm)+
  labs(title = "", x = "Telomere length", y = "X-Ra") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 36),
        axis.title = element_text(size = 38),
        plot.title = element_text(size = 38),
        legend.text = element_text(size = 34))


ggsave("./Figures/Figure5C.png", plot1, width = 12, height = 10, dpi = 300)


### Association with methylation age
### CpG data to compute methylation age in TruDiagnostic 
### (Data request varun@trudiagnostic.com) 
########################################################

load("./Data/metAgeTrue.RData")
mm1 <- cbind(rownames(metck), data.frame(metck))
colnames(mm1)[1] <- "ProbeID"

load("./Data/metAgeTrueH.RData")
mm2 <- cbind(rownames(metck), data.frame(metck))
colnames(mm2)[1] <- "ProbeID"

mm <- rbind(mm1,mm2)

##load annotation for CpGs
load("./Data/cpgid.RData")

cpgx <- rownames(cpgid)[which(cpgid$Chromosome=="X")]

##load annotation for CpGs in each clock
load("./Data/coefHorvath.rda")
load("./Data/coefLevine.rda")

metage <- DNAmAge(mm, normalize = FALSE)

metage <- as.data.frame(metage)
rownames(metage) <- rownames(pheno)

phenofemale$levine <- metage[rownames(phenofemale),]$Levine
phenofemale$horvath <- metage[rownames(phenofemale),]$Horvath

sel <-  rowSums(is.na(phenofemale[,c("age", "xra", "TelomereValues", "Bcell","CD4T", "CD8T", "Eos","Mono", "Neu" , "NK", "Slide")]))==0


mod <- glm(xra ~ levine + cannum + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, data=phenofemale[sel,])
summary(mod)

mod <- glm(xra ~ horvath + cannum + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, data=phenofemale[sel,])
summary(mod)


#####################################################
##Association between XRa and Cancer risk (Figure 6A)
##Validation in GSE142536_BLOOD.R
#####################################################

phenofemale$cannum <- as.numeric(factor(phenofemale$Cancer_Diagnosis_any))-1
phenofemale$xraesc <- log(phenofemale$xra*100)

mod <- glm(cannum ~ xraesc+age, 
           family="binomial", data=phenofemale)

mod <- glm(cannum ~ xraesc + age + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide, 
           family="binomial", data=phenofemale)

summary(mod)


png("./Figures/Figure6A.png", width = 6, height = 6, units = 'in', res = 300)
visreg(mod, "xraesc",ylab="Probability of cancer", main="", 
       scale='response', jitter = TRUE, xlab="log(X-Ra*100)", overlay=TRUE)
dev.off()

