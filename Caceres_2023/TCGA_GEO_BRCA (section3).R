## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##Frequency and distribution of X-Ra in breast cancer
## Data obtained from TCGA_XRa_Calling (phcancer.RData)

#### libraries
library(survival)
library(survminer)
library(GEOquery)
library(chrXRa)
library(contsurvplot)
library(riskRegression)
library(survival)
library(ggplot2)
library(pammtools)
library(smplot2)
library(data.table)
library(irr)
library(lmerTest)
library(curatedTCGAData)
library(TCGAbiolinks)

###X-Ra in BRCA
load("./Data/phcancer.RData")

#association of XRa with cancer risk and survival (Figure 4A)
#We select cancer samples only (caco=1)
########################################

cancerdat <- phcancer$BRCA[phcancer$BRCA["caco"]==1,]
rownames(cancerdat) <- cancerdat$ID

summary(lm(xra ~ caco + age, data=phcancer$BRCA))$coefficients["caco",c(1,2)]

cancerdat$months <- cancerdat$times/30

mod1 <- coxph(Surv(months, event) ~ xra + age, data=cancerdat, x=TRUE)
summary(mod1)


plot2 <- plot_surv_area(time="months",
               status="event",
               variable="xra",
               data=cancerdat,
               model=mod1,
               horizon=seq(0, .7, 0.05))+
        theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))

png("./Figures/Figure4A1.png", width = 6, height = 6, units = 'in', res = 300)
plot2
dev.off()



###Validation of less survival with X-Ra in TNBC GSE78754
#########################################################

gsm <- getGEO("GSE78754",destdir ="./Data", AnnotGPL =TRUE)

gsm <- gsm[[1]]

phenobb <- pData(phenoData(gsm))
met <- exprs(gsm)

feat <- fData(gsm)

smmet <- feat$CHR=="X"
  
metX <- met[smmet,]

sel <- which(phenobb$`gender:ch1`=="female")

#obtain beta values and phenotypes for females 
metfemale <- t(metX)[sel, ]
phenofemale <- phenobb[sel, ]

phenofemale$xra <- XRa(metfemale, colnames(metfemale))


####Associations with phenotypes
phenofemale$age <- as.numeric(phenofemale$"age at diagnosis:ch1") 
phenofemale$event <- as.numeric(factor(phenofemale$"survival:ch1"))-1
phenofemale$event[phenofemale$event==2]<- NA
phenofemale$time <- as.numeric(phenofemale$"month of follow up/to death:ch1")


mod1 <- coxph(Surv(time, event) ~ xra + age, data=phenofemale, x=TRUE)
summary(mod1)


plot2 <- plot_surv_area(time="time",
                        status="event",
                        variable="xra",
                        data=phenofemale,
                        model=mod1,
                        horizon=seq(0, .7, 0.05))+
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))

png("./Figures/Figure4A2.png", width = 6, height = 6, units = 'in', res = 300)
plot2
dev.off()


##association with relapse in TNBC GSE141441 (Figure S9)
###########################################################

gsm <- getGEO("GSE141441",destdir ="./Data", AnnotGPL =TRUE)
gsm <- gsm[[1]]

phenobb <- pData(phenoData(gsm))
met <- exprs(gsm)

load("./Data/cpgid.RData")

cpgx <- rownames(cpgid)[which(cpgid$Chromosome=="X")]
smmet <- intersect(rownames(met), cpgx)

metX <- met[smmet,]
metfemale <- t(metX)
phenofemale <- phenobb

#compute XRa
phenofemale$xra <- XRa(metfemale, colnames(metfemale))

#phenotypes
dat <- data.frame(relapse=as.numeric(phenofemale$"relapse:ch1"),
                  XRa=phenofemale$xra*10,
                  chemo=phenofemale$`study_arm`,
                  surg=as.numeric(phenofemale$`treatment:ch1`=="SurgOnly")
)


summary(glm(relapse ~ XRa + chemo, data=dat[dat$surg==0,], family="binomial"))

df <- dat[dat$surg==0,]

plot4 <- ggplot(data=df, mapping=aes(x=XRa,y=relapse)) + geom_point() +
  stat_smooth(method="glm", method.args=list(family=binomial)) +
  labs(title = "TNBC", x = "X-Ra (1e-1)", y = "Prob. of relapse") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))+
  geom_text(aes(x = 2.5, y = 0.25, label = "OR=2.20, P=0.04"), size=8)

ggsave("./Figures/FigureS9.png", plot4, width = 12, height = 10, dpi = 300)



##Association with visit for TNBC under chemotherapy treatment (GSE184159)
####Figure S10

##Download data locally and select CpGs in chr X, save methylation data in metXGSE184159.RData

##Run only once###############################################
met <- read.delim("./Data/GSE184159_beta_detp.csv.gz", sep=",")

load("./Data/cpgid.RData")

met <- met[, -grep("p_", names(met))]
ids <- met[,1]

met <- met[, -1]
rownames(met) <- ids

cpgx <- rownames(cpgid)[which(cpgid$Chromosome=="X")]

smmet <- intersect(rownames(met), cpgx)

metX <- met[smmet,]

save(metX, file="./Data/metXGSE184159.RData")
######################################################


load("./Data/metXGSE184159.RData")

#load phenotypes
gsm <- getGEO("GSE184159",destdir ="./Data", AnnotGPL =TRUE)[[1]]
phenobb <- pData(phenoData(gsm))

##merge data
subidsmet <- sapply(strsplit(colnames(metX), "_"), function(x) paste(x[c(2,3)], collapse=""))
subidspheno <- sapply(strsplit(phenobb$title, "_"), function(x) paste(x[c(1,3)], collapse=""))

identical(subidsmet, subidspheno)


##identify subject and visit
phenobb$sub <- sapply(strsplit(phenobb$title, "_"), function(x) x[1])
phenobb$time <- sapply(strsplit(phenobb$title, "_"), function(x) x[3])

table(phenobb$`timepoint:ch1`, phenobb$time)

###select data for repated visits
sel1 <- phenobb$time%in%c("A","B","C")
sel2 <-names(table(phenobb$sub[sel1]))[table(phenobb$sub[sel1])>1]
sel <- sel1 & phenobb$sub%in%sel2

metfemale <- as.matrix(metX[,sel])
phenofemale <- phenobb[sel, ]

table(phenofemale$sub)

#compute Xra
cpgids <- rownames(metfemale)
metfemale <- t(metfemale)

phenofemale$xra <- XRa(metfemale, cpgids)

##phenotypes
phenofemale$resp <- as.numeric(factor(phenofemale$"response:ch1", labels = c(1,0,1)))
phenofemale$time <- as.numeric(as.factor(phenofemale$time))

phenofemale2visits <- phenofemale[phenofemale$time%in%c(1,2),]


##arrange data for test-retest ICC calculation
sel1 <- phenofemale$time%in%c(1,2)
sel2 <-names(table(phenofemale$sub[sel1]))[table(phenofemale$sub[sel1])>1]
sel <- sel1 & phenofemale$sub%in%sel2

phenofemale2visits <- phenofemale[sel,]

oo <- order(phenofemale2visits$sub)
opheno <- phenofemale2visits[oo,]

data <- matrix(opheno$xra[opheno$time%in%c(1,2)], ncol=2, byrow = TRUE)

##check
matrix(opheno$title[opheno$time%in%c(1,2)], ncol=2, byrow = TRUE)

icc(data, model="twoway", type = "consistency", unit = "single")

##Test association with visit adjusted by response
mod <- lmerTest::lmer(xra ~ time + resp  + (1 |sub), 
                      data=phenofemale2visits)

summary(mod)

# Create a new data frame with predicted values 
dat <- data.frame(time = phenofemale2visits$time,
                  sub = phenofemale2visits$sub,
                  xra = phenofemale2visits$xra)



mod <- lmer(xra ~ time  + (1|sub),  data=dat)

summary(mod)

pred <- data.frame(dat, predicted=predict(mod))

plot3 <- ggplot(pred, aes(x = time, y = xra)) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = "blue")  +
  labs(x = "Time", y = "X-Ra", color = " ", title = "TNBC Chemotherapy (GSE184159)") +
  guides(color = guide_legend(title = " ")) +
  scale_x_continuous(breaks = c(1, 2), labels =c("A", "B")) +
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))+
  geom_line(data = dat, aes(x = time, y = xra, group = sub),size = 0.15, color ="blue")

ggsave("./Figures/FigureS10.png", plot3, width = 6, height = 6, dpi = 300, bg="white")



###################################
##association with TNBC (Figure 4B)
###################################

query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)


# ER content 
er.cols <- grep("^er",colnames(clinical.BCRtab.all$clinical_patient_brca))

dclin <- data.frame(clinical.BCRtab.all$clinical_patient_brca)[-(1:2),]
rownames(dclin) <- dclin$bcr_patient_barcode
dclin <- dclin[rownames(cancerdat),]

## ER status
er <- dclin$er_status_by_ihc 
er <- factor(er, labels = c("NA", "NA", "0", "1"))
cancerdat$er <- as.numeric(as.character(er))

## HER status
her <- dclin$her2_status_by_ihc 
her <- factor(her, labels = c("NA", "NA", "NA", "NA","0", "1"))
cancerdat$her <- as.numeric(as.character(her))

## PER status
per <- dclin$pr_status_by_ihc 
per <- factor(per, labels = c("NA", "NA", "0", "1"))
cancerdat$per <- as.numeric(as.character(per))

summary(glm(er ~ xra + age, family="binomial", data=cancerdat))
summary(glm(per ~ xra + age, family="binomial", data=cancerdat))
summary(glm(her ~ xra + age, family="binomial", data=cancerdat))

cancerdat$TN <- cancerdat$er==0 & cancerdat$per==0 & cancerdat$her==0
summary(glm(TN ~ xra + age, family="binomial", data=cancerdat))

mod <- glm(TN ~ xra + age, family="binomial", data=cancerdat)
summary(mod)

png("./Figures/Figure4b1.png", width = 6, height = 6, units = 'in', res = 300)
visreg(mod, "xra", ylab="Probability of TNBC", main="", 
       scale='response', jitter = TRUE, xlab="X-Ra", overlay=TRUE,
       cex.lab=1.4, cex.axis=1.4)
dev.off()


coxph(Surv(times, event) ~  xra + TN + age, data=cancerdat)
coxph(Surv(times, event) ~   TN + age, data=cancerdat)


####Validation of association between X-Ra and TNBC in GSE225845 (Figure4B2)
gsm <- getGEO("GSE225845",destdir ="./Data", AnnotGPL =TRUE)
gsm <- gsm[[1]]
phenobb <- pData(phenoData(gsm))

###download methylation data and select X cromsome CpGs save them in metXTNBC 
###run only once##############################################################

met <- fread("./Data/GSE225845_tumors_normalized_betas.txt.gz")
selout <- is.na(met[1,]) | (met[1,]==0)

met <- met[, !selout, with=FALSE]

ids <- met[,1]
ids <- unlist(ids)
met <- met[, -1]
subids <- colnames(met)

met <- as.matrix(met)
rownames(met) <- ids

load("./Data/cpgid.RData")

cpgx <- rownames(cpgid)[which(cpgid$Chromosome=="X")]
smmet <- intersect(colnames(met), cpgx)

metX <- t(met[,smmet])
metX <- metX
metX <- as.data.frame(metX)

save(metX, file="./Data/metXTNBC.RData")
####################################################

load(file="./Data/metXTNBC.RData")

ids <- colnames(metX)

phenobb <- pData(phenoData(gsm))
rownames(phenobb) <- phenobb$`methylation id (basenames):ch1`
phenobb <- phenobb[ids,]

metfemale <- t(metX)
phenofemale <- phenobb

metfemale <- metfemale[phenobb$"Sex:ch1"=="F", ]
phenofemale <- phenofemale[phenobb$"Sex:ch1"=="F", ]


#compute XRa
phenofemale$xra <- XRa(metfemale, colnames(metfemale))
phenofemale$type <- factor(phenofemale$`molecular subtype:ch1`, labels = c("HER2+/HER+", "HER2+/HER+", "TNBC"))
phenofemale$TNBC <- as.numeric(phenofemale$type)-1

mod <- glm(TNBC ~ xra , family="binomial", data=phenofemale)
summary(mod)

png("./Figures/Figure4B2.png", width = 6, height = 6, units = 'in', res = 300)
visreg(mod, "xra", ylab="Probability of TNBC", main="", 
       scale='response', jitter = TRUE, xlab="X-Ra", overlay=TRUE, 
       cex.lab=1.4, cex.axis=1.4)
dev.off()


