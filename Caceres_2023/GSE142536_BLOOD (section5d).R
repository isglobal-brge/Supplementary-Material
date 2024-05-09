## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##Validation of association between XRa in blood and cancer 
#from a longitudinal study GSE142536
########################################################## 

#### libraries
library(chrXRa)

#### libraries
library(SummarizedExperiment)
library(arm)
library(GEOquery)
library(irr)
library(methylclock)
library(lmerTest)
library(data.table)
library(ggpubr)
library(visreg)


###Methylation data was downloaded from 
###https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142536/suppl/GSE142536%5FRyoichi%5Fprocessed%5Fsignals.csv.gz
###CpGs from chromosome X were selected and saved in metXGSE142536.RData
###Methylation age was computed and saved in metageGSE142536.RData
############################################################################


##Run only 
#once#################################################################
met <- read.delim("./Data/GSE142536_Ryoichi_processed_signals.csv.gz", sep=",")

load("./Data/cpgid.RData")

met <- met[, -grep("Pval", names(met))]
ids <- met[,1]

met <- met[, -1]
rownames(met) <- ids

cpgx <- rownames(cpgid)[which(cpgid$Chromosome=="X")]

smmet <- intersect(rownames(met), cpgx)

metX <- met[smmet,]

save(metX, file="./Data/metXGSE142536.RData")

#compute methylation age
mm <- cbind(rownames(met), met)
colnames(mm)[1] <- "ProbeID"

metage <- DNAmAge(mm, normalize = FALSE)

save(metage, file="./Data/metageGSE142536.RData")
##################################################################

#load phenos
gsm <- getGEO("GSE142536",destdir ="./Data", AnnotGPL =TRUE)[[1]]
phenobb <- pData(phenoData(gsm))

#load methylation
load("./Data/metXGSE142536.RData")

#check IDS
#identical(phenobb$label_ch1, colnames(metX))

metfemale <- as.matrix(metX[,phenobb$`gender:ch1`=="F"])
phenofemale <- phenobb[phenobb$`gender:ch1`=="F", ]

cpgids <- rownames(metfemale)
metfemale <- t(metfemale)


### Comute XRa and arrage phenotypes for associations
phenofemale$xra <- XRa(metfemale, cpgids)
phenofemale$age <- as.numeric(phenofemale$`age:ch1`)
phenofemale$cancer <- as.numeric(phenofemale$`cancer:ch1`=="Y")
phenofemale$cd4t <- as.numeric(phenofemale$`proportion of  cd4t cell:ch1`)
phenofemale$cd8t <- as.numeric(phenofemale$`proportion of cd8t cell:ch1`)
phenofemale$bcell <- as.numeric(phenofemale$`proportion of bcell:ch1`)
phenofemale$gran <- as.numeric(phenofemale$`proportion of gran:ch1`)
phenofemale$nkcell <- as.numeric(phenofemale$`proportion of nk cell:ch1`)
phenofemale$mono <- as.numeric(phenofemale$`proportion of mono:ch1`)
phenofemale$time <- as.numeric(as.factor(phenofemale$`time point:ch1`))


  
#Compute intraclass correlatioon coefficient
data <- matrix(phenofemale$xra, ncol=3, byrow = TRUE)
icc(data, model = "twoway", type = "agreement", unit = "single")


phenofemale$timesq <- phenofemale$time^2
phenofemale$surgery <- phenofemale$`surgery:ch1`
#adjustements by age and surgery do no change the estimates

phenofemale$sub <- sapply(strsplit(phenofemale$label_ch1, "\\."), function(x) x[[1]])

mod <- lmerTest::lmer(xra ~ time + cancer + age + cd4t + cd8t +  bcell + gran + nkcell + mono +(1|sub), 
          data=phenofemale)

summary(mod)



# Create a new data frame with predicted values for each level of 'cancer'
dat <- data.frame(time = phenofemale$time,
                  cancer = factor(phenofemale$cancer, labels = c("No Cancer", "Cancer")),
                  sub = phenofemale$sub,
                  xra = phenofemale$xra,
                  age = phenofemale$age)
  
              

mod <- lmer(xra ~ time + cancer + age + (1|sub), 
                   data=dat)


pred <- data.frame(dat, predicted=predict(mod))

plot1 <- ggplot(pred, aes(x = time, y = xra, color = cancer)) +
  geom_smooth(aes(fill = cancer), method = "lm", se = TRUE) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"), guide = "none") +
  labs(x = "Time", y = "X-Ra", color = " ") +
  guides(color = guide_legend(title = " ")) +
  scale_x_continuous(breaks = c(1, 2, 3), labels =c("BL", "PoD1", "PoD4/7")) +
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))+
  geom_line(data = dat, aes(x = time, y = xra, color = cancer, group = sub),size = 0.15)
            
          
# Display the plot
ggsave("./Figures/Figure6B.png", plot1, width = 8, height = 6, dpi = 300, bg="white")


###########################################
#compare with biological aging (Figure S13)
###########################################

load("./Data/metageGSE142536.RData")

dfmetage <- data.frame(metage)[phenobb$`gender:ch1`=="F",]

iid <- sapply(strsplit(dfmetage$id, "\\."), function(x) x[[1]])
identical(phenofemale$sub, iid)

phenofemaleage <- data.frame(phenofemale, dfmetage)

plot(phenofemaleage$xra, phenofemaleage$Levine)
cor.test(phenofemaleage$xra, phenofemaleage$Levine)


mod <- lmerTest::lmer(Levine ~ time + age + cancer  + cd4t + cd8t +  bcell + gran + nkcell + mono + (1|sub), 
                      data=phenofemaleage)

summary(mod)



# Create a new data frame with predicted values for each level of 'cancer'
dat <- data.frame(time = phenofemaleage$time,
                  cancer = factor(phenofemaleage$cancer, labels = c("No Cancer", "Cancer")),
                  sub = phenofemaleage$sub,
                  Levine = phenofemaleage$Levine,
                  age = phenofemaleage$age)



mod <- lmer(Levine ~ time + cancer + age + (1|sub), 
            data=dat)


pred <- data.frame(dat, predicted=predict(mod))

plot <- ggplot(pred, aes(x = time, y = Levine, color = cancer)) +
  geom_smooth(aes(fill = cancer), method = "lm", se = TRUE) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"), guide = "none") +
  labs(x = "Time", y = "Levine's methylation age", color = " ") +
  guides(color = guide_legend(title = " ")) +
  scale_x_continuous(breaks = c(1, 2, 3), labels =c("BL", "PoD1", "PoD4/7")) +
  theme_minimal()+
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))+
  geom_line(data = dat, aes(x = time, y = Levine, color = cancer, group = sub),size = 0.15)


# Display the plot
ggsave("./Figures/FigureS13.png", plot, width = 8, height = 6, dpi = 300, bg=TRUE)


##Association of xRa in blood with cancer from the GSE237036 study
##################################################################


###Methylation data was downloaded from 
###https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237036/suppl/GSE237036%5Fmatrix%5Fprocessed.txt.gz
###CpGs from chromosome X were selected and saved in metXGSE142536.RData
###Methylation age was computed and saved in metGSE237036.RData
############################################################################

##Run only 
#once#################################################################
met <- read.delim("./Data/GSE237036_matrix_processed.txt.gz")

load("./Data/cpgid.RData")
cpgx <- rownames(cpgid)[which(cpgid$Chromosome=="X")]

ids <- met[,1]
met <- met[, -1]
rownames(met) <- ids

smmet <- intersect(rownames(met), cpgx)
metX <- met[smmet,]
save(metX, file="./Data/metGSE237036.RData")

#Immunce cell count.
library(FlowSorted.Blood.EPIC)
common <- intersect(rownames(met), IDOLOptimizedCpGs)
met <- met[common, ]
common <- IDOLOptimizedCpGs %in% rownames(met) 

propEPIC2<-projectCellType_CP (
  met[IDOLOptimizedCpGs[common],],
  IDOLOptimizedCpGs.compTable, contrastWBC=NULL, nonnegative=TRUE,
  lessThanOne=FALSE)

percEPIC2<-round(propEPIC2*100,1)

save(percEPIC2, file="./Data/percEPIC2GSE237036.RData")

############################################

#load phenos
gsm <- getGEO("GSE237036",destdir ="./Data", AnnotGPL =TRUE)[[1]]
phenobb <- pData(phenoData(gsm))

#load immune cell count
load("./Data/percEPIC2GSE237036.RData")

#load biological age
load("./Data/metageGSE237036.RData")
metage <- as.data.frame(metage)


#load methylation
load("./Data/metGSE237036.RData")

#check IDS
#identical(phenobb$title, colnames(metX))

idssmet <- gsub("_", "", colnames(metX))
idssmet <- gsub("c", "C", idssmet)
idphenos <- sapply(strsplit(phenobb$title, " "), function(x) x[[4]])

comn <- intersect(idssmet, idphenos)

rownames(phenobb) <- idphenos
colnames(metX) <- idssmet
rownames(percEPIC2) <- idssmet
rownames(metage) <- idssmet

phenobb <- phenobb[comn, ]
metX  <- metX[, comn]
percEPIC2 <- percEPIC2[comn, ]
metage <- metage[comn, ]

metfemale <- as.matrix(metX)
phenofemale <- phenobb

cpgids <- rownames(metfemale)
metfemale <- t(metfemale)


### Comute XRa and arrage phenotypes for associations
phenofemale$xra <- log(XRa(metfemale, cpgids)*100)

phenofemale$cancer <- as.numeric(phenofemale$`disease state:ch1`=="BC")
phenofemale$CD8T <- percEPIC2[,"CD8T"]
phenofemale$CD4T <- percEPIC2[,"CD4T"]
phenofemale$NK <- percEPIC2[,"NK"]
phenofemale$Bcell <- percEPIC2[,"Bcell"]
phenofemale$Mono <- percEPIC2[,"Mono"]
phenofemale$age <- metage$Levine

mod <- glm(cancer ~ xra+age, 
           data=phenofemale, family="binomial")

summary(mod)

visreg(mod, "xra", ylab="Probability of cancer", main="", 
       scale='response', jitter = TRUE, xlab="X-Ra (10e-2)", overlay=TRUE)


mod <- glm(cancer ~ xra+CD8T+CD4T+NK+Bcell+Neu+age, 
           data=phenofemale, family="binomial")

summary(mod)



formeta <- data.frame(ES=c(0.7656834,  1.78896), SD=c(0.3178651, 1.30142))


matdat <- meta::metagen(formeta[,1], formeta[,2],
                        studlab= c("TrueDiag", "GSE237036"),
                        level.ci = 0.95, sm="OR")


matdat
