## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##X-Ra in blood 
#from the MESA study GSE56046
############################# 

#### libraries
library(XRa)
library(SummarizedExperiment)
library(arm)
library(GEOquery)
library(limma)
library(sva)
library(EnhancedVolcano)
library(ggplot2)



#MONOCYTES data from chrX in 
#GSE56046_methylome_normalized.txt.gz
######################################

#load phenos
gsm <- getGEO("GSE56046",destdir ="./Data", AnnotGPL =TRUE)[[1]]
phenobb <- pData(phenoData(gsm))
subidspheno <- substr(phenobb$title, 1, 6)
rownames(phenobb) <-  subidspheno


#Transcription data
gsm <- getGEO("GSE56047",destdir ="./data", AnnotGPL =TRUE)[[1]]

expr <- exprs(gsm)
genesIDs <- fData(gsm)
rownames(expr) <- genesIDs$'Gene symbol'

phenotrans <- pData(phenoData(gsm))
subidsphenotrans <- substr(phenotrans$title, 1, 6)
cmsubs <- intersect(subidsphenotrans, subidspheno)
nmscom <- rownames(phenotrans)
names(nmscom) <- subidspheno
nmscom <- nmscom[cmsubs]

phenobb <- phenobb[names(nmscom), ]

#phenotypes do not include sex, but a combination of sex with cohort
#we infered sex with transcriptomic data from Y chromosome
####################################################################

load("./Data/geneids.RData")
iny <- unlist(geneids$symbol[geneids$seqnames=="chrY"])
iny <- iny[complete.cases(iny)]

#function to compute relative expression from chr Y
Rexp <- function(expr, genenm)
{
  selChr <- rownames(expr)%in%genenm
  selnoChr <- !rownames(expr)%in%genenm
  
  #remove transcipts with zero counts across all individuals
  subdataChr <- expr[selChr,]
  selin <-rowMeans(subdataChr!=0)!=0
  subdataChr <- subdataChr[selin,]
  
  subdatanoChr <- expr[selnoChr,]
  selin <- rowMeans(subdatanoChr!=0)!=0
  subdatanoChr <- subdatanoChr[selin,]
  
  #compute mean of log2 counts
  datChrMean <- colMeans(log2(subdataChr-min(subdataChr, na.rm=TRUE)+1), na.rm=TRUE)
  datnoChrMean <- colMeans(log2(subdatanoChr-min(subdatanoChr, na.rm=TRUE)+1), na.rm = TRUE)
  
  #estimate the relative expression of Y with respect to genome
  datChrMean-datnoChrMean
}


#computes Ry 
phenobb$Ry <- Rexp(expr[,nmscom], iny)

clus <- letters[as.numeric(phenobb$"racegendersite:ch1")]
clas <- tapply(phenobb$Ry, factor(clus),  mean)
sel <- clas< -0.35
phenobb$sex <- clus%in%names(clas)[sel]

#confirm sex classification
plot(phenobb$age, phenobb$Ry, col=as.numeric(phenobb$sex)+1)



###Methylation data was downloaded from 
###https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56046/suppl/GSE56046%5Fmethylome%5Fnormalized.txt.gz
###CpGs from chromosome X selected and saved in metXMonocytes.RData
###################################################################

load(file="./Data/metXMonocytes.RData")
met <- metX
met <- (2^met)/((2^met+1))

cpgmono <- rownames(met)

metfemale <- t(met)[phenobb$sex, ]
phenofemale <- phenobb[phenobb$sex, ]


#compute XRa
phenofemale$xra <- XRa(metfemale, cpgmono)

# Data for Plot
phenofemale$age <- as.numeric(phenofemale$`age:ch1`)
phenofemale$cohort <- phenofemale$`racegendersite:ch1`
phenofemale$well <- phenofemale$`well:ch1`
phenofemale$tcell <- as.numeric(phenofemale$`tcell:ch1`)
phenofemale$neutro <- as.numeric(phenofemale$`neutro:ch1`)
phenofemale$nkcell <- as.numeric(phenofemale$`nkcell:ch1`)
phenofemale$bcell <- as.numeric(phenofemale$`bcell:ch1`)

mod <- lm(xra ~ age + cohort + tcell + neutro + nkcell + bcell, 
          data=phenofemale)
summary(mod)

phenofemale$age.dic <-  phenofemale$age > 65

mod3 <- lm(xra ~ age.dic + cohort + tcell + neutro + nkcell + bcell, 
          data=phenofemale)
summary(mod)

save(mod3, file="./Data/modMonocites.RData")


##########################################################
### Differential expression on XRa (Figure S11, Table S7)
##########################################################

exprfemale <- expr[,nmscom]
exprfemale <- exprfemale[,phenobb$sex]
exprfemale <- exprfemale[!is.na(rowSums(exprfemale)),]

mod0 <- model.matrix( ~ age + cohort +tcell+ neutro+nkcell+bcell  , data = phenofemale)
mod <- model.matrix( ~ xra + age + cohort +tcell+ neutro+nkcell+bcell, data = phenofemale)

#this gives 0 svas 
#ns <- num.sv(exprfemale, mod, method="be")
#ss <- sva(exprfemales, mod, mod0, n.sv=ns)$sv

modss <- mod

fit <- lmFit(exprfemale, modss)
fit <- eBayes(fit)

tt <- topTable(fit, coef="xra", number=Inf)


lab <- tt[1:7,"ID"]
lab <- lab[lab!=""]

tt$genes <- rep("", nrow(tt))
tt$genes[1:7] <-  tt$ID[1:7]

png("./Figures/FigureS11.png", width = 6, height = 6, units = 'in', res = 300)
EnhancedVolcano(tt, lab = tt$genes, 
                selectLab  = lab, 
                x = 'logFC', y = 'P.Value', 
                xlim=c(-15, 15),
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                legendPosition = 'bottom',
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'black',
                title = "",
                subtitle = "Differential expression",
                hline = 0.01,
                legendLabSize  =10)

dev.off()     

### Enrichment of XRa associations genes under XCI 
##################################################

#list of genes that  are always inactive and always escape 
esc <- read.delim("./Data/Suppl.Table.1.csv", as.is=TRUE, sep=";", header=TRUE, skip=1)

#genes that are inactive
inactive <- esc$Gene.name[esc$Combined.XCI.status=="inactive"]

#genes that always escape
always <- esc$Gene.name[esc$Combined.XCI.status=="escape"]

#associations in inactivated genes
ttinac <- tt[which(tt$ID%in%inactive),]
n <- length(unique(ttinac$ID))
k <- length(unique(ttinac$ID[ttinac[, "P.Value"]<1e-3]))
N <- length(unique(tt$ID))
K <- length(unique(tt$ID[tt[, "P.Value"]<1e-3]))

write.table(ttinac[ttinac[, "P.Value"]<1e-3,], file="./Data/TableS7.txt", row=FALSE, quote=FALSE, sep="\t")

#enrichment test
fisher.test(matrix(c(N,K,n,k), ncol=2))

#significant asociation of XIST
tt[tt$ID%in%"XIST", c(1,2,5)]

