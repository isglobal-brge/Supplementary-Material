## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##Frequency and distribution of X-Ra in 12 tumors 
## Data obtained from TCGA_XRa_Calling (phcancer.RData)

#### libraries
library(survival)
library(cvAUC)
library(pROC)
library(AF)
library(limma)
library(sva)
library(EnhancedVolcano)
library(clusterProfiler)
library(curatedTCGAData)
library(TCGAbiolinks)
library(org.Hs.eg.db)
library(biomaRt)
library(visreg)
library(arm)


###X-Ra calling
load("./Data/phcancer.RData")

cancersel <- names(phcancer)

########Multiple cancers

####Transcriptomic data downloaded with curatedTCGAData
## DE expression analysis of X-Ra in each cancer type
## adjusted by SVA and age.
## associations for each tumor are saved in assocs and saved in assocsSVAns.RData
#run only once###################################################################

assocs <- list()

for(type in names(phcancer)){
  
  canrna  <- curatedTCGAData(type, "RNASeq2GeneNorm", FALSE, version="2.0.1")
  
  mergerna <- lapply(1:length(canrna), function(n){
    x <- canrna[[n]]
    data.frame(assay(x), check.names = FALSE)
  })
  
  
  selgenes <- lapply(mergerna, function(x) rownames(x))
  if(length(selgenes)==1){ 
    selgenes <- unlist(selgenes)
  }else{
    selgenes <- do.call(intersect, selgenes)}
  
  mergerna <- lapply(mergerna, function(x) x[selgenes,])
  counts <- do.call(cbind, mergerna)
  
  #We select cancer samples only 
  counts <- counts[,grep("01A",colnames(counts))]
  colnames(counts) <- substr(colnames(counts), 1, 12)
  
  filter <- apply(counts, 1, function(x) mean(x>15)>0.25)
  counts <- counts[filter,]
  
  cancerdat <- phcancer[[type]][phcancer[[type]]["caco"]==1,]
  rownames(cancerdat) <- cancerdat$ID
  
  commonsubs <- intersect(rownames(cancerdat), colnames(counts))
  
  genesIDs <- sapply(strsplit(rownames(counts), "\\|") , function(x)  x[[1]] )
  
  nmsgenes <- sapply(strsplit(genesIDs, "\\|"), function(x) x[1])
  
  rownames(counts) <- nmsgenes
  
  cancercounts <- counts[,commonsubs]
  cancerdat <-  cancerdat[commonsubs, ]
  
  cancerdat$eff <- as.numeric(cancerdat$xra)
  cancerdat$age <- as.numeric(cancerdat$age)
  
  cc <- complete.cases(cancerdat[c("eff", "age")])
  cancerdat <- cancerdat[cc, ]
  cancercounts<-cancercounts[,cc]
  cancercounts <- as.matrix(cancercounts)
  
  mod0 <- model.matrix(~ age, data=cancerdat)
  mod <- model.matrix(~ eff + age, data=cancerdat)
  
  ns <- num.sv(cancercounts, mod, method="be")
  ss <- svaseq(cancercounts, mod, mod0, n.sv=ns)$sv
  colnames(ss) <- paste("cov", 1:ncol(ss), sep="")
  
  
  modss <- cbind(mod, ss)
  
  design <- model.matrix(~  eff + age, data=cancerdat)
  v <- voom(cancercounts, design=design)
  
  fit <- lmFit(v, modss)
  fit <- eBayes(fit)
  
  assocs[[type]] <- fit
  save(assocs, file="./Data/assocsSVAns.RData")
}
################################################

###Volcano plots (Figure S3)
#############################

load("./Data/assocsSVAns.RData")

for(x in names(assocs)){
  
  #  x <- names(assocs)[1]  
  tt <- topTable(assocs[[x]], coef="eff", number=Inf)
  
  tt$genes <- rep("", nrow(tt))
  tt$genes[1:7] <-  tt$ID[1:7]
  lab <- tt[1:7,"ID"]
  lab <- lab[lab!=""]
  
  fl <- paste("./Figures/","Volcano.png", sep=x)
  png(fl, width = 6, height = 6, units = 'in', res = 300)
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
                  title = x,
                  subtitle = "Differential expression",
                  hline = 0.01,
                  legendLabSize  =10)
  dev.off()
  
}

###Perform meta analysis accross cancer types.
## Save results in matdatSVAns and select significant results(Table S2)
##run only once#########################################################

bse <- lapply(assocs, function(x){ 
  beta <- x$coefficients[, "eff"]
  se <- (sqrt(x$s2.post)*x$stdev.unscaled)[, "eff"]
  cbind(beta,se)
})

nms <- lapply(bse, function(x) rownames(x))
nms <- table(unlist(nms))
nms <- names(nms)[nms==length(names(phcancer))]

bb <- lapply(bse, function(x) x[nms,1])
bb <- do.call(cbind, bb)

ss <- lapply(bse, function(x) x[nms,2])
ss <- do.call(cbind, ss)

metaa <- lapply(1:nrow(bb), function(x){
  matdat <- try(meta::metagen(bb[x,], 
                              ss[x,],
                              studlab= colnames(bb),
                              level.ci = 0.95), silent = TRUE)
  matdat
})

gnms <- rownames(bb)
sel <- sapply(metaa, function(x) class(x)[1]=="metagen")

metaa <- metaa[sel]
gnms <- gnms[sel]
names(metaa) <- gnms

save(metaa, file="./Data/matdatSVAns.RData")

##############################################

load("./Data/matdatSVAns.RData")

gres <- lapply(metaa, function(x) c(x$statistic.random, x$pval.random))
gres <- do.call(rbind, gres)
colnames(gres) <- c("beta", "pval")

gres <- gres[order(gres[,2]),]
padj <- p.adjust(gres[,2])
gres <- cbind(gres, padj)

sigmat <- gres[gres[,"padj"]<0.05,]

head(sigmat)

write.table(sigmat, file="./DATA/TableS2.txt", quote=FALSE, sep="\t")


##Extract meta-analysis for XIST (Figure 3B)
############################################

png("./Figures/Figure3B.png", width = 8, height = 8, units = 'in', res = 300)

meta::forest(metaa[["XIST"]], layout="JAMA",  
             leftlabs=c("TCGA", "Log2FC 95% CI"), xlab="XIST Log2FC Differential Expression",
             title="")

dev.off()

#####Enrichment analysis (Figure S4)
####################################

selgenes <- rownames(sigmat)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneids<- getBM(c('hgnc_symbol','ensembl_gene_id','entrezgene_id','chromosome_name','start_position','end_position'), filters = 'hgnc_symbol',  values=selgenes, mart=ensembl)
selgenes <- geneids[geneids$hgnc_symbol%in%selgenes, "entrezgene_id"]

#run enrichment in GO
GOmet1 <- enrichGO(gene = selgenes, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01, pAdjustMethod="fdr")

png("./Figures/FigureS4.png", width = 6, height = 6, units = 'in', res = 300)
clusterProfiler::dotplot(GOmet1) 
dev.off()

###############################################
####Association between cancer status and X-Ra
####Cancer status on X-Ra (Figure S5)
###############################################


formeta <-lapply(phcancer, function(x){ 
                 ##reescale X-Ra
                 x$xraes <- x$xra*10^2
                 summary(glm(caco ~ xraes  + age,
                             family="binomial",
                             data=x))$coefficients["xraes",c(1,2)]
                 })

formeta <- do.call(rbind, formeta)

matdat <- meta::metagen(formeta[,1], formeta[,2],
                        studlab= cancersel,
                        level.ci = 0.95, sm="OR")


png("./Figures/FigureS5.png", width = 8, height = 8, units = 'in', res = 300)
meta::forest(matdat, layout="JAMA",  
       leftlabs=c("TCGA", "OR 95% CI"), xlab="Odds Ratio for Cancer risk",
       main="Effect of cancer status on X-Ra")
dev.off()




####compare permutation test of XRa between cancer and undiseased tissues
####Figure S6############################################################

####Permutation test saved in cpglevelspermCancer.RData##################
#Run only once###########################################################

cpglevelsperm <- lapply(cancersel, function(cc){
  print(cc)
  met_cancer <- met0[[cc]][,grep("01A", colnames(met0[[cc]]))]
  colnames(met_cancer) <- substr(colnames(met_cancer), 1, 12)
  met_cancer <- met_cancer[,colnames(met_cancer)%in%phcancer[[cc]]$ID]
  
  metfemale <- t(met_cancer)
  cpgmono  <- rownames(met_cancer)
  
  ig <- IGlevels(metfemale, cpgmono)
  obs <- mean(ig <0.2, na.rm=TRUE)
  
  
  permtest <- sapply(1:1000, function(i){
    sn <- sample(1:ncol(metfemale))
    spgsam <- cpgmono[sn]
    igsam <- IGlevels(metfemale, spgsam)
    mean(igsam<0.2, na.rm=TRUE)})
  
  list(obs=obs, pemtest=permtest)
})

names(cpglevelsperm) <- cancersel

save(cpglevelsperm, file="./Data/cpglevelspermCancer.RData")
##################################################################
##################################################################

##permuation in cancer
load("./Data/cpglevelspermCancer.RData")
cpglevelspermCancer <- cpglevelsperm

obscancer <- unlist(lapply(cancersel, function(cc){ cpglevelspermCancer[[cc]]$obs}))
names(obscancer) <- cancersel

pval <- unlist(lapply(cancersel, function(cc){ mean(cpglevelspermCancer[[cc]]$pemtest< obscancer[cc], na.rm=TRUE)}))

##permuation in undiseased tissue
load("./Data/cpglevelsperm.RData")
obs <- unlist(lapply(cancersel, function(cc){ cpglevelsperm[[cc]]$obs}))
names(obs) <- cancersel

pval <- unlist(lapply(cancersel, function(cc){ mean(cpglevelsperm[[cc]]$pemtest< obs[cc], na.rm=TRUE)}))


tissues <- c("Breast", "Bladder", "Colon", "Head and Neck", "Kidney",
             "Liver", "Lung (LUAD)", "Lung (LUSC)", "Pancreas", "Rectum",
             "Thyroid", "Uterus") 

##Differences between X-Ra to the null in undiseased tissues and cancer
difXra <- sapply(1:12, function(cc){
  tt <- t.test(cpglevelsperm[[cc]]$pemtest-obs[cc], cpglevelspermCancer[[cc]]$pemtest-obscancer[cc])
  c(meandif=tt$estimate[1]-tt$estimate[2])})


##correlation between association estimates of XRa and cancer and differences between
##Xra and null distributions

cor.test(formeta[,1], difXra)

png("./Figures/FigureS6.png", width = 6, height = 6, units = 'in', res = 300)

par(mfrow=c(3,4))

for(cc in 1:12){
  
  hist(cpglevelsperm[[cc]]$pemtest, 
       xlim=c(0.05, 0.25), ylim=c(0,300), col="lightblue", border = "transparent", 
       main=tissues[cc],xlab="X-Ra", cex.lab=1.3, cex.axis=1.2, cex.main=1.5)
  hist(cpglevelspermCancer[[cc]]$pemtest, add=TRUE, col="salmon", border = "transparent")
  lines(x=c(obs[cc], obs[cc]),y=c(0,200), lwd=1.5)
  lines(x=c(obscancer[cc], obscancer[cc]),y=c(0,200), lwd=1.5)
  points(obs[cc], 200, col="blue", pch=16)
  points(obscancer[cc], 200, col="red", pch=16)
  
  lines(c(mean(cpglevelsperm[[cc]]$pemtest), obs[cc]), c(300,300), col="blue")
  lines(c(mean(cpglevelsperm[[cc]]$pemtest), mean(cpglevelsperm[[cc]]$pemtest)), c(300,290), col="blue")
  lines(c(obs[cc], obs[cc]), c(300,290), col="blue")
  
  
  lines(c(mean(cpglevelspermCancer[[cc]]$pemtest), obscancer[cc]), c(275,275), col="red")
  lines(c(mean(cpglevelspermCancer[[cc]]$pemtest), mean(cpglevelspermCancer[[cc]]$pemtest)), c(275,265), col="red")
  lines(c(obscancer[cc], obscancer[cc]), c(275,265), col="red")
  
}

dev.off()









##effect on cancer status of samples of the proportion of 2-allele demethylation 
##in CpGs of escapees (Xesc) (Figure S7)
#######################################

formeta <-lapply(phcancer, function(x){ 
  x$xesces <- x$xesc*10^2
  summary(glm(caco ~ xesces  + age,
              family="binomial",
              data=x))$coefficients["xesces",c(1,2)]
})

formeta <- do.call(rbind, formeta)

matdat <- meta::metagen(formeta[,1], formeta[,2],
                        studlab= cancersel,
                        level.ci = 0.95, sm="OR")


png("./Figures/FigureS7.png", width = 8, height = 8, units = 'in', res = 300)
meta::forest(matdat, layout="JAMA",  
             leftlabs=c("TCGA", "OR 95% CI"), xlab="Odds Ratio for Cancer risk",
             main="Effect of cancer status on Xesc")
dev.off()


#Correlation between the proportion of 1-allele methylation in CpGs of escapees 
#and the proportion of 2-allele methylation of XCI CpGs(Figure S8)
##############################################################################

corxraxep <-lapply(names(phcancer), function(x){ 
  ##reescale X-Ra
  dat <- phcancer[[x]]
  cbind(dat$xra*10^2, dat$xescinterm*10^2)

})

corxraxep <- do.call(rbind,corxraxep)

cor.test(corxraxep[,1],corxraxep[,2])


png("./Figures/FigureS8.png", width = 8, height = 8, units = 'in', res = 300)
plot(corxraxep[,1],corxraxep[,2], 
     xlab="2-allele demethylation of XCI CpGs (X-Ra)", 
     ylab="1-allele methylation of CpGs in escapees")

dev.off()

##########################################################################################
## Cancer risk models with proportion of 2-allele demethylation in XCI CpGs, proportion of
## 1-allele methylation of CpGs in escapees and age (Figure 3C)
###############################################################

formeta <-lapply(names(phcancer), function(x){ 
  ##reescale X-Ra
  dat <- phcancer[[x]]
  dat$xraes <- dat$xra*10^2
  dat$xescintermes <- dat$xescinterm*10^2
  model <- glm(caco ~ xescintermes  + xraes + age,
              data=dat[complete.cases(dat),], family = "binomial")
  
  c(summary(model)$coefficients["xraes",c(1,2)],
    summary(model)$coefficients["xescintermes",c(1,2)])

  })


formeta <- do.call(rbind, formeta)

matdat <- meta::metagen(formeta[,1], formeta[,2],
                        studlab= cancersel,
                        level.ci = 0.95, sm="OR")


png("./Figures/Figure3Ca.png", width = 8, height = 8, units = 'in', res = 300)
meta::forest(matdat, layout="JAMA",  
             leftlabs=c("TCGA", "OR 95% CI"), xlab="Odds Ratio for Cancer risk",
             main="Effect of cancer status on X-Ra")
dev.off()


matdat <- meta::metagen(formeta[,3], formeta[,4],
                        studlab= cancersel,
                        level.ci = 0.95, sm="OR")


png("./Figures/Figure3Cb.png", width = 8, height = 8, units = 'in', res = 300)
meta::forest(matdat, layout="JAMA",  
             leftlabs=c("TCGA", "OR 95% CI"), xlab="Odds Ratio for Cancer risk",
             main="Effect of cancer status on X-Ra")
dev.off()



#which Cpgs in escapees contribute to 1-allele methylation  of CpGs in escapees in undiseased tissue? 
####################################################################################################

load("./Data/cpgid.RData") #annotation data

#Supplementary data from Tukiainen T, et al. Nat Publ Gr. 2017;550.
esc <- read.delim("./Data/Suppl.Table.1.csv", as.is=TRUE, sep=";", header=TRUE, skip=1)

always <- esc$Gene.name[esc$Combined.XCI.status=="escape"]
cpgidalways <- rownames(cpgid)[cpgid$Gene_Symbol%in%always]

##control
freq1allelecont <- lapply(cancersel, function(cc){
  met_control <- met0[[cc]][,grep("11A", colnames(met0[[cc]]))]
  colnames(met_control) <- substr(colnames(met_control), 1, 12)
  met_control <- met_control[,colnames(met_control)%in%phcancer[[cc]]$ID]
  
  metfemale <- t(met_control)
  cpgmono  <- rownames(met_control)
  mm <-   metfemale[,colnames(metfemale)%in%cpgidalways]
  
  distcpgs <- lapply(1:nrow(mm), function (i){
    out <- mm[i,] 
    colnames(mm)[(out > 0.4 & out <0.6)]
  })
  
  distcpgs
})

names(freq1allelecont) <- cancersel

freq1allelecont <- unlist(freq1allelecont)
freq1allelecont <- freq1allelecont[!is.na(freq1allelecont)]


nsb <- lapply(cancersel, function(cc){
  met_control <- met0[[cc]][,grep("11A", colnames(met0[[cc]]))]
  colnames(met_control) <- substr(colnames(met_control), 1, 12)
  met_control <- met_control[,colnames(met_control)%in%phcancer[[cc]]$ID]
  
  metfemale <- t(met_control)
  cpgmono  <- rownames(met_control)
  mm <-   metfemale[,colnames(metfemale)%in%cpgidalways]
  nrow(mm)})

nsubscont <- sum(unlist(nsb))

mean(table(freq1allelecont)/nsubscont<0.2)

tbcont <- table(freq1allelecont)/nsubscont

#which Cpgs in escapees contribute to 1-allele methylation of CpGs in escapees in cancer tissue? 
###############################################################################################

freq1allele <- lapply(cancersel, function(cc){
  met_cancer <- met0[[cc]][,grep("01A", colnames(met0[[cc]]))]
  colnames(met_cancer) <- substr(colnames(met_cancer), 1, 12)
  met_cancer <- met_cancer[,colnames(met_cancer)%in%phcancer[[cc]]$ID]
  
  metfemale <- t(met_cancer)
  cpgmono  <- rownames(met_cancer)
  mm <-   metfemale[,colnames(metfemale)%in%cpgidalways]
  
  distcpgs <- lapply(1:nrow(mm), function (i){
    out <- mm[i,] 
    colnames(mm)[(out > 0.2 &  out < 0.6)]
  })
  
  distcpgs
})

names(freq1allele) <- cancersel

freq1allele <- unlist(freq1allele)
freq1allele <- freq1allele[!is.na(freq1allele)]


nsb <- lapply(cancersel, function(cc){
  met_cancer <- met0[[cc]][,grep("01A", colnames(met0[[cc]]))]
  colnames(met_cancer) <- substr(colnames(met_cancer), 1, 12)
  met_cancer <- met_cancer[,colnames(met_cancer)%in%phcancer[[cc]]$ID]
  
  metfemale <- t(met_cancer)
  cpgmono  <- rownames(met_cancer)
  mm <-   metfemale[,colnames(metfemale)%in%cpgidalways]
  nrow(mm)})

nsubs <- sum(unlist(nsb))

mean(table(freq1allele)/nsubs<0.2)

tb <- table(freq1allele)/nsubs

hist(tb, br=100, col="blue")

selcpgs <- names(tb[tb < 0.05])

##############################################################################
###is the frequency of 1-allele methylation in gemetolog Cpgs higher in cancer 
###than in undiseased tissue?
#############################

#homologs activity measures (PMID: 32790864)
XYgametologs<- matrix(c("RPS4Y1","RPS4X",
                   "KDM5D","KDM5C",
                   "CYorf15A","TXLNG",
                   "DDX3Y","DDX3X",
                   "EIF1AY", "EIF1AX",
                   "USP9Y", "USP9X",
                   "UTY", "KDM6A",
                   "ZFY", "ZFX"), nrow=2)


cpghomol <- rownames(cpgid)[cpgid$Gene_Symbol%in%XYgametologs[2,]]

N <- length(tbcont)
K <- sum(cpghomol%in%names(tbcont))            
n <- length(tb)
x <- sum(cpghomol%in%names(tb))
prop.table(matrix(c(N,K,n,x), ncol=2), margin = 2)
fisher.test(matrix(c(N,K,n,x), ncol=2))


####################################################################################
#are loss of function mutations in XCI associated with 2-allele demethylation in XCI
#and muations in escapees with 1-allele methylation in escapees in cancer?
####################################################################################

####LOF
load("./Data/geneids.RData")

esc <- read.delim("./Data/Suppl.Table.1.csv", as.is=TRUE, sep=";", header=TRUE, skip=1)
always <- esc$Gene.name[esc$Combined.XCI.status=="escape"]
inactive <- esc$Gene.name[esc$Combined.XCI.status=="inactive"]

#mutation data for the TCGA downloaded from
#https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.somaticsniper_snv.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
maf <- read.delim("./Data/GDC-PANCAN.somaticsniper_snv.tsv.gz")

maf$Variant_Classification <- maf$effect

#create a single dataframe with XRa and Xesc calling across all tumors
Cancer <- lapply(names(phcancer), function(cc){
  data.frame(phcancer[[cc]][phcancer[[cc]]$caco==1,], type=cc)}) 

Cancer <- do.call(rbind, Cancer)       
rownames(Cancer) <- Cancer$ID


##create a LOF mutations matrix subjectsXgenes
mutations <- names(table(maf$effect))

#select LOF mutations
sel <- mutations %in% c("coding_sequence_variant", 
                        "incomplete_terminal_codon_variant;coding_sequence_variant", 
                        "splice_acceptor_variant", 
                        "splice_acceptor_variant;NMD_transcript_variant", 
                        "splice_donor_variant", 
                        "splice_donor_variant;NMD_transcript_variant", 
                        "start_lost", 
                        "start_lost;splice_region_variant", 
                        "stop_gained", 
                        "stop_gained;NMD_transcript_variant", 
                        "stop_gained;splice_region_variant", 
                        "stop_lost")

mafLOF <- maf[sel,]

#select common subjects in mutation data and phenotype data, and genes in XCI and escapees
subnames <- unique(substr(mafLOF$Sample_ID, 1, 12)) 
cmnames <- subnames[subnames%in%rownames(Cancer)]
setgenes <- sort(unique(mafLOF$gene))
setgenes <- setgenes[setgenes %in% c(always, inactive)] 

#create mutation matrix
mutmat <- lapply(cmnames, function(sub){
  wichsub <- grep(sub, mafLOF$Sample_ID)
  as.numeric(setgenes%in%mafLOF$gene[wichsub])})

mutmat <- do.call(rbind, mutmat)
colnames(mutmat) <- setgenes
rownames(mutmat) <- cmnames

#which subjects have at least one mutation in XCI or in escapees
lofinactive <- rowSums(mutmat[,colnames(mutmat)%in%inactive])>0
lofalways <- rowSums(mutmat[,colnames(mutmat)%in%always])>0

Cancermutlof <- Cancer[cmnames, ]
xra <- Cancermutlof$xra
xescinterm <- Cancermutlof$xescinterm
type <- Cancermutlof$type

datmut <-   data.frame(xra, xescinterm, lofinactive, lofalways, type=factor(type))  

mod1 <- bayesglm(lofinactive~xra+type, family="binomial", data=datmut)
mod2 <- bayesglm(lofalways~xescinterm+type, family="binomial", data=datmut)


png("./Figures/Figure3D.png", width = 16, height = 8, units = 'in', res = 300)
par(mfrow=c(1,2))
visreg(mod1, "xra",main="Probability of LOF in any XCI gene", ylab="", 
       scale='response', jitter = TRUE, xlab="Proportion of 2-allele demetylation in XCI (X-Ra)", overlay=TRUE)
visreg(mod2, "xescinterm",main="Probability of LOF in any escapee", ylab="", 
       scale='response', jitter = TRUE, xlab="Propoortion of 1-allele metylation in escapees (X-esc)", overlay=TRUE)

dev.off()

########################################
###Association between X-Ra and survival
########################################


#remove cancer studies with only 0.05 individuals with events
selin <- sapply(phcancer, function(x) prop.table(table(x[,"event"]))[2])
selin <- selin > 0.05

formeta <-lapply(phcancer[selin], function(x){ 
  summary(coxph(Surv(times, event) ~ xra + age, data = x[x$caco==1,]))$coefficients["xra",c(1,3)]
})


formeta <- do.call(rbind, formeta)

matdat <- meta::metagen(formeta[,1], formeta[,2],
                        studlab= cancersel[selin],
                        level.ci = 0.95, sm="OR")


matdat

################################################
#####Attributable risk of X-Ra on cancer status.
################################################

##pooled AUC to dichotomize X-Ra at the Youden distance
phpooled <- do.call(rbind, phcancer)
datroc <- data.frame(test =phpooled$xra, 
                     response=phpooled$caco, 
                     type=as.factor(rep(1:12,sapply(phcancer, nrow))),
                     age = phpooled$age)  
                                                          

out <- cvAUC(datroc$test, datroc$response) #compute ROC
ci.cvAUC(datroc$test, datroc$response)

rocobj <- roc(datroc$response, datroc$test)
youden <- coords(rocobj, "best", best.method="youden")
youden

#########################################
##Attributable risk for dichotomous X-Ra
#########################################

datroc$dictest <- datroc$test > youden$threshold 

tb <- table(datroc$dictest, datroc$type)
tb <- round(prop.table(tb,2),2)
colnames(tb) <- cancersel
tb

model <- glm(response ~ dictest +age + type, dat=datroc, family = "binomial")
summary(model)

or <- exp(summary(model)$coeff["dictestTRUE",1])
condtab <- prop.table(table(datroc$response, datroc$dictest),1)

afad <- condtab[2,2]*(or-1)/or

af <- AFglm(model, data= datroc, exposure="dictest", 
            clusterid="type", case.control = TRUE)

int <- c(af$AF.est + qnorm(0.975)*sqrt(af$AF.var),
  af$AF.est - qnorm(0.975)*sqrt(af$AF.var))

c(af$AF.est, int)

