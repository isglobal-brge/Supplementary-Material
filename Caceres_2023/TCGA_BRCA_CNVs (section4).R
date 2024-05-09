## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

## Association of X-Ra in breast cancer with somatic mutations
## Data obtained from TCGA_XRa_Calling (phcancer.RData)

#### libraries
library(XRa)
library(curatedTCGAData)
library(TCGAbiolinks)
library(GenomicRanges)
library(arm)
library(survival)
library(qqman)
library(readxl)
library(RTCGA.mutations)


###X-Ra in BRCA
load("./Data/phcancer.RData")

#association of XRa with cancer risk and survival (Figure 4B)
#We select cancer samples only (caco=1)
#######################################

Cancer <- phcancer$BRCA[phcancer$BRCA$caco==1,]
rownames(Cancer) <- Cancer$ID

####CNVs in BRCA for TCGA 
##run only once and save in cnv.BCRtab.all.RData
################################################

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")


GDCdownload(query)
cnv.BCRtab.all <- GDCprepare(query)

save(cnv.BCRtab.all, file="./Data/cnv.BCRtab.all.RData")
##########################################################

load(file="./Data/cnv.BCRtab.all.RData")

cnvs  <- as.data.frame(cnv.BCRtab.all)
seltumor <-  grep("-01A-",cnvs$Sample)
cnvs <- cnvs[seltumor,]
cnvs$Chromosome[cnvs$Chromosome%in%"X"] <- 23 


#define gains
#############

ss <- cnvs$Segment_Mean > log2(3/2)
dups <- cnvs[ss,]

selsubs <- unique(substr(cnvs$Sample, 1, 12))
Cancer <- Cancer[selsubs,]

### Run association between XRa and frequency of gains in 
#sliding windows of 0.5Mb across the genome (Figure 4C)
########################################################

phenoname <- "xra"
gwasdups <- lapply(1:23, function(chr){
  print(chr)
    ss <- dups$Chromosome == chr
    grchr <- GRanges(paste("chr", dups[ss,"Chromosome"], sep=""),
                     IRanges(start=dups[ss,"Start"], 
                             end=dups[ss,"End"]),
                     sub = dups[ss,"Sample"])
    
    mn <- min(start(grchr))
    mx <- max(start(grchr))
    
    #define window
    blocks <- round(seq(mn,mx, by=500000))
    
    #test association between XRa and gain presence in each window
    assocchr <- lapply(1:(length(blocks)-1), function(bb){
      
      block <- GRanges(paste("chr", chr, sep=""),
                       IRanges(start=blocks[bb],
                               end=blocks[bb+1]))
      
      selover <- data.frame(findOverlaps(block,grchr))[,2]
      
      selsubscnv <- unique(as.character(grchr$sub[selover]))
      selsubscnv <- substr(selsubscnv,1, 12)
      
      dat <- Cancer
      dat$cnv <- rep(0,nrow(dat))
      sel <- rownames(dat)%in%selsubscnv
      if(sum(sel)>0) dat[sel,]$cnv <- 1
      
      dat$pheno <- log2(dat[,phenoname])

      assocblock <- summary(bayesglm(
        pheno ~ cnv+age, 
        data=dat))$coeff["cnv",c(1,4)]
      
      data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
                 beta=exp(assocblock[1]), P=assocblock[2], 
                 freq=sum(sel), N=sum(complete.cases(dat)))
      
      
    })
  
  
    assocchr <- do.call(rbind, assocchr)
})
  
gwasdups <- do.call(rbind, gwasdups)

##Compute frequencies
fr <- gwasdups$freq/gwasdups$N
gwasdups.fr <- gwasdups

gwasdups.fr$chr[gwasdups.fr$chr==23] <- "X"
gwasdups.fr$chr <- paste("chr",gwasdups.fr$chr, sep="") 

#select significat results
Pad <- p.adjust(gwasdups.fr$P)
sigdups <- gwasdups.fr[Pad <0.05,]

##map results
load("./Data/geneids.RData")
GRgeneids <- GRanges(geneids)

sigchr <- names(table(sigdups$chr))

sres <- lapply(1:length(sigchr), function(i){
  sc <- sigdups[sigdups$chr==sigchr[i],]
  sigblock <- GRanges(sc[which(min(sc$P)==sc$P),])
  mapg <- as.vector(na.omit(unlist(GRgeneids[subjectHits(findOverlaps(sigblock,GRgeneids)),]$symbol)))
  sigblock$genes  <- paste(mapg, collapse=",")
  sigblock
})

antomost <- as.data.frame(do.call(c, sres))

#write resutls in tables 
write.table(sigdups[order(sigdups$P),], file="./Data/TableS3_dups.txt", quote=FALSE, row=FALSE, col=TRUE, sep="\t")
write.table(antomost[order(antomost$P),], file="./Data/TableS4_dups.txt", quote=FALSE, row=FALSE, col=TRUE, sep="\t")

# Make the Manhattan plot on the gwasResults dataset (Fig 4C gains)
gwasdups$Genes <- rep(NA, length=nrow(gwasdups))
gwasdups$Genes[gwasdups$chr==20 & gwasdups$start==63580664] <- "TNFRSF6B" #
gwasdups$Genes[gwasdups$chr==8 & gwasdups$start==127581254] <- "MYC"
gwasdups$Genes[gwasdups$chr==10 & gwasdups$start==2045792] <- "MIR6072"
gwasdups$Genes[gwasdups$chr==7 & gwasdups$start==38043355] <- "STARD3NL"
gwasdups$Genes[gwasdups$chr==16 & gwasdups$start==9010777] <- "C16orf72"
gwasdups$Genes[gwasdups$chr==5 & gwasdups$start==20515532] <- "CDH18"
gwasdups$Genes[gwasdups$chr==1 & gwasdups$start==156562920] <- "NTRK1"
gwasdups$Genes[gwasdups$chr==17 & gwasdups$start==65650733] <- "CEP112"
gwasdups$Genes[gwasdups$chr==2 & gwasdups$start==39012784 ] <- "MAP4K3"
gwasdups$Genes[gwasdups$chr==12 & gwasdups$start==69051460] <- "FRS2"
gwasdups$Genes[gwasdups$chr==14 & gwasdups$start==22239283] <- "ABHD4"
gwasdups$Genes[gwasdups$chr==21 & gwasdups$start==45854168] <- "COL6A1"
gwasdups$Genes[gwasdups$chr==6 & gwasdups$start==23649661] <- "NRSN1"
gwasdups$Genes[gwasdups$chr==11 & gwasdups$start==47698510] <- "PTPRJ"
gwasdups$Genes[gwasdups$chr==19 & gwasdups$start==28090910] <- "LOC100420587"
gwasdups$Genes[gwasdups$chr==3 & gwasdups$start==5020930 ] <- "ARL8B"


png("./Figures/Figure4C_dups.png", width = 12, height = 6, units = 'in', res = 300)
mh(gwasdups, chr="chr", bp="start", suggestiveline = F, p="P" , 
          genomewideline = -log10(0.05/nrow(gwasdups)), 
          snp = "Genes", annotatePval = 1e-4, ylim=c(0,20), annotateTop=FALSE, 
          cex.axis=1.25, cex.lab=1.3, ylab="-log P", cex=1.5)
dev.off()


#define deletions
#################

ss<- cnvs$Segment_Mean < -1
dels <- cnvs[ss,]

### Run association between XRa and frequency of deletions in 
#sliding windows of 0.5Mb across the genome 
#############################################################

gwasdels <- lapply(1:23, function(chr){
  print(chr)
  ss <- dels$Chromosome == chr
  grchr <- GRanges(paste("chr", dels[ss,"Chromosome"], sep=""),
                   IRanges(start=dels[ss,"Start"], 
                           end=dels[ss,"End"]),
                   sub = dels[ss,"Sample"])
  
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  
  blocks <- round(seq(mn,mx, by=500000))
  
  assocchr <- lapply(1:(length(blocks)-1), function(bb){
    
    block <- GRanges(paste("chr", chr, sep=""),
                     IRanges(start=blocks[bb],
                             end=blocks[bb+1]))
    
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    
    selsubscnv <- unique(as.character(grchr$sub[selover]))
    selsubscnv <- substr(selsubscnv,1, 12)
    
    dat <- Cancer
    dat$cnv <- rep(0,nrow(dat))
    sel <- rownames(dat)%in%selsubscnv
    if(sum(sel)>0) dat[sel,]$cnv <- 1
    
    dat$pheno <- log2(dat[,phenoname])
    
    assocblock <- summary(bayesglm(
      pheno ~ cnv+age, 
      data=dat))$coeff["cnv",c(1,4)]
    
    data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
               beta=exp(assocblock[1]), P=assocblock[2], 
               freq=sum(sel), N=sum(complete.cases(dat)))
    
    
  })
  
  
  assocchr <- do.call(rbind, assocchr)
})


gwasdels <- do.call(rbind, gwasdels)

##compute frequencies
fr <- gwasdels$freq/gwasdels$N
gwasdels.fr <- gwasdels

gwasdels.fr$chr[gwasdels.fr$chr==23] <- "X"
gwasdels.fr$chr <- paste("chr",gwasdels.fr$chr, sep="") 

#select significant results
Pad <- p.adjust(gwasdels.fr$P)
sigdels <- gwasdels.fr[Pad <0.05,]

##Annotate genes
GRgeneids <- GRanges(geneids)

sigchr <- names(table(sigdels$chr))

sres <- lapply(1:length(sigchr), function(i){
  sc <- sigdels[sigdels$chr==sigchr[i],]
  sigblock <- GRanges(sc[which(min(sc$P)==sc$P),])
  mapg <- as.vector(na.omit(unlist(GRgeneids[subjectHits(findOverlaps(sigblock,GRgeneids)),]$symbol)))
  sigblock$genes  <- paste(mapg, collapse=",")
  sigblock
})

antomost <- as.data.frame(do.call(c, sres))

#write resutls in tables 
write.table(sigdels[order(sigdels$P),], file="./Data/TableS3_dels.txt", quote=FALSE, row=FALSE, col=TRUE, sep="\t")
write.table(antomost[order(antomost$P),], file="./Data/TableS4-dels.txt", quote=FALSE, row=FALSE, col=TRUE, sep="\t")


# Make the Manhattan plot on the gwasResults dataset (Fig 4C dels)
gwasdels$Genes <- rep(NA, length=nrow(gwasdels))
gwasdels$Genes[gwasdels$chr==8 & gwasdups$start==3581254] <- "CSMD1" #


# Make the Manhattan plot on the gwasResults dataset
png("./Figures/Figure4C_dels.png", width = 12, height = 6, units = 'in', res = 300)
mh(gwasdels, chr="chr", bp="start", suggestiveline = F, p="P" , 
          genomewideline = -log10(0.05/nrow(gwasdels)), 
          snp = "Genes", annotatePval = 1e-4, ylim=c(0,20), annotateTop=FALSE, 
          cex.axis=1.25, cex.lab=1.3, ylab="-log P")
dev.off()

#######Associations with telomere length data obtained from 
## Barthel, et al 2017 PMID: 28135248 
###########################################################

telomere <- read_excel("./Data/NIHMS893058-supplement-Supplementary_Table_1.xlsx")

telomere$PatientID

TL <- telomere$TL
names(TL) <- telomere$PatientID
TL <- TL[rownames(Cancer)]

dat <- Cancer
dat$TL <- TL

summary(glm(xra ~ TL + age, data=dat))

##Association telomere length and gains in RTEL1 
phenoname <- "TL"
gwasdupsTL <- lapply(8, function(chr){
  print(chr)
  ss <- dups$Chromosome == chr
  grchr <- GRanges(paste("chr", dups[ss,"Chromosome"], sep=""),
                   IRanges(start=dups[ss,"Start"], 
                           end=dups[ss,"End"]),
                   sub = dups[ss,"Sample"])
  
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  
  blocks <- round(seq(mn,mx, by=500000))
  
  assocchr <- lapply(1:(length(blocks)-1), function(bb){
    
    block <- GRanges(paste("chr", chr, sep=""),
                     IRanges(start=blocks[bb],
                             end=blocks[bb+1]))
    
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    
    selsubscnv <- unique(as.character(grchr$sub[selover]))
    selsubscnv <- substr(selsubscnv,1, 12)
    
    dat$cnv <- rep(0,nrow(dat))
    sel <- rownames(dat)%in%selsubscnv
    if(sum(sel)>0) dat[sel,]$cnv <- 1
    
    dat$pheno <- log2(dat[,phenoname])
    
    assocblock <- summary(bayesglm(
      pheno ~ cnv+age, 
      data=dat))$coeff["cnv",c(1,4)]
    
    data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
               beta=exp(assocblock[1]), P=assocblock[2], 
               freq=sum(sel), N=sum(complete.cases(dat)))
    
    
  })
  
  
  assocchr <- do.call(rbind, assocchr)
})

gwasdupsTL <- do.call(rbind, gwasdupsTL)

Pad <- p.adjust(gwasdupsTL$P)
sigdups <- gwasdupsTL[Pad <0.05,]
sigdups

####Somatic mutations data obtained from RTCGA
##############################################

maf <- mutationsTCGA(BRCA.mutations)

setgenes <- sort(unique(maf$Hugo_Symbol))
subnames <- unique(substr(maf$bcr_patient_barcode, 1, 12)) 

mutmat <- lapply(subnames, function(sub){
  wichsub <- grep(sub, maf$bcr_patient_barcode)
  as.numeric(setgenes%in%maf$Hugo_Symbol[wichsub])})

mutmat <- do.call(rbind, mutmat)

colnames(mutmat) <- setgenes
rownames(mutmat) <- subnames

Cancermut <- Cancer[subnames, ]

asmut <- lapply(1:ncol(mutmat), function(x){
  print(x)
  
  dat <- data.frame(xra=Cancermut$xra, mut=mutmat[,x],
                    age =Cancermut$age)
  summary(bayesglm(xra~mut+age, data=dat))$coef[2,c(1,4)]})

asmut <- do.call(rbind,asmut)
rownames(asmut) <- colnames(mutmat)

#select frequent mutations
selmut <- colMeans(mutmat)>0.05
colMeans(mutmat)[selmut]

#correct for multiple comparisons
asmutsel <- asmut[selmut,]
asmutsel[p.adjust(asmutsel[,2])<0.05,]

colMeans(mutmat)[c("CDH1","TP53")]

