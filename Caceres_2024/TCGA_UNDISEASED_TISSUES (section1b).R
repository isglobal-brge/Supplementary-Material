## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##Frequency and distribution of X-Ra in 12 undiseased tissues 
## Data obtained from TCGA_XRa_Calling (phcancer.RData)

#libraries
library(chrXRa)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(ggplot2)
library(GenomicRanges)
library(ggforestplot)


###X-Ra calling
load("./Data/phcancer.RData")

#methylation data from the TCGA were downloaded using curatedTCGAData and TCGAbiolinks, 
#and CpGs from chr X were filtered. 
#Data was stored in the metTCGA object.

load("./Data/metTCGA.RData")
#case control status of samples
splitcaco <- lapply(met0, function(x) 
  list(cancer=grep("01A", colnames(x)), 
       healthy=grep("11A", colnames(x))))


#select cancer data with >5 normal tissue samples (12 selected tissues in total)
ln <- lapply(splitcaco, function(x)  length(x$healthy))
sel <- unlist(ln)>5 

cancertypes <- names(met0)
cancersel <- cancertypes[sel]


####Permutation test (Figure 2D) saved in cpglevelsperm.RData
#Run onl once################################################

cpglevelsperm <- lapply(cancersel, function(cc){
  print(cc)
  met_control <- met0[[cc]][,grep("11A", colnames(met0[[cc]]))]
  colnames(met_control) <- substr(colnames(met_control), 1, 12)
  met_control <- met_control[,colnames(met_control)%in%phcancer[[cc]]$ID]

  metfemale <- t(met_control)
  cpgmono  <- rownames(met_control)
  
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

save(cpglevelsperm, file="./Data/cpglevelsperm.RData")
##################################################################

load("./Data/cpglevelsperm.RData")

obs <- unlist(lapply(cancersel, function(cc){ cpglevelsperm[[cc]]$obs}))
names(obs) <- cancersel
              
pval <- unlist(lapply(cancersel, function(cc){ mean(cpglevelsperm[[cc]]$pemtest< obs[cc], na.rm=TRUE)}))


tissues <- c("Breast", "Bladder", "Colon", "Head and Neck", "Kidney",
             "Liver", "Lung (LUAD)", "Lung (LUSC)", "Pancreas", "Rectum",
             "Thyroid", "Uterus") 

png("./Figures/Ficure2D.png", width = 6, height = 6, units = 'in', res = 300)

par(mfrow=c(3,4))
for(cc in 1:12){
hist(cpglevelsperm[[cc]]$pemtest, 
     xlim=c(0.05, 0.20), ylim=c(0,300), col="blue", border = "transparent", 
     main=tissues[cc],xlab="X-Ra", cex.lab=1.3, cex.axis=1.2, cex.main=1.5)
lines(x=c(obs[cc], obs[cc]),y=c(0,200), col="red", lwd=1.5)
points(obs[cc], 200, col="red", pch=16)
}
  
dev.off()

##Population frequency of 2-allele demethylation calls in XCI CpGs (Figure 2A)
################################################################################

#gather data needed for X-Ra
load("./Data/cpgid.RData") #annotation data

#Supplementary data from Tukiainen T, et al. Nat Publ Gr. 2017;550.
esc <- read.delim("./Data/Suppl.Table.1.csv", as.is=TRUE, sep=";", header=TRUE, skip=1)
inactive <- esc$Gene.name[esc$Combined.XCI.status=="inactive"]
cpgidinactive <- rownames(cpgid)[cpgid$Gene_Symbol%in%inactive]

always <- esc$Gene.name[esc$Combined.XCI.status=="escape"]
cpgidalways <- rownames(cpgid)[cpgid$Gene_Symbol%in%always]


freqHypo <- lapply(cancersel, function(cc){
  met_control <- met0[[cc]][,grep("11A", colnames(met0[[cc]]))]
  colnames(met_control) <- substr(colnames(met_control), 1, 12)
  met_control <- met_control[,colnames(met_control)%in%phcancer[[cc]]$ID]
  
  metfemale <- t(met_control)
  cpgmono  <- rownames(met_control)
  mm <-   metfemale[,colnames(metfemale)%in%cpgidinactive]

  distcpgs <- lapply(1:nrow(mm), function (i){
    out <- mm[i,] 
    colnames(mm)[(out < 0.2)]
  })

  distcpgs
})

names(freqHypo) <- cancersel

freqHypoCpG <- unlist(freqHypo)
freqHypoCpG <- freqHypoCpG[!is.na(freqHypoCpG)]


nsb <- lapply(cancersel, function(cc){
  met_control <- met0[[cc]][,grep("11A", colnames(met0[[cc]]))]
  colnames(met_control) <- substr(colnames(met_control), 1, 12)
  met_control <- met_control[,colnames(met_control)%in%phcancer[[cc]]$ID]
  
  metfemale <- t(met_control)
  cpgmono  <- rownames(met_control)
  mm <-   metfemale[,colnames(metfemale)%in%cpgidinactive]
  nrow(mm)})

nsubs <- sum(unlist(nsb))

mean(table(freqHypoCpG)/nsubs<0.2)


png("./Figures/Figure2A.png", width = 6, height = 6, units = 'in', res = 300)

hist(table(freqHypoCpG)/nsubs, br=200, col="blue", border = "transparent",
     xlab="Population frequency",
     ylab="Number of CpGs",
     main="2-allele demethylation call in CpGs under XCI \n (12 undiseased tissues)", 
     cex.lab=1.3, cex.axis=1.5, cex.main=1.5)

dev.off()


##where do the CpGs with 2-allele demathylation contributing to X-Ra come from? (Figure S1)
###########################################################################################

tb <- table(freqHypoCpG)

selcpgs <- names(tb[tb/nsubs < 0.2])

cpgscoords <- as.numeric(cpgid[names(tb[selcpgs]),]$Genomic_Coordinate)/10^6
popfreq <- as.numeric(tb[selcpgs])/nsubs

png("./Figures/FigureS1.png", width = 6, height = 6, units = 'in', res = 300)
plot(cpgscoords,
     popfreq, pch=16, 
     ylab="Population frequency",
     xlab="ChrX coordinates (Mb)", 
     cex=0.5,
     main="2-allele demethylation call in CpGs under XCI \n (12 healthy tissues)" ) 

dev.off()

##Enrichment of XCI CpGs with 2-allele demethylation in chromatin states (Figure 2B)
####################################################################################

#annotation 15-chromatin states released from the ROADMAP Epigenomics Mapping Consortium (ChromHMM v1.10)
load("./Data/crom15.RData")

rownames(crom15)<- crom15$HT12v4.ArrayAddress
selcpgs <- names(tb)
gr <- crom15[cpgidinactive,]

nmcr <- names(gr)[8:22]

cromatine <- lapply(nmcr, function(x){
crstate <- gr[[x]]
names(crstate) <- rownames(gr)

crstate <- crstate[!is.na(crstate)]
N <- length(crstate)
K <- sum(crstate)
n <- length(selcpgs)
x <- sum(selcpgs%in%names(crstate[crstate]))

unlist(c(fisher.test(matrix(c(N,K,n,x), ncol=2))[c("estimate", "conf.int", "p.value")], K=K, x=x))

})

cromatine <- as.data.frame(do.call(rbind, cromatine))
cromatine$p.adjust <- p.adjust(cromatine[,4])
rownames(cromatine) <- nmcr

# Check that nomenclature nmcr coincides with description 
# https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html

desc <- c("Active TSS",
          "Flanking Acive TTS",
          "Transcr. at gene 5' and 3'",
          "Weak transcription",
          "Strong transcription",
          "Genic enhancers",
          "Enhancers",
          "ZNF genes & repeats",
          "Heterochromatin",
          "Flanking Bivalent TSS/Enh",
          "Bivalent/Poised TSS",
          "Bivalent Enhancer",
          "Repressed PolyComb",
          "Weak Repressed PolyComb",
          "Quiescent/Low")



#Enrichment in DNase hypersensitivity states
#data downloaded from 
#https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
###############################################################################

files <- c("./Data/E028-DNase.hotspot.fdr0.01.peaks.v2.bed.gz",
           "./Data/E029-DNase.hotspot.fdr0.01.peaks.v2.bed.gz",
           "./Data/E094-DNase.hotspot.fdr0.01.peaks.v2.bed.gz",
           "./Data/E097-DNase.hotspot.fdr0.01.peaks.v2.bed.gz",
           "./Data/E100-DNase.hotspot.fdr0.01.peaks.v2.bed.gz",
           "./Data/E109-DNase.hotspot.fdr0.01.peaks.v2.bed.gz")

names(files) <- c("Breast", "Monocytes", "Gastric", "Ovary", "Muscle", "Small Intestine")


DNaseenrich <- lapply(1:6, function(i){
DNasa <- read.delim(files[i], header=FALSE)
DNasa <- DNasa[DNasa[,1]=="chrX",]

grchr <- GRanges(DNasa[,1],
                 IRanges(start=DNasa[,2], 
                         end=DNasa[,3]))

cpgidinactive <- rownames(cpgid)[cpgid$Gene_Symbol%in%inactive]

bb <- cpgidinactive[1]

cpgrg <- GRanges("chrX",
                 IRanges(start=as.numeric(cpgid[cpgidinactive, 3]),
                         end=as.numeric(cpgid[cpgidinactive, 3])))

selover <- data.frame(findOverlaps(cpgrg,grchr))

cpgDNAse <- cpgidinactive[selover[,1]]

selcpgs <- names(tb)

N <- length(cpgidinactive)
K <- length(cpgDNAse)
n <- length(selcpgs)
x <- sum(selcpgs%in%cpgDNAse)
  
unlist(c(fisher.test(matrix(c(N,K,n,x), ncol=2))[c("estimate", "conf.int", "p.value")], K=K, x=x))
})


DNaseenrich <- as.data.frame(do.call(rbind, DNaseenrich))
DNaseenrich$p.adjust <- p.adjust(DNaseenrich[,2])
rownames(DNaseenrich) <- names(files)

desc2 <- paste(paste0("DNase (", names(files) , sep=""), ")", sep="")


#join 15 chromatin enrichment with DNase enrichment
enrich <- data.frame(name=desc, beta=log(cromatine[,1]), 
                     se=(log(cromatine[,2])-log(cromatine[,1]))/qnorm(0.975), 
                     pvalue=cromatine[,4])

enrich <- rbind(enrich, data.frame(name=desc2, beta=log(DNaseenrich[,1]), 
                                   se=(log(DNaseenrich[,2])-log(DNaseenrich[,1]))/qnorm(0.975), 
                                   pvalue=DNaseenrich[,4]))

enrich <- enrich[order(enrich[,2], decreasing = TRUE),]

# Forestplot
p1 <- forestplot(
  df = enrich,
  estimate = beta,
  logodds = TRUE,
  title = "Enrichment of 2-allele demathylation in XCI CpGs",
  xlab = "Odds ratio (95% CI)")+
  theme_minimal() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12))+
  geom_point(aes(x = exp(enrich$beta)), size = 4)


ggsave("./Figures/Figure2B.png", p1, width = 6, height = 6, dpi = 300, bg="white")

#####Are low population frequency (<0.05) of 2-allele demethylation informative on X-Ra? (Figure 2E)
####################################################################################################

distcpgs <- freqHypo[[cc]]

met_control <- met0[[cc]][,grep("11A", colnames(met0[[cc]]))]
colnames(met_control) <- substr(colnames(met_control), 1, 12)
met_control <- met_control[,colnames(met_control)%in%phcancer[[cc]]$ID]

metfemale <- t(met_control)
cpgmono  <- rownames(met_control)
mm <-   metfemale[,colnames(metfemale)%in%cpgidinactive]
mm <- ncol(mm)

tb <- table(freqHypo[[cc]])


#X-Ra based on rare hypomethylation events
xralow <- sapply(distcpgs, function(x){
  selcpgs <- names(tb)[which((tb/nsubs<0.05))]
  xx <- unlist(x)[!is.na(x)]
  sum(xx%in%selcpgs)/mm})

xra <- XRa(metfemale, colnames(metfemale))

plot(xra, xralow, pch=16)

df <- data.frame(x = xra, y = xralow)

# Create dataframe
df <- data.frame(x = xra, y = xralow)

# Create the plot with regression line and confidence bands
plot <- ggplot(df, aes(x = x, y = y)) +
  geom_point(shape = 16, size=3) +  # plot points
  theme_minimal() + # Optional: customize the theme
  xlab("X-Ra") +  # Label for the x-axis
  ylab("X-Ra (fr < 0.05)") + # Label for the y-axis
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))

ggsave("./Figures/Figure2E.png", plot, width = 6, height = 8, dpi = 300, bg="white")



