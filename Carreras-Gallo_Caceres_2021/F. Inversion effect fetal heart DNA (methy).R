###########################################################################
## DIFFERENTIAL ANALYSIS BY INVERSION IN FETAL HEART DNA (meta-analysis) ##
###########################################################################

# Load libraries 
library(scoreInvHap)
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(cowplot)
library(limma)
library(SummarizedExperiment)
library(Gviz)
library(stringr)

# Load methylation
load("results/methylation/finalQC_files/2021-03-08/gset.autosomic.Rdata")

#######################
## Define inversions ##
#######################

vcf_paths <- dir("results/VariantCalling/SNV/", pattern = ".gz$", full.names = TRUE)

## inv-8p23.1
range <- inversionGR["inv8_001"]
seqlevelsStyle(range) <- "NCBI"

## Subset VCFs to inv-8p23.1 region
a <- lapply(vcf_paths, function(x) {
  vcf <- readVcf(x, genome = "hg19", 
                 param = ScanVcfParam(which = range))
  id <- strsplit(basename(x), "_")[[1]][2]
  writeVcf(vcf, paste0("scripts/Natalia_ISGlobal/tempVCFs/", id, ".inv86.vcf.gz"), index = TRUE)
})

## Merge inv-8p23.1 vcfs
system("bcftools merge --force-samples scripts/Natalia_ISGlobal/tempVCFs/*inv8.vcf.bgz -m none -o scripts/Natalia_ISGlobal/tempVCFs/merged.inv8.vcf.gz -O z")
vcf <- readVcf("scripts/Natalia_ISGlobal/tempVCFs/merged.inv8vcf.gz", genome = "hg19")

check <- checkSNPs(vcf)
vcf.filt <- check$genos
inv16sc <- scoreInvHap(SNPlist = vcf.filt, inv = "inv8_001")
save(inv8sc, file = "scripts/Natalia_ISGlobal/inv8Genos.scoreInvHap.Rdata")


## inv-16p11.2
range <- inversionGR["inv16_009"]
seqlevelsStyle(range) <- "NCBI"

## Subset VCFs to inv-16p11.2 region
a <- lapply(vcf_paths, function(x) {
  vcf <- readVcf(x, genome = "hg19", 
                 param = ScanVcfParam(which = range))
  id <- strsplit(basename(x), "_")[[1]][2]
  writeVcf(vcf, paste0("scripts/Natalia_ISGlobal/tempVCFs/", id, ".inv16.vcf.gz"), index = TRUE)
})

## Merge inv-16p11.2 vcfs
system("bcftools merge --force-samples scripts/Natalia_ISGlobal/tempVCFs/*inv16.vcf.bgz -m none -o scripts/Natalia_ISGlobal/tempVCFs/merged.inv16.vcf.gz -O z")
vcf <- readVcf("scripts/Natalia_ISGlobal/tempVCFs/merged.inv16.vcf.gz", genome = "hg19")

check <- checkSNPs(vcf)
vcf.filt <- check$genos
inv16sc <- scoreInvHap(SNPlist = vcf.filt, inv = "inv16_009")
save(inv16sc, file = "scripts/Natalia_ISGlobal/inv16Genos.scoreInvHap.Rdata")

## inv-17q21.31
range <- inversionGR["inv17_007"]
seqlevelsStyle(range) <- "NCBI"

## Subset VCFs to inv-17q21.31 region
a <- lapply(vcf_paths, function(x) {
  vcf <- readVcf(x, genome = "hg19", 
                 param = ScanVcfParam(which = range))
  id <- strsplit(basename(x), "_")[[1]][2]
  writeVcf(vcf, paste0("scripts/Natalia_ISGlobal/tempVCFs/", id, ".inv17.vcf.gz"), index = TRUE)
})

## Merge inv-17q21.31 vcfs
system("bcftools merge --force-samples scripts/Natalia_ISGlobal/tempVCFs/*inv17.vcf.bgz -m none -o scripts/Natalia_ISGlobal/tempVCFs/merged.inv17.vcf.gz -O z")
vcf <- readVcf("scripts/Natalia_ISGlobal/tempVCFs/merged.inv17.vcf.gz", genome = "hg19")

check <- checkSNPs(vcf)
vcf.filt <- check$genos
inv17sc <- scoreInvHap(SNPlist = vcf.filt, inv = "inv17_007")
save(inv17sc, file = "scripts/Natalia_ISGlobal/inv17Genos.scoreInvHap.Rdata")


####################################
## Differential analysis in cases ##
####################################

cases <- gset[, gset$Status != "Control"]

## inv-8p23.1
load("scripts/Natalia_ISGlobal/inv8Genos.scoreInvHap.Rdata")

inv8class <- classification(inv8sc)
cases$inv8 <- factor(inv8class[colnames(cases)], levels = c("NN", "NI", "II"))
cases$inv8add <- as.numeric(cases$inv8) - 1
modelcases <- model.matrix(~ inv8add + Sex, colData(cases))
lmfitcas <- lmFit(getBeta(cases[, rownames(modelcases)]), design = modelcases)
lmFitecas <- eBayes(lmfitcas)
tabcas <- topTable(lmFitecas, n = Inf, coef = 2)
tabcas <- cbind(tabcas, data.frame(rowRanges(cases)[rownames(tabcas)]))
tabcas8 <- subset(tabcas, seqnames == "chr8")
save(tabcas8, file = "scripts/Natalia_ISGlobal/inv8_cpgAssoc.Rdata")

## inv-16p11.2
load("scripts/Natalia_ISGlobal/inv16Genos.scoreInvHap.Rdata")

inv16class <- classification(inv16sc)
cases$inv16 <- factor(inv16class[colnames(cases)], levels = c("NN", "NI", "II"))
cases$inv16add <- as.numeric(cases$inv16) - 1
modelcases <- model.matrix(~ inv16add + Sex, colData(cases))
lmfitcas <- lmFit(getBeta(cases[, rownames(modelcases)]), design = modelcases)
lmFitecas <- eBayes(lmfitcas)
tabcas <- topTable(lmFitecas, n = Inf, coef = 2)
tabcas <- cbind(tabcas, data.frame(rowRanges(cases)[rownames(tabcas)]))
tabcas16 <- subset(tabcas, seqnames == "chr16")
save(tabcas16, file = "scripts/Natalia_ISGlobal/inv16_cpgAssoc.Rdata")



## inv-17q21.31
load("scripts/Natalia_ISGlobal/inv17Genos.scoreInvHap.Rdata")

inv17class <- classification(inv17sc)
cases$inv17 <- factor(inv17class[colnames(cases)], levels = c("NN", "NI", "II"))
cases$inv17add <- as.numeric(cases$inv17) - 1
modelcases <- model.matrix(~ inv17add + Sex, colData(cases))
lmfitcas <- lmFit(getBeta(cases[, rownames(modelcases)]), design = modelcases)
lmFitecas <- eBayes(lmfitcas)
tabcas <- topTable(lmFitecas, n = Inf, coef = 2)
tabcas <- cbind(tabcas, data.frame(rowRanges(cases)[rownames(tabcas)]))
tabcas17 <- subset(tabcas, seqnames == "chr17")
save(tabcas17, file = "scripts/Natalia_ISGlobal/inv17_cpgAssoc.Rdata")


#################################################
## Select only CpG sites from inversion region ##
#################################################

#Define the start and the end of the inversion region
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")
coordinates <- as.data.frame(inversionGR[c("inv8_001","inv16_009","inv17_007"),])
rownames(coordinates) <- coordinates$Cytogenetic.location

#Define the start and the end of the inversion region +/- 1Mb
inv_regions <- as.data.frame(inversionGR[c("inv8_001","inv16_009","inv17_007"),])
inv_regions$start <- inv_regions$start-1000000
inv_regions$end <- inv_regions$end+1000000
rownames(inv_regions) <- inv_regions$Cytogenetic.location

#Dataset with CpG sites that overlap with SNPs
snps_cpgs <- read.table("/home/isglobal.lan/ncarreras/homews/TFM/Data_Tables/450Khg19SNP.tsv", header=TRUE, sep = "\t")
cpgs_unmask <- snps_cpgs[which(!snps_cpgs$MASK_general),]$probeID

#Selection of CpG sites
select_cpg <- function(tabcas, inv){
  tabcas$CpG <- rownames(tabcas)
  tabcas$Inversion <- inv
  
  inv_methy <- tabcas[which(tabcas$start>=inv_regions[inv,"start"] &
                               tabcas$start<=inv_regions[inv,"end"]),]
  #Select CpG sites withuot SNPs
  inv_methy <- inv_methy[rownames(inv_methy) %in% cpgs_unmask,]
  
  #Adjust by Bonferroni
  inv_methy$adj.P.Val <- p.adjust(inv_methy$P.Value, method="bonferroni")
  
  #Sort by adj.P.Val column
  inv_methy <- inv_methy[order(inv_methy$adj.P.Val),]
}

inv_methy8 <- select_cpg(tabcas8,"8p23.1")
dim(inv_methy8) #898 CpG sites

inv_methy16 <- select_cpg(tabcas16,"16p11.2")
dim(inv_methy16) #409 CpG sites

inv_methy17 <- select_cpg(tabcas17,"17q21.31")
dim(inv_methy17) #698 CpG sites

inv_methy <- rbind(inv_methy8,inv_methy16,inv_methy17)
save(inv_methy,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/heart_Carlos/inv_methy_heart.Rdata")

inv_methy_sig <- inv_methy[which(inv_methy$P.Value<0.005),]

inv_methy_sig <- inv_methy_sig[order(inv_methy_sig$P.Value),]
writexl::write_xlsx(inv_methy_sig,
                    "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/heart_Carlos/inv_methy_heart.xlsx")

####################
## Manhattan plot ##
####################

manhattan <- function(dataset,inv, analysis){
  
  chr <- str_extract(inv,"[\\d]+")
  
  title_plot <- paste0("Differentially methylated CpG sites for ",inv)
  
  #Define the plot title and path
  path_plot_fam <- paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/",
                          analysis, "/Summary_plots/Manhattan_inv_methy",
                          substr(inv,1,2),
                          "_fam.png")
  path_plot_ggviz <-paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/",
                           analysis, "/Summary_plots/Manhattan_inv_methy",
                           substr(inv,1,2),
                           "_ggviz.png")
  
  ############
  ## ggplot ##
  ############
  
  e <- ggplot(dataset, aes(x = start, y = -log10(P.Value))) + 
    geom_point(col="grey", size = 1) +
    geom_point(data=dataset[which(dataset$P.Value<(0.05/nrow(dataset))),],aes(x = start, 
                                                                         y = -log10(P.Value),), 
               col="#2B6A6C", size=1.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    labs(x="Chromosome Position", y="-log10(p-value)", title=title_plot) +
    geom_hline(yintercept = -log10(0.05/nrow(dataset)), color= "#B80D48", linetype="dashed")
  
  png(path_plot_fam, res = 1200, width = 6, height = 5, units = "in")
  par(mar=c(5,5,2,2))
  
  print(e)
  
  dev.off()
  
  ##########
  ## Gviz ##
  ##########
  
  #Define the tracks
  
  #Track ideogram (chromosome)
  itrack <- IdeogramTrack(genome = "hg19", chromosome = paste0("chr",chr), 
                          cex = 1.3, fontsize = 14)
  
  #Track for the chr positions
  gtrack <- GenomeAxisTrack(cex = 1.3, fontsize = 11)
  
  #Track manhattan plot
  dtrack <- DataTrack(data = -log10(dataset$P.Value), start = dataset$start,
                      end = dataset$end, 
                      chromosome = paste0("chr",chr), genome = "hg19", 
                      name = "-log10(p-value)",
                      size = 13, 
                      background.title = "#525252",
                      fontsize=24)
  
  #Track inversion region
  invtrack <- GeneRegionTrack(coordinates[which(coordinates$seqnames==paste0("chr",chr)),], genome = "hg19",
                              chromosome = paste0("chr",chr), 
                              name = "Inv",
                              background.title = "#F27F24", fill="#F29724", col="#F29724")
  
  png(path_plot_ggviz, res = 1200, width = 6, height = 6.3, units = "in")
  par(mar=c(5,5,2,2))
  plotTracks(list(itrack,gtrack,dtrack,invtrack), 
             cex.title = 1, sizes=c(1.4,1.5,7,0.7),
             title.width = 1.8)
  dev.off()
}

manhattan(inv_methy8, "8p23.1", analysis = "heart_Carlos")
manhattan(inv_methy16, "16p11.2", analysis = "heart_Carlos")
manhattan(inv_methy17, "17q21.31", analysis = "heart_Carlos")
