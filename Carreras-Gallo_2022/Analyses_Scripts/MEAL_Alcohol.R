##############################
## ALCOHOL UNIFORM ANALYSES ##
##############################

#Select model
model <- "5levels"

#Set working directory
setwd(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Uniform_analysis/Alcohol/",model))

#Load libraries
library(minfi)
library(MEAL)
library(ggplot2)
library(plyr)
library(qqman)

#Set dataset
dataset <- "SVA"
load(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_",dataset,".Rdata"))

#Set covariates
covs <- c("Biological_Sex","age", "Ethnicity",
          "BMI","Level_of_Education","Slide",
          "Bcell","CD4T","CD8T","Eos","Mono","Neu","NK",
          "Tobacco_Use")

#Set variable

##### Model A: Never (0) vs Ever (1)
if (model=="Never_vs_Ever"){
  GRset$Alcohol_ever <- factor(ifelse(GRset$Alcohol_Use_per_week=="Never","Never","Ever"),
                             levels=c("Never","Ever"))
  var <- "Alcohol_ever"
  coef <- 2
}

##### Model B: Never (0) vs Minimum 3 times per week
if (model=="Never_vs_min3week"){
  GRset$Alcohol_min3 <- revalue(GRset$Alcohol_Use_per_week, c("Never"="Never", 
                                                 "3-5 times per week"="Min 3 times per week",
                                                 "Regularly"="Min 3 times per week"))
  GRset$Alcohol_min3[GRset$Alcohol_min3!="Never" & GRset$Alcohol_min3!="Min 3 times per week"] <- NA
  GRset$Alcohol_min3 <- droplevels(GRset$Alcohol_min3)
  var <- "Alcohol_min3"
  coef <- 2
}

#### Model C: Factor
if (model=="5levels"){
  var <- "Alcohol_Use_per_week"
  coef <- 2:5
}

table(GRset[[var]])

gc(reset=TRUE)

#Perform differential analysis
resmean <- runPipeline(set = GRset, 
                       variable_names = var,
                       covariable_names	= covs,
                       analyses = "DiffMean",
                       verbose = TRUE,
                       method="robust")

#Save result from the differential analysis
save(resmean,file=paste0("resmean_alcohol_",model,".Rdata"))

#Get topcpgs
topcpgs <- getAssociation(resmean, coef=coef, 
                          fNames=c("chr","pos","HGNC_GeneSymbol","UCSC_RefGene_Group"), 
                          rid="DiffMean")

topcpgs$CpG <- rownames(topcpgs)

#Save results
save(topcpgs,file=paste0("topcpgs_alcohol_",model,".Rdata"))
writexl::write_xlsx(topcpgs[which(topcpgs$adj.P.Val<0.05),],paste0("topcpgs_alcohol_",model,".xlsx"))


############
## QQPLOT ##
############

png("qqplot.png", 
    res = 1200, width = 10, height = 10, units = "in")
par(mar=c(5,5,2,2))

plot(resmean, rid = "DiffMean", type = "qq") + ggtitle("QQ-plot for Alcohol effect")

dev.off()


####################
## MANHATTAN PLOT ##
####################

#Choose line of significance
last_sig <- rev(topcpgs[which(topcpgs$adj.P.Val<0.05),]$P.Value)[1]
first_nosig <- topcpgs[which(topcpgs$adj.P.Val>0.05),]$P.Value[1]

line_manha <- -log10(mean(last_sig,first_nosig))

if (is.na(line_manha)){
  line_manha <- -log10(first_nosig-first_nosig/2)
}
#Chromosome column as numeric
topcpgs$chr_num <- sapply(topcpgs$chr, function(x) strsplit(x,"chr")[[1]][2])
topcpgs$chr_num[topcpgs$chr_num=="X"] <- "23"
topcpgs$chr_num[topcpgs$chr_num=="Y"] <- "24"
topcpgs$chr_num <- as.numeric(topcpgs$chr_num)

#Select only 1 gene per CpG and NA for intergenics
topcpgs$Gene_1 <- NA

for (i in 1:15){
  if (topcpgs$HGNC_GeneSymbol[i]==""){
    topcpgs$Gene_1[i] <- NA
  }
  else{
    topcpgs$Gene_1[i] <- strsplit(topcpgs$HGNC_GeneSymbol[i],";")[[1]][1]
  }
}

#Create Manhattan plot
png("manhattan.png", 
    res = 900, width = 8, height = 4, units = "in")
par(mar=c(5,5,2,2))

manhattan(topcpgs, chr = "chr_num", bp = "pos", p = "P.Value",  snp = "Gene_1",
          suggestiveline = -log10(0.0001),
          genomewideline = line_manha,
          annotatePval = 9e-18,
          col = c("#97c3d7", "#3478a9","#a4d38f","#54a353","#f19192","#d72d2d","#f4b275"), 
          main = "Alcohol EWAS (N=3424)",
          chrlabs = c(1:22, "X", "Y"),
          cex.axis=0.6)

dev.off()
