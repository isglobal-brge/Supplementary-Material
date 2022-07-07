##############################
## TOBACCO UNIFORM ANALYSES ##
##############################

#Select model
model <- "7levels"

#Set working directory
setwd(paste0("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Uniform_analysis/Tobacco/",model))

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
covs <- c("Biological_Sex","age","Ethnicity",
          "BMI","Level_of_Education","Slide",
          "Bcell","CD4T","CD8T","Eos","Mono","Neu","NK",
          "Alcohol_Use_per_week")

#Set variable

##### Model A: None (0) vs Any (1)
if(model=="None_vs_Any"){
  GRset$Tobacco_any <- factor(ifelse(GRset$Tobacco_Use=="None","None","Any"),
                              levels=c("None","Any"))
  var <- "Tobacco_any"
  coef <- 2
}

##### Model B: None (0) vs Minimum 1 per day
if(model=="None_vs_min1day"){
  GRset$Tobacco_min1 <- revalue(GRset$Tobacco_Use, c("None"="None", 
                                                 "1-5 cigarettes per day"="Min 1 per day",
                                                 "6-10 cigarettes per day"="Min 1 per day",
                                                 "11-20 cigarettes per day"="Min 1 per day",
                                                 "More than 20 cigarettes per day"="Min 1 per day"))
  GRset$Tobacco_min1[GRset$Tobacco_min1!="None" & GRset$Tobacco_min1!="Min 1 per day"] <- NA
  GRset$Tobacco_min1 <- droplevels(GRset$Tobacco_min1)
  var <- "Tobacco_min1"
  coef <- 2
}

#### Model C: Factor
if(model=="7levels"){
  var <- "Tobacco_Use"
  coef <- 2:7
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
save(resmean,file=paste0("resmean_tobacco_",model,".Rdata"))

#Get topcpgs
topcpgs <- getAssociation(resmean, coef=coef, 
                          fNames=c("chr","pos","HGNC_GeneSymbol","UCSC_RefGene_Group"), 
                          rid="DiffMean")

topcpgs$CpG <- rownames(topcpgs)

#Save results
save(topcpgs,file=paste0("topcpgs_tobacco_",model,".Rdata"))
writexl::write_xlsx(topcpgs[which(topcpgs$adj.P.Val<0.05),],paste0("topcpgs_tobacco_",model,".xlsx"))


############
## QQPLOT ##
############

png("qqplot.png", 
    res = 1200, width = 10, height = 10, units = "in")
par(mar=c(5,5,2,2))

plot(resmean, rid = "DiffMean", type = "qq") + ggtitle("QQ-plot for Tobacco effect")

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

#Create Manhattan plot
png("manhattan.png", 
    res = 1200, width = 10, height = 6, units = "in")
par(mar=c(5,5,2,2))

manhattan(topcpgs, chr = "chr_num", bp = "pos", p = "P.Value",  snp = "CpG",
          suggestiveline = -log10(0.001),
          genomewideline = line_manha,
          col = c("#97c3d7", "#3478a9","#a4d38f","#54a353","#f19192","#d72d2d","#f4b275"), 
          main = "EWAS of Tobacco Consumption (N=3424)",
          chrlabs = c(1:22, "X", "Y"))

dev.off()
