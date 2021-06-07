########################################
## Principal Component Analysis (PCA) ##
########################################

#Load libraries

library(ggpubr)
library(SummarizedExperiment)

#Load datasets

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy17.Rdata")

#PCA analysis and PCA plot (scores)

pca_methy <- function(methy,inv, name_inv){
  cpgs_inv <- t(assay(methy))
  invgenotypes_inv <- methy[[inv]]
  pc <- princomp(cpgs_inv)
  data <- cbind(as.data.frame(pc$scores[,1]),
                as.data.frame(pc$scores[,2]),
                invgenotypes_inv)
  PoV <- pc$sdev^2/sum(pc$sdev^2)*100
  colnames(data) <- c("PCA1","PCA2","Genotype")
  png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/PCA_",inv,"methy.png"))
  ggscatterhist(
    data, x = "PCA1", y = "PCA2",
    xlab = paste0("PCA1 (",round(PoV[1],1),"%)"),
    ylab = paste0("PCA2 (",round(PoV[2],1),"%)"),
    color = "Genotype", size = 3, alpha = 0.6,
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    margin.params = list(fill = "Genotype", color = "black", size = 0.2)
  )
  dev.off()
  return(pc)
}

pc_methy8 <- pca_methy(methy8,"inv8_001", "8p23.1")
pc_methy16 <- pca_methy(methy16,"inv16_009", "16p11.2")
pc_methy17 <- pca_methy(methy17,"inv17_007", "17q21.31")

#Loadings plot
loading_plot <- function(pc,name){
  PoV <- pc$sdev^2/sum(pc$sdev^2)*100
  data <- data.frame(Load1=pc$loadings[,1],
                     Load2=pc$loadings[,2])
  png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/Loading_",name,".png"))
  ggscatterhist(
    data, x = "Load1", y = "Load2",
    xlab = paste0("PCA1 (",round(PoV[1],1),"%)"),
    ylab = paste0("PCA2 (",round(PoV[2],1),"%)"),
    size = 3, alpha = 0.6
  )
  dev.off()
}

loading_plot(pc_methy8,"inv8_001_methy")
loading_plot(pc_methy16,"inv16_009_methy")
loading_plot(pc_methy17,"inv17_007_methy")

#Biplot
png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/biplot8.png"))
biplot(pc_methy8)
dev.off()

png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/biplot16.png"))
biplot(pc_methy16)
dev.off()

png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/biplot17.png"))
biplot(pc_methy17)
dev.off()

#Correlation PCA1 ~ genotype
corr <- function(pc,inv,dataset){
  genotype <- as.numeric(dataset[[inv]])-1
  pca1 <- pc$scores[,1]
  res <- summary(lm(pca1 ~ genotype))
  return(res)
}

res8 <- corr(pc_methy8,"inv8_001",methy8)
res16 <- corr(pc_methy16,"inv16_009",methy16)
res17 <- corr(pc_methy17,"inv17_007",methy17)