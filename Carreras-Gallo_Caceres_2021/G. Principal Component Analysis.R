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
  colnames(data) <- c("PC1","PC2","Genotype")
  png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/PCA_",inv,"methy.png"),
      res = 1200,
      width = 3.25, height = 3.25, units = "in", pointsize = 3)
  par(mar=c(5,5,2,2))
  ggscatterhist(
    data, x = "PC1", y = "PC2",
    xlab = paste0("PC1 (",round(PoV[1],1),"%)"),
    ylab = paste0("PC2 (",round(PoV[2],1),"%)"),
    color = "Genotype", size = 1, alpha = 0.6,
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    margin.params = list(fill = "Genotype", color = "black", size = 0.2),
    font.x = c(12, "bold"),
    font.y = c(12, "bold"),
    font.legend = c(10),
    font.xtickslab = 8,
    font.ytickslab = 8,
    ggtheme = theme(axis.title.y=element_text(margin=margin(r=10)),
                    axis.title.x=element_text(margin=margin(t=10)),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour="black", size = 0.3),
                    legend.key = element_rect(fill=NA))
  )
  dev.off()
  return(pc)
}

pc_methy8 <- pca_methy(methy8,"inv8_001", "8p23.1")
pc_methy16 <- pca_methy(methy16,"inv16_009", "16p11.2")
pc_methy17 <- pca_methy(methy17,"inv17_007", "17q21.31")

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