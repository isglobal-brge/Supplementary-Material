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


###########################################
## expression genes ~ PCA patterns methy ##
###########################################

#Load datasets

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans17.Rdata")

#Select the IDs shared between gene expression and DNA methylation
IDs <- intersect(colnames(methy8),colnames(trans8))

#Discard the IDs that are not present in the DNA methylation dataset
trans8 <- trans8[,IDs]
trans16 <- trans16[,IDs]
trans17 <- trans17[,IDs]

#Discard the IDs from the PCA calculated
scores_methy8 <- pc_methy8$scores[IDs,1:2]
scores_methy16 <- pc_methy16$scores[IDs,1:2]
scores_methy17 <- pc_methy17$scores[IDs,1:2]

sex <- trans8$sex
age <- trans8$age

#Study the correlation between the expression of the genes and the PC1 and PC2
corr_PC_expr <- function(inv, inv_cov){
  if (inv=="8p23.1"){
    trans <- trans8
    scores <- scores_methy8
  }
  if (inv=="16p11.2"){
    trans <- trans16
    scores <- scores_methy16
  }
  if (inv=="17q21.31"){
    trans <- trans17
    scores <- scores_methy17
  }
  
  df <- data.frame(Transcript=as.character(),
                   Gene_Symbol=as.character(),
                   Inversion=as.character,
                   r2_PCA1=as.numeric(),
                   pval_PCA1=as.numeric(),
                   r2_PCA2=as.numeric(),
                   pval_PCA2=as.numeric(),
                   stringsAsFactors = F)
  
  for (i in 1:nrow(trans)){
    pca1 <- summary(lm(assays(trans)[["exprs"]][i,] ~ scores[,1] + sex + age + inv_cov))
    r2_pca1 <- pca1$adj.r.squared
    pval_pca1 <- pca1$coefficients[2,4]
    
    pca2 <- summary(lm(assays(trans)[["exprs"]][i,] ~ scores[,2] + sex + age + inv_cov))
    r2_pca2 <- pca2$adj.r.squared
    pval_pca2 <- pca2$coefficients[2,4]
    
    df_i <- data.frame(Transcript=rownames(trans)[i],
                       Gene_Symbol=rowData(trans)$GeneSymbol_Affy[i],
                       Inversion=inv,
                       r2_PCA1=r2_pca1,
                       pval_PCA1=pval_pca1,
                       r2_PCA2=r2_pca2,
                       pval_PCA2=pval_pca2)
    df <- rbind(df,df_i)
  }
  
  return(df)
}

df_corr8 <- corr_PC_expr("8p23.1", trans8$inv8_001)
df_corr16 <- corr_PC_expr("16p11.2", trans8$inv16_009)
df_corr17 <- corr_PC_expr("17q21.31", trans8$inv17_007)

#Sort the results by pval
df_corr8 <- df_corr8[order(df_corr8$pval_PCA1),]
df_corr16 <- df_corr16[order(df_corr16$pval_PCA1),]
df_corr17 <- df_corr17[order(df_corr17$pval_PCA1),]

#Look at the significant correlations
df_corr8[which(df_corr8$pval_PCA1<0.05/nrow(df_corr8)),]
df_corr16[which(df_corr16$pval_PCA1<0.05/nrow(df_corr16)),]
df_corr17[which(df_corr17$pval_PCA1<0.05/nrow(df_corr17)),]



#Save the results

writexl::write_xlsx(df_corr8, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/correlation_expression/df_corr8_adj.xlsx")
writexl::write_xlsx(df_corr16, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/correlation_expression/df_corr16_adj.xlsx")
writexl::write_xlsx(df_corr17, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/PCA/correlation_expression/df_corr17_adj.xlsx")

