##########################################################################################
######### VENN DIAGRAMS TO COMPARE GENES DIFFERENTIALLY EXPRESSED AND METHYLATED #########
##########################################################################################

#Load libraries

library(ggplot2)
library(VennDiagram)

#Create a funtion to separate the gene symbols

split_symbol <- function(top) {
  a <- paste(top,collapse=";")
  b <- unique(strsplit(a,";")[[1]])
  c <- b[which(b!="")]
  return (c)
}

#Load results from the transcriptome and select the significant ones

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/inv_trans.Rdata")
inv_trans_sig <- inv_trans_mean[which(inv_trans_mean$p.adj<0.05),]

topmean_trans_8 <- split_symbol(inv_trans_sig[which(inv_trans_sig$Inversion=="8p23.1"),"Gene_Symbol"])
topmean_trans_16 <- split_symbol(inv_trans_sig[which(inv_trans_sig$Inversion=="16p11.2"),"Gene_Symbol"])
topmean_trans_17 <- split_symbol(inv_trans_sig[which(inv_trans_sig$Inversion=="17q21.31"),"Gene_Symbol"])

#Load results from the methylome and select the significant ones

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/inv_methy.Rdata")
inv_methy_sig <- inv_methy_mean[which(inv_methy_mean$p.adj<0.05),]

topmean_methy_8 <- split_symbol(inv_methy_sig[which(inv_methy_sig$Inversion=="8p23.1"),"Gene_Symbol"])
topmean_methy_16 <- split_symbol(inv_methy_sig[which(inv_methy_sig$Inversion=="16p11.2"),"Gene_Symbol"])
topmean_methy_17 <- split_symbol(inv_methy_sig[which(inv_methy_sig$Inversion=="17q21.31"),"Gene_Symbol"])

#See the intersent gene symbols

intersect(topmean_trans_8,topmean_methy_8)
intersect(topmean_trans_16,topmean_methy_16)
intersect(topmean_trans_17,topmean_methy_17)

#Create the Venn diagrams

#Inversion 8p23.1
venn.diagram(
  x = list(topmean_trans_8, topmean_methy_8),
  category.names = c("eQTL", "mQTL"),
  filename = "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/venn_8_transmethy_mean.png",
  output=TRUE,
  
  # Circles
  lwd = 1.5,
  lty = 'blank',
  fill = c("#B3E2CD","#FDCDAC"),
  
  # Numbers
  cex = 3,
  fontface = "bold",
  fontfamily = "serif",
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(-0.4, -0.45),
  cat.fontfamily = "serif"
)

#Inversion 16p11.2
venn.diagram(
  x = list(topmean_trans_16, topmean_methy_16),
  category.names = c("eQTL", "mQTL"),
  filename = "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/venn_16_transmethy_mean.png",
  output=TRUE,
  
  # Circles
  lwd = 1.5,
  lty = 'blank',
  fill = c("#B3E2CD","#FDCDAC"),
  
  # Numbers
  cex = 3,
  fontface = "bold",
  fontfamily = "serif",
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(-0.45, -0.45),
  cat.fontfamily = "serif"
)

#Inversion 17q21.31
venn.diagram(
  x = list(topmean_trans_17, topmean_methy_17),
  category.names = c("eQTL", "mQTL"),
  filename = "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/venn_17_transmethy_mean.png",
  output=TRUE,
  
  # Circles
  lwd = 1.5,
  lty = 'blank',
  fill = c("#B3E2CD","#FDCDAC"),
  
  # Numbers
  cex = 3,
  fontface = "bold",
  fontfamily = "serif",
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(0.1, 0.05),
  cat.fontfamily = "serif"
)
