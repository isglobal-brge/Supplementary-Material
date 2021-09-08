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

analysis="metanalysis"

#Load results from the transcriptome and select the significant ones

load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/",analysis,"/inv_trans.Rdata"))
inv_trans_sig <- inv_trans[which(inv_trans$p.adj<0.05),]

topmean_trans_8 <- split_symbol(inv_trans_sig[which(inv_trans_sig$Inversion=="8p23.1"),"Gene_Symbol"])
topmean_trans_16 <- split_symbol(inv_trans_sig[which(inv_trans_sig$Inversion=="16p11.2"),"Gene_Symbol"])
topmean_trans_17 <- split_symbol(inv_trans_sig[which(inv_trans_sig$Inversion=="17q21.31"),"Gene_Symbol"])

#Load results from the methylome and select the significant ones

load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/",analysis,"/inv_methy.Rdata"))
inv_methy_sig <- inv_methy[which(inv_methy$p.adj<0.05),]

topmean_methy_8 <- split_symbol(inv_methy_sig[which(inv_methy_sig$Inversion=="8p23.1"),"Gene_Symbol"])
topmean_methy_16 <- split_symbol(inv_methy_sig[which(inv_methy_sig$Inversion=="16p11.2"),"Gene_Symbol"])
topmean_methy_17 <- split_symbol(inv_methy_sig[which(inv_methy_sig$Inversion=="17q21.31"),"Gene_Symbol"])

#Load results from the methylome and select the significant ones
load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/",analysis,"/it_sig_expr.Rdata"))

allsymbols_to_symbol <- function(allsymbols, join=T) {
  allsymbols <- paste(allsymbols,collapse=";")
  symbols<-strsplit(allsymbols,";")[[1]]
  symbols[symbols=="NA" | symbols==""]=NA
  symbols<-unique(symbols[!is.na(symbols)])
  if (join){
    return (paste(symbols, collapse=";"))
  }
  return(symbols)
}

cys_8 <- allsymbols_to_symbol(it_sig[which(it_sig$Inversion=="8p23.1"),]$cis.eQTM, join=F)
cys_16 <- allsymbols_to_symbol(it_sig[which(it_sig$Inversion=="16p11.2"),]$cis.eQTM, join=F)
cys_17 <- allsymbols_to_symbol(it_sig[which(it_sig$Inversion=="17q21.31"),]$cis.eQTM, join=F)

#See the intersect gene symbols

intersect(topmean_trans_8,topmean_methy_8)
intersect(topmean_trans_16,topmean_methy_16)
intersect(topmean_trans_17,topmean_methy_17)

######################
## eQTL vs cis-eQTM ##
######################

venn.diagram(
  x = list(topmean_trans_8, cys_8),
  category.names = c("eQTL", "cis-eQTM"),
  filename = paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/",analysis,"/venn_8_expr.png"),
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
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0.05,0.05),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "serif"
)

venn.diagram(
  x = list(topmean_trans_16, cys_16),
  category.names = c("eQTL", "cis-eQTM"),
  filename = paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/",analysis,"/venn_16_expr.png"),
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
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0.05,0.05),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "serif"
)


venn.diagram(
  x = list(topmean_trans_17, cys_17),
  category.names = c("eQTL", "cis-eQTM"),
  filename = paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/",analysis,"/venn_17_expr.png"),
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
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0.05,0.05),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "serif"
)



##################
## eQTL vs mQTL ##
##################

#Inversion 8p23.1
venn.diagram(
  x = list(topmean_trans_8, topmean_methy_8),
  category.names = c("eQTL", "mQTL"),
  filename = paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/",analysis,"/venn_8_transmethy_mean.png"),
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
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(-0.4, -0.45),
  cat.fontfamily = "serif"
)

#Inversion 16p11.2
venn.diagram(
  x = list( topmean_trans_16, topmean_methy_16),
  category.names = c("eQTL","mQTL"),
  filename = paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/",analysis,"/venn_16_transmethy_mean.png"),
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
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(0.05, 0.05),
  cat.fontfamily = "serif"
)

#Inversion 17q21.31
venn.diagram(
  x = list(topmean_trans_17, topmean_methy_17),
  category.names = c("eQTL", "mQTL"),
  filename = paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Venn_trans_vs_methy/",analysis,"/venn_17_transmethy_mean.png"),
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
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(0.1, 0.05),
  cat.fontfamily = "serif"
)
