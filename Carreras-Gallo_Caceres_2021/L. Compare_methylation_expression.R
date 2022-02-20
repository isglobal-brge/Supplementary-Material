###########################################
## ADD GENE GROUP TO METHYLATION RESULTS ##
###########################################

#Load libraries
library("readxl")
library(VennDiagram)

#Load methylation data
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy17.Rdata")

#Load table with the effects of inversion on CpG methylation
inv_methy <- data.frame(read_excel("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy.xlsx",
                      sheet = 1))

inv_methy$Gene_group <- NA

for (i in 1:nrow(inv_methy)){
  if (inv_methy$Inversion[i]=="8p23.1"){
    methy <- rowData(methy8)
  }
  if (inv_methy$Inversion[i]=="16p11.2"){
    methy <- rowData(methy16)
  }
  if (inv_methy$Inversion[i]=="17q21.31"){
    methy <- rowData(methy17)
  }
  cpg <- inv_methy$CpG[i]
  inv_methy$Gene_group[i] <- methy[cpg,"UCSC_RefGene_Group"]
}

writexl::write_xlsx(inv_methy, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy_gene_group.xlsx")


#Load table with the effects of inversion*exposure on CpG methylation
inter_methy <- data.frame(read_excel("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy.xlsx"))

inter_methy$Gene_Group <- NA

for (i in 1:nrow(inter_methy)){
  if (inter_methy$Inversion[i]=="8p23.1"){
    methy <- rowData(methy8)
  }
  if (inter_methy$Inversion[i]=="16p11.2"){
    methy <- rowData(methy16)
  }
  if (inter_methy$Inversion[i]=="17q21.31"){
    methy <- rowData(methy17)
  }
  cpg <- inter_methy$CpG[i]
  inter_methy$Gene_Group[i] <- methy[cpg,"UCSC_RefGene_Group"]
}

writexl::write_xlsx(inter_methy, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy_GeneGroup.xlsx")



###################
## VENN DIAGRAMS ##
###################

#Create a funtion to separate the gene symbols
split_symbol <- function(top) {
  a <- paste(top,collapse=";")
  b <- unique(strsplit(a,";")[[1]])
  c <- b[which(b!="")]
  return (c)
}

#Load table with the effects of inversion on gene expression
inv_trans <- data.frame(read_excel("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/metanalysis/inv_trans.xlsx",
                                   sheet = 1))

#Store the gene symbols
topmean_trans_8 <- split_symbol(inv_trans[which(inv_trans$Inversion=="8p23.1"),"Gene_Symbol"])
topmean_trans_16 <- split_symbol(inv_trans[which(inv_trans$Inversion=="16p11.2"),"Gene_Symbol"])
topmean_trans_17 <- split_symbol(inv_trans[which(inv_trans$Inversion=="17q21.31"),"Gene_Symbol"])

topmean_methy_8 <- split_symbol(inv_methy[which(inv_methy$Inversion=="8p23.1"),"Gene_Symbol"])
topmean_methy_16 <- split_symbol(inv_methy[which(inv_methy$Inversion=="16p11.2"),"Gene_Symbol"])
topmean_methy_17 <- split_symbol(inv_methy[which(inv_methy$Inversion=="17q21.31"),"Gene_Symbol"])

#See the intersect gene symbols
intersect(topmean_trans_8,topmean_methy_8)
intersect(topmean_trans_16,topmean_methy_16)
intersect(topmean_trans_17,topmean_methy_17)

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
