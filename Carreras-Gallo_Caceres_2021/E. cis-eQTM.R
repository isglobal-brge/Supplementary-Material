######################################################
## cis-eQTM for differentially methylated CpG sites ##
######################################################

eQTM <- read.table("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/eQTM/eQTM_autosome_adj.cells_SIG.txt", sep="\t", header = T)

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy.Rdata")

it_sig <- inv_methy[which(inv_methy$p.adj<0.05),]

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

#Add the genes they have an effect

it_sig$cis.eQTM <- NA

for (i in 1:nrow(it_sig)){
  associations <- eQTM[which(eQTM$CpG==it_sig$CpG[i]),]
  genes <- allsymbols_to_symbol(associations$TC_gene)
  it_sig$cis.eQTM[i] <- genes
}

save(it_sig,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/it_sig_expr.Rdata")
writexl::write_xlsx(it_sig, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/it_sig_expr.xlsx")


inv8 <- it_sig[which(it_sig$Inversion=="8p23.1"),]
allsymbols_to_symbol(inv8$cis.eQTM, join=F)

inv16 <- it_sig[which(it_sig$Inversion=="16p11.2"),]
allsymbols_to_symbol(inv16$cis.eQTM, join=F)

inv17 <- it_sig[which(it_sig$Inversion=="17q21.31"),]
allsymbols_to_symbol(inv17$cis.eQTM, join=F)

