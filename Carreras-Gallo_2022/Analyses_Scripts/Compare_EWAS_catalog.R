###############################
## COMPARE WITH EWAS CATALOG ##
###############################

library(VennDiagram)

setwd("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Uniform_analysis/")

#Load EWAS catalogues
cat <- c("alcohol_consumption","lifetime_cannabis_use","smoking")
catalogues <- list()
for (i in 1:length(cat)){
  catalogues[[i]] <- read.table(file=paste0("EWAS_catalog/",cat[i],".tsv"),sep = '\t',header=T)
}

#Load results from TruDiagnostic
truD <- list()
subst <- c("Alcohol","Marijuana","Tobacco")

for (i in 1:length(subst)){
  if (subst[i]=="Tobacco"){
    load(list.files(path=paste0(subst[i],"/7levels"),pattern="topcpgs.+Rdata",full.names=TRUE))
    truD[[i]] <- topcpgs
  }
  else{
    load(list.files(path=paste0(subst[i],"/5levels"),pattern="topcpgs.+Rdata",full.names=TRUE))
    truD[[i]] <- topcpgs
  }
}

#Change names from both lists
names(truD) <- subst
names(catalogues) <- subst

#Create a data.frame with the number of publications and the most significant p-value for each CpG in the catalog
publ_sig_list <- list()

for (i in 1:length(catalogues)){
  
  #Count how many times each CpG appears as significant
  publ_sig <- data.frame(table(catalogues[[i]]$cpg))
  names(publ_sig)[1] <- "CpG"
  
  #Include the most significant P.Value for each CpG
  publ_sig$Min.P.Val <- apply(publ_sig,1,function(x){
    min(catalogues[[i]][which(catalogues[[i]]$cpg==x[1]),"p"])}
  )

  #Include the sum of the beta value signs for each CpG
  publ_sig$Sum.Beta <- apply(publ_sig,1,function(x){
    if (all(is.na(catalogues[[i]][which(catalogues[[i]]$cpg==x[1]),"beta"]))){
      NA
    }
    else{
      sum(sign(catalogues[[i]][which(catalogues[[i]]$cpg==x[1]),"beta"]),na.rm=TRUE)}
    }
  )
  
  #Include the P.Value from TruDiagnostic EWAS
  publ_sig$P.Val.TruD <- apply(publ_sig,1,function(x){
    if(x[1]%in%truD[[i]]$CpG){
      truD[[i]][which(truD[[i]]$CpG==x[1]),"P.Value"]
    }
    else{
      NA
    }
  }
  )
  
  #Include the Regularly Beta value from TruDiagnostic EWAS
  publ_sig$Beta.TruD.Regularly <- apply(publ_sig,1,function(x){
    if(x[1]%in%truD[[i]]$CpG){
      truD[[i]][which(truD[[i]]$CpG==x[1]),4]
    }
    else{
      NA
    }
  }
  )
  
  #Include whether the P.Value is significant or not according to the P.Val.Adj
  publ_sig$Sig.TruD <- apply(publ_sig,1,function(x){
    if(x[1]%in%truD[[i]]$CpG){
      truD[[i]][which(truD[[i]]$CpG==x[1]),"adj.P.Val"]<0.05
    }
    else{
      NA
    }
  }
  )
  
  #Include whether the beta sign is equal or not
  publ_sig$Equal.Sign <- sign(publ_sig$Sum.Beta)==sign(publ_sig$Beta.TruD.Regularly)
  
  #Order the results by P.Value and afterwards by Freq
  publ_sig <- publ_sig[order(publ_sig$Min.P.Val),]
  publ_sig <- publ_sig[order(publ_sig$Freq, decreasing = TRUE),]
  
  #Save the result in the list
  publ_sig_list[[i]] <- publ_sig
}

names(publ_sig_list) <- subst

#Save results
save(publ_sig_list, file="EWAS_catalog/publ_sig_list.Rdata")
writexl::write_xlsx(publ_sig_list[[1]][1:200,],path="EWAS_catalog/catalog_alcohol_compare.xlsx")
writexl::write_xlsx(publ_sig_list[[2]],path="EWAS_catalog/catalog_marijuana_compare.xlsx")
writexl::write_xlsx(publ_sig_list[[3]][1:200,],path="EWAS_catalog/catalog_tobacco_compare.xlsx")

#Venn diagrams

venn <- function(substance){
  #Create a data frame with the publ_sig for this substance
  df <- publ_sig_list[[substance]]
  
  #Remove the CpG sites that are not present in TruDiagnostic
  df_clean <- df[which(!is.na(df$P.Val.TruD)),]
  
  #Create a vector with all the CpG sites that are significant in at least half of the studies
  min_freq <- df[1,"Freq"]/3
  if (min_freq>5){
    min_freq <- 5
  }
  cpg_sig <- df_clean[which(df_clean$Freq>=min_freq),]$CpG
  
  
  #Select the significant CpG sites in TruDiagnostic
  df_truD <- truD[[substance]]
  df_sig <- df_truD[which(df_truD$P.Value<0.0001),]
  cpg_tru <- df_sig$CpG
  
  #Create Venn diagram
  venn.diagram(
    x = list(cpg_sig, cpg_tru),
    category.names = c("EWAS catalog","TruDiagnostic"),
    filename = paste0("EWAS_catalog/venn_diagram_EWAScatalog_",substance,".png"),
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 2000 , 
    width = 2000 , 
    resolution = 600,
    compression = "lzw",
    
    # Circles
    lwd = 1.5,
    lty = 'blank',
    fill = c("#B3E2CD","#FDCDAC"),
    
    # Numbers
    cex = 0.9,
    fontface = "plain",
    fontfamily = "serif",
    
    # Set names
    cat.cex = 0.9,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(0,0),
    cat.dist = c(0.05,0.05),
    cat.fontfamily = "serif"
  )
  
}

venn(substance="Alcohol")
venn(substance="Marijuana")
venn(substance="Tobacco")
