########################################################
## DIFFERENTIAL ANALYSIS BY INVERSION (meta-analysis) ##
########################################################

#Load libraries
library(SummarizedExperiment)
library(rexposome)
library(MEAL)
library(MetaDE)
library(readxl)
library(meta)

#Load inversion information
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")

#Define a function to perform differential analysis
diff_inv <- function(dataset, inversion, omics, covs){
  print(paste(inversion, "-", omics))
  if (omics=="Transcriptome"){
    mat <- assays(dataset)$expr
  }
  if (omics=="Methylome"){
    mat <- assay(dataset)
  }
  mod <- model.matrix(formula(paste("~", inversion, "+", paste(colnames(covs)[4:ncol(covs)], collapse = "+"))), data=covs)
  mat <- mat[,rownames(mod)]
  resmean <- runDiffMeanAnalysis(set = mat, model = mod, 
                                 resultSet = TRUE, method="robust")
  topfeatures <- getAssociation(resmean, coef=inversion, 
                                fNames=NULL, rid="DiffMean")
  
  return(topfeatures)
}

############################################
## DIFFERENTIAL ANALYSIS IN TRANSCRIPTOME ##
############################################

#Load trancriptome datasets for the three inversions
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans17.Rdata")

#Select adjustment factors data for transcriptome
covs_trans <- data.frame(inv8_001=as.numeric(trans8$inv8_001)-1,
                         inv16_009=as.numeric(trans8$inv16_009)-1,
                         inv17_007=as.numeric(trans8$inv17_007)-1,
                         sex=trans8$sex,
                         age=trans8$age)
covs_trans <- cbind(covs_trans, colData(trans8)[,c(18:23,54:63)])

#Define a function to obtain a data frame with the results of the differential analysis
diff_inv_trans <- function(dataset, inversion, covs_trans){
  #Create an empty data frame
  inv_trans <- data.frame(Inversion=character(),
                            Transcript=character(),
                            Location=character(),
                            Gene_Symbol=character(),
                            LogFC=numeric(),
                            CI.L=numeric(),
                            CI.R=numeric(),
                            AveExpr=numeric(),
                            t=numeric(),
                            pval=numeric(),
                            B=numeric(),
                            SE=numeric(),
                            stringsAsFactors=FALSE)
  
  #Differential analysis for all the features in the inversion region
  topgenes <- diff_inv(dataset, inversion, omics="Transcriptome", covs_trans)
  for (transcript in rownames(topgenes)){
    df <- data.frame(Inversion=inversionGR[inversion,]$Cytogenetic.location,
                     Transcript = transcript,
                     Location=paste(as.character(seqnames(dataset[transcript,])),":",as.character(start(dataset[transcript,])),"-",as.character(end(dataset[transcript,])),sep=""),
                     Gene_Symbol=rowData(dataset)[transcript,"GeneSymbol_Affy"],
                     LogFC=topgenes[transcript,"logFC"],
                     CI.L=topgenes[transcript,"CI.L"],
                     CI.R=topgenes[transcript,"CI.R"],
                     AveExpr=topgenes[transcript,"AveExpr"],
                     t=topgenes[transcript,"t"],
                     pval=topgenes[transcript,"P.Value"],
                     B=topgenes[transcript,"B"],
                     SE=topgenes[transcript,"SE"])
    
    inv_trans <- rbind(inv_trans,df)
  }
  return(inv_trans)
}


cohorts <- levels(factor(trans8$cohort))

rundiff_cohort <- function(coh){
  print(paste("Calculating differential analysis for",coh,"cohort"))
  trans8 <- trans8[,which(trans8$cohort==coh)]
  trans16 <- trans16[,which(trans16$cohort==coh)]
  trans17 <- trans17[,which(trans17$cohort==coh)]
  covs_trans <- covs_trans[colnames(trans8),]
  
  it8 <- diff_inv_trans(trans8, "inv8_001", covs_trans)
  it16 <- diff_inv_trans(trans16, "inv16_009", covs_trans)
  it17 <- diff_inv_trans(trans17, "inv17_007", covs_trans)
  inv_trans <- rbind(it8,it16,it17)

  #To order inversions as 8 - 16 - 17
  inv_trans$Inversion <- factor(inv_trans$Inversion, levels=c("8p23.1","16p11.2","17q21.31"))
  
  
  #Calculate padj adjusted by cpgs and exposures (all pval together)
  inv_trans$padj <- p.adjust(inv_trans$pval, method="bonferroni")
  
  #Sort by padj column
  inv_trans <- inv_trans[order(inv_trans$padj),]
  
  #Save the results in tables
  save(inv_trans,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/metanalysis/inv_trans_",coh,".Rdata"))
  
  it_sig <- inv_trans[which(inv_trans$padj<0.1),]
  writexl::write_xlsx(it_sig, path=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/metanalysis/inv_trans_",coh,".xlsx"))
}

#Perform the same steps but separately for each cohort
res_trans <- lapply(cohorts,rundiff_cohort)



#Create a function to load the inv_methy table for each cohort
load_dataset <- function(coh){
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/metanalysis/inv_trans_",coh,".Rdata"))
  return(inv_trans)
}

#Join all the inv_methy_coh in a list
allit <- parallel::mclapply(cohorts,load_dataset,mc.cores=6)
names(allit) <- cohorts


#Create a dataframe with one transcript in each comparing the 6 cohorts
subset_df <- function(transcript){
  df <- allit[["BIB"]][1,]
  df$cohort <- NA
  df <- df[0,]
  for (coh in cohorts){
    it <- allit[[coh]]
    dd <- it[which(it$Transcript==transcript),]
    dd$cohort <- coh
    df <- rbind(df,dd)
  }
  return(df)
}
transcripts <- allit[[1]]$Transcript

list <- lapply(transcripts, subset_df)
names(list) <- transcripts

meta_function <- function(transcript){
    met <- metagen(TE=list[[transcript]][,"LogFC"], 
                   seTE=list[[transcript]][,"SE"], sm="RR",
                   studlab=cohorts)
    res_met <- do.call(cbind,summary(met)$random)
    
    #Define the direction in each cohort
    dirs <- sign(list[[transcript]][,"LogFC"])
    pvals <- list[[transcript]][,"pval"]<0.05
    dirstudies <- dirs*pvals
    
    dirstudies[dirstudies==-1] <- "-"
    dirstudies[dirstudies==1] <- "+"
    dirstudies[dirstudies==0] <- "ns"
    
    names(dirstudies) <- paste0("DIR_",cohorts)
    df_inv <- cbind(list[[transcript]][1,c(1:4)], res_met, t(as.data.frame(dirstudies)))
    return(df_inv)
  }
  
list_inv <- lapply(transcripts, meta_function)
inv_trans <- do.call(rbind, list_inv)

# Adjust by Bonferroni
inv_trans$p.adj <- p.adjust(inv_trans$p, method="bonferroni")
inv_trans <- inv_trans[order(inv_trans$p.adj),]

#Change p column for pval column
colnames(inv_trans)[10] <- "pval"

#Remove df column
inv_trans <- inv_trans[,-12]

save(inv_trans,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/metanalysis/inv_trans.Rdata"))

it_sig <- inv_trans[which(inv_trans$p.adj<0.05),]
writexl::write_xlsx(it_sig, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/metanalysis/inv_trans.xlsx")




########################################
## DIFFERENTIAL ANALYSIS IN METHYLOME ##
########################################

#Load methylome datasets
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy17.Rdata")


#Select adjustment factors data for methylome
covs_methy <- data.frame(inv8_001=as.numeric(methy8$inv8_001)-1,
                         inv16_009=as.numeric(methy8$inv16_009)-1,
                         inv17_007=as.numeric(methy8$inv17_007)-1,
                         sex=methy8$sex,
                         age=methy8$age)
covs_methy <- cbind(covs_methy, pData(methy8)[,c(40:45,56:65)])

#Define a function to obtain a data frame with the results of the differential analysis
diff_inv_methy <- function(dataset, inversion,covs_methy){
  #Create an empty data frame
  inv_methy <- data.frame(Inversion=character(),
                            CpG=character(),
                            Location=character(),
                            Gene_Symbol=character(),
                            LogFC=numeric(),
                            CI.L=numeric(),
                            CI.R=numeric(),
                            AveExpr=numeric(),
                            t=numeric(),
                            pval=numeric(),
                            B=numeric(),
                            SE=numeric(),
                            stringsAsFactors=FALSE)
  
  #Differential analysis for all the features in the inversion region
  topcpgs <- diff_inv(dataset, inversion, omics="Methylome", covs_methy)
  for (cpg in rownames(topcpgs)){
    df <- data.frame(Inversion=inversionGR[inversion,]$Cytogenetic.location,
                     CpG = cpg,
                     Location=paste0(as.character(seqnames(dataset[cpg,])),":",as.character(start(dataset[cpg,]))),
                     Gene_Symbol=rowData(dataset)[cpg,"UCSC_RefGene_Name"],
                     LogFC=topcpgs[cpg,"logFC"],
                     CI.L=topcpgs[cpg,"CI.L"],
                     CI.R=topcpgs[cpg,"CI.R"],
                     AveExpr=topcpgs[cpg,"AveExpr"],
                     t=topcpgs[cpg,"t"],
                     pval=topcpgs[cpg,"P.Value"],
                     B=topcpgs[cpg,"B"],
                     SE=topcpgs[cpg,"SE"])
    
    inv_methy <- rbind(inv_methy,df)
  }
  return(inv_methy)
}

cohorts <- levels(factor(methy8$cohort))

rundiff_cohort <- function(coh){
  print(paste("Calculating differential analysis for",coh,"cohort"))
  methy8 <- methy8[,which(methy8$cohort==coh)]
  methy16 <- methy16[,which(methy16$cohort==coh)]
  methy17 <- methy17[,which(methy17$cohort==coh)]
  covs_methy <- covs_methy[colnames(methy8),]
  
  it8 <- diff_inv_methy(methy8, "inv8_001", covs_methy)
  it16 <- diff_inv_methy(methy16, "inv16_009", covs_methy)
  it17 <- diff_inv_methy(methy17, "inv17_007", covs_methy)
  inv_methy <- rbind(it8,it16,it17)
  
  #To order inversions as 8 - 16 - 17
  inv_methy$Inversion <- factor(inv_methy$Inversion, levels=c("8p23.1","16p11.2","17q21.31"))
  
  
  #Calculate padj adjusted by cpgs and exposures (all pval together)
  inv_methy$padj <- p.adjust(inv_methy$pval, method="bonferroni")
  
  #Sort by padj column
  inv_methy <- inv_methy[order(inv_methy$padj),]
  
  #Save the results in tables
  save(inv_methy,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy_",coh,".Rdata"))
  
  it_sig <- inv_methy[which(inv_methy$padj<0.1),]
  writexl::write_xlsx(it_sig, path=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy_",coh,".xlsx"))
}

#Perform the same steps but separately for each cohort
res_methy <- lapply(cohorts,rundiff_cohort)



#Create a function to load the inv_methy table for each cohort 
load_dataset <- function(coh){
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy_",coh,".Rdata"))
  return(inv_methy)
}

#Join all the inv_methy_coh in a list
allit <- parallel::mclapply(cohorts,load_dataset,mc.cores=6)
names(allit) <- cohorts


#Create a dataframe with one cpg in each comparing the 6 cohorts
subset_df <- function(cpg){
  df <- allit[["BIB"]][1,]
  df$cohort <- NA
  df <- df[0,]
  for (coh in cohorts){
    it <- allit[[coh]]
    dd <- it[which(it$CpG==cpg),]
    dd$cohort <- coh
    df <- rbind(df,dd)
  }
  return(df)
}
cpgs <- allit[[1]]$CpG

list <- lapply(cpgs, subset_df)
names(list) <- cpgs

meta_function <- function(cpg){
  met <- metagen(TE=list[[cpg]][,"LogFC"], 
                 seTE=list[[cpg]][,"SE"], sm="RR",
                 studlab=cohorts)
  res_met <- do.call(cbind,summary(met)$random)
  
  #Define the direction in each cohort
  dirs <- sign(list[[cpg]][,"LogFC"])
  pvals <- list[[cpg]][,"pval"]<0.05
  dirstudies <- dirs*pvals
  
  dirstudies[dirstudies==-1] <- "-"
  dirstudies[dirstudies==1] <- "+"
  dirstudies[dirstudies==0] <- "ns"
  
  names(dirstudies) <- paste0("DIR_",cohorts)
  df_inv <- cbind(list[[cpg]][1,c(1:4)], res_met, t(as.data.frame(dirstudies)))
  return(df_inv)
}

list_inv <- lapply(cpgs, meta_function)
inv_methy <- do.call(rbind, list_inv)

# Adjust by Bonferroni
inv_methy$p.adj <- p.adjust(inv_methy$p, method="bonferroni")
inv_methy <- inv_methy[order(inv_methy$p.adj),]

#Change p column for pval column
colnames(inv_methy)[10] <- "pval"

#Remove df column
inv_methy <- inv_methy[,-12]

save(inv_methy,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy.Rdata"))

it_sig <- inv_methy[which(inv_methy$p.adj<0.05),]
writexl::write_xlsx(it_sig, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy.xlsx")








