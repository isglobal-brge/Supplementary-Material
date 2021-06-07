########################################
## DIFFERENTIAL ANALYSIS BY INVERSION ##
########################################

#Load libraries
library(MEAL)
library(SummarizedExperiment)
library(rexposome)

#Load inversion information
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")

#Load covariate dataset
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppreg_final.Rdata")

#Define a function to perform differential analysis
diff_inv <- function(dataset, inversion, omics, covs){
  print(paste(inversion, "-", omics))
  if (omics=="Transcriptome"){
    mat <- assays(dataset)$expr
  }
  if (omics=="Methylome"){
    mat <- assay(dataset)
  }
  mod <- model.matrix(formula(paste("~", inversion,"+ sex + age + cohort + trim_conception + parity + maternal_education + maternal_bmi + maternal_age + maternal_smoke")), data=covs)
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

#Select samples that overlap between covariates dataset and transcriptome
ids_trans <- intersect(colnames(imppreg_final),colnames(trans8))

trans8 <- trans8[,ids_trans]
trans16 <- trans16[,ids_trans]
trans17 <- trans17[,ids_trans]

imppreg_trans <- imppreg_final[,ids_trans]

#Select adjustment factors data for transcriptome
covs_trans <- data.frame(sex=imppreg_trans$sex,
                         age=imppreg_trans$hs_child_age_days_None,
                         cohort=imppreg_trans$cohort,
                         trim_conception=imppreg_trans$h_trimcon_None,
                         parity=imppreg_trans$h_parity_None,
                         maternal_education=imppreg_trans$h_edumc_None,
                         maternal_bmi=imppreg_trans$h_mbmi_None,
                         maternal_age=imppreg_trans$h_age_None,
                         maternal_smoke=expos(imppreg_trans)$e3_asmokyn_p_None,
                         inv8_001=as.numeric(imppreg_trans$inv8_001)-1,
                         inv16_009=as.numeric(imppreg_trans$inv16_009)-1,
                         inv17_007=as.numeric(imppreg_trans$inv17_007)-1)
rownames(covs_trans) <- colnames(imppreg_trans)

#Define a function to obtain a data frame with the results of the differential analysis
diff_inv_trans <- function(dataset, inversion){
  #Create an empty data frame
  inter_trans <- data.frame(Inversion=character(),
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
    
    inter_trans <- rbind(inter_trans,df)
  }
  return(inter_trans)
}

#Differential analysis separately for each inversion
it8 <- diff_inv_trans(trans8, "inv8_001")
it16 <- diff_inv_trans(trans16, "inv16_009")
it17 <- diff_inv_trans(trans17, "inv17_007")

#Bind all the datasets
inv_trans <- rbind(it8,it16,it17)

#Sort inversions
inv_trans$Inversion <- factor(inv_trans$Inversion, levels=c("8p23.1","16p11.2","17q21.31"))

#Calculate p-value adjusted by Bonferroni
inv_trans$p.adj <- p.adjust(inv_trans$pval, method="bonferroni")

#Sort by p.adj column
inv_trans <- inv_trans[order(inv_trans$p.adj),]

#Save the results in tables
save(inv_trans,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/inv_trans.Rdata"))

inv_trans_sig <- inv_trans[which(inv_trans$p.adj<0.05),]
writexl::write_xlsx(inv_trans_sig, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/inv_trans.xlsx")


########################################
## DIFFERENTIAL ANALYSIS IN METHYLOME ##
########################################

#Load methylome datasets
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy17.Rdata")

#Select samples that overlap between covariates dataset and methylome
ids_methy <- intersect(colnames(imppreg_final),colnames(methy8))

methy8 <- methy8[,ids_methy]
methy16 <- methy16[,ids_methy]
methy17 <- methy17[,ids_methy]

imppreg_methy <- imppreg_final[,ids_methy]

#Select adjustment factors data for methylome
covs_methy <- data.frame(sex=imppreg_methy$sex,
                         age=imppreg_methy$hs_child_age_days_None,
                         cohort=imppreg_methy$cohort,
                         trim_conception=imppreg_methy$h_trimcon_None,
                         parity=imppreg_methy$h_parity_None,
                         maternal_education=imppreg_methy$h_edumc_None,
                         maternal_bmi=imppreg_methy$h_mbmi_None,
                         maternal_age=imppreg_methy$h_age_None,
                         maternal_smoke=expos(imppreg_methy)$e3_asmokyn_p_None,
                         inv8_001=as.numeric(imppreg_methy$inv8_001)-1,
                         inv16_009=as.numeric(imppreg_methy$inv16_009)-1,
                         inv17_007=as.numeric(imppreg_methy$inv17_007)-1)
rownames(covs_methy) <- colnames(imppreg_methy)


#Define a function to obtain a data frame with the results of the differential analysis
diff_inv_methy <- function(dataset, inversion){
  #Create an empty data frame
  inter_methy <- data.frame(Inversion=character(),
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
    
    inter_methy <- rbind(inter_methy,df)
  }
  return(inter_methy)
}

#Differential analysis separately for each inversion
it8 <- diff_inv_methy(methy8, "inv8_001")
it16 <- diff_inv_methy(methy16, "inv16_009")
it17 <- diff_inv_methy(methy17, "inv17_007")

#Bind datasets
inv_methy <- rbind(it8,it16,it17)

#Sort inversions
inv_methy$Inversion <- factor(inv_methy$Inversion, levels=c("8p23.1","16p11.2","17q21.31"))

#Calculate p-value adjusted by Bonferroni
inv_methy$p.adj <- p.adjust(inv_methy$pval, method="bonferroni")

#Sort by p.adj column
inv_methy <- inv_methy[order(inv_methy$p.adj),]

#Save the results in tables
save(inv_methy,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/inv_methy.Rdata"))

inv_methy_sig <- inv_methy[which(inv_methy$p.adj<0.05),]
writexl::write_xlsx(inv_methy_sig, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/inv_methy.xlsx")
