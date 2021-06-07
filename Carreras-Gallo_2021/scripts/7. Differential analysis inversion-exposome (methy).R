############################################################################
########### DIFFERENTIAL ANALYSIS INVERSION*EXPOSOME ON METHYLOME ##########
############################################################################

#Load libraries
library(SummarizedExperiment)
library(rexposome)
library(MEAL)
library(MetaDE)
library(readxl)
library(meta)

#Load data
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy17.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppost_final.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppreg_final.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")

#Import the excel table to get the exposure abreviations
excel_table <- read_excel("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Expodata.xlsx",sheet = "Exposure_Covariate")
exp_abrev <- cbind(excel_table[,"Variable_name_TRANS"],excel_table[,"Label for tables"],excel_table[,"Group"],excel_table[,"Subgroup"])
rownames(exp_abrev) <- exp_abrev[,"Variable_name_TRANS"]
# source("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/aux_functions.R")

#Discard exposure facility richness
imppost_final <- imppost_final[-38,]

#Select samples that intersect between exposome and methylome datasets
intersect_samples<-intersect(colnames(methy8),colnames(imppost_final))
methy8 <- methy8[,intersect_samples]
methy16 <- methy16[,intersect_samples]
methy17 <- methy17[,intersect_samples]

#Select adjustment factors data
imppreg_methy <- imppreg_final[,intersect_samples]

phenos <- data.frame(sex=imppreg_methy$sex,
                     age=imppreg_methy$hs_child_age_days_None,
                     cohort=imppreg_methy$cohort,
                     trim_conception=imppreg_methy$h_trimcon_None,
                     parity=imppreg_methy$h_parity_None,
                     maternal_education=as.numeric(imppreg_methy$h_edumc_None),
                     maternal_bmi=imppreg_methy$h_mbmi_None,
                     maternal_age=imppreg_methy$h_age_None,
                     maternal_smoke=expos(imppreg_methy)$e3_asmokyn_p_None,
                     inv8_001=as.numeric(imppreg_methy$inv8_001)-1,
                     inv16_009=as.numeric(imppreg_methy$inv16_009)-1,
                     inv17_007=as.numeric(imppreg_methy$inv17_007)-1)
rownames(phenos) <- colnames(imppreg_methy)

#Store the name of the postnatal exposures
exp_variables_post<-colnames(expos(imppost_final))

#Select the exposures and the phenotypes for the intersect_samples
exps <- expos(imppost_final[,intersect_samples])

#Create the functions needed for the interaction differential analysis
diff_inv_expo <- function(dataset, exps, phenos, inversion, expo){
  print(paste(inversion,"-",expo))
  mat <- assay(dataset)
  phenoModel <- data.frame(expi=as.numeric(exps[[expo]]), 
                           inv=phenos[[inversion]])
  phenoModel <- cbind(phenoModel,phenos[,1:9])
  
  rownames(phenoModel) <- phenos$HelixID
  mod <- model.matrix(~ expi:inv + expi + inv + sex + age + 
                        trim_conception + parity + maternal_education + 
                        maternal_bmi + maternal_age + maternal_smoke, 
                      data=phenoModel)
  resmean <- try(runDiffMeanAnalysis(set = mat, model = mod, 
                                     resultSet = TRUE, method="robust"),
                 silent=TRUE)
  if (class(resmean)=="try-error"){
    if(resmean[1]=="Error in rlm.default(x = X, y = y, weights = w, ...) : \n  'x' is singular: singular fits are not implemented in 'rlm'\n"){
      return("")
    }
  }
  topfeatures <- getAssociation(resmean, coef="expi:inv", 
                                fNames=NULL, rid="DiffMean")
  return(topfeatures)
}

diff_inter <- function(dataset, exps, phenos, inversion, expo){
  #Create an empty data frame
  inter_methy <- data.frame(Exposure=character(),
                            Exposure_abrev=character(),
                            Family=character(),
                            Inversion=character(),
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
  topcpgs <- diff_inv_expo(dataset, exps, phenos, inversion, expo)
  for (cpg in rownames(topcpgs)){
    df <- data.frame(Exposure=expo,
                     Exposure_abrev=exp_abrev[expo,][["Label for tables"]],
                     Family=fData(imppost_final)[expo,"Group"],
                     Inversion=inversionGR[inversion,]$Cytogenetic.location,
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

call <- function(expo, exps, phenos, methy8, methy16, methy17){
  it8 <- diff_inter(methy8, exps, phenos, "inv8_001", expo)
  it16 <- diff_inter(methy16, exps, phenos, "inv16_009", expo)
  it17 <- diff_inter(methy17, exps, phenos, "inv17_007", expo)
  it <- rbind(it8,it16,it17)
  return(it)
}

#Create functions to adjust only by the number of CpG sites and the number of exposures
padj.cpgs <- function(it){
  it$padj.cpgs <- NA
  for (expo in unique(it$Exposure)){
    pvals <- p.adjust(it[which(it$Exposure==expo),]$pval, method="bonferroni")
    it[which(it$Exposure==expo),]$padj.cpgs <- pvals
  }
  return(it)
}

padj.expos <- function(it){
  it$padj.expos <- NA
  for (methy in unique(it$CpG)){
    pvals <- p.adjust(it[which(it$CpG==methy),]$pval, method="bonferroni")
    it[which(it$CpG==methy),]$padj.expos <- pvals
  }
  return(it)
}

#In this analysis, we are going to perform the differential analysis separately 
#by cohorts and, after that, perform a meta-analysis to combine the results

#########################################
## 1. DIFFERENTIAL ANALYSIS BY COHORTS ##
#########################################

cohorts <- levels(factor(methy8$cohort))

rundiff_cohort <- function(coh){
  print(paste("Calculating differential analysis for",coh,"cohort"))
  methy8 <- methy8[,which(methy8$cohort==coh)]
  methy16 <- methy16[,which(methy16$cohort==coh)]
  methy17 <- methy17[,which(methy17$cohort==coh)]
  exps <- exps[colnames(methy8),]
  phenos <- phenos[colnames(methy8),]
  
  it <- parallel::mclapply(exp_variables_post, call, mc.cores=10, exps, phenos, methy8, methy16, methy17)
  inter_methy <- do.call(rbind,it)

  #To order inversions as 8 - 16 - 17
  inter_methy$Inversion <- factor(inter_methy$Inversion, levels=c("8p23.1","16p11.2","17q21.31"))
  
  
  #Calculate padj adjusted by cpgs and exposures (all pval together)
  inter_methy$padj <- p.adjust(inter_methy$pval, method="bonferroni")

  #Calculate padj only adjusted by cpgs
  inter_methy <- padj.cpgs(inter_methy)
  
  #Calculate padj only adjusted by expos
  inter_methy <- padj.expos(inter_methy)
  
  #Sort by padj column
  inter_methy <- inter_methy[order(inter_methy$padj),]
  
  #Save the results in tables
  save(inter_methy,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy_",coh,".Rdata"))
  
  it_sig <- inter_methy[which(inter_methy$padj<0.1),]
  writexl::write_xlsx(it_sig, path=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy_",coh,".xlsx"))
}

#Perform the same steps but separately for each cohort
res_methy <- lapply(cohorts,rundiff_cohort)

######################################
## 2. META-ANALYSIS BY LOGFC AND SE ##
######################################

#Create a function to load the inter_methy table for each cohort and add an interaction column
load_dataset <- function(coh){
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy_",coh,".Rdata"))
  inter_methy$Interaction <- as.character(interaction(inter_methy$Exposure, 
                                                      inter_methy$Inversion))
  return(inter_methy)
}

#Join all the inter_methy_coh in a list
allit <- parallel::mclapply(cohorts,load_dataset,mc.cores=6)
names(allit) <- cohorts

#Find the interactions that are present in all the cohorts (inters)
allinters <- lapply(allit, function(it) unique(it$Interaction) )
num_allinters <- table(do.call(c,allinters))
inters <- names(num_allinters[which(num_allinters==6)])

#Create a function to combine the results of the different cohorts
meta_interaction <- function(interaction, cohorts, cpgs){
  
  print(interaction)
  
  #Create an empty list
  alldd <- vector(mode = "list", length = length(cohorts))
  names(alldd) <- cohorts
  
  #Subset a data frame for the specific interaction and store it in alldd
  subset_df <- function(cpg){
    df <- allit[["BIB"]][1,]
    df$cohort <- NA
    df <- df[0,]
    for (coh in cohorts){
      it <- allit[[coh]]
      dd <- it[which(it$Interaction==interaction &
                       it$CpG==cpg),]
      dd$cohort <- coh
      df <- rbind(df,dd)
    }
    return(df)
  }
  
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
    df_inter <- cbind(list[[cpg]][1,c(1:7)], res_met, t(as.data.frame(dirstudies)))
    return(df_inter)
  }
  
  list_inter <- lapply(cpgs, meta_function)
  df_final <- do.call(rbind, list_inter)
  return(df_final)
}

call_meta <- function(inter){
  cpgs <- allit[[1]][which(allit[[1]]$Interaction==inter),]$CpG
  df <- meta_interaction(inter,cohorts, cpgs=cpgs)
  return(df)
}

it_metanalysis <- parallel::mclapply(inters[1:3], call_meta, mc.cores=10)

inter_methy <- do.call(rbind, it_metanalysis)

# Adjust by Bonferroni
inter_methy$p.adj <- p.adjust(inter_methy$p, method="bonferroni")
inter_methy <- inter_methy[order(inter_methy$p.adj),]

#Change p column for pval column
colnames(inter_methy)[13] <- "pval"

save(inter_methy,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy.Rdata"))

it_sig <- inter_methy[which(inter_methy$p.adj<0.05),]
writexl::write_xlsx(it_sig, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy.xlsx")
