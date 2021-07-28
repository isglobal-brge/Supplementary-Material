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
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/exps.Rdata")

#Import the excel table to get the exposure abbreviations
excel_table <- read_excel("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Expodata.xlsx",sheet = "Exposure_Covariate")
excel_table <- data.frame(excel_table)
rownames(excel_table) <- excel_table$Variable_name_TRANS

exp_abrev<- excel_table[colnames(exps),]

#Select samples that intersect between exposome and methylome datasets
intersect_samples<-intersect(colnames(methy8),rownames(exps))
methy8 <- methy8[,intersect_samples]
methy16 <- methy16[,intersect_samples]
methy17 <- methy17[,intersect_samples]

imppost_final <- imppost_final[,intersect_samples]

#Select adjustment factors data
phenos <- data.frame(inv8_001=as.numeric(methy8$inv8_001)-1,
                     inv16_009=as.numeric(methy8$inv16_009)-1,
                     inv17_007=as.numeric(methy8$inv17_007)-1,
                     sex=methy8$sex,
                     age=methy8$age_sample_months,
                     maternal_education=as.numeric(imppost_final$h_edumc_None),
                     bmi=imppost_final$hs_c_bmi_None)

phenos <- cbind(phenos, pData(methy8)[,c(40:45,56:65)])

#Store the name of the exposures
exp_variables <- colnames(exps)

#Select the exposures from the intersect_samples
exps <- exps[intersect_samples,]

#Create the functions needed for the interaction differential analysis
diff_inv_expo <- function(dataset, exps, phenos, inversion, expo){
  print(paste(inversion,"-",expo))
  mat <- assay(dataset)
  phenoModel <- data.frame(expi=exps[,expo],
                           inv=phenos[[inversion]])
  phenoModel <- cbind(phenoModel,phenos[,4:ncol(phenos)])
  
  rownames(phenoModel) <- phenos$HelixID
  mod <- model.matrix(formula(paste("~ expi:inv +",paste(colnames(phenoModel),collapse=" + "))),
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
                            Period=character(),
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
                     Exposure_abrev=exp_abrev[expo,"Label.for.tables"],
                     Family=exp_abrev[expo,"Group"],
                     Period=exp_abrev[expo,"Period"],
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
  
  it <- parallel::mclapply(exp_variables, call, mc.cores=10, exps, phenos, methy8, methy16, methy17)
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
  save(inter_methy,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy_",coh,".Rdata"))
  
  it_sig <- inter_methy[which(inter_methy$padj<0.1),]
  writexl::write_xlsx(it_sig, path=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy_",coh,".xlsx"))
}

#Perform the same steps but separately for each cohort
res_methy <- lapply(cohorts,rundiff_cohort)

######################################
## 2. META-ANALYSIS BY LOGFC AND SE ##
######################################

#cohorts <- c("BIB","EDEN","KANC","MOBA","RHEA","SAB")

#Create a function to load the inter_methy table for each cohort and add an interaction column
load_dataset <- function(coh){
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy_",coh,".Rdata"))
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
    pvals <- list[[cpg]][,"pval"]<0.5
    dirstudies <- dirs*pvals
    
    dirstudies[dirstudies==-1] <- "-"
    dirstudies[dirstudies==1] <- "+"
    dirstudies[dirstudies==0] <- "ns"
    
    names(dirstudies) <- paste0("DIR_",cohorts)
    df_inter <- cbind(list[[cpg]][1,c(1:8)], res_met, t(as.data.frame(dirstudies)))
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

it_metanalysis <- parallel::mclapply(inters, call_meta, mc.cores=10)

inter_methy <- do.call(rbind, it_metanalysis)

# Adjust by Bonferroni
inter_methy$p.adj <- p.adjust(inter_methy$p, method="bonferroni")
inter_methy <- inter_methy[order(inter_methy$p.adj),]

#Change p column for pval column
colnames(inter_methy)[14] <- "pval"

#Remove level and df columns
inter_methy <- inter_methy[,-c(15,16)]

#Select only one gene symbol per CpG
allsymbols_to_symbol <- function(allsymbols) {
  symbols <- strsplit(allsymbols,";")
  symbols <- lapply(symbols, function(x) unique(x[!is.na(x)]))
  symbols <- lapply(symbols, function(x) paste(x,collapse=";"))
  vector <- do.call(c,symbols)
  return (vector)
}

inter_methy$Gene_Symbol <- allsymbols_to_symbol(inter_methy$Gene_Symbol)

save(inter_methy,file=paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy.Rdata"))

it_sig <- inter_methy[which(inter_methy$p.adj<0.05),]
writexl::write_xlsx(it_sig, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy.xlsx")

#Create a table counting the number of interactions per exposure
select_df <- read_excel("/home/isglobal.lan/ncarreras/homews/TFM/Data_tables/Expodata_select.xlsx")
select_df <- data.frame(select_df)
select_df <- select_df[which(select_df$Selection=="yes"),c(1:2,7,4)]
rownames(select_df) <- select_df$Variable_name_TRANS

select_df$Num_inters <- NA
select_df$Inversion_inters <- NA

for (i in 1:nrow(select_df)){
  expo <- rownames(select_df)[i]
  num_inters <- nrow(it_sig[which(it_sig$Exposure==expo),])
  select_df$Num_inters[i] <- num_inters
  inv <- paste0(sort(unique(it_sig[which(it_sig$Exposure==expo),]$Inversion)),collapse=";")
  select_df$Inversion_inters[i] <- inv
}

writexl::write_xlsx(select_df, "/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_by_expo.xlsx")
