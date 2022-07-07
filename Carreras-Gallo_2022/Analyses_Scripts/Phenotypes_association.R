######################################################
## PHENOTYPES ASSOCIATED WITH SUBSTANCE CONSUMPTION ##
######################################################

library(minfi)
library(tidyverse)
library(gridExtra)
library(ggplot2)

#Set working directory
setwd("/PROJECTES/GENOMICS/TruDiagnostic/Analysis/Uniform_analysis/Phenotypes")

#Load GenomicRatioSet
load("/PROJECTES/GENOMICS/TruDiagnostic/Final_datasets/GRset_SVA.Rdata")

GRset$Obesity <- GRset$BMI_class=="obesity"

#Remove individuals with Marijuana_frequency == NA (individuals who use other drugs)
GRset_mari <- GRset[,which(!is.na(GRset$Marijuana_frequency))]

#Create a function for creating an ODDs plot
pheno_subst <- function(substance, var, phenos, covs, levels,GRset){
  data <- pData(GRset)[,c(var,covs,phenos)]
  list <- list()
  for (i in 1:length(phenos)){
    pheno <- phenos[i]
    
    #Skip phenotype if there are more than 2 pheno-var with 0 individuals
    if(sum(table(GRset[[var]],GRset[[pheno]])<3)>=3){
      print(paste0(pheno," has not enough individuals"))
      next
    }
    
    if(pheno=="Obesity"){
      covs <- covs[!covs %in% "BMI"]
    }
    
    formula <- formula(paste0(pheno,"~",var," + ",paste0(covs,collapse="+")))
    res <- glm(formula = formula, data = data, family="binomial")
    sumres <- data.frame(summary(res)$coefficients[2:length(levels(GRset[[var]])),])
    
    #Remove rows were standard error is bigger than 2 bc it means that there are not enough participants to evaluate the association
    pos <- sumres$Std..Error<2
    sumres <- sumres[pos,]
    levels_a <- levels[pos]
    
    #Calculate odds ratio based on estimate value
    sumres$Odds <- exp(sumres$Estimate)
    
    #Calculate confidence interval [low,high] based on the standard error
    sumres$CIlow <- exp(sumres$Estimate - sumres$Std..Error*2)
    sumres$CIhigh <- exp(sumres$Estimate + sumres$Std..Error*2)
    sumres$Level <- levels_a
    sumres$Inter <- paste0("[",round(sumres$CIlow,digits=2),",",round(sumres$CIhigh,digits=2),"]")
    
    
    #Create the plot for the odds ratio if significant
    if(any(sumres$Pr...z..<0.05)){
      p <- ggplot(sumres, aes(x = Odds, y = 1:nrow(sumres))) +
        geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
        geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = .2, color = "gray50") +
        geom_point(size = 3.5, color = "orange") +
        theme_bw() +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(breaks = 1:nrow(sumres), labels = levels_a) +
        scale_x_continuous(breaks = seq(-2,3,0.2)) +
        ylab(paste0(substance," Consumption")) +
        xlab("Odds ratio (95% interval)") +
        ggtitle(paste0(pheno," Odds Ratio"))
      
      #Create the tables next to the plot
      ## Create the table-base pallete
      table_base <- ggplot(sumres, aes(y=Level)) +
        ylab(NULL) + xlab("  ") + 
        theme(plot.title = element_text(hjust = 0.5, size=12), 
              axis.text.x = element_text(color="white", hjust = -3, size = 25), ## This is used to help with alignment
              axis.line = element_blank(),
              axis.text.y = element_blank(), 
              axis.ticks = element_blank(),
              axis.title.y = element_blank(), 
              legend.position = "none",
              panel.background = element_blank(), 
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              plot.background = element_blank())
      
      ## ODDs ratio table
      tab1 <- table_base + 
        labs(title = "space") +
        geom_text(aes(y = 1:nrow(sumres), x = 1, label = sprintf("%0.2f", round(Odds, digits = 2))), size = 4) + ## decimal places
        ggtitle("OR")
      
      ## 95% CI table
      tab2 <- table_base +
        geom_text(aes(y = 1:nrow(sumres), x = 1, label = Inter), size = 4) + 
        ggtitle("95% CI")
      
      ## p-value
      tab3 <- table_base +
        geom_text(aes(y = 1:nrow(sumres), x = 1, label = round(Pr...z..,digits=5)), size = 4) + 
        ggtitle("P-value")
      
      ## Merge tables with plot
      lay <-  matrix(c(1,1,1,1,1,1,1,1,1,1,2,3,3,4,4), nrow = 1)
      png(paste0(substance,"/",substance,"_",pheno,".png"),
          res = 1200, width = 10, height = 5, units = "in")
      
      grid.arrange(p, tab1, tab2, tab3, layout_matrix = lay)
      dev.off()
    }
    sumres$Phenotype <- pheno
    rownames(sumres) <- NULL
    list[[i]] <- sumres
  }
  df <- do.call(rbind,list)
  
  #Change columns order
  df <- df[,c(10,8,1:7,9)]
  
  #Change column names
  colnames(df)[6] <- "P.Value"
  
  #Sort by p-value
  df <- df[order(df$P.Value),]
  
  #Adjust p-value by multiple comparisons (dividing by the number of phenotypes)
  df$P.Value.Adj <- df$P.Value*length(phenos)
  
  save(df, file=paste0(substance,"/",substance,"_phenos.Rdata"))
  writexl::write_xlsx(df,path=paste0(substance,"/",substance,"_phenos.xlsx"))
}

#Covariates
covs <- c("Biological_Sex","age", "Ethnicity",
          "BMI", "Level_of_Education","Tobacco_Use","Alcohol_Use_per_week")

#Factor levels (A for alcohol and marijuana; B for tobacco)
levelsA <- c("On special\noccasions","Once per\nweek","3-5 times\nper week","Regularly")
levelsB <- c("Less than 1\nper week","Less than 1\nper day",
             "1-5\nper day","6-10\nper day","11-20\nper day",
             "More than 20\nper day")

#############
## ALCOHOL ##
#############

phenos_alco <- c("Cardiovascular_any", "High_Blood_Pressure",
                 "Respiratory_Disease_any", "Pneumonia",
                 "Endocrine_Disease_any","Diabetes_2",
                 "Gastrointestinal_any",
                 "Musculoskeletal_any","Gout","Rheumatoid_Arthritis",
                 "Neuropsychological_any","Depression","Anxiety",
                 "Cancer_Diagnosis_any","Obesity")

pheno_subst(substance="Alcohol",var="Alcohol_Use_per_week",phenos=phenos_alco,covs=covs[-7],levels=levelsA,GRset=GRset)

###############
## MARIJUANA ##
###############

phenos_mari <- c("Cardiovascular_any","Heart_Attack","High_Blood_Pressure",
                 "Respiratory_Disease_any","Chronic_Bronchitis","Pneumonia",
                 "Neuropsychological_any","Depression","Anxiety",
                 "Cancer_Diagnosis_any")

pheno_subst(substance="Marijuana",var="Marijuana_frequency",phenos=phenos_mari,covs=covs,levels=levelsA, GRset=GRset_mari)

#############
## TOBACCO ##
#############

phenos_toba <- c("Cardiovascular_any","High_Blood_Pressure",
                 "Respiratory_Disease_any","Chronic_Bronchitis","Pneumonia",
                 "Endocrine_Disease_any", "Diabetes_2",
                 "Musculoskeletal_any",
                 "Neuropsychological_any","Depression","Anxiety",
                 "Cancer_Diagnosis_any")

pheno_subst(substance="Tobacco",var="Tobacco_Use",phenos=phenos_toba,covs=covs[-6],levels=levelsB,GRset=GRset)
