###################################################################
## Interaction Cadmium*8p23.1 and cg16810626 (GATA4) methylation ##
###################################################################

## Load libraries

library(ggplot2)
library(meta)

## Load datasets

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppost_final.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppreg_final.Rdata")

## Select participants that intersect between exposome and methylome
ids_expo_methy <- intersect(colnames(imppost_final),colnames(methy8))
imppost_final <- imppost_final[,ids_expo_methy]
imppreg_methy <- imppreg_final[,ids_expo_methy]

## Create a data.frame all the variables (including covariates)
covs <- data.frame(sex=imppreg_methy$sex,
                   age=imppreg_methy$hs_child_age_days_None,
                   cohort=imppreg_methy$cohort,
                   trim_conception=imppreg_methy$h_trimcon_None,
                   parity=imppreg_methy$h_parity_None,
                   maternal_education=as.numeric(imppreg_methy$h_edumc_None),
                   maternal_bmi=imppreg_methy$h_mbmi_None,
                   maternal_age=imppreg_methy$h_age_None,
                   maternal_smoke=expos(imppreg_methy)$e3_asmokyn_p_None,
                   inv8_001=imppreg_methy$inv8_001,
                   cadmium=expos(imppost_final)$hs_cd_c_Log2)

rownames(covs) <- colnames(imppreg_methy)
data <- cbind(covs, t(assay(methy8["cg16810626",ids_expo_methy])))

## Calculate residual values for cg16810626 methylation 
data$cg16810626_residual <- residuals(lm(formula=paste0("cg16810626 ~ ",paste(colnames(covs),collapse=" + ")), data=data))

theme_plots <- theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
                     axis.title.x = element_text(size = 18),
                     axis.title.y = element_text(size = 18),
                     axis.text.x = element_text(size = 12),
                     axis.text.y = element_text(size = 12),
                     strip.text.y = element_text(size = 18))

#####################################################
## INTERACTION PLOT (Cadmium*8p23.1 VS cg16810626) ##
#####################################################

interaction_plot <- function(data,path){
  png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Specific_cases/GATA4_SOX7/cg16810626_cadmium/interaction_methy",path,".png"))
  ggplot(data, aes(x=cadmium, y=cg16810626_residual, col=inv8_001))+
    geom_point()+
    labs(x="log2(Cadmium)", y="cg16810626 methylation (residuals)", col="Genotype")+
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~inv8_001, ncol=1, strip.position = "right") +
    scale_color_manual(values = c("N/N" = "#00AFBB", "N/I" = "#E7B800", "I/I" = "#FC4E07")) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text.y = element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  dev.off()
  res<- lm(formula=cg16810626_residual ~ inv8_001*cadmium, data=data)
  return(summary(res))
}

#Interaction plot with all samples
all_samples <- interaction_plot (data=data, path="")

#Interaction plot without outliers
noOutliers <- interaction_plot(data=data[which(!data$cadmium>-1),], path="_noOutliers")

###################################################
#### FOREST PLOT FOR CADMIUM*8p23.1 - cg16810626 ##
###################################################

cohorts <- c("BIB", "EDEN", "KANC", "MOBA", "RHEA", "SAB")

#Load results by cohorts
#Create a function to load the inter_methy table and add an Interaction column
load_dataset <- function(coh){
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy_",coh,".Rdata"))
  inter_methy$Interaction <- as.character(interaction(inter_methy$Exposure, 
                                                      inter_methy$Inversion))
  return(inter_methy)
}

#Join all the inter_methy_coh in a list
allit <- parallel::mclapply(cohorts,load_dataset,mc.cores=6)
names(allit) <- cohorts

subset_function <- function(coh,inter,cpg){
  it_all <- allit[[coh]]
  it <- it_all[which(it_all$Interaction==inter&
                       it_all$CpG==cpg),]
}

gata4_cohorts <- t(sapply(cohorts,subset_function,inter="hs_cd_c_Log2.8p23.1",cpg="cg16810626"))

met <- metagen(TE=do.call(rbind,gata4_cohorts[,"LogFC"]), 
               seTE=do.call(rbind,gata4_cohorts[,"SE"]), sm="MD",
               studlab=rownames(gata4_cohorts), comb.fixed = FALSE)

png("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Specific_cases/GATA4_SOX7/forest_gata4_cadmium.png",
    height = 600, width = 1200, units = "px")
forest.meta(met, studlab=TRUE)
dev.off()

#######################
## CADMIUM HISTOGRAM ##
#######################

png("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Specific_cases/GATA4_SOX7/cg16810626_cadmium/hist_cadmium.png")
hist(data$cadmium, xlab="Cadmium", main = "Cadmium exposure histogram")
dev.off()

#Plot by cohorts
data <- data[order(data$cohort),]

png("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Specific_cases/GATA4_SOX7/cg16810626_cadmium/cadmium_cohorts.png")
ggplot(data=data, aes(x=1:nrow(data),y=cadmium, col=cohort))+
  geom_point()+
  labs(y="Cadmium Exposure", x="Samples")+
  theme_plots
  
dev.off()

################################################
## GATA4 EXPRESSION vs cg16810626 METHYLATION ##
################################################

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/trans8.Rdata")
transcript_GATA4 <- rownames(rowData(trans8)[which(rowData(trans8)$GeneSymbol_Affy=="GATA4"),])

IDs_3 <- intersect(rownames(data),colnames(trans8))
data_expr <- data[IDs_3,]
data_expr$GATA4 <- assay(trans8)[transcript_GATA4,IDs_3]

png("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Specific_cases/GATA4_SOX7/cg16810626_cadmium/expr_methy.png")
ggplot(data=data_expr, aes(x=cg16810626_residual, y=GATA4))+
  geom_point()+
  labs(x="cg16810626 methylation (residual)", y="GATA4 expression")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_plots
dev.off()

