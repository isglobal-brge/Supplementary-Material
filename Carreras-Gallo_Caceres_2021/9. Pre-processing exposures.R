#Load libraries
library(readxl)
library(rexposome)
library(ggplot2)

#Load data
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppost_final.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Imppreg_final.Rdata")

#Import the excel table to get the exposure abbreviations
excel_table <- read_excel("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Expodata.xlsx",sheet = "Exposure_Covariate")
excel_table <- data.frame(excel_table)
rownames(excel_table) <- excel_table$Variable_name_TRANS

exp_abrev<- excel_table[c(rownames(imppost_final),rownames(imppreg_final)),]
dim(exp_abrev)

################################################################################
################### REMOVE EXPOSURES WITH A LOT OF MISSINGS ####################
################################################################################

#Missings percentages from character to numeric
for (i in 24:30){
  exp_abrev[,i] <- as.numeric(gsub(",", ".", exp_abrev[,i]))
}

table(exp_abrev[,24]<10)

#Discard exposures with too much missing
miss_10_no_coh_20 <- exp_abrev[which(exp_abrev$X..missing.HELIX.Subcohort<10 &
                                       exp_abrev$BIB<20 &
                                       exp_abrev$EDEN<20 &
                                       exp_abrev$INMA<20 &
                                       exp_abrev$KANC<20 &
                                       exp_abrev$MOBA<20 &
                                       exp_abrev$RHEA<20),]
dim(miss_10_no_coh_20)

miss_10_no_coh_20 <- miss_10_no_coh_20[order(miss_10_no_coh_20$Group),]
miss_10_no_coh_20 <- miss_10_no_coh_20[order(miss_10_no_coh_20$Period),]

# writexl::write_xlsx(miss_10_no_coh_20[,c("Period","Group","Subgroup",
#                                          "Variable_name_TRANS","Label.for.tables",
#                                          "Transformation","Type","ref.category",
#                                          "Categories","Units","description",
#                                          "X..missing.HELIX.Subcohort", "BIB",
#                                          "EDEN","INMA","KANC","MOBA","RHEA")], 
#                     "/home/isglobal.lan/ncarreras/homews/TFM/Data_tables/Expodata_no_missings.xlsx")

################################################################################
############# SELECT ONLY THE MOST IMPORTANT EXPOSURES PER GROUP ###############
################################################################################

select_df <- read_excel("/home/isglobal.lan/ncarreras/homews/TFM/Data_tables/Expodata_select.xlsx")
select_df <- data.frame(select_df)
select_df <- select_df[which(select_df$Selection=="yes"),]
exp_variables <- rownames(select_df) <- select_df$Variable_name_TRANS

#Select the exposures from the intersect_samples
post_preg <- cbind(expos(imppost_final),expos(imppreg_final))
post_preg <- post_preg[,exp_variables]

################################################################################
######################### WINSORIZE NUMERIC CATEGORIES #########################
################################################################################

#Select exposures annotated as numeric
nums <- rownames(select_df[which(select_df$Type=="Numeric"),])

#Function to winsorize outliers
winsorize <- function(x, probs = NULL, cutpoints = NULL , replace = c(cutpoints[1], cutpoints[2]), verbose = TRUE){
  dummy = is.integer(x)
  if (!is.null(probs)){
    stopifnot(is.null(cutpoints))
    stopifnot(length(probs)==2)
    cutpoints <- quantile(x, probs, type = 1, na.rm = TRUE)
  } else if (is.null(cutpoints)){
    l <- quantile(x, c(0.25, 0.50, 0.75), type = 1, na.rm = TRUE) 
    cutpoints <- c(l[2]-3*(l[3]-l[1]), l[2]+3*(l[3]-l[1])) 
  } else{
    stopifnot(length(cutpoints)==2)
  }
  if (is.integer(x)) cutpoints <- round(cutpoints)
  bottom <- x < cutpoints[1]
  top <- x > cutpoints[2]
  if (verbose){
    length <- length(x)
    message(paste(100*sum(bottom, na.rm = TRUE)/length,"% observations replaced at the bottom"))
    message(paste(100*sum(top, na.rm = TRUE)/length,"% observations replaced at the top"))
  }
  x[bottom] <- replace[1]
  x[top] <- replace[2]
  if (dummy){
    x <- as.integer(x)
  }
  x
}

for (i in nums){
  post_preg[,i] <- winsorize(post_preg[,i])
}

################################################################################
######################### TRANSFORM FACTOR TO NUMERIC ##########################
################################################################################

#Select exposures annotated as factors
factors <- rownames(select_df[which(select_df$Type=="Factor"),])

#Relevel factors to get the correct reference category
post_preg[,"hs_dmdtp_cdich_None"] <- factor(post_preg[,"hs_dmdtp_cdich_None"], levels=c("Undetected","Detected"))
post_preg[,"hs_smk_parents_None"] <- factor(post_preg[,"hs_smk_parents_None"], levels=c("neither","one","both"))
post_preg[,"e3_asmokyn_p_None"] <- factor(post_preg[,"e3_asmokyn_p_None"], levels=c("no","yes"))

#Transform factors to numeric
post_preg_num <- data.frame(sapply(post_preg, as.numeric, MARGIN=2))
rownames(post_preg_num) <- rownames(post_preg)

#Save the final data frame
exps <- post_preg_num
save(exps,file="/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/exps.Rdata")

########################################################################
######################### EXPOSOME HISTOGRAMS ##########################
########################################################################

hist_expo <- function(data){
  a <- ggplot(data = data) +
    geom_bar(aes(x=Group, fill=Group)) +
    theme(axis.text.x = element_text(size=12.5),
          axis.text=element_text(size=16),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
          strip.text.x = element_text(size=16),
          legend.position = 0,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "#F9F6F6", colour = NA)) +
    scale_fill_manual(values = c("Lifestyle" = "#FF4500", "OCs" = "#DEB887", "Phthalates" = "#CD853F",
                                 "Metals" = "#8FBC8F", "Phenols" = "#EE82EE", "Indoor air"="#9ACD32",
                                 "Building Environment" = "#B22222", "Tobacco Smoke" = "#696969",
                                 "Social and economic capital" = "#000080", "PFASs" = "#B0C4DE", "PBDEs" = "#20B2AA", 
                                 "OP Pesticides" = "#008080", "Natural Spaces" = "#228B22",
                                 "Meteorological" = "#FFD700", "Essential minerals" = "#FFB90F",
                                 "Air Pollution" = "#7B68EE", "Traffic" = "#F4A460", "Others"= "#8B8682")) +
    labs(x="Exposure family", y="Number of exposures", title="Early-life exposome")+
    facet_grid(~Period) +
    coord_flip()
  
  png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Data_description/Hist_Exposome_horitzontal.png"), 
      res = 1200, width = 7.5, height = 6, units = "in")
  par(mar=c(5,5,2,2))
  print(a)
  dev.off()
}

hist_expo(select_df)

########################################################################
######################### EXPOSOME DESCRIPTION #########################
########################################################################

data <- select_df[,c("Label.for.tables", "Group", "Period", "Units", 
                     "Transformation", "description")]

colnames(data) <- c("Exposure","Exposure Family","Period","Units","Transformation","Description")
writexl::write_xlsx(data, "/home/isglobal.lan/ncarreras/homews/TFM/Data_tables/Expos_list.xlsx")
