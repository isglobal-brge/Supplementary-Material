########################
## PLOTS HITS INV*EXP ##
########################

## Load libraries

library(ggplot2)
library(meta)
library(ggpubr)
library(rexposome)

## Load datasets

load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy16.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy17.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/exps.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")
names(inversionGR) <- inversionGR$Cytogenetic.location
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/metanalysis/inter_methy.Rdata")
inter_methy <- inter_methy[order(inter_methy$p.adj),]

#Select samples that intersect between exposome and methylome datasets
intersect_samples<-intersect(colnames(methy8),rownames(exps))
methy8 <- methy8[,intersect_samples]
methy16 <- methy16[,intersect_samples]
methy17 <- methy17[,intersect_samples]

exps <- exps[intersect_samples,]


#Function for the factor variables
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Create the interaction plot and the forest plot
plots_inter <- function(i,inter_plot=T,forest_plot=T, factor=F, levels=NULL){
  cpg <- inter_methy[i,"CpG"]
  expo <- inter_methy[i,"Exposure"]
  expo_abrev <- inter_methy[i,"Exposure_abrev"]
  inv <- as.character(inter_methy[i,"Inversion"])
  inter <- as.character(interaction(expo, inv))
  
  
  if (cpg %in% rownames(methy8)){
    dataset <- methy8
  } else if (cpg %in% rownames(methy16)){
    dataset <- methy16
  } else {
    dataset <- methy17
  }
  
  ######################
  ## Interaction plot ##
  ######################
  
  if(inter_plot){
    if (factor){
      covs <- data.frame(Inv_numeric=as.numeric(dataset[[inversionGR[inv,]$scoreInvHap.name]])-1,
                               Exposure=factor(exps[,expo]),
                               CpG=c(assay(dataset[cpg,])),
                               Inversion=dataset[[inversionGR[inv,]$scoreInvHap.name]])
      levels(covs$Exposure) <- levels
      covs_se <- summarySE(covs, measurevar="CpG", groupvars=c("Exposure","Inversion"))
      plotA <- ggplot(covs_se, aes(x=Exposure, y=CpG, col=Inversion))+
        geom_errorbar(aes(ymin=CpG-se, ymax=CpG+se), width=.1)+
        geom_line(aes(group=Inversion))+
        geom_point()
      
    } else {
      covs <- data.frame(Inv_numeric=as.numeric(dataset[[inversionGR[inv,]$scoreInvHap.name]])-1,
                               Exposure=exps[,expo],
                               CpG=c(assay(dataset[cpg,])),
                               Inversion=dataset[[inversionGR[inv,]$scoreInvHap.name]])
      plotA <- ggplot(covs, aes(x=Exposure, y=CpG, col=Inversion))+
        geom_point(size=1)+
        geom_smooth(method = "lm", se = FALSE)
    }
    
    plotB <- plotA +
      labs(x=expo_abrev, y=paste0("Methylation"), col="Genotype")+
      facet_wrap(~Inversion, ncol=1, strip.position = "right") +
      scale_color_manual(values = c("N/N" = "#00AFBB", "N/I" = "#E7B800", "I/I" = "#FC4E07")) +
      theme(axis.title.x = element_text(size = 17, face="bold",),
            axis.title.y = element_text(size = 17, face="bold",),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            strip.text.y = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position= "top",
            legend.title=element_text(size=15), 
            legend.text=element_text(size=15))
    png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Specific_cases/Hits_inv_exp/Hit",i,"_InterPlot_nores.png"),
        res = 1200, width = 5, height = 5, units = "in")
    par(mar=c(5,5,2,2))
    print(plotB)
    dev.off()
  }
  
  #################
  ## Forest plot ##
  #################
  
  if(forest_plot){
    cohorts <- c("BIB", "EDEN", "KANC", "MOBA", "RHEA", "SAB")
    
    #Load results by cohorts
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
    
    subset_function <- function(coh,inter,cpg){
      it_all <- allit[[coh]]
      it <- it_all[which(it_all$Interaction==inter&
                           it_all$CpG==cpg),]
    }
    
    subset <- t(sapply(cohorts,subset_function,inter=inter,cpg=cpg))
    
    met <- metagen(TE=do.call(rbind,subset[,"LogFC"]), 
                   seTE=do.call(rbind,subset[,"SE"]), sm="MD",
                   studlab=rownames(subset), comb.fixed = FALSE)
    
    png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Specific_cases/Hits_inv_exp/Hit",i,"_ForestPlot.png"), 
        res = 1200, width = 9, height = 4.5, units = "in")
    forest.meta(met, studlab=TRUE)
    dev.off()
  }
}

#TDH
plots_inter(2, factor=T, levels=c("<6","6-9",">9"))

#GATA4
plots_inter(9)

#C8orf79
plots_inter(12, factor=T, levels=c("Neither","One","Both"))

