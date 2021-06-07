########################################################################
## MANHATTAN PLOT FOR DIFFERENTIAL ANALYSIS DUE TO INVERSION*EXPOSOME ##
########################################################################

library(SummarizedExperiment)
library(ggplot2)
library(Gviz)
library(stringr)


#Define the start and the end of the inversion region
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")
coordinates <- as.data.frame(inversionGR[c("inv8_001","inv16_009","inv17_007"),])
rownames(coordinates) <- coordinates$Cytogenetic.location

#Define the start and the end of the inversion region +/- 1Mb
inv_regions <- as.data.frame(inversionGR[c("inv8_001","inv16_009","inv17_007"),])
inv_regions$start <- inv_regions$start-1000000
inv_regions$end <- inv_regions$end+1000000
rownames(inv_regions) <- inv_regions$Cytogenetic.location



#Create the function to plot the Manhattan plots (with ggplot2 and Gviz)
manhattan <- function(dataset,name_dataset,inv, analysis, coh="", lim_y){
  
  chr <- str_extract(inv,"[\\d]+")
  
  #Store the position to locate the points from the Manhattan plot
  #For methylome, use the CpG site position
  if (name_dataset=="Methylome"){
    positions <- sapply(strsplit(as.character(dataset$Location),":"),"[",2)
    dataset$position <- as.numeric(positions)
  }
  #For transcriptome, use the middle point of the transcript
  if (name_dataset=="Transcriptome"){
    trans <- function(s) {
      cc <- strsplit(s,":")
      range <- unlist(cc)[2]
      cc2 <- strsplit(range,"-")
      first <- as.numeric(unlist(cc2)[1])
      last <- as.numeric(unlist(cc2)[2])
      return(mean(first,last))
    }
    positions <- sapply(as.character(dataset$Location),trans)
    dataset$position <- as.numeric(positions)
  }
  
  #Subset the interactions of interest
  mydata <- dataset[which(dataset$Inversion==inv),]
  
  #Define the plot title and path
  title_plot <- paste0("Significant interactions for ",inv," - ",name_dataset)
  path_plot_fam <- paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_",
                          tolower(substr(name_dataset,1,5)),
                          "_interaction/",analysis,"/Summary_plots/Manhattan_inter_",
                          tolower(substr(name_dataset,1,5)),
                          substr(inv,1,2),
                          "_fam",coh,".png")
  path_plot_ggviz <- paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_",
                            tolower(substr(name_dataset,1,5)),
                            "_interaction/",analysis,"/Summary_plots/Manhattan_inter_",
                            tolower(substr(name_dataset,1,5)),
                            substr(inv,1,2),
                            "_ggviz",coh,".png")
  
  ############
  ## ggplot ##
  ############
  
  e <- ggplot(mydata, aes(x = position, y = -log10(pval))) + 
    geom_point(col="grey") +
    geom_point(size=2, data=mydata[which(mydata$pval<(0.05/nrow(dataset))),],aes(x = position, 
                                                                                 y = -log10(pval), 
                                                                                 col= Family)) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    labs(x="Chromosome Position", y="-log10(p-value)", title=title_plot) +
    ylim(0,lim_y) +
    geom_hline(yintercept = -log10(0.05/nrow(dataset)), color= "red", linetype="dashed") +
    scale_color_manual(values = c("Lifestyle" = "#FF4500", "OCs" = "#DEB887", "Phthalates" = "#CD853F",
                                  "Metals" = "#8FBC8F", "Phenols" = "#EE82EE", "Indoor air"="#9ACD32",
                                  "Built Environment" = "#B22222", "Tobacco Smoke" = "#696969",
                                  "Socio-eco capital" = "#000080", "PFASs" = "#B0C4DE", "PBDEs" = "#20B2AA", 
                                  "OP Pesticides" = "#008080", "Natural Spaces" = "#228B22",
                                  "Meteorological" = "#FFD700", "Essential minerals" = "#FFB90F",
                                  "Air Pollution" = "#7B68EE", "Traffic" = "#F4A460", "Noise"= "#8B8682"),
                       limits=c("Lifestyle", "OCs", "Phthalates",
                                "Metals", "Phenols", "Indoor air",
                                "Built Environment", "Tobacco Smoke",
                                "Socio-eco capital", "PFASs", "PBDEs", 
                                "OP Pesticides", "Natural Spaces",
                                "Meteorological", "Essential minerals",
                                "Air Pollution", "Traffic", "Noise"))
  
  ggsave(plot = e, width = 10, height = 5, units = "in",
         filename = path_plot_fam)
  
  ##########
  ## Gviz ##
  ##########
  
  #Define the tracks
  
  #Track ideogram (chromosome)
  itrack <- IdeogramTrack(genome = "hg19", chromosome = paste0("chr",chr), cex = 1.3)
  
  #Track for the chr positions
  gtrack <- GenomeAxisTrack(cex = 1.3)
  
  #Track manhattan plot
  dtrack <- DataTrack(data = -log10(mydata$pval), start = mydata$position,
                      end = mydata$position,
                      chromosome = paste0("chr",chr), genome = "hg19",
                      name = "-log10(p-value)",
                      size = 13, ylim=c(0,lim_y),
                      background.title = "#20B2AA")
  
  #Track inversion region
  invtrack <- GeneRegionTrack(coordinates[which(coordinates$seqnames==paste0("chr",chr)),], genome = "hg19",
                              chromosome = paste0("chr",chr),
                              name = "Inv",
                              background.title = "#981218", fill="#D64F55", col="#D64F55")
  
  png(path_plot_ggviz, width=6, height=5.5, units = "in", res=300)
  plotTracks(list(itrack,gtrack,dtrack,invtrack),
             cex.title = 1, sizes=c(1.4,1.5,7,0.7))
  dev.off()
}


#META-ANALYSIS
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy.Rdata")

#Decreasing order to see the most associated points on top
inter_methy <- inter_methy[order(inter_methy$p.adj, decreasing=T),]

for (inv in c("8p23.1","16p11.2","17q21.31")){
  manhattan(inter_methy,"Methylome",inv, analysis="NN_ref", lim_y = 16)
}

#ANALYSIS BY COHORTS
plots <- function(coh){
  load(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Results_methy_interaction/NN_ref/inter_methy_",coh,".Rdata"))
  for (inv in c("8p23.1","16p11.2","17q21.31")){
    for (type in c("Mean","Variance")){
      manhattan(inter_methy,"Methylome",inv,type,coh=coh, analysis="NN_ref")
    }
  }
}

cohorts <- c("BIB", "EDEN", "KANC", "MOBA", "RHEA", "SAB")
mclapply(cohorts,plots, mc.cores=6)

