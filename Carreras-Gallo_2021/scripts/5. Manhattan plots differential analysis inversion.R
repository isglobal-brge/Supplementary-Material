#####################################################################
## MANHATTAN PLOT FOR DIFFERENTIAL ANALYSIS DUE TO INVERSION STATE ##
#####################################################################

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
manhattan <- function(dataset,name_dataset,inv, lim_y){
  
  chr <- str_extract(inv,"[\\d]+")
  
  #Store the position to locate the points from the Manhattan plot
  #For methylome, use the CpG site position
  if (name_dataset=="Methylome"){
    positions <- sapply(strsplit(as.character(dataset$Location),":"),"[",2)
    dataset$position <- as.numeric(positions)
    title_plot <- paste0("Differentially methylated CpG sites for ",inv)
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
    title_plot <- paste0("Differentially expressed genes for ",inv)
  }
  
  #Subset the rows of interest
  mydata <- dataset[which(dataset$Inversion==inv),]
  
  #Define the plot title and path
  path_plot_fam <- paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_",
                          tolower(substr(name_dataset,1,5)),
                         "/Summary_plots/Manhattan_inv_",
                          tolower(substr(name_dataset,1,5)),
                          substr(inv,1,2),
                          "_fam.png")
  path_plot_ggviz <- paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_",
                            tolower(substr(name_dataset,1,5)),
                            "/Summary_plots/Manhattan_inv_",
                            tolower(substr(name_dataset,1,5)),
                            substr(inv,1,2),
                            "_ggviz.png")
  
  ############
  ## ggplot ##
  ############
  
  e <- ggplot(mydata, aes(x = position, y = -log10(pval))) + 
    geom_point(col="grey") +
    geom_point(data=mydata[which(mydata$pval<(0.05/nrow(dataset))),],aes(x = position, 
                                                        y = -log10(pval),), col="orange") +
    geom_text(data=mydata[which(mydata$pval<(0.05/nrow(dataset))),],
             aes(label=Gene_Symbol),hjust=-0.1, vjust=-0.1) +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    labs(x="Chromosome Position", y="-log10(p-value)", title=title_plot) +
    geom_hline(yintercept = -log10(0.05/nrow(dataset)), color= "red", linetype="dashed")
  
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
                      size = 13,
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

#Load results from the differential analysis
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/inv_methy.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_trans/inv_trans.Rdata")

#Call function for methylome
for (inv in c("8p23.1","16p11.2","17q21.31")){
    manhattan(inv_methy,"Methylome",inv, lim_y=250)
}

#Call function for transcriptome
#First, select only one gene symbol for each transcript
for (i in 1:nrow(inv_trans)){
  inv_trans$Gene_Symbol[i] <- strsplit(inv_trans$Gene_Symbol[i],";")[[1]][1]
}
for (inv in c("8p23.1","16p11.2","17q21.31")){
    manhattan(inv_trans,"Transcriptome",inv, lim_y=110)
}
