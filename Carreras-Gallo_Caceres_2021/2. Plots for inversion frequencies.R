####################################################################
########## FREQUENCIES OF INVERSIONS IN DIFFERENT COHORTS ##########
####################################################################

#Load libraries
library(SummarizedExperiment)
library(ggplot2)

#Load data
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/Methy_final.Rdata")
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/inversionGR.rda")


#Function to create plots
plot_frec<-function(inv, x_lab){
  tit<-paste("Frecuency of inv-",inversionGR[inv,]$Cytogenetic.location," depending on ",tolower(x_lab),sep="")
  plot <- ggplot(data = as.data.frame(pData(methy_final)), aes_string(x=tolower(x_lab), fill=inv)) +
    geom_bar(position="fill", alpha=0.7) +
    labs(x=x_lab, y="Percentage",title=tit) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 13),
          plot.title = element_text(hjust=0.5,size = 16, face = "bold"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_fill_manual("Genotype", values = c("N/N" = "#00AFBB", "N/I" = "#E7B800", "I/I" = "#FC4E07")) +
    scale_y_continuous(labels = scales::percent)
  path <- paste("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Data_description/Plots_frec_inversions/freq_",tolower(x_lab),inv,".png",sep="")
  png(path, res = 1200, width = 6, height = 3.25, units = "in")
  print(plot)
  dev.off()
}

plot_frec("inv8_001", x_lab="Cohort")
plot_frec("inv16_009", x_lab="Cohort")
plot_frec("inv17_007", x_lab="Cohort")

plot_frec("inv8_001", x_lab="Sex")
plot_frec("inv16_009", x_lab="Sex")
plot_frec("inv17_007", x_lab="Sex")