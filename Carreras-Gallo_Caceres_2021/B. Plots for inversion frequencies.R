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


########################
## Inversion gradient ##
########################

#Load libraries
library(ggplot2)

#Load data
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Final_datasets/methy8.Rdata")

#Define phenotypes
phenos_methy <- data.frame(cohort=methy8$cohort,
                           Inversion.8p23.1=methy8$inv8_001,
                           Inversion.16p11.2=methy8$inv16_009,
                           Inversion.17q21.31=methy8$inv17_007)

#Create a function to plot the gradient of each inversion

grad.plot <- function(inversion,gradient,direction,x,y1,y2){
  info.inv <- table(phenos_methy$cohort,phenos_methy[[inversion]])
  info.inv<- as.data.frame.matrix(info.inv)
  
  info.inv$propI <- (info.inv[,2]+info.inv[,3]*2)/(rowSums(info.inv)*2)
  info.inv$cohort <- rownames(info.inv)
  coord <- data.frame(cohort=c("RHEA","SAB","EDEN","BIB","KANC","MOBA"),
                      Latitude=c(39.51, 42.17, 46.20, 54.22, 55.34, 60.64),
                      Longitude=c(21.77, 1.91, 1.99, -2.05, 23.70, 8.41))
  
  info.inv <- merge(info.inv,coord,by="cohort")
  info.inv <- info.inv[order(info.inv[[gradient]]),]
  
  
  plot <- ggplot(data = info.inv, aes_string(x=gradient, y="propI", label="cohort")) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x) +
    geom_text(aes(label=cohort),hjust=0.5, vjust=-0.5) +
    labs(x=gradient, y="Inverted Allele proportion", title=paste0("Inversion ",substr(inversion,11,nchar(inversion))," ",direction," gradient")) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 13),
          plot.title = element_text(hjust=0.5,size = 16, face = "bold"))
  
  corr <- cor.test(info.inv$propI,info.inv[[gradient]])
  plot <- plot +  annotate("text", x = x, y = y1, label = paste0("P = ",round(corr[[3]],3))) +
    annotate("text", x = x, y = y2, label = paste0("r = ",round(corr[[4]],3)))
  png(paste0("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Data_description/Gradient_inversions/",
             gradient,"_",inversion,".png"), res = 1200, width = 6, height = 3.25, units = "in")
  print(plot)
  dev.off()
  return(cor.test(info.inv$propI,info.inv[[gradient]]))
}

corr8lat <- grad.plot(inversion="Inversion.8p23.1", gradient="Latitude", direction="south-north", x=59, y1=0.7, y2=0.68)
corr16lat <- grad.plot(inversion="Inversion.16p11.2", gradient="Latitude", direction="south-north", x=41.5, y1=0.47, y2=0.45)
corr17lat <- grad.plot(inversion="Inversion.17q21.31", gradient="Latitude", direction="south-north", x=59, y1=0.33, y2=0.31)

corr8lon <- grad.plot(inversion="Inversion.8p23.1", gradient="Longitude", direction="east-west")
corr16lon <- grad.plot(inversion="Inversion.16p11.2", gradient="Longitude", direction="east-west")
corr17lon <- grad.plot(inversion="Inversion.17q21.31", gradient="Longitude", direction="east-west")
