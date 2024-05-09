## Copyright © 2024 FUNDACIÓN PRIVADA INSTITUTO DE SALUD GLOBAL BARCELONA. All rights reserved.

##Computation of X-Ra in 12 different tissues with cancer 
#and constitutional biopsies 

#Merged data with clinical data is stored in phcancer object and 
#TCGA_XRa.txt (TableS1 of the manuscript).

#### libraries
library(chrXRa)
library(curatedTCGAData)
library(TCGAbiolinks)
library(ggplot2)


#methylation data from the TCGA were downloaded and CpGs from chr X filtered. 
#Data was stored in the metTCGA object.
#run only once##############################################################

load("./Data/metTCGA.RData")

#case control status of samples
splitcaco <- lapply(met0, function(x) 
  list(cancer=grep("01A", colnames(x)), 
       healthy=grep("11A", colnames(x))))


#select cancer data with >5 normal tissue samples (12 slected tissues in total)
ln <- lapply(splitcaco, function(x)  length(x$healthy))
sel <- unlist(ln)>5 

cancertypes <- names(met0)
cancersel <- cancertypes[sel]


##clinical variables were downloaded for all tissues
query <- GDCquery(project = paste("TCGA-",cancersel, sep=""), 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")


GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)

#Clinical data by tissue
namesclin <- grep("patient",names(clinical.BCRtab.all))
namesclin <- names(clinical.BCRtab.all)[namesclin]

#Fomat clinical data into a single data.frame survInfo

survInfo <- lapply(namesclin, function(x){ 
  clincan <- clinical.BCRtab.all[[x]][-c(1:2),]
  times <- sapply(1:nrow(clincan), function(i) sum(as.numeric(clincan[i,"death_days_to"]),as.numeric(clincan[i, "last_contact_days_to"]), na.rm=TRUE))
  gender <- clincan$gender
  event <- as.numeric(clincan$vital_status=="Dead")
  age <- -as.numeric(clincan$birth_days_to)/365
  admin.disease_code <- rep(strsplit(x, "_")[[1]][3], length=nrow(clincan))
  subnames <- clincan$bcr_patient_barcode
  out <- data.frame(times, event, gender, age, admin.disease_code)  
  rownames(out) <- subnames
  out
  })

survInfo <- do.call(rbind, survInfo)

#call X-Ra in each cancer study and merge with clinical data
phcancer <- lapply(cancersel, function(x){

  #split methylation data into tumor and healthy samples
  met_cancer <- met0[[x]][,grep("01A", colnames(met0[[x]]))]
  met_control <- met0[[x]][,grep("11A", colnames(met0[[x]]))]
 
  colnames(met_cancer) <- substr(colnames(met_cancer), 1, 12)
  colnames(met_control) <- substr(colnames(met_control), 1, 12)
  
  #match with clinical data
  selcancer <- survInfo$admin.disease_code==tolower(x)
  clindat <- survInfo[selcancer, ]

  #methlyation data for female cancer samples
  nmssubs <- intersect(colnames(met_cancer), rownames(clindat))
  met_cancer <- met_cancer[,nmssubs]
  ph_cancer <- clindat[nmssubs,]
  fnms <- rownames(ph_cancer)[which(ph_cancer$gender=="FEMALE")]
  met_cancer <- met_cancer[,fnms]
  met_cancer <- as.matrix(met_cancer)
  ph_cancer <- ph_cancer[fnms,]

  cpgidtcga <- rownames(met_cancer)
  
  #proportion of 2-allele demethylation in XCI cpgs
  ph_cancer$xra <- XRa(t(met_cancer), cpgidtcga)
  
  #proportion of 2-allele demethylation in escapees
  ph_cancer$xesc <- Xesc(t(met_cancer), cpgidtcga)
  
  #proportion of 1-allele methylation in escapees
  ph_cancer$xescinterm <- Xesc(t(met_cancer), cpgidtcga, threshold0=0.4, threshold=0.6)
  
  
  #methlyation data for female normal tissue
  nmssubs <- intersect(colnames(met_control), rownames(clindat))
  met_control <- met_control[,nmssubs]
  ph_control <- clindat[nmssubs,]
  fnms <- rownames(ph_control)[which(ph_control$gender=="FEMALE")]
  met_control <- met_control[,fnms]
  met_control <- as.matrix(met_control)
  ph_control <- ph_control[fnms,]
  
  cpgidtcga <- rownames(met_control)
  
  #proportion of 2-allele demethylation in XCI cpgs
  ph_control$xra <- XRa(t(met_control), cpgidtcga)
  
  #proportion of 2-allele demethylation in escapees
  ph_control$xesc <- Xesc(t(met_control), cpgidtcga)

  #proportion of 1-allele methylation in escapees
  ph_control$xescinterm <- Xesc(t(met_control), cpgidtcga, threshold0=0.4, threshold=0.6)
  
  #bind cancer and normal tissue data
  ph <- rbind(ph_cancer, ph_control)
  ID <- c(rownames(ph_cancer), rownames(ph_control))
  nmage <- grep("age", names(ph))
  age <- as.numeric(ph[,"age"])
  times <- as.numeric(ph[,"times"]) 
  event <- as.numeric(ph[,"event"]) 
  xra <- as.numeric(ph$xra)
  xesc <- as.numeric(ph$xesc)
  xescinterm <- as.numeric(ph$xescinterm)
  
  caco <- c(rep(1, nrow(ph_cancer)), rep(0, nrow(ph_control)))
  
  data.frame(ID, xra, xesc, xescinterm, caco, age, times, event)
  
})

names(phcancer) <- cancersel

save(phcancer, file="./Data/phcancer.RData")
write.table(do.call(rbind, phcancer), file="./Data/TCGA_XRa.txt", sep="\t", quote=FALSE)


#2-allele demethylation in XCI CpGs(X-Ra) and CpGs in escapees in normal tissue (Figure 2C)
##########################################################################################

load("./Data/phcancer.RData")

cpglevels <- lapply(cancersel, function(cc){
  
  
  met_cancer <- met0[[cc]][,grep("01A", colnames(met0[[cc]]))]
  met_control <- met0[[cc]][,grep("11A", colnames(met0[[cc]]))]
  
  colnames(met_control) <- substr(colnames(met_control), 1, 12)
  colnames(met_cancer) <- substr(colnames(met_cancer), 1, 12)
  
  met_control <- met_control[,colnames(met_control)%in%phcancer[[cc]]$ID]
  met_cancer <- met_cancer[,colnames(met_cancer)%in%phcancer[[cc]]$ID]
  
  
  list(eg_control=EGlevels(t(met_control), rownames(met_control)),
       ig_control=IGlevels(t(met_control), rownames(met_control)),
       eg_cancer=EGlevels(t(met_cancer), rownames(met_cancer)),
       ig_cancer=IGlevels(t(met_cancer), rownames(met_cancer)))
})

names(cpglevels) <- cancersel


ig <- unlist(lapply(cancersel, function(cc){ cpglevels[[cc]]$ig_control}))
eg <- unlist(lapply(cancersel, function(cc){ cpglevels[[cc]]$eg_control}))


# Data for Plot
metlevels2 <- rep(NA, length(ig))
metlevels2[1:length(eg)] <- eg
met_data <- data.frame(metlevels = ig, metlevels2 = metlevels2)

table(unlist(lapply(phcancer, function(x) x$caco)))

# Create the density plot
plot <- ggplot(data = met_data)

plot <- plot + labs(title = "Chr X methlylation in 446 samples \n (12 undiseased tissues)", 
                    x = "CpG methylation levels", y = "Density") + theme_bw()



plot <- plot + geom_density(aes(x = metlevels2), color = "blue", size = 0.75) + xlim(-0.05, 1.05) + ylim(0, 5)+
  scale_color_manual(values = c("CpGs in escapees" = "blue"),
                     labels = c("CpGs in escapees")) 



d <- ggplot_build(plot)$data[[1]]
plot <- plot + geom_area(data = subset(d, x >= -0.2 & x <= 0.2),
                         aes(x = x, y = y, fill = "CpGs in escapees"),
                         alpha = 0.3) +   
  scale_fill_manual(values = c("CpGs in escapees" = "grey"),
                    labels = c("CpGs in escapees")) 


plot <- plot + geom_density(aes(x = metlevels), color = "red", size = 1.25) 
plot <- plot + geom_vline(xintercept = 0.2, linetype = "dashed")


d <- ggplot_build(plot)$data[[3]]
plot <- plot + geom_area(data = subset(d, x >= 0 & x <= 0.2),
                         aes(x = x, y = y, fill = "CpGs in inactives"),
                         alpha = 0.4)


plot <- plot +   scale_fill_manual(values = c("CpGs in inactives" = "red"),
                    labels = c("CpGs in inactives")) +
  
  
  guides(fill = guide_legend(override.aes = list(color = "red")),
         color = guide_legend(override.aes = list(fill = "grey"))) +
  theme_minimal() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))

ggsave("./Figures/Figure2C.png", plot, width = 8, height = 6, dpi = 300, bg="white")


#2-allele demethylation in XCI CpGs(X-Ra) and 1-allele methylation in CpGs of escapees 
#in tumors (Figure 3A)
######################################################################################


ig <- unlist(lapply(cancersel, function(cc){ cpglevels[[cc]]$ig_cancer}))
eg <- unlist(lapply(cancersel, function(cc){ cpglevels[[cc]]$eg_cancer}))

ln <- lapply(splitcaco, function(x)  length(x$cancer))


table(unlist(lapply(phcancer, function(x) x$caco)))


# Data for Plot
metlevels2 <- rep(NA, length(ig))
metlevels2[1:length(eg)] <- eg
met_data <- data.frame(metlevels = ig, metlevels2 = metlevels2)

sig <- sort(eg[!is.na(eg)])
f <- density(sig)
froot <- abs(f$x-0.4)
xx <- f$x[which(froot==min(froot))]
yy <- f$y[which(froot==min(froot))]

lines_data1 <- data.frame(x1 = xx, x2 = xx, y1 = 0, y2 = yy)

sig <- sort(eg[!is.na(eg)])
f <- density(sig)
froot <- abs(f$x-0.6)
xx <- f$x[which(froot==min(froot))]
yy <- f$y[which(froot==min(froot))]

lines_data2 <- data.frame(x1 = xx, x2 = xx, y1 = 0, y2 = yy)


sig <- sort(ig[!is.na(ig)])
f <- density(sig)
froot <- abs(f$x-0.2)
xx <- f$x[which(froot==min(froot))]
yy <- f$y[which(froot==min(froot))]

lines_data3 <- data.frame(x1 = xx, x2 = xx, y1 = 0, y2 = yy)

# Create the density plot
plot <- ggplot(data = met_data)

plot <- plot + labs(title = "Chr X methlylation in 3,307 samples \n (12 tumours)", 
                    x = "CpG methylation levels", y = "Density") + theme_bw()



plot <- plot + geom_density(aes(x = metlevels2), color = "blue", size = 0.75) + xlim(-0.05, 1.05) + ylim(0, 5)+
  scale_color_manual(values = c("CpGs in escapees" = "blue"),
                     labels = c("CpGs in escapees")) 



d <- ggplot_build(plot)$data[[1]]
plot <- plot + geom_area(data = subset(d, x >= -0.2 & x <= 0.2),
                         aes(x = x, y = y, fill = "CpGs in escapees"),
                         alpha = 0.3) 

d2 <- ggplot_build(plot)$data[[1]]

plot <- plot + geom_area(data = subset(d, x >= 0.4 & x <= 0.6),
                         aes(x = x, y = y, fill = "CpGs2"),
                         alpha = 0.3) 

plot <- plot + geom_density(aes(x = metlevels), color = "red", size = 1.25) 


d3 <- ggplot_build(plot)$data[[4]]
plot <- plot + geom_area(data = subset(d3, x >= 0 & x <= 0.2),
                         aes(x = x, y = y, fill = "CpGs in inactives"),
                         alpha = 0.4)


plot <- plot +  scale_fill_manual(values = c("CpGs in escapees" = "grey", "CpGs2"="blue", "CpGs in inactives" = "red"),
                    labels = c("CpGs in escapees", "CpGs2", "CpGs in inactives"))+ 
  theme_minimal() +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14))

plot <- plot +  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_data1, 
                             linetype = "dashed")
plot <- plot +  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_data2, 
                             linetype = "dashed")
plot <- plot +  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_data3, 
                             linetype = "dashed")


ggsave("./Figures/Figure3A.png", plot, width = 8, height = 6, dpi = 300, bg="white")

