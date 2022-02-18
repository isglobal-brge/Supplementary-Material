#######################################
## COMPARE CpGs FROM BLOOD AND HEART ##
#######################################

#Load inv_methy from HELIX
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/metanalysis/inv_methy.Rdata")
helix <- inv_methy
helix_sig <- helix[which(helix$p.adj<0.05),]

#Load inv_methy from Heart tissue
load("/home/isglobal.lan/ncarreras/data/NataliaCarreras/paper/Differential_analysis_inv/Results_methy/heart_Carlos/inv_methy_heart.Rdata")
heart <- inv_methy
heart_sig <- heart[which(heart$P.Value<0.05),]

#Compare the CpG sites
nrow(helix)
nrow(heart)
length(intersect(helix$CpG,heart$CpG))

nrow(helix_sig)
nrow(heart_sig)
cpgs_equal <- intersect(helix_sig$CpG,heart_sig$CpG)
length(cpgs_equal)

#Look at the CpG sites that intersect (whether the direction is the same)
heart[cpgs_equal,c(1,4,13,14)]

rownames(helix) <- helix$CpG

helix[cpgs_equal,c(1,4,5,10)]

#Create a table for comparing both
compare <- cbind(helix[cpgs_equal,c(1,4,5,10)],heart[cpgs_equal,c(1,4)])
compare$direction <- (compare$TE<0)==(compare$logFC<0)