plot_samples <- function(GRSet, epi){
  ## 1.2. Create the GenomicRatioSet
  case_samples <- which(colnames(GRSet) == epi$sample)
  control_samples  <- which(GRSet$status == "control")
  methy <- GRSet[,c(case_samples, control_samples)]
  ## 1.3. Data frame with epimutation result
  dmr <- epi[1,c("sample", "chromosome", "start", "end", "cpg_ids")]
  ## 1.4. Plot control and case samples
  dev.new(width = 1080, height = 1350, unit = "px")
  epimutacions::plot_epimutations(dmr, methy)
  ggplot2::ggsave(paste0("vignette/fig/", dmr$sample,"_",dmr$chromosome,"-", dmr$start,"-", dmr$end,".png"))
}
#1. GSE51032 dataset
## 1.1. Load data
load("C:/Users/nla94/Desktop/GSE51032.rda")
load("hpc/outputs/GSE51032-TPR/outputfile-TPR-GSE51032-20-1.rda")
##1.2. Create plot
plot_samples(GSE51032, results[1,])

#2. GSE111629 dataset
## 2.1. Load data
load("C:/Users/nla94/Desktop/GSE111629.rda")
load("hpc/outputs/GSE111629-TPR/outputfile-TPR-GSE111629-20-1.rda")
##2.2. Create plot
plot_samples(GSE111629, results[1,])
plot_samples(GSE111629, results[2,])
plot_samples(GSE111629, results[3,])
plot_samples(GSE111629, results[4,])
plot_samples(GSE111629, results[5,])
