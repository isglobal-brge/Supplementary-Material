# 1. Create a list with all simulation results and the validated epimutations

create_result_list<- function(folder_name, rate,  variable, case_sample = NULL, save_name, epi_validated = NULL){
  RES.n <- list()
  for(n in seq(20, 100, 10)){
    RES.i <- list()
    for(i in 1:100){
      load(paste0("hpc/outputs/", folder_name,"/outputfile-", rate, "-", variable, "-", n, "-", i,".rda"))
      RES.i[[length(RES.i)+1]] <- results
    }
    res <- as.data.frame(data.table::rbindlist(RES.i))
    if(!is.null(case_sample)){
    res <- res[res$sample %in% case_sample,]
    }
    if(!is.null(epi_validated)){
      query <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(res), seqnames.field = "chromosome", start.field = "start", end.field= "end")
      subject  <- epi_validated
      hits <- GenomicRanges::findOverlaps(query, subject, minoverlap = 0, maxgap = 1000,  type ="equal")
      overlaps <- GenomicRanges::pintersect(query[S4Vectors::queryHits(hits)], subject[S4Vectors::subjectHits(hits)])
      percentOverlap <- BiocGenerics::width(overlaps) / BiocGenerics::width(subject[S4Vectors::subjectHits(hits)])
      res$accuracy[S4Vectors::queryHits(hits)] <- percentOverlap 
      validated <- as.data.frame(subject)
      res$chromosome_validated[S4Vectors::queryHits(hits)] <- as.character(validated[S4Vectors::subjectHits(hits),1])
      res$start_validated[S4Vectors::queryHits(hits)] <- validated[S4Vectors::subjectHits(hits),2]
      res$end_validated[S4Vectors::queryHits(hits)] <- validated[S4Vectors::subjectHits(hits),3]
      res$sz_validated[S4Vectors::queryHits(hits)] <- validated[S4Vectors::subjectHits(hits),4]
      #res <- res[,c(1:6, 18:21, 7:17)] 
      res <- res[,c(1:6, 19:22, 7:18)] 
    }
    RES.n[[length(RES.n)+1]] <- res 
    
  }
  
  names(RES.n) <- paste0("n", seq(20, 100, 10))
  
  result_table <- RES.n
  
  save(result_table, file = paste0("result_files/1-Result_list/", save_name,".rda"))
}


# 2. Function to calculate the TPR and accuracy
TPR <- function(file_name, epi_validated){
  #2.1. Load the data file containing the results list  
  load(paste0("result_files/1-Result_list/", file_name, ".rda"))
  #2.2. Set the parameters
  control_n <- seq(20, 100, 10)
  #methods <- c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta")
  methods <- c("beta", "wbeta", "IQR")
  #2.3. Calculate TPR and accuracy
  regions <-  do.call(list, lapply(seq(nrow(epi_validated)), function(i){
    TPR <- do.call(list, lapply(seq(length(result_table)), function(n){
      rst <-  do.call(rbind, lapply(seq_len(length(methods)), function(j){
        df_n <- result_table[[n]]
        keep_region <- which(df_n$chromosome_validated == epi_validated$seqnames[i] & df_n$start_validated == epi_validated$start[i])
        region <- df_n[keep_region,]
        keep_method <- grepl(methods[j],unlist(region[,"epi_id"]))
        out <- region[keep_method,]
        percent <- round(mean(na.omit(out[,"accuracy"])), digits = 2)
          TPR <- round(length(which(out[,"chromosome"] != 0))/nrow(df_n), digits = 2)
        table <- c(methods[j], as.numeric(control_n[n]), as.numeric(TPR), round(as.numeric(percent), digits = 3))
      }))
      rst <- as.data.frame(rst)
      rownames(rst)
      colnames(rst) <- c("method", "n", "TPR", "accuracy")
      rst
    }))
    names(TPR) <- paste0("n", seq(from = 20, to = 100, by = 10))
    TPR
  }))
  names(regions) <- paste0(epi_validated$seqnames, "_", epi_validated$start, "_", epi_validated$end)
  TPR <- regions
  #2.4. Save results
  save(TPR, file = paste0("result_files/2-TPR_accuracy/", file_name, ".rda"))
}

#3. Function to calculate the FPR 
FPR <- function(file_name){
  #3.1. Load the data file containing the results list
  load(paste0("result_files/1-Result_list/", file_name, ".rda"))
  #3.2. Set the parameters  
  control_n <- seq(20, 100, 10)
  #methods <- c("manova", "mlm", "mahdistmcd", "isoforest", "quantile", "beta")
  methods <- c("beta", "wbeta", "IQR")
  #3.3. Calculate the FPR 
  FPR <- do.call(list, lapply(seq(length(result_table)), function(n){
    rst <- do.call(rbind, lapply(seq_len(length(methods)), function(i){
      results <- grepl(methods[i],unlist(result_table[[n]]["epi_id"]))
      results <- result_table[[n]][results,]
      if(all(is.na(results$cpg_n))){
        FPR <- 0
      }else{
        keep <- which(!is.na(results$cpg_n))
        results <- results[keep,]
        FPR <- nrow(unique(results))/400
      }
      
      table <- c(methods[i], as.numeric(control_n[n]), as.numeric(FPR))
      table
      
    }))
    rst <- as.data.frame(rst)
    colnames(rst) <- c("method", "n", "FPR")
    rst
  }))
  names(FPR) <-  paste0("n", seq(from = 20, to = 100, by = 10))
  
  #3.4. Save results
  save(FPR, file = paste0("C:result_files/2-FPR/", file_name, ".rda"))
}

#4. Function to calculate the table contaning the TPR, FPR and accuracy

TPR_FPR_accuracy <- function(FPR_file_name, TPR_file_name, final_file_name){
  # 4.1. Load the data files containing TPR, FPR and accuracy list 
  load(paste0("result_files/2-TPR_accuracy/", TPR_file_name, ".rda"))
  load(paste0("result_files/2-FPR/", FPR_file_name, ".rda"))
  # 4.2. Create the table
  df <- do.call(list, lapply(seq_len(length(TPR)), function(i) {
    TPR_n <- as.data.frame(data.table::rbindlist(TPR[[i]]))
    FPR <- as.data.frame(data.table::rbindlist(FPR))
    out <- merge(TPR_n, FPR, by = c("method", "n"), sort = FALSE)
    out[,3:5] <- apply(out[,3:5], 2, function(x) as.numeric(x))
    out[,3:5] <- out[,3:5] * 100
    out
  }))
  names(df) <- names(TPR)
  # 4.3. Save results  
  save(df, file = paste0("result_files/3-TPR_FPR_accuracy/",final_file_name, ".rda"))
}

# 4. Function to create the results graphics

library(ggplot2)


graph <- function(file_name){
  
  load(paste0("result_files/3-TPR_FPR_accuracy/", file_name,".rda"))
  
  names <- strsplit(names(df), "_")
  dataset_name <- sub(".*_", "", file_name)
  
  do.call(list, lapply(seq(length(df)), function(i){
    df_region <- df[[i]]
    df_region <- reshape2::melt(df_region, id = c("method", "n"))
    
    dev.new(width = 1080, height = 1350, unit = "px")
    graph <- ggplot(df_region, aes(x = as.numeric(n), y = as.numeric(value), colour = factor(method))) +
      #geom_point()+
      geom_jitter(width = 1) + geom_smooth(method = "loess", se = FALSE, size = 0.75) +
      #geom_smooth(method = lm) +
      facet_grid( rows =vars(variable)) 
      #facet_grid( variable ~ as.numeric(n))
    style <- graph + theme_bw() + theme(strip.background =element_rect(fill="white")) +
      #theme(axis.text.x = element_blank(), 
            #axis.ticks.x = element_blank(), 
            #axis.title.x = element_blank(),
            #strip.background = element_rect(color="Black", size=0.2, linetype="solid")) +
      #ggtitle(paste0(dataset_name,": ", names[[i]][1], ":", names[[i]][2], " - ", names[[i]][3])) +
      #ggtitle(dataset_name) +
      ggtitle("ramr simulations") +
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_colour_manual("Method", values = c("#E6C570", "#4F766F", "#0182B0", "#8D4925", "#F3505A", "#6C117B")) + 
      ylab("%")   + xlab("Sample size")
    style
    #ggsave(paste0( "result_files/4-Graph/", dataset_name, "_", names[[i]][1], "_", names[[i]][2], "_", names[[i]][3], ".png"))
    ggsave("result_files/4-Graph/simulations.png")
  }))
}

#4.1. Graph version 2

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

library(ggplot2)
library(gridExtra)

setwd("C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021")
file_name <- "ramr-TPR_FPR_accuracy_GSE111629"
file_name <- "TPR_FPR_accuracy_GSE111629"
load(paste0("result_files/3-TPR_FPR_accuracy/", file_name,".rda"))
df_bind <- rbind(df[[1]], df[[2]], df[[3]], df[[4]])
file_name <- "ramr-TPR_FPR_accuracy_GSE51032"
file_name <- "TPR_FPR_accuracy_GSE51032"
load(paste0("result_files/3-TPR_FPR_accuracy/", file_name,".rda"))
df <- rbind(df_bind, df[[1]])
rm("df_bind", "file_name")


df_TPR <- df[, c("method", "n", "TPR")]
df_accuracy<- df[, c("method", "n", "accuracy")]
df_FPR <- df[, c("method", "n", "FPR")]

TPR <- ggplot(df_TPR, 
              aes(x = as.numeric(n), 
                  y = as.numeric(TPR), 
                  colour = factor(method))) +
  geom_jitter(width = 1) + 
  geom_smooth(method = "loess", se = FALSE, size = 0.75) +
  ylim(0, 100)


FPR <- ggplot(df_FPR, 
              aes(x = as.numeric(n), 
                  y = as.numeric(FPR), 
                  colour = factor(method))) +
  geom_jitter(width = 0) + 
  geom_smooth(method = "loess", se = FALSE, size = 0.75) +
  ylim(0,60)


dev.new(width = 1080, height = 1350, unit = "px")


#ramr
grid.arrange(TPR + 
               theme_Publication() + 
               scale_colour_manual("Method", values = c("#386cb0","#fdb462","#7fc97f")) +
               ylab("TPR (%)") + 
               xlab("Sample size"),
             FPR + 
               geom_hline(yintercept = 5,
                          linetype="dashed",
                          color = "black", 
                          size=1) + 
               theme_Publication() + 
               scale_colour_manual("Method", values = c("#386cb0","#fdb462","#7fc97f")) +
               ylab("FPR (%)") + 
               xlab("Sample size"),nrow = 1, top=textGrob("ramr",gp=gpar(fontsize=20,font=3)))

#epimutacions
grid.arrange(TPR + 
               theme_Publication() + 
               scale_colour_manual("Method", values = c("#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")) +
               ylab("TPR (%)") + 
               xlab("Sample size"),
             FPR + 
               geom_hline(yintercept = 5,
                          linetype="dashed",
                          color = "black", 
                          size=1) + 
             theme_Publication() + 
               scale_colour_manual("Method", values = c("#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")) +
               ylab("FPR (%)") + 
               xlab("Sample size"),nrow=1, top=textGrob("epimutacions",gp=gpar(fontsize=20,font=3)))







