library(ggplot2)


graph <- function(file_name){
  
  load(paste0("result_files/3-TPR_FPR_accuracy/", file_name,".rda"))

  names <- strsplit(names(df), "_")
  dataset_name <- sub(".*_", "", file_name)

  do.call(list, lapply(seq(length(df)), function(i){
    df_region <- df[[i]]
    df_region <- reshape2::melt(df_region, id = c("method", "n"))
  
    dev.new(width = 1080, height = 1350, unit = "px")
    graph <- ggplot(df_region, aes(x = method, y = as.numeric(value), fill = factor(method))) +
                    geom_bar(stat="identity") + 
                    facet_grid(variable ~ as.numeric(n),
                             scales = "free_y")
    style <- graph + theme_minimal() + 
             theme(axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank(), 
                   axis.title.x = element_blank(),
                   strip.background = element_rect(color="Black", size=0.2, linetype="solid")) +
             ggtitle(paste0(dataset_name,": ", names[[i]][1], ":", names[[i]][2], " - ", names[[i]][3])) +
             theme(plot.title = element_text(hjust = 0.5)) + 
             scale_fill_manual("Method", values = c("#8F9CDA", "#DE9BB0", "#CE6D96", "#3972AA", "8E95D9", "#C5A18C")) + 
             ylab("%") 
    style
    ggsave(paste0( "result_files/4-Graph/", dataset_name, "_", names[[i]][1], "_", names[[i]][2], "_", names[[i]][3], ".png"))
  }))
}


graph("TPR_FPR_accuracy_GSE111629")
graph("TPR_FPR_accuracy_GSE51032")

           


