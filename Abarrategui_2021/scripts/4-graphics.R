#1. Load files
load.dir <- "result_files/TPR_FPR_accuracy"
file_name <- "GSE111629_chr_start_end"
load(paste0(load.dir, "/", file_name,".rda"))
#2. Set the parameters
dataset <- "GSE111629"
chr <- "chr5"
start <- "67583971"
end <- "67584381"

##2.1. Only for GSE111629 dataset
df <- df$chr5_67583971_67584381

#3. Generate the graphics
#delete <- which(is.nan(as.numeric(df$accuracy)))
#df <- df[-delete,]

library(ggplot2)

##3.1. graph: TPR, Accuracy and FPR graph

df <- reshape2::melt(df, id = c("method", "n"))
dev.new(width = 1080, height = 1350, unit = "px")
graph <- ggplot(df, aes(x = method, y = as.numeric(value), fill = factor(method))) +
  geom_bar(stat="identity") + 
  facet_grid(variable ~ as.numeric(n))

style <- graph + theme_minimal() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        strip.background = element_rect(color="Black", size=0.2, linetype="solid")) +
  ggtitle(paste0(dataset,": ", chr, ":", start, " - ", end)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual("Method", values = c("#8F9CDA", "#DE9BB0", "#CE6D96", "#3972AA", "8E95D9", "#C5A18C")) + ylab("%") +
  ylim(0,1) 
style

save.dir<- "result_files/Graphic/"
ggsave(paste0(save.dir,dataset, "_", chr, "_", start, "_", end, ".png"))




           


