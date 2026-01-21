library(ggplot2)
library(ggthemes)
library(gtable)
library(grid)

library(dplyr)
library(stringr)
library(readr)





benchfiles <- list.files("csv_results/", pattern = "*SVD_H*", full.names = TRUE) 

# c("win", "macOS", "linux-hpc")

finalBenchmark <- do.call(rbind, sapply(benchfiles, read_table, simplify = FALSE)) %>%
    filter(system == "macOS" & !is.na(mean) & date> "2025-06-1" & rows<2500) %>%
    mutate(msize = paste0(rows,"x",cols))


#   Plot només multiple Cores - SVD by Blocks
# --------------------------------------------


p <- finalBenchmark %>%
    filter(( funct !="R" & funct != "BDSM-hdf5-full") & nCores > 1) %>%
    group_by( funct, system, rows, cols, expr, nCores) %>%
    dplyr::summarise( mean_g = mean(mean, na.rm = TRUE)) %>%
    mutate(levels = sub(".*_(q[0-9])_.*", "\\1", expr),
           blocs = sub(".*_(k[0-9]{1,2})$", "\\1", expr)) %>%
    ggplot() + 
    aes( x = cols, y= mean_g) +
    geom_line(aes(color = blocs, linetype = levels ), size = 1 ) +
    facet_grid( rows ~ nCores, scale = "free") +
    ylab("mean execution time (seconds)") +
    xlab("number of columns")  + 
    ggplot2::theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"),
          axis.title = element_text(size=19,face="bold"),
          axis.text =  element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=17, face = "bold"))


p <- p + scale_y_continuous(sec.axis = sec_axis(~ . , name = "number of rows", breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = "number of cores", breaks = NULL, labels = NULL)) 



png(paste0("plots/FinalPlotsPaper/svd_BigDataStatMeth_kq.png"),   width = 1024, height = 790)
    print(p) 
dev.off()
