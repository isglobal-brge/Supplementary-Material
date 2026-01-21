library(ggplot2)
library(ggthemes)
library(gtable)
library(grid)

library(dplyr)
library(readr)


benchfiles <- list.files("csv_results", pattern = "*Mult*", full.names = TRUE, recursive = TRUE) 

# c("win", "macOS", "linux-hpc")

finalBenchmark <- do.call(rbind, sapply(benchfiles, read_table, simplify = FALSE)) %>%
    filter(!is.na(mean)) %>%
    # filter(system == "linux-hpc" & !is.na(mean)) %>%
    mutate(msize = paste0(rows,"x",cols),
           system = as.factor(system),
           blockSize = as.factor(blockSize))

#   Plot només multiple Cores - Multiplication by Blocks
# ------------------------------------------------------

data_analisis <- finalBenchmark %>%
    filter(( funct !="R" & funct != "BDSM-mem" & funct != "BDSM-mem-Block") & nCores > 1 & system == "macOS") %>%
    group_by( funct, system, rows, cols, expr, nCores, blockSize) %>%
    dplyr::summarise( mean_g = mean(mean, na.rm = TRUE)) %>%
    write.csv2("summary_data/analisis_dades_matmult_block_hdf5.csv", quote = FALSE, row.names = FALSE)

p <- finalBenchmark %>%
    filter(( funct !="R" & funct != "BDSM-mem" & funct != "BDSM-mem-Block") & nCores > 1 & system == "macOS") %>%
    group_by( funct, system, rows, cols, expr, nCores, blockSize) %>%
    dplyr::summarise( mean_g = mean(mean, na.rm = TRUE)) %>%
    # mutate(levels = sub(".*_(q[0-9])_.*", "\\1", expr),
    #        blocs = sub(".*_(k[0-9]{1,2})$", "\\1", expr)) %>%
    ggplot() + 
    aes( x = cols, y= mean_g) +
    geom_line(aes(color = blockSize ), size = 1 ) +
    facet_grid(  ~ nCores, scale = "free") +
    ylab("mean execution time (seconds)") +
    xlab("number of columns and rows")  + 
    ggplot2::theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"),
          axis.title = element_text(size=19,face="bold"),
          axis.text =  element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=17, face = "bold"))


 p + scale_x_continuous(sec.axis = sec_axis(~ . , name = "number of cores", breaks = NULL, labels = NULL)) 
    


png(paste0("plots/FinalPlotsPaper/matmult_BigDataStatMeth_BlockSize.png"), width = 1024, height = 400)
    p + scale_x_continuous(sec.axis = sec_axis(~ . , name = "number of cores", breaks = NULL, labels = NULL)) 
dev.off()


#   Plot només multiple R vs HDF5 
# ------------------------------------------------------


data_analisis <- finalBenchmark %>%
    filter(  funct != "BDSM-mem-Block" & system == "macOS" & ( blockSize != 1024 & blockSize != 2048 ) & rows<30000 & cols<30000) %>%
    group_by( funct, rows, cols, expr, nCores, blockSize) %>%
    dplyr::summarise( mean_g = mean(mean, na.rm = TRUE)) %>%
    write.csv2("summary_data/analisis_dades_matmult.csv", quote = FALSE, row.names = FALSE)



p <- finalBenchmark %>%
    filter(  funct != "BDSM-mem-Block" & system == "macOS" & ( blockSize != 1024 & blockSize != 2048 ) & rows<30000 & cols<30000) %>%
    group_by( funct, rows, cols, expr, nCores, blockSize) %>%
    dplyr::summarise( mean_g = mean(mean, na.rm = TRUE)) %>%
    ggplot() + 
    aes( x = cols, y= mean_g) +
    geom_line(aes(color = funct ), size = 1 ) +
    facet_grid( blockSize ~ nCores, scale = "free") +
    ylab("mean execution time (seconds)") +
    xlab("number of columns and rows")  + 
    ggplot2::theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"),
          axis.title = element_text(size=19,face="bold"),
          axis.text =  element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=17, face = "bold"))


p + scale_y_continuous(sec.axis = sec_axis(~ . , name = "block size", breaks = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = "number of cores", breaks = NULL, labels = NULL)) 

png(paste0("plots/FinalPlotsPaper/matmult_R_BigDataStatMeth.png"), width = 1024, height = 790)
    p + scale_x_continuous(sec.axis = sec_axis(~ . , name = "number of cores", breaks = NULL, labels = NULL)) 
dev.off()


