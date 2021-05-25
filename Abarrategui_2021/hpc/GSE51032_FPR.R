library(epimutacions)

#1. Load the data
n <- Sys.getenv("N")
i <- Sys.getenv("I")

print(paste0("n=", n, ", i=", i)) 

input.filename <- paste0("inputs/GSE51032",".rda")
variable <- load(file=input.filename)
load("inputs/epimutations_found_GSE51032.rda")
gr_found <- GenomicRanges::makeGRangesFromDataFrame(found_GSE51032)
load("inputs/names_GSE51032.rda")


#2. Filtering inputs
found <- which(colnames(GSE51032) %in% unique(names))
GSE51032 <- GSE51032[,-found]


#3. Functions
##3.1 Random region generation
random_regions <- function(object, span = 20000, min_cpg = 3){
  my_chr <- paste0("chr",1:22)
  fd<- as.data.frame(SummarizedExperiment::rowRanges(object))
  
  #initialise list to store chromosome sizes
  my_chr_max <- list()
  my_chr_min <- list()
  for (i in 1:length(my_chr)){
    a <- fd[fd$seqnames == my_chr[i],]
    my_chr_max[[i]] <- max(a$start)
    my_chr_min[[i]] <- min(a$start)
  }
  #checkout my_chr_size
  names(my_chr_max) <- my_chr
  names(my_chr_min) <- my_chr
  
  #initialise some vectors for storing random coordinates
  my_random_start  <- vector()
  my_random_end    <- vector()
  my_random_chr    <- vector()
  regions <- list()
  for(i in 1:1){
    my_random_chr[i] <- sample(x = my_chr,size = 1)
    my_max <- my_chr_max[[my_random_chr[i]]]  - span
    my_random_start[i] <- round(runif(n = 1, min = my_chr_min[[my_chr[i]]], max = my_max))
    my_random_end[i] <- my_random_start[i] + span
    region <- fd[fd$seqnames %in% my_random_chr[i] & fd$start >= my_random_start[i] & fd$end <= my_random_end[i],]
    if(nrow(region) >= min_cpg){
      regions[[i]] <- region
    }else{
      while(nrow(region) < min_cpg){
        my_random_chr[i] <- sample(x = my_chr,size = 1)
        my_max <- my_chr_max[[my_random_chr[i]]]  - span
        my_random_start[i] <- round(runif(n = 1, min = my_chr_min[[my_chr[i]]], max = my_max))
        my_random_end[i] <- my_random_start[i] + span
        region <- fd[fd$seqnames %in% my_random_chr[i] & fd$start >= my_random_start[i] & fd$end <= my_random_end[i],]
      }
      regions[[i]] <- region
    }
  }
  return(regions)
}
##3.2. Simulation function
simulations_FPR <- function(methy, n = 100, methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "beta"), chr = NULL , start = NULL, end = NULL, epi_params = epimutacions::epi_parameters(),bump_cutoff = 0.1, min_cpg = 3,  verbose = TRUE){
  
  #select "n" control samples randomly
  control_samples  <- which(GSE51032$status == "control")
  control_samples_ncol <- sample(control_samples, size = n, replace = FALSE)
  control_panel <- GSE51032[,control_samples_ncol]
  #select case sample randomly
  
  #total_case_samples <- ncol(case_samples)
  case_samples <- which(GSE51032$status == "case")
  case_samples_ncol <- sample(case_samples, size = 4, replace = FALSE)
  case_samples <- GSE51032[,case_samples_ncol]
  
  regions <- random_regions(GSE51032[,c(colnames(control_samples), colnames(case_samples))], span = 20000, min_cpg = 3)
  
  
  rst <- do.call(rbind, lapply(seq_len(length(methods)), function(j) {
    rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii) {
      epimutacions::epimutations(case_samples, control_panel, 
                                 method = methods[j], 
                                 chr = unique(regions[[1]][,"seqnames"]),
                                 start = min(regions[[1]][,"start"]) - 10000,
                                 end = max(regions[[1]][,"end"]) + 10000,   
                                 epi_params = epi_params,
                                 min_cpg = min_cpg,
                                 bump_cutoff = bump_cutoff,
                                 verbose = verbose)
      
    }))
  }))
  if(is.null(rst)){
    return(NA)
  }
  if(nrow(rst) == 0){
    return(NA)
  }else{
    rst <- as.data.frame(rst)
    gr_results <- GenomicRanges::makeGRangesFromDataFrame(rst)
    library(GenomicRanges)
    rst <- rst[gr_results %outside% gr_found,]
    return(rst)
  }
}

#4. Simulation
out <- simulations_FPR(GSE51032, n = n) 

#5. Save the results
output.filename <- paste0("outputs/outputfile-","FPR-", variable,"-", n, "-", i, ".rda")
save(out, file=output.filename)
