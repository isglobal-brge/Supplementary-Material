library(epimutacions)

#1. Load the data
n <- Sys.getenv("N")
i <- Sys.getenv("I")

print(paste0("n=", n, ", i=", i)) 

input.filename <- paste0("inputs/GSE111629",".rda")
variable <- load(file=input.filename)
#object   <- eval(parse(text=variable))
#print(object)
load("inputs/epimutations_found_GSE111629.rda")
gr_found <- GenomicRanges::makeGRangesFromDataFrame(found_GSE111629)
load("inputs/names_GSE111629.rda")

#2. Filtering inputs
found <- which(colnames(GSE111629) %in% unique(names_GSE111629))
GSE111629 <- GSE111629[,-found]

control_samples <- GSE111629[,GSE111629$disease_state == "control"]
case_samples <- GSE111629[,GSE111629$disease_state == "case"]

case_samples <- case_samples[,-which(colnames(case_samples) == "GSM3035497")]

#3. Functions
##3.1. Simulation function
simulations_FPR <- function(case_samples, control_samples, n = 100, methods = c("manova", "mlm", "mahdistmcd", "isoforest", "barbosa", "qn"), chr = NULL , start = NULL, end = NULL, epi_params = epimutacions::epi_parameters(),bump_cutoff = 0.1, min_cpg = 3,  verbose = TRUE){
  
  if(is.null(case_samples)){
    stop("Please provide a valid 'case_samples'")
  }
  if(is.null(control_samples)){
    stop("Please provide a valid 'control_samples'")
  }
  print(1)
  #select "n" control samples randomly
  total_control_samples <- ncol(control_samples)
  samples_ncol <- sample.int(total_control_samples, size = n, replace = FALSE)
  control_samples <- control_samples[,samples_ncol] 
  #select case sample randomly
  
  #total_case_samples <- ncol(case_samples)
  print(2)
  samples_ncol <- sample.int(ncol(case_samples), size = 4, replace = FALSE)
  case_samples <- case_samples[,samples_ncol] 
  print(3)
  print(colnames(case_samples))
  print(colnames(control_samples))
  rst <- do.call(rbind, lapply(seq_len(length(methods)), function(j) {
    print(4)
    rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii) {
      print(5)
      epimutacions::epimutations(case_samples[,ii], 
                   control_panel = control_samples, 
                   method = methods[j], 
                   chr = chr, 
                   start = start, 
                   end = end,  
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
##3.2 Random region generation
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

#4. Simulation
regions <- random_regions(GSE111629, span = 20000, min_cpg = 3)
out <- simulations_FPR(case_samples, control_samples,
                       n = n,
                       chr = unique(regions[[1]][,"seqnames"]),
                       start = min(regions[[1]][,"start"]) - 10000,
                       end = max(regions[[1]][,"end"]) + 10000) 

#5. Save the results
output.filename <- paste0("outputs/outputfile-", variable,"-out-", n, "-", i, ".rda")
save(out, file=output.filename)
