# Simulation functions

# 1. TPR
TPR <- function(file_name, case_sample, methods, n,  chr, start, end){
  #1.1 Load the the data
  input.filename <- paste0("inputs/",file_name,".rda")
  variable <- load(file=input.filename)
  object   <- eval(parse(text=variable))
  #1.2 Filtering inputs
  control_samples  <- which(object$status == "control")
  case_samples <- object[,case_sample]
  #1.3 Epimutations detection
  results <- do.call(rbind, lapply(seq_len(length(methods)), function(j) {
    rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii){
      samples_ncol <- sample(control_samples, size = n, replace = FALSE)
      control_panel <- GSE51032[,samples_ncol]
      time <- system.time(res <- epimutacions::epimutations(case_samples, control_panel,
                                                            method = methods[j], 
                                                            chr = chr[ii],
                                                            start = start[ii],
                                                            end = end[ii]))
      res$time <- time[[2]]
      res
      
    }))
  }))
  
  #1.4. Save the results
  output.filename <- paste0("outputs/outputfile-","TPR-", variable,"-", n, "-", i, ".rda")
  save(results, file=output.filename)
  print(paste0("Written to file '", output.filename, "'"))
}

# 2. Random region generation
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

# 3. FPR 
FPR <- function(file_name, epi_found, methods, n){
  #1. Load the data
  input.filename <- paste0("inputs/",file_name,".rda")
  variable <- load(file=input.filename)
  object <- eval(parse(text = variable))
  load(paste0("inputs/",epi_found,".rda"))
  
  #2. Filtering inputs
  found <- which(colnames(object) %in% unique(names))
  object <- object[,-found]
  
  ##2.1. select "n" control samples randomly
  control_samples  <- which(object$status == "control")
  control_samples_ncol <- sample(control_samples, size = n, replace = FALSE)
  control_panel <- object[,control_samples_ncol]
  
  ##2.2. select 1 case sample randomly
  #total_case_samples <- ncol(case_samples)
  case_samples <- which(object$status == "case")
  case_samples_ncol <- sample(case_samples, size = 4, replace = FALSE)
  case_samples <- object[,case_samples_ncol]
  
  regions <- random_regions(object[,c(colnames(control_samples), colnames(case_samples))], span = 20000, min_cpg = 3)
  
  
  results <- do.call(rbind,lapply(seq_len(length(methods)), function(j) {
    rst_case <- do.call(rbind, lapply(seq_len(ncol(case_samples)), function(ii) {
      time  <- system.time(res <- epimutacions::epimutations(case_samples, control_panel, 
                                                             method = methods[j], 
                                                             chr = unique(regions[[1]][,"seqnames"]),
                                                             start = min(regions[[1]][,"start"]) - 10000,
                                                             end = max(regions[[1]][,"end"]) + 10000))
      
      res$time <- time[[2]]
      res            
    }))
  }))
  
  #3. Save the results
  output.filename <- paste0("outputs/outputfile-","FPR-", variable,"-", n, "-", i, ".rda")
  save(results, file=output.filename)
  print(paste0("Written to file '", output.filename, "'"))
}


