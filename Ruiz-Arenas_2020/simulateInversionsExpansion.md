# Summary

Evalute recombClust accuracy in simulated invesions with expansion. This script analyses the inversion generated with nextflow:

```{bash}
~/nextflow run simulateInversionsExpansion.nf --initialxml in.xml --paramFile par.xml \
 --lens 500000 --freqs 0.5 --sim 1000 -resume
```

invertFREGENE is buggy to generate inversions with expansions. Generate more simulations
to ensure we have at least 100 simulations in total.

# Compute accuracy based on 100 first simulations 

```{r}
library(dplyr)

recombFiles <- dir("results/recombClustClassifications", full.names = TRUE)[1:100] 

mapDF <- data.frame(recomb = recombFiles, stringsAsFactors = FALSE) %>%
  mutate(idFile = gsub("recombClustClassifications", "", recomb),
         idFile = gsub("recombClass.Rdata", "log_ids.txt", idFile))

getAccuracy <- function(recombFile, idsFile){
  
  load(recombFile)
  inv_ids <- read.table(idsFile)
  ids <- rep(0, 4000)
  ids[inv_ids$V1 + 1] <- 1
  
  max(mean(ids == (class == 2)), mean(ids == (class == 1)))
}

sapply(seq_len(nrow(mapDF)), function(x) getAccuracy(mapDF[x, 1], mapDF[x, 2]))
```
