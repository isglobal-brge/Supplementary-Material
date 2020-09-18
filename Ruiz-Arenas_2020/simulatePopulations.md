# Summary

Evalute recombClust accuracy in synthetic populations with different number of chromosomes or SNP-blocks. Simulate two populations with different combinations of blocks. The algorithm is described in the section *Subpopulations with different combinations of blocks* from generateSyntheticDatasets.Rmd.

## Load libraries

```{r}
library(recombClust)
library(parallel)
library(cluster)
```

# Define functions

```{r}
makeRecombTables <- function(nchr, i){
  probA <- runif(1, min = 0.55, max = 0.95)
  probB <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sample(c("0", "1"), nchr, prob = c(probA, 1-probA), replace = TRUE),
                  sample(c("0", "1"), nchr, prob = c(probB, 1-probB), replace = TRUE)), 
                nrow = nchr)
  mat  
}

makeLinkTables <- function(mat){
  vec <- bindMat(mat)
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "11", 
                     ifelse(mat[, 2] == "10", "10", 
                            ifelse(mat[, 2] == "01", "01", "00")))
  mat
}
bindMat <- function(mat) apply(mat, 1, paste, collapse = "")

getPopulation <- function(nchr, nBlocks, odd = TRUE){
  
  ### Create Recomb pairs datasets ####
  recTablesL <- lapply(seq_len(nBlocks), makeRecombTables, nchr = nchr)
  recTablesR <- lapply(seq_len(nBlocks), makeRecombTables, nchr = nchr)
  
  recBlocks <- lapply(seq_len(nBlocks), function(x){
    cbind(sample(bindMat(recTablesL[[x]])),  sample(bindMat(recTablesR[[x]])))
  })
  
  ## Create linkage Blocks
  linkBlocks <- lapply(recTablesL, makeLinkTables)
  
  if (odd){
    pop <- c(recBlocks[seq(1, nBlocks, 2)], linkBlocks[seq(2, nBlocks, 2)])
  } else {
    pop <- c(linkBlocks[seq(1, nBlocks, 2)], recBlocks[seq(2, nBlocks, 2)])
  }
  pop <- pop[unlist(lapply(seq_len(nBlocks/2), function(x) c(x, x + nBlocks/2)))]
  pop
}

makeMixture <- function(nchr, nBlocks){
  pop1 <- getPopulation(nchr = nchr, nBlocks = nBlocks)
  pop2 <- getPopulation(nchr = nchr, nBlocks = nBlocks, odd = FALSE)
  mixpop <- mapply(rbind, pop1, pop2, SIMPLIFY = FALSE)
  mixpop
}
```

# Number of individuals

Generate scenarios with 100 SNP-blocks per population and varying number of individuals. Each scenario is simulated 1000 times.

```{r}
nBlocks <- 100 ## SNP block pairs per populations
numSim <- 1000 ## Number simulations
nInds <- c(10, 15, 20, 25, 30)
```

## Generate mixture populations

```{r}
set.seed(0)
mixPopulations <- lapply(nInds, function(nchr){
  lapply(seq_len(numSim), function(x) makeMixture(nchr = nchr, nBlocks = nBlocks))
})
```

## Generate single populations

```{r}
set.seed(0)
singlePopulations <- lapply(nInds*2, function(nchr){
  lapply(seq_len(numSim), function(x) getPopulation(nchr = nchr, nBlocks = nBlocks))
})
save(mixPopulations, singlePopulations, file = "simulations/SimBlockPops.Rdata")
```

## Apply recombClust

```{r}
runRecomb <- function(blockList){
  lapply(blockList, recombClust:::LDmixtureModel, maxSteps = 1e3, prob0 = 0.5)
}

mixRecombRes <- lapply(mixPopulations, function(l){
  mclapply(l, runRecomb, mc.cores = 30)
})
singleRecombRes <- lapply(singlePopulations, function(l){
  mclapply(l, runRecomb, mc.cores = 30)
})
save(mixRecombRes, singleRecombRes, file = "simulations/SimBlockRecombRes.Rdata")
```


## Get PC and run clustering

```{r}
getClust <- function(recObj){
  indsmat <- do.call(cbind, lapply(recObj, `[[`, "r1"))
  pc <- prcomp(indsmat, rank. = 2)$x
  class <- kmeans(pc, centers = 2, nstart = 100)$cluster
  list(pc = pc, class = class)
}
mixPCs <- lapply(mixRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})
singlePCs <- lapply(singleRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})

```

## Compute silhouette values

```{r}
computeSilhouette <- function(l){
  summary(silhouette(l$class, dist(l$pc)))$si.summary[4]
}

mixSilhouete <- lapply(mixPCs, function(x) sapply(x, computeSilhouette))
singleSilhouete <- lapply(singlePCs, function(x) sapply(x, computeSilhouette))
```

# Number of SNP-block

Generate scenarios with 100 individuals per population and varying number of SNP-blocks. Each scenario is simulated 1000 times.

```{r}
nBlocks <- c(10, 15, 20, 40, 70, 100) ## SNP block pairs per populations
numSim <- 1000 ## Number simulations
nInds <- 100
```

## Generate mixture populations

```{r}
set.seed(0)
mixPopulationsBlock <- lapply(nBlocks, function(nBlocks){
  lapply(seq_len(numSim), function(x) makeMixture(nchr = nInds, nBlocks = nBlocks))
})
```

## Generate single populations

```{r}
set.seed(0)
singlePopulationsBlock <- lapply(nBlocks, function(nBlocks){
  lapply(seq_len(numSim), function(x) getPopulation(nchr = nInds*2, nBlocks = nBlocks))
})
save(mixPopulationsBlock, singlePopulationsBlock, file = "simulations/SimBlockPopsBlocks.Rdata")
```

## Apply recombClust

```{r}
mixBlockRecombRes <- lapply(mixPopulationsBlock, function(l){
  mclapply(l, runRecomb, mc.cores = 30)
})
singleBlockRecombRes <- lapply(singlePopulationsBlock, function(l){
  mclapply(l, runRecomb, mc.cores = 30)
})
save(mixBlockRecombRes, singleBlockRecombRes, file = "simulations/SimBlockRecombResBlocks.Rdata")
```

## Get PC and run clustering

```{r}
mixBlockPCs <- lapply(mixBlockRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})
singleBlockPCs <- lapply(singleBlockRecombRes, function(l){
  mclapply(l, getClust, mc.cores = 30)
})
```
