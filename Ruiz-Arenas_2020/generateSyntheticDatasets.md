# Summary

This document contains the code to generate the syntethic datasets used to evaluate recombClust performance. 

# Load libraries

```{r}
library(recombClust)
library(parallel)
```

# Initialization and intrablock linkage datasets

We generated datasets to test initialization of mixture parameter and with different linkage disequilibrium between the SNPs in the blocks. Each simulation was generated 200 times. This datasets consist of a population of two SNP-blocks, each comprised of two SNPs. Inside a population, two subpopulation of chromosomes exist.

```{r}
## Define general parameters
numSim <- 200
set.seed(0)
```


## Create SNP blocks 

Blocks without linkage between SNPs (D' < 1, R2 < 1). Each blocks contains 1000 individuals. The SNPs are represented with 0 or 1, where 0 is the major allele (AF between 0.55 and 0.95).

```{r}
NoLDtablesL <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.55, max = 0.95)
  probB <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), replace = TRUE),
                  sample(c("0", "1"), 1000, prob = c(probB, 1-probB), replace = TRUE)), 
                nrow = 1000)
  mat  
})

NoLDtablesR <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.55, max = 0.95)
  probD <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), replace = TRUE),
                  sample(c("0", "1"), 1000, prob = c(probD, 1-probD), replace = TRUE)), 
                nrow = 1000)
  mat  
})
```

Blocks with high LD (D' = 1, R2 < 1)

```{r}
LDtablesL <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.55, max = 0.95)
  probB <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), replace = TRUE)),
                  sort(sample(c("0", "1"), 1000, prob = c(probB, 1-probB), replace = TRUE))), 
                nrow = 1000)
  mat  
})

LDtablesR <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.55, max = 0.95)
  probD <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), replace = TRUE)),
                  sort(sample(c("0", "1"), 1000, prob = c(probD, 1-probD), replace = TRUE))), 
                nrow = 1000)
  mat  
})
```

Blocks with complete LD (D' = 1, R2 = 1)

```{r}
R2tablesL <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(0, nrow = 1000, ncol = 2)
  mat[, 1] <- c(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), replace = TRUE))
  A <- mat[, 1] == "0"
  mat[A, 2] <- "0"
  mat[!A, 2] <- "1"
  mat
})

R2tablesR <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(0, nrow = 1000, ncol = 2)
  mat[, 1] <- c(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), replace = TRUE))
  C <- mat[, 1] == "0"
  mat[C, 2] <- "0"
  mat[!C, 2] <- "1"
  mat
})
```

## Create populations by binding two blocks 

We generated populations with two SNP-blocks, by binding two syntethic blocks. SNP-blocks could be merged using two possiblities: No linkage (D' < 1, R2 < 1) or complete linkage (D' = 1, R2 = 1). 

### Define binding function

```{r}
bindMat <- function(mat) apply(mat, 1, paste, collapse = "")
```

### Populations without linkage between SNP-blocks (D' < 1, R2 < 1)

Subpopulations are named with NoLD (populations with full linkage between the blocks) and the linkage between the SNPs (N: no linkage, L: high linkage, R: complete linkage), in the first and second blocks.


```{r}
NoLD_N_N_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(NoLDtablesL[[x]])),  sample(bindMat(NoLDtablesR[[x]])))
})

NoLD_N_L_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(NoLDtablesL[[x]])),  sample(bindMat(LDtablesR[[x]])))
})

NoLD_N_R_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(NoLDtablesL[[x]])), sample(bindMat(R2tablesR[[x]])))
})

NoLD_L_N_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(LDtablesL[[x]])),  sample(bindMat(NoLDtablesR[[x]])))
})

NoLD_L_L_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(LDtablesL[[x]])), sample(bindMat(LDtablesR[[x]])))
})

NoLD_L_R_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(LDtablesL[[x]])),  sample(bindMat(R2tablesR[[x]])))
})

NoLD_R_N_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(R2tablesL[[x]])), sample(bindMat(NoLDtablesR[[x]])))
})

NoLD_R_L_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(R2tablesL[[x]])), sample(bindMat(LDtablesR[[x]])))
})

NoLD_R_R_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(R2tablesL[[x]])),  sample(bindMat(R2tablesR[[x]])))
})

NoLDTables <- list(NoLD_N_N = NoLD_N_N_Tables,
                   NoLD_N_L = NoLD_N_L_Tables,
                   NoLD_N_R = NoLD_N_R_Tables,
                   NoLD_L_N = NoLD_L_N_Tables,
                   NoLD_L_L = NoLD_L_L_Tables, 
                   NoLD_L_R = NoLD_L_R_Tables,
                   NoLD_R_N = NoLD_R_N_Tables,
                   NoLD_R_L = NoLD_R_L_Tables,
                   NoLD_R_R = NoLD_R_R_Tables)
```

### Populations with full linkage between SNP-blocks (D' = 1, R2 = 1)

Subpopulations are named with R2 (populations with full linkage between the blocks) and the linkage between the SNPs (N: no linkage, L: high linkage, R: complete linkage). 

```{r}
R2_N_Tables <- lapply(1:numSim, function(x){
  vec <- bindMat(NoLDtablesL[[x]])
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "11", 
                     ifelse(mat[, 2] == "10", "10", 
                            ifelse(mat[, 2] == "01", "01", "00")))
  mat
})

R2_L_Tables <- lapply(1:numSim, function(x){
  vec <- bindMat(LDtablesL[[x]])
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "11", 
                     ifelse(mat[, 2] == "10", "10", 
                            ifelse(mat[, 2] == "01", "01", "00")))
  mat
})

R2_R_Tables <- lapply(1:numSim, function(x){
  vec <- bindMat(R2tablesL[[x]])
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "11", 
                     ifelse(mat[, 2] == "10", "10", 
                            ifelse(mat[, 2] == "01", "01", "00")))
  mat
})

R2tables <- list(R2_N = R2_N_Tables, R2_L = R2_L_Tables, R2_R = R2_R_Tables )
```

## Create mixture of subpopulations

Merge previous populations by concatening the individuals. We generated all possible combination of previous subpopulations.

### Define function

```{r}
mixTabs <- function(ind, list1, list2) {
  rbind(list1[[ind]], list2[[ind]])
}
```

### Create mixture populations

Elements name contains the names of the two subpopulations.

```{r}
mixedTables <- list(
  NoLD_N_N_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_N_N_Tables, list2 = R2_N_Tables),
  NoLD_N_N_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_N_N_Tables, list2 = R2_L_Tables),
  NoLD_N_N_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_N_N_Tables, list2 = R2_R_Tables),
  NoLD_N_L_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_N_L_Tables, list2 = R2_N_Tables),
  NoLD_N_L_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_N_L_Tables, list2 = R2_L_Tables),
  NoLD_N_L_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_N_L_Tables, list2 = R2_R_Tables),
  NoLD_N_R_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_N_R_Tables, list2 = R2_N_Tables),
  NoLD_N_R_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_N_R_Tables, list2 = R2_L_Tables),
  NoLD_N_R_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_N_R_Tables, list2 = R2_R_Tables),
  NoLD_L_N_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_L_N_Tables, list2 = R2_N_Tables),
  NoLD_L_N_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_L_N_Tables, list2 = R2_L_Tables),
  NoLD_L_N_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_L_N_Tables, list2 = R2_R_Tables),
  NoLD_L_L_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_L_L_Tables, list2 = R2_N_Tables),
  NoLD_L_L_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_L_L_Tables, list2 = R2_L_Tables),
  NoLD_L_L_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_L_L_Tables, list2 = R2_R_Tables),
  NoLD_L_R_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_L_R_Tables, list2 = R2_N_Tables),
  NoLD_L_R_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_L_R_Tables, list2 = R2_L_Tables),
  NoLD_L_R_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_L_R_Tables, list2 = R2_R_Tables),
  NoLD_R_N_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_R_N_Tables, list2 = R2_N_Tables),
  NoLD_R_N_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_R_N_Tables, list2 = R2_L_Tables),
  NoLD_R_N_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_R_N_Tables, list2 = R2_R_Tables),
  NoLD_R_L_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_R_L_Tables, list2 = R2_N_Tables),
  NoLD_R_L_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_R_L_Tables, list2 = R2_L_Tables),
  NoLD_R_L_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_R_L_Tables, list2 = R2_R_Tables),
  NoLD_R_R_R2_N = lapply(1:numSim, mixTabs, list1 = NoLD_R_R_Tables, list2 = R2_N_Tables),
  NoLD_R_R_R2_L = lapply(1:numSim, mixTabs, list1 = NoLD_R_R_Tables, list2 = R2_L_Tables),
  NoLD_R_R_R2_R = lapply(1:numSim, mixTabs, list1 = NoLD_R_R_Tables, list2 = R2_R_Tables))
```


## Apply recombClust

We run recombClust initializating mixture parameter (prob0) at a low or high value:

```{r}
allresHigh <- mclapply(c(NoLDTables, R2tables, mixedTables), function(x) 
  lapply(x, function(y) LDmixtureModel(y, maxSteps = 1e3, prob0 = 0.95)), 
  mc.cores = 20)

allresLow  <- mclapply(c(NoLDTables, R2tables, mixedTables), function(x) 
  lapply(x, function(y) LDmixtureModel(y, maxSteps = 1e3, prob0 = 0.05)), 
  mc.cores = 20)

save(allresHigh, allresLow, file = "testProb0.Rdata")
```

# Nucleotide divergence

We generated datasets to test whether nucleotide divergence between subpopulations affect recombClust performance. Each simulation was generated 200 times.

```{r}
numSim <- 200
set.seed(0)
```


## Create blocks

Blocks will be created with SNPs in high LD (D' = 1, R2 < 1, L of previous datasets). Left blocks are named with a and b and right blocks with c and d. We generated blocks with different major allele. If letter is upper case (e.g. A), the major allele is 0; if the letter is lower case (e.g. a), the major allele is 1.

```{r}
LDtablesAB <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.55, max = 0.95)
  probB <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), replace = TRUE)),
                  sort(sample(c("0", "1"), 1000, prob = c(probB, 1-probB), replace = TRUE))), 
                nrow = 1000)
  mat  
})

LDtablesAb <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.55, max = 0.95)
  probB <- runif(1, min = 0.05, max = 0.45)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), 
                              replace = TRUE)),
                  sort(sample(c("0", "1"), 1000, prob = c(probB, 1-probB), 
                              replace = TRUE), decreasing = TRUE)), 
                nrow = 1000)
  mat  
})

LDtablesaB <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.05, max = 0.45)
  probB <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), 
                              replace = TRUE), decreasing = TRUE),
                  sort(sample(c("0", "1"), 1000, prob = c(probB, 1-probB), replace = TRUE))), 
                nrow = 1000)
  mat  
})

LDtablesab <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.05, max = 0.45)
  probB <- runif(1, min = 0.05, max = 0.45)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), 
                              replace = TRUE), decreasing = TRUE),
                  sort(sample(c("0", "1"), 1000, prob = c(probB, 1-probB), 
                              replace = TRUE), decreasing = TRUE)), 
                nrow = 1000)
  mat  
})

LDtablesCD <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.55, max = 0.95)
  probD <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), replace = TRUE)),
                  sort(sample(c("0", "1"), 1000, prob = c(probD, 1-probD), replace = TRUE))), 
                nrow = 1000)
  mat  
})

LDtablesCd <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.55, max = 0.95)
  probD <- runif(1, min = 0.05, max = 0.45)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), replace = TRUE)),
                  sort(sample(c("0", "1"), 1000, prob = c(probD, 1-probD), 
                              replace = TRUE), decreasing = TRUE)), 
                nrow = 1000)
  mat  
})

LDtablescD <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.05, max = 0.45)
  probD <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), 
                              replace = TRUE), decreasing = TRUE),
                  sort(sample(c("0", "1"), 1000, prob = c(probD, 1-probD), replace = TRUE))), 
                nrow = 1000)
  mat  
})

LDtablescd <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.05, max = 0.45)
  probD <- runif(1, min = 0.05, max = 0.45)
  mat <- matrix(c(sort(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), 
                              replace = TRUE), decreasing = TRUE),
                  sort(sample(c("0", "1"), 1000, prob = c(probD, 1-probD), 
                              replace = TRUE), decreasing = TRUE)), 
                nrow = 1000)
  mat  
})
leftTables <- list(AB = LDtablesAB, Ab = LDtablesAb, aB = LDtablesaB, ab = LDtablesab)
rightTables <- list(CD = LDtablesCD, Cd = LDtablesCd, cD = LDtablescD, cd = LDtablescd)
```

## Create populations by binding two blocks 

We generated populations with two SNP-blocks, by binding two syntethic blocks. SNP-blocks could be merged using two possiblities: No linkage (D' < 1, R2 < 1) or complete linkage (D' = 1, R2 = 1). 

### Populations without linkage between SNP-blocks (D' < 1, R2 < 1)

```{r}
bindMat <- function(mat) apply(mat, 1, paste, collapse = "")
makeNoLDTables <- function(index, list1, list2){
  cbind(sample(bindMat(list1[[index]])),  sample(bindMat(list2[[index]])))
}

NoLDTablesAlleles <- lapply(leftTables, function(leftList) {
  lapply(rightTables, function(rightList){
    lapply(1:numSim, function(ind) makeNoLDTables(ind, leftList, rightList))
  })
})
NoLDTablesAlleles <- unlist(NoLDTablesAlleles, recursive = FALSE)
```

### Populations with full linkage between SNP-blocks (D' = 1, R2 = 1)

```{r}
R2TablesAlleles <- lapply(1:numSim, function(x){
  vec <- bindMat(LDtablesAB[[x]])
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "11", 
                     ifelse(mat[, 2] == "10", "10", 
                            ifelse(mat[, 2] == "01", "01", "00")))
  mat
})
```

## Create mixture populations

Combine subpopulations without linkage with subpopulations with linkage.

```{r}
mixTabs <- function(ind, list1, list2) {
  rbind(list1[[ind]], list2[[ind]])
}

mixedTablesAlleles <- lapply(NoLDTablesAlleles, function(list1) {
  lapply(1:numSim, function(ind) mixTabs(ind, list1, R2TablesAlleles))
})
```

## Apply recombClust

```{r}
allresAllele  <- mclapply(c(NoLDTablesAlleles, list(R2 = R2TablesAlleles), mixedTablesAlleles), function(x) 
  lapply(x, function(y) LDmixtureModel(y, maxSteps = 1e3, prob0 = 0.05)), 
  mc.cores = 20)
save(allresAllele, file = "testAllele.Rdata")
```

# Mixture proportion

We evaluated recombClust performance under different subpopulations proportion. We reuse the blocks generated in the previous step.

## Define functions

```{r}
bindMat <- function(mat) apply(mat, 1, paste, collapse = "")
makeNoLDTables <- function(index, list1, list2){
  cbind(sample(bindMat(list1[[index]])),  sample(bindMat(list2[[index]])))
}

mixTabsProps <- function(ind, list1, list2, props) {
  tab1 <- list1[[ind]]
  tab2 <- list2[[ind]]
  if (props > 0.5){
    props <- 1 - props
    rows1 <- props*nrow(tab2)/(1 - props)
    tab1 <- tab1[sample(nrow(tab1), rows1), ]
  } else {
    rows2 <- props*nrow(tab2)/(1 - props)
    tab2 <- tab2[sample(nrow(tab2), rows2), ]
  }
  rbind(tab1, tab2)
}
```

## Create populations by binding two blocks 

We generated populations with two SNP-blocks, by binding two syntethic blocks. SNP-blocks could be merged using two possiblities: No linkage (D' < 1, R2 < 1) or complete linkage (D' = 1, R2 = 1). 

### Populations without linkage between SNP-blocks (D' < 1, R2 < 1)

```{r}
NoLDTables1 <- lapply(1:numSim, function(ind) 
  makeNoLDTables(ind, LDtablesAB, LDtablesCD))

NoLDTables2 <- lapply(1:numSim, function(ind) 
  makeNoLDTables(ind, LDtablesab, LDtablescd))
```

### Populations with full linkage between SNP-blocks (D' = 1, R2 = 1)

```{r}
R2Tables1 <- lapply(1:numSim, function(x){
  vec <- bindMat(LDtablesAB[[x]])
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "11", 
                     ifelse(mat[, 2] == "10", "10", 
                            ifelse(mat[, 2] == "01", "01", "00")))
  mat
})

R2Tables2 <- lapply(1:numSim, function(x){
  vec <- bindMat(LDtablesab[[x]])
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "00", 
                     ifelse(mat[, 2] == "10", "01", 
                            ifelse(mat[, 2] == "01", "10", "11")))
  mat
})
```

## Create mixture populations

Combine subpopulations without linkage with subpopulations with linkage. We generated all combinations of subpopulations without LD (NoLD) and with full linkage (LD). Each combination was performed 9 times, with different subpopulation proportions (from 1:9 to 9:1).

```{r}
mixTabs <- function(ind, list1, list2) {
  rbind(list1[[ind]], list2[[ind]])
}

props <- seq(0.1, 0.9, 0.1)
mixedTablesCombsProps <- list(
  NoLD_NoLD = lapply(props, function(prop) lapply(1:numSim, function(ind) 
    mixTabsProps(ind, NoLDTables1, NoLDTables2, prop))),
  NoLD_LD =  lapply(props, function(prop) lapply(1:numSim, function(ind) 
    mixTabsProps(ind, NoLDTables1, R2Tables2, prop))),
  LD_NoLD = lapply(props, function(prop) lapply(1:numSim, function(ind) 
    mixTabsProps(ind, R2Tables1, NoLDTables2, prop))),
  LD_LD = lapply(props, function(prop) lapply(1:numSim, function(ind) 
    mixTabsProps(ind, R2Tables1, R2Tables2, prop)))
)
```

## Apply recombClust

```{r}
allresCombsProps  <- lapply(mixedTablesCombsProps, function(z) lapply(z, function(x)
  mclapply(x, function(y) LDmixtureModel(y, maxSteps = 1e3, prob0 = 0.05), 
           mc.cores = 20)))

save(allresCombs, allresCombsProps, file = "testCombs.Rdata")
```


# Subpopulations with different combinations of blocks 

In this last segment, we generated two subpopulations consisting of 10 SNP-blocks of 2 SNPs. One subpopulation was in high linkage in odd blocks and no linkage in even blocks. For the other subpopulation, the linkage between the blocks was inversed: low linkage in odd blocks and high linkage in even blocks.

## Create block datasets 

Each subpopulation will contain 500 chromosomes.

```{r}
numSim <- 10 ## Generate population with ten pairs of blocks

# Each population has 500 individuals
set.seed(3)
```

We generated all blocks without linkage between the SNPs:

```{r}
NoLDtablesL <- lapply(1:numSim, function(x){
  probA <- runif(1, min = 0.55, max = 0.95)
  probB <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sample(c("0", "1"), 1000, prob = c(probA, 1-probA), replace = TRUE),
                  sample(c("0", "1"), 1000, prob = c(probB, 1-probB), replace = TRUE)), 
                nrow = 1000)
  mat  
})

NoLDtablesR <- lapply(1:numSim, function(x){
  probC <- runif(1, min = 0.55, max = 0.95)
  probD <- runif(1, min = 0.55, max = 0.95)
  mat <- matrix(c(sample(c("0", "1"), 1000, prob = c(probC, 1-probC), replace = TRUE),
                  sample(c("0", "1"), 1000, prob = c(probD, 1-probD), replace = TRUE)), 
                nrow = 1000)
  mat  
})


```

## Create populations by binding two blocks 

We generated two types of blocks combinations, one without linkage (NoLD_Tables), and another with full linkage (R2_Tables).

```{r}
### Create NoLD pairs datasets ####
bindMat <- function(mat) apply(mat, 1, paste, collapse = "")

## Blocks not in LD
NoLD_Tables <- lapply(1:numSim, function(x){
  cbind(sample(bindMat(NoLDtablesL[[x]])),  sample(bindMat(NoLDtablesR[[x]])))
})


## Blocks in LD
R2_Tables <- lapply(1:numSim, function(x){
  vec <- bindMat(NoLDtablesL[[x]])
  mat <- cbind(vec, vec)
  mat[, 2] <- ifelse(mat[, 2] == "11", "11", 
                     ifelse(mat[, 2] == "10", "10", 
                            ifelse(mat[, 2] == "01", "01", "00")))
  mat
})
```

## Create two subpopulations and mix

```{r}
pop1 <- c(NoLD_Tables[seq(1, 10, 2)], R2_Tables[seq(2, 10, 2)])
pop1 <- pop1[unlist(lapply(1:5, function(x) c(x, x + 5)))]

pop2 <- c(R2_Tables[seq(1, 10, 2)], NoLD_Tables[seq(2, 10, 2)])
pop2 <- pop2[unlist(lapply(1:5, function(x) c(x, x + 5)))]

mixpop <- mapply(rbind, pop1, pop2, SIMPLIFY = FALSE)
```

## Apply recombClust

```{r}
modelRes <- lapply(mixpop, LDmixtureModel, maxSteps = 1e3, prob0 = 0.5)
```

## Apply Genotypes clustering

We compared recombClust with invClust, a method to cluster individuals based on genotypes.

### Load libraries

```{r}
library(snpStats)
library(invClust)
```

### Preprocess

Convert data to SNPmatrix

```{r}
mixpopChr <- Reduce(cbind, mixpop)
#### Split SNP pairs into SNPs
mixpopChr <- t(apply(mixpopChr, 1, function(x) unlist(strsplit(x, ""))))
## Convert SNPs to numeric
class(mixpopChr) <- "numeric"
rownames(mixpopChr) <- 1:2000
colnames(mixpopChr) <- c(letters, LETTERS)[1:40]

### Test genotypes
mixpopGeno <- t(sapply(seq(1, 2000, 2), function(i) mixpopChr[i, ] + mixpopChr[i + 1, ]))
hapGeno <- new("SnpMatrix", mixpopGeno + 1)
```

### Apply invClust

```{r}
invClustGeno <- invClust(roi = invdf, wh = 1, geno = hapGeno,
                        annot = map, dim = 2, tol = 0.5)


save(pop1, pop2, invClustGeno, modelRes, file = "simulatedPopulations.Rdata")

```

