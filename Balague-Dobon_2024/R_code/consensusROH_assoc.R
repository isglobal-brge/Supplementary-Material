library(GenomicRanges)
library(parallel)
library(arm)

load("/data/dataloh.RData")

for(i in rownames(phenoloh)){
  phenoloh[i,"NLOH"] <- nrow(getloh[getloh$Sample == i,])
}

males <- phenoloh$IID[phenoloh$SEX == 1]
females <- phenoloh$IID[phenoloh$SEX == 2]
phenolohm <- phenoloh[phenoloh$SEX == 1,]
phenolohf <- phenoloh[phenoloh$SEX == 2,]
getlohm <- getloh[getloh$Sample %in% males,]
getlohf <- getloh[getloh$Sample %in% females,]


### CONSENSUS ROH REGIONS - ASSOCIATIONS

phenoname <- "C>GMLA"
phenoname <- "vivo"

gwasloh <- lapply(1:22, function(chr){
  ss <- getloh[,2] == chr
  grchr <- GRanges(paste("chr", getloh[ss,2], sep=""),
                   IRanges(start=getloh[ss,3], end=getloh[ss,4]),
                   sub = getloh[ss,1])
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  blocks <- round(seq(mn,mx, by=500000))
  assocchr <- mclapply(1:(length(blocks)-1), function(bb){
    block <- GRanges(paste("chr", chr, sep=""),
                     IRanges(start=blocks[bb],
                             end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenoloh
    dat$loh <- rep(0,nrow(dat))
    dat[selsubsloh,]$loh <- 1
    dat$pheno <- dat[,phenoname]
    
    assocblock <- summary(bayesglm(
      pheno ~ loh+SEX*AGE+
        euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
      family="binomial", data=dat))$coeff["loh",c(1,4)]
    data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
               OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat[selsubsloh,]$loh==1))
  })
  assocchr <- do.call(rbind, assocchr)
})
gwasloh <- do.call(rbind, gwasloh)
rownames(gwasloh) <- NULL
fr005 <- sum(table(phenoloh[,phenoname]))*0.05
tbm <- gwasloh[gwasloh$freq > fr005, ]
rownames(tb) <- NULL

# males

phenoname <- "C>GMLA"
gwaslohm <- lapply(1:22, function(chr){
  ss <- getlohm[,2] == chr
  grchr <- GRanges(paste("chr", getlohm[ss,2], sep=""),
                   IRanges(start=getlohm[ss,3], end=getlohm[ss,4]),
                   sub = getlohm[ss,1])
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  blocks <- round(seq(mn,mx, by=500000))
  assocchrm <- mclapply(1:(length(blocks)-1), function(bb){
    block <- GRanges(paste("chr", chr, sep=""),
                     IRanges(start=blocks[bb],
                             end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenolohm
    dat$loh <- rep(0,nrow(dat))
    dat[selsubsloh,]$loh <- 1
    dat$pheno <- dat[,phenoname]
    
    assocblock <- summary(bayesglm(
      pheno ~ loh+AGE+
        euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
      family="binomial", data=dat))$coeff["loh",c(1,4)]
    data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
               OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat[selsubsloh,]$loh==1))
  })
  assocchrm <- do.call(rbind, assocchrm)
})
gwaslohm <- do.call(rbind, gwaslohm)
rownames(gwaslohm) <- NULL
fr005m <- sum(table(phenolohm[,phenoname]))*0.05
tbm <- gwaslohm[gwaslohm$freq > fr005m, ]
rownames(tbm) <- NULL

# females

phenoname <- "C>GMLA"
gwaslohf <- lapply(1:22, function(chr){
  ss <- getlohf[,2] == chr
  grchr <- GRanges(paste("chr", getlohf[ss,2], sep=""),
                   IRanges(start=getlohf[ss,3], end=getlohf[ss,4]),
                   sub = getlohf[ss,1])
  mn <- min(start(grchr))
  mx <- max(start(grchr))
  blocks <- round(seq(mn,mx, by=500000))
  assocchrf <- mclapply(1:(length(blocks)-1), function(bb){
    block <- GRanges(paste("chr", chr, sep=""),
                     IRanges(start=blocks[bb],
                             end=blocks[bb+1]))
    selover <- data.frame(findOverlaps(block,grchr))[,2]
    selsubsloh <- unique(as.character(grchr$sub[selover]))
    dat <- phenolohf
    dat$loh <- rep(0,nrow(dat))
    dat[selsubsloh,]$loh <- 1
    dat$pheno <- dat[,phenoname]
    
    assocblock <- summary(bayesglm(
      pheno ~ loh+AGE+
        euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
      family="binomial", data=dat))$coeff["loh",c(1,4)]
    data.frame(chr=chr, start=blocks[bb], end=blocks[bb+1],
               OR=exp(assocblock[1]), P=assocblock[2], freq=sum(dat[selsubsloh,]$loh==1))
  })
  
  assocchrf <- do.call(rbind, assocchrf)
})
gwaslohf <- do.call(rbind, gwaslohf)
rownames(gwaslohf) <- NULL
fr005f <- sum(table(phenolohf[,phenoname]))*0.05
tbf <- gwaslohf[gwaslohf$freq > fr005f, ]
rownames(tbf) <- NULL
