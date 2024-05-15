library(snpStats)
library(haplo.stats)
library(genetics)
library(SNPassoc)
library(parallel)

# EXAMPLE FOR CHROMOSOME 2 HIT

# HAPLOTYPE
load("/data/dataloh.RData")

snpmat <- read.plink("/data/SCOURGE.bed",
                     "/data/SCOURGE.bim",
                     "/data/SCOURGE.fam")

#Annotation data
annot <- read.delim("Z:/SCOURGE/Axiom_SpainBA.na36.r2.a2.annot.csv", sep=",", comment.char ="#")
rsannot <- as.character(annot[,"Extended.RSID"])
names(rsannot) <- annot[,"Probe.Set.ID"]

selsubs <- rownames(phenoloh)
selsubs <- intersect(rownames(snpmat$genotypes), rownames(phenoloh))
phenohit <- phenoloh[selsubs, c(5:16, 55, 74)]
genos <- snpmat$genotypes[selsubs,]
map<- snpmat$map
dim(genos)


### chromosome 2 
hit <- tb[3,c(1:6)]

selsnps <- snpmat$map$chromosome == hit$chr & 
  snpmat$map$position >= hit$start & snpmat$map$position <= hit$end

genohit <- genos[,selsnps]
genomat <- as.matrix(genohit)
maphit <- map[selsnps,]

img <- data.frame(lapply(1:ncol(genomat), function(i) as.numeric(genomat[,i])))
img[img==0] <- NA
img <- img - 1
tbh <- prop.table(table(as.matrix(img)))
tbh

#expected heterozygocity
2*sqrt(tbh[1])*sqrt(tbh[3])
#observed heterozygocity
tbh[2]


getalleles <- lapply(1: ncol(genomat), function(nr){
  nm <- as.numeric(genomat[,nr])
  rs <- rep(NA, nrow(genohit))
  allele <- maphit[nr,]
  rs[which(nm == 1)] <- paste(allele$allele.1, allele$allele.1,sep="")
  rs[which(nm == 2)] <- paste(allele$allele.1, allele$allele.2,sep="")
  rs[which(nm == 3)] <- paste(allele$allele.2, allele$allele.2,sep="")
  factor(rs)})

getalleles <- data.frame(getalleles)
colnames(getalleles) <- colnames(genomat)
datahaplo <- cbind(phenohit, getalleles)
data.hap <- setupSNP(datahaplo, (ncol(phenohit)+1):ncol(datahaplo), sep="")
nmsnps <- colnames(data.hap)[(ncol(phenohit)+1):ncol(datahaplo)]


assoc <- mclapply(5:length(nmsnpsm), function(i){
  nmsnpssel <- nmsnpsm[(i-4):i]
  genoH <- make.geno(data.hapm, nmsnpssel)
  mod.adj <- haplo.glm(`C>GMLA` ~ genoH + AGE*SEX +
                         euPC1 + euPC2 + euPC3 + euPC4 + euPC5 + euPC6 + euPC7 + euPC8 + euPC9 + euPC10,
                       family="binomial",
                       haplo.effect="recessive",
                       locus.label=nmsnpssel,
                       allele.lev=attributes(genoH)$unique.alleles,
                       control = haplo.glm.control(haplo.freq.min=0.05))
  sm <- summary(mod.adj)$coef
  min(sm[grep("geno",rownames(sm)), "p"])
})

p <- unlist(assoc)
plot(maphit[5:nrow(maphit),4]/10^6,-log10(p),
     xlab="Coordinates (Mb)", ylab="-log10(p)",
     pch="", main="Haplotype association")
for(i in 5:length(nmsnps))
  lines(c(maphit[(i-4),4], maphit[i,4])/10^6, c(-log10(p[i-4]), -log10(p[i-4])))




### ANALYZE SIGNIFICANT HAPLOTYPE vs SEVERITY

wh <- which(unlist(assoc)==min(unlist(assoc))) + 4
nmsnpssel <- nmsnps[(wh-4):wh]
genoHm <- make.geno(data.hapm, nmsnpssel)
mod.adj <- haplo.glm(`C>GMLA` ~ genoH + AGE*SEX +
                       euPC1 + euPC2 + euPC3 + euPC4 + euPC5 + euPC6 + euPC7 + euPC8 + euPC9 + euPC10,
                     family="binomial",
                     haplo.effect="recessive",
                     locus.label=nmsnpssel,
                     allele.lev=attributes(genoHm)$unique.alleles,
                     control = haplo.glm.control(haplo.freq.min=0.025))
intervals(mod.adj)

roi <- maphit[gsub("\\.", "-", nmsnpssel), ]
roi <- data.frame(roi, rs=rsannot[rownames(roi)])




### ANALYZE SIGNIFICANT HAPLOTYPE VS ROH

chr <- hit$chr
ss <- getloh[,2] == chr
grchr <- GRanges(paste("chr", getloh[ss,2], sep=""),
                 IRanges(start=getloh[ss,3],
                         end=getloh[ss,4]), sub = getloh[ss,1])
block <- GRanges(paste("chr", chr, sep=""),
                 IRanges(start=min(roi$position),
                         end=max(roi$position)))

selover <- data.frame(findOverlaps(block,grchr))[,2]
selsubsloh <- unique(as.character(grchr$sub[selover]))
selsubsloh <- intersect(selsubsloh, rownames(phenohit))
dat <- phenohit
dat$loh <- rep(0,nrow(dat))
dat[selsubsloh,]$loh <- 1
dat$pheno <- dat[,phenoname]
dat <- dat[dat$SEX == 1,]
assocblock <- summary(bayesglm(pheno~loh+SEX*AGE+
                                 euPC1+euPC2+euPC3+euPC4+euPC5+
                                 euPC6+euPC7+euPC8+euPC9+euPC10,
                               family="binomial",
                               data=dat))$coeff["loh",c(1,4)]
res <- c(exp(assocblock[1]), assocblock)
names(res)[1] <- "OR"
res


loh <- dat$loh
mod.loh <- haplo.glm(loh ~ genoHm,
                     family="binomial",
                     haplo.effect="recessive",
                     locus.label=nmsnpssel,
                     allele.lev=attributes(genoH)$unique.alleles,
                     control = haplo.glm.control(haplo.freq.min=0.025))
intervals(mod.loh)



