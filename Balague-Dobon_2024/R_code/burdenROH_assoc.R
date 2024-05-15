
library(ggplot2)

load("/data/dataloh.RData")

getloh$length <- getloh$End - getloh$Start
df <- merge(getloh, phenoloh[,c(2,5,6)], by.x = 'Sample', by.y = 'IID')

metrics <- data.frame(sample = unique(phenoloh$IID[!is.na(phenoloh$IID)]), NROH = NA, SROH = NA)
rownames(metrics) <- metrics$sample
for (id in rownames(metrics)){
  metrics[id, 'NROH'] <- nrow(df[df$Sample == id,])
  metrics[id, 'SROH'] <- sum(df$length[df$Sample == id])/1000000
}

phenoloh2 <- merge(phenoloh, metrics, by.x = 'IID', by.y = 'sample')
phenoloh2$SEX <- as.factor(phenoloh2$SEX)
phenoloh2$SEX[phenoloh2$SEX == 0] <- NA


summary(bayesglm(vivo ~ NROH+AGE*SEX+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2))
summary(bayesglm(vivo ~ NROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 1,]))
summary(bayesglm(vivo ~ NROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 2,]))



summary(bayesglm(vivo ~ SROH+AGE*SEX+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2))
summary(bayesglm(vivo ~ SROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 1,]))
summary(bayesglm(vivo ~ SROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 2,]))




summary(bayesglm(`C>GMLA` ~ NROH+AGE*SEX+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2))
summary(bayesglm(`C>GMLA` ~ NROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 1,]))
summary(bayesglm(`C>GMLA` ~ NROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 2,]))



summary(bayesglm(`C>GMLA` ~ SROH+AGE*SEX+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2))
summary(bayesglm(`C>GMLA` ~ SROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 1,]))
summary(bayesglm(`C>GMLA` ~ SROH+AGE+
                   euPC1+euPC2+euPC3+euPC4+euPC5+euPC6+euPC7+euPC8+euPC9+euPC10,
                 family="binomial", data=phenoloh2[phenoloh2$SEX == 2,]))


library(vtable)
sumtable(phenoloh2[,c('SEX', 'NROH', 'SROH')], group = 'SEX')
