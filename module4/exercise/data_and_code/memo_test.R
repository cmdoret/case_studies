setwd(dir = "~/Documents/Master/sem_1/Case_study/module4/exercise/data_and_code/") 
source('memo_test_methods.R', chdir=TRUE)

gam = read.table("data/GBM_gene_GAM.txt",header=T,row.names=1)

sorted.gene.alterations = sort(rowSums(gam))
sorted.sample.alterations = sort(colSums(gam))

# Q: For how many samples and genes do you have data?

gbm <- read.table("data/GBM_gene_GAM.txt",header=T, row.names = 1)
freq.gene <- rowSums(gbm)
library("car")
pca.genes<-prcomp(gbm)
PC.genes<-pca.genes[2]
PC.genes[[1]]

# Plot the alteration distributions
# per gene:
plot(sorted.gene.alterations,type='b',xaxt="n",xlab="",ylab="Number of Altered Cases",main="Alteration Distribution per Gene")
axis(1,at=c(1:nrow(gam)), labels=rownames(gam)[order(rowSums(gam))] ,las=2,cex.axis=0.4)
# per sample:
plot(sorted.sample.alterations,type='b',xaxt="n",xlab="Tumor Samples",ylab="Number of Altered Genes",main="Alteration Distribution per Gene")

# Q: Try a different way to plot this data (e.g. plotting density, barplot, percentage of alteration rather than actual number, etc.)

free_permut <- function(M, Q=100, events=NULL) {
  #	M <- as.matrix(M == TRUE);
  M <- as.matrix(M);
  R = sample(M)
  dim(R) = dim(M)
  rownames(R) = rownames(M)
  colnames(R) = colnames(M)
  return(R)
}
#image(free_permut(as.matrix(gbm)))
for(i in 1:100){points(rowSums(free_permut(gbm)),type="b",col=i)}


# Test the effect of different null models (i.e. expectations) on sample and gene distributions
# 1) unconstrained permutation
data.preserve.none = permute.preserve.none(gam)

# 2) Preserve gene distribution
data.preserve.gene = permute.preserve.gene(gam)

row_permut <- function(M){
  R <-matrix(nrow=nrow(M),ncol=ncol(M)) 
  R <-apply(M, FUN=sample,MARGIN = 1)
  R <- t(R)
  colnames(R)<-colnames(M)
  rownames(R)<-rownames(M)
  return(R)
}
for(i in 1:1){points(colSums(row_permut(gbm)),type="b",col=i)}


# 3) Preserve gene & sample distribution
data.preserve.all = permute.preserve.all(gam)

# Plot the real/observed and random/expected alteration distribution using different models:

plot(sort(rowSums(gam)),type='b',xaxt="n",xlab="",ylab="Number of Altered Cases")
axis(1,at=c(1:nrow(gam)), labels=rownames(gam)[order(rowSums(gam))] ,las=2,cex.axis=0.4)

points(rowSums(data.preserve.none)[order(rowSums(gam))],type='b',col="red")

plot(sort(colSums(gam)),type='b',xaxt="n",xlab="Tumor Samples",ylab="Number of Altered Genes")
points(colSums(data.preserve.none)[order(colSums(gam))],type='b',col="red")

# Q: Plot the alteration distributions according to the other models


# Define the modules to test for mutual exclusivity

# GBM modules
m1 = c("CDKN2A","CDK4","RB1")
m2 = c("PIK3CA","PTEN","PIK3R1")
m3 = c("EGFR","ERBB2","PDGFRA")
m4 = c("PIK3CA","PIK3R1","NF1")
m5 = c("RB1","PTEN","TP53")
m6 = c("CDKN2A","MDM2","MDM4","TP53")

modules = list(m1,m2,m3,m4,m5,m6)

# Q: based on the manuscript, define which modules to test (you should test positive and negative example)

# If you want to test if the data includes genes in your modules

(m1 %in% rownames(gam))


# TEST MODULES
# Here we are generating 1000 random matrices according to each null model (memoData)
# and each time we test our modules (testModules) and look at the results (summary)

data.preserve.none <- memoData(gam, Q=100, N=1000, permuteFun=permute.preserve.none, verbose=TRUE)
results = testModules(data.preserve.none, modules)
summary(results);

data.preserve.gene <- memoData(gam, Q=100, N=1000, permuteFun=permute.preserve.gene, verbose=TRUE)
results = testModules(data.preserve.gene, modules)
summary(results);

data.preserve.both <- memoData(gam, Q=100, N=1000, permuteFun=permute.preserve.all, verbose=TRUE)
results = testModules(data.preserve.both, modules)
summary(results);

# Q: make a table with modules as rows and p-value (-log10) as columns and plot these data

# UCEC modules
m1 = c("PIK3CA","PTEN","PIK3R1")
m2 = c("CTNNB1","KRAS","SOX17")

# CRC modules
m1 = c("BRAF","NRAS","KRAS")
m2 = c("IGF2","PIK3CA","PTEN","PIK3R1")

# STAD modules
m1 = c("RHOA","RASA1","KRAS")
m2 = c("PIK3CA","PTEN","PIK3R1")

# LGG modules
m1 = c("IDH1","IDH2")
m2 = c("CDKN2A","MDM4","TP53")
m3 = c("PIK3CA","PTEN","PIK3R1")

###
# Exo

STAD <- read.table("data/STAD_gene_GAM.txt",header=T,row.names=1)
freq.stad <- rowSums(STAD) 
plot(sort(freq.stad),type='b',xaxt="n",xlab="",ylab="Number of Altered Cases",main="Alteration Distribution per Gene")
axis(1,at=c(1:nrow(STAD)), labels=rownames(STAD)[order(rowSums(STAD))] ,las=2,cex.axis=0.4)
sample.stad <- colSums(STAD) 
plot(sort(sample.stad),type='b',xaxt="n",xlab="",ylab="Number of Altered genes",main="Alteration Distribution per Case")
axis(1,at=c(1:ncol(STAD)), labels=colnames(STAD)[order(colSums(STAD))] ,las=2,cex.axis=0.4)


mEBV <- c("PIK3CA","PDCD1LG2", "CDKN2A","CD274")
mMSI <- c("RASA1","KRAS","PIK3CA","PTEN")
mGS <- c("RHOA","CDH1")
mCIN <- c("TP53","ERBB2","KRAS")
mCTRL <- c("TP53","PTEN","CDH1","CDK6")
m1 = c("RHOA","RASA1","KRAS")
m2 = c("PIK3CA","PTEN","PIK3R1")
modules <- list(mEBV,mMSI,mGS,mCIN,mCTRL,m1,m2)


data.preserve.both <- memoData(STAD, Q=100, N=1000, permuteFun=permute.preserve.sample, verbose=TRUE)
results = testModules(data.preserve.both, modules)

out <-summary(results)
propre_sample <-cbind(moduleNames=c("Control","m1","Chrom. instability","Microsat. instability","EBV","Genomically stable","m2"),out)

data.preserve.both <- memoData(STAD, Q=100, N=1000, permuteFun=permute.preserve.sample, verbose=TRUE)
results = testModules(data.preserve.both, modules)

out <-summary(results)
propre_sample <-cbind(moduleNames=c("Control","m1","Chrom. instability","Microsat. instability","EBV","Genomically stable","m2"),out)

data.preserve.both <- memoData(STAD, Q=100, N=1000, permuteFun=permute.preserve.sample, verbose=TRUE)
results = testModules(data.preserve.both, modules)

out <-summary(results)
propre_sample <-cbind(moduleNames=c("Control","m1","Chrom. instability","Microsat. instability","EBV","Genomically stable","m2"),out)


