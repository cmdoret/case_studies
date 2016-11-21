#' @description This script contains all exercises for the second part of module 2 for the course "case studies in bioinformatics"
#' @author Cyril Matthey-Doret
# Mon Nov 21 09:33:44 2016 ------------------------------



# GENERATE ARTIFICIAL DATA , first random and equal for groups then add a shift for informative genes 

SEED <- 444 
set.seed(SEED)

# parameters
Ngroup <- 10; Nlearn <- 2*Ngroup; 
Y <- c( rep(0,Ngroup), rep(1,Ngroup) ) 
Ngenes <- 2000;
# VALUE DISTRIBUTION
mean.x <- 0; sd.x <- 1 
effect.mean <- 3

# SIMULATED GENE EXPRESSION VALUES IN A MATRIX X
# variables /FEATURES IN COLUMNS, observations / SAMPLES IN ROWS
XS <- matrix( rnorm(( Ngenes * Nlearn), mean=mean.x, sd=sd.x), ncol= Ngenes, nrow= Nlearn)  
dim(XS)

#  GENERATE CLASS LABELS (JUST FOR PRINCIPLE RANDOM ORDER LABELS AS IN PAPER)
YS <- rep(1,Nlearn); 
group1 <- sample(x=1:Nlearn, Ngroup)
YS[group1] <- 0
table(YS)

# LEARNING - SELECTION OF FEATURES AND COEFFICIENTS FOR THE CLASSIFIER 
# SELECTION OF FEATURES : top genes in a t-test

Tvalues <- apply(X=XS, MARGIN=2, FUN=testfunction.t, Z=YS)  
# where YS are the group identifiers and testfunction.t is a function we write
# and that returns the t-value for each gene 
# ***EXERCISE E3:  WRITE THE FUNCTION testfunction.t ****


boxplot(Tvalues, abs(Tvalues) )



#########################################################################################
#########################################################################################
# EXAMPLE SOLUTION E3

testfunction.t <- function(u,Z)  
  # receiving from caller two parameters, a vector of values u and a vector of labels Z
  # these can be provided by "apply", the first by its main mechanism, the second as a named parameter 
{
  # this can be run for each gene in a matrix using the function call inside the apply function 
  u1 <- u[Z==0]     #  values for group "0" 
  u2 <- u[Z==1]   	#  values for group "1"
  st <- t.test (x=u1, y=u2)  #  object returned from the t-test function, positive when mean u1 HIGHER 
  return (st$statistic)  #  return to caller the t-stat value from the object st 
}

# ADD DISCRIMINATIVE GENE 

effect.group <- 1 
effect.mean <- 3
NBgenes.selected <- 1

# add  "Effect" always to group 1, so can reproduce in indep test set without having 
# to store the information on this direction of the effect 
# for simplicity take always the gene "1" 
sel <- (YS == effect.group);
XS[sel, 1 ] <- XS[sel, 1 ] + effect.mean;

Tvalues <- apply(X=XS, MARGIN=2, FUN=testfunction.t, Z=YS)  
boxplot(Tvalues)
boxplot( list(Tvalues=Tvalues, absTvalues=abs(Tvalues) ) )
points(1, Tvalues[1], pch=16, col="red", cex=2) 
points(2, abs(Tvalues[1]), pch=16, col="red", cex=2) 


# ***
selgenes <- (order(abs(Tvalues), decreasing = TRUE))[1] 
# *** EXERCISE 4: generalize for more than one gene ***
model.coefficients <- Tvalues[selgenes]
scores <- XS[, selgenes] * model.coefficients  
# *** EXERCISE 4: generalize for more than one gene ***
# ***
Nselected <- 2
selgenes <- (order(abs(Tvalues), decreasing = TRUE))[1:Nselected] 
# *** EXERCISE 4: generalize for more than one gene ***
model.coefficients <- Tvalues[selgenes]
scores <- XS[, selgenes] * model.coefficients  

cutpoint <- ( mean( scores[YS==0] ) + mean( scores[YS==1] ) ) /  2 
pred <- sign(scores < cutpoint)
correct <- sum(pred == Y)
# return(correct / (2*Ngroup) )


# ***EXERCISE E4a:  GENERALIZE THE CODE FOR MORE THAN ONE GENE  ****


# ***EXERCISE E4b:  WRITE A FUNCTION FOR DATA GENERATION   ****

datagen <- function(Ngroup, Ngenes, meanvalue, sdvalue, effect.mean, Nselected){
  Nlearn <- Ngroup*2
  XS <- matrix(rnorm(mean = meanvalue, sd = sdvalue,n = (Ngenes*Nlearn)),nrow = Nlearn,byrow = T)
  YS <- rep(1, Nlearn);group1 <- sample(x=1:Nlearn,Ngroup);YS[group1] <- 0
  sel <- (YS == 1)
  XS[sel,1:Nselected] <- XS[sel,1:Nselected]+ effect.mean
  data <- list(X=XS,Y=YS)
  return(data)
}

# that returns generated data as a list , ans is called like this: 
data <- datagen(localSEED=4, Ngroup=200, Ngenes=9, meanvalue=2 , sdvalue=0.3 , effect.mean=3, Nselected=1 )

# ***EXERCISE E4c:  WRITE A FUNCTION FOR CONSTRUCTION OF THE CLASSIFIER   ****
# that returns a classifier as a list , ans is called like this: 


testfunction.t <- function(u,Z)  
  # receiving from caller two parameters, a vector of values u and a vector of labels Z
  # these can be provided by "apply", the first by its main mechanism, the second as a named parameter 
{
  # this can be run for each gene in a matrix using the function call inside the apply function 
  u1 <- u[Z==0]     #  values for group "0" 
  u2 <- u[Z==1]   	#  values for group "1"
  st <- t.test (x=u1, y=u2)  #  object returned from the t-test function, positive when mean u1 HIGHER 
  return (st$statistic)  #  return to caller the t-stat value from the object st 
}

classifier.function <- function(data, Ngenes, Nobs){
  Tvalues <- apply(X=data[["X"]], MARGIN=2, FUN=testfunction.t, Z=data[["Y"]])  
  #par(mfrow=c(1,2))
  #boxplot(Tvalues)
  #boxplot( list(Tvalues=Tvalues, absTvalues=abs(Tvalues) ) )
  #points(rep(1,Ngenes), Tvalues[1:Ngenes], pch=16, col="red", cex=2) 
  #points(rep(2,Ngenes), abs(Tvalues[1:Ngenes]), pch=16, col="red", cex=2) 
  selgenes <- (order(abs(Tvalues), decreasing = TRUE))[Ngenes] 
  # *** EXERCISE 4: generalize for more than one gene ***
  model.coefficients <- Tvalues[selgenes]
  scores <- data[["XS"]][, selgenes] * model.coefficients
  # *** EXERCISE 4: generalize for more than one gene ***
  model.coefficients <- Tvalues[selgenes]
  scores <- data[["X"]][, selgenes] * model.coefficients  
  cutpoint <- ( mean( scores[data[["Y"]]==0] ) + mean( scores[data[["Y"]]==1] ) ) /  2 # Point between distributions
  pred <- sign(scores < cutpoint)
  correct <- sum(pred == data[["Y"]])/ (nrow(data[["X"]]))
  return(list(correct=correct,coef=model.coefficients,scores=scores,features=selgenes,cutoff=cutpoint))
}

classifier  <-  classifier.function (data, Ngenes=1, Nobs=400)
# where the data is a list returned by the datageneration.function of E4b

# ***EXERCISE E4d:  PUT THINGS TOGETHER TO COMPUTE CORRECT CLASSIFICATION PROPORTION   ****

  # test it for the case that the classifier uses two genes 

# ***EXERCISE E5:  REPRODUCE THE ANALYIS DONE IN THE PAPER AND GENERATE A VERSION OF THEIR FIGURE  ****
# (If the simulation is too slow for 2'000 repetitions and 6'000 genes reduce
# the number of simulations in the lab and run the fuller version later at home)
# ***EXERCISE E5a: write code for one full cross-validation 



benchmark <- data.frame()
for(t in 10:300){
  corr.prop.loo <- rep(0,length(sim[1,1,]))
  xsim <- replicate(n = t,datagen(Ngroup=10, Ngenes=60, meanvalue=2 , sdvalue=0.3 , effect.mean=3, Nselected=1 )[["X"]])
  ysim <- replicate(n = t,datagen(Ngroup=10, Ngenes=60, meanvalue=2 , sdvalue=0.3 , effect.mean=3, Nselected=1 )[["Y"]])
  for(tmp in 1:length(sim[1,1,])){
    pred.learn <- rep(0,19)
    tmp.data <- list(X=xsim[,,tmp],Y=ysim[,tmp])
    for(i in 1:length(tmp.data[["X"]][,1])){
      learn.data <- list(X=tmp.data[["X"]][-i,],Y=tmp.data[["Y"]][-i])
      classifier.learn  <-  classifier.function (learn.data, Ngenes=2, Nobs=length(learn.data$X[,1]))
      score.learn <- as.numeric(data$X[i,classifier.learn$features] %*% classifier.learn$coef)
      pred.learn[i] <- sign(score.learn <  classifier.learn$cutoff)
      print(paste(sep=" ", "We are at itereation nr", i, "of simulation nr", tmp))
    }
    corr.prop.loo[tmp] <- sum(pred.learn == tmp.data$Y)/length(xsim[,1,1])
  }
  benchmark$as.character(t) <- corr.prop.loo
}

benchmark <- data.frame()
for(t in 10:300){
  corr.prop.loo <- rep(0,length(sim[1,1,]))
  xsim <- replicate(n = t,datagen(Ngroup=10, Ngenes=60, meanvalue=2 , sdvalue=0.3 , effect.mean=3, Nselected=1 )[["X"]])
  ysim <- replicate(n = t,datagen(Ngroup=10, Ngenes=60, meanvalue=2 , sdvalue=0.3 , effect.mean=3, Nselected=1 )[["Y"]])
  for(tmp in 1:length(sim[1,1,])){
    pred.learn <- rep(0,19)
    tmp.data <- list(X=xsim[,,tmp],Y=ysim[,tmp])
    for(i in 1:length(tmp.data[["X"]][,1])){
      learn.data <- list(X=tmp.data[["X"]][-i,],Y=tmp.data[["Y"]][-i])
      lapply(1:length(data[["X"]][,1]), function(d){
        return(list(X=data[["X"]][-d,], Y=data[["Y"]][-d]))})
      classifier.learn  <-  classifier.function (learn.data, Ngenes=2, Nobs=length(learn.data$X[,1]))
      score.learn <- as.numeric(data$X[i,classifier.learn$features] %*% classifier.learn$coef)
      pred.learn[i] <- sign(score.learn <  classifier.learn$cutoff)
      print(paste(sep=" ", "We are at itereation nr", i, "of simulation nr", tmp))
    }
    corr.prop.loo[tmp] <- sum(pred.learn == tmp.data$Y)/length(xsim[,1,1])
  }
  benchmark$as.character(t) <- corr.prop.loo
}

# ***EXERCISE E5b: write code for one incomplete cross-validation 
# ***EXERCISE E5c: do an appropriate number of repetitions to collect results for
# resubstitution, full and incomplete cross-validation  

