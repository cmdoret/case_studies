

#########################################################################################
#########################################################################################


# ***EXERCISE E5:  REPRODUCE THE ANALYIS DONE IN THE PAPER AND 
		GENERATE A VERSION OF THEIR FIGURE  ****
(If the simulation is too slow for 2'000 repetitions and 6'000 genes reduce
the number of simulations in the lab and run the fuller version later at home)
***EXERCISE E5a: write code for one full cross-validation 
***EXERCISE E5b: write code for one incomplete cross-validation 
***EXERCISE E5c: do an appropriate number of repetitions to collect results for
resubstitution, full and incomplete cross-validation  


# ***EXERCISE E6:  Scenario: there is one informative gene and only one
gene is selected for use in the classifier: Generate and interpret 
the analogous figure (set effect.mean = 2.5) ****
(An illustration of "Winner's bias") 
Extension E6b: test the classifiers also on an independent second (large)
dataset and compare the performances in selected cases.


# ***EXERCISE E7:  In the situation of the exercise E6: ****
How large should the sample size of the dataset be so that one can build a
classifier that reaches at least 90% of the maximal performance in at
least 90% of the cases ?
(An illustration of a "Power Study") 



#########################################################################################
#########################################################################################
# EXAMPLE SOLUTION E5a

# What we did before was learning and testing in one set
# in crossvalidation we separate the two things 

Ngroup <- 10; Nlearn <- 2*Ngroup; 
NBgenesMeasured <- 6000; # genes in data
NBclassifierFeatures <- 10 
mean.x <- 0; sd.x <- 1; effect.mean <- 0

SEED <- 444
data <- datageneration.function(localSEED=SEED, Ngroup=Ngroup, Ngenes=NBgenesMeasured, meanvalue=mean.x, sdvalue=sd.x, effect=effect.mean )
str(data)

# learn each time after removing one case and test on this one alone 
# use for example a for loop for this leave-one-out CV

# one full LOO crossvalidation :
pred.loo <- rep(NA, Nlearn)
for ( loo in c(1:Nlearn) )
{
learndata <- data
learndata$X <- learndata$X[-loo, ]
learndata$Y <- learndata$Y[-loo]
classifier.loo <- classifier.function(mydata=learndata, NBgenes.selected=NBclassifierFeatures, NBobservations=(Nlearn-1) )
score.loo <- as.numeric( data$X[ loo, classifier.loo$features] %*% classifier.loo$coefficients )
pred.loo[loo] <- sign(score.loo <  classifier.loo$cutoff)
}
corr.prop.loo <- sum(pred.loo == data$Y) / Nlearn

# measure how long does this take ?: can repeat it many times without a long waiting time ? 


#########################
# ***EXERCISE E5b:

# one incompelte LOO crossvalidation : gene selected on all data
# then coefficients with LOO using exactly only those genes 

SEED <- 444
data <- datageneration.function(localSEED=SEED, Ngroup=Ngroup, Ngenes=NBgenesMeasured, meanvalue=mean.x, sdvalue=sd.x, effect=effect.mean )

# select genes on all data :
classifier <- classifier.function(mydata=data, NBgenes.selected=NBclassifierFeatures, NBobservations=Nlearn) # this is the original classifier as for resubstitution, trained on all data
selectedgenes <- classifier$features

# one LOO crossvalidation loop:
pred.iCV <- rep(NA, Nlearn)
for ( loo in c(1:Nlearn) )
{
learndata <- data
testdata.loo <- learndata$X[loo, selectedgenes]
learndata$X <- learndata$X[-loo, selectedgenes]
learndata$Y <- learndata$Y[-loo]
classifier.icv <- classifier.function(mydata=learndata, NBgenes.selected=NBclassifierFeatures, NBobservations=(Nlearn-1) )
score.iCV <- as.numeric( testdata.loo[classifier.icv$features] %*% classifier.icv$coefficients)
pred.iCV[loo] <- sign(score.iCV <  classifier.icv$cutoff)
}
corr.prop.iCV <- sum(pred.iCV == data$Y) / Nlearn

# measure how long does this take ?: can repeat it many times without a long waiting time ? 



#########################
# ***EXERCISE E5c:


# parameters
Ngroup <- 10; Nlearn <- 2*Ngroup; 
NBgenesMeasured <- 6000; # genes in data, can use about 1000 for an alternate faster analysis 
NBclassifierFeatures <- 10
mean.x <- 0; sd.x <- 1 
effect.mean <- 0; effect.group <- 1 

# other variables
repetitions  <- 20  # 200 
corr.prop.resub <- rep(NA, repetitions)
corr.prop.loo <- rep(NA, repetitions)
corr.prop.iCV <- rep(NA, repetitions)


SEED <- 443
for (i in  c(1:repetitions) )
{
SEED <- SEED+1; set.seed(SEED);

cat("\n data")
data <- datageneration.function(localSEED=SEED, Ngroup=Ngroup, Ngenes=NBgenesMeasured, meanvalue=mean.x , sdvalue=sd.x , effect=effect.mean )

cat("\n resubstitution case")
classifier <- classifier.function(mydata=data, NBgenes.selected=NBclassifierFeatures, NBobservations=Nlearn)
corr.prop.resub[i] <- classifier$correct.proportion
selectedgenes <- classifier$features


cat("\n incomplete LOO crossvalidation")
# incomplete LOO crossvalidation : # could put this section into a function 
pred.iCV <- rep(NA, Nlearn)
for ( loo in c(1:Nlearn) )
{
learndata <- data
testdata.loo <- learndata$X[loo, selectedgenes]
learndata$X <- learndata$X[-loo, selectedgenes]
learndata$Y <- learndata$Y[-loo]
classifier.icv <- classifier.function(mydata=learndata, NBgenes.selected=NBclassifierFeatures, NBobservations=(Nlearn-1) )
score.iCV <- as.numeric( testdata.loo[classifier.icv$features] %*% classifier.icv$coefficients)
pred.iCV[loo] <- sign(score.iCV <  classifier.icv$cutoff)
}
corr.prop.iCV[i] <- sum(pred.iCV == data$Y) / Nlearn


cat("\n full LOO crossvalidation")
# full LOO crossvalidation : # could put this section into a function 
pred.loo <- rep(NA, Nlearn)
for ( loo in c(1:Nlearn) )
{
learndata <- data
learndata$X <- learndata$X[-loo, ]
learndata$Y <- learndata$Y[-loo]
classifier.loo <- classifier.function(mydata=learndata, NBgenes.selected=NBclassifierFeatures, NBobservations=(Nlearn-1) )
score.loo <- as.numeric( data$X[ loo, classifier.loo$features] %*% classifier.loo$coefficients )
pred.loo[loo] <- sign(score.loo <  classifier.loo$cutoff)
if ( (loo %% 5) == 0)  {cat("\n loo=",loo)}
}
corr.prop.loo[i] <- sum(pred.loo == data$Y) / Nlearn


if ( (i %% 10) == 0)  {cat("\n\n repetitions i=",i)}
}


table(corr.prop.resub) 
table(corr.prop.loo) 
table(corr.prop.iCV) 


# 6000 genes, 200 repetitions 

> table(corr.prop.resub) 
corr.prop.resub
0.95    1 
   6  194 
> table(corr.prop.loo) 
corr.prop.loo
0.05  0.1 0.15  0.2 0.25  0.3 0.35  0.4 0.45  0.5 0.55  0.6 0.65  0.7 0.75  0.8 0.85  0.9 
   3    9    7   10   19   14   13   15   13   19   15   15   14   12   11    7    2    2 
> table(corr.prop.iCV) 
corr.prop.iCV
 0.9 0.95    1 
   1   25  174 



summary(corr.prop.resub) 
summary(corr.prop.loo) 
summary(corr.prop.iCV) 


> summary(corr.prop.resub) 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.9500  1.0000  1.0000  0.9925  1.0000  1.0000 
> summary(corr.prop.loo) 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1000  0.3000  0.5000  0.4325  0.5000  0.8500 
> summary(corr.prop.iCV) 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.95    0.95    1.00    0.98    1.00    1.00 
> 


boxplot( list ( corr.prop.resub=corr.prop.resub, 
corr.prop.loo=corr.prop.loo, corr.prop.iCV=corr.prop.iCV ))

results.1 <- list ( corr.prop.resub=corr.prop.resub, 
corr.prop.loo=corr.prop.loo, corr.prop.iCV=corr.prop.iCV )




Myname <- "XY"
pdf(paste (Myname , "exercise5c.pdf"))

hist( results.1$corr.prop.resub, ylim=c(0,repetitions), xlim=c(-0.1,1.1), breaks=-0.025+c(0:21)*0.05, col=rgb(0,0,0,0.7), main="", xlab="proportion correct classifications")
par(new=TRUE)
hist( results.1$corr.prop.iCV, ylim=c(0,repetitions), xlim=c(-0.1,1.1), breaks=-0.025+c(0:21)*0.05 , 
col=rgb(1,1,0,0.5), main="", xlab="")
par(new=TRUE)
hist( results.1$corr.prop.loo, ylim=c(0,repetitions), xlim=c(-0.1,1.1), breaks=-0.025+c(0:21)*0.05 , 
col=rgb(1,0,0,1), main="", xlab="")
title("Histogram proportions correct")

shifts <- 10*jitter(rep(0,repetitions))
boxplot( results.1, col= c(rgb(0,0,0,1), rgb(1,0,0,1), rgb(1,1,0,1) ) )
points( shifts + rep(1, repetitions ), results.1$corr.prop.resub, pch=15, col="blue", cex=0.3)
points( shifts + rep(2, repetitions ), results.1$corr.prop.loo, pch=15, col="blue", cex=0.3)
points( shifts + rep(3, repetitions ), results.1$corr.prop.iCV, pch=15, col="blue", cex=0.3)
title("Boxplot proportions correct")

dev.off()



##########

results.1 <- list ( corr.prop.resub=corr.prop.resub, 
corr.prop.loo=corr.prop.loo, corr.prop.iCV=corr.prop.iCV )

boxplot(results.1)

#########################################################################################
#########################################################################################









