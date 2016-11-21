

#########################################################################################

# This was a preparation 
# now with more genes say 2'000 and generalize further
# of which most or all will not be different btw the two groups 

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



#########################################################################################
#########################################################################################
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
boxplot(Tvalues, abs(Tvalues) )
points(1, Tvalues[1], pch=16, col="red", cex=2) 
points(2, abs(Tvalues[1]), pch=16, col="red", cex=2) 


# ***
selgenes <- (order(abs(Tvalues), decreasing = TRUE))[1] 
# *** EXERCISE 4: generalize for more than one gene ***
model.coefficients <- Tvalues[selgenes]
scores <- XS[, selgenes] * model.coefficients  
# *** EXERCISE 4: generalize for more than one gene ***
# ***

cutpoint <- ( mean( scores[YS==0] ) + mean( scores[YS==1] ) ) /  2 
pred <- sign(scores < cutpoint)
correct <- sum(pred == Y)
# return(correct / (2*Ngroup) )


# ***EXERCISE E4a:  GENERALIZE THE CODE FOR MORE THAN ONE GENE  ****


# ***EXERCISE E4b:  WRITE A FUNCTION FOR DATA GENERATION   ****
# that returns generated data as a list , ans is called like this: 
data <- datageneration.function(localSEED=SEED, Ngroup=Ngroup, Ngenes=NBclassifierFeatures, meanvalue=mean.x , sdvalue=sd.x , effect=effect.mean )

# ***EXERCISE E4c:  WRITE A FUNCTION FOR CONSTRUCTION OF THE CLASSIFIER   ****
# that returns a classifier as a list , ans is called like this: 
classifier  <-  classifier.function (mydata, NBgenes.selected, NBobservations)
# where the data is a list returned by the datageneration.function of E4b

# ***EXERCISE E4d:  PUT THINGS TOGETHER TO COMPUTE CORRECT CLASSIFICATION PROPORTION   ****
# test it for the case that the classifier uses two genes 


#########################################################################################
#########################################################################################

# EXAMPLE SOLUTION E4

# ***EXERCISE E4a:

NBgenes.selected <- 2
selgenes <- (order(abs(Tvalues), decreasing = TRUE))[1:NBgenes.selected ]  
model.coefficients <- Tvalues[selgenes]
scores <- XS[, selgenes] %*% model.coefficients  # test for 1 gene and multiple genes 


#########################
# ***EXERCISE E4b:

datageneration.function <-  function (localSEED, Ngroup=10, Ngenes=2000, meanvalue=0, sdvalue=1, effect=3, shiftgroup=1)
{
set.seed(localSEED)
Nlearn <- 2*Ngroup; 
XS <- matrix( rnorm(( Ngenes * Nlearn), mean=meanvalue, sd=sdvalue), ncol= Ngenes, nrow= Nlearn)  
YS <- rep(1,Nlearn); group1 <- sample(x=1:Nlearn, Ngroup); YS[group1] <- 0
sel <- (YS == shiftgroup;
XS[sel, 1 ] <- XS[sel, 1 ] + effect.mean;
# remark: always only the values of gene 1 are shifted, to easy keep track or reproduce in new data
datalist <- list (X=XS, Y=YS)
return (datalist)
}


#########################
# ***EXERCISE E4c:


classifier.function <-  function (mydata, NBgenes.selected, NBobservations)
{
Tvalues <- apply(X=mydata$X, MARGIN=2, FUN=testfunction.t, Z=mydata$Y)  
selgenes <- (order(abs(Tvalues), decreasing = TRUE))[1:NBgenes.selected ]   
model.coefficients <- Tvalues[selgenes]
if(NBgenes.selected > 1) {scores <- mydata$X[, selgenes] %*% model.coefficients}
	else {scores <- mydata$X[, selgenes] * model.coefficients}
cutpoint <- ( mean( scores[mydata$Y==0] ) + mean( scores[mydata$Y==1] ) ) /  2 
pred <- sign(scores < cutpoint)
corr.prop <- sum(pred == mydata$Y) / NBobservations
classifier <- list(features=selgenes, coefficients=model.coefficients, cutoff=cutpoint, predicted=pred, correct.proportion=corr.prop)
return (classifier)
}




#########################
# ***EXERCISE E4d

# parameters
Ngroup <- 10; 
NBgenesMeasured <- 2000; # genes in data
mean.x <- 0; sd.x <- 1; effect.mean <- 3; effect.group <- 1 
SEED <- 444

set.seed(SEED)
data <- datageneration.function(localSEED=SEED, Ngroup=Ngroup, Ngenes=NBgenesMeasured, meanvalue=mean.x, sdvalue=sd.x, effect=effect.mean )
str(data)

# test multiple situations :

NBclassifierFeatures <-  1 
classifier <- classifier.function(mydata=data, NBgenes=NBclassifierFeatures, NBobservations=2*Ngroup)
str(classifier)
classifier$correct.proportion

NBclassifierFeatures <-  1 
classifier <- classifier.function(mydata=data, NBgenes=NBclassifierFeatures, NBobservations=2*Ngroup)
str(classifier)
classifier$correct.proportion


effect.mean <- 0
SEED <- 444
set.seed(SEED)
data <- datageneration.function(localSEED=SEED, Ngroup=Ngroup, Ngenes=NBgenesMeasured, meanvalue=mean.x, sdvalue=sd.x, effect=effect.mean )
str(data)

NBclassifierFeatures <-  1 
classifier <- classifier.function(mydata=data, NBgenes=NBclassifierFeatures, NBobservations=2*Ngroup)
str(classifier)
classifier$correct.proportion

NBclassifierFeatures <-  2 
classifier <- classifier.function(mydata=data, NBgenes.selected=NBclassifierFeatures, NBobservations=2*Ngroup)
str(classifier)
classifier$correct.proportion

NBclassifierFeatures <-  10 
classifier <- classifier.function(mydata=data, NBgenes.selected=NBclassifierFeatures, NBobservations=2*Ngroup)
str(classifier)
classifier$correct.proportion





#########################################################################################
#########################################################################################
# needed for the report 

# ***EXERCISE E5:  REPRODUCE THE ANALYIS DONE IN THE PAPER AND 
		GENERATE A VERSION OF THEIR FIGURE  ****
(If the simulation is too slow for 2'000 repetitions and 6'000 genes reduce
the number of simulations in the lab and run the fuller version later at home)
***EXERCISE E5a: write code for one full cross-validation 
***EXERCISE E5b: write code for one incomplete cross-validation 
***EXERCISE E5c: do an appropriate number of repetitions to collect results for
resubstitution, full and incomplete cross-validation  



#########################################################################################
#########################################################################################
# supplementary suggested additional exercises
   

# ***EXERCISE E6:  Scenario: there is one informative gene and only one
gene is selected for use in the classifier: Generate and interpret 
the analogous figure (set effect.mean = 2) ****
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
