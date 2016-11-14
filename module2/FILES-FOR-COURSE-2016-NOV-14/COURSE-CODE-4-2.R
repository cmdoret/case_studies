

SEED <- 444
set.seed(SEED)
Ngroup <- 10
effect.mean <- 3
y0 <- rnorm(Ngroup , 0, 1) 
y1 <- rnorm(Ngroup , effect.mean, 1) 
# how strong the separation ? How to visualize ? 


boxplot(y1,y2)
# other graphical help ?

points( rep(1, Ngroup ), y0, pch=16, col="blue", cex=1.5)
points( rep(2, Ngroup ), y1, pch=15, col="red", cex=1.5)


# show the distribution of the simulated data
x <- seq(1,12, by=0.01) 
plot( x, dnorm(x, 5, 1), pch=".", cex=3)
points( x, dnorm(x, 8, 1) , pch=".", col="red", cex=3)

# how discriminative is this gene ?
# how good is a classifier that uses this gene ?

t.test(y0, y1)
tt <- t.test(y0, y1)
str(tt)
Tstat <- tt$statistic
Tstat
# negative when higher values in second group y1
print(tt)


# Correct classification rate ?
mean(y0)
mean(y1)
cutpoint <- (mean(y0) + mean(y1)) / 2 
cutpoint

y <- c(y0, y1)
# group labels
Y <- c( rep(0,Ngroup), rep(1,Ngroup) ) 
pred <- sign(y > cutpoint )
pred  
correct <- sum(pred == Y)
correct
correct.ratio <- correct / (2*Ngroup)
correct.ratio 

# How to generalize ?

# generalize in relation to gene "direction": 
# what if gene would be higher in group 0 ?

# discriminatory score with direction
# direction information for example in the t-statistics 

scores <- Tstat * y
scores
scores[Y==0]
scores[Y==1]
mean(scores[Y==0])
mean(scores[Y==1])

cutpoint <- (mean(scores[Y==0]) + mean(scores[Y==1])) / 2 
cutpoint
pred <- sign(scores < cutpoint )
pred 
correct <- sum(pred == Y)
correct
correct.ratio <- correct / (2*Ngroup)
correct.ratio 


#########################################################################################
#########################################################################################


# EXERCISE 1+2

# E1a. How un-stable is this 90% correct classification rate across example datasets ?
Take multiple such datasets, collect results and compare results in a table  

# E1b. ADDITIONS to E1. 
# Show results in a histogram and store as a pdf file 
# Write a function that does the computation for a (single case), so that the “core-code" for the 	# simulation is simplified, most "details" are handled by the function, and "key parameters"
# are passed to the function but defined in the “core-code"


# E2a. How accurate is this 90% correct classification rate obtained in the example ?
Take a large dataset  
How large ?

# E2b. ADDITIONS to E2. 
# Use the function (modified if needed) from the previous example
# DO 3 simulations IN PARALLEL PLOT RUNNING PARTIAL RESULT  


#########################################################################################

# EXAMPLE SOLUTION E1A

# BY SIMULATION, TAKING REPEATED DATA SETS 

# parameters
Ngroup <- 10; effect.mean <- 3
Y <- c( rep(0,Ngroup), rep(1,Ngroup) ) 
repetitions <- 1000

# variable / preparations 
results <- rep(NA,repetitions)

# do the simulation
SEED <- 444
for (i in  c(1:repetitions) )
{
SEED <- SEED+1; set.seed(SEED);
y0 <- rnorm(Ngroup , 0, 1) 
y1 <- rnorm(Ngroup , effect.mean, 1) 
y <- c(y0,y1);
# how strong the separation ? How to visualize ? 
tt <- t.test(y0, y1)
Tstat <- tt$statistic
scores <- Tstat * y
cutpoint <- (mean(scores[Y==0]) + mean(scores[Y==1])) / 2 
pred <- sign(scores < cutpoint)
correct <- sum(pred == Y)
results[i] <- correct
if ( (i %% 100) == 0)  {cat("\n i=",i)}
}

# analyze the simulation results
summary(results)
sd(results)
results.proportion <- results / (2*Ngroup)
summary(results.proportion )
table(results.proportion )
boxplot(results.proportion)
hist(results.proportion, xlim=c(-0.1,1.1), breaks=-0.025+c(0:21)*0.05 )
Myname <- "XY"
pdf(paste (Myname , "exercise1b.pdf"))
hist(results.proportion, xlim=c(-0.1,1.1), breaks=-0.025+c(0:21)*0.05 )
dev.off()


######################
# EXAMPLE SOLUTION E1B

# write a function to simplify the loop etc.
determine.perf.fu <- function(localSEED, Ngroup=10, effect.mean=3)
{
set.seed(localSEED);
Y <- c( rep(0,Ngroup), rep(1,Ngroup) );
y0 <- rnorm(Ngroup , 0, 1); y1 <- rnorm(Ngroup , effect.mean, 1); y <- c(y0,y1);
# how strong the separation ? How to visualize ? 
tt <- t.test(y0, y1); Tstat <- tt$statistic
scores <- Tstat * y; 
cutpoint <- (mean(scores[Y==0]) + mean(scores[Y==1])) / 2 
pred <- sign(scores < cutpoint)
correct <- sum(pred == Y)
return(correct)
}


# rewrite the simplified main loop etc.

# parameters
repetitions <- 1000

# variable / preparations 
results2 <- rep(NA,repetitions)

# do the simulation
SEED <- 444
for (i in  c(1:repetitions) )
{
SEED <- SEED+1; 
results2[i]  <- determine.perf.fu (SEED)
if ( (i %% 100) == 0)  {cat("\n i=",i)}
}



# analyze the simulation results
results.proportion <- results2 / (2*Ngroup)
summary(results.proportion )
table(results.proportion )


#########################################################################################

# EXAMPLE SOLUTION E2a

# BY SIMULATION, TAKING REPEATED DATA SETS, USING A FUNCTION 
# ADDITIONS   
# 3 IN PARALLEL
# PLOT RUNNING PARTIAL RESULTS DURING THE SIMULATION  



# MODIFY THE PREVIOUS FUNCTION  

determine.perf.fu2 <- function(localSEED, Ngroup=10, effect.mean=3)
{
set.seed(localSEED);
Y <- c( rep(0,Ngroup), rep(1,Ngroup) );
y0 <- rnorm(Ngroup , 0, 1); y1 <- rnorm(Ngroup , effect.mean, 1); y <- c(y0,y1);
# how strong the separation ? How to visualize ? 
tt <- t.test(y0, y1); Tstat <- tt$statistic
scores <- Tstat * y; 
cutpoint <- (mean(scores[Y==0]) + mean(scores[Y==1])) / 2 
pred <- sign(scores < cutpoint)
correct <- sum(pred == Y)
return(correct / (2*Ngroup) )
}

 
 
 
######################
# EXAMPLE SOLUTION E2a

# parameters
effect.mean <- 3
NBparallel <- 5

# variable / preparations 
group.sizes <- c(10,20,40,80,160,320,640,1280,2560,5120,10240,20480)
resultsmat <- matrix (NA, ncol=length(group.sizes), nrow=NBparallel)
plot( x=0, y=0, type="n", xlim=c(0, length(group.sizes)), ylim=c(0,1), xlab= "sample size case", ylab= "proportion correct") 


# do the simulation
SEED <- 444
i <- 0
for (Ngroup in  group.sizes )
{
i <- i+1; 
for ( k in 1:NBparallel)
{
SEED <- SEED+1;
resultsmat[k,i] <- determine.perf.fu2 (SEED, Ngroup=Ngroup)
points(x=i+k/20, y=resultsmat[k,i], pch=16, col=k)
}
if ( (i %% 1) == 0)  {cat("\n i=",i)}
}

summary(resultsmat)
apply (resultsmat,2,range)


plot( x=0, y=0, type="n", xlim=c(0, length(group.sizes)), ylim=c(0.8,1), 
xlab= "sample size case", ylab= "proportion correct") 

for (i in c(1:NBparallel) )
{
points(x= c(1: length(group.sizes)), y=resultsmat[i,] , pch=16, col=i+1)
}

for (i in c(1:NBparallel) )
{
lines (x= c(1: length(group.sizes)), y=resultsmat[i,], col=i+1)
}

# convergence ? 



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
classifier.function <-  function (mydata, NBgenes, NBobservations)
# where the data is a list returned by the datageneration.function of E4b

# ***EXERCISE E4d:  PUT THINGS TOGETHER TO COMPUTE CORRECT CLASSIFICATION PROPORTION   ****
# test it for the case that the classifier uses two genes 

#########################################################################################
#########################################################################################





