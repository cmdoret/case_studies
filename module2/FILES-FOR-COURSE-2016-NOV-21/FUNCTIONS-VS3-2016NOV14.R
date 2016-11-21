
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




datageneration.function <-  function (localSEED, Ngroup=10, Ngenes=2000, meanvalue=0, sdvalue=1, effect=3, shiftgroup=1)
{
set.seed(localSEED)
Nlearn <- 2*Ngroup; 
XS <- matrix( rnorm(( Ngenes * Nlearn), mean=meanvalue, sd=sdvalue), ncol= Ngenes, nrow= Nlearn)  
YS <- rep(1,Nlearn); group1 <- sample(x=1:Nlearn, Ngroup); YS[group1] <- 0
sel <- (YS == shiftgroup);
XS[sel, 1 ] <- XS[sel, 1 ] + effect;
# remark: always only the values of gene 1 are shifted, to easy keep track or reproduce in new data
datalist <- list (X=XS, Y=YS)
return (datalist)
}



classifier.function <-  function (mydata, NBgenes.selected, NBobservations)
{
Tvalues <- apply(X=mydata$X, MARGIN=2, FUN=testfunction.t, Z=mydata$Y)  
selgenes <-  (order(abs(Tvalues), decreasing = TRUE))[1:NBgenes.selected ] 
# selgenes <- sort( (order(abs(Tvalues), decreasing = TRUE))[1:NBgenes.selected ]) 
# sorting for use in iCV  ? 
model.coefficients <- Tvalues[selgenes]
if(NBgenes.selected > 1) {scores <- mydata$X[, selgenes] %*% model.coefficients}
	else {scores <- mydata$X[, selgenes] * model.coefficients}
cutpoint <- ( mean( scores[mydata$Y==0] ) + mean( scores[mydata$Y==1] ) ) /  2 
pred <- sign(scores < cutpoint)
corr.prop <- sum(pred == mydata$Y) / NBobservations
classifier <- list(features=selgenes, coefficients=model.coefficients, cutoff=cutpoint, predicted=pred, correct.proportion=corr.prop)
return (classifier)
}

























