

SEED <- 444
set.seed(SEED)
Ngroup <- 10
effect.group <- 1 
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

scores <- Tstat * c(y0, y1)
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


# EXERCISE 1 & 2

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

















