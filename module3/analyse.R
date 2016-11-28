library(seqinr)
library(stats)
library(ape)

#################
#Install packages
#################
#Run  the commands below if you need to install these packages

#install.packages("seqinr")
#install.packages("ape")


##############
#Load the data
##############
setwd("/home/cyril/Documents/Master/sem_1/Case_study/module3/")
files.human <- list.files("PDZLigands/human", pattern="*.fa", full.names=TRUE)
files.worm <- list.files("PDZLigands/worm", pattern="*.fa", full.names=TRUE)


#########################
#Get the PDZ domain names
#########################

#Remove the extension (".pep.fa") by splitting the string according to the delimiter ".pep"
nm <- c()
for (f in files.human){
  nm <- c(nm, strsplit(f, split=".pep", fixed=T )[[1]][1])
}
#Remove the path ("PDZLigands/human/")
PDZ.human <- c()
for (x in nm){
  PDZ.human <- c(PDZ.human, strsplit(x, split="/", fixed=T )[[1]][3])
}

#Do the same for the worm data. Use a more concise approach with the lapply() function
nm <- lapply(files.worm, function(x){strsplit(x, split=".pep", fixed=T )[[1]][1] })
PDZ.worm <- lapply(nm, function(x){strsplit(x, split="/", fixed=T )[[1]][3] })
PDZ.worm <- unlist(PDZ.worm)   #This is because the output of lapply() is a list and not a vector

#This would be an even more concise version for loading the PDZ names in Worm
#PDZ.worm <- unlist(lapply(files.worm, function(x){strsplit(  strsplit(x, ".pep")[[1]][1], "/", fixed=T )[[1]][3] }))

PDZ.all <- c(PDZ.human, PDZ.worm)
le <- length(PDZ.all)

PDZ.org <- c(rep("H", length(PDZ.human)), rep("W", length(PDZ.worm)))
names(PDZ.org) <- PDZ.all


####################
#Get the PDZ classes
####################
PDZ.class <- read.table("PDZclass.txt", skip=1)
class2 <- as.numeric(PDZ.class[,2])
names(class2) <- PDZ.all
class7 <- as.numeric(PDZ.class[,3])
names(class7) <- PDZ.all
class16 <- as.character(PDZ.class[,4])
names(class16) <- PDZ.all

#########################
#Load the background freq 
#########################
#This is only useful for the clustering in part 6
b <- read.table("phageLibraryNNKTheoreticalCodonBias.txt")
bg.fr <- b[,2]


########################################################
#Get the sequence of the peptide ligands for each domain
########################################################
loadFile <- function(x){
  f <- read.fasta(file=x, as.string = TRUE, seqtype = "AA")
  nm <- names(f)
  seq <- unlist(lapply(nm, function(x){f[[x]][1]}))
  return(seq)
}

seq.human <- lapply(files.human, loadFile)
seq.worm <- lapply(files.worm, loadFile)
seq.all <- c(seq.human, seq.worm)
names(seq.all) <- PDZ.all


#######################################
#Load the PDZ domain sequence alignment
#######################################
#f <- read.fasta(file="PDZ_SMART_MUSCLE_sub.fa", as.string = TRUE, seqtype = "AA")
#f <- read.fasta(file="PDZ_SMART_CLUSTAL_sub.fa", as.string = TRUE, seqtype = "AA")
f <- read.fasta(file="PDZ_phage_MUSCLE.fa", as.string = TRUE, seqtype = "AA")

PDZ.seq <- unlist(lapply(PDZ.all, function(x){f[[x]][1]})) # Extract the sequence of each PDZ domain
bs <- which(unlist(strsplit(f[[1]][1], ""))=="B")   # Find the position of the binding site

bs.seq <- unlist(lapply(PDZ.seq, function(x){ paste(unlist(strsplit(x, ""))[bs], sep="", collapse="") }))
names(bs.seq) <- PDZ.all



####################
#Create PWM matrices - part 2b
####################

#This is useful to assign a integer to each amino acid, including gaps ('X' or '-')
Naa <- 20 #Number of amino acids
gap <- 1/Naa #Contribution of gaps in the PWM
map <- c(1:(Naa+1), Naa+1)
names(map) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X", "-")

#This is the map for grouping amino acids based on their physico-chemical properties (only useful for part 6)
reduced_map <- c(4, 4, 3, 3, 4, 5, 2, 4, 2, 4, 4, 1, 4, 1, 2, 1, 1, 4, 4, 4)
names(reduced_map) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

#Choose between the two options for the two Figures to be reproduced (7 and 2 in the paper)
#Fig=2
Fig=7

compCount <- function(x){
  
  #Build a 20x10 matrix with 1 at position (i,j) if the sequence contains amino acid 'i' at position 'j'. 
  #Use the map vector to know the position (matrix row) corresponding to each amino acid. 
  #For 'X' or '-', put 'gap=1/20' in the whole column (non-specific).
  M <- matrix(data = rep(0,200),nrow = 20)
  x <- unname(strsplit(x,"")[[1]])
  for(j in 1:length(x)){
    if(unname(map[x[j]]) <= 20){M[unname(map[x[j]]),j] = 1}
    else{M[,j] = 1/20}
  }
  return(M)
}

compCount2 <- function(x){
  
  #Build a 20x10 matrix with 1 at position (i,j) if the sequence contains amino acid 'i' at position 'j'. 
  #Use the map vector to know the position (matrix row) corresponding to each amino acid. 
  #For 'X' or '-', put 'gap=1/20' in the whole column (non-specific).
  M <- matrix(data = rep(0,10*max(reduced_map)),ncol=10)
  x <- unname(strsplit(x,"")[[1]])
  for(j in 1:length(x)){
    if(x[j] %in% names(reduced_map)){M[unname(reduced_map[x[j]]),j] = 1}
    else{M[,j] = 1/20}
  }
  return(M)
}

compFreq <- function(x){
  
  x <- unlist(x)
  m <- matrix(0, nrow=20, ncol=10)
  
  #For each sequence, build a matrix with 1 at position (i,j) if the sequence contains amino acid 'i' at position 'j' 
  ml <- lapply(x, function(y){ compCount(y) })

  #Sum all the matrices from each individual sequences
  m <- Reduce('+', ml)
  m <- m/length(ml)
  if(Fig==2){
    #Only when doing the clustering (part 6), add the codon bias
    m <- t(sapply(1:length(m[,1]),function(r){m[r,]/bg.fr[r]}))
    m <- sapply(1:length(m[1,]),function(c){m[,c]*(1/sum(m[,c]))
      })
  
    #Only when doing the clustering (part 6), group the lines corresponding to similar amino acids (based on reduced_map)
    redmat <- matrix(nrow=5,ncol=10)
    m <- t(sapply(1:5,function(n){
      colSums(matrix(m[which(reduced_map==n),],ncol=10))
    }))
  }
  
  #Normalize each column of m
  
  return(m)
}
Fig=7
pwm.all <- lapply(seq.all, compFreq)

####################################
#Compute the similarity between PWMs - part 2c
####################################

compPWMSim <- function(x,y){
  #Compute the similarity between two pwm matrices as defined in the article
  w <- ncol(x)
  D <- (1/(w)) * sum(sapply(1:w,function(p){
    (1/2) * sum(sapply(1:nrow(x),function(r){
      (x[r,p]-y[r,p])^2
      }))
    }))
  S <- 1-D
  return(S)
}

#Build a le x le matrix with similarity values between every pair of PWMs


sim.pwm <- matrix(nrow=le, ncol=le)
for(i in 1:(le-1)){
  sim.pwm[i,i] <- 1
  for(j in (i+1):le){
    sim.pwm[i,j] <- compPWMSim(pwm.all[[i]], pwm.all[[j]]) 
    sim.pwm[j,i] <- sim.pwm[i,j]
  }
}
sim.pwm[le,le] <- 1
rownames(sim.pwm) <- PDZ.all
colnames(sim.pwm) <- PDZ.all


############################################
#Compute the distances in the sequence space - 3
############################################

compSeqSim <- function(x,y){
  s.x <- unlist(strsplit(unname(x),""))
  s.y <- unlist(strsplit(unname(y),""))
  #Compute the precentage of identity between two PDZ domain binding sites (string x and y)
  S <- sum(s.x==s.y)/length(s.x)
  return(S)
}

#Build a le x le matrix  with similarity values between every pair of PDZ binding site (bs.seq)
sim.seq <- matrix(nrow=le, ncol=le)
for(i in 1:(le-1)){
  sim.seq[i,i] <- 1
  for(j in (i+1):le){
    sim.seq[i,j] <- compSeqSim(bs.seq[[i]], bs.seq[[j]]) 
    sim.seq[j,i] <- sim.seq[i,j]
  }
}
sim.pwm[le,le] <- 1
rownames(sim.seq) <- PDZ.all
colnames(sim.seq) <- PDZ.all


#################################
#Plot the two kinds of similarity as in Fig 7 - 4
#################################

#Build a le x le matrix  with 1 if the two PDZ domains are in the same class, 0 if not, and -1 in the diagonal, based on the class16.
sim.class <- matrix(nrow=le, ncol=le)
for(i in 1:le){
  sim.class[i,i] <- -1
  for(j in (i+1):le){
    if(unname(class16[PDZ.all[[i]]])==unname(class16[PDZ.all[[j]]])){sim.class[i,j] <- 1} 
    else{sim.class[i,j] <- 0}  
    sim.class[j,i] <- sim.class[i,j]
  }
}
rownames(sim.class) <- PDZ.all
colnames(sim.class) <- PDZ.all

#Plot the comparison of similarity values, with different colors whether the pairs of PDZ domains are in the same class (ind1) or not (ind2)
ind1 <- which(sim.class==1)
ind2 <- which(sim.class==0)
if(Fig==7){
  pdf("output/Fig7_compSim.pdf")
  plot(sim.seq[ind2], sim.pwm[ind2], col="blue", xlim=c(0,1), ylim=c(0.4,1), pch=0, cex=0.5) #This is similar to Figure 2
  points(sim.seq[ind1], sim.pwm[ind1], col="red") 
  dev.off()
}

##################
#Do the clustering - 6
##################                     

#Go back to the compFreq function to build the PWM

#Include the codon bias (basically divide the counts observed for each amino acid by the bias in bg.fr)
#.... (in compFreq)

#Include the merging of similar amino acids in compFreq, based on the 'reduced_map' vector (basically merge the rows corresponding to the same values in reduced_map). 
#This will create PWM of size 5x10, instead of 20x10
#... (in compFreq)

#Define a distance object as 1-sim.pwm and the as.dist function
if(Fig==2){
  dis <- 1-sim.pwm
  rownames(dis) <- unlist(lapply(1:le, function(x){paste(class2[x], class7[x], class16[x], PDZ.all[x], sep="  ")})) #This will be useful to show the classes on the dendrogram
  d <- as.dist(dis)

  #Do the hierarchical clustering
  hc <- hclust(d, method = "average")

  #Plot the data with plot(as.dendrogram(hc), horiz=T, axes=F). Use par(cex=0.5, mar=c(5, 2, 2, 8)) to define the margin and text size
  par(cex=0.5, mar=c(5, 2, 2, 8))
  plot(as.dendrogram(hc), horiz=T, axes=F)
  
}


