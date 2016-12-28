require(igraph);
require(foreach);
require(BiRewire);
require(reshape);

# MEMO DATA SETTING

memoData <- function(M, events=NULL, Q=100, N=10000, permuteFun=permute.preserve.both, ...) {
	# Force convert the matrix into binary
	M <- as.matrix(M == TRUE);
	
	if( getDoParRegistered() ) {
		randomMs <- foreach(i=1:N, ...) %dopar% permuteFun(M, Q, events);
	} else {
		randomMs <- foreach(i=1:N) %do% permuteFun(M, Q, events);
	}

	if(is.null(events)) {
		eids <- rownames(M);
		events <- memoEvents(eids);
	}
	
	M <- M*1
	randomMs = lapply(randomMs, function(x) x*1)
	res <- list(M=M, R=randomMs, events=events);
	class(res) <- "memoData";
	return(res);
}

memoEvents <- function(eventIds, eventGraphIds=NULL, eventDescriptions=NULL, eventTypes=NULL) {
	if( is.null(eventGraphIds) ) {
		eventGraphIds <- as.list(eventIds);
		names(eventGraphIds) <- eventIds;
	}

	if( is.null(eventDescriptions) ) {
		eventDescriptions <- as.list(eventIds);
		names(eventDescriptions) <- eventIds;
	}

	if( is.null(eventTypes) ) {
		eventTypes <- as.list(rep("ALTERATION", length(eventIds)));
		names(eventTypes) <- eventIds;
	}

	result <- list(eventIds=eventIds, eventGraphIds=eventGraphIds, eventDescriptions=eventDescriptions, eventTypes=eventTypes);
	class(result) <- "memoEvents";
	return(result);
}

# MEMO TEST PROCEDURE

testModules <- function(memo.data, moduleIndices, FUN=memo.statistics, p.threshold=.05, samples=NULL, ...) {
	
	M = memo.data$M
	sids <- colnames(M);
	if(is.null(samples)) { samples <- sids; }

	results <- list();
	if(is.null(names(moduleIndices))) {
		names(moduleIndices) <- 1:length(moduleIndices)
	}
	
	for(aModule_name in names(moduleIndices)) {
		
		aModule = moduleIndices[[aModule_name]]
        
        if(all(aModule %in% rownames(M))){
		
            observedValue <- FUN(M[aModule, samples]);
		
            # Now pool values from the random matrices
            if( getDoParRegistered() ) {
                permutedValues <- foreach(singleR=memo.data$R, .combine="cbind", ...) %dopar% { testR <- singleR[aModule, samples]; FUN(testR) };
            } else {
                permutedValues <- sapply(memo.data$R, FUN=function(singleR) { testR <- singleR[aModule, samples];FUN(testR) });
            }
		
            # Estimate the p-val, +1s to correct for zero p-vals.
            pVal <- (sum(permutedValues >= observedValue) + 1) / (length(permutedValues) + 1);
            coverage <- length( which( colSums(memo.data$M[aModule, samples]) > 0 ) );

            aResult <- list( moduleName=aModule_name, module=aModule, pvalue=pVal, observedValue=observedValue, permutedValues=permutedValues, coverage=coverage);
                class(aResult) <- "modResult";
            results[[length(results)+1]] <- aResult;
        
        }
    }
	
	return(results);
	
}

memo.statistics <- function(M) {
	cnt <- sum(colSums(M) > 0);
	return(cnt);
}

summary <- function(modResults, events=NULL, adjust.method="fdr") {
	
	pvals <- c(NULL);
	setName <- c(NULL);
	modNames <- c(NULL);
	modIds <- c(NULL);
	sizes <- c(NULL);
	coverages <- c(NULL);
	expected <- c(NULL);
	
	for(result in modResults) {
		# result <- modResults[[resultn]]
		modId <- paste(result$module, collapse=",");
		modIds <- c(modIds, modId);
		setName <- c(setName, result$moduleName)
		if( is.null(events) ) {
			modNames <- c(modNames, modId);
		} else {
			modNames <- c(modNames, paste(events$eventDescriptions[result$module], collapse=","));
		}
		pvals <- c(pvals, result$pvalue);
		sizes <- c(sizes, length(result$module));
		coverages <- c(coverages, result$coverage);
		expected <- c(expected, mean(result$permutedValues))
	}
	
	results <- data.frame(Set=setName, 
			Module=modNames, moduleSize=sizes, observedAltered=coverages, expectedAltered=expected,
			pVal=pvals, pAdjVal=p.adjust(pvals, method=adjust.method),
			stringsAsFactors=FALSE);
	
	if(nrow(results) > 0) {
		results <- results[with(results, order(pVal)), ];
	}	
	
	return(results);
}

# PERMUTATION STRATEGIES

permute.preserve.none <- function(M, Q=100, events=NULL) {
#	M <- as.matrix(M == TRUE);
	M <- as.matrix(M);
	R = sample(M)
	dim(R) = dim(M)
	rownames(R) = rownames(M)
	colnames(R) = colnames(M)
	return(R)
}

permute.preserve.gene <- function(M, Q=100, events=NULL) {
	M <- as.matrix(M == TRUE);
	R = apply(M,1,function(x) sample(x))
	R = t(R)
	rownames(R) = rownames(M)
	colnames(R) = colnames(M)
	return(R)
}
permute.preserve.sample <- function(M, Q=100, events=NULL) {
  M <- as.matrix(M == TRUE);
  R = apply(M,2,function(x) sample(x))
  R = R
  rownames(R) = rownames(M)
  colnames(R) = colnames(M)
  return(R)
}


permute.preserve.all <- function(M, Q=100, events=NULL) {
	numOfEdges <- sum(M);
	R <- birewire.rewire.bipartite(M, max.iter=numOfEdges*Q, verbose=FALSE);
	return(R);
}