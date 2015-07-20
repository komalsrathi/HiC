#Mouse Hi-C
#call interactions within a sample
#
#Copyright: Robert Sugar (robert.sugar@ebi.ac.uk), EMBL-EBI, 2012

library(plyr)

# midPointDistance <- function(loci1, loci2, absolute=TRUE)
# {
# 	absornot <- if(absolute) abs else identity
# 	as.vector(ifelse(seqnames(loci1) == seqnames(loci2), absornot(start(loci1) + end(loci1) - (start(loci2) + end(loci2))) / 2, Inf))
# }

midPointDistanceDataFrame <- function(interactions, absolute=TRUE)
{
	absornot <- if(absolute) abs else identity
	as.vector(ifelse(interactions$chr1 == interactions$chr2, absornot(interactions$locus1 - interactions$locus2), Inf))
}

midPointDistance <- function(x, y=NULL, absolute=TRUE)
{
	absornot <- if(absolute) abs else identity
	if(is(x, "GenomicInteraction"))
	{
		loci1 <- x@source
		loci2 <- x@target
	} else if(is(x, "data.frame"))
	{
		return(as.vector(ifelse(x$chr1 == x$chr2, absornot(x$locus1 - x$locus2), Inf)))
	} else
	{
		loci1 <- x
		loci2 <- y
	}

	as.vector(ifelse(seqnames(loci1) == seqnames(loci2), absornot(start(loci1) + end(loci1) - (start(loci2) + end(loci2))) / 2, Inf))
}

distanceCorrectionEstimate <- function(interactions, binby=c("quantiles", "breaks"), quantiles=100, breaks=NULL, dosample=FALSE, sampleSize=1e8, sites=calculateCoverage(interactions))
{

	binby <- match.arg(binby)

	#estimate distance curve
	if(is(interactions, "GenomicInteraction"))
	{
		interactionDistance <- midPointDistance(interactions@source, interactions@target)	
	} else if(is(interactions, "data.frame"))
	{
		interactionDistance <- midPointDistanceDataFrame(interactions)	
	} else stop("interactions is of wrong type")
	#interactions <- interactions[order(interactionDistance)]
	#interactions$cumulative <- cumsum(interactions$frequencies)
	
	#sam<-rep(interactionDistance, interactions$contact)
	#breaks <- quantile(sam, 1:quantiles/quantiles)
	#breaks<-unique(breaks)
	
#	N <- sum(interactions$contact)
#	weights <- (seq(0, quantiles) / quantiles) * N #breaking the 0-N interval to quantiles pieces
#	breaks <- interactionDistance[findInterval(weights, interactions$cumulative), all.inside=TRUE] #TODO: rightmost? make sure it works at the two edges
#	interactionCategories <- findInterval(interactionDistance, breaks)
	
	breaks <- quantile(interactionDistance, 1:quantiles/quantiles)
	breaks <- unique(breaks)
	interactionCategories <- findInterval(interactionDistance, breaks)	
	sitesAtDistances <- distancesOfBins(sites, breaks, dosample, sampleSize)

	# TEST
	# library(plyr)
	# f <- function(x) sum(x@frequency)
	# result <- ddply(interactions, .(categories), "f")

	library(data.table)
	distancesByCategories <- data.table(reads=ifelse(is(interactions, "GenomicInteraction"), interactions@frequency, interactions$frequencies),
		categories=interactionCategories)
	setkey(distancesByCategories, "categories")
	result <- as.data.frame(distancesByCategories[, sum(reads), by="categories"])
	names(result) <- c("categories", "reads")
	merged <- merge(result, sitesAtDistances, by="categories")
	merged$expected <- merged$reads / merged$binpairs
	merged$distanceModifier <- (merged$reads / sum(merged$reads)) / (merged$binpairs / sum(merged$binpairs)) #how much more contacts do we expect in this distance compared to the average

	interactions$expectedInteractionsByDistance <- merged$expected[interactionCategories + 1]
	interactions$distanceModifier <- merged$distanceMod[interactionCategories + 1]

	return(interactions)
}

#returns the distribution pairwise of bin-bin distances 
distancesOfBins <- function(bins, breaks, dosample=FALSE, sampleSize=1e7)
{
	library(data.table)
	if(is.data.frame(bins)) bins <- GRanges(seqnames=bins$chr, ranges=bins$locus)
	allPossible <- (length(bins) ^ 2 -  length(bins)) %/% 2

	if(dosample) #TODO sample cis only?
	{
		binpairs <- GenomicInteraction(sample(bins, sampleSize, replace=TRUE), sample(bins, sampleSize, replace=TRUE))
		distances <- midPointDistance(binpairs@source, binpairs@target)
		distances <- distances[distances > 0] #no self-interaction
		distanceTable <- table(findInterval(distances, breaks))
		distanceTable <- (distanceTable / sum(distanceTable)) * allPossible #scale up
		summed <- as.data.frame(distanceTable, stringsAsFactors=FALSE)
		names(summed) <- c("categories", "binpairs")
		return(summed)
	} else
	{
		#calculate binpairs of a given distance
		binsByChromosome <- split(bins, seqnames(bins))
		distancesPerChromosome = lapply(binsByChromosome, function(x)
				{
					stopifnot(length(x)^2 < 2^31) #check int overflow

					midpoints = (end(x) + start(x)) / 2
					distanceMatrix = outer(midpoints, midpoints, function(a,b) findInterval(abs(a-b), breaks))
					#distanceMatrix[lower.tri(distanceMatrix)] <- 0 #make sure we do not count interactions twice (a-b and b-a)
					distanceTable = table(distanceMatrix[upper.tri(distanceMatrix, diag=FALSE)])
				})

		#trans
		cis <- sum(sapply(distancesPerChromosome, sum))
		trans <- allPossible - cis
		distancesPerChromosome$trans <- data.frame(Var1=length(breaks), Freq=trans) #trans in the last bin

		#sum it all up
		dpc <- data.table(do.call("rbind", lapply(distancesPerChromosome, 
			function(tbl) as.data.frame(tbl, stringsAsFactors=FALSE))))
		setkey(dpc, "Var1")
		summed <- as.data.frame(dpc[,sum(Freq), by="Var1"])
		names(summed) <- c("categories", "binpairs")
		return(summed)
	}
}

#for(i in 2^(1:15)){ print(paste(i, system.time((distancesOfBins(baitGR[1:min(i, length(baitGR))], 1e7)))))}

#estimates the density function given a set of distances
estimateDensityFromDistances <- function(distances, applyLog = TRUE)
{
	distances <- distances[!is.infinite(distances)]
	if(!applyLog)
	{
		densfunc <- density(distances)
		#plot(densfunc, main = "estimate of distances", xlab = "log(distance)", ylab = "density function")

		#ifun <- approxfun(2^(densfunc$x), densfunc$y / binwidths)
		ifun <- splinefun(densfunc$x, densfunc$y)
	} else
	{
		distances <- log2(distances + 1)

		densfunc <- density(distances)
		#plot(densfunc, main = "estimate of distances", xlab = "log(distance)", ylab = "density function")
		
		binwidths <- (diff(exp(densfunc$x)))
		binwidths <- c(binwidths, binwidths[length(binwidths)])
		
		#ifun <- approxfun(2^(densfunc$x), densfunc$y / binwidths)
		ifun <- splinefun(2^(densfunc$x), densfunc$y / binwidths)
	}

	#plot(2^(seq(0,32,0.01)), ifun(2^(seq(0,32,0.01))),main = "estimate of distances", xlab = "log(distance)", ylab = "density function", log = 'xy')
	
	return(ifun)

}

numberOfBinpairsWithDistanceBetweenTwoLargerRegions <- function(x1, x2, res, d)
{
#	if we estimated distance with a given resolution, but have regions of different size - how to calculate expected contacts?
#	for example if bins are adjacent megabase bins, we underestimate contacts assuming that the average contact frequency is 
#	the one 1Mb away - as some regions are much closer than that
	
#	x1: region1Length
#	x2: region2Length, resolution, distanceBetweenRegions
#	res: granularity with which we estimated contacts
#	d: distance of the two regions (end of x1 and start of x2)
#returns: a vector containing the number of binpairs with a given distance
	if(x1 %/% res * res != x1 || x2 %/% res * res != x2) stop("lengths should be a multiple of the resolution")	
	
	smaller <- min(x1, x2)
	larger <- max(x1, x2)
	
	binPairsPerDistance <- list()
	i = 0
	for(dist in seq(d + res, smaller + d, by = res))
	{
		i = i + 1
		binPairsPerDistance[[as.character(dist)]] <- i
	}
	for(dist in seq(smaller + d, larger + d, by = res))
	{
		binPairsPerDistance[[as.character(dist)]] <- i
	}
	for(dist in seq(larger + d + res, larger + smaller + d - res, by = res))
	{
		i = i - 1
		binPairsPerDistance[[as.character(dist)]] <- i
	}
	
	binPairsPerDistance
}

numberOfBinpairsPerDistance <- function(chromosomeLengths = NULL, bins = NULL, resolution)
{	
	if(!is.null(chromosomeLengths))
	{
		bins <- getBins(chromosomeLengths, binsize = resolution)
		possibleDistances <- seq(0, max(chromosomeLengths), resolution)

		#cis interactions
		#if all bins are of equal size, for a given distance, every bin further away than "distance" from the end of the chromosome will have a partner in that distance 
		allDistanceOccurences <- lapply(possibleDistances, function(distance)
				{
					occurencesPerChromosome <- lapply(chromosomeLengths, function(chromosome)
							{
								max(chromosome - distance, 0) %/% resolution
							})
					
					sum(unlist(occurencesPerChromosome))
				})
		
		names(allDistanceOccurences) <- possibleDistances

		#trans interactions
		numberOfBinsByChromosome <- chromosomeLengths %/% resolution + 1
		numberOfBins = sum(numberOfBinsByChromosome)
		transInteractions <- lapply(numberOfBinsByChromosome, function(x)
				{
					x * (numberOfBins - x) #Inf is for trans hits
				})
		
		allDistanceOccurences[['Inf']] <- sum(unlist(transInteractions))
		
		return(allDistanceOccurences)		
	}
	
	if(!is.null(bins))
	{
		#calculate binpairs of a given distance
		binsByChromosome <- split(bins, seqnames(bins))
		distancesPerChromosome = lapply(binsByChromosome, function(x)
				{
					midpoints = (end(x) + start(x)) / 2
					distanceMatrix = outer(midpoints, midpoints, function(x,y) abs(x-y))
					#distanceMatrix[lower.tri(distanceMatrix)] <- 0 #make sure we do not count interactions twice (a-b and b-a)
					distanceTable = table(distanceMatrix)
				})
		allDistanceOccurences = list()
		notUsed <- lapply(distancesPerChromosome, function(chromosome) vapply(names(chromosome), function(x)
							{
								allDistanceOccurences[[x]] <<- max(allDistanceOccurences[[x]] + chromosome[[x]], chromosome[[x]]) #add to existing or create
								return(0)
							}, 0))
		#trans bins
		allDistanceOccurences[['Inf']] <- 0
		numberOfBinsByChromosome <- lapply(binsByChromosome, length)
		numberOfBins = sum(unlist(numberOfBinsByChromosome))
		notUsed <- lapply(numberOfBinsByChromosome, function(x)
				{
					allDistanceOccurences[['Inf']] <<- allDistanceOccurences[['Inf']] + x * (numberOfBins - x) #Inf is for trans hits
					return(0)
				}
		)
		allDistanceOccurences = lapply(allDistanceOccurences, function(x) x / 2) #we counted each interaction twice
		allDistanceOccurences[['0']] = allDistanceOccurences[['0']] * 2 #self-interactions only appear once
		notUsed <- lapply(allDistanceOccurences, function (x)
				{
					if(x != as.integer(x)) 
						stop('assert: trans interactions should be an integer value')
				})
		return(allDistanceOccurences)
	}
}

estimateContactsByDistance <- function(interactions, applyBinning = FALSE, resolution = 1e7, chromosomeLengths = NULL)
{
	#chromosomeLengths: whole genome will be used for the estimate
	
#	#bin it
#	lastBin = max(distancesNotInf)	
#	binStarts = seq(0, lastBin, binSize)
	
	#if(resolution > 0) interactions$distance = (interactions$distance %/% resolution) * resolution
	
	if(applyBinning) interactions = binInteractions(interactions, resolution)
	
	interactions$distance <- ifelse(interactions$chr1 == interactions$chr2, (abs(interactions$locus2 - interactions$locus1)), Inf)		
	
	#sums by distance
	sumsByDistance <- vapply(split(interactions$frequencies, interactions$distance), sum, 0)
	
	if(is.null(chromosomeLengths))
	{
		bins = getBinsFromInteractions(interactions, asGranges = TRUE)
		allDistanceOccurences <- numberOfBinpairsPerDistance(chromosomeLengths = NULL, bins, resolution)
	} else
	{
		allDistanceOccurences <- numberOfBinpairsPerDistance(chromosomeLengths, bins = NULL, resolution)		
	}
	
	observedDf = data.frame(distance = names(sumsByDistance), observed = unlist(sumsByDistance), stringsAsFactors = FALSE)	
	expectedDf = data.frame(distance = names(allDistanceOccurences), expected = unlist(allDistanceOccurences), stringsAsFactors = FALSE)	
	
	distanceDf = merge(observedDf, expectedDf, by = 'distance', all = TRUE, stringsAsFactors = FALSE)
	
	#replace NAs with 0
	distanceDf$expected[is.na(distanceDf$expected)] <- 0
	distanceDf$observed[is.na(distanceDf$observed)] <- 0
	
	distanceDf$expected <- distanceDf$observed / distanceDf$expected
	
	#return(distanceDf[, c('distance', 'expected')])
	
	expectedInteractionsByDistance <- distanceDf$expected
	names(expectedInteractionsByDistance) <- distanceDf$distance
	estimateContactsByDistanceResolution <- resolution
	infIndex <- which(distanceDf$distance == Inf)
	linApprox <- approxfun(distanceDf$distance[-infIndex], distanceDf$expected[-infIndex])
	#loessApprox <- loess(as.numeric(distanceDf$distance[-infIndex]) ~ distanceDf$expected[-infIndex])
	
	estimateContactsByDistance = function(distance, binsize1 = estimateContactsByDistanceResolution, binsize2 = binsize1,approximation = c("linear", "loglinear", "spline", "logspline"))
	{	
		if(approximation[1] == "linear")
		{
			baseEstimate <- ifelse(distance == Inf, expectedInteractionsByDistance['Inf'], linApprox(distance)) 
#		} else if(approximation[1] == "loess")
#		{
#			baseEstimate <- ifelse(distance == Inf, expectedInteractionsByDistance['Inf'], loessApprox(distance))
		} else stop("approximation method not supported")
		
		sizeAdjustmentFactor <- (binsize1 / resolution) * (binsize2 / resolution)
		return(baseEstimate * sizeAdjustmentFactor)		
	}
	
	return(estimateContactsByDistance)
	
	#normalise expected values
#	distanceDf$expected = (distanceDf$expected /  sum(distanceDf$expected)) * sum(distanceDf$observed)
#calculate pvalues based on poisson distribution
}

callInteractions  <- function(gi, contactsEstimate, minDistance = 0, callIndividualInteractions = TRUE)
{
	#interactions: in a data frame
	#targetLoci: interactions to call
	#bins: bins where the interactions should be called

	if(!is(gi, "GenomicInteraction")) stop("gi must be a GenomicInteraction")
	
	distance <- midPointDistance(gi@source, gi@target)	
	expected <- ifelse(distance <= minDistance, NA, contactsEstimate(distance, width(gi@source), width(gi@target)))
	
	iframe <- data.frame(frequency = gi@frequency, expected = expected)
	gi@elementMetadata$distance = distance
	gi@elementMetadata$expected = expected
	
	if(callIndividualInteractions)
	{
		gi@elementMetadata$pvalues = apply(iframe, 1, function(x) if(is.na(x['expected'])) NA else poisson.test(x['frequency'], x['expected'])$p.value)
		gi@elementMetadata$adjustedPvalues <- p.adjust(gi@elementMetadata$pvalues, method = "BH")
	}
	#interactions$padj <- adjustedPvalues
	return(gi)
}

callInteractionsByBins <- function(gi, contactEstimate, targets, bins, excludeDiagonal = FALSE)
{
	#diag: include intra-bin interactions
	
	#by standard binning
	#calculate number of targets in bins
	targetsPerBins <- countOverlaps(bins, targets)

	matches <- findOverlaps(gi@source, bins)@matchMatrix
	gi@elementMetadata$sourceTargetCount <- 0
	#make sure source/targe falls only in one bin
	if(length(matches[,1]) != unique(length(matches[,1]))) stop("source range falls into multiple bins")	
	gi@elementMetadata$sourceTargetCount[matches[,1]] <- targetsPerBins[matches[,2]]
	
	matches <- findOverlaps(gi@target, bins)@matchMatrix
	gi@elementMetadata$targetTargetCount <- 0
	#make sure target/targe falls only in one bin
	if(length(matches[,1]) != unique(length(matches[,1]))) stop("target range falls into multiple bins")	
	gi@elementMetadata$targetTargetCount[matches[,1]] <- targetsPerBins[matches[,2]]
	
	#pick the interactions between targets
	interactionsBetweenTargets <- which((gi@elementMetadata$sourceTargetCount > 0) & (gi@elementMetadata$targetTargetCount > 0))

	# expected:
	binsWithTargets <- which(targetsPerBins > 0)
	#if a bin contains n targets and the other m, there are n*m possible targetInteractions
	numberOfTargetInteractions <- outer(binsWithTargets, binsWithTargets, function(x,y) targetsPerBins[x] * targetsPerBins[y])
	numberOfTargetInteractions[lower.tri(numberOfTargetInteractions, diag = excludeDiagonal)] <- 0
	distances <- outer(binsWithTargets, binsWithTargets, function(x,y) midPointDistance(bins[x], bins[y]))
	expectedContacts <- contactsEstimate(distances)
	expectedInteractions <- expectedContacts * numberOfTargetInteractions
	expectedSum <- sum(expectedInteractions)
		
	observed <- gi@frequency * gi@elementMetadata$targetTargetCount * gi@elementMetadata$sourceTargetCount
	if(excludeDiagonal) observed <- subset(observed, midPointDistance(gi@source, gi@target) != 0)
	observedSum <- sum(observed)

	#perform Poisson test
	#the event count in the Poisson test should not be influenced by the number of targets - they influence only the weights 
	eventCountMultiplier <- observedSum / sum(gi@frequency[gi@elementMetadata$targetTargetCount * gi@elementMetadata$sourceTargetCount > 0])
	significance <- poisson.test(round(observedSum / eventCountMultiplier), expectedSum / eventCountMultiplier)
	
	#high-resolution interactions of targets are connected more strongly than lower-res or neighboring same-res regions 
	#(e.g. the interactions are focused on the actual targets)
}

#callInteractions <- function(interactions, estimateContactsByDistance)
#{
#	#interactions: in a data frame
#	#targetLoci: interactions to call
#	#bins: bins where the interactions should be called
#		
#	#calculate distances
#	#returns a function that calculates estimated contacts based on distance
#		
#	
#
#	#select bins
#	
#	
#	#interactions$expected  <- estimateContactsByDistance(interactions$distance)
#	distances <- 
#	expected <- 
#	
#	#call all
#	interactions$pvalues = apply(interactions, 1, function(x) poisson.test(as.numeric(x['frequencies']), as.numeric(x['expected']))$p.value)
#	adjustedPvalues <- p.adjust(interactions$pvalues, method = "BH")
#	interactions$padj <- adjustedPvalues
#	
#	#interactions <- merge(interactions, expectedContactsByDistance, by = 'distance', all.x = TRUE)
#
#	return(interactions)
#}

estimateRandomInteractions <- function(coverage, totalNumberOfInteractions, clearDiagonal = TRUE, applysqrt = FALSE)
{	
	#create matrix
	size = length(coverage)
	estimate <- matrix(0, nrow = size, ncol = size) #matrix full of 0s
	colnames(estimate) = rownames(estimate) = paste(seqnames(coverage), start(coverage), sep = "_")
	
	coverageVector = coverage@elementMetadata$coverage
	
	for(i in seq(along = coverageVector))
	{
		for(j in seq(along = coverageVector))
		{
			if(!clearDiagonal || i!= j)
			{
				estimate[i, j] = ifelse(applysqrt, sqrt(coverageVector[i] * coverageVector[j]), coverageVector[i] * coverageVector[j])
			}
		}
	}
	
	#normalise by totalNumberOfInteractions
	totalInEstimate = sum(estimate)
	correctionFactor = totalNumberOfInteractions / totalInEstimate
	#if(clearDiagonal) correctionFactor = correctionFactor * ((size * size - size) / (size * size)) #take into account that the diagonal is missing
	
	estimate = estimate * correctionFactor
	
	return(estimate)
}

#------- Binomial Estimate -------

#filters interactions by regions
#
#bothInRegions: where both ends of the interactions are in the regions. Result will be a smaller matrix
#anyInRegions: where any (one or both) of the the ends are in the region. Will result in a "swiss cheese" matrix with holes
#oneInRegions: where exactly one end is in the regions. 
#
#betweenRegionsAndRegions2: where one end is in regions and the other is in regions2
#
#keep: keep the selected interactions
#remove: remove the selected interactions
#both: return both "keep" and "remove" in a list

subsetInteractionsByRegion <- function(interactions, regions, type=c("bothInRegions", "anyInRegions", "oneInRegions", "betweenRegionsAndRegions2"), action=c("keep", "remove", "both"), regions2=NULL)
{
	if(is.null(regions) || length(regions) == 0) stop("regions should be filled")

	if(class(interactions)=="data.frame")
	{
		interactions <- as(interactions, "GenomicInteraction")
	}
	
	type = match.arg(type)
	action = match.arg(action)
	
	#figure out ends
	firstEndInRegions <- !is.na(findOverlaps(interactions@source, regions, select="first"))
	secondEndInRegions <- !is.na(findOverlaps(interactions@target, regions, select="first"))
	
	if(type=="bothInRegions")
	{
		chosenOnes <- firstEndInRegions & secondEndInRegions
	} else
	if(type=="anyInRegions")
	{
		chosenOnes <- firstEndInRegions | secondEndInRegions
	} else
	if(type=="oneInRegions")
	{
		chosenOnes <- xor(firstEndInRegions, secondEndInRegions)
	} else
	if(type=="betweenRegionsAndRegions2")
	{
		if(is.null(regions2) || length(regions2) == 0) stop("regions2 should be filled when using betweenRegionsAndRegions2")
		firstEndInRegions2 <- !is.na(findOverlaps(interactions@source, regions2, select="first"))
		secondEndInRegions2 <- !is.na(findOverlaps(interactions@target, regions2, select="first"))
		
		chosenOnes <- (firstEndInRegions & secondEndInRegions2) | (firstEndInRegions2 & secondEndInRegions)
	}
	
	if(action=="remove")
	{
		return(giRange(interactions, !chosenOnes))
	} else if(action=="keep")
	{
		return(giRange(interactions, chosenOnes))		
	} else #both
	{
		kept <- giRange(interactions, chosenOnes)
		removed <- giRange(interactions, !chosenOnes)
		return(list(kept=kept, removed=removed))
	}
}

#calculates the expected probability of regions based on coverages
#
#bothInRegions: where both ends of the interactions are in the regions. Result will be a smaller matrix
#anyInRegions: where any (one or both) of the the ends are in the region. Will result in a "swiss cheese" matrix with holes
#oneInRegions: where exactly one end is in the regions. 
#
#betweenRegionsAndRegions2: where one end is in regions and the other is in regions2
#
#keep: keep the selected interactions
#remove: remove the selected interactions
calculateExpectedWeightByRegion <- function(coverages, regions, type=c("bothInRegions", "anyInRegions", "oneInRegions", "betweenRegionsAndRegions2"), action=c("remove", "keep"), regions2=NULL)
{
	
}	

binomialTest <- function(interactions, considerFullMatrix = TRUE, noiseEstimate = 1, noDiagonal = TRUE, distanceCorrection=FALSE)
{
	if(noDiagonal)
		interactions <- subset(interactions, !((locus1 == locus2) & (as.character(chr1) == as.character(chr2))))
	
	#noiseEstimate: value between 0 and 1
	coverages <- calculateCoverage(interactions, relativeCoverage = FALSE)
	coverages$coverage <- coverages$coverage / sum(coverages$coverage)
		
	#add source and target coverage to the interactions
	names(coverages) <- c('chr1', 'locus1', 'coverage_source')
	interactions = merge(interactions, coverages, by = c('chr1', 'locus1'), all.x = TRUE)
	names(coverages) <- c('chr2', 'locus2', 'coverage_target')
	interactions = merge(interactions, coverages, by = c('chr2', 'locus2'), all.x = TRUE)
	names(coverages) <- c('chr', 'locus', 'coverage')

	N = sum(interactions$frequencies)
		
	#predict using coverages
	interactions$probabiltyOfInteraction <- (interactions$coverage_source) * (interactions$coverage_target)  
	if(distanceCorrection)
	{
		if(nrow(coverages) ^ 2 > 2^31)
		{
			dosample = TRUE
			sampleSize = 1e8
		} else
		{
			dosample = FALSE
			sampleSize = Inf
		}

		interactions <- distanceCorrectionEstimate(interactions, binby="quantiles", quantiles=100, breaks=NULL, dosample=dosample, sampleSize=sampleSize, sites=coverages)
		interactions$probabiltyOfInteraction <- interactions$probabiltyOfInteraction * interactions$distanceModifier
	}	
	
	#we multiply non-self interactions by two because 
	#we consider 1->2 and 2->1 interaction the same (diagonal is expected to be taken out)
	if(!noDiagonal)
	{
		attach(interactions)
		nonSelfInteractions <- !((locus1 == locus2) & (as.character(chr1) == as.character(chr2)))
		interactions$probabiltyOfInteraction[nonSelfInteractions] <- interactions$probabiltyOfInteraction[selfInteractions] * 2
		detach(interactions)
	}		
	
	#correct for the fact that the diagonal has been removed -> interactions probabilities will be adjusted
	if(noDiagonal)
	{
		sumOfDiagonal <- sum(coverages$coverage ^ 2)
		interactions$probabiltyOfInteraction <- interactions$probabiltyOfInteraction / (1 - sumOfDiagonal)
	}
		
	#if we have a good estimate for noise levels, we could increase power by 
	if(noiseEstimate < 1)
		interactions$probabiltyOfInteraction <- interactions$probabiltyOfInteraction * noiseEstimate
	
	interactions$predicted <- interactions$probabiltyOfInteraction * N #knowing the probability of a single interaction, we can calculate how many we expect
	
	interactions$logFoldChange <- log2(interactions$frequencies / interactions$predicted)
	
	#binomial/poisson test
	
	interactions$pvalue <- apply(interactions, 1, function(x)
		{
			binom.test(as.numeric(x[["frequencies"]]), N, as.numeric(x[["probabiltyOfInteraction"]]), alternative = "greater")$p.value
			#interactions$pvalue[i] <- poisson.test(interactions$frequencies[i], N * probabiltyOfInteraction[i], alternative = "greater")$p.value
			#print(interactions$pvalue[i])
		}	
	)

	#multiple testing correction
	if(considerFullMatrix)
	{
		numberOfInteractions <- ((nrow(coverages) ^ 2)-nrow(coverages))/2
		interactions$qvalue <- p.adjust(interactions$pvalue, method = "BH", n = numberOfInteractions)
	} else
	{
		interactions$qvalue <- p.adjust(interactions$pvalue, method = "BH")
	}
	return(interactions)
}

callBinomialFromLogit <- function(interactions, considerFullMatrix, parallel=FALSE, cores=4)
{
	library("data.table")
	if(parallel)
	{
		print("running garbage collector before parallel fork")
		gc()
	}

	interactions$logFoldChange <- log2(interactions@frequency / interactions$predicted_logit)
	
	#binomial/poisson test
	N <- sum(interactions@frequency)
	df <- data.frame(frequency=as.numeric(interactions@frequency), probability_logit=as.numeric(interactions$probability_logit))

	if(parallel)
	{
		library(multicore)
		binomParams <- as.data.frame(t(cbind(
			as.numeric(df[["frequency"]]), 
			as.numeric(df[["probability_logit"]]
				))))
		
		interactions$pvalue <- unlist(mclapply(binomParams, function(x)
			{
			   binom.test(x[1], N, x[2], alternative = "greater")$p.value
			}, 
			mc.cores=cores))
	} else
	{
		interactions$pvalue <- apply(df, 1, function(x)
			{
				binom.test(x[["frequency"]], N, x[["probability_logit"]], alternative = "greater")$p.value
			}	
		)
	}

	#multiple testing correction
	if(considerFullMatrix)
	{
		numberOfSites <- length(unique(c(unique(interactions@source), unique(interactions@target))))
		numberOfInteractions <- ((numberOfSites ^ 2) - numberOfSites) / 2
		interactions$qvalue <- p.adjust(interactions$pvalue, method = "BH", n = numberOfInteractions)
	} else
	{
		interactions$qvalue <- p.adjust(interactions$pvalue, method = "BH")
	}

	return(interactions)
}

estimateNoiseFromCisTrans <- function(interactions, chromosomeWeights = NULL)
{
	#estimates noise levels by looking at the number of cis vs. trans interactions
	#it is a conservative estimates that assumes that all trans reads are due to noise
	#chromosomeWeights: this could be simply the length of the chromosomes - if not given it will be calculated from the coverage
	
	if(!('frequencies' %in% names(interactions)))
		interactions$frequencies <- rep(1, nrow(interactions))
	transRatio <- sum(interactions$frequencies[as.character(interactions$chr1) != as.character(interactions$chr2)]) / (sum(interactions$frequencies))

	if(is.null(chromosomeWeights))
	{
		chromosomeWeightsSource <- ddply(interactions, "chr1", function(x) sum(x[['frequencies']])); names(chromosomeWeightsSource) <- c("chr", "weightSource")
		chromosomeWeightsTarget <- ddply(interactions, "chr2", function(x) sum(x[['frequencies']])); names(chromosomeWeightsTarget) <- c("chr", "weightTarget")
		merged <- merge(chromosomeWeightsSource, chromosomeWeightsTarget, by = "chr", all = TRUE)
		merged[is.na(merged)] <- 0
		chromosomeWeights <- merged$weightSource + merged$weightTarget
		names(chromosomeWeights) <- merged$chr
	}
	
	chrWeightMatrix <- outer(chromosomeWeights, chromosomeWeights) #multiply weights pairwise
	cis = sum(diag(chrWeightMatrix))
	total = sum(chrWeightMatrix)
	trans = total - cis
	
	#assuming that all trans reads are noisy - that also means that some cis will randomly be noise
	noiseLevelEstimate <- min(transRatio * (total / trans), 1)
}

#separates interactions by distance
#Inf means a trans interaction
separateInteractionsByDistance <- function(interactions, distances=c(10000, 100000, 1000000, Inf), printSummary=TRUE, sampleName="")
{
	if(is.data.frame(interactions))
		interactions <- as(interactions, "GenomicInteraction")
	d <- midPointDistance(interactions@source, interactions@target)
	category <- findInterval(d, distances)
	
	result <- list()	
	for(i in 0:length(distances))
	{
		if(i == 0)
		{
			label <- paste("-", distances[1], sep="")
		} else
		if(i == length(distances) - 1 && is.infinite(distances[length(distances)]))
		{
			label <- paste(distances[length(distances) - 1], "-", sep="")
		} else
		if(i == length(distances))
		{
			label <- paste(distances[length(distances)], "-", sep="")
			if(is.infinite(distances[length(distances)]))
				label= "trans"
		} else
		{
			label = paste(distances[i], distances[i+1], sep='-')
		}
		
#		item <- giRange(interactions, category == i)
#		if(is.nonEmpty(item)) result[[label]] <- item 
		
		result[[label]] <- giRange(interactions, category == i)
	}
	
	if(printSummary)
	{
		print(paste("Separating ", length(interactions), " interactions in ", sampleName, sep=""))
		for(i in seq_along(result))
		{
			print(paste(names(result[i]), ": ", length(result[[i]]), sep=""))
		}
	}
	
	return(result)	
}

topSignificant <- function(interactions, filterBy=c("qvalue", "pvalue"), cutoff=0.05, maxItems=Inf, minLogFoldChange =-Inf)
{
	filterBy <- match.arg(filterBy)
	interactions <- giRange(interactions, interactions@elementMetadata[[filterBy]] <= cutoff)
	interactions <- giRange(interactions, interactions@elementMetadata$logFoldChange >= minLogFoldChange)
		
	if(length(interactions) == 0)
		return(interactions)
	ranks <- order(interactions@elementMetadata[[filterBy]])
	interactions <- giRange(interactions, ranks[seq(min(length(ranks),maxItems))])
}

saveInteractionsAsCsv <- function(interactions, filename, allFields=FALSE)
{
	ep <- as.data.frame(interactions) 
	if(allFields == FALSE) ep <- ep[c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "frequency", "logFoldChange", "pvalue", "qvalue")]
	names(ep)<-sub("seqnames", "chr", names(ep))
	write.csv(ep, quote=FALSE, row.names = FALSE, file=filename)	
}

saveInteractionsAsTxt <- function(interactions, filename, allFields=FALSE)
{
	ep <- as.data.frame(interactions) 
	if(allFields == FALSE) ep <- ep[c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "frequency", "logFoldChange", "pvalue", "qvalue")]
	names(ep)<-sub("seqnames", "chr", names(ep))
	write.table(ep, quote=FALSE, row.names = FALSE, sep='\t', file=filename)	
}

saveTopSignificantByRange <- function(interactions, sampleName, distances=c(10000, 100000, 1000000, Inf), filterBy=c("qvalue", "pvalue"), cutoff=0.05, maxItems=10000)
{
	interactions <- topSignificant(interactions, filterBy, cutoff)
	
	save(interactions, file=paste(sampleName, "_significant", sep=""))
	saveInteractionsAsCsv(interactions, paste(sampleName, "_significant", ".csv", sep=""))
	
	interactionsByRange <- separateInteractionsByDistance(interactions, distances, printSummary=TRUE, sampleName=sampleName)
	for(i in seq(along=interactionsByRange))
	{
		interactions <- interactionsByRange[[i]]
		save(interactions, file=paste(sampleName, "_", names(interactionsByRange)[i], sep=""))
		interactions <- topSignificant(interactionsByRange[[i]], filterBy, cutoff, maxItems)
		saveInteractionsAsCsv(interactions, paste(sampleName, "_", names(interactionsByRange)[i], ".csv", sep=""))
	}
}

#BED file extended with interaction data could be read by
#the Wash U Human Epigenome Browser. File format specification:
# http://washugb.blogspot.co.uk/2012/09/prepare-custom-long-range-interaction.html
#
#extraFields: FALSE (only necessary ones), TRUE (all fields), or specific fields, for example c("pvalue", "logFoldChange") 

saveAsBed <- function(interactions, filename, scoreField=c("frequency", "frequencies"), extraFields = FALSE, seqmonk=FALSE)
{
	if(!is(interactions, "GenomicInteraction"))
		interactions <- as(interactions, "GenomicInteraction")

	# if(!seqmonk)
	# {
	# 	#add mirror interactions
	# 	interactions <- GenomicInteraction(c(interactions@source, interactions@target),
	# 									   c(interactions@target, interactions@source),
	# 									   c(interactions@frequency, interactions@frequency),
	# 									   rbind(interactions@elementMetadata, interactions@elementMetadata))
	# }

	interactions <- as.data.frame(interactions)
	
	#for genomic interactions
	names(interactions) <- sub("seqnames", "chr", names(interactions))
	#for chr-locus
	if(		all(c("locus1", "locus2") %in% names(interactions)) &&
			all(!(c("start1", "start2", "end1", "end2") %in% names(interactions))))
	{
		names(interactions) <- sub("locus", "start", names(interactions))
		interactions$end1 <- interactions$start1	
		interactions$end2 <- interactions$start2	
	}
	
	#pick right scoreFiled
	scoreField <- scoreField[scoreField %in% names(interactions)][1]
	
	#add bed extension if needed
	if(is.empty(grep(".bed$", filename)))
		filename <- sub("$", ".bed", filename)
	
	if(!seqmonk) #WashU format (http://washugb.blogspot.co.uk/2012/09/prepare-custom-long-range-interaction.html)
	{
        #create a reverse interaction for all entries
        #they have to be there separately		mirrorInteractions <- interactions
        #pick fields to flip:
        mirrorInteractions <- interactions
        fieldsToFlip <- sapply(names(interactions), function(x) sub("other1", "2", sub("2", "1", sub("1", "other1", x))) %in% names(interactions))
        names(mirrorInteractions)[fieldsToFlip] <- sub("other1", "2", sub("2", "1", sub("1", "other1", names(mirrorInteractions)[fieldsToFlip]))) #flip 1 and 2 in names
        interactions <- rbind(interactions, mirrorInteractions)

		#create interaction field in format: chrX:123-456,3.14	
		interactions$interactionPartner <- paste(interactions$chr2, ":", interactions$start2, "-", interactions$end2, ",", interactions[,scoreField] ,sep="")
		interactions$ID <- row.names(interactions)
		interactions$relativeDirection <- "."
		
		basicFields <- c("chr1", "start1", "end1", "interactionPartner", "ID", "relativeDirection")
		if(extraFields == FALSE)
		{
			fields <- basicFields
		} else
		{
			if(extraFields == TRUE) extraFields <- names(interactions)
			
			fields <- c(basicFields, extraFields[!(extraFields %in% basicFields)])
		}
		
		write.table(interactions[fields], quote=FALSE, row.names = FALSE, , sep = "\t", col.names=FALSE, file=filename)		
	}	else #write interactions in two consecutive rows
	{
		if(!any(grepl("names1", names(interactions)))) interactions$names1 <- "NA"
		if(!any(grepl("names2", names(interactions)))) interactions$names2 <- "NA"

		oddRows <- interactions[, c("chr1", "start1", "end1", "names1", scoreField, "strand1")]			
		names(oddRows) <- c("chr", "start", "end", "names", scoreField, "strand")
		evenRows <- interactions[, c("chr2", "start2", "end2", "names2", scoreField, "strand2")]			
		names(evenRows) <- c("chr", "start", "end", "names", scoreField, "strand")
		#merge odd and even rows
		halfLine <- nrow(oddRows)
		merged <- rbind(oddRows, evenRows)
		toWrite <- merged[rep(seq(1,halfLine), rep(2, halfLine)) + rep(c(0,halfLine), halfLine) , ] #zig-zag: if 3 rows the seq is: 1,4,2,5,3,6
		# toWrite$thickStart <- toWrite$start
		# toWrite$thickEnd <- toWrite$end
		# toWrite$rgb <- "255,0,0"
		write.table(toWrite, quote=FALSE, row.names = FALSE, sep = "\t", file=filename, col.names=FALSE)		
	}
}

saveInteractionsAsBed <- function(gr, name)
{
	if(length(gr) > 0)
	{
		filename <- paste(name, "_", gr[1]$names, ".bed", sep="")
		library(rtracklayer)
		print(paste("saving", filename))
		export.bed(gr, filename)
		return(NULL)
	}
}


#pools interactions from different samples
#samples: filenames or data frames
#aggregate: update frequencies to sum up all samples
#outputFilename: save results in outputFilename under the name interactions
poolInteractions <- function(samples, aggregate=TRUE, deduplicate=FALSE)
{
	if(all(is.character(samples)) || all("character" == unlist(lapply(samples, class))))
	{
		sampleInteractions <- list()
		for(i in samples)
		{
			loaded <- load(i)
			interactions <- get(loaded)
			if(!is.data.frame(interactions))
				stop("poolInteractions: expecting data frame in file")
			sampleInteractions[[i]] <- interactions
		}
	} else
	if(all("data.frame" == unlist(lapply(samples, class))))
	{
		sampleInteractions <- samples
	} else
	{
		stop("poolInteractions: expecting data frames or filenames")
	}
	
	interactions <- do.call(rbind, sampleInteractions)
	
	print(paste("pooling:", paste(names(samples), collapse=',')))
	print(paste("raw:", nrow(interactions)))
	
	if(aggregate)
	{
		interactions <- countDuplicates(interactions)
		print(paste("unique:", nrow(interactions)))
	}
	
	if(deduplicate)
	{
		#set unmapped frequencies to 1
		print(paste("deduplication: removing", sum(interactions$frequencies) - nrow(interactions), "from", sum(interactions$frequencies)))		
		interactions$frequencies <- rep(1, length(interactions$frequencies))
		
	}
	
	return(interactions)
}

deduplicate <- function(interactions)
{
	#sort it
	interactions <- interactions[order(interactions$chr1, interactions$locus1, interactions$chr2, interactions$locus2), ]	
	#bin it
	interactions <- countDuplicates(interactions)

	#set unmapped frequencies to 1
	print(paste("deduplication: removing", sum(interactions$frequencies) - nrow(interactions), "from", sum(interactions$frequencies)))		
	interactions$frequencies <- rep(1, length(interactions$frequencies))

	return(interactions)
}

#------- tests -------
testEstimateContactsByDistance <- function()
{
	library(GenomicRanges)
	RESOLUTION = 1e7
	load("test/testInteractionMatrix")
	interactions = binInteractions(testInteractionMatrix, RESOLUTION)
	bins = getBinsFromInteractions(interactions, asGranges = TRUE)
	
	contactsEstimate <- estimateContactsByDistance(interactions, applyBinning = FALSE, resolution = 1e7)
	gi <- as(testInteractionMatrix, "GenomicInteraction")
	interactionsByTarget <- binInteractionsByTarget(gi, targets)
	gi <- callInteractions(gi, contactsEstimate)
	
	#by standard binning
	#calculate number of targets in bins
	#targetsPerBin <- ...
	#interactions@elementMetadata$ <- contactsEstimate()
	#interactions@elementMetadata$expected <- contactsEstimate()
	#observedSum <- sum(gi@frequency)
	#expectedSum <- sum(gi@elementMetadata$expected)
	#significance <- poisson.test(observedSum, expectedSum)
	
	#high-resolution interactions of targets are connected more strongly than lower-res or neighboring same-res regions 
	#(e.g. the interactions are focused on the actual targets)
}

testBinomial <- function()
{
	interactions <- data.frame(chr1=c("chr1", "chr1","chr2"), 
			locus1=c(1,1,2), 
			chr2=c("chr1", "chr2", "chr2"), 
			locus2=c(1,2,2),
			frequencies=c(1,2,1))
	binomialTest(interactions)
	
	interactions <- data.frame(chr1=c("chr1"), 
			locus1=c(1), 
			chr2=c("chr2"), 
			locus2=c(2),
			frequencies=c(2))
	binomialTest(interactions)

	interactions <- data.frame(chr1=c("chr1", "chr1","chr2"), 
			locus1=c(1,1,2), 
			chr2=c("chr1", "chr2", "chr2"), 
			locus2=c(1,2,2),
			frequencies=c(1,5,4))
	binomialTest(interactions)	

	interactions <- data.frame(chr1=c("chr1", "chr1","chr1", "chr2", "chr2", "chr3"), 
			locus1=c(1,1,1,2,2,3), 
			chr2=c("chr1", "chr2", "chr3", "chr2", "chr3", "chr3"), 
			locus2=c(1,2,3,2,3,3),
			frequencies=c(1,6,10,5,14,9))
	binomialTest(interactions)	
	
	interactions <- data.frame(chr1=c("chr1","chr1", "chr2"), 
			locus1=c(1,1,2), 
			chr2=c("chr2", "chr3","chr3"), 
			locus2=c(2,3,3),
			frequencies=c(6,10,14))
	binomialTest(interactions)	
}
