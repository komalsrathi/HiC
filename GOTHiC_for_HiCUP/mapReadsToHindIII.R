# TODO: Add comment
# 
# Author: borbalagerle
###############################################################################

processMappabilityScores <- function(chromosomeNames = mouseChromosomes, bins =  getHindIIIsites(noOverlap=TRUE, pointRanges=FALSE), directory = ".", maxFragmentLength = 600, maxlength = 250e6)
{
	#take +/- maxFragmentLength from restriction sites - make sure the ranges don't overlap with the previous/next restriction site boundary
#	if(is.finite(maxFragmentLength) && maxFragmentLength > 0)
#	{
#		prevstart <- c(0, start(bins)[-length(start(bins))])
#		prevstart[(diff(prevstart) < 0)] = 0 #deal with chromosome boundaries
##		nextstart <- c(start(bins)[-1], Inf)
##		nextstart[(which(diff(nextstart) < 0)) + 1] = Inf #deal with chromosome boundaries
#		bins@ranges = IRanges(start = pmax(start(bins) - 600, prevstart, 1), end = pmin(start(bins) + 600, end(bins)))		
#	}
	
	bins@elementMetadata$mappability = NA
	for(chrom in chromosomeNames)
	{
		f <- file(paste(directory, "/", chrom, "b.out", sep = ""), "rb")
		mappability <- readBin(f, "raw", n = maxlength)
		mappability <- as.integer(as.integer(mappability) == 01)
		
		#bin it
		sitesOnChrom <- as.character(seqnames(bins)) == chrom
		cummap <- cumsum(mappability)
		
		if(is.finite(maxFragmentLength) && maxFragmentLength > 0)
		{
			#calculate (sum) the mappability maxFragmentLength inwards from start to end
			from_start_inwards <- pmin(start(bins[sitesOnChrom]) + maxFragmentLength, end(bins[sitesOnChrom]))
			mappability_start <-  cummap[from_start_inwards] - cummap[start(bins[sitesOnChrom])]
			from_end_inwards <- pmax(end(bins[sitesOnChrom]) - maxFragmentLength, start(bins[sitesOnChrom]))
			mappability_end <- cummap[end(bins[sitesOnChrom])] - cummap[from_end_inwards]
			
			bins@elementMetadata$mappability[sitesOnChrom] <- mappability_start + mappability_end
		} else
		{
			bins@elementMetadata$mappability[sitesOnChrom] <- cummap[end(bins[sitesOnChrom])] - cummap[start(bins[sitesOnChrom])]
		}
		
		close(f)
	}
	#bins[which(is.na(bins@elementMetadata$mappability))]@elementMetadata$mappability = 0
	return(bins)
}

processGcContent <- function(chromosomeNames = mouseChromosomes, organism = 'BSgenome.Mmusculus.UCSC.mm9', species = 'Mmusculus', bins =  getHindIIIsites(noOverlap=TRUE, pointRanges=FALSE), directory = ".", windowsize = 200)
{
	library(BSgenome)
	library(organism, character.only=TRUE)
	genome= get(species)
	
	bins@elementMetadata$gc = NA
	for(chrom in chromosomeNames)
	{	
		#bin it
		sitesOnChrom <- as.character(seqnames(bins)) == chrom
		
		if(is.finite(windowsize) && windowsize > 0)
		{
			#calculate the gc content windowsize inwards from start to end
			from_start_inwards <- pmin(start(bins[sitesOnChrom]) + windowsize, end(bins[sitesOnChrom]))
			view_start = Views(genome[[chrom]], start = start(bins[sitesOnChrom]), end = from_start_inwards)
			gc_start <-  alphabetFrequency(view_start)
			from_end_inwards <- pmax(end(bins[sitesOnChrom]) - windowsize, start(bins[sitesOnChrom]))
			view_end = Views(genome[[chrom]], start = from_end_inwards, end = end(bins[sitesOnChrom]))
			gc_end <-  alphabetFrequency(view_end)
			gc_all <- (gc_start + gc_end)
			gccontent <- (gc_all[, 'G'] + gc_all[, 'C']) / (gc_all[, 'G'] + gc_all[, 'C'] + gc_all[, 'A'] + gc_all[, 'T'])
			
			bins@elementMetadata$gc[sitesOnChrom] <- gccontent
		} else stop("invalid window size")
	}
	#bins[which(is.na(bins@elementMetadata$mappability))]@elementMetadata$mappability = 0
	return(bins)
}

getHindIIIsites <- function(organism = "mouse", pointRanges = TRUE, noOverlap = TRUE)
{
	#organism: mouse and human 
	#pointRanges: if TRUE the range start and end are the same. If FALSE, the ranges between restriction sites are returned. Note: first and last restriction fragments of the chromosomes are truncated/ignored respectively
	
	#load sorted paired alignments in GRanges
	load(paste(BASEDATADIR, paste(organism, "HindIIIsites.Rd", sep = ""), sep = "/"))
	
	if(pointRanges) 
	{
		return(GRanges(seqnames=allHindIIIsites$chr, ranges=IRanges(start=allHindIIIsites$locus, end=allHindIIIsites$locus)))
	} else
	{
		seqnames = c(allHindIIIsites$chr[1], allHindIIIsites$chr)
		start = c(allHindIIIsites$locus[1], allHindIIIsites$locus)
		end = c(allHindIIIsites$locus, allHindIIIsites$locus[length(allHindIIIsites$locus)])
		
		if(noOverlap) end <- end - 1 #ranges are [100, 199] [200, 299] instead of [100, 200] [200 300]
		
		#fix ranges that lie between two chromosomes
		chromosomeBeginningsOrEnds = which(start > end)
		end[chromosomeBeginningsOrEnds] = start[chromosomeBeginningsOrEnds]
		
		ranges = GRanges(seqnames, ranges=IRanges(start, end))
		
		#remove zero-length fragments (e.g. beginning/end of chromosomes). TODO: should we include these fragments?
		return(ranges[-which(start(ranges) == end(ranges))])		
	}
}

getHindIIIsitesFromHicup <- function(filename, asRanges = TRUE)
{
	sites = read.table(filename, stringsAsFactor=FALSE, header=TRUE, skip=1)
	sites$chr <- sub("^", "chr", sites$Chromosome)
	sites$chr <- sub("chrCHR", "chr", sites$chr) #sometimes there is CHR
	sites$start <- sites$Fragment_Start_Position
	sites$locus <- sites$Fragment_Start_Position
	sites$end <- sites$Fragment_End_Position
	sites <- sites[, c("chr", "locus", "start", "end")]	
	
	if(asRanges)
	{
		library(GenomicRanges)
		return(GRanges(seqnames=sites$chr, ranges=IRanges(start=sites$start, end=sites$end)))
	}
}

checkRestrictionSiteBias <- function(paired_reads_1, paired_reads_2, sampleName, BINSIZES = 25000, removeOutliersOverPercentile = 0.99)
{
	library("GenomicRanges")
	
	load(paste(BASEDATADIR, "chromosomeLengths.Rd", sep = "/"))
	ranges <- getHindIIIsites()
	
	allReads <- c(paired_reads_1, paired_reads_2)

	pdf(paste(FIGURESDIR, "/restrictionSiteBias_", sampleName, ".pdf", sep = ""))	
	
	for(BINSIZE in BINSIZES)
	{
		resolution <- resolutionString(BINSIZE)
		
		#create bins as Granges
		binsPerChromosome = (chromosomeLengths %/% BINSIZE) + 1
		binLoci =  unlist(sapply(binsPerChromosome, function(x) seq(0, length = x, by = BINSIZE)))	
		bins = GRanges(seqnames =
						Rle(mouseChromosomes, binsPerChromosome),
				ranges = IRanges(c(binLoci), width = BINSIZE))
		
		rangesInBins = countOverlaps(bins, ranges)
		readsInBins = countOverlaps(bins, allReads)
		
		readsVsHindIII = data.frame(reads2 = rangesInBins, coverage = readsInBins)
		
		#remove 0s
		readsVsHindIII <- readsVsHindIII[readsVsHindIII$reads2 > 0 & readsVsHindIII$coverage > 0, ]	
	
		#remove outliers
		readsVsHindIII <- readsVsHindIII[readsVsHindIII$reads2 <= quantile(readsVsHindIII$reads2, removeOutliersOverPercentile) &
										 readsVsHindIII$coverage <= quantile(readsVsHindIII$coverage, removeOutliersOverPercentile), ]	
		
		rho = cor(readsVsHindIII$reads2, readsVsHindIII$coverage, method = "spearman")
		
		#smoothScatterPlot with regression:	
		smoothScatter(readsVsHindIII, main = paste("Restriction site bias in ", sampleName, "(resolution:", resolution, ")"))
		fit<-lm(readsVsHindIII$coverage ~ readsVsHindIII$reads2)
		abline(fit$coefficients[1], fit$coefficients[2])	
		lines(lowess(readsVsHindIII))	
		legend("topright", legend = paste("Spearman's rho =", round(rho, digits = 2)))	
		
		#plot(readsVsHindIII, main = paste("Restriction site bias (logarithmic axes) in ", sampleName), log = "xy")
		#fit<-lm(readsVsHindIII$coverage ~ readsVsHindIII$reads2)
		#abline(fit$coefficients[1], fit$coefficients[2])	
		#lines(lowess(readsVsHindIII))	
	}
	
	dev.off()
}

checkFragmentLengthBias <- function(paired_reads_1, paired_reads_2, sampleName, trimOutliers = c(0.05, 0.95))
{
	#trimOutliers: if the fragment length is below the first and above the second percentile value the fragments are not considered
	
	library("GenomicRanges")
	
	ranges <- getHindIIIsites(pointRanges = FALSE)
	fragmentLengths = end(ranges) - start(ranges)
	
	trimValues = quantile(fragmentLengths, trimOutliers)
	
	normalFragments = (fragmentLengths >= trimValues[1] & fragmentLengths <= trimValues[2])
	ranges = ranges[normalFragments]
	fragmentLengths = fragmentLengths[normalFragments]
	
	allReads <- c(paired_reads_1, paired_reads_2)
	
	overlaps = countOverlaps(ranges, allReads)
	
	pdf(paste(FIGURESDIR, "/restrictionFragmentLengthBias_", sampleName, ".pdf", sep = ""))	
	
	#---
	#fragment length vs. coverage
	#---
	rho = cor(fragmentLengths, overlaps, method = "spearman")
	#smoothScatterPlot with regression:	
	smoothScatter(fragmentLengths, overlaps, main = paste("Fragment length site bias in ", sampleName))
	fit<-lm(overlaps ~ fragmentLengths)
	abline(fit$coefficients[1], fit$coefficients[2])	
	lines(lowess(fragmentLengths, overlaps))	
	legend("topright", legend = paste("Spearman's rho =", round(rho, digits = 2)))
	
	#---
	#fragment length vs. coverage corrected for fragment length
	#---	
	overlaps = overlaps / fragmentLengths
	rho = cor(fragmentLengths, overlaps, method = "spearman")	
	#smoothScatterPlot with regression:	
	smoothScatter(fragmentLengths, overlaps, main = paste("Fragment length site bias in ", sampleName))
	fit<-lm(overlaps ~ fragmentLengths)
	abline(fit$coefficients[1], fit$coefficients[2])	
	lines(lowess(fragmentLengths, overlaps))	
	legend("topright", legend = paste("Spearman's rho =", round(rho, digits = 2)))
	
	
	dev.off()
}

toGranges <- function(reads)
{
	if(!all(c("chr1", "locus1", "chr2", "locus2")) %in% names(reads))
		
		reads_source <- reads[, c("chr1", "locus1")]
	names(reads_source) <- c("chr", "locus")	
	reads_target <- reads[, c("chr2", "locus2")]
	names(reads_target) <- c("chr", "locus")		
	all_reads <- rbind(reads_source, reads_target)
	all_reads$chr <- gsub("0(.)","\\1", all_reads$chr) #replace chr01 with chr1 to match original
	
	GRanges(seqnames =	Rle(all_reads$chr), 
			ranges = IRanges(all_reads$locus, width = 1))
	
}

checkCoverageCorrelation <- function(reads1, reads2, sampleName, sample1Name = "sample1", sample2Name = "sample2", BINSIZES = 1e6)#c(10000, 25000, 50000, 100000, 250000, 500000, 1e6, 2.5e6, 5e6, 1e7), trimOutliers = c(0, 1))
{
	library("GenomicRanges")	
	
	if(class(reads1) == "data.frame") reads1 = toGranges(reads1)
	if(class(reads2) == "data.frame") reads2 = toGranges(reads2)	
	
	load(paste(BASEDATADIR, "chromosomeLengths.Rd", sep = "/"))
	
	pdf(paste(FIGURESDIR, "/coverageCorrelation_", sampleName, "-", sample1Name, "_vs_", sample2Name, ".pdf", sep = ""))	
	
	correlations_pearson = list()
	correlations_spearman = list()
	
	for(BINSIZE in BINSIZES)
	{
		resolution <- resolutionString(BINSIZE)
		
		#create bins as Granges
		bins = getBins(chromosomeLengths, mouseChromosomes, BINSIZE)
		
		reads1InBins = countOverlaps(bins, reads1)
		reads2InBins = countOverlaps(bins, reads2)
		
		reads1Vsreads2 = data.frame(reads1 = reads1InBins, reads2 = reads2InBins)
		
		#remove 0s
		reads1Vsreads2 <- reads1Vsreads2[reads1Vsreads2$reads1 > 0 & reads1Vsreads2$reads2 > 0, ]	
			
		rho = cor(reads1Vsreads2$reads1, reads1Vsreads2$reads2, method = "spearman")	
		correlations_spearman[paste(sampleName, sample1Name, "vs", sample2Name, BINSIZE, sep = "_")] = rho
		print(paste("Spearman's rho for read correlation between", sample1Name, "and", sample2Name, ":", rho))
		pearson = cor(reads1Vsreads2$reads1, reads1Vsreads2$reads2, method = "pearson")
		correlations_pearson[paste(sampleName, sample1Name, "vs", sample2Name, BINSIZE, sep = "_")] = pearson
		print(paste("Pearson's r for read correlation between", sample1Name, "and", sample2Name, ":", pearson))
		
		#smoothScatterPlot with regression:	
		smoothScatter(reads1Vsreads2, main = paste("Correlation of coverages ", sampleName, "(resolution:", resolution, ")"))
		fit<-lm(reads1Vsreads2$reads1 ~ reads1Vsreads2$reads2)
		abline(fit$coefficients[1], fit$coefficients[2])	
		lines(lowess(reads1Vsreads2))	
		legend("topright", legend = paste("Spearman's rho =", round(rho, digits = 3), "\nPearson's r =", round(pearson, digits = 3)))	
		
		#plot(reads1Vsreads2, main = paste("Restriction site bias (logarithmic axes) in ", sampleName), log = "xy")
		#fit<-lm(reads1Vsreads2$coverage ~ reads1Vsreads2$reads2)
		#abline(fit$coefficients[1], fit$coefficients[2])	
		#lines(lowess(reads1Vsreads2))	
	}
	
	dev.off()
	
	return(list(pearson = correlations_pearson, spearman = correlations_spearman))
}


mapReadsToHindIII <- function(pairedReadsFile, sampleName = sub("_paired","", pairedReadsFile))
{
	library("ShortRead")

	#load sorted paired alignments in GRanges
	load(pairedReadsFile)
	hindIIIRanges <- getHindIIIsites()
	
	checkRestrictionSiteBias(paired_reads_1, paired_reads_2, sampleName, RESOLUTIONS)
	
	#find HindIII sites	
	hindIII_1 <- ifelse(as.vector(paired_reads_1%in%hindIIIRanges == TRUE),match(paired_reads_1,hindIIIRanges),precede(paired_reads_1, hindIIIRanges))  
	hindIII_2 <- ifelse(as.vector(paired_reads_2%in%hindIIIRanges == TRUE),match(paired_reads_2,hindIIIRanges),precede(paired_reads_2, hindIIIRanges))
	#take out reads where we can't assign hindIII site
	sort_reads_c1=paired_reads_1[is.na(hindIII_1)==FALSE]
	sort_reads_c2=paired_reads_2[is.na(hindIII_2)==FALSE]
	id1=as.character(elementMetadata(sort_reads_c1)$id)
	ids1 = sapply(strsplit(id1, '[/ ]'), '[[', 1) 
	rm(id1)
	id2=as.character(elementMetadata(sort_reads_c2)$id)
	ids2 = sapply(strsplit(id2, '[/ ]'), '[[', 1) 
	rm(id2)
	s1=which(ids1%in%ids2)
	s2=which(ids2%in%ids1)
	hindIII1=hindIII_1[is.na(hindIII_1)==FALSE]
	hindIII2=hindIII_2[is.na(hindIII_2)==FALSE]
	hindIII1=hindIII1[s1]
	hindIII2=hindIII2[s2]
			
	#fragment length bias
	checkFragmentLengthBias(paired_reads_1, paired_reads_2, sampleName)
	
	#TODO - remove sites that do not lie within X (~500) bps of a restriction site	- MAXFRAGMENTSIZE > d1 + d2 (distances between read and nearest restriction site)
	
	#assign hindIII restriction sites to reads
	loci1 <- hindIIIRanges[hindIII1]
	loci2 <- hindIIIRanges[hindIII2]
	interactingLoci <- GRangesList(locus1 = loci1, locus2 = loci2)
	outputfilename <- paste("interactingLoci", sampleName, sep = "_") #remove "paired", add "interactingLoci"
	save(interactingLoci, file = outputfilename)
	
	#calculate coverage
	all_reads <- c(paired_reads_1, paired_reads_2)
	exportCoverage(all_reads, sampleName)
	
	rm(paired_reads_1,paired_reads_2)
	rm(hindIII_1,hindIII_2,loci1,loci2,ids1,ids2,s1,s2)	
}

#maps Hicup reads to restriction sites
#interactions: data frame or GenomicInteraction
#fragments: GenomicRanges or a filename
mapHicupToRestrictionFragment <- function(interactions, fragments)
{
	if(is.data.frame(interactions))
	{
		interactions <- as(interactions, "GenomicInteraction")
	}
	
	if(is.character(fragments))
	{
		fragments <- getHindIIIsitesFromHicup(fragments)
	}
	
	firstFragment <- findOverlaps(interactions@source, fragments, select="first")
	secondFragment <- findOverlaps(interactions@target, fragments, select="first")
	
	#filter out non-mapping loci
	validInteractions <- !is.na(firstFragment) & !is.na(secondFragment)
	firstFragment <- firstFragment[validInteractions]
	secondFragment <- secondFragment[validInteractions]
	
	#make sure A-B and B-A are treated together -> put fragment with smaller ID in front
	sources <- pmin(firstFragment, secondFragment)
	targets <- pmax(firstFragment, secondFragment)
	
	gi <- GenomicInteraction(fragments[sources], fragments[targets])
	return(as.data.frame(gi, chr_locus_format = TRUE))
}

exportCoverage <- function(reads, sampleName)
{
	library("ShortRead")
	library(rtracklayer)
	coverage_reads = coverage(reads)
	export.bedGraph(as(coverage_reads,"RangedData"), paste("coverage_", sampleName, ".bed", sep =""), name=paste("coverage_", sampleName, sep =""))
}
