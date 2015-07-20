#Mouse Hi-C
#imports HiCUP output
#
#Copyright: Robert Sugar (robert.sugar@ebi.ac.uk), EMBL-EBI, 2012

fixChromosomeNames <- function(chrnames)
{
	#capital to small
	chrnames <- sub("CHR", "chr", chrnames)
	chrnames <- sub("", "chr", chrnames)
	chrnames <- sub("chrchr", "chr", chrnames)
}

importHicup <- function(filename, checkConsistency=TRUE, filetype=ifelse(grepl("\\.bam$", filename), "bam", "table"))
{
	#the output of hicup is a sam file, that looks like uniques_ORIGINALFILE_trunc.sam
	#this has to be converted using the hicupToTable tool
	
	library(GenomicRanges)

	if(filetype=="table")
	{
		tbl <- read.table(filename)
		#in the sam file the two ends of the interactions are consecutive rows	
		odd <- tbl[seq(1,nrow(tbl),2), ]
		names(odd) <- c("id1", "flag1", "chr1", "locus1")	
		levels(odd$chr1) <- fixChromosomeNames(levels(odd$chr1))
		even <- tbl[seq(2,nrow(tbl),2), ]
		names(even) <- c("id2", "flag2", "chr2", "locus2")	
		levels(even$chr2) <- fixChromosomeNames(levels(even$chr2))
	} else if(filetype=="bam")
	{
		library(ShortRead)
		tbl <- readAligned(filename, type="BAM")
		tbl <- as(tbl, "GRanges")

		#in the bam file the two ends of the interactions are consecutive rows	
		odd <- tbl[seq(1,length(tbl),2)]
		odd <- data.frame(chr1=as.vector(odd@seqnames), locus1=as.integer(odd@ranges@start), id1=as.character(odd$id))
		levels(odd$chr1) <- fixChromosomeNames(levels(odd$chr1))
		even <- tbl[seq(2,length(tbl),2)]
		even <- data.frame(chr2=as.vector(even@seqnames), locus2=as.integer(even@ranges@start), id2=as.character(even$id))
		levels(even$chr2) <- fixChromosomeNames(levels(even$chr2))
	}

	joined <- cbind(odd, even)
	
	if(nrow(odd) != nrow(even)) stop('importHicup: reads must be paired in consecutive rows')

	if(checkConsistency)
	{
		if(!all(joined$id1==joined$id2)) stop('importHicup: reads must be paired in consecutive rows')
	}
	
	return(joined[, c("chr1", "locus1", "chr2", "locus2")])
}

aggregateInteractions <- function(interactions, deduplicate=FALSE)
{
	#sort it
	interactions <- interactions[order(interactions$chr1, interactions$locus1, interactions$chr2, interactions$locus2), ]	
	#bin it
	interactions <- countDuplicates(interactions)

	if(deduplicate)
	{
		#set unmapped frequencies to 1
		print(paste("deduplication: removing", sum(interactions$frequencies) - nrow(interactions), "from", sum(interactions$frequencies)))		
		interactions$frequencies <- rep(1, length(interactions$frequencies))
	}

	return(interactions)
}
