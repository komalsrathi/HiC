# prepare input file
# filename.bam is the output of HiCUP
# samtools view -h filename.bam | grep -v "^@" | cut -f -4 > filename.txt
# This will be imported as HICUPFILENAME

# import scripts
source('GenomicInteraction.R')
source('callInteractions.R')
source('binomial_method_HiC.R')
source('importHicup.R')
source('mapReadsToHindIII.R')
library(GOTHiC)

###################################### function definitions ###########################
.binInteractions <- GOTHiC:::.binInteractions
.countDuplicates <- GOTHiC:::.countDuplicates
getHindIIIsitesFromHicup <- function(filename, asRanges = TRUE)
{
  sites = read.table(filename, stringsAsFactor=FALSE, header=TRUE, skip=1)
  sites$chr <- sub("^", "chr", sites$Chromosome)
  sites$chr <- sub("chrCHR", "chr", sites$chr) #sometimes there is CHR
  sites$chr <- sub("chrchr", "chr", sites$chr) #sometimes there is chr***
  sites$start <- sites$Fragment_Start_Position
  sites$locus <- sites$Fragment_Start_Position
  sites$end <- sites$Fragment_End_Position
  sites <- sites[, c("chr", "locus", "start", "end")]
  
  if(asRanges)
  {
    library(GenomicRanges)
    return(GRanges(seqnames=sites$chr, ranges=IRanges(start=sites$start,
                                                      end=sites$end)))
  }
}
###################################### function definitions ###########################

HICUPFILENAME <- 'filename.txt'
fragmentFile <- 'filename of HiCUP digest file'
interactions <- importHicup(HICUPFILENAME)
interactions <- mapHicupToRestrictionFragment(interactions, fragmentFile)
RESOLUTION <- 1
interactions <- .binInteractions(interactions, RESOLUTION)
ncores <- 10
interactions <- binomialHiCagg(interactions, parallel=TRUE, cores=ncores)

