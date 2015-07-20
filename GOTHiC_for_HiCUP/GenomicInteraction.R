#S4 class to store an interaction matrix between genomic loci
#
#Copyright: Robert Sugar (robert.sugar@ebi.ac.uk), EMBL-EBI, 2012

#Handling elementMetadata:
#GenomicInteraction contains its own elementMetadata that is information about the interaction
#on top of that source and target contain their own elementMetadata that refers to the metadata about the source or target only

library(GenomicRanges)

setClass(Class="GenomicInteraction",
		  representation(
		  source = "GRanges", #one side of the interaction
		  target = "GRanges", #the other side of the interaction
		  frequency = "numeric",	#how many time it is seen
		  elementMetadata = "DataFrame"
		  )
 		)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###



setValidity("GenomicInteraction",  function(object)
		{
			#check lengths
			sourceLength = length(object@source)
			otherLengths = c(length(object@target), length(object@frequency), nrow(object@elementMetadata))	
			if(!all(sourceLength == otherLengths))
				return("fields have different lengths")
			
			#frequencies should be non-negative
			if(min(object@frequency) < 0)
				return("frequencies have to be non-negative")
			
#			if (!is.null(rownames(object@elementMetadata)))
#				return("slot 'elementMetadata' cannot contain row names")
			
			TRUE
		}
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###		
GenomicInteraction <- function(source, target, frequency = rep(1, length(source)), elementMetadata = NULL, ...)  
{
	#"..." provides additional columns in elementMetadata - but only if the elementMetadata itself is NULL
	if (is.null(elementMetadata)) 
	{
		elementMetadata <- DataFrame(...)
		if (ncol(elementMetadata) == 0L)
			elementMetadata <- new("DataFrame", nrows = length(source))
	}
	
	if(is.data.frame(elementMetadata))
		elementMetadata <- DataFrame(elementMetadata)

	#TODO remark: GRanges does not allow rownames in elementMetadata - we do here, because there is no "names" field
#	if (!is.null(rownames(elementMetadata))) {
#		rownames(elementMetadata) <- NULL	
	
	new("GenomicInteraction", source=source, target=target, frequency=frequency, elementMetadata=elementMetadata)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters.
###

setGeneric("seqnames1", function(x) standardGeneric("seqnames1"))
setMethod("seqnames1", "GenomicInteraction", function(x) x@source@seqnames)
setGeneric("chr1", function(x) standardGeneric("chr1"))
setMethod("chr1", "GenomicInteraction", function(x) x@source@seqnames)
setGeneric("ranges1", function(x, ...) standardGeneric("ranges1", ...))
setMethod("ranges1", "GenomicInteraction", function(x, ...) x@source@ranges)
setGeneric("strand1", function(x) standardGeneric("strand1"))
setMethod("strand1", "GenomicInteraction", function(x) x@source@strand)
setGeneric("seqinfo1", function(x) standardGeneric("seqinfo1"))
setMethod("seqinfo1", "GenomicInteraction", function(x) x@source@seqinfo)

setGeneric("start1", function(x, ...) standardGeneric("start1", ...))
setMethod("start1", "GenomicInteraction", function(x, ...) start(ranges(x@source)))
setGeneric("locus1", function(x, ...) standardGeneric("locus1", ...))
setMethod("locus1", "GenomicInteraction", function(x, ...) start(ranges(x@source)))
setGeneric("end1", function(x, ...) standardGeneric("end1", ...))
setMethod("end1", "GenomicInteraction", function(x, ...) end(ranges(x@source)))
setGeneric("width1", function(x) standardGeneric("width1"))
setMethod("width1", "GenomicInteraction", function(x) width(ranges(x@source)))
setGeneric("names1", function(x) standardGeneric("names1"))
setMethod("names1", "GenomicInteraction", function(x) names(x@source))

setGeneric("seqnames2", function(x) standardGeneric("seqnames2"))
setMethod("seqnames2", "GenomicInteraction", function(x) x@target@seqnames)
setGeneric("chr2", function(x) standardGeneric("chr2"))
setMethod("chr2", "GenomicInteraction", function(x) x@target@seqnames)
setGeneric("ranges2", function(x, ...) standardGeneric("ranges2", ...))
setMethod("ranges2", "GenomicInteraction", function(x, ...) x@target@ranges)
setGeneric("strand2", function(x) standardGeneric("strand2"))
setMethod("strand2", "GenomicInteraction", function(x) x@target@strand)
setGeneric("seqinfo2", function(x) standardGeneric("seqinfo2"))
setMethod("seqinfo2", "GenomicInteraction", function(x) x@target@seqinfo)


setGeneric("start2", function(x, ...) standardGeneric("start2", ...))
setMethod("start2", "GenomicInteraction", function(x, ...) start(ranges(x@target)))
setGeneric("locus2", function(x, ...) standardGeneric("locus2", ...))
setMethod("locus2", "GenomicInteraction", function(x, ...) start(ranges(x@target)))
setGeneric("end2", function(x, ...) standardGeneric("end2", ...))
setMethod("end2", "GenomicInteraction", function(x, ...) end(ranges(x@target)))
setGeneric("width2", function(x) standardGeneric("width2"))
setMethod("width2", "GenomicInteraction", function(x) width(ranges(x@target)))
setGeneric("names2", function(x) standardGeneric("names2"))
setMethod("names2", "GenomicInteraction", function(x) names(x@target))


setGeneric("from", function(x) standardGeneric("from"))
setMethod("from", "GenomicInteraction", function(x) x@source)
setGeneric("to", function(x) standardGeneric("to"))
setMethod("to", "GenomicInteraction", function(x) x@target)
setGeneric("frequency", function(x) standardGeneric("frequency"))
setMethod("frequency", "GenomicInteraction", function(x) x@frequency)
#setGeneric("elementMetadata", function(x) standardGeneric("elementMetadata"))
setMethod("elementMetadata", "GenomicInteraction", function(x) x@elementMetadata)
setGeneric("values", function(x) standardGeneric("values"))
setMethod("values", "GenomicInteraction", function(x) x@elementMetadata)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###
setAs("data.frame", "GenomicInteraction", function(from) 
		{		
			if(nrow(from) == 0)
				return(GenomicInteraction(GRanges(), GRanges()))
			
			#strand
			strand1 = from[['strand1']]
			if(is.null(strand1)) strand1 <- Rle("*", nrow(from))
			strand2 = from[['strand2']]
			if(is.null(strand2)) strand2 <- Rle("*", nrow(from))
			
			if(all(c("chr1", "locus1", "chr2", "locus2") %in% names(from)))
			{				
				source = GRanges(seqnames =	Rle(from$chr1), 
						ranges = IRanges(from$locus1, width = 1, names = from[['names1']]),
						strand = strand1)
				target = GRanges(seqnames =	Rle(from$chr2), 
						ranges = IRanges(from$locus2, width = 1, names = from[['names2']]),
						strand = strand2)
				
				names(from) <- sub("frequencies", "frequency", names(from))
			} else if(all(c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2") %in% names(from)))
			{
				#put all fields ending with "1" in source, all ending with "2" in target, frequency in "frequency" and the rest in elementMetadata
				
				#source
				source <- GRanges(seqnames = Rle(from$seqnames1),
						ranges = IRanges(start = from$start1, end =  from$end1, names = from[['names1']]),
						strand = strand1)												
				
				#target				
				target <- GRanges(seqnames = Rle(from$seqnames2),
						ranges = IRanges(start = from$start2, end = from$end2, names = from[['names2']]),
						strand = strand2
				)				
			} else
				stop('data frame should contain fields: "chr1", "locus1", "chr2", "locus2" OR "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"')
			
			#frequency
			if('frequency' %in% names(from))
			{
				frequency <- from$frequency
			} else
			{
				frequency <- rep(1, nrow(from))
			}
			
			#metadata1
			fieldsEndingIn1 <- names(from)[grep("1$", names(from))]
			metadataFields1 <- fieldsEndingIn1[!(fieldsEndingIn1 %in% #reserved fields are not metadata
								c("seqnames1", "start1", "end1", "width1", "strand1", "names1", "element1", "elementMetadata1", "chr1", "locus1"))]							
			sourceMetadata <- DataFrame(from[, metadataFields1])
			names(sourceMetadata) <- sub('1$', '', metadataFields1)
			values(source) <- sourceMetadata 
			gi <- GenomicInteraction(source, target, frequency)
			
			#metadata2
			fieldsEndingIn2 <- names(from)[grep("2$", names(from))]
			#reserved fields are not metadata
			metadataFields2 <- fieldsEndingIn2[!(fieldsEndingIn2 %in% 
								c("seqnames2", "start2", "end2", "width2", "strand2", "names2", "element2", "elementMetadata2", "chr2", "locus2"))]				
			targetMetadata <- DataFrame(from[, metadataFields2])
			names(targetMetadata) <- sub('2$', '', metadataFields2)
			values(target) <- targetMetadata 
			
			#elementMetadata - all fields left out so far: fields not ending with 1 or 2 AND not the field 'frequency'
			elementMetadataFields <- names(from)[grep('[12]$|^frequency$', names(from), , invert = TRUE)]
			elementMetadata <- DataFrame(from[, elementMetadataFields])
			names(elementMetadata) <- elementMetadataFields
			
			gi <- GenomicInteraction(source, target, frequency, elementMetadata)
			
			return(gi)				
		})

setMethod("as.data.frame", "GenomicInteraction",
		function(x, row.names=NULL, elementMetadata=TRUE, chr_locus_format = FALSE, ...)
		{
			
			if(chr_locus_format)
			{
				result = data.frame(chr1=as.factor(seqnames(x@source)),
						locus1=start(x@source),
						chr2=as.factor(seqnames(x@target)),
						locus2=start(x@target),
						frequencies = x@frequency,
						row.names=row.names,
						stringsAsFactors=FALSE)
				return(result)
			} else
			{
				elementMetadata1 <- as.data.frame(values(x@source))
				if(ncol(elementMetadata1) > 0)	names(elementMetadata1) <- paste(names(elementMetadata1), '1', sep = "")
				elementMetadata2 <- as.data.frame(values(x@target))
				if(ncol(elementMetadata2) > 0) names(elementMetadata2) <- paste(names(elementMetadata2), '2', sep = "")
				
				result <- data.frame(seqnames1=as.factor(seqnames(x@source)),
						start1=start(x@source),
						end1=end(x@source),
						width1=width(x@source),
						strand1=as.factor(strand(x@source)),
						elementMetadata1,
						seqnames2=as.factor(seqnames(x@target)),
						start2=start(x@target),
						end2=end(x@target),
						width2=width(x@target),
						strand2=as.factor(strand(x@target)),
						elementMetadata2,
						frequency = x@frequency,
						as.data.frame(x@elementMetadata),
						row.names=row.names,
						stringsAsFactors=FALSE)
				
				names1 <- names(x@source)
				if(!is.null(names1)) result$names1 <- names1
				names2 <- names(x@source)
				if(!is.null(names2)) result$names2 <- names2	
				
				return(result)
			}
		}
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.
### For an object with a pure S4 slot representation, these both map to
### initialize. Reference classes will want to override 'update'. Other
### external representations need further customization.

#setMethod("update", "GenomicRanges",
#		function(object, ...)
#		{
#			initialize(object, ...)
#		}
#)

#setGeneric("clone", function(x, ...) standardGeneric("clone"))

#setMethod("clone", "ANY",
#		function(x, ...)#
#		{
#			initialize(x, ...)
#		}
#)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Methods
###
setMethod("length", "GenomicInteraction", function(x) {
			#all lengths should be equal, as checked by setValidity
			length(x@source)
		})

#returns fields
setMethod("$", "GenomicInteraction", function(x, name)
		{
			#known fields
			knownFields <- c("seqnames1", "start1", "end1", "width1", "strand1", "names1", "element1", "elementMetadata1", "chr1", "locus1",
					"seqnames2", "start2", "end2", "width2", "strand2", "names2", "element2", "elementMetadata2", "chr2", "locus2",
					"from", "to", "frequency", "frequencies", "elementMetadata", "values")
			
			if(name %in% knownFields)
			{
				return(eval(call(name, x)))
			}
			
			#metadata for source
			if(length(grep('1$', name)) > 0)
			{
				metadataField <- sub('1$', '', name)
				return(x@source@elementMetadata[, metadataField])
			}
			
			#metadata for target
			if(length(grep('2$', name)) > 0)
			{
				metadataField <- sub('2$', '', name)
				return(x@target@elementMetadata[, metadataField])
			}
			
			#elementMetadata
			return(x@elementMetadata[, name])
			
			return(NULL)
		})

		
setMethod("$<-", "GenomicInteraction", function(x, name, value)
{
	#known fields
	knownFields <- c("seqnames1", "start1", "end1", "width1", "strand1", "names1", "element1", "elementMetadata1", "chr1", "locus1",
			"seqnames2", "start2", "end2", "width2", "strand2", "names2", "element2", "elementMetadata2", "chr2", "locus2",
			"from", "to", "frequency", "frequencies", "elementMetadata", "values")	

	#TODO: write replacement methods
	if(name %in% knownFields)
	{
		return(eval(call(paste(name, '<-', sep = ''), x)))
	} else
	
	#metadata for source
	if(length(grep('1$', name)) > 0)
	{
		metadataField <- sub('1$', '', name)
		x@source@elementMetadata[, metadataField] <- value
	} else
	
	#metadata for target
	if(length(grep('2$', name)) > 0)
	{
		metadataField <- sub('2$', '', name)
		x@target@elementMetadata[, metadataField] <- value
	} else
	#elementMetadata
	{

		x@elementMetadata[, name] <- value
	}
	return(x)
})

giRange <- function(gi, range) GenomicInteraction(gi@source[range], gi@target[range], gi@frequency[range], gi@elementMetadata[range, ])

#giApply <- function(gi, fun) GenomicInteraction(gi@source[range], gi@target[range], gi@frequency[range], gi@elementMetadata[range, ])

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test
###

testGenomicInteractions <- function()
{
	gr1 <- GRanges(seqnames =
							Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
					ranges =
							IRanges(1:10, width = 10:1, names = head(letters,10)),
					strand =
							Rle(strand(c("-", "+", "*", "+", "-")),
									c(1, 2, 2, 3, 2)),
					meta = 10:1
			)
			
	gr2 <- GRanges(seqnames =
					Rle(c("chr1", "chr1", "chr3", "chr5"), c(1, 3, 2, 4)),
			ranges =
					IRanges(seq(1, 20, 2), width = seq(20, 1, -2), names = tail(letters,10)),
			strand =
					Rle(strand(c("-", "+", "*", "+", "-")),
							c(1, 2, 2, 3, 2)),
			meta = 1:10
	)
	
	metadata = letters[10:1]
	
	gi <- GenomicInteraction(gr1, gr2, 0:9, meta = metadata)
	gi_with_dataFrame <- GenomicInteraction(gr1, gr2, 0:9, data.frame(meta = metadata))
	
	#error cases
	gi_wrong_class <- GenomicInteraction(gr1, 0:9, 0:9)
	gi_too_many_frequencies <- GenomicInteraction(gr1, gr2, 0:10)
	gi_negative_frequencies <- GenomicInteraction(gr1, gr2, -1:8)
	gi_no_target <- GenomicInteraction(gr1, NULL, 0:9)
	gi_wrong_metadata <- GenomicInteraction(gr1, gr2, 0:9, data.frame(meta = metadata[1:5]))
	gi_wrong_metadata2 <- GenomicInteraction(gr1, gr2, 0:9, meta = metadata[1:5])
	gi_wrong_metadata3 <- GenomicInteraction(gr1, gr2, 0:9, list(meta = metadata))
	#TODO: handle when elementMetadata in the two GRanges have the same name... maybe source has preference?
	
	#coercion
	gi_dataframe <- as.data.frame(gi)
	gi_dataframe
	
	gi_chr_locus <- as.data.frame(gi, chr_locus_format = TRUE)
	gi_chr_locus
	
	gi2 <- as(gi_chr_locus, "GenomicInteraction")
	gi2
	
	gi2 <- as(gi_dataframe, "GenomicInteraction")
	gi2
	
	#TODO: check that the conversion is ok
	
	#getters
	f.list <- showMethods(classes="GenomicInteraction", printTo=FALSE)
	functions <- sapply(
			strsplit( f.list[grep("Function",f.list)],' '),
			function(l){gsub('\"|:','',l[[2]])}
	)
	functions = functions[grep("*[12]$", functions)]
	functions = c(functions, "from", "to", "frequency", "elementMetadata", "values")
	
	result <- sapply(functions, function(x)
			{ 
				print(paste("calling", x))
				print(eval(call(x, gi)))
			})
}
