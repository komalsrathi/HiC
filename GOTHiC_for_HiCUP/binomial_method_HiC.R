# TODO: Add comment
# 
# Author: borbalagerle
###############################################################################


# take GenomicRangesList with interactions in locus1 and locus2
binomialHiC=function(interactingLoci, res, removeDiagonal=TRUE,cistrans='all',filterdist=10000)
{
	library(data.table)
	#filter for reads that are closer than 10kb to remove self ligations 
	df_int <- data.frame(as.vector(seqnames(interactingLoci[[1]])), start(ranges(interactingLoci[[1]])), as.vector(seqnames(interactingLoci[[2]])), start(ranges(interactingLoci[[2]])))
	colnames(df_int) <- c("chr1", "locus1", "chr2", "locus2")
	df_filtered <- df_int
	df_filtered$dist <-as.vector(ifelse(df_filtered$chr1 == df_filtered$chr2, abs(df_filtered$locus1 - df_filtered$locus2), Inf))
	df_filtered <-df_filtered[df_filtered$dist>filterdist,]
	df_filtered <-df_filtered[,1:4]
	#bin interactions according to resolution
	binned_df_filtered <- binInteractions(df_filtered, res)
	binned_df_filtered$int1 <-paste(binned_df_filtered$chr1,binned_df_filtered$locus1,sep='_')
	binned_df_filtered$int2 <-paste(binned_df_filtered$chr2,binned_df_filtered$locus2,sep='_')
	#diagonal removal
	if(removeDiagonal)
	{
		subs <-which(binned_df_filtered$int1!=binned_df_filtered$int2)
		binned_df_filtered <- binned_df_filtered[subs,]	
	}
	if(cistrans=='cis'){
		subs <-which(binned_df_filtered$chr1==binned_df_filtered$chr2)
		binned_df_filtered <- binned_df_filtered[subs,]	
	}
	if(cistrans=='trans'){
		subs <-which(binned_df_filtered$chr1!=binned_df_filtered$chr2)
		binned_df_filtered <- binned_df_filtered[subs,]	
	}

	#all read pairs used in binomial
	numberOfReadPairs <- sum(binned_df_filtered$frequencies)
	#calculate coverage 
	all_bins <- unique(c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2)))
	#all_bins <- c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2))
	all_bins <- sort(all_bins)
	
	binned_dt=data.table(binned_df_filtered)
	covA <- binned_dt[,sum(c(frequencies)),by=int1]	
	covB <- binned_dt[,sum(c(frequencies)),by=int2]
	covA <- setkey(covA,key='int1')
	colnames(covB)[1] <- 'int1'
	covB <- setkey(covB,key='int1')
	cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
	cov$V1.x[is.na(cov$V1.x)]=0
	cov$V1.y[is.na(cov$V1.y)]=0
	cov$coverage=cov$V1.x+cov$V1.y
	coverage=cov$coverage
	names(coverage)=cov$int1
	
	sumcov <- sum(coverage)
	relative_coverage <- coverage/sumcov
	names(relative_coverage)=names(coverage)
	binned_df_filtered$cov1 <- relative_coverage[binned_df_filtered$int1]
	binned_df_filtered$cov2 <- relative_coverage[binned_df_filtered$int2]
	#probability correction assuming on average equal probabilities for all interactions
	numberOfAllInteractions <- length(all_bins)^2
	upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
	chromos <- unique(binned_df_filtered$chr1)
	chrlens <- c()
	for(cr in chromos){ 
		chrlens[cr] <- max(length(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),length(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
	}
	cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
	transBinNumber <- upperhalfBinNumber-cisBinNumber
	
	diagonalProb <- sum(relative_coverage^2)
	if(cistrans=='all'){
		probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
	}
	if(cistrans=='cis'){
		probabilityCorrection <- upperhalfBinNumber/cisBinNumber
	}
	if(cistrans=='trans'){
		probabilityCorrection <- upperhalfBinNumber/transBinNumber
	}
	
	
	binned_df_filtered$probability <- binned_df_filtered$cov1*binned_df_filtered$cov2*2*probabilityCorrection

	binned_df_filtered$predicted <- binned_df_filtered$probability * numberOfReadPairs
	
	binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
			{
				binom.test(as.numeric(x[["frequencies"]]), numberOfReadPairs, as.numeric(x[["probability"]]), alternative = "greater")$p.value
			}	
	)
	
	#observed over expected log ratio
	binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
	#multiple testing correction separately for matrices with all interactions/only cis/only transs
	
	if(cistrans=='all'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber+length(all_bins))}
	}
	if(cistrans=='cis'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber+length(all_bins))}	
	}
	if(cistrans=='trans'){
		binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=transBinNumber)	
	}
	
	return(binned_df_filtered)
}

###############################################################################################
# take hicupfile mapped to HindIII sites
res=1
binomialHiChicup=function(hicupinteraction, res, removeDiagonal=TRUE,cistrans='all',filterdist=10000)
{
	library("data.table")
	#filter for reads that are closer than 10kb to remove self ligations 
	df_int <- hicupinteraction[order(hicupinteraction$chr1, hicupinteraction$locus1, hicupinteraction$chr2, hicupinteraction$locus2), ]
	df_filtered <- df_int
	df_filtered$dist <-as.vector(ifelse(df_filtered$chr1 == df_filtered$chr2, abs(df_filtered$locus1 - df_filtered$locus2), Inf))
	df_filtered <-df_filtered[df_filtered$dist>filterdist,]
	df_filtered <-df_filtered[,1:4]
	#bin interactions according to resolution
	if(res>1){
		binned_df_filtered <- binInteractions(df_filtered, res)
	}
	if(res==1){
		binned_df_filtered <- countDuplicates(df_filtered)
	}
	
	binned_df_filtered$int1 <-paste(binned_df_filtered$chr1,binned_df_filtered$locus1,sep='_')
	binned_df_filtered$int2 <-paste(binned_df_filtered$chr2,binned_df_filtered$locus2,sep='_')
	#diagonal removal
	if(removeDiagonal)
	{
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$int1!=binned_df_filtered$int2,]	
	}
#	if(promoters!=NULL){
#		binned_df_filteredGR=GRangesList()
#		binned_df_filteredGR$locus1=GRanges(seqnames=binned_df_filtered$chr1,ranges=IRanges(start=binned_df_filtered$locus1,end=binned_df_filtered$locus1))
#		binned_df_filteredGR$locus2=GRanges(seqnames=binned_df_filtered$chr2,ranges=IRanges(start=binned_df_filtered$locus2,end=binned_df_filtered$locus2))
#		subs <- which(binned_df_filtered)
#	}
	
	if(cistrans=='cis'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1==binned_df_filtered$chr2,]	
	}
	if(cistrans=='trans'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1!=binned_df_filtered$chr2,]	
	}
	
	#all read pairs used in binomial
	numberOfReadPairs <- sum(binned_df_filtered$frequencies)
	#calculate coverage 
	all_bins <- unique(c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2)))
	all_bins <- sort(all_bins)
	binned_dt=data.table(binned_df_filtered)
	covA <- binned_dt[,sum(c(frequencies)),by=int1]	
	covB <- binned_dt[,sum(c(frequencies)),by=int2]
	covA <- setkey(covA,key='int1')
	colnames(covB)[1] <- 'int1'
	covB <- setkey(covB,key='int1')
	cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
	cov$V1.x[is.na(cov$V1.x)]=0
	cov$V1.y[is.na(cov$V1.y)]=0
	cov$coverage=cov$V1.x+cov$V1.y
	coverage=cov$coverage
	names(coverage)=cov$int1
	
	sumcov <- sum(coverage)
	relative_coverage <- coverage/sumcov
	names(relative_coverage)=names(coverage)
	binned_df_filtered$cov1 <- relative_coverage[binned_df_filtered$int1]
	binned_df_filtered$cov2 <- relative_coverage[binned_df_filtered$int2]
	#probability correction assuming on average equal probabilities for all interactions
	numberOfAllInteractions <- length(all_bins)^2
	upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
	
	if(cistrans!='all'){
	chromos <- unique(binned_df_filtered$chr1)
	chrlens <- c()
	for(cr in chromos){ 
		chrlens[cr] <- max(length(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),length(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
	}
	cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
	transBinNumber <- upperhalfBinNumber-cisBinNumber
	}
	
	diagonalProb <- sum(relative_coverage^2)
	if(cistrans=='all'){
		probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
	}
	if(cistrans=='cis'){
		probabilityCorrection <- upperhalfBinNumber/cisBinNumber
	}
	if(cistrans=='trans'){
		probabilityCorrection <- upperhalfBinNumber/transBinNumber
	}
	
	
	binned_df_filtered$probability <- binned_df_filtered$cov1*binned_df_filtered$cov2*2*probabilityCorrection
	
	binned_df_filtered$predicted <- binned_df_filtered$probability * numberOfReadPairs
	
	binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
			{
				binom.test(as.numeric(x[["frequencies"]]), numberOfReadPairs, as.numeric(x[["probability"]]), alternative = "greater")$p.value
			}	
	)
	#observed over expected log ratio
	binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
	#multiple testing correction separately for matrices with all interactions/only cis/only transs
	
	if(cistrans=='all'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber+length(all_bins))}
	}
	if(cistrans=='cis'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber+length(all_bins))}	
	}
	if(cistrans=='trans'){
		binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=transBinNumber)	
	}
	
	return(binned_df_filtered)
}

################# coverage from random ligations ################

binomialHiChicupRL=function(hicupinteraction, RLinteraction, res, removeDiagonal=TRUE,cistrans='all',filterdist=10000)
{
	library("data.table")
#filter for reads that are closer than 10kb to remove self ligations 
	df_int <- hicupinteraction[order(hicupinteraction$chr1, hicupinteraction$locus1, hicupinteraction$chr2, hicupinteraction$locus2), ]
	df_filtered <- df_int
	df_filtered$dist <-as.vector(ifelse(df_filtered$chr1 == df_filtered$chr2, abs(df_filtered$locus1 - df_filtered$locus2), Inf))
	df_filtered <-df_filtered[df_filtered$dist>filterdist,]
	df_filtered <-df_filtered[,1:4]
	
	rl_int <- RLinteraction[order(RLinteraction$chr1, RLinteraction$locus1, RLinteraction$chr2, RLinteraction$locus2), ]
	rl_filtered <- rl_int
	rl_filtered$dist <-as.vector(ifelse(rl_filtered$chr1 == rl_filtered$chr2, abs(rl_filtered$locus1 - rl_filtered$locus2), Inf))
	rl_filtered <-rl_filtered[rl_filtered$dist>filterdist,]
	rl_filtered <-rl_filtered[,1:4]
	
#bin interactions according to resolution
	if(res>1){
		binned_df_filtered <- binInteractions(df_filtered, res)
	}
	if(res==1){
		binned_df_filtered <- countDuplicates(df_filtered)
	}
	binned_df_filtered$int1 <-paste(binned_df_filtered$chr1,binned_df_filtered$locus1,sep='_')
	binned_df_filtered$int2 <-paste(binned_df_filtered$chr2,binned_df_filtered$locus2,sep='_')
	
	if(res>1){
		binned_rl_filtered <- binInteractions(rl_filtered, res)
	}
	if(res==1){
		binned_rl_filtered <- countDuplicates(rl_filtered)
	}
	binned_rl_filtered$int1 <-paste(binned_rl_filtered$chr1,binned_rl_filtered$locus1,sep='_')
	binned_rl_filtered$int2 <-paste(binned_rl_filtered$chr2,binned_rl_filtered$locus2,sep='_')
#diagonal removal
	if(removeDiagonal)
	{
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$int1!=binned_df_filtered$int2,]	
	
		binned_rl_filtered <- binned_rl_filtered[binned_rl_filtered$int1!=binned_rl_filtered$int2,]
	}
#	if(promoters!=NULL){
#		binned_df_filteredGR=GRangesList()
#		binned_df_filteredGR$locus1=GRanges(seqnames=binned_df_filtered$chr1,ranges=IRanges(start=binned_df_filtered$locus1,end=binned_df_filtered$locus1))
#		binned_df_filteredGR$locus2=GRanges(seqnames=binned_df_filtered$chr2,ranges=IRanges(start=binned_df_filtered$locus2,end=binned_df_filtered$locus2))
#		subs <- which(binned_df_filtered)
#	}
	
	if(cistrans=='cis'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1==binned_df_filtered$chr2,]	
		
		binned_rl_filtered <- binned_rl_filtered[binned_rl_filtered$chr1==binned_rl_filtered$chr2,]
	}
	if(cistrans=='trans'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1!=binned_df_filtered$chr2,]	
	
		binned_rl_filtered <- binned_rl_filtered[binned_rl_filtered$chr1!=binned_rl_filtered$chr2,]	
	}
	
#all read pairs used in binomial
	numberOfReadPairs <- sum(binned_df_filtered$frequencies)
	numberOfReadPairsRL <- sum(binned_rl_filtered$frequencies)
#calculate coverage 
	all_bins <- unique(c(binned_rl_filtered$int1,binned_rl_filtered$int2))
	all_bins <- sort(all_bins)
	
	binned_dt=data.table(binned_rl_filtered)
	covA <- binned_dt[,sum(c(frequencies)),by=int1]	
	covB <- binned_dt[,sum(c(frequencies)),by=int2]
	covA <- setkey(covA,key='int1')
	setnames(covB, 1, 'int1')
	covB <- setkey(covB,key='int1')
	cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
	cov$V1.x[is.na(cov$V1.x)]=0
	cov$V1.y[is.na(cov$V1.y)]=0
	cov$coverage=cov$V1.x+cov$V1.y
	coverage=cov$coverage
	names(coverage)=cov$int1
	
	sumcov <- sum(coverage)
	relative_coverage <- coverage/sumcov
	names(relative_coverage)=names(coverage)

#take only those interactions for which we have coverage information
	binned_df_filtered <- binned_df_filtered[binned_df_filtered$int1%in%all_bins&binned_df_filtered$int2%in%all_bins,]
	
	binned_df_filtered$cov1 <- relative_coverage[binned_df_filtered$int1]
	binned_df_filtered$cov2 <- relative_coverage[binned_df_filtered$int2]
#probability correction assuming on average equal probabilities for all interactions
	numberOfAllInteractions <- length(all_bins)^2
	upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
	
	if(cistrans!='all'){
		chromos <- unique(binned_df_filtered$chr1)
		chrlens <- c()
		for(cr in chromos){ 
			chrlens[cr] <- max(length(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),length(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
		}
		cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
		transBinNumber <- upperhalfBinNumber-cisBinNumber
	}
	
	diagonalProb <- sum(relative_coverage^2)
	if(cistrans=='all'){
		probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
	}
	if(cistrans=='cis'){
		probabilityCorrection <- upperhalfBinNumber/cisBinNumber
	}
	if(cistrans=='trans'){
		probabilityCorrection <- upperhalfBinNumber/transBinNumber
	}
	
	
	binned_df_filtered$probability <- binned_df_filtered$cov1*binned_df_filtered$cov2*2*probabilityCorrection
	
	binned_df_filtered$predicted <- binned_df_filtered$probability * numberOfReadPairs
	
	binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
									   {
									   binom.test(as.numeric(x[["frequencies"]]), numberOfReadPairs, as.numeric(x[["probability"]]), alternative = "greater")$p.value
									   }	
									   )
#observed over expected log ratio
	binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
#multiple testing correction separately for matrices with all interactions/only cis/only transs
	
	if(cistrans=='all'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber+length(all_bins))}
	}
	if(cistrans=='cis'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber+length(all_bins))}	
	}
	if(cistrans=='trans'){
		binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=transBinNumber)	
	}
	
	return(binned_df_filtered)
}

#############binomial test from aggregated reads############
binomialHiCagg=function(hicupinteraction, hindGR, removeDiagonal=TRUE, cistrans='all', parallel=FALSE, cores=8, dc=FALSE, dcstep=5000, sName=NULL)
{
	library("data.table")
	if(parallel)
	{
		library(parallel)
		print("running garbage collector before parallel fork")
		gc()
	}

	binned_df_filtered <-hicupinteraction
	binned_df_filtered$int1 <-paste(binned_df_filtered$chr1,binned_df_filtered$locus1,sep='_')
	binned_df_filtered$int2 <-paste(binned_df_filtered$chr2,binned_df_filtered$locus2,sep='_')
#diagonal removal
	if(removeDiagonal)
	{
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$int1!=binned_df_filtered$int2,]	
	}
	if(cistrans=='cis'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1==binned_df_filtered$chr2,]	
	}
	if(cistrans=='trans'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1!=binned_df_filtered$chr2,]	
	}
	
#all read pairs used in binomial
	numberOfReadPairs <- sum(binned_df_filtered$frequencies)
#calculate coverage 
	all_bins <- unique(c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2)))
	all_bins <- sort(all_bins)
	if(nrow(binned_df_filtered)>1e8)
	{
		t <- ceiling(nrow(binned_df_filtered)/1e8)
		dfList <- list()
		dfList[[1]] <- binned_df_filtered[1:1e8,]
		for(i in 2:t){
			dfList[[i]] <- binned_df_filtered[(((i-1)*1e8)+1):min((i*1e8),nrow(binned_df_filtered)),]
		}
		dtList <- lapply(dfList, data.table)
		covAs <- lapply(dtList, function(x) x[,sum(frequencies), by=int1])
		covBs <- lapply(dtList, function(x) x[,sum(frequencies), by=int2])
		covAm <- do.call(rbind, covAs)
		covBm <- do.call(rbind, covBs)
		covA <- covAm[,sum(V1),by=int1]
		covB <- covBm[,sum(V1),by=int2]	
	}else{
	binned_dt=data.table(binned_df_filtered)
	covA <- binned_dt[,sum(frequencies),by=int1]	
	covB <- binned_dt[,sum(frequencies),by=int2]
	}
	covA <- setkey(covA,key='int1')
	setnames(covB, 1,'int1')
	covB <- setkey(covB,key='int1')
#check consistency due to data.table bug after ~100 million entries
#uiques, and coverage values
	if(length(covA$int1) != length(unique(covA$int1))) stop(paste0("Error: coverage inconsistencies. State saved in ", "binomialState_", nrow(binned_df_filtered)))
	if(any(covA$V1 > 0 & covA$V1 < 1e-10)) stop(paste0("Error: coverage inconsistencies. State saved in ", "binomialState_", nrow(binned_df_filtered)))
	if(length(covB$int2) != length(unique(covB$int2))) stop(paste0("Error: coverage inconsistencies. State saved in ", "binomialState_", nrow(binned_df_filtered)))
	if(any(covB$V1 > 0 & covB$V1 < 1e-10)) stop(paste0("Error: coverage inconsistencies. State saved in ", "binomialState_", nrow(binned_df_filtered)))
	   
	cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
	cov$V1.x[is.na(cov$V1.x)]=0
	cov$V1.y[is.na(cov$V1.y)]=0
	cov$coverage=cov$V1.x+cov$V1.y
	coverage=cov$coverage
	names(coverage)=cov$int1
	sumcov <- sum(coverage)
	relative_coverage <- coverage/sumcov
	names(relative_coverage)=names(coverage)
	binned_df_filtered$coverage_source <- relative_coverage[binned_df_filtered$int1]
	binned_df_filtered$coverage_target <- relative_coverage[binned_df_filtered$int2]
#distance correction
	if(dc){
#put end of fragment
		binned_df_filtered$probabilityOfInteraction <- binned_df_filtered$coverage_source*binned_df_filtered$coverage_target*2
		sourceGR=GRanges(seqnames=binned_df_filtered$chr1, ranges=IRanges(start=binned_df_filtered$locus1, end=binned_df_filtered$locus1))
		targetGR=GRanges(seqnames=binned_df_filtered$chr2, ranges=IRanges(start=binned_df_filtered$locus2, end=binned_df_filtered$locus2))
		ss=findOverlaps(sourceGR, hindGR, type='any')
		st=findOverlaps(targetGR, hindGR, type='any')
		binned_df_filtered$end_source=end(ranges(hindGR[subjectHits(ss)]))
		binned_df_filtered$end_target=end(ranges(hindGR[subjectHits(st)]))
		print("end put")
#calculate distance
		distances=pmin(abs(binned_df_filtered$locus1-binned_df_filtered$locus2), abs(binned_df_filtered$locus1-binned_df_filtered$end_target), abs(binned_df_filtered$end_source-binned_df_filtered$locus2), abs(binned_df_filtered$end_source-binned_df_filtered$end_target))
		binned_df_filtered$distance=distances
		binned_df_filtered$distance[binned_df_filtered$chr1!=binned_df_filtered$chr2]=Inf
#make distance bins
		chromosomes=unique(binned_df_filtered$chr1)
		maxbin=max(binned_df_filtered$locus1[binned_df_filtered$chr1=='chrY'])
		steps=seq(0,maxbin,dcstep)
		print("distances calculated")
#calculate total readcount/distance bin to find out until where to do the distance correction
		totals=c()
		for(k in 1:(length(steps)-1)){
			totals[k]=sum(binned_df_filtered$frequencies[binned_df_filtered$distance>=steps[k] & binned_df_filtered$distance<steps[k+1]])
		}
#threshold
		threshold=max(totals)/1e6
		distTresh=steps[1:(max(which(totals>threshold))+1)]
#calculate the number of potential interactions per distance
		chromosomes=chromosomes[-20]
		chrlens <- c()
		for(cr in chromosomes){ 
			chrlens[cr] <- max(max(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),max(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
		}
	
	potentialInt=c()
		for(k in 1:(length(distTresh)-1)){
			potPerChr=c()
			for(l in 1:length(chrlens)){
				numberOfBinsPerChr=ceiling(chrlens[l]/distTresh[k+1])
				potPerChr[l]=ceiling((chrlens[l]-distTresh[k+1])/dcstep)
				potPerChr[l]=potPerChr[l]+1
			}
			potentialInt[k]=sum(potPerChr)
		}

#calculate the median readcount for a given distance
		means=c()
		for(k in 1:(length(distTresh)-1)){
			nonzero=length(binned_df_filtered$frequencies[binned_df_filtered$distance>=steps[k] & binned_df_filtered$distance<steps[k+1]])
			means[k]=mean(c(binned_df_filtered$frequencies[binned_df_filtered$distance>=steps[k] & binned_df_filtered$distance<steps[k+1]], rep(0, times=ifelse((potentialInt[k]-nonzero)>0,(potentialInt[k]-nonzero),0))))
#			means[k]=mean(binned_df_filtered$frequencies[binned_df_filtered$distance>=steps[k] & binned_df_filtered$distance<steps[k+1]])		
		}
		print("means calculated")
		ss=which(binned_df_filtered$distance!=Inf)
		ranges1=IRanges(start=binned_df_filtered$distance[ss], end=binned_df_filtered$distance[ss])
		ranges2=IRanges(start=distTresh[1:(length(distTresh)-1)], end=distTresh[2:length(distTresh)])
		ssr=findOverlaps(ranges1, ranges2, type='any')
		
		binned_df_filtered$distmean=1
		ssdr=ss[queryHits(ssr)]
		binned_df_filtered$distmean[ssdr]=means[subjectHits(ssr)]
		
		binned_df_filtered$correctedexp=binned_df_filtered$probabilityOfInteraction*binned_df_filtered$distmean
		binned_df_filtered$correctedexp[ssdr]=binned_df_filtered$correctedexp[ssdr]/sum(binned_df_filtered$distmean[ssdr])
		totalProb=sum(binned_df_filtered$probabilityOfInteraction[ssdr])
		totalNew=sum(binned_df_filtered$correctedexp[ssdr])
		correc=totalProb/totalNew
		binned_df_filtered$correctedexp[ssdr]=binned_df_filtered$correctedexp[ssdr]*correc
		binned_df_filtered$probabilityOfInteraction <-binned_df_filtered$correctedexp	
#saving state
		save(binned_df_filtered, relative_coverage, all_bins, numberOfReadPairs, removeDiagonal, cistrans, parallel, cores, file=paste0("binomState_", nrow(binned_df_filtered)))
		
	}
	

			
#probability correction assuming on average equal probabilities for all interactions
	numberOfAllInteractions <- length(all_bins)^2
	upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
	
	if(cistrans!='all'){
		chromos <- unique(binned_df_filtered$chr1)
		chrlens <- c()
		for(cr in chromos){ 
			chrlens[cr] <- max(length(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),length(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
		}
		cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
		transBinNumber <- upperhalfBinNumber-cisBinNumber
	}
	
	diagonalProb <- sum(relative_coverage^2)
	if(cistrans=='all'){
		probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
	}
	if(cistrans=='cis'){
		probabilityCorrection <- upperhalfBinNumber/cisBinNumber
	}
	if(cistrans=='trans'){
		probabilityCorrection <- upperhalfBinNumber/transBinNumber
	}
	
	if(dc){
		binned_df_filtered$probabilityOfInteraction <- binned_df_filtered$probabilityOfInteraction * probabilityCorrection
	}else
	{
		
	binned_df_filtered$probabilityOfInteraction <- binned_df_filtered$coverage_source*binned_df_filtered$coverage_target*2*probabilityCorrection
	}
	
	binned_df_filtered$predicted <- binned_df_filtered$probabilityOfInteraction * numberOfReadPairs
	
	print(summary(binned_df_filtered$probabilityOfInteraction))
	
	if(parallel)
	{
		library(parallel)
		if(nrow(binned_df_filtered)>1e8)
		{
			t <- ceiling(nrow(binned_df_filtered)/1e8)
			dfList <- list()
			dfList[[1]] <- binned_df_filtered[1:1e8,]
			for(i in 2:t){
				dfList[[i]] <- binned_df_filtered[(((i-1)*1e8)+1):min((i*1e8),nrow(binned_df_filtered)),]
			}
			dtList <- lapply(dfList, function(x) as.data.frame(t(cbind(as.numeric(x[["frequencies"]]), 
																	   as.numeric(x[["probabilityOfInteraction"]])))))
			pvalues=list()
			for(i in 1:length(dtList)){
				pvalues[[i]] <-unlist(mclapply(dtList[[i]], function(x)
													 {
													 binom.test(x[1], numberOfReadPairs, x[2], alternative = "greater")$p.value
													 }, 
													 mc.cores=cores))
			}
			save(pvalues, file=paste0('pvalues_', sName))
			pvals=unlist(pvalues)
			binned_df_filtered$pvalue <- pvals
		}else{
			
			binomParams <- as.data.frame(t(cbind(
												 as.numeric(binned_df_filtered[["frequencies"]]), 
												 as.numeric(binned_df_filtered[["probabilityOfInteraction"]]
															))))
			
# pvalue <- lapply(binomParams, function(x)
# 	{
# 	   binom.test(x[1], numberOfReadPairs, x[2], alternative = "greater")$p.value
# 	})
			binned_df_filtered$pvalue <- unlist(mclapply(binomParams, function(x)
														 {
														 binom.test(x[1], numberOfReadPairs, x[2], alternative = "greater")$p.value
														 }, 
														 mc.cores=cores))
			}
		}else{
		if(nrow(binned_df_filtered)>1e8)
		{
			t <- ceiling(nrow(binned_df_filtered)/1e8)
			dfList <- list()
			dfList[[1]] <- binned_df_filtered[1:1e8,]
			for(i in 2:t){
				dfList[[i]] <- binned_df_filtered[(((i-1)*1e8)+1):min((i*1e8),nrow(binned_df_filtered)),]
			}
			pvalues=list()
			for(i in 1:length(dfList)){
				pvalues[[i]] <-apply(dfList[[i]], 1, function(x)
									  {
									  binom.test(as.numeric(x[["frequencies"]]), numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
									  }	
									  )
			}
			pvals=unlist(pvalues)
			binned_df_filtered$pvalue <- pvals
		}else{
			
		binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
									   {
									   binom.test(as.numeric(x[["frequencies"]]), numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
									   }	
									   )
		}
	}
	
#observed over expected log ratio
	binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
#multiple testing correction separately for matrices with all interactions/only cis/only transs
	
	if(cistrans=='all'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber+length(all_bins))}
	}
	if(cistrans=='cis'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber+length(all_bins))}	
	}
	if(cistrans=='trans'){
		binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=transBinNumber)	
	}
	
	return(binned_df_filtered)
}

binomialHiCagg5C=function(hicupinteraction, cistrans='all')
{
	if(cistrans != 'all') stop("cistrans not implemented")

	library("data.table")
	binned_df_filtered <-hicupinteraction
	binned_df_filtered$int1 <-paste(binned_df_filtered$chr1,binned_df_filtered$locus1,sep='_')
	binned_df_filtered$int2 <-paste(binned_df_filtered$chr2,binned_df_filtered$locus2,sep='_')
	
	if(cistrans=='cis'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1==binned_df_filtered$chr2,]	
	}
	if(cistrans=='trans'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1!=binned_df_filtered$chr2,]	
	}
	
#all read pairs used in binomial
	numberOfReadPairs <- sum(binned_df_filtered$frequencies)
#calculate coverage 
	all_bins <- unique(c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2)))
	all_bins <- sort(all_bins)
	binned_dt=data.table(binned_df_filtered)
	covA <- binned_dt[,sum(c(frequencies)),by=int1]
	setnames(covA, 1:2, c("id", "coverage"))
	covA <- setkey(covA,key='id')

	covB <- binned_dt[,sum(c(frequencies)),by=int2]
	setnames(covB, 1:2, c("id", "coverage"))
	covB <- setkey(covB,key='id')

	covA$relative_coverage <- covA$coverage / numberOfReadPairs
	covB$relative_coverage <- covB$coverage / numberOfReadPairs

	binned_df_filtered$coverage_source <- covA[binned_df_filtered$int1]$relative_coverage
	binned_df_filtered$coverage_target <- covB[binned_df_filtered$int2]$relative_coverage

	#probability correction assuming on average equal probabilities for all interactions
	numberOfAllInteractions <- nrow(covA) * nrow(covB)
	
	# if(cistrans!='all'){
	# 	chromos <- unique(binned_df_filtered$chr1)
	# 	chrlens <- c()
	# 	for(cr in chromos){ 
	# 		chrlens[cr] <- max(length(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),length(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
	# 	}
	# 	cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
	# 	transBinNumber <- upperhalfBinNumber-cisBinNumber
	# }
	
	# if(cistrans=='cis'){
	# 	probabilityCorrection <- upperhalfBinNumber/cisBinNumber
	# }
	# if(cistrans=='trans'){
	# 	probabilityCorrection <- upperhalfBinNumber/transBinNumber
	# }
	
	
	binned_df_filtered$probabilityOfInteraction <- binned_df_filtered$coverage_source*binned_df_filtered$coverage_target
	
	binned_df_filtered$predicted <- binned_df_filtered$probabilityOfInteraction * numberOfReadPairs
	
	binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
									   {
									   binom.test(as.numeric(x[["frequencies"]]), numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
									   }	
									   )
#observed over expected log ratio
	binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
#multiple testing correction separately for matrices with all interactions/only cis/only transs
	
	if(cistrans=='all'){
		binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=numberOfAllInteractions)
	}
	# if(cistrans=='cis'){
	# 	binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber+length(all_bins))}	
	# }
	# if(cistrans=='trans'){
	# 	binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=transBinNumber)	
	# }
	
	return(binned_df_filtered)
}

test_multicore <- function()
{
	# small <- get(load("../hicuptest/interactions_aggregated_pooledHicuptest_1frag"))
	small <- get(load("../hicuptest/interactions_binomial_pooledHicuptest_1frag"))

	sampleSizes <- c(1e3)
	cores <- c(1)

	sampleSizes <- c(1e5, 1e6, 1e7)
	cores <- c(2, 4, 8)

	binomialTimes <- list()
	for(i in sampleSizes)
	{
		interactions <- small[sample(nrow(small), i, replace=TRUE), ]
		print(paste("standard binomial with size", i))
		biT <- system.time(binomialHiCagg(interactions))
		binomialTimes[[paste("binomial", i, 1, sep="_")]] <- c(i, 1, biT[['elapsed']])
		print(biT)
		for(j in cores)
		{
			print(paste("parallel binomial with size", i, "cores:", j))
			biT <- system.time(binomialHiCagg(interactions, parallel=TRUE, cores=j))
			binomialTimes[[paste("binomial", i, j, sep="_")]] <- c(i, j, biT[['elapsed']])
			print(biT)
		}
		save(binomialTimes, file=paste("binomialTimes", i, sep="_"))

	}

	df <- as.data.frame(t(data.frame(binomialTimes)))
	names(df) <- c("samples", "cores", "seconds")

	save(binomialTimes, df, file="parallelBinomial")

	library(ggplot2)
	pdffig("parallelBinomial")
	df$cores <- as.character(df$cores)
	ggplot(data=df, aes(x=samples, y=seconds, group=(cores), colour=cores, fill=cores)) +
		scale_x_log10() + 
	 	# geom_line() +
		geom_bar(position="dodge", stat="identity")
		# scale_y_log10()

	dev.off()
}

coverageForbinomialHiCagg <- function(hicupinteraction, removeDiagonal=TRUE, cistrans='all')
{
	library("data.table")

	binned_df_filtered <-hicupinteraction
	binned_df_filtered$int1 <-paste(binned_df_filtered$chr1,binned_df_filtered$locus1,sep='_')
	binned_df_filtered$int2 <-paste(binned_df_filtered$chr2,binned_df_filtered$locus2,sep='_')
#diagonal removal
	if(removeDiagonal)
	{
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$int1!=binned_df_filtered$int2,]	
	}
	
	if(cistrans=='cis'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1==binned_df_filtered$chr2,]	
	}
	if(cistrans=='trans'){
		binned_df_filtered <- binned_df_filtered[binned_df_filtered$chr1!=binned_df_filtered$chr2,]	
	}
	
#all read pairs used in binomial
	numberOfReadPairs <- sum(binned_df_filtered$frequencies)
#calculate coverage 
	# all_bins <- unique(c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2)))
	all_bins <- unique(c(unique(binned_df_filtered$int1), unique(binned_df_filtered$int2)))
	all_bins <- sort(all_bins)
	binned_dt=data.table(binned_df_filtered)
	covA <- binned_dt[,sum(c(frequencies)),by=int1]	
	covB <- binned_dt[,sum(c(frequencies)),by=int2]
	covA <- setkey(covA,key='int1')
	setnames(covB, 1,'int1')
	covB <- setkey(covB,key='int1')
	cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
	cov$V1.x[is.na(cov$V1.x)]=0
	cov$V1.y[is.na(cov$V1.y)]=0
	cov$coverage=cov$V1.x+cov$V1.y
	coverage=cov$coverage
	names(coverage)=cov$int1
	
	sumcov <- sum(coverage)
	relative_coverage <- coverage/sumcov
	names(relative_coverage)=names(coverage)
	binned_df_filtered$coverage_source <- relative_coverage[binned_df_filtered$int1]
	binned_df_filtered$coverage_target <- relative_coverage[binned_df_filtered$int2]
#probability correction assuming on average equal probabilities for all interactions
	numberOfAllInteractions <- length(all_bins)^2
	upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
	
	if(cistrans!='all'){
		chromos <- unique(binned_df_filtered$chr1)
		chrlens <- c()
		for(cr in chromos){ 
			chrlens[cr] <- max(length(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),length(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
		}
		cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
		transBinNumber <- upperhalfBinNumber-cisBinNumber
	}
	
	diagonalProb <- sum(relative_coverage^2)
	if(cistrans=='all'){
		probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
	}
	if(cistrans=='cis'){
		probabilityCorrection <- upperhalfBinNumber/cisBinNumber
	}
	if(cistrans=='trans'){
		probabilityCorrection <- upperhalfBinNumber/transBinNumber
	}
	
	
	binned_df_filtered$probabilityOfInteraction <- binned_df_filtered$coverage_source*binned_df_filtered$coverage_target*2*probabilityCorrection
	
	binned_df_filtered$predicted <- binned_df_filtered$probabilityOfInteraction * numberOfReadPairs

	binned_df_filtered

	# save(binned_df_filtered, file="binned_df_filtered_final_g")
}
