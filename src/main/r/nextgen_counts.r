loadPileupData <- function(config, sample)
{
	samplename <- ifelse(config@filter,concat(sample,'.filtered'), sample)
	filename <- concat(config@pileup.dir,'/',samplename,'.txt'); print(filename) # load the corresponding data file
	data <- loadDataFrame(filename)
	print(concat('loaded file ',filename,'. contains ',nrow(data),' reads'))
	ref <- get_ref_for_sample(sample)
	startnt <- config@refs[ref,'startnt']
	data$ntnum <- data$position + startnt
	return(data)
}
#data <- loadPileupData(config,'110719-4.p29.NS3aa36@NS3-36')

#getCodonPositionsForRegion <- function(config, region)
#{
#	ref <- config@regions[region,'ref']
#	startnt <- config@refs[ref,'startnt']
#	startntrel <- config@refs[ref,'startntrel']
#	sequence <- config@refs[ref,'sequence']
#	endntrel <- startntrel + nchar(sequence)
#	#print(sequence)
#	#print(translateCodon(sequence))
#	positions <- data.frame()
#	print(startntrel)
#	offset <- (startntrel-1) %% 3
#	print(offset)
#	if (offset!=0)
#	{
#		print(concat('adjusting to the start of a codon: offset=',offset))
#		startntrel <- startntrel + 3 - offset
#		startnt <- startnt + 3 - offset
#		if ((startntrel-1) %% 3!=0)
#			stop(concat('start number is not a multiple of 3: ',startntrel,' in region ',region))
#	}
#	codon <- (startntrel-1)/3 + 1
#	print(concat('ref: ',ref,' startnt: ',startnt,', startntrel: ',startntrel,', codon: ',codon))
#	
#	for (position in seq(startntrel,endntrel,3))
#	{
#		index <- position-startntrel+1
#		if (index+2 > nchar(sequence))
#			break
#		refcodon <- extractSequence(sequence, index, index+2)
#		refaa <- translateCodon(refcodon)
#		positions <- rbind(positions, data.frame(codon=codon, ntnum=startnt+position, relpos=position, refcodon=refcodon, refaa=refaa))
#		codon <- codon + 1
#	}
#	startaa <- config@regions[region,'startaa']
#	endaa <- config@regions[region,'endaa']
#	#if (!is.na(startaa))
#	#	positions <- positions[which(positions$codon >= startaa),]
#	#if (!is.na(endaa))
#	#	positions <- positions[which(positions$codon <= endaa),]
#	print(concat('sequence: ',seqinr::c2s(positions$refaa)))
#	return(positions)
#}
#

#################################################################

getNtCounts <- function(data, ntnum)
{
	counts <- data.frame()
	nts <- data[which(data$ntnum==ntnum),'nt']
	freqs <- sort(xtabs(as.data.frame(nts)), decreasing=TRUE)
	total <- sum(freqs)				
	rank <- 1
	for (nt in names(freqs))
	{
		count <- freqs[nt]
		freq <- count/total
		row <- data.frame(ntnum=ntnum, aanum=NA, nt=nt, rank=rank, count=count, freq=freq) 
		counts <- rbind(counts,row)
		rank <- rank+1
	}
	return(counts)
}
#counts <- getNtCounts(data,3885)

getCodonCounts <- function(data, ntnum)
{
	counts <- data.frame()
	codons <- extractCodonData(data,ntnum, FALSE)
	if (length(codons)==0)
		return(counts)
	#print(codons)
	#aanum <- positions[which(positions$ntnum==ntnum),'codon']
	#print(paste('ntnum=',ntnum,'numcodons=',length(codons)))
	freqs <- sort(xtabs(as.data.frame(codons)), decreasing=TRUE)
	total <- sum(freqs)		
	rank <- 1
	for (codon in names(freqs))
	{
		count <- freqs[codon]
		freq <- count/total
		row <- data.frame(ntnum=ntnum, aanum=NA, codon=codon, rank=rank, count=count, freq=freq) #sample=sample, 
		counts <- rbind(counts,row)
		rank <- rank +1
	}
	counts$aa <- sapply(counts$codon,translateCodon)
	return(counts)
}
#counts <- getCodonCounts(data,3885)

getAaCounts <- function(data, ntnum)
{
	counts <- data.frame()
	codons <- extractCodonData(data,ntnum,FALSE)
	if (length(codons)==0)
		return(counts)
	#aanum <- positions[which(positions$ntnum==ntnum),'codon']
	aa <- sapply(codons,translateCodon)
	freqs <- sort(xtabs(as.data.frame(aa)), decreasing=TRUE)
	total <- sum(freqs)
	rank <- 1
	for (aa in names(freqs))
	{
		count <- freqs[aa]
		freq <- count/total
		row <- data.frame(ntnum=ntnum, aanum=NA, aa=aa, rank=rank, count=count, freq=freq)
		counts <- rbind(counts,row)
		rank <- rank +1
	}
	return(counts)
}
#counts <- getAaCounts(data,3885)

########################################################################

createNtCountTable <- function(config, data, params)
{
	positions <- getCodonPositionsForRegion(config,params@region)
	counts <- data.frame()
	for (ntnum in positions$ntnum)
	{
		for (offset in 0:2)
		{
			ntcounts <- getNtCounts(data,ntnum + offset)			
			if (nrow(ntcounts)>0)
			{
				ntcounts$aanum <- positions[which(positions$ntnum==ntnum),'codon']
				counts <- rbind(counts,ntcounts)		
			}
		}
	}
	counts$nt <- as.character(counts$nt)
	counts$ntnum <- factor(counts$ntnum)
	counts$aanum <- factor(counts$aanum)
	counts <- .appendSampleParams(counts, params)
	return(counts)
}
#params <- new('sampleparams',subject='10201689', sample='110719-4.p29.NS3aa36@NS3-36', region='NS3-36@NS3aa36', replicate=1)
#data <- loadPileupData(config,params@sample)
#counts <- createNtCountTable(config,data,params)

createCodonCountTable <- function(config, data, params)
{
	counts <- data.frame()
	positions <- getCodonPositionsForRegion(config,params@region)
	for (ntnum in positions$ntnum)
	{
		codoncounts <- getCodonCounts(data,ntnum)
		if (nrow(codoncounts)>0)
		{
			codoncounts$aanum <- positions[which(positions$ntnum==ntnum),'codon']
			counts <- rbind(counts,codoncounts)
		}
	}
	counts$codon <- as.character(counts$codon)
	counts$ntnum <- factor(counts$ntnum)
	counts$aanum <- factor(counts$aanum)
	counts <- .appendSampleParams(counts, params)
	return(counts)
}
#counts <- createCodonCountTable(config,data,params)

createAminoAcidCountTable <- function(config, data, params)
{
	counts <- data.frame()
	positions <- getCodonPositionsForRegion(config,params@region)
	for (ntnum in positions$ntnum)
	{	
		aacounts <- getAaCounts(data,ntnum)
		if (nrow(aacounts)>0)
		{
			aacounts$aanum <- positions[which(positions$ntnum==ntnum),'codon']
			counts <- rbind(counts,aacounts)
		}
	}
	counts$ntnum <- factor(counts$ntnum)
	counts$aanum <- factor(counts$aanum)
	counts <- .appendSampleParams(counts, params)
	return(counts)
}
#counts <- createAminoAcidCountTable(config,data,params@region)

.appendSampleParams <- function(counts, params)
{
	try({
		counts$subject <- as.character(params@subject)
		counts$sample <- as.character(params@sample)
		counts$replicate <- as.character(params@replicate)
		counts$region <- as.character(params@region)
	}, silent=FALSE)
	cols <- names(counts)
	cols <- c('subject','sample','replicate','region',cols[1:(length(cols)-4)])
	counts <- counts[,cols]
	return(counts)
}

count_codons_for_region <- function(config, data, params, variantdata)
{
	print(concat('count_codons_for_region: ',as.character(params@region)))
	variantdata@nt <- rbind(variantdata@nt, createNtCountTable(config, data, params))
	variantdata@codons <- rbind(variantdata@codons, createCodonCountTable(config, data, params))
	variantdata@aa <- rbind(variantdata@aa, createAminoAcidCountTable(config, data, params))
	return(variantdata)
}
#params@region <- 'NS3-36@NS3aa36'
#variantdata <- count_codons_for_region(config, data, params, variantdata)

count_codons_for_sample <- function(config, params, variantdata)
{
	print(concat('count_codons_for_sample: ',params@sample))
	data <- loadPileupData(config,params@sample)
	for (region in get_regions_for_sample(config,params@sample))
	{
		params@region <- region
		variantdata <- count_codons_for_region(config, data, params, variantdata)
	}
	return(variantdata)
}
#count_codons_for_sample(config,params,variantdata)

#params <- new('sampleparams', subject='8538159', sample='110719-2.p39.NS3aa156@NS3-156-R', region='NS3-156-R@NS3aa156', replicate=1)
count_codons_for_subject <- function(config, params)
{
	if (is.character(params))
		params <- new('sampleparams', subject=params)
	print(concat('count_codons_for_subject: ',as.character(params@subject)))
	variantdata <- new('variantdata')
	samples <- config@runs[which(config@runs$subject==params@subject),]
	if (nrow(samples)==0)
		stop(concat('cannot find any samples for subject ',params@subject))
	for (sample in samples[,'sample'])
	{
		params@sample <- sample
		params@replicate <- samples[sample,'replicate']
		try({variantdata <- count_codons_for_sample(config, params, variantdata)}, silent=FALSE)
	}
	counts.dir <- concat(config@counts.dir,'/')
	writeTable(variantdata@nt, concat(counts.dir,params@subject,'.nt.txt'), row.names=FALSE)
	writeTable(variantdata@codons, concat(counts.dir,params@subject,'.codons.txt'), row.names=FALSE)
	writeTable(variantdata@aa, concat(counts.dir,params@subject,'.aa.txt'), row.names=FALSE)
	return(variantdata)
}
#counts <- count_codons_for_subject(config,'8538159')
#counts <- count_codons_for_subject(config,'10348001')
#counts <- count_codons_for_subject(config, '10201689')

