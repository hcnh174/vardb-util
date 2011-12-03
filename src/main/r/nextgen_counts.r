loadPileupData <- function(config, sample)
{
	samplename <- ifelse(config@filter,concat(sample,'.filtered'), sample)
	filename <- concat(config@pileup.dir,'/',samplename,'.txt'); print(filename) # load the corresponding data file
	data <- loadDataFrame(filename)
	printcat('loaded file ',filename,'. contains ',nrow(data),' reads')
	ref <- getRefForSample(sample)
	startnt <- getField(config@refs,ref,'start')
	data$ntnum <- data$position + startnt
	return(data)
}
#data <- loadPileupData(config,'PXB0220-0002.wk08__HCV-KT9')
#
#extractCodonData <- function(data, ntnum, drop.ambig=FALSE)
#{
#	#ntnum <- ntnum - 1 #hack
#	nt1 <- data[which(data$ntnum==ntnum),'nt']
#	nt2 <- data[which(data$ntnum==ntnum+1),'nt']
#	nt3 <- data[which(data$ntnum==ntnum+2),'nt']
#	printcat('nt1=',length(nt1),' nt2=',length(nt2),' nt3=',length(nt3))
#	codons <- paste(nt1, nt2, nt3, sep='')
#	#if (drop.ambig)
#	#	codons <- removeAmbiguousCodons(codons)
#	return(codons)
#}
##extractCodonData(data,6348)

extractCodonData <- function(data, ntnum)
{
	cols <- c('read','nt')
	nts1 <- data[which(data$ntnum==ntnum),cols]; nts1 <- nts1[order(nts1$read),]
	nts2 <- data[which(data$ntnum==ntnum+1),cols]; nts2 <- nts2[order(nts2$read),]
	nts3 <- data[which(data$ntnum==ntnum+2),cols]; nts3 <- nts3[order(nts3$read),]
	
	reads <- nts1$read
	reads <- reads[which(reads %in% nts2$read)]
	reads <- reads[which(reads %in% nts3$read)]
	reads <- sort(reads)
	
	nt1 <- nts1[which(nts1$read %in% reads),'nt']
	nt2 <- nts2[which(nts2$read %in% reads),'nt']
	nt3 <- nts3[which(nts3$read %in% reads),'nt']
	
	codons <- paste(nt1, nt2, nt3, sep='')
	codons <- removeAmbiguousCodons(codons)
	return(codons)
}
#extractCodonData(data,6348)

#################################################################

getNtCounts <- function(data, ntnum)
{
	nts <- data[which(data$ntnum==ntnum & data$nt!='N'),'nt']
	freqs <- sort(xtabs(as.data.frame(nts)), decreasing=TRUE)
	total <- sum(freqs)		
	rank <- 1
	counts <- data.frame()
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
#counts <- getNtCounts(data,6348)

getCodonCounts <- function(data, ntnum)
{
	counts <- data.frame()
	codons <- extractCodonData(data,ntnum)
	if (length(codons)==0)
		return(counts)
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
	codons <- extractCodonData(data,ntnum)
	if (length(codons)==0)
		return(counts)
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

createNtCountTable <- function(config, data, params, positions)
{
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

createCodonCountTable <- function(config, data, params, positions)
{
	counts <- data.frame()
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

createAminoAcidCountTable <- function(config, data, params, positions)
{
	counts <- data.frame()	
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
		counts$group <- as.character(params@group)
		counts$subject <- as.character(params@subject)
		counts$sample <- as.character(params@sample)
		counts$column <- as.character(params@column)
		counts$ref <- as.character(params@ref)
		counts$region <- as.character(params@region)
	}, silent=FALSE)
	cols <- names(counts)
	cols <- c('group','subject','sample','column','ref','region',cols[1:(length(cols)-6)])
	counts <- counts[,cols]
	return(counts)
}

countCodonsForRegion <- function(config, data, params, variantdata)
{
	printcat('countCodonsForRegion: ',as.character(params@region))
	positions <- getCodonPositionsForRegion(config,params@region,params@ref)
	variantdata@nt <- rbind(variantdata@nt, createNtCountTable(config, data, params, positions))
	variantdata@codons <- rbind(variantdata@codons, createCodonCountTable(config, data, params, positions))
	variantdata@aa <- rbind(variantdata@aa, createAminoAcidCountTable(config, data, params, positions))
	return(variantdata)
}
#params@region <- 'NS3aa156'
#variantdata <- countCodonsForRegion(config, data, params, variantdata)

getCodonCountFilename <- function(config, sample, type)
{
	return(concat(config@counts.dir,'/',sample,'.',type,'.txt'))
}
#getCodonCountFilename(config,'KT9.random__HCV-KT9','codons')

countCodonsForSample <- function(config, sample)
{
	printcat('countCodonsForSample: ',sample)
	params <- new('sampleparams', sample=sample)
	row <-config@data[which(config@data$sample==sample),]
	params@group <- unique(row$group)
	params@subject <- unique(row$subject)
	params@column <- unique(row$column)
	params@ref <- getRefForSample(sample)
	data <- loadPileupData(config,sample)
	variantdata <- new('variantdata')	
	for (region in getRegionsForSample(config,sample))
	{
		params@region <- region
		try(variantdata <- countCodonsForRegion(config, data, params, variantdata))
	}
	
	writeTable(variantdata@nt, getCodonCountFilename(config,sample,'nt'), row.names=FALSE)
	writeTable(variantdata@codons, getCodonCountFilename(config,sample,'codons'), row.names=FALSE)
	writeTable(variantdata@aa, getCodonCountFilename(config,sample,'aa'), row.names=FALSE)
	return(variantdata)
}
#countCodonsForSample(config,'PXB0220-0002.wk08__HCV-KT9_PXB0220-0002')#KT9.random__HCV-KT9

#uses GATK custom walker to count all the bases at each position
countBasesForSample <- function(config, sample, bam.dir=config@bam.dir, out.dir=config@basecounts.dir)
{
	reffile <- getRefFile(config,getRefForSample(sample))
	bamfile <- concat(bam.dir,'/',sample,'.bam')
	outfile <- concat(out.dir,'/',sample,'.txt')
	checkFileExists(reffile)
	checkFileExists(bamfile)	
	str <- 'java -Xmx8g'
	str <- concat(str,' -cp $VARDB_UTIL_HOME/target/gatk-walkers.jar:$GATK_HOME/GenomeAnalysisTK.jar')
	str <- concat(str,' org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants')
	str <- concat(str,' -dt NONE')
	str <- concat(str,' -R ',reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -o ',outfile)
	runCommand(str)
	checkFileExists(outfile)
}

countBases <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		print(sample)
		try(countBasesForSample(config,sample))
	}
}
#countBases(config)