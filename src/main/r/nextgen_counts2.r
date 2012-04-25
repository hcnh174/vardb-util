getCodonCountFilename <- function(config, sample, type)
{
	return(concat(config@counts.dir,'/',sample,'.',type,'.txt'))
}
#getCodonCountFilename(config,'KT9.random__HCV-KT9','codons')
#getRefForSample('HCV-KT9.random__HCV-KT9')

#uses GATK custom walker to count all the bases at each position
countCodonsForId <- function(config, id, bam.dir=config@bam.dir, filter=config@filter)
{
	region <- getField(config@data,id,'region') #region <- config@data[id,'region']
	ref <- getField(config@data,id,'ref') #ref <- config@data[id,'ref']
	sample <- concat(id,'__',ref)
	reffile <- getRefFile(config,ref)
	bamfile <- ifelse(filter, concat(bam.dir,'/',sample,'.filtered.bam'), concat(bam.dir,'/',sample,'.bam'))
	ntcountsfile <- getCodonCountFilename(config,id,'nt')
	codoncountsfile <- getCodonCountFilename(config,id,'codons')
	aacountsfile <- getCodonCountFilename(config,id,'aa')
	checkFileExists(reffile)
	checkFileExists(bamfile)
	
	if (config@force | !file.exists(ntcountsfile) | !file.exists(codoncountsfile) | !file.exists(aacountsfile))
	{
		str <- 'java -Xmx8g'
		str <- concat(str,' -cp $VARDB_UTIL_HOME/target/gatk-walkers.jar:$GATK_HOME/GenomeAnalysisTK.jar')
		str <- concat(str,' org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants')
		str <- concat(str,' -dt NONE')
		str <- concat(str,' -et NO_ET')
		str <- concat(str,' -R ',reffile)
		str <- concat(str,' -I ',bamfile)
		str <- concat(str,' --validation_strictness strict')
		str <- concat(str,' --ntcounts ',ntcountsfile)
		str <- concat(str,' --codoncounts ',codoncountsfile)
		str <- concat(str,' --aacounts ',aacountsfile)
		str <- concat(str,' -L ',getIntervalForRegion(config,ref,region))
		runCommand(str)
	}
	else printcat('skipping CountVariants command because count files already exist for sample: ',id)
	checkFileExists(ntcountsfile)
	checkFileExists(codoncountsfile)
	checkFileExists(aacountsfile)
}
#countCodonsForId(config,'nextgen1-1E',filter=TRUE)
#countCodons(config)


getCodonCountSubsetForSample <- function(config, sample, region, filetype)
{
	printcat('sample=',sample,', region=',region)
	row <- config@data[which(config@data$sample==sample & config@data$region==region),]
	if (nrow(row)==0)
		return(data.frame())
	if (nrow(row)>1)
		throw2('more than one row returned for sample=',sample,' region=',region,': ',row)
	rowname <- row$id
	printcat('rowname=',rowname)
	filename <- getCodonCountFilename(config,rowname,filetype)
	printcat('filename=',filename)
	checkFileExists(filename)
	data <- loadDataFrame(filename)
	data <- renameColumn(data,'position','ntnum')
#	if (filetype=='nt')
#	{
#		library(reshape)
#		data <- excludeColumns(data,'depth')
#		data.melted <- melt(data, id=c('ntnum','refnt'))
#		colnames(data.melted) <- splitFields('ntnum,refnt,nt,count')
#		data.melted$nt <- toupper(data.melted$nt)
#		data.melted <- data.melted[order(data.melted$count, decreasing=TRUE),]
#		data.melted <- data.melted[order(data.melted$ntnum),]
#		data <- data.melted
#	}
	loc <- getLocationForRegion(config,row$ref,region)
	data$group <- row$group
	data$column <- row$column
	data$ref <- row$ref
	data$sample <- sample
	data$rowname <- rowname
	data$region <- region
	data$aanum <- sapply(data$ntnum,function(ntnum)
	{
		aanum <- (ntnum-loc$start)/3 + 1
		return(aanum)
	})
	return(data)
}
#head(getCodonCountSubsetForSample(config,'KT9.plasmid__HCV-KT9','NS3aa36','nt'))

#getCodonCountSubset <- function(config, samples, region, filetype, start, end=start, minreads=0)
#{
#	data <- NULL
#	for (sample in samples)
#	{
#		data.sample <- getCodonCountSubsetForSample(config,sample,region,filetype)
#		if (nrow(data.sample)==0)
#			next
#		if (is.null(data))
#			data <- data.sample
#		else data <- rbind(data,data.sample)
#	}
#	data.subset <- data
#	data.subset <- data[which(data$aanum>=start & data$aanum<=end & data$count>=minreads),]
#	#data.subset <- data.subset[which(data.subset$region %in% splitFields(region)),]
#	data.subset$column <- factor(data.subset$column)
#	data.subset$aanum <- factor(data.subset$aanum)
#	return(data.subset)
#}
##getCodonCountSubset(config,getSamplesForSubGroup(config,'hcv_infection','hcv_infection'), 'NS5Aaa31', 'aa', 31)
#

getCodonCountSubset <- function(config, samples, region, filetype, start, end=start, minreads=0)
{
	data <- NULL
	for (sample in samples)
	{
		try({
			data.sample <- getCodonCountSubsetForSample(config,sample,region,filetype)
			if (nrow(data.sample)==0)
				next
			if (is.null(data))
				data <- data.sample
			else data <- rbind(data,data.sample)
		})
	}
	if (is.null(data))
		throw2('No data for in region ',region, ' for samples: ',samples)
	data.subset <- data[which(data$aanum>=start & data$aanum<=end & data$count>=minreads),]	
	if (nrow(data.subset)==0)
		throw2('No data for in region ',region, ' for samples: ',samples)
	
	data.subset$column <- factor(data.subset$column)
	data.subset$aanum <- factor(data.subset$aanum)
	
	#indicate whether the aa is the same as the ref or not
	ref <- getRefForSamples(config,samples)
	data.subset$isrefaa <- apply(data.subset, 1, function(row)
	{ 
		try({
			refaa <- getReferenceAminoAcid(config, ref, region, row['aanum'])
			return(row['aa']==refaa)
		})
	})
	data.subset$isrefcodon <- apply(data.subset, 1, function(row)
	{ 
		refcodon <- getReferenceCodon(config, ref, region, row['aanum'])
		return(row['codon']==refcodon)
	})

	data.subset$subtype <- ifelse(!data.subset$isrefaa, 'non-synonymous', ifelse(!data.subset$isrefcodon, 'synonymous', 'reference'))
	
	return(data.subset)
}
#data.subset <- getCodonCountSubset(config,getSamplesForSubject(config,'PXB0220-0030'), 'NS3aa36', 'codons', 36)
