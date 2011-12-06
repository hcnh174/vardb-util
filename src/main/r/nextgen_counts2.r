getCodonCountFilename <- function(config, sample, type)
{
	return(concat(config@counts.dir,'/',sample,'.',type,'.txt'))
}
#getCodonCountFilename(config,'KT9.random__HCV-KT9','codons')
#getRefForSample('HCV-KT9.random__HCV-KT9')

#uses GATK custom walker to count all the bases at each position
countBasesForSample <- function(config, id, bam.dir=config@bam.dir)
{
	region <- config@data[id,'region']
	ref <- config@data[id,'ref']
	sample <- concat(id,'__',ref)
	reffile <- getRefFile(config,ref)
	bamfile <- concat(bam.dir,'/',sample,'.bam')
	ntcountsfile <- getCodonCountFilename(config,id,'nt')
	codoncountsfile <- getCodonCountFilename(config,id,'codons')
	aacountsfile <- getCodonCountFilename(config,id,'aa')
	checkFileExists(reffile)
	checkFileExists(bamfile)
	
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
	checkFileExists(ntcountsfile)
	checkFileExists(codoncountsfile)
	checkFileExists(aacountsfile)
}
#countBasesForSample(config,'nextgen1-7E')

countBases <- function(config, ids=rownames(config@data))
{
	for (id in ids)
	{
		try(countBasesForSample(config,id))
	}
}
#countBases(config)


getCodonCountSubsetForSample <- function(config, sample, region, filetype)
{
	printcat('sample=',sample,', region=',region)
	row <- config@data[which(config@data$sample==sample & config@data$region==region),]
	if (nrow(row)==0)
		return(data.frame())
	if (nrow(row)>1)
		throw('more than one row returned for sample=',sample,' region=',region,': ',row)
	rowname <- row$id
	printcat('rowname=',rowname)
	filename <- getCodonCountFilename(config,rowname,filetype)
	printcat('filename=',filename)
	checkFileExists(filename)
	data <- loadDataFrame(filename)
	loc <- getLocationForRegion(config,row$ref,region)
	data$ntnum <- data$position #hack - change this column name in walker output files
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

getCodonCountSubset <- function(config, samples, region, filetype, start, end=start, minreads=0)
{
	data <- NULL
	for (sample in samples)
	{
		data.sample <- getCodonCountSubsetForSample(config,sample,region,filetype)
		if (nrow(data.sample)==0)
			next
		if (is.null(data))
			data <- data.sample
		else data <- rbind(data,data.sample)
	}
	data.subset <- data
	data.subset <- data[which(data$aanum>=start & data$aanum<=end & data$count>=minreads),]
	#data.subset <- data.subset[which(data.subset$region %in% splitFields(region)),]
	data.subset$column <- factor(data.subset$column)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,getSamplesForSubGroup(config,'hcv_infection','hcv_infection'), 'NS5Aaa31', 'aa', 31)
