source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'common.r',sep=''))
loadUtilFiles('nextgen')

setCurDir('nextgen2')
config <- new('nextgenconfig')


#setwd('n:/vcf/')
#setwd('f:/analysis/nextgen2/vcf/')
setwd('C:/Documents and Settings/nhayes/My Documents/My Dropbox/vcf')
print(getwd())


isFieldIncluded <- function(fields,field)
{
	if (is.null(fields))
		return(TRUE)
	fields <- splitFields(fields)
	return(containsElement(fields,field))
}
#isFieldIncluded(NULL,'GT') #TRUE
#isFieldIncluded('GT,AC','GT') #TRUE
#isFieldIncluded('GT,AC','AB') #FALSE

getNumberOfLinesToSkip <- function(filename)
{
	lines <- readLines(filename)
	for (index in 1:length(lines))
	{
		line <- lines[index]
		if (substr(line,1,2)!='##')
			return(index-1)
	}
	throw(concat('cannot find end of metadata header in file ',filename))
}
#getNumberOfLinesToSkip('merged.filtered.vcf')

getInfoFields <- function(data)
{
	rowname <- rownames(data)[1] # get the first row
	pairs <- strsplit(data[rowname,'INFO'],';')[[1]]
	fields <- c()
	for (pair in pairs)
	{
		fields <- c(fields,strsplit(pair,'=')[[1]][1])
	}
	return(fields)
}

extractFieldData <- function(fields,values,field)
{
	fields <- splitFields(fields,':')
	values <- splitFields(values,':')
	for (i in 1:length(fields))
	{
		if (fields[i]==field)
			return(values[i])
	}
}
#extractFieldData('GT:AB:AD:DP:FA:GQ:MQ0:PL', '0/0:.:68,0:68:0.000:99:0:0,156,1446', 'AD')

loadVcf <- function(filename)
{
	numskip <- getNumberOfLinesToSkip(filename)	
	data <- read.table(filename, header=TRUE, encoding='UTF-8', sep = '\t', skip=numskip,
			comment.char='', stringsAsFactors=FALSE, check.names=FALSE)
	cols <- colnames(data)
	cols[1] <- 'CHROM'
	colnames(data) <- cols
	return(data)
}
data <- loadVcf('merged.annotated.filtered.vcf')
#data <- loadVcf('merged.all.vcf')
data

plotDepthForSample <- function(config, data, sample)
{
	results <- data.frame()
	for (rowname in rownames(data))
	{
		position <- as.numeric(data[rowname,'POS'])
		values <- data[rowname,sample]
		fields <- strsplit(data[rowname,'FORMAT'],':')[[1]]
		dp <- as.numeric(extractFieldData(fields,values,'DP'))
		results <- rbind(results, data.frame(position=position, dp=dp))
	}
	plot(dp ~ position, results, type='l')
	#return(results)
}
#plotDepthForSample(config,data,'KT9.plasmid.KT9')


extractReplicateFromSampleName <- function(sample)
{
	replicate <- strsplit(sample,'\\.')[[1]][2]
	if (substr(replicate,1,2)=='wk')
		replicate <- as.numeric(substr(c(replicate),3,100)) #drop off the wk prefix and convert to numbers
	return(replicate)
}
#extractReplicateFromSampleName('KT9.plasmid.KT9') # returns plasmid
#extractReplicateFromSampleName('KT9.wk10.KT9') # returns 10


extractVcf <- function(config, data, subject, position)
{
	#get list of samples for given subject
	samples <- getSamplesForSubject(config,subject)
	#row <- data[which(data$POS==position),c('REF','ALT','QUAL','FILTER',samples)]
	fields <- strsplit(data[which(data$POS==position),'FORMAT'],':')[[1]]
	results <- data.frame()
	for (sample in samples)
	{
		replicate <- extractReplicateFromSampleName(sample)
		values <- data[which(data$POS==position),sample]
		ad <- extractFieldData(fields,values,'AD')
		dp <- extractFieldData(fields,values,'DP')
		counts <- as.numeric(splitFields(ad,','))
		refnum <- counts[1]
		altnum <- counts[2]
		total <- refnum + altnum
		freq <- altnum/total#format(altnum/total, digits=3)
		row <- data.frame(replicate=replicate, refnum=refnum, altnum=altnum, total=total, freq=freq, depth=dp)
		results <- rbind(results,row)
	}
	return(results)
	
}
#data.subset <- extractVcf(config,data,'PXB0218-0007',3526)
#data.subset <- extractVcf(config,data,'KT9',3526)


plotVariantFrequenciesForSubjectAndPosition <- function(config, data, subject, position, xlab='Weeks after inoculation', ylab='Variant frequency')
{
	try({
		data.subset <- extractVcf(config,data,subject,position)
		data.subset <- data.subset[complete.cases(data.subset),]
		title <- concat('Subject: ',subject,', nt',position)
		ref <- data[which(data$POS==position),'REF']
		alt <- data[which(data$POS==position),'ALT']
		ylab <- concat(ylab,' (',ref,'/',alt,')')
		if (nrow(data.subset)>0)
		{
			print('')
			print(concat('Subject: ',subject,' position: ',position))
			print(data.subset)
			plot(freq ~ replicate, data.subset, type='b', main=title, ylab=ylab, xlab=xlab, ylim=c(0,1))
		}
	}, silent=FALSE)
}
#plotVariantFrequenciesForSubjectAndPosition(config,data,'PXB0218-0007',3526)

plotVariantFrequenciesForSubject <- function(config, data, subject)
{
	for (position in data[,'POS'])
	{
		plotVariantFrequenciesForSubjectAndPosition(config,data,subject,position)
	}

}
#plotVariantFrequenciesForSubject(config,data,'PXB0218-0007')

plotVariantFrequencies <- function(config,data,position=NULL,filename='d:/temp/variantplots.pdf')
{
	pdf(filename)
	for (subject in config@subjects$subject)
	{
		if (!is.null(position))
			plotVariantFrequenciesForSubjectAndPosition(config,data,subject,position)
		else plotVariantFrequenciesForSubject(config,data,subject)
	}
	dev.off()
}
plotVariantFrequencies(config,data)



#
#for (subject in config@subjects)
#{
#	#data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq',subject)
#	data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq',subject)
#	print(data)
#}
#
#
#
#data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq','PXB0220-0002')
#y <- as.numeric(data[4,splitFields('PXB0220-0002.wk08.KT9:freq,PXB0220-0002.wk09.KT9:freq,PXB0220-0002.wk10.KT9:freq,PXB0220-0002.wk11.KT9:freq,PXB0220-0002.wk12.KT9:freq,PXB0220-0002.wk13.KT9:freq')])
#x <- c(8,9,10,11,12,13)
#plot(y ~ x)
#
#
#data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq')
#
##
##plotVariantFrequenciesForSubject <- function(config, data, subject)
##{
##	samples <- getSampleForSubject(config,subject)
##	fields <- paste(samples,'freq',sep=':')
##	replicates <- getReplicatesForSubject(config,subject)
##	x <- as.numeric(substr(c(replicates),3,100)) #drop off the wk prefix and convert to numbers
##	par(ask=TRUE)
##	for (no in rownames(data))
##	{
##		try({
##				y <- as.numeric(data[no,fields])
##				plot(y ~ x, type='b', main=concat('NT position: ',data[no,'POS']), xlab='Weeks after inoculation',
##						ylab='Variant frequency', ylim=c(0,1))
##			}, silent=FALSE)
##	}
##	par(ask=FALSE)
##}
##plotVariantFrequenciesForSubject(config,data,subject)
#
#plotVariantFrequencies <- function(config, data)
#{
#	for (subject in config@subjects)
#	{
#		plotVariantFrequenciesForSubject(config,data,subject)
#	}
#}
#plotVariantFrequencies(config,data)
#
#convertVcf <- function(filename)
#{
#	data <- loadVcf(filename)
#	data[,]
#	
#}
#data <- convertVcf('merged.filtered.vcf')
#head(data)



#loadVcf <- function(filename, format.fields=NULL, subject=NULL)
#{
#	numskip <- getNumberOfLinesToSkip(filename)	
#	data <- read.table(filename, header=TRUE, encoding='UTF-8', sep = '\t', skip=numskip,
#			comment.char='', stringsAsFactors=FALSE, check.names=FALSE)
#	cols <- colnames(data)
#	cols[1] <- 'CHROM'
#	colnames(data) <- cols
#
#	use.samples <- cols[10:length(cols)]
#	if (!is.null(subject))
#	{
#		subj.samples <- apply(config@samples[which(config@samples$subject %in% splitFields(subject)),c('sample','ref')],1,function(values){concat(values[1],'.',values[2])})
#		use.samples <- intersectValues(use.samples,subj.samples)
#	}
#	info.fields <- getInfoFields(data)
#	
#	#expand the info column into separate columns
#	for (no in rownames(data))
#	{
#		pairs <- strsplit(data[no,'INFO'],';')[[1]]
#		for (pair in pairs)
#		{
#			arr <- strsplit(pair,'=')[[1]]
#			field <- arr[1]
#			value <- arr[2]
#			if (isFieldIncluded(info.fields,field))
#				data[no,field] <- value
#		}
#	}
#	
#	# use the FORMAT field to expand the sample fields
#	for (no in rownames(data))
#	{
#		fields <- strsplit(data[no,'FORMAT'],':')[[1]]
#		if (is.null(format.fields))
#			format.fields <- fields
#		for (sample in use.samples)
#		{
#			values <- strsplit(data[no,sample],':')[[1]]
#			if (length(values)==length(fields))
#			{
#				for (i in 1:length(fields))
#				{
#					field <- fields[i]
#					value <- values[i]
#					if (isFieldIncluded(format.fields,field))
#						data[no,concat(sample,':',field)] <- value
#					if (field=='AD')
#					{
#						if (!is.na(value))
#						{
#							counts=as.numeric(strsplit(value,',')[[1]])
#							data[no,concat(sample,':refnum')] <- counts[1]
#							data[no,concat(sample,':altnum')] <- counts[2]
#							data[no,concat(sample,':total')] <- counts[1] + counts[2]
#							data[no,concat(sample,':freq')] <- counts[2]/(counts[1] + counts[2])
#						}
#					}
#				}
#			}
#		}
#	}
#	
#	fields <- splitFields('CHROM,POS,REF,ALT,QUAL,FILTER')
#	fields <- appendValues(fields, info.fields)
#	for (field in splitFields(format.fields))
#	{
#		for (sample in use.samples)
#		{
#			fields <- appendValues(fields, concat(sample,':',field))
#		}
#	}
#	fields <- removeElements(fields,'Samples') #hack!
#	#print(fields)
#	return(data[,fields])
#	#return(data)
#}
#data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq','PXB0218-0007')
#data

