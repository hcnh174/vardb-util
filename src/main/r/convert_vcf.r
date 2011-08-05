source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'common.r',sep=''))
loadUtilFiles('nextgen')

setCurDir('nextgen2')
config <- new('nextgenconfig')


setwd('n:/vcf/')
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


loadVcf <- function(filename, format.fields=NULL, subject=NULL)
{
	numskip <- getNumberOfLinesToSkip(filename)	
	data <- read.table(filename, header=TRUE, encoding='UTF-8', sep = '\t', skip=numskip,
			comment.char='', stringsAsFactors=FALSE, check.names=FALSE)
	cols <- colnames(data)
	cols[1] <- 'CHROM'
	colnames(data) <- cols

	use.samples <- cols[10:length(cols)]
	if (!is.null(subject))
	{
		subj.samples <- apply(config@samples[which(config@samples$subject %in% splitFields(subject)),c('sample','ref')],1,function(values){concat(values[1],'.',values[2])})
		use.samples <- intersectValues(use.samples,subj.samples)
	}
	info.fields <- getInfoFields(data)
	
	#expand the info column into separate columns
	for (no in rownames(data))
	{
		pairs <- strsplit(data[no,'INFO'],';')[[1]]
		for (pair in pairs)
		{
			arr <- strsplit(pair,'=')[[1]]
			field <- arr[1]
			value <- arr[2]
			if (isFieldIncluded(info.fields,field))
				data[no,field] <- value
		}
	}
	
	# use the FORMAT field to expand the sample fields
	for (no in rownames(data))
	{
		fields <- strsplit(data[no,'FORMAT'],':')[[1]]
		if (is.null(format.fields))
			format.fields <- fields
		for (sample in use.samples)
		{
			values <- strsplit(data[no,sample],':')[[1]]
			if (length(values)==length(fields))
			{
				for (i in 1:length(fields))
				{
					field <- fields[i]
					value <- values[i]
					if (isFieldIncluded(format.fields,field))
						data[no,concat(sample,':',field)] <- value
					if (field=='AD')
					{
						if (!is.na(value))
						{
							counts=as.numeric(strsplit(value,',')[[1]])
							data[no,concat(sample,':refnum')] <- counts[1]
							data[no,concat(sample,':altnum')] <- counts[2]
							data[no,concat(sample,':total')] <- counts[1] + counts[2]
							data[no,concat(sample,':freq')] <- counts[2]/(counts[1] + counts[2])
						}
					}
				}
			}
		}
	}
	
	fields <- splitFields('CHROM,POS,REF,ALT,QUAL,FILTER')
	fields <- appendValues(fields, info.fields)
	for (field in splitFields(format.fields))
	{
		for (sample in use.samples)
		{
			fields <- appendValues(fields, concat(sample,':',field))
		}
	}
	fields <- removeElements(fields,'Samples') #hack!
	#print(fields)
	return(data[,fields])
	#return(data)
}
data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq','PXB0218-0007')
data

for (subject in config@subjects)
{
	#data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq',subject)
	data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq',subject)
	print(data)
}



data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq','PXB0220-0002')
y <- as.numeric(data[4,splitFields('PXB0220-0002.wk08.KT9:freq,PXB0220-0002.wk09.KT9:freq,PXB0220-0002.wk10.KT9:freq,PXB0220-0002.wk11.KT9:freq,PXB0220-0002.wk12.KT9:freq,PXB0220-0002.wk13.KT9:freq')])
x <- c(8,9,10,11,12,13)
plot(y ~ x)


data <- loadVcf('merged.annotated.filtered.vcf','AD,DP,freq')


plotVariantFrequenciesForSubject <- function(config, data, subject)
{
	samples <- getSampleForSubject(config,subject)
	fields <- paste(samples,'freq',sep=':')
	replicates <- getReplicatesForSubject(config,subject)
	x <- as.numeric(substr(c(replicates),3,100)) #drop off the wk prefix and convert to numbers
	par(ask=TRUE)
	for (no in rownames(data))
	{
		try({
				y <- as.numeric(data[no,fields])
				plot(y ~ x, type='b', main=concat('NT position: ',data[no,'POS']), xlab='Weeks after inoculation',
						ylab='Variant frequency', ylim=c(0,1))
			}, silent=FALSE)
	}
	par(ask=FALSE)
}
#plotVariantFrequenciesForSubject(config,data,subject)

plotVariantFrequencies <- function(config, data)
{
	for (subject in config@subjects)
	{
		plotVariantFrequenciesForSubject(config,data,subject)
	}
}
plotVariantFrequencies(config,data)

convertVcf <- function(filename)
{
	data <- loadVcf(filename)
	data[,]
	
}
data <- convertVcf('merged.filtered.vcf')
head(data)

