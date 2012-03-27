
setClass("nextgenconfig",
		representation(			
				config.dir='character',
				data.dir='character',
				out.dir='character',
				ref.dir='character',
				index.dir='character',
				data='data.frame',
				refs='data.frame',
				regions='data.frame',
				genes='data.frame',
				positions='list',
				subjects='vector',
				samples='vector',
				groups='vector',				
				fastq.dir='character',
				bam.dir='character',
				vcf.dir='character',
				qc.dir='character',
				pileup.dir='character',
				basecounts.dir='character',
				tables.dir='character',
				charts.dir='character',
				counts.dir='character',
				consensus.dir='character',
				coverage.dir='character',
				tmp.dir='character',
				check.dir='character',
				profile='character',
				force='logical',
				trim='logical',
				mask='logical',
				dedup='logical',
				filter='logical',
				map.quality='character',
				minfreq='numeric',
				minreads='numeric',
				minlength='numeric'
		),
		prototype(
				index.dir='indexes',
				profile='default',
				force=FALSE,
				trim=TRUE,
				mask=TRUE,
				dedup=TRUE,
				filter=FALSE,
				map.quality='>30',
				minfreq=0.0,
				minreads=0,
				minlength=34
		)
)

setMethod("initialize", "nextgenconfig", function(.Object, config.dir='.', data.dir='~/nextgen/data', out.dir='out')
{	
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	if (!file.exists(config.dir))
		throw('config directory does not exist: ',config.dir)
	.Object@config.dir <- config.dir
	
	if (!file.exists(data.dir))
		throw('data directory does not exist: ',data.dir)
	.Object@data.dir <- data.dir
	
	if (!file.exists(out.dir))
		throw('out directory does not exist: ',out.dir)
	.Object@out.dir <- out.dir
	
	params <- loadDataFrame(concat(.Object@config.dir,'/params.txt'), idcol='name')
	for (name in row.names(params))
	{
		type <- class(slot(.Object,name))
		value <- getField(params,name,'value')
		if (type=='logical')
			value <- as.logical(value)
		if (type=='numeric')
			value <- as.numeric(value)
		slot(.Object,name) <- value
	}	
	.Object@regions <- loadDataFrame(concat(.Object@config.dir,'/regions.txt'), idcol='id')
	.Object@genes <- loadDataFrame(concat(.Object@config.dir,'/genes.txt'), idcol='id')
	
	.Object@data <- loadDataFrame(concat(.Object@config.dir,'/data.txt'), idcol='id')
	.Object@data <- excludeColumns(.Object@data, 'OLDgroup,OLDtable,OLDcolumn')
	if (is.null(.Object@data$profile))
		.Object@data$profile <- .Object@data$run
	if (.Object@profile!='default')
	{
		printcat('subsetting samples from selected profile: ',.Object@profile)
		.Object@profile <- splitFields(.Object@profile)
		print(.Object@profile)
		.Object@data <- .Object@data[which(.Object@data$profile %in% .Object@profile),]
	}

	#print(head(.Object@data))
	# add run information to the data table
	#return(.Object)
	runs <- loadDataFrame(concat(.Object@config.dir,'/runs.txt'), idcol='run')
	for (id in row.names(.Object@data))
	{
		run <- .Object@data[id,'run']
		checkRowExists(runs,run)
		.Object@data[id,'rundata'] <- runs[run,'folder']
	}

	# add region information to the data table
	for (id in row.names(.Object@data))
	{
		region <- .Object@data[id,'region']
		checkRowExists(.Object@regions,region)
	}
	
	# add subject information to the data table
	subjects <- loadDataFrame(concat(.Object@config.dir,'/subjects.txt'), idcol='id')
	for (id in row.names(.Object@data))
	{
		subject <- .Object@data[id,'subject']
		checkRowExists(subjects,subject)
		.Object@data[id,'group'] <- subjects[subject,'group']
		.Object@data[id,'ref'] <- subjects[subject,'ref']
	}
	
	# add group information to the data table
	groups <- loadDataFrame(concat(.Object@config.dir,'/groups.txt'), idcol='id')
	for (id in row.names(.Object@data))
	{
		group <- .Object@data[id,'group']
		checkRowExists(groups,group)
		.Object@data[id,'description'] <- groups[group,'description']
	}
	
	# add repeat information to the data table
	repeats <- loadDataFrame(concat(.Object@config.dir,'/repeats.txt'))
	for (id in row.names(.Object@data))
	{
		subject <- .Object@data[id,'subject']
		column <- .Object@data[id,'column']
		if (!is.na(column))
		{
			row <- repeats[which(repeats$subject==subject & repeats$name==column),]
			#printcat('subject=',subject,', column=',column)
			#print(row)
			if (is.na(row)[1])
			{
				printcat('subject: ',subject,'; column: ',column)
				throw('cannot find repeated measure with subject ',subject,' and name ',column)
			}
		}
	}
	#print(head(.Object@data))
	
	#set table names - use the group if no repeats; if repeats, use the group plus subject
	for (id in row.names(.Object@data))
	{		
		group <- .Object@data[id,'group']
		subject <- .Object@data[id,'subject']
		column <- .Object@data[id,'column']
		ref <- .Object@data[id,'ref']
		stem <- subject
		if (!is.na(column))
		{
			.Object@data[id,'table'] <- concat(group,'_',subject)
			stem <- concat(subject,'.',column)
		}
		else 
		{
			.Object@data[id,'column'] <- subject
			.Object@data[id,'table'] <- group
		}
		.Object@data[id,'stem'] <- stem
		.Object@data[id,'sample'] <- concat(stem,'__',ref)		
	}
	#return(.Object)
	
	.Object@samples <- unique(.Object@data$sample)
	.Object@subjects <- unique(.Object@data$subject)
	.Object@groups <- unique(.Object@data$group)

	.Object@ref.dir <- concat(.Object@out.dir,'/ref')
	.Object@fastq.dir <- concat(.Object@out.dir,'/fastq')
	.Object@bam.dir <- concat(.Object@out.dir,'/bam')
	.Object@vcf.dir <- concat(.Object@out.dir,'/vcf')
	.Object@qc.dir <- concat(.Object@out.dir,'/qc')
	.Object@counts.dir <- concat(.Object@out.dir,'/counts')
	.Object@pileup.dir <- concat(.Object@out.dir,'/pileup')
	.Object@basecounts.dir <- concat(.Object@out.dir,'/basecounts')
	.Object@tables.dir <- concat(.Object@out.dir,'/tables')
	.Object@charts.dir <- concat(.Object@out.dir,'/charts')
	.Object@consensus.dir <- concat(.Object@out.dir,'/consensus')
	.Object@coverage.dir <- concat(.Object@out.dir,'/coverage')
	.Object@tmp.dir <- concat(.Object@out.dir,'/tmp')
	.Object@check.dir <- concat(.Object@out.dir,'/check')
	
	reffilename <- concat(.Object@config.dir,'/refs.txt')
	data <- loadDataFrame(reffilename, idcol='ref')
	sequences <- readFastaFiles(.Object@config.dir,'^refs.*\\.fasta$')
	for (id in names(sequences))
	{
		data[id,'sequence'] <- sequences[[id]]
	}
	.Object@refs <- data
	
	#needed <- unique(config@data$ref)
	#present <- rownames(config@refs)
	#if (length(removeElements(needed,present)))
	unmatched <- uniqueValues(unique(.Object@data$ref),rownames(.Object@refs))
	if (length(unmatched)>0)
		throw('ref(s) missing: ',unmatched)

	.Object@data <- .Object@data[order(.Object@data$order),]
	
	return(.Object)
})


##############################################################

setClass("variantdata",
	representation(nt="data.frame",
		codons="data.frame",
		aa="data.frame"))

##############################################################

setClass("sampleparams",
	representation(
		group='character',
		subject='character',
		sample='character',
		column='character',
		ref='character',
		region='character'
	)
)

##############################################################


#	
#	for (row in row.names(.Object@data))
#	{
#		# if there is no table listed, use the group name
#		group <- .Object@data[row,'group']
#		if (is.na(.Object@data[row,'table']))
#			.Object@data[row,'table'] <- group
#		
#		subject <- .Object@data[row,'subject']
#		col <- .Object@data[row,'column']
#		stem <- subject
#		if (is.na(col) | col=='')
#		{
#			col <- subject
#			.Object@data[row,'column'] <- col
#		}
#		else 
#		{
#			stem <- concat(subject,'.',col)
#		}
#		
#		.Object@data[row,'stem'] <- stem
#		ref <- subjects[subject,'ref']
#		if (is.na(ref))#check to make sure each subject is defined
#			throw('subject ',subject,' is not defined in subjects.txt')
#		sample <- concat(stem,'__',ref)
#		.Object@data[row,'ref'] <- ref
#		.Object@data[row,'sample'] <- sample
#	}