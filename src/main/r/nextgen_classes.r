
setClass("nextgenconfig",
		representation(			
				config.dir='character',
				out.dir='character',
				ref.dir='character',
				index.dir='character',
				data='data.frame',
				refs='data.frame',
				regions='data.frame',
				treatments='data.frame',
				titers='data.frame',
				goals='data.frame',
				genes='data.frame',
				positions='list',
				subjects='vector',
				samples='vector',
				groups='vector',
				illumina.dir='character',
				fastq.dir='character',
				bam.dir='character',
				vcf.dir='character',
				qc.dir='character',
				pileup.dir='character',
				tables.dir='character',
				counts.dir='character',
				consensus.dir='character',
				coverage.dir='character',
				tmp.dir='character',
				check.dir='character',
				profile='character',
				trim='logical',
				filter='logical',
				map.quality='character'
		),
		prototype(
				out.dir='out',
				index.dir='indexes',
				goals=data.frame(),
				profile='default',
				trim=TRUE,
				filter=TRUE,
				map.quality='>30'
		)
)

setMethod("initialize", "nextgenconfig", function(.Object, config.dir='.')
{	
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	if (!file.exists(config.dir))
		throw('config directory does not exist: ',config.dir)
	.Object@config.dir <- config.dir
	
	params <- loadDataFrame(concat(.Object@config.dir,'/params.txt'), idcol='name')
	for (name in row.names(params))
	{
		type <- class(slot(.Object,name))
		value <- getField(params,name,'value')
		if (type=='logical')
			value <- as.logical(value)
		slot(.Object,name) <- value
	}	
	.Object@regions <- loadDataFrame(concat(.Object@config.dir,'/regions.txt'), idcol='id')
	.Object@titers <- loadDataFrame(concat(.Object@config.dir,'/titers.txt'))
	.Object@genes <- loadDataFrame(concat(.Object@config.dir,'/genes.txt'), idcol='id')
	.Object@data <- loadDataFrame(concat(.Object@config.dir,'/data.txt'), idcol='id')#
	if (.Object@profile!='default')
	{
		print(concat('subsetting samples from selected profile: ',.Object@profile))
		.Object@profile <- splitFields(.Object@profile)
		print(.Object@profile)
		.Object@data <- .Object@data[which(.Object@data$profile %in% .Object@profile),]
	}
	
	subjects <- loadDataFrame(concat(.Object@config.dir,'/subjects.txt'), idcol='id')
	for (row in row.names(.Object@data))
	{
		subject <- .Object@data[row,'subject']
		col <- .Object@data[row,'column']
		stem <- subject
		if (is.na(col) | col=='')
		{
			col <- subject
			.Object@data[row,'column'] <- col
		}
		else 
		{
			stem <- concat(subject,'.',col)
		}
		
		.Object@data[row,'stem'] <- stem
		ref <- subjects[subject,'ref']
		sample <- concat(stem,'__',ref)
		.Object@data[row,'ref'] <- ref
		.Object@data[row,'sample'] <- sample
		#print(sample)
	}	
	
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
	.Object@tables.dir <- concat(.Object@out.dir,'/tables')
	.Object@consensus.dir <- concat(.Object@out.dir,'/consensus')
	.Object@coverage.dir <- concat(.Object@out.dir,'/coverage')
	.Object@tmp.dir <- concat(.Object@out.dir,'/tmp')
	.Object@check.dir <- concat(.Object@out.dir,'/check')
	
	reffilename <- concat(.Object@config.dir,'/refs.txt')
	data <- loadDataFrame(reffilename, idcol='ref')
	sequences <- readFastaFiles(.Object@config.dir,'^refs.*\\.fasta')
	for (id in names(sequences))
	{
		data[id,'sequence'] <- sequences[[id]]
	}
	.Object@refs <- data
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
