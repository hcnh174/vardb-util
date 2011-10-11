
setClass("nextgenconfig",
		representation(			
				config.dir='character',
				out.dir='character',
				ref.dir='character',
				index.dir='character',
				runs='data.frame',
				refs='data.frame',
				regions='data.frame',
				treatments='data.frame',
				titers='data.frame',
				goals='data.frame',
				features='data.frame',
				positions='list',
				subjects='vector',
				samples='vector',
				illumina.dir='character',
				fastq.dir='character',
				bam.dir='character',
				vcf.dir='character',
				qc.dir='character',
				pileup.dir='character',
				tables.dir='character',
				counts.dir='character',
				tmp.dir='character',
				trim='logical',
				filter='logical',
				map.quality='character'
		),
		prototype(
				#config.dir='config',	
				out.dir='out',
				index.dir='indexes',
				goals=data.frame(),
				trim=TRUE,
				filter=TRUE,
				map.quality='>30'
		)
)

setMethod("initialize", "nextgenconfig", function(.Object, config.dir='config')
{	
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	.Object@config.dir <- config.dir
	.Object@runs <- loadDataFrame(concat(.Object@config.dir,'/runs.txt'), idcol='run')
	.Object@regions <- loadDataFrame(concat(.Object@config.dir,'/regions.txt'), idcol='region')
	.Object@titers <- loadDataFrame(concat(.Object@config.dir,'/titers.txt'))
	.Object@goals <- loadDataFrame(concat(.Object@config.dir,'/goals.txt'), idcol='goal')
	.Object@features <- loadDataFrame(concat(.Object@config.dir,'/features.txt'), idcol='feature')
	.Object@samples <- unique(.Object@runs$sample)
	.Object@subjects <- unique(.Object@runs$subject)
	.Object@ref.dir <- concat(.Object@out.dir,'/ref')
	.Object@fastq.dir <- concat(.Object@out.dir,'/fastq')
	.Object@bam.dir <- concat(.Object@out.dir,'/bam')
	.Object@vcf.dir <- concat(.Object@out.dir,'/vcf')
	.Object@qc.dir <- concat(.Object@out.dir,'/qc')
	.Object@counts.dir <- concat(.Object@out.dir,'/counts')
	.Object@pileup.dir <- concat(.Object@out.dir,'/pileup')
	.Object@tables.dir <- concat(.Object@out.dir,'/tables')
	.Object@tmp.dir <- concat(.Object@out.dir,'/tmp')
	
	params <- loadDataFrame(concat(.Object@config.dir,'/params.txt'), idcol='name')
	for (name in row.names(params))
	{
		type <- class(slot(.Object,name))
		value <- getField(params,name,'value')
		if (type=='logical')
			value <- as.logical(value)
		slot(.Object,name) <- value
	}
	
	reffilename <- concat(.Object@config.dir,'/refs.txt')
	fastafilename <- concat(.Object@config.dir,'/refs.fasta')
	data <- loadDataFrame(reffilename, idcol='ref')
	sequences <- read.fasta(file = fastafilename, as.string = TRUE, seqtype = "DNA", forceDNAtolower=TRUE)
	#remove everything after the | in the sequence ID
	ids <- names(sequences)
	ids <- sapply(ids,function(value){return(strsplit(value,'\\|')[[1]][1])}, USE.NAMES=FALSE)
	names(sequences) <- ids
	for (id in names(sequences))
	{
		seq <- sequences[[id]][1]
		data[id,'sequence'] <- seq
	}
	.Object@refs <- data
	.Object
})

##############################################################

setClass("variantdata",
		representation(nt="data.frame",
			codons="data.frame",
			aa="data.frame"))

##############################################################

setClass("sampleparams",
		representation(subject='character',
				sample='character',
				replicate='numeric',
				ref='character',
				region='character',
				drop.ambig='logical',
				nt.cutoff='numeric'),
		prototype(drop.ambig=TRUE,
				nt.cutoff=0))
