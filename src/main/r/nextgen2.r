#library(gsubfn)
#library(seqinr)
#library(methods)
#library(reshape)
#library(R2wd)

setClass("variantdata", representation(nt="data.frame", codons="data.frame", aa="data.frame"))

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
			subjects='vector',
			samples='vector',
			illumina.dir='character',
			fastq.dir='character',
			bam.dir='character',
			vcf.dir='character',
			qc.dir='character',
			pileup.dir='character',
			counts.dir='character',
			tmp.dir='character',
			trim='logical',
			filter='logical'),
	prototype(
			config.dir='config',	
			out.dir='out',
			index.dir='indexes',
			trim=FALSE,
			filter=FALSE))

setMethod("initialize", "nextgenconfig", function(.Object)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	
	.Object@runs <- loadDataFrame(concat(.Object@config.dir,'/runs.txt'), idcol='run')
	.Object@regions <- loadDataFrame(concat(.Object@config.dir,'/regions.txt'), idcol='region')
	.Object@titers <- loadDataFrame(concat(.Object@config.dir,'/titers.txt'))
	.Object@goals <- loadDataFrame(concat(.Object@config.dir,'/goals.txt'), idcol='goal')
	.Object@samples <- unique(.Object@runs$sample)
	.Object@subjects <- unique(.Object@runs$subject)
	.Object@ref.dir <- concat(.Object@out.dir,'/ref')
	.Object@fastq.dir <- concat(.Object@out.dir,'/fastq')
	.Object@bam.dir <- concat(.Object@out.dir,'/bam')
	.Object@vcf.dir <- concat(.Object@out.dir,'/vcf')
	.Object@qc.dir <- concat(.Object@out.dir,'/qc')
	.Object@counts.dir <- concat(.Object@out.dir,'/counts')
	.Object@pileup.dir <- concat(.Object@out.dir,'/pileup')
	.Object@tmp.dir <- concat(.Object@out.dir,'/tmp')
	
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
#		reffile <- concat(.Object@ref.dir,'/',id,'.fasta')
#		if (!file.exists(reffile))
#		{
#			print(concat('writing ref file ',reffile))
#			write.fasta(s2c(seq), id, file.out=reffile)
#			checkFileExists(reffile)
#		}
	}
	.Object@refs <- data
	.Object
})

writeRefs <- function(config)
{
	for (ref in rownames(config@refs))
	{
		reffile <- concat(config@ref.dir,'/',ref,'.fasta')
		if (!file.exists(reffile))
		{
			print(concat('writing ref file ',reffile))
			seq <- config@refs[ref,'sequence']
			write.fasta(s2c(seq), ref, file.out=reffile)
			checkFileExists(reffile)
		}
	}
}
#writeRefs(config)

#.loadRefs <- function(filename='config/refs.txt', fasta='config/refs.fasta')
#{
#	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
#	data <- loadDataFrame(filename, idcol='ref')
#	sequences <- read.fasta(file = fasta, as.string = TRUE, seqtype = "DNA", forceDNAtolower=TRUE)
#	#remove everything after the | in the sequence ID
#	ids <- names(sequences)
#	ids <- sapply(ids,function(value){return(strsplit(value,'\\|')[[1]][1])}, USE.NAMES=FALSE)
#	names(sequences) <- ids
#	for (id in names(sequences))
#	{
#		seq <- sequences[[id]][1]
#		data[id,'sequence'] <- seq
#		reffile <- concat(config@ref.dir,'/',id,'.fasta')
#		print(concat('writing ref file ',reffile))
#		write.fasta(s2c(seq), ref, file.out=reffile)
#		checkFileExists(reffile)
#	}
#	return(data)
#}
# .loadRefs()

getSamplesForSubject <- function(config, subject)
{
	samples <- config@runs[which(config@runs$subject==subject),'sample']
	return(samples)
}
#getSamplesForSubject(config,'10348001')

getReplicatesForSubject <- function(config, subject)
{
	replicates <- config@samples[which(config@samples$subject==subject),c('replicate')]	
	return(replicates)
}
#getReplicatesForSubject(config,'PXB0220-0002')

get_ref_for_sample <- function(sample)
{
	return(strsplit(sample,'__', fixed=TRUE)[[1]][2])
}
#get_ref_for_sample('110628-2.CTE247-21.1__HCV-NS3-156')

#get_reffile <- function(config, ref)
#{
#	ref.dir <- config@ref.dir
#	reffile <- concat(config@ref.dir,'/',ref,'.fasta')
#	if (!file.exists(reffile))
#	{
#		seq <- config@refs[ref,'sequence']
#		print(concat('writing ref file ',reffile))
#		write.fasta(s2c(seq), ref, file.out=reffile)
#	}
#	checkFileExists(reffile)
#	return(reffile)
#}
##get_reffile(config,'HBV-RT')

get_reffile <- function(config, ref)
{
	reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	checkFileExists(reffile)
	return(reffile)
}
#get_reffile(config,'HBV-RT')

get_regions_for_sample <- function(config, sample)
{
	return(config@runs[config@runs[,'sample']==sample,'region'])
}

getRegionsForSubject <- function(config, subject)
{
	return(unique(config@runs[which(config@runs$subject==subject),'region']))
}
#getRegionsForSubject(config,'10348001')

#############################################################################


#getCodonPositionsForRef <- function(config, ref)
#{
#	startnt <- config@refs[ref,'startnt']
#	startntrel <- config@refs[ref,'startntrel']
#	sequence <- config@refs[ref,'sequence']
#	if (is.na(startnt)) throw('cannot find startnt for ref: ',ref)
#	if (is.na(startntrel)) throw('cannot find startntrel for ref: ',ref)
#	if (is.na(sequence)) throw('cannot find sequence for ref: ',ref)
#	positions <- data.frame()
#	#print(concat('startnt=',startnt,' startntrel=',startntrel))
#	codon <- (startntrel-(startntrel%%3))/3 + 1	
#	if ((startntrel)%%3!=0) # if the startntrel is in the middle of a codon, increment the codon number by 1
#		codon <- codon + 1
#	for (index in 1:nchar(sequence))
#	{
#		relpos <- startntrel + index - 1
#		if (relpos%%3!=0)
#			next
#		ntnum <- startnt + index
#		refcodon <- extractSequence(sequence, index, index+2)
#		if (nchar(refcodon)<3)
#			break
#		refaa <- translateCodon(refcodon)
#		#print(concat('index=',index,' relpos=',relpos,', ntnum=',ntnum,', refcodon=',refcodon,', refaa=',refaa))
#		positions <- rbind(positions, data.frame(codon=codon, ntnum=ntnum, relpos=relpos, refcodon=refcodon, refaa=refaa))#position
#		codon <- codon + 1
#	}
#	return(positions)
#}

getCodonPositionsForRegion <- function(config, region)
{
	ref <- config@regions[region,'ref']
	if (is.na(ref)) throw('cannot find ref for region: ',region)
	startnt <- config@refs[ref,'startnt']
	if (is.na(startnt)) throw('cannot find startnt for ref: ',ref)
	startntrel <- config@regions[region,'startntrel']
	if (is.na(startntrel)) throw('cannot find startntrel for region: ',region)
	sequence <- config@refs[ref,'sequence']	
	if (is.na(sequence)) throw('cannot find sequence for ref: ',ref)
	positions <- data.frame()
	#print(concat('startnt=',startnt,' startntrel=',startntrel))
	codon <- (startntrel-(startntrel%%3))/3 + 1	
	if ((startntrel)%%3!=0) # if the startntrel is in the middle of a codon, increment the codon number by 1
		codon <- codon + 1
	for (index in 1:nchar(sequence))
	{
		relpos <- startntrel + index - 1
		if (relpos%%3!=0)
			next
		ntnum <- startnt + index
		refcodon <- extractSequence(sequence, index, index+2)
		if (nchar(refcodon)<3)
			break
		refaa <- translateCodon(refcodon)
		#print(concat('index=',index,' relpos=',relpos,', ntnum=',ntnum,', refcodon=',refcodon,', refaa=',refaa))
		positions <- rbind(positions, data.frame(codon=codon, ntnum=ntnum, relpos=relpos, refcodon=refcodon, refaa=refaa))#position
		codon <- codon + 1
	}
	return(positions)
}
#getCodonPositionsForRegion(config,'HBVRT')
#getCodonPositionsForRegion(config,'NS3aa36')
#getCodonPositionsForRegion(config,'NS3aa156')
#getCodonPositionsForRegion(config,'NS3aa156R')
#getCodonPositionsForRegion(config,'NS5Aaa31')
#getCodonPositionsForRegion(config,'NS5Aaa93')

###########################################################################

translateSequence <- function(sequence)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	try(return(c2s(translate(s2c(sequence)))))
}
#translateCodon('GGG')

translateCodon <- function(codon)
{
	if (nchar(codon)!=3)
		return(NA)
	return(translateSequence(codon))
}
#translateCodon('GGG')

extractSequence <- function(sequence, start, end)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	if (end>nchar(sequence))
		end <- nchar(sequence)
	sequence <- s2c(sequence)
	return(c2s(sequence[start:end]))
}
#extractSequence(refs['KT9','sequence'],3420,5312)

################################################################################

plotTiter <- function(config, subject)
{
	titers <- config@titers[which(config@titers$subject==subject),]
	titers <- titers[complete.cases(titers),]
	if (nrow(titers)==0)
		return()
	plt <- plot(titer ~ week, titers, type='b', xlim=c(min(titers$week),max(titers$week)), 
			main=concat('Viral titer: ',subject[1]), ylab='Titer (log10)', xlab='Weeks after inoculation')
	treatments <- config@treatments[which(config@treatments$subject==subject),]
	mtext(paste(unique(treatments$drug), collapse=', '))
	for (rowname in rownames(treatments))
	{
		week <- treatments[rowname,'week']
		drug <- treatments[rowname,'drug']
		amount <- treatments[rowname,'amount']
		abline(v=week, col='red', lty=2)
		#mtext(x=week, labels=drug)
	}
	
	# plot the sampling dates as points
	#config@samples[which(config@samples$subject==subject),]
	replicates <- config@runs[which(config@runs$subject==subject),'replicate']
	for (week in replicates)
	{
		titer <- titers[which(titers$week==week),'titer']
		points(y=titer, x=week, cex=2, pch=16)
	}
	return(plt)
}
#plotTiter(config,'PXB0220-0002')

plotTiters <- function(config)
{
	subjects1 <- unique(config@titers$subject)
	subjects2 <- unique(config@subjects$subject)
	subjects <- intersectValues(subjects1, subjects2)
	
	oldpar <- par(mfrow=calcMfrowLayout(length(subjects),2))
	#oldpar <- par(mfrow=c(2,3))
	for (subject in subjects)
	{
		try(plotTiter(config,subject), silent=FALSE)
	}
	par(oldpar)
}
#plotTiters(config)

