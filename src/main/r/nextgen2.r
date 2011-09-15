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
			subjects='vector',
			samples='vector',
			illumina.dir='character',
			fastq.dir='character',
			bam.dir='character',
			vcf.dir='character',
			qc.dir='character',
			pileup.dir='character',
			counts.dir='character',
			tmp.dir='character'),
	prototype(
			config.dir='config',	
			out.dir='out',
			index.dir='indexes'))

setMethod("initialize", "nextgenconfig", function(.Object)
{
	.Object@refs <- .loadRefs(concat(.Object@config.dir,'/refs.txt'), concat(.Object@config.dir,'/refs.fasta'))
	.Object@runs <- loadDataFrame(concat(.Object@config.dir,'/runs.txt'), idcol='run')
	.Object@regions <- loadDataFrame(concat(.Object@config.dir,'/regions.txt'), idcol='region')
	.Object@titers <- loadDataFrame(concat(.Object@config.dir,'/titers.txt'))	
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
	.Object
})

.loadRefs <- function(filename='config/refs.txt', fasta='config/refs.fasta')
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	data <- loadDataFrame(filename, idcol='ref')
	sequences <- read.fasta(file = fasta, as.string = TRUE, seqtype = "DNA", forceDNAtolower=TRUE)
	#remove everything after the | in the sequence ID
	ids <- names(sequences)
	ids <- sapply(ids,function(value){return(strsplit(value,'\\|')[[1]][1])}, USE.NAMES=FALSE)
	names(sequences) <- ids
	for (id in names(sequences))
	{
		data[id,'sequence'] <- sequences[[id]][1]
	}
	return(data)
}
# .loadRefs()

getSamplesForSubject <- function(config, subject)
{
	samples <- config@samples[which(config@samples$subject==subject),c('sample','ref')]
	samples <- unique(paste(samples$sample,samples$ref,sep='.'))
	return(samples)
}
#getSamplesForSubject(config,'PXB0220-0002')

getReplicatesForSubject <- function(config, subject)
{
	replicates <- config@samples[which(config@samples$subject==subject),c('replicate')]	
	return(replicates)
}
#getReplicatesForSubject(config,'PXB0220-0002')

get_ref_for_sample <- function(sample)
{
	return(strsplit(sample,'@')[[1]][2])
}
#get_ref_for_sample('110617HBV.HBV07@HBV-RT')

get_reffile <- function(config, ref)
{
	ref.dir <- config@ref.dir
	reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	if (!file.exists(reffile))
	{
		seq <- config@refs[ref,1]	
		print(concat('writing ref file ',reffile))
		write.fasta(s2c(seq), ref, file.out=reffile)
	}
	checkFileExists(reffile)
	return(reffile)
}
#get_reffile(config,'IL28B-70')

get_regions_for_sample <- function(config, sample)
{
	return(config@runs[config@runs[,'sample']==sample,'region'])
}

#########################################################################3

getVarRefName <- function(refid, index)
{
	if (index==1)
		name <- refid
	else name <- paste(refid,'var',(index-1), sep='')
	return(name)
}

# makes variant reference sequences using all possible combinations of specified variants codons
makeVariantsForRef <- function(config,refid)
{
	#print(concat('makeVariantsForRef.refid=',refid))
	refs <- config@refs
	variants <- config@variants
	# get variants for the current reference sequence
	ref.variants <- subset(variants,ref==refid)
	if (nrow(ref.variants)==0)
		return(data.frame())
	# extract reference sequence
	sequence <- tolower(refs[refid,'sequence'])
	#aastart <- as.numeric(refs[refid,'aastart'])
	codons <- strapply(sequence, "...")[[1]] # split into triplets
	sets <- list()
	positions <- getCodonPositionsForRegion(config,refid)
	aastart <- positions[1,'codon']
	for (no in 1:length(codons))
	{
		aa <- no + aastart - 1
		sets[[paste('aa',aa,sep='')]] <- codons[no]
	}
	#print(sets)
	
	#set the codon position as the rowname
	rownames(ref.variants) <- ref.variants$codon
	for (codon in rownames(ref.variants))
	{	
		aa <- paste('aa',codon,sep='')
		sets[[aa]] <- appendUniqueValues(sets[[aa]], tolower(ref.variants[codon,'variants']))
	}
	#print(sets)	
	# expand the grid to try every combination
	varseqs <- expand.grid(sets)
	# convert the grid to strings instead of factors so can be concatenated
	varseqs <- data.frame(lapply(varseqs, as.character), stringsAsFactors=FALSE)
	return(varseqs)
}
#makeVariantsForRef(config,'NS3aa156')
#makeVariantsForRef(config,'KT9')

writeVariantsForRef <- function(config, refid, varseqs)
{
	ref.dir <- config@ref.dir
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		varseq <- paste(varseqs[index,], collapse='')
		filename <- concat(ref.dir,name,'.fasta')
		print(concat('writing variant file ',filename))
		write.fasta(s2c(varseq), name, file.out=filename)
	}
}

# determine varref names by counting how many variants are present for a particular sample
getVarRefNames <- function(config, ref)
{
	#print(concat('getVarRefNames.ref=',ref))
	varseqs <- makeVariantsForRef(config,ref)
	if (nrow(varseqs)==0)
		return(ref)
	refnames <- c()
	for (index in 1:nrow(varseqs))
	{
		refnames <- c(refnames,getVarRefName(ref,index))
	}
	return(refnames)
}
#getVarRefNames(config, 'KT9')

indexReferences <- function(config, refid, varseqs, file)#, ref.dir='ref', index.dir='indexes')
{	
	ref.dir <- config@ref.dir
	index.dir <- config@index.dir
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		appendFile('bowtie-build ',ref.dir,name,'.fasta',' ',index.dir,name, file=file)
	}
	
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		appendFile('bwa index ',ref.dir,name,'.fasta', file=file)
	}
}

makeVariants <- function(config, file='index_refs.txt')
{
	file <- concat(config@out.dir, file)
	system(concat('mkdir ',config@ref.dir))
	cat('', file=file)
	for (refid in rownames(refs))
	{
		varseqs <- makeVariantsForRef(config, refid)
		writeVariantsForRef(config, refid, varseqs)
		indexReferences(config, refid, varseqs, file)
	}
}

####################################################################

#map_sample_reads <- function(config, sample, ref, file)
#{
#	#then loop through each reference
#	for (varref in getVarRefNames(config, ref))
#	{
#		appendFile('python map_reads.py ',sample,' ',varref, file=file)
#	}
#}
#
#map_reads <- function(config, file='out/map_reads.txt')
#{
#	createFile(file=file)
#	appendFile('mkdir sam; rm sam/*', file=file)
#	appendFile('mkdir bam; rm bam/*', file=file)
#	appendFile('mkdir sai; rm sai/*', file=file)
#	appendFile('mkdir unmapped; rm unmapped/*', file=file)
#	appendFile('mkdir variants; rm variants/*', file=file)
#	
#	#samples <- loadSamples()
#	samples <- config@samples
#	for (sample in rownames(samples))
#	{
#		ref <- samples[sample,'ref']
#		map_sample_reads(config,sample,ref,file)
#	}
#}

##################################################################

loadVariantData <- function(sample, ref, merged=FALSE)#, nums)
{
	if (merged==TRUE)
		filename <- concat('variants/',sample,'.',ref,'.merged.txt')
	else filename <- concat('variants/',sample,'.',ref,'.txt')
	#print(paste('Loading file',filename))
	data <- loadDataFrame(filename)
	#print(paste('Loaded',filename))
#	
#	positions <- c()
#	for (num in nums)
#	{
#		position <- (num-1)*3
#		positions <- c(positions, position, position+1, position+2)
#	}
#	
#	data <- subset(data, position %in% positions)
	#data <- subset(data, position >= startpos & position <= endpos)
	return(data)
}


removeAmbiguousCodons <- function(codons)
{
	newvalues <- c()
	for (codon in codons)
	{
		if (regexpr("N", codon) == -1)
			newvalues <- appendValues(newvalues,codon)
	}
	return(newvalues)
}
#removeAmbiguousCodons(splitFields('GCT,GTT,NGT'))

###########################################################################

translateCodon <- function(codon)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	#print(codon)
	if (nchar(codon)!=3)
		return(NA)
	try(return(c2s(translate(s2c(codon)))))
}
#translateCodon('GGG')

extractSequence <- function(sequence, start, end)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	sequence <- s2c(sequence)
	return(c2s(sequence[start:end]))
}
#extractSequence(refs['KT9','sequence'],3420,5312)

extractRefSequence <- function(config, ref, start, end)
{
	return(extractSequence(config@refs[ref,'sequence'],start,end))
}
#extractRefSequence(config,'KT9',3420,5312)

getReferenceCodon <- function(config, subject, region)
{
	ref <- config@ref
	#ref <- min(config@samples[which(config@samples$subject==subject),'ref'])
	coords <- as.numeric(strsplit(config@regions[region,'ntfocus'],'-')[[1]])
	refcodon <- toupper(extractRefSequence(config,ref,coords[1],coords[2]))
	return(refcodon)
}
#getReferenceCodon(config,'PXB0219-0011','NS3aa36')

plotReadDistributions <- function(filename="histograms.pdf")
{
	pdf(filename)
	#par(mfrow=c(2,2))
	#par(ask=TRUE)
	for (smpl in rownames(samples))
	{
		try({
			ref <- samples[smpl,'ref']
			filename <- concat('variants/',smpl,'.',ref,'.txt')
			data <- loadDataFrame(filename)
			table <- createNtCountTable(data, cutoff=0)
			runs.subset <- subset(runs, sample==smpl)
			#par(mfrow=c(1,nrow(runs.subset)))
			par(mfrow=c(2,nrow(runs.subset)))
			for (run in rownames(runs.subset))
			{
				region <- runs.subset[run,'region']
				start <- regions[region,'start']
				end <- regions[region,'end']
				xmin <- min(start)-10
				xmax <- max(end) + 10
				ymax <- max(table$total)+10
				xlab <- concat('position (',xmin,'-',xmax,')')
				#main <- concat(smpl,': ',run,': ',region,': ',start,':',end)
				plot(top1 ~ as.numeric(position), table, ylim=c(0,ymax), xlim=c(xmin,xmax), type='h',
						xlab=xlab, ylab='read coverage', main=run, sub=region)
				plot(top1 ~ as.numeric(position), table, ylim=c(0,500), xlim=c(xmin,xmax), type='h',
						xlab=xlab, ylab='read coverage', main=run, sub=region)
			}
			#par(mfrow=c(1,1))
		}, silent=FALSE)
	}
	#par(ask=FALSE)
	par(mfrow=c(1,1))
	dev.off()
}

####################################################################


extractCodonData <- function(data, ntnum, drop.ambig=FALSE)
{
	#hack!!!
	ntnum <- ntnum - 1
	nt1 <- data[which(data$ntnum==ntnum),'nt']
	nt2 <- data[which(data$ntnum==ntnum+1),'nt']
	nt3 <- data[which(data$ntnum==ntnum+2),'nt']
	codons <- paste(nt1, nt2, nt3, sep='')
	if (drop.ambig)
		codons <- removeAmbiguousCodons(codons)
	return(codons)
}
#extractCodonData(data,3495)

################################################################################

appendSampleParams <- function(counts, params)
{
	counts$subject <- as.character(params@subject)
	counts$sample <- as.character(params@sample)
	counts$replicate <- as.character(params@replicate)
	counts$region <- as.character(params@region)
	cols <- names(counts)
	cols <- c('subject','sample','replicate','region',cols[1:(length(cols)-4)])
	counts <- counts[,cols]
	return(counts)
}

#
#createNtCountTable <- function(config, data, params)
#{
#	ref <- config@regions[which(config@regions==params@region),'ref']
#	startnt <- config@regions[params@region,'startnt']
#	counts <- data.frame()
#	positions <- getCodonPositionsForRegion(config,params@region)
#	positions <- subset(positions,ntnum < max(data$position+startnt))
#	for (ntnum in positions$ntnum)
#	{
#		#try({
#		aanum <- positions[which(positions$ntnum==ntnum),'codon']
#		for (offset in 0:2)
#		{
#			ntnum2 <- ntnum + offset
#			#if (ntnum2 < max(data$position+startnt))
#			nt <- data[which(data$position+startnt==ntnum2),'nt']
#			freqs <- sort(xtabs(as.data.frame(nt)), decreasing=TRUE)
#			total <- sum(freqs)				
#			rank <- 1
#			for (base in names(freqs))
#			{
#				count <- freqs[base]
#				freq <- count/total
#				row <- data.frame(ntnum=ntnum2, aanum=aanum, nt=base, rank=rank, count=count, freq=freq) 
#				counts <- rbind(counts,row)
#				rank <- rank+1
#			}
#		}
#		#}, silent=FALSE)
#	}
#	counts$nt <- as.character(counts$nt)
#	counts$ntnum <- factor(counts$ntnum)
#	counts$aanum <- factor(counts$aanum)
#	counts <- appendSampleParams(counts, params)
#	return(counts)
#}

#getCodonPositionsForRegion <- function(config, region)
#{
#	#region <- 'NS3-156@NS3aa156'
#	ref <- config@regions[which(config@regions==region),'ref']
#	startnt <- config@regions[region,'startnt']
#	startntrel <- config@regions[region,'startntrel']
#	#print(concat(start,':',end))
#	sequence <- config@refs[ref,'sequence']#, regions[region,'start'], regions[region,'end'])
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
#	for (position in seq(startntrel,endntrel-1,3))
#	{
#		try({
#					index <- position-startntrel+1
#					refcodon <- extractSequence(sequence, index, index+2)
#					refaa <- translateCodon(refcodon)
#					positions <- rbind(positions, data.frame(codon=codon, ntnum=startnt+position, relpos=position, refcodon=refcodon, refaa=refaa))
#					codon <- codon + 1
#				})
#	}	
#	return(positions)
#}
##getCodonPositionsForRegion(config,'HBV-RT@HBV-RT')
##getCodonPositionsForRegion(config,'NS3-156@NS3aa156')

getCodonPositionsForRegion <- function(config, region)
{
	ref <- config@regions[region,'ref']
	startnt <- config@refs[ref,'startnt']
	startntrel <- config@refs[ref,'startntrel']
	sequence <- config@refs[ref,'sequence']#, regions[region,'start'], regions[region,'end'])
	endntrel <- startntrel + nchar(sequence)
	#print(sequence)
	#print(translateCodon(sequence))
	positions <- data.frame()
	print(startntrel)
	offset <- (startntrel-1) %% 3
	print(offset)
	if (offset!=0)
	{
		print(concat('adjusting to the start of a codon: offset=',offset))
		startntrel <- startntrel + 3 - offset
		startnt <- startnt + 3 - offset
		if ((startntrel-1) %% 3!=0)
			stop(concat('start number is not a multiple of 3: ',startntrel,' in region ',region))
	}
	codon <- (startntrel-1)/3 + 1
	print(concat('ref: ',ref,' startnt: ',startnt,', startntrel: ',startntrel,', codon: ',codon))
	for (position in seq(startntrel,endntrel-1,3))
	{
		#try({
		index <- position-startntrel+1
		refcodon <- extractSequence(sequence, index, index+2)
		refaa <- translateCodon(refcodon)
		positions <- rbind(positions, data.frame(codon=codon, ntnum=startnt+position, relpos=position, refcodon=refcodon, refaa=refaa))
		codon <- codon + 1
		#})
	}
	startaa <- config@regions[region,'startaa']
	endaa <- config@regions[region,'endaa']
	if (!is.na(startaa))
		positions <- positions[which(positions$codon >= startaa),]
	if (!is.na(endaa))
		positions <- positions[which(positions$codon <= endaa),]
	return(positions)
}
#positions <- getCodonPositionsForRegion(config,'NS3-156@NS3aa156')
#positions <- getCodonPositionsForRegion(config,'NS3-36@NS3aa36')
#getCodonPositionsForRegion(config,'HBV-RT@HBV-RT')

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
#counts <- getNtCounts(data,3731)

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
#counts <- getCodonCounts(data,3602)

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
#counts <- getAaCounts(data,3731)

########################################################################

createNtCountTable <- function(config, data, params)
{
	counts <- data.frame()
	positions <- getCodonPositionsForRegion(config,params@region)
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
	counts <- appendSampleParams(counts, params)
	return(counts)
}
#params <- new('sampleparams',subject='10201689', sample='110719-4.p29.NS3aa36@NS3-36', region='NS3-36@NS3aa36', replicate=1)
#data <- loadPileupData(config,params@sample)
#counts <- createNtCountTable(config,data,params@region)

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
	#counts <- appendSampleParams(counts, params)
	return(counts)
}
#counts <- createCodonCountTable(config,data,params@region)

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
	counts <- appendSampleParams(counts, params)
	return(counts)
}
#counts <- createAminoAcidCountTable(config,data,params@region)

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

loadPileupData <- function(config, sample)
{
	filename <- concat(config@pileup.dir,'/',sample,'.filtered.txt'); print(filename) # load the corresponding data file
	data <- loadDataFrame(filename)
	print(concat('loaded file ',filename,'. contains ',nrow(data),' reads'))
	ref <- get_ref_for_sample(sample)
	startnt <- config@refs[ref,'startnt']
	data$ntnum <- data$position + startnt
	return(data)
}
#data <- loadPileupData(config,'110719-4.p29.NS3aa36@NS3-36')

count_codons_for_sample <- function(config, params, variantdata)
{
	print(concat('count_codons_for_sample: ',params@sample))
	#ref <- get_ref_for_sample(params@sample)
	#filename <- concat(config@pileup.dir,'/',params@sample,'.filtered.txt'); print(filename) # load the corresponding data file
	#data <- loadDataFrame(filename)
	#print(concat('loaded file ',filename,'. contains ',nrow(data),' reads'))
	# each sample has one or more runs targeting different regions
	data <- loadPileupData(config,params@sample)
	for (region in get_regions_for_sample(config,params@sample))
	{
		params@region <- region
		variantdata <- count_codons_for_region(config, data, params, variantdata)
	}
	return(variantdata)
}
#count_codons_for_sample(config,params,variantdata)

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
#counts <- count_codons_for_subject(config, '10201689')
#count_codons_for_subject(config,'10348001')
#
#count_codons <- function(config, params=NULL)
#{
#	#samples <- config@samples
#	print('count_codons')
#	if (is.null(params))
#		params <- new('sampleparams')
#	for (subject in config@subjects)
#	{
#		params@subject <- subject
#		count_codons_for_subject(config, params)
#	}
#}
#count_codons()



################################################################################

# estimates the number of reads based on the size of the file, assuming a ratio of 7757 for uncompressed fastq
estimateReadCount <- function(mb)
{
	return(round(mb*7757))
}
#estimateReadCount(112.5)
#reportAminoAcidChange(config,'PXB0218-0007','NS3aa156',156)

#####################################################################################3

getCodonCountSubset <- function(config, subject, region, filetype, start, end=start, cutoff=0)
{
	filename <- concat(config@counts.dir,'/',subject,'.',filetype,'.txt')
	data <- loadDataFrame(filename)
	data.subset <- data[which(data$region==region & data$aanum>=start &data$aanum<=end & data$count>=cutoff),]
	data.subset$replicate <- factor(data.subset$replicate)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,'PXB0218-0007','NS3aa156','codons',156)

reportAminoAcidChange <- function(config, subject, region, log=FALSE, updown=5)
{
	aanum <- config@regions[region,'aafocus']
	if (!is.integer(aanum))
		stop(concat('cannot find aafocus aanum for region: ',region,' (',aanum,')'))
	aanum <- as.integer(aanum)	
	start <- aanum - updown
	end <- aanum + updown
	
	data.subset <- getCodonCountSubset(config,subject,region,'aa', start, end)
	
	print(data.subset[which(data.subset$aanum==aanum),splitFields('replicate,rank,aa,count,freq')])
	#print(head(data.subset))
	numcol <- length(unique(data.subset$rank))
	col <- gray(numcol:0 / numcol)
	if (log)
		frmla <- as.formula(log10(count) ~ aanum | replicate)
	else frmla <- as.formula(count ~ aanum | replicate)	
	chrt <- barchart(frmla, data.subset, group=rank, horizontal=FALSE, stack=TRUE,
			main=subject, xlab='Amino acid number', ylab='Amino acid count', sub=region,
			col=col, strip=FALSE, strip.left=TRUE, #strip.text = list(cex = 0.75),
			#auto.key = list(space = "right"),
			layout = c(1,length(unique(data.subset$replicate))))
	print(chrt)
	addLine(v=aanum - start + 1 - 0.5, col='red', lty=2)
	addLine(v=aanum - start + 1 + 0.5, col='red', lty=2)
	return(data.subset)
}
#reportAminoAcidChange(config,'PXB0218-0007','NS3aa156')

reportAminoAcidChanges <- function(config, subject=NULL, ...)
{
	if (is.null(subject))
		subjects <- unique(config@samples$subject)
	else subjects <- c(subject)
	for (subject in subjects)
	{
		filename <- concat(config@out.dir,subject,'.aa.pdf')
		pdf(filename)
		samples <- config@samples[which(config@samples$subject==subject),'sample']
		regions <- unique(config@runs[which(config@runs$sample %in% samples),'region'])
		for (region in regions)
		{
			reportAminoAcidChange(config, subject, region, ...)
		}
		dev.off()
	}
}
#reportAminoAcidChanges(config, 'PXB0218-0007', log=FALSE)

########################################################################

reportCodonChange <- function(config, subject, region, log=FALSE, updown=5)
{
	aanum <- config@regions[region,'aafocus']
	if (!is.integer(aanum))
		stop(concat('cannot find aafocus aanum for region: ',region,' (',aanum,')'))
	aanum <- as.integer(aanum)	
	start <- aanum - updown
	end <- aanum + updown
	
	data.subset <- getCodonCountSubset(config,subject,region,'codons', start, end)
	
	print(data.subset[which(data.subset$aanum==aanum),splitFields('replicate,rank,aa,codon,count,freq')])
	#print(head(data.subset))
	numcol <- length(unique(data.subset$rank))
	col <- gray(numcol:0 / numcol)
	if (log)
		frmla <- as.formula(log10(count) ~ aanum | replicate)
	else frmla <- as.formula(count ~ aanum | replicate)
	chrt <- barchart(frmla, data.subset, group=rank, horizontal=FALSE, stack=TRUE,
			main=subject, xlab='Codon number', ylab='Codon count', sub=region,
			col=col, strip=FALSE, strip.left=TRUE, #strip.text = list(cex = 0.75),
			#auto.key = list(space = "right"),
			layout = c(1,length(unique(data.subset$replicate))))
	print(chrt)
	addLine(v=aanum - start + 1 - 0.5, col='red', lty=2)
	addLine(v=aanum - start + 1 + 0.5, col='red', lty=2)
	return(data.subset)
}
#reportCodonChange(config,'PXB0218-0007','NS3aa156')

reportCodonChanges <- function(config, subject=NULL, ...)
{
	if (is.null(subject))
		subjects <- unique(config@samples$subject)
	else subjects <- c(subject)
	for (subject in subjects)
	{
		filename <- concat(config@out.dir,subject,'.codons.pdf')
		pdf(filename)
		samples <- config@samples[which(config@samples$subject==subject),'sample']
		regions <- unique(config@runs[which(config@runs$sample %in% samples),'region'])
		for (region in regions)
		{
			reportCodonChange(config, subject, region, ...)
		}
		dev.off()
	}
}
#reportCodonChanges(config, 'PXB0218-0007')

#######################################################################3

makeCodonBarchart <- function(config, subject, region)
{
	aanum <- as.integer(config@regions[region,'aafocus'])
	data.subset <- getCodonCountSubset(config,subject,region,'codons',aanum)
	chrt <- barchart(count ~ replicate, data.subset, group=codon, horizontal=FALSE, stack=TRUE,
			main=subject, sub=concat('codon ',aanum), xlab='Replicate', ylab='Count',
			auto.key = list(space = "right"))
	print(chrt)
}
#makeCodonBarchart(config,'PXB0218-0007','NS3aa156')

makeAminoAcidBarchart <- function(config, subject, region)
{
	aanum <- as.integer(config@regions[region,'aafocus'])
	data.subset <- getCodonCountSubset(config,subject,region,'aa',aanum)
	chrt <- barchart(count ~ replicate, data.subset, group=aa, horizontal=FALSE, stack=TRUE,
			main=subject, sub=concat('amino acid ',aanum), xlab='Replicate', ylab='Count',
			auto.key = list(space = "right"))
	print(chrt)
}
#makeAminoAcidBarchart(config,'PXB0218-0007','NS3aa156')

makeAminoAcidBarcharts <- function(config, subject=NULL, ...)
{
	if (is.null(subject))
		subjects <- unique(config@samples$subject)
	else subjects <- splitFields(subject)
	if (length(subjects)==1)
		filename <- concat(config@out.dir,concat('barcharts.',subjects[1],'.aa.pdf'))
	else filename <- concat(config@out.dir,'barcharts.aa.pdf')
	pdf(filename)
	for (subject in subjects)
	{
		#filename <- concat(config@out.dir,subject,'.aa2.pdf')
		#pdf(filename)
		samples <- config@samples[which(config@samples$subject==subject),'sample']
		regions <- unique(config@runs[which(config@runs$sample %in% samples),'region'])
		for (region in regions)
		{
			makeAminoAcidBarchart(config, subject, region, ...)
		}
		#dev.off()
	}
	dev.off()
}
#makeAminoAcidBarcharts(config, 'PXB0218-0007')
#makeAminoAcidBarcharts(config)

makeCodonBarcharts <- function(config, subject=NULL, ...)
{
	if (is.null(subject))
		subjects <- unique(config@samples$subject)
	else subjects <- splitFields(subject)
	if (length(subjects)==1)
		filename <- concat(config@out.dir,concat('barcharts.',subjects[1],'.codons.pdf'))
	else filename <- concat(config@out.dir,'barcharts.codons.pdf')
	pdf(filename)
	for (subject in subjects)
	{
		#filename <- concat(config@out.dir,subject,'.codons2.pdf')
		#pdf(filename)
		samples <- config@samples[which(config@samples$subject==subject),'sample']
		regions <- unique(config@runs[which(config@runs$sample %in% samples),'region'])
		for (region in regions)
		{
			makeCodonBarchart(config, subject, region, ...)
		}
		#dev.off()
	}
	dev.off()
}
#makeCodonBarcharts(config, 'PXB0218-0007')
#makeCodonBarcharts(config)

#############################################################3

makeVariantTable <- function(config, type, subject, region, cutoff=0)
{
	aanum <- as.integer(config@regions[region,'aafocus'])
	data.subset <- getCodonCountSubset(config,subject,region,type,aanum,cutoff=cutoff)
	if (type=='codons')
		frmla <- as.formula(codon ~ replicate)
	else frmla <- as.formula(aa ~ replicate)
	counts <- cast(data.subset, frmla, value='count', fun.aggregate=function(x) return(x[1])); counts
	counts <- counts[order(counts[,2], decreasing=TRUE),]
	row.names(counts) <- seq(nrow(counts))
	#row <- data.frame(codon="Total", t(colSums(counts[,-1], na.rm=TRUE)))
	#colnames(row) <- colnames(counts)
	#counts <- rbind(counts,row)
	
	# add an asterisk to indicate the reference codon
	refcodon <- getReferenceCodon(config,subject,region)
	if (type=='codons')
	{
		if (length(which(counts$codon==refcodon))==1)
			counts[which(counts$codon==refcodon),'codon'] <- concat(refcodon,'*')
		else
		{
			row <- counts[1,]
			row[,'codon'] <- concat(refcodon,'*')
			row[,-1] <- NA
			counts <- rbind(counts,row)
		}			
	}
	if (type=='aa')
	{
		refaa <- translateCodon(refcodon)
		if (length(which(counts$aa==refaa))==1)
			counts[which(counts$aa==refaa),'aa'] <- concat(refaa,'*')
		else
		{
			row <- counts[1,]
			row[,'aa'] <- concat(refaa,'*')
			row[,-1] <- NA
			counts <- rbind(counts,row)
		}	
	}
	return(counts)
}
#makeVariantTable(config,'codons','PXB0219-0011','NS3aa36')
#makeVariantTable(config,'codons','PXB0219-0018','NS5Aaa31')

makeVariantTables <- function(config, type, subject=NULL, ...)
{
	if (is.null(subject))
		subjects <- unique(config@samples$subject)
	else subjects <- splitFields(subject)
	tables <- list()
	for (subject in subjects)
	{
		#filename <- concat(config@out.dir,subject,',',type,'.txt')
		samples <- config@samples[which(config@samples$subject==subject),'sample']
		regions <- unique(config@runs[which(config@runs$sample %in% samples),'region'])
		tables[[subject]] <- list()
		for (region in regions)
		{
			tbl <- makeVariantTable(config, type, subject, region, ...)
			tables[[subject]][[region]] <- tbl
		}
	}
	return(tables)
}
#tables <- makeVariantTables(config, 'codons', 'KT9')#PXB0218-0007

appendVariantTablesToWord <- function(tables)
{
	for (subject in names(tables))
	{
		wdHeading(level=2,concat('Subject: ',subject))
		for (region in names(tables[[subject]]))
		{
			wdHeading(level=2,concat('Region: ',region))
			tbl <- tables[[subject]][[region]]
			tbl <- replaceNAs(tbl)
			tbl <- format(tbl)			
			wdTable(tbl)
		}
		wdPageBreak()
	}
}


outputVariantTablesToWord <- function(subjects=NULL, filename='tables.doc',...)
{
	filename <- concat(getwd(),'/',config@out.dir,filename)
	
	#codon.tables <- makeCodonTables(config,subjects,...)
	#aa.tables <- makeAminoAcidTables(config,subjects,...)
	
	codon.tables <- makeVariantTables(config,'codons',subjects,...)
	aa.tables <- makeVariantTables(config,'aa',subjects,...)
	
	wdGet(visible=FALSE)
	wdNewDoc(filename)
	wdSection('Codon tables', newpage=FALSE)
	appendVariantTablesToWord(codon.tables)
	wdSection('Amino acid tables', newpage=TRUE)
	appendVariantTablesToWord(aa.tables)
	wdSave(filename)
	wdQuit()
}
#outputVariantTablesToWord()

appendVariantTablesToLatex <- function(config,tables)
{
	for (subject in names(tables))
	{
		#print(concat('Mutation: ',as.character(config@subjects[subject,'mutation'][[1]])))
		for (region in names(tables[[subject]]))
		{
			tbl <- tables[[subject]][[region]]
			caption <- concat('Subject: ',as.character(subject[1]),', Region: ',as.character(region[1]))
			caption <- concat(caption,'\\newline  Description: ',as.character(config@subjects[subject,'description']))
			caption <- concat(caption,'\\newline')
			xtbl <- xtable(tbl, caption=caption, digits=0)
			print(xtbl, include.rownames=FALSE, caption.placement='top', latex.environments='flushleft')
		}
	}
}

getVariantCounts <- function(config, subject, replicate, ranges=NULL)
{
	#print(concat(subject,'.',replicate))
	filename <- concat(config@counts.dir,subject,'.nt.txt')
	#print(filename)
	data <- loadDataFrame(filename)
	data.subset <- data[which(data$replicate==replicate),]
	if (!is.null(ranges))
	{
		ntnums <- parseRanges(ranges)
		data.subset <- data.subset[which(data.subset$ntnum %in% ntnums),]
	}
	if (nrow(data.subset)==0)
		throw(concat('cannot find data any rows for sample ',subject,'.',replicate))
	counts <- cast(data.subset, ntnum ~ rank, value='count', fun.aggregate=function(x) return(x[1]))#; counts
	counts <- replaceNAs(counts, replacestr=0)
	if (ncol(counts)<4)
	{
		print(counts)
		throw(concat('not enough cols for sample ',subject,'.',replicate,': ',ncol(counts)))
	}
	counts$variants <- apply(counts[,3:ncol(counts)], 1, sum)
	counts$total <- apply(counts[,2:ncol(counts)], 1, sum)
	counts$freq <- counts$variants/counts$total
	return(counts)
}
#counts <- getVariantCounts(config,'PXB0218-0007','wk15','6300-6800')
#counts <- getVariantCounts(config,'PXB0218-0007','wk14','3490-4100'); head(counts)

estimateSequencingError <- function(config, subject, ranges=NULL)#, refsample='KT9', refreplicate='plasmid')
{
	replicates <- config@samples[which(config@samples$subject==subject),'replicate']
	
	replicate <- replicates[1]
	counts <- getVariantCounts(config,subject,replicate,ranges)
	
	sample <- concat(subject,'.',replicate)
	freq.summary <- summary(counts$freq)
	
	colors <- c('blue','red','green','orange','brown')
	main <- concat('Variant frequency by position for ',subject)
	if (!is.null(ranges))
		main <- concat(main,' (',ranges,')')
	print(xyplot(freq ~ ntnum, counts, type='l', ylim=c(0,0.1), main=main, col=colors[1],
					ylab='Variant frequency', xlab='NT position', sub=concat('Median: ',freq.summary[3])))
	#add lines to show the median error rate
	addLine(h=freq.summary[3], col='lightgrey')
	addLine(h=freq.summary[2], col='lightgrey', lty=2)
	addLine(h=freq.summary[5], col='lightgrey', lty=2)
	
	for (region in config@runs[which(config@runs$sample==sample),'region'])
	{
		start <- as.integer(config@regions[region,'start'])
		end <- as.integer(config@regions[region,'end'])
		# add gray lines around the aafocus regions
		addLine(v=start, col='lightgrey')
		addLine(v=end, col='lightgrey')
		aanum <- as.integer(config@regions[region,'aafocus'])
		focusnts <- unique(data.subset[which(data.subset$aanum==aanum),'ntnum'])
		#print(focusnts)
		if (length(focusnts)>0)
		{
			# add yellow lines to show the aafocus positions 
			addLine(v=min(focusnts, na.rm=TRUE), col='yellow', lty=2)
			addLine(v=max(focusnts, na.rm=TRUE), col='yellow', lty=2)
		}
	}
	
	trellis.focus("panel",1,1,highlight = FALSE)
	for (i in 2:length(replicates))
	{
		try({
					replicate <- replicates[i]
					counts2 <- getVariantCounts(config,subject,replicate,ranges)
					panel.lines(x=counts2$ntnum, y=counts2$freq, col=colors[i])
				}, silent=FALSE)
	}
	if (subject!='KT9')
	{
		counts2 <- getVariantCounts(config,'KT9','plasmid',ranges)
		panel.lines(x=counts2$ntnum, y=counts2$freq, col='lightgrey')
	}
	trellis.unfocus()
	
	return(counts)
}
#counts <- estimateSequencingError(config,'PXB0218-0007','3490-4100')#'6300-6800'

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
	samples <- config@samples[which(config@samples$subject==subject),]
	for (week in samples$value)
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


#
#createCodonCountTable <- function(config, data, params)
#{
#	ref <- config@regions[which(config@regions==params@region),'ref']
#	startnt <- config@regions[params@region,'startnt']
#	
#	counts <- data.frame()
#	positions <- getCodonPositionsForRegion(config,params@region)
#	positions <- subset(positions,ntnum < max(data$position+startnt))
#	for (ntnum in positions$ntnum)
#	{
#		#print(concat('ntnum=',ntnum,', startnt=',startnt,', diff=',(ntnum-startnt)))
#		codons <- extractCodonData(data,ntnum-startnt,params@drop.ambig)
#		if (length(codons)==0)
#			next
#		#print(codons)
#		aanum <- positions[which(positions$ntnum==ntnum),'codon']
#		print(paste('aanum=',aanum,'ntnum=',ntnum,'numcodons=',length(codons)))
#		freqs <- sort(xtabs(as.data.frame(codons)), decreasing=TRUE)
#		total <- sum(freqs)		
#		rank <- 1
#		for (codon in names(freqs))
#		{
#			count <- freqs[codon]
#			freq <- count/total
#			row <- data.frame(ntnum=ntnum, aanum=aanum, codon=codon, rank=rank, count=count, freq=freq) #sample=sample, 
#			counts <- rbind(counts,row)
#			rank <- rank +1
#		}
#	}
#	counts$codon <- as.character(counts$codon)
#	counts$aa <- sapply(counts$codon,translateCodon)
#	counts$ntnum <- factor(counts$ntnum)
#	counts$aanum <- factor(counts$aanum)
#	counts <- appendSampleParams(counts, params)
#	return(counts)
#}
##counts <- createCodonCountTable(config,data,params)
#
#createAminoAcidCountTable <- function(config, data, params)
#{
#	ref <- config@regions[which(config@regions==params@region),'ref']
#	startnt <- config@regions[params@region,'startnt']
#	
#	counts <- data.frame()
#	positions <- getCodonPositionsForRegion(config,params@region)
#	positions <- subset(positions,ntnum < max(data$position+startnt))
#	for (ntnum in positions$ntnum)
#	{		
#		codons <- extractCodonData(data,ntnum-startnt,params@drop.ambig)
#		if (length(codons)==0)
#			next
#		aanum <- positions[which(positions$ntnum==ntnum),'codon']
#		aa <- sapply(codons,translateCodon)
#		freqs <- sort(xtabs(as.data.frame(aa)), decreasing=TRUE)
#		total <- sum(freqs)
#		rank <- 1
#		for (aa in names(freqs))
#		{
#			count <- freqs[aa]
#			freq <- count/total
#			row <- data.frame(ntnum=ntnum, aanum=aanum, aa=aa, rank=rank, count=count, freq=freq)
#			counts <- rbind(counts,row)
#			rank <- rank +1
#		}
#	}
#	counts$ntnum <- factor(counts$ntnum)
#	counts$aanum <- factor(counts$aanum)
#	counts <- appendSampleParams(counts, params)
#	return(counts)
#}
##counts <- createAminoAcidCountTable(config,data,params)
