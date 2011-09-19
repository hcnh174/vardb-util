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
	samples <- config@runs[which(config@runs$subject==subject),'sample']
	#samples <- config@samples[which(config@samples$subject==subject),c('sample','ref')]
	#samples <- unique(paste(samples$sample,samples$ref,sep='.'))
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

getRegionsForSubject <- function(config, subject)
{
	return(unique(config@runs[which(config@runs$subject==subject),'region']))
}
#getRegionsForSubject(config,'10348001')

#############################################################################

#getCodonPositionsForRegion <- function(config, region)
#{
#	ref <- config@regions[region,'ref']	
#	startnt <- config@refs[ref,'startnt']
#	startntrel <- config@refs[ref,'startntrel']
#	print(concat('startnt=',startnt,' startntrel=',startntrel))
#	sequence <- config@refs[ref,'sequence']
#	if (startntrel<0)
#	{
#		diff <- abs(startntrel)
#		startntrel <- startntrel + diff
#		startnt <- startnt + diff
#		sequence <- c2s(s2c(sequence)[(diff+1):nchar(sequence)])
#	}
#	endntrel <- startntrel + nchar(sequence)
#	#print(concat('startnt=',startnt,' startntrel=',startntrel))
#	positions <- data.frame()
#	offset <- (startntrel) %% 3
#	if (offset!=0)
#	{
#		print(concat('adjusting to the start of a codon: offset=',offset))
#		startntrel <- startntrel + 3 - offset
#		startnt <- startnt + 3 - offset
#		if ((startntrel) %% 3!=0)
#			stop(concat('start number is not a multiple of 3: ',startntrel,' in region ',region))
#	}
#	codon <- (startntrel)/3 + 1
#	#print(concat('ref: ',ref,' startnt: ',startnt,', startntrel: ',startntrel,', codon: ',codon))
#	for (position in seq(startntrel,endntrel,3))
#	{
#		index <- position-startntrel+1
#		if (index+2 > nchar(sequence))
#			break
#		refcodon <- extractSequence(sequence, index, index+2)
#		refaa <- translateCodon(refcodon)
#		ntnum <- startnt + position - startntrel
#		positions <- rbind(positions, data.frame(codon=codon, ntnum=ntnum, relpos=position, refcodon=refcodon, refaa=refaa))#position
#		codon <- codon + 1
#	}
#	startaa <- config@regions[region,'startaa']
#	endaa <- config@regions[region,'endaa']
#	if (!is.na(startaa))
#		positions <- positions[which(positions$codon >= startaa),]
#	if (!is.na(endaa))
#		positions <- positions[which(positions$codon <= endaa),]
#	#print(concat('sequence: ',seqinr::c2s(positions$refaa)))
#	return(positions)
#}


getCodonPositionsForRegion <- function(config, region)
{
	ref <- config@regions[region,'ref']	
	startnt <- config@refs[ref,'startnt']
	startntrel <- config@refs[ref,'startntrel']
	print(concat('startnt=',startnt,' startntrel=',startntrel))
	sequence <- config@refs[ref,'sequence']
	positions <- data.frame()
	#codon <- (startntrel-1-(startntrel-1)%%3)/3 + 1	
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
		print(concat('index=',index,' relpos=',relpos,', ntnum=',ntnum,', refcodon=',refcodon,', refaa=',refaa))
		positions <- rbind(positions, data.frame(codon=codon, ntnum=ntnum, relpos=relpos, refcodon=refcodon, refaa=refaa))#position
		codon <- codon + 1
	}
#	startaa <- config@regions[region,'startaa']
#	endaa <- config@regions[region,'endaa']
#	if (!is.na(startaa))
#		positions <- positions[which(positions$codon >= startaa),]
#	if (!is.na(endaa))
#		positions <- positions[which(positions$codon <= endaa),]
	print(concat('sequence: ',seqinr::c2s(positions$refaa)))
	return(positions)
}
#getCodonPositionsForRegion(config,'NS3-36@NS3aa36')
#getCodonPositionsForRegion(config,'NS3-156@NS3aa156')
#positions <- getCodonPositionsForRegion(config,'HBV-RT@HBV-RT')
#positions <- getCodonPositionsForRegion(config,'NS3-156-R@NS3aa156')
#getCodonPositionsForRegion(config,'NS5A-31-R@NS5Aaa31')

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

##################################################################

loadVariantData <- function(sample, ref, merged=FALSE)#, nums)
{
	if (merged==TRUE)
		filename <- concat('variants/',sample,'.',ref,'.merged.txt')
	else filename <- concat('variants/',sample,'.',ref,'.txt')
	#print(paste('Loading file',filename))
	data <- loadDataFrame(filename)
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

extractRefSequence <- function(config, ref, start, end)
{
	return(extractSequence(config@refs[ref,'sequence'],start,end))
}
#extractRefSequence(config,'KT9',3420,5312)

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

#getReferenceCodon <- function(config, region)
#{
#	positions <- getCodonPositionsForRegion(config,region)
#	coords <- as.numeric(strsplit(config@regions[region,'ntfocus'],'-')[[1]])
#	refcodon <- toupper(extractRefSequence(config,ref,coords[1],coords[2]))
#	return(refcodon)
#}
#getReferenceCodon(config,'PXB0219-0011','NS3aa36')

################################################################################

getCodonCountSubset <- function(config, subject, region, filetype, start, end=start, cutoff=0)
{
	filename <- concat(config@counts.dir,'/',subject,'.',filetype,'.txt')
	data <- loadDataFrame(filename)
	data.subset <- data[which(data$region==region & data$aanum>=start & data$aanum<=end & data$count>=cutoff),]
	data.subset$replicate <- factor(data.subset$replicate)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,'10348001','HBV-RT@HBV-RT','codons',200,250)

makeVariantTable <- function(config, type, subject, region, cutoff=0)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	aanum <- as.integer(config@regions[region,'aafocus'])
	positions <- getCodonPositionsForRegion(config,region)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	data.subset <- getCodonCountSubset(config,subject,region,type,aanum,cutoff=cutoff)
	#frmla <- ifelse(type=='codons', as.formula(codon ~ replicate), as.formula(aa ~ replicate))
	if (type=='codons')
		frmla <- as.formula(codon ~ replicate)
	else frmla <- as.formula(aa ~ replicate)
	counts <- cast(data.subset, frmla, value='count', fun.aggregate=function(x) return(x[1])); counts
	counts <- counts[order(counts[,2], decreasing=TRUE),]
	row.names(counts) <- seq(nrow(counts))
	# add an asterisk to indicate the reference codon
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
#makeVariantTable(config,'codons','8538159','NS3-156-R@NS3aa156')

makeVariantTables <- function(config, type, subject=NULL, ...)
{
	if (is.null(subject))
		subjects <- config@subjects
	else subjects <- splitFields(subject)
	tables <- list()
	for (subject in subjects)
	{
		tables[[subject]] <- list()
		for (region in getRegionsForSubject(config,subject))
		{
			try({
				print(concat('subject=',subject,', region=',region))
				tbl <- makeVariantTable(config, type, subject, region, ...)
				tables[[subject]][[region]] <- tbl
			}, silent=FALSE)
		}
	}
	return(tables)
}
tables <- makeVariantTables(config, 'codons')
#tables <- makeVariantTables(config, 'codons', '8538159')
#tables <- makeVariantTables(config, 'codons', '10348001')

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

