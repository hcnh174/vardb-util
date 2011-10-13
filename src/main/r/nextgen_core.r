loadConfig <- function(dir=NULL)
{
	if (is.null(dir))
		dir <- concat('../config/',getCurDir())
	config <- new('nextgenconfig',config.dir=dir)
	try({config <- preloadCodonPositions(config)})
	return(config)
}

writeRefs <- function(config)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
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

getSamplesForSubject <- function(config, subject)
{
	samples <- config@data[which(config@data$subject==subject),'sample']
	return(samples)
}
#getSamplesForSubject(config,'10348001')

getReplicatesForSubject <- function(config, subject)
{
	replicates <- config@samples[which(config@samples$subject==subject),c('replicate')]	
	return(replicates)
}
#getReplicatesForSubject(config,'PXB0220-0002')

getRefForSample <- function(sample)
{
	ref <- strsplit(sample,'__', fixed=TRUE)[[1]][2] #use the part after the delimiter (__)
	ref <- strsplit(ref,'.', fixed=TRUE)[[1]][1] # remove any extensions
	return(ref)
}
#getRefForSample('10348001.20040315__HBV-RT.filtered')

#get_ref_for_subject <- function(config, subject, region)
#{
#	ref <- unique(config@data[which(config@data$subject==subject & config@data$region==region),'ref'])
#	if (length(ref)>1)
#		throw('multiple refs found for subject+region: ',joinFields(ref,','))
#	return(ref)
#}
#get_ref_for_subject(config,'PXB0218-0007','NS3aa156')

getRefForSubject <- function(config, subject)
{
	ref <- unique(config@data[which(config@data$subject==subject),'ref'])
	if (length(ref)>1)
		throw('multiple refs found for subject+region: ',joinFields(ref,','))
	return(ref)
}
#getRefForSubject(config,'CTE247-21')

getRefForGroup <- function(config, group)
{
	ref <- unique(config@data[which(config@data$group==group),'ref'])
	if (length(ref)>1)
		throw('multiple refs found for group+region: ',joinFields(ref,','))
	return(ref)
}
#getRefForGroup(config,'G9')

getRefFile <- function(config, ref)
{
	reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	checkFileExists(reffile)
	return(reffile)
}
#getRefFile(config,'HBV-RT')

getRegionsForSample <- function(config, sample)
{
	return(unique(config@data[which(config@data$sample==sample),'region']))
}
#getRegionsForSample(config,'PXB0220-0030.8__HCV-KT9')

getRegionsForSubject <- function(config, subject)
{
	return(unique(config@data[which(config@data$subject==subject),'region']))
}
#getRegionsForSubject(config,'10348001')

getRegionsForGroup <- function(config, group)
{
	return(unique(config@data[which(config@data$group==group),'region']))
}
#getRegionsForGroup(config,'G9')


#############################################################################

getCodonPositionsForRegion <- function(config, region, ref)
{
	gene <- getField(config@regions,region,'gene')
	region.start <- getField(config@regions,region,'start')
	region.end <- getField(config@regions,region,'end')
	region.focus <- getField(config@regions,region,'focus')
	positions <- getCodonPositionsForGene(config,gene,ref)
	positions <- positions[which(positions$codon>=region.start & positions$codon<=region.end),]
	positions$focus <- ifelse(positions$codon==region.focus,'*','')
	return(positions)
}
#getCodonPositionsForRegion(config,'NS3aa156','HCV-KT9')

#getCodonPositionsForRegion <- function(config, region)
#{
#	#first check if it is cached in the config
#	positions <- config@positions[[region]]
#	if (is.null(positions))
#		positions <- calculateCodonPositionsForRegion(config,region)
#	return(positions)
#}
#getCodonPositionsForRegion(config,'NS3aa156')

#calculateCodonPositionsForGene <- function(config, gene)
#{
#	feature <- getField(config@regions,region,'feature')
#	ref <- getField(config@regions,region,'ref')
#	ref.start <- getField(config@refs,ref,'start')
#	feature.start <- getField(config@features,feature,'start')
#	feature.end <- getField(config@features,feature,'end')
#	region.start <- getField(config@regions,region,'startaa')
#	region.end <- getField(config@regions,region,'endaa')
#	region.focus <- getField(config@regions,region,'focusaa')	
#	sequence <- getField(config@refs,ref,'sequence')
#	
#	startntrel <- ref.start-feature.start
#	printParams(feature,ref,ref.start,feature.start,feature.end,startntrel,region.start,region.end,region.focus)
#	
#	codon <- (startntrel-(startntrel%%3))/3 + 1	
#	if ((startntrel)%%3!=0) # if the startntrel is in the middle of a codon, increment the codon number by 1
#		codon <- codon + 1
#	positions <- data.frame()
#	for (index in 1:nchar(sequence))
#	{
#		relpos <- startntrel + index - 1
#		if (relpos%%3!=0)
#			next
#		ntnum <- feature.start+((codon-1)*3)
#		if (index+2>nchar(sequence))
#			break
#		refcodon <- extractSequence(sequence, index, index+2)
#		refaa <- translateCodon(refcodon)
#		positions <- rbind(positions, data.frame(codon=codon, ntnum=ntnum, refcodon=refcodon, refaa=refaa))#position,relpos=relpos, 
#		codon <- codon + 1
#	}
#	positions <- positions[which(positions$codon>=region.start & positions$codon<=region.end),]
#	positions$focus <- ifelse(positions$codon==region.focus,'*','')
#	return(positions)
#}
#

getCodonPositionsForGene <- function(config, gene, ref)
{
	#first check if it is cached in the config
	positions <- config@positions[[gene]][[ref]]
	if (is.null(positions))
		positions <- calculateCodonPositionsForGene(config,gene,ref)
	return(positions)
}
#getCodonPositionsForGene(config,gene,ref)

calculateCodonPositionsForGene <- function(config, gene, ref)
{
	id <- config@genes[which(config@genes$gene==gene & config@genes$ref==ref),'id']
	gene.start <- config@genes[id,'start']
	gene.end <- config@genes[id,'end']
	ref.start <- getField(config@refs,ref,'start')
	refseq <- getField(config@refs,ref,'sequence')
	sequence <- extractSequence(refseq, gene.start, gene.end)
	positions <- splitCodons(sequence, start=gene.start)
	return(positions)
}
#gene <- 'HCV-NS3'; ref <- 'HCV-KT9'
#positions <- calculateCodonPositionsForGene(config,gene,ref)

#calculateCodonPositionsForRegion(config,'NS3aa156')
#calculateCodonPositionsForRegion(config,'NS3aa36')
#calculateCodonPositionsForRegion(config,'HBVRT')
#calculateCodonPositionsForRegion(config,'NS3aa156R')
#calculateCodonPositionsForRegion(config,'NS5Aaa31')
#calculateCodonPositionsForRegion(config,'NS5Aaa93')

#preloadCodonPositions <- function(config)
#{
#	for (region in unique(config@data$region))
#	{
#		config@positions[[region]] <- calculateCodonPositionsForRegion(config,region)
#	}
#	return(config)
#}
preloadCodonPositions <- function(config)
{
	for (id in config@genes$id)
	{
		gene <- config@genes[id,'gene']
		ref <- config@genes[id,'ref']
		config@positions[[gene]][[ref]] <- calculateCodonPositionsForGene(config,gene,ref)
	}
	return(config)
}
#config <- preloadCodonPositions(config)

###########################################################################

cleanSequence <- function(sequence)
{
	#library(gsubfn) 
	sequence <- gsub ("[0-9]+", "\\1", sequence, ignore.case=T, perl=T)
	sequence <- gsub ("\\s+", "\\1", sequence, ignore.case=T, perl=T)
	return(sequence)
}

translateSequence <- function(sequence)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	try(return(c2s(translate(s2c(sequence)))))
}
#translateCodon('GGG')

translateCodon <- function(codon)
{
	codon <- as.character(codon)
	if (nchar(codon)!=3)
		return(NA)
	return(translateSequence(codon))
}
#translateCodon('GGG')

extractSequence <- function(sequence, start, end, pad=FALSE)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	if (start<1)
		throw('start position is less than 1: ',start)
	if (end<1)
		throw('end position is less than 1: ',start)
	if (end>nchar(sequence))
		throw('end is beyond length of sequence: end=',end,' length=',nchar(sequence))#end <- nchar(sequence)
	sequence <- s2c(sequence)
	return(c2s(sequence[start:end]))
}
#extractSequence(refs['KT9','sequence'],3420,5312)


#
#extractSequence <- function(sequence, start, end, pad=TRUE)
#{
#	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
#	offset <- 0
#	if (start<1)
#	{
#		if (pad)
#		{
#			offset <- abs(start)+1
#			start <- 1
#		}
#		else throw('start position is less than 1: ',start)
#	}
#	if (end<1)
#	{
#		throw('end position is less than 1: ',start)
#	}
#	#if (end>nchar(sequence))
#	#	end <- nchar(sequence)
#	sequence <- s2c(sequence)
#	sequence <- c2s(sequence[start:end])
#	sequence <- concat(paste(rep('_',offset),collapse=''),sequence)
#	return(sequence)
#}
#extractSequence('gcgcctatcacagcatactcccaa',1,3)

extractCodon <- function(sequence, codon)
{
	start <- (codon-1)*3+1
	end <- start+2
	printParams(start,end)
	return(extractSequence(sequence,start,end))
}
#extractCodon('gcgcctatcacagcatactcccaa',1)
#extractCodon('gcg cct atc aca gca tac tcc caa',1)

splitCodons <- function(sequence,start=1)
{
	require(gsubfn) 
	#sequence <- cleanSequence(sequence)
	codons <- strapply(sequence, "...")[[1]]
	ntnum <- start
	data <- data.frame()
	for (codonnum in 1:length(codons))
	{
		refcodon <- codons[codonnum]
		refaa <- translateCodon(refcodon)
		data[codonnum,'ntnum'] <- ntnum
		data[codonnum,'codon'] <- codonnum
		data[codonnum,'refcodon'] <- refcodon
		data[codonnum,'refaa'] <- refaa
		ntnum <- ntnum+3
	}
	return(data)
}
#splitCodons(sequence)

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
	replicates <- config@data[which(config@data$subject==subject),'replicate']
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

