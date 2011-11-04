loadConfig <- function(dir=NULL)
{
	if (is.null(dir))
		dir <- concat('../config/',getCurDir())
	config <- new('nextgenconfig',config.dir=dir)
	try({config <- preloadCodonPositions(config)})
	return(config)
}

fastqToFasta <- function(infile,outfile=NULL)
{
	if (is.null(outfile))
		outfile <- concat(stripExtension(infile),'.fasta')	
	str <- 'fastq_to_fasta'
	str <- concat(str,' -n')#keep sequences with unknown (N) nucleotides.
	str <- concat(str,' -v')#Verbose - report number of sequences.
	str <- concat(str,' -i ',infile)
	str <- concat(str,' -o ',outfile)
	runCommand(str)
}
#fastqToFasta(concat(config@consensus.dir,'/PXB0218-0007.wk10__HCV-KT9.consensus.fastq'))

fastq2fasta <- function(infile,outfile=NULL)
{
	if (is.null(outfile))
		outfile <- concat(stripExtension(infile),'.fasta')	
	str <- concat('prinseq-lite.pl -fastq ',infile,' -out_format 1 -out_good ',outfile)
	runCommand(str)
}
#fastq2fasta(concat(config@consensus.dir,'/PXB0218-0007.wk10__HCV-KT9.consensus.fastq'))

readFastaFile <- function(filename)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	print(filename)
	seqs <- read.fasta(file=filename, as.string=TRUE, seqtyp='DNA', forceDNAtolower=TRUE)
	data <- list()
	for (id in names(seqs))
	{
		newid <- strsplit(id,'\\|')[[1]][1]
		data[[newid]] <- trim(seqs[[id]][1])
	}
	return(data)
}
#readFastaFile('../config/merged/fragments.fasta')

# reads multiple fasta files and returns results in a single list
readFastaFiles <- function(dir, pattern='*.fasta')
{
	filenames <- list.files(dir,pattern)
	#print(filenames)
	data <- list()
	for (filename in filenames)
	{
		data.file <- readFastaFile(concat(dir,'/',filename))
		#print(head(data.file))
		data <- append(data, data.file)
	}
	return(data)
}
#readFastaFiles('../config/merged','^refs.*\\.fasta')


writeFastaFile <- function(filename, seqs, ids=names(seqs))
{
	if (class(seqs)=='character')
	{
		seq <- seqs
		seqs <- list()
		seqs[[ids]] <- seq
	}
	for (id in names(seqs))
	{
		seqs[[id]] <- s2c(seqs[[id]])
	}
	write.fasta(seqs, ids, file.out=filename)
}
#writeFastaFile('','aaaccccgggttt','HCV1')
#seqs <- list(); seqs[['HCV1']] <- 'aaaccccgggttt'; writeFastaFile('',seqs)

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

getSamplesForGroup <- function(config, group)
{
	samples <- unique(config@data[which(config@data$group==group),'sample'])
	return(samples)
}
#getSamplesForGroup(config,'KT9')

getSamplesForSubGroup <- function(config, group, subgroup)
{
	samples <- unique(config@data[which(config@data$group==group & config@data$table==subgroup),'sample'])
	return(samples)
}
#getSamplesForSubGroup(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy')

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
	refs <- unique(config@data[which(config@data$group==group),'ref'])
	#print(refs)
	mapref <- unique(config@refs[which(config@refs$ref %in% refs),'mapping'])
	if (length(mapref)>1)
		throw('multiple mapping refs found for group+region: ',joinFields(mapref,','))
	return(mapref)
}
#getRefForGroup(config,'NS3_NS5A_inhibitors')

getRefFile <- function(config, ref)
{
	reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	checkFileExists(reffile)
	return(reffile)
}
#getRefFile(config,'HBV-RT')

getRefSequence <- function(config, ref)
{
	return(getField(config@refs,ref,'sequence'))
}
#getRefSequence(config,'HCV-KT9')

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

getRegionsForSubGroup <- function(config, group, subgroup)
{
	return(unique(config@data[which(config@data$group==group & config@data$table==subgroup),'region']))
}
#getRegionsForSubGroup(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy')

getTablesForGroup <- function(config, group)
{
	return(unique(config@data[which(config@data$group==group),'table']))
}
#getTablesForGroup(config,'BMS-790052_BMS-650032')

getStemsForSamples <- function(config, samples)
{
	return(unique(config@data[which(config@data$sample %in% samples),'stem']))
}
#getStemsForSamples(config, config@samples)

#############################################################################

getCodonPositionsForRegion <- function(config, region, ref)
{
	gene <- getField(config@regions,region,'gene')
	region.start <- getField(config@regions,region,'start')
	region.end <- getField(config@regions,region,'end')
	region.focus <- splitFields(getField(config@regions,region,'focus'))
	positions <- getCodonPositionsForGene(config,gene,ref)
	positions <- positions[which(positions$codon>=region.start & positions$codon<=region.end),]
	positions$focus <- ifelse(positions$codon %in% region.focus,'*','')
	return(positions)
}
#getCodonPositionsForRegion(config,'NS3aa156','HCV-KT9')

getFociForRegion <- function(config, region)
{
	return(splitFields(getField(config@regions,region,'focus')))
}
#getFociForRegion(config,'NS3aa156')


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
	#ref.start <- getField(config@refs,ref,'start')
	refseq <- getField(config@refs,ref,'sequence')
	mapref <- getField(config@refs,ref,'mapping')
	row <- config@genes[which(config@genes$gene==gene & config@genes$ref==mapref),]
	gene.start <- row$start #config@genes[id,'start']
	gene.end <- row$end #config@genes[id,'end']	
	sequence <- extractSequence(refseq, gene.start, gene.end)
	positions <- splitCodons(sequence, start=gene.start)
	return(positions)
}
#positions <- calculateCodonPositionsForGene(config,'HCV-NS3','HCV-KT9')

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
	sequence <- gsub("[0-9]+", "\\1", sequence, ignore.case=T, perl=T)
	sequence <- gsub("\\s+", "\\1", sequence, ignore.case=T, perl=T)
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

