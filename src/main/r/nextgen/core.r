loadConfig <- function(config.dir=NULL, data.dir=NULL, out.dir=NULL)
{
	if (is.null(config.dir))
		config.dir <- concat('~/nextgen/config/',getCurDir())
	if (is.null(data.dir))
		data.dir <-'~/nextgen/data/'
	if (is.null(out.dir))
		out.dir <-'out'
	config <- new('nextgenconfig',config.dir=config.dir, data.dir=data.dir, out.dir=out.dir)
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

convertConsensusToFasta <- function(config, sample)
{
	filename <- concat(config@consensus.dir,'/',sample,'.consensus.fastq')
	lines <- readLines(filename)
	str <- ''
	for (line in lines[-1])
	{
		if (substring(line,1,1)=='+')
			break
		str <- concat(str,line)
	}
	outfile <- concat(config@consensus.dir,'/',sample,'.consensus.fasta')
	writeFastaFile(outfile,str,sample)
	return(outfile)
}
#convertConsensusToFasta(config,'CTE247-21__HCV-KT9')

writeConsensusForBam <- function(config,sample,bam.dir=config@bam.dir, out.dir=config@consensus.dir)
{
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	bamfile <- concat(bam.dir,'/',sample,'.bam')
	fastqfile <- concat(out.dir,'/',sample,'.consensus.fastq')
	runCommand('samtools mpileup -uf ',reffile,' ',bamfile,' | bcftools view -cg - | vcfutils.pl vcf2fq > ',fastqfile)
	checkFileExists(fastqfile)
	convertConsensusToFasta(config,sample)
}
#writeConsensusForBam(config,'10464592.1__HCV-NS3-156')

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
	if (length(filenames)==0)
		throw2('could not find any fasta files in directory: ',dir)
	print(filenames)
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
		#if (!file.exists(reffile))
		#{
			printcat('writing ref file ',reffile)
			seq <- config@refs[ref,'sequence']
			write.fasta(s2c(seq), ref, file.out=reffile)
			checkFileExists(reffile)
			runCommand('bwa index ',reffile)
		#}
	}
}
#writeRefs(config)

########################################3

checkGroupExists <- function(config, group)
{
	if (!(group %in% config@groups))
		throw2('group "',group,'" does not exist')
}
#checkGroupExists(config,'HCV')

checkSubjectExists <- function(config, subject)
{
	if (!(subject %in% config@subjects))
		throw2('subject "',subject,'" does not exist')
}
#checkSubjectExists(config,'HCV')

checkSampleExists <- function(config, sample)
{
	if (!(sample %in% config@samples))
		throw2('sample "',sample,'" does not exist')
}
#checkSampleExists(config,'KT9.plasmid__HCV-KT9')#'notexist')


#####################################################

getIdsForGroup <- function(config, group)
{
	checkGroupExists(config,group)
	ids <- sort(unique(config@data[which(config@data$group==group),'id']))
	if (length(ids)==0)
		throw2('could not find any ids for group: ',group)
	return(ids)
}
#getIdsForGroup(config,'KT9')

getSamplesForSubject <- function(config, subject)
{
	checkSubjectExists(config,subject)
	samples <- sort(unique(config@data[which(config@data$subject==subject),'sample']))
	if (length(samples)==0)
		throw2('could not find any samples for subject: ',subject)
	return(samples)
}
#getSamplesForSubject(config,'10348001')

getSamplesForGroup <- function(config, group)
{
	checkGroupExists(config,group)
	samples <- sort(unique(config@data[which(config@data$group==group),'sample']))
	if (length(samples)==0)
		throw2('could not find any samples for group: ',group)
	return(samples)
}
#getSamplesForGroup(config,'KT9')


getSubjectsForGroup <- function(config, group)
{
	checkGroupExists(config,group)
	subjects <- sort(unique(config@data[which(config@data$group==group),'subject']))
	if (length(subjects)==0)
		throw2('could not find any subjects for group: ',group)
	return(subjects)
}
#getSubjectsForGroup(config,'KT9')

getSamplesForSubGroup <- function(config, group, subgroup)
{
	checkGroupExists(config,group)
	samples <- sort(unique(config@data[which(config@data$group==group & config@data$table==subgroup),'sample']))
	if (length(samples)==0)
		throw2('could not find any samples for group: ',group,' and subgroup ',subgroup)
	return(samples)
}
#getSamplesForSubGroup(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy')
#
getStemForSample <- function(sample)
{
	stem <- stripPath(sample)
	stem <- strsplit(stem,'__', fixed=TRUE)[[1]][1] #use the part before the delimiter (__)
	return(stem)
}
#getStemForSample('out/bam/nextgen3-2G__HCV-KT9.bam')

getRefForSample <- function(sample)
{
	ref <- strsplit(sample,'__', fixed=TRUE)[[1]][2] #use the part after the delimiter (__)
	ref <- strsplit(ref,'.', fixed=TRUE)[[1]][1] # remove any extensions
	return(ref)
}
#getRefForSample('10348001.20040315__HBV-RT.filtered')

getRefForSubject <- function(config, subject)
{
	checkSubjectExists(config,subject)
	ref <- unique(config@data[which(config@data$subject==subject),'ref'])
	if (length(ref)>1)
		throw2('multiple refs found for subject+region: ',joinFields(ref,','))
	return(ref)
}
#getRefForSubject(config,'CTE247-21')

getRefForGroup <- function(config, group)
{
	checkGroupExists(config,group)
	refs <- unique(config@data[which(config@data$group==group),'ref'])
	#print(refs)
	mapref <- unique(config@refs[which(config@refs$ref %in% refs),'mapping'])
	if (length(mapref)>1)
		throw2('multiple mapping refs found for group+region: ',joinFields(mapref,','))
	return(mapref)
}
#getRefForGroup(config,'NS3_NS5A_inhibitors')

getRefForSamples <- function(config, samples)
{
	refs <- unique(config@data[which(config@data$sample %in% samples),'ref'])
	mapref <- unique(config@refs[which(config@refs$ref %in% refs),'mapping'])
	if (length(mapref)>1)
		throw2('multiple mapping refs found for samples: ',joinFields(samples,','))
	return(mapref)
}
#getRefForSamples(config,getSamplesForSubGroup(config,'hcv_infection','hcv_infection'))


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

sortRegions <- function(config, regions)
{
	sorted <- c()
	for (region in config@regions$id)
	{
		if (region %in% regions)
			sorted <- c(sorted,region)
	}
	return(sorted)
}
#sortRegions(config,splitFields('NS3aa156,NS3aa36'))

getRegionsForSample <- function(config, sample)
{
	checkSampleExists(config,group)
	return(sortRegions(config,unique(config@data[which(config@data$sample==sample),'region'])))
}
#getRegionsForSample(config,'PXB0220-0030.8__HCV-KT9')

getRegionsForSubject <- function(config, subject)
{
	checkSubjectExists(config,group)
	return(sortRegions(config,unique(config@data[which(config@data$subject==subject),'region'])))
}
#getRegionsForSubject(config,'10348001')

getRegionsForGroup <- function(config, group)
{
	checkGroupExists(config,group)
	return(sortRegions(config,unique(config@data[which(config@data$group==group),'region'])))
}
#getRegionsForGroup(config,'G9')

getRegionsForSubGroup <- function(config, group, subgroup)
{
	checkGroupExists(config,group)
	return(sortRegions(config,unique(config@data[which(config@data$group==group & config@data$table==subgroup),'region'])))
}
#getRegionsForSubGroup(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy')

getTablesForGroup <- function(config, group)
{
	checkGroupExists(config,group)
	return(unique(config@data[which(config@data$group==group),'table']))
}
#getTablesForGroup(config,'BMS-790052_BMS-650032')

getStemsForSamples <- function(config, samples)
{
	return(unique(config@data[which(config@data$sample %in% samples),'stem']))
}
#getStemsForSamples(config, config@samples)

getGroupForSubject <- function(config, subject)
{
	checkSubjectExists(config,group)
	return(unique(config@data[which(config@data$subject==subject),'group']))
}
##getGroupForSubject(config, 'PXB0220-0002')
#
#getIntervalsForSample <- function(config, sample)
#{
#	ref <- getRefForSample(sample)
#	print(ref)
#	regions <- getRegionsForSample(config, sample)
#	if (length(regions)==0)
#	{
#		stem <- getStemForSample(sample)
#		regions <- config@data[stem,'region']
#	}
#	print(regions)
#	genes <- unique(config@regions[which(config@regions$id %in% regions),'gene'])
#	print(genes)
#	intervals <- c()
#	for (gene in genes)
#	{
#		start <- config@genes[gene,'start']
#		end <- config@genes[gene,'end']
#		interval <- concat(ref,':',start,'-',end)
#		intervals <- c(intervals,interval)
#	}
#	return(intervals)	
#}
##getIntervalsForSample(config,'KT9.random__HCV-KT9')

#getLocationForRegion <- function(config, ref, region)
#{
#	mapref <- unique(config@refs[which(config@refs$ref==ref),'mapping'])
#	#print(mapref)
#	gene <- unique(config@regions[region,'gene'])
#	#print(gene)
#	row <- config@genes[which(config@genes$gene==gene & config@genes$ref==mapref),]
#	#print(row)
#	location <- list(ref=ref, start=row$start, end=row$end)
#	return(location)	
#}
##getLocationForRegion(config,'HCV-KT9','NS3aa156')

getLocationForRegion <- function(config, ref, region)
{
	mapref <- getField(config@refs,which(config@refs$ref==ref),'mapping')#mapref <- unique(config@refs[which(config@refs$ref==ref),'mapping'])
	gene <- getField(config@regions,region,'gene')#gene <- unique(config@regions[region,'gene'])
	#printcat('mapref=',mapref,' gene=',gene)
	start <- getField(config@genes,which(config@genes$gene==gene & config@genes$ref==mapref),'start')
	end <- getField(config@genes,which(config@genes$gene==gene & config@genes$ref==mapref),'end')
	location <- list(ref=ref, start=start, end=end)
	return(location)	
}
#getLocationForRegion(config,'HCV-KT9','NS3aa156')

getIntervalForRegion <- function(config, ref, region)
{
	loc <- getLocationForRegion(config,ref,region)
	interval <- concat(loc$ref,':',loc$start,'-',loc$end)
	return(interval)	
}
#getIntervalForRegion(config,'HCV-KT9','NS3aa156')


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
#getCodonPositionsForRegion(config,'NS5Aaa31','HCV-KT9')
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

getReferenceCodon <- function(config, ref, region, aanum)
{
	gene <- config@regions[region,'gene']
	#ntstart <- config@genes[which(config@genes$gene==gene & config@genes$ref==ref),'start']
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	return(refcodon)
}

getReferenceAminoAcid <- function(config, ref, region, aanum)
{
	refcodon <- getReferenceCodon(config,ref,region,aanum)
	refaa <- translateCodon(refcodon)
	return(refaa)
}
#getReferenceAminoAcid(config,'HCV-KT9','NS3aa36',36)

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
	try(return(c2s(seqinr::translate(s2c(sequence)))))
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

extractSequence <- function(sequence, start, end)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	if (start<1)
		throw2('start position is less than 1: ',start)
	if (end<1)
		throw2('end position is less than 1: ',start)
	if (end>nchar(sequence))
		throw2('end is beyond length of sequence: end=',end,' length=',nchar(sequence))#end <- nchar(sequence)
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

reverseComplement <- function(seq)
{
	require(seqinr, quietly=TRUE, warn.conflicts=FALSE)
	return(c2s(rev(comp(s2c(seq)))))
}
#reverseComplement('aaaattttggggcccc')


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


indexBam <- function(bamfile)
{
	checkFileExists(bamfile)
	runCommand('samtools index ',bamfile)
	checkFileExists(concat(bamfile,'.bai'))
}

sam2bam <- function(samfile,bamfile)
{
	checkFileExists(samfile)
	tmpbamfile <- concat(bamfile,'.sort.bam')
	print(tmpbamfile)
	runCommand('samtools view -bS ',samfile,' > ',tmpbamfile)
	checkFileExists(tmpbamfile)
	runCommand('samtools sort ',tmpbamfile,' ',stripExtension(bamfile))
	checkFileExists(bamfile)
	indexBam(bamfile)	
	deleteFile(tmpbamfile)
	return(bamfile)
}

mergeBamFiles <- function(config, filenames, outfile)
{
	if (!config@force & file.exists(outfile))
		return(outfile)
	filenames <- splitFields(filenames)
	if (length(filenames)==0)
		throw2('no bam files to merge')
	#throws an error if there is only one file, so copy it using the new name
	if (length(filenames)==1)		
	{
		copyFile(filenames,outfile)
		copyFile(concat(filenames,'.bai'),concat(outfile,'.bai'))
		return(outfile)
	}	
	str <- concat('samtools merge')
	if (config@force)
		str <- concat(str,' -f')
	str <- concat(str,' ',outfile)
	for (filename in filenames)
	{
		checkFileExists(filename)
		str <- concat(str,' ',filename)
	}
	runCommand(str)
	checkFileExists(outfile)
	indexBam(outfile)
	return(outfile)
}
#mergeBamFiles(config, c('out/NS5Aaa31.bam','out/NS5Aaa93.bam'), 'out/NS5A.bam')
#mergeBamFiles(config, c('out/test.bam'), 'out/test-merged.bam')


trimBamToRegion <- function(config, sample, start, end=start, ref=getRefForSample(sample), bam.dir=config@bam.dir)
{
	bamfile <-concat(bam.dir,'/',sample,'.bam')
	newbamfile <-concat(bam.dir,'/',sample,'.region.bam')
	checkFileExists(bamfile)
	region <- concat(ref,':',start,'-',end)
	runCommand('samtools view -b ',bamfile,' ',region,' > ',newbamfile)
	checkFileExists(newbamfile)
	indexBam(newbamfile)
	return(newbamfile)
}
#trimBamToRegion(config,'nextgen3-2G__HCV-KT9',6348)
#trimBamToRegion(config,'nextgen3-5G__HCV-KT9',6534)

findVariants <- function(config,sample)
{
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	bamfile <- concat(config@bam.dir,'/',sample,'.bam')
	bcffile <- concat(config@vcf.dir,'/',sample,'.bcf')
	vcffile <- concat(config@vcf.dir,'/',sample,'.vcf')
	checkFileExists(bamfile)
	checkFileExists(reffile)
	runCommand('samtools mpileup -uD -f ',reffile,' ',bamfile,' | bcftools view -vcg - > ',vcffile)
	checkFileExists(vcffile)
	#runCommand('bcftools view ',bcffile,' | vcfutils.pl varFilter > ',vcffile)
	#checkFileExists(vcffile)
	return(vcffile)
}

hasRegion <- function(config, subject, region)
{
	return(nrow(config@data[which(config@data$subject==subject & config@data$region==region),])>0)
}
#hasRegion(config,'subject','NS5Aaa31')


getTrimmedExtension <- function(config, trim=config@trim)
{
	return(ifelse(trim,'.trimmed',''))
}
#getTrimmedExtension(config,TRUE)

#getMaskedExtension <- function(config, mask=config@mask)
#{
#	return(ifelse(mask,'.masked',''))
#}
#
#getDedupExtension <- function(config, dedup=config@dedup)
#{
#	return(ifelse(dedup,'.dedup',''))
#}

getFastqExtension <- function(config, trim=config@trim)#, mask=config@mask, dedup=config@dedup)
{
	return(concat(getTrimmedExtension(config,trim),'.fastq'))
	#return(concat(getTrimmedExtension(config,trim),getMaskedExtension(config,mask),getDedupExtension(config,dedup),'.fastq'))
}
#getFastqExtension(config,TRUE)

getFastqFilename <- function(config, stem, fastq.dir=config@fastq.dir, trim=config@trim)#, dedup=config@dedup)
{
	ext <- getFastqExtension(config,trim=trim)#,dedup)
	return(concat(config@fastq.dir,'/',stem,ext))
}
#getFastqFilename(config, 'nextgen3-4C')

addReadGroup <- function(config, sample, bam.dir=config@bam.dir)
{
	bamfile <- concat(bam.dir,'/',sample,'.bam')
	tmpfile <- concat(bam.dir,'/tmp-',sample,'.bam')
	str <- 'java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar'
	str <- concat(str,' INPUT=',bamfile)
	str <- concat(str,' OUTPUT=',tmpfile)
	str <- concat(str,' RGSM="',sample,'"')
	str <- concat(str,' RGLB="',sample,'"')
	str <- concat(str,' RGID="',sample,'"')
	str <- concat(str,' RGPL=illumina')
	str <- concat(str,' RGPU=barcode')
	str <- concat(str,' SORT_ORDER=coordinate')
	#str <- concat(str,' CREATE_INDEX=true')
	runCommand(str)
	checkFileExists(tmpfile)
	deleteFile(bamfile)
	renameFile(tmpfile,bamfile)
	indexBam(bamfile) #bamfile
	return(bamfile)
}

stripSubjectFromSample <- function(sample)
{
	return(strsplit(sample,'.', fixed=TRUE)[[1]][2])
}
#stripSubjectFromSample('PXB0220-0030.wk17')

