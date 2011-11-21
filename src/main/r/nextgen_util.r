# estimates the number of reads based on the size of the file, assuming a ratio of 7757 for uncompressed fastq
estimateReadCount <- function(mb)
{
	return(round(mb*7757))
}
#estimateReadCount(112.5)

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

plotReadDistributions <- function(config, outfile="histograms.pdf")
{
	pdf(outfile)
	for (sample in rownames(config@samples))
	{
		try({
			filename <- concat(config@pileup.dir,'/',sample,'.txt')
			data <- loadDataFrame(filename)
			table <- createNtCountTable(data, cutoff=0)
			runs.subset <- subset(runs, sample==smpl)
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
	par(mfrow=c(1,1))
	dev.off()
}
#plotReadDistributions(config)

####################################################

displayCodons <- function(sequence,start=1)
{
	library(gsubfn) 
	sequence <- cleanSequence(sequence)
	codons <- strapply(sequence, "...")[[1]]
	ntnum <- start
	data <- data.frame()
	for (num in 1:length(codons))
	{
		codon <- codons[num]
		aa <- translateCodon(codon)
		data[num,'ntnum'] <- ntnum
		data[num,'codon'] <- codon
		data[num,'aa'] <- aa
		ntnum <- ntnum+3
	}
	return(t(data))
}
#data <- displayCodons('gcgcctatcac agcatactcc caacagacgc ggggcttact tggctgcatc atcactagcc', start=3420)

######################################################

getPileupConsensusSequence <- function(config,sample)
{
	data <- loadPileupData(config,sample)
	data.subset <- as.data.frame(data[,c('position','nt')])
	xtab <- xtabs(~nt + position, data.subset)
	sequence <- c()
	for (col in colnames(xtab))
	{
		nt <- names(sort(xtab[,col], decreasing=TRUE)[1])
		sequence <- c(sequence,nt)
	}
	return(tolower(joinFields(sequence,'')))
}
#getPileupConsensusSequence(config,'PXB0220-0030.8__HCV-KT9')
#getPileupConsensusSequence(config,'11551793.1__HCV-NS3-36')

######################################################

#diagnostics

viewBam <- function(config, sample, alignment_status='Aligned', pf_status='All')
{
	bamfile <- concat(config@bam.dir,'/',sample,'.bam')
	checkFileExists(bamfile)
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/ViewSam.jar'
	str <- concat(str,' INPUT=',bamfile)
	str <- concat(str,' ALIGNMENT_STATUS=',alignment_status)#{Aligned, Unaligned, All}
	str <- concat(str,' PF_STATUS=',pf_status) #{PF, NonPF, All}
	runCommand(str)
}
#viewBam(config,'10464592.1__HCV-NS3-156', alignment_status='Unaligned')

#######################################################################

getRawReadCounts <- function(config)
{
	for (rowname in rownames(config@data))
	{
		row <-  config@data[rowname,]
		dir <- concat(config@tmp.dir,'/',row$sample)
		filename <- concat(row$folder,'_',row$barcode,'_L00',row$lane,'_R1_001.fastq')
		kb <- getFileSize(concat(dir,'/',filename))
		reads <- estimateReadCount(kb/1000)
		print(concat(filename,': ',reads))
		config@data[rowname,'reads'] <- reads
	}
	runcounts <- config@data[,splitFields('sample,region,reads')]
	outfile <- concat(config@tmp.dir,'/','readcounts-runs.txt')
	writeTable(runcounts,outfile,row.names=FALSE)
	
	samplecounts <- reshape::cast(runcounts, sample~., value='reads', fun.aggregate=sum)
	outfile <- concat(config@tmp.dir,'/','readcounts-samples.txt')
	writeTable(samplecounts,outfile,row.names=FALSE)
	return(list(runcounts=runcounts, samplecounts=samplecounts))
}
#getRawReadCounts(config)

getMapStats <- function(bamfile)
{
	checkFileExists(bamfile)
	total <- as.numeric(system(concat('samtools flagstat ',bamfile,' | awk \'{print $1}\' | head -n 1'), intern=TRUE))
	mapped <- as.numeric(system(concat('samtools flagstat ',bamfile,' | awk \'{print $1}\' | head -n 3 | tail -n 1'), intern=TRUE))
	prop <- mapped/total
	#print(concat('total=',total,', mapped=',mapped,' proportion=',prop))
	return(list(total=total, mapped=mapped, prop=prop))
}
#getMapStats(concat(config@bam.dir,'/','etsuko_yamada__HCV-HCJ4.bam'))
#
#getMappingReport <- function(config, bam.dir=config@bam.dir, samples=config@samples, outfile=concat(config@tmp.dir,'/','readcounts-mapped.txt'))
#{
#	counts <- data.frame()
#	for (sample in samples)
#	{
#		try({
#			stats <- getMapStats(concat(bam.dir,'/',sample,'.bam'))
#			print(concat('mapstats for ',sample,': ',stats))
#			counts[sample,'sample'] <- sample
#			counts[sample,'total'] <- stats$total
#			counts[sample,'mapped'] <- stats$mapped
#			counts[sample,'prop'] <- stats$prop
#		})
#	}
#	#outfile <- concat(out.dir,'/','readcounts-mapped.txt')
#	writeTable(counts,outfile,row.names=FALSE)
#	return(counts)
#}
##getMappingReport(config)

getMappingReport <- function(config, bam.dir=config@bam.dir, samples=config@samples, outfile=concat(config@tmp.dir,'/','readcounts-mapped.txt'))
{
	counts <- data.frame()
	for (filename in list.files(bam.dir, pattern='\\.bam$'))
	{
		try({
			stats <- getMapStats(concat(bam.dir,'/',filename))
			sample <- stripExtension(filename)
			print(concat('mapstats for ',sample,': ',stats))
			counts[sample,'sample'] <- filename
			counts[sample,'total'] <- stats$total
			counts[sample,'mapped'] <- stats$mapped
			counts[sample,'prop'] <- stats$prop
		})
	}
	writeTable(counts,outfile,row.names=FALSE)
	return(counts)
}
#getMappingReport(config)


bambino <- function(config, sample)
{
	bamfile <- concat(config@bam.dir,'/',sample,'.bam')
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	runCommand('java -Xmx3072m -cp $BAMBINO_HOME/bambino_bundle_1.03.jar Ace2.AceViewer -bam ',bamfile,' -fasta ',reffile)
}
#bambino(config,'KT9.plasmid__KT9')

#########################################################################

loadConsensusFile <- function(config, sample)
{
	filename <- concat(config@consensus.dir,'/',sample,'.consensus.fastq')
	lines <- readLines(filename)
	str <- ''
	for (line in lines[-1])
	{
		if (substring(line,1,1)=='+')
			break
		#print(line)
		str <- concat(str,line)
	}
	return(str)
}
#loadConsensusFile(config,'CTE247-21__HCV-KT9')

loadConsensusFiles <- function(config, samples=config@samples)
{
	# read each fasta file
	seqs <- list()
	for (sample in samples)
	{
		try(seqs[[sample]] <- loadConsensusFile(config,sample))
	}
	return(seqs)
}
#loadConsensusFiles(config)

makeConsensusFasta <- function(config, samples=config@samples)
{
	seqs <- loadConsensusFiles(config,samples)
	writeFastaFile(concat(config@tmp.dir,'/consensus.fasta'),seqs)
}
#makeConsensusFasta(config)

###############################################################################


extractSampleFile <- function(config,id)
{
	tmp.dir <- concat(config@check.dir,'/tmp')
	fastq.dir <- concat(config@check.dir,'/fastq')
	row <-  config@data[id,]
	folder <- row$folder
	barcode <- row$barcode
	lane <- row$lane
	dir.from <- concat(config@data.dir,'/',row$rundata,'/Unaligned/Project_',folder,'/Sample_',folder,'/')
	pattern <- concat('^',folder,'_',barcode,'_L00',lane,'_R1_.*\\.fastq.*')
	for (filename in list.files(dir.from,pattern))
	{
		ext <- ifelse(getFileExtension(filename)=='gz','.fastq.gz','.fastq')
		newfilename <- concat(tmp.dir,'/',id,ext)
		runCommand('cp ', dir.from, filename,' ',newfilename)
	}
	runCommand('gunzip ',tmp.dir,'/*')

	dir.from <- concat(tmp.dir,'/',id)
	fastqfile <- concat(fastq.dir,'/',id,'.fastq')
	runCommand('cat ',tmp.dir,'/* > ',fastqfile)
	checkFileExists(fastqfile)
	
	runCommand('rm -r ',tmp.dir,'/*')
}
#extractSampleFile(config, 'nextgen4-3D')

############################################################

checkSampleMappingForRun <- function(config, id, refs)
{
	#check a specific sample using several metrics
	dir <- config@check.dir
	tmp.dir <- concat(dir,'/tmp')
	fastq.dir <- concat(dir,'/fastq')
	bam.dir <- concat(dir,'/bam')	
	extractSampleFile(config,id)
	trimSolexaqa(config,id,fastq.dir)
	fqfile <- concat(fastq.dir,'/',id,'.fastq')
	fqtrfile <- concat(fastq.dir,'/',id,'.trimmed.fastq')
	samples <- c()
	for (ref in splitFields(refs))
	{
		if (is.na(config@refs[ref,'start']))
		{
			print(concat('skipping nonexistent ref: ',ref))
			next
		}
		print(ref)
		for (sample in c(concat(id,'__',ref)))#,concat(id,'.trimmed','__',ref)))
		{
			print(sample)
			samples <- c(samples,sample)
			reffile <- getRefFile(config,ref)			
			bwa(fqfile, reffile, bam.dir, outstem=sample)
			writeConsensusForBam(config, sample, bam.dir=bam.dir, out.dir=dir)
		}
	}
	runCommand('rm -r ',bam.dir,'/*.sam')
	runCommand('rm -r ',bam.dir,'/*.sai')
	fixBaiFiles(config,bam.dir)
	print(samples)
	return(samples)
}

checkSampleMappingForSubject <- function(config, subject, refs)
{
	makeSubDirs(config@check.dir,'tmp,fastq,bam')
	#writeRefs(config)
	
	samples <- c()
	for (id in config@data[which(config@data$subject==subject),'id'])
	{
		try({
			smpls <- checkSampleMappingForRun(config,id,refs)
			samples <- c(samples,smpls)
		})
	}
	dir <- config@check.dir
	bam.dir <- concat(dir,'/bam')
	outfile <- concat(dir,'/','readcounts-',subject,'.txt')
	report <- getMappingReport(config, bam.dir=bam.dir, samples,  outfile=outfile)
	return(report)
}

##########################################################################

writeConsensusForBam <- function(config,sample,bam.dir=config@bam.dir, out.dir=config@consensus.dir)
{
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	bamfile <- concat(bam.dir,'/',sample,'.bam')
	fastqfile <- concat(out.dir,'/',sample,'.consensus.fastq')
	runCommand('samtools mpileup -uf ',reffile,' ',bamfile,' | bcftools view -cg - | vcfutils.pl vcf2fq > ',fastqfile)
	checkFileExists(fastqfile)
	#fastq2fasta(fastqfile)
}
#writeConsensusForBam(config,'10464592.1__HCV-NS3-156')

writeCoverageForBam <- function(config,sample)
{
	bamfile <- concat(config@bam.dir,'/',sample,'.bam')
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	outfilebase <- concat(config@coverage.dir,'/',sample)
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T DepthOfCoverage'
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -o ',outfilebase)
	str <- concat(str,' -R ',reffile)
	runCommand(str)
}
#writeCoverageForBam(config,'PXB0220-0002.wk09__HCV-KT9')


###############################################

displayConfig <- function(config)
{
	data <- loadDataFrame(concat(config@config.dir,'/data.txt'), idcol='id')
	#readcounts <- loadDataFrame(concat(config@config.dir,'/readcounts.txt'), idcol='id')
	indent <- '  '
	for (group in unique(data$group))
	{
		print(concat('group: ',group))
		grouprows <- data[which(data$group==group),]
		for (subgroup in unique(grouprows$table))
		{
			print(concat(indent,'subgroup: ',subgroup))
			subgrouprows <- grouprows[which(grouprows$table==subgroup),]
			for (replicate in unique(subgrouprows$column))
			{
				replicaterows <- subgrouprows[which(subgrouprows$column==replicate),]
				stem <- unique(replicaterows$stem)
				print(stem)
				#counts <- readcounts[which(readcounts$id==stem),]
				counts <- list(mapped=0,total=0)
				print(concat(indent,indent,'replicate: ',replicate,', counts=',counts$mapped,'/',counts$total,' (',counts$prop,')'))				
				for (run in unique(replicaterows$id))
				{
					row <- replicaterows[which(replicaterows$id==run),]
					print(concat(indent,indent,indent,'run: ',run,', region=',row$region))
				}
			}
		}
	}
}
#displayConfig(config)

###########################################################################

exportUnmappedReads <- function(config,stem)
{
	bamfile <- concat(config@bam.dir,'/',stem,'.bam')
	checkFileExists(bamfile)
	unmapped.dir <- concat(config@out.dir,'/unmapped')	
	fastqfile <- concat(unmapped.dir,'/',stem,'.unmapped.fastq')
	fastafile <- concat(unmapped.dir,'/',stem,'.unmapped.fasta')
	tablefile <- concat(unmapped.dir,'/',stem,'.table.txt')
	countsfile <- concat(unmapped.dir,'/',stem,'.counts.txt')
	uniquefile <- concat(unmapped.dir,'/',stem,'.unique.txt')
	runCommand('bam2fastq --force --no-aligned --unaligned --no-filtered -o ',fastqfile,' ',bamfile)
	runCommand('perl $VARDB_RUTIL_HOME/fq_all2std.pl fq2fa ',fastqfile,' > ',fastafile)
	runCommand('perl $VARDB_RUTIL_HOME/fastq2table.pl -i ',fastqfile,' -o ',tablefile)
	runCommand('sort ',tablefile,' | uniq -dc > ',countsfile)
	runCommand('sort --numeric-sort --reverse ',countsfile,' > ',uniquefile)
	
	runCommand('rm ',fastqfile)
	runCommand('rm ',fastafile)
	runCommand('rm ',tablefile)
	runCommand('rm ',countsfile)
}
#exportUnmappedReads(config,'nextgen4-8F__HCV-KT9')



convertUniqueReadsToFasta <- function(config,sample,mincount=1000,maxdiff=5)#100000
{
	#refid <- config@data[stem,'ref']
	#ref <- getRefSequence(config,refid)
	filename <- concat(config@out.dir,'/unmapped/',sample,'.unique.txt')
	data <- read.table(filename, header=FALSE, encoding='UTF-8', comment.char='#', stringsAsFactors=FALSE)
	colnames(data) <- c('count','read')
	data <- subset(data,count>=mincount)
	data <- data[order(data$count,decreasing=TRUE),]
	head(data)
	seqs <- list()
	for (rowname in rownames(data))
	{
		count <- data[rowname,'count']
		seq <- tolower(data[rowname,'read'])
		readname <- concat('read',rowname,'_',count)
		seqs[readname] <- seq
	}
	#head(seqs)
	outfile <- concat(stripExtension(filename),'.fasta')
	writeFastaFile(outfile,seqs)
	return(list(filename=outfile, sequences=seqs))
}
#convertUniqueReadsToFasta(config,'nextgen1-2A__HCV-KT9')

analyzeUnmappedReads <- function(config,stem)
{
	refid <- config@data[stem,'ref']
	ref <- getRefSequence(config,refid)
	reffile <- getRefFile(config,refid)
	sample <- concat(stem,'__',refid)
	exportUnmappedReads(config,sample)
	res <- convertUniqueReadsToFasta(config,sample)
	infile <- res$filename
	reads <- res$sequences
	outfile <- concat(sample,'.txt')
	db <- concat(refid,'_db')
	runCommand('makeblastdb -in ',reffile,' -dbtype nucl -out ',db)
	#runCommand('blastn -db ',db,' -query ',infile,' -out ',outfile,' -word_size 9 -dust no -evalue 1000 -ungapped -max_target_seqs 1 -outfmt 7')#6
	#runCommand('blastn -db ',db,' -query ',infile,' -out ',outfile,' -word_size 7 -dust no -evalue 1000 -gapopen 5 -gapextend 3 -reward 2 -penalty -3 -max_target_seqs 1 -outfmt 7')#6
	runCommand('blastn -db ',db,' -query ',infile,' -out ',outfile,' -word_size 4 -dust no -evalue 1000 -ungapped -max_target_seqs 1 -outfmt 7')#6
	checkFileExists(outfile)
	data <- read.table(outfile, header=FALSE, encoding='UTF-8', comment.char='#', stringsAsFactors=FALSE)
	colnames(data) <- splitFields('query_id,subject_id,identity,align_length,mismatches,gap_opens,q_start,q_end,s_start,s_end,evalue,bitscore')
	rownames(data) <- data$query_id
	newseqs <- list()
	newseqs[[refid]] <- ref
	for (readname in rownames(data)) #rev(
	{
		#print(data[readname,])
		q_start <- data[readname,'q_start']
		q_end <- data[readname,'q_end']
		s_start <- data[readname,'s_start']
		s_end <- data[readname,'s_end']
		read <- reads[[readname]]
		seqname <- readname
		# assume no insertions and expand the range to include the whole read
		if (s_start>s_end) #use reverse complement
		{
			read <- reverseComplement(read)
			tmp <- s_start; s_start <- s_end; s_end <- tmp
			seqname <- concat(readname,'_rev')
		}
		s_start <- s_start-q_start+1
		s_end <- s_end+nchar(read)-q_end
		#s_start <- s_start-10
		#s_end <- s_end+10
		#oldread <- extractSequence(ref,s_start,s_end)
		#printcat('start: ',s_start,', end: ',s_end,', length: ',nchar(oldread))
		#print(oldread)
		#print(read)
		newseqs[[seqname]] <- read
	}
	outfile <- concat(stripExtension(infile),'.strand.fasta')
	writeFastaFile(outfile,newseqs)
}
#analyzeUnmappedReads(config,'nextgen1-2A')
#analyzeUnmappedReads(config,'nextgen2-5I')
