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

getMappingReport <- function(config, bam.dir=config@bam.dir, samples=config@samples, outfile=concat(config@tmp.dir,'/','readcounts-mapped.txt'))
{
	counts <- data.frame()
	for (sample in samples)
	{
		try({
			stats <- getMapStats(concat(bam.dir,'/',sample,'.bam'))
			print(concat('mapstats for ',sample,': ',stats))
			counts[sample,'sample'] <- sample
			counts[sample,'total'] <- stats$total
			counts[sample,'mapped'] <- stats$mapped
			counts[sample,'prop'] <- stats$prop
		})
	}
	#outfile <- concat(out.dir,'/','readcounts-mapped.txt')
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
	dir.from <- concat('../data/',row$rundata,'/Unaligned/Project_',folder,'/Sample_',folder,'/')
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


checkSampleMapping <- function(config, id, refs)
{
	#check a specific sample using several metrics
	dir <- config@check.dir
	tmp.dir <- concat(dir,'/tmp')
	fastq.dir <- concat(dir,'/fastq')
	bam.dir <- concat(dir,'/bam')
	makeSubDirs(dir,'tmp,fastq,bam')
	extractSampleFile(config,id)
	trimSolexaqa(config,id,fastq.dir)
	fqfile <- concat(fastq.dir,'/',id,'.fastq')
	fqtrfile <- concat(fastq.dir,'/',id,'.trimmed.fastq')
	samples <- c()
	for (ref in splitFields(refs))
	{
		for (sample in c(concat(id,'__',ref),concat(id,'.trimmed','__',ref)))
		{
			samples <- c(samples,sample)
			reffile <- getRefFile(config,ref)			
			bwa(fqfile, reffile, bam.dir, outstem=sample)
			filterBam(config, sample, bam.dir=bam.dir)
			sample.filtered <- concat(sample,'.filtered')
			samples <- c(samples, sample.filtered)
			writeConsensusForBam(config, sample, bam.dir=bam.dir, out.dir=dir)
			writeConsensusForBam(config, sample.filtered, bam.dir=bam.dir, out.dir=dir)
		}
	}
	runCommand('rm -r ',bam.dir,'/*.sam')
	runCommand('rm -r ',bam.dir,'/*.sai')
	fixBaiFiles(config,bam.dir)
	print(samples)
	outfile <- concat(dir,'/','readcounts-',id,'.txt')
	report <- getMappingReport(config, bam.dir=bam.dir, samples,  outfile=outfile)
	return(report)
}
#checkSampleMapping(config,'nextgen4-3D','HCV-KT9,HCV-HCJ4')
checkSampleMapping(config,'nextgen4-3B','HCV-KT9,HCV-HCJ4')#NS5Aaa31
checkSampleMapping(config,'nextgen4-5H','HCV-KT9,HCV-HCJ4')#NS5Aaa93


