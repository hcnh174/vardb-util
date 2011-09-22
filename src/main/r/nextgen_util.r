#diagnostics

show_coverage <- function(config,stem)
{
	#tmp.dir <- config@tmp.dir
	bamfile <- concat(config@tmp.dir,'/',stem,'.bam')
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T DepthOfCoverage'
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -R ',config@reffile)
	run_command(str)
	#run_command('java -jar GenomeAnalysisTK.jar -I aln.bam -R hsRef.fa -T DepthOfCoverage -L intervals.txt -U -S SILENT')
}


solexa_qa <- function(config,sample)
{
	run_command('cd ',config@qc.dir,'; SolexaQA.pl ../fastq/',sample,'.fastq')
}
#solexa_qa(config,'KT9.plasmid')

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

########################################################

calcOffsets <- function(config,ref)
{
	#ref <- 'HCV-NS3-156'
	#sequence <- config@refs[ref,'sequence']
}

refseq <- config@refs['HCV-KT9','sequence']
seq <- config@refs['HCV-NS3-156','sequence']

library(Biostrings)

psa1 <- pairwiseAlignment(pattern = refseq, subject = seq)
startnt <- pattern(psa1)@range@start #3861
#NS3	3420	5312
startntrel <- startnt-3420
#441



