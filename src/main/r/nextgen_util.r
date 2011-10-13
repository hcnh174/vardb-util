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
#data <- displayCodons('g
#				3421 cgcctatcac agcatactcc caacagacgc ggggcttact tggctgcatc atcactagcc
#				3481 ttacgggccg ggacaagaac caggtcgagg gagaggttca aatagtctcc accgcaacac
#				3541 aaaccttcct ggcaacctgc gtcaacggcg tgtgctggac tgtctttcac ggcgccggct
#				3601 cgaagaccct agctggccca aagggtccca tcacccaaat gtacaccaat gtagaccaag
#				3661 atcttgttgg ctggcaggcg ccccctggag cgcgctccat gacgccatgc acctgcggca
#				3721 gctcggacct ctacttggtc acgagacatg ctgatgtcat cccggtgcgc cggcggggag
#				3781 acagtagggg gagcctgctc tcccccaggc ccgtctccta cctgaagggc tcttcgggtg
#				3841 gcccactgct ctgcccttcg gggcacgttg tgggcatctt ccgggccgct gtatgcaccc
#				3901 ggggggtcgc aaaagcggtg gacttcgtac ccgttgagtc tatggaaact acaatgcggt
#				3961 ctccggtctt cacagataac tcatcccccc cggccgtacc gcagacattc caagtggcac
#				4021 atctacacgc ccccactggc agcggcaaga gtactaaagt gccagctgca tacgcagccc
#				4081 aagggtacaa ggtgctcgtc ctgaacccgt ccgttgccgc caccttaggg tttggagcgt
#				4141 acatgtccaa ggcacatggt gtagacccta acatcagaac tggggtaagg accatcacca
#				4201 cgggcgcccc catcacgtac tccacctacg gcaagttcct cgccgacggt ggttgctctg
#				4261 ggggcgccta tgatatcata atatgtgatg agtgccactc aactgactcg actaccatct
#				4321 tgggcattgg cacagttctg gaccaagcgg agacggctgg agcgcgactc gtcgtgctcg
#				4381 ccaccgctac gcctccagga tcagtcaccg tgccacaccc taatattgag gaggtggccc
#				4441 tgtccaccac tggagagatc cccttctatg gcaaggccat ccccattgag gccatcaagg
#				4501 gggggaggca tctcattttc tgccattcaa aaaaaaagtg tgatgagctc gccgcaaagc
#				4561 tgtcaaacct cggaatcaac gctgtagcgt attaccgggg tctcgatgtg tccgtcatac
#				4621 caactggcgg ggacgtcgtt gtcgtggcaa cagacgcttt aatgacgggc tttaccggcg
#				4681 actttgactc agtgatcgac tgtaacacgt gtgtcaccca aacagtcgat ttcagcttgg
#				4741 atcccacctt caccattgag acgacgaccg tgccccaaga cgcggtgtcg cgctcgcagc
#				4801 ggcggggtag gactggtaga ggtaggagag gcatctacag gtttgtgact ccaggagaac
#				4861 ggccctcggg catgttcgat tcctcggtcc tgtgtgagtg ctatgacgcg ggctgtgctt
#				4921 ggtacgagct cacgcctgct gaaacctcgg ttaggttacg ggcttaccta aatacaccag
#				4981 ggttgcccgt ttgccaggac catctggagt tctgggagag cgtcttcaca ggcctcaccc
#				5041 atatagatgc ccatttccta tcccagacca agcaggcagg agataacttc ccctatctgg
#				5101 tagcatacca ggctacagtg tgcgccaggg cccaagctcc acctccatca tgggatcaaa
#				5161 tgtggaagtg tctcatacgg ctgaaaccta cactgcacgg gcagacgccc ctgctgtata
#				5221 ggctaggagc cgttcaaaat gaggtcaccc tcacacaccc tataaccaaa tacatcatgg
#				5281 catgcatgtc ggctgacctg gaggtcgtca c', start=3420)

		########################################################
		
#		refseq <- config@refs['HCV-HCJ4','sequence']
#		seq <- config@refs['HCV-NS3-156','sequence']
#		
#		library(Biostrings)
#		
#		psa1 <- pairwiseAlignment(pattern = refseq, subject = seq, type='local')
#		startnt <- pattern(psa1)@range@start #3861
##NS3	3420	5312

#refseq <- config@refs['HCV-KT9','sequence']
#seq <- config@refs['HCV-KT9-NS3','sequence']

findFragmentStartPosition <- function(config, refseqid, seqid)
{
	refseq <-  getField(config@refs,refseqid,'sequence')
	seq <-  getField(config@refs,seqid,'sequence')
	psa1 <- Biostrings::pairwiseAlignment(pattern = refseq, subject = seq, type='local', gapOpening = -1000000)
	print(psa1)
	startnt <- Biostrings::pattern(psa1)@range@start
	#return(Biostrings::pattern(psa1)@range)
	endnt <- startnt + nchar(seq) -1 
	print(concat(startnt,'..',endnt))
	return(startnt)
}
#findFragmentStartPosition(config,'HCV-HCJ4','HCV-KT9-NS3') #3410..5302
#findFragmentStartPosition(config,'HCV-HCJ4','HCV-KT9-NS5A') #6248..7588


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

solexa_qa <- function(config,sample)
{
	run_command('cd ',config@qc.dir,'; SolexaQA.pl ../fastq/',sample,'.fastq -sanger')
}
#solexa_qa(config,'KT9.plasmid__KT9')

viewBam <- function(config, sample, alignment_status='Aligned', pf_status='All')
{
	bamfile <- concat(config@bam.dir,'/',sample,'.bam')
	checkFileExists(bamfile)
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/ViewSam.jar'
	str <- concat(str,' INPUT=',bamfile)
	str <- concat(str,' ALIGNMENT_STATUS=',alignment_status)#{Aligned, Unaligned, All}
	str <- concat(str,' PF_STATUS=',pf_status) #{PF, NonPF, All}
	run_command(str)
}
#viewBam(config,'10464592.1__HCV-NS3-156', alignment_status='Unaligned')

consensus <- function(config,sample)
{
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	bamfile <- concat(config@bam.dir,'/',sample,'.bam')	
	outfile <- concat(config@tmp.dir,'/',sample,'.consensus.txt')
	run_command('samtools mpileup -uf ',reffile,' ',bamfile,' | bcftools view -cg - | vcfutils.pl vcf2fq > ',outfile)
	checkFileExists(outfile)
}
#consensus(config,'10464592.1__HCV-NS3-156')

