#
#
#reportAminoAcidChange <- function(config, subject, region, log=FALSE, updown=5)
#{
#	#aanum <- config@regions[region,'focus']
#	if (!is.integer(aanum))
#		stop(concat('cannot find aafocus aanum for region: ',region,' (',aanum,')'))
#	aanum <- as.integer(aanum)
#	start <- aanum - updown
#	end <- aanum + updown
#	
#	data.subset <- getCodonCountSubset(config,subject,region,'aa', start, end)
#	
#	print(data.subset[which(data.subset$aanum==aanum),splitFields('replicate,rank,aa,count,freq')])
#	#print(head(data.subset))
#	numcol <- length(unique(data.subset$rank))
#	col <- gray(numcol:0 / numcol)
#	if (log)
#		frmla <- as.formula(log10(count) ~ aanum | replicate)
#	else frmla <- as.formula(count ~ aanum | replicate)	
#	chrt <- barchart(frmla, data.subset, group=rank, horizontal=FALSE, stack=TRUE,
#			main=subject, xlab='Amino acid number', ylab='Amino acid count', sub=region,
#			col=col, strip=FALSE, strip.left=TRUE, #strip.text = list(cex = 0.75),
#			#auto.key = list(space = "right"),
#			layout = c(1,length(unique(data.subset$replicate))))
#	print(chrt)
#	addLine(v=aanum - start + 1 - 0.5, col='red', lty=2)
#	addLine(v=aanum - start + 1 + 0.5, col='red', lty=2)
#	return(data.subset)
#}
##reportAminoAcidChange(config,'','NS3aa156')
#
#
#
#reportAminoAcidChanges <- function(config, subject=NULL, ...)
#{
#	if (is.null(subject))
#		subjects <- unique(config@samples$subject)
#	else subjects <- c(subject)
#	for (subject in subjects)
#	{
#		filename <- concat(config@out.dir,subject,'.aa.pdf')
#		pdf(filename)
#		samples <- config@samples[which(config@samples$subject==subject),'sample']
#		regions <- unique(config@runs[which(config@runs$sample %in% samples),'region'])
#		for (region in regions)
#		{
#			reportAminoAcidChange(config, subject, region, ...)
#		}
#		dev.off()
#	}
#}
##reportAminoAcidChanges(config, 'PXB0218-0007', log=FALSE)
#
##########################################################################
#
#reportCodonChange <- function(config, subject, region, log=FALSE, updown=5)
#{
#	aanum <- config@regions[region,'aafocus']
#	if (!is.integer(aanum))
#		stop(concat('cannot find aafocus aanum for region: ',region,' (',aanum,')'))
#	aanum <- as.integer(aanum)	
#	start <- aanum - updown
#	end <- aanum + updown
#	
#	data.subset <- getCodonCountSubset(config,subject,region,'codons', start, end)
#	
#	print(data.subset[which(data.subset$aanum==aanum),splitFields('replicate,rank,aa,codon,count,freq')])
#	#print(head(data.subset))
#	numcol <- length(unique(data.subset$rank))
#	col <- gray(numcol:0 / numcol)
#	if (log)
#		frmla <- as.formula(log10(count) ~ aanum | replicate)
#	else frmla <- as.formula(count ~ aanum | replicate)
#	chrt <- barchart(frmla, data.subset, group=rank, horizontal=FALSE, stack=TRUE,
#			main=subject, xlab='Codon number', ylab='Codon count', sub=region,
#			col=col, strip=FALSE, strip.left=TRUE, #strip.text = list(cex = 0.75),
#			#auto.key = list(space = "right"),
#			layout = c(1,length(unique(data.subset$replicate))))
#	print(chrt)
#	addLine(v=aanum - start + 1 - 0.5, col='red', lty=2)
#	addLine(v=aanum - start + 1 + 0.5, col='red', lty=2)
#	return(data.subset)
#}
##reportCodonChange(config,'PXB0218-0007','NS3aa156')
#
#reportCodonChanges <- function(config, subject=NULL, ...)
#{
#	if (is.null(subject))
#		subjects <- unique(config@samples$subject)
#	else subjects <- c(subject)
#	for (subject in subjects)
#	{
#		filename <- concat(config@out.dir,subject,'.codons.pdf')
#		pdf(filename)
#		samples <- config@samples[which(config@samples$subject==subject),'sample']
#		regions <- unique(config@runs[which(config@runs$sample %in% samples),'region'])
#		for (region in regions)
#		{
#			reportCodonChange(config, subject, region, ...)
#		}
#		dev.off()
#	}
#}
##reportCodonChanges(config, 'PXB0218-0007')

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

#######################################################################3

#makeAminoAcidBarchart <- function(config, sample, region, ...)
#{
#	filename <- getCodonCountFilename(config,sample,'aa')
#	data <- loadDataFrame(filename)
#	data <- data[which(data$region==region),]
#	data$aanum <- factor(data$aanum)
#	chrt <- barchart(count ~ aanum, data, group=aa, horizontal=FALSE, stack=TRUE,
#			main=sample, xlab=concat(region,' region'), ylab='Count',
#			auto.key = list(space = "right"))
#	print(chrt)
#}
##makeAminoAcidBarchart(config,'')

makeAminoAcidBarchart <- function(config, sample, region, ...)
{
	filename <- getCodonCountFilename(config,sample,'aa')
	data <- loadDataFrame(filename)
	data <- data[which(data$region==region),]
	data$aanum <- factor(data$aanum)
	chrt <- barchart(count ~ aanum, data, group=aa, horizontal=FALSE, stack=TRUE,
			main=sample, xlab=concat(region,' region'), ylab='Count',
			auto.key = list(space = "right")
	)
	print(chrt)
}
#makeAminoAcidBarchart(config,'PXB0220-0002.wk08__HCV-KT9_PXB0220-0002', 'NS5Aaa93')

makeAminoAcidBarcharts <- function(config, group, ...)
{
	samples <- getSamplesForGroup(config,group)
	pdffile <- concat(config@charts.dir,'/barcharts-',group,'-aa.pdf')
	pdf(pdffile, width = 10, height = 3)
	regions <- getRegionsForGroup(config,group)
	for (region in regions)
	{
		for (sample in samples)
		{
			makeAminoAcidBarchart(config, sample, region, ...)
		}
	}
	dev.off()
	#openPdf(pdffile)
}
#makeAminoAcidBarcharts(config, 'NS5A_L31V_Y93H_mutations_maintained')

#########################################################

reportAminoAcidChangeBarChart <- function(config, subject, region, usecounts=TRUE)
{
	group <- getGroupForSubject(config, subject)
	start <- config@regions[region,'start']
	end <- config@regions[region,'end']
	samples <- getSamplesForSubject(config,subject)
	data.subset <- 	getCodonCountSubset(config, samples, region, 'aa', start, end, minreads=10)

	#update the start and end values
	start <- min(as.numeric(levels(data.subset$aanum)[data.subset$aanum]))
	end <- max(as.numeric(levels(data.subset$aanum)[data.subset$aanum]))
	
	frmla <- freq ~ aanum | column
	main <- concat(group,': ',subject); main <- gsub ('_',' ',main,ignore.case=T,perl=T)
	xlab <- 'Amino acid number'
	ylab <- 'Amino acid count'
	ylim <- c(0,1)
	if (usecounts)
	{
		frmla <- count ~ aanum | column
		ylim <- c(0,max(data.subset$count))
		ylab <- 'Amino acid frequency'
	}
	chrt <- barchart(frmla, data.subset, group=aa, horizontal=FALSE, stack=TRUE,
			main=main, xlab='Amino acid number', sub=region, ylab=ylab,# ylim=ylim,
			par.settings=list(axis.text=list(cex=0.7), fontsize=list(text=10)),
			auto.key = list(space = "right"),
			strip = FALSE,
			strip.left = strip.custom(bg='#F2F2F2',fg='#F2F2F2', horizontal = FALSE),#,style = 4, 
			layout = c(1,length(unique(data.subset$column))))
	print(chrt)#, more=TRUE)
	for (aanum in as.numeric(splitFields(config@regions[region,'focus'])))
	{
		addLine(v = aanum - start + 1 - 0.5, col='red', lty=2)
		addLine(v = aanum - start + 1 + 0.5, col='red', lty=2)
	}
	return(data.subset)
}
#data.subset <- reportAminoAcidChangeBarChart(config,'PXB0220-0002','NS5Aaa93')
#data.subset <- reportAminoAcidChangeBarChart(config,'junko_ogawa','NS3aa36',TRUE)


reportAminoAcidChangeXyPlot <- function(config, subject, region, minfreq=0.01, minreads=10)
{
	start <- config@regions[region,'start']
	end <- config@regions[region,'end']
	samples <- getSamplesForSubject(config,subject)
	data.subset <- 	getCodonCountSubset(config, samples, region, 'aa', start, end, minreads=minreads)
	data.subset$aanum <- as.numeric(levels(data.subset$aanum)[data.subset$aanum])
	data.subset <- subset(data.subset, aa!='*')
	#data.subset <- data.subset[order(data.subset$sample,data.subset$aanum,data.subset$rank),]
	data.subset <- subset(data.subset, rank==2)
	print(xyplot(freq ~ aanum | column, data.subset,
			main=subject, xlab='Amino acid number', ylab='Frequency of minor variant', sub=region,
			strip=FALSE, strip.left=TRUE, type='h', xlim=c(start,end), ylim=c(0,1),
			par.settings=list(axis.text=list(cex=0.6), fontsize=list(text=8)),
			scales = list(x = list(at = start:end)),
			layout = c(1,length(unique(data.subset$column)))))

	addLine(h = minfreq, col='gray', lty=2)
	for (aanum in as.numeric(splitFields(config@regions[region,'focus'])))
	{
		addLine(v = aanum - 0.5, col='red', lty=2)
		addLine(v = aanum + 0.5, col='red', lty=2)
	}
}
#reportAminoAcidChangeXyPlot(config,'PXB0220-0002','NS5Aaa93')


reportAminoAcidChanges <- function(config, group)
{
	subjects <- getSubjectsForGroup(config,group)
	pdffile <- concat(config@charts.dir,'/multibarcharts-',group,'-aa.pdf')
	pdf(pdffile)#, width = 10, height = 3)
	regions <- getRegionsForGroup(config,group)
	for (region in regions)
	{
		for (subject in subjects)
		{
			printcat('subject: ',subject,', region: ',region)
			try(reportAminoAcidChangeBarChart(config, subject, region, usecounts=TRUE))
			#try(reportAminoAcidChangeBarChart(config, subject, region, usecounts=FALSE))
			#try(reportAminoAcidChangeXyPlot(config, subject, region))
		}
	}
	dev.off()
	#openPdf(pdffile)
}
#reportAminoAcidChanges(config, 'NS5A_L31V_Y93H_mutations_maintained')
