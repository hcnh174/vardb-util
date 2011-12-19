
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


reportAminoAcidChangeBarChart <- function(config, subject, region, usecounts=TRUE)
{
	if (!hasRegion(config,subject,region))
	{
		printcat('region ',region,' is not available for subject ',subject)
		return()
	}
	group <- getGroupForSubject(config, subject)
	start <- config@regions[region,'start']
	end <- config@regions[region,'end']
	samples <- getSamplesForSubject(config,subject)
	data.subset <- 	getCodonCountSubset(config, samples, region, 'aa', start, end, minreads=10)
	data.subset <- subset(data.subset, aa!='*')
	if (nrow(data.subset)==0)
	{
		printcat('group=',group,' start=',start,' end=',end,' samples=',samples)
		print(head(data.subset))
		throw('No data for subject=',subject,' in region ',region)
	}
	start <- min(as.numeric(levels(data.subset$aanum)[data.subset$aanum]))
	end <- max(as.numeric(levels(data.subset$aanum)[data.subset$aanum]))
	col <- getAminoAcidColors(data.subset$aa)	
	frmla <- freq ~ aanum | column
	main <- concat(group,': ',subject); main <- simpleCap(gsub ('_',' ',main,ignore.case=T,perl=T))
	xlab <- 'Amino acid number'; ylab <- 'Amino acid count'
	#ylim <- c(0,1)
	if (usecounts)
	{
		frmla <- count ~ aanum | column
		#ylim <- c(0,8500)
		ylab <- 'Amino acid frequency'
	}
	chrt <- barchart(frmla, data.subset, group=aa, horizontal=FALSE, stack=TRUE, box.width = 1,
			main=main, xlab=xlab, sub=region, ylab=ylab, #ylim=ylim,
			par.settings=list(axis.text=list(cex=0.7), fontsize=list(text=10), 
					superpose.polygon = list(col=col)),
			auto.key = list(space = "right"),
			strip = FALSE, strip.left = strip.custom(bg='#F2F2F2',fg='#F2F2F2', horizontal = FALSE),
			layout = c(1,length(unique(data.subset$column))))
	print(chrt)
	for (aanum in as.numeric(splitFields(config@regions[region,'focus'])))
	{
		addLine(v = aanum - start + 1 - 0.5, col='red', lty=2)
		addLine(v = aanum - start + 1 + 0.5, col='red', lty=2)
	}
	return(data.subset)
}
#data.subset <- reportAminoAcidChangeBarChart(config,'PXB0220-0002','NS5Aaa93')

reportAminoAcidChangesForGroup <- function(config, group, make.pdf=FALSE)
{
	subjects <- getSubjectsForGroup(config,group)
	if (make.pdf)
	{
		pdffile <- concat(config@charts.dir,'/charts-',group,'-aa.pdf')
		pdf(pdffile)#, width = 10, height = 3)
	}
	regions <- getRegionsForGroup(config,group)
	for (region in regions)
	{
		for (subject in subjects)
		{
			printcat('subject: ',subject,', region: ',region)
			reportAminoAcidChangeBarChart(config, subject, region, usecounts=TRUE)
		}
	}
	if (make.pdf)
		dev.off()
	#openPdf(pdffile)
}
#reportAminoAcidChangesForGroup(config, 'hcv_infection')
#reportAminoAcidChangesForGroup(config, 'NS5A_L31V_Y93H_mutations_maintained')

reportAminoAcidChanges <- function(config, groups=config@groups, single.pdf=TRUE)
{
	if (single.pdf)
	{
		pdffile <- concat(config@charts.dir,'/charts.pdf')
		pdf(pdffile)#, width = 10, height = 3)
	}
	for (group in groups)
	{
		reportAminoAcidChangesForGroup(config, group, make.pdf=!single.pdf)
	}
	if (single.pdf)
		dev.off()
}
#reportAminoAcidChanges(config)
#reportAminoAcidChanges(config, single.pdf=FALSE)


######################################################################################################
#
#reportCodonChangeForSubject <- function(config, subject, region)
#{
#	if (!hasRegion(config,subject,region))
#	{
#		printcat('region ',region,' is not available for subject ',subject)
#		return()
#	}
#	group <- getGroupForSubject(config, subject)
#	start <- config@regions[region,'start']
#	end <- config@regions[region,'end']
#	samples <- getSamplesForSubject(config,subject)
#	data.subset <- 	getCodonCountSubset(config, samples, region, 'codons', start, end, minreads=0)
#	#data.subset <- subset(data.subset, aa!='*')
#	if (nrow(data.subset)==0)
#	{
#		printcat('group=',group,' start=',start,' end=',end,' samples=',samples)
#		print(head(data.subset))
#		throw('No data for subject=',subject,' in region ',region)
#	}
#	#aanums <- as.numeric(levels(data.subset$aanum)[data.subset$aanum])
#	#start <- min(aanums)
#	#end <- max(aanums)
#	#printcat('start=',start,' end=',end)
#	
#	#indicate whether the aa is the same as the ref or not
#	ref <- getRefForSamples(config,samples)
#	data.subset$isref <- 'synonymous'
#	oldwarn <- options('warn')
#	options(warn=-1) #temporarily disable warnings
#	for (aanum in data.subset$aanum)
#	{
#		refaa <- getReferenceAminoAcid(config, ref, region, aanum)
#		data.subset[which(data.subset$aanum==aanum & data.subset$aa==refaa),'isref'] <- 'non-synonymous'
#	}
#	options(oldwarn)
#
#	main <- concat(group,': ',subject); main <- simpleCap(gsub ('_',' ',main,ignore.case=T,perl=T))
#	xlab <- 'Amino acid number'; ylab <- 'Amino acid count'
#	col <- c('lightgrey','red') 
#	chrt <- barchart(freq ~ aanum | column, data.subset, group=isref, horizontal=FALSE, stack=TRUE, box.width = 1,
#			main=main, xlab=xlab, sub=region, ylab=ylab, ylim=c(0,1),
#			par.settings=list(axis.text=list(cex=0.7), fontsize=list(text=10), superpose.polygon = list(col=col)),
#			auto.key = list(space = "right"),
#			strip = FALSE, strip.left = strip.custom(bg='#F2F2F2',fg='#F2F2F2', horizontal = FALSE),
#			layout = c(1,length(unique(data.subset$column))))
#	print(chrt)
#	for (aanum in as.numeric(splitFields(config@regions[region,'focus'])))
#	{
#		addLine(v = aanum - start + 1 - 0.5, col='orange')#, lty=2)
#		addLine(v = aanum - start + 1 + 0.5, col='orange')#, lty=2)
#	}
#	print(data.subset)
#	return(data.subset)
#}
##data.subset <- reportCodonChangeForSubject(config,'PXB0220-0030','NS5Aaa93')
#

reportCodonChangeForSubject <- function(config, subject, region)
{
	if (!hasRegion(config,subject,region))
	{
		printcat('region ',region,' is not available for subject ',subject)
		return()
	}
	group <- getGroupForSubject(config, subject)
	start <- config@regions[region,'start']
	end <- config@regions[region,'end']
	samples <- getSamplesForSubject(config,subject)
	data.subset <- 	getCodonCountSubset(config, samples, region, 'codons', start, end, minreads=0)
	#data.subset <- subset(data.subset, aa!='*')

	main <- concat(group,': ',subject); main <- simpleCap(gsub ('_',' ',main,ignore.case=T,perl=T))
	xlab <- 'Amino acid number'; ylab <- 'Amino acid count'
	col <- c('red','lightgrey','pink') 
	chrt <- barchart(count ~ aanum | column, data.subset, group=subtype, horizontal=FALSE, stack=TRUE, box.width = 1,
			main=main, xlab=xlab, sub=region, ylab=ylab, ylim=c(0,1),
			par.settings=list(axis.text=list(cex=0.7), fontsize=list(text=10), superpose.polygon = list(col=col)),
			#auto.key = list(space = "right"),
			strip = FALSE, strip.left = strip.custom(bg='#F2F2F2',fg='#F2F2F2', horizontal = FALSE),
			layout = c(1,length(unique(data.subset$column))))
	print(chrt)
	for (aanum in as.numeric(splitFields(config@regions[region,'focus'])))
	{
		addLine(v = aanum - start + 1 - 0.5, col='orange')#, lty=2)
		addLine(v = aanum - start + 1 + 0.5, col='orange')#, lty=2)
	}
	return(data.subset)
}
#data.subset <- reportCodonChangeForSubject(config,'PXB0220-0030','NS5Aaa93')


reportCodonChangesForSubjects <- function(config, subjects, regions)
{
	for (region in regions)
	{
		pdffile <- concat(config@charts.dir,'/codon-changes-',region,'.pdf')
		pdf(pdffile)
		for (subject in subjects)
		{
			reportCodonChangeForSubject(config, subject, region)
		}
		dev.off()		
	}
}
#reportCodonChangesForSubjects(config,subjects,splitFields('NS3aa36,NS3aa156,NS5Aaa31,NS5Aaa93'))

