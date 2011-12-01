
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
			reportAminoAcidChangeBarChart(config, subject, region, usecounts=TRUE)
#			samples <- getSamplesForSubject(config,subject)
#			for (aanum in getFociForRegion(config,region))
#			{
#				printcat(' aanum: ',aanum)
#				tbl <- makeAminoAcidVariantTable(config, samples, region, aanum, minreads=2)[,-1]
#				#gplots::textplot(tbl, show.rownames=FALSE, show.colnames=TRUE, cex=1, halign='left', valign='top')
#				gplots::textplot(capture.output(tbl), cex=1, halign='left', valign='top') #, show.rownames=FALSE
#				title(concat(group,'-',subject,'-',region,'-',aanum))
#			}
		}
	}
	dev.off()
	#openPdf(pdffile)
}
#reportAminoAcidChanges(config, 'hcv_infection')
#reportAminoAcidChanges(config, 'NS5A_L31V_Y93H_mutations_maintained')
