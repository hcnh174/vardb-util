#install.packages('GenABEL', dep=TRUE)
#library(nlme)
#library(LDheatmap)
#library(genetics)
library(GenABEL)

qualityControl <- function(data0, maf=0.05, make.plot=FALSE)
{
	ids <- data0@phdata[which(data0@phdata$ibs==0),]$id
	qc <- check.marker(data0, p.level=0, imphetasmissing=F, ibs.mrk=0, fdr=0.2, maf=maf, idsubset=ids)
	print(summary(qc))
	print(paste('is male:',qc$ismale))
	print(paste('is female:',qc$isfemale))
	
	if (make.plot==TRUE)
		plot(qc)
	
	#make a subset of the data removing flagged ids and snps
	data1 <- data0[qc$idok, qc$snpok]
	
	#recheck HWE inflation 
	descriptives.marker(data1)[2]
	return(data1)
}


analyzeSnps <- function(model, data, data.summary, gtmode='additive', top=10, trait.type='guess', sort='Pc1df', show.plot=T)
{
	response <- as.character(model[2])
	print(check.trait(response,data=data, graph=FALSE))
	data.qt <- mlreg(model, data, trait.type=trait.type, gtmode=gtmode)
	topsnps <- descriptives.scan(data.qt, top=top, sort=sort)
	snptable <- extractSnpTable(topsnps,data.summary)
	snptable <- addTraitSummaryData(data, snptable, response)
	print(lambda(data.qt))
	print(snptable)
	if (show.plot)
	{
		oldpar <- par(mfrow=c(1,2)); on.exit(par(oldpar))
		createManhattanPlot(data.qt)
		createQQPlot(data.qt)
		#par(mfrow=c(1,1))
	}
	return(list(qt=data.qt, topsnps=topsnps, snptable=snptable))
}
#results <- analyzeSnps(plt01 ~ 1, data1maf5, data.summary, show.plot=T)

extractSnpTable <- function(topsnps, data.summary) #, cutoff=0.05)
{
	snptable <- data.frame(row.names=row.names(topsnps),
			Chromosome=topsnps$Chromosome,
			Position=topsnps$Position,
			N=topsnps$N,
			P=topsnps$P1df,
			Padj=topsnps$Pc1df,
			stringsAsFactors=FALSE)
	
	for (snp in row.names(snptable))
	{
		snptable[snp,"Alleles"] <- paste(data.summary[snp,c('A1','A2')],collapse='')
		snptable[snp,"MAF"] <- data.summary[snp,"Q.2"]
		snptable[snp,"AA"] <- data.summary[snp,"P.11"]
		snptable[snp,"AB"] <- data.summary[snp,"P.12"]
		snptable[snp,"BB"] <- data.summary[snp,"P.22"]
		snptable[snp,"HW"] <- data.summary[snp,"Pexact"]
	}
	#snptable <- subset(snptable, P<cutoff)
	return(snptable)
}

addTraitSummaryData <- function(data, snptable, response)
{
	snps <- rownames(snptable)
	phdata <- addSnpsToPhdata(data,snps)
	for (snp in snps)
	{
		phdata.sub <- subsetNA(phdata,c(response,snp))
		temp <- by(phdata.sub[[response]], phdata.sub[[snp]], function(x)
		{
			c(n=length(x), mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE))
		})
		
		temp11 <- temp[1][[1]]
		temp12 <- temp[2][[1]]
		temp22 <- temp[3][[1]]
		
		snptable[snp,"n11"] <- temp11[['n']]
		snptable[snp,"n12"] <- temp12[['n']]
		snptable[snp,"n22"] <- temp22[['n']]
		
		snptable[snp,"ave11"] <- temp11[['mean']]
		snptable[snp,"ave12"] <- temp12[['mean']]
		snptable[snp,"ave22"] <- temp22[['mean']]
		
		snptable[snp,"sd11"] <- temp11[['sd']]
		snptable[snp,"sd12"] <- temp12[['sd']]
		snptable[snp,"sd22"] <- temp22[['sd']]
		
		# add Pearson correlation
		try({
			snpnum <- paste(snp,'num',sep='')
			cor.res <- cor.test(as.formula(paste('~',response,'+',snpnum)), phdata)
			#print(cor.res)
			snptable[snp,'Pcor'] <- cor.res$p.value
		}, silent=FALSE)
	}
	print(snptable)
	return(snptable)
}

createQQPlot <- function(data.qt)
{
	print(estlambda(data.qt@results$Pc1df))
}
#createQQPlot(results$qt)

createManhattanPlot <- function(data.qt, df="Pc1df")
{
	plot(data.qt, df=df)
	#add.plot(an, df="Pc1df", col=c("lightblue","lightgreen"))
}

showTopSnps <- function(an)
{
	descriptives.scan(an0, sort="Pc1df")
}

addSnpsToPhdata <- function(data,snps)
{
	phdata <- data@phdata
	for (snp in splitFields(snps))
	{
		if (is.null(phdata[[snp]]))
		{
			print(snp)
			phdata[[snp]] <- factor(sapply(as.character(data[,snp]),function(value){sub('/', '', value)}))
			phdata[[paste(snp,'num',sep='')]] <- as.numeric(data[,snp])
		}
	}
	return(phdata)
}

#library(genetics)
calculateLinkageDisequilibrium <- function(data,snps)
{
	# r2s using package genetics
	ld <- LD(as.genotype(data[,splitFields(snps)])) #$"R^2"
	print(ld)
	print(ld$"R^2")
	return(ld)
}
#calculateLinkageDisequilibrium(data0,c("rs8099917", "rs12979860"))

makeLinkageMap <- function(data, snp, numbases=1000, snp.names=snp)
{
	position <- map(data)[snp]
	snps <- snpnames(data)[(map(data) > position-numbases & map(data) < position+numbases)]
	data.subset <- data[, snps]
	genotypes <- as.genotype(data.subset@gtdata)
	positions <- map(data.subset)
	if (snp.names[1]=='all')
		snp.names <- snps
	heatmap <- LDheatmap(genotypes, genetic.distances = positions, SNP.name = snp.names)#, color = grey.colors(20))	
}
#makeLinkageMap(data4,'rs8099917',1000)# c('rs8099917','rs12979860'))

makeSnpBoxplots <- function(data, response, snp, predictors='')
{
	oldpar <- par(ask=T); on.exit(par(oldpar))
	data@phdata <- addSnpsToPhdata(data,snps)
	for (snp in splitFields(snps))
	{
		makeSnpBoxplot(data,response,snp,predictors)
	}
	#par(ask=F)
}

makeSnpBoxplot <- function(data, response, snp, predictors='')
{
	phdata <- addSnpsToPhdata(data,snp)
	snpnum <- paste(snp,'num',sep='')
	model <- paste(response,'~',paste(c(snpnum,splitFields(predictors)),collapse=' + '))
	fit <- lm(model, phdata)
	fit.sum <- summary(fit)
	print(fit.sum)
	p.value <- fit.sum$coefficients[snpnum,'Pr(>|t|)']
	boxplot(as.formula(paste(response,'~',snp)), phdata, xlab=snp, ylab=response, main=model, sub=p.value)
}
#makeSnpBoxplot(data1maf5, 'hbdiff', 'rs6051639')

displayMaf <- function(data)
{
	sumgt <- summary(data@gtdata)
	afr <- sumgt[,"Q.2"]
	maf <- pmin(afr,(1-afr))
	oldpar <- par(mfcol=c(2,1)); on.exit(par(oldpar))
	hist(afr)
	hist(maf)
	#par(oldpar)
	
	print(catable(afr,c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99)))
	print(catable(maf,c(0,0.01, 0.05, 0.1, 0.2), cum=T))
}

##############################################################################3

plotIbs <- function(data1)
{
	#detecting genetic substructure
	data1.gkin <- ibs(data1[, data1@gtdata@chromosome != "X"], weight="freq")
	
	#create a distance matrix
	data1.dist <- as.dist(0.5 - data1.gkin)
	data1.mds <- cmdscale(data1.dist)
	plot(data1.mds)
	return(data1.mds)
}
#data1maf5.mds <- plotIbs(data1maf5)

extractCluster <- function(data1, data1.mds=NULL)
{
	if (is.null(data1.mds))
		data1.mds <- plotIbs(data1)
	
	cat('\n','Enter the number of clusters to create based on the image','\n')
	numclusters <- scan(n=1) 
	
	# just return the original data set if no value is entered
	if (length(numclusters)==0 | numclusters==1)
	{
		return(data1)
	}
	
	#identify points by cluster
	km <- kmeans(data1.mds, centers=numclusters, nstart=1000)
	for (clusternum in 1:numclusters)
	{
		cluster <- names(which(km$cluster==clusternum))
		print(paste('cluster ',clusternum,': ',length(cluster), separator=''))
	}
	
	cat('\n','Enter the number(s) of the clusters to use','\n')
	useclusters <- scan()
	
	ids <- c()
	for (clusternum in useclusters)
	{
		ids <- c(ids,names(which(km$cluster==clusternum)))
	}
	data2 <- data1[ids, ]
	return(data2)
}
#data2maf5 <- extractCluster(data1maf5, data1maf5.mds)

##########################################################################

#
#qualityControlStep1 <- function(data0)
#{
#	#quality control steps
#	#first round, don't check for HWE errors
#	qc1 <- check.marker(data0, p.level=0, imphetasmissing=F); #ibs.mrk=0
#	print(summary(qc1))
#	print(paste('is male:',qc1$ismale))
#	print(paste('is female:',qc1$isfemale))
#	
#	#make a subset of the data removing flagged ids and snps
#	data1 <- data0[qc1$idok, qc1$snpok]
#	#data1 <- Xfix(data1)
#	
#	#recheck HWE inflation 
#	descriptives.marker(data1)[2]
#	return(data1)
#}
#
#qualityControlStep2 <- function(data1, cutoff=0.5)
#{
#	#detecting genetic substructure
#	data1.gkin <- ibs(data1[, data1@gtdata@chromosome != "X"], weight="freq")
#	#data1.gkin[1:5, 1:5]
#	
#	#create a distance matrix
#	data1.dist <- as.dist(0.5 - data1.gkin)
#	data1.mds <- cmdscale(data1.dist)
#	plot(data1.mds)
#	
#	cat('\n','Enter the number of clusters to create based on the image','\n')
#	numclusters <- scan(n=1) 
#	
#	# just return the original data set if no value is entered
#	if (length(numclusters)==0 | numclusters==1)
#	{
#		return(data1)
#	}
#	
#	#identify points by cluster
#	km <- kmeans(data1.mds, centers=numclusters, nstart=1000)
#	for (clusternum in 1:numclusters)
#	{
#		cluster <- names(which(km$cluster==clusternum))
#		print(paste('cluster ',clusternum,': ',length(cluster), separator=''))
#	}
#	
#	cat('\n','Enter the numbers of the clusters to use','\n')
#	useclusters <- scan()
#	
#	ids <- c()
#	for (clusternum in useclusters)
#	{
#		ids <- c(ids,names(which(km$cluster==clusternum)))
#	}
#	data2 <- data1[ids, ]
#	return(data2)
#}
#
#qualityControlStep3 <- function(data2)
#{
#	#repeat QC with HWE checks
#	qc2 <- check.marker(data2, fdr=0.2, imphetasmissing=F)
#	summary(qc2)
#	plot(qc2)
#	data3 <- data2[qc2$idok, qc2$snpok]
#	return(data3)
#}
#
#filterByMaf <- function(data4, maf=0.05)
#{
#	qc1 <- check.marker(data4, p.level=0, imphetasmissing=F, maf=maf);
#	print(summary(qc1))
#	#make a subset of the data removing flagged ids and snps
#	data5 <- data4[qc1$idok, qc1$snpok]
#	return(data5)
#}
#
#analyzeSnps <- function(model, data, data.summary, times=1, top=10, trait.type='guess', sort='Pc1df', show.plot=T, snps=NULL)
#{
#	data.qt <- NA
#	if (!is.null(snps))
#		data.qt <- qtscore(model, data, times=times, trait.type=trait.type, snps=snps, quiet=T)
#	else data.qt <- qtscore(model, data, times=times, trait.type=trait.type, quiet=T)
#	topsnps <- descriptives.scan(data.qt, top=top, sort=sort)
#	snptable <- extractSnpTable(topsnps,data.summary)#,p.value.col=sort)
#	#print(snptable)
#	if (show.plot)
#		plot(data.qt, df="Pc1df")
#	return(snptable)
#}
#
#analyzeSnps2 <- function(model, data, data.summary, gtmode='additive', top=10, trait.type='guess', sort='Pc1df', show.plot=T)
#{
#	data.qt <- mlreg(model, data, trait.type=trait.type, gtmode=gtmode)
#	topsnps <- descriptives.scan(data.qt, top=top, sort=sort)
#	snptable <- extractSnpTable(topsnps,data.summary)#,p.value.col=sort)
#	#print(snptable)
#	if (show.plot)
#		plot(data.qt, df="Pc1df")
#	return(snptable)
#}
#
#analyzeEmpSnps <- function(model, data, data.summary, times=50, top=20, trait.type='guess', sort='Pc1df', show.plot=T, gtmode='additive')
#{	
#	op <- options()
#	options(warn=-1)
#	data.qt <- mlreg(model, data, trait.type=trait.type, gtmode=gtmode)
#	topsnps <- descriptives.scan(data.qt, top=top, sort=sort)
#	snptable <- extractSnpTable(topsnps,data.summary)#,p.value.col=sort)
#	
#	if (show.plot)
#	{
#		#par(ask=T)
#		plot(data.qt, df="Pc1df")
#		#createQQPlot(data.qt)
#		#par(ask=F)
#	}
#	
#	data.qt <- emp.qtscore(model, data, times=times, trait.type=trait.type, quiet=F)
#	topsnps <- descriptives.scan(data.qt, top=top, sort=sort)
#	snptable2 <- extractSnpTable(topsnps,data.summary,cutoff=1)#,p.value.col=sort)
#	
#	snptable$Pemp <- 1
#	for (snp in row.names(snptable))
#	{
#		snptable[snp,'Pemp'] <- snptable2[snp,'P']
#	}
#	
#	snptable <- subset(snptable, !is.na(Pemp))
#	
#	#print(snptable)
#	options(op)
#	return(snptable)
#}
##analyzeEmpSnps(plt04 ~ age + sex, data0, data.summary, times=10)

#
#plotIndividualTrendLines <- function(data.stacked, iids, color, weeks, by.individual=T, by.group=F)
#{
#	for (iid in iids)
#	{
#		#print(paste("iid=",iid))
#		individual <- data.stacked[which(data.stacked$iid==iid),]
#		#if (nrow(individual)==numweeks)
#		#{
#			#print(paste(iid,genotype))
#			if (by.individual)
#				lines(x=individual$week, y=individual$response, col=color)
#		#}
#	}
#	#gray is too hard to see for a best fit line
#	if (color=='gray')
#		color='darkgray'
#	if (by.group)
#	{
#		#try(abline(lm(response ~ week, data=subset(data.stacked,iid %in% iids)), col=color, lwd=2), silent=TRUE)
#		#try(abline(lm(response ~ as.numeric(levels(week)), data=subset(data.stacked,iid %in% iids)), col=color, lwd=2), silent=TRUE)
#		try(abline(lm(response ~ as.numeric(week), data=subset(data.stacked, iid %in% iids)), col=color, lwd=2), silent=TRUE)
#		#lines(lowess(data024$week, data024$response), col=color, lwd=2, lty=2)
#	}
#}
##plotIndividualTrendData(data5@phdata,'taqman','rs7079002', weeks='0,4', by.group=T, by.individual=F)
#
#stackTrendData <- function(phdata,responsetype, weeks='1,2,4', snp='id', omit.na=F, transform=F)
#{
#	data.stacked <- data.frame(iid=phdata$id, week=0, response=phdata[[paste(responsetype,"0",sep="")]], snp=phdata[[snp]], stringsAsFactors=F)
#	for (week in splitFields(weeks))
#	{
#		if (week!='0')
#		{
#			data.week <- data.frame(iid=phdata$id, week=week, response=phdata[[paste(responsetype,week,sep="")]], snp=phdata[[snp]])
#			data.stacked <- rbind(data.stacked,data.week)
#		}
#	}
#	#data.stacked <- data.stacked[!is.na(data.stacked$response),]
#	#data.stacked$snp <- data.stacked[[snp]]
#	if (omit.na)
#		data.stacked <- na.omit(data.stacked)
#	if (transform)
#	{
#		power <- powerTransform(data.stacked$response)
#		data.stacked$response <- data.stacked$response^power$start
#	}
#	return(data.stacked)
#}
##stackTrendData(phdata,'plt')
#
#plotIndividualTrendData <- function(phdata,responsetype,snp,weeks='0,2,4', by.individual=T, by.group=F)
#{
#	data.stacked <- stackTrendData(phdata,responsetype,weeks)
#	
#	plot(data.stacked$response ~ data.stacked$week, main=paste(responsetype,snp), ylab=responsetype, xlab="week")
#	
#	colors <- c('gray','blue','red','green')
#	colorindex <- 1
#	
#	genotypes <- phdata[[snp]]
#	names(genotypes) <- row.names(phdata)
#	uniquegenotypes <- unique(genotypes)
#	uniquegenotypecounts <- data.frame(genotype=uniquegenotypes, count=0)
#	names(uniquegenotypecounts) <- c('genotype','count')
#	
#	for (genotype in uniquegenotypes)
#	{
#		uniquegenotypecounts[which(uniquegenotypecounts$genotype==genotype),'count'] <- length(which(genotypes==genotype))
#	}
#	
#	uniquegenotypecounts <- uniquegenotypecounts[order(-uniquegenotypecounts$count),]
#	
#	for (genotype in uniquegenotypecounts$genotype)
#	{
#		color <- colors[colorindex]
#		iids <- names(genotypes[genotypes==genotype])
#		plotIndividualTrendLines(data.stacked,iids,color,weeks,by.individual,by.group)
#		colorindex <- colorindex+1
#	}
#	#xy=locator(1)
#	legend('topleft', legend=uniquegenotypes, fill=colors[1:length(uniquegenotypes)])
#}
##plotIndividualTrendData(phdata,'plt','rs4811002')
#
#plotIndividualTrendDataByVariable <- function(phdata, responsetype, variable, weeks='0,2,4', by.individual=T, by.group=T, colors=c('gray','blue','red','green'))
#{
#	data.stacked <- stackTrendData(phdata,responsetype,weeks)
#
#	plot(data.stacked$response ~ data.stacked$week, main=paste(responsetype,variable), ylab=responsetype, xlab="week")
#	
#	#colors <- c('gray','blue','red','green')
#	colorindex <- 1
#	
#	variables <- phdata[,variable]
#	uniquevariables <- unique(phdata[[variable]])
#	uniquevariablecounts <- data.frame(variable=uniquevariables, count=0)
#	names(uniquevariablecounts) <- c('variable','count')
#	
#	print(uniquevariables)
#	print(uniquevariablecounts)
#	
#	for (uniquevariable in uniquevariables)
#	{
#		count <- length(which(variables==uniquevariable))
#		uniquevariablecounts[which(uniquevariablecounts$variable==uniquevariable),'count'] <- count
#	}
#	
#	print(uniquevariablecounts)
#	uniquevariablecounts <- uniquevariablecounts[order(-uniquevariablecounts$count),]
#	print(uniquevariablecounts)
#	
#	for (uniquevariable in uniquevariablecounts$variable)
#	{
#		color <- colors[colorindex]
#		iids=row.names(phdata[which(phdata[[variable]]==uniquevariable),])	
#		plotIndividualTrendLines(data.stacked,iids,color,numweeks,by.individual,by.group)
#		colorindex <- colorindex+1
#	}
#	#xy=locator(1)
#	legend('topleft', legend=uniquevariables, fill=colors[1:length(uniquevariables)])
#}
#
#plotIndividualTrendDataByBinaryVariable <- function(phdata, responsetype, variable, weeks='0,2,4', by.individual=T, by.group=T, colors=c('gray','red'), main=paste(responsetype,variable))
#{
#	data.stacked <- stackTrendData(phdata,responsetype,weeks)
#	plot(data.stacked$response ~ data.stacked$week, main=main, ylab=responsetype, xlab="week")
#	colorindex <- 1
#	uniquevariables <- c(0,1)
#	for (uniquevariable in uniquevariables)
#	{
#		color <- colors[colorindex]
#		iids=row.names(phdata[which(phdata[[variable]]==uniquevariable),])	
#		plotIndividualTrendLines(data.stacked,iids,color,numweeks,by.individual,by.group)
#		colorindex <- colorindex+1
#	}
#	legend('topleft', legend=uniquevariables, fill=colors[1:length(uniquevariables)])
#}
#
#plotIndividualTrendDataBySnp <- function(data,responsetype,snps=NULL,numweeks=2)
#{
#	if (is.null(snps))
#		snps=getSnpIds(data)
#	oldpar <- par(ask=T)
#	for (snp in snps)
#	{
#		plotIndividualTrendData(data,responsetype,snp,numweeks)
#	}
#	par(oldpar)
#}
#
#plotSnpAgainstResiduals <- function(phdata, response, predictors, snp)
#{
#	model <- paste(response,' ~ ',predictors)
#	fit <- lm(as.formula(model), data=phdata)
#	boxplot(fit$residuals ~ phdata[names(fit$residuals),snp], main=snp, ylab="residuals", sub=model)
#}
##makePhdataPlots(data.1b,'rs4811002', 'plt', weeks='0,1,2')
#
#plotSnpAgainstResidualsStacked <- function(phdata, response, weeks, snp)
#{
#	data.stacked <- stackTrendData(phdata,response,weeks,snp)
#	#fit <- lm(response ~ as.numeric(week), data=data.stacked)
#	fit <- lm(response ~ as.factor(week), data=data.stacked)
#	boxplot(fit$residuals ~ data.stacked[names(fit$residuals),'snp'], main=snp, ylab="residuals")
#}
##plotSnpAgainstResidualsStacked(data,'plt','0,1,2,3,4','rs13289749')
#
#plotSnpAgainstFittedStacked <- function(phdata, response, weeks, snp)
#{
#	data.stacked <- stackTrendData(phdata,response,weeks,snp)
#	fit <- lm(response ~ as.factor(snp) + as.numeric(week), data=data.stacked)
#	boxplot(fit$fitted.values ~ data.stacked[names(fit$fitted.values),'snp'], ylab='fitted')
#}
##plotSnpAgainstFittedStacked(data,'plt','0,1,2,3,4','rs13289749')
#
#plotQuantitiveTraitBySnpGenotype <- function(phdata, response, predictor, snp)
#{
#	combined <- data.frame(response=phdata[[response]], predictor=phdata[[predictor]], genotype=phdata[[snp]])
#	combined <- combined[!is.na(combined$response),]
#	combined <- combined[!is.na(combined$predictor),]
#	names(combined) <- c('response','predictor','genotype')
#	
#	plot(combined$predictor,combined$response, type='n', main=snp, xlab=predictor, ylab=response)
#	
#	colors <- c('gray','blue','red','green')
#	colorindex <- 1
#	
#	uniquegenotypes <- levels(combined$genotype)
#	uniquegenotypecounts <- data.frame(genotype=uniquegenotypes, count=0)
#	
#	for (genotype in uniquegenotypes)
#	{
#		uniquegenotypecounts[which(uniquegenotypecounts$genotype==genotype),'count'] <- length(which(combined$genotype==genotype))
#	}
#	
#	#print(uniquegenotypecounts)
#	uniquegenotypecounts <- uniquegenotypecounts[order(-uniquegenotypecounts$count),]
#	#print(uniquegenotypecounts)
#	
#	#for (genotype in uniquegenotypes)
#	for (genotype in uniquegenotypecounts$genotype)
#	{
#		color <- colors[colorindex]
#		#subset <- combined[which(combined[[snp]]==genotype),]
#		subset <- combined[which(combined$genotype==genotype),]
#		points(subset$predictor, subset$response, col=color, pch=20, cex=0.5)
#		try(abline(lm(response ~ predictor, data=subset), col=color, lwd=1), silent=TRUE)
#		colorindex <- colorindex+1
#	}
#	
#	legend('topleft', legend=uniquegenotypes, fill=colors[1:length(uniquegenotypes)])
#}
#
#makePlots <- function(data, snp, field, weeks='0,2,4', ...)
#{
#	data@phdata <- addSnpsToPhdata(data,snp)
#	phdata <- data@phdata
#	fit <- makePhdataPlots(phdata,snp,field,weeks,...)
#	return(fit)
#}
#
#makePhdataPlots <- function(phdata, snp, field, weeks='0,2,4', predictors='sex,age', logtransform=F)
#{	
#	if (is.null(phdata[[snp]]))
#		throw("snp has not been set in phdata: ", snp) 
#	oldpar <- par(mfrow=c(2,2))
#	week0 <- paste(field,'0',sep='')
#	arr <- strsplit(weeks,',')[[1]]
#	weekX <- arr[length(arr)]
#	response <- paste(field,weekX,sep='')
#	
#	if (logtransform==T)
#	{
#		phdata[[week0]] <- log(phdata[[week0]])
#		phdata[[response]] <- log(phdata[[response]])
#	}
#	
#	plotQuantitiveTraitBySnpGenotype(phdata,response,week0,snp)
#	plotSnpAgainstResiduals(phdata,response,week0,snp)
#	#plotSnpAgainstResidualsStacked(phdata,field,weeks,snp)
#	#plotSnpAgainstFittedStacked(phdata,field,weeks,snp)
#	plotIndividualTrendData(phdata,field,snp,weeks,by.individual=T,by.group=F)
#	plotIndividualTrendData(phdata,field,snp,weeks,by.individual=F,by.group=T)
#	par(mfrow=c(1,1))
#	#par(oldpar)
#	
##	if (logtransform==T)
##	{
##		week0 <- makeLogField(week0)
##		response <- makeLogField(response)
##	}
#	
#	preds <- paste(c(week0,snp,splitFields(predictors)),collapse=' + ')
#	print(preds)
#	model <- paste(response,' ~ ',preds)
#	print(model)
#	fit <- lm(as.formula(model), data=phdata)
#	print(anova(fit))
#	#print(summary(fit))
#	#crPlots(fit, one.page=TRUE, ask=FALSE)
#	#regressionDiagnostics(fit)
#	#print(stepAIC(fit, direction="backward"))
#	return(fit)
#}
##makePlots(data4,'rs7565424', 'plt', weeks='0,4', logtransform=T)
##makePlots(data4,'rs1691580', 'plt', weeks='0,1')
#
#summarizeTrendData <- function(phdata,responsetype,snp,weeks='0,2,4')
#{
#	weeksarr <- strsplit(weeks,',')[[1]]
#	week0 <- paste(responsetype,'0',sep='')	
#	weekX <- paste(responsetype,weeksarr[length(weeksarr)],sep='')
#	model <- paste(weekX,' ~ ',week0,' + ',snp)
#	fit <- lm(as.formula(model), data=phdata)
#	print(summary(fit))
#	
#	data.stacked <- stackTrendData(phdata,responsetype,weeks,snp)
#	print(summaryBy(response ~ snp + week, data=data.stacked, FUN=function(x)
#						(c(median=median(x), mean=format(mean(x),digits=2), stdev=format(sd(x),digits=2)))))
#	#aggregate(data.stacked$response, FUN=mean, by=list(week=data.stacked$week, snp=data.stacked$rs4745466))
#	#summary.formula(response ~ week + rs4745466, data=data.stacked)
#	#by(data.stacked$response, list(week=data.stacked$week,snp=data.stacked$rs4745466), function(x)(c(mean=mean(x), sd=sd(x))))
#	
#	#densityplot(~response | snp, group=week, data=data.stacked, auto.key=TRUE)
#	bwplot(response ~ week | snp, data=data.stacked, auto.key=TRUE, layout=c(3,1), main=snp)
#}
##summarizeTrendData(data5@phdata,'taqman',snp='rs4745466',weeks='0,2,4')
#
#
#getSharedSnps <- function(data, formula, top, trait.type='guess')
#{
#	group1.qt <- qtscore(formula, data, times=1, idsubset=data@phdata[data@phdata$group==1,'id'], trait.type=trait.type)
#	group2.qt <- qtscore(formula, data, times=1, idsubset=data@phdata[data@phdata$group==2,'id'], trait.type=trait.type)
#	group3.qt <- qtscore(formula, data, times=1, idsubset=data@phdata[data@phdata$group==3,'id'], trait.type=trait.type)
#	
#	group1.table <- descriptives.scan(group1.qt, sort="Pc1df", top=top)
#	group2.table <- descriptives.scan(group2.qt, sort="Pc1df", top=top)
#	group3.table <- descriptives.scan(group3.qt, sort="Pc1df", top=top)
#	
#	group1.dataframe <- data.frame(row.names=row.names(group1.table), snp=row.names(group1.table), p=group1.table$Pc1df, stringsAsFactors=FALSE)
#	group2.dataframe <- data.frame(row.names=row.names(group2.table), snp=row.names(group2.table), p=group2.table$Pc1df, stringsAsFactors=FALSE)
#	group3.dataframe <- data.frame(row.names=row.names(group3.table), snp=row.names(group3.table), p=group3.table$Pc1df, stringsAsFactors=FALSE)
#	
#	combined.dataframe <- data.frame()
#	for (snp in row.names(group1.dataframe))
#	{
#		combined.dataframe[snp,'p1'] <- group1.dataframe[snp,'p']
#	}
#	
#	for (snp in row.names(group2.dataframe))
#	{
#		combined.dataframe[snp,'p2'] <- group2.dataframe[snp,'p']
#	}
#	
#	for (snp in row.names(group3.dataframe))
#	{
#		combined.dataframe[snp,'p3'] <- group3.dataframe[snp,'p']
#	}
#	combined.dataframe <- combined.dataframe[which(!is.na(combined.dataframe$p1) & (!is.na(combined.dataframe$p2) | !is.na(combined.dataframe$p3))),]
#	return(combined.dataframe)
#}
#
#addSnpGenotypeField <- function(data, snp, refgenotype, values1, values2)
#{
#	field <- paste(snp,refgenotype,sep='')
#	data[[field]] <- NA
#	data[which(data[[snp]] %in% splitFields(values1)),field] <- 1
#	data[which(data[[snp]] %in% splitFields(values2)),field] <- 0
#	return(data)
#}
#
#
#createMinimumField <- function(data,newfield,fields)
#{
#	data[,newfield] <- NA
#	for (rowname in row.names(data))
#	{
#		values <- c()
#		for (field  in fields)
#		{
#			value <- data[rowname,field]
#			if (!is.na(value))
#				values <- c(values,value)
#		}
#		#print(values)
#		if (length(values)>0)
#			data[rowname,newfield]<- min(values, na.rm=T)
#	}
#	return(data)
#}
#
#############################################################################
#
#getSnpByWeekInteraction <- function(data, responsetype, snp, weeks='0,1,2,3,4', make.plot=T)
#{
#	data.stacked <- stackTrendData(data,responsetype,weeks,snp)
#	#data.stacked <- na.omit(data.stacked)
#	if (make.plot)
#	{
#		title <- snp
#		with(data.stacked,
#		{
#			interaction.plot(week, factor(snp), response, lwd=3, main=title, ylab=responsetype, xlab="week", trace.label="group") #ylim=c(10, 60), lty=c(1, 12),
#		})
#	}
#	data.stacked.aov <- aov(data.stacked$response ~ factor(data.stacked$snp)*factor(data.stacked$week) + Error(factor(data.stacked$iid)))
#	data.stacked.aov.sum <- summary(data.stacked.aov)
#	#data.stacked.aov.sum.iid <- data.stacked.aov.sum[["Error: factor(data.stacked$iid)"]]
#	p.value <- data.stacked.aov.sum[['Error: Within']][[1]]['factor(data.stacked$snp):factor(data.stacked$week)','Pr(>F)']
#	#p.value <- data.stacked.aov.sum.iid[[1]]["factor(data.stacked$snp):factor(data.stacked$week)","Pr(>F)"]
#	return(p.value)
#}
##getSnpByWeekInteraction(data.1b,responsetype='plt',snp='rs4811002')
#
#getSnpByWeekInteractions <- function(data, responsetype, snps, weeks='0,1,2,3,4')
#{
#	pvalues <- NULL
#	for (snp in snps)
#	{
#		p.value <- getSnpByWeekInteraction(data,responsetype=responsetype,snp=snp, weeks=weeks, make.plot=F)
#		marker <- createSignificantMarker(p.value)
#		#print(paste(snp,'=',p.value))
#		if (is.null(pvalues))
#			pvalues <- data.frame(snp=snp, p=p.value, marker=marker)	
#		else pvalues <- rbind(pvalues, data.frame(snp=snp, p=p.value, marker=marker))
#	}
#	pvalues <- pvalues[order(pvalues$p),]
#	return(pvalues)
#}
##getSnpByWeekInteractions(data.1b,responsetype='plt',snps=pltfields)
#
#makeSnpByWeekInteractionPlot <- function(data, responsetype, snp, weeks)
#{
#	mypanel1 <- function(x, y, subscripts, groups)
#	{
#		panel.superpose(x,y,subscripts,groups)
#		#panel.xyplot(x,y,pch=19)
#		#panel.rug(x,y)
#		panel.grid(h=-1, v=-1)
#		panel.lmline(x,y,col='red', lwd=1, lty=2)
#	}
#
#	data.stacked <- stackTrendData(data,responsetype,weeks,snp)
#	#data.stacked <- na.omit(data.stacked)
#	graph1 <- xyplot(response ~ as.numeric(week) | snp, groups=iid, type="o", data.stacked, 
#				layout=c(3,1), main=snp, xlab='Weeks', ylab=responsetype, panel=mypanel1)#panel.superpose
#	graph2 <- bwplot(response ~ week | snp, data.stacked, 
#				layout=c(3,1), main=snp, xlab='Weeks', ylab=responsetype)
#	plot(graph1, split=c(1,1,1,2))
#	plot(graph2, split=c(1,2,1,2), newpage=FALSE)
#	#plot(graph1, split=c(1,1,2,1))
#	#plot(graph2, split=c(2,1,2,1), newpage=FALSE)
#}
##makeSnpByWeekInteractionPlot(data.1b,'plt','rs4811002','0,0.5,1,2,3,4')
#
#######################################################
#
#getLinearSnpByWeekInteraction <- function(data, responsetype, snp, weeks='0,1,2,3,4', make.plot=F)
#{
#	data.stacked <- stackTrendData(data,responsetype,weeks,snp)
#	try({
#		time.linear <- lme(response ~ factor(snp)*as.numeric(week), random=list(iid = pdDiag(~as.numeric(week))), data.stacked)
#		#summary(time.linear)
#		time.linear.aov <- anova(time.linear)
#		p.value <- time.linear.aov['factor(snp):as.numeric(week)','p-value']
#		return(p.value)},silent=F)
#	return(1.0)
#}
##getLinearSnpByWeekInteraction(data,responsetype='plt',snp='rs4811002')
##getLinearSnpByWeekInteraction(data.1b,responsetype='plt',snp='rs4811002')
#
#getLinearSnpByWeekInteractions <- function(data, responsetype, snps, weeks='0,1,2,3,4')
#{
#	pvalues <- NULL
#	for (snp in snps)
#	{
#		p.value <- getLinearSnpByWeekInteraction(data,responsetype=responsetype,snp=snp, weeks=weeks)
#		marker <- createSignificantMarker(p.value)
#		#print(paste(snp,'=',p.value))
#		if (is.null(pvalues))
#			pvalues <- data.frame(snp=snp, p=p.value, marker=marker)	
#		else pvalues <- rbind(pvalues, data.frame(snp=snp, p=p.value, marker=marker))
#	}
#	pvalues <- pvalues[order(pvalues$p),]
#	return(pvalues)
#}
##getLinearSnpByWeekInteractions(data.1b,responsetype='plt',snps=pltfields)
#
##################################################################################33
#
#getGlsSnpByWeekInteraction <- function(data, responsetype, snp, weeks='0,1,2,3,4', make.plot=F)
#{
#	data.stacked <- stackTrendData(data,responsetype,weeks,snp,omit.na=T)
#	try({
#		data.stacked.grouped <- groupedData(response ~ as.numeric(factor(snp))*as.numeric(week) | iid, data=data.stacked)
#		#fit.cs <- gls(response ~ factor(snp)*factor(week), data=data.stacked.grouped, corr=corCompSymm(, form= ~ 1 | iid) )
#		fit.un <- gls(response ~ factor(snp)*factor(week), data=data.stacked.grouped, corr=corSymm(form = ~ 1 | iid), weights = varIdent(form = ~ 1 | week))
#		#fit.ar1 <- gls(response ~ factor(snp)*factor(week), data=data.stacked.grouped, corr=corAR1(, form= ~ 1 | iid))
#		#fit.arh1 <- gls(response ~ factor(snp)*factor(week), data=data.stacked.grouped, corr=corAR1(, form= ~ 1 | iid), weight=varIdent(form = ~ 1 | week))
#				
#		#print(snp)
#		#aicscores <- data.frame(type='cs', aic=summary(fit.cs)$AIC, p=anova(fit.cs)['factor(snp):factor(week)','p-value'])
#		#aicscores <- rbind(aicscores, data.frame(type='un', aic=summary(fit.un)$AIC, p=anova(fit.un)['factor(snp):factor(week)','p-value']))
#		#aicscores <- rbind(aicscores, data.frame(type='ar1', aic=summary(fit.ar1)$AIC, p=anova(fit.ar1)['factor(snp):factor(week)','p-value']))
#		#aicscores <- rbind(aicscores, data.frame(type='arh1', aic=summary(fit.arh1)$AIC, p=anova(fit.arh1)['factor(snp):factor(week)','p-value']))
#		#aicscores <- aicscores[order(aicscores$aic),]
#		#print(aicscores)
#		#method <- aicscores[1,'type']
#		#pvalue <- aicscores[1,'p']
#		
#		pvalue <- anova(fit.un)['factor(snp):factor(week)','p-value']
#		print(paste(snp,pvalue))
#		return(pvalue)
#	},silent=T)
#	return(1.0)
#}
##getGlsSnpByWeekInteraction(data,responsetype='plt',snp='rs13289749')
#
#getGlsSnpByWeekInteractions <- function(data, responsetype, snps, weeks='0,1,2,3,4')
#{
#	pvalues <- NULL
#	for (snp in snps)
#	{
#		p.value <- getGlsSnpByWeekInteraction(data,responsetype=responsetype,snp=snp, weeks=weeks)
#		marker <- createSignificantMarker(p.value)
#		#print(paste(snp,'=',p.value))
#		if (is.null(pvalues))
#			pvalues <- data.frame(snp=snp, p=p.value, marker=marker)	
#		else pvalues <- rbind(pvalues, data.frame(snp=snp, p=p.value, marker=marker))
#	}
#	pvalues <- pvalues[order(pvalues$p),]
#	return(pvalues)
#}
##getGlsSnpByWeekInteractions(data,responsetype='plt',snps=pltfields)

