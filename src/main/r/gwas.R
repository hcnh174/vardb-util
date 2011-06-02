library(Design,T)

#plink.dir <- "c:/projects/analysis/gwas/plink/"
#setwd(plink.dir)

createPlinkQQPlot <- function(data)
{
	#data <- read.table(file=filename, header=TRUE)
	print(data[1:10,])
	plot(-log(data$QQ,10), -log(data$UNADJ,10), main="QQ-plot", xlab="QQ", ylab="UNADJ")
	abline(a=0, b=1)
	identify(-log(data$QQ,10), -log(data$UNADJ,10), labels=data$SNP)
}

plotQuantitiveTraitBySnpGenotype <- function(data, response, predictor, snp)
{
	plot(makeFormula(response,predictor), type="n", main=snp, data=data)
	AA <- data[data[[snp]]=='AA',]
	AB <- data[data[[snp]]=='AB',]
	BB <- data[data[[snp]]=='BB',]
	points(AA[[predictor]], AA[[response]], col="black")
	points(AB[[predictor]], AB[[response]], col="blue")
	points(BB[[predictor]], BB[[response]], col="red")	
}

plotSnpAgainstResiduals <- function(data, response, predictor, snp)
{
	fit <- lm(makeFormula(response,predictor), data=data)
	boxplot(fit$residuals ~ data[names(fit$residuals),snp], main=snp, ylab="residuals")
}

plotIndividualTrendDataBySnp <- function(data,responsetype,snps=NULL,numweeks=2)
{
	if (is.null(snps))
		snps=getSnpIds(data)
	oldpar <- par(ask=T)
	for (snp in snps)
	{
		plotIndividualTrendData(data,responsetype,snp,numweeks)
	}
	par(oldpar)
}

plotIndividualTrendLines <- function(data024,iids,color,numweeks)
{
	for (iid in iids)
	{
		#print(paste("iid=",iid))
		individual <- data024[which(data024$iid==iid),]
		if (nrow(individual)==numweeks)
		{
			#print(paste(iid,genotype))
			lines(x=individual$week, y=individual$response, col=color)
		}
	}
}

plotIndividualTrendData <- function(data,responsetype,snp,numweeks)
{
	data0 <- data.frame(iid=data$IID, week=0, response=data[[paste(responsetype,"0",sep="")]])
	data2 <- data.frame(iid=data$IID, week=2, response=data[[paste(responsetype,"2",sep="")]])
	data4 <- data.frame(iid=data$IID, week=4, response=data[[paste(responsetype,"4",sep="")]])
	if (numweeks==2)
		data024 <- rbind(data0,data4)
	else data024 <- rbind(data0,data2,data4)
	
	data024 <- data024[!is.na(data024$response),]
	
	aa <- data[data[[snp]]=='AA','IID']
	ab <- data[data[[snp]]=='AB','IID']
	bb <- data[data[[snp]]=='BB','IID']
	
	#print(paste("aa=",length(aa)))
	#print(paste("ab=",length(ab)))
	#print(paste("bb=",length(bb)))

	plot(data024$response ~ data024$week, main=paste(responsetype,snp), ylab=responsetype, xlab="week")
	
	if (length(aa) > length(bb))
	{
		plotIndividualTrendLines(data024,aa,'gray',numweeks)
		plotIndividualTrendLines(data024,ab,'green',numweeks)
		plotIndividualTrendLines(data024,bb,'red',numweeks)	
	}
	else
	{
		plotIndividualTrendLines(data024,bb,'gray',numweeks)
		plotIndividualTrendLines(data024,ab,'green',numweeks)
		plotIndividualTrendLines(data024,aa,'red',numweeks)
	}	
}
		
getSnpIds <- function(data)
{
	return(names(pltdata)[grep(pattern="^rs[0-9]+", names(pltdata))])
}

getPolymorphicSnpIds <- function(data)
{
	snps <- getSnpIds(data)
	for (snp in snps)
	{
		print(paste(snp,length(unique(data[[snp]]))))
	}
}
