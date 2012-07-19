library(ggplot2)
library(gridExtra)
library(reshape)
library(limma)

cutoff <- 1.5

alpha <- 90
col1 <- rgb(0,255,0,alpha,maxColorValue=255)
col2 <- rgb(255,0,0,alpha,maxColorValue=255) 

databases <- splitFields('miranda,mirbase,targetscan'); at.least <- 3

loadHccData <- function(filename)
{
	deadentries <- loadDataFrame('dead-entries.txt')
	data <- loadDataFrame(filename, idcol='mirna')
	data <- subset(data, !(mirna %in% deadentries$mirna))
	#data <- cbind(data[,c(1,2)],normalizeBetweenArrays(as.matrix(data[,c(-1,-2)], method='scale')))
	return(data)
}

getPairedData <- function(data)
{
	#data.paired <- data[,c(splitFields('mirna,id'),concat('tumor',1:6),concat('nontumor',1:6),concat('serum',1:6))]
	data.paired <- data[,c(splitFields('mirna,id'),concat('tumor',1:6),concat('nontumor',1:6))]
	return(data.paired)
}

getSerumData <- function(data)
{
	data.serum <- data[,c(splitFields('mirna,id'),concat('healthy',1:10),concat('hcc',1:16),concat('hbv',c(2:8,10:12)))]
	#data.serum <- data[,c(splitFields('mirna,id'),concat('healthy',1:10),concat('hcc',1:16),concat('serum',1:6),concat('hbv',c(2:8,10:12)))]
	#data.serum <- renameColumns(data.serum,concat('serum',1:6),concat('hcc',11:16))
	return(data.serum)
}

appendSubjectData <- function(data.long, filename='subjects.txt')
{
	subjects <- loadDataFrame(filename,id='subject')
	
	for (field in splitFields('healthy,hcc,chbv,advanced,hbv,hcv,nbnc,hbvhcc,hcvhcc,nbnchcc,virus,age,BH,BW,WBC,NE,HGB,PLT,AST,ALT,gGTP,TB,DB,ALPH,ALB,CRE,AFP,PT,BS,A1C,hdl,cholesterol,ha,ft4,ft3,TSH,Fe,TG,FERRITIN,HBsAg,HBsAb,HBeAg,HBeAb,BUN,PIVKA,PCR'))#PCR
	{
		data.long[[field]] <- sapply(data.long$subject, function(subject)
		{
			return(subjects[which(subjects$subject==subject),field][1])
		})
	}
	return(data.long)
}

meltPairedArrayData <- function(data)
{
	data.long <- melt(data, id=c('mirna','id'))
	data.long$variable <- as.character(data.long$variable) 
	data.long$trt <- sub("[0-9]+","",data.long$variable)
	data.long$subject <- concat('hcc',as.numeric(sub("[a-zA-Z]+","",data.long$variable))+10)
	data.long <- appendSubjectData(data.long)
	return(data.long)
}

meltSerumArrayData <- function(data)
{
	#data.long <- meltArrayData(data)
	data.long <- melt(data, id=c('mirna','id'))
	data.long$variable <- as.character(data.long$variable) 
	data.long$trt <- sub("[0-9]+","",data.long$variable)
	data.long <- renameColumn(data.long, 'variable', 'subject')
	data.long$trt <- factor(data.long$trt)#, levels=c('nontumor','tumor','serum'))#,'diff'
	data.long <- appendSubjectData(data.long)
	return(data.long)
}

as.yesno(healthy) <- function(value)
{
	if (value %in% c('y','n'))
		return(value)
	return(ifelse(value==TRUE,'y','n'))
}
#as.yesno('y')

getPatientSubset <- function(data.long, trt=NULL, healthy=NULL, hcc=NULL, chbv=NULL, hbv=NULL, hcv=NULL, nbnc=NULL, advanced=NULL)
{
	tmp <- data.long
	if (!is.null(trt)) tmp <- tmp[which(tmp$trt %in% splitFields(trt)),]
	if (!is.null(healthy)) tmp <- tmp[which(tmp$healthy==as.yesno(healthy)),]
	if (!is.null(hcc)) tmp <- tmp[which(tmp$hcc==as.yesno(hcc)),]
	if (!is.null(chbv)) tmp <- tmp[which(tmp$chbv==as.yesno(chbv)),]
	if (!is.null(hbv)) tmp <- tmp[which(tmp$hbv==as.yesno(hbv)),]
	if (!is.null(hcv)) tmp <- tmp[which(tmp$hcv==as.yesno(healthy)),]
	if (!is.null(advanced)) tmp <- tmp[which(tmp$advanced==as.yesno(advanced)),]
	if (!is.null(nbnc)) tmp <- tmp[which(tmp$nbnc==as.yesno(nbnc)),]
	print(nrow(tmp))
	print(table(tmp$trt))
	return(tmp)
}
#data.tmp <- getPatientSubset(data.serum.long, chbv=FALSE)

plotMirna <- function(data.long, mirna, field, xlab=field, geom=c("boxplot","point"), ylim=NULL)
{
	data.mirna <- data.long[which(data.long$mirna==mirna),]
	data.mirna <- data.mirna[!is.na(data.mirna$value),]
	data.mirna <- data.mirna[!is.na(data.mirna[[field]]),]
	
	fit <- kruskal.test(as.formula(concat('value ~ as.factor(',field,')')), data.mirna)
	fit$table <- aggregate(as.formula(concat('value ~ ',field)), data = data.mirna, median)
	main <- concat('P=',format(fit$p.value, digits=5))
	if (!is.null(ylim))
		fit$plot <- qplot(data.mirna[[field]], data.mirna$value, geom=geom, main=main, xlab=xlab, ylab=mirna, ylim=ylim)
	else fit$plot <- qplot(data.mirna[[field]], data.mirna$value, geom=geom, main=main, xlab=xlab, ylab=mirna)
	if (nrow(fit$table==2))
		fit$ratio <- log2(fit$table[2,2])/log2(fit$table[1,2])
	return(fit)
}
#fit <- plotMirna(data.serum.long, 'hsa-miR-122', 'hcc')

analyzeMirna <- function(data.paired.long, data.serum.long, mirna, geom=c("boxplot","point"), make.plots=TRUE)
{
	data.tmp <- data.serum.long
	fit.hcc <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy/HBV vs. HCC', )
	
	data.tmp <- getPatientSubset(data.serum.long, trt='healthy,hcc')
	fit.hcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy controls vs. HCC', )

	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='n')),]
	fit.early_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy controls vs. early HCC', )
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='y')),]
	fit.advanced_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy controls vs. advanced HCC', )
	
	data.tmp <- getPatientSubset(data.paired.long, trt='tumor,nontumor')	
	fit.paired <- plotMirna(data.tmp, mirna, 'trt', xlab='Paired liver samples', )	
	
	data.tmp <- getPatientSubset(data.serum.long, trt='healthy,hcc', hbv=FALSE)
	fit.nonhbvhcc <- plotMirna(data.tmp, mirna, 'hcc', xlab='healthy vs. non-HBV-HCC', )
	
	data.tmp <- getPatientSubset(data.serum.long, trt='hcc')
	fit.hbvhcc <- plotMirna(data.tmp, mirna, 'hbvhcc', xlab='Non-HBV-HCC vs. HBV-HCC', )
	
	data.tmp <- getPatientSubset(data.serum.long, trt='hcc')
	fit.hcvhcc <- plotMirna(data.tmp, mirna, 'hcvhcc', xlab='Non-HCV-HCC vs. HCV-HCC', )
	
	data.tmp <- getPatientSubset(data.serum.long, trt='healthy,hbv')
	fit.hbv_vs_healthy <- plotMirna(data.tmp, mirna, 'hbv', xlab='healthy vs. HBV', )
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='n' & data.serum.long$hbv=='y')),]
	fit.early_hbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy controls vs. early HBV-related HCC', )
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='y' & data.serum.long$hbv=='y')),]
	fit.advanced_hbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy controls vs. advanced HBV-related HCC', )
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='n' & data.serum.long$hbv=='n')),]
	fit.early_nonhbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy controls vs. early non-HBV-related HCC', )
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='y' & data.serum.long$hbv=='n')),]
	fit.advanced_nonhbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy controls vs. advanced non-HBV-related HCC', )
	
	if (make.plots)
	{
		pdf(concat('out/plots-',make.names(mirna),'.pdf'))
		#age.plot <- plot(value ~ age, data.serum.long[which(data.serum.long$mirna==mirna),], main=mirna)
		print(grid.arrange(fit.hcc$plot, fit.hcc_vs_healthy$plot, #age.plot, 
				fit.early_vs_healthy$plot, fit.advanced_vs_healthy$plot, fit.paired$plot, nrow=2))
		dev.off()
	}
	return(list(paired=fit.paired, hcc=fit.hcc, hcc_vs_healthy=fit.hcc_vs_healthy, 
					nonhbvhcc=fit.nonhbvhcc, hbvhcc=fit.hbvhcc,	hcvhcc=fit.hcvhcc, hbv_vs_healthy=fit.hbv_vs_healthy,
					early_vs_healthy=fit.early_vs_healthy, advanced_vs_healthy=fit.advanced_vs_healthy, 
					early_hbvhcc_vs_healthy=fit.early_hbvhcc_vs_healthy, advanced_hbvhcc_vs_healthy=fit.advanced_hbvhcc_vs_healthy,
					early_nonhbvhcc_vs_healthy=fit.early_nonhbvhcc_vs_healthy, advanced_nonhbvhcc_vs_healthy=fit.advanced_nonhbvhcc_vs_healthy))
}
#analyzeMirna(data.paired.long, data.serum.long, mirna='hsa-miR-122', make.plots=TRUE)


getMirnaList <- function(filename='list.txt')
{
	data <- loadDataFrame(filename)
	return(data$mirna)
}
#getMirnaList()

analyzeMirnas <- function(data.paired.long, data.serum.long, mirnas=unique(data.serum.long$mirna), make.plots=FALSE)
{
	results <- data.frame(mirna=mirnas)
	rownames(results) <- results$mirna
	for (mirna in splitFields(mirnas))
	{
		try({
			fits <- analyzeMirna(data.paired.long, data.serum.long, mirna=mirna, make.plots=make.plots)	
			results[mirna,'hcc'] <- fits$hcc$p.value
			results[mirna,'hcc_vs_healthy'] <- fits$hcc_vs_healthy$p.value
			results[mirna,'paired'] <- fits$paired$p.value
			results[mirna,'nonhbvhcc'] <- fits$nonhbvhcc$p.value
			results[mirna,'hbvhcc'] <- fits$hbvhcc$p.value
			results[mirna,'hcvhcc'] <- fits$hcvhcc$p.value
			results[mirna,'hbv_vs_healthy'] <- fits$hbv_vs_healthy$p.value		
			results[mirna,'early_vs_healthy'] <- fits$early_vs_healthy$p.value
			results[mirna,'advanced_vs_healthy'] <- fits$advanced_vs_healthy$p.value
			results[mirna,'early_hbvhcc_vs_healthy'] <- fits$early_hbvhcc_vs_healthy$p.value
			results[mirna,'advanced_hbvhcc_vs_healthy'] <- fits$advanced_hbvhcc_vs_healthy$p.value
			results[mirna,'early_nonhbvhcc_vs_healthy'] <- fits$early_nonhbvhcc_vs_healthy$p.value
			results[mirna,'advanced_nonhbvhcc_vs_healthy'] <- fits$advanced_nonhbvhcc_vs_healthy$p.value
		})
	}
	results <- results[order(results$hcc, results$hcc_vs_healthy, results$paired, results$nonhbvhcc, results$hbvhcc, results$hcvhcc, results$hbv_vs_healthy),]
	writeTable(results,concat('out/results.txt'), row.names=FALSE)
	return(results)
}
#results <- analyzeMirnas(data.paired.long, data.serum.long)
#results <- analyzeMirnas(data.paired.long, data.serum.long, unique(c(results.hcc$mirna, results.hcc_vs_healthy$mirna)))
#results <- analyzeMirnas(data.paired.long, data.serum.long, mirnas=getMirnaList(), make.plots=TRUE)


#############################################################################

compareTreatments <- function(data.long, trt1, trt2, field='trt', pvalue=0.05, cutoff=NULL, minlength=2, make.plots=FALSE, name=NULL, adjust='fdr', p.field='p.kruskal', p.adj=NULL, keep.na=TRUE)
{
	data.long <- data.long[which(data.long[[field]] %in% c(trt1,trt2)),]
	data.long <- data.long[!is.na(data.long[[field]]),]
	data <- cast(data.long, as.formula(concat('mirna ~ ',field)), value='value', fun.aggregate=median)
	data$ratio <- log2(data[[trt2]]/data[[trt1]])
	rownames(data) <- data$mirna
	for (mirna in rownames(data))
	{
		try({
			data.mirna <- data.long[which(data.long$mirna==mirna),]
			data.mirna <- data.mirna[!is.na(data.mirna$value),]
			values1 <- data.mirna[which(data.mirna[[field]]==trt1),'value']
			values2 <- data.mirna[which(data.mirna[[field]]==trt2),'value']
			#values1 <- subset(data.mirna, trt==trt1)$value
			#values2 <- subset(data.mirna, trt==trt2)$value
			#fit <- kruskal.test(value ~ as.factor(trt), data.mirna)
			fit <- kruskal.test(as.formula(concat('value ~ as.factor(',field,')')), data.mirna)
			data[mirna,'n1'] <- length(values1)
			data[mirna,'n2'] <- length(values2)
			data[mirna,'p.kruskal'] <- fit$p.value	
			fit <- t.test(values1,values2)
			data[mirna,'p.t'] <- fit$p.value
		}, silent=TRUE)
	}
	#print(head(data))
	data$p.value <- data[[p.field]]
	data$p.adjusted <- p.adjust(data$p.value, method=adjust)#, n=length(data$p.value))
	data <- data[order(-data$ratio, data$p.adjusted, data$p.value),]
	#data <- data[order(data$t,-data$ratio),]
	if (!is.null(pvalue))
		data <- data[which(data$p.value<=pvalue),]	#data <- data[which(data$t<=pvalue),]
	if (!is.null(p.adj))
		data <- data[which(data$p.adjusted<=p.adj),]	#data <- data[which(data$t<=pvalue),]
	if (!is.null(cutoff))
		data <- data[which(abs(data$ratio) >= cutoff),]
	if (!is.null(name))
		writeTable(data,concat('out/compare-',name,'.txt'), row.names=FALSE)	
	if (!is.null(name) & make.plots)
	{
		pdffile <- concat('out/plots-',name,'.pdf')
		geom <- c("boxplot","point")
		pdf(pdffile)
		for (mirna in data$mirna)
		{
			#print(qplot(trt, value, data=data.long[which(data.long$mirna==mirna),], geom=c("boxplot","point"), main=concat(mirna,' (',format(data[mirna,'p.value'], digits=5),')')))
			#boxplot(as.formula(concat('value ~ ',field)), data.long[which(data.long$mirna==mirna),], main=mirna, sub=data[mirna,'p.value'])
			#fit <- kruskal.test(as.formula(concat('value ~ as.factor(',field,')')), data.mirna)
			data.mirna <- data.long[which(data.long$mirna==mirna),]
			data.mirna <- data.mirna[!is.na(data.mirna$value),]
			#data.mirna <- data.mirna[!is.na(data.mirna[[field]]),]
			#print(qplot(field, value, data=data.mirna, geom=c("boxplot","point"), main=concat(mirna,' (',format(data[mirna,'p.value'], digits=5),')')))
			main <- concat(mirna,' (',format(data[mirna,'p.value'], digits=5),')')
			print(qplot(data.mirna[[field]], data.mirna$value, geom=geom, main=main, xlab=field, ylab=mirna))
		}
		dev.off()
	}
	createWekaFile(data.long, mirnas=data$mirna, class=field, name=name)
	return(data)
}
#compareTreatments(data.paired.long,'nontumor','tumor', cutoff=cutoff, make.plots=TRUE, name='paired'); results.paired


compareCovariates <- function(data.long, field='trt', pvalue=0.05, cutoff=NULL, minlength=2, make.plots=FALSE, name=field, adjust='fdr', p.field='p.kruskal', p.adj=NULL, keep.na=TRUE)
{
	#data <- cast(data.long, as.formula(concat('mirna ~ ',field)), value='value', fun.aggregate=median)
	data <- data.frame(mirna=unique(data.long$mirna))
	rownames(data) <- data$mirna
	for (mirna in rownames(data))
	{
		try({
					data.mirna <- data.long[which(data.long$mirna==mirna),]
					data.mirna <- data.mirna[!is.na(data.mirna$value),]
					fit <- lm(as.formula(concat('value ~ ',field)), data.mirna)
					fit.sum <- summary(fit)
					data[mirna,'p.value'] <- fit.sum$coefficients[field,'Pr(>|t|)']
				}, silent=TRUE)
	}
	data$p.adjusted <- p.adjust(data$p.value, method=adjust)#, n=length(data$p.value))
	data <- data[order(data$p.adjusted, data$p.value),]
	if (!is.null(pvalue))
		data <- data[which(data$p.value<=pvalue),]	#data <- data[which(data$t<=pvalue),]
	if (!is.null(p.adj))
		data <- data[which(data$p.adjusted<=p.adj),]	#data <- data[which(data$t<=pvalue),]
	if (!is.null(name))
		writeTable(data,concat('out/compare-',name,'.txt'), row.names=FALSE)	
	if (!is.null(name) & make.plots)
	{
		pdffile <- concat('out/plots-',name,'.pdf')
		pdf(pdffile)
		for (mirna in data$mirna)
		{
			data.mirna <- data.long[which(data.long$mirna==mirna),]
			plot(as.formula(concat('value ~ ',field)), data.mirna, main=mirna[1], sub=data[mirna,'p.value'])
			abline(lm(as.formula(concat('value ~ ',field)), data.mirna))
		}
		dev.off()
	}
	return(data)
}
#compareCovariates(data.serum.long, field='ALT', make.plots=TRUE)

createWekaFile <- function(data.long, mirnas, class='trt', name=class)
{
	mirnas <- splitFields(mirnas)
	data.long$class <- data.long[[class]]
	#data.long <- subset(data.long, !is.na(class))
	tmp <- subset(data.long, mirna %in% mirnas)[,c(splitFields('mirna,subject,value,class,age,ALT,hbv,hcv'))]
	data.pred <- cast(tmp, class + subject + age + ALT + hbv + hcv ~ mirna, value='value', fun.aggregate=median) 
	colnames(data.pred) <- make.names(colnames(data.pred))
	data.pred$class <- factor(data.pred$class)
	data.weka <- data.pred[,c(make.names(mirnas),'class')]
	filename <- concat('out/',name,'.arff')
	write.arff(data.weka, file = filename, class_col='class')
	printcat('wrote arff table to file: ',filename)
}
