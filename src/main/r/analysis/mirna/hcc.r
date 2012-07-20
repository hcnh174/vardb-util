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
	
	for (field in splitFields('healthy,hcc,chbv,early,advanced,hbv,hcv,nbnc,hbvhcc,hcvhcc,nbnchcc,virus,age,BH,BW,WBC,NE,HGB,PLT,AST,ALT,gGTP,TB,DB,ALPH,ALB,CRE,AFP,PT,BS,A1C,hdl,cholesterol,ha,ft4,ft3,TSH,Fe,TG,FERRITIN,HBsAg,HBsAb,HBeAg,HBeAb,BUN,PIVKA,PCR'))#PCR
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

as.yesno <- function(value)
{
	if (value %in% c('y','n'))
		return(value)
	return(ifelse(value==TRUE,'y','n'))
}
#as.yesno('y')

getPatientSubset <- function(data.long, trt=NULL, healthy=NULL, hcc=NULL, chbv=NULL, hbv=NULL, hcv=NULL, nbnc=NULL, advanced=NULL, verbose=FALSE)
{
	tmp <- data.long
	if (!is.null(trt)) tmp <- tmp[which(tmp$trt %in% splitFields(trt)),]
	if (!is.null(healthy)) tmp <- tmp[which(tmp$healthy==as.yesno(healthy)),]
	if (!is.null(hcc)) tmp <- tmp[which(tmp$hcc==as.yesno(hcc)),]
	if (!is.null(chbv)) tmp <- tmp[which(tmp$chbv==as.yesno(chbv)),]
	if (!is.null(hbv)) tmp <- tmp[which(tmp$hbv==as.yesno(hbv)),]
	if (!is.null(hcv)) tmp <- tmp[which(tmp$hcv==as.yesno(hcv)),]
	if (!is.null(advanced)) tmp <- tmp[which(tmp$advanced==as.yesno(advanced)),]
	if (!is.null(nbnc)) tmp <- tmp[which(tmp$nbnc==as.yesno(nbnc)),]
	if (verbose)
	{
		print(nrow(tmp))
		print(table(tmp$trt))
		print(unique(tmp$subject))
	}
	return(tmp)
}
#data.tmp <- getPatientSubset(data.serum.long, chbv=FALSE)

plotMirna <- function(data.long, mirna, field, xlab=field, geom=c("boxplot","point"), ylim=NULL, paired=FALSE, fontsize=10, color='black')
{
	data.mirna <- data.long[which(data.long$mirna==mirna),]
	data.mirna <- data.mirna[!is.na(data.mirna$value),]
	data.mirna <- data.mirna[!is.na(data.mirna[[field]]),]
	
	#fit <- kruskal.test(as.formula(concat('value ~ as.factor(',field,')')), data.mirna)
	fit <- wilcox.test(as.formula(concat('value ~ as.factor(',field,')')), data.mirna, paired=paired)
	fit$table <- aggregate(as.formula(concat('value ~ ',field)), data = data.mirna, median)
	main <- concat('P=',format(fit$p.value, digits=5))
	if (fit$p.value<0.05) color='red'
	if (!is.null(ylim))
		fit$plot <- qplot(data.mirna[[field]], data.mirna$value, geom=geom, main=main, xlab=xlab, ylab=mirna, ylim=ylim)
	else fit$plot <- qplot(data.mirna[[field]], data.mirna$value, geom=geom, main=main, xlab=xlab, ylab=mirna)
	fit$plot <- fit$plot + opts(
							plot.title=theme_text(size=fontsize, colour=color),
							axis.title.x=theme_text(size=fontsize),
							axis.title.y=theme_text(size=fontsize, angle = 90)
							) 
	if (nrow(fit$table==2))
		fit$ratio <- log2(fit$table[2,2])/log2(fit$table[1,2])
	return(fit)
}
#fit <- plotMirna(data.serum.long, 'hsa-miR-122', 'hcc'); fit$plot


analyzeMirna <- function(data.paired.long, data.serum.long, mirna, geom=c("boxplot","point"), make.plots=TRUE)
{
	ylim <- c(0,ceiling(max(data.serum.long[which(data.serum.long$mirna==mirna),'value'])))
	data.tmp <- data.serum.long
	fit.hcc_vs_healthyhbv <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy/HBV v HCC', ylim=ylim)
	
	data.tmp <- getPatientSubset(data.serum.long, trt='healthy,hcc')
	fit.hcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy v HCC', ylim=ylim)

	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='n')),]
	fit.early_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy v early HCC', ylim=ylim)
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='y')),]
	fit.advanced_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy v advanced HCC', ylim=ylim)
	
	data.tmp <- getPatientSubset(data.paired.long, trt='tumor,nontumor')
	fit.tumor_vs_nontumor <- plotMirna(data.tmp, mirna, 'trt', xlab='Paired liver samples', paired=TRUE)
	
	###	
	data.tmp <- getPatientSubset(data.serum.long, trt='healthy,hcc', hcv=FALSE, nbnc=FALSE)
	fit.hbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='healthy v HBV-HCC', ylim=ylim)

	data.tmp <- getPatientSubset(data.serum.long, trt='healthy,hcc', hbv=FALSE)
	fit.nonhbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='healthy v non-HBV-HCC', ylim=ylim)
	
	data.tmp <- getPatientSubset(data.serum.long, trt='hcc')
	fit.nonhbvhcc_vs_hbvhcc <- plotMirna(data.tmp, mirna, 'hbvhcc', xlab='Non-HBV-HCC v HBV-HCC', ylim=ylim)
	
	data.tmp <- getPatientSubset(data.serum.long, trt='hcc')
	fit.nonhcvhcc_vs_hcvhcc <- plotMirna(data.tmp, mirna, 'hcvhcc', xlab='Non-HCV-HCC v HCV-HCC', ylim=ylim)
	
	data.tmp <- getPatientSubset(data.serum.long, trt='healthy,hbv')
	fit.hbv_vs_healthy <- plotMirna(data.tmp, mirna, 'hbv', xlab='healthy v HBV', ylim=ylim)
	
	data.tmp <- getPatientSubset(data.serum.long, trt='hcc,hbv', hbv=TRUE)
	fit.hbvhcc_vs_hbv <- plotMirna(data.tmp, mirna, 'hcc', xlab='HBV v HBV-HCC', ylim=ylim)
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='n' & data.serum.long$hbv=='y')),]
	fit.early_hbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy v early HBV-HCC', ylim=ylim)
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='y' & data.serum.long$hbv=='y')),]
	fit.advanced_hbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy v advanced HBV-HCC', ylim=ylim)
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='n' & data.serum.long$hbv=='n')),]
	fit.early_nonhbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy v early non-HBV-HCC', ylim=ylim)
	
	data.tmp <- data.serum.long[which(data.serum.long$trt=='healthy' | (data.serum.long$trt=='hcc' & data.serum.long$advanced=='y' & data.serum.long$hbv=='n')),]
	fit.advanced_nonhbvhcc_vs_healthy <- plotMirna(data.tmp, mirna, 'hcc', xlab='Healthy v advanced non-HBV-HCC', ylim=ylim)
	
	if (make.plots)
	{
		pdf(concat('plots/plots-',make.names(mirna),'.pdf'))
		#age.plot <- plot(value ~ age, data.serum.long[which(data.serum.long$mirna==mirna),], main=mirna)
		print(grid.arrange(
				fit.hcc_vs_healthyhbv$plot,
				fit.hcc_vs_healthy$plot,
				fit.tumor_vs_nontumor$plot, 
				fit.early_vs_healthy$plot,
				fit.advanced_vs_healthy$plot,
				fit.hbvhcc_vs_hbv$plot, nrow=2, main='Early and advanced HCC'))
		print(grid.arrange(
				fit.hbvhcc_vs_healthy$plot,
				fit.nonhbvhcc_vs_healthy$plot,		
				fit.early_hbvhcc_vs_healthy$plot,
				fit.early_nonhbvhcc_vs_healthy$plot,
				fit.advanced_hbvhcc_vs_healthy$plot,
				fit.advanced_nonhbvhcc_vs_healthy$plot, nrow=2, main='HBV- and non-HBV-related HCC'))
		print(grid.arrange(	
				fit.nonhbvhcc_vs_hbvhcc$plot,
				fit.nonhcvhcc_vs_hcvhcc$plot,
				fit.hbv_vs_healthy$plot, nrow=2, main='Differences attributable to HBV or HCV'))
		dev.off()
	}
	return(list(
			tumor_vs_nontumor=fit.tumor_vs_nontumor,
			hcc_vs_healthyhbv=fit.hcc_vs_healthyhbv,
			hcc_vs_healthy=fit.hcc_vs_healthy,
			hbvhcc_vs_healthy=fit.hbvhcc_vs_healthy,
			hbvhcc_vs_hbv=fit.hbvhcc_vs_hbv,
			nonhbvhcc_vs_healthy=fit.nonhbvhcc_vs_healthy,
			nonhbvhcc_vs_hbvhcc=fit.nonhbvhcc_vs_hbvhcc,
			nonhcvhcc_vs_hcvhcc=fit.nonhcvhcc_vs_hcvhcc,
			hbv_vs_healthy=fit.hbv_vs_healthy,
			early_vs_healthy=fit.early_vs_healthy,
			advanced_vs_healthy=fit.advanced_vs_healthy, 
			early_hbvhcc_vs_healthy=fit.early_hbvhcc_vs_healthy,
			advanced_hbvhcc_vs_healthy=fit.advanced_hbvhcc_vs_healthy,
			early_nonhbvhcc_vs_healthy=fit.early_nonhbvhcc_vs_healthy,
			advanced_nonhbvhcc_vs_healthy=fit.advanced_nonhbvhcc_vs_healthy))
}
#analyzeMirna(data.paired.long, data.serum.long, mirna='hsa-miR-122', make.plots=TRUE)

getMirnaList <- function(filename='list.txt')
{
	data <- loadDataFrame(filename)
	return(unique(data$mirna))
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
			results[mirna,'hcc_vs_healthyhbv'] <- fits$hcc_vs_healthyhbv$p.value
			results[mirna,'R_hcc_vs_healthyhbv'] <- fits$hcc_vs_healthyhbv$ratio
			results[mirna,'hcc_vs_healthy'] <- fits$hcc_vs_healthy$p.value
			results[mirna,'R_hcc_vs_healthy'] <- fits$hcc_vs_healthy$ratio
			results[mirna,'early_vs_healthy'] <- fits$early_vs_healthy$p.value
			results[mirna,'R_early_vs_healthy'] <- fits$early_vs_healthy$ratio
			results[mirna,'advanced_vs_healthy'] <- fits$advanced_vs_healthy$p.value
			results[mirna,'R_advanced_vs_healthy'] <- fits$advanced_vs_healthy$ratio
			results[mirna,'tumor_vs_nontumor'] <- fits$tumor_vs_nontumor$p.value
			results[mirna,'R_tumor_vs_nontumor'] <- fits$tumor_vs_nontumor$ratio
			results[mirna,'tumor_vs_nontumor'] <- fits$tumor_vs_nontumor$p.value
			results[mirna,'R_tumor_vs_nontumor'] <- fits$tumor_vs_nontumor$ratio
			results[mirna,'hbvhcc_vs_hbv'] <- fits$hbvhcc_vs_hbv$p.value
			results[mirna,'R_hbvhcc_vs_hbv'] <- fits$hbvhcc_vs_hbv$ratio
			
			results[mirna,'hbvhcc_vs_healthy'] <- fits$hbvhcc_vs_healthy$p.value
			results[mirna,'nonhbvhcc_vs_healthy'] <- fits$nonhbvhcc_vs_healthy$p.value
			results[mirna,'early_hbvhcc_vs_healthy'] <- fits$early_hbvhcc_vs_healthy$p.value
			results[mirna,'early_nonhbvhcc_vs_healthy'] <- fits$early_nonhbvhcc_vs_healthy$p.value
			results[mirna,'advanced_hbvhcc_vs_healthy'] <- fits$advanced_hbvhcc_vs_healthy$p.value			
			results[mirna,'advanced_nonhbvhcc_vs_healthy'] <- fits$advanced_nonhbvhcc_vs_healthy$p.value			
			
			results[mirna,'hbv_vs_healthy'] <- fits$hbv_vs_healthy$p.value
			results[mirna,'nonhbvhcc_vs_hbvhcc'] <- fits$nonhbvhcc_vs_hbvhcc$p.value
			results[mirna,'nonhcvhcc_vs_hcvhcc'] <- fits$nonhcvhcc_vs_hcvhcc$p.value
		})
	}
	results <- results[order(
					results$hcc_vs_healthyhbv,
					results$hcc_vs_healthy,
					results$early_vs_healthy,
					results$advanced_vs_healthy,
					results$tumor_vs_nontumor),]
#					results$nonhbvhcc_vs_hbvhcc,
#					results$nonhbvhcc_vs_hbvhcc,
#					results$nonhcvhcc_vs_hcvhcc,
#					results$hbv_vs_healthy
				
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
