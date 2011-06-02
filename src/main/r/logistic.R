library(popbio)
library(coin)
library(rms,T)

#source(paste(libdir,'external.r',sep=''))

univariateLogisticTest <- function(data, response, field, verbose=FALSE, use.wilcox=TRUE)
{
	print(paste('use.wilcox =',use.wilcox))
	cls <- class(data[[field]])
	print(paste(field,cls))
	n <- length(na.omit(data[,field]))
	
	print(paste('univariate test for',response,field))
	lrm.fit <- lrm(makeFormula(response,field), data=data)
	lrm.fit.summary <- summary(lrm.fit)
	
	ci.lower <- format(lrm.fit.summary[' Odds Ratio','Lower 0.95'], digits=3)
	ci.upper <- format(lrm.fit.summary[' Odds Ratio','Upper 0.95'], digits=3)
	or <- lrm.fit.summary[' Odds Ratio','Effect']
	
	if (is.numeric(data[[field]]))
	{
		if (use.wilcox)
		{
			fit <- wilcox.test(as.formula(paste(field," ~ ",response)), data)
			pvalue <- fit$p.value
			return(list(fit=fit, p.value=pvalue, n=n, or=or, ci.lower=ci.lower, ci.upper=ci.upper))
		}
		else
		{
			fit <- anova(lrm.fit)
			pvalue <- fit[field,'P']
			return(list(fit=fit, p.value=pvalue, n=n, or=or, ci.lower=ci.lower, ci.upper=ci.upper))
		}
	}
	else if (is.factor(data[[field]]))
	{
		#fit <- chisq.test(xtabs(as.formula(paste('~', response, '+', field)), data=data))
		#return(list(fit=fit, p.value=fit$p.value, n=n))
		fit <- fisher.test(xtabs(as.formula(paste('~', response, '+', field)), data=data))
		return(list(fit=fit, p.value=fit$p.value, n=n, or=or, ci.lower=ci.lower, ci.upper=ci.upper))
	}
	else throw('no handler for data type',cls)
}
#univariateLogisticTest(data, 'svr', 'ISG15', verbose=T, use.wilcox=T)

univariateLogisticTests <- function(data, response, fields, verbose=TRUE, use.wilcox=T)
{
	fits <- list()
	use.wlcx <- use.wilcox
	vrbs <- verbose
	for (field in splitFields(fields))
	{
		try({
			fit <- univariateLogisticTest(data, response, field, verbose=vrbs, use.wilcox=use.wlcx)
			fits[[field]] <- fit
		}, silent=F)
	}
	return(fits)
}
#uniresults <- univariateLogisticTests(data,'nvr',fields,T)

findSignificantUnivariateLogisticFields <- function(uniresults, cutoff=0.05, verbose=TRUE)
{
	sigfields <- c()
	for (field in names(uniresults))
	{
		uniresult <- uniresults[[field]]
		if (uniresult$p.value<=cutoff)
		{
			sigfields <- append(sigfields,field)
			print(uniresult$fit)
		}
	}
	if (length(sigfields)==0)
		throw(paste("no fields significant at P=",cutoff,'using all fields'))
	if (verbose) print(paste("significant univariate effects: ",sigfields))
	return(sigfields)
}
#findSignificantUnivariateLogisticFields(uniresults)
#
#findBestFields <- function(data,response,sigfields,verbose=TRUE)
#{
#	#data.sigfields <<- na.omit(data[,c(response,sigfields)])
#	#fit <- glm(makeFormula(response,sigfields), data=data.sigfields, na.action=na.fail, family=binomial())
#	fit <- glm(makeFormula(response,sigfields), data, na.action=na.fail, family=binomial())
#	step <- stepAIC(fit, direction='both')
#	bestfields <- attr(step$terms,'term.labels')
#	if (verbose)
#	{
#		print("Model using only significant fields")
#		print(step)		
#		print(bestfields)
#	}
#	if (length(bestfields)==0)
#		throw('No fields retained in final model - using all fields')
#	return(bestfields)
#}
#
#multivariateLogisticTests <- function(data,response,bestfields,penalized=FALSE,verbose=TRUE)
#{
#	bestfit <- lrm(makeFormula(response,bestfields), data=data, x=T, y=T)
#	adjbestfit <-bestfit
#	if (penalized)
#		adjbestfit <- penalizedLogisticRegression(data, response, bestfields)
#	
#	if (verbose) print("Model using only aic fields")
#	if (verbose) print(bestfields)
#	if (verbose) print(summary(adjbestfit))
#	if (verbose) print(anova(adjbestfit))
#	
#	makeLogisticDiagnosticPlots(adjbestfit)
#	
#	#try with glm to get diagnostics
#	glmfit <- glm(makeFormula(response,bestfields), data=data, family=binomial())
#	print('***********************************')
#	print('GLM analysis')
#	print('***********************************')
#	print(summary(glmfit))
#	print(anova(glmfit))
#	print(paste('Overdispersion:',testOverdispersion(glmfit)))
#	print('***********************************')
#	plot(glmfit)
#	return(list(fit=adjbestfit, bestfit=bestfit, adjbestfit=adjbestfit, glmfit=glmfit, bestfields=bestfields))
#}
##


#multivariateLogisticTests <- function(data,response,sigfields,penalized=FALSE,verbose=TRUE, vif.cutoff=NULL)
#{
#	data.sigfields <<- na.omit(data[,c(response,sigfields)])
#	fit <- glm(makeFormula(response,sigfields), data=data.sigfields, na.action=na.fail, family=binomial())
#	step <- stepAIC(fit, direction='both')
#	bestfields <- attr(step$terms,'term.labels')
#	if (verbose)
#	{
#		print("Model using only significant fields")
#		print(step)		
#		print(bestfields)
#	}
#	
#	if (length(bestfields)==0)
#		throw('No fields retained in final model')
#	
#	bestfit <- lrm(makeFormula(response,bestfields), data=data, x=T, y=T)
#	if (!is.null(vif.cutoff))
#	{
#		bestfields <- c()
#		fit.vif <- vif(bestfit)
#		for (field in names(fit.vif))
#		{
#			print(paste(field,'=',fit.vif[field]))
#			if (fit.vif[field] < vif.cutoff)
#				bestfields <- c(bestfields,field)
#		}
#		print(bestfields)
#		bestfit <- lrm(makeFormula(response,bestfields), data=data, x=T, y=T)
#	}
#
#	adjbestfit <-bestfit
#	if (penalized)
#		try(adjbestfit <- penalizedLogisticRegression(data, response, bestfields),silent=F)
#	
#	if (verbose)
#	{
#		print("Model using only aic fields")
#		print(bestfields)
#		print(summary(adjbestfit))
#		print(anova(adjbestfit))
#		makeLogisticDiagnosticPlots(adjbestfit,data)
#	}	
#
#	#try with glm to get diagnostics
#	glmfit <- glm(makeFormula(response,bestfields), data=data, family=binomial())
#	if (verbose)
#	{
#		print('***********************************')
#		print('GLM analysis')
#		print('***********************************')
#		print(summary(glmfit))
#		print(anova(glmfit))
#		print(paste('Overdispersion:',testOverdispersion(glmfit)))
#		print('***********************************')
#		plot(glmfit)
#	}
#	return(list(fit=adjbestfit, bestfit=bestfit, adjbestfit=adjbestfit, glmfit=glmfit, bestfields=bestfields))
#}


findBestFields <- function(data,response,sigfields,verbose=TRUE)
{
	#data.sigfields <<- na.omit(data[,c(response,sigfields)])
	data.sigfields <<- na.omit(data[,c(response,sigfields)])
	fit <- glm(makeFormula(response,sigfields), data=data.sigfields, na.action=na.fail, family=binomial())
	step <- stepAIC(fit, direction='both')
	bestfields <- attr(step$terms,'term.labels')
	if (verbose)
	{
		print("Model using only significant fields")
		print(step)		
		print(bestfields)
	}
	if (length(bestfields)==0)
		throw('No fields retained in final model - using all fields')
	return(bestfields)
}
#findBestFields(data.naive,'svr',fields.genotype1b)

multivariateLogisticTests <- function(data,response,sigfields,usefields=NULL, penalized=FALSE,verbose=TRUE, vif.cutoff=NULL)
{
	if (is.null(usefields))
		bestfields <- findBestFields(data,response,sigfields,verbose)
	else bestfields <- usefields
	
	bestfit <- lrm(makeFormula(response,bestfields), data=data, x=T, y=T)
	if (!is.null(vif.cutoff))
	{
		bestfields <- c()
		fit.vif <- vif(bestfit)
		for (field in names(fit.vif))
		{
			print(paste(field,'=',fit.vif[field]))
			if (fit.vif[field] < vif.cutoff)
				bestfields <- c(bestfields,field)
		}
		print(bestfields)
		bestfit <- lrm(makeFormula(response,bestfields), data=data, x=T, y=T)
	}
	
	adjbestfit <-bestfit
	if (penalized)
		try(adjbestfit <- penalizedLogisticRegression(data, response, bestfields),silent=F)
	
	if (verbose)
	{
		print("Model using only aic fields")
		print(bestfields)
		print(summary(adjbestfit))
		print(anova(adjbestfit))
		makeLogisticDiagnosticPlots(adjbestfit,data)
	}	
	
	#try with glm to get diagnostics
	glmfit <- glm(makeFormula(response,bestfields), data=data, family=binomial())
	if (verbose)
	{
		print('***********************************')
		print('GLM analysis')
		print('***********************************')
		print(summary(glmfit))
		print(anova(glmfit))
		print(paste('Overdispersion:',testOverdispersion(glmfit)))
		print('***********************************')
		plot(glmfit)
	}
	return(list(fit=adjbestfit, bestfit=bestfit, adjbestfit=adjbestfit, glmfit=glmfit, bestfields=bestfields))
}

createLogisticRegressionTable <- function(uniresults, multiresults, response)
{	
	table <- data.frame(variable=names(uniresults), row.names=names(uniresults), stringsAsFactors=FALSE)
	for (field in names(uniresults))
	{
		uniresult <- uniresults[[field]]		
		table[field,'N1'] <- uniresult$n
		table[field,'OR1'] <- format(uniresult$or, digits=3)
		table[field,'CI1'] <- paste('(',uniresult$ci.lower,'-',uniresult$ci.upper,')', sep='') 
		table[field,'P1'] <- format(uniresult$p.value, digits=4)
		table[field,'sig1'] <- createSignificantMarker(uniresult$p.value)
	}
	
	multiaov <- anova(multiresults$fit);
	numrows <- nrow(multiaov)-1
	chisquares <- multiaov[1:numrows,'Chi-Square']
	pvalues <- multiaov[1:numrows,'P']
	#takes care of case with only one variable
	rownames <- row.names(multiaov)[1:numrows]
	names(chisquares) <- rownames
	names(pvalues) <- rownames
	n2 <- multiresults$fit$stats['Obs'][1]	
	multiresults.summary <- summary(multiresults$fit)
	print(multiresults.summary)
	
	for (field in names(chisquares))
	{
		try({
		summary.rownum <- grep(field,rownames(multiresults.summary))+1
		ci.lower <- format(multiresults.summary[summary.rownum,'Lower 0.95'], digits=3)
		ci.upper <- format(multiresults.summary[summary.rownum,'Upper 0.95'], digits=3)
		or2 <- multiresults.summary[summary.rownum,'Effect']
		print(paste('odds ratio for field',field,'=',or2))
		table[field,'N2'] <- n2
		table[field,'OR2'] <- format(or2, digits=3)
		table[field,'CI2'] <- paste('(',ci.lower,'-',ci.upper,')', sep='') 
		table[field,'P2'] <- format(pvalues[field], digits=4)
		table[field,'sig2'] <- createSignificantMarker(pvalues[field])
		},silent=T)
	}
	table <- fixVariableNames(table)
	table <- replaceNAs(table)
	return(table)
}

writeLogisticRegressionTable <- function(dataframe, response, filename=NULL)
{
	if (is.null(filename))
		filename <- paste('out/table-',response,'.txt',sep='')
	#rownames(dataframe) <- rep('',length(rownames(dataframe)))
	print(dataframe)
	print(paste('Writing dataframe to ',getwd(),'/',filename, sep=''))
	colnames <- c('Variable','n','OR','(95% CI)','P','','n','OR','(95% CI)','P','')
	write.table(dataframe, filename, quote=FALSE, row.names = FALSE, col.names=colnames, sep='\t', na = '')
}

multipleLogisticRegression <- function(data, response, fields, usefields=NULL, cutoff=0.05, vif.cutoff=NULL, penalized=TRUE, verbose=TRUE, bestfields=NULL, filename='', use.wilcox=T)
{
	#for each factor do a univariate non-parametric test
	uniresults <- univariateLogisticTests(data, response, fields, verbose=verbose, use.wilcox=use.wilcox)
	sigfields <- findSignificantUnivariateLogisticFields(uniresults,cutoff)
	#if (is.null(bestfields))
	#	bestfields <- findBestFields(data,response,sigfields,verbose)
	multiresults <- multivariateLogisticTests(data,response=response,sigfields=sigfields,usefields=usefields,penalized=penalized,vif.cutoff,verbose=verbose)
	table <- createLogisticRegressionTable(uniresults,multiresults)
	results <- list(uniresults=uniresults, sigfields=sigfields, bestfields=multiresults$bestfields, multiresults=multiresults, table=table)
	writeLogisticRegressionTable(table,response,filename)
	return(results)
}
#results <- multipleLogisticRegression(data,'vr',fields, cutoff=cutoff, penalized=T, verbose=F, 'out/table-vr-penalized.txt')


makeLogisticPlot <- function(data, response, predictor)
{
	#library(popbio)
	data.subset <- na.omit(data[,c(response,predictor)])
	logi.hist.plot(data.subset[[predictor]], data.subset[[response]], type = "dit", boxp = TRUE, rug = TRUE, mainlabel=response, xlabel=predictor)	
}
#makeLogisticPlot(data,'nvr','ISG15_trimmed')

getR2shrinkage <- function(fit)
{
	v <- validate(fit, B=140)
	optimism <- v['R2','optimism']
	return(optimism)
}

makeLogisticDiagnosticPlots <- function(fit, data)
{
	#variance inflation factor
	print(vif(fit))
	
	#influential points
	try({
		inf <- which.influence(fit, cutoff=.2)
		show.influence(inf, data, report='no')
	})
	
	#r2 shrinkage
	print(paste('R2 shrinkage:',getR2shrinkage(fit)))
	
	#fit
	print(resid(fit, 'gof'))
	
	par(ask=T)	
	try({resid(fit, 'partial', pl=TRUE)}) #pl=TRUE pl='loess'
	
	plot(summary(fit), log=T)
	plot(anova(fit))
	
	try(plot.lrm.partial(fit),silent=T)
	
	#plot(calibrate(fit))
	
	par(ask=F)
}

penalizedLogisticRegression <- function(data, response, bestfields, penlty=NULL)
{
	if (is.null(penlty))
		penlty <- findOptimalPenalty(data,response,bestfields)
	fit <- lrm(makeFormula(response,bestfields), data, x=T, y=T)
	penalty.matrix <- diag(diag(var(fit$x)))
	fit2 <- lrm(makeFormula(response,bestfields), data, penalty=penlty, penalty.matrix=penalty.matrix, x=T,y=T)
	val2 <- validate(fit2, method="boot", B=300)
	print(paste('penalty:',penlty))
	print(fit)
	print(summary(fit))
	print(fit2)
	print(summary(fit2))
	print(val2)
	return(fit2)
}
#penalizedLogisticRegression(data,'vr', splitFields('IL28,MxA,OAS1,ISG15'), 0.1)

findOptimalPenalty <- function(data, response, fields, verbose=F)
{
	data <- na.omit(data[,c(response,fields)])
	penlty <- 0
	penalties <- c()
	rocs <- c()
	briers <- c()
	deviances <- c()
	for(penlty in seq(0,.15, by=.005))
	{
		if (penlty==0)
		{
			fit <- lrm(makeFormula(response,fields), data, x=T, y=T)
			#model <- as.formula(paste(response," ~ as.numeric(",paste(fields, collapse=") + as.numeric("),')'))
			#print(model)
			#fit <- lrm(model, data, x=T, y=T)
			
			X <- fit$x
			Y <- fit$y
			penalty.matrix <- diag(diag(var(X)))
			Xnew <- predict(fit, data, type="x", incl.non.slopes=FALSE)
			Ynew <- data$vr
		}
		else
		{
			fit <- lrm.fit(X,Y, penalty.matrix=penlty*penalty.matrix)
			pred.logit <- fit$coef[1] + (Xnew %*% fit$coef[-1])
			pred <- plogis(pred.logit)
			C.index <- somers2(pred, Ynew)["C"]
			Brier <- mean((pred-Ynew)^2)
			Deviance <- -2*sum( Ynew*log(pred) + (1-Ynew)*log(1-pred) )
			#cat("\nPenalty :",penlty,"\n")
			#cat("ROC area:",format(C.index)," Brier score:",format(Brier), " -2 Log L:",format(Deviance),"\n")
			
			penalties <- c(penalties,penlty)
			rocs <- c(rocs,C.index)
			briers <- c(briers,Brier)
			deviances <- c(deviances,Deviance)
		}
	}
	results <- data.frame(penalty=penalties, roc=rocs, brier=briers, deviance=deviances)
	results <- results[order(results$brier),]
	if (verbose)
	{
		print(results)
		#plot(results$brier ~ 1:(nrow(results)-1))
	}
	return(results[1,'penalty'])
}
#findOptimalPenalty(data,'vr',splitFields('IL28,MxA,OAS1,ISG15,fibrosis'),verbose=T)



plot.summary.rms2 <- function (x, at, log = FALSE, q = c(0.7, 0.8, 0.9, 0.95, 0.99), 
		xlim, nbar, cex = 1, nint = 10, cex.c = 0.5, cex.t = 1, clip = c(-1e+30, 
				1e+30), main, fieldnames=NULL, ...) 
{
	scale <- attr(x, "scale")
	adjust <- attr(x, "adjust")
	Type <- x[, "Type"]
	x <- x[Type == 1, , drop = FALSE]
	lab <- dimnames(x)[[1]]
	if (!is.null(fieldnames))
		lab=fieldnames
	effect <- x[, "Effect"]
	se <- x[, "S.E."]
	if (!log && any(Type == 2)) {
		fun <- exp
		tlab <- scale[2]
	}
	else {
		fun <- function(x) x
		if (log) {
			if (length(scale) == 2) 
				tlab <- scale[2]
			else tlab <- paste("exp(", scale[1], ")", sep = "")
		}
		else tlab <- scale[1]
	}
	if (!length(scale)) 
		tlab <- ""
	if (!missing(main)) 
		tlab <- main
	augment <- if (log | any(Type == 2)) 
				c(0.1, 0.5, 0.75, 1)
			else 0
	n <- length(effect)
	out <- qnorm((max(q) + 1)/2)
	if (missing(xlim) && !missing(at)) 
		xlim <- range(if (log) logb(at) else at)
	else if (missing(xlim)) {
		xlim <- fun(range(c(effect - out * se, effect + out * 
										se)))
		xlim[1] <- max(xlim[1], clip[1])
		xlim[2] <- min(xlim[2], clip[2])
	}
	else augment <- c(augment, if (log) exp(xlim) else xlim)
	fmt <- function(k) {
		m <- length(k)
		f <- character(m)
		for (i in 1:m) f[i] <- format(k[i])
		f
	}
	lb <- ifelse(is.na(x[, "Diff."]), lab, lab)
	#lb <- ifelse(is.na(x[, "Diff."]), lab, paste(lab, " - ", 
	#    fmt(x[, "High"]), ":", fmt(x[, "Low"]), sep = ""))
	plot.new()
	par(new = TRUE)
	mxlb <- 0.1 + max(strwidth(lb, units = "inches", cex = cex))
	tmai <- par("mai")
	on.exit(par(mai = tmai))
	par(mai = c(tmai[1], mxlb, 1.5 * tmai[3], tmai[4]))
	outer.widths <- fun(effect + out * se) - fun(effect - out * 
					se)
	if (missing(nbar)) 
		nbar <- n
	npage <- ceiling(n/nbar)
	is <- 1
	for (p in 1:npage) {
		ie <- min(is + nbar - 1, n)
		plot(1:nbar, rep(0, nbar), xlim = xlim, ylim = c(1, nbar), 
				type = "n", axes = FALSE, xlab = "", ylab = "")
		if (cex.t > 0) 
			title(tlab, cex = cex.t)
		lines(fun(c(0, 0)), c(nbar - (ie - is), nbar), lty = 2)
		if (log) {
			pxlim <- pretty(exp(xlim), n = nint)
			pxlim <- sort(unique(c(pxlim, augment)))
			pxlim <- pxlim[pxlim >= exp(xlim[1])]
			if (!missing(at)) 
				pxlim <- at
			axis(3, logb(pxlim), lab = format(pxlim))
		}
		else {
			pxlim <- pretty(xlim, n = nint)
			pxlim <- sort(unique(c(pxlim, augment)))
			pxlim <- pxlim[pxlim >= xlim[1]]
			if (!missing(at)) 
				pxlim <- at
			axis(3, pxlim)
		}
		imax <- (is:ie)[outer.widths[is:ie] == max(outer.widths[is:ie])][1]
		for (i in is:ie) {
			confbar(nbar - (i - is + 1) + 1, effect[i], se[i], 
					q = q, type = "h", fun = fun, cex = cex.c, labels = i == 
							imax, clip = clip, ...)
			mtext(lb[i], 2, 0, at = nbar - (i - is + 1) + 1, 
					cex = cex, adj = 1, las = 1)
		}
		if (adjust != "") {
			adjto <- paste("Adjusted to:", adjust, sep = "")
			xx <- par("usr")[2]
			if (nbar > ie) 
				text(xx, nbar - (ie - is + 1), adjto, adj = 1, 
						cex = cex)
			else title(sub = adjto, adj = 1, cex = cex)
		}
		is <- ie + 1
	}
	invisible()
}
#plot.summary.rms2

########################################################33

predictY <- function(data, coeff, fields, hold=c())
{
	y <- coeff[1]
	for (field in fields)
	{
		if (containsElement(hold,field))
		{
			print(paste('Skipping field',field))
			y <- y + coeff[field]
		}
		else y <- y + coeff[field]*data[[field]]
	}
	y <- 1/(1+exp(-1*y))
	return(y)
}
#y2 <- predictY2(data.naive,coeff,fields.svr)
#subsetNA(data.frame(y=y2, svr=data.naive$svr))


getPerformanceMeasure <- function(pred,name)
{
	perf <- performance(pred,name)
	value <- perf@y.values[[1]]
	print(paste(name,'=',value))
	return(value)
}

plotROC <- function(data, response, fields, ...)
{
	fields <- splitFields(fields)
	fit <- glmLogisticRegression(data,response,fields)
	fit.sum <- summary(fit)
	coeff <-fit.sum$coefficients[,1]
	
	data.rocr <- data.frame(predictions=predictY(data,coeff,fields), labels=data[[response]])
	data.rocr <- subsetNA(data.rocr)
	#print(head(data.rocr))
	pred <- prediction(data.rocr$predictions, data.rocr$labels)
	#perf <- performance(pred,"auc")
	#auc <- perf@y.values[[1]]
	
	auc <- getPerformanceMeasure(pred,'auc')
#	tpr <- getPerformanceMeasure(pred,'tpr')
#	tpr <- getPerformanceMeasure(pred,'fpr')
#	ppv <- getPerformanceMeasure(pred,'ppv')
#	npv <- getPerformanceMeasure(pred,'npv')
#	acc <- getPerformanceMeasure(pred,'acc')
	
	perf <- performance(pred,"tpr","fpr")
	plot(perf,colorize=F, ...) #sub=paste('AUC =',format(auc,digits=2))
	return(auc)
	#return(data.rocr)
}
#plotROC(data,'svr','age,bmi,X1B,logviralload,logGPT,rs8099917,rs1127354,rbv,pegifn', title='SVR - all patients')

#########################################################################


simpleUnivariateLogisticTest <- function(data, response, field, alpha=0.05)
{
	cls <- class(data[[field]])
	#print(paste(field,cls))
	n <- length(na.omit(data[,field]))
	
	#print(paste('univariate test for',response,field))
	
	fit <- glm(makeFormula(response,field), family=binomial, data=data)
	fit.sum <- summary(fit)
	coeff <- fit.sum$coefficients[,1]
	or <- exp(coeff[[field]])
	conf.int <- confint(fit, level=1-alpha)
	ci.lower <- format(exp(conf.int[2,1]), digits=3)
	ci.upper <- format(exp(conf.int[2,2]), digits=3)
	
	if (is.numeric(data[[field]]))
	{
		fit <- wilcox.test(as.formula(paste(field," ~ ",response)), data)
		return(list(fit=fit, p.value=fit$p.value, n=n, or=or, ci.lower=ci.lower, ci.upper=ci.upper))
	}
	else if (is.factor(data[[field]]))
	{
		fit <- fisher.test(xtabs(as.formula(paste('~', response, '+', field)), data=data))
		return(list(fit=fit, p.value=fit$p.value, n=n, or=or, ci.lower=ci.lower, ci.upper=ci.upper))
	}
	else throw('no handler for data type',cls)
}
#simpleUnivariateLogisticTest(data, 'svr', 'male1')

simpleUnivariateLogisticTests <- function(data, response, fields)
{
	fits <- list()
	for (field in fields)
	{
		try({
					fit <- simpleUnivariateLogisticTest(data, response, field)
					fits[[field]] <- fit
				}, silent=F)
	}
	return(fits)
}
#uniresults <- simpleUnivariateLogisticTests(data.naive,'svr',fields.svr)

glmLogisticRegression <- function(data, response, fields, filename=NULL, alpha=0.05)
{
	uniresults <- simpleUnivariateLogisticTests(data, response, fields)
	
	table <- data.frame(variable=fields, row.names=fields, stringsAsFactors=FALSE)
	for (field in fields)
	{
		print(field)
		try({
					uniresult <- uniresults[[field]]		
					table[field,'N1'] <- uniresult$n
					table[field,'OR1'] <- format(uniresult$or, digits=3)
					table[field,'CI1'] <- paste('(',uniresult$ci.lower,'-',uniresult$ci.upper,')', sep='') 
					table[field,'P1'] <- format(uniresult$p.value, digits=4)
					table[field,'sig1'] <- createSignificantMarker(uniresult$p.value)
				}, silent=FALSE)
	}
	
	fit <- glm(makeFormula(response,fields), data=data, family=binomial)
	fit.sum <- summary(fit)
	#fit.aov <- anova(fit)
	
	coeff <- fit.sum$coefficients[,1]
	ors <- exp(coeff)
	conf.int <- confint(fit, level=1-alpha)
	ci.lowers <- exp(conf.int[,1])
	ci.uppers <- exp(conf.int[,2])	
	n <- length(fit$fitted.values)
	
	for (field in fields)
	{
		try({
					or <- ors[[field]]
					ci.lower <- format(ci.lowers[[field]], digits=3)
					ci.upper <- format(ci.uppers[[field]], digits=3)
					p.value <- fit.sum$coefficients[field,'Pr(>|z|)']
					est <- fit.sum$coefficients[field,'Estimate']
					
					#print(paste('odds ratio for field',field,'=',or))
					
					table[field,'N2'] <- n
					table[field,'Est2'] <- format(est, digits=4)
					table[field,'OR2'] <- format(or, digits=3)
					table[field,'CI2'] <- paste('(',ci.lower,'-',ci.upper,')', sep='') 
					table[field,'P2'] <- format(p.value, digits=4)
					table[field,'sig2'] <- createSignificantMarker(p.value)
				},silent=F)
	}
	table <- fixVariableNames(table)
	table <- replaceNAs(table)
	
	if (is.null(filename))
		filename <- paste('out/table-',response,'.txt',sep='')
	
	if (!is.null(filename))
	{
		print(table)
		print(paste('Writing dataframe to ',getwd(),'/',filename, sep=''))
		colnames <- c('Variable','n','OR','(95% CI)','P','','n','Coeff','OR','(95% CI)','P','')
		write.table(table, filename, quote=FALSE, row.names = FALSE, col.names=colnames, sep='\t', na = '')
	}
	return(fit)
}

