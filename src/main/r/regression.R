library(gvlma)

findSignificantUnivariateFields <- function(uniresults, cutoff=0.05)
{
	sigfields <- c()
	for (field in names(uniresults))
	{
		fit <- uniresults[[field]]
		sum <- summary(fit)
		pvalue <- sum$coefficients[,'Pr(>|t|)'][2]
		if (pvalue<=cutoff)
		{
			sigfields <- append(sigfields,field)
			print(summary(fit))
		}
	}
	print(paste('significant univariate effects: ',sigfields))
	return(sigfields)
}

univariateTest <- function(data, model, showSummary=FALSE)
{
	frmla <- as.formula(model)
	print(frmla)
	fit <- lm(frmla, data=data)
	#try(fit <- lm(as.formula(model), data=data), silent=TRUE)
	if (showSummary)
		print(summary(fit))
	return(fit)
}

univariateTests <- function(data, response, fields, showSummary=FALSE)
{
	fits <- list()
	for (field in fields)
	{
		try(fits[[field]] <- univariateTest(data, paste(response,' ~ ',field), showSummary))
	}
	return(fits)
}

multivariateTests <- function(data,response,sigfields,trace=0)
{
	#data.subset <- data[,c(response,sigfields)]
	#data.subset <- na.omit(data.subset)
	data.sigfields <<- na.omit(data[,c(response,sigfields)])
	#fit <- lm(makeFormula(response,sigfields), data=data) #data.subset)
	fit <- lm(makeFormula(response,sigfields), data=data.sigfields) #data.subset)
	print('Model using only significant fields')
	print(summary(fit))
	
	#subset <- remove.NA(subset)
	#fit <- lm(makeFormula(response,sigfields), data=subset)
	#aic <- stepAIC(fit, direction='backward', trace=trace)
	aic <- stepAIC(fit, direction='both', trace=0)#, trace=trace)
	#bestfields <- names(aic$coefficients)[-1]
	bestfields <- attr(aic$terms,'term.labels')
	print(bestfields)
	
	bestfit <- lm(makeFormula(response,bestfields), data=data) #data.subset
	print('Model using only aic fields')
	print(summary(bestfit))
	return(bestfit)	
}
#multivariateTests(kawaoka.2b,'titerlog4',splitFields('age,sex,rs12979860CC,rs8099917TT,rs12980275AA,titerlog0'))

createRegressionTable <- function(uniresults, multiresults, response)
{	
	table <- data.frame(variable=names(uniresults),row.names=names(uniresults), stringsAsFactors=FALSE)
	for (field in names(uniresults))
	{
		unisum <- summary(uniresults[[field]])
		pvalue <- unisum$coefficients[,'Pr(>|t|)'][2]
		table[field,'N1'] <- length(uniresults[[field]]$fitted.values)
		table[field,'Estimate1'] <- format(unisum$coefficients[,'Estimate'][2], digits=3)
		table[field,'P1'] <- format(pvalue, digits=4)
		table[field,'sig1'] <- createSignificantMarker(pvalue)
	}
	
	sum <- summary(multiresults);
	numrows <- nrow(sum$coefficients)-1
	estimates <- sum$coefficients[,'Estimate'][1:numrows+1]
	pvalues <- sum$coefficients[,'Pr(>|t|)'][1:numrows+1]
	n <- length(multiresults$fitted.values)
	
	for (field in names(estimates))
	{
		table[field,'N2'] <- n
		table[field,'Estimate2'] <- format(estimates[field], digits=3)
		table[field,'P2'] <- format(pvalues[field], digits=4)
		table[field,'sig2'] <- createSignificantMarker(pvalues[field])
	}
	table <- fixVariableNames(table)
	table <- replaceNAs(table,'N2')
	table <- replaceNAs(table,'Estimate2')
	table <- replaceNAs(table,'P2')
	table <- replaceNAs(table,'sig2')
	return(table)
}

writeRegressionTable <- function(dataframe, response, filename='')
{
	dataframe <- dataframe[!is.na(dataframe$Estimate1),];
	dataframe <- format(dataframe, scientific=FALSE, digits=1, nsmall=3)
	if (filename=='')
		filename <- paste('table-',response,'.txt',sep='')	
	#filename <- paste('table-',response,'.txt',sep='')	
	print(dataframe)
	print(paste('Writing dataframe to ',getwd(),'/',filename, sep=''))	
	colnames <- c('Variable','N','Coefficient','P','','N','Coefficient','P','')
	write.table(dataframe, filename, quote=FALSE, row.names = FALSE, col.names=colnames, sep='\t', na = '')
}

multipleRegression <- function(data,response,fields, cutoff=0.05, filename='')
{
	uniresults <- univariateTests(data, response, fields)
	sigfields <- findSignificantUnivariateFields(uniresults,cutoff)
	multiresults <- multivariateTests(data,response,sigfields)
	table <- createRegressionTable(uniresults,multiresults)
	results <- list(uniresults=uniresults, sigfields=sigfields, multiresults=multiresults, table=table)
	writeRegressionTable(table,response,filename)
	return(results)
}


findSignificantNonparametricFields <- function(data, response, fields, cutoff=0.05)
{
	fits <- list()
	for (field in fields)
	{
		fit <- kruskal.test(as.formula(paste(response,' ~ ',field)), data=data)
		fits[[field]] <- fit
		if (fit$p.value <= cutoff)
			print(fit)
	}
	#sigfields <- c()
	#for (field in names(fits))
	#{
	#fit <- fits[[field]]
	#pvalue <- fit$p.value
	#if (pvalue<=cutoff)
	#{
	#	sigfields <- append(sigfields,field)
	#	print(fit)
	#}
	#}
	#print(paste('significant univariate effects: ',sigfields))
	#return(sigfields)
	#return(fits)
}
#findSignificantNonparametricFields(data.1b, 'wbc01', wbcfields)

#from R in Action
##########################################################
# This function determines the relative importance of each independent #
# variable to the dependent variable in an OLS regression. The code is #
# adapted from an SPSS program generously provided by Dr. Johnson. #
# See Johnson (2000, Multivariate Behavioral Research, 35, 1-19) for #
# an explanation of how the relative weights are derived. #
########################################################################
relweights <- function(fit,...)
{
	R <- cor(fit$model) # correlation matrix with criterion in
	nvar <- ncol(R) # number of variables
	rxx <- R[2:nvar, 2:nvar] # correlations among predictors
	rxy <- R[2:nvar, 1] # correlations between predictors and criterion
	svd <- eigen(rxx) # singular value decomposition
	evec <- svd$vectors # eigenvectors
	ev <- svd$values # eigenvalues
	delta <- diag(sqrt(ev)) # diag matrix with sqrts of eigenvalues
	# correlations between original predictors and new orthogonal variables
	lambda <- evec %*% delta %*% t(evec)
	lambdasq <- lambda ^ 2 # square of the correlations
	# regression coefficients of Y on orthogonal variables
	beta <- solve(lambda) %*% rxy
	rsquare <- colSums(beta^2)
	rawwgt <- lambdasq %*% beta ^ 2	# raw relative weights
	import <- (rawwgt / rsquare) * 100 # rescale to % of R-square
	lbls <- names(fit$model[2:nvar]) # predictor labels
	rownames(import) <- lbls
	colnames(import) <- "Weights"
	# plot results
	barplot(t(import),names.arg=lbls, ylab="% of R-Square", xlab="Predictor Variables", main="Relative Importance of Predictor Variables", sub=paste("R-Square = ", round(rsquare, digits=3)), ...)
	return(import)
}

residplot <- function(fit, nbreaks=10)
{
	x <- rstudent(fit)
	h <- hist(x, breaks=nbreaks, freq=FALSE, xlab='Studentized Residual', main='Distribution of Errors')
	xfit <- seq(min(x), max(x), length=40)
	lines(xfit, dnorm(xfit), col='blue', lwd=2)
	lines(density(x)$x, density(x)$y, col='red', lwd=2, lty=2)
	legend('topright', legend=c('Normal Curve', 'Density Curve'), lty=1:2, col=c('blue','red'), cex=.7)
}

findInfluentialObservations <- function(data, fit)
{
	cutoff <- 4/(nrow(data)-length(fit$coefficients)-2)
	plot(fit, which=4, cook.levels=cutoff)
	abline(h=cutoff, lty=2, col='red')
}

plotHatValues <- function(data,fit)
{
	p <- length(coefficients(fit))
	n <- length(fitted(fit))
	plot(hatvalues(fit), main='Index plot of hat values')
	abline(h=c(2,3)*p/n, col='red', lty=2)
	identify(1:n, hatvalues(fit), row.names(data))
}

regressionDiagnostics <- function(fit)
{
	print(summary(fit))
	oldpar <- par(mfrow=c(2,2), ask=T)
	plot(fit)
	par(mfrow=c(1,1))
	qqPlot(fit, simulate=TRUE)
	#plot studentized residuals
	residplot(fit)
	#test auto-correlation
	print(durbinWatsonTest(fit))
	#test linearity
	print(crPlots(fit, one.page=TRUE, ask=FALSE))
	#test homoscedasticity - non-constant error variance
	print(ncvTest(fit))
	spreadLevelPlot(fit)
	gvmodel <- gvlma(fit)
	print(summary(gvmodel))
	print(sqrt(vif(fit)) > 2) #problem?
	#test outliers
	print(outlierTest(fit))
	par(oldpar)
}

makeScatterplotMatrix <- function(data,fields,groups='')
{
	#data.subset <- na.omit(data[,splitFields('nvr,ISG15_trimmed,IL28,MxA')])
	data.subset <- na.omit(data[,c(splitFields(fields),splitFields(groups))])
	model <- paste('~', paste(splitFields(fields),collapse='+'))
	if (groups!='')
		model <- paste(model, '|', paste(splitFields(groups),collapse='+'))
	scatterplotMatrix(as.formula(model), data=data.subset)
	#scatterplotMatrix(as.formula(~ ISG15_trimmed + IL28 + MxA | nvr, data=data.subset)
}
#makeScatterplotMatrix(data,'ISG15_trimmed,IL28,MxA','nvr')

#test overdispersion
testOverdispersion <- function(fit)
{
	pvalue <- pchisq(summary(fit)$dispersion * fit$df.residual, fit$df.residual, lower=F)
}

### The following is a function adapted from http//www.math.mcmaster.capeters4f03s4f03_0607index.html
### roc.plot() will plot the ROC curve given two vectors of scores, the first for the treatment group (y==1) and the second for the control group (y==0).
roc.plot <- function (sd, sdc, newplot = TRUE, ...) 
{
	sall <- sort(c(sd, sdc))
	sens <- 0
	specc <- 0
	for (i in length(sall):1) {
		sens <- c(sens, mean(sd >= sall[i], na.rm = T))
		specc <- c(specc, mean(sdc >= sall[i], na.rm = T))
	}
	if (newplot) {
		plot(specc, sens, xlim = c(0, 1), ylim = c(0, 1), type = "l", 
				xlab = "1-specificity", ylab = "sensitivity", ...)
		abline(0, 1)
	}
	else lines(specc, sens, ...)
	npoints <- length(sens)
	area <- sum(0.5 * (sens[-1] + sens[-npoints]) * (specc[-1] - 
						specc[-npoints]))
	lift <- (sens - specc)[-1]
	cutoff <- sall[lift == max(lift)][1]
	sensopt <- sens[-1][lift == max(lift)][1]
	specopt <- 1 - specc[-1][lift == max(lift)][1]
	list(area = area, cutoff = cutoff, sensopt = sensopt, specopt = specopt)
}
