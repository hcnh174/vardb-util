#install.packages('mice')

#library(Design)
library(car)
library(MASS)
library(lattice)
#library(rcom)
#library(agce)
#library(leaps)
#library(corrgram)
#library(vcd)
library(R.oo)
#library(reshape)
#library(doBy)

#library(gplot)
#library(HH)
#library(rrcov)

options(contrasts=c("contr.sum","contr.poly"))

loadUtilFiles <- function(filenames)
{
	libdir <- Sys.getenv("VARDB_RUTIL_HOME")
	for (filename in splitFields(filenames))
	{
		source(paste(libdir,filename,'.r',sep=''))
	}
}

setCurDir <- function(dir)
{
	setwd(paste(Sys.getenv("ANALYSIS_HOME"),'/',dir,'/',sep=''))
}

loadDataframe <- function(filename, stringsAsFactors=default.stringsAsFactors())
{
	dataframe <- read.table(filename, header=TRUE, encoding='UTF-8', sep = '\t', comment.char='', stringsAsFactors=stringsAsFactors)
	return(dataframe)
}

#plotDistribution <- function(values, ylab='values')
#{	
#	values <- values[!is.na(values)];
#	print(outlier.test(lm(values ~ 1)))
#	oldpar <- par(mfrow=c(1,2))
#	qqnorm(values, ylab=ylab)
#	qqline(values)
#	hist(values)
#	par(oldpar)
#	#if (min(values>0))
##		print(powerTransform(values))
#	print(sort(values))
#}

# check normality
plotDistribution <- function(values, ylab='values')
{	
	values <- values[!is.na(values)]
	print(sort(values))
	log.values <- log(values)
	print(outlierTest(lm(values ~ 1)))
	oldpar <- par(mfrow=c(1,3))
	qqnorm(values, ylab=ylab)
	try(qqline(values),silent=T)
	hist(values)
	qqnorm(log.values, ylab=ylab)
	try(qqline(log.values),silent=T)
	par(mfrow=c(1,1))
	if (min(values>0))
		print(powerTransform(values))
}

plotDistributions <- function(dataframe, fields)
{
	oldpar <- par(ask=T)
	for (field in splitFields(fields))
	{
		print('******************************')
		print(field)
		plotDistribution(dataframe[,field], field)
	}
	par(oldpar)
}

# automate regression tests

splitFields <- function(str)
{
	if (length(str)==0)
		return(c())
	if (length(str)>1)
		return(str)
	strsplit(str,",")[[1]]
}

makeFormula <- function(responses,fields)
{
	as.formula(paste(paste(responses, collapse="+")," ~ ",paste(fields, collapse="+")))
}


#traceback()

countObservations <- function(data)
{
	print('Number of observations by variable')
	for (name in names(data))
	{
		obs <- length(which(!is.na(data[[name]])))
		#num.rows <- nrow(data)
		print(paste(name,': ',obs,sep=''))
	}	
}

countObservations <- function(data)
{
	table <- data.frame(variable=names(data),row.names=names(data), stringsAsFactors=FALSE) #
	for (name in names(data))
	{
		obs <- length(which(!is.na(data[[name]])))
		table[name,'Count'] <- obs
	}
	
	table <- table[order(table$Count),]
	print('Number of observations by variable')
	#print(paste(name,': ',obs,sep=''))
	print(table)
}
#countObservations(data)

#replaceNAs <- function(table,colname,replacestr='')
#{
#	ids <- which(is.na(table[[colname]]))
#	if (length(ids)>0)
#		table[ids,colname] <- replacestr
#	return(table)
#}

replaceNAs <- function(table,colnames=names(table),replacestr='')
{
	for (colname in colnames)
	{
		ids <- which(is.na(table[[colname]]))
		if (length(ids)>0)
			table[ids,colname] <- replacestr
	}
	return(table)
}

addFieldName <- function(fieldnames, fields, fieldname)
{
	for (field in splitFields(fields))
	{
		fieldnames[[field]] <- fieldname
	}
	return(fieldnames)
}

addFieldSuffix <- function(field, suffix)
{
	return(paste(field,suffix,sep='_'))
}

stripFieldSuffix <- function(field)
{
	return(strsplit(field,"_")[[1]][1])
}

stripFieldSuffixes <- function(fields)
{
	newfields <- c()
	for (field in fields)
	{
		newfields <- c(newfields,stripFieldSuffix(field))
	}
	return(newfields)
}

wrapText <- function(str, prefix, suffix)
{
	return(paste(prefix,str,suffix,sep=''))
}

makeLogField <- function(field)
{
	return(wrapText(field,'log(',')'))
}

fixVariableName <- function(table, rowname, newname, colname='variable')
{
	rownames <- c(rowname)
	rownames <- c(rownames,addFieldSuffix(rowname,'trans'))
	rownames <- c(rownames,addFieldSuffix(rowname,'midcut'))
	rownames <- c(rownames,addFieldSuffix(rowname,'bestcut'))
	rownames <- c(rownames,addFieldSuffix(rowname,'adj'))#
	rownames <- c(rownames,addFieldSuffix(rowname,'adj2'))#	
	ids <- which(table[[colname]] %in% rownames)
	if (length(ids)>0)
		table[ids,colname] <- newname
	return(table)
}

#empty placeholder funtion - redefine later in script to add functionality 
fixVariableNames <- function(table)
{
	return(table)
}

subsetNA <- function(dataframe, fields=NULL)
{
	if (is.null(fields))
		fields <- colnames(dataframe)
	for (field in fields)
	{
		dataframe <- dataframe[which(!is.na(dataframe[[field]])),]
	}
	dataframe <- dataframe[,fields]
	return (dataframe)
}

createSplitField <- function(data, field, newfield, cutoff)
{
	data[[newfield]] <- NA
	data[which(data[[field]] < cutoff),newfield] <- 0
	data[which(data[[field]] >= cutoff),newfield] <- 1
	data[[newfield]] <- factor(data[[newfield]])
	#levels(data[[newfield]]) <- c(paste(newfield,'<',cutoff,sep=''), paste(newfield,'>=',cutoff,sep=''))
	return(data)
}

createSplitFieldByMedian <- function(data, field)
{
	cutoff <- median(data[[field]], na.rm=T)
	print(paste('median cutoff for field',field,':',cutoff))
	newfield <- addFieldSuffix(field,'midcut')# paste(field,'_midcut',sep='')
	#newfield <- paste(field,'_cut',sep='')
	#print(paste('cutoff=',cutoff))
	#print(paste('newfield',newfield))
	data <- createSplitField(data,field,newfield,cutoff)
	return(data)
}
#
#findBestCutpoint <- function(data, response, field)
#{
#	data.subset=subsetNA(data,c(response,field))
#	newfield <- paste(field,'_cut2',sep='')
#	data.subset[,newfield] <- NA
#	for (rowname in row.names(data.subset))
#	{
#		cutpoint <- data.subset[rowname,field]
#		data.subset[,newfield] <- NA
#		data.subset[,newfield] <- ifelse(data.subset[[field]] <= cutpoint, 0, 1)
#		if (length(unique(data.subset[,newfield]))>1)
#		{
#			counts <- xtabs(as.formula(paste('~', response, '+', newfield)), data=data.subset)
#			fit <- chisq.test(counts)
#			pvalue <- fit['p.value']
#			#print(paste('cutpoint=',cutpoint,', p.value=',pvalue))
#			data.subset[rowname,'pvalue'] <- pvalue
#		}
#	}
#	data.subset <- data.subset[order(data.subset$pvalue),]
#	bestcutpoint <- data.subset[1,field]  # skip the first one
#	data.subset <- createSplitField(data, field, newfield, bestcutpoint)
#	xprops(as.formula(paste('~', response, '+', newfield)), data=data.subset)
#	#counts <- xtabs(as.formula(paste('~', response, '+', newfield)), data=data.subset)
#	#print(chisq.test(counts))
#	return(bestcutpoint)
#	#return(data.subset)
#}
##findBestCutpoint(data,'nvr','MxA')

findBestCutpoint <- function(data, response, field, show.plot=F, show.table=F)
{
	data.subset <- subsetNA(data,c(response,field))
	#data.subset <- na.omit(data[,c(response,field)])
	newfield <- addFieldSuffix(field,'bestcut') #paste(field,'_bestcut',sep='')
	data.subset[,newfield] <- NA
	data.subset$index <- NA
	data.subset$cutpoint <- NA
	data.subset$pvalue <- NA
	index <- 0
	for (rowname in row.names(data.subset))
	{
		cutpoint <- data.subset[rowname,field]
		data.subset[,newfield] <- NA
		data.subset[,newfield] <- ifelse(data.subset[[field]] <= cutpoint, 0, 1)
		data.subset[rowname,'index'] <- index
		index <- index + 1
		if (length(unique(data.subset[,newfield]))>1)
		{
			counts <- xtabs(as.formula(paste('~', response, '+', newfield)), data=data.subset)
			fit <- fisher.test(counts)
			#fit <- chisq.test(counts)
			pvalue <- fit['p.value']
			#print(paste('cutpoint=',cutpoint,', p.value=',pvalue))
			data.subset[rowname,'pvalue'] <- pvalue
			data.subset[rowname,'cutpoint'] <- cutpoint
		}
	}
	data.subset <- data.subset[order(data.subset$pvalue),]
	bestcutpoint <- data.subset[1,field]
	#bestindex <- data.subset[1,'index']
	if (show.table)
	{
		data.subset <- createSplitField(data, field, newfield, bestcutpoint)
		print(data.subset)		
		print(xprops(as.formula(paste('~', response, '+', newfield)), data=data.subset, fisher=T))
	}
	if (show.plot)
	{
		#plotDistribution(data.subset$pvalue)
		plot(-log(pvalue)~cutpoint,data.subset[order(data.subset$cutpoint),], type='b', main=field)
		abline(v = bestcutpoint, col = "red")
	}	
	return(list(bestcutpoint=bestcutpoint, pvalue=data.subset[1,'pvalue']))
}
#findBestCutpoint(data,'nvr','ISG15',T)

findFieldsWithBestCutpoint <- function(data,response,fields,show.plot=F)
{
	if (show.plot)
		par(ask=T)
	pvalues <- data.frame(field=fields, row.names=fields)
	for (field in fields)
	{
		field <- stripFieldSuffix(field) #strsplit(field,"_")[[1]][1]
		results <- findBestCutpoint(data,response,field,show.plot=show.plot)
		#print(paste(field,results$bestcutpoint, results$pvalue))
		pvalues[field,'bestcutpoint'] <- results$bestcutpoint
		pvalues[field,'pvalue'] <- results$pvalue
		#print(paste(field,bestcutpoint))
	}
	print(pvalues[order(pvalues$pvalue),])
	par(ask=F)
}
#findFieldsWithBestCutpoint(data,'nvr',isgupfields)

createSplitFieldByBestCutpoint <- function(data, response, field)
{
	results <- findBestCutpoint(data,response,field)
	newfield <- paste(field,'_bestcut',sep='')
	data <- createSplitField(data,field,newfield,results$bestcutpoint)
	attr(data,paste(field,'bestcutpoint',sep='_')) <- results$bestcutpoint
	return(data)
}

createSplitByAge <- function(data)
{
	data$age10 <- floor(data$age/10)*10
	return(data)
}

chisquare<-function(counts, fisher=F)
{
	print(counts)
	print(prop.table(counts))
	print(chisq.test(counts))
	if (fisher)
		print(fisher.test(counts))
}

xprops <- function(formula, data, fisher=F)
{
	temp <- xtabs(formula, data=data)
	chisquare(temp, fisher)
	total0 <- temp[1,1]+temp[2,1]
	success0 <- temp[2,1]
	total1 <- temp[1,2]+temp[2,2]
	success1 <- temp[2,2]
	print(paste('percent for level 0:',100*(success0/total0),'%'))
	print(paste('percent for level 1:',100*(success1/total1),'%'))
}

calculateMaf <- function(aa, ab, bb)
{
	total <- aa + ab + bb
	pqtotal <- total*2
	qnum <- bb*2+ab
	pnum <- pqtotal-qnum
	q <- qnum/pqtotal
	p <- 1-q
	#q=(bb*2+ab)/(total*2)
	#p=1-q
	print(paste('aa=',aa,'ab=',ab,'bb=',bb,'total=',total))
	print(paste('aa=',aa/total,'ab=',ab/total,'bb=',bb/total))
	print(paste('p=',pnum,'q=',qnum,'pqtotal=',pqtotal))
	print(paste('p=',p,'q=',q))
	return(q)
}
#calculateMaf(630,151,14)

getMedian <- function(data.subset, field)
{
	smmry <- summary(data.subset[[field]])
	return(paste(smmry[3],' (',smmry[2],'-',smmry[5],')', sep=''))
}


wilcoxonTest <- function(data, model, verbose=FALSE)
{
	wilcox.test(titer ~ nvr, data) 
	fit <- wilcox.test(as.formula(model), data=data)
	if (verbose)
		print(fit)
	return(fit)
}

wilcoxonTests <- function(data, predictor, fields, verbose=FALSE)
{
	fits <- list()
	for (field in fields)
	{
		fit <- wilcoxonTest(data, paste(field," ~ ",predictor), verbose)
		fits[[field]] <- fit
	}
	return(fits)
}

findSignificantWilcoxonFields <- function(uniresults, cutoff=0.05, verbose=TRUE)
{
	sigfields <- c()
	for (field in names(uniresults))
	{
		fit <- uniresults[[field]]
		pvalue <- fit['p.value']
		if (pvalue<=cutoff)
		{
			sigfields <- append(sigfields,field)
			print(fit)
		}
	}
	if (length(sigfields)==0)
	{
		print(paste("no fields significant at P=",cutoff,'using all fields'))
		return(names(uniresults))
	}
	if (verbose) print(paste("significant univariate effects: ",sigfields))
	return(sigfields)
}

#kruskal.test(formula, data, subset, na.action, ...)


createSignificantMarker <- function(value)
{
	#print(paste('value=',value))
	if (is.na(value))
		return('NA')
	if (value < 0.001)
		return('***')
	else if (value < 0.01)
		return('**')
	else if (value < 0.05)
		return('*')
	else if (value < 0.1)
		return('.')
	else return('')
}

skew <- function(x)
{
	(sum((x-mean(x))^3/sqrt(var(x))^3)/length(x))
}

kurtosis <- function(x)
{
	(sum((x-mean(x))^4/var(x)^2)/length(x) - 3)
}

.ls.objects <- function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 10) {
	# based on postings by Petr Pikal and David Hinds to the r-help list in 2004
	# modified by: Dirk Eddelbuettel (http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session) 
	# I then gave it a few tweaks (show size as megabytes and use defaults that I like)
	# a data frame of the objects and their associated storage needs.
	napply <- function(names, fn) sapply(names, function(x)
					fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.size <- napply(names, object.size) / 10^6 # megabytes
	obj.dim <- t(napply(names, function(x)
						as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.dim)
	names(out) <- c("Type", "Size", "Rows", "Columns")
	out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}

containsElement <- function(x, value)
{
	for (curval in splitFields(x))
	{
		if (curval==value)
			return(TRUE)		
	}
	return(FALSE)
}
#containsElement(splitFields('a,b,c,d,e,f,g'),'d')
#containsElement(splitFields('a,b,c,d,e,f,g'),'q')

appendValues <- function(values1,values2)
{
	return(c(splitFields(values1),splitFields(values2)))
}
#appendValues('a,b,c,d,e,f,g','h,i,j,k,l')
#appendValues(c(),'q')

removeElements <- function(x, y)
{
	newvalues <- c()
	for (curval in splitFields(x))
	{
		if (!containsElement(y,curval))
			newvalues <- appendValues(newvalues,curval)
	}
	return(newvalues)
}
#removeElements('a,b,c,d,e,f,g','b,f')

######################################################3

getCounts <- function(data, field, values, separator='|')
{
	counts <- c()
	for (value in values)
	{
		counts <- c(counts, length(which(data[[field]]==value)))
	}
	return(counts)
}

addCountRow <- function(table, data, field, label, values, subsets=c(), separator='|')
{
	#if (label=='') label=field
	table[field,'variable'] <- label
	table[field,'all'] <- paste(getCounts(data,field,values),collapse=separator)
	if (length(subsets)>0)
	{
		for (num in 1:length(subsets))
		{
			table[field,paste('subset',num,sep='')] <- paste(getCounts(subsets[num][[1]],field,values),collapse=separator)
		}
	}
	return(table)	
}

getMedianAndRange <- function(data.subset, field)
{
	smmry <- summary(data.subset[[field]])
	return(paste(smmry[3],' (',format(smmry[1], digits=2),'-',format(smmry[6], digits=2),')', sep=''))
}

addMedianRow <- function(table, data, field, label, subsets=c())
{
	#if (label=='') label=field
	table[field,'variable'] <- label
	table[field,'all'] <- getMedianAndRange(data,field)
	if (length(subsets)>0)
	{
		for (num in 1:length(subsets))
		{
			table[field,paste('subset',num,sep='')] <- getMedianAndRange(subsets[num][[1]],field)
		}
	}
	return(table)	
}


pad <- function(num,ch='0')
{
	if (num < 10)
		num <- paste('0', num, sep='')
	return(as.character(num))
}

# http://tolstoy.newcastle.edu.au/R/help/05/12/17151.html
list2ascii <- function(x,file=paste(deparse(substitute(x)),".txt",sep="")) {
	
	# MHP July 7, 2004
	# R or S function to write an R list to an ASCII file.
	# This can be used to create files for those who want to use
	# a spreadsheet or other program on the data.
	#
	tmp.wid = getOption("width")  # save current width
	options(width=10000)          # increase output width
	sink(file)                    # redirect output to file
	print(x)                      # print the object
	sink()                        # cancel redirection
	options(width=tmp.wid)        # restore linewidth
	return(invisible(NULL))       # return (nothing) from function
	
}

writeTable <- function(table, filename, verbose=TRUE, row.names=TRUE)
{
	#col.names <- c('id',names(table))
	write.table(table, filename, quote=FALSE, row.names = row.names, sep='\t', na = '')
	if (verbose)
		print(paste('wrote table to file:',filename))
}


