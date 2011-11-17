library(lattice)
library(R.oo)
#library(car)
#library(MASS)
#library(xtable)
#library(R.utils)

#library(rcom)
#library(agce)
#library(leaps)
#library(corrgram)
#library(vcd)
#library(reshape)
#library(doBy)
#library(gplot)
#library(HH)
#library(rrcov)

lattice.options(default.args = list(as.table = TRUE))
options(contrasts=c("contr.sum","contr.poly"))
options("width"=200)
#options(warnPartialMatchArgs=TRUE)
#options(warnPartialMatchAttr=TRUE)
#options(warnPartialMatchDollar=TRUE)

#################################################################
		
#concat <- function(...)
#{
#	return(paste(..., sep=''))
#}

#http://jermdemo.blogspot.com/2011/10/making-rs-paste-act-more-like-concat.html
concat <- function(..., sep='', collapse=NULL)
{
	strings<-list(...)
	#catch NULLs, NAs
	if (all(unlist(plyr::llply(strings,length))>0) &&	all(!is.na(unlist(strings))))
	{
		do.call("paste", c(strings, list(sep = sep, collapse = collapse)))
	}
	else
	{
		NULL
	}
}
#concat('abc','def','ghi')

printcat <- function(...)
{
	print(concat(...))
}
#printcat('abc','def','ghi')

loadUtilFiles <- function(filenames)
{
	libdir <- Sys.getenv("VARDB_RUTIL_HOME")
	filenames <- splitFields(filenames)
	filenames <- c(filenames,'external')
	for (filename in filenames)
	{
		#source(paste(libdir,filename,'.r',sep=''))
		source(concat(libdir,'/',filename,'.r'))
	}
}

setCurDir <- function(dir)
{
	setwd(paste(Sys.getenv("ANALYSIS_HOME"),'/',dir,'/',sep=''))
}

loadLibrary <- function(name)
{
	require(name, quietly=TRUE, warn.conflicts=FALSE, character.only=TRUE)	
}

loadDataFrame <- function(filename, idcol=NULL, stringsAsFactors=FALSE)#default.stringsAsFactors())
{
	dataframe <- read.table(filename, header=TRUE, encoding='UTF-8', sep = '\t', comment.char='#', stringsAsFactors=stringsAsFactors)
	if (!is.null(idcol))
		rownames(dataframe) <- dataframe[[idcol]]
	return(dataframe)
}

#temporarily alias the old version
#loadDataframe <- loadDataFrame

appendLog <- function(str, logfile='commands.log')
{
	str <- concat(str,'\n')
	cat(str,file=logfile,append=TRUE)
}

runCommand <- function(...)
{
	command <- concat(...)
	#tstmp <- concat('[',format(Sys.time(),'%b %d %H:%M:%S'),']')
	tstmp <- concat('[',Sys.time(),']')
	#timestamp()
	print(tstmp)
	print(command)
	if (command!='')
	{
		system(command)
		#appendLog(command)
	}
}
#runCommand('ls')

testCommand <- function(...)
{
	print(concat(...))
}

getField <- function(data, id, col, msg=concat('cannot find col [',col,'] for row [',id,'] in dataframe'))
{
	value <- data[id,col]
	if (is.null(value) || is.na(value))
	{
		if (nrow(data[id,])==0)
			throw('no rows found with id ',id)
		if (is.null(data[[col]]))
			throw('column not found: ',col)
		throw(msg)
	}
	return(value)
}
#getField(config@regions,'NS3aa156','feature')

# check normality
plotDistribution <- function(values, ylab='values')
{	
	values <- values[!is.na(values)]
	print(sort(values))
	log.values <- log(values)
	print(outlierTest(lm(values ~ 1)))
	oldpar <- par(mfrow=c(1,3)); on.exit(par(oldpar))	
	qqnorm(values, ylab=ylab)
	try(qqline(values),silent=T)
	hist(values)
	qqnorm(log.values, ylab=ylab)
	try(qqline(log.values),silent=T)
	#par(mfrow=c(1,1))
	if (min(values>0))
		print(powerTransform(values))
}

plotDistributions <- function(dataframe, fields)
{
	oldpar <- par(ask=T); on.exit(par(oldpar))
	for (field in splitFields(fields))
	{
		print('******************************')
		print(field)
		plotDistribution(dataframe[,field], field)
	}
	#par(oldpar)
}

# automate regression tests

splitFields <- function(str, delimiter=',')
{
	if (length(str)==0)
		return(c())
	if (length(str)>1)
		return(str)
	strsplit(str,delimiter)[[1]]
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
		{oldpar <- par(ask=T); on.exit(par(oldpar))}
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
	#par(ask=F)
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

appendUniqueValues <- function(values1,values2)
{
	return(unique(appendValues(values1,values2)))
}

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

intersectValues <- function(values1, values2)
{
	values <- c()
	for (value in splitFields(values1))
	{
		if (containsElement(values2,value))
			values <- appendValues(values,value)
	}
	return(values)
}
#intersectValues('A,B,C','B,D')

excludeColumns <- function(data, cols)
{
	cols <- splitFields(cols)
	keep <- removeElements(colnames(data),cols)
	return(data[,keep])
}
#data <- excludeColumns(data,'constant1,constant2')

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

writeTable <- function(table, filename, verbose=TRUE, row.names=TRUE, col.names=TRUE, eol='\n')
{
	#col.names <- c('id',names(table))
	write.table(table, filename, quote=FALSE, row.names=row.names, col.names=col.names, sep='\t', na = '', eol=eol)
	if (verbose)
		print(paste('wrote table to file:',filename))
}

createFile <- function(file)
{
	cat('', file=file)
}

appendFile <- function(file, ...)
{
	print(concat(...))
	cat(...,'\n', sep='', file=file, append=TRUE)
}

getFileSize <- function(filename)
{
	return(file.info(filename)$size/1000)
}
#getFileSize(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))

stripPath <- function(filenames)
{
	return(sub("(.*)[/\\](.*)$", "\\2", filenames))
}
#stripPath('C:/workspace/vardb-util/src/main/r/tables.rnw')
#stripPath('C:\\workspace\\vardb-util\\src\\main\\r\\tables.rnw')

charAt <- function(str, index)
{
	return(substr(str, index, index))
}
#charAt('abcdefghijk',5)

#returns the current subdirectory; use getwd() for the full path
getCurDir <- function(dir=getwd())
{
	if (charAt(dir, nchar(dir))=='/')
		dir <- substr(dir,1,nchar(dir)-1)
	subdirs <- splitFields(dir,'/')
	return(subdirs[length(subdirs)])
}
#getCurDir('C:/workspace/vardb-util/src/main/r/')

getFileExtension <- function (filenames) 
{
	filenames <- as.character(filenames)
	n <- length(filenames)
	if (length(grep("\\.", filenames)) < n) 
		return(filenames)
	ext <- sub("(.*)\\.(.*)$", "\\2", filenames)
	return(ext)
}
#getFileExtension('L_CTTGTA_L007_R1_001.fastq.gz')

openPdf <- function(filename)
{
	system(concat('evince ',filename)) #system(concat('open ',pdffile))
}

sweaveToPdf <- function(filename, clean=TRUE)
{
	stem <- stripPath(filename)
	stem <- stripExtension(stem)
	auxfile <- concat(stem,'.aux')
	logfile <- concat(stem,'.log')
	texfile <- concat(stem,'.tex')
	pdffile <- concat(stem,'.pdf')
	
	#if (isOpen(pdffile))
	#	throw(concat('Pdf file is already open: ',pdffile))
	
	require(R.utils)
	
	if (isFile(auxfile))
		system(concat('rm ',logfile))
	if (isFile(auxfile))
		system(concat('rm ',auxfile))
	if (isFile(texfile))
		system(concat('rm ',texfile))
	if (isFile(pdffile))
		system(concat('rm ',pdffile))
	
	#Sweave("tables.Rnw", driver=cacheSweaveDriver());
	Sweave(filename)
	if (!isFile(texfile))
		throw(concat('Tex file not created for Sweave file: ',filename))
	system(concat('pdflatex ',texfile)) #-quiet
	if (!isFile(pdffile))
		throw(concat('Pdf file not created for tex file: ',texfile))
	if (clean)
	{
		system(concat('rm ',auxfile))
		system(concat('rm ',logfile))
		system(concat('rm ',texfile))
		system(concat('rm --f Rplots.pdf'))
		system(concat('rm --f ',stem,'-*.pdf'))
	}
	openPdf(pdffile)
}
#sweaveToPdf('tables.Rnw')


parseRanges <- function(str)
{
	items <- strsplit(str,',')[[1]]
	nums <- c()
	for (item in items)
	{
		pair <- strsplit(item,'-')[[1]]
		if (length(pair)==1)
			nums <- c(nums,as.integer(pair))
		else nums <- c(nums,seq(as.integer(pair[1]), as.integer(pair[2]), 1))
	}
	return(nums)
}
#parseRanges('3495-3736,3861-4097,6321-6506,6510-6792')

#########################################################################3

calcMfrowLayout <- function(num, maxperrow=2)
{
	numrows <- floor(num/maxperrow)
	if (num %% maxperrow >0)
		numrows <- numrows+1
	mfrow <- c(numrows,maxperrow)
	return(mfrow)
}
#calcMfrowLayout(10,3)


checkFileExists <- function(filename)
{
	if (!file.exists(filename))
		throw(concat('file does not exist: ',filename))
}
#checkFileExists('test.txt')

checkFilesExist <- function(...)
{
	filenames <- c(...)
	for (filename in filenames)
	{
		checkFileExists(filename)
	}
}
#checkFilesExist('abc.txt','def.txt','ghi.txt')

startsWith <- function(str, target)
{
	#str <- 'fsdfsf'; target <- 'HBV'
	str <- substr(str,1,nchar(target))
	return(str==target)
	#return (charmatch(target, str)==1)
}
#startsWith('HBV.TransHigh_6','HBV.TransHigh')
#startsWith('global normalization','HBV.TransHigh')

joinFields <- function(fields, delimiter=',')
{
	return(paste(fields,collapse=delimiter))
}
#joinFields(c('A','B','C'))

printParams <- function(...)
{
	args <- as.list(substitute(list(...)))[-1L]
	values <- list(...)
	arr <- c()
	for (i in 1:length(args))
	{
		name <- args[[i]]
		value <- values[[i]] 
		arr <- c(arr,concat(name,'=',value))
	}
	print(joinFields(arr,', '))
}
#printParams(region,ref)

setLatticePropertyField <- function(property,field,value)
{
	oldprop <- trellis.par.get(property)
	prop <- oldprop
	prop[[field]] <- value
	trellis.par.set(property,prop)
	return(oldprop)
}
#setLatticePropertyField('strip.background','col','red')

makeBackupFile <- function(filename)
{
	if (file.exists(filename))
	{
		bakfilename <- concat(filename,'.bak')
		file.rename(filename,bakfilename)
		print(bakfilename)
	}		
}

# copied directly from limma removeExt function
stripExtension <- function (filenames) 
{
	filenames <- as.character(filenames)
	n <- length(filenames)
	if (length(grep("\\.", filenames)) < n) 
		return(filenames)
	ext <- sub("(.*)\\.(.*)$", "\\2", filenames)
	if (all(ext[1] == ext)) 
		return(sub("(.*)\\.(.*)$", "\\1", filenames))
	else return(filenames)
}
#stripExtension('tables.Rnw')

changeExtension <- function(filename, newext)
{
	stem <- stripExtension(filename)
	return(concat(stem,'.',newext))
}
#changeExtension('test.fastq','bam')

######################################################

makeRow <- function(data, vals=NULL)
{
	row <- data[1,]
	row[,-1] <- NA
	for (col in names(vals))
	{
		row[,col] <- vals[[col]]
	}	
	return(row)
}
#makeRow(counts, list(codon='ABC'))

addRow <- function(data, vals=NULL)
{
	row <- makeRow(data,vals)
	data <- rbind(data,row)
	return(data)
}
#addRow(counts, list(codon='ABC'))

makeSubDirs <- function(dir, subdirs)
{
	runCommand('mkdir ',dir,' -p')
	for (subdir in splitFields(subdirs))
	{
		subdir <- concat(dir,'/',subdir)
		runCommand('mkdir ',subdir,' -p')
	}
}

