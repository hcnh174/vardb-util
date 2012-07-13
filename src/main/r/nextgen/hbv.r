#analyzeG2A <- function(config, group='HBV_nucleoside_analogues', exclude.ntnums=NULL, include.ntnums=NULL, minvariant=100)
#{
#	data <- data.frame()
#	subjects <- getSubjectsForGroup(config,group)
#	for (subject in subjects)
#	{
#		ids <- config@data[which(config@data$subject==subject),]$id
#		if (length(ids)==1)
#			next()
#		#throw('only one sample for subject ',subject)
#		id1 <- ids[1]
#		id2 <- ids[length(ids)]
#		
#		row <- config@data[id1,]
#		filename1 <- concat(config@counts.dir,'/',id1,'.nt.txt')
#		filename2 <- concat(config@counts.dir,'/',id2,'.nt.txt')
#		
#		data.subset1 <- loadDataFrame(filename1)
#		data.subset2 <- loadDataFrame(filename2)
#		rownames(data.subset1) <- data.subset1$position
#		rownames(data.subset2) <- data.subset2$position
#		
#		for (ntnum in rownames(data.subset1))
#		{
#			if (!rowExists(data.subset2, ntnum))
#				next()
#			row1 <- data.subset1[ntnum,]
#			row2 <- data.subset2[ntnum,]
#			data <- rbind(data, data.frame(position=ntnum, subject=row$subject, refnt=row1$refnt,
#							a1=row1$a, c1=row1$c, g1=row1$g, t1=row1$t, depth1=row1$depth,
#							a2=row2$a, c2=row2$c, g2=row2$g, t2=row2$t, depth2=row2$depth))
#		}
#	}
#	
#	data$ag1 <- data$a1+data$g1
#	data$ag2 <- data$a2+data$g2
#	data <- subset(data, ag1>minvariant & ag2>minvariant)
#	data$ratio1 <- data$a1/(data$a1+data$g1)
#	data$ratio2 <- data$a2/(data$a2+data$g2)
#	data$diff <- data$ratio2-data$ratio1
#	data <- data[!is.nan(data$diff),]
#	
#	filename <- concat(config@counts.dir,'/ntcounts.txt')
#	writeTable(data,filename,row.names=FALSE)
#	
#	#counts$dir <- ifelse(counts$diff<0,'down',ifelse(counts$diff>0,'up',NA))
#	
#	
#	return(data)
#}
##counts <- analyzeG2A(config)

analyzeG2A <- function(config, group='HBV_nucleoside_analogues', exclude.ntnums=NULL, include.ntnums=NULL, minvariant=100)
{
	library(reshape)
	
	data <- data.frame(position=integer(), refnt=character())
	subjects <- getSubjectsForGroup(config,group)
	# find the max number of time points
	maxlength <- 0
	for (subject in subjects)
	{
		ids <- config@data[which(config@data$subject==subject),]$id
		if (length(ids)>maxlength)
			maxlength <- length(ids)
	}
	print(maxlength)
	
	for (num in 1:maxlength)
	{
		data[[concat('a',num)]] <- integer()
		data[[concat('c',num)]] <- integer()
		data[[concat('g',num)]] <- integer()
		data[[concat('t',num)]] <- integer()
		data[[concat('depth',num)]] <- integer()
	}
	
	
	for (subject in subjects)
	{
		ids <- config@data[which(config@data$subject==subject),]$id
		if (length(ids)==1)
			next()
		#throw('only one sample for subject ',subject)
		
		id <- ids[1]
		filename <- concat(config@counts.dir,'/',id,'.nt.txt')
		data.subset <- loadDataFrame(filename)
		rownames(data.subset) <- data.subset$position
		colnames(data.subset) <- splitFields('position,refnt,a1,c1,g1,t1,depth1')
		#data.subset$ratio1 <- data.subset$a1/(data.subset$a1+data.subset$g1)
		
		for (num in 2:maxlength)
		{
			data.subset[[concat('a',num)]] <- NA
			data.subset[[concat('c',num)]] <- NA
			data.subset[[concat('g',num)]] <- NA
			data.subset[[concat('t',num)]] <- NA
			data.subset[[concat('depth',num)]] <- NA
		}

		for (num in 2:length(ids))
		{
			id <- ids[num]
			filename <- concat(config@counts.dir,'/',id,'.nt.txt')
			data.subset2 <- loadDataFrame(filename)
			rownames(data.subset2) <- data.subset2$position
			
			for (ntnum in rownames(data.subset))
			{
				if (!rowExists(data.subset2, ntnum))
					next()
				row <- data.subset2[ntnum,]
				data.subset[ntnum,concat('refnt',num)] <- row$refnt
				data.subset[ntnum,concat('a',num)] <- row$a
				data.subset[ntnum,concat('c',num)] <- row$c
				data.subset[ntnum,concat('g',num)] <- row$g
				data.subset[ntnum,concat('t',num)] <- row$t
				data.subset[ntnum,concat('depth',num)] <- row$depth
			}
			#data.subset[[concat('ratio',num)]] <- data.subset[[concat('a',num)]]/(data.subset[[concat('a',num)]] + data.subset[[concat('g',num)]])
		}

		for (num in 1:length(ids))
		{
			data.subset[[concat('ratio',num)]] <- data.subset[[concat('a',num)]]/(data.subset[[concat('a',num)]] + data.subset[[concat('g',num)]])
		}
		
		filename <- concat(config@counts.dir,'/ntcounts-',subject,'.txt')
		writeTable(data.subset,filename,row.names=FALSE)
		#print(head(data.subset))
		#data <- rbind(data, data.subset)
		
	}
	#data <- subset(data, ag1>minvariant & ag2>minvariant)
	#data <- data[!is.nan(data$diff),]
	
	#filename <- concat(config@counts.dir,'/ntcounts.txt')
	#writeTable(data,filename,row.names=FALSE)
	
	#counts$dir <- ifelse(counts$diff<0,'down',ifelse(counts$diff>0,'up',NA)
	
	return(data)
}
#counts <- analyzeG2A(config)

#
#analyzeG2AforSubject <- function(config, subject, minvariant=100)
#{	
#	ids <- config@data[which(config@data$subject==subject),]$id
#	if (length(ids)==1)
#		next()
#	#throw('only one sample for subject ',subject)
#	
#	id <- ids[1]
#	filename <- concat(config@counts.dir,'/',id,'.nt.txt')
#	data <- loadDataFrame(filename)
#	rownames(data) <- data$position
#	colnames(data) <- splitFields('position,refnt,a1,c1,g1,t1,depth1')
#	#data$ratio1 <- data$a1/(data$a1+data$g1)
#	
#	
#	for (num in 2:length(ids))
#	{
#		id <- ids[num]
#		filename <- concat(config@counts.dir,'/',id,'.nt.txt')
#		data2 <- loadDataFrame(filename)
#		rownames(data2) <- data2$position
#		
#		for (ntnum in rownames(data))
#		{
#			if (!rowExists(data2, ntnum))
#				next()
#			row <- data2[ntnum,]
#			data[ntnum,concat('refnt',num)] <- row$refnt
#			data[ntnum,concat('a',num)] <- row$a
#			data[ntnum,concat('c',num)] <- row$c
#			data[ntnum,concat('g',num)] <- row$g
#			data[ntnum,concat('t',num)] <- row$t
#			data[ntnum,concat('depth',num)] <- row$depth
#		}
#		#data[[concat('ratio',num)]] <- data[[concat('a',num)]]/(data[[concat('a',num)]] + data[[concat('g',num)]])
#	}
#	
#	for (num in 1:length(ids))
#	{
#		data[[concat('ratio',num)]] <- data[[concat('a',num)]]/(data[[concat('a',num)]] + data[[concat('g',num)]])
#	}
#	
#	data <- subset(data, (a1+g1)>minvariant & depth1>minvariant)
#	
#	filename <- concat(config@counts.dir,'/ntcounts-',subject,'.txt')
#	writeTable(data,filename,row.names=FALSE)
#	return(data)
#}
##counts <- analyzeG2AforSubject(config,subject)
#
##boxplot(counts$ratio1, counts$ratio2, counts$ratio3, counts$ratio4, counts$ratio5, counts$ratio6)


analyzeG2AforSubject <- function(config, subject, minvariant=100)
{	
	ids <- config@data[which(config@data$subject==subject),]$id
	if (length(ids)==1)
		return()
	#throw('only one sample for subject ',subject)
	
	id <- ids[1]
	filename <- concat(config@counts.dir,'/',id,'.nt.txt')
	data <- loadDataFrame(filename)
	rownames(data) <- data$position
	colnames(data) <- splitFields('position,refnt,a1,c1,g1,t1,depth1')
	#data$ratio1 <- data$a1/(data$a1+data$g1)
	
	
	for (num in 2:length(ids))
	{
		id <- ids[num]
		filename <- concat(config@counts.dir,'/',id,'.nt.txt')
		data2 <- loadDataFrame(filename)
		rownames(data2) <- data2$position
		
		for (ntnum in rownames(data))
		{
			if (!rowExists(data2, ntnum))
				next()
			row <- data2[ntnum,]
			data[ntnum,concat('refnt',num)] <- row$refnt
			data[ntnum,concat('a',num)] <- row$a
			data[ntnum,concat('c',num)] <- row$c
			data[ntnum,concat('g',num)] <- row$g
			data[ntnum,concat('t',num)] <- row$t
			data[ntnum,concat('depth',num)] <- row$depth
		}
		#data[[concat('ratio',num)]] <- data[[concat('a',num)]]/(data[[concat('a',num)]] + data[[concat('g',num)]])
	}
	
	for (num in 1:length(ids))
	{
		data[[concat('ratio',num)]] <- data[[concat('a',num)]]/(data[[concat('a',num)]] + data[[concat('g',num)]])
	}
	
	data <- subset(data, (a1+g1)>minvariant & depth1>minvariant)
	
	#filename <- concat(config@counts.dir,'/ntcounts-',subject,'.txt')
	#writeTable(data,filename,row.names=FALSE)
	return(data)
}
#counts <- analyzeG2AforSubject(config,subject)


