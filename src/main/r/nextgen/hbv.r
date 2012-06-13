
analyzeG2A <- function(config, group='HBV_nucleoside_analogues', exclude.ntnums=NULL, include.ntnums=NULL, minvariant=100)
{
	data <- data.frame()
	subjects <- getSubjectsForGroup(config,group)	
	for (subject in subjects)
	{
		ids <- config@data[which(config@data$subject==subject),]$id
		if (length(ids)==1)
			next()
		#throw('only one sample for subject ',subject)
		id1 <- ids[1]
		id2 <- ids[length(ids)]
		
		row <- config@data[id1,]
		filename1 <- concat(config@counts.dir,'/',id1,'.nt.txt')
		filename2 <- concat(config@counts.dir,'/',id2,'.nt.txt')
		
		data.subset1 <- loadDataFrame(filename1)
		data.subset2 <- loadDataFrame(filename2)
		rownames(data.subset1) <- data.subset1$position
		rownames(data.subset2) <- data.subset2$position
		
		for (ntnum in rownames(data.subset1))
		{
			if (!rowExists(data.subset2, ntnum))
				next()
			row1 <- data.subset1[ntnum,]
			row2 <- data.subset2[ntnum,]
			data <- rbind(data, data.frame(position=ntnum, subject=row$subject, refnt=row1$refnt,
							a1=row1$a, c1=row1$c, g1=row1$g, t1=row1$t, depth1=row1$depth,
							a2=row2$a, c2=row2$c, g2=row2$g, t2=row2$t, depth2=row2$depth))
		}
	}
	
	data$ag1 <- data$a1+data$g1
	data$ag2 <- data$a2+data$g2
	data <- subset(data, ag1>minvariant & ag2>minvariant)
	data$ratio1 <- data$a1/(data$a1+data$g1)
	data$ratio2 <- data$a2/(data$a2+data$g2)
	data$diff <- data$ratio2-data$ratio1
	data <- data[!is.nan(data$diff),]
	
	filename <- concat(config@counts.dir,'/ntcounts.txt')
	writeTable(data,filename,row.names=FALSE)
	
	#counts$dir <- ifelse(counts$diff<0,'down',ifelse(counts$diff>0,'up',NA))
	
	
	return(data)
}
counts <- analyzeG2A(config)

