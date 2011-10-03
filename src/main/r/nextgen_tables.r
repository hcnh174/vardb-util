
getCodonCountSubset <- function(config, subject, region, filetype, start, end=start, cutoff=0)
{
	filename <- concat(config@counts.dir,'/',subject,'.',filetype,'.txt')
	data <- loadDataFrame(filename)
	data.subset <- data[which(data$region==region & data$aanum>=start & data$aanum<=end & data$count>=cutoff),]
	data.subset$replicate <- factor(data.subset$replicate)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,'KT9','NS3aa156','codons',156)

makeVariantTable <- function(config, type, subject, region, cutoff=0)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	aanum <- as.integer(config@regions[region,'focusaa'])
	positions <- getCodonPositionsForRegion(config,region)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	data.subset <- getCodonCountSubset(config,subject,region,type,aanum,cutoff=cutoff)
	if (type=='codons')
		frmla <- as.formula(codon ~ replicate)
	else frmla <- as.formula(aa ~ replicate)
	counts <- cast(data.subset, frmla, value='count', fun.aggregate=function(x) return(x[1])); counts
	counts <- counts[order(counts[,2], decreasing=TRUE),]
	row.names(counts) <- seq(nrow(counts))
	# add an asterisk to indicate the reference codon
	if (type=='codons')
	{
		if (length(which(counts$codon==refcodon))==0)
		{
			row <- counts[1,]
			row[,'codon'] <- refcodon
			row[,-1] <- NA
			counts <- rbind(counts,row)
		}
		counts$aa <- sapply(counts$codon,translateCodon)
		counts[which(counts$codon==refcodon),'codon'] <- concat(refcodon,'*')
		cols <- c(1,length(colnames(counts)),2:(length(colnames(counts))-1))
		counts <- counts[,cols]
	}
	if (type=='aa')
	{
		refaa <- translateCodon(refcodon)
		if (length(which(counts$aa==refaa))==0)
		{
			row <- counts[1,]
			row[,'aa'] <- refaa
			row[,-1] <- NA
			counts <- rbind(counts,row)
		}
		counts[which(counts$aa==refaa),'aa'] <- concat(refaa,'*')		
	}
	return(counts)
}
#makeVariantTable(config,'aa','KT9','NS3aa156')
#makeVariantTable(config,'codons','KT9','NS3aa156')
#makeVariantTable(config,'codons','8538159','NS3-156-R@NS3aa156')

makeVariantTables <- function(config, type, subject=NULL, ...)
{
	if (is.null(subject))
		subjects <- config@subjects
	else subjects <- splitFields(subject)
	tables <- list()
	for (subject in subjects)
	{
		tables[[subject]] <- list()
		for (region in getRegionsForSubject(config,subject))
		{
			#try({
			#print(concat('subject=',subject,', region=',region))
			tbl <- makeVariantTable(config, type, subject, region, ...)
			tables[[subject]][[region]] <- tbl
			#}, silent=FALSE)
		}
	}
	return(tables)
}
#tables <- makeVariantTables(config, 'codons')
#tables <- makeVariantTables(config, 'codons', '8538159')
#tables <- makeVariantTables(config, 'codons', '10348001')

getSubjectsByGoal <- function(config, goal)
{
	return(unique(config@runs[which(config@runs$goal==goal),'subject']))
}
#getSubjectsByGoal(config,'goal1')
#
#makeVariantTablesByGoal <- function(config, type, goal=NULL, ...)
#{
#	if (is.null(goals))
#		goals <- rownames(config@goals)
#	else goals <- splitFields(goal)
#	tables <- list()
#	for (goal in goals)
#	{
#		for (subject in getSubjectsByGoal(config,goal))
#		{
#			tables[[subject]] <- list()
#			for (region in getRegionsForSubject(config,subject))
#			{
#				#try({
#				#print(concat('subject=',subject,', region=',region))
#				tbl <- makeVariantTable(config, type, subject, region, ...)
#				tables[[subject]][[region]] <- tbl
#				#}, silent=FALSE)
#			}
#		}
#	}
#	return(tables)
#}
#tables <- makeVariantTables(config, 'codons')
#tables <- makeVariantTables(config, 'codons', '8538159')
#tables <- makeVariantTables(config, 'codons', '10348001')

appendVariantTablesToWord <- function(tables)
{
	for (subject in names(tables))
	{
		wdHeading(level=2,concat('Subject: ',subject))
		for (region in names(tables[[subject]]))
		{
			wdHeading(level=2,concat('Region: ',region))
			tbl <- tables[[subject]][[region]]
			tbl <- replaceNAs(tbl)
			tbl <- format(tbl)			
			wdTable(tbl)
		}
		wdPageBreak()
	}
}


outputVariantTablesToWord <- function(subjects=NULL, filename='tables.doc',...)
{
	filename <- concat(getwd(),'/',config@out.dir,filename)
	
	#codon.tables <- makeCodonTables(config,subjects,...)
	#aa.tables <- makeAminoAcidTables(config,subjects,...)
	
	codon.tables <- makeVariantTables(config,'codons',subjects,...)
	aa.tables <- makeVariantTables(config,'aa',subjects,...)
	
	wdGet(visible=FALSE)
	wdNewDoc(filename)
	wdSection('Codon tables', newpage=FALSE)
	appendVariantTablesToWord(codon.tables)
	wdSection('Amino acid tables', newpage=TRUE)
	appendVariantTablesToWord(aa.tables)
	wdSave(filename)
	wdQuit()
}
#outputVariantTablesToWord()

getLabelForReplicate <- function(config, subject, replicate)
{
	value <- unique(config@runs[which(config@runs$subject==subject & config@runs$replicate==replicate),'label'])
	if (length(value)>1)
		throw('more than one label found for replicate: subject=',subject,', replicate=',replicate,': ',joinFields(value))
	return(value)
}
#getLabelForReplicate(config,'KT9',1)

applyReplicateLabels <- function(config, tbl, subject)
{
	cols <- c()
	for (col in colnames(tbl))
	{
		if (col %in% c('codon','aa'))
			cols <- c(cols,col)
		else
		{
			col <- getLabelForReplicate(config,subject,col)
			cols <- c(cols,col)
		}
	}
	colnames(tbl) <- cols
	return(tbl)
}
#applyReplicateLabels(config, codon.tables[['PXB0220-0002']][['NS5Aaa31']],'PXB0220-0002')

appendVariantTablesToLatex <- function(config,tables)
{
	for (subject in names(tables))
	{
		#print(concat('Mutation: ',as.character(config@subjects[subject,'mutation'][[1]])))
		for (region in names(tables[[subject]]))
		{
			tbl <- tables[[subject]][[region]]
			caption <- concat('Subject: ',as.character(subject[1]),', Region: ',as.character(region[1]))
			caption <- concat(caption,'\\newline  Description: ',as.character(config@subjects[subject,'description']))
			caption <- concat(caption,'\\newline')
			xtbl <- xtable(tbl, caption=caption, digits=0)
			print(xtbl, include.rownames=FALSE, caption.placement='top', latex.environments='flushleft')
		}
	}
}

###############################################

write_table <- function(config, codon.tables, goal, subject, region)
{
	tbl <- codon.tables[[subject]][[region]]
	tbl <- applyReplicateLabels(config,tbl,subject)
	ref <- get_ref_for_subject(config,subject,region)
	filename <- concat(config@out.dir,'/tables/table-',goal,'-',subject,'-',region,'-',ref,'.txt')
	writeTable(tbl,filename,row.names=FALSE)
}

write_tables <- function(config)
{ 
	codon.tables <- makeVariantTables(config,'codons',cutoff=2)
	for (goal in rownames(config@goals))
	{
		subjects <- getSubjectsByGoal(config,goal)
		for (subject in subjects)
		{
			for (region in getRegionsForSubject(config,subject))
			{
				write_table(config, codon.tables, goal, subject, region)
			}
		}
	}
}
#write_tables(config)
