getCodonCountSubset <- function(config, group, region, filetype, start, end=start, cutoff=0)
{
	filename <- concat(config@counts.dir,'/',group,'.',filetype,'.txt')
	data <- loadDataFrame(filename)
	data.subset <- data[which(data$region==region & data$aanum>=start & data$aanum<=end & data$count>=cutoff),]
	#data.subset$replicate <- factor(data.subset$replicate)
	data.subset$label <- factor(data.subset$label)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,'G9','NS3aa36','codons',36)

makeVariantTable <- function(config, type, group, region, cutoff=0)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	gene <- config@regions[region,'gene']
	aanum <- as.integer(config@regions[region,'focus'])
	ref <- getRefForGroup(config,group)
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	data.subset <- getCodonCountSubset(config,group,region,type,aanum,cutoff=cutoff)
	if (type=='codons')
		frmla <- as.formula(codon ~ label)
	else frmla <- as.formula(aa ~ label)
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
#makeVariantTable(config,'aa','G9','NS3aa36')
#makeVariantTable(config,'aa','KT9','NS3aa156')
#makeVariantTable(config,'aa','KT9','NS3aa156')
#makeVariantTable(config,'codons','KT9','NS3aa156')
#makeVariantTable(config,'codons','8538159','NS3-156-R@NS3aa156')


makeVariantTables <- function(config, type, group=NULL, ...)
{
	if (is.null(group))
		groups <- config@groups
	else groups <- splitFields(group)
	tables <- list()
	for (group in groups)
	{
		tables[[group]] <- list()
		for (region in getRegionsForGroup(config,group))
		{
			try({
			print(concat('group=',group,', region=',region))
			tbl <- makeVariantTable(config, type, group, region, ...)
			tables[[group]][[region]] <- tbl
			}, silent=FALSE)
		}
	}
	return(tables)
}
#tables <- makeVariantTables(config, 'codons', 'PXB0220-0030')

getSubjectsByGoal <- function(config, goal)
{
	return(unique(config@runs[which(config@runs$goal==goal),'subject']))
}
#getSubjectsByGoal(config,'goal1')


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

writeCodonTable <- function(config, codon.tables, group, region)
{
	tbl <- codon.tables[[group]][[region]]
	ref <- getRefForGroup(config,group)
	identifier <- concat(group,'-',region,'-',ref)
	filename <- concat(config@out.dir,'/tables/table-',identifier,'.txt')
	try(colnames(tbl)[1] <-identifier, silent=FALSE)
	writeTable(tbl,filename,row.names=FALSE)
}

writeCodonTables <- function(config, groups=config@groups, cutoff=2)
{ 
	codon.tables <- makeVariantTables(config,'codons',groups,cutoff=cutoff)
	for (group in groups)
	{
		for (region in getRegionsForGroup(config,group))
		{
			writeCodonTable(config, codon.tables, group, region)
		}
	}
	return(codon.tables)
}
#codon.tables <- writeCodonTables(config)
#writeCodonTables(config,'confirm_with_new_reagents')
#writeCodonTables(config,'PXB0220-0030')

