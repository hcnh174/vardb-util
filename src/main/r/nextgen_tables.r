getCodonCountSubset <- function(config, group, subgroup, region, filetype, start, end=start, cutoff=0)
{
	samples <- getSamplesForSubGroup(config,group,subgroup)
	data <- NULL
	for (sample in samples)
	{
		filename <- getCodonCountFilename(config,sample,filetype)
		data.sample <- loadDataFrame(filename)
		#print(head(data.sample))
		if (is.null(data))
			data <- data.sample
		else data <- rbind(data,data.sample)
	}
	print(head(data))
	data.subset <- data[which(data$region==region & data$aanum>=start & data$aanum<=end & data$count>=cutoff),]
	print(head(data.subset))
	data.subset$column <- factor(data.subset$column)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy','NS3aa36','codons',)

makeVariantTable <- function(config, type, group, subgroup, region, cutoff=0)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	gene <- config@regions[region,'gene']
	aanum <- as.integer(config@regions[region,'focus'])
	ref <- getRefForGroup(config,group)
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	data.subset <- getCodonCountSubset(config,group,subgroup,region,type,aanum,cutoff=cutoff)
	if (type=='codons')
		frmla <- as.formula(codon ~ column)
	else frmla <- as.formula(aa ~ column)
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
#tbl <- makeVariantTable(config, 'codons','MP-424', 'MP-424', 'NS3aa156')
#makeVariantTable(config,'aa','BMS-790052_BMS-650032','undetectable_in_absence_of_therapy','NS3aa36')

#makeVariantTables <- function(config, type, group=NULL, ...)
#{
#	if (is.null(group))
#		groups <- config@groups
#	else groups <- splitFields(group)
#	tables <- list()
#	for (group in groups)
#	{
#		tables[[group]] <- list()
#		for (region in getRegionsForGroup(config,group))
#		{
#			try({
#			print(concat('group=',group,', region=',region))
#			tbl <- makeVariantTable(config, type, group, region, ...)
#			tables[[group]][[region]] <- tbl
#			}, silent=FALSE)
#		}
#	}
#	return(tables)
#}
##tables <- makeVariantTables(config, 'codons', 'PXB0220-0030')

makeVariantTables <- function(config, type, groups=config@groups, ...)
{
	tables <- list()
	for (group in groups)
	{
		tables[[group]] <- list()
		for (subgroup in getTablesForGroup(config,group))
		{
			tables[[group]][[subgroup]] <- list()
			for (region in getRegionsForSubGroup(config,group,subgroup))
			{
				#try({
					print(concat('group=',group,', subgroup=',subgroup,', region=',region))
					tbl <- makeVariantTable(config, type, group, subgroup, region, ...)
					tables[[group]][[subgroup]][[region]] <- tbl
				#}, silent=FALSE)
			}
		}
	}
	return(tables)
}
#codon.tables <- makeVariantTables(config,'codons','MP-424',cutoff=2)
#tables <- makeVariantTables(config, 'codons', 'BMS-790052_BMS-650032')

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

writeCodonTable <- function(config, codon.tables, group, subgroup, region)
{
	tbl <- codon.tables[[group]][[subgroup]][[region]]
	ref <- getRefForGroup(config,group)
	identifier <- concat(group,'-',subgroup,'-',region,'-',ref)
	filename <- concat(config@tables.dir,'/table-',identifier,'.txt')
	try(colnames(tbl)[1] <-identifier, silent=FALSE)
	writeTable(tbl,filename,row.names=FALSE)
}

writeCodonTables <- function(config, groups=config@groups, cutoff=2)
{ 
	codon.tables <- makeVariantTables(config,'codons',groups,cutoff=cutoff)
	for (group in groups)
	{
		for (subgroup in getTablesForGroup(config,group))
		{
			print(subgroup)
			for (region in getRegionsForSubGroup(config,group,subgroup))
			{
				writeCodonTable(config, codon.tables, group, subgroup, region)
			}
		}
	}
	return(codon.tables)
}
#codon.tables <- writeCodonTables(config,'MP-424')
#writeCodonTables(config,'confirm_with_new_reagents')
#writeCodonTables(config,'PXB0220-0030')

#################################################################

#require(XLConnect)

