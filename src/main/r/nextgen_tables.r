getCodonCountSubset <- function(config, group, subgroup, region, filetype, start, end=start, cutoff=0)
{
	samples <- getSamplesForSubGroup(config,group,subgroup)
	data <- NULL
	for (sample in samples)
	{
		filename <- getCodonCountFilename(config,sample,filetype)
		data.sample <- loadDataFrame(filename)
		if (is.null(data))
			data <- data.sample
		else data <- rbind(data,data.sample)
	}
	data.subset <- data[which(data$region==region & data$aanum>=start & data$aanum<=end & data$count>=cutoff),]
	data.subset$column <- factor(data.subset$column)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy','NS3aa36','codons',)

makeCodonVariantTable <- function(config, group, subgroup, region, aanum, cutoff=0)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	gene <- config@regions[region,'gene']
	ref <- getRefForGroup(config,group)
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	data.subset <- getCodonCountSubset(config,group,subgroup,region,'codons',aanum,cutoff=cutoff)
	counts <- cast(data.subset, codon ~ column, value='count', fun.aggregate=function(x) return(x[1])); counts
	if (length(which(counts$codon==refcodon))==0)
		counts <- addRow(counts, list(codon=refcodon))
	counts$aa <- sapply(counts$codon,translateCodon)
	counts$isref <- ifelse(counts$codon==refcodon,1,0)
	counts[which(counts$codon==refcodon),'codon'] <- concat(refcodon,'*')
	counts <- counts[order(counts$isref,counts$aa,counts$codon, decreasing=TRUE),]
	countcols <- removeElements(colnames(counts),'aa,codon,isref')
	counts <- counts[,c('codon','aa',countcols)]
	
	mafs <- c(); sums <- c()
	for (col in countcols)
	{
		colsum <- sum(counts[[col]], na.rm=TRUE)
		if (length(colsum)==0) colsum <- 0
		sums <- c(sums,colsum)
		
		maf <- sort(as.numeric(names(xtabs(~counts[[col]]))), decreasing=TRUE)[2]
		if (is.na(maf)) maf <- 0
		mafs <- c(mafs,maf)
	}

	totalrow <- makeRow(counts,list(codon='total'))
	totalrow[,countcols] <- sums
	counts <- rbind(counts,totalrow)
	
	freqrow <- makeRow(counts,list(codon='freq'))
	freqrow[,countcols] <- format(mafs/sums, digits=3)
	counts <- rbind(counts,freqrow)

	return(counts)
}
#group <- 'MP-424'; subgroup <- 'undetectable_in_absence_of_therapy'; region <- 'NS3aa36'; aanum <- 36
#makeCodonVariantTable(config,group,subgroup,region,aanum,50)

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

writeCodonTables <- function(config, groups=config@groups, cutoff=2)
{
	tbls <- list()
	for (group in groups)
	{
		ref <- getRefForGroup(config,group)
		for (subgroup in getTablesForGroup(config,group))
		{
			#print(subgroup)
			for (region in getRegionsForSubGroup(config,group,subgroup))
			{
				for (aanum in getFociForRegion(config,region))
				{
					tbl <- makeCodonVariantTable(config, group, subgroup, region, aanum, cutoff)
					identifier <- concat(group,'-',subgroup,'-',region,'-',aanum,'-',ref)
					filename <- concat(config@tables.dir,'/table-',identifier,'.txt')
					writeTable(tbl,filename,row.names=FALSE)
					tbls[[identifier]] <- tbl
				}
			}
		}
	}
	return(tbls)
}
#writeCodonTables(config,'MP-424')
#writeCodonTables(config,'confirm_with_new_reagents')
#writeCodonTables(config,'PXB0220-0030')

#################################################################

#require(XLConnect)

