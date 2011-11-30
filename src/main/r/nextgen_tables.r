#getCodonCountSubset <- function(config, group, subgroup, region, filetype, start, end=start, minreads=0)
#{
#	samples <- getSamplesForSubGroup(config,group,subgroup)
#	data <- NULL
#	for (sample in samples)
#	{
#		filename <- getCodonCountFilename(config,sample,filetype)
#		data.sample <- loadDataFrame(filename)
#		if (is.null(data))
#			data <- data.sample
#		else data <- rbind(data,data.sample)
#	}
#	data.subset <- data[which(data$region==region & data$aanum>=start & data$aanum<=end & data$count>=minreads),]
#	data.subset$column <- factor(data.subset$column)
#	data.subset$aanum <- factor(data.subset$aanum)
#	return(data.subset)
#}
##getCodonCountSubset(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy','NS3aa36','codons',)
#
#getCodonCountSubset <- function(config, samples, region, filetype, start, end=start, minreads=0)
#{
#	#samples <- getSamplesForSubGroup(config,group,subgroup)
#	data <- NULL
#	for (sample in samples)
#	{
#		filename <- getCodonCountFilename(config,sample,filetype)
#		data.sample <- loadDataFrame(filename)
#		if (is.null(data))
#			data <- data.sample
#		else data <- rbind(data,data.sample)
#	}
#	data.subset <- data[which(data$aanum>=start & data$aanum<=end & data$count>=minreads),]
#	data.subset <- data[which(data$region %in% splitFields(region)),]
#	data.subset$column <- factor(data.subset$column)
#	data.subset$aanum <- factor(data.subset$aanum)
#	return(data.subset)
#}
##getCodonCountSubset(config,'BMS-790052_BMS-650032','undetectable_in_absence_of_therapy','NS3aa36','codons',)

getCodonCountSubset <- function(config, samples, region, filetype, start, end=start, minreads=0)
{
	data <- NULL
	for (sample in samples)
	{
		filename <- getCodonCountFilename(config,sample,filetype)
		data.sample <- loadDataFrame(filename)
		if (is.null(data))
			data <- data.sample
		else data <- rbind(data,data.sample)
	}
	data.subset <- data[which(data$aanum>=start & data$aanum<=end & data$count>=minreads),]
	data.subset <- data.subset[which(data.subset$region %in% splitFields(region)),]
	data.subset$column <- factor(data.subset$column)
	data.subset$aanum <- factor(data.subset$aanum)
	return(data.subset)
}
#getCodonCountSubset(config,getSamplesForSubGroup(config,'hcv_infection','hcv_infection'), 'NS5Aaa31', 'aa', 31)

makeCodonVariantTable <- function(config, group, subgroup, region, aanum, minreads=0)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	gene <- config@regions[region,'gene']
	ref <- getRefForGroup(config,group)
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	samples <- getSamplesForSubGroup(config,group,subgroup)
	#data.subset <- getCodonCountSubset(config,group,subgroup,region,'codons',aanum,minreads=minreads)
	data.subset <- getCodonCountSubset(config,samples,region,'codons',aanum,minreads=minreads)
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
	
	freqrow <- makeRow(counts,list(codon='maf'))
	freqrow[,countcols] <- format(mafs/sums, digits=3)
	counts <- rbind(counts,freqrow)

	return(counts)
}
#group <- 'MP-424'; subgroup <- 'undetectable_in_absence_of_therapy'; region <- 'NS3aa36'; aanum <- 36
#makeCodonVariantTable(config,group,subgroup,region,aanum,50)


###############################################

writeCodonTables <- function(config, groups=config@groups, minreads=config@minreads)
{
	printcat('using minreads: ',minreads)
	tbls <- list()
	for (group in groups)
	{
		ref <- getRefForGroup(config,group)
		printcat('ref for group: ',ref)
		for (subgroup in getTablesForGroup(config,group))
		{
			printcat(' subgroup: ',subgroup)
			for (region in getRegionsForSubGroup(config,group,subgroup))
			{			
				printcat('  region: ',region)
				for (aanum in getFociForRegion(config,region))
				{
					printcat('   aanum: ',aanum)
					tbl <- makeCodonVariantTable(config, group, subgroup, region, aanum, minreads)
					identifier <- concat(group,'-',subgroup,'-',region,'-',aanum,'-',ref)
					filename <- concat(config@tables.dir,'/table-codons-',identifier,'.txt')
					writeTable(tbl,filename,row.names=FALSE)
					tbls[[identifier]] <- tbl
				}
			}
		}
	}
	return(tbls)
}
#writeCodonTables(config,'hcv_infection')
#writeCodonTables(config,'MP-424')
#writeCodonTables(config,'confirm_with_new_reagents')
#writeCodonTables(config,'PXB0220-0030')


###############################################
#
#makeAminoAcidVariantTable <- function(config, group, subgroup, region, aanum, minreads=0)
#{
#	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
#	gene <- config@regions[region,'gene']
#	ref <- getRefForGroup(config,group)
#	positions <- getCodonPositionsForGene(config,gene,ref)
#	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
#	refaa <- translateCodon(refcodon)
#	samples <- getSamplesForSubGroup(config,group,subgroup)
#	#data.subset <- getCodonCountSubset(config,group,subgroup,region,'aa',aanum,minreads=minreads)
#	data.subset <- getCodonCountSubset(config,samples,region,'aa',aanum,minreads=minreads)
#	return(data.subset)
#	
#	counts <- cast(data.subset, aa ~ column, value='count', fun.aggregate=function(x) return(x[1]))
#	if (length(which(counts$aa==refaa))==0)
#		counts <- addRow(counts, list(aa=refaa))
#	counts$isref <- ifelse(counts$aa==refaa,1,0)
#	counts[which(counts$aa==refaa),'aa'] <- concat(refaa,'*')
#	counts <- counts[order(counts$isref,counts$aa,decreasing=TRUE),]
#	print(colnames(counts))
#	#return(counts)
#	countcols <- removeElements(colnames(counts),'aa,isref')
#	counts <- counts[,c('aa',countcols)]
#
#	mafs <- c(); sums <- c()
#	for (col in countcols)
#	{
#		colsum <- sum(counts[[col]], na.rm=TRUE)
#		if (length(colsum)==0) colsum <- 0
#		sums <- c(sums,colsum)
#		
#		maf <- sort(as.numeric(names(xtabs(~counts[[col]]))), decreasing=TRUE)[2]
#		if (is.na(maf)) maf <- 0
#		mafs <- c(mafs,maf)
#	}
#	
#	totalrow <- makeRow(counts,list(aa='total'))
#	totalrow[,countcols] <- sums
#	counts <- rbind(counts,totalrow)
#	
#	freqrow <- makeRow(counts,list(aa='maf'))
#	freqrow[,countcols] <- format(mafs/sums, digits=3)
#	counts <- rbind(counts,freqrow)
#	
#	return(counts)
#}
##group <- 'hcv_infection'; subgroup <- 'hcv_infection'; region <- 'NS3'; aanum <- 36
##group <- 'confirm_plasmid_with_new_reagents'; subgroup <- 'confirm_plasmid_with_new_reagents'; region <- 'NS3'; aanum <- 36
##counts <- makeAminoAcidVariantTable(config,group,subgroup,region,aanum)
##counts


makeAminoAcidVariantTable <- function(config, group, subgroup, region, aanum, minreads=0)
{	
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	gene <- config@regions[region,'gene']
	ref <- getRefForGroup(config,group)
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	refaa <- translateCodon(refcodon)
	samples <- getSamplesForSubGroup(config,group,subgroup)
	data.subset <- getCodonCountSubset(config,samples,region,'aa',aanum,minreads=minreads)	
	counts <- cast(data.subset, aa ~ column, value='count', fun.aggregate=function(x) return(x[1]))
	if (length(which(counts$aa==refaa))==0)
		counts <- addRow(counts, list(aa=refaa))
	counts$isref <- ifelse(counts$aa==refaa,1,0)
	counts[which(counts$aa==refaa),'aa'] <- concat(refaa,'*')
	counts <- counts[order(counts$isref,counts$aa,decreasing=TRUE),]
	countcols <- removeElements(colnames(counts),'aa,isref')
	counts <- counts[,c('aa',countcols)]	
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
	
	totalrow <- makeRow(counts,list(aa='total'))
	totalrow[,countcols] <- sums
	counts <- rbind(counts,totalrow)
	
	freqrow <- makeRow(counts,list(aa='maf'))
	freqrow[,countcols] <- format(mafs/sums, digits=3)
	counts <- rbind(counts,freqrow)
	
	return(counts)
}
#counts <- makeAminoAcidVariantTable(config, 'hcv_infection', 'hcv_infection',  'NS5Aaa31', 31)



writeAminoAcidTables <- function(config, groups=config@groups, minreads=config@minreads)
{
	printcat('using minreads: ',minreads)
	tbls <- list()
	for (group in groups)
	{
		ref <- getRefForGroup(config,group)
		printcat('ref: ',ref)
		for (subgroup in getTablesForGroup(config,group))
		{
			printcat(' subgroup: ',subgroup)
			for (region in getRegionsForSubGroup(config,group,subgroup))
			{
				printcat('  region: ',subgroup)
				for (aanum in getFociForRegion(config,region))
				{
					printcat('   aanum: ',aanum)
					tbl <- makeAminoAcidVariantTable(config, group, subgroup, region, aanum, minreads)
					identifier <- concat(group,'-',subgroup,'-',region,'-',aanum,'-',ref)
					filename <- concat(config@tables.dir,'/table-aa-',identifier,'.txt')
					writeTable(tbl,filename,row.names=FALSE)
					tbls[[identifier]] <- tbl
				}
			}
		}
	}
	return(tbls)
}
#tables <- writeAminoAcidTables(config,'hcv_infection')

#################################################################

#require(XLConnect)

