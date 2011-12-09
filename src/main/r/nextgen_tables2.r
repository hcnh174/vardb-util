makeCodonVariantTable <- function(config, samples, region, aanum, minreads=0, show.total=TRUE, show.freq=TRUE)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	ref <- getRefForSamples(config,samples)
	#gene <- config@regions[region,'gene']
	#positions <- getCodonPositionsForGene(config,gene,ref)
	#refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	refcodon <- getReferenceCodon(config, ref, region, aanum)
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
	
	if (show.total)
	{
		totalrow <- makeRow(counts,list(codon='total'))
		totalrow[,countcols] <- sums
		counts <- rbind(counts,totalrow)
	}
	
	if (show.freq)
	{
		freqrow <- makeRow(counts,list(codon='maf'))
		freqrow[,countcols] <- format(mafs/sums, digits=3)
		counts <- rbind(counts,freqrow)
	}
	
	return(counts)
}
#counts <- makeCodonVariantTable(config, getSamplesForSubGroup(config,'hcv_infection','hcv_infection'),  'NS5Aaa31', 31)
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
					samples <- getSamplesForSubGroup(config,group,subgroup)
					tbl <- makeCodonVariantTable(config, samples, region, aanum, minreads)
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

makeAminoAcidVariantTable <- function(config, samples, region, aanum, minreads=0, show.total=TRUE, show.freq=TRUE)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	ref <- getRefForSamples(config,samples)
	#gene <- config@regions[region,'gene']
	#positions <- getCodonPositionsForGene(config,gene,ref)
	#refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	#refaa <- translateCodon(refcodon)
	refaa <- getReferenceAminoAcid(config, ref, region, aanum)
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
	
	if (show.total)
	{
		totalrow <- makeRow(counts,list(aa='total'))
		totalrow[,countcols] <- sums
		counts <- rbind(counts,totalrow)
	}
	
	if (show.freq)
	{
		freqrow <- makeRow(counts,list(aa='maf'))
		freqrow[,countcols] <- format(mafs/sums, digits=3)
		counts <- rbind(counts,freqrow)
	}
	
	try(rownames(counts) <- counts$aa)
	return(counts)
}
#counts <- makeAminoAcidVariantTable(config, getSamplesForSubGroup(config,'hcv_infection','hcv_infection'),  'NS5Aaa31', 31)

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
					samples <- getSamplesForSubGroup(config,group,subgroup)
					tbl <- makeAminoAcidVariantTable(config, samples, region, aanum, minreads)
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


outputTablesToSpreadsheet <- function(config, groups=config@groups, minreads=config@minreads)
{
	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
	filename <- 'tables.xlsx'
	deleteFile(filename)
	wb <- loadWorkbook(filename, create = TRUE)
	for (group in groups)
	{
		sheet <- fixSheetName(group)
		createSheet(wb, name = sheet)
		setCellText(wb,sheet,group)
		ref <- getRefForGroup(config,group)
		#printcat('ref: ',ref)
		for (subgroup in getTablesForGroup(config,group))
		{
			#printcat(' subgroup: ',subgroup)
			for (region in getRegionsForSubGroup(config,group,subgroup))
			{
				#printcat('  region: ',subgroup)
				gene <- strsplit(region,'aa', fixed=TRUE)[[1]][1]
				for (aanum in getFociForRegion(config,region))
				{
					#printcat('   aanum: ',aanum)
					samples <- getSamplesForSubGroup(config,group,subgroup)
					tbl <- makeAminoAcidVariantTable(config, samples, region, aanum, minreads, show.total=FALSE, show.freq=TRUE)
					title <- concat(subgroup,' ',gene,'aa',aanum)
					writeTableToWorksheet(wb, sheet, tbl, title=title, footer=TRUE)
				}
			}
		}
	}
	saveWorkbook(wb)
}
#outputTablesToSpreadsheet(config,'hcv_infection',minreads=100)
#outputTablesToSpreadsheet(config,minreads=100)


makeReferenceVsVariantTable <- function(config, group, region, aanum, minreads=0)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	samples <- getSamplesForGroup(config,group)
	ref <- getRefForSamples(config,samples)
	data <- getCodonCountSubset(config, samples, region, 'aa', aanum,minreads=minreads)
	data$sample <- sapply(data$sample,function(sample){return(getStemForSample(sample))})
	refaa <- getReferenceAminoAcid(config, ref, region, aanum)
	data[which(data$aa==refaa),'aa'] <- concat(refaa,'*')
	refaa <- concat(refaa,'*')
	counts <- cast(data, sample ~ aa, value='count', fun.aggregate=function(x) return(x[1]))
	if (is.null(counts[[refaa]]))
		counts[[refaa]] <- 0
	colsums <- data.frame()
	for (col in colnames(counts)[-1])
	{
		colsums[col,'aa'] <- col
		colsums[col,'sum'] <- sum(counts[[col]], na.rm=TRUE)
	}
	#print(colsums)
	cols <- colsums[order(colsums$sum, decreasing=TRUE),'aa']
	counts <- counts[,c('sample',cols)]
	filename <- concat(config@tables.dir,'/',group,'-bysubject.',region,'-',aanum,'.txt')
	writeTable(counts,filename)
	return(counts)
}
#counts <- makeReferenceVsVariantTable(config,'MP-424','NS3aa36', 36)

makeReferenceVsVariantTables <- function(config, groups=config@groups, minreads=0)
{
	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
	filename <- 'tables-by-subject.xlsx'
	deleteFile(filename)
	wb <- loadWorkbook(filename, create = TRUE)
	for (group in groups)
	{
		sheet <- fixSheetName(group)
		createSheet(wb, name = sheet)
		setCellText(wb,sheet,group)
		for (region in getRegionsForGroup(config,group))
		{
			gene <- strsplit(region,'aa', fixed=TRUE)[[1]][1]
			for (aanum in getFociForRegion(config,region))
			{
				tbl <- makeReferenceVsVariantTable(config,group,region,aanum,minreads=minreads)
				title <- concat(group,' ',gene,'aa',aanum)
				writeTableToWorksheet(wb, sheet, tbl, title=title)
			}
		}
	}
	saveWorkbook(wb)
}
#makeReferenceVsVariantTables(config, minreads=100)
#makeReferenceVsVariantTables(config, groups='MP-424', minreads=100)
