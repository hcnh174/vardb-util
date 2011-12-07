#getCodonCountSubset <- function(config, samples, region, filetype, start, end=start, minreads=0)
#{
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
#	data.subset <- data.subset[which(data.subset$region %in% splitFields(region)),]
#	data.subset$column <- factor(data.subset$column)
#	data.subset$aanum <- factor(data.subset$aanum)
#	return(data.subset)
#}
##getCodonCountSubset(config,getSamplesForSubGroup(config,'hcv_infection','hcv_infection'), 'NS5Aaa31', 'aa', 31)

makeCodonVariantTable <- function(config, samples, region, aanum, minreads=0, show.total=TRUE, show.freq=TRUE)
{
	require(reshape, quietly=TRUE, warn.conflicts=FALSE)
	gene <- config@regions[region,'gene']
	ref <- getRefForSamples(config,samples)
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
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
	gene <- config@regions[region,'gene']
	ref <- getRefForSamples(config,samples)
	#ref <- getRefForGroup(config,group)
	positions <- getCodonPositionsForGene(config,gene,ref)
	refcodon <- toupper(as.character(positions[which(positions$codon==aanum),'refcodon']))
	refaa <- translateCodon(refcodon)
	#samples <- getSamplesForSubGroup(config,group,subgroup)
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


createHeaderStyle <- function(wb)
{
	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
	style <- createCellStyle(wb)
	setFillPattern(style, fill = XLC$"FILL.SOLID_FOREGROUND")
	setFillBackgroundColor(style, color = XLC$"COLOR.WHITE")
	setFillForegroundColor(style, color = XLC$"COLOR.GREY_25_PERCENT")
	setBorder(style, side = c("all"), type = XLC$"BORDER.THIN", color = c(XLC$"COLOR.AUTOMATIC"))
	return(style)
}

createBorderStyle <- function(wb)
{
	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
	style <- createCellStyle(wb)
	setBorder(style, side = c("all"), type = XLC$"BORDER.THIN", color = c(XLC$"COLOR.AUTOMATIC"))
	return(style)
}

fixSheetName <- function(name)
{
	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
	if (nchar(name)>31)
		name <- substr(name,1,31)
	return(name)
}

setCellText <- function(wb, sheet, text, row=1, col=1)
{
	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
	writeWorksheet(wb,text,sheet,startRow=row,startCol=col,header=FALSE)
}
#
#writeTableToWorksheet <- function(wb, sheet, tbl, startRow=NULL, startCol=1, title=NULL,
#		style=createBorderStyle(wb),
#		headerstyle=createHeaderStyle(wb),
#		footerstyle=createHeaderStyle(wb))
#{
#	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
#	if (is.null(startRow))
#		startRow <- getLastRow(wb, sheet)+2
#	if (!is.null(title))
#	{
#		setCellText(wb,sheet,title,startRow)
#		startRow <- getLastRow(wb, sheet)+1
#	}
#	writeWorksheet(wb,tbl,sheet,startRow)
#	endRow <- getLastRow(wb, sheet)
#	startCol <- 1
#	endCol <- ncol(tbl)
#	#printcat('startRow=',startRow,', endRow=',endRow)
#	colnums <- c()
#	rownums <- c()
#	for (col in seq(startCol, endCol))
#	{
#		for (row in seq(startRow,endRow))
#		{
#			#printcat('row=',row,' col=',col)
#			if (row==startRow)
#				setCellStyle(wb,sheet,row,col,headerstyle)
#			else if (row==endRow)
#				setCellStyle(wb,sheet,row,col,footerstyle)
#			else setCellStyle(wb,sheet,row,col,style)
#		}
#	}
#}

styleCells <- function(wb, sheet, startRow, startCol, endRow, endCol, style)
{
	colnums <- c()
	rownums <- c()
	for (col in seq(startCol, endCol))
	{
		for (row in seq(startRow,endRow))
		{
			rownums <- c(rownums, row)
			colnums <- c(colnums, col)			
		}
	}
	setCellStyle(wb,sheet,rownums,colnums,style)
}

writeTableToWorksheet <- function(wb, sheet, tbl, startRow=NULL, startCol=1, title=NULL,
		style=createBorderStyle(wb),
		headerstyle=createHeaderStyle(wb),
		footerstyle=createHeaderStyle(wb))
{
	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
	if (is.null(startRow))
		startRow <- getLastRow(wb, sheet)+2
	if (!is.null(title))
	{
		setCellText(wb,sheet,title,startRow)
		startRow <- getLastRow(wb, sheet)+1
	}
	writeWorksheet(wb,tbl,sheet,startRow)
	endRow <- getLastRow(wb, sheet)
	startCol <- 1
	endCol <- ncol(tbl)
	
	styleCells(wb,sheet,startRow,startCol,endRow,endCol,style)
	styleCells(wb,sheet,startRow,startCol,startRow,endCol,headerstyle)
	styleCells(wb,sheet,endRow,startCol,endRow,endCol,footerstyle)
	
	
	#printcat('startRow=',startRow,', endRow=',endRow)
#	colnums <- c()
#	rownums <- c()
#	for (col in seq(startCol, endCol))
#	{
#		for (row in seq(startRow,endRow))
#		{
#			rownums <- c(rownums, row)
#			colnums <- c(colnums, col)			
#		}
#	}
#	printcat('colnums: ',colnums)
#	printcat('rownums: ',rownums)
#	
#	setCellStyle(wb,sheet,rownums,colnums,style)
	
#	for (col in seq(startCol, endCol))
#	{
#		for (row in seq(startRow,endRow))
#		{
#			#printcat('row=',row,' col=',col)
#			if (row==startRow)
#				setCellStyle(wb,sheet,row,col,headerstyle)
#			else if (row==endRow)
#				setCellStyle(wb,sheet,row,col,footerstyle)
#			else setCellStyle(wb,sheet,row,col,style)
#		}
#	}
}

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
					writeTableToWorksheet(wb, sheet, tbl, title=title)
				}
			}
		}
	}
	saveWorkbook(wb)
}
#outputTablesToSpreadsheet(config,'hcv_infection',minreads=100)
#outputTablesToSpreadsheet(config,minreads=100)

#
#outputTablesToSpreadsheet <- function(config, groups=config@groups, minreads=config@minreads)
#{
#	require(XLConnect, quietly=TRUE, warn.conflicts=FALSE)
#	filename <- 'tables.xlsx'
#	deleteFile(filename)
#	wb <- loadWorkbook(filename, create = TRUE)
#	
#	headerstyle <- createCellStyle(wb)
#	setFillPattern(headerstyle, fill = XLC$"FILL.SOLID_FOREGROUND")
#	setFillBackgroundColor(headerstyle, color = XLC$"COLOR.WHITE")
#	setFillForegroundColor(headerstyle, color = XLC$"COLOR.GREY_25_PERCENT")
#	setBorder(headerstyle, side = c("all"), type = XLC$"BORDER.THIN", color = c(XLC$"COLOR.AUTOMATIC"))
#	
#	footerstyle <- createCellStyle(wb)
#	setFillPattern(footerstyle, fill = XLC$"FILL.SOLID_FOREGROUND")
#	setFillBackgroundColor(footerstyle, color = XLC$"COLOR.WHITE")
#	setFillForegroundColor(footerstyle, color = XLC$"COLOR.GREY_25_PERCENT")
#	setBorder(footerstyle, side = c("all"), type = XLC$"BORDER.THIN", color = c(XLC$"COLOR.AUTOMATIC"))
#
#	borderstyle <- createCellStyle(wb)
#	setBorder(borderstyle, side = c("all"), type = XLC$"BORDER.THIN", color = c(XLC$"COLOR.AUTOMATIC"))
#
#	for (group in groups)
#	{
#		sheet <- group
#		if (nchar(sheet)>31)
#			sheet <- substr(sheet,1,31)
#		createSheet(wb, name = sheet)
#		writeWorksheet(wb,group,sheet,startRow=1,startCol=1,header=FALSE)
#		ref <- getRefForGroup(config,group)
#		printcat('ref: ',ref)
#		for (subgroup in getTablesForGroup(config,group))
#		{
#			printcat(' subgroup: ',subgroup)
#			for (region in getRegionsForSubGroup(config,group,subgroup))
#			{
#				printcat('  region: ',subgroup)
#				gene <- strsplit(region,'aa', fixed=TRUE)[[1]][1]
#				for (aanum in getFociForRegion(config,region))
#				{
#					printcat('   aanum: ',aanum)
#					samples <- getSamplesForSubGroup(config,group,subgroup)
#					tbl <- makeAminoAcidVariantTable(config, samples, region, aanum, minreads, show.total=FALSE, show.freq=TRUE)
#					identifier <- concat(subgroup,' ',gene,'aa',aanum)
#					
#					startRow <- getLastRow(wb, sheet)+2
#					writeWorksheet(wb,identifier,sheet,startRow=startRow,startCol=1,header=FALSE)
#					startRow <- getLastRow(wb, sheet)+1
#					writeWorksheet(wb,tbl,sheet,startRow)
#					endRow <- getLastRow(wb, sheet)
#					startCol <- 1
#					endCol <- ncol(tbl)
#					#printcat('startRow=',startRow,', endRow=',endRow)
#					for (col in seq(startCol, endCol))
#					{
#						for (row in seq(startRow,endRow))
#						{
#							#printcat('row=',row,' col=',col)
#							if (row==startRow)
#								setCellStyle(wb,sheet,row,col,headerstyle)
#							else if (row==endRow)
#								setCellStyle(wb,sheet,row,col,footerstyle)
#							else setCellStyle(wb,sheet,row,col,borderstyle)
#						}
#					}
#				}
#			}
#		}
#	}
#	
#	saveWorkbook(wb)
#}
##outputTablesToSpreadsheet(config,'hcv_infection',minreads=100)



