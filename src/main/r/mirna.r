library(reshape)
library(RmiR.Hs.miRNA)
library(org.Hs.eg.db)
library(topGO)

makeHeatmapMatrix <- function(data, use.log2=TRUE)
{
	#data <- data[3:length(colnames(data))]
	data <- data[complete.cases(data),]
	if (use.log2)
	{
		for (col in colnames(data))
		{
			data[[col]] <- log2(data[[col]])
		}
	}
	data.matrix <- as.matrix(data)
	#data.matrix <- data.matrix[complete.cases(data.matrix),]
	return(data.matrix)
}

drawHeatmap <- function(data, cexRow=0.2, use.log2=FALSE, main='Heatmap')
{
	data.matrix <- makeHeatmapMatrix(data, use.log2)
	hv <- heatmap(data.matrix, na.rm=T, ylab='miRNAs', main=main, Colv=TRUE, Rowv=TRUE, cexRow=cexRow, margins=c(10,5))
}
#drawHeatmap(data.human)

drawHeatmap2 <- function(data, cexRow=0.2, use.log2=FALSE, main='Heatmap', distmethod='euclidean', clustmethod='ward', ...)
{
	require(gplots)
	data.matrix <- makeHeatmapMatrix(data, use.log2)
	heatmap.2(data.matrix, col=redgreen(75), scale='row', main=main, #ColSideColors=patientcolors,
			key=TRUE, symkey=FALSE, density.info='none', trace='none', cexRow=cexRow, margins=c(10,5),
			distfun=function(x){dist(x, method=distmethod)},# 'euclidean,maximum,manhattan,canberra,binary,minkowski'			
			hclustfun=function(x){hclust(x, method=clustmethod)}, #'ward,single,complete,average,mcquitty,median,centroid'
			...
	)
}
#drawHeatmap2(data.human)

getGeneSymbols <- function()
{
	mapped_genes <- mappedkeys(org.Hs.egSYMBOL)
	symbols <- as.list(org.Hs.egSYMBOL[mapped_genes])
	return(symbols)
}
#symbols <- getGeneSymbols()

getGeneNames <- function()
{
	mapped_genes <- mappedkeys(org.Hs.egGENENAME)
	genenames <- as.list(org.Hs.egGENENAME[mapped_genes])
	return(genenames)
}

annotateTargets <- function(targets)
{
	symbols <- getGeneSymbols()
	#genenames <- getGeneNames()
	for (gene_id in targets$gene_id)
	{
		#print(gene_id)
		try({
			symbol <- symbols[gene_id]
			if (!is.null(symbol))
			{
				print(concat(gene_id,': ',symbol))
				targets[which(targets$gene_id==gene_id),'gene'] <- symbol
			}
		})
		#targets[which(targets$gene_id==gene_id),'description'] <- genenames[gene_id]
	}
	return(targets)
}

compareTreatments <- function(data.melted, trt1, trt2, pvalue=0.05, cutoff=NULL)
{
	data.subset <- subset(data.melted, trt==trt1 | trt==trt2)
	#data.subset <- data.subset[which(data.subset$value!=0),]
	#data.subset$value <- log2(data.subset$value)
	data <- cast(data.subset, mirna ~ trt, value='value', fun.aggregate=median)#mean
	data$ratio <- log2(data[[trt1]]/data[[trt2]])
	rownames(data) <- data$mirna
	for (mirna in rownames(data))
	{
		values1 <- data.subset[which(data.subset$mirna==mirna & data.subset$trt==trt1),'value']
		values2 <- data.subset[which(data.subset$mirna==mirna & data.subset$trt==trt2),'value']
		fit <- t.test(values1,values2)
		data[mirna,'mean1'] <- mean(values1)
		data[mirna,'mean2'] <- mean(values2)
		data[mirna,'pvalue'] <- fit$p.value
	}
	data <- data[order(data$pvalue),]
	data <- data[which(data$pvalue<=pvalue),]
	if (!is.null(cutoff))
	{
		if (cutoff<1)
			data <- data[which(data$ratio <= cutoff),]
		else data <- data[which(data$ratio >= cutoff),]
	}
	return(data)
}
#compareTreatments(data.melted,'healthy','hcc', cutoff=1.5)
#
#getMiRNATargets <- function(data, name, at.least=3, databases=splitFields('miranda,mirbase,mirtarget2,pictar,tarbase,targetscan'))
#{
#	#data <- makeContrast(data.melted, trt1, trt2)
#	#data <- data[which(data$pvalue<=pvalue),]
#	targets <- NULL
#	for (db in databases)
#	{
#		rows <- dbGetQuery(RmiR.Hs.miRNA_dbconn(), 
#				concat("SELECT * FROM ",db," WHERE mature_miRNA='",
#						joinFields(data$mirna, delimiter="' OR mature_miRNA='"),"'"))
#		if (nrow(rows)>0)
#		{
#			rows$db <- db
#			rows <- rows[,splitFields('mature_miRNA,gene_id,db')]
#			rows <- unique(rows)
#			print(head(rows))
#			if (is.null(targets))
#				targets <- rows
#			else targets <- rbind(targets,rows)
#		}
#	}
#	targets$pair <- concat(targets$mature_miRNA,':',targets$gene_id)
#	filename <- concat('out/targets-',name,'-all.txt')
#	writeTable(targets, filename, row.names=FALSE)
#	counts <- xtabs(~pair, targets)
#	counts <- sort(counts, decreasing=TRUE)
#	pairs <- names(counts[which(counts>=at.least)])
#	targets2 <- unique(targets[which(targets$pair %in% pairs), splitFields('mature_miRNA,gene_id')])
#	targets2 <- annotateTargets(targets2)
#	filename <- concat('out/targets-',name,'.txt')
#	writeTable(targets2, filename, row.names=FALSE)
#	
#	targets3 <- targets2[which(!is.na(targets2$gene)),]
#	attributes <- data.frame()
#	attributes <- rbind(attributes,data.frame(node_id=targets3$gene, node_type='gene'))
#	attributes <- rbind(attributes,data.frame(node_id=targets3$mature_miRNA, node_type='mirna'))
#	filename <- concat('out/attributes-',name,'.txt')
#	writeTable(attributes, filename, row.names=FALSE)
#	
#	return(targets2)
#}

getMiRNATargets <- function(mirnas, name, at.least=3, databases=splitFields('miranda,mirbase,mirtarget2,pictar,tarbase,targetscan'))
{
	targets <- NULL
	for (db in databases)
	{
		rows <- dbGetQuery(RmiR.Hs.miRNA_dbconn(), 
				concat("SELECT * FROM ",db," WHERE mature_miRNA='",
						joinFields(mirnas, delimiter="' OR mature_miRNA='"),"'"))
		if (nrow(rows)>0)
		{
			rows$db <- db
			rows <- rows[,splitFields('mature_miRNA,gene_id,db')]
			rows <- unique(rows)
			print(head(rows))
			if (is.null(targets))
				targets <- rows
			else targets <- rbind(targets,rows)
		}
	}
	if (is.null(targets))
		throw('no targets found for any database: Are miRNA names in the form hsa-miR-122?')
	targets$pair <- concat(targets$mature_miRNA,':',targets$gene_id)
	filename <- concat('out/targets-',name,'-all.txt')
	writeTable(targets, filename, row.names=FALSE)
	counts <- xtabs(~pair, targets)
	counts <- sort(counts, decreasing=TRUE)
	pairs <- names(counts[which(counts>=at.least)])
	targets2 <- unique(targets[which(targets$pair %in% pairs), splitFields('mature_miRNA,gene_id')])
	targets2 <- annotateTargets(targets2)
	filename <- concat('out/targets-',name,'.txt')
	writeTable(targets2, filename, row.names=FALSE)
	
	genes <- unique(targets2$gene)
	filename <- concat('out/genes-',name,'.txt')
	writeTable(genes, filename, row.names=FALSE)
	
	targets3 <- targets2[which(!is.na(targets2$gene)),]
	attributes <- data.frame()
	attributes <- rbind(attributes,data.frame(node_id=targets3$gene, node_type='gene'))
	attributes <- rbind(attributes,data.frame(node_id=targets3$mature_miRNA, node_type='mirna'))
	filename <- concat('out/attributes-',name,'.txt')
	writeTable(attributes, filename, row.names=FALSE)
	
	return(targets2)
}

loadGOdata <- function(targets, ontology='BP')
{
	allgenes <- names(revmap(getGeneSymbols()))
	genes <- unique(targets$gene)
	geneList <- factor(as.integer(allgenes %in% genes))
	names(geneList) <- genes
	#str(geneList)
	
	GOdata <- new("topGOdata",
		ontology = ontology,
		allGenes = geneList,
		nodeSize = 5,
		annot = annFUN.org, 
		mapping = "org.Hs.eg.db",
		ID = "symbol")
	return(GOdata)
}
#GOdata <- loadGOdata(targets.hcv_vs_nash)


getTopGoTerms <- function(targets, ontology='BP', filename=NULL, show.plot=FALSE)
{
	GOdata <- loadGOdata(targets,ontology)
	resultFisher <- getSigGroups(GOdata, new("classicCount", testStatistic = GOFisherTest, name = "Fisher test"))
	resultKS <- getSigGroups(GOdata, new("classicScore", testStatistic = GOKSTest, name = "KS tests"))
	resultElim <- getSigGroups(GOdata, new("elimCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 0.01))
	resultWeight <- getSigGroups(GOdata, new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio"))
	
	sig.tab <- GenTable(GOdata, P.value = resultFisher, #KS=resultKS, elim=resultElim, weight=resultWeight, 
			orderBy = "classic", ranksOf = "classic", topNodes = 20, numChar=1000)
	
	if (!is.null(filename))
		writeTable(sig.tab, filename, row.names=FALSE)
	
	if (show.plot)
	{
		l <- list(classic = score(resultFisher), KS = score(resultKS), elim = score(resultElim), weight = score(resultWeight))	
		oldPar <- par(mfrow = c(2, 2))
		for (nn in names(l))
		{
			p.val <- l[[nn]]
			hist(p.val[p.val < 1], br = 50, xlab = "p values", main = paste("Histogram for method:", nn, sep = " "))
		}
		par(mfrow = c(1, 1))
	}
	
	return(sig.tab)
}
#getTopGoTerms(targets.hcv_vs_nash)

