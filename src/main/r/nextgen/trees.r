library(ape)

#http://the-praise-of-insects.blogspot.com/2010/01/transitions-and-transversions-in-r.html
#http://the-praise-of-insects.blogspot.com/2010/04/transitions-in-r-redux.html

titv <- function(dat)
{
	mat <- as.matrix(dat)
	res <- matrix(NA, ncol=dim(mat)[1], nrow=dim(mat)[1], dimnames=list(x=names(dat), y=names(dat)))
	for (i in 1:(dim(mat)[1] - 1))
	{
		for (j in (i+1):dim(mat)[1])
		{
			vec <- as.numeric(mat[i,])+as.numeric(mat[j,])-8
			res[j,i] <- sum(!is.na(match(vec,c(200,56))))#Transitions
			res[i,j] <- sum(!is.na(match(vec,c(152,168,88,104))))#Transversions
		}
	}
	res
}

#filename <- 'out/consensus/consensus-PXB0219-0018.6321-6810.fasta'
#seqs <- read.dna(filename, format = "fasta")
#ti <- titv(seqs)
#tv[lower.tri(tv)] #Number of transversions
#ti[lower.tri(ti)] #Number of transitions

mergeCodonCountsForSample <- function(config, stem)
{
	data <- NULL
	for (id in rownames(config@data[which(config@data$stem==stem),]))
	{
		filename <- concat(config@counts.dir,'/',id,'.codons.txt')
		#print(filename)
		data.subset <- loadDataFrame(filename)
		#print(head(data.subset))
		if (is.null(data))
			data <- data.subset
		else data <- rbind(data,data.subset)
	}
	#filename <- concat(config@counts.dir,'/',stem,'.codons.txt')
	#writeTable(data,filename,row.names=FALSE)
	return(data)
}
#data <- mergeCodonCountsForSample(config,'PXB0219-0018.wk09')


getVariantSequencesForSample <- function(config, stem, start, end, minfreq=0.001, mincount=50)
{
	data <- mergeCodonCountsForSample(config,stem)
	data <- data[which(data$freq>=minfreq & data$count>=mincount),]
	data <- subset(data, rank<=2)
	maxranks <- max(data$rank)
	varseqs <- list()
	for (rank in 1:maxranks)
	{
		codons <- c()
		for (ntnum in seq(start,end,3))
		{
			refcodon <- unique(data[which(data$position==ntnum & data$rank==1),'codon'])[1]
			if (is.na(refcodon))
				refcodon <- '---'
			codon <- unique(data[which(data$position==ntnum & data$rank==rank),'codon'])[1]
			if (is.na(codon))
				codon <- refcodon
			codons <- c(codons,codon)
		}
		#print(length(codons))
		replicate <- strsplit(stem,'.', fixed=TRUE)[[1]][2]
		name <- concat(replicate,'-',rank)
		varseqs[[name]] <- codons
	}
	#outfile <- concat(config@tmp.dir,'/variants-',stem,'.fasta')
	#writeFastaFile(outfile,varseqs)
	return(varseqs)
}
#seqs <- getVariantSequencesForSample(config, stem='PXB0220-0030.wk29', start=6321, end=6810)
#seqs <- getVariantSequencesForSample(config, stem='PXB0219-0018.wk09', start=6321, end=6810)

#assumes a list of arrays of codon triplets: list(seq1=c('GCT','CGT','---','AGT'), seq2=c()), etc.
removeGappedColumns <- function(seqs)
{
	skipcols <- c()
	for (col in 1:length(seqs[[1]]))#1:((end-start+3)/3))
	{
		for (seqname in names(seqs))
		{
			codon <- seqs[[seqname]][col]
			if (codon=='---')
			{
				skipcols <- c(skipcols, col)
				break
			}
		}
	}
	#print(skipcols)
	newseqs <- list()
	for (seqname in names(seqs))
	{
		newseqs[[seqname]] <- seqs[[seqname]][-skipcols]
	}
	return(newseqs)
}
#seqs.gap <- removeGappedColumns(seqs)

makeTreeForSubject <- function(config, subject, start, end)
{
	seqs <- list()
	samples <- unique(config@data[which(config@data$subject==subject), 'sample'])
	for (sample in samples)
	{
		stem <- getStemForSample(sample)
		seqs <- append(seqs,getVariantSequencesForSample(config, stem, start, end))
	}
	seqs <- removeGappedColumns(seqs)
	seqs <- lapply(seqs,function(seq){paste(seq,collapse='')})
	filename <- concat(config@tmp.dir,'/variants-',subject,'.fasta')
	writeFastaFile(filename,seqs)
	#filename <- makeConsensusFasta(config, samples, name=subject, start=pair$start, end=pair$end)
	seqs <- read.dna(filename, format = "fasta")
	dist <- dist.dna(seqs)
	tree <- nj(dist)
	plot(tree)
	#write.tree(tree, file = 'out/tree.newick', append = FALSE, digits = 10, tree.names = FALSE)
	ti <- titv(seqs)#woodmouse
	tv <- t(ti)
	#print('transversions')
	#print(tv[lower.tri(tv)]) #Number of transversions
	#print('transitions')
	#print(ti[lower.tri(ti)]) #Number of transitions
}
#makeTreeForSubject(config, 'PXB0220-0030', start=6321, end=6810)

