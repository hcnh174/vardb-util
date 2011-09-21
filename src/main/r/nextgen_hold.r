loadVariantData <- function(sample, ref, merged=FALSE)#, nums)
{
	if (merged==TRUE)
		filename <- concat('variants/',sample,'.',ref,'.merged.txt')
	else filename <- concat('variants/',sample,'.',ref,'.txt')
	#print(paste('Loading file',filename))
	data <- loadDataFrame(filename)
	return(data)
}

getVarRefName <- function(refid, index)
{
	if (index==1)
		name <- refid
	else name <- paste(refid,'var',(index-1), sep='')
	return(name)
}

# makes variant reference sequences using all possible combinations of specified variants codons
makeVariantsForRef <- function(config,refid)
{
	#print(concat('makeVariantsForRef.refid=',refid))
	refs <- config@refs
	variants <- config@variants
	# get variants for the current reference sequence
	ref.variants <- subset(variants,ref==refid)
	if (nrow(ref.variants)==0)
		return(data.frame())
	# extract reference sequence
	sequence <- tolower(refs[refid,'sequence'])
	#aastart <- as.numeric(refs[refid,'aastart'])
	codons <- strapply(sequence, "...")[[1]] # split into triplets
	sets <- list()
	positions <- getCodonPositionsForRegion(config,refid)
	aastart <- positions[1,'codon']
	for (no in 1:length(codons))
	{
		aa <- no + aastart - 1
		sets[[paste('aa',aa,sep='')]] <- codons[no]
	}
	#print(sets)
	
	#set the codon position as the rowname
	rownames(ref.variants) <- ref.variants$codon
	for (codon in rownames(ref.variants))
	{	
		aa <- paste('aa',codon,sep='')
		sets[[aa]] <- appendUniqueValues(sets[[aa]], tolower(ref.variants[codon,'variants']))
	}
	#print(sets)	
	# expand the grid to try every combination
	varseqs <- expand.grid(sets)
	# convert the grid to strings instead of factors so can be concatenated
	varseqs <- data.frame(lapply(varseqs, as.character), stringsAsFactors=FALSE)
	return(varseqs)
}
#makeVariantsForRef(config,'NS3aa156')
#makeVariantsForRef(config,'KT9')

writeVariantsForRef <- function(config, refid, varseqs)
{
	ref.dir <- config@ref.dir
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		varseq <- paste(varseqs[index,], collapse='')
		filename <- concat(ref.dir,name,'.fasta')
		print(concat('writing variant file ',filename))
		write.fasta(s2c(varseq), name, file.out=filename)
	}
}

# determine varref names by counting how many variants are present for a particular sample
getVarRefNames <- function(config, ref)
{
	#print(concat('getVarRefNames.ref=',ref))
	varseqs <- makeVariantsForRef(config,ref)
	if (nrow(varseqs)==0)
		return(ref)
	refnames <- c()
	for (index in 1:nrow(varseqs))
	{
		refnames <- c(refnames,getVarRefName(ref,index))
	}
	return(refnames)
}
#getVarRefNames(config, 'KT9')

indexReferences <- function(config, refid, varseqs, file)#, ref.dir='ref', index.dir='indexes')
{	
	ref.dir <- config@ref.dir
	index.dir <- config@index.dir
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		appendFile('bowtie-build ',ref.dir,name,'.fasta',' ',index.dir,name, file=file)
	}
	
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		appendFile('bwa index ',ref.dir,name,'.fasta', file=file)
	}
}

makeVariants <- function(config, file='index_refs.txt')
{
	file <- concat(config@out.dir, file)
	system(concat('mkdir ',config@ref.dir))
	cat('', file=file)
	for (refid in rownames(refs))
	{
		varseqs <- makeVariantsForRef(config, refid)
		writeVariantsForRef(config, refid, varseqs)
		indexReferences(config, refid, varseqs, file)
	}
}

##################################################################

removeAmbiguousCodons <- function(codons)
{
	newvalues <- c()
	for (codon in codons)
	{
		if (regexpr("N", codon) == -1)
			newvalues <- appendValues(newvalues,codon)
	}
	return(newvalues)
}
#removeAmbiguousCodons(splitFields('GCT,GTT,NGT'))

##########################################################


getVariantCounts <- function(config, subject, replicate, ranges=NULL)
{
	#print(concat(subject,'.',replicate))
	filename <- concat(config@counts.dir,subject,'.nt.txt')
	#print(filename)
	data <- loadDataFrame(filename)
	data.subset <- data[which(data$replicate==replicate),]
	if (!is.null(ranges))
	{
		ntnums <- parseRanges(ranges)
		data.subset <- data.subset[which(data.subset$ntnum %in% ntnums),]
	}
	if (nrow(data.subset)==0)
		throw(concat('cannot find data any rows for sample ',subject,'.',replicate))
	counts <- cast(data.subset, ntnum ~ rank, value='count', fun.aggregate=function(x) return(x[1]))#; counts
	counts <- replaceNAs(counts, replacestr=0)
	if (ncol(counts)<4)
	{
		print(counts)
		throw(concat('not enough cols for sample ',subject,'.',replicate,': ',ncol(counts)))
	}
	counts$variants <- apply(counts[,3:ncol(counts)], 1, sum)
	counts$total <- apply(counts[,2:ncol(counts)], 1, sum)
	counts$freq <- counts$variants/counts$total
	return(counts)
}
#counts <- getVariantCounts(config,'PXB0218-0007','wk15','6300-6800')
#counts <- getVariantCounts(config,'PXB0218-0007','wk14','3490-4100'); head(counts)


#extractRefSequence <- function(config, ref, start, end)
#{
#	return(extractSequence(config@refs[ref,'sequence'],start,end))
#}
#extractRefSequence(config,'KT9',3420,5312)

