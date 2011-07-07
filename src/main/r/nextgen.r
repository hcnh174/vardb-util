library(gsubfn) 

loadRefs <- function(filename='refs.txt')
{
	data <- loadDataframe(filename, stringsAsFactors=FALSE)
	rownames(data) <- data$id
	return(data)
}

loadVariants <- function(filename='variants.txt')
{
	data <- loadDataframe(filename, stringsAsFactors=FALSE)
	return(data)
}

# makes variant reference sequences using all possible combinations of specified variants codons
makeVariantsForRef <- function(refid,refs,variants, out.dir='out/')
{
	# extract reference sequence
	sequence <- tolower(refs[refid,'sequence'])
	aastart <- as.numeric(refs[refid,'aastart'])
	# divide into codons
	codons <- strapply(sequence, "...")[[1]]
	sets <- list()
	for (no in 1:length(codons))
	{
		aa <- no + aastart - 1
		sets[[paste('aa',aa,sep='')]] <- codons[no]
	}
	#print(sets)	
	
	# get variants for the current reference sequence
	ref.variants <- subset(variants,ref==refid)
	#set the codon position as the rowname
	rownames(ref.variants) <- ref.variants$codon
	for (codon in rownames(ref.variants))
	{	
		aa <- paste('aa',codon,sep='')
		sets[[aa]] <- appendUniqueValues(sets[[aa]], tolower(ref.variants[codon,'variants']))
	}
	# expand the grid to try every combination
	varseqs <- expand.grid(sets)
	# convert the grid to strings instead of factors so can be concatenated
	varseqs <- data.frame(lapply(varseqs, as.character), stringsAsFactors=FALSE)
	
	for (index in 1:nrow(varseqs))
	{
		if (index==1)
			name <- refid
		else name <- paste(refid,'var',(index-1), sep='')
		varseq <- paste(varseqs[index,], collapse='')
		filename <- paste(out.dir,name,'.fasta', sep='')
		print(paste('writing variant file',filename))
		cat('>',name,'\n',varseq,'\n', sep='', file=filename)
	}
}
#refs <- loadRefs()
#variants <- loadVariants()
#makeVariantsForRef('NS3aa156',refs,variants)

makeVariants <- function(refs,variants, out.dir='out/')
{
	for (refid in rownames(refs))
	{
		makeVariantsForRef(refid,refs,variants)
	}
}
