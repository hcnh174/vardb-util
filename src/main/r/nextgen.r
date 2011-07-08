library(gsubfn) 


loadSamples <- function(filename='samples.txt')
{
	data <- loadDataframe('samples.txt', stringsAsFactors=FALSE)
	rownames(data) <- data$identifier
	return(data)
}


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

#############################################################################

preprocess <- function(dir='/home/nelson/nextgen2/', path='GA_RunData/110624_HWUSI-EAS1611_00063_FC639J3AAXX/Unaligned',
		temp.dir='tmp', file='out/preprocess.txt', fastq.dir='fastq')
{
	samples <- loadSamples()
	cat('', file=file)
	appendFile('cd ',dir, file=file)
	appendFile('mkdir ',temp.dir, file=file)
	appendFile('mkdir ',fastq.dir, file=file)
	appendFile('rm -r ',temp.dir,'/*', file=file)
	appendFile('rm ',fastq.dir,'/*', file=file)
	#cat('mkdir ',fastq.dir,'\n', sep='', file=file, append=TRUE)
	for (identifier in samples$identifier)
	{
		row <- samples[identifier,]
		dir.to <- concat(dir,temp.dir,'/',identifier)
		project <- row$project
		barcode <- row$barcode
		lane <- row$lane
		dir.from <- concat(dir,path,'/Project_',project,'/Sample_',project,'/')
		filename <- concat(project,'_',barcode,'_L00',lane,'_R1_*.fastq.gz')
		
		appendFile('mkdir ',dir.to, file=file)
		appendFile('cd ',dir.from, file=file)
		appendFile('cp ',filename,' ',dir.to, file=file)
		appendFile('gunzip ',dir.to,'/*', file=file)
		appendFile('cat ',dir.to,'/* > ',dir,fastq.dir,'/',identifier,'.fq', file=file)
		appendFile('cd ',dir, file=file)
		appendFile('', file=file)
	}
}

#########################################################################3

getVarRefName <- function(refid, index)
{
	if (index==1)
		name <- refid
	else name <- paste(refid,'var',(index-1), sep='')
	return(name)
}

# makes variant reference sequences using all possible combinations of specified variants codons
makeVariantsForRef <- function(refid,refs,variants)
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
	return(varseqs)
}

writeVariantsForRef <- function(refid, varseqs, ref.dir) #='out/')
{
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		#if (index==1)
		#	name <- refid
		#else name <- paste(refid,'var',(index-1), sep='')
		varseq <- paste(varseqs[index,], collapse='')
		filename <- concat(ref.dir,'/',name,'.fasta')
		print(concat('writing variant file ',filename))
		cat('>',name,'\n',varseq,'\n', sep='', file=filename)
	}
}
#refs <- loadRefs()
#variants <- loadVariants()
#makeVariantsForRef('NS3aa156',refs,variants)

# determine varref names by counting how many variants are present for a particular sample
getVarRefNames <- function(ref)
{
	refs <- loadRefs()
	variants <- loadVariants()
	varseqs <- makeVariantsForRef(ref,refs,variants)
	refnames <- c()
	for (index in 1:nrow(varseqs))
	{
		refnames <- c(refnames,getVarRefName(ref,index))
	}
	return(refnames)
}

indexReferences <- function(refid, varseqs, file, ref.dir='ref', index.dir='indexes')
{	
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		cat('bowtie-build ',ref.dir,'/',name,'.fasta',' ',index.dir,'/',name,'\n', sep='', file=file, append=TRUE)	
	}
	
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		cat('bwa index ',ref.dir,'/',name,'.fasta','\n', sep='', file=file, append=TRUE)
	}
}

makeVariants <- function(ref.dir='ref', index.dir='indexes', file='out/index_refs.txt')
{
	refs <- loadRefs()
	variants <- loadVariants()
	system(concat('mkdir ',ref.dir))
	system(concat('mkdir ',index.dir))
	cat('', file=file)
	for (refid in rownames(refs))
	{
		varseqs <- makeVariantsForRef(refid, refs, variants)
		writeVariantsForRef(refid, varseqs, ref.dir)
		indexReferences(refid, varseqs, file)
	}
}

##################################################################

loadVariantData <- function(sample, ref, merged=FALSE)#, nums)
{
	if (merged==TRUE)
		filename <- concat('variants/',sample,'.',ref,'.merged.txt')
	else filename <- concat('variants/',sample,'.',ref,'.txt')
	#print(paste('Loading file',filename))
	data <- loadDataframe(filename)
	#print(paste('Loaded',filename))
#	
#	positions <- c()
#	for (num in nums)
#	{
#		position <- (num-1)*3
#		positions <- c(positions, position, position+1, position+2)
#	}
#	
#	data <- subset(data, position %in% positions)
	#data <- subset(data, position >= startpos & position <= endpos)
	return(data)
}


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

###########################################################################

getCodonCodeTable <- function()
{
	codes <- data.frame()
	codes['TTT','aa'] <- 'F'
	codes['TTC','aa'] <- 'F'
	codes['TTA','aa'] <- 'L'
	codes['TTG','aa'] <- 'L'
	codes['CTT','aa'] <- 'L'
	codes['CTC','aa'] <- 'L'
	codes['CTA','aa'] <- 'L'
	codes['CTG','aa'] <- 'L'
	codes['ATT','aa'] <- 'I'
	codes['ATC','aa'] <- 'I'
	codes['ATA','aa'] <- 'I'
	codes['ATG','aa'] <- 'M'
	codes['GTT','aa'] <- 'V'
	codes['GTC','aa'] <- 'V'
	codes['GTA','aa'] <- 'V'
	codes['GTG','aa'] <- 'V'
	codes['TCT','aa'] <- 'S'
	codes['TCC','aa'] <- 'S'
	codes['TCA','aa'] <- 'S'
	codes['TCG','aa'] <- 'S'
	codes['CCT','aa'] <- 'P'
	codes['CCC','aa'] <- 'P'
	codes['CCA','aa'] <- 'P'
	codes['CCG','aa'] <- 'P'
	codes['ACT','aa'] <- 'T'
	codes['ACC','aa'] <- 'T'
	codes['ACA','aa'] <- 'T'
	codes['ACG','aa'] <- 'T'
	codes['GCT','aa'] <- 'A'
	codes['GCC','aa'] <- 'A'
	codes['GCA','aa'] <- 'A'
	codes['GCG','aa'] <- 'A'
	codes['TAT','aa'] <- 'Y'
	codes['TAC','aa'] <- 'Y'
	codes['TAA','aa'] <- 'X'
	codes['TAG','aa'] <- 'X'
	codes['CAT','aa'] <- 'H'
	codes['CAC','aa'] <- 'H'
	codes['CAA','aa'] <- 'Q'
	codes['CAG','aa'] <- 'Q'
	codes['AAT','aa'] <- 'N'
	codes['AAC','aa'] <- 'N'
	codes['AAA','aa'] <- 'K'
	codes['AAG','aa'] <- 'K'
	codes['GAT','aa'] <- 'D'
	codes['GAC','aa'] <- 'D'
	codes['GAA','aa'] <- 'E'
	codes['GAG','aa'] <- 'E'
	codes['TGT','aa'] <- 'C'
	codes['TGC','aa'] <- 'C'
	codes['TGA','aa'] <- 'X'
	codes['TGG','aa'] <- 'W'
	codes['CGT','aa'] <- 'R'
	codes['CGC','aa'] <- 'R'
	codes['CGA','aa'] <- 'R'
	codes['CGG','aa'] <- 'R'
	codes['AGT','aa'] <- 'S'
	codes['AGC','aa'] <- 'S'
	codes['AGA','aa'] <- 'R'
	codes['AGG','aa'] <- 'R'
	codes['GGT','aa'] <- 'G'
	codes['GGC','aa'] <- 'G'
	codes['GGA','aa'] <- 'G'
	codes['GGG','aa'] <- 'G'
	return(codes)
}
codes <- getCodonCodeTable()

translateCodon <- function(codon)
{
	aa <- codes[codon,'aa']
	if (is.na(aa))
		aa <- 'X'
	return(aa)
}
#translateCodon('GGG')


