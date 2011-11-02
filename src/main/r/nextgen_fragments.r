
findFragmentStartPosition <- function(refseq, seq)
{
	#refseq <-  getField(config@refs,refseqid,'sequence')
	#seq <-  getField(config@refs,seqid,'sequence')
	psa1 <- Biostrings::pairwiseAlignment(pattern = refseq, subject = seq, type='local', gapOpening = -1000000)
	print(psa1)
	startnt <- Biostrings::pattern(psa1)@range@start
	#return(Biostrings::pattern(psa1)@range)
	endnt <- startnt + nchar(seq) -1 
	print(concat(startnt,'..',endnt))
	return(list(start=startnt,end=endnt))
}
#findFragmentStartPosition(config,'HCV-HCJ4','HCV-KT9-NS3') #3410..5302
#findFragmentStartPosition(config,'HCV-HCJ4','HCV-KT9-NS5A') #6248..7588

findFragmentStartPositions <- function(refseq, suffix, filename='../config/merged/fragments.fasta', outfile='../config/merged/fragments.txt')
{
	seqs <- readFastaFile(filename)
	starts <- data.frame()
	for (id in names(seqs))
	{
		region <- strsplit(id,'__', fixed=TRUE)[[1]][2]
		pos <- findFragmentStartPosition(refseq, seqs[[id]])
		starts[id,'id'] <- id
		starts[id,'region'] <- region
		starts[id,'start'] <- pos$start
	}
	makeBackupFile(outfile)
	writeTable(starts,outfile,row.names=FALSE)
	return(starts)
}
#refseq <-  getRefSequence(config,'HCV-HCJ4'); starts <- findFragmentStartPositions(refseq,'NS3-36')

findBestFragmentStartPositions <- function(refseq, filename='../config/merged/fragments.fasta', startsfile='../config/merged/fragments.txt')
{
	seqs <- readFastaFile(filename)
	data <- loadDataFrame(startsfile, idcol='id')
	data <- data[order(data$region,data$start),]
	for (region in unique(data$region))
	{
		print(region)
		rows <- data[which(data$region==region),]
		starts <- rows$start
		beststart <- as.numeric(names(sort(table(starts), decreasing=TRUE)[1]))
		samplelen <- 60
		for (id in rownames(rows))
		{
			start <- rows[id,'start']
			print(concat(id,': ',start))
			print(tolower(substring(refseq,start,start+samplelen-1)))
			print(tolower(substring(seqs[[id]],1,samplelen)))
		}
	}
	return(starts)
}
#refseq <-  getRefSequence(config,'HCV-HCJ4'); starts <- findBestFragmentStartPositions(refseq)

getSharedNameForFragment <- function(refid,id)
{
	sharedid <- strsplit(id,'__', fixed=TRUE)[[1]][1]
	return(concat(refid,'_',sharedid))
}
#getSharedNameForFragment('HCV-KT9','CTE247-21__NS3-156')

#mergeFragmentWithReferenceSequence <- function(refseq, seq, sep='')
#{
#	loc <- findFragmentStartPosition(refseq,seq)
#	prefix <- tolower(substring(refseq,1,loc$start-1))
#	suffix <- tolower(substring(refseq,loc$end+1))
#	hybrid <- concat(prefix,sep,toupper(seq),sep,suffix)
#	return(hybrid)
#}

mergeFragmentWithReferenceSequence <- function(refseq, seq, start, sep='')
{
	prefix <- tolower(substring(refseq,1,start-1))
	suffix <- tolower(substring(refseq,start+nchar(seq)))
	hybrid <- concat(prefix,sep,toupper(seq),sep,suffix)
	return(hybrid)
}
#refseq <-  getField(config@refs,'HCV-HCJ4','sequence')
#seq <- getField(config@refs,'HCV-KT9-NS3','sequence')
#mergeFragmentWithReferenceSequence(refseq,seq)
#mergeFragmentWithReferenceSequence('aabbbcccc','BBB',3)

createHybridReferenceSequences <- function(refid, filename='../config/merged/fragments.fasta', startsfile='../config/merged/fragments.txt',
		outfile='../config/merged/hybrid.fasta')
{
	refseq <-  getRefSequence(config,refid)
	seqs <- readFastaFile(filename)
	data <- list()
	# extract the shared id for each fragment (part before the separator)
	for (id in names(seqs))
	{
		sharedid <- getSharedNameForFragment(refid,id)
		#set the initial ref sequence for each one - overwrites if all ready present
		data[[sharedid]] <- refseq
	}
	startdata <- loadDataFrame(startsfile, idcol='id')
	for (id in names(seqs))
	{
		seq <- seqs[[id]]
		sharedid <- getSharedNameForFragment(refid,id)
		refseq <- data[[sharedid]]
		start <- startdata[id,'start']
		# if there are multiple fragments sharing the same prefix, then they will be updated sequentially
		data[[sharedid]] <- mergeFragmentWithReferenceSequence(refseq,seq,start)
	}
	writeFastaFile(outfile,data)
	return(data)
}
#createHybridReferenceSequences('HCV-HCJ4')
#
#
#createHybridReferenceSequences <- function(refid, filename='../config/merged/fragments.fasta', outfile='../config/merged/hybrid.fasta')
#{
#	refseq <-  getRefSequence(config,refid)
#	seqs <- readFastaFile(filename)
#	data <- list()
#	# extract the shared id for each fragment (part before the separator)
#	for (id in names(seqs))
#	{
#		sharedid <- getSharedNameForFragment(refid,id)
#		#set the initial ref sequence for each one - overwrites if all ready present
#		data[[sharedid]] <- refseq
#	}
#	for (id in names(seqs))
#	{
#		seq <- seqs[[id]]
#		sharedid <- getSharedNameForFragment(refid,id)
#		refseq <- data[[sharedid]]
#		# if there are multiple fragments sharing the same prefix, then they will be updated sequentially
#		data[[sharedid]] <- mergeFragmentWithReferenceSequence(refseq,seq)
#	}
#	writeFastaFile(outfile,data)
#	return(data)
#}
##createHybridReferenceSequences('HCV-HCJ4')
