
findFragmentStartPosition <- function(refseq, seq)
{
	#refseq <-  getField(config@refs,refseqid,'sequence')
	#seq <-  getField(config@refs,seqid,'sequence')
	psa1 <- Biostrings::pairwiseAlignment(pattern = refseq, subject = seq, type='local', gapOpening = -1000000)
	print(psa1)
	startnt <- Biostrings::pattern(psa1)@range@start
	#return(Biostrings::pattern(psa1)@range)
	endnt <- startnt + nchar(seq) -1 
	printcat(startnt,'..',endnt)
	return(list(start=startnt,end=endnt))
}
#findFragmentStartPosition(config,'HCV-HCJ4','HCV-KT9-NS3') #3410..5302
#findFragmentStartPosition(config,'HCV-HCJ4','HCV-KT9-NS5A') #6248..7588

findFragmentStartPositions <- function(config, refid, region, filename=concat(config@config.dir,'/fragments-',region,'.fasta'), outfile=concat(config@config.dir,'/tmp/fragments-',refid,'-',region,'.txt'))
{
	refseq <-  getRefSequence(config,refid)
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
#findFragmentStartPositions(config,'HCV-HCJ4','NS5A-93')

findBestFragmentStartPositions <- function(config, refid, region, 
		filename=concat(config@config.dir,'/fragments-',region,'.fasta'), 
		startsfile=concat(config@config.dir,'/tmp/fragments-',refid,'-',region,'.txt'),
		cutoff=0.4)
{
	refseq <-  getRefSequence(config,refid)
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
			printcat(id,': ',start)
			seq1 <- tolower(substring(refseq,start,start+samplelen-1))
			seq2 <- tolower(substring(seqs[[id]],1,samplelen))
			mismatches <- MiscPsycho::stringMatch(seq1, seq2, normalize = 'no')
			freq <- mismatches/nchar(seq2)
			if (freq > cutoff)
			{
				#aln <- Biostrings::pairwiseAlignment(seq1,seq2,type='local')#gapOpening = -1000000)#,)#,
				#offset <- aln@subject@range@start
				#newstart <- start+offset
				#data[id,'oldstart'] <- start
				#data[id,'start'] <- newstart
				print('')
				#print(aln)
				print(seq1)
				print(seq2)
				printcat(freq,' [',mismatches,'/',nchar(seq2),']')
				#printcat('offset=',offset,', oldstart=',start,', newstart=',newstart))
				print('')
				print('###########################################################')
				print('')
			}
		}
	}
	#writeTable(data,startsfile,row.names=FALSE)
	return(starts)
}
#findBestFragmentStartPositions(config,'HCV-KT9','NS3-36')

getSharedNameForFragment <- function(refid,id)
{
	sharedid <- strsplit(id,'__', fixed=TRUE)[[1]][1]
	return(concat(refid,'_',sharedid))
}
#getSharedNameForFragment('HCV-KT9','CTE247-21__NS3-156')

mergeFragmentWithReferenceSequence <- function(refseq, seq, start, sep='')
{
	prefix <- substring(refseq,1,start-1)
	suffix <- substring(refseq,start+nchar(seq))
	hybrid <- concat(prefix,sep,seq,sep,suffix)
	return(hybrid)
}
#refseq <-  getField(config@refs,'HCV-HCJ4','sequence')
#seq <- getField(config@refs,'HCV-KT9-NS3','sequence')
#mergeFragmentWithReferenceSequence(refseq,seq)

createHybridReferenceSequences <- function(config, refid,
	fragments.dir=concat(config@config.dir,'/fragments'),
	startsfile=concat(fragments.dir,'/fragments-',refid,'.txt'), #concat(config@config.dir,'/fragments/fragments-',refid,'.txt'),
	outfile=concat(config@config.dir,'/tmp/hybrid-',refid,'.fasta'))
{
	origrefseq <-  tolower(getRefSequence(config,refid))
	seqs <- readFastaFiles(fragments.dir,'^fragments.*\\.fasta')
	data <- list()
	# extract the shared id for each fragment (part before the separator)
	#set the initial ref sequence for each one - overwrites if all ready present
	for (id in names(seqs))
	{
		sharedid <- getSharedNameForFragment(refid,id)
		print(sharedid)
		data[[sharedid]] <- origrefseq
	}
	startdata <- loadDataFrame(startsfile, idcol='id')
	print(head(startdata))
	# if there are multiple fragments sharing the same prefix, then they will be updated sequentially
	for (id in names(seqs))
	{
		seq <- toupper(seqs[[id]])
		sharedid <- getSharedNameForFragment(refid,id)
		refseq <- data[[sharedid]]
		start <- startdata[id,'start']
		printcat('id=',id,', sharedid=',sharedid)
		hybrid <- mergeFragmentWithReferenceSequence(refseq,seq,start)
		data[[sharedid]] <- hybrid
	}
	writeFastaFile(outfile,data)
	return(data)
}
#data <- createHybridReferenceSequences(config, refid='HCV-HCJ4')
#writeRefs(config)