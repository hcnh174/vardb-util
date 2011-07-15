library(gsubfn)
library(seqinr)
library(methods)

setClass("variantdata", representation(nt="data.frame", codons="data.frame", aa="data.frame"))

setClass("sampleparams",
		representation(subject='character',
				sample='character',
				week='numeric',
				region='character',
				drop.ambig='logical',
				nt.cutoff='numeric',
				out.dir='character'),
		prototype(drop.ambig=TRUE,
				nt.cutoff=0,
				out.dir='out/'))

counts.dir <- 'counts/' 
config.dir <- 'config/'
variants.dir <- 'variants/'

loadRuns <- function(filename=concat(config.dir,'runs.txt'))
{
	data <- loadDataFrame(filename, stringsAsFactors=FALSE)
	rownames(data) <- data$run
	return(data)
}

loadSamples <- function(filename=concat(config.dir,'samples.txt'))
{
	data <- loadDataFrame(filename, stringsAsFactors=FALSE)
	rownames(data) <- data$sample
	return(data)
}

loadRefs <- function(filename=concat(config.dir,'refs.fasta'))
{
	data <- data.frame()
	sequences <- read.fasta(file = filename, as.string = TRUE, seqtype = "DNA", forceDNAtolower=TRUE)
	for (id in names(sequences))
	{
		data[id,'sequence'] <- sequences[[id]][1]
	}
	return(data)
}

loadRegions <- function(filename=concat(config.dir,'regions.txt'))
{
	data <- loadDataFrame(filename, stringsAsFactors=FALSE)
	rownames(data) <- data$region
	return(data)
}

loadVariants <- function(filename=concat(config.dir,'variants.txt'))
{
	data <- loadDataFrame(filename, stringsAsFactors=FALSE)
	return(data)
}

loadTiters <- function(filename=concat(config.dir,'titers.txt'))
{
	data <- loadDataFrame(filename)
	return(data)
}

loadTreatments <- function(filename=concat(config.dir,'treatments.txt'))
{
	data <- loadDataFrame(filename)
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

# copy all files from the same sample (subject + date) to the same folder
preprocess2 <- function(dir='/home/nelson/nextgen2/', path='GA_RunData/110624_HWUSI-EAS1611_00063_FC639J3AAXX/Unaligned',
		temp.dir='tmp', file='out/preprocess.txt', fastq.dir='fastq')
{
	runs <- loadRuns()
	createFile(file=file)
	appendFile('cd ',dir, file=file)
	appendFile('mkdir ',temp.dir, file=file)
	appendFile('mkdir ',fastq.dir, file=file)
	appendFile('rm -r ',temp.dir,'/*', file=file)
	appendFile('rm ',fastq.dir,'/*', file=file)
	
	for (sample in unique(runs$sample))
	{
		dir.to <- concat(temp.dir,'/',sample)
		appendFile('mkdir ',dir.to, file=file)
	}
	appendFile('', file=file)
	
	for (run in rownames(runs))
	{
		row <- runs[run,]
		sample <- row$sample
		dir.to <- concat(temp.dir,'/',sample)
		project <- row$project
		barcode <- row$barcode
		lane <- row$lane		
		dir.from <- concat(path,'/Project_',project,'/Sample_',project,'/')
		filename <- concat(project,'_',barcode,'_L00',lane,'_R1_*.fastq.gz')
		
		appendFile('cp ',dir.from,filename,' ',dir.to, file=file)
	}
	appendFile('', file=file)
	
	for (sample in unique(runs$sample))
	{
		dir.from <- concat(temp.dir,'/',sample)
		appendFile('gunzip ',dir.from,'/*', file=file)
	}
	appendFile('', file=file)
	
	for (sample in unique(runs$sample))
	{
		dir.from <- concat(temp.dir,'/',sample)
		appendFile('cat ',dir.from,'/* > ',fastq.dir,'/',sample,'.fastq', file=file)
		
	}
	appendFile('', file=file)
	#appendFile('rm -r ',temp.dir,'/*', file=file)
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
#makeVariantsForRef('KT9',refs,variants)

writeVariantsForRef <- function(refid, varseqs, ref.dir) #='out/')
{
	for (index in 1:nrow(varseqs))
	{
		name <- getVarRefName(refid,index)
		varseq <- paste(varseqs[index,], collapse='')
		filename <- concat(ref.dir,'/',name,'.fasta')
		print(concat('writing variant file ',filename))
		write.fasta(s2c(varseq), name, file.out=filename)
		#write.fasta(strsplit(varseq,"")[[1]], name, file.out=filename)
	}
}

# determine varref names by counting how many variants are present for a particular sample
getVarRefNames <- function(ref)
{
	refs <- loadRefs()
	variants <- loadVariants()
	varseqs <- makeVariantsForRef(ref,refs,variants)
	if (nrow(varseqs)==0)
		return(ref)
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
	#system(concat('mkdir ',index.dir))
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
	data <- loadDataFrame(filename)
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

translateCodon <- function(codon)
{
	return(c2s(translate(s2c(codon))))
}
#translateCodon('GGG')

extractSequence <- function(sequence, start, end)
{
	sequence <- s2c(sequence)
	return(c2s(sequence[start:end]))
}
#extractSequence(refs['KT9','sequence'],3420,5312)


plotReadDistributions <- function(filename="histograms.pdf")
{
	pdf(filename)
	#par(mfrow=c(2,2))
	#par(ask=TRUE)
	for (smpl in rownames(samples))
	{
		try({
			ref <- samples[smpl,'ref']
			filename <- concat('variants/',smpl,'.',ref,'.txt')
			data <- loadDataFrame(filename)
			table <- createNtCountTable(data, cutoff=0)
			runs.subset <- subset(runs, sample==smpl)
			#par(mfrow=c(1,nrow(runs.subset)))
			par(mfrow=c(2,nrow(runs.subset)))
			for (run in rownames(runs.subset))
			{
				region <- runs.subset[run,'region']
				start <- regions[region,'start']
				end <- regions[region,'end']
				xmin <- min(start)-10
				xmax <- max(end) + 10
				ymax <- max(table$total)+10
				xlab <- concat('position (',xmin,'-',xmax,')')
				#main <- concat(smpl,': ',run,': ',region,': ',start,':',end)
				plot(top1 ~ as.numeric(position), table, ylim=c(0,ymax), xlim=c(xmin,xmax), type='h',
						xlab=xlab, ylab='read coverage', main=run, sub=region)
				plot(top1 ~ as.numeric(position), table, ylim=c(0,500), xlim=c(xmin,xmax), type='h',
						xlab=xlab, ylab='read coverage', main=run, sub=region)
			}
			#par(mfrow=c(1,1))
		}, silent=FALSE)
	}
	#par(ask=FALSE)
	par(mfrow=c(1,1))
	dev.off()
}

####################################################################

getCodonPositionsForRegion <- function(region)
{
	offset <- 0
	parent <- regions[region,'parent']
	if (!is.na(parent))
		offset <- regions[parent,'start']
	start <- regions[region,'start'] - offset
	end <- regions[region,'end'] - offset
	#print(concat(start,':',end))
	sequence <- extractSequence(refs['KT9','sequence'], regions[region,'start'], regions[region,'end'])
	#print(sequence)
	#print(translateCodon(sequence))
	positions <- data.frame()
	codon <- start/3 + 1
	if (start %% 3 !=0)
		stop(concat('start number is not a multiple of 3: ',start,' in region ',region))
	for (position in seq(start,end,3))
	{
		refcodon <- extractSequence(refs['KT9','sequence'], offset+position, offset+position+2)
		refaa <- translateCodon(refcodon)
		positions <- rbind(positions, data.frame(codon=codon, ntnum=offset+position, relpos=position, refcodon=refcodon, refaa=refaa))
		codon <- codon + 1
	}	
	return(positions)
}
#head(getCodonPositionsForRegion('NS3aa156'))
#head(getCodonPositionsForRegion('NS3aa36'))
#head(getCodonPositionsForRegion('NS5Aaa31'))
#head(getCodonPositionsForRegion('NS5Aaa93'))


extractCodonData <- function(data, ntnum, drop.ambig=FALSE)
{
	#hack!!!
	ntnum <- ntnum - 1
	nt1 <- data[which(data$position==ntnum),'nt']
	nt2 <- data[which(data$position==ntnum+1),'nt']
	nt3 <- data[which(data$position==ntnum+2),'nt']
	codons <- paste(nt1, nt2, nt3, sep='')		
	if (drop.ambig)
		codons <- removeAmbiguousCodons(codons)
	return(codons)
}
#extractCodonData(data,3495)

################################################################################

appendSampleParams <- function(counts, params)
{
	counts$subject <- as.character(params@subject)
	counts$sample <- as.character(params@sample)
	counts$week <- as.character(params@week)
	counts$region <- as.character(params@region)
	cols <- names(counts)
	cols <- c('subject','sample','week','region',cols[1:(length(cols)-4)])
	counts <- counts[,cols]
	return(counts)
}
#
#createNtCountTableAlt <- function(data, params)
#{
#	region <- params@region
#	start <- regions[region,'start']
#	end <- regions[region,'end']	
#	data <- subset(data, position >= start & position <= end)
#	
#	data$nt <- factor(data$nt, levels=c('A','C','G','T','N'))
#	table <- xtabs(~position + nt, data)
#	
#	for (nt in c('A','C','G','T','N'))
#	{
#		table[,nt] <- ifelse(table[,nt] > params@nt.cutoff, table[,nt], 0)
#	}
#	position <- as.numeric(rownames(table))
#	total <- apply(table, 1, sum)
#	
#	top1 <- apply(table, 1, function(values){sort(values, decreasing=TRUE)[[1]]})
#	top2 <- apply(table, 1, function(values){sort(values, decreasing=TRUE)[[2]]})
#	top3 <- apply(table, 1, function(values){sort(values, decreasing=TRUE)[[3]]})
#	top4 <- apply(table, 1, function(values){sort(values, decreasing=TRUE)[[4]]})
#	top5 <- apply(table, 1, function(values){sort(values, decreasing=TRUE)[[5]]})
#	
#	table <- cbind(table,position)
#	table <- cbind(table,top1)
#	table <- cbind(table,top2)
#	table <- cbind(table,top3)
#	table <- cbind(table,top4)
#	table <- cbind(table,top5)
#	table <- cbind(table,total)
#	table <- data.frame(table)
#	
#	positions <- getCodonPositionsForRegion(region)
#	table$codon <- sapply(table$position, function(ntnum){
#				return(positions[which(positions$ntnum <= ntnum & positions$ntnum+2 >= ntnum),'codon'])
#			})
#	table <- appendSampleParams(table, params)
#	return(table)
#}
##table.nt <- createNtCountTableAlt(data, 'NS3aa36')

createNtCountTable <- function(data, params)
{
	region <- params@region
	start <- regions[region,'start']
	end <- regions[region,'end']	
	data <- subset(data, position >= start & position <= end)
	
	counts <- data.frame()
	positions <- getCodonPositionsForRegion(params@region)
	for (ntnum in positions$ntnum)
	{
		try({
			aanum <- positions[which(positions$ntnum==ntnum),'codon']
			for (offset in 0:2)
			{
				ntnum2 <- ntnum + offset
				nt <- data[which(data$position==ntnum2),'nt']
				#print(head(nt))
				#if (params@drop.ambig)
				#	nt <- removeAmbiguousCodons(nt)
				#print(concat('nt vector length: ',length(nt)))
				freqs <- sort(xtabs(as.data.frame(nt)), decreasing=TRUE)
				total <- sum(freqs)				
				rank <- 1
				for (base in names(freqs))
				{
					count <- freqs[base]
					freq <- count/total
					row <- data.frame(ntnum=ntnum, aanum=aanum, nt=base, rank=rank, count=count, freq=freq) 
					counts <- rbind(counts,row)
					rank <- rank+1
				}
			}
		}, silent=FALSE)
	}
	counts$nt <- as.character(counts$nt)
	counts$ntnum <- factor(counts$ntnum)
	counts$aanum <- factor(counts$aanum)
	counts <- appendSampleParams(counts, params)
	return(counts)
}

createCodonCountTable <- function(data, params)
{
	counts <- data.frame()
	positions <- getCodonPositionsForRegion(params@region)
	for (ntnum in positions$ntnum)
	{
		codons <- extractCodonData(data,ntnum,params@drop.ambig)
		aanum <- positions[which(positions$ntnum==ntnum),'codon']
		#print(paste('aanum=',aanum,'ntnum=',ntnum,'numcodons=',length(codons)))
		freqs <- sort(xtabs(as.data.frame(codons)), decreasing=TRUE)
		total <- sum(freqs)
		
		rank <- 1
		for (codon in names(freqs))
		{
			count <- freqs[codon]
			freq <- count/total
			row <- data.frame(ntnum=ntnum, aanum=aanum, codon=codon, rank=rank, count=count, freq=freq) #sample=sample, 
			counts <- rbind(counts,row)
			rank <- rank +1
		}
	}
	counts$codon <- as.character(counts$codon)
	counts$aa <- sapply(counts$codon,translateCodon)
	counts$ntnum <- factor(counts$ntnum)
	counts$aanum <- factor(counts$aanum)
	counts <- appendSampleParams(counts, params)
	return(counts)
}
#data <- loadDataFrame('variants/218-7_03-01.KT9.txt')
#params <- new('sampleparams',subject='218-7', sample='218-7_03-01', region='NS3aa36')
#table.codons <- createCodonCountTable(data,params)
#head(table.codons, n=20)

createAminoAcidCountTable <- function(data, params)
{
	counts <- data.frame()
	positions <- getCodonPositionsForRegion(params@region)
	for (ntnum in positions$ntnum)
	{		
		codons <- extractCodonData(data,ntnum,params@drop.ambig)
		aanum <- positions[which(positions$ntnum==ntnum),'codon']
		aa <- sapply(codons,translateCodon)
		freqs <- sort(xtabs(as.data.frame(aa)), decreasing=TRUE)
		total <- sum(freqs)
		rank <- 1
		for (aa in names(freqs))
		{
			count <- freqs[aa]
			#print(attributes(count))
			freq <- count/total
			row <- data.frame(ntnum=ntnum, aanum=aanum, aa=aa, rank=rank, count=count, freq=freq)
			counts <- rbind(counts,row)
			rank <- rank +1
		}
	}
	counts$ntnum <- factor(counts$ntnum)
	counts$aanum <- factor(counts$aanum)
	counts <- appendSampleParams(counts, params)
	return(counts)
}
#data <- loadDataFrame('variants/218-7_03-01.KT9.txt')
#params <- new('sampleparams',subject='218-7', sample='218-7_03-01', region='NS3aa36')
#table.aa <- createAminoAcidCountTable(data,params)


count_codons_for_region <- function(data, params, variantdata)
{
	print(concat('count_codons_for_region: ',as.character(params@region)))
	variantdata@nt <- rbind(variantdata@nt, createNtCountTable(data, params))
	variantdata@codons <- rbind(variantdata@codons, createCodonCountTable(data, params))
	variantdata@aa <- rbind(variantdata@aa, createAminoAcidCountTable(data, params))	
	return(variantdata)
}

count_codons_for_sample <- function(params, variantdata)
{
	sample <- as.character(params@sample)
	print(concat('count_codons_for_sample: ',sample))
	ref <- samples[sample,'ref'] # look up the ref for the sample
	filename <- concat(variants.dir,sample,'.',ref,'.txt'); print(filename) # load the corresponding data file
	data <- loadDataFrame(filename)
	print(concat('loaded file ',filename,'. contains ',nrow(data),' reads'))
	# each sample has several runs targeting different regions
	for (region in runs[runs[,'sample']==sample,'region'])
	{
		params@region <- region
		variantdata <- count_codons_for_region(data, params, variantdata)
	}
	return(variantdata)
}
#tables <- count_codons_for_sample('218-7_03-01')

count_codons_for_subject <- function(params)
{
	if (is.character(params))
		params <- new('sampleparams', subject=params)
	print(concat('count_codons_for_subject: ',as.character(params@subject)))
	variantdata <- new('variantdata')
	for (sample in samples[which(samples$subject == params@subject),'sample'])
	{
		params@sample <- sample
		params@week <- samples[sample,'week']
		try({variantdata <- count_codons_for_sample(params, variantdata)}, silent=FALSE)
	}
	out.dir <- params@out.dir
	#out.dir <- counts.dir
	writeTable(variantdata@nt, concat(out.dir,params@subject,'.nt.txt'), row.names=FALSE)
	writeTable(variantdata@codons, concat(out.dir,params@subject,'.codons.txt'), row.names=FALSE)
	writeTable(variantdata@aa, concat(out.dir,params@subject,'.aa.txt'), row.names=FALSE)
	return(variantdata)
}
#variantdata <- count_codons_for_subject('218-7')

count_codons <- function(params=NULL)
{
	print('count_codons')
	if (is.null(params))
		params <- new('sampleparams')	
	for (subject in unique(samples$subject))
	{
		params@subject <- subject
		count_codons_for_subject(params)
	}
}
#count_codons()

# estimates the number of reads based on the size of the file, assuming a ratio of 7757 for uncompressed fastq
estimateReadCount <- function(mb)
{
	return(round(mb*7757))
}
#estimateReadCount(112.5)

