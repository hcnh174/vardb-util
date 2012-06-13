preprocess <- function(config, samples=config@samples, subdirs='ref,fastq,tmp,bam,unmapped,vcf,basecounts,qc,counts,tables,charts,consensus,coverage')#pileup,tmp
{
	makeSubDirs(config@out.dir,subdirs)
	writeRefs(config)
	
	fastq.dir <- config@fastq.dir
	temp.dir <- config@tmp.dir
	fastq.tmp.dir <- concat(fastq.dir,'/tmp')
	#runCommand('mkdir ',fastq.dir,' -p; rm -r ',fastq.dir,'/*')
	#runCommand('mkdir ',fastq.tmp.dir,' -p; rm ',fastq.tmp.dir,'/*')
	
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		fastqfile <- concat(fastq.dir,'/',rowname,'.fastq')
		#if (file.exists(fastqfile))
		#	next
		runCommand('rm ',fastq.tmp.dir,'/*')		
		row <-  config@data[rowname,]
		folder <- row$folder
		barcode <- row$barcode
		lane <- row$lane
		dir.from <- concat(config@data.dir,row$rundata,'/Unaligned/Project_',folder,'/Sample_',folder,'/')
		dir.to <- concat(fastq.tmp.dir)
		#pattern <- concat(folder,'_',barcode,'_L00',lane,'_R1_*')
		pattern <- concat(folder,'_',barcode,'_L00',lane,'_R1_001*')
		runCommand('cp ', dir.from, pattern,' ',fastq.tmp.dir)
		runCommand('gunzip ',fastq.tmp.dir,'/*')		
		runCommand('cat ',fastq.tmp.dir,'/* > ',fastqfile)
		checkFileExists(fastqfile)
	}
}
#preprocess(config)
#preprocess(config, getSamplesForGroup(config,'KT9'))

#############################################################

solexaqa <- function(config,stem)
{
	runCommand('cd ',config@qc.dir,'; SolexaQA.pl ../fastq/',stem,'.fastq -sanger')
}
#solexaqa(config,'nextgen1-1E')
#solexaqa(config,'nextgen3-7L')

trimSolexaqa <- function(config,stem, fastq.dir=config@fastq.dir, minlength=config@minlength)
{
	outfile <- concat(stem,'.trimmed.fastq')
	if (config@force==FALSE & file.exists(concat(fastq.dir,'/',outfile)))
		return()
	olddir <- getwd()
	tryCatch({
		setwd(fastq.dir)
		runCommand('DynamicTrim.pl ',stem,'.fastq',' -phredcutoff 30')
		runCommand('LengthSort.pl ',stem,'.fastq.trimmed',' -length 36')
		checkFileExists(concat(stem,'.fastq.trimmed.single'))
		runCommand('mv ',stem,'.fastq.trimmed.single',' ',outfile)
		runCommand('rm ',stem,'.fastq.trimmed')
		runCommand('rm ',stem,'.fastq.trimmed.discard')
	}, finally=setwd(olddir))
}
#trimSolexaqa(config,'nextgen1-7A')

trimSamplesByProfile <- function(config, profile)
{
	for (rowname in rownames(config@data[which(config@data$profile==profile),]))
	{
		trimSolexaqa(config,rowname)
	}
}
#trimSamplesByProfile(config,'nextgen1')

trimSamples <- function(config, samples=config@samples)
{
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		trimSolexaqa(config,rowname)
	}
}
#trimSamples(config, getSamplesForGroup(config,'KT9'))

##########################################


maskReads <- function(config, stem, q=5,
		infile=concat(config@fastq.dir,'/',stem,getFastqExtension(config,mask=FALSE,dedup=FALSE)),
		outfile=concat(config@fastq.dir,'/',stem,getFastqExtension(config,mask=TRUE,dedup=FALSE)))
{
	printcat(infile)
	printcat(outfile)
	checkFileExists(infile)
	runCommand('fastq_masker -q ',q,' -v -i ',infile,' -o ',outfile)
	checkFileExists(outfile)
	runCommand('head ',outfile)
}
#maskReads(config, 'nextgen3-2G')

maskSamples <- function(config, samples=config@samples)
{
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		maskReads(config,rowname)
	}
}

##########################################

collapseDuplicates <- function(config, stem, 
		infile=concat(config@fastq.dir,'/',stem,getFastqExtension(config,dedup=FALSE)),
		outfile=concat(config@fastq.dir,'/',stem,getFastqExtension(config,dedup=TRUE)))
{
	printcat(infile)
	printcat(outfile)
	checkFileExists(infile)
	runCommand('fastx_collapser -v -i ',infile,' -o ',outfile)
	checkFileExists(outfile)
}
#collapseDuplicates(config, 'nextgen3-2G')
#collapseDuplicates(config, 'nextgen3-5G')

dedupSamples <- function(config, samples=config@samples)
{
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		collapseDuplicates(config,rowname)
	}
}

######################################################################

bwa <- function(fqfile, reffile, outdir, outstem=NULL)
{
	checkFileExists(fqfile)
	checkFileExists(reffile)
	sample <- stripPath(fqfile)
	sample <- stripExtension(sample)
	if (is.null(outstem)) outstem=sample
	
	samfile <- concat(outdir,'/',outstem,'.sam')
	saifile <- concat(outdir,'/',outstem,'.sai')
	bamfile <- concat(outdir,'/',outstem,'.bam')
	baifile <- concat(outdir,'/',outstem,'.bai')
	if (config@force | !file.exists(bamfile))
	{
		if (config@force | !file.exists(concat(reffile,'.amb')))
			runCommand('bwa index ',reffile)
		#delete any existing files
		runCommand('rm ',outdir,'/',outstem,'.*')
		runCommand('bwa aln -n 4 ',reffile,' ',fqfile,' > ',saifile) # -k 3 -n 4
		#runCommand('bwa aln ',reffile,' ',fqfile,' > ',saifile) # -k 3 -n 4
		runCommand('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
		checkFileExists(samfile)
		sam2bam(samfile,bamfile)
	}
	else printcat('skipping bwa command because bamfile already exists: ',bamfile)
	runCommand('rm ',samfile)
	runCommand('rm ',saifile)
	return(bamfile)
}

runBwa <- function(config, stem, ref=config@data[stem,'ref'], trim=config@trim, dedup=TRUE)
{
	reffile <- getRefFile(config,ref)
	fqfile <- getFastqFilename(config,stem)
	outstem <- concat(stem,'__',ref)
	bamfile <- bwa(fqfile,reffile,config@bam.dir,outstem)
	addReadGroup(config,outstem)
	print(getMapStats(config,bamfile))
	return(bamfile)
}
#runBwa(config,'nextgen3-2H')

#############################################################################################

tmap <- function(fqfile, reffile, outdir, outstem=NULL)
{
	checkFileExists(fqfile)
	checkFileExists(reffile)
	sample <- stripPath(fqfile)
	sample <- stripExtension(sample)
	if (is.null(outstem)) outstem=sample
	
	samfile <- concat(outdir,'/',outstem,'.sam')
	saifile <- concat(outdir,'/',outstem,'.sai')
	bamfile <- concat(outdir,'/',outstem,'.bam')
	baifile <- concat(outdir,'/',outstem,'.bai')
	if (config@force | !file.exists(bamfile))
	{
		if (config@force | !file.exists(concat(reffile,'.tmap.anno')))
			runCommand('tmap index -f ',reffile)
		#delete any existing files
		runCommand('rm ',outdir,'/',outstem,'.*')
		runCommand('tmap map1 --fn-fasta ',reffile,' --fn-reads ',fqfile,' --reads-format fastq --duplicate-window -1 -s ',samfile)
		checkFileExists(samfile)
		sam2bam(samfile,bamfile)
	}
	else printcat('skipping bwa command because bamfile already exists: ',bamfile)
	runCommand('rm ',samfile)
	runCommand('rm ',saifile)
	return(bamfile)
}

runTmap <- function(config, stem, ref=config@data[stem,'ref'], trim=config@trim, dedup=TRUE)
{
	reffile <- getRefFile(config,ref)
	fqfile <- getFastqFilename(config,stem)
	outstem <- concat(stem,'__',ref)
	bamfile <- tmap(fqfile,reffile,config@bam.dir,outstem)
	addReadGroup(config,outstem)
	print(getMapStats(config,bamfile))
	return(bamfile)
}
#runTmap(config,'nextgen3-2H')

###################################################################################################

mapReadsByProfile <- function(config, profile)
{
	for (rowname in rownames(config@data[which(config@data$profile==profile),]))
	{
		printcat('runTmap: ',rowname)
		try(runTmap(config,rowname))
		#try(runBwa(config,rowname))
	}
	#outputBams(config, samples)
}
#mapReadsByProfile(config,'nextgen1')

mapReads <- function(config, samples=config@samples)
{
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		#printcat('runBwa: ',rowname)
		#try(runBwa(config,rowname))
		try(runTmap(config,rowname))
	}
	#outputBams(config, samples)
}
#mapReads(config)
#mapReads(config, getSamplesForGroup(config,'MP-424'))
#mapReads(config, getSamplesForSubject(config,'10464592'))


writeConsensusForBams <- function(config,ids=rownames(config@data))
{
	for (id in ids)
	{
		ref <- config@data[id,'ref']
		sample <- concat(id,'__',ref)
		try(writeConsensusForBam(config,sample))
		#writeCoverageForBam(config,sample)
	}
}
#writeConsensusForBams(config)

#################################################################

mergeBamsForSample <- function(config, sample, bam.dir=config@bam.dir, out.dir=config@bam.dir)
{
	outfile <- concat(out.dir,'/',sample,'.bam')
	filenames <- c()
	for (rowname in rownames(config@data[which(config@data$sample==sample),]))
	{
		print(rowname)
		ref <- config@data[rowname,'ref']
		filename <- concat(bam.dir,'/',rowname,'__',ref,'.bam')
		filenames <- c(filenames,filename)
	}
	print(filenames)
	mergeBamFiles(config,filenames,outfile)
	return(outfile)
}

mergeBamsForSamples <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		try(mergeBamsForSample(config,sample))
	}
}
#mergeBamsForSamples(config)

################################################################3

fixBaiFiles <- function(config, bam.dir=config@bam.dir)
{
	for (filename in list.files(bam.dir, pattern='\\.bam$'))
	{
		fixBaiFile(config,filename)
	}
}
#fixBaiFiles(config)

##################################################################3

filterBam <- function(config, sample, bam.dir=config@bam.dir, map.quality=config@map.quality)
{
	#bamtools filter -in out/bam/KT9.plasmid.bam -out out/bam/KT9.plasmid.filtered.bam -mapQuality ">30" 
	infile <- concat(bam.dir,'/',sample,'.bam')
	outfile <- concat(bam.dir,'/',sample,'.filtered.bam')
	checkFileExists(infile)
	str <- 'bamtools filter'
	str <- concat(str,' -in ',infile)
	str <- concat(str,' -out ',outfile)
	str <- concat(str,' -mapQuality "',map.quality,'"')
	runCommand(str)
	checkFileExists(outfile)
	runCommand('samtools index ',outfile)
	return(outfile)
}
#filterBam(config,'nextgen1-1E__HCV-KT9')

filterBams <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		printcat('filter_bam: ',sample)
		try(filterBam(config,sample))
	}
}
#filterBams(config)

countCodonsForSample <- function(config, sample)
{
	for (id in config@data[which(config@data$sample==sample),'id'])
	{
		#try(countCodonsForSample(config,id))
		countCodonsForId(config,id)
	}
}

countCodons <- function(config, ids=rownames(config@data))
{
	for (id in ids)
	{
		#try(countCodonsForSample(config,id))
		countCodonsForId(config,id)
	}
}
#countCodons(config)


makePiecharts <- function(config)
{
	sweaveToPdf(concat(Sys.getenv("VARDB_RUTIL_HOME"),'/nextgen_piecharts.rnw'))
}
#makePiecharts(config)

concatTablesByGroup <- function(config, groups=config@groups)
{
	for (group in groups)
	{
		for (type in splitFields('codons,aa'))
		{
			infiles <- concat(config@tables.dir,'/table-',type,'-',group,'*.txt /dev/null')
			outfile <- concat(config@tables.dir,'/group-',type,'-',group,'.txt')
			runCommand('tail -n +1 ',infiles,' > ',outfile)
		}
	}
}
#concatTablesByGroup(config)
#concatTablesByGroup(config,'BMS-605339')

###############################################################

analyzeReadsForSample <- function(config,sample)
{
	mapReads(config,sample)
	#mergeBamsForSamples(config,sample)
	if (config@filter) filterBams(config,sample)
	#writeConsensusForBam(config,sample)
	#findVariants(config,sample)
	#exportPileup(config,sample)
	countCodonsForSample(config,sample)
}
#analyzeReadsForSample(config,'KT9.random__HCV-KT9')

analyzeReadsForGroup <- function(config,group,make.pdf=FALSE)
{
	for (sample in getSamplesForGroup(config,group))
	{
		try(analyzeReadsForSample(config,sample))
	}
	#writeCodonTables(config,group)
	#writeAminoAcidTables(config,group)
	#concatTablesByGroup(config,group)
	reportAminoAcidChangesForGroup(config,group,make.pdf=make.pdf)
	makeReferenceVsVariantTables(config,group=group,minreads=100)
}
#analyzeReadsForGroup(config,'MP-424')
#analyzeReadsForGroup(config,'hcv_infection')
#analyzeReadsForGroup(config,'KT9')


analyzeReads <- function(config, groups=config@groups)
{
	for (group in groups)
	{
		try(analyzeReadsForGroup(config,group))
	}
}
#analyzeReads(config)

#analyzeReads<- function(config)
#{
#	preprocess(config)
#	trim_all(config)
#	##remove_exact_duplicates_for_all_samples(config)
#	map_reads(config)
#	#merge_bams(config)
#	##realign_indels(config,stem)
#	#recalibrate(config,stem) #concat(stem,'.realigned'))
#	#output_bam(config,stem,'recal')
#	
#	#call_variants(config,stem)
#	#filter_variants(config,stem)
#	#unmerge_bams(config)
#	output_bams(config)
#	fix_bai_files(config)
#	filter_bams(config)
#	export_pileup(config)
#	count_codons(config)
#	make_piecharts(config)
#	#export_unmapped_reads(config,concat(stem,'.nodup'))
#}


