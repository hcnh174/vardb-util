#preprocess <- function(config, samples=config@samples, subdirs='tmp,ref,fastq,tmp,bam,unmapped,vcf,pileup,qc,counts,tables,charts,consensus,coverage')
#{
#	makeSubDirs(config@out.dir,subdirs)
#	writeRefs(config)
#	
#	fastq.dir <- config@fastq.dir
#	temp.dir <- config@tmp.dir
#	sort.dir <- concat(config@tmp.dir,'/sort')
#	runCommand('mkdir ',sort.dir,' -p')
#	runCommand('rm -r ',sort.dir,'/*')
#	
#	stems <- getStemsForSamples(config,samples)
#	#stems <-unique(config@data[which(config@data$sample %in% samples),'stem'])
#	for (stem in stems)
#	{
#		dir.to <- concat(sort.dir,'/',stem)
#		runCommand('mkdir ',dir.to)
#	}
#	
#	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
#	{
#		row <-  config@data[rowname,]
#		folder <- row$folder
#		barcode <- row$barcode
#		lane <- row$lane
#		#dir.from <- concat('../data/',row$rundata,'/Unaligned/Project_',folder,'/Sample_',folder,'/')
#		dir.from <- concat(config@data.dir,'/',row$rundata,'/Unaligned/Project_',folder,'/Sample_',folder,'/')
#		dir.to <- concat(sort.dir,'/',row$stem)
#		pattern <- concat('^',folder,'_',barcode,'_L00',lane,'_R1_.*\\.fastq.*')
#		for (filename in list.files(dir.from,pattern))
#		{
#			ext <- ifelse(getFileExtension(filename)=='gz','.fastq.gz','.fastq')
#			newfilename <- concat(dir.to,'/',stem,'.',row$region,'.',rowname,ext)
#			#printcat('cp ', dir.from, filename,' ',newfilename)
#			runCommand('cp ', dir.from, filename,' ',newfilename)
#		}
#	}
#	
#	for (stem in stems)
#	{
#		dir.from <- concat(sort.dir,'/',stem)
#		runCommand('gunzip ',dir.from,'/*')
#	}
#	
#	for (stem in stems)
#	{
#		dir.from <- concat(sort.dir,'/',stem)
#		fastqfile <- concat(fastq.dir,'/',stem,'.fastq')
#		runCommand('cat ',dir.from,'/* > ',fastqfile)
#		checkFileExists(fastqfile)
#	}
#	
#	#config <- getRawReadCounts(config)
#	#runCommand('rm -r ',temp.dir,'/*')
#}
##preprocess(config)

preprocess <- function(config, samples=config@samples, subdirs='tmp,ref,fastq,tmp,bam,unmapped,vcf,pileup,qc,counts,tables,charts,consensus,coverage')
{
	makeSubDirs(config@out.dir,subdirs)
	writeRefs(config)
	
	fastq.dir <- config@fastq.dir
	temp.dir <- config@tmp.dir
	#fastq.dir <- concat(config@tmp.dir,'/fastq')
	fastq.tmp.dir <- concat(fastq.dir,'/tmp')
	runCommand('mkdir ',fastq.dir,' -p; rm -r ',fastq.dir,'/*')
	runCommand('mkdir ',fastq.tmp.dir,' -p; rm ',fastq.tmp.dir,'/*')
	
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
		pattern <- concat(folder,'_',barcode,'_L00',lane,'_R1_*')
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
#solexa_qa(config,'KT9.plasmid__KT9')

trimSolexaqa <- function(config,stem, fastq.dir=config@fastq.dir)
{
	outfile <- concat(stem,'.trimmed.fastq')
	if (file.exists(concat(fastq.dir,'/',outfile)))
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

#trimSamples <- function(config, samples=config@samples)
#{
#	for (stem in getStemsForSamples(config,samples))
#	{
#		trimSolexaqa(config,stem)
#	}
#}
trimSamples <- function(config, samples=config@samples)
{
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		trimSolexaqa(config,rowname)
	}
}
#trimSamples(config)
#
#mergeTrimmedFiles <- function(config, samples=config@samples)
#{
#	for (sample in samples)
#	{
#		stems <- getStemsForSamples(config,sample)
#		filenames <- c()
#		for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
#		{
#			filenames <- c(filenames,concat(config@fastq.dir,'/',rowname,'.trimmed.fastq'))
#		}
#		printcat(filenames)
#	}
#}

######################################################################

addReadGroups <- function(config, sample, tmp.dir=config@tmp.dir)
{
	#tmp.dir <- config@tmp.dir
	samfile <- concat(tmp.dir,'/',sample,'.sam')
	bamfile <- concat(tmp.dir,'/',sample,'.bam')
	str <- 'java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar'
	str <- concat(str,' INPUT=',samfile)
	str <- concat(str,' OUTPUT=',bamfile)
	str <- concat(str,' RGSM="',sample,'"')
	str <- concat(str,' RGLB="',sample,'"')
	str <- concat(str,' RGID="',sample,'"')
	str <- concat(str,' RGPL=illumina')
	str <- concat(str,' RGPU=barcode')
	str <- concat(str,' SORT_ORDER=coordinate')
	str <- concat(str,' CREATE_INDEX=true')
	runCommand(str)
	checkFileExists(bamfile)
	baifile <- concat(tmp.dir,'/',sample,'.bai')
	checkFileExists(baifile)
	return(bamfile)
}
#addReadGroups(config,'110617HBV.HBV07@HBV-RT')

#bwa <- function(fqfile, reffile, outdir, outstem=NULL)
#{
#	sample <- stripPath(fqfile)
#	sample <- stripExtension(sample)
#	if (is.null(outstem)) outstem=sample
#	
#	samfile <- concat(outdir,'/',outstem,'.sam')
#	saifile <- concat(outdir,'/',outstem,'.sai')
#	
#	if (!file.exists(concat(reffile,'.amb')))
#		runCommand('bwa index ',reffile)
#	runCommand('bwa aln ',reffile,' ',fqfile,' > ',saifile)
#	runCommand('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
#	checkFileExists(samfile)
#	
#	bamfile <- addReadGroups(config,outstem,outdir)
#	return(bamfile)
#}
##bwa('fastq/test.fastq','ref/hcv.fasta','tmp')

#runBwa <- function(config, sample, ref=getRefForSample(sample), trim=config@trim)
#{
#	print(ref)
#	stem <- getStemsForSamples(config,sample)
#	fastq.ext <- ifelse(trim,'.trimmed.fastq','.fastq')
#	reffile <- getRefFile(config,ref)
#	print(reffile)
#	fqfile <- concat(config@fastq.dir,'/',stem,fastq.ext)	
#	outstem <- concat(stem,'__',ref)
#	bwa(fqfile,reffile,config@bam.dir,outstem)
#	#addReadGroups(config,sample)	
#}
##runBwa(config,getSamplesForGroup(config,'MP-424')[1])


#mapReadsByProfile <- function(config, profile)
#{
#	for (rowname in rownames(config@data[which(config@data$profile==profile),]))
#	{
#		ref <- config@data[rowname,'ref']
#		try(runBwa(config,rowname,ref))
#		return()
#	}
#}
#mapReadsByProfile(config,'nextgen1')

#mapReads <- function(config, samples=config@samples)
#{
#	for (sample in samples)
#	{
#		printcat('runBwa: ',sample)
#		try(runBwa(config,sample))
#	}
#	outputBams(config, samples)
#}

fixBaiFile <- function(config, bamfile)
{
	stem <- stripExtension(bamfile)
	oldbaifile <- concat(stem,'.bai')
	baifile <- concat(stem,'.bam.bai')
	print(oldbaifile)
	print(baifile)
	#oldbaifile <- concat(bam.dir,'/',stem,'.bai')
	#baifile <- concat(bam.dir,'/',stem,'.bam.bai')
	if (file.exists(oldbaifile) & !file.exists(baifile))
		runCommand('mv "',oldbaifile,'" "',baifile,'"')
}
#fixBaiFile(config)

bwa <- function(fqfile, reffile, outdir, outstem=NULL)
{
	sample <- stripPath(fqfile)
	sample <- stripExtension(sample)
	if (is.null(outstem)) outstem=sample
	
	samfile <- concat(outdir,'/',outstem,'.sam')
	saifile <- concat(outdir,'/',outstem,'.sai')
	bamfile <- concat(outdir,'/',outstem,'.bam')
	if (config@force | !file.exists(bamfile))
	{		
		if (config@force | !file.exists(concat(reffile,'.amb')))
			runCommand('bwa index ',reffile)
		runCommand('bwa aln ',reffile,' ',fqfile,' > ',saifile) # -k 3 -n 4
		runCommand('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
		checkFileExists(samfile)
	
		addReadGroups(config,outstem,outdir)
	}
	runCommand('rm ',samfile)
	runCommand('rm ',saifile)
	fixBaiFile(config,bamfile)
}
#bwa('fastq/test.fastq','ref/hcv.fasta','tmp')

runBwa <- function(config, stem, ref=config@data[stem,'ref'], trim=config@trim)#ref=getRefForSample(sample)
{
	print(ref)
	fastq.ext <- ifelse(trim,'.trimmed.fastq','.fastq')
	reffile <- getRefFile(config,ref)
	print(reffile)
	fqfile <- concat(config@fastq.dir,'/',stem,fastq.ext)	
	outstem <- concat(stem,'__',ref)
	bwa(fqfile,reffile,config@bam.dir,outstem)
	#addReadGroups(config,sample)	
}
#runBwa(config,'nextgen4-7A')

mapReadsByProfile <- function(config, profile)
{
	for (rowname in rownames(config@data[which(config@data$profile==profile),]))
	{
		printcat('runBwa: ',rowname)
		try(runBwa(config,rowname))
	}
	#outputBams(config, samples)
}
#mapReadsByProfile(config,'nextgen1')

mapReads <- function(config, samples=config@samples)
{
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		printcat('runBwa: ',rowname)
		try(runBwa(config,rowname))
	}
	#outputBams(config, samples)
}
#mapReads(config)
#mapReads(config, getSamplesForGroup(config,'MP-424'))
#mapReads(config, getSamplesForSubject(config,'10464592'))

#################################################################

mergeBamsForSample <- function(config, sample, bam.dir=config@bam.dir, out.dir=config@bam.dir)
{
	outfile <- concat(out.dir,'/',sample,'.bam')
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
	str <- concat(str,' CREATE_INDEX=true')
	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
	
	for (rowname in rownames(config@data[which(config@data$sample==sample),]))
	{
		ref <- config@data[rowname,'ref']
		infile <- concat(bam.dir,'/',rowname,'__',ref,'.bam')
		str <- concat(str,' INPUT=',infile)
	}
	str <- concat(str,' OUTPUT=',outfile)
	#print(str)
	runCommand(str)
	checkFileExists(outfile)
	fixBaiFile(config,outfile)
}
#mergeBamsForSample(config,'mitsuto_fujita__HCV-KT9_mitsuto_fujita')

mergeBamsForSamples <- function(config, samples=config@samples)
{
	#stems <- getStemsForSamples(config,samples)
	#3for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	for (sample in samples)
	{
		mergeBamsForSample(config,sample)
	}
}
#mergeBamsForSamples(config)

############################################
#
#outputBam <- function(config,sample,suffix='')
#{
#	infile <- concat(config@tmp.dir,'/',sample,suffix,'.bam')
#	outfile <- concat(config@bam.dir,'/',sample,'.bam')
#	runCommand('cp ',infile,' ',outfile)
#	checkFileExists(outfile)
#	
#	infile <- concat(config@tmp.dir,'/',sample,suffix,'.bai')
#	outfile <- concat(config@bam.dir,'/',sample,'.bai')
#	runCommand('cp ',infile,' ',outfile)
#	checkFileExists(outfile)
#	#runCommand('samtools index ',outfile)
#}
##outputBam(config,'merged','.realigned.recal')
#
#outputBams <- function(config,samples=config@samples,suffix='') #ifelse(config@trim,'.trimmed',''))
#{
#	for (sample in samples)
#	{
#		outputBam(config,sample,suffix)
#	}
#	fixBaiFiles(config)
#}
##outputBams(config,getSamplesForGroup(config,'NS5A_L31V_Y93H_mutations_maintained'))
##outputBams(config,getSamplesForGroup(config,'MP-424'))
##outputBams(config,getSamplesForSubject(config,'10464592'))

writeConsensusForBams <- function(config,samples=config@samples)
{
	for (sample in samples)
	{
		try(writeConsensusForBam(config,sample))
		#writeCoverageForBam(config,sample)
	}
}
#writeConsensusForBams(config)

################################################################3

fixBaiFiles <- function(config, bam.dir=config@bam.dir)
{
	for (filename in list.files(bam.dir, pattern='\\.bam$'))
	{
		fixBaiFile(config,filename)
	}
}
#fixBaiFiles(config)

#fixBaiFiles <- function(config, bam.dir=config@bam.dir)
#{
#	for (filename in list.files(bam.dir, pattern='\\.bam$'))
#	{
#		stem <- stripExtension(filename)
#		oldbaifile <- concat(bam.dir,'/',stem,'.bai')
#		baifile <- concat(bam.dir,'/',stem,'.bam.bai')
#		if (file.exists(oldbaifile) & !file.exists(baifile))
#			runCommand('mv "',oldbaifile,'" "',baifile,'"')
#	}
#}
##fixBaiFiles(config)

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

filterBams <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		printcat('filter_bam: ',sample)
		try(filterBam(config,sample))
	}
}
#filterBams(config)

####################################################

exportPileupForSample <- function(config,sample,filtered=TRUE)
{
	ref <- getRefForSample(sample)
	samplename <- ifelse(filtered,concat(sample,'.filtered'), sample)
	runCommand('python $VARDB_RUTIL_HOME/export_pileup.py ', samplename,' ',ref,' ',config@bam.dir,' ',config@pileup.dir)
}
#exportPileupForSample(config,'110617HBV-1.10348001.20020530__HBV-RT',filtered=FALSE)

exportPileup <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		exportPileupForSample(config,sample,filtered=FALSE)
		if (config@filter)
			exportPileupForSample(config,sample,filtered=TRUE)
	}
}
#exportPileup(config,getSamplesForSubject(config,'10464592'))

################################################

countCodons <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		try(countCodonsForSample(config,sample))
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
			infiles <- concat(config@tables.dir,'/table-',type,'-',group,'*.txt')
			outfile <- concat(config@tables.dir,'/group-',type,'-',group,'.txt')
			runCommand('tail -n +1 ',infiles,' > ',outfile)
		}
	}
}
#concatTablesByGroup(config)
#concatTablesByGroup(config,'PXB0219-0011')

###############################################################

analyzeReadsForSample <- function(config,sample)
{
	mapReads(config,sample)
	if (config@filter) filterBams(config,sample)
	writeConsensusForBams(config,sample)
	exportPileup(config,sample)
	countCodons(config,sample)
}
#analyzeReadsForSample(config,'KT9.random__HCV-KT9')

analyzeReadsForGroup <- function(config,group)
{
	samples <- getSamplesForGroup(config,group)
	analyzeReadsForSample(config,samples)
	writeCodonTables(config,group)
	concatTablesByGroup(config,group)
	makeAminoAcidBarcharts(config,group)
}
#analyzeReadsForGroup(config,'MP-424')
#analyzeReadsForGroup(config,'hcv_infection')
#analyzeReadsForGroup(config,'KT9')

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
