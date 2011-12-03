preprocess <- function(config, samples=config@samples, subdirs='tmp,ref,fastq,tmp,bam,unmapped,vcf,pileup,basecounts,qc,counts,tables,charts,consensus,coverage')
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
#solexaqa(config,'nextgen3-7L')

trimSolexaqa <- function(config,stem, fastq.dir=config@fastq.dir, minlength=config@minlength)
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

trimSamples <- function(config, samples=config@samples)
{
	stems <- getStemsForSamples(config,samples)
	for (rowname in rownames(config@data[which(config@data$stem %in% stems),]))
	{
		trimSolexaqa(config,rowname)
	}
}
#trimSamples(config)

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

#
#addReadGroups <- function(config, sample, tmp.dir=config@tmp.dir)
#{
#	#tmp.dir <- config@tmp.dir
#	samfile <- concat(tmp.dir,'/',sample,'.sam')
#	bamfile <- concat(tmp.dir,'/',sample,'.bam')
#	str <- 'java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar'
#	str <- concat(str,' INPUT=',samfile)
#	str <- concat(str,' OUTPUT=',bamfile)
#	str <- concat(str,' RGSM="',sample,'"')
#	str <- concat(str,' RGLB="',sample,'"')
#	str <- concat(str,' RGID="',sample,'"')
#	str <- concat(str,' RGPL=illumina')
#	str <- concat(str,' RGPU=barcode')
#	str <- concat(str,' SORT_ORDER=coordinate')
#	str <- concat(str,' CREATE_INDEX=true')
#	runCommand(str)
#	checkFileExists(bamfile)
#	baifile <- concat(tmp.dir,'/',sample,'.bai')
#	checkFileExists(baifile)
#	return(bamfile)
#}
##addReadGroups(config,'110617HBV.HBV07@HBV-RT')
#
#fixBaiFile <- function(config, bamfile)
#{
#	stem <- stripExtension(bamfile)
#	oldbaifile <- concat(stem,'.bai')
#	baifile <- concat(stem,'.bam.bai')
#	print(oldbaifile)
#	print(baifile)
#	if (file.exists(oldbaifile) & !file.exists(baifile))
#		runCommand('mv "',oldbaifile,'" "',baifile,'"')
#}
##fixBaiFile(config)
#
#bwa <- function(fqfile, reffile, outdir, outstem=NULL)
#{
#	sample <- stripPath(fqfile)
#	sample <- stripExtension(sample)
#	if (is.null(outstem)) outstem=sample
#	
#	samfile <- concat(outdir,'/',outstem,'.sam')
#	saifile <- concat(outdir,'/',outstem,'.sai')
#	bamfile <- concat(outdir,'/',outstem,'.bam')
#	baifile <- concat(outdir,'/',outstem,'.bai')
#	baifile2 <- concat(outdir,'/',outstem,'.bam.bai')
#	if (config@force | !file.exists(bamfile))
#	{
#		if (config@force | !file.exists(concat(reffile,'.amb')))
#			runCommand('bwa index ',reffile)
#		#delete any existing files
#		runCommand('rm ',outdir,'/',outstem,'.*')
#		runCommand('bwa aln -n 4 ',reffile,' ',fqfile,' > ',saifile) # -k 3 -n 4
#		#runCommand('bwa aln ',reffile,' ',fqfile,' > ',saifile) # -k 3 -n 4
#		runCommand('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
#		checkFileExists(samfile)
#	
#		addReadGroups(config,outstem,outdir)
#	}
#	runCommand('rm ',samfile)
#	runCommand('rm ',saifile)
#	runCommand('mv "',baifile,'" "',baifile2,'"')
#	return(bamfile)
#}
##bwa('fastq/test.fastq','ref/hcv.fasta','tmp')

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
	runCommand('rm ',samfile)
	runCommand('rm ',saifile)
	return(bamfile)
}

#runBwa <- function(config, stem, ref=config@data[stem,'ref'], trim=config@trim)
#{
#	print(ref)
#	fastq.ext <- ifelse(trim,'.trimmed.fastq','.fastq')
#	reffile <- getRefFile(config,ref)
#	print(reffile)
#	fqfile <- concat(config@fastq.dir,'/',stem,fastq.ext)	
#	outstem <- concat(stem,'__',ref)
#	bamfile <- bwa(fqfile,reffile,config@bam.dir,outstem)
#	print(getMapStats(config,bamfile))
#	return(bamfile)
#}
##runBwa(config,'nextgen1-3F') #'nextgen1-2E')#'nextgen2-5I')
#
#runBwa <- function(config, stem, ref=unique(config@data[which(config@data$id==stem | config@data$stem==stem),'ref']), trim=config@trim)
#{
#	printcat('stem=',stem,', ref=',ref)
#	fastq.ext <- ifelse(trim,'.trimmed.fastq','.fastq')
#	reffile <- getRefFile(config,ref)
#	print(reffile)
#	fqfile <- concat(config@fastq.dir,'/',stem,fastq.ext)
#	print(fqfile)
#	outstem <- concat(stem,'__',ref)
#	bamfile <- bwa(fqfile,reffile,config@bam.dir,outstem)
#	print(getMapStats(config,bamfile))
#	return(bamfile)
#}
##runBwa(config,'nextgen4-7A')

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
#runBwa(config,'nextgen3-2G')

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

#mergeBamsForSample <- function(config, sample, bam.dir=config@bam.dir, out.dir=config@bam.dir)
#{
#	tmpoutfile <- concat(out.dir,'/tmp-',sample,'.bam')
#	
#	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
#	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
#	str <- concat(str,' CREATE_INDEX=true')
#	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
#	
#	for (rowname in rownames(config@data[which(config@data$sample==sample),]))
#	{
#		ref <- config@data[rowname,'ref']
#		infile <- concat(bam.dir,'/',rowname,'__',ref,'.bam')
#		str <- concat(str,' INPUT=',infile)
#	}
#	str <- concat(str,' OUTPUT=',tmpoutfile)
#	runCommand(str)
#	checkFileExists(tmpoutfile)
#	
#	outfile <- concat(out.dir,'/',sample,'.bam')
#	str <- 'java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar'
#	str <- concat(str,' INPUT=',tmpoutfile)
#	str <- concat(str,' OUTPUT=',outfile)
#	str <- concat(str,' RGSM="',sample,'"')
#	str <- concat(str,' RGLB="',sample,'"')
#	str <- concat(str,' RGID="',sample,'"')
#	str <- concat(str,' RGPL=illumina')
#	str <- concat(str,' RGPU=barcode')
#	str <- concat(str,' SORT_ORDER=coordinate')
#	str <- concat(str,' CREATE_INDEX=true')
#	
#	runCommand(str)
#	checkFileExists(outfile)
#	baifile <- concat(out.dir,'/',sample,'.bai')
#	checkFileExists(baifile)
#	fixBaiFile(config,outfile)
#	runCommand('rm ',out.dir,'/tmp-',sample,'.*')
#	return(outfile)	
#}
##mergeBamsForSample(config,'sample__HCV-KT9')

#mergeBamsForSample <- function(config, sample, bam.dir=config@bam.dir, out.dir=config@bam.dir)
#{
#	outfile <- concat(out.dir,'/',sample,'.bam')
#	
#	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
#	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
#	str <- concat(str,' CREATE_INDEX=true')
#	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
#	
#	for (rowname in rownames(config@data[which(config@data$sample==sample),]))
#	{
#		ref <- config@data[rowname,'ref']
#		infile <- concat(bam.dir,'/',rowname,'__',ref,'.bam')
#		str <- concat(str,' INPUT=',infile)
#	}
#	str <- concat(str,' OUTPUT=',outfile)
#	#print(str)
#	runCommand(str)
#	checkFileExists(outfile)
#	fixBaiFile(config,outfile)
#}
##mergeBamsForSample(config,'sample__HCV-KT9_sample')

mergeBamsForSample <- function(config, sample, bam.dir=config@bam.dir, out.dir=config@bam.dir)
{
	outfile <- concat(out.dir,'/',sample,'.bam')
	filenames <- c()
	for (rowname in rownames(config@data[which(config@data$sample==sample),]))
	{
		print(rowname)
		ref <- config@data[rowname,'ref']
		filename <- concat(bam.dir,'/',rowname,'__',ref,'.bam')
		#clipregion <- config@data[rowname,'clipregion']
		#if (!is.na(clipregion))
		#	filename <- trimBamToRegion(config,concat(rowname,'__',ref),as.numeric(clipregion))
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


############################################

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

#exportPileupForSample <- function(config,sample,ref=getRefForSample(sample), bam.dir=config@bam.dir, out.dir=config@pileup.dir, filtered=TRUE)
#{
#	#ref <- getRefForSample(sample)
#	samplename <- ifelse(filtered,concat(sample,'.filtered'), sample)
#	#runCommand('python $VARDB_RUTIL_HOME/export_pileup.py ', samplename,' ',ref,' ',config@bam.dir,' ',config@pileup.dir)
#	runCommand('python $VARDB_RUTIL_HOME/export_pileup.py ', samplename,' ',ref,' ',bam.dir,' ',out.dir)
#}
##exportPileupForSample(config,'110617HBV-1.10348001.20020530__HBV-RT',filtered=FALSE)

exportPileupForSample <- function(config,sample,ref=getRefForSample(sample), bam.dir=config@bam.dir, out.dir=config@pileup.dir, filtered=FALSE)
{
	samplename <- ifelse(filtered,concat(sample,'.filtered'), sample)
	runCommand('python $VARDB_RUTIL_HOME/export_pileup.py ', samplename,' ',ref,' ',bam.dir,' ',out.dir)
}
#exportPileupForSample(config,'sample','ref','out','out')

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
	mergeBamsForSamples(config,sample)
	if (config@filter) filterBams(config,sample)
	writeConsensusForBams(config,sample)
	findVariants(config,sample)
	exportPileup(config,sample)
	countCodons(config,sample)
}
#analyzeReadsForSample(config,'KT9.random__HCV-KT9')

analyzeReadsForGroup <- function(config,group)
{
	for (sample in getSamplesForGroup(config,group))
	{
		try(analyzeReadsForSample(config,sample))
	}
	writeCodonTables(config,group)
	writeAminoAcidTables(config,group)
	concatTablesByGroup(config,group)
	reportAminoAcidChanges(config,group)
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


