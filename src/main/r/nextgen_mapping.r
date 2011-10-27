preprocess <- function(config, samples=config@samples, subdirs='tmp,ref,fastq,tmp,bam,unmapped,vcf,pileup,qc,counts,tables,consensus')
{
	runCommand('mkdir ',config@out.dir,' -p')
	for (subdir in splitFields(subdirs))
	{
		subdir <- concat(config@out.dir,'/',subdir)
		#runCommand('mkdir ',subdir,' -p; rm ',subdir,'/*')
		runCommand('mkdir ',subdir,' -p')
	}
	writeRefs(config)
	
	fastq.dir <- config@fastq.dir
	temp.dir <- config@tmp.dir
	runCommand('rm -r ',temp.dir,'/*')
	for (sample in samples)
	{
		dir.to <- concat(temp.dir,'/',sample)
		runCommand('mkdir ',dir.to)
	}
	
	for (rowname in rownames(config@data[which(config@data$sample %in% samples),]))
	{
		row <-  config@data[rowname,]
		sample <- row$sample
		dir.to <- concat(temp.dir,'/',sample)
		folder <- row$folder
		barcode <- row$barcode
		lane <- row$lane
		dir.from <- concat('../data/',row$rundata,'/Unaligned/Project_',folder,'/Sample_',folder,'/')
		filename <- concat(folder,'_',barcode,'_L00',lane,'_R1_*.fastq.gz')
		runCommand('cp ', dir.from, filename,' ',dir.to)
		filename <- concat(folder,'_',barcode,'_L00',lane,'_R1_*.fastq')
		runCommand('cp ', dir.from, filename,' ',dir.to)
	}	
	
	for (sample in samples)
	{
		dir.from <- concat(temp.dir,'/',sample)
		runCommand('gunzip ',dir.from,'/*')
	}
	
	for (sample in samples)
	{
		dir.from <- concat(temp.dir,'/',sample)
		fastqfile <- concat(fastq.dir,'/',sample,'.fastq')
		runCommand('cat ',dir.from,'/* > ',fastqfile)
		checkFileExists(fastqfile)
	}
	
	#config <- getRawReadCounts(config) 
	#runCommand('rm -r ',temp.dir,'/*')
}
#preprocess(config)


#############################################################

trimSolexaqa <- function(config,sample)
{
	olddir <- getwd()
	setwd(config@fastq.dir)
	runCommand('DynamicTrim.pl ',sample,'.fastq',' -phredcutoff 30')
	runCommand('LengthSort.pl ',sample,'.fastq.trimmed',' -length 36')
	checkFileExists(concat(sample,'.fastq.trimmed.single'))
	runCommand('mv ',sample,'.fastq.trimmed.single',' ',sample,'.trimmed.fastq')
	runCommand('rm ',sample,'.fastq.trimmed')
	runCommand('rm ',sample,'.fastq.trimmed.discard')
	setwd(olddir)
}
#trimSolexaqa(config,'PXB0220-0002.wk13')

trimSamples <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		trimSolexaqa(config,sample)
	}
}
#trimSamples(config)

######################################################################

addReadGroups <- function(config, sample)
{
	tmp.dir <- config@tmp.dir
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
}
#addReadGroups(config,'110617HBV.HBV07@HBV-RT')

runBwa <- function(config, sample)
{
	fastq.ext <- ifelse(config@trim,'.trimmed.fastq','.fastq')
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	fqfile <- concat(config@fastq.dir,'/',sample,fastq.ext)
	tmp.dir <- config@tmp.dir
	
	samfile <- concat(tmp.dir,'/',sample,'.sam')
	saifile <- concat(tmp.dir,'/',sample,'.sai')
	
	runCommand('bwa index ',reffile)
	runCommand('bwa aln ',reffile,' ',fqfile,' > ',saifile)
	runCommand('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
	checkFileExists(samfile)
	
	addReadGroups(config,sample)
	
	#runCommand('rm ',samfile)
	#runCommand('rm ',saifile)
}
#runBwa(config,'110617HBV-1.10348001.20020530__HBV-RT')

mapReads <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		print(concat('runBwa: ',sample))
		try(runBwa(config,sample))
	}
	outputBams(config, samples)
}
#mapReads(config)
#mapReads(config, getSamplesForSubject(config,'10464592'))

############################################

mergeBamsForRef <- function(config,ref)
{
	tmp.dir <- config@tmp.dir
	outfile <- concat(tmp.dir,'/merged@',ref,'.bam')
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
	str <- concat(str,' CREATE_INDEX=true')
	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
	
	for (sample in unique(config@data[which(config@data$ref==ref),'sample']))
	{
		ref <- getRefForSample(sample)
		infile <- concat(tmp.dir,'/',sample,'.bam')
		str <- concat(str,' INPUT=',infile)
	}
	str <- concat(str,' OUTPUT=',outfile)
	runCommand(str)
}
#mergeBamsForRef(config,ref)

mergeBams <- function(config)
{
	for (ref in rownames(config@refs))
	{
		mergeBamsForRef(config,ref)
	}
}
#mergeBams(config)

###################################################################

analyzeCovariates <- function(config,stem,suffix)
{
	tmp.dir <- config@tmp.dir
	bamfile <- concat(tmp.dir,'/',stem,'.bam')
	recalfile <- concat(tmp.dir,'/',stem,'.recal.csv')
	output.dir <- concat(config@qc.dir,'/',stem,'.',suffix)
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T CountCovariates'
	str <- concat(str,' -l INFO')
	str <- concat(str,' -R ',config@reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -B:mask,VCF ',config@maskfile)
	str <- concat(str,' --standard_covs')
	#str <- concat(str,' -cov ReadGroupCovariate')
	#str <- concat(str,' -cov QualityScoreCovariate')
	#str <- concat(str,' -cov CycleCovariate')
	#str <- concat(str,' -cov DinucCovariate')
	#str <- concat(str,' -cov HomopolymerCovariate')
	#str <- concat(str,' -cov MappingQualityCovariate')
	#str <- concat(str,' -cov MinimumNQSCovariate')
	str <- concat(str,' -recalFile ',recalfile)
	runCommand(str)
	
	checkFileExists(recalfile)
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/AnalyzeCovariates.jar'
	str <- concat(str,' -resources $GTAK_HOME/resources')
	str <- concat(str,' -recalFile ',recalfile)
	str <- concat(str,' -outputDir ',output.dir)
	runCommand(str)
}
#analyzeCovariates(config,'merged','before')

recalibrate <- function(config,stem) 
{
	tmp.dir <- config@tmp.dir
	newstem <- concat(stem,'.recal')
	bamfile <- concat(tmp.dir,'/',stem,'.bam')
	recalfile <- concat(tmp.dir,'/',newstem,'.csv')
	outfile <- concat(tmp.dir,'/',newstem,'.bam')
	
	analyzeCovariates(config,stem,'before')
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T TableRecalibration'
	str <- concat(str,' -l INFO')
	str <- concat(str,' -R ',config@reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -recalFile ',recalfile)	
	str <- concat(str,' -o ',outfile)
	runCommand(str)
	
	checkFileExists(outfile)
	
	analyzeCovariates(config,stem,'after')
	return(newstem)
}
#recalibrate(config,'merged')

outputBam <- function(config,sample,suffix='')
{
	infile <- concat(config@tmp.dir,'/',sample,suffix,'.bam')
	outfile <- concat(config@bam.dir,'/',sample,'.bam')
	runCommand('cp ',infile,' ',outfile)
	checkFileExists(outfile)
	
	infile <- concat(config@tmp.dir,'/',sample,suffix,'.bai')
	outfile <- concat(config@bam.dir,'/',sample,'.bai')
	runCommand('cp ',infile,' ',outfile)
	checkFileExists(outfile)
	#runCommand('samtools index ',outfile)
}
#outputBam(config,'merged','.realigned.recal')

outputBams <- function(config,samples=config@samples,suffix='')
{
	for (sample in samples)
	{
		outputBam(config,sample,suffix)
	}
	fixBaiFiles(config)
}
#outputBams(config)
#outputBams(config,getSamplesForSubject(config,'10464592'))

writeConsensusForBam <- function(config,sample)
{
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	bamfile <- concat(config@bam.dir,'/',sample,'.bam')
	fastqfile <- concat(config@consensus.dir,'/',sample,'.consensus.fastq')
	runCommand('samtools mpileup -uf ',reffile,' ',bamfile,' | bcftools view -cg - | vcfutils.pl vcf2fq > ',fastqfile)
	checkFileExists(fastqfile)
	fastq2fasta(fastqfile)
}
#writeConsensusForBam(config,'10464592.1__HCV-NS3-156')

writeConsensusForBams <- function(config,samples=config@samples)
{
	for (sample in samples)
	{
		writeConsensusForBam(config,sample)
	}
}
#writeConsensusForBams(config)

################################################################3

callVariants <- function(config,stem)
{	
	bamfile <- concat(config@bam.dir,'/',stem,'.bam')
	outfile <- concat(config@vcf.dir,'/',stem,'.vcf')
	
	str <- 'java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper'
	str <- concat(str,' -R ',config@reffile)
	str <- concat(str,' -I ',bamfile)
	#str <- concat(str,' -stand_call_conf 10.0')	#30.0' #50.0
	#str <- concat(str,' -stand_emit_conf 10.0')
	str <- concat(str,' -L config/',config@ref,'.interval_list')
	#str <- concat(str,' -dcov 50')
	str <- concat(str,' -o ',outfile)
	runCommand(str)
	
	checkFileExists(outfile)
}
#callVariants(config,'merged')

filterVariants <- function(config, stem)
{
	#ref <- config@ref
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	infile <- concat(config@vcf.dir,'/',stem,'.vcf')
	outfile <- concat(config@vcf.dir,'/',stem,'.filtered.vcf')
	
	str <- 'java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T VariantFiltration'
	str <- concat(str,' -R ',config@reffile)
	str <- concat(str,' -B:variant,VCF ',infile)
	str <- concat(str,' -o ',outfile)
	#str <- concat(str,' --clusterWindowSize 10')
	
	#basic indel filtering
	#str <- concat(str,' --filterExpression \'MQ0 ><- 4 && ((MQ0 / (1.0 * DP)) > 0.1)\'')  #expression to match 10% of reads with MAPQ0
	#str <- concat(str,' --filterName \'HARD_TO_VALIDATE\'')
	#str <- concat(str,' --filterExpression \'SB ><- -1.0\'')
	#str <- concat(str,' --filterName \'StrandBiasFilter\'')
	#str <- concat(str,' --filterExpression \'QUAL < 10\'')
	#str <- concat(str,' --filterName \'QualFilter\'')
	
	#snp filtering
	#str <- concat(str,' --filterExpression \'MQ0 ><- 4 && ((MQ0 / (1.0 * DP)) > 0.1)\'')
	#str <- concat(str,' --filterName \'HARD_TO_VALIDATE\'')
	#str <- concat(str,' --filterExpression \'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10\'')
	#str <- concat(str,' --filterName GATKStandard')
	
	#hard filtering
	str <- concat(str,' --filterExpression \'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10\'') #For exomes with deep coverage per sample
	#str <- concat(str,' --filterExpression \'DP > 100 || MQ0 > 40 || SB > -0.10\'') #whole genomes with deep coverage 
	str <- concat(str,' --filterName GATKStandard')
	
	runCommand(str)
	checkFileExists(outfile)
}
#filterVariants(config,'merged')

mpileupVcf <- function(config,stem)
{
	#ref <- config@ref
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	bamfile <- concat(config@bam.dir,'/',stem,'.bam')
	bcffile <- concat(config@tmp.dir,'/',stem,'.bcf')
	vcffile <- concat(config@vcf.dir,'/',stem,'.mpileup.vcf')
	
	runCommand('samtools mpileup -u -f ',config@reffile,' ',bamfile,' > ',bcffile)
	runCommand('bcftools view ',bcffile,' > ',vcffile)
	checkFileExists(vcffile)
}
#mpileupVcf(config,'merged')

unmergeBamsForRef <- function(config,stem,readgroup)
{
	infile <- concat(config@bam.dir,'/',stem,'.bam')
	outfile <- concat(config@bam.dir,'/',readgroup,'.bam')
	runCommand('samtools view -bhu -o ',outfile,' -r ',readgroup,' ',infile)
	checkFileExists(outfile)
	runCommand('samtools index ',outfile)
}
#unmergeBamsForRef(config,stem,readgroup)

unmergeBams <- function(config,stem)
{
	for (ref in rownames(config@refs))
	{
		unmergeBamsForRef(config,ref)
	}
}
#unmergeBams(config,'merged')

################################################################

fixBaiFiles <- function(config)
{
	for (filename in list.files(config@bam.dir, pattern='\\.bam$'))
	{
		stem <- stripExtension(filename)
		baifile <- concat(config@bam.dir,'/',stem,'.bam.bai')
		if (!file.exists(baifile))
		{
			oldbaifile <- concat(config@bam.dir,'/',stem,'.bai')
			runCommand('mv "',oldbaifile,'" "',baifile,'"')
		}
	}
}
#fixBaiFiles(config)

##################################################################3

filterBam <- function(config, sample)
{
	#bamtools filter -in out/bam/KT9.plasmid.bam -out out/bam/KT9.plasmid.filtered.bam -mapQuality ">30" 
	infile <- concat(config@bam.dir,'/',sample,'.bam')
	outfile <- concat(config@bam.dir,'/',sample,'.filtered.bam')
	checkFileExists(infile)
	str <- 'bamtools filter'
	str <- concat(str,' -in ',infile)
	str <- concat(str,' -out ',outfile)
	str <- concat(str,' -mapQuality "',config@map.quality,'"')
	runCommand(str)
	checkFileExists(outfile)
	runCommand('samtools index ',outfile)
}
#
#filterBams <- function(config, samples=NULL)
#{
#	require(foreach, quietly=TRUE, warn.conflicts=FALSE)
#	require(doMC, quietly=TRUE, warn.conflicts=FALSE)
#	registerDoMC()
#	if (is.null(samples))
#		samples <- config@samples
#	foreach(i=1:length(samples), .combine = cbind) %dopar% {
#		sample <- as.character(samples[i])
#		print(concat('filter_bam: ',sample))
#		try(filterBam(config,sample))
#		sample
#	}
#}
##filterBams(config)

filterBams <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		print(concat('filter_bam: ',sample))
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
	for (sample in config@samples)
	{
		countCodonsForSample(config,sample)
	}
}
#countCodons(config)

#countCodons <- function(config, groups=config@groups)
#{	
#	for (group in groups)
#	{
#		print(concat('countCodonsForGroup: ',group))
#		try(countCodonsForGroup(config, group))
#	}
#}
##countCodons(config,'confirm_with_new_reagents')

makePiecharts <- function(config)
{
	sweaveToPdf(concat(Sys.getenv("VARDB_RUTIL_HOME"),'/nextgen_piecharts.rnw'))
}
#makePiecharts(config)

concatTablesByGroup <- function(config, groups=config@groups)
{
	for (group in groups)
	{
		infiles <- concat(config@tables.dir,'/table-',group,'*.txt')
		outfile <- concat(config@tables.dir,'/group-',group,'.txt')
		#runCommand('cat ',infiles,' > ',outfile)
		runCommand('tail -n +1 ',infiles,' > ',outfile)
	}
}
#concatTablesByGroup(config)
#concatTablesByGroup(config,'PXB0219-0011')

###############################################################

analyzeReadsForSample <- function(config,sample)
{
	preprocess(config,sample)
	trimSamples(config,sample)
	mapReads(config,sample)
	filterBams(config,sample)
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
}
#analyzeReadsForGroup(config,'KT9')


#
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
