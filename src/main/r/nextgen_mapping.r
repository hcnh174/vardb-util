preprocess <- function(config, subdirs='tmp,ref,fastq,tmp,bam,unmapped,vcf,pileup,qc,counts,tables')
{
	run_command('mkdir ',config@out.dir,' -p')
	for (subdir in splitFields(subdirs))
	{
		subdir <- concat(config@out.dir,'/',subdir)
		run_command('mkdir ',subdir,' -p; rm ',subdir,'/*')
	}
	writeRefs(config)
	
	fastq.dir <- config@fastq.dir
	temp.dir <- config@tmp.dir
	#runs <- config@runs
	run_command('rm -r ',temp.dir,'/*')
	for (sample in config@samples)
	{
		dir.to <- concat(temp.dir,'/',sample)
		run_command('mkdir ',dir.to)
	}
	run_command('')
	
	for (run in rownames(config@runs))
	{
		row <-  config@runs[run,]
		sample <- row$sample
		dir.to <- concat(temp.dir,'/',sample)
		folder <- row$folder
		barcode <- row$barcode
		lane <- row$lane
		dir.from <- concat(config@illumina.dir,'/Project_',folder,'/Sample_',folder,'/')
		filename <- concat(folder,'_',barcode,'_L00',lane,'_R1_*.fastq.gz')
		run_command('cp ', dir.from, filename,' ',dir.to)
		filename <- concat(folder,'_',barcode,'_L00',lane,'_R1_*.fastq')
		run_command('cp ', dir.from, filename,' ',dir.to)
	}	
	run_command('')
	
	for (sample in config@samples)
	{
		dir.from <- concat(temp.dir,'/',sample)
		run_command('gunzip ',dir.from,'/*')
	}
	run_command('')
	
	for (sample in config@samples)
	{
		dir.from <- concat(temp.dir,'/',sample)
		fastqfile <- concat(fastq.dir,'/',sample,'.fastq')
		run_command('cat ',dir.from,'/* > ',fastqfile)
		checkFileExists(fastqfile)
	}
	run_command('')
	run_command('rm -r ',temp.dir,'/*')
}
#preprocess(config)

#############################################################

trim_solexaqa <- function(config,sample)
{
	olddir <- getwd()
	setwd(config@fastq.dir)
	run_command('DynamicTrim.pl ',sample,'.fastq',' -phredcutoff 30')
	run_command('LengthSort.pl ',sample,'.fastq.trimmed',' -length 36')
	checkFileExists(concat(sample,'.fastq.trimmed.single'))
	run_command('mv ',sample,'.fastq.trimmed.single',' ',sample,'.trimmed.fastq')
	run_command('rm ',sample,'.fastq.trimmed')
	run_command('rm ',sample,'.fastq.trimmed.discard')
	setwd(olddir)
}
#trim_solexaqa(config,'PXB0220-0002.wk13')
#cd out/fastq; DynamicTrim.pl PXB0220-0002.wk13.fastq
#LengthSort.pl PXB0220-0002.wk13.fastq.trimmed

trim_all <- function(config)
{
	#if (!config@trim)
	#	return()
	for (sample in config@samples)
	{
		trim_solexaqa(config,sample)
		#trim_prinseq(config,sample)
	}
}
#trim_all(config)

#trim_all <- function(config)
#{
#	if (!config@trim)
#		return()
#	require(foreach, quietly=TRUE, warn.conflicts=FALSE)
#	require(doMC, quietly=TRUE, warn.conflicts=FALSE)
#	registerDoMC()	
#	foreach(i=1:length(config@samples), .combine = cbind) %dopar% {
#		sample <- as.character(config@samples[i])
#		print(concat('trim_solexaqa: ',sample))
#		try(trim_solexaqa(config,sample))
#		sample
#	}
#}
#trim_all(config)

######################################################################

add_read_groups <- function(config, sample)
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
	run_command(str)
	checkFileExists(bamfile)
	baifile <- concat(tmp.dir,'/',sample,'.bai')
	checkFileExists(baifile)
}
#add_read_groups(config,'110617HBV.HBV07@HBV-RT')

run_bwa <- function(config, sample)#, .trimmed.fastq #, q=20
{
	fastq.ext <- ifelse(config@trim,'.trimmed.fastq','.fastq')
	ref <- get_ref_for_sample(sample)
	reffile <- get_reffile(config,ref)
	fqfile <- concat(config@fastq.dir,'/',sample,fastq.ext)
	tmp.dir <- config@tmp.dir
	
	samfile <- concat(tmp.dir,'/',sample,'.sam')
	saifile <- concat(tmp.dir,'/',sample,'.sai')
	
	run_command('bwa index ',reffile)
	run_command('bwa aln ',reffile,' ',fqfile,' > ',saifile)
	run_command('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
	checkFileExists(samfile)
	
	add_read_groups(config,sample)
	
	#run_command('rm ',samfile)
	#run_command('rm ',saifile)
}
#run_bwa(config,'110617HBV-1.10348001.20020530__HBV-RT')
#run_bwa(config,'110617HBV-1.10348001.20020624__HBV-RT')
#run_bwa(config,'110617HBV-1.10348001.20040128__HBV-RT')
#run_bwa(config,'110617HBV-1.10348001.20040315__HBV-RT')
#run_bwa(config,'110617HBV-1.10348001.20040804__HBV-RT')
#run_bwa(config,'110628-1.PXB0220-0030.13__HCV-NS5A-31')

map_reads <- function(config)
{
	#require(foreach, quietly=TRUE, warn.conflicts=FALSE)
	#require(doMC, quietly=TRUE, warn.conflicts=FALSE)
	#registerDoMC()
	#foreach(i=1:length(config@samples), .combine = cbind) %do% { #dopar
	#	sample <- as.character(config@samples[i])
	for (sample in config@samples)
	{
		print(concat('run_bwa: ',sample))
		try(run_bwa(config,sample))
		#sample
	}
}
#map_reads(config)

############################################

merge_bams_for_ref <- function(config,ref)
{
	tmp.dir <- config@tmp.dir
	outfile <- concat(tmp.dir,'/merged@',ref,'.bam')
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
	str <- concat(str,' CREATE_INDEX=true')
	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
	
	for (sample in unique(config@runs[which(config@runs$ref==ref),'sample']))
	{
		ref <- get_ref_for_sample(sample)
		infile <- concat(tmp.dir,'/',sample,'.bam')
		str <- concat(str,' INPUT=',infile)
	}
	str <- concat(str,' OUTPUT=',outfile)
	run_command(str)
}
#merge_bams(config)

merge_bams <- function(config)
{
	for (ref in rownames(config@refs))
	{
		merge_bams_for_ref(config,ref)
	}
}
#merge_bams(config)

###################################################################

analyze_covariates <- function(config,stem,suffix)
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
	run_command(str)
	
	checkFileExists(recalfile)
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/AnalyzeCovariates.jar'
	str <- concat(str,' -resources $GTAK_HOME/resources')
	str <- concat(str,' -recalFile ',recalfile)
	str <- concat(str,' -outputDir ',output.dir)
	run_command(str)
}
#analyze_covariates(config,'merged','before')

recalibrate <- function(config,stem) 
{
	tmp.dir <- config@tmp.dir
	newstem <- concat(stem,'.recal')
	bamfile <- concat(tmp.dir,'/',stem,'.bam')
	recalfile <- concat(tmp.dir,'/',newstem,'.csv')
	outfile <- concat(tmp.dir,'/',newstem,'.bam')
	
	analyze_covariates(config,stem,'before')
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T TableRecalibration'
	str <- concat(str,' -l INFO')
	str <- concat(str,' -R ',config@reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -recalFile ',recalfile)	
	str <- concat(str,' -o ',outfile)
	run_command(str)
	
	checkFileExists(outfile)
	
	analyze_covariates(config,stem,'after')
	return(newstem)
}
#recalibrate(config,'merged')

output_bam <- function(config,sample,suffix='')
{
	infile <- concat(config@tmp.dir,'/',sample,suffix,'.bam')
	outfile <- concat(config@bam.dir,'/',sample,'.bam')
	run_command('cp ',infile,' ',outfile)
	checkFileExists(outfile)
	
	infile <- concat(config@tmp.dir,'/',sample,suffix,'.bai')
	outfile <- concat(config@bam.dir,'/',sample,'.bai')
	run_command('cp ',infile,' ',outfile)
	checkFileExists(outfile)
	#run_command('samtools index ',outfile)
}
#output_bam(config,'merged','.realigned.recal')

output_bams <- function(config,stem,suffix='')
{
	for (sample in config@samples)
	{
		output_bam(config,sample,suffix)
	}
}
#output_bams(config)

################################################################3

call_variants <- function(config,stem)
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
	run_command(str)
	
	checkFileExists(outfile)
}
#call_variants(config,'merged')

filter_variants <- function(config, stem)
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
	
	run_command(str)
	checkFileExists(outfile)
}
#filter_variants(config,'merged')

mpileup_vcf <- function(config,stem)
{
	#ref <- config@ref
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	bamfile <- concat(config@bam.dir,'/',stem,'.bam')
	bcffile <- concat(config@tmp.dir,'/',stem,'.bcf')
	vcffile <- concat(config@vcf.dir,'/',stem,'.mpileup.vcf')
	
	run_command('samtools mpileup -u -f ',config@reffile,' ',bamfile,' > ',bcffile)
	run_command('bcftools view ',bcffile,' > ',vcffile)
	checkFileExists(vcffile)
}
#mpileup_vcf(config,'merged')


#export_read_group('merged','KT9.plasmid')

unmerge_bams_for_ref <- function(config,stem,readgroup)
{
	infile <- concat(config@bam.dir,'/',stem,'.bam')
	outfile <- concat(config@bam.dir,'/',readgroup,'.bam')
	run_command('samtools view -bhu -o ',outfile,' -r ',readgroup,' ',infile)
	checkFileExists(outfile)
	run_command('samtools index ',outfile)
}
#merge_bams(config)

unmerge_bams <- function(config,stem)
{
	for (ref in rownames(config@refs))
	{
		unmerge_bams_for_ref(config,ref)
	}
	
#	for (sample in config@samples)
#	{
#		#export_read_group(config,stem,sample)
#		ref <- get_ref_for_sample(sample)
#		sample <- concat(sample,'.',ref,'.rg') #hack!
#		filter_bam(config,sample)
#		##remove_duplicates(sample,ref)
#		export_pileup(config,concat(sample,'.filtered'),ref)
#	}
}
#unmerge_bams(config,'merged')

################################################################

fix_bai_files <- function(config)
{
	for (filename in list.files(config@bam.dir, pattern='\\.bam$'))
	{
		stem <- stripExtension(filename)
		baifile <- concat(config@bam.dir,'/',stem,'.bam.bai')
		if (!file.exists(baifile))
		{
			oldbaifile <- concat(config@bam.dir,'/',stem,'.bai')
			run_command('mv "',oldbaifile,'" "',baifile,'"')
		}
	}
}
#fix_bai_files(config)

##################################################################3

filter_bam <- function(config, sample)
{
	#bamtools filter -in out/bam/KT9.plasmid.bam -out out/bam/KT9.plasmid.filtered.bam -mapQuality ">30" 
	infile <- concat(config@bam.dir,'/',sample,'.bam')
	outfile <- concat(config@bam.dir,'/',sample,'.filtered.bam')
	checkFileExists(infile)
	str <- 'bamtools filter'
	str <- concat(str,' -in ',infile)
	str <- concat(str,' -out ',outfile)
	str <- concat(str,' -mapQuality ">30"')
	run_command(str)
	checkFileExists(outfile)
	run_command('samtools index ',outfile)
}

#filter_bams <- function(config)
#{
#	for (sample in config@samples)
#	{
#		filter_bam(config,sample)
#	}
#}

filter_bams <- function(config)
{
	require(foreach, quietly=TRUE, warn.conflicts=FALSE)
	require(doMC, quietly=TRUE, warn.conflicts=FALSE)
	registerDoMC()
	foreach(i=1:length(config@samples), .combine = cbind) %dopar% {
		sample <- as.character(config@samples[i])
		print(concat('filter_bam: ',sample))
		try(filter_bam(config,sample))
		sample
	}
}
#filter_bams(config)

####################################################

export_pileup_for_sample <- function(config,sample,filtered=TRUE)
{
	ref <- get_ref_for_sample(sample)
	samplename <- ifelse(filtered,concat(sample,'.filtered'), sample)
	run_command('python $VARDB_RUTIL_HOME/export_pileup.py ', samplename,' ',ref,' ',config@bam.dir,' ',config@pileup.dir)
}
#export_pileup_for_sample(config,'110617HBV-1.10348001.20020530__HBV-RT',filtered=FALSE)

export_pileup <- function(config)
{
	for (sample in config@samples)
	{
		export_pileup_for_sample(config,sample,filtered=FALSE)
		export_pileup_for_sample(config,sample,filtered=TRUE)
	}
}
#export_pileup(config)

#export_pileup <- function(config)
#{
#	require(foreach, quietly=TRUE, warn.conflicts=FALSE)
#	require(doMC, quietly=TRUE, warn.conflicts=FALSE)
#	registerDoMC() #cores=4
#	#print(concat('workers: ',getDoParWorkers(),' (',getDoParName(),' ',getDoParVersion(),')'))
#	foreach(i=1:length(config@samples), .combine = cbind) %dopar% {#dopar
#		sample <- as.character(config@samples[i])
#		print(concat('export_pileup: ',sample))
#		#try({
#			export_pileup_for_sample(config,sample,filtered=FALSE)
#			export_pileup_for_sample(config,sample,filtered=TRUE)
#		#}, silent=FALSE)
#		sample
#	}
#}

################################################

count_codons <- function(config)
{
	require(foreach, quietly=TRUE, warn.conflicts=FALSE)
	require(doMC, quietly=TRUE, warn.conflicts=FALSE)
	registerDoMC() #cores=4
	#print(concat('workers: ',getDoParWorkers(),' (',getDoParName(),' ',getDoParVersion(),')'))
	foreach(i=1:length(config@subjects), .combine = cbind) %dopar% {
		subject <- as.character(config@subjects[i])
		print(concat('count_codons_for_subject: ',subject))
		try(count_codons_for_subject(config, subject))
		subject
	}
}
#count_codons(config)

#make_tables <- function(config)
#{
#	sweaveToPdf(concat(Sys.getenv("VARDB_RUTIL_HOME"),'/tables.rnw'))
#	#run_command('Rscript $VARDB_RUTIL_HOME/make_tables.r dir=',config@count.dir)
#}

make_piecharts <- function(config)
{
	sweaveToPdf(concat(Sys.getenv("VARDB_RUTIL_HOME"),'/nextgen_piecharts.rnw'))
}
#make_piecharts(config)
