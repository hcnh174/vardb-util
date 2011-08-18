source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen')
print(getwd())
testrun <- FALSE

config <- new('nextgenconfig')

if (FALSE)
{
	testrun <- TRUE
	setCurDir('nextgen2')
	config <- new('nextgenconfig')
	#config@config.dir <- 'n:/config'
	#config <- initialize(config)	
}

run_command <- function(..., dir=NULL)
{
	command <- concat(...)
	olddir <- getwd()
	if (!is.null(dir) & testrun==FALSE)
		setwd(dir)
	print(command)
	if (command!='' & testrun==FALSE)
		system(command)
	if (!is.null(dir) & testrun==FALSE)
		setwd(olddir)
}
#run_command('ls')
#run_command('ls', dir='c:/temp')

checkFileExists <- function(filename)
{
	if (!file.exists(filename))
		throw(concat('file does not exist: ',filename))
}
#checkFileExists('test.txt')

make_folders <- function(config, subdirs='fastq,tmp,bam,unmapped,vcf,pileup,qc,counts')
{
	dir <- config@out.dir
	run_command('mkdir ',dir,' -p')
	for (subdir in splitFields(subdirs))
	{
		subdir <- concat(dir,'/',subdir)
		run_command('mkdir ',subdir,' -p; rm ',subdir,'/*')
	}
}
#make_folders(config)

preprocess <- function(config)# path='GA_RunData/110624_HWUSI-EAS1611_00063_FC639J3AAXX/Unaligned', temp.dir='tmp', fastq.dir='fastq')
{
	fastq.dir <- config@fastq.dir
	temp.dir <- config@tmp.dir
	runs <- config@runs
	run_command('rm -r ',temp.dir,'/*')
	for (sample in unique(runs$sample))
	{
		dir.to <- concat(temp.dir,'/',sample)
		run_command('mkdir ',dir.to)
	}
	run_command('')
	
	for (run in rownames(runs))
	{
		row <- runs[run,]
		sample <- row$sample
		dir.to <- concat(temp.dir,'/',sample)
		project <- row$project
		barcode <- row$barcode
		lane <- row$lane		
		dir.from <- concat(config@illumina.dir,'/Project_',project,'/Sample_',project,'/')
		filename <- concat(project,'_',barcode,'_L00',lane,'_R1_*.fastq.gz')
		run_command('cp ',dir.from,filename,' ',dir.to)
	}
	run_command('')
	
	for (sample in unique(runs$sample))
	{
		dir.from <- concat(temp.dir,'/',sample)
		run_command('gunzip ',dir.from,'/*')
	}
	run_command('')
	
	for (sample in unique(runs$sample))
	{
		dir.from <- concat(temp.dir,'/',sample)
		fastqfile <- concat(fastq.dir,'/',sample,'.fastq')
		run_command('cat ',dir.from,'/* > ',fastqfile)
		checkFileExists(fastqfile)
	}
	run_command('')
	run_command('rm -r ',temp.dir)
}
#preprocess(config)

#qseq2fastq <- function(sample)
#{
#	run_command('perl qseq2fastq.pl -a qseq/',sample,'.qseq -v T')
#	run_command('mv ',sample,'.fastq fastq/',sample,'.fastq')
#}

solexa_qa <- function(config,sample)
{
	run_command('cd ',config@qc.dir,'; SolexaQA.pl ../fastq/',sample,'.fastq')
}
#solexa_qa(config,'KT9.plasmid')

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
#trim(config,'PXB0220-0002.wk13')
#cd out/fastq; DynamicTrim.pl PXB0220-0002.wk13.fastq
#LengthSort.pl PXB0220-0002.wk13.fastq.trimmed

trim_prinseq <- function(config,sample)
{
	fastq.dir <- config@fastq.dir
	infile <- concat(fastq.dir,'/',sample,'.fastq')
	outfile <- concat(fastq.dir,'/',sample,'.trimmed') #program adds .fastq extension automatically
	str <- 'prinseq-lite.pl'
	str <- concat(str,' -fastq ',infile)
	str <- concat(str,' -out_good ',outfile)
	str <- concat(str,' -out_bad null')
	str <- concat(str,' -min_len 36')
	str <- concat(str,' -min_qual_score 30')
	str <- concat(str,' -ns_max_n 0')
	#str <- concat(str,' -stats_all')
	run_command(str)
	checkFileExists(concat(outfile,'.fastq'))
}

trim_all <- function(config)
{
	for (sample in rownames(config@samples))
	{
		trim_solexaqa(config,sample)
		#trim_prinseq(config,sample)
	}
}

remove_exact_duplicates <- function(config,sample)
{
	fastq.dir <- config@fastq.dir
	infile <- concat(fastq.dir,'/',sample,'.trimmed.fastq')
	outfile <- concat(fastq.dir,'/',sample,'.trimmed.dedup') #program adds .fastq extension automatically
	run_command('prinseq-lite.pl -fastq ',infile,' -out_good ',outfile,' -out_bad null -derep 1')#4
	checkFileExists(concat(outfile,'.fastq'))
}
#remove_exact_duplicates('PXB0220-0002.wk12')

#testrun <- TRUE
remove_exact_duplicates_for_all_samples <- function(config)
{
	for (sample in rownames(config@samples))
	{
		remove_exact_duplicates(config,sample)
	}
}
#remove_exact_duplicates_for_all_samples(config)

######################################################################

run_bwa <- function(config, sample, fastq.ext='.trimmed.fastq', q=20)#, .trimmed.fastq
{
	#ref <- config@ref
	stem <- concat(sample,'.',config@ref)
	fqfile <- concat(config@fastq.dir,'/',sample,fastq.ext)
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	tmp.dir <- config@tmp.dir

	saifile <- concat(tmp.dir,'/',stem,'.sai')
	samfile <- concat(tmp.dir,'/',stem,'.sam')
	
	run_command('bwa aln -q ',q,' ',config@reffile,' ',fqfile,' > ',saifile)
	run_command('bwa samse ',config@reffile,' ',saifile,' ',fqfile,' > ',samfile)
	checkFileExists(samfile)
}
#run_bwa(config,'PXB0220-0002.wk13')

run_bowtie <- function(config, sample, fastq.ext='.trimmed.fastq')
{
	#ref <- config@ref
	stem <- concat(sample,'.',config@ref)
	fqfile <- concat(config@fastq.dir,'/',sample,fastq.ext)
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	tmp.dir <- config@tmp.dir
	
	samfile <- concat(tmp.dir,'/',stem,'.sam')
	#unmappedfile <- concat(tmp.dir,'/unmapped/',stem,'.un.txt')
	run_command('bowtie -S ',config@reffile,' ',fqfile,' ',samfile) #--un unmappedfile
	checkFileExists(samfile)
}
#run_bowtie(config,'PXB0220-0002.wk13')

add_read_groups <- function(config,sample)
{
	#ref <- config@ref
	stem <- concat(sample,'.',config@ref)
	newstem <- concat(stem,'.rg')
	tmp.dir <- config@tmp.dir
	infile <- concat(tmp.dir,'/',stem,'.sam')
	outfile <- concat(tmp.dir,'/',newstem,'.bam')

	str <- 'java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar'
	str <- concat(str,' INPUT=',infile)
	str <- concat(str,' OUTPUT=',outfile)
	str <- concat(str,' RGSM="',sample,'"')
	str <- concat(str,' RGLB="',sample,'"')
	str <- concat(str,' RGID="',sample,'"')
	str <- concat(str,' RGPL=illumina')
	str <- concat(str,' RGPU=barcode')
	str <- concat(str,' SORT_ORDER=coordinate')
	str <- concat(str,' CREATE_INDEX=true')	
	run_command(str)
	return(newstem)
}
#add_read_groups(config,'PXB0220-0002.wk13')

#mark_duplicates <- function(config,stem) #sample,'.',ref,'.rg'
#{
#	newstem <- concat(stem,'.dedup')
#	tmp.dir <- concat(config@out.dir,'/tmp')
#	metricsfile <- concat(tmp.dir,'/qc/',newstem,'.metrics')
#	infile <- concat(tmp.dir,'/',stem,'.bam')
#	outfile <-  concat(tmp.dir,'/',newstem,'.bam')
#	
#	str <- 'java -Xmx2g -jar $PICARD_HOME/MarkDuplicates.jar'
#	str <- concat(str,' INPUT=',infile)
#	str <- concat(str,' OUTPUT=',outfile)
#	str <- concat(str,' METRICS_FILE=',metricsfile)
#	run_command(str)
#	run_command('samtools index ',outfile)
#	return(newstem)
#}
#mark_duplicates(config,'PXB0220-0002.wk13.KT9.rg')

map_reads_for_sample <- function(config, sample)
{
	stem <- run_bwa(config,sample)
	#stem <- run_bowtie(config,sample)
	stem.rg <- add_read_groups(config,sample)
	#stem.rg.dedup <- mark_duplicates(config,stem.rg)
}
#map_reads_for_sample(config,'PXB0220-0002.wk13')

map_reads_for_all_samples <- function(config)
{
	for (sample in rownames(config@samples))
	{
		map_reads_for_sample(config,sample)
	}
}
#map_reads_for_all_samples(config)

merge_bams <- function(config)
{
	#ref <- config@ref
	tmp.dir <- config@tmp.dir
	outfile <- concat(tmp.dir,'/merged.bam')
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
	str <- concat(str,' CREATE_INDEX=true')
	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
	for (sample in rownames(config@samples))
	{
		infile <- concat(tmp.dir,'/',sample,'.',config@ref,'.rg.bam') #.dedup
		str <- concat(str,' INPUT=',infile)
	}
	str <- concat(str,' OUTPUT=',outfile)
	run_command(str)
}


show_coverage <- function(config,stem)
{
	#tmp.dir <- config@tmp.dir
	bamfile <- concat(config@tmp.dir,'/',stem,'.bam')
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T DepthOfCoverage'
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -R ',config@reffile)
	run_command(str)
	#run_command('java -jar GenomeAnalysisTK.jar -I aln.bam -R hsRef.fa -T DepthOfCoverage -L intervals.txt -U -S SILENT')
}

###################################################################

realign_indels <- function(config, stem)
{
	#ref <- config@ref
	tmp.dir <- config@tmp.dir
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	bamfile <- concat(tmp.dir,'/',stem,'.bam')
	intervalfile <- concat(tmp.dir,'/',stem,'.intervals')
	outfile <- concat(tmp.dir,'/',stem,'.realigned.bam')
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator'
	str <- concat(str,' -R ',config@reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -o ',intervalfile)
	run_command(str)
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T IndelRealigner'
	str <- concat(str,' -R ',config@reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -targetIntervals ',intervalfile)
	str <- concat(str,' -o ',outfile)
	run_command(str)
	
	checkFileExists(outfile)
}
#realign_indels(config,'merged')

analyze_covariates <- function(config,stem,suffix)
{
	#ref <- config@ref
	tmp.dir <- config@tmp.dir
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	bamfile <- concat(tmp.dir,'/',stem,'.bam')
	#maskfile <- concat(config@config.dir,'/',ref,'.mask.vcf')
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
	#ref <- config@ref
	tmp.dir <- config@tmp.dir
	newstem <- concat(stem,'.recal')
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	#maskfile <- concat(config@config.dir,'/',ref,'.mask.vcf')
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

output_bam <- function(config,stem,suffix)
{
	infile <- concat(config@tmp.dir,'/',stem,'.',suffix,'.bam')
	outfile <- concat(config@bam.dir,'/',stem,'.bam')
	run_command('cp ',infile,' ',outfile)
	checkFileExists(outfile)
	run_command('samtools index ',outfile)
}
#output_bam(config,'merged','realigned.recal')

################################################################3

call_variants <- function(config,stem)
{	
	#ref <- config@ref
	#reffile <- concat(config@ref.dir,'/',ref,'.fasta')
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

export_read_group <- function(config,stem,readgroup)
{
	infile <- concat(config@bam.dir,'/',stem,'.bam')
	outfile <- concat(config@bam.dir,'/',readgroup,'.bam')
	run_command('samtools view -bhu -o ',outfile,' -r ',readgroup,' ',infile)
	checkFileExists(outfile)
	run_command('samtools index ',outfile)
}
#export_read_group('merged','KT9.plasmid')

export_pileup <- function(config,sample)
{
	print(concat('ref: ',config@ref))
	run_command('python $VARDB_RUTIL_HOME/export_pileup.py ', sample,' ',config@ref,' ',config@bam.dir,' ',config@pileup.dir)
}
#export_pileup(config,'merged','KT9')

filter_bam <- function(config, stem)
{
	#bamtools filter -in out/bam/KT9.plasmid.bam -out out/bam/KT9.plasmid.filtered.bam -mapQuality ">30" 
	infile <- concat(config@bam.dir,'/',stem,'.bam')
	outfile <- concat(config@bam.dir,'/',stem,'.filtered.bam')
	checkFileExists(infile)
	str <- 'bamtools filter'
	str <- concat(str,' -in ',infile)
	str <- concat(str,' -out ',outfile)
	str <- concat(str,' -mapQuality ">30"')
	run_command(str)
	checkFileExists(outfile)
	run_command('samtools index ',outfile)
}

export_read_groups <- function(config,stem)
{
	for (sample in rownames(config@samples))
	{
		#export_read_group(config,stem,sample)
		filter_bam(config,sample)
		##remove_duplicates(sample,ref)
		export_pileup(config,concat(sample,'.filtered'))
	}
}
#export_read_groups(config,'merged')

##################################################################3
#
#export_unmapped_reads <- function(stem)
#{
#	bamfile <- concat('tmp/',stem,'.bam')
#	fastqfile <- concat('unmapped/',stem,'.unmapped.fastq')
#	trimmedfastqfile <- concat('unmapped/',stem,'.unmapped.fastq.trimmed')
#	fastafile <- concat('unmapped/',stem,'.unmapped.fasta')
#	#run_command('bam2fastq --force --no-aligned --unaligned --no-filtered -o ',fastqfile,' ',bamfile)
#	#run_command('cd unmapped; DynamicTrim.pl ',concat(stem,'.unmapped.fastq'))
#	run_command('perl $VARDB_RUTIL_HOME/fq_all2std.pl fq2fa ',trimmedfastqfile,' > ',fastafile)
#	run_command('perl $VARDB_RUTIL_HOME/fastq2table.pl -a ',trimmedfastqfile)
#	run_command('sort merged.nodup.unmapped.fastq.trimmed.table.txt | uniq -dc > unmapped/merged.unique.txt')
#}
##export_unmapped_reads('bam/PXB0220-0002.wk11.bam')
#
#remove_duplicates <- function(stem,ref)
#{
#	metricsfile <- concat('tmp/',stem,'.nodup.metrics')
#	outfile <-  concat('tmp/',stem,'.nodup.bam')
#	
#	str <- 'java -Xmx2g -jar $PICARD_HOME/MarkDuplicates.jar'
#	str <- concat(str,' INPUT=bam/',stem,'.bam')
#	str <- concat(str,' OUTPUT=',outfile)
#	str <- concat(str,' METRICS_FILE=',metricsfile)
#	str <- concat(str,' REMOVE_DUPLICATES=true')
#	run_command(str)
#	run_command('samtools index ',outfile)
#}

count_codons <- function(config)
{
	for (subject in config@subjects$subject)
	{
		count_codons_for_subject(config, subject)
		#run_command('Rscript $VARDB_RUTIL_HOME/count_codons.r subject=',subject,' countdir=',countdir)
	}
}

make_tables <- function(config)
{
	sweaveToPdf(concat(Sys.getenv("VARDB_RUTIL_HOME"),'/tables.rnw'))
	#run_command('Rscript $VARDB_RUTIL_HOME/make_tables.r dir=',config@count.dir)
}

analyze_reads_merged <- function(config)
{
	stem <- 'merged'
	#make_folders(config)
	#preprocess(config)
	#trim_all(config)
	##remove_exact_duplicates_for_all_samples(config)
	#map_reads_for_all_samples(config)
	#merge_bams(config)
	#realign_indels(config,stem)
	#recalibrate(config,stem) #concat(stem,'.realigned'))
	#output_bam(config,stem,'recal') #'realigned.recal')
	
	#call_variants(config,stem)
	#filter_variants(config,stem)
	#export_read_groups(config,stem)
	#count_codons(config)
	make_tables(config)
	#export_unmapped_reads(config,concat(stem,'.nodup'))
}

args <- commandArgs(TRUE) # from R.utils package
config@ref <- args$ref
config@out.dir <- args$out

analyze_reads_merged(config)
#test2
#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r ref=KT9 out=out
#Rscript ~/workspace/vardb-util/src/main/r/analyze_reads_merged.r ref=KT9 out=out

