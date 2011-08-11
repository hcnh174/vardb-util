source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'common.r',sep=''))
loadUtilFiles('nextgen')
#setCurDir('nextgen2')
print(getwd())

config <- new('nextgenconfig')

run_command <- function(...)
{
	command <- concat(...)
	print(command)	
	if (command!='') system(command)
}

preprocess <- function(config, path='GA_RunData/110624_HWUSI-EAS1611_00063_FC639J3AAXX/Unaligned',
		temp.dir='tmp', fastq.dir='fastq')
{
	runs <- config@runs
	run_command('rm -r tmp/*')	
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
		dir.from <- concat(path,'/Project_',project,'/Sample_',project,'/')
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
		run_command('cat ',dir.from,'/* > ',fastq.dir,'/',sample,'.fastq')
		
	}
	run_command('')
	#run_command('rm -r ',temp.dir,'/*')
}

qseq2fastq <- function(sample)
{
	run_command('perl qseq2fastq.pl -a qseq/',sample,'.qseq -v T')
	run_command('mv ',sample,'.fastq fastq/',sample,'.fastq')
}

trim <- function(sample)
{
	run_command('cd trimmed; DynamicTrim.pl ../fastq/',sample,'.fastq')
}

solexa_qa <- function(sample)
{
	run_command('cd quality; SolexaQA.pl ../fastq/',sample,'.fastq')
}

remove_exact_duplicates <- function(sample)
{
	infile <- concat('fastq/',sample,'.fastq')
	outfile <- concat('fastq/',sample,'.dedup') #program adds .fastq extension automatically
	run_command('prinseq-lite.pl -fastq ',infile,' -out_good ',outfile,' -derep 14')
}
#remove_exact_duplicates('PXB0220-0002.wk12')

remove_exact_duplicates_for_all_samples <- function(config)
{
	for (sample in rownames(config@samples))
	{
		remove_exact_duplicates(sample)
	}
}

run_bwa <- function(sample,ref)
{
	stem <- concat(sample,'.',ref)
	fqfile <- concat('fastq/',sample,'.fastq')
	reffile <- concat('ref/',ref,'.fasta')
	saifile <- concat('tmp/',stem,'.sai')
	samfile <- concat('tmp/',stem,'.sam')
	bamfile <- concat('tmp/',stem,'.bam')
	
	run_command('bwa aln ',reffile,' ',fqfile,' > ',saifile)
	run_command('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
	run_command('samtools view -bS -o ',bamfile,' ',samfile)
	run_command('samtools sort ',bamfile,' tmp/',stem)
	run_command('samtools index ',bamfile)
	return(stem)
}
#run_bwa('PXB0220-0002.wk13.lane8','KT9')

add_read_groups <- function(sample,ref)
{
	stem <- concat(sample,'.',ref)
	newstem <- concat(stem,'.rg')
	str <- 'java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar'
	str <- concat(str,' INPUT=tmp/',stem,'.bam')
	str <- concat(str,' OUTPUT=tmp/',newstem,'.bam')
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
#add_read_groups('PXB0220-0002.wk13','8','KT9')

mark_duplicates <- function(stem) #sample,'.',ref,'.rg'
{
	newstem <- concat(stem,'.dedup')
	metricsfile <- concat('qc/',newstem,'.metrics')
	outfile <-  concat('tmp/',newstem,'.bam')
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MarkDuplicates.jar'
	str <- concat(str,' INPUT=tmp/',stem,'.bam')
	str <- concat(str,' OUTPUT=',outfile)
	str <- concat(str,' METRICS_FILE=',metricsfile)
	run_command(str)
	run_command('samtools index ',outfile)
	return(newstem)
}
#mark_duplicates('PXB0220-0002.wk13.KT9.rg')

map_reads_for_sample <- function(sample, ref)
{
	stem <- run_bwa(sample,ref)
	stem.rg <- add_read_groups(sample, ref)
	stem.rg.dedup <- mark_duplicates(stem.rg)
	return(stem.rg.dedup)
}
#map_reads_for_sample('PXB0220-0002.wk13','KT9')

extract_read_group <- function(infile,readgroup)
{
	outfile <- concat('bam/',readgroup,'.bam')
	print(concat('sample=',infile))
	run_command('samtools view -hu -o ',outfile,' -r ',readgroup,' ',infile) #b
	run_command('samtools index ',outfile)
}
#extract_read_group('tmp/KT9.plasmid.lane1.KT9.rg.bam','KT9.plasmid')

realign_indels <- function(stem,ref)
{
	reffile <- concat('ref/',ref,'.fasta')
	bamfile <- concat('tmp/',stem,'.bam')
	intervalfile <- concat('tmp/',stem,'.intervals')
	outfile <- concat('tmp/',stem,'.realigned.bam')
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator'
	str <- concat(str,' -R ',reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -o ',intervalfile)
	run_command(str)
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T IndelRealigner'
	str <- concat(str,' -R ',reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -targetIntervals ',intervalfile)
	str <- concat(str,' -o ',outfile)
	run_command(str)
}
#realign_indels('merged','KT9')

analyze_covariates <- function(stem,ref,suffix)
{
	reffile <- concat('ref/',ref,'.fasta')
	bamfile <- concat('tmp/',stem,'.bam')
	maskfile <- concat('config/',ref,'.mask.vcf')
	recalfile <- concat('tmp/',stem,'.recal.csv')
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T CountCovariates'
	str <- concat(str,' -l INFO')
	str <- concat(str,' -R ',reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -B:mask,VCF ',maskfile)
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
	
	str <- 'java -Xmx2g -jar $GTAK_HOME/AnalyzeCovariates.jar'
	str <- concat(str,' -resources $GTAK_HOME/resources')
	str <- concat(str,' -recalFile ',recalfile)
	str <- concat(str,' -outputDir ./qc/',stem,'.',suffix)
	run_command(str)
}

recalibrate <- function(stem,ref) 
{
	newstem <- concat(stem,'.recal')
	reffile <- concat('ref/',ref,'.fasta')
	bamfile <- concat('tmp/',stem,'.bam')
	maskfile <- concat('config/',ref,'.mask.vcf')
	recalfile <- concat('tmp/',newstem,'.csv')
	outfile <- concat('tmp/',newstem,'.bam')
	
	#analyze_covariates(stem,ref,'before')

	str <- 'java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T TableRecalibration'
	str <- concat(str,' -l INFO')
	str <- concat(str,' -R ',reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' -recalFile ',recalfile)	
	str <- concat(str,' -o ',outfile)
	#run_command(str)

	analyze_covariates(newstem,ref,'after')
	return(newstem)
}
#recalibrate('merged','KT9')

call_variants <- function(stem,ref)
{	
	reffile <- concat('ref/',ref,'.fasta')
	bamfile <- concat('bam/',stem,'.bam')
	outfile <- concat('vcf/',stem,'.vcf')
	
	str <- 'java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper'
	str <- concat(str,' -R ',reffile)
	str <- concat(str,' -I ',bamfile)
	#str <- concat(str,' -stand_call_conf 10.0')	#30.0' #50.0
	#str <- concat(str,' -stand_emit_conf 10.0')
	str <- concat(str,' -L config/',ref,'.interval_list')
	#str <- concat(str,' -dcov 50')
	str <- concat(str,' -o ',outfile)
	run_command(str)
}
#call_variants('merged','KT9')

filter_variants <- function(stem, ref)
{
	reffile <- concat('ref/',ref,'.fasta')
	infile <- concat('vcf/',stem,'.vcf')
	outfile <- concat('vcf/',stem,'.filtered.vcf')
	
	str <- 'java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T VariantFiltration'
	str <- concat(str,' -R ',reffile)
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
}
#filter_variants('merged','KT9')

mpileup_vcf <- function(stem,ref)
{
	reffile <- concat('ref/',ref,'.fasta')
	run_command('samtools mpileup -u -f ',reffile,' bam/',stem,'.bam > tmp/',stem,'.bcf')
	run_command('bcftools view tmp/',stem,'.bcf > vcf/',stem,'.mpileup.vcf')
}
#mpileup_vcf('merged','KT9')

export_pileup <- function(sample,ref)
{
	run_command('python $VARDB_RUTIL_HOME/export_pileup.py ',sample,' ',ref)
}
#export_pileup('merged','KT9')

output_bam <- function(stem,suffix)
{
	infile <- concat('tmp/',stem,'.',suffix,'.bam')
	outfile <- concat('bam/',stem,'.bam')
	run_command('cp ',infile,' ',outfile)
	run_command('samtools index ',outfile)
}
#output_bam('merged','realigned.recal')

###########################################################

map_reads_for_all_samples <- function(config)
{
	for (sample in rownames(config@samples))
	{
		ref <- config@samples[sample,'ref']
		map_reads_for_sample(sample,ref)
	}
}

merge_bams <- function(config)
{
	outfile <- 'tmp/merged.bam'
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
	str <- concat(str,' CREATE_INDEX=true')
	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
	for (sample in rownames(config@samples))
	{
		ref <- config@samples[sample,'ref']
		str <- concat(str,' INPUT=tmp/',sample,'.',ref,'.rg.dedup.bam')
	}
	str <- concat(str,' OUTPUT=',outfile)
	run_command(str)
	#run_command('samtools index ',outfile)
}
#merge_bams(config)

export_read_group <- function(stem,readgroup)
{
	infile <- concat('bam/',stem,'.bam')
	outfile <- concat('bam/',readgroup,'.bam')
	run_command('samtools view -bhu -o ',outfile,' -r ',readgroup,' ',infile)
	run_command('samtools index ',outfile)
}
#export_read_group('merged','KT9.plasmid')

export_read_groups <- function(config,stem,ref)
{
	for (sample in rownames(config@samples))
	{
		#export_read_group(stem,sample)
		remove_duplicates(sample,ref)
		#export_pileup(sample,ref)
	}
}
#export_read_groups(config,'merged')

export_unmapped_reads <- function(stem)
{
	bamfile <- concat('tmp/',stem,'.bam')
	fastqfile <- concat('unmapped/',stem,'.unmapped.fastq')
	trimmedfastqfile <- concat('unmapped/',stem,'.unmapped.fastq.trimmed')
	fastafile <- concat('unmapped/',stem,'.unmapped.fasta')
	#run_command('bam2fastq --force --no-aligned --unaligned --no-filtered -o ',fastqfile,' ',bamfile)
	#run_command('cd unmapped; DynamicTrim.pl ',concat(stem,'.unmapped.fastq'))
	run_command('perl $VARDB_RUTIL_HOME/fq_all2std.pl fq2fa ',trimmedfastqfile,' > ',fastafile)
	run_command('perl $VARDB_RUTIL_HOME/fastq2table.pl -a ',trimmedfastqfile)
	run_command('sort merged.nodup.unmapped.fastq.trimmed.table.txt | uniq -dc > unmapped/merged.unique.txt')
}
#export_unmapped_reads('bam/PXB0220-0002.wk11.bam')

remove_duplicates <- function(stem,ref)
{
	metricsfile <- concat('tmp/',stem,'.nodup.metrics')
	outfile <-  concat('tmp/',stem,'.nodup.bam')
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MarkDuplicates.jar'
	str <- concat(str,' INPUT=bam/',stem,'.bam')
	str <- concat(str,' OUTPUT=',outfile)
	str <- concat(str,' METRICS_FILE=',metricsfile)
	str <- concat(str,' REMOVE_DUPLICATES=true')
	run_command(str)
	run_command('samtools index ',outfile)
}

count_codons <- function(config)
{
	for (subject in config@subjects$subject)
	{
		run_command('Rscript $VARDB_RUTIL_HOME/count_codons.r subject=',subject)
	}
}

make_tables <- function()
{
	run_command('Rscript $VARDB_RUTIL_HOME/make_tables.r')
}

analyze_reads_merged <- function(config)
{
	stem <- 'merged'
	ref <- 'KT9'
	#preprocess(config)
	remove_exact_duplicates_for_all_samples(config)
	#map_reads_for_all_samples(config)
	#merge_bams(config)
	#realign_indels(stem,ref)
	#recalibrate(concat(stem,'.realigned'),ref)
	#output_bam(stem,'realigned.recal')
	
	#call_variants(stem,ref)
	#filter_variants(stem,ref)
	#export_read_groups(config,stem,ref)
	#count_codons(config)
	#make_tables()
	#export_unmapped_reads(concat(stem,'.nodup'))
}

#map_reads_for_sample('KT9.specific','KT9')

analyze_reads_merged(config)
#solexa_qa('KT9.plasmid')

#Rscript ~/workspace/vardb-util/src/main/r/analyze_reads_merged.r
