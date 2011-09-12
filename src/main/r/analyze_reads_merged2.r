source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen2')
testrun <- FALSE
config <- new('nextgenconfig')
config@illumina.dir <- 'GA_RunData/110802_HWUSI-EAS1611_00068_FC634PPAAXX/Unaligned'

make_folders <- function(config, subdirs='ref,fastq,tmp,bam,unmapped,vcf,pileup,qc,counts')
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

preprocess <- function(config)
{
	fastq.dir <- config@fastq.dir
	temp.dir <- config@tmp.dir
	#runs <- config@runs
	run_command('rm -r ',temp.dir,'/*')
	for (sample in config@samplenames)
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
	
	for (sample in config@samplenames)
	{
		dir.from <- concat(temp.dir,'/',sample)
		run_command('gunzip ',dir.from,'/*')
	}
	run_command('')
	
	for (sample in config@samplenames)
	{
		dir.from <- concat(temp.dir,'/',sample)
		fastqfile <- concat(fastq.dir,'/',sample,'.fastq')
		run_command('cat ',dir.from,'/* > ',fastqfile)
		checkFileExists(fastqfile)
	}
	run_command('')
	#run_command('rm -r ',temp.dir)
}
#preprocess(config)

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
#trim_solexaqa(config,'PXB0220-0002.wk13')
#cd out/fastq; DynamicTrim.pl PXB0220-0002.wk13.fastq
#LengthSort.pl PXB0220-0002.wk13.fastq.trimmed

trim_all <- function(config)
{
	for (sample in config@samplenames)
	{
		trim_solexaqa(config,sample)
		#trim_prinseq(config,sample)
	}
}
#trim_all(config)

######################################################################

find_ref_for_sample <- function(config, sample)
{
	refs <- unique(config@runs[which(config@runs$sample==sample),'ref'])
	if (length(refs)>1)
		throw('more than 1 unique ref per sample: ',sample)
	return(refs[1])
}
#find_ref_for_sample(config,'110617HBV.HBV01')

get_reffile <- function(config, ref)
{
	ref.dir <- config@ref.dir
	reffile <- concat(config@ref.dir,'/',ref,'.fasta')
	seq <- config@refs[ref,1]	
	print(concat('writing ref file ',reffile))
	write.fasta(s2c(seq), ref, file.out=reffile)
	return(reffile)
}
#get_reffile(config,'IL28B-70')

run_bwa <- function(config, sample, ref, fastq.ext='.trimmed.fastq', q=20)#, .trimmed.fastq
{
	reffile <- get_reffile(config,ref)
	print(concat('ref for sample ',sample,': ',ref,' (',reffile,')'))
	stem <- concat(sample,'.',ref)
	fqfile <- concat(config@fastq.dir,'/',sample,fastq.ext)
	tmp.dir <- config@tmp.dir

	saifile <- concat(tmp.dir,'/',stem,'.sai')
	samfile <- concat(tmp.dir,'/',stem,'.sam')
	
	run_command('bwa index ',reffile)
	run_command('bwa aln -q ',q,' ',reffile,' ',fqfile,' > ',saifile)
	run_command('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
	checkFileExists(samfile)
}
#run_bwa(config,'110617HBV.HBV03') #'PXB0220-0002.wk13')

add_read_groups <- function(config,sample,ref)
{
	stem <- concat(sample,'.',ref)
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
#add_read_groups(config,'110617HBV.HBV03')

map_reads_for_sample <- function(config, sample)
{
	try({
		print(concat('mapping reads for sample: ',sample))
		ref <- find_ref_for_sample(config,sample)
		run_bwa(config,sample,ref)
		add_read_groups(config,sample,ref)
	}, silent=FALSE)
}
#map_reads_for_sample(config,'110617HBV.HBV03')

map_reads_for_all_samples <- function(config)
{
	for (sample in config@samplenames)
	{
		map_reads_for_sample(config,sample)
	}
}
#map_reads_for_all_samples(config)

merge_bams <- function(config)
{
	tmp.dir <- config@tmp.dir
	outfile <- concat(tmp.dir,'/merged.bam')
	
	str <- 'java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar'
	str <- concat(str,' MERGE_SEQUENCE_DICTIONARIES=true')
	str <- concat(str,' CREATE_INDEX=true')
	str <- concat(str,' VALIDATION_STRINGENCY=LENIENT')
	for (sample in config@samplenames)
	{
		ref <- find_ref_for_sample(config,sample)
		infile <- concat(tmp.dir,'/',sample,'.',ref,'.rg.bam') #.dedup
		str <- concat(str,' INPUT=',infile)
	}
	str <- concat(str,' OUTPUT=',outfile)
	run_command(str)
}
#merge_bams(config)

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
	tmp.dir <- config@tmp.dir
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

export_pileup <- function(config,sample,ref)
{
	print(concat('ref: ',ref))
	run_command('python $VARDB_RUTIL_HOME/export_pileup.py ', sample,' ',ref,' ',config@bam.dir,' ',config@pileup.dir)
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
	for (sample in config@samplenames)
	{
		#export_read_group(config,stem,sample)
		ref <- find_ref_for_sample(config,sample)
		sample <- concat(sample,'.',ref,'.rg') #hack!
		filter_bam(config,sample)
		##remove_duplicates(sample,ref)
		export_pileup(config,concat(sample,'.filtered'),ref)
	}
}
#export_read_groups(config,'merged')

##################################################################3

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
	make_folders(config)
	preprocess(config)
	trim_all(config)
	##remove_exact_duplicates_for_all_samples(config)
	map_reads_for_all_samples(config)
	merge_bams(config)
	##realign_indels(config,stem)
	#recalibrate(config,stem) #concat(stem,'.realigned'))
	#output_bam(config,stem,'recal')
	
	#call_variants(config,stem)
	#filter_variants(config,stem)
	#export_read_groups(config,stem)
	#count_codons(config)
	#make_tables(config)
	#export_unmapped_reads(config,concat(stem,'.nodup'))
}

args <- commandArgs(TRUE) # from R.utils package
config@out.dir <- args$out

analyze_reads_merged(config)

#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r out=out

