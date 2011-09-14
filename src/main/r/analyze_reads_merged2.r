source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen2')
testrun <- FALSE
config <- new('nextgenconfig')
config@illumina.dir <- 'GA_RunData/110802_HWUSI-EAS1611_00068_FC634PPAAXX/Unaligned'

.make_folders <- function(config, subdirs='ref,fastq,tmp,bam,unmapped,vcf,pileup,qc,counts')
{
	dir <- config@out.dir
	run_command('mkdir ',dir,' -p')
	for (subdir in splitFields(subdirs))
	{
		subdir <- concat(dir,'/',subdir)
		run_command('mkdir ',subdir,' -p; rm ',subdir,'/*')
	}
}
# .make_folders(config)

preprocess <- function(config)
{
	.make_folders(config)
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

add_read_groups <- function(config,sample)
{
	tmp.dir <- config@tmp.dir
	samfile <- concat(tmp.dir,'/',sample,'.sam')
	bamfile <- concat(tmp.dir,'/',sample,'.bam')
	baifile <- concat(tmp.dir,'/',sample,'.bai')
	
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
	checkFileExists(baifile)	
}
#add_read_groups(config,'110617HBV.HBV07@HBV-RT')

run_bwa <- function(config, sample, fastq.ext='.trimmed.fastq', q=20)#, .trimmed.fastq
{
	ref <- get_ref_for_sample(sample)
	reffile <- get_reffile(config,ref)
	fqfile <- concat(config@fastq.dir,'/',sample,fastq.ext)
	tmp.dir <- config@tmp.dir

	samfile <- concat(tmp.dir,'/',sample,'.sam')
	saifile <- concat(tmp.dir,'/',sample,'.sai')
	
	run_command('bwa index ',reffile)
	run_command('bwa aln -q ',q,' ',reffile,' ',fqfile,' > ',saifile)
	run_command('bwa samse ',reffile,' ',saifile,' ',fqfile,' > ',samfile)
	checkFileExists(samfile)
	
	add_read_groups(config,sample)
	
	run_command('rm ',samfile)
	run_command('rm ',saifile)
}
#run_bwa(config,'110617HBV.HBV07@HBV-RT')

map_reads <- function(config)
{
	for (sample in config@samplenames)
	{
		run_bwa(config,sample)
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
	
	for (sample in unique(config@runs[which(config@runs$ref=='HBV-RT'),'sample']))
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
	
#	for (sample in config@samplenames)
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

filter_bams <- function(config)
{
	for (sample in config@samplenames)
	{
		filter_bam(config,sample)
	}
}
#filter_bams(config)

####################################################

export_pileup_for_sample <- function(config,sample)
{
	ref <- get_ref_for_sample(sample)
	run_command('python $VARDB_RUTIL_HOME/export_pileup.py ', sample,'.filtered ',ref,' ',config@bam.dir,' ',config@pileup.dir)
}
#export_pileup(config,'merged','KT9')

export_pileup <- function(config)
{
	for (sample in config@samplenames)
	{
		export_pileup_for_sample(config,sample)
	}
}
export_pileup(config)

################################################

count_codons <- function(config)
{
	#for (subject in config@subjects$subject)
	for (subject in unique(config@runs$subject))
	{
		try(count_codons_for_subject(config, subject))
	}
}
#count_codons(config)

make_tables <- function(config)
{
	sweaveToPdf(concat(Sys.getenv("VARDB_RUTIL_HOME"),'/tables.rnw'))
	#run_command('Rscript $VARDB_RUTIL_HOME/make_tables.r dir=',config@count.dir)
}

analyze_reads<- function(config)
{
	preprocess(config)
	trim_all(config)
	##remove_exact_duplicates_for_all_samples(config)
	map_reads(config)
	merge_bams(config)
	##realign_indels(config,stem)
	#recalibrate(config,stem) #concat(stem,'.realigned'))
	#output_bam(config,stem,'recal')
	
	#call_variants(config,stem)
	#filter_variants(config,stem)
	#unmerge_bams(config)
	#count_codons(config)
	#make_tables(config)
	#export_unmapped_reads(config,concat(stem,'.nodup'))
}

args <- commandArgs(TRUE) # from R.utils package
config@out.dir <- args$out

analyze_reads(config)

#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r out=out



variants <- count_codons_for_subject(config, '10348001')








#
#
#getCodonPositionsForRegion <- function(config, region)
#{
#	#region <- 'NS3-156@NS3aa156'
#	ref <- config@regions[which(config@regions==region),'ref']
#	startnt <- config@regions[region,'startnt']
#	startntrel <- config@regions[region,'startntrel']
#	#print(concat(start,':',end))
#	sequence <- config@refs[ref,'sequence']#, regions[region,'start'], regions[region,'end'])
#	endntrel <- startntrel + nchar(sequence)
#	#print(sequence)
#	#print(translateCodon(sequence))
#	positions <- data.frame()
#	codon <- startntrel/3 + 1
#	if (startntrel %% 3 !=0)
#		stop(concat('start number is not a multiple of 3: ',startntrel,' in region ',region))
#	print(concat('ref: ',ref,' startnt: ',startnt,', startntrel: ',startntrel,', codon: ',codon))
#	for (position in seq(startntrel,endntrel,3))
#	{
#		try({
#		index <- position-startntrel+1
#		print(index)
#		refcodon <- extractSequence(sequence, index, index+2)
#		refaa <- translateCodon(refcodon)
#		positions <- rbind(positions, data.frame(codon=codon, ntnum=startnt+position, relpos=position, refcodon=refcodon, refaa=refaa))
#		codon <- codon + 1
#		})
#	}	
#	return(positions)
#}
##getCodonPositionsForRegion(config,'NS3-156@NS3aa156')#'HBV-RT@HBV-RT')
#
#createNtCountTable <- function(config, data, params)
#{	
#	region <- params@region
#	offset <- config@regions[region,'offset']
#	start <- config@regions[region,'start']
#	end <- config@regions[region,'end']
#	data <- subset(data, position >= start-offset & position <= end-offset)
#	
#	counts <- data.frame()
#	positions <- getCodonPositionsForRegion(config,params@region)
#	for (ntnum in positions$ntnum)
#	{
#		try({
#			aanum <- positions[which(positions$ntnum==ntnum),'codon']
#			for (offset in 0:2)
#			{
#				ntnum2 <- ntnum + offset
#				nt <- data[which(data$position==ntnum2),'nt']
#				#print(head(nt))
#				#if (params@drop.ambig)
#				#	nt <- removeAmbiguousCodons(nt)
#				#print(concat('nt vector length: ',length(nt)))
#				freqs <- sort(xtabs(as.data.frame(nt)), decreasing=TRUE)
#				total <- sum(freqs)				
#				rank <- 1
#				for (base in names(freqs))
#				{
#					count <- freqs[base]
#					freq <- count/total
#					row <- data.frame(ntnum=ntnum2, aanum=aanum, nt=base, rank=rank, count=count, freq=freq) 
#					counts <- rbind(counts,row)
#					rank <- rank+1
#				}
#			}
#		}, silent=FALSE)
#	}
#	counts$nt <- as.character(counts$nt)
#	counts$ntnum <- factor(counts$ntnum)
#	counts$aanum <- factor(counts$aanum)
#	counts <- appendSampleParams(counts, params)
#	return(counts)
#}
#
#count_codons_for_region <- function(config, data, params, variantdata)
#{
#	print(concat('count_codons_for_region: ',as.character(params@region)))
#	variantdata@nt <- rbind(variantdata@nt, createNtCountTable(config, data, params))
#	variantdata@codons <- rbind(variantdata@codons, createCodonCountTable(config, data, params))
#	variantdata@aa <- rbind(variantdata@aa, createAminoAcidCountTable(config, data, params))	
#	return(variantdata)
#}
#
#count_codons_for_sample <- function(config, params, variantdata)
#{
#	sample <- as.character(params@sample)
#	print(concat('count_codons_for_sample: ',sample))
#	#ref <- config@samples[sample,'ref'] # look up the ref for the sample
#	ref <- get_ref_for_sample(sample)
#	filename <- concat(config@pileup.dir,'/',sample,'.filtered.txt'); print(filename) # load the corresponding data file
#	data <- loadDataFrame(filename)
#	print(concat('loaded file ',filename,'. contains ',nrow(data),' reads'))
#	# each sample has one or more runs targeting different regions
#	for (region in config@runs[ config@runs[,'sample']==sample,'region'])
#	{
#		params@region <- region
#		variantdata <- count_codons_for_region(config, data, params, variantdata)
#	}
#	return(variantdata)
#}
##tables <- count_codons_for_sample('218-7_03-01')
#
#params <- '10348001'
#
#count_codons_for_subject <- function(config, params)
#{
#	if (is.character(params))
#		params <- new('sampleparams', subject=params)
#	print(concat('count_codons_for_subject: ',as.character(params@subject)))
#	variantdata <- new('variantdata')
#	samples <- config@runs[which(config@runs$subject==params@subject),]
#	if (nrow(samples)==0)
#		stop(concat('cannot find any samples for subject ',params@subject))
#	for (sample in samples[,'sample'])
#	{
#		params@sample <- sample
#		params@replicate <- samples[sample,'replicate']
#		try({variantdata <- count_codons_for_sample(config, params, variantdata)}, silent=FALSE)
#	}
#	counts.dir <- concat(config@counts.dir,'/')
#	writeTable(variantdata@nt, concat(counts.dir,params@subject,'.nt.txt'), row.names=FALSE)
#	writeTable(variantdata@codons, concat(counts.dir,params@subject,'.codons.txt'), row.names=FALSE)
#	writeTable(variantdata@aa, concat(counts.dir,params@subject,'.aa.txt'), row.names=FALSE)
#	return(variantdata)
#}
##variants <- count_codons_for_subject(config, '110719-4.p30.NS3aa36@NS3-36')
