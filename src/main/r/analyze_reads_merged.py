import sys
import os

def run_command(command):
	print command
	os.system(command)

def run_bwa(sample,ref):
	fqfile = "fastq/"+sample+".fastq"
	reffile = "ref/"+ref+".fasta"
	saifile = "tmp/"+sample+"."+ref+".sai"
	samfile = "tmp/"+sample+"."+ref+".sam"
	run_command("bwa aln "+reffile+" "+fqfile+" > "+saifile)
	run_command("bwa samse "+reffile+" "+saifile+" "+fqfile+" > "+samfile)

def sam2bam(sample,ref):
	stem = sample+"."+ref
	sam="tmp/"+stem+".sam"
	bam="tmp/"+stem+".bam"
	run_command("samtools view -bS -o "+bam+" "+sam)
	run_command("samtools sort "+bam+" tmp/"+stem)
	run_command("samtools index "+bam)
	
def add_read_groups(sample,ref):
	stem = sample+"."+ref
	str = "java -Xmx2g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar"
	str = str+" INPUT=tmp/"+stem+".bam"
	str = str+" OUTPUT=tmp/"+stem+".rg.bam"
	str = str+" RGSM="+stem
	str = str+" SORT_ORDER=coordinate RGID=1 RGLB=HCV RGPL=illumina RGPU=barcode"
	str = str+" CREATE_INDEX=true"
	run_command(str)

def mark_duplicates(sample,ref):
	stem = sample+"."+ref+".rg"
	metricsfile = "qc/"+stem+".dedup.metrics"
	outfile =  "tmp/"+stem+".dedup.bam"

	str = "java -Xmx2g -jar $PICARD_HOME/MarkDuplicates.jar"
	str = str+" INPUT=tmp/"+stem+".bam"
	str = str+" OUTPUT="+outfile
	str = str+" METRICS_FILE="+metricsfile
	run_command(str)
	run_command("samtools index "+outfile)

def realign_indels(stem,ref):
	
	reffile = "ref/"+ref+".fasta"
	bamfile = "tmp/"+stem+".bam"
	intervalfile = "tmp/"+stem+".intervals"
	outfile = "tmp/"+stem+".realigned.bam"

	str = "java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator"
	str = str+" -R "+reffile
	str = str+" -I "+bamfile
	str = str+" -o "+intervalfile
	run_command(str)
	
	str = "java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T IndelRealigner"
	str = str+" -R "+reffile
	str = str+" -I "+bamfile
	str = str+" -targetIntervals "+intervalfile
	str = str+" -o "+outfile
	run_command(str)

def recalibrate(stem,ref):  
	reffile = "ref/"+ref+".fasta"
	bamfile = "tmp/"+stem+".bam"
	maskfile = "config/"+ref+".mask.vcf"
	recalfile = "tmp/"+stem+".recal.csv"
	outfile = "tmp/"+stem+".recal.bam"
	
	str = "java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T CountCovariates"
	str = str+" -l INFO"
	str = str+" -R "+reffile
	str = str+" -I "+bamfile
	str = str+" -B:mask,VCF "+maskfile
	str = str+" --standard_covs"
	str = str+" -cov ReadGroupCovariate"
	str = str+" -cov QualityScoreCovariate"
	str = str+" -cov CycleCovariate"
	str = str+" -cov DinucCovariate"
	str = str+" -cov HomopolymerCovariate"
	str = str+" -cov MappingQualityCovariate"
	str = str+" -cov MinimumNQSCovariate"
	str = str+" -recalFile "+recalfile
	run_command(str)
	
	str = "java -Xmx2g -jar $GTAK_HOME/AnalyzeCovariates.jar"
	str = str+" -resources $GTAK_HOME/resources"
	str = str+" -recalFile "+recalfile
	str = str+" -outputDir ./qc/"+stem
	run_command(str) 

 	str = "java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T TableRecalibration"
 	str = str+" -l INFO"
 	str = str+" -R "+reffile
 	str = str+" -I "+bamfile
 	str = str+" -recalFile "+recalfile	
	str = str+" -o "+outfile	
	run_command(str)

def output_bam(stem,suffix):
	infile = "tmp/"+stem+"."+suffix+".bam"
	outfile = "bam/"+stem+".bam"
	run_command("cp "+infile+" "+outfile)
	run_command("samtools index "+outfile)

def cleanup(sample,ref):
	run_command("rm tmp/"+sample+"."+ref+".*")

def call_variants(stem,ref):

	reffile = "ref/"+ref+".fasta"
	bamfile = "bam/"+stem+".bam"
	outfile = "vcf/"+stem+".vcf"
	
	str = "java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper"
	str = str+" -R "+reffile
	str = str+" -I "+bamfile
	#str = str+" -B:mask,VCF config/"+ref+".mask.vcf"
	#str = str+" -stand_call_conf 10.0"	#30.0" #50.0
	#str = str+" -stand_emit_conf 10.0"
	str = str+" -L config/"+ref+".interval_list"
	#str = str+" -dcov 50"
	str = str+" -o "+outfile
	#str = str+" --output_mode EMIT_ALL_SITES"
	run_command(str)

def annotate_variants(stem,ref):

	reffile = "ref/"+ref+".fasta"
	bamfile = "bam/"+stem+".bam"
	infile = "vcf/"+stem+".vcf"
	outfile = "vcf/"+stem+".annotated.vcf"
	
	str = "java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T VariantAnnotator"
	str = str+" -l INFO"
	str = str+" -R "+reffile
	str = str+" -I "+bamfile
	str = str+" -B:variant,VCF "+infile
	str = str+" -o "+outfile
	str = str+" --useAllAnnotations"	
	run_command(str)

def filter_variants(stem, ref):
	reffile = "ref/"+ref+".fasta"
	infile = "vcf/"+stem+".annotated.vcf"
	outfile = "vcf/"+stem+".annotated.filtered.vcf"

	str = "java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T VariantFiltration"
	str = str+" -R "+reffile
	str = str+" -B:variant,VCF "+infile
	str = str+" -o "+outfile
	#str = str+" --clusterWindowSize 10"
	
	#basic indel filtering
	#str = str+" --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\""  #expression to match 10% of reads with MAPQ0
	#str = str+" --filterName \"HARD_TO_VALIDATE\""	
	#str = str+" --filterExpression \"SB >= -1.0\""
	#str = str+" --filterName \"StrandBiasFilter\""
	#str = str+" --filterExpression \"QUAL < 10\""
	#str = str+" --filterName \"QualFilter\""
	
	#snp filtering
	#str = str+" --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\""
	#str = str+" --filterName \"HARD_TO_VALIDATE\""
	#str = str+" --filterExpression \"QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10\""
	#str = str+" --filterName GATKStandard"
	
	#hard filtering
	str = str+" --filterExpression \"QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10\"" #For exomes with deep coverage per sample
	#str = str+" --filterExpression \"DP > 100 || MQ0 > 40 || SB > -0.10\"" #whole genomes with deep coverage 
	str = str+" --filterName GATKStandard"
	
	run_command(str)

def variants_to_table(stem, ref):
	reffile = "ref/"+ref+".fasta"
	vcffile = "vcf/"+stem+".vcf"
	str = "java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T VariantsToTable"
	str = str+" -R "+reffile
	str = str+" -B:variant,VCF "+vcffile
	str = str+" -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F AF"	
	str = str+" -o vcf/"+stem+".table"
	run_command(str)

###################################################3

def analyze_reads_for_sample(sample, ref):
	run_bwa(sample,ref)
	sam2bam(sample,ref)
	add_read_groups(sample,ref)
	mark_duplicates(sample,ref)

def analyze_reads_for_subject(subject,replicates,ref):
	for replicate in replicates:
		analyze_reads_for_sample(subject+"."+replicate, ref)

def analyze_reads_for_all_samples():
	analyze_reads_for_subject('KT9',['plasmid','random','specific'],'KT9')
	analyze_reads_for_subject('PXB0218-0007',['wk10','wk11','wk12','wk13','wk15'],'KT9')
	analyze_reads_for_subject('PXB0219-0011',['wk08','wk09','wk10','wk11','wk12'],'KT9')
	analyze_reads_for_subject('PXB0219-0018',['wk08','wk09','wk10','wk12','wk14','wk15'],'KT9')
	analyze_reads_for_subject('PXB0220-0002',['wk08','wk09','wk10','wk11','wk12','wk13'],'KT9')

def merge_bams():
	ref = 'KT9'
	prefix = " INPUT=tmp/"
	suffix = '.'+ref+'.rg.dedup.bam'
	outfile = "tmp/merged.bam"

	str = "java -Xmx2g -jar $PICARD_HOME/MergeSamFiles.jar"
	for replicate in ['plasmid','random','specific']:
		str = str+prefix+"KT9."+replicate+suffix
		
	for replicate in ['wk10','wk11','wk12','wk13','wk15']:
		str = str+prefix+"PXB0218-0007."+replicate+suffix
		
	for replicate in ['wk08','wk09','wk10','wk11','wk12']:
		str = str+prefix+"PXB0219-0011."+replicate+suffix
		
	for replicate in ['wk08','wk09','wk10','wk12','wk14','wk15']:
		str = str+prefix+"PXB0219-0018."+replicate+suffix
		
	for replicate in ['wk08','wk09','wk10','wk11','wk12','wk13']:
		str = str+prefix+"PXB0220-0002."+replicate+suffix
	
	str = str+" OUTPUT="+outfile
	run_command(str)
	run_command("samtools index "+outfile)
	
	
def analyze_reads_merged(ref):
	stem = 'merged'
	#analyze_reads_for_all_samples()
	#merge_bams()
	#realign_indels(stem,ref)
	#recalibrate(stem+'.realigned',ref)
	#output_bam(stem,'realigned.recal')
	
	#call_variants(stem,ref)
	#annotate_variants(stem,ref)
	filter_variants(stem,ref)
	#variants_to_table(stem,ref)
	#variants_to_table(stem+'.filtered',ref)
	#export_pileup(sample,ref)
	#cleanup(sample,ref)

analyze_reads_merged('KT9')

#NS3aa156	3885-3887
#NS3aa36	3525-3527
#NS5Aaa31	6348-6351
#NS5Aaa93	6534-6537


