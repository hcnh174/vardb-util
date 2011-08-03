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

	str = "java -Xmx2g -jar $PICARD_HOME/MarkDuplicates.jar"
	str = str+" INPUT=tmp/"+stem+".bam"
	str = str+" OUTPUT=tmp/"+stem+".dedup.bam"
	str = str+" METRICS_FILE="+metricsfile
	run_command(str)
	run_command("samtools index tmp/"+stem+".dedup.bam")

def realign_indels(sample,ref):
	stem = sample+"."+ref+".rg.dedup" #assume duplicates already marked
	reffile = "ref/"+ref+".fasta"
	bamfile = "tmp/"+stem+".bam"
	intervalfile = "tmp/"+stem+".intervals"

	str = "java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator"
	str = str+" -R "+reffile
	str = str+" -I "+bamfile
	str = str+" -o "+intervalfile
	run_command(str)
	
	str = "java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T IndelRealigner"
	str = str+" -R "+reffile
	str = str+" -I "+bamfile
	str = str+" -targetIntervals "+intervalfile
	str = str+" -o tmp/"+stem+".realigned.bam"
	run_command(str)

def recalibrate(sample,ref): 
	stem = sample+"."+ref+".rg.dedup.realigned" #assume duplicates marked and realigned 
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
	str = str+" -outputDir ./qc/"+sample+"."+ref
	run_command(str) 

 	str = "java -Xmx2g -jar $GTAK_HOME/GenomeAnalysisTK.jar -T TableRecalibration"
 	str = str+" -l INFO"
 	str = str+" -R "+reffile
 	str = str+" -I "+bamfile
 	str = str+" -recalFile "+recalfile	
	str = str+" -o "+outfile	
	run_command(str)

def output_bam(sample,ref,suffix):
	infile = "tmp/"+sample+"."+ref+"."+suffix+".bam"
	outfile = "bam/"+sample+"."+ref+".bam"
	run_command("cp "+infile+" "+outfile)
	run_command("samtools index "+outfile)

def cleanup(sample,ref):
	run_command("rm tmp/"+sample+"."+ref+".*")

def call_variants(sample,ref):
	stem = sample+"."+ref
	reffile = "ref/"+ref+".fasta"
	vcffile = "vcf/"+stem+".vcf"
	
	str = "java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper"
	str = str+" -R "+reffile
	str = str+" -I bam/"+stem+".bam"
	str = str+" -B:mask,VCF config/"+ref+".mask.vcf"
	str = str+" -o "+vcffile
	str = str+" -stand_call_conf 10.0"	#30.0" #50.0
	str = str+" -stand_emit_conf 10.0"
	str = str+" -L config/"+ref+".interval_list"
	#str = str+" -dcov 50"
	run_command(str)

def filter_variants(sample,ref):
	stem = sample+"."+ref
	reffile = "ref/"+ref+".fasta"
	vcffile = "vcf/"+stem+".vcf"

	str = "java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T VariantFiltration"
	str = str+" -R "+reffile
	str = str+" -B:variant,VCF "+vcffile
	str = str+" -o vcf/"+stem+".filtered.vcf"
	str = str+" --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\""
	str = str+" --filterName \"HARD_TO_VALIDATE\""
	str = str+" --filterExpression \"SB >= -1.0\""
	str = str+" --filterName \"StrandBiasFilter\""
	str = str+" --filterExpression \"QUAL < 10\""
	str = str+" --filterName \"QualFilter\""
	str = str+" --clusterWindowSize 10"
	str = str+" --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\""
	str = str+" --filterName \"HARD_TO_VALIDATE\""
	str = str+" --filterExpression \"QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10\""
	str = str+" --filterName GATKStandard"
	run_command(str)

def export_pileup(sample,ref):
	run_command("python $VARDB_RUTIL_HOME/export_pileup.py "+sample+" "+ref)

def analyze_reads(sample, ref):
	run_bwa(sample,ref)
	sam2bam(sample,ref)
	add_read_groups(sample,ref)
	mark_duplicates(sample,ref)
	realign_indels(sample,ref)
	recalibrate(sample,ref)
	output_bam(sample,ref,"rg.dedup.realigned.recal")
	call_variants(sample,ref)
	filter_variants(sample,ref)
	export_pileup(sample,ref)
	cleanup(sample,ref)

sample = sys.argv[1]
ref = sys.argv[2]

#analyze_reads(sample,ref)
#export_pileup(sample,ref)
#call_variants(sample,ref)
#filter_variants(sample,ref)

#call_variants(sample,ref)
#filter_variants(sample,ref)

def call_variants_combined(subject,replicates,ref):

	reffile = "ref/"+ref+".fasta"
	vcffile = "vcf/"+subject+".vcf"
	
	str = "java -jar $GTAK_HOME/GenomeAnalysisTK.jar -T UnifiedGenotyper"
	str = str+" -R "+reffile
	for replicate in replicates:
		str = str+" -I bam/"+subject+"."+replicate+"."+ref+".bam"
	#str = str+" -B:mask,VCF config/"+ref+".mask.vcf"
	#str = str+" -stand_call_conf 10.0"	#30.0" #50.0
	#str = str+" -stand_emit_conf 10.0"
	str = str+" -L config/"+ref+".interval_list"
	#str = str+" -dcov 50"
	str = str+" -o "+vcffile
	#str = str+" --output_mode EMIT_ALL_SITES"
	run_command(str)

call_variants_combined('KT9',['plasmid','random','specific'],'KT9')
call_variants_combined('PXB0218-0007',['wk10','wk11','wk12','wk13','wk15'],'KT9')
call_variants_combined('PXB0219-0011',['wk08','wk09','wk10','wk11','wk12'],'KT9')
call_variants_combined('PXB0219-0018',['wk08','wk09','wk10','wk12','wk14','wk15'],'KT9')
call_variants_combined('PXB0220-0002',['wk08','wk09','wk10','wk11','wk12','wk13'],'KT9')




