import sys
import os


def run_command(command):
	print command
	os.system(command)

def qseq2fastq(sample):
	run_command("perl qseq2fastq.pl -a qseq/"+sample+".qseq -v T")
	run_command("mv "+sample+".fastq fastq/"+sample+".fastq")

def btrim(sample,ref):
	#run_command("perl fastq_btrim.pl -a fastq/"+sample+".fastq -t T -v T")
	run_command("cat fastq/"+sample+".fastq | perl ../tools/bin/trim_bwa_style.pl -q 20 > trimmed/"+sample+".trimmed.fastq")

def trim(sample,ref):
	run_command("cd trimmed; DynamicTrim.pl ../fastq/"+sample+".fastq")
	run_command("cd quality; SolexaQA.pl ../fastq/"+sample+".fastq")


def run_bowtie(sample,ref):
	fqfile = "trimmed/"+sample+".fastq.trimmed"
	#fqfile = "fastq/"+sample+".fastq"
	samfile = "sam/"+sample+"."+ref+".sam"
	#command = "bowtie -v 3 -S "+ref+" trimmed/"+sample+"_Btrim.fq sam/"+sample+"."+ref+".sam"
	#run_command("bowtie -v 2 -S "+ref+" fastq/"+sample+".fq sam/"+sample+"."+ref+".sam")
	#run_command("bowtie -v 3 -S "+ref+" fastq/"+sample+".fq sam/"+sample+"."+ref+".sam")
	#run_command("bowtie -S "+ref+" fastq/"+sample+".fq sam/"+sample+"."+ref+".sam")
	#run_command("bowtie -S --seedmms 3 --tryhard --un unmapped/"+sample+"."+ref+".un.txt "+ref+" fastq/"+sample+".fastq sam/"+sample+"."+ref+".sam")
	#run_command("bowtie -S --un unmapped/"+sample+"."+ref+".un.txt "+ref+" "+fqfile+" sam/"+sample+"."+ref+".sam")
	run_command("bowtie -S --un unmapped/"+sample+"."+ref+".un.txt "+ref+" "+fqfile+" "+samfile)

def run_bwa(sample,ref):
	#fqfile = "trimmed/"+sample+"_Btrim.fq"
	fqfile = "trimmed/"+sample+".fastq.trimmed"
	#fqfile = "fastq/"+sample+".fastq"
	reffile = "ref/"+ref+".fasta"
	saifile = "sai/"+sample+"."+ref+".sai"
	samfile = "sam/"+sample+"."+ref+".sam"
	run_command("bwa aln "+reffile+" "+fqfile+" > "+saifile) # -q 20
	run_command("bwa samse "+reffile+" "+saifile+" "+fqfile+" > "+samfile)

#def run_smalt(sample,ref):
#	run_command("smalt_x86_64 "+reffile+" "+fqfile+" > "+saifile) # -q 20

def run_mosaik(sample,ref):
	stem = sample+"."+ref
	run_command("MosaikBuild -fr ref/"+ref+".fasta -oa mosaik/"+ref+".dat")
	run_command("MosaikBuild -q fastq/"+sample+".fastq -st illumina -out mosaik/"+stem+".mkb")
	run_command("MosaikAligner -mm 2 -act 20 -bw 13 -mhp 100 -in mosaik/"+stem+".mkb -ia mosaik/"+ref+".dat -out mosaik/"+stem+".mka -rur mosaik/"+stem+".un")#4
	run_command("MosaikSort -in mosaik/"+stem+".mka -out mosaik/"+stem+".mks")
	run_command("MosaikText -in mosaik/"+stem+".mks -sam sam/"+stem+".sam")
	run_command("gunzip sam/"+stem+".sam.gz")

def sam2bam(sample,ref):
	run_command("./sam2bam.sh "+sample+" "+ref)

def export_pileup(sample,ref):
	run_command("python export_pileup.py "+sample+" "+ref)

def vcf(sample,ref):
	reffile = "ref/"+ref+".fasta"
	stem = sample+"."+ref
	run_command("samtools mpileup -u -f "+reffile+" bam/"+stem+".bam > bcf/"+stem+".bcf")
	run_command("bcftools view bcf/"+stem+".bcf > vcf/"+stem+".vcf")

def map_reads(sample,ref):
	try:
		print 'sample = %s, ref = %s' % (sample, ref)
		#qseq2fastq(sample)
		#trim(sample,ref)
		run_bowtie(sample,ref)
		#run_bwa(sample,ref)
		#run_mosaik(sample,ref)
		#run_smalt(sample,ref)
		sam2bam(sample,ref)
		vcf(sample,ref)
		export_pileup(sample,ref)
	except:
	    print "Unexpected error:", sys.exc_info()[0]
	    raise

def count_codons(sample,ref,numrefs):
	run_command("Rscript count_codons.r "+sample+" "+ref+" "+str(numrefs)+" \"24-47,144-167\"")

def count_lines(sample):
	os.system("wc --lines qseq/" + sample+".qseq")


sample = sys.argv[1]
ref = sys.argv[2]

map_reads(sample,ref)
#trim(sample,ref)
#sam2bam(sample,ref)
#vcf(sample,ref)

#vcf(sample,ref)
#export_pileup(sample,ref)

