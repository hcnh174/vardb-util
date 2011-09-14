#diagnostics

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


solexa_qa <- function(config,sample)
{
	run_command('cd ',config@qc.dir,'; SolexaQA.pl ../fastq/',sample,'.fastq')
}
#solexa_qa(config,'KT9.plasmid')
