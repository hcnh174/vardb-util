getCodonCountFilename <- function(config, sample, type)
{
	return(concat(config@counts.dir,'/',sample,'.',type,'.txt'))
}
#getCodonCountFilename(config,'KT9.random__HCV-KT9','codons')
#getRefForSample('HCV-KT9.random__HCV-KT9')

#uses GATK custom walker to count all the bases at each position
countBasesForSample <- function(config, sample, bam.dir=config@bam.dir)#, out.dir=config@basecounts.dir)
{
	ref <- getRefForSample(sample)
	reffile <- getRefFile(config,ref)
	bamfile <- concat(bam.dir,'/',sample,'.bam')
	ntcountsfile <- getCodonCountFilename(config,sample,'nt')
	codoncountsfile <- getCodonCountFilename(config,sample,'codons')
	aacountsfile <- getCodonCountFilename(config,sample,'aa')
	checkFileExists(reffile)
	checkFileExists(bamfile)

	str <- 'java -Xmx8g'
	str <- concat(str,' -cp $VARDB_UTIL_HOME/target/gatk-walkers.jar:$GATK_HOME/GenomeAnalysisTK.jar')
	str <- concat(str,' org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants')
	str <- concat(str,' -dt NONE')
	str <- concat(str,' -et NO_ET')	
	str <- concat(str,' -R ',reffile)
	str <- concat(str,' -I ',bamfile)
	str <- concat(str,' --validation_strictness strict')
	str <- concat(str,' --ntcounts ',ntcountsfile)
	str <- concat(str,' --codoncounts ',codoncountsfile)
	str <- concat(str,' --aacounts ',aacountsfile)
	for (interval in getIntervalsForSample(config,sample))
	{
		str <- concat(str,' -L ',interval)
	}
	runCommand(str)
	#checkFileExists(ntcountsfile)
	#checkFileExists(codoncountsfile)
	#checkFileExists(aacountsfile)
}
#countBasesForSample(config,'nextgen1-1A__HCV-KT9')
#countBasesForSample(config,'KT9.random__HCV-KT9')


countBases <- function(config, samples=config@samples)
{
	for (sample in samples)
	{
		print(sample)
		try(countBasesForSample(config,sample))
	}
}
#countBases(config)