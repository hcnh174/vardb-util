source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen_classes,nextgen_core,nextgen_util,nextgen_mapping,nextgen_counts,nextgen_tables,nextgen_fragments,nextgen_report')
config <- loadConfig()

writeRefs(config)

#analyzeReads<- function(config)
#{
	preprocess(config)
	trimSamples(config)
	mapReads(config)
	outputBams(config)
	#fixBaiFiles(config)
	filterBams(config)
	writeConsensusForBams(config)
	exportPileup(config)
	countCodons(config)
	writeCodonTables(config)
	concatTablesByGroup(config)
	#makePiecharts(config)
#}

args <- commandArgs(TRUE) # from R.utils package
config@out.dir <- args$out

analyzeReads(config)


#analyzeReadsForGroup(config,'KT9')

#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r out=out

for (group in removeElements(config@groups,'KT9'))
{
	analyzeReadsForGroup(config,group)
}

#clearNextgenOutput(config)

analyzeReadsForGroup(config,'KT9')
analyzeReadsForGroup(config,'confirm_plasmid_with_new_reagents')
analyzeReadsForGroup(config,'HBV_nucleoside_analogues')
analyzeReadsForGroup(config,'hcv_infection')

analyzeReadsForGroup(config,'MP-424')
analyzeReadsForGroup(config,'BMS-790052_BMS-650032')

analyzeReadsForGroup(config,'BMS-790052_MP-424')
analyzeReadsForGroup(config,'NS3_V36A_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_Y93H_mutations_maintained')
analyzeReadsForGroup(config,'BMS-605339')
analyzeReadsForGroup(config,'BMS-788329')
analyzeReadsForGroup(config,'BMS-821095')

writeConsensusForBam(config,'nextgen4-8F__HCV-KT9')


doForSample <- function(config,sample)
{
	mapReads(config,sample)	
	#if (config@filter) filterBams(config,sample)
	#writeConsensusForBams(config,sample)
	#exportPileup(config,sample)
	#countCodons(config,sample)
}
#doForSample(config,config@samples)

doForGroup <- function(config,group)
{
	samples <- getSamplesForGroup(config,group)
	#trimSamples(config,samples)
	doForSample(config,samples)
	#analyzeReadsForSample(config,samples)
	#writeCodonTables(config,group,minreads=1)
	#writeAminoAcidTables(config,group,minreads=1)
	#concatTablesByGroup(config,group)
	#makeAminoAcidBarcharts(config,group)
}

for (rowname in rownames(config@data))
{
	ref <- config@data[rowname,'ref']
	stem <- concat(rowname,'__',ref)
	print(stem)
	try(exportUnmappedReads(config,stem))
}


#stem <- 'nextgen2-6I'#nextgen2-5I'
runBwa(config,'nextgen2-6I')
runBwa(config,'nextgen2-6J')
runBwa(config,'nextgen3-1L')
runBwa(config,'nextgen3-8J')
#analyzeUnmappedReads(config,stem)
getMapStats('out/bam/nextgen2_6I__HCV-KT9.bam')



getMapStats('out/bam/nextgen3-7L__HCV-KT9.bam')
getMapStats('out/bam/PXB0220-0002.wk08__HCV-KT9.bam')
exportPileup(config,'PXB0220-0002.wk08__HCV-KT9')
countCodons(config,'PXB0220-0002.wk08__HCV-KT9')



doForGroup(config,'BMS-605339')
getMapStats('out/bam/nextgen4-8E__HCV-KT9.bam')
getMapStats('out/bam/nextgen4-8F__HCV-KT9.bam')
exportUnmappedReads(config,'nextgen4-8E__HCV-KT9')
exportUnmappedReads(config,'nextgen4-8F__HCV-KT9')


doForGroup(config,'BMS-821095')
doForGroup(config,'confirm_plasmid_with_new_reagents')
doForGroup(config,'HBV_nucleoside_analogues')

doForGroup(config,'MP-424')
doForGroup(config,'BMS-790052_BMS-650032')
doForGroup(config,'hcv_infection')

doForGroup(config,'BMS-790052_MP-424')
doForGroup(config,'NS3_V36A_mutation_maintained')
doForGroup(config,'NS5A_L31V_mutation_maintained')

doForGroup(config,'NS5A_L31V_Y93H_mutations_maintained')

doForGroup(config,'BMS-788329')

doForGroup(config,'KT9')

####################################################################

mergeBamsForSamples(config)
#getMappingReport(config)
exportPileup(config)
countCodons(config)
writeCodonTables(config)
writeAminoAcidTables(config)
concatTablesByGroup(config)


