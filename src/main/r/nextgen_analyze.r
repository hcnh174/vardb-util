source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen_classes,nextgen_core,nextgen_util,nextgen_mapping,nextgen_counts,nextgen_tables,nextgen_fragments,nextgen_report')
config <- loadConfig()

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

analyzeReadsForGroup(config,'KT9')
analyzeReadsForGroup(config,'confirm_plasmid_with_new_reagents')
analyzeReadsForGroup(config,'HBV_nucleoside_analogues')
analyzeReadsForGroup(config,'MP-424')
analyzeReadsForGroup(config,'BMS-790052_BMS-650032')
analyzeReadsForGroup(config,'hcv_infection')
analyzeReadsForGroup(config,'BMS-790052_MP-424')
analyzeReadsForGroup(config,'NS3_V36A_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_Y93H_mutations_maintained')

doForSample <- function(config,sample)
{
	#mapReads(config,sample)
	#if (config@filter) filterBams(config,sample)
	#writeConsensusForBams(config,sample)
	#exportPileup(config,sample)
	#countCodons(config,sample)
}
#doForSample(config,config@samples)

doForGroup <- function(config,group)
{
	samples <- getSamplesForGroup(config,group)
	#doForSample(config,samples)
	#analyzeReadsForSample(config,samples)
	writeCodonTables(config,group,minreads=40)
	writeAminoAcidTables(config,group,minreads=40)
	concatTablesByGroup(config,group)
	makeAminoAcidBarcharts(config,group)
}
doForGroup(config,'KT9')
doForGroup(config,'confirm_plasmid_with_new_reagents')
doForGroup(config,'HBV_nucleoside_analogues')

doForGroup(config,'MP-424')
doForGroup(config,'BMS-790052_BMS-650032')
doForGroup(config,'hcv_infection')

doForGroup(config,'BMS-790052_MP-424')
doForGroup(config,'NS3_V36A_mutation_maintained')
doForGroup(config,'NS5A_L31V_mutation_maintained')

doForGroup(config,'NS5A_L31V_Y93H_mutations_maintained')
