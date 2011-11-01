source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen_classes,nextgen_core,nextgen_util,nextgen_mapping,nextgen_counts,nextgen_tables')
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

analyzeReadsForGroup(config,'MP-424')
analyzeReadsForGroup(config,'BMS-790052_BMS-650032')
analyzeReadsForGroup(config,'hcv_infection')
analyzeReadsForGroup(config,'BMS-790052_MP-4242')
analyzeReadsForGroup(config,'NS3_V36A_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_Y93H_mutations_maintained')

#KT9
#confirm_plasmid_with_new_reagents
#HBV_nucleoside_analogues
