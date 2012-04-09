source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen_classes,nextgen_core,nextgen_util,nextgen_mapping,nextgen_counts2,nextgen_tables2,nextgen_fragments,nextgen_charts')
config <- loadConfig()

writeRefs(config)

#analyzeReads<- function(config)
#{
	preprocess(config)
	trimSamples(config)
	#maskSamples(config)
	#dedupSamples(config)
	
	mapReads(config)
	mergeBamsForSamples(config)
	#countBases(config)
	
	#outputBams(config)
	#fixBaiFiles(config)
	#filterBams(config)
	writeConsensusForBams(config)
	exportPileup(config)
	countCodons(config)
	writeCodonTables(config)
	writeAminoAcidTables(config)
	concatTablesByGroup(config)
	#makePiecharts(config)
#}

args <- commandArgs(TRUE) # from R.utils package
config@out.dir <- args$out

analyzeReads(config)

#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r out=out

clearNextgenOutput(config)

analyzeReadsForGroup(config,'KT9')
analyzeReadsForGroup(config,'confirm_plasmid_with_new_reagents')
analyzeReadsForGroup(config,'HBV_nucleoside_analogues')
analyzeReadsForGroup(config,'hcv_infection')

analyzeReadsForGroup(config,'MP-424')
analyzeReadsForGroup(config,'MK-7009')
analyzeReadsForGroup(config,'BMS-790052_BMS-650032')

analyzeReadsForGroup(config,'BMS-790052_MP-424')
analyzeReadsForGroup(config,'NS3_V36A_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_mutation_maintained')
analyzeReadsForGroup(config,'NS5A_L31V_Y93H_mutations_maintained')
analyzeReadsForGroup(config,'BMS-605339')
analyzeReadsForGroup(config,'BMS-788329')
analyzeReadsForGroup(config,'BMS-821095')
analyzeReadsForGroup(config,'BMS-790052')
analyzeReadsForGroup(config,'TMC435')

#"TMC435" ,"MK-7009","BMS-790052_BMS-650032","MP-424" ,"unknown","NS3_V36A_mutation_maintained"


#nextgen_ppt.r - only runs on Win32 version
minreads <- 100
makeReferenceVsVariantPresentationByGroup(config,'KT9', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'confirm_plasmid_with_new_reagents', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'HBV_nucleoside_analogues', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'hcv_infection', minreads=minreads)

makeReferenceVsVariantPresentationByGroup(config,'MP-424', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'MK-7009', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'BMS-790052_BMS-650032', minreads=minreads)

makeReferenceVsVariantPresentationByGroup(config,'BMS-790052_MP-424', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'NS3_V36A_mutation_maintained', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'NS5A_L31V_mutation_maintained', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'NS5A_L31V_Y93H_mutations_maintained', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'BMS-605339', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'BMS-788329', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'BMS-821095', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'BMS-790052', minreads=minreads)
makeReferenceVsVariantPresentationByGroup(config,'TMC435', minreads=minreads)

