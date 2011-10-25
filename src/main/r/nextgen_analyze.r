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
	makePiecharts(config)
#}

args <- commandArgs(TRUE) # from R.utils package
config@out.dir <- args$out

analyzeReads(config)

#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r out=out

analyzeReadsForSample <- function(config,sample)
{
	#sample <- 'KT9.random__HCV-KT9'
	preprocess(config,sample)
	trimSamples(config,sample)
	mapReads(config,sample)
	filterBams(config,sample)
	writeConsensusForBams(config,sample)
	exportPileup(config,sample)
}

analyzeReadsForGroup <- function(config,group)
{
	samples <- config@data[which(config@data$group==group),'sample']
	analyzeReadsForSample(config,samples)
	countCodons(config,group)
	writeCodonTables(config,group)
}