source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen_classes,nextgen_core,nextgen_mapping,nextgen_counts,nextgen_tables')
config <- loadConfig()

analyze_reads<- function(config)
{
	preprocess(config)
	trim_all(config)
	map_reads(config)
	output_bams(config)
	fix_bai_files(config)
	filter_bams(config)
	export_pileup(config)
	count_codons(config)
	make_piecharts(config)
	#export_unmapped_reads(config,concat(stem,'.nodup'))
}

args <- commandArgs(TRUE) # from R.utils package
config@out.dir <- args$out

analyze_reads(config)

#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r out=out
