source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen2,nextgen_mapping,nextgen_counts,nextgen_tables')
config <- new('nextgenconfig')
config@illumina.dir <- 'GA_RunData/110802_HWUSI-EAS1611_00068_FC634PPAAXX/Unaligned'

analyze_reads<- function(config)
{
	preprocess(config)
	#trim_all(config)
	##remove_exact_duplicates_for_all_samples(config)
	map_reads(config)
	#merge_bams(config)
	##realign_indels(config,stem)
	#recalibrate(config,stem) #concat(stem,'.realigned'))
	#output_bam(config,stem,'recal')
	
	#call_variants(config,stem)
	#filter_variants(config,stem)
	#unmerge_bams(config)
	output_bams(config)
	#count_codons(config)
	#make_tables(config)
	make_piecharts(config)
	#export_unmapped_reads(config,concat(stem,'.nodup'))
}

args <- commandArgs(TRUE) # from R.utils package
config@out.dir <- args$out

analyze_reads(config)

#Rscript $VARDB_RUTIL_HOME/analyze_reads_merged.r out=out



counts <- count_codons_for_subject(config, '10201689')
