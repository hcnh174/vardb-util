source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('classes,core,util,mapping,counts,tables,fragments,charts', subdir='nextgen')
#setCurDir('nextgen/merged')
setwd('\\\\ubuntu\\merged')
config <- loadConfig('d:/projects/analysis/nextgen/merged','\\\\TS-XELF44\\share\\nextgendata')#

#makeReferenceVsVariantPresentation(config,'MP-424', minreads=100)
#makeReferenceVsVariantPresentation(config,minreads=100)