source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'common.r',sep=''))
loadUtilFiles('nextgen')
#setCurDir('nextgen2')
print(getwd())

config <- new('nextgenconfig')

preprocess(config)

