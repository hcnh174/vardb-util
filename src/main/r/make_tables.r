source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'common.r',sep=''))
loadUtilFiles('nextgen')
setCurDir('nextgen2')

config <<- new('nextgenconfig')
#config@counts.dir <- 'n:/counts/'
config@counts.dir <- 'c:/Documents and Settings/nhayes/My Documents/My Dropbox/counts/'

sweaveToPdf(concat(Sys.getenv("VARDB_RUTIL_HOME"),'tables.rnw'))
