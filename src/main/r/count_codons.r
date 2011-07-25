source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'common.r',sep=''))
loadUtilFiles('nextgen')
#setCurDir('nextgen2')

print(getwd())

config <- new('nextgenconfig')
#config@variants.dir <- 'c:/temp/variants-bwa/'
#config@variants.dir <- 'f:/variants/'

args <- commandArgs(TRUE)
subject <- args[1]

#subject <- 'PXB0219-0011'
#subject <- 'KT9'

params <- new('sampleparams', subject=subject)
print(concat('subject=',subject))

#Rscript c:/workspace/vardb-util/src/main/r/count_codons.r KT9
#Rscript count_codons.r KT9
#Rscript count_codons.r PXB0218-0007
#Rscript count_codons.r PXB0219-0011
#Rscript count_codons.r PXB0219-0018
#Rscript count_codons.r PXB0220-0002
variants <- count_codons_for_subject(config, subject)

