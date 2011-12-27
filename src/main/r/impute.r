source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
#loadUtilFiles('nextgen2,nextgen_mapping,nextgen_counts,nextgen_tables')

setClass("imputeconfig",
	representation(
		Ne='numeric',
		chunksize='numeric',
		ref.dir='character',
		in.dir='character',
		out.dir='character',
		tmp.dir='character',
		chr.lengths='data.frame'
	),
	prototype(
		Ne=20000,
		chunksize=5000000,
		ref.dir='~/mnt/impute/ref',
		in.dir='~/mnt/impute/in',
		out.dir='~/mnt/impute/out',
		tmp.dir='~/mnt/impute/tmp'
	)
)

setMethod("initialize", "imputeconfig", function(.Object)
{	
	.Object@chr.lengths <- loadDataFrame(concat(Sys.getenv("VARDB_RUTIL_HOME"),'/chrlengths.txt'), idcol='chr')
	.Object@chr.lengths$length <- as.numeric(.Object@chr.lengths$length)
	.Object
})

padChr <- function(n)
{
	chr <- as.character(n)
	if (n<10)
		chr <- concat('0',chr)
	return(chr)
}
#padChr(2);padChr(20)

createSubfolders <- function(config, dir)
{
	for (chr in 1:22)
	{
		chrstr <- padChr(chr)
		subdir <- concat(dir,'/chr',chrstr)
		runCommand('mkdir ',subdir)
	}
}
#createSubfolders(config, config@tmp.dir)
#createSubfolders(config, config@out.dir)

############################################################

exportChromosome <- function(config, data, chr)
{
	chrstr <- padChr(chr)
	print(concat('exporting chr', chr))
	dir <- 'd:/temp'
	genofile <- concat(dir,'/gwas.chr',chrstr,'.gen')
	samplefile <- concat(dir,'/gwas.chr',chrstr,'.sample')
	strandfile <- concat(dir,'/gwas.chr',chrstr,'.strand')
	
	export.impute(data[,data@gtdata@chromosome==chr], genofile=genofile,
			samplefile=samplefile, strandfile=strandfile, cachesizeMb=128)
}
#exportChromosome(config, data0, 19)

exportChromosomes <- function(config, data)
{
	for (chr in 1:22) #	for (chr in c(2,4,6,8,10,15,17))
	{
		exportChromosome(config,data,chr)
	}
}
#exportChromosomes(config,data1)

#exportChromosomes <- function(config, data)
#{
#	for (chr in c(2,4,6,7,8,10,15,17))
#	{
#		exportChromosome(config,data,chr)
#	}
#}
##exportChromosomes(config,data1)

#####################################################################

liftOverChromosome <- function(config, chr, chainfile=concat(config@ref.dir,'/hg18ToHg19.over.chain'))
{
	chrstr <- padChr(chr)	
	strandfile <- concat(config@in.dir,'/gwas.chr',chrstr,'.strand')
	print(strandfile)
	
	data <- read.table(strandfile, header=FALSE, sep=' ', stringsAsFactors=FALSE)
	colnames(data) <- splitFields('snp,pos,strand')
	data$chr <- concat('chr',chr)
	data$pos2 <- data$pos+1
	data <- data.frame(chr=data$chr, start=data$pos, end=data$pos2, name=data$snp)
	
	oldbedfile <- concat(config@tmp.dir,'/gwas.chr',chrstr,'.strand.hg18.bed')
	newbedfile <- concat(config@tmp.dir,'/gwas.chr',chrstr,'.strand.hg19.bed')
	unmappedfile <- concat(config@tmp.dir,'/unmapped.chr',chrstr,'.txt')
	
	writeTable(data, oldbedfile, row.names=FALSE, col.names=FALSE)
	
	str <- 'liftOver'	
	#str <- concat(str,' -multiple')
	str <- concat(str,' ',oldbedfile)
	str <- concat(str,' ',chainfile)
	str <- concat(str,' ',newbedfile)
	str <- concat(str,' ',unmappedfile)
	#print(str)
	runCommand(str)
	
	checkFileExists(newbedfile)	
	data <- read.table(newbedfile, header=FALSE, sep='\t', stringsAsFactors=FALSE)
	colnames(data) <- splitFields('chr,pos,pos2,snp')
	#print(head(data))
	
	data <- data.frame(snp=data$snp, pos=data$pos)
	data$strand <- '+'
	#print(head(data))
	
	strandfile2 <- concat(config@in.dir,'/gwas.chr',chrstr,'.hg19.strand')
	writeTable(data, strandfile2, row.names=FALSE, col.names=FALSE)
	
	#runCommand('rm ',oldbedfile)
	#runCommand('rm ',newbedfile)
	#runCommand('rm ',unmappedfile)
}
#liftOverChromosome(config,2)

liftOver <- function()
{
	for (chr in 1:22)
	{
		liftOverChromosome(config,chr)
	}
}
#liftOver()

#####################################################################

partitionChromosome <- function(config,chr)
{
	total <- as.numeric(config@chr.lengths[as.character(chr),'length'])
	#print(concat('chr',chr,' has ',total,' bases'))
	partitions <- data.frame()
	for (start in seq(1,total,config@chunksize))
	{
		end <- start+config@chunksize-1
		if (end > total)
			end <- total
		#print(concat('chunk start=',start,', end=',end))
		partitions <- rbind(partitions, data.frame(start=start, end=end))
	}
	return(partitions)
}
#partitions <- partitionChromosome(config,1)

#####################################################################

formatInterval <- function(chr,start,end)
{
	return(concat('chr',padChr(chr),'.',format(start,scientific=FALSE),'.',format(end,scientific=FALSE)))
}
#formatInterval(22, 20.4e6, 20.5e6)

prephase_interval <- function(config, chr, start, end)
{
	print(concat('prephase_interval: ',chr,' ',start,':',end))
	chrstr <- padChr(chr)
	
	ref.sub.dir <- concat(config@ref.dir,'/ALL_1000G_phase1interim_jun2011_impute')
	
	mapfile <- concat(ref.sub.dir,'/genetic_map_chr',chr,'_combined_b37.txt')
	genotypefile <- concat(config@in.dir,'/gwas.chr',chrstr,'.gen')
	strandfile <- concat(config@in.dir,'/gwas.chr',chrstr,'.hg19.strand')
	outfile <- concat(config@tmp.dir,'/chr',chrstr,'/prephase.',formatInterval(chr,start,end),'.impute2')
	
	checkFilesExist(mapfile,genotypefile,strandfile)
	
	str <- 'impute2'
	str <- concat(str,' -phase')
	str <- concat(str,' -include_buffer_in_output')
	str <- concat(str,' -m ',mapfile)
	str <- concat(str,' -g ',genotypefile)
	str <- concat(str,' -strand_g ',strandfile)
	str <- concat(str,' -int ',start,' ',end)
	str <- concat(str,' -Ne ',config@Ne)
	str <- concat(str,' -o ',outfile)
	#print(str)
	runCommand(str)
	
	checkFileExists(outfile)
}
#prephase_interval(config, 22, 20.4e6, 20.5e6)

prephase_chromosome <- function(config, chr)
{
	print(concat('prephase_chromosome: ',chr))
	partitions <- partitionChromosome(config,chr)
	for (rowname in rownames(partitions))
	{
		start <- partitions[rowname,'start']
		end <- partitions[rowname,'end']
		prephase_interval(config,chr,start,end)
	}
}
#prephase_chromosome(config,22)

prephase <- function(config)
{
	print('prephase')
	for (chr in 1:22)
	{		
		prephase_chromosome(config,chr)
	}
}
#prephase(config)

#########################################

impute_interval <- function(config, chr, start, end)
{
	print(concat('impute_interval: ',chr,' ',start,':',end))
	
	chrstr <- padChr(chr)
	ref.dir <- concat(config@ref.dir,'/ALL_1000G_phase1interim_jun2011_impute')
	mapfile <- concat(ref.dir,'/genetic_map_chr',chr,'_combined_b37.txt')
	stem <- concat(ref.dir,'/ALL_1000G_phase1interim_jun2011_chr',chr,'_impute')
	refhapfile <- concat(stem,'.hap.gz')
	legendfile <- concat(stem,'.legend')
	
	knownhapsfile <- concat(config@tmp.dir,'/chr',chrstr,'/prephase.',formatInterval(chr,start,end),'.impute2')
	outfile <- concat(config@out.dir,'/chr',chrstr,'/impute.',formatInterval(chr,start,end),'.impute2')
	
	checkFilesExist(mapfile,refhapfile,legendfile,knownhapsfile)
	
	str <- 'impute2'
	str <- concat(str,' -m ',mapfile)
	str <- concat(str,' -h ',refhapfile)
	str <- concat(str,' -l ',legendfile)
	str <- concat(str,' -known_haps_g ',knownhapsfile)
	str <- concat(str,' -int ',start,' ',end)
	str <- concat(str,' -Ne ',config@Ne)
	str <- concat(str,' -o ',outfile)
	#print(str)
	runCommand(str)
	
	checkFileExists(outfile)
}
#impute(config, 22, 20.4e6, 20.5e6)

impute_chromosome <- function(config, chr, partitions=NULL)
{
	print(concat('impute_chromosome: ',chr))
	if (is.null(partitions))
		partitions <- partitionChromosome(config,chr)
	for (rowname in rownames(partitions))
	{
		start <- partitions[rowname,'start']
		end <- partitions[rowname,'end']
		try(impute_interval(config,chr,start,end))
	}
}
#impute_chromosome(config,19)

impute <- function(config)
{
	print('prephase')
	for (chr in 1:23)
	{		
		try(impute_chromosome(config,chr))
	}
}
#impute(config)

##########################################

snptest_chromosome <- function(config, phenofile, pheno, chr)
{
	chrstr <- padChr(chr)
	infile <- concat('genome/chr',chrstr,'.impute2')
	outfile <- concat(config@out.dir,'/',pheno,'-chr',chrstr,'.out')
			
	str <- 'snptest'
	str <- concat(str,' -pheno ',pheno)
	str <- concat(str,' -frequentist 1')
	str <- concat(str,' -method threshold')
	str <- concat(str,' -nowarn')
	str <- concat(str,' -hwe')
	str <- concat(str,' -data ',infile,' ',phenofile)
	str <- concat(str,' -o out/',pheno,'-chr',chrstr,'.out')
	print(str)
	#runCommand(str)
	
	#runCommand("tr -s \"NA\" \""+str(n)+"\" < out/"+pheno+"-chr"+str(n)+".out > out/"+pheno+"-chr"+str(n)+".tmp")
}
#snptest_chromosome(config,'patients.txt','svr',1)





concatenate_results <- function(pheno)
{
	runCommand("cat out/"+pheno+"-chr1.tmp > out/"+pheno+".out")	
	for (n in 2:23)
	{
		print(concat("concatening results for chr ",n))
		write_chromosome_number(pheno,n)
		runCommand("tail -n +2 out/",pheno,"-chr",n,".tmp >> out/",pheno,".out")
	}
	runCommand("cut -f 1- --delimiter=' ' --output-delimiter='\t' out/",pheno,".out > out/",pheno,".txt")
	#runCommand("rm out/*.out")
	#runCommand("rm out/*.tmp")
}

snptest <- function(config, phenofile, pheno)
{
	for (n in 1:22)
	{
		snptest_chromosome(config,phenofile,pheno)
		
	}
	
	#add chrosome numbers
	
}
#snptest('patients.txt','svr')

#snptest(phenofile,pheno)
#concatenate_results(pheno)
#runCommand("head out/"+pheno+".txt")
#runCommand("tail out/"+pheno+".txt")
#cp out/*.txt ~/analysis/imputed

# python snptest.py ./patients.txt svr
# python snptest.py ./patients.txt hcv0
# python snptest.py ./patients.txt hbdiff
# python snptest.py ./patients.txt plt01
# python snptest.py ./patients.txt plt04
# python snptest.py ./patients.txt wbc01
# python snptest.py ./patients.txt wbc04
# python snptest.py ./patients.txt hcc



