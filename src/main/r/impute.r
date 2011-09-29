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
			ref.dir='~/impute/ref',
			in.dir='~/impute/in',
			out.dir='~/impute/out',
			tmp.dir='~/impute/tmp'
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
	for (chr in 1:23)
	{
		#if (chr==19)
		#	next
		exportChromosome(config,data,chr)
	}
}
#exportChromosomes(config,data0)

#####################################################################

liftOver <- function(config, chr, chainfile)
{
	#dir <- '~/out.dir'
	dir <- 'd:/temp'
	chrstr <- padChr(chr)	
	strandfile <- concat(dir,'/gwas.chr',chrstr,'.strand')
	print(strandfile)
	
	dataframe <- read.table(strandfile, header=FALSE, sep=' ', stringsAsFactors=FALSE)
	colnames(dataframe) <- splitFields('snp,pos,strand')
	dataframe$chr <- concat('chr',chr)
	print(head(dataframe))
	#data <- data.frame(chr=data$chr, start=data$pos, end=data$pos, name=data$snp)
	return(data)
}
#data <- liftOver(config,19,'')

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

prephase_interval <- function(config, chr, start, end)
{
	print(concat('prephase_interval: ',chr,' ',start,':',end))
	chrstr <- padChr(chr)
	mapfile <- concat(config@in.dir,'/example.chr',chrstr,'.map')
	genotypefile <- concat(config@in.dir,'/example.chr',chrstr,'.study.gens')
	strandfile <- concat(config@in.dir,'/example.chr',chrstr,'.study.strand')
	outfile <- concat(config@tmp.dir,'/prephase.chr',chrstr,'.impute2')
	
	str <- 'impute2'
	str <- concat(str,' -phase')
	str <- concat(str,' -include_buffer_in_output')
	str <- concat(str,' -m ',mapfile)
	str <- concat(str,' -g ',genotypefile)
	str <- concat(str,' -strand_g ',strandfile)
	str <- concat(str,' -int ',start,' ',end)
	str <- concat(str,' -Ne ',config@Ne)
	str <- concat(str,' -o ',outfile)
	print(str)
	#run_command(str)
}
#prephase_chromosome(config, 22, 20.4e6, 20.5e6)

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

prephase <- function(config)
{
	print('prephase')
	for (chr in 1:23)
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
	
	mapfile <- concat(config@ref.dir,'/genetic_map_chr',chr,'_combined_b37.txt')
	stem <- concat(config@ref.dir,'/ALL_1000G_phase1interim_jun2011_chr',chrstr,'_impute')
	refhapfile <- concat(stem,'.hap.gz')
	legendfile <- concat(stem,'.legend')
	
	knownhapsfile <- concat(config@tmp.dir,'/prephase.chr',chrstr,'.impute2')
	outfile <- concat(config@out.dir,'/impute.chr',chrstr,'.impute2')
	
	str <- 'impute2'
	str <- concat(str,' -m ',mapfile)
	str <- concat(str,' -h ',refhapfile)
	str <- concat(str,' -l ',legendfile)
	str <- concat(str,' -known_haps_g ',knownhapsfile)
	str <- concat(str,' -int ',start,' ',end)
	str <- concat(str,' -Ne ',config@Ne)
	str <- concat(str,' -o ',outfile)
	print(str)
	#run_command(str)
}
#impute(config, 22, 20.4e6, 20.5e6)

impute_chromosome <- function(config, chr)
{
	print(concat('impute_chromosome: ',chr))
	partitions <- partitionChromosome(config,chr)
	for (rowname in rownames(partitions))
	{
		start <- partitions[rowname,'start']
		end <- partitions[rowname,'end']
		impute_interval(config,chr,start,end)
	}
}

impute <- function(config)
{
	print('prephase')
	for (chr in 1:23)
	{		
		impute_chromosome(config,chr)
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
	#run_command(str)
	
	#run_command("tr -s \"NA\" \""+str(n)+"\" < out/"+pheno+"-chr"+str(n)+".out > out/"+pheno+"-chr"+str(n)+".tmp")
}
#snptest_chromosome(config,'patients.txt','svr',1)





concatenate_results <- function(pheno)
{
	run_command("cat out/"+pheno+"-chr1.tmp > out/"+pheno+".out")	
	for (n in 2:23)
	{
		print(concat("concatening results for chr ",n))
		write_chromosome_number(pheno,n)
		run_command("tail -n +2 out/",pheno,"-chr",n,".tmp >> out/",pheno,".out")
	}
	run_command("cut -f 1- --delimiter=' ' --output-delimiter='\t' out/",pheno,".out > out/",pheno,".txt")
	#run_command("rm out/*.out")
	#run_command("rm out/*.tmp")
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
#run_command("head out/"+pheno+".txt")
#run_command("tail out/"+pheno+".txt")
#cp out/*.txt ~/analysis/imputed

# python snptest.py ./patients.txt svr
# python snptest.py ./patients.txt hcv0
# python snptest.py ./patients.txt hbdiff
# python snptest.py ./patients.txt plt01
# python snptest.py ./patients.txt plt04
# python snptest.py ./patients.txt wbc01
# python snptest.py ./patients.txt wbc04
# python snptest.py ./patients.txt hcc



