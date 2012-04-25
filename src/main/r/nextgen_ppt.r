Sys.setenv(JAVA_HOME='C:\\Program Files (x86)\\Java\\jdk1.6.0_31')
library(XLConnect)
library(R2PPT)

source(paste(Sys.getenv("VARDB_RUTIL_HOME"),'/common.r',sep=''))
loadUtilFiles('nextgen_classes,nextgen_core,nextgen_util,nextgen_mapping,nextgen_counts2,nextgen_tables2,nextgen_fragments,nextgen_charts')
#setCurDir('nextgen/merged')
setwd('\\\\ubuntu\\merged')
config <- loadConfig('d:/projects/analysis/nextgen/merged','\\\\TS-XELF44\\share\\nextgendata')#

makeReferenceVsVariantPresentationByGroup <- function(config, group, ppt=NULL, minreads=0, rowheight=20, colwidth=50, vspacer=5)
{
	ppt <- PPT.AddTitleSlide(ppt, title=concat('Results for ',group), title.fontsize=25, subtitle.fontsize=20)		
	subjects <- getSubjectsForGroup(config,group)
	for (subject in subjects)
	{
		title <- simpleCap(gsub ('_',' ',subject,ignore.case=T,perl=T))
		ppt <- PPT.AddTitleOnlySlide(ppt,title=title, title.fontsize=20)
		left <- 30;	top <- 100
		samples <- getSamplesForSubject(config,subject)
		for (region in getRegionsForSubject(config,subject))
		{
			gene <- strsplit(region,'aa', fixed=TRUE)[[1]][1]
			for (aanum in getFociForRegion(config,region))
			{
				try({
					tbl <- makeReferenceVsVariantTable(config,samples,region,aanum,minreads=minreads)			
					if (nrow(tbl)==1)
						tbl[1,1] <- 'sample'
					else tbl$sample <- sapply(tbl$sample,stripSubjectFromSample)
					tbl <- renameColumn(tbl,'sample',concat(gene,'aa',aanum))
					tbl <- replaceNAs(tbl)
					width <- (ncol(tbl)+1)*colwidth
					height <- (nrow(tbl)+1)*rowheight
					printcat('width=',width,', height=',height)
					ppt<-PPT.AddDataFrame(ppt, df=tbl, row.names=FALSE, col.names=TRUE, size=c(left,top,width,height))
					top <- top+height+vspacer
				})
			}
		}
	}
	return(ppt)
}
#makeReferenceVsVariantPresentationByGroup(config,'MP-424')

makeReferenceVsVariantPresentation <- function(config, groups=config@groups, ...)
{
	filename <- 'tables-ref_vs_variants.ppt'
	if (length(groups)==1)
		filename <- concat('tables-ref_vs_variants-',groups,'.ppt')
	deleteFile(filename)
	ppt <- PPT.Init(visible=FALSE, method='RDCOMClient')
	ppt <- PPT.AddTitleSlide(ppt, title=concat('Deep sequencing results'), subtitle=Sys.Date())	
	for (group in groups)
	{
		ppt <- makeReferenceVsVariantPresentationByGroup(config,group,ppt,...)
	}
	PPT.SaveAs(ppt,filename)
}
#makeReferenceVsVariantPresentation(config,'MP-424', minreads=100)
#makeReferenceVsVariantPresentation(config,minreads=100)