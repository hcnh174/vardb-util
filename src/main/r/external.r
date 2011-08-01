

# copied directly from limma removeExt function
stripExtension <- function (filenames) 
{
	filenames <- as.character(filenames)
	n <- length(filenames)
	if (length(grep("\\.", filenames)) < n) 
		return(filenames)
	ext <- sub("(.*)\\.(.*)$", "\\2", filenames)
	if (all(ext[1] == ext)) 
		return(sub("(.*)\\.(.*)$", "\\1", filenames))
	else return(filenames)
}
#stripExtension('tables.Rnw')


#http://r.789695.n4.nabble.com/Lattice-histogram-with-vertical-lines-td876969.html
#only works with lattice plots
addLine<- function(a=NULL, b=NULL, v = NULL, h = NULL, ..., once=F)
{
	tcL <- trellis.currentLayout()
	k<-0
	for(i in 1:nrow(tcL))
	{
		for(j in 1:ncol(tcL))
		{
			if (tcL[i,j] > 0)
			{
				k<-k+1
				trellis.focus("panel", j, i, highlight = FALSE)
				if (once)
					panel.abline(a=a[k], b=b[k], v=v[k], h=h[k], ...)
				else panel.abline(a=a, b=b, v=v, h=h, ...)
				trellis.unfocus()
			}
		}
	}
}
#addLine(v=aanum - start + 1 - 0.5, col='red', lty=2)

# http://tolstoy.newcastle.edu.au/R/help/05/12/17151.html
list2ascii <- function(x,file=paste(deparse(substitute(x)),".txt",sep="")) {
	
	# MHP July 7, 2004
	# R or S function to write an R list to an ASCII file.
	# This can be used to create files for those who want to use
	# a spreadsheet or other program on the data.
	#
	tmp.wid = getOption("width")  # save current width
	options(width=10000)          # increase output width
	sink(file)                    # redirect output to file
	print(x)                      # print the object
	sink()                        # cancel redirection
	options(width=tmp.wid)        # restore linewidth
	return(invisible(NULL))       # return (nothing) from function
	
}

skew <- function(x)
{
	(sum((x-mean(x))^3/sqrt(var(x))^3)/length(x))
}

kurtosis <- function(x)
{
	(sum((x-mean(x))^4/var(x)^2)/length(x) - 3)
}

.ls.objects <- function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 10) {
	# based on postings by Petr Pikal and David Hinds to the r-help list in 2004
	# modified by: Dirk Eddelbuettel (http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session) 
	# I then gave it a few tweaks (show size as megabytes and use defaults that I like)
	# a data frame of the objects and their associated storage needs.
	napply <- function(names, fn) sapply(names, function(x)
					fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.size <- napply(names, object.size) / 10^6 # megabytes
	obj.dim <- t(napply(names, function(x)
						as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.dim)
	names(out) <- c("Type", "Size", "Rows", "Columns")
	out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}


