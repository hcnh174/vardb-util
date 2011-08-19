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


plot.summary.rms2 <- function (x, at, log = FALSE, q = c(0.7, 0.8, 0.9, 0.95, 0.99), 
		xlim, nbar, cex = 1, nint = 10, cex.c = 0.5, cex.t = 1, clip = c(-1e+30, 
				1e+30), main, fieldnames=NULL, ...) 
{
	scale <- attr(x, "scale")
	adjust <- attr(x, "adjust")
	Type <- x[, "Type"]
	x <- x[Type == 1, , drop = FALSE]
	lab <- dimnames(x)[[1]]
	if (!is.null(fieldnames))
		lab=fieldnames
	effect <- x[, "Effect"]
	se <- x[, "S.E."]
	if (!log && any(Type == 2)) {
		fun <- exp
		tlab <- scale[2]
	}
	else {
		fun <- function(x) x
		if (log) {
			if (length(scale) == 2) 
				tlab <- scale[2]
			else tlab <- paste("exp(", scale[1], ")", sep = "")
		}
		else tlab <- scale[1]
	}
	if (!length(scale)) 
		tlab <- ""
	if (!missing(main)) 
		tlab <- main
	augment <- if (log | any(Type == 2)) 
				c(0.1, 0.5, 0.75, 1)
			else 0
	n <- length(effect)
	out <- qnorm((max(q) + 1)/2)
	if (missing(xlim) && !missing(at)) 
		xlim <- range(if (log) logb(at) else at)
	else if (missing(xlim)) {
		xlim <- fun(range(c(effect - out * se, effect + out * 
										se)))
		xlim[1] <- max(xlim[1], clip[1])
		xlim[2] <- min(xlim[2], clip[2])
	}
	else augment <- c(augment, if (log) exp(xlim) else xlim)
	fmt <- function(k) {
		m <- length(k)
		f <- character(m)
		for (i in 1:m) f[i] <- format(k[i])
		f
	}
	lb <- ifelse(is.na(x[, "Diff."]), lab, lab)
	#lb <- ifelse(is.na(x[, "Diff."]), lab, paste(lab, " - ", 
	#    fmt(x[, "High"]), ":", fmt(x[, "Low"]), sep = ""))
	plot.new()
	par(new = TRUE)
	mxlb <- 0.1 + max(strwidth(lb, units = "inches", cex = cex))
	tmai <- par("mai")
	on.exit(par(mai = tmai))
	par(mai = c(tmai[1], mxlb, 1.5 * tmai[3], tmai[4]))
	outer.widths <- fun(effect + out * se) - fun(effect - out * 
					se)
	if (missing(nbar)) 
		nbar <- n
	npage <- ceiling(n/nbar)
	is <- 1
	for (p in 1:npage) {
		ie <- min(is + nbar - 1, n)
		plot(1:nbar, rep(0, nbar), xlim = xlim, ylim = c(1, nbar), 
				type = "n", axes = FALSE, xlab = "", ylab = "")
		if (cex.t > 0) 
			title(tlab, cex = cex.t)
		lines(fun(c(0, 0)), c(nbar - (ie - is), nbar), lty = 2)
		if (log) {
			pxlim <- pretty(exp(xlim), n = nint)
			pxlim <- sort(unique(c(pxlim, augment)))
			pxlim <- pxlim[pxlim >= exp(xlim[1])]
			if (!missing(at)) 
				pxlim <- at
			axis(3, logb(pxlim), lab = format(pxlim))
		}
		else {
			pxlim <- pretty(xlim, n = nint)
			pxlim <- sort(unique(c(pxlim, augment)))
			pxlim <- pxlim[pxlim >= xlim[1]]
			if (!missing(at)) 
				pxlim <- at
			axis(3, pxlim)
		}
		imax <- (is:ie)[outer.widths[is:ie] == max(outer.widths[is:ie])][1]
		for (i in is:ie) {
			confbar(nbar - (i - is + 1) + 1, effect[i], se[i], 
					q = q, type = "h", fun = fun, cex = cex.c, labels = i == 
							imax, clip = clip, ...)
			mtext(lb[i], 2, 0, at = nbar - (i - is + 1) + 1, 
					cex = cex, adj = 1, las = 1)
		}
		if (adjust != "") {
			adjto <- paste("Adjusted to:", adjust, sep = "")
			xx <- par("usr")[2]
			if (nbar > ie) 
				text(xx, nbar - (ie - is + 1), adjto, adj = 1, 
						cex = cex)
			else title(sub = adjto, adj = 1, cex = cex)
		}
		is <- ie + 1
	}
	invisible()
}
#plot.summary.rms2

#http://www.r-bloggers.com/exporting-r-output-to-ms-word-with-r2wd-an-example-session/
wdBody.anything <- function(output)
{
	# This function takes the output of an object and prints it line by line into the word document
	# Notice that in many cases you will need to change the text font into courier new roman...
	a <- capture.output(output)
	for(i in seq_along(a))
	{
		wdBody(format(a[i]))
	}
}

