#http://factbased.blogspot.jp/2012/04/working-with-strings.html#more
strhead <- function(s,n=1)
{
	if(n<0) 
		substr(s,1,nchar(s)+n) 
	else 
		substr(s,1,n)
}

strtail <- function(s,n=1) 
{
	if(n<0) 
		substring(s,1-n) 
	else 
		substring(s,nchar(s)-n+1)
}

strsubst <- function(template, map, verbose=getOption("verbose")) {
	pat <- "\\$\\([^\\)]+\\)"
	res <- template
	map[["$"]] <- "$"
	m <- gregexpr(pat, template)
	idx <- which(sapply(m, function(x) x[[1]]!=-1)) # faster than 1:length(template)?
	for (i in idx) {
		line <- template[[i]]
		if(verbose) cat("input: |", template[[i]], "|\n")
		starts <- m[[i]]
		ml <- attr(m[[i]], "match.length")
		sym <- substring(line, starts+2, starts+ml-2)
		repl <- map[sym]
		idx1 <- is.null(repl)
		repl[idx1] <- sym[idx1]
		norepl <- substring(line, c(1, starts+ml), c(starts-1, nchar(line)))
		res[[i]] <- paste(norepl, c(repl, ""), sep="", collapse="") # more elegant?
		if (verbose) cat("output: |", res[[i]], "|\n")
	}
	return(res)
}

strparse <- function(pat, x) {
	parsed <- regexpr(pat, x, perl=TRUE)
	if (length(x)==1) {
		if(parsed[1]==-1) return(NULL)
		st <- attr(parsed, "capture.start")[1,]
		m <- substring(x, st, st + attr(parsed, "capture.length")[1,]-1)
		names(m) <- attr(parsed, "capture.names")
	} else {
		m <- do.call(rbind, lapply(seq_along(parsed), function(i) {
							if(parsed[i] == -1) return("")
							st <- attr(parsed, "capture.start")[i, ]
							substring(x[i], st, st + attr(parsed, "capture.length")[i, ] - 1)
						}))
		colnames(m) <- attr(parsed, "capture.names")
	}
	return(m)
}

strrecode <- function(pats, repls, x, ...) {
	res <- rep(NA, length(x))
	hits <- rep(FALSE, length(x))
	for (i in seq_along(pats)) {
#browser()
		new_hits <- grepl(pats[[i]],x[!hits],...)
		res[!hits][new_hits] <- repls[[i]]
		hits[!hits][new_hits] <- TRUE
		if(all(hits)) break
	}
	return(res)
}
