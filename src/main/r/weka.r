################################################################################
## File: 
##    readarff.r
##
## Description:
##    Provides a function for reading ARFF files into the R statistical
##    package. Currently supports dense file formats without support
##    for timestamps.
##
## Author:
##    Craig A. Struble, Ph.D.
##
## Created:
##    Mon Aug 18 15:16:19 CDT 2003
##
## Last Modified:
##    Wed Sep 10 14:39:44 CDT 2003 - modified lastChar variable to attempt to
##        support S-PLUS.
################################################################################
#
#read.arff <- function (file, header = TRUE, sep = ",", quote="\"", dec=".",
#		fill = TRUE, ...)
#{
#	comment.char <- "%"                          # ARFF file comment char
#	enum.open <- "{"
#	enum.close <- "}"
#	sep <- ","                                   # value separator
#	
#	###########################################################################
#	# Get the next token from the ARFF file
#	# file - an open file or connection
#	###########################################################################
#	nextArffToken <- function(file, lastChar) {
#		aToken <- FALSE # something implying end of file
#		if (lastChar == "") {
#			aChar <- readChar(file, 1)
#		} else {
#			aChar <- lastChar
#		}
#		
#		while (aToken == FALSE && aChar != "") {
#			# skip whitespace
#			while (regexpr("[[:space:]]", aChar) != -1 || aChar == sep) {
#				aChar <- readChar(file, 1)
#			}
#			# Check for comment
#			if (aChar == comment.char) {
#				# skip until the end of the line
#				readLines(file, n=1,ok=TRUE)
#				# Get the next character
#				aChar <- readChar(file,1)
#			} else {
#				aToken <- aChar
#			}
#		}
#		
#		if (aToken == "\"" || aToken == "'") { # Quoted string
#			quote <- aToken
#			aToken <- ""
#			aChar <- readChar(file, 1)
#			while(aChar != "" && aChar != quote) {
#				aToken <- paste(aToken,aChar,sep="")
#				if (aChar == "\\") { # Escaped character in string
#					aChar <- readChar(file, 1)
#					aToken <- paste(aToken,aChar,sep="")
#				}
#				aChar <- readChar(file, 1)
#			}
#			if (aChar == quote) {
#				aChar <- readChar(file, 1)
#			}
#		} else { # Unquoted string
#			aChar <- readChar(file, 1)
#			# Until a brace, whitespace, comma, or EOF
#			while (aChar != "" 
#					&& regexpr("[[:space:]]", aChar) == -1
#					&& !aChar %in% c(enum.open, enum.close, "\"", "'", sep)
#					&& aToken != enum.open && aToken != enum.close)
#			{
#				aToken <- paste(aToken,aChar,sep="")
#				aChar <- readChar(file, 1)
#			}
#		}
#		
#		return(list(token=aToken, lastChar=aChar))
#	}
#	
#	# Check function parameters and open the file if necessary
#	if(is.character(file)) {
#		file <- file(file,"r")
#		on.exit(close(file))
#	}
#	if (!inherits(file,"connection"))
#		stop("argument `file' must be a character string or connection")
#	
#	if (!isOpen(file)) {
#		open(file,"r")
#		on.exit(close(file))
#	}
#	
#	# Read in the header information
#	# Read the @relation line, first token should be @relation,
#	# otherwise an error
#	res <- nextArffToken(file, "")
#	if (tolower(res$token) != "@relation") {
#		stop("@relation is the first token expected, got: ", aToken)
#	}
#	res <- nextArffToken(file, res$lastChar)
#	relation <- res$token
#	
#	# Read in attributes
#	res <- nextArffToken(file, res$lastChar)
#	colClasses <- c()
#	col.names <- c()
#	while(tolower(res$token) == "@attribute") {
#		res <- nextArffToken(file, res$lastChar)
#		aName <- res$token
#		col.names <- append(col.names, aName)
#		res <- nextArffToken(file, res$lastChar)
#		aType <- tolower(res$token)
#		if (aType %in% c("numeric", "real", "integer")) {
#			# The value is a numeric type
#			colClasses <- append(colClasses, "numeric")
#		} else if (aType %in% c("string")) {
#			# The value is a character string type
#			colClasses <- append(colClasses, "character")
#		} else if (aType %in% c("date")) {
#			stop("ARFF date type is not supported")
#		} else if (aType == enum.open) {
#			# nominal type read until enum.close
#			colClasses <- append(colClasses, NA)
#			res <- nextArffToken(file, res$lastChar)
#			while (res$token != enum.close && res$token != FALSE) {
#				res <- nextArffToken(file, res$lastChar)
#			}
#		} else {
#			stop("invalid data type")
#		}
#		res <- nextArffToken(file, res$lastChar)
#	}
#	
#	# Read the data
#	if (tolower(res$token) != "@data") {
#		stop("@data expected, got: ", aToken)
#	}
#	
#	table <- read.table(file, na.strings="?", sep=",", 
#			comment.char=comment.char, col.names=col.names,
#			colClasses=colClasses)
#	
#	return(table)
#}
#
#
## Output to ARFF file format for use with Weka
## Version 1 2005-05-05
## Nigel Sim <nigel.sim@jcu.edu.au>
## 
## Karl Young <Karl.Young at radiology dot ucsf dot edu>
## added support for missing values, 2006-03-16
##
#write.arff <- function (x, file = "data", ncolumns = if (is.charactor(x)) 1 else 5, append = FALSE, class_col = NULL){
#	if (append == TRUE){
#		stop("Append not yet supported")
#	}
#	
#	if(is.character(file)) {
#		file <- file(file,"w")
#		on.exit(close(file))
#	}
#	
#	if (!inherits(file,"connection"))
#		stop("argument `file' must be a character string or connection")
#	
#	if (!isOpen(file)) {
#		open(file,"r")
#		on.exit(close(file))
#	}
#	
#	# Write header
#	writeLines("% Output from R",file)
#	writeLines("@relation ROutput",file)
#	
#	colNames <- dimnames(x)[[2]]
#	for (colName in colNames){
#		if (is.factor(x[1,colName])){
#			values <- unique(x[colName])[,1]
#			writeLines(paste("@attribute '", colName, "' {", paste(values,collapse=","), "}", sep=""),file)
#		} else {
#			colType <- "REAL"
#			writeLines(paste("@attribute '", colName, "' ", colType, sep=""),file)
#		}
#	}
#	
#	writeLines("@data", file)
#	# Write data
#	rows <- c()
#	for (i in 1:nrow(x)) {
#		# FIXME This messy code is to ensure that the factor types stay as factors
#		row <- ""
#		for (c in x[i,]){
#			if(is.na(c)) { c <- "?" } # Added Line For Handling Missing Values
#			if (row == "")
#				row <- c
#			else
#				row <- paste(row,c,sep=",")
#		}
#		rows <- c(rows, row)
#		if (length(rows) > 50){
#			writeLines(rows, file)
#			rows <- c()
#		}
#	}
#	
#	if (length(rows) > 0){
#		writeLines(rows, file)
#		rows <- c()
#	}
#	
#}
