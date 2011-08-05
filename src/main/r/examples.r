#installation
install.packages("ggplot2") 

library()

#dataframes
attributes(someframe)

#plotting

source("http://research.stowers-institute.org/efg/R/Color/Chart/ColorChart.R")

lines(lowess(wt,mpg), col='blue', lwd=2, lty=2)

library(car)
scatterplot(mpg~wt)
scatterplot(mpg~wt | cyl, data=mtcars, lwd=2, legend.plot=T, labels=row.names(mtcars))
scatterplot(response ~ predictor, data=data, xlab="X", ylab="Y")


pairs(~mpg+disp+drat+wt,data=mtcars)
scatterplot.matrix(~mpg+disp+drat+wt,data=mtcars)
scatterplot.matrix(~mpg+disp+drat+wt|cyl,data=mtcars, diagonal='histogram')
scatterplot.matrix(~mpg+disp+drat+wt|cyl,data=mtcars, by.groups=T, diagonal='histogram')



rainbow(4)
library(corrgram)
corrgram(mtcars,order=TRUE,lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt)

library(vcd)
mosaic(Titanic, shade=T, legend=T)


capture.output(anova(data.glm),file="output.txt")


library(arm)
plot(data$svr, fitted(fit), type='n')
curve(invlogit(coef(fit)[1]+coef(fit)[2]*x), add=TRUE)

library(ggplot2)
qplot(trait,data=mydat,facets=geno~.)

sessionInfo()

#regression diagnostics
library(car)
fit <- lm(mpg~disp+hp+wt+drat, data=mtcars) 
qqPlot(fit, main="QQ Plot")
leveragePlots(fit, ask=FALSE) # leverage plots 
av.plots(fit, one.page=TRUE, ask=FALSE)
cutoff <- 4/((nrow(mtcars)-length(fit$coefficients)-2))
plot(fit, which=4, cook.levels=cutoff)
# Influence Plot
influencePlot(fit, main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )
# Evaluate homoscedasticity
# non-constant error variance test
ncvTest(fit)
# plot studentized residuals vs. fitted values
spreadLevelPlot(fit)
vif(fit) # variance inflation factors
sqrt(vif(fit)) > 2 # problem?
# Evaluate Nonlinearity
# component + residual plot
crPlots(fit, one.page=TRUE, ask=FALSE)
# Ceres plots
ceresPlots(fit, one.page=TRUE, ask=FALSE)
# Test for Autocorrelated Errors
durbinWatsonTest(fit)
# Global test of model assumptions
library(gvlma)
gvmodel <- gvlma(fit)
summary(gvmodel)


#R2 shrinkage
library(bootstrap)
theta.fit <- function(x,y){lsfit(x,y)}
theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}
kawaoka.subset <- na.omit(kawaoka.2b)
# matrix of predictors
X <- as.matrix(kawaoka.subset[c('age','sex','rs12979860CC','rs8099917TT','rs12980275AA','titerlog0')])
# vector of predicted values
y <- as.matrix(kawaoka.subset[c('titerlog4')])
results <- crossval(X,y,theta.fit,theta.predict,ngroup=10)
cor(y, fit$fitted.values)^2 # raw R2
cor(y,results$cv.fit)^2 # cross-validated R2



#relative importance
fit <- lm(titerlog4 ~ age + rs8099917TT + titerlog0, data=kawaoka.2b.subset)
calc.relimp(fit, type=c('lmg','last','first','pratt'),rela=T)
boot <- boot.relimp(fit,b=1000, type=c('lmg','last','first','pratt'), rank=T, diff=T, rela=T)
booteval.relimp(boot)
plot(booteval.relimp(boot,sort=T))


#memory usage
gc(verbose=T)# garbage collection
rm(list=ls()) # remove all objects from memory
.ls.objects() # list all objects by memory
round(memory.limit()/1048576.0, 2) # see the current memory limit in Mb
memory.profile() 
objects()
object.size(data0)
storage.mode(data0)


#R generic methods
methods('summary')
getAnywhere(summary) # finds where a function is defined regardless of package

#correspondence analysis
library(ca)
mytable <- with(data4@phdata, table(svr,sex)) # create a 2 way table
prop.table(mytable, 1) # row percentages
prop.table(mytable, 2) # column percentages
fit <- ca(mytable)
print(fit) # basic results
summary(fit) # extended results
plot(fit) # symmetric map
plot(fit, mass = TRUE, contrib = "absolute", map ="rowgreen", arrows = c(FALSE, TRUE)) # asymmetric map 

#rename something
names(Z)[names(Z)=='smoking'] <- 'smoke'

#display pch and color codes
show.pch()
show.col()
character.table()

#help - show overview information
help(Overview, library='Hmisc')
apropos("foo", mode="function")


with(mtcars, {
	summary(mpg, disp, wt)
	plot(mpg, disp)
	plot(mpg, wt)
})


status <- factor(status, ordered=TRUE)

class(data)
methods("print")
methods(class='lm')

define classes: setClass()
create objects: new()
define generics: setGeneric()
define methods: setMethods()
convert objects: as(), setAs()
check object validity: setValidity(), validObject()
access registry: showClass(), showMethods(), getMethod()


showClass("pixmap")
showMethods("show")
getMethod("show", "pixmap")

df[grep("trna", df$common_name, ignore.case=T),]


##########################################

# from R.bloggers 7/15/2011
library(affydata)
data(Dilution)
Dilution
phenoData(Dilution)
pData(Dilution)

# first plot
boxplot(Dilution,col=c(2,2,3,3))
##pick only a few genes to reduce calculation time
gn <- sample(geneNames(Dilution),100)
pms <- pm(Dilution[,3:4], gn)
mva.pairs(pms)

#normalization
normalized.Dilution &lt;- Biobase::combine(normalize(Dilution[, 1:2]),
		normalize(Dilution[, 3:4]))
normalize.methods(Dilution)

#second plot
boxplot(normalized.Dilution, col=c(2,2,3,3), main="Normalized Arrays")
pms <- pm(normalized.Dilution[, 3:4],gn)
mva.pairs(pms)

#######################################################################

if (!require('RJSONIO')) install.packages('RJSONIO', repos = 'http://cran.r-project.org')
if (!require('png')) install.packages('png', repos = 'http://cran.r-project.org')
if (!require('ReadImages')) install.packages('ReadImages', repos = 'http://cran.r-project.org')
install.packages("RXKCD", repos="http://R-Forge.R-project.org")


suppressWarnings(as.integer(c(input.list, recursive=TRUE)))
sprintf("%s and %s","value1", "value2")