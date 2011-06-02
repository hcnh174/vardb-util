
predictSvr <- function(data, fields, fit)
{
	data <- subsetNA(data,c(fields,'response','svr'))
	for (id in row.names(data))
	{
		y <- predict(fit, type='fitted.ind', newdata=data[id,])
		#y2 <- predictVr(data[id,],fit) #
		#print(y)
		svr <- data[id,'svr']
		predSVR <- ifelse(y<0.5,0,1)
		correct <- ifelse(svr==predSVR, 1, 0)
		print(paste('svr',svr,'y',y,'predSVR',predSVR,'correct',correct,'response',data[id,'response']))
		data[id,'y'] <- y
		#data[id,'y2'] <- y2
		data[id,'predSVR'] <- predSVR 
		data[id,'correct'] <- correct
	}
	print(length(which(data$correct==1))/length(data$correct))
	#data[,splitFields('response,correct')]
	data <- data[order(data$correct,data$response),splitFields('response,predSVR,correct')]
	
	return(data)
}