library(qcc)

#
#COMPUTE LAMBDA
#

getLambda = function(ELAI, W){
	n = length(ELAI)
	if( n>W ){
        	#compute sum of squared forecasting errors     
        	ssError = function(L, y, W){
			#compute ewma
        	        ewmaOut = ewmaSmooth(1:W, y[1:W], lambda=L)
			#sum squared forecasting error
			return( sum((ewmaOut$y[2:W] - y[1:(W-1)])^2) )
        	}
		#minimize sum of squared forcasting errors
        	out = optimize(ssError, c(0, 1), ELAI, W)$minimum
	} else{ out = 0.2 } #return default for lambda until control window filled
        return( out )
}

#
#CHECK CONVERGENCE
#

isConverge = function(ELAI, L, W){	
	n = length(ELAI)
	#first collect control window
	if( n>W ){
		#compute ewma and limits
        	ewmaOut = ewma(ELAI[1:W], lambda=L, newdata=ELAI[W:n], plot=F)
		#find points outside of control limits
		isOut = ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1]
		#no violations in window & violations beyond the control window
		cFlag = !any(isOut[1:W]) & any(isOut[W+1:length(isOut)])
	} else{ cFlag = FALSE } #no convergence until control window is filled
	return( cFlag )
}
