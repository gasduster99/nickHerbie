rm(list=ls())
#
#EWMA CONVERGENCE OF ROSENBROCK
#

#
library(tgp)
library(qcc)

#rosenbrock test function
f = function(x){ 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2 }

#
#EWMA CONVERGENCE CHART
#

#window size
W = 30
#a flag to monitor convergence
notConverge = TRUE
#a function to compute lambda
getLambda = function(ELAI, W){
	#minimize sum of squared forecasting errors	
	ssError = function(L, y, W){
		yLen = length(y)
		ewmaOut = ewma(y[1:W], lambda=L, newdata=y[W:yLen], plot=F)	
		#sum((ewmaOut$y[1:yLen]-y)^2)
		return( sum((ewmaOut$y[1:yLen+1]-y)^2) )
	}
	out = optimize(ssError, c(0, 1), rev(ELAI), W)
	return( out$minimum )
}
#a function for checking convergence
ewmaCC = function(ELAI, L, W){	
	#
	ewmaOut = ewma(rev(EI)[1:W], lambda=L, newdata=rev(EI)[W:length(EI)], plot=F)
	#return a list
		#1)convergence state
			#out = ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1]
		#2)L
}


#
#TGP SURROGATE
#

#define search domain
rect = cbind(c(-2, -3), c(2, 5))
#1) Collect an initial set, X, from the domain.
X = lhs(40, rect)
#2) Compute f(X).
Z = f(X)
Zmax = c(min(Z))
#init
elai = c()
tgpOut = NULL
#
#while( notConverge ){
	#3-5) Fit Surrogate; Collect candidate set; Compute EI; 
	tgpOut = optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=tgpOut, improv=c(1,1), trace=T, verb=0)	
	ex = matrix(tgpOut$X, ncol=2)	
	#6) Add argmax EI to X.
        X = rbind(X, ex)#
        Z = c(Z, f(ex))
        Zmax = c(Zmax, min(Z))
	#
	#maxI = which( EimprovAll$rank==1 )
	samp = tgpOut$obj$trace$preds$improv[ tgpOut$obj$improv$rank==1 ][[1]]
	ei = mean(samp)
	vi = var(samp)
	#
	elai = c(elai, log((ei^2)/(vi+ei^2)^0.5))
	
        #EimprovAll = tgpOut$obj$improv 
        #maxSamples = unlist( improvSamples[maxI] )
        #m = mean(maxSamples)
        maxes[it] = m
#	#7) Check Convergence.
#	notConverge = ewmaCC(ei, W)
#}


#
#EM = dim(rect)[1]
##sweeps of optimization
#M = 60#150#120
##intializing
#Zmax = c(min(Z))
#wSamples = c()
#samples = c()
#sLen = NULL
#maxes = matrix(NaN, nrow=M, ncol=1)
#out = NULL





