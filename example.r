rm(list=ls())
library(tgp)
source('ewmaConvChart.r')

#rosenbrock test function
f = function(x){ 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2 }
#window size
W = 30

#
#OPTIMIZATION
#

#define search domain
rect = cbind(c(-2, -3), c(2, 5))
#1) Collect an initial set, X, from the domain.
X = lhs(40, rect)
#2) Compute f(X).
Z = f(X); Zmin = min(Z)
#3-5) Fit Surrogate; Collect candidate set; Compute EI; 
tgpOut = optim.step.tgp(f, X=X, Z=Z, rect=rect, improv=c(1,1), trace=T, verb=0)
#6) Add argmax EI to X.
X = rbind(X, tgpOut$X)

#prepare for next iteration
Z = c(Z, f(tgpOut$X)); Zmin = c(Zmin, min(Z))
#maximum improvement samples
samp = tgpOut$obj$trace$preds$improv[ tgpOut$obj$improv$rank==1 ][[1]]
samples = list(); samples = c(samples, samp)
#EI & ELAI
ei   = mean(samp)
vi   = var(samp)
elai = log((ei^2)/(vi+ei^2)^0.5)

#continue with additional iterations until convergence is reached
while( !isConverge(rev(elai), getLambda(rev(elai), W), W) ){
	#3-5) Fit Surrogate; Collect candidate set; Compute EI; 
	tgpOut = optim.step.tgp(f, X=X, Z=Z, rect=rect, improv=c(1,1), trace=T, verb=0)
	#6) Add argmax EI to X.
	X = rbind(X, tgpOut$X)
	
	#prepare for next iteration
	Z = c(Z, f(tgpOut$X)); Zmin = c(Zmin, min(Z))
	#maximum improvement samples
	samp = tgpOut$obj$trace$preds$improv[ tgpOut$obj$improv$rank==1 ][[1]]
	samples = c(samples, samp)
	#EI & ELAI
	ei   = mean(samp)
	vi   = var(samp)
	elai = c( elai, log((ei^2)/(vi+ei^2)^0.5) )
	#realtime track convergence
	n = length(elai)	
	if(n>W){
		layout(t(c(1,2)))
		plot(Zmin, type='l')
	 	ewma(rev(elai)[1:W], lambda=getLambda(rev(elai), W), newdata=rev(elai)[W:n])
	}
}




#        Zmax = c(Zmax, min(Z))
#
#elai = c()
#tgpOut = NULL
#
##rec
#Zmax = c(min(Z))
##
##while( !isConverge ){
#	#3-5) Fit Surrogate; Collect candidate set; Compute EI; 
#	tgpOut = optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=tgpOut, improv=c(1,1), trace=T, verb=0)	
#	ex = matrix(tgpOut$X, ncol=2)	
#	#6) Add argmax EI to X.
#        X = rbind(X, ex)#
#        Z = c(Z, f(ex))
#        Zmax = c(Zmax, min(Z))
#	#
#	#maxI = which( EimprovAll$rank==1 )
#	samp = tgpOut$obj$trace$preds$improv[ tgpOut$obj$improv$rank==1 ][[1]]
#	ei = mean(samp)
#	vi = var(samp)
#	#
#	elai = c(elai, log((ei^2)/(vi+ei^2)^0.5))
#	
#        #EimprovAll = tgpOut$obj$improv 
#        #maxSamples = unlist( improvSamples[maxI] )
#        #m = mean(maxSamples)
#        #maxes[it] = m
##	#7) Check Convergence.
##	notConverge = ewmaCC(ei, W)
##}


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




