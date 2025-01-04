rm(list=ls())

#Using the gaussian process models that fit using RJ-mcmc in the tgp package this file make run for an extended period of time.  

#
#EWMA CONVERGENCE OF ROSENBROCK
#

#
library(tgp)
library(qcc)

#rosenbrock test function
rosenbrock = function(x){ 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2 }

#create and objective function
f = rosenbrock 

#
#EWMA CONVERGENCE CHART
#

#window size
W = 30

#a function to compute lambda
getLambda = function(ELAI, W){
	#minimize sum of squared forecasting errors	
	ssError = function(L, y, W){
		yLen = length(y)
		ewmaOut = ewma(y[1:W], lambda=L, newdata=y[W:yLen], plot=F)	
		return( sum((ewmaOut$y[1:yLen+1]-y)^2) )
	}
	out = optimize(ssError, c(0, 1), rev(ELAI), W)
	return( out$minimum )
}

#a function for checking convergence
ewmaCC = function(ELAI, L, W){	
	#compute the ewma chart on the control window
	ewmaOut = ewma(rev(ELAI)[1:W], lambda=L, newdata=rev(ELAI)[(W+1):length(ELAI)], plot=F)
	#check elai against the limits
        allLimitCheck = rev(ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1])
        A = length(allLimitCheck)
        outWindow = allLimitCheck[1:(A-W)]   #outside of control window
        innWindow = allLimitCheck[(A-W+1):A] #inside the control window
        #
        return( any(outWindow) & !any(innWindow) )
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
Zmin = c(min(Z))

#initialize
i = 1
elai = c()
tgpOut = NULL
isConverged = F
while( !isConverged ){
	#3-5) Fit Surrogate; Collect candidate set; Compute EI; 
	tgpOut = optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=tgpOut, improv=c(1,1), trace=T, verb=0)	
	ex = matrix(tgpOut$X, ncol=2)	
	#6) Add argmax EI to X.
        X = rbind(X, ex)#
        Z = c(Z, f(ex))
        Zmin = c(Zmin, min(Z))
	
	#compute ELAI
	samp = tgpOut$obj$trace$preds$improv[ tgpOut$obj$improv$rank==1 ][[1]]
	ei = mean(samp)
	vi = var(samp)
	elai = c(elai, log((ei^2)/(vi+ei^2)^0.5) )

	#7) Check Convergence	
	#only check after window is filled
	if(i>W){	
		#compute lambda
		L = getLambda(elai, W)
		#convergence test
		isConverged = ewmaCC(elai, L, W)
	}
	
	#index iteration	
	i = i+1
}


#
#CONVERGENCE PLOTS
#

#
#a function for plotting the EWMA chart
ewmaPlot = function(ELAI, L, W){
        #compute the ewma chart on the control window
        ewmaOut = ewma(rev(ELAI)[1:W], lambda=L, newdata=rev(ELAI)[(W+1):length(ELAI)], plot=F)
        #check elai against the limits
        allLimitCheck = rev(ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1])
        A = length(allLimitCheck)
        outWindow = allLimitCheck[1:(A-W)]   #outside of control window
        innWindow = allLimitCheck[(A-W+1):A] #inside the control window
        #	
	ys  = rev(c(ewmaOut$statistics, ewmaOut$newstats))
	its = 1:length(ys)
	#
	plot(its, ys,
	        pch=3,
	        xlab="Iteration Number",
	        ylab="Summary Statistics",
	        main="EWMA Convergence Chart"
	)
        points(ewmaOut$x[!allLimitCheck], rev(ewmaOut$y)[!allLimitCheck], pch=20)
	points(ewmaOut$x[allLimitCheck], rev(ewmaOut$y)[allLimitCheck], col='red')
	lines(ewmaOut$x, rev(ewmaOut$y))	
	matplot(matrix(rev(ewmaOut$limits), ncol=2), type='l', lty=2, col='black', add=T)
	abline(v=(length(ys)+0.5)-W, lty=3)
	legend( "topright", c(expression(Y[i]), expression(Z[i]*" Violation"), expression(Z[i]*" Control"), "Control Limit", NA,"\nControl Window\nBoundary"),
		pch=c(3 , 1, 20, NA, NA, NA),
		lty=c(NA, NA, NA, 2, 3, NA),
		col=c("black", "red", "black", "black", "black", NA)
	)
}

#
ewmaPlot(elai, L, W)

#
dev.new()
plot( seq(1, length(Zmin)), Zmin,
	"l",
	xlab="Iteration Number",
	ylab="Objective Value",
	main="Best Z Value",
	ylim=c(0, max(Zmin))
)
abline( h=0, lty=3 )
legend("left", c("\nCurrent Best\n Z Value", "\nTheoretical\n Minimum", NA),
        lty=c(1, 3, NA)
)


