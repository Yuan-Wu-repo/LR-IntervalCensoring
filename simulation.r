library(MASS)
library(pracma)
library(splines)
library(survival)
library(ICsurv)

source("functions.R")

n.variable <- 3
size <- 50
n.spline <- 6
run <- 1000
lambda<- 1
spline.ord <- 4
spline.degree <- 3

etamatrix <- matrix(0,(n.variable+n.spline),run)
estnaive <- matrix(0, n.variable, run)

senaive<-matrix(0, n.variable, run)
knotmatrix<-matrix(0,(n.spline + spline.degree),run)
sematrix <- matrix(0, n.variable, run)


cover.wald1 <- 0
cover.wald2 <- 0
cover.wald3 <- 0
cover.wald.whole <- 0

cover.lr1 <- 0
cover.lr2 <- 0
cover.lr3 <- 0
cover.lr.whole <- 0

cover1n <- 0
cover2n <- 0
cover3n <- 0
coverwholen <- 0

cover.lr1n <- 0
cover.lr2n <- 0
cover.lr3n <- 0
cover.lrn.whole <- 0


icsurvse <- matrix(0, n.variable, run)

cover.icsurv1 <- 0
cover.icsurv2 <- 0
cover.icsurv3 <- 0
cover.icsurv.whole <- 0


for (run_index in 1 : run)	{
	
	cat("\n", "\n", "\n", "run_index=",run_index) 


	z1<-runif(size, 0, 1)		  
	z2<-rnorm(size, 0, 1)
	z3<-rbinom(size,1,0.5)
	z<-rbind(z1,z2,z3)
	
	beta0<-matrix(c(-1, 0.5, 1.5), n.variable, 1)

	u<-runif(size, 0 ,1)

	# simulating proportional hazard survival data with exponential baseline survival function
	# \Lamda(t)=t^{1/2} S(t)=exp(-t^{1/2}exp(\beta z))
	t<-(log(u) * (-1 / lambda) * c(exp(-t(beta0) %*% z)))^2
	
	delta1 <- rep(0, size)		

	delta2 <- rep(0, size)

	delta3 <- rep(0, size)		  

	ctu <- rep(0, size)

	ctv <- rep(0 ,size)

	for (ind in 1 : size) {
	
	preuv <- rexp(2, 0.5)
	if (t[ind] <= 5)
	{
		if (t[ind] <= min(preuv)) {
			delta1[ind] <- 1
			ctu[ind] <- min(min(preuv), 5)
			ctv[ind] <- ctu[ind] + 10^(-8)
		} else {
			if (t[ind] > min(preuv) & t[ind] <= max(preuv)) {
			delta2[ind] <- 1
			ctu[ind] <- max(preuv[preuv < t[ind]])
			ctv[ind] <- min(min(preuv[preuv >= t[ind]]), 5)
			} else {										
				if (t[ind] > max(preuv)) {
					delta3[ind] <- 1
					ctv[ind] <- max(preuv)
					ctu[ind] <- ctv[ind] - 10^(-8)
				}
			}			
		}	
	} else {
		delta3[ind] <- 1
		ctv[ind] <- min(max(preuv), 5)
		ctu[ind] <- ctv[ind] - 10^(-8)
	}	

	}	


	
	
	###################################naive method using survival package###########################
	newt <- rep(0, size)
	newt[delta1 == 1] <- ctu[delta1 == 1] / 2
	newt[delta2 == 1] <- (ctu[delta2 == 1] + ctv[delta2 == 1]) / 2
	newt[delta3 == 1] <- ctv[delta3 == 1]
	newdelta <- rep(0, size)
	newdelta[delta1 == 1 | delta2 == 1] <- 1

	fitnaive <- coxph(Surv(newt, newdelta) ~ z[1,] + z[2,] + z[3,], method="breslow", init = c(0, 0, 0))
	### to avoid diag(fitnaive$var) has negative components, when it is negative let it be 0
	diagnaive <- apply(rbind(diag(fitnaive$var), rep(0, n.variable)), 2, max)
			
	fitcoef <- c(fitnaive$coef)

	if ((fitcoef[1] + qnorm(0.975, 0, 1) * sqrt(diagnaive)[1] >= 0) &
		(fitcoef[1] - qnorm(0.975, 0, 1) * sqrt(diagnaive)[1] <= 0))				
		cover1n <- cover1n + 1
		
	if ((fitcoef[2] + qnorm(0.975, 0, 1) * sqrt(diagnaive)[2] >= 0) &
		(fitcoef[2] - qnorm(0.975, 0, 1) * sqrt(diagnaive)[2] <= 0))				
		cover2n <- cover2n + 1
		
	if ((fitcoef[3] + qnorm(0.975, 0, 1) * sqrt(diagnaive)[3] >= 0) &
		(fitcoef[3] - qnorm(0.975, 0, 1) * sqrt(diagnaive)[3] <= 0))				
		cover3n <- cover3n + 1	

	# whole Wald test based on naive converting
	if(pchisq(q = t(fitcoef) %*% ginv(fitnaive$var) %*% 
				  fitcoef, df = 3, lower.tail = FALSE) > 0.05)
		coverwholen <- coverwholen + 1 

	estnaive[, run_index] <- fitcoef
	senaive[, run_index] <- sqrt(diagnaive)

	### likelihood ratio test for naive method
	# likelihood ratio test based on naive converting for beta=0	
	if(summary(fitnaive)$logtest[3] > 0.05)
		cover.lrn.whole <- cover.lrn.whole + 1 	
		
	### note that loglik[2] is the loglikelihood for final coefficients	
	loglikefull.naive <- fitnaive$loglik[2]

	### partial coeffcients cox regression
	fitnaive <- coxph(Surv(newt, newdelta) ~  z[2,] + z[3,], method="breslow")

	loglikepart1.naive <- fitnaive$loglik[2]

	fitnaive <- coxph(Surv(newt, newdelta) ~  z[1,] + z[3,], method="breslow")

	loglikepart2.naive <- fitnaive$loglik[2]

	fitnaive <- coxph(Surv(newt, newdelta) ~  z[1,] + z[2,], method="breslow")

	loglikepart3.naive <- fitnaive$loglik[2]

	## likelihood ratio test for three parameters for naive method
	if(pchisq(q = 2 * (loglikefull.naive - loglikepart1.naive), df = 1, lower.tail = FALSE) > 0.05)
		cover.lr1n <- cover.lr1n + 1
		
	if(pchisq(q = 2 * (loglikefull.naive - loglikepart2.naive), df = 1, lower.tail = FALSE) > 0.05)
		cover.lr2n <- cover.lr2n + 1
		
	if(pchisq(q = 2 * (loglikefull.naive - loglikepart3.naive), df = 1, lower.tail = FALSE) > 0.05)
		cover.lr3n <- cover.lr3n + 1
	############## end naive method ################################################



	#### start sieve interval censoring method #############
	knotb <- get_knots(ctu, ctv, delta1, delta2, delta3)

	### bspline results		 
	bspu <- t(splineDesign(knots = knotb, x = ctu, ord = 4))
	### from bspline to ispline  
	pre.ispu <- apply(bspu, 2, rev)

	ispu <- apply(pre.ispu, 2, cumsum)

	ispu <- apply(ispu, 2 ,rev)

	ispu <- ispu[-1,]

	### bspline results
	bspv <- t(splineDesign(knots = knotb, x = ctv, ord = 4))
	### from bspline to ispline  
	pre.ispv <- apply(bspv, 2, rev)

	ispv <- apply(pre.ispv, 2, cumsum)

	ispv <- apply(ispv, 2, rev)

	ispv <- ispv[-1,]

	## create vecters of 1s with length of the numbers of left, interval and right censored###
	number.ones1 <- matrix(1, 1, sum(delta1))
	number.ones2 <- matrix(1, 1, sum(delta2))
	number.ones3 <- matrix(1, 1, sum(delta3))
	
	oldz <- z	  
	  
	index <- 0			
	  
	  if (index == 0) {
		 z <- oldz
		 fixsubx <- matrix(0, 1, size)
		 n.total <- n.spline + n.variable		 
	  }	 

	############################################################################
	############## unconstrained R function for finding minimum point ##########
	############################################################################		
	out <- optim(par = rep(0, (n.variable + n.spline)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))

	etamatrix[,run_index]<-out$par
			
	###############################################################################################################
	#### variance based on theory from Huang et al. (2008) or Zhang et al. (2010)'s least square approach ##########
	hat_O <- observ.beta(out$par) - observ.beta.lambda(out$par) %*% ginv(observ.lambda(out$par)) %*% t(observ.beta.lambda(out$par))


	### to avoid diag(ginv(hat_O)) has negative components, when it is negative let it be 0
	diag.inv <- apply(rbind(diag(ginv(hat_O)), rep(0, n.variable)), 2, max)
	se<-sqrt(diag.inv / size)	
	sematrix[, run_index] <- se

	if ((etamatrix[(n.spline + 1), run_index] + qnorm(0.975, 0, 1) * sematrix[1, run_index] >= 0)&
		(etamatrix[(n.spline + 1), run_index] - qnorm(0.975, 0, 1) * sematrix[1, run_index] <= 0))
		cover.wald1 <- cover.wald1 + 1	

	if ((etamatrix[(n.spline + 2), run_index] + qnorm(0.975, 0, 1) * sematrix[2, run_index] >= 0)&
		(etamatrix[(n.spline + 2), run_index] - qnorm(0.975, 0, 1) * sematrix[2, run_index] <= 0))
		cover.wald2 <- cover.wald2 + 1	
		
	if ((etamatrix[(n.spline + 3), run_index] + qnorm(0.975, 0, 1) * sematrix[3, run_index] >= 0)&
		(etamatrix[(n.spline + 3), run_index] - qnorm(0.975, 0, 1) * sematrix[3, run_index] <= 0))
		cover.wald3 <- cover.wald3 + 1	
		
	## wald test for beta=0
	if(pchisq(q = t(etamatrix[(n.spline + 1) : (n.spline + 3), run_index]) %*% hat_O %*% 
				  etamatrix[(n.spline + 1) : (n.spline + 3), run_index] * size, df = 3, lower.tail = FALSE) > 0.05)
		cover.wald.whole <- cover.wald.whole + 1  	

	#### end for wald test based sieve interval censoring variance method 1 (Huang et al. 2008 or Zhang et al 2010) #############

	############## starting using ICsurv package for variance estimation ##################################
	d1<-rep(0, size)

	d2<-rep(0, size)

	d3<-rep(0, size)

	d1[delta1 == 1] <- 1

	d2[delta2 == 1] <- 1

	d3[delta3 == 1] <- 1

	Li <- rep(0, size)

	Ri <- rep(0, size)

	Li[delta1 == 1] <- 0

	Li[delta2 == 1] <- ctu[delta2 == 1]

	Li[delta3 == 1] <- ctv[delta3 == 1]

	Ri[delta1 == 1] <- ctu[delta1 == 1]

	Ri[delta2 == 1] <- ctv[delta2 == 1]

	Ri[delta3 == 1] <- 0

	Xp <- t(z)

	isLi <- t(splineDesign(knots = knotb, x = Li, ord = 4))

	isLi <- apply(isLi, 2, rev)

	isLi <- apply(isLi, 2, cumsum)

	isLi <- apply(isLi, 2 ,rev)

	bLi <- t(isLi[-1,])

	bLi[bLi == 0] <- 10^(-10)

	isRi <- t(splineDesign(knots = knotb, x = Ri, ord = 4))

	isRi <- apply(isRi, 2, rev)

	isRi <- apply(isRi, 2, cumsum)

	isRi <- apply(isRi, 2 ,rev)

	bRi <- t(isRi[-1,])

	bRi[bRi == 0] <- 10^(-10)

	b1 <- etamatrix[(n.spline + 1) : (n.spline + 3), run_index]

	g1 <- exp(etamatrix[1 : n.spline, run_index])

	v <- PH.Louis.ICsurv(b1, g1, bLi, bRi, d1, d2, d3, Xp)

	A <- v[1 : n.variable, 1 : n.variable]
		
	B <- v[1 : n.variable, (n.variable + 1) : (n.variable + n.spline)]

	C <- v[(n.variable + 1) : (n.variable + n.spline), 1 : n.variable]

	D <- v[(n.variable + 1) : (n.variable + n.spline), (n.variable + 1) : (n.variable + n.spline)]
					
	var.b = ginv(A - B %*% ginv(D) %*% C)

	### to avoid diag(fitnaive$var) has negative components, when it is negative let it be 0
	icsurvdiag <- apply(rbind(diag(var.b), rep(0, n.variable)), 2, max)

	icsurvse[, run_index] <- sqrt(icsurvdiag)

		
	if ((etamatrix[(n.spline + 1), run_index] + qnorm(0.975, 0, 1) * icsurvse[1, run_index] >= 0)&
		(etamatrix[(n.spline + 1), run_index] - qnorm(0.975, 0, 1) * icsurvse[1, run_index] <= 0))
		cover.icsurv1 <- cover.icsurv1 + 1	

	if ((etamatrix[(n.spline + 2), run_index] + qnorm(0.975, 0, 1) * icsurvse[2, run_index] >= 0)&
		(etamatrix[(n.spline + 2), run_index] - qnorm(0.975, 0, 1) * icsurvse[2, run_index] <= 0))
		cover.icsurv2 <- cover.icsurv2 + 1	
		
	if ((etamatrix[(n.spline + 3), run_index] + qnorm(0.975, 0, 1) * icsurvse[3, run_index] >= 0)&
		(etamatrix[(n.spline + 3), run_index] - qnorm(0.975, 0, 1) * icsurvse[3, run_index] <= 0))
		cover.icsurv3 <- cover.icsurv3 + 1	

	## wald test for beta=0
	if(pchisq(q = t(etamatrix[(n.spline + 1) : (n.spline + 3), run_index]) %*% (A - B %*% ginv(D) %*% C) %*% 
				  etamatrix[(n.spline + 1) : (n.spline + 3), run_index], df = 3, lower.tail = FALSE) > 0.05)
		cover.icsurv.whole <- cover.icsurv.whole + 1  	
	######### end for wald test approach with icsurv package variance #################################################


	######################### likelihood ratio approach for sieve estimation ###############################
	## loglikelihood for full model
	loglikefull <- loglike(out$par)
	 
	index <- 1
	  
	  if (index != 0) {
		 z <- oldz[-index,]
		 fixsubx <- matrix(c(0*oldz[index,]), 1, size)
		 n.total <- n.spline + n.variable -1			 
	  }	 

	out <- optim(par = rep(0, (n.variable + n.spline -1)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))

	## loglikelihood for model for test beta1
	loglikepart1 <- loglike(out$par)

	index <- 2

	  if (index != 0) {
		 z <- oldz[-index,]
		 fixsubx <- matrix(c(0*oldz[index,]), 1, size)
		 n.total <- n.spline + n.variable -1			 
	  }	 

	out <- optim(par = rep(0, (n.variable + n.spline -1)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))

	## loglikelihood for model for test beta2
	loglikepart2 <- loglike(out$par)

	index <- 3
	  
	  if (index != 0) {
		 z <- oldz[-index,]
		 fixsubx <- matrix(c(0*oldz[index,]), 1, size)
		 n.total <- n.spline + n.variable -1			 
	  }	 

	out <- optim(par = rep(0, (n.variable + n.spline -1)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))

	## loglikelihood for model for test beta2
	loglikepart3 <- loglike(out$par)

	## likelihood ratio test for three parameters 
	if(pchisq(q = -2 * (loglikefull - loglikepart1), df = 1, lower.tail = FALSE) > 0.05)
		cover.lr1 <- cover.lr1 + 1
		
	if(pchisq(q = -2 * (loglikefull - loglikepart2), df = 1, lower.tail = FALSE) > 0.05)
		cover.lr2 <- cover.lr2 + 1
		
	if(pchisq(q = -2 * (loglikefull - loglikepart3), df = 1, lower.tail = FALSE) > 0.05)
		cover.lr3 <- cover.lr3 + 1	

	# loglikelihood for model all regression coefficients equal to 0
	fixsubx <- matrix(0, 1, 3) %*% oldz
	n.total <- n.spline		 
			
	out <- optim(par = rep(0, n.spline), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))

	loglikenull <- loglike(out$par)

	## likelihood ratio test for beta=0
	if(pchisq(q = -2 * (loglikefull - loglikenull), df = 3, lower.tail = FALSE) > 0.05)
		cover.lr.whole <- cover.lr.whole + 1 

	######### end for likelihood ratio test for sieve interval censoring #################################################	
			
}

eta.mean<-c(rowMeans(exp(etamatrix[1:n.spline,])),rowMeans(etamatrix[(n.spline+1):(n.spline+n.variable),]))
eta.mean

se.mean <- rowMeans(sematrix)
se.mean

se.ICsurv.mean <- rowMeans(icsurvse)
se.ICsurv.mean

estnaivemean <- rowMeans(estnaive)
estnaivemean

senaivemean <- rowMeans(senaive)
senaivemean

c(sd(etamatrix[(n.spline + 1),]),
  sd(etamatrix[(n.spline + 2),]),
  sd(etamatrix[(n.spline + 3),]))
  


run - c(cover.wald1, cover.wald2, cover.wald3)
run - cover.wald.whole

run - c(cover.icsurv1, cover.icsurv2, cover.icsurv3)
run - cover.icsurv.whole

run - c(cover.lr1, cover.lr2, cover.lr3)
run - cover.lr.whole

run - c(cover1n, cover2n, cover3n)
run - coverwholen

run - c(cover.lr1n, cover.lr2n, cover.lr3n)
run - cover.lrn.whole
