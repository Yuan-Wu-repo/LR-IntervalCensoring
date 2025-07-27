


############################ the log likelihood function #############
loglike <- function(eta)
{
	eta <- matrix(c(eta), (n.total),1)
	alpha<-eta[(1):(n.spline),]
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	
	if (n.total > n.spline)	{
		beta<-eta[(n.spline+1):(n.total),]
		subx<-exp(c(t(beta) %*% z + fixsubx))	
	} else
		subx <- exp(c(fixsubx))
		
	log1 <- sum(log(1 - exp(-subu[delta1==1] * subx[delta1==1])))
	
	log2 <- sum(log(exp(-subu[delta2==1] * subx[delta2==1])
				-exp(-subv[delta2==1] * subx[delta2==1])))
				
	log3 <- -sum(subv[delta3==1] * subx[delta3==1])
	
	res <- -log1 - log2 - log3
	
	res
}
########################################################################

################ derivative of the log likelihood function  ##############################################

gradient <- function(eta) 
{
    
	eta <- matrix(c(eta), (n.total),1)
	alpha<-eta[(1):(n.spline),]
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	
	if (n.total > n.spline)	{
		beta<-eta[(n.spline+1):(n.total),]
		subx<-exp(c(t(beta) %*% z + fixsubx))	
	} else
		subx <- exp(c(fixsubx))
		
	diag_alpha<-diag(c(exp(alpha)),n.spline)
	
		
	### derivative for alpha
	ggalpha11 <- t(diag_alpha %*% ispu[,delta1==1]) * 
	exp(-subu[delta1==1] * subx[delta1==1]) * subx[delta1==1]
	
	ggalpha12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	ggalpha1 <- number.ones1 %*% (ggalpha11 / ggalpha12)
		
	ggalpha21 <- -t(diag_alpha %*% ispu[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * subx[delta2==1]
	
	ggalpha22 <- t(diag_alpha %*% ispv[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * subx[delta2==1]
	
	ggalpha23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * subx[delta2==1]) 
	
	ggalpha2 <- number.ones2 %*% ((ggalpha21 + ggalpha22) / ggalpha23)
	
	ggalpha3 <- -number.ones3 %*% (t(diag_alpha %*% ispv[,delta3==1]) * subx[delta3==1])
	
	ggalpha <- ggalpha1 + ggalpha2 + ggalpha3
	
	### derivarive for beta
	
	if (n.total > n.spline) {
		ggbeta11 <- t(z[,delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) * 
		subu[delta1==1] * subx[delta1==1]
		
		ggbeta12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
		
		ggbeta1 <- number.ones1 %*% (ggbeta11 / ggbeta12)
		
		ggbeta21 <- -t(z[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * 
		subu[delta2==1] * subx[delta2==1]
		
		ggbeta22 <- t(z[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * 
		subv[delta2==1] * subx[delta2==1]
		
		ggbeta23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * subx[delta2==1])
		
		ggbeta2 <- number.ones2 %*% ((ggbeta21 + ggbeta22) / ggbeta23)
		
		ggbeta3 <- -number.ones3 %*% (t(z[,delta3==1]) * subv[delta3==1] * subx[delta3==1])
		
		ggbeta <- ggbeta1 + ggbeta2 + ggbeta3
		
		res <- -c(c(ggalpha), c(ggbeta))
	} else
		res <- -c(ggalpha)
		
	res
		
}
##############################################################################################################

##################### for variance estimation based on Zhang et al. (2010) ##############################
observ.beta <- function(eta)
{
	eta <- matrix(c(eta), (n.variable+n.spline),1)
	alpha<-eta[(1):(n.spline),]
	beta<-eta[(n.spline+1):(n.variable+n.spline),]
	
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	subx<-exp(c(t(beta)%*%z))	
	
	#part 1 beta observation
	
	obeta11 <- t(z[,delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) * 
	subu[delta1==1] * subx[delta1==1]
	
	obeta12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	obeta1 <- t((obeta11 / obeta12)) %*% (obeta11 / obeta12)
	
	#part 2 beta observation
		
	obeta21 <- -t(z[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * 
	subu[delta2==1] * subx[delta2==1]
	
	obeta22 <- t(z[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * 
	subv[delta2==1] * subx[delta2==1]
	
	obeta23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * subx[delta2==1])
	
	obeta2 <- t((obeta21 + obeta22) / obeta23) %*% ((obeta21 + obeta22) / obeta23)	
	
	# part 3 beta observation
	
	obeta31 <- -t(z[,delta3==1]) * subv[delta3==1] * subx[delta3==1]
	
	obeta3 <- t(obeta31) %*% obeta31
	
	# total beta observation
	
	res <- (obeta1 + obeta2 + obeta3) / size
	
	res
	
}


observ.beta.lambda <- function(eta)
{
	
	eta <- matrix(c(eta), (n.variable+n.spline),1)
	alpha<-eta[(1):(n.spline),]
	beta<-eta[(n.spline+1):(n.variable+n.spline),]
	
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	subx<-exp(c(t(beta)%*%z))	
	
	#part 1 beta cross Lambda observation
	
	obeta11 <- t(z[,delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) * 
	subu[delta1==1] * subx[delta1==1]
	
	obeta12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	olambda11 <- t(bspu[, delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) *
	subx[delta1==1]
	
	olambda12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	obeta.lambda1 <- t((obeta11 / obeta12)) %*% (olambda11 / olambda12)

	#part 2 beta cross Lambda observation

	obeta21 <- -t(z[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * 
	subu[delta2==1] * subx[delta2==1]
	
	obeta22 <- t(z[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * 
	subv[delta2==1] * subx[delta2==1]
	
	obeta23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * 
	subx[delta2==1])
	
	olambda21 <- -t(bspu[, delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda22 <- t(bspv[, delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * 
	subx[delta2==1])
	
	obeta.lambda2 <- t((obeta21 + obeta22) / obeta23) %*% ((olambda21 + olambda22) / olambda23)
	
	# part 3 beta cross Lambda observation
	
	obeta31 <- -t(z[,delta3==1]) * subv[delta3==1] * subx[delta3==1]
	
	olambda31 <- -t(bspv[,delta3==1]) * subx[delta3==1]
	
	obeta.lambda3 <- t(obeta31) %*% olambda31
	
	# total beta cross Lambda observation
	
	res <- (obeta.lambda1 + obeta.lambda2 + obeta.lambda3) / size
	
	res
	
}		

observ.lambda <- function(eta)
{
	
	eta <- matrix(c(eta), (n.variable+n.spline),1)
	alpha<-eta[(1):(n.spline),]
	beta<-eta[(n.spline+1):(n.variable+n.spline),]
	
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	subx<-exp(c(t(beta)%*%z))	
	
	#part 1  Lambda observation
	
	olambda11 <- t(bspu[, delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) *
	subx[delta1==1]
	
	olambda12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	olambda1 <- t((olambda11 / olambda12)) %*% (olambda11 / olambda12)

	#part 2  Lambda observation

	olambda21 <- -t(bspu[, delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda22 <- t(bspv[, delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * 
	subx[delta2==1])
	
	olambda2 <- t((olambda21 + olambda22) / olambda23) %*% ((olambda21 + olambda22) / olambda23)
	
	# part 3 Lambda observation
	
	olambda31 <- -t(bspv[,delta3==1]) * subx[delta3==1]
	
	olambda3 <- t(olambda31) %*% olambda31
	
	# total beta observation
	
	res <- (olambda1 + olambda2 + olambda3) / size
	
	res
	
}	
########################################################################################################	

	
################ for variance estimation based on ICsurv (function obtained from ICsurv package) ########
PH.Louis.ICsurv <- function (b, g, bLi, bRi, d1, d2, d3, Xp) 
{
    xb <- Xp %*% b
    A <- matrix(0, length(b), length(b))
    B <- matrix(0, length(g), length(b))
    C <- matrix(0, length(g), length(g))
    for (i in 1:length(b)) {
        for (j in 1:length(b)) {
            A[i, j] <- -sum(((d1 + d2) * bRi %*% g + d3 * bLi %*% 
                g) * exp(xb) * (Xp[, i] * Xp[, j]))
        }
    }
    for (i in 1:length(b)) {
        for (j in 1:length(g)) {
            B[j, i] <- -sum(((d1 + d2) * bRi[, j] + d3 * bLi[, 
                j]) * exp(xb) * Xp[, i])
        }
    }
    ci <- 1 - exp(-(bRi %*% g) * exp(xb))
    di <- 1 - exp(-(bRi %*% g - bLi %*% g) * exp(xb))
    di[d2 == 0] = 1
    hi <- d1 * (bRi %*% g) * exp(xb)
    ti <- d2 * (bRi %*% g - bLi %*% g) * exp(xb)
    D <- rbind(cbind(A, t(B)), cbind(B, C))
    E <- matrix(0, length(b), length(b))
    F <- matrix(0, length(g), length(b))
    G <- matrix(0, length(g), length(g))
    VZi <- (hi^2/ci) * (1 - 1/ci) + hi/ci
    VWi <- (ti^2/di) * (1 - 1/di) + ti/di
    for (i in 1:length(b)) {
        for (j in 1:length(b)) {
            E[i, j] <- sum((VZi + (d2 + d3) * VWi) * Xp[, i] * 
                Xp[, j])
        }
    }
    for (i in 1:length(b)) {
        for (j in 1:length(g)) {
            zpil <- bRi[, j] * g[j]/(bRi %*% g)
            num <- (bRi %*% g - bLi %*% g)
            num[d2 == 0] <- 1
            wpil <- (bRi[, j] - bLi[, j]) * g[j]/num
            CovZilZi <- zpil * (hi/ci) * (1 + hi - hi/ci)
            CovWilWi <- wpil * (ti/di) * (1 + ti - ti/di)
            F[j, i] <- 1/(g[j]) * sum((CovZilZi + (d2 + d3) * 
                CovWilWi) * Xp[, i])
        }
    }
    for (i in 1:length(g)) {
        for (j in 1:length(g)) {
            zpi <- bRi[, i] * g[i]/(bRi %*% g)
            num <- (bRi %*% g - bLi %*% g)
            num[d2 == 0] <- 1
            wpi <- (bRi[, i] - bLi[, i]) * g[i]/num
            zpj <- bRi[, j] * g[j]/(bRi %*% g)
            num <- (bRi %*% g - bLi %*% g)
            num[d2 == 0] <- 1
            wpj <- (bRi[, j] - bLi[, j]) * g[j]/num
            Covz <- zpi * zpj * (hi^2/ci) * (1 - 1/ci)
            Covw <- wpi * wpj * (ti^2/di) * (1 - 1/di)
            G[i, j] <- 1/(g[i] * g[j]) * sum(Covz + (d2 + d3) * 
                Covw)
        }
    }
    H <- rbind(cbind(E, t(F)), cbind(F, G))
    hess <- -(D + H)
    return(hess)
}
###########################################################################################
	
################### knot sequence of spline functions ####################################
get_knots <- function(ctu, ctv, delta1, delta2, delta3)
{
	knot_pre<-c(ctu[delta1==1], ctu[delta1==1], ctu[delta2==1], ctv[delta2==1], ctv[delta3==1], ctv[delta3==1])

	### qurtile knots
	qmin <- 0
	q1 <- quantile(knot_pre, (1 : (n.spline + 1 - spline.ord))/(n.spline + 2 - spline.ord))
	qmax <- max(knot_pre) + 0.001
	knot <- c(rep(qmin, spline.ord), q1, rep(qmax, spline.ord))

	knot
}
###########################################################################################

library(MASS)
library(pracma)
library(splines)
library(survival)
library(ICsurv)

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

 

 