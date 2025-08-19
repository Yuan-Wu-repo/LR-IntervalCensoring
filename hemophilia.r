library(MASS)
library(pracma)
library(splines)
library(survival)
library(ICsurv)

source("functions.R")

data(Hemophilia)

delta1 <- Hemophilia$d1
delta2 <- Hemophilia$d2
delta3 <- Hemophilia$d3
size <- dim(Hemophilia)[1]
ctu <- rep(0, size)
ctv <- rep(0, size)
ctu[delta1 == 1] <- Hemophilia$R[delta1 == 1]
ctv[delta1 == 1] <- ctu[delta1 == 1] + 0.001
ctu[delta2 == 1] <- Hemophilia$L[delta2 == 1]
ctv[delta2 == 1] <- Hemophilia$R[delta2 == 1]
ctv[delta3 == 1] <- Hemophilia$L[delta3 == 1]
ctu[delta3 == 1] <- ctv[delta3 == 1] - 0.001
z <- rbind(Hemophilia$Low, Hemophilia$Medium, Hemophilia$High)

newt <- rep(0, size)
newt[delta1 == 1] <- ctu[delta1 == 1] / 2
newt[delta2 == 1] <- (ctu[delta2 == 1] + ctv[delta2 == 1]) / 2
newt[delta3 == 1] <- ctv[delta3 == 1]
newdelta <- rep(0, size)
newdelta[delta1 == 1 | delta2 == 1] <- 1

d1 <- Hemophilia$d1
d2 <- Hemophilia$d2
d3 <- Hemophilia$d3
Li <- Hemophilia$L
Ri <- Hemophilia$R
Xp <- cbind(Hemophilia$Low, Hemophilia$Medium, Hemophilia$High)

n.variable <- 3
n.spline <- 11
spline.ord <- 4

###################################
##### naive convering method ######
###################################
fitnaive <- coxph(Surv(newt, newdelta) ~ z[1,] + z[2,] + z[3,], method="breslow", init = c(0, 0, 0))
### to avoid diag(fitnaive$var) has negative components, when it is negative let it be 0
diagnaive <- apply(rbind(diag(fitnaive$var), rep(0, n.variable)), 2, max)

fitcoef <- c(fitnaive$coef)
fitvar <- c(fitnaive$var)

#### note that loglik[2] is the loglikelihood for final coefficients
loglikefull.naive <- fitnaive$loglik[2]

### Wald test p-value for each variable among Low, Medium, High being 0
2*pnorm(abs(fitcoef[1]/sqrt(diagnaive)[1]), mean=0, sd=1, lower.tail = FALSE)
2*pnorm(abs(fitcoef[2]/sqrt(diagnaive)[2]), mean=0, sd=1, lower.tail = FALSE)
2*pnorm(abs(fitcoef[3]/sqrt(diagnaive)[3]), mean=0, sd=1, lower.tail = FALSE)
	
### Wald test p-value for all variables being 0
pchisq(q = t(fitcoef) %*% ginv(fitnaive$var) %*% fitcoef, df = 3, lower.tail = FALSE)

### partial coeffcients cox regression	
fitnaive <- coxph(Surv(newt, newdelta) ~  z[2,] + z[3,], method="breslow")
loglikepart1.naive <- fitnaive$loglik[2]
fitnaive <- coxph(Surv(newt, newdelta) ~  z[1,] + z[3,], method="breslow")
loglikepart2.naive <- fitnaive$loglik[2]
fitnaive <- coxph(Surv(newt, newdelta) ~  z[1,] + z[2,], method="breslow")
loglikepart3.naive <- fitnaive$loglik[2]

## likelihood ratio test p-value for each variable among Low, Medium, High being 0
pchisq(q = 2 * (loglikefull.naive - loglikepart1.naive), df = 1, lower.tail = FALSE) 
pchisq(q = 2 * (loglikefull.naive - loglikepart2.naive), df = 1, lower.tail = FALSE)
pchisq(q = 2 * (loglikefull.naive - loglikepart3.naive), df = 1, lower.tail = FALSE)

### Likelihood ratio test p-value for all variable being 0
summary(fitnaive)$logtest[3]
############## end naive converting method #################

#######################################
#### interval censoring method ########
#######################################

############## unconstrained R function for finding optimization point
	
knotb <- get_knots(ctu, ctv, delta1, delta2, delta3)
### bspline results		 
bspu <- t(splineDesign(knots = knotb, x = ctu, ord = spline.ord))
### from bspline to ispline  
pre.ispu <- apply(bspu, 2, rev)
ispu <- apply(pre.ispu, 2, cumsum)
ispu <- apply(ispu, 2 ,rev)
ispu <- ispu[-1,]

### bspline results
bspv <- t(splineDesign(knots = knotb, x = ctv, ord = spline.ord))
### from bspline to ispline  
pre.ispv <- apply(bspv, 2, rev)
ispv <- apply(pre.ispv, 2, cumsum)
ispv <- apply(ispv, 2, rev)
ispv <- ispv[-1,]

## create vecters of 1s with length of the numbers of left, interval and right censored
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
out <- optim(par = rep(0, (n.variable + n.spline)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))
betaest <- out$par[(n.spline+1) : (n.spline+n.variable)]
		
#### variance based on theory from Huang et al. (2008) or Zhang et al. (2010)'s least square approach 
hat_O <- observ.beta(out$par) - observ.beta.lambda(out$par) %*% ginv(observ.lambda(out$par)) %*% t(observ.beta.lambda(out$par))

### to avoid diag(ginv(hat_O)) has negative components, when it is negative let it be 0
diag.inv <- apply(rbind(diag(ginv(hat_O)), rep(0, n.variable)), 2, max)
se<-sqrt(diag.inv / size)	


### Wald test p-value for each variable among Low, Medium, High being 0
2*pnorm(abs(betaest[1]/se[1]), mean=0, sd=1, lower.tail = FALSE)
2*pnorm(abs(betaest[2]/se[2]), mean=0, sd=1, lower.tail = FALSE)
2*pnorm(abs(betaest[3]/se[3]), mean=0, sd=1, lower.tail = FALSE)
	
### Wald test p-value for all variables being 0
pchisq(q = t(betaest) %*% hat_O %*% betaest * size, df = 3, lower.tail = FALSE)
	
#### beta estimation and wald tests using ICsurv package 
fitsemi <- PH.ICsurv.EM(d1, d2, d3, Li, Ri, Xp, n.int = n.spline-3, order = spline.ord-1, g0 = rep(1,n.spline), b0 = rep(0,n.variable), t.seq = seq (0,57,1), tol = 0.001)

#### to avoid diag(fitnaive$var) has negative components, when it is negative let it be 0
diagvar.b <- apply(rbind(diag(fitsemi$var.b), rep(0, n.variable)), 2, max)

### Wald test p-value for each variable among Low, Medium, High being 0
2*pnorm(abs(fitsemi$b[1]/sqrt(diagvar.b[1])), mean=0, sd=1, lower.tail = FALSE)
2*pnorm(abs(fitsemi$b[2]/sqrt(diagvar.b[2])), mean=0, sd=1, lower.tail = FALSE)
2*pnorm(abs(fitsemi$b[3]/sqrt(diagvar.b[3])), mean=0, sd=1, lower.tail = FALSE)

### Wald test p-value for all variables being 0
pchisq(q = t(fitsemi$b) %*% ginv(fitsemi$var.b) %*% fitsemi$b, df = 3, lower.tail = FALSE)

##### likelihood ratio approach for sieve estimation 
## loglikelihood for full model
loglikefull <- loglike(out$par)
 
index <- 1
if (index != 0) {
 z <- oldz[-index,]
 fixsubx <- matrix(c(0*oldz[index,]), 1, size)
 n.total <- n.spline + n.variable -1			 
}	 
out <- optim(par = rep(0, (n.variable + n.spline -1)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))
## loglikelihood for model without variable Low
loglikepart1 <- loglike(out$par)

index <- 2
if (index != 0) {
 z <- oldz[-index,]
 fixsubx <- matrix(c(0*oldz[index,]), 1, size)
 n.total <- n.spline + n.variable -1			 
}	 
out <- optim(par = rep(0, (n.variable + n.spline -1)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))
## loglikelihood for model without variable Medium
loglikepart2 <- loglike(out$par)

index <- 3  
if (index != 0) {
 z <- oldz[-index,]
 fixsubx <- matrix(c(0*oldz[index,]), 1, size)
 n.total <- n.spline + n.variable -1			 
}	 
out <- optim(par = rep(0, (n.variable + n.spline -1)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))
## loglikelihood for model without variable High
loglikepart3 <- loglike(out$par)

## likelihood ratio test for each variable among Low, Medium, High being 0
pchisq(q = -2 * (loglikefull - loglikepart1), df = 1, lower.tail = FALSE)
pchisq(q = -2 * (loglikefull - loglikepart2), df = 1, lower.tail = FALSE) 
pchisq(q = -2 * (loglikefull - loglikepart3), df = 1, lower.tail = FALSE) 

fixsubx <- matrix(0, 1, 3) %*% oldz
n.total <- n.spline		 
out <- optim(par = rep(0, n.spline), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))
## loglikelihood for model with no variables left
loglikenull <- loglike(out$par)

## likelihood ratio test for all variables being 0
pchisq(q = -2 * (loglikefull - loglikenull), df = 3, lower.tail = FALSE) 
#### end of interval censoring approach ####
			