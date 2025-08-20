This project comprehensively studies the analysis of the interval-censored data under the <b>Cox proportional hazards model</b>. 


The following tasks are facilitated by the functions included in <b>functions.r</b>.

- Estimating the <b>sieve spline-based MLE</b> of the Cox regression parameters and the hazard function

- Estimating the variance of the sieve MLE of the Cox regression parameters using two existing approaches

- Generating the knot sequence for <b>B-splines (or I-splines)</b>


The script <b>simulaiton.r</b> aims to compare <b>our proposed likelihood ratio (LR) testing approach</b> and the existing <b>Wald testing approaches</b>. 

The script <b>hemophilia.r</b> applies different testing approaches analyzing this public available data set.

The script <b>UMvsICsurv.r</b> compares the computing time and estimation bias between <b>the proposed unconstrained maximization (UM)</b> and <b>the EM algorithm adopted in ICsurv package</b>.