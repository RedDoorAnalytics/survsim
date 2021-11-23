adopath ++ "H:\Survival Analysis\survsim\survsim_version_1"

/* Standard Weibull */

	clear 
	pr drop _all
	set obs 1000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,10)
	gen sex = rbinomial(1,0.2)
	survsim stime, n(1000) lambdas(0.1) gammas(1.5) covariates(trt 1 age 0.01 sex -1)
	gen died=1
	stset stime, f(died)
	streg trt age sex, nohr dist(w)

/* Mixture Weibull */

	clear 
	pr drop _all
	set obs 1000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,10)
	survsim stime, n(1000) mixture lambdas(0.1 0.1) gammas(1 1.5) cov(trt 0.5) 
	
/* Competing risks */

	clear
	pr drop _all
	set obs 100000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,2)
	survsim stime status, cr ncr(3) lambdas(0.1 0.1 0.1) ///
	gammas(0.6 0.5 0.7) covariates(trt 0.5 0.1 0.2)
	

