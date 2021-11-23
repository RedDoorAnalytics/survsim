adopath ++ "C:\Users\Michael\Documents\Survival Analysis\Simulation\survsim\SSC\version_2_0_4_editing"
pr drop _all

mata mata clear

clear all

do ""
//twoway function y = exp(-1 +0.02*x - 0.03*x^2 + 0.005*x^3), range(0 5)

clear all
set obs 1000

gen trt = runiform()>0.5
gen enter = rnormal(5,0.01)
set seed 28762
survsim test1 event1, loghazard(trt :+ log(0.1) :+ log(1.2) :+ (1.2-1):*log(#t)) cov(trt -0.5) maxt(10) nodes(50)
stset test1, f(event1) //enter(enter)
streg trt, dist(w) nohr


// tde interacted with time
clear all
set obs 10000
set seed 28762
gen trt = runiform()>0.5
survsim test1, loghazard(log(0.1) :+ log(1.2) :+ (1.2-1):*log(#t)) maxt(10) cov(trt -0.5) tde(trt 0.02)
gen event1 = test1<10
stset test1, f(event1)
streg trt, dist(w) nohr
stpm2 trt, scale(h) df(3) tvc(trt) dftvc(2)
predict hr, hrnumerator(trt 1) ci
line hr* _t,sort name(h1,replace)


// tde interacted with time^2
clear all
set obs 10000
set seed 28762
gen trt = runiform()>0.5
survsim test1, loghazard(log(0.1) :+ log(1.2) :+ (1.2-1):*log(#t)) maxt(10) cov(trt -0.5) tde(trt 0.02) tdefunc(#t:^2)
gen event1 = test1<10
stset test1, f(event1)
streg trt, dist(w) nohr
stpm2 trt, scale(h) df(3) tvc(trt) dftvc(2)
predict hr, hrnumerator(trt 1) ci
twoway line hr* _t,sort name(h2,replace)

//FP baseline
clear all
set obs 10000
set seed 28762
gen trt = runiform()>0.5
survsim test1, loghazard(#t :- 0.2:*#t:^2 :+ 0.7:*#t:^3) maxt(10) cov(trt -0.5)
gen event1 = test1<10
stset test1, f(event1)
stpm2 trt, scale(h) df(3) tvc(trt) dftvc(2)
predict h1, hazard
scatter h1 _t




/* Standard Weibull */

	clear 
	pr drop _all
	set obs 10000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,10)
	gen sex = rbinomial(1,0.2)
	survsim stime, lambdas(0.1) cov(trt -0.5) d(exp) tde(trt 0.5)
	//survsim stime, lambdas(0.1) gammas(1.5) covariates(trt 1 age 0.01 sex -1) tde(trt 0.2)
	gen died=stime<5
	replace stime = 5 if died==0
	stset stime, f(died)
	//streg trt, dist(e) nohr 
	stpm2 trt, scale(hazard) df(1) tvc(trt) dftvc(1)
	predict hr, hrnumerator(trt 1) ci
	line hr* _t, sort

/* Standard Gompertz */

	clear 
	pr drop _all
	set obs 10000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,10)
	gen sex = rbinomial(1,0.2)
	survsim stime, lambdas(0.1) dist(exp) covariates(trt 1 age 0.01 sex -1) tde(trt 0.4)
	gen died=stime<5
	replace stime = 5 if died==0
	stset stime, f(died)
	streg trt age sex, dist(gom) nohr
	stpm2 trt age sex, scale(hazard) df(3) tvc(trt) dftvc(1)
	predict hr, hrnumerator(trt 1) ci
	line hr* _t, sort
	
/* Mixture */
	local pmix=0.5
	local l1 = 0.1
	local l2 = 0.1
	local g1 = 1.5
	local g2 = 0.5

	clear 
	pr drop _all
	set obs 1000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,10)
	survsim stime, mixture lambdas(0.1 0.1) gammas(1 1.1) d(w) pmix(0.5) maxt(5)
	gen died=stime<5
	replace stime = 5 if died==0
	stset stime, f(died)
	//streg trt age sex, dist(gom) nohr
	stpm2 trt, scale(hazard) df(3)
	
	
	clear 
	pr drop _all
	set obs 1000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,10)
	survsim stime,  mixture lambdas(0.1 0.1) gammas(0.1 2) cov(trt 0.5) dist(gompertz) maxt(5)
	
	local pmix=0.8
	local l1 = 10
	local l2 = 1
	local g1 = 2
	local g2 = 15
	//twoway function y=`pmix'*exp((`l1'/`g1')*(1-exp(`g1'*x))) +  (1-`pmix')*exp((`l2'/`g2')*(1-exp(`g2'*x))) , range(0 5)
	
/* Competing risks */

	clear
	set seed 3476147
	pr drop _all
	set obs 10000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,2)
	survsim stime status, cr ncr(3) lambdas(0.1 0.1 0.1) ///
	covariates(trt 0.5 0.1 0.2) dist(exp) tde(trt 0.02 0.05 0.05)

	clear
	set seed 3476147
	pr drop _all
	set obs 10000
	gen trt = rbinomial(1,0.5)
	gen age = rnormal(45,2)
	survsim stime status, cr ncr(3) lambdas(0.1 0.1 0.1) ///
	gammas(0.6 0.5 0.7) covariates(trt 0.5 0.1 0.2) dist(w)
	
	
//======================================================================================================================================================//
//SiM paper examples


//Example 1: standard parametric distribution with a binary and continuous covariate

clear
//sjlog using example1, replace
//Simulate 1000 survival times
set obs 1000
//Generate a binary treatment group indicator and a continuous age covariate
gen treatment = rbinomial(1,0.5)
gen age = rnormal(65,12)
//Simulate times from an exponential distribution
survsim time1, dist(exp) lambda(0.1) covariates(treatment -0.5 age 0.01)
//Simulate times from a Gompertz distribution
survsim time2, dist(gompertz) lambda(0.1) gamma(1.2) covariates(treatment -0.5 age 0.01)
//Simulate times from a Weibull distribution
survsim time3, dist(weibull) lambda(0.1) gamma(1.2) covariates(treatment -0.5 age 0.01)
//sjlog close, replace


//Example 2: 2-component mixture Weibull with proportional hazards

clear
//sjlog using example2, replace
//Simulate 1000 survival times
set obs 1000
//Generate a binary treatment group indicator and a continuous age covariate
gen treatment = rbinomial(1,0.5)
gen age = rnormal(65,12)
//Simulate times from a 2-component mixture Weibull with proportional hazards
survsim time1, mixture lambdas(0.1 0.05) gammas(1 1.5) pmix(0.5) covariates(treatment -0.5 age 0.01) maxtime(5)
//sjlog close, replace

//Example 3: FP PH

clear
//sjlog using example3, replace
//Simulate 1000 survival times
set obs 1000
//Generate a binary treatment group indicator and a continuous age covariate
gen trt = rbinomial(1,0.5)
gen age = rnormal(65,12)
//Simulate times from a user-defined log hazard function
survsim stime event, loghazard(-18 :+ 7.3:*#t:-11.5:*#t:^(0.5):*log(#t) :+ 9.5:*#t:^0.5) maxt(5) nodes(30) covariates(trt -0.5 age 0.02)
//sjlog close, replace

//Example 4: FP NPH

clear
//sjlog using example4, replace
//Simulate 1000 survival times
set obs 1000
//Generate a binary treatment group indicator and a continuous age covariate
gen trt = rbinomial(1,0.5)
gen age = rnormal(65,12)
//Simulate times from a user-defined log hazard function
survsim stime event, loghazard(-18 :+ 7.3:*#t:-11.5:*#t:^(0.5):*log(#t) :+ 9.5:*#t:^0.5) maxt(5) nodes(30) cov(trt -0.7 age 0.02) tde(trt 1) tdefunc(0.01:*#t :+ 0.4:*log(#t))
//sjlog close, replace


clear 
//sjlog using example5, replace
//Simulate 1000 survival times
set obs 1000
//Generate a binary treatment group indicator, a continuous age covariate and a bad prognosis indicator variable
gen trt = runiform()>0.5
gen badprog = runiform()<0.4
gen age = rnormal(65,12)
//Generate an indicator variable for patients who swap treatment group, dependent on bad prognosis
gen sp = cond(badprog==1,runiform()<0.4,runiform()<0.2)
//Generate the time of swapping treatment group
gen tswap = cond(sp==1,4.5*runiform() + 0.5,5)
//Define the treatment effect
local trteffect = -0.5
survsim test1 event1, loghazard(-2.3:+2:*#t:-#t:^(2):+0.12:*#t:^3 :`trteffect':*((#t:<=tswap):*trt :+ (#t:>tswap):*(1:-trt))) maxt(5) nodes(30) cov(badprog `=log(1.5)' age 0.02)
//sjlog close, replace


clear all
//sjlog using example6, replace
//Simulate 1000 survival times
set obs 1000
//Generate a binary treatment group indicator and a continuous age covariate
gen trt = runiform()>0.5
gen age = rnormal(65,12)
//Define the association between the biomarker and survival
local alpha = 0.25
//Generate the random intercept and random slopes for the longitudinal submodel
gen b0 = rnormal(0,1)
gen b1 = rnormal(1,0.5)
survsim stime event, loghazard(-2.3:+2:*#t:-#t:^(2):+0.12:*#t:^3 :+ `alpha' :* (b0 :+ b1 :* #t)) maxt(5) nodes(30) mint(0.03) cov(trt -0.5 age 0.02)
//Generate observed biomarker values at times 0, 1, 2, 3 , 4 years
gen id = _n					
expand 5
bys id: gen meastime = _n-1
//Remove observations after event or censoring time
bys id: drop if meastime>=stime
//Generate observed biomarker values incorporating measurement error
gen response = b0 + b1*meastime + rnormal(0,0.5)
//sjlog close, replace



//Example ??: Delayed entry

clear
//Simulate 1000 survival times
set obs 1000
//Generate a binary treatment group indicator and age at entry
gen treatment = rbinomial(1,0.5)
gen age_entry = rnormal(30,3)
//Simulate times from a user defined log-hazard function, with delayed entry
survsim time1 event, loghazard(0.01:*#t:^(-2) :- 8:*#t:^(-0.5)) covariates(treatment -0.5) enter(age_entry) maxtime(50)

