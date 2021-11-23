//local drive Z:\
local drive /Users/Michael/Documents

//merlin and galahad
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/merlin/galahad/galahad"
clear all
do ./build/buildmlib.do
mata mata clear

//survsim
adopath ++ "`drive'/survsim/survsim"

clear
matrix tmat = (.,1,2\.,.,3\.,.,.)
mat list tmat
set obs 1000
set seed 9865
gen trt = runiform()>0.5
gen age = rnormal(50,3)

//2->3 depends on attained age 
cap drop time* state* event*
survsim time state event, transmatrix(tmat)                                     	///
		hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t})) covariates(trt 0.1))    ///
		hazard2(dist(weibull) lambda(0.01) gamma(1.3) covariates(trt -0.5))         ///
		hazard3(user(0.1 :* ({t}:+age) :^ 1.5) covariates(trt -0.5) tde(trt 0.1))  	///
		maxtime(3)

//2->3 depends on attained age and time since start	
cap drop time* state* event*
survsim time state event, transmatrix(tmat)                                  				///
		hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*({t}:+age))) covariates(trt 0.1))   	///
		hazard2(dist(weibull) lambda(0.01) gamma(1.3) covariates(trt -0.5))      			///
		hazard3(user(0.1 :* ({t}:+age) :^ 1.5 :* exp(0.1 :* log({t}))))  					///
		maxtime(3) 

//2->3 depends on attained age and main timescale, and time at entry
cap drop time* state* event*
survsim time state event, transmatrix(tmat)                                  				///
		hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*({t}:+age))) covariates(trt 0.1))   	///
		hazard2(dist(weibull) lambda(0.01) gamma(1.3) covariates(trt -0.5))      			///
		hazard3(user(0.1 :* ({t}:+age) :^ 1.5 :* exp(0.1 :* log({t}) :- 0.1 :* {t01})))  	///
		maxtime(3) 
