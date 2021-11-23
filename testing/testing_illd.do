//local drive Z:\
local drive /Users/Michael/Documents

//merlin
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
clear all
cd "`drive'/merlin"
do ./build/buildmlib.do
mata mata clear

//predictms
adopath ++ "`drive'/multistate/multistate/predictms"
clear all
cd "`drive'/multistate/multistate"
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

timer clear 
timer on 1
survsim time state event, transmatrix(tmat)                                     ///
	hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t})) covariates(trt 0.1))      ///
	hazard2(dist(weibull) lambda(0.01) gamma(1.3) covariates(trt -0.5))         ///
	hazard3(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1) 			///
			tdefunction(log({t})))  											///
	maxtime(3)
timer off 1
timer list
list if inlist(_n,1,4,16,112)

cap drop time* state* event*

gen sstate = (runiform()>0.3) + 1
	
survsim time state event, transmatrix(tmat)                                  ///
	hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t})) covariates(trt 0.1))   ///
	hazard2(dist(weibull) lambda(0.01) gamma(1.3) covariates(trt -0.5))      ///
	hazard3(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1)) tdefunc(log({t}))  ///
	maxtime(3) startstate(sstate)
	
list if inlist(_n,1,4,16,112)
	
