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
set obs 1000000
set seed 9865
gen trt = runiform()>0.5

survsim time state event, transmatrix(tmat)                                     ///
	hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t})) covariates(trt 0.1))      ///
	hazard2(dist(weibull) lambda(0.01) gamma(1.6) covariates(trt -0.5))         ///
	hazard3(dist(weibull) lambda(0.01) gamma(1.3) covariates(trt -0.5) reset)  	///
	maxtime(3)
	
	
stset time2 if state1==2, f(event2) enter(time1)
streg trt, dist(weib) nohr

gen rt =  time2 - time1
stset rt if state1==2, f(event2)
streg trt, dist(weib) nohr

gen test3 = event1==1 & state1==3
stset time1 , f(test3)

streg trt, dist(weib) nohr
