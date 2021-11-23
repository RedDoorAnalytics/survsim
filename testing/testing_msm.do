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
set obs 1000
set seed 98732
gen trt = runiform()>0.5

mat tmat = (.,1,2\.,.,.\.,.,.)

gen t0 = runiform() * 2

set seed 9898
survsim stime state event , 		transmatrix(tmat) 									///
									hazard1(dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5))			///
									hazard2(dist(weib) lambda(0.1) gamma(0.8) cov(trt 0.1))			///
									maxtime(10) ltruncated(t0)

cap drop stime*
cap drop state*
cap drop event*
set seed 9898

survsim stime state event , 	transmatrix(tmat) 																		///
								hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t})) cov(trt -0.5))						///
								hazard2(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1) tdefunction(log({t})))	///
								maxtime(10) ltruncated(t0)

su stime* state*
cap drop stime*
cap drop state*
cap drop event*
set seed 9898
					
mata:
real matrix testf(t)
{
resq = exp(-2 :+ 0.2:* log(t) :+ 0.1:*t)
return(resq)
}
end


survsim stime state event , 	transmatrix(tmat) 																		///
								hazard1(user(testf({t})) cov(trt -0.5))						///
								hazard2(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1) tdefunction(log({t})))	///
								maxtime(10) ltruncated(t0)
su stime* state*
