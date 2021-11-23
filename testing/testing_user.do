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

clear
set obs 100000
set seed 98732
gen trt = runiform()>0.5

gen lt0 = runiform()*2
gen b0 = 1 
survsim stime died, hazard(b0:*0.1:*1.2:*{t}:^(1.2:-1)) maxtime(10) covariates(trt -0.5) ltruncated(lt0)
merlin (stime trt, family(weib, failure(died) ltruncated(lt0)))

cap drop stime
cap drop died
survsim stime died, chazard(exp(b0):*0.1:*{t}:^(1.2)) maxtime(10) covariates(trt -0.5) ltruncated(lt0)
merlin (stime trt, family(weib, failure(died) ltruncated(lt0)))

cap drop stime
cap drop died
survsim stime died, mixture lambdas(0.1 0.1) gammas(1.2 0.8) dist(weib) maxt(10) cov(trt -0.5) ltruncated(lt0)
merlin (stime trt, family(rp, failure(died) df(2) ltruncated(lt0)))
