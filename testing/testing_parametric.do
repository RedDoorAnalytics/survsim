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
set obs 1000000
set seed 98732
gen trt = runiform()>0.5

gen lt0 = runiform()*2

survsim stime, dist(exp) lambdas(0.1) covariates(trt -0.5)

cap drop stime

survsim stime died, dist(gomp) lambdas(0.1) gammas(1.2) maxt(10) covariates(trt -0.5) ltruncated(lt0)

stset stime, f(died) enter(lt0)
streg trt, dist(gomp) nohr

// stpm2 trt, df(1) scale(h) noorthog 
merlin (_t  trt, family(gomp, failure(_d) ltruncated(_t0)))
