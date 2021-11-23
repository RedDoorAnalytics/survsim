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
gen age = rnormal()
gen lt0 = runiform()*2

survsim stime died, 	hazard(0.1:*1.2:*{t}:^(0.2))	 		///
						covariates(trt -0.5 age 0.1) 			///
						maxt(10) marginal


stset stime, f(died)
streg trt age, dist(weib) nohr
