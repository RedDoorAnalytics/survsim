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

gen lt0 = runiform()*2

survsim time state event , 	hazard1(user(0.1 :* 1.5 :* {t} :^ (0.5)))			///
							hazard2(dist(weib) lambda(0.1) gamma(1.2))			///
							maxtime(10)			

