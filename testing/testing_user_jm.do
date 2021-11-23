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
set obs 2000
set seed 98732
gen id = _n
gen trt = runiform()>0.5

gen b0 = rnormal()
gen b1 = rnormal()


survsim stime died, hazard(0.1:*1.2:*{t}:^(1.2:-1) :* exp(0.2 :* (b0 :+ b1:*{t})) ) 	///
					maxtime(10) covariates(trt -0.5)

expand 5
bys id : gen time = _n-1
drop if stime < time
bys id: gen y = b0 + b1 * time + rnormal(0,0.5)

bys id: replace stime = . if _n>1
bys id: replace died = . if _n>1

merlin 	(stime trt EV[y], family(weib, failure(died)) timevar(stime))	///
		(y time M1[id]@1 time#M2[id]@1, family(gaussian) timevar(time))
