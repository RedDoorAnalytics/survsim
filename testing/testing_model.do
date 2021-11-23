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

webuse brcancer,clear
stset rectime, f(censrec=1) scale(365)

// stpm2 hormon, scale(h) df(3) //tvc(hormon) dftvc(2)
//streg hormon, dist(weib) nohr
merlin (_t hormon, family(pwe, knots(1 3) failure(_d)))
predict s1, surv
est store m1
expand 100

survsim stime died, model(m1) maxt(7)

merlin (stime hormon, family(pwe, knots(1 3) failure(died)))
predict s2, surv

twoway (scatter s1 _t)(scatter s2 stime)
