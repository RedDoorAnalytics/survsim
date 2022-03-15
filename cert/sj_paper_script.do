
//local drive Z:\
local drive /Users/Michael/Documents/reddooranalytics/products

//merlin and galahad
// cd "`drive'/merlin"
// adopath ++ "`drive'/merlin"
// adopath ++ "`drive'/merlin/merlin"
// adopath ++ "`drive'/merlin/galahad/galahad"
// clear all
// do ./build/buildmlib.do
// mata mata clear

//survsim
adopath ++ "`drive'/survsim/survsim"


//Simulating survival times from standard parametric distributions

clear
set obs 300
set seed 134987
gen trt = runiform()>0.5
gen age = rnormal(50,3)

survsim stime, distribution(weibull) lambda(0.1) gamma(1.2)  ///
               covariates(trt -0.5 age 0.01)

gen censtime = runiform() * 5

survsim stime2 died2, distribution(weibull) lambda(0.1) gamma(1.2)     ///
                      covariates(trt -0.5 age 0.01) maxtime(censtime)

//Simulating survival times from a user-defined (log) (cumulative) hazard function

clear
set obs 500
set seed 134987
gen trt = runiform()>0.5

survsim stime1 died1, loghazard(-1:+0.02:*{t}:-0.03:*{t}:^2:+0.005:*{t}:^3)  ///
                      covariates(trt -0.5) maxtime(1.5)

survsim stime2 died2, loghazard(-1:+0.02:*{t}:-0.03:*{t}:^2:+0.005:*{t}:^3)    ///
                      covariates(trt -0.5) tde(trt 0.03) tdefunction(log({t})) ///
                      maxtime(1.5)

//Simulating survival times from a fitted \texttt{merlin} survival model

webuse brcancer, clear
stset rectime, f(censrec=1) scale(365)
merlin (_t hormon , family(weibull, failure(_d)))

estimates store m1

survsim stime5 died5, model(m1) maxtime(7)

stset stime5, failure(died5)
merlin (_t hormon , family(weibull, failure(_d)))

//Simulating competing risks data from specified cause-specific hazard functions

clear
set seed 398
set obs 1000
gen trt = runiform()>0.5
survsim time state event , hazard1(dist(weibull) lambda(0.1) gamma(0.8))  ///
                           hazard2(dist(exponential) lambda(0.02)         /// 
                           covariates(trt -0.5)) maxtime(10)

list if _n<=5

tabulate state1 event1

cap drop time* state* event*
set seed 32984575
survsim time state event,                                    ///
        hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t}))  ///
                covariates(trt 0.1))                         ///
        hazard2(dist(weibull) lambda(0.01) gamma(1.3)        ///
                covariates(trt -0.5))                        ///
        hazard3(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) ///
                tde(trt 0.1) tdefunction(log({t})))          ///
        maxtime(3)
tabulate state1 event1

//Simulating from an illness-death model

matrix tmat = (.,1,2\.,.,3\.,.,.)

mat colnames tmat = "healthy" "ill" "dead"
mat rownames tmat = "healthy" "ill" "dead"
mat list tmat

clear
set obs 1000
set seed 9865
gen trt = runiform()>0.5

survsim time state event, transmatrix(tmat)                  ///
        hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t}))  ///
                covariates(trt 0.1))                         ///
        hazard2(dist(weibull) lambda(0.01) gamma(1.3)        ///
                covariates(trt -0.5))                        ///
        hazard3(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) ///
                tde(trt 0.1) tdefunction(log({t})))          ///
        maxtime(3)

list if inlist(_n,1,4,16,112), compress

//Simulating from an illness-death model with multiple timescales

capture drop time* state* event*
set seed 98798
survsim time state event, transmatrix(tmat)                           ///
        hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t}))           ///
                covariates(trt 0.1))                                  ///
        hazard2(dist(weibull) lambda(0.01) gamma(1.3)                 ///
                covariates(trt -0.5))                                 ///
        hazard3(user(0.1 :* {t} :^ 1.5 :* exp(-0.05 :* ({t}:-{t0})))  ///
                covariates(trt -0.5)                                  ///
                tde(trt 0.1) tdefunction(log({t})))                   ///
        maxtime(3)
