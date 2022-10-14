
local drive /Users/Michael/Documents/reddooranalytics/products/survsim
cd `drive'

//survsim
adopath ++ "`drive'/survsim"
pr drop _all

clear
set obs 1000
set seed 9865
gen trt = runiform()>0.5
matrix tmat = (.,1,2\3,.,4\5,6,.)

survsim time state event, transmatrix(tmat)                                        ///
                hazard1(user(exp(-2 :+ 0.2:* log({t}) :+ 0.1:*{t})) covariates(trt 0.1))   ///
                hazard2(dist(weibull) lambda(0.01) gamma(1.3) covariates(trt -0.5))        ///
                hazard3(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1) tdefunction(log({t}))) ///
                hazard4(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1) tdefunction(log({t}))) ///
                hazard5(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1) tdefunction(log({t}))) ///
                hazard6(user(0.1 :* {t} :^ 1.5) covariates(trt -0.5) tde(trt 0.1) tdefunction(log({t}))) ///
                maxtime(3)
