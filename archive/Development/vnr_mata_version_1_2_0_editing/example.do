pr drop _all
adopath ++ "E:\Survival Analysis\survsim\SSC\version_1_2_0_editing"

clear all
set obs 100000
gen trt = runiform()>0.5
cap drop t_new
set seed 87578769
timer clear
timer on 1
tr:survsim t_new, basehaz(0.1*1.2*#t^(1.2-1)) cov(trt -0.2) tde(trt 0.1) tdefunc(log(#t)) centol(0.00000001)
timer off 1
gen died = t_new<5
replace t_new = 5 if died==0
stset t_new, f(died)
streg trt, dist(w) nohr
est store s1

set seed 87578769
cap drop t_old
timer on 2
pr drop _all
survsim t_old, lambdas(0.1) gammas(1.2) dist(weibull) cov(trt -0.2) tde(trt 0.1)
timer off 2
timer list
gen died2 = t_old<5
replace t_old = 5 if died2==0
stset t_old, f(died2)
streg trt, dist(w) nohr
est store s2
est tab s*


/* Time dependent effects */

clear
pr drop _all
set obs 10000
gen trt = runiform()>0.5
gen bin = runiform()>0.5
cap drop t_new
//survsim t_new, basehaz(0.1*1.5*#t^(1.5-1)) cov(trt -0.5) tde(trt 0.1) tdefunction(log(#t))
survsim t_new6, basehaz(exp(log(#t))) cov(trt -0.5)

//survsim t_new, basehaz(0.1*1.5*@x^(1.5-1)) cov(trt -0.5 bin 0.5) tde(trt 0.1 bin -0.1) logtime
//survsim t_new, lambdas(0.1) gammas(1.5) dist(weibull) cov(trt -0.5) tde(trt 0.15)
gen died = t_new<5
replace t_new = 5 if died==0
stset t_new, f(died)
streg trt, dist(w) nohr
stcox trt
estat phtest,detail
stpm2 trt bin, df(3) scale(h) tvc(trt) dftvc(1)
predict hr2, hrnumerator(bin 1) ci
twoway (line hr2* _t, sort)

. clear
. set obs 1000
. gen trt = rbinomial(1,0.5)
. tr:survsim stime, lambdas(0.1) gammas(1.5) cov(trt -0.5) dist(w)
. gen died = stime <= 5
. replace stime = 5 if died == 0
. stset stime, f(died = 1)
. streg trt, dist(weibull) nohr

. clear
. set obs 1000
. gen trt = rbinomial(1,0.5)
. survsim stime, n(1000) lambdas(0.1) gammas(0.05) dist(gompertz)

. clear
. set obs 1000
. gen trt = rbinomial(1,0.5)
. survsim stime, n(1000) mixture lambdas(0.1 0.05) gammas(1 1.5) pmix(0.5) dist(w)

. clear
. set obs 1000
. gen trt = rbinomial(1,0.5)
. survsim stime status, n(1000) cr ncr(4) lambdas(0.1 0.05 0.1 0.05) gammas(0.5 1.5 1 1.2) cov(trt 0.2 0.1 -0.1 0.4) dist(w)

. clear
. set obs 1000
. gen trt = rbinomial(1,0.5)
. survsim stime, n(1000) lambdas(0.1) gammas(1.5) cov(trt -0.5) tde(trt 0.05) dist(w)
