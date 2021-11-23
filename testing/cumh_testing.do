adopath ++ "C:\Users\Michael\Documents\Survival Analysis\Simulation\survsim\SSC\version_2_0_3_editing"
pr drop _all

clear
set obs 1000
gen s0 = runiform()>0.5

//twoway function y = exp(-0.2+0.5*log(x)-0.04*log(x)^(2)+0.4*x^(0.5)), range(0 7)


survsim stime died, logcumh(-0.2:+0.5:*log(#t):-0.04:*log(#t):^(2):+0.4:*#t:^(0.5)) maxt(7) nodes(30) ///
         covariates(s0 -2.4) tde(s0 -0.22) tdefunction(log(#t))
stset stime, f(died)
stpm2 s0, scale(h) df(5) tvc(s0) dftvc(1)

gen trt = runiform()>0.5
survsim stime8 died3, logcumhazard(-1:+0.02:*#t:-0.03:*#t:^2:+0.005:*#t:^3) cov(trt -0.5) tde(trt 0.03) maxt(1.5)

twoway function y = (0.02 -0.06*x + 0.015*x^2)* exp(-1+0.02*x-0.03*x^2+0.005*x^3), range(0 7)
