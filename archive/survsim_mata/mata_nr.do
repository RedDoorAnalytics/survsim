clear all
local lambda1 0.01
local lambda2 0.01
local gamma1 0.1
local gamma2 1.1
local pmix 0.5
local N 50000
mata
lambda1 = strtoreal(st_local("lambda1"))
lambda2 = strtoreal(st_local("lambda2"))
gamma1 = strtoreal(st_local("gamma1"))
gamma2 = strtoreal(st_local("gamma2"))
pmix = strtoreal(st_local("pmix"))
N = strtoreal(st_local("N"))

U = runiform(N,1)
function ww(x,a) {
 fn = `pmix' *exp(-`lambda1'*x^(`gamma1'))   + (1-`pmix')*exp(-`lambda2'*x^(`gamma2')) - a
 df = -`pmix'*`lambda1'*`gamma1'*x^(`gamma1' - 1)*exp(-`lambda1'*x^(`gamma1')) - (1-`pmix')*`lambda2'*`gamma2'*x^(`gamma2' - 1)*exp(-`lambda2'*x^(`gamma2'))
 return(fn,df)
}

a = 0.3
mm_nrroot(x=1, &ww(), 0, 1000, a)
Y = J(N,1,.)
for(i=1;i<=N;i++){
 bob = mm_nrroot(x=1, &ww(), 0, 1000, a=U[i,1])
 Y[i,1] = x
}
Y
end
set obs `N'
gen t =  .
mata st_store(.,"t",Y)
gen died = t<=5
replace t = 5 if died==0
stset t, f(died=1)
stpm2 , df(5) scale(h)
predict s, s
twoway (line s _t, sort) ///
  (function `pmix' *exp(-`lambda1'*x^(`gamma1'))   + (1-`pmix')*exp(-`lambda2'*x^(`gamma2')), range(0 5))

  
  
  