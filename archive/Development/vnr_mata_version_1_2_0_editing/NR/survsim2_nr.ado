program define survsim2_nr
	syntax newvarname(min=1 max=2), 										///
										Lambdas(numlist min=1)				///	-Scale parameters-
																			///
									[										/// -Options-
																			///
										Gammas(numlist min=1)				///	-Shape parameters-
										Distribution(string)				/// -Parametric distribution-
																			///
										N(string)							///	-Number of survival times to simulate-
										COVariates(string)					///	-Baseline covariates, e.g, (sex 0.5 race -0.4)-
										TDE(string)							///	-Time dependent effects to interact with log(time)-
																			///
									/* Mixture */							///
										MIXture								///	-Simulate 2-component mixture-
										PMix(real 0.5)						///	-Mixture parameter-
																			///
									/* Competing risks */					///
										CR									///	-Simulate competing risks-
										NCR(string)							///	-Number of competing risks-
																			///
									/* Newton-Raphson scheme */				///
										CENTOL(real 0.0001)					///	-Tolerance of Newton-Raphson iterations-
										SHOWdiff							///	-Dislay error of N-R iterations-
																			///
									]


	tokenize `varlist'
		
	gen double `1' = .	
	
	forvalues i=1/2 {
		local l`i' : word `i' of `lambdas'
		local inputs "`inputs',`l`i''"
	}
	if "`distribution'"!="exp" {
		forvalues i=1/2 {
			local g`i' : word `i' of `gammas'				
			local inputs "`inputs',`g`i''"
		}				
	}

	mata: sim("`1'"`inputs',`pmix',obs,&ee(),linpred,U)
	
end	


mata:
	function ee(x,a$function_params,pmix,lin) {
		fn = $eqn_xb - a
		df = $eqn_dxb
		//fn = (pmix*exp(-l1*x^(g1)) + (1-pmix)*exp(-l2*x^(g2)))^(exp(lin)) - a
		//df = exp(lin)*((pmix*exp(-l1*x^(g1)) + (1-pmix)*exp(-l2*x^(g2)))^(exp(lin)-1))*(-pmix*l1*g1*x^(g1-1)*exp(-l1*x^(g1)) - (1-pmix)*l2'*g2*x^(g2-1)*exp(-l2*x^(g2)))
		return(fn,df)
	}
end


mata:
mata set matastrict off
  void sim(string scalar name
			$mataprogsyntax,
			real scalar pmix,
			real scalar obs,
			pointer (function) scalar ee,
			numeric matrix linpred,
			numeric matrix U)
{
    st_view(final=.,.,name)
	
	Y = J(obs,1,.)
	for(i=1;i<=obs;i++){
		bob=mm_nrroot(x=1,&ee(),0,1000000,a=U[i,1]$function_params,pmix,lin=linpred[i,1])
		if (x==.){
			bob=mm_nrroot(x=0.0000000000000001,&ee(),0,1000000,a=U[i,1]$function_params,pmix,lin=linpred[i,1])
		}
		Y[i,1] = x
		
	}
     
	final[,] = Y
	
}
end
