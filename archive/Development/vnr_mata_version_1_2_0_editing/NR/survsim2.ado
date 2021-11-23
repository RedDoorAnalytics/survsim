*! version 1.2.0 Feb2012

/*
History
MJC Feb2012
MJC 15Nov2011 v1.1.2: Fixed bug when generating covariate tempvars with competing risks.
MJC 20sep2011 v1.1.1: Exponential distribution added.
MJC 10sep2011 v1.1.0: Added Gompertz distribution. Time-dependent effects available for all models except mixture. showerror option added.
MJC 09sep2011 v1.0.1: Time dependent effects now allowed for standard Weibull.
*/
program define survsim2

	version 11.2
	
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

	local newsyntax `0'
	cap pr drop survsim2_nr

	local newvarname `varlist'
	local nvars : word count `varlist'

	/********************************************************************************************************************************************************************************************/
	/* Distribution */
	
	if "`distribution'"=="" {
		local dist "weibull"
	}
	else {
		local l = length("`distribution'")
		if substr("exponential",1,max(1,`l')) == "`distribution'" {
			local dist "exp"
		}
		else if substr("gompertz",1,max(3,`l')) == "`distribution'" {
			local dist "gompertz"
		}
		else if substr("weibull",1,max(1,`l')) == "`distribution'" {
			local dist "weibull"
		}
		else {
			di as error "Unknown distribution"
			exit 198
		}
	}
	
	/********************************************************************************************************************************************************************************************/
	/* Error checks */
		
		foreach l of numlist `lambdas' {
			if `l'<0 {
				di as error "lambdas must be > 0"
				exit 198
			}
		}
		
		if "`dist'"=="exp" & "`gammas'"!="" {
			di as error "gammas cannot be specified with distribution(exponential)"
			exit 198
		}
		
		if "`dist'"=="weibull" {
			foreach g of numlist `gammas' {
				if `g'<0 {
					di as error "gammas must be > 0"
					exit 198
				}
			}
		}
		
		if "`n'"!="" {
			cap confirm integer number `n'
			if _rc>0 {
				di as error "n must be an integer"
				exit 198
			}
		}
		else {
			local n = _N
		}
		
		cap set obs `n'
		mata: obs = `n'
		tempvar use
		gen `use'= _n <= `n'
		
		if "`mixture'"!="" & "`cr'"!="" {
			di as error "Can only specify one of mixture/cr"
			exit 198
		} 
		
		if `nvars'>1 & "`cr'"=="" {
			di as error "Two variables can only be specified when cr is used"
			exit 198
		}
		
		local nlambdas : word count `lambdas'
		local ngammas  : word count `gammas'
				
		if "`cr'"=="" & "`mixture'"=="" {
			if "`nlambdas'"!="1"{
				di as error "Number of lambda's must be 1 under a standard parametric model"
				exit 198
			}
			if "`ngammas'"!="1" & "`gammas'"!="" {
				di as error "Number of gamma's must be 1 under a standard parametric model"
				exit 198
			}				
		}
		
		if "`mixture'"!="" {
			if "`nlambdas'"!="2" {
				di as error "Number of lambdas must be 2 under a mixture model"
				exit 198
			}
			if "`ngammas'"!="2" & "`gammas'"!="" {
				di as error "Number of gamma must be 2 under a mixture model"
				exit 198
			}		
		}
		
		if "`cr'"!="" & "`nvars'"!="2" {
			di as error "2 variables must be specified for cr"
			exit 198
		}
		
		if "`cr'"!="" {
		
			if "`ncr'"=="" {
				di as error "ncr must be specified"
				exit 198
			}
			
			cap confirm integer number `ncr'
			if _rc>0 {
				di as error "ncr must be an integer"
				exit 198
			}
			
			if `ncr'<2 {
				di as error "ncr must be >1"
				exit 198
			}
			
			if "`nlambdas'"!="`ncr'"{
				di as error "Number of lambdas must equal ncr"
				exit 198
			}
			if "`ngammas'"!="`ncr'" & "`dist'"!="exp" {
				di as error "Number of gammas must equal ncr"
				exit 198
			}
			
		}
		
		if "`tde'"!="" & "`mixture'"!="" {
			di as error "time-dependent effects not allowed under a mixture model"
			exit 198
		}
		
		if ("`mixture'"=="" & "`cr'"=="") & "`showdiff'"!="" {
			di as error "showdiff only available when mixture or cr is specified"
			exit 198
		}
				
	/********************************************************************************************************************************************************************************************/
	/* Defaults */
	
		if wordcount(`"`mixture' `cr'"')>0 {
			local model = trim("`mixture' `cr'")
		}
		else {
			local model "`dist'"
		}
		
		if "`showdiff'"!="" {
			local show "noisily"
		}
		else {
			local show "quietly"
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Baseline covariates and time-dependent effects */	
	
		if "`covariates'"!="" {
			
			/* Standard parametric or mixture*/
			if "`model'"!="cr" {
				tokenize `covariates'
				local ncovlist : word count `covariates'
				local ncovvars = `ncovlist'/2
				cap confirm integer number `ncovvars'
				if _rc>0 {
					di as error "Variable/number missing in covariates"
					exit 198
				}
				local ind = 1
				forvalues i=1/`ncovvars' {
					cap confirm var ``ind''
					if _rc {
						local errortxt "invalid covariates(... ``ind'' ``=`ind'+1'' ...)"
						local error = 1
					}
					cap confirm num ``=`ind'+1''
					if _rc {
						local errortxt "invalid covariates(... ``ind'' ``=`ind'+1'' ...)"
						local error = 1
					}
					tempvar vareffect`i'
					gen double `vareffect`i'' = ``ind''*``=`ind'+1'' if `use'
		
					local ind = `ind' + 2
				}
				if "`error'"=="1" {
					di as error "`errortxt'"
					exit 198
				}
				local cov_linpred "`vareffect1'"
				if `ncovvars'>1 {
					forvalues k=2/`ncovvars' {
						local cov_linpred "`cov_linpred' + `vareffect`k''"
					}
				}
				local cov_linpred "* exp(`cov_linpred')"
			}
			
			/* Competing risks */
			else {
				tokenize `covariates'
				local ncovlist : word count `covariates'	
				local ncovvars = `ncovlist'/`=`ncr'+1'
				cap confirm integer number `ncovvars'
				if _rc>0 {
					di as error "Variable/number missing in covariates"
					exit 198
				}

				local ind = 1
				forvalues i=1/`ncovvars' {
					cap confirm var ``ind''
					if _rc {
						local errortxt "invalid covariates(... ``ind'' ``=`ind'+1'' ...)"
						local error = 1
					}
					forvalues j=1/`ncr' {
						cap confirm num ``=`ind'+`j'''
						if _rc {
							local errortxt "invalid covariates(... ``ind'' ``=`ind'+`j''' ...)"
							local error = 1
						}
						/* Create effect for ith variable and jth risk */
						tempvar vareffect_`i'_`j'
						gen double `vareffect_`i'_`j'' = ``ind''*``=`ind'+`j''' if `use'
					}
					local ind = `ind' + `ncr' + 1
				}
				if "`error'"=="1" {
					di as error "`errortxt'"
					exit 198
				}
				forvalues k=1/`ncr' {
					local cov_linpred_`k' "`vareffect_1_`k''"
				}
				if `ncovvars'>1 {
					forvalues p=2/`ncovvars' {
						forvalues m=1/`ncr' {
							local cov_linpred_`m' "`cov_linpred_`m'' + `vareffect_`p'_`m''"
						}
					}
				}
				forvalues k=1/`ncr' {
					local cov_linpred_`k' "* exp(`cov_linpred_`k'')"
				}			
					
			}
		}
		
		if "`tde'"!="" {
			/* Standard parametric or mixture*/
			if "`model'"!="cr" {
				tokenize `tde'
				local ntde : word count `tde'	
				local ntdevars = `ntde'/2
				cap confirm integer number `ntdevars'
				if _rc>0 {
					di as error "Variable/number missing in tde"
					exit 198
				}

				local ind = 1
				forvalues i=1/`ntdevars' {
					cap confirm var ``ind''
					if _rc {
						local errortxt "invalid tde(... ``ind'' ``=`ind'+1'' ...)"
						local error = 1
					}
					cap confirm num ``=`ind'+1''
					if _rc {
						local errortxt "invalid tde(... ``ind'' ``=`ind'+1'' ...)"
						local error = 1
					}
					tempvar tdeeffect`i'
					gen double `tdeeffect`i'' = ``ind''*``=`ind'+1'' if `use'

					local ind = `ind' + 2
				}
				if "`error'"=="1" {
					di as error "`errortxt'"
					exit 198
				}
				local tde_linpred "`tdeeffect1'"
				if `ntdevars'>1 {
					forvalues k=2/`ntdevars' {
						local tde_linpred "`tde_linpred' + `tdeeffect`k''"
					}
				}
				local tde_linpred "+ `tde_linpred'"
			}
			
			/* Competing risks */
			else {
				tokenize `tde'
				local ntdelist : word count `tde'	
				local ntdevars = `ntdelist'/`=`ncr'+1'
				cap confirm integer number `ntdevars'
				if _rc>0 {
					di as error "Variable/number missing in tde"
					exit 198
				}
				local ind = 1
				forvalues i=1/`ntdevars' {
					cap confirm var ``ind''
					if _rc {
						local errortxt "invalid tde(... ``ind'' ``=`ind'+1'' ...)"
						local error = 1
					}
					forvalues j=1/`ncr' {
						cap confirm num ``=`ind'+`j'''
						if _rc {
							local errortxt "invalid tde(... ``ind'' ``=`ind'+`j''' ...)"
							local error = 1
						}
						/* Create effect for ith variable and jth risk */
						tempvar tdeeffect_`i'_`j'
						gen double `tdeeffect_`i'_`j'' = ``ind''*``=`ind'+`j''' if `use'
					}
					local ind = `ind' + `ncr' + 1
				}
				if "`error'"=="1" {
					di as error "`errortxt'"
					exit 198
				}
				forvalues k=1/`ncr' {
					local tde_linpred_`k' "`tdeeffect_1_`k''"
				}
				if `ntdevars'>1 {
					forvalues p=2/`ntdevars' {
						forvalues m=1/`ncr' {
							local tde_linpred_`m' "`tde_linpred_`m'' + `tdeeffect_`p'_`m''"
						}
					}
				}
				forvalues k=1/`ncr' {
					local tde_linpred_`k' "+ `tde_linpred_`k''"
				}			
			}
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Preliminaries */

		tempvar u
		qui gen double `u' = runiform() if `use'
		if "`mixture'"!="" {	
			mata: U = st_data(.,"`u'","`use'")
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Equations for N-R */
	
		tempvar nr_time nr_time_old
	
		/* Standard exponential */
			if "`model'"=="exp" {
				local lambdastart : word 1 of `lambdas'
			}

		/* Standard Weibull/Gompertz */
			else if "`model'"=="weibull" | "`model'"=="gompertz" {
				local lambdastart : word 1 of `lambdas'
				local gammastart  : word 1 of `gammas'	
			}
			
		/* Mixture Weibull */
			else if "`model'"=="mixture" {
				forvalues i=1/2 {
					local l`i' : word `i' of `lambdas'
				}
				local lambdastart = `pmix'*`l1'+(1-`pmix')*`l2'
				if "`dist'"!="exp" {
					forvalues i=1/2 {
						local g`i' : word `i' of `gammas'				
					}				
					local gammastart = `pmix'*`g1'+(1-`pmix')*`g2'			
				}

				/* Equations */
				if "`dist'"=="exp" {
					local base_surv "( pmix*exp(-l1*x) + (1-pmix)*exp(-l2*x) )"
					global eqn_xb exp(-((-log(`base_surv')) `cov_linpred'))
					local numer "(-l1*pmix*exp(-l1*x) -l2*(1-pmix)*exp(-l2*x))"
					global eqn_dxb $eqn_xb * `numer' `cov_linpred' / `base_surv'
					
					global function_params ",l1,l2"
					global mataprogsyntax ",real scalar l1, real scalar l2"
				}
				else if "`dist'"=="weibull" {
					local base_surv "(pmix*exp(-l1*x^(g1))+(1-pmix)*exp(-l2*x^(g2)))"
					global eqn_xb exp(-((-log(`base_surv')) `cov_linpred'))
					local numer "(-l1*g1*pmix*x^(g1-1)*exp(-l1*x^(g1))-l2*g2*(1-pmix)*x^(g2-1)*exp(-l2*x^(g1)))"
					global eqn_dxb $eqn_xb * `numer' `cov_linpred' / `base_surv'
					
					global function_params ",l1,l2,g1,g2"
					global mataprogsyntax ",real scalar l1, real scalar l2,real scalar g1,real scalar g2"
				}
				else {
					local base_surv "(pmix*exp((l1/g1)*(1-exp(g1*x))) + (1-pmix)*exp((l2/g2)*(1-exp(g2*x))))"
					global eqn_xb exp(-((-log(`base_surv')) `cov_linpred'))
					local numer "pmix*exp((l1/g1)*(1-exp(g1*x)))*(-l1*exp(g1*x)) + (1-pmix)*exp((l2/g2)*(1-exp(g2*x))) * (-l2*exp(g2*x))"
					global eqn_dxb $eqn_xb `cov_linpred' *(`numer') / `base_surv'
				}
			}
			
			mata: linpred = J(`n',1,0)
			
			survsim2_nr `newsyntax'
		
	
end

/* Mata program to generate status indicator under a competing risks model */
mata:
mata set matastrict off
  void genstatus(string scalar name,		///
				 numeric matrix pmatrix)
  {
    st_view(final=.,.,name)
    n=strtoreal(st_local("n"))
	
	for (i=1; i<=n; i++) {
		p = pmatrix[i,.]
		final[i,.] = rdiscrete(1,1,p)
	}
  }
end



