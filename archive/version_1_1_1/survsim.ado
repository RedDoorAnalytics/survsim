*! version 1.1.1 20sep2011 MJC

/*
History
MJC 20sep2011 v1.1.1: Exponential distribution added.
MJC 10sep2011 v1.1.0: Added Gompertz distribution. Time-dependent effects available for all models except mixture. showerror option added.
MJC 09sep2011 v1.0.1: Time dependent effects now allowed for standard Weibull.
*/

program define survsim
	version 11.0
	
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
									/* Frailty */							///
									/*	Frailty(string)						///
										FParameters(numlist min=2 max=2)*/	///
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
		
	/*	if "`frailty'"!="" & "`frailty'"!="gamma" {
			di as error "Unknown frailty"
			exit 198
		}*/
		
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
					local ind = `ind' + 5
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
					local ind = `ind' + 5
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
	/* Generate frailty's */
	
	/*	if "`frailty'"!="" {
			tempvar fvar
			local f1 : word 1 of `fparameters'
			local f2 : word 2 of `fparameters'
			gen `fvar' = rgamma(`f1',`f2')
			local fadd "+ `fvar'"
		}*/
	
	/********************************************************************************************************************************************************************************************/
	/* Preliminaries */

		tempvar u
		qui gen double `u' = runiform() if `use'

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
					local base_surv "( `pmix'*exp(-`l1'*`nr_time_old') + (1-`pmix')*exp(-`l2'*`nr_time_old') )"
					local eqn_xb "exp(-((-log(`base_surv')) `cov_linpred'))"
					local numer "(-`l1'*`pmix'*exp(-`l1'*`nr_time_old') -`l2'*(1-`pmix')*exp(-`l2'*`nr_time_old'))"
					local eqn_dxb "`eqn_xb' * `numer' `cov_linpred' / `base_surv'"
				}
				else if "`dist'"=="weibull" {
					local base_surv "( `pmix'*exp(-`l1'*`nr_time_old'^(`g1')) + (1-`pmix')*exp(-`l2'*`nr_time_old'^(`g2')) )"
					local eqn_xb "exp(-((-log(`base_surv')) `cov_linpred'))"
					local numer "(-`l1'*`g1'*`pmix'*`nr_time_old'^(`g1'-1)*exp(-`l1'*`nr_time_old'^(`g1')) -`l2'*`g2'*(1-`pmix')*`nr_time_old'^(`g2'-1)*exp(-`l2'*`nr_time_old'^(`g1')))"
					local eqn_dxb "`eqn_xb' * `numer' `cov_linpred' / `base_surv'"
				}
				else {
					local base_surv "( `pmix'*exp((`l1'/`g1')*(1-exp(`g1'*`nr_time_old'))) +  (1-`pmix')*exp((`l2'/`g2')*(1-exp(`g2'*`nr_time_old'))) )"
					local eqn_xb "exp(-((-log(`base_surv')) `cov_linpred'))"
					local numer "`pmix'*exp((`l1'/`g1')*(1-exp(`g1'*`nr_time_old'))) * (-`l1'*exp(`g1'*`nr_time_old'))  +  (1-`pmix')*exp((`l2'/`g2')*(1-exp(`g2'*`nr_time_old'))) * (-`l2'*exp(`g2'*`nr_time_old'))"
					local eqn_dxb "`eqn_xb' `cov_linpred' *(`numer') / `base_surv'"
				}
			}
			
		/* Competing risks */
			else if "`model'"=="cr"{
				local lambdastart = 0
				forvalues i=1/`ncr' {
					local l`i' : word `i' of `lambdas'
					local lambdastart = `lambdastart' + `l`i''
				}
				if "`dist'"!="exp" {
					local gammastart = 1
					forvalues i=1/`ncr' {
						local g`i' : word `i' of `gammas'				
						local gammastart = `gammastart' * `g`i''
					}				
				}
				
				/* Equations */
				if "`dist'"=="exp" {
					local eqn_xb "exp((-`l1'*`nr_time_old'^(1 `tde_linpred_1')) `cov_linpred_1' /(1 `tde_linpred_1'))"
					local eqn_dxb "(-`l1'*`nr_time_old'^((1 `tde_linpred_1')-1)`cov_linpred_1')"
					forvalues j=2/`ncr' {
						local eqn_xb "`eqn_xb' * exp((-`l`j''*`nr_time_old'^(1 `tde_linpred_`j'')) `cov_linpred_`j'' /(1 `tde_linpred_`j''))"			
						local eqn_dxb "`eqn_dxb' - (`l`j''*`nr_time_old'^((1 `tde_linpred_`j'')-1)`cov_linpred_`j'')"
					}
				}
				else if "`dist'"=="weibull" {
					local eqn_xb "exp((-`l1'*`g1'*`nr_time_old'^(`g1' `tde_linpred_1')) `cov_linpred_1' /(`g1' `tde_linpred_1'))"
					local eqn_dxb "(-`l1'*`g1'*`nr_time_old'^((`g1' `tde_linpred_1')-1)`cov_linpred_1')"
					forvalues j=2/`ncr' {
						local eqn_xb "`eqn_xb' * exp((-`l`j''*`g`j''*`nr_time_old'^(`g`j'' `tde_linpred_`j'')) `cov_linpred_`j'' /(`g`j'' `tde_linpred_`j''))"			
						local eqn_dxb "`eqn_dxb' - (`l`j''*`g`j''*`nr_time_old'^((`g`j'' `tde_linpred_`j'')-1)`cov_linpred_`j'')"
					}
				}
				else {
					local eqn_xb "exp((`l1'/(`g1' `tde_linpred_1'))`cov_linpred_1'*(1-exp((`g1' `tde_linpred_1')*`nr_time_old')))"
					local eqn_dxb "(`l1' `cov_linpred_1' *(-1)*exp((`g1' `tde_linpred_1')*`nr_time_old'))"
					forvalues j=2/`ncr' {
						local eqn_xb "`eqn_xb' * exp((`l`j''/(`g`j'' `tde_linpred_`j''))`cov_linpred_`j''*(1-exp((`g`j'' `tde_linpred_`j'')*`nr_time_old')))"			
						local eqn_dxb "`eqn_dxb' + (`l`j'' `cov_linpred_`j'' *(-1)*exp((`g`j'' `tde_linpred_`j'')*`nr_time_old'))"
					}
				}
				local eqn_dxb "(`eqn_xb') * (`eqn_dxb')"
			}
		
		if "`gammas'"!="" {
			if `gammastart'==0 {
				local gammastart = 0.0001
			}
		}
		if `lambdastart'==0 {
			local lambdastart = 0.0001
		}

	/********************************************************************************************************************************************************************************************/
	/* Standard Weibull/Gompertz calculations OR Starting values */
	
		if "`model'"!="cr" {
			if "`dist'"=="exp" {
				qui gen double `nr_time' 		= (-ln(`u')*(1 `tde_linpred')/(`lambdastart' `cov_linpred' `fadd'))^(1/(1 `tde_linpred')) if `use'
				qui gen double `nr_time_old' 	= `nr_time' if `use'
			}
			else if "`dist'"=="weibull" {
				qui gen double `nr_time' 		= (-ln(`u')*(`gammastart' `tde_linpred')/(`lambdastart'*`gammastart' `cov_linpred' `fadd'))^(1/(`gammastart' `tde_linpred')) if `use'
				qui gen double `nr_time_old' 	= `nr_time' if `use'
			}
			else{
				qui gen double `nr_time' 		= (1/(`gammastart' `tde_linpred'))*log(1-(((`gammastart' `tde_linpred')*log(`u'))/(`lambdastart' `cov_linpred'))) if `use'
				qui gen double `nr_time_old'	= `nr_time' if `use'
			}
		}
		else if "`model'"=="cr"{
			if "`dist'"=="exp" {
				qui gen double `nr_time' 		= -ln(`u')/(`lambdastart') if `use'
				qui gen double `nr_time_old' 	= `nr_time'	 if `use'
			}
			else if "`dist'"=="weibull" {
				qui gen double `nr_time' 		= (-ln(`u')/(`lambdastart'))^(1/`gammastart') if `use'
				qui gen double `nr_time_old' 	= `nr_time'	 if `use'
			}
			else {
				qui gen double `nr_time' 		= (1/(`gammastart'))*log(1-(`gammastart'*log(`u')/(`lambdastart'))) if `use'
				qui gen double `nr_time_old' 	= `nr_time'	 if `use'
			}
		}

	/********************************************************************************************************************************************************************************************/
	/* Newton-Raphson */

		if "`model'"=="mixture" | "`model'"=="cr" {
			
			/* Based on Therese Andersson's centile prediction option of stpm2_pred */
			local done 0
			while !`done' {
				qui gen double _nr_xb  = `eqn_xb' - `u' if `use'
				qui gen double _nr_dxb = `eqn_dxb' if `use'
				qui replace `nr_time'  = max(`nr_time_old' - _nr_xb/_nr_dxb,0.0000000000000001) if `use'
				qui gen double _error  = abs(`nr_time' - `nr_time_old') if `use'
				qui su _error if `use'
				`show' di in green `r(max)'
				if `r(max)'<`centol' {
					local done 1
				}
				else {
					drop _nr_xb _nr_dxb _error
					qui replace `nr_time_old' = `nr_time' if `use'
				}
			}
			cap drop _nr_xb _nr_dxb _error		
			
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Final variables */
	
		if "`model'"!="cr" {
			qui gen double `newvarname' = `nr_time'
		}
		else {
			local var1 : word 1 of `varlist'
			local var2 : word 2 of `varlist'
			qui gen double `var1' = `nr_time'
			
			/* Generate cause specific hazards at survival times */
			if "`dist'"=="exp" {
				local totalhaz "`l1'*(`nr_time')^(1 `tde_linpred_1'-1) `cov_linpred_1'"
				forvalues i=1/`ncr' {
					tempvar haz_`i'
					qui gen double `haz_`i'' = `l`i''*(`nr_time')^(1  `tde_linpred_`i''-1) `cov_linpred_`i''		
					if "`i'"!="1" {
						local totalhaz "`totalhaz' + `haz_`i''"
					}
				}			
			}
			else if "`dist'"=="weibull" {
				local totalhaz "`l1'*`g1'*(`nr_time')^(`g1' `tde_linpred_1'-1) `cov_linpred_1'"
				forvalues i=1/`ncr' {
					tempvar haz_`i'
					qui gen double `haz_`i'' = `l`i''*`g`i''*(`nr_time')^(`g`i''  `tde_linpred_`i''-1) `cov_linpred_`i''		
					if "`i'"!="1" {
						local totalhaz "`totalhaz' + `haz_`i''"
					}
				}
			}
			else {
				local totalhaz "`l1'*exp((`g1' `tde_linpred_1')*`nr_time') `cov_linpred_1'"
				forvalues i=1/`ncr' {
					tempvar haz_`i'
					qui gen double `haz_`i'' = `l`i''*exp((`g`i'' `tde_linpred_`i'')*`nr_time') `cov_linpred_`i''		
					if "`i'"!="1" {
						local totalhaz "`totalhaz' + `haz_`i''"
					}
				}
			}
			tempvar haz_all
			qui gen double `haz_all' = `totalhaz'
			forvalues i=1/`ncr' {
				tempvar p`i'
				qui gen double `p`i'' = `haz_`i''/`haz_all'
				local pvars "`pvars' `p`i''"
			}
			/* Event code */
			mata: pmatrix = st_data(.,tokens(st_local("pvars")))
			
			tempvar status
			qui gen `status'=.
			mata mata: genstatus("`status'",pmatrix)		
			qui gen `var2' = `status'
		}
	
	
	
	
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
			
	
