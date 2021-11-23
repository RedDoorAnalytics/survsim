/**********************************************************************/
/***																***/
/*** 	Program to simulate survival data 							***/
/***																***/
/*** 	Models:														***/
/***		Standard Weibull										***/
/***		2-component mixture Weibull								***/
/***		Competing risks with Weibull cause-specific hazards 	***/
/***																***/
/**********************************************************************/

// MJC 09sep2011 version 1.0.2 - Time dependent effects now allowed for standard weibull
// MJC 26Aug2011 version 1.0.1 - version 11 update.

program define survsim
	version 11.0
	
	syntax newvarname(min=1 max=2), 										///
										Lambdas(numlist min=1)				///	-Scale parameters for Weibull's-
										Gammas(numlist min=1)				///	-Shape parameters for Weibull's-
																			///
									[										/// -Options-
																			///
										N(string)							///	-Number of survival times to simulate-
										COVariates(string)					///	-Baseline covariates, e.g, (sex 0.5 race -0.4)-
										TDE(string)							///	-Time dependent effects to interact with log(time)-
										CENTOL(real 0.0001)					///	-Tolerance of Newton-Raphson iterations-
																			///
									/* Mixture Weib-Weib */					///
										MIXTURE								///	-Simulate 2-mixture Weibull data-
										PMix(real 0.5)						///	-Mixture parameter-
																			///
									/* Competing risks */					///
										CR									///	-Simulate competing risks data-
										NCR(string)							///	-Number of competing riskd-
																			///
									]

	
	local newvarname `varlist'
	local nvars : word count `varlist'

	/********************************************************************************************************************************************************************************************/
	/* Error checks */
		
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
		
		if "`mixture'"!="" & "`cr'"!="" & "`jm'"!=""  {
			di as error "Can only specify one of mixture/cr/jm"
			exit 198
		} 
		
		if `nvars'>1 & "`cr'"=="" {
			di as error "Two variables can only be specified when cr is used"
			exit 198
		}
		
		local nlambdas : word count `lambdas'
		local ngammas : word count `gammas'
		
		if "`cr'"=="" & "`jm'"=="" & "`mixture'"=="" {
			if "`nlambdas'"!="1"{
				di as error "Number of lambdas must be 1 under a standard Weibull model"
				exit 198
			}
			if "`ngammas'"!="1"{
				di as error "Number of gamma must be 2 under a standard Weibull model"
				exit 198
			}				
		}
		
		if "`mixture'"!="" {
			if "`nlambdas'"!="2" {
				di as error "Number of lambdas must be 2 under a mixture model"
				exit 198
			}
			if "`ngammas'"!="2" {
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
			if "`ngammas'"!="`ncr'" {
				di as error "Number of gammas must equal ncr"
				exit 198
			}
			
		}
		
		if "`tde'"!="" & ("`mixture'"!="" | "`cr'"!="") {
			di as error "option tde only allowed for standard Weibull"
			exit 198
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Defaults */
	
	if wordcount(`"`mixture' `cr'"')>0 {
		local model = trim("`mixture' `cr'")
	}
	else {
		local model "weibull"
	}
		
	/********************************************************************************************************************************************************************************************/
	/* Baseline covariates */	
	
	if "`covariates'" != "" {
		
		/* Standard Weibull or mixture Weibull*/
		if "`model'"=="weibull" | "`model'"=="mixture" {
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
		else if "`model'"=="cr" {
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
	
	/********************************************************************************************************************************************************************************************/
	/* Preliminaries */

		cap set obs `n'

		tempvar lhs u
		qui gen `u' = runiform() if `use'
		qui gen double `lhs' = 1-`u' if `use'				

	/********************************************************************************************************************************************************************************************/
	/* Equations for N-R */
	
		tempvar nr_time nr_time_old
	
		/* Standard Weibull */
			if "`model'"=="weibull" {
				local lambdastart : word 1 of `lambdas'
				local gammastart : word 1 of `gammas'	
			}
		
		/* Mixture Weibull */
			else if "`model'"=="mixture" {
				forvalues i=1/2 {
					local l`i' : word `i' of `lambdas'
					local g`i' : word `i' of `gammas'				
				}
				local lambdastart = `pmix'*`l1'+(1-`pmix')*`l2'
				local gammastart = `pmix'*`g1'+(1-`pmix')*`g2'			

				/* Equations */
				local base_surv "(`pmix'*exp(-`l1'*`nr_time_old'^(`g1')) + (1-`pmix')*exp(-`l2'*`nr_time_old'^(`g2'))          )"
				local eqn_xb "exp(-((-log(`base_surv')) `cov_linpred'))"
				local numer "(-`l1'*`g1'*`pmix'*`nr_time_old'^(`g1'-1)*exp(-`l1'*`nr_time_old'^(`g1')) -`l2'*`g2'*(1-`pmix')*`nr_time_old'^(`g2'-1)*exp(-`l2'*`nr_time_old'^(`g1')))"
				local eqn_dxb "`eqn_xb' * `numer' `cov_linpred' / `base_surv'"
			}
			
		/* Competing risks */
			else if "`model'"=="cr"{
				local lambdastart = 0
				local gammastart = 1
				forvalues i=1/`ncr' {
					local l`i' : word `i' of `lambdas'
					local g`i' : word `i' of `gammas'				
					local lambdastart = `lambdastart' + `l`i''
					local gammastart = `gammastart' * `g`i''
				}
				
				/* Equations */
				local eqn_xb "exp(-`l1'*`nr_time_old'^(`g1') `cov_linpred_1')"
				local eqn_dxb "(-`l1'*`g1'*`nr_time_old'^(`g1'-1)`cov_linpred_1')"
				forvalues j=2/`ncr' {
					local eqn_xb "`eqn_xb' * exp(-`l`j''*`nr_time_old'^(`g`j'') `cov_linpred_`j'')"			
					local eqn_dxb "`eqn_dxb' - (`l`j''*`g`j''*`nr_time_old'^(`g`j''-1) `cov_linpred_`j'')"
				}
				
				local eqn_dxb "(`eqn_xb') * (`eqn_dxb')"
			}

	/**********************************************************************************************************/
	/* Starting values */
	
		if "`model'"=="weibull" | "`model'"=="mixture" {
			qui gen double `nr_time' 		= (-ln(`lhs')*(`gammastart' `tde_linpred')/(`lambdastart'*`gammastart' `cov_linpred'))^(1/(`gammastart' `tde_linpred')) if `use'
			qui gen double `nr_time_old' 	= `nr_time' if `use'
		}
		else if "`model'"=="cr"{
			qui gen double `nr_time' 		= (-ln(`lhs')/(`lambdastart'))^(1/`gammastart') if `use'
			qui gen double `nr_time_old' 	= `nr_time'	 if `use'
		}

	/**********************************************************************************************************/
	/* Newton-Raphson */

		if "`model'"=="mixture" | "`model'"=="cr" {
			
			/* Based on Therese Andersson's centile prediction option of stpm2_pred */
			local done 0
			while !`done' {
				qui gen double _nr_xb = `eqn_xb' - `lhs' if `use'
				qui gen double _nr_dxb = `eqn_dxb' if `use'
				qui replace `nr_time' = max(`nr_time_old' - _nr_xb/_nr_dxb,0.0000000000000001) if `use'
				qui gen double _error = abs(`nr_time' - `nr_time_old') if `use'
				qui su _error if `use'
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
		
	/**********************************************************************************************************/
	/* Final variables */
	
		if "`model'"=="weibull" | "`model'"=="mixture" {
			qui gen double `newvarname' = `nr_time'
		}
		else if "`model'"=="cr" {
			local var1 : word 1 of `varlist'
			local var2 : word 2 of `varlist'
			qui gen double `var1' = `nr_time'
			
			/* Generate cause specific hazards at survival times */
			local totalhaz "`l1'*`g1'*(`nr_time')^(`g1'-1) `cov_linpred_1'"
			forvalues i=1/`ncr' {
				tempvar haz_`i'
				qui gen double `haz_`i'' = `l`i''*`g`i''*(`nr_time')^(`g`i''-1) `cov_linpred_`i''		
				if "`i'"!="1" {
					local totalhaz "`totalhaz' + `haz_`i''"
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
			
	
