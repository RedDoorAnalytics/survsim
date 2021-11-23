/**********************************************************************/
/***																***/
/*** 	Program to simulate survival data 							***/
/***																***/
/*** 	Models:														***/
/***		Standard Weibull										***/
/***		2-component mixture Weibull								***/
/***		Competeing risks with Weibull cause-specific hazards 	***/
/*** 		Joint model data										***/
/***																***/
/**********************************************************************/

program define survsim

	syntax newvarname(min=1 max=2), 	N(string)							///	-Number of survival times to simulate-
										LAMBDAS(numlist)					///	-Scale parameters for Weibull's-
										GAMMAS(numlist)						///	-Shape parameters for Weibull's-
																			///
									[										/// -Options-
																			///
										SCOVariates(string)					///	-Baseline covariates, e.g, (sex 0.5 race -0.4)-
										CENTOL(real 0.0001)					///	-Tolerance of Newton-Raphson iterations-
																			///
									/* Mixture Weib-Weib */					///
										MIXTURE								///	-Simulate 2-mixture Weibull data-
										PMix(real 0.5)						///	-Mixture parameter-
																			///
									/* Competing risks */					///
										CR									///	-Simulate competeing risks data-
										NCR(string)							///	-Number of competing riskd-
																			///
									/* Joint model - standard Weibull */	///
										JM									/// -Simulate joint model data-
										BETAS(numlist min=2 max=2)			///	-Intercept and slope values-
										SDS(numlist min=2 max=2)			///	-Standard deviations of random effects-
										CORRelation(real 0.5)				///	-Correlation between random effects
										ALPHA(real 0)						///
										LCOVariates(string)					///
										NLS									///
																			///
									]

	
	local newvarname `varlist'
	local nvars : word count `varlist'

	/********************************************************************************************************************************************************************************************/
	/* Error checks */
		
		if "`mixture'"!="" & "`cr'"!="" & "`jm'"!=""  {
			di as error "Can only specify one of mixture/cr/jm"
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
		
		if "`jm'"!="" {
		
			if `nlambdas'>1 {
				di as error "Number of lambdas must be 1 under a joint model"
				exit 198
			}
			if `ngammas'>1 {
				di as error "Number of gamma must be 2 under a joint model"
				exit 198
			}		
		
		
		
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Defaults */
	
	if wordcount(`"`mixture' `cr' `jm'"')>0 {
		local model = trim("`mixture' `cr' `jm'")
	}
	else {
		local model "weibull"
	}
		
	/********************************************************************************************************************************************************************************************/
	/* Baseline covariates */	
	
	if "`scovariates'" != "" {
		
		/* Standard Weibull or mixture Weibull*/
		if "`model'"=="weibull" | "`model'"=="mixture" {
			tokenize `scovariates'
			local ncovlist : word count `scovariates'
			local ncovvars = `ncovlist'/2
			cap confirm integer number `ncovvars'
			if _rc>0 {
				di as error "Variable/number missing in scovariates"
				exit 198
			}
			local ind = 1
			forvalues i=1/`ncovvars' {
				cap confirm var ``ind''
				if _rc {
					local errortxt "invalid scovariates(... ``ind'' ``=`ind'+1'' ...)"
					local error = 1
				}
				cap confirm num ``=`ind'+1''
				if _rc {
					local errortxt "invalid scovariates(... ``ind'' ``=`ind'+1'' ...)"
					local error = 1
				}
				tempvar vareffect`i'
				gen double `vareffect`i'' = ``ind''*``=`ind'+1''
	
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
			tokenize `scovariates'
			local ncovlist : word count `scovariates'	
			local ncovvars = `ncovlist'/`=`ncr'+1'
			cap confirm integer number `ncovvars'
			if _rc>0 {
				di as error "Variable/number missing in scovariates"
				exit 198
			}

			local ind = 1
			forvalues i=1/`ncovvars' {
				cap confirm var ``ind''
				if _rc {
					local errortxt "invalid scovariates(... ``ind'' ``=`ind'+1'' ...)"
					local error = 1
				}
				forvalues j=1/`ncr' {
					cap confirm num ``=`ind'+`j'''
					if _rc {
						local errortxt "invalid scovariates(... ``ind'' ``=`ind'+`j''' ...)"
						local error = 1
					}
					/* Create effect for ith variable and jth risk */
					tempvar vareffect_`i'_`j'
					gen double `vareffect_`i'_`j'' = ``ind''*``=`ind'+`j'''
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

	/*** Baseline covariates for joint model ***/
	if "`jm'"!="" & ("`scovariates'"!="" | "`lcovariates'"!="") {
	
		/* Longitudinal covariates */
		if "`lcovariates'"!="" {
			tokenize `lcovariates'
			local nlongcovlist : word count `lcovariates'
			local nlongcovvars = `nlongcovlist'/2
			cap confirm integer number `nlongcovvars'
			if _rc>0 {
				di as error "Variable/number missing in covariates"
				exit 198
			}						

			local ind = 1
			forvalues i=1/`nlongcovvars' {
				cap confirm var ``ind''
				if _rc {
					local errortxt "invalid lcovariates(... ``ind'' ``=`ind'+1'' ...)"
					local error = 1
				}
				cap confirm num ``=`ind'+1''
				if _rc {
					local errortxt "invalid lcovariates(... ``ind'' ``=`ind'+1'' ...)"
					local error = 1
				}
				tempvar longvareffect`i'
				gen double `longvareffect`i'' = ``ind''*``=`ind'+1''
	
				local ind = `ind' + 2
			}
			if "`error'"=="1" {
				di as error "`errortxt'"
				exit 198
			}
			local longcov_linpred "`longvareffect1'"
			if `nlongcovvars'>1 {
				forvalues k=2/`nlongcovvars' {
					local longcov_linpred "`longcov_linpred' + `longvareffect`k''"
				}
			}
			local longcov_linpred "+ `longcov_linpred'"
		}
		
		/* Survival covariates */
		if "`scovariates'"!="" {
			tokenize `scovariates'
			local nsurvcovlist : word count `lcovariates'
			local nsurvcovvars = `nsurvcovlist'/2
			cap confirm integer number `nsurvcovvars'
			if _rc>0 {
				di as error "Variable/number missing in scovariates"
				exit 198
			}						

			local ind = 1
			forvalues i=1/`nsurvcovvars' {
				cap confirm var ``ind''
				if _rc {
					local errortxt "invalid scovariates(... ``ind'' ``=`ind'+1'' ...)"
					local error = 1
				}
				cap confirm num ``=`ind'+1''
				if _rc {
					local errortxt "invalid scovariates(... ``ind'' ``=`ind'+1'' ...)"
					local error = 1
				}
				tempvar survvareffect`i'
				gen double `survvareffect`i'' = ``ind''*``=`ind'+1''
	
				local ind = `ind' + 2
			}
			if "`error'"=="1" {
				di as error "`errortxt'"
				exit 198
			}
			local survcov_linpred "`survvareffect1'"
			if `nsurvcovvars'>1 {
				forvalues k=2/`nsurvcovvars' {
					local survcov_linpred "`survcov_linpred' + `survvareffect`k''"
				}
			}
			local survcov_linpred "+ `survcov_linpred'"
		}
	}
	
	/********************************************************************************************************************************************************************************************/
	/* Preliminaries */

		cap set obs `n'

		tempvar lhs u
		qui gen `u' = runiform()
		qui gen double `lhs' = 1-`u'				

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

		/* Joint model */
			else if "`model'"=="jm"{
			
				local lambdastart : word 1 of `lambdas'
				local gammastart : word 1 of `gammas'			
			
				tempname betasmat sdsmat corrmat
				mat `betasmat' = J(1,2,.)
				mat `sdsmat' = J(1,2,.)
				forvalues i=1/2 {
					local b`i' : word `i' of `betas'
					local sd`i' : word `i' of `sds'
					mat `betasmat'[1,`i'] = `b`i''
					mat `sdsmat'[1,`i'] = `sd`i''
				}
				mat `corrmat' =  (1, `correlation', 1)
				
				tempvar int slope
				drawnorm `int' `slope', means(`betasmat') sds(`sdsmat') corr(`corrmat') cstorage(lower)
		
				local longlinpred "(`int') + (`slope') * `nr_time_old' `longcov_linpred'"

				local cumhazard "`lambdastart'*`nr_time_old'^(`gammastart') * exp(`alpha'*(`longlinpred') `survcov_linpred')"
			
				local eqn_xb "exp(-(`cumhazard'))"
				
				local eqn_dxb "(-1)*(`eqn_xb')*((`cumhazard')*(`alpha')*(`slope') + (`cumhazard')*(`gammastart')*(`nr_time_old'^(-1)))"
			
			/*	tempname betasmat sdsmat corrmat
				mat `betasmat' = J(1,2,.)
				mat `sdsmat' = J(1,2,.)
				forvalues i=1/2 {
					local b`i' : word `i' of `betas'
					local sd`i' : word `i' of `sds'
					mat `betasmat'[1,`i'] = `b`i''
					mat `sdsmat'[1,`i'] = `sd`i''
				}
				mat `corrmat' =  (1, `correlation', 1)
				
				tempvar int slope
				drawnorm `int' `slope', means(`betasmat') sds(`sdsmat') corr(`corrmat') cstorage(lower)

				local longlinpred "(`int') + (`slope') * `nr_time_old' `longcov_linpred'"
			
				local hazard "(`lambdastart')*exp(`alpha'*(`longlinpred') `survcov_linpred')"
		
				local cumhazard "(`hazard')/(`alpha'*(`slope')) - (`lambdastart'*exp(`alpha'*(`int' `longcov_linpred') `survcov_linpred')/(`alpha'*(`slope')))"
				
				local eqn_xb "exp((-1)*(`cumhazard'))"
			
				local eqn_dxb "(`eqn_xb')*((-1)*(`hazard'))"
			*/
			
			}
			
	/**********************************************************************************************************/
	/* Starting values */
	
		if "`model'"=="weibull" | "`model'"=="mixture" {
			qui gen double `nr_time' 		= (-ln(`lhs')/(`lambdastart' `cov_linpred'))^(1/`gammastart')
			qui gen double `nr_time_old' 	= `nr_time'
		}
		else if "`model'"=="cr"{
			qui gen double `nr_time' 		= (-ln(`lhs')/(`lambdastart'))^(1/`gammastart')
			qui gen double `nr_time_old' 	= `nr_time'	
		}
		else if "`model'"=="jm" {
			qui gen double `nr_time' 		= (-ln(`lhs')/(`lambdastart'))^(1/`gammastart')
			qui gen double `nr_time_old' 	= `nr_time'
		}
/*		gen time = `nr_time_old'
	gen haz = `hazard'
	gen cumhaz1 = (`hazard')/(`alpha'*(`slope'))*/
	/**********************************************************************************************************/
	/* Newton-Raphson */

		if "`model'"=="mixture" | "`model'"=="cr" | ("`model'"=="jm" & "`nls'"=="") {
		
			local done 0
			while !`done' {
				qui gen double nr_xb = `eqn_xb' - `lhs'
				qui gen double nr_dxb = `eqn_dxb'
				qui replace `nr_time' = max(`nr_time_old' - nr_xb/nr_dxb,0.0000000000000001)
				qui gen double error = abs(`nr_time' - `nr_time_old')
				 su error
				if `r(max)'<`centol' {
					local done 1
				}
				else {
					drop nr_xb nr_dxb error
					qui replace `nr_time_old' = `nr_time'
				}
			*local done 1
			}
			cap drop nr_xb nr_dxb error		
			
		}
		else if "`model'"=="jm" & "`nls'"!="" {

			forvalues i=1/`n' {
			
				tempvar y
				gen `y'=0 in 1/2
				replace `y' = 1 in 1
			
				replace `int' = `int' `longcov_linpred'
				tempvar scov
				gen `scov' = 0 `survcov_linpred'
				
				nl jmww @ `y' in 1/2, parameters(A B) initial(A 1 B 1) 	///
						unif(`=`u'[`i']') lambda1(`lambdastart') gamma1(`gammastart') 	///
						alpha(`alpha') beta0(`=`int'[`i']') beta1(`=`slope'[`i']') scov(`=`scov'[`i']')

				replace `nr_time' = _b[_cons] in `i'
				drop `y'
		
			}
		
		}
		
	/**********************************************************************************************************/
	/* Final variables */
	
		if "`model'"=="weibull" | "`model'"=="mixture" | "`model'"=="jm"{
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
			
	
