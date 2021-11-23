*! version 1.2.0

/*
Notes
-> add quadrature and root finding
-> add maxt() option
*/

/*
History
MJC ?????2012
MJC 15Nov2011 v1.1.2: Fixed bug when generating covariate tempvars with competing risks.
MJC 20sep2011 v1.1.1: Exponential distribution added.
MJC 10sep2011 v1.1.0: Added Gompertz distribution. Time-dependent effects available for all models except mixture. showerror option added.
MJC 09sep2011 v1.0.1: Time dependent effects now allowed for standard Weibull.
*/

program define survsim
	version 11.2
	syntax newvarname(min=1 max=2), 										///
																			///
									[										/// -Options-
										BASEHazard(string)					///	-Define the baseline hazard function-
										Lambdas(numlist min=1)				///	-Scale parameters-
										Gammas(numlist min=1)				///	-Shape parameters-
										Distribution(string)				/// -Parametric distribution-
																			///
										N(string)							///	-Number of survival times to simulate-
										COVariates(string)					///	-Baseline covariates, e.g, (sex 0.5 race -0.4)-
										TDE(string)							///	-Time dependent effects to interact with log(time)-
										TDEFUNCtion(string)					/// -Function of time to interact with covariates defined in tde-
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
									/* Quadrature */						///
										NODES(real 15)						///
										GL									///
										ALPHA(real 0)						///
										BETA(real 0)						///
									]
									
	if "`basehazard'"!="" {
		cap pr drop survsim_nr
	}
	tokenize varlist
	local newvarname `varlist'
	local nvars : word count `varlist'

	/********************************************************************************************************************************************************/
	/* Error checks */
		
		if "`basehazard'"!="" & wordcount(`"`lambdas' `gammas' `distribution' `cr' `mixture'"')>0 {
			di as error "basehazard cannot be specified"
			exit 198
		}
		
		if "`basehazard'"=="" & "`distribution'"=="" {
			di as error "distribution must be specified"
			exit 198
		}
		
		if "`distribution'"!="" {
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

		if "`lambdas'"!="" {
			foreach l of numlist `lambdas' {
				if `l'<0 {
					di as error "lambdas must be > 0"
					exit 198
				}
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
				
		if "`cr'"=="" & "`mixture'"=="" & "`basehazard'"=="" {
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

	/********************************************************************************************************************************************************/
	/* Prelims */
		
		/* Define model choice */	
		if wordcount(`"`mixture' `cr'"')>0 {
			local model = trim("`mixture' `cr'")
		}
		else if "`basehazard'"!="" {
			local model user
		}
		else {
			local model "`dist'"
		}
		
		/* Define whether NR scheme is to be used */
		if ("`basehazard'"!="" | "`tdefunction'"!="") & "`cr'"=="" {
			local nr yes
		}
		
		/* Sample size */		
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
		
		/* Show difference in NR iterations - need to add mata pass */
		if "`showdiff'"!="" {
			local show "noisily"
		}
		else {
			local show "quietly"
		}
		
		/* Marksample */	
		tempvar use
		qui gen `use' = _n <= `n'
		
		/* Generate random draws */
		tempvar u
		qui gen double `u' = runiform() if `use'
		
		if "`nr'"!="" {
			mata: n = `n'									/* Sample size */
			mata: centol = `centol'							/* Tolerance for NR scheme */
			qui gen double `newvarname' = .
			global newvarname `newvarname'
			mata: u = st_data(.,"`u'","`use'")
		}

	/********************************************************************************************************************************************************/
	/* Time-independent and time-dependent covariates */

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
						di as error "invalid covariates(... ``ind'' ``=`ind'+1'' ...)"
						exit 198
					}
					cap confirm num ``=`ind'+1''
					if _rc {
						di as error "invalid covariates(... ``ind'' ``=`ind'+1'' ...)"
						exit 198
					}
					tempvar vareffect`i'
					gen double `vareffect`i'' = ``ind''*``=`ind'+1'' if `use'
		
					local ind = `ind' + 2
				}
				local cov_linpred "`vareffect1'"
				if `ncovvars'>1 {
					forvalues k=2/`ncovvars' {
						local cov_linpred "`cov_linpred' + `vareffect`k''"
					}
				}
				if "`nr'"=="" {	
					local cov_linpred "* exp(`cov_linpred')"
				}
				else {
					tempvar linpred
					gen double `linpred' = `cov_linpred'
					mata: linpred = st_data(.,"`linpred'","`use'")
				}
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
						di as error "invalid covariates(... ``ind'' ``=`ind'+1'' ...)"
						exit 198
					}
					forvalues j=1/`ncr' {
						cap confirm num ``=`ind'+`j'''
						if _rc {
							di as error "invalid covariates(... ``ind'' ``=`ind'+`j''' ...)"
							exit 198
						}
						/* Create effect for ith variable and jth risk */
						tempvar vareffect_`i'_`j'
						gen double `vareffect_`i'_`j'' = ``ind''*``=`ind'+`j''' if `use'
					}
					local ind = `ind' + `ncr' + 1
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

			/* Standard parametric or mixture */
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
						di as error "invalid tde(... ``ind'' ``=`ind'+1'' ...)"
						exit 198
					}
					cap confirm num ``=`ind'+1''
					if _rc {
						di as error  "invalid tde(... ``ind'' ``=`ind'+1'' ...)"
						exit 198
					}
					tempvar tdeeffect`i'
					gen double `tdeeffect`i'' = ``ind''*``=`ind'+1'' if `use'

					local ind = `ind' + 2
				}
				local tde_linpred "`tdeeffect1'"
				if `ntdevars'>1 {
					forvalues k=2/`ntdevars' {
						local tde_linpred "`tde_linpred' + `tdeeffect`k''"
					}
				}
				
				if "`nr'"=="" {
					local tde_linpred "+ `tde_linpred'"
				}
				else {
					tempvar tdevar_linpred
					gen double `tdevar_linpred' = `tde_linpred'
					mata: tdelinpred = st_data(.,"`tdevar_linpred'","`use'")	
				}
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
						di as error "invalid tde(... ``ind'' ``=`ind'+1'' ...)"
						exit 198
					}
					forvalues j=1/`ncr' {
						cap confirm num ``=`ind'+`j'''
						if _rc {
							di as error "invalid tde(... ``ind'' ``=`ind'+`j''' ...)"
							exit 198
						}
						/* Create effect for ith variable and jth risk */
						tempvar tdeeffect_`i'_`j'
						gen double `tdeeffect_`i'_`j'' = ``ind''*``=`ind'+`j''' if `use'
					}
					local ind = `ind' + `ncr' + 1
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

	/********************************************************************************************************************************************************/
	/* Lambda and Gamma values from syntax */
		
		/* Standard parametric */
		if "`model'"=="exp" {
			local l1 : word 1 of `lambdas'
			local lambdastart = `l1'
		}
		else if ("`model'"=="weibull" | "`model'"=="gompertz") {
			local l1 : word 1 of `lambdas'
			local g1 : word 1 of `gammas'	
			local lambdastart = `l1'
			local gammastart = `g1'
		}
		
		/* 2-component mixture */
		if "`model'"=="mixture" {
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
		}
		
		/* CR */
		if "`model'"=="cr" {
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
		}
		
	/********************************************************************************************************************************************************/
	/* Newton-Raphson set-up: needed if user-defined hazard function, tdefunction specified, mixture or cr model */
		
		if "`nr'"!="" {
			
			/* Baseline hazard function */
				
				/* User specified */
				if "`basehazard'"!="" {	
					mata: st_global("newhaz", subinstr("`basehazard'","#t","x",.))
					mata: st_global("newhaz", subinstr("$newhaz","*",":*",.))			//could avoid these next lines if defined a mata function and used mm_nrroot, depends on syntax
					mata: st_global("newhaz", subinstr("$newhaz","+",":+",.))
					mata: st_global("newhaz", subinstr("$newhaz","-",":-",.))
					mata: st_global("newhaz", subinstr("$newhaz","^",":^",.))
					mata: st_global("newhaz", subinstr("$newhaz","/",":/",.))
				}
				
				/* 2-component mixture */
				if "`model'"=="mixture" {
					if ""=="exp" {
						global newhaz (`l1':*`pmix':*exp(-`l1':*x) :+ `l2':*(1:-`pmix'):*exp(-`l2':*x)) :/ (`pmix':*exp(-`l1':*x) :+ (1:-`pmix'):*exp(-`l2':*x))
					}
					else if ""=="weibull" {
						global newhaz (`l1':*`g1':*x:^(`g1':-1):*`pmix':*exp(-`l1':*x:^(`g1')) :+ `l2':*`g2':*x:^(`g2':-1):*(1:-`pmix'):*exp(-`l2':*x:^(`g2'))) :/ (`pmix':*exp(-`l1':*x:^(`g1')) :+ (1:-`pmix'):*exp(-`l2':*x:^(`g2')))
					}
					else {
						local numer `l1':*exp(`g1':*x):*`pmix':*exp((-`l1'):/(`g1'):*(exp(`g1':*x):-1)) :+ `l2':*exp(`g2':*x):*(1:-`pmix'):*exp((-`l2'):/(`g2'):*(exp(`g2':*x):-1))
						local demon `pmix':*exp((-`l1'):/(`g1'):*(exp(`g1':*x):-1)) :+ (1:-`pmix'):*exp((-`l2'):/(`g2'):*(exp(`g2':*x):-1))
						global newhaz (`numer'):/(`denom')
					}				
				}

				/* Standard parametric model */
				if "`model'"=="exp" {
					global newhaz `l1'
				}
				else if "`model'"=="weibull" {
					global newhaz `l1':*`g1':*x:^(`g1':-1)
				}
				else if "`model'"=="gompertz"{
					global newhaz `l1':*exp(`g1':*x)
				}
				
			/* Time-dependent effects */
				
				/* If time-dependent effects then need to add interaction with generated coefficients x covariates x time */
				if "`tde'"!="" & "`tdefunction'"!="" {
					mata: st_local("tdenewfunc", subinstr("`tdefunction'","#t","x",.))
					mata: st_global("tdenewfunc", subinstr("$tdenewfunc","*",":*",.))
					mata: st_global("tdenewfunc", subinstr("$tdenewfunc","+",":+",.))
					mata: st_global("tdenewfunc", subinstr("$tdenewfunc","-",":-",.))
					mata: st_global("tdenewfunc", subinstr("$tdenewfunc","^",":^",.))
					mata: st_global("tdenewfunc", subinstr("$tdenewfunc","/",":/",.))
					global newhaz ($newhaz) :* exp(tdelinpred :* (`tdenewfunc'))
				}
				else mata: tdelinpred = -99
			
			/* Gauss-Kronrod */
				tempname knodes kweights
				mat `knodes' 	= J(1,15,.)
				mat `kweights' 	= J(1,15,.)

				local i=1
				foreach node of numlist 0.991455371120813 -0.991455371120813 0.949107912342759 -0.949107912342759 0.864864423359769 -0.864864423359769 0.741531185599394 -0.741531185599394 0.586087235467691 -0.586087235467691 0.405845151377397 -0.405845151377397 0.207784955007898 -0.207784955007898 0 {
					mat `knodes'[1,`i'] = `node'
					local `++i'
				}
				local i=1
				foreach weight of numlist 0.022935322010529 0.022935322010529 0.063092092629979 0.063092092629979 0.104790010322250 0.104790010322250 0.140653259715525 0.140653259715525 0.169004726639267 0.169004726639267 0.190350578064785 0.190350578064785 0.204432940075298 0.204432940075298 0.209482141084728  {
					mat `kweights'[1,`i'] = `weight'
					local `++i'
				}
						
			/* Pass to Mata GK nodes and weights */
				mata: nodes = J(`n',1,st_matrix("`knodes'"))
				mata: weights = J(`n',1,st_matrix("`kweights'"))
				
			/* Starting values - could improve */
				tempvar nr_time_old nr_time
				gen `nr_time_old' 	= 1 if `use'
				gen `nr_time' 		= `nr_time_old' if `use'
				mata: nr_time_old 	= st_data(.,"`nr_time_old'","`use'")
				mata: nr_time		= st_data(.,"`nr_time'","`use'")
		
			di "$newhaz"
			di "$tdenewfunc"
				survsim_nr

		}

		
		
		
		
		
		
		
		
		
	/********************************************************************************************************************************************************************************************/
	/* Stata implementation of N-R scheme */
	
		tempvar nr_time nr_time_old

		/* Mixture Weibull */
		if "`model'"=="mixture" & "`nr'"=="" {
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
		else if "`model'"=="cr" & "`nr'"==""{
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

		if "`nr'"=="" {	
			if "`gammas'"!="" {
				if `gammastart'==0 {
					local gammastart = 0.0001
				}
			}
			if `lambdastart'==0 {
				local lambdastart = 0.0001
			}
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Standard exp/Weibull/Gompertz calculations OR Starting values */
	
		if "`nr'"=="" {
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
		}
		
	/********************************************************************************************************************************************************************************************/
	/* Newton-Raphson for standard mixture and competing risks */

		if ("`model'"=="mixture" | "`model'"=="cr") & "`nr'"=="" {
			
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
	
		if "`nr'"=="" {
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

program define gaussquad_survsim, rclass
	syntax [, N(integer -99) LEGendre CHEB1 CHEB2 HERmite JACobi LAGuerre alpha(real 0) beta(real 0)]
	
    if `n' < 0 {
        display as err "need non-negative number of nodes"
		exit 198
	}
	if wordcount(`"`legendre' `cheb1' `cheb2' `hermite' `jacobi' `laguerre'"') > 1 {
		display as error "You have specified more than one integration option"
		exit 198
	}
	local inttype `legendre'`cheb1'`cheb2'`hermite'`jacobi'`laguerre' 
	if "`inttype'" == "" {
		display as error "You must specify one of the integration type options"
		exit 198
	}

	tempname weights nodes
	mata gq_survsim("`weights'","`nodes'")
	return matrix weights = `weights'
	return matrix nodes = `nodes'
end

mata:
	void gq_survsim(string scalar weightsname, string scalar nodesname)
{
	n =  strtoreal(st_local("n"))
	inttype = st_local("inttype")
	i = range(1,n,1)'
	i1 = range(1,n-1,1)'
	alpha = strtoreal(st_local("alpha"))
	beta = strtoreal(st_local("beta"))
		
	if(inttype == "legendre") {
		muzero = 2
		a = J(1,n,0)
		b = i1:/sqrt(4 :* i1:^2 :- 1)
	}
	else if(inttype == "cheb1") {
		muzero = pi()
		a = J(1,n,0)
		b = J(1,n-1,0.5)
		b[1] = sqrt(0.5)
    }
	else if(inttype == "cheb2") {
		muzero = pi()/2
		a = J(1,n,0)
		b = J(1,n-1,0.5)
	}
	else if(inttype == "hermite") {
		muzero = sqrt(pi())
		a = J(1,n,0)
		b = sqrt(i1:/2)
	}
	else if(inttype == "jacobi") {
		ab = alpha + beta
		muzero = 2:^(ab :+ 1) :* gamma(alpha + 1) * gamma(beta + 1):/gamma(ab :+ 2)
		a = i
		a[1] = (beta - alpha):/(ab :+ 2)
		i2 = range(2,n,1)'
		abi = ab :+ 2 :* i2
		a[i2] = (beta:^2 :- alpha^2):/(abi :- 2):/abi
		b = i1
        b[1] = sqrt(4 * (alpha + 1) * (beta + 1):/(ab :+ 2):^2:/(ab :+ 3))
        i2 = i1[2..n-1]
        abi = ab :+ 2 :* i2
        b[i2] = sqrt(4 :* i2 :* (i2 :+ alpha) :* (i2 :+ beta) :* (i2 :+ ab):/(abi:^2 :- 1):/abi:^2)
	}
	else if(inttype == "laguerre") {
		a = 2 :* i :- 1 :+ alpha
		b = sqrt(i1 :* (i1 :+ alpha))
		muzero = gamma(alpha :+ 1)
    }

	A= diag(a)
	for(j=1;j<=n-1;j++){
		A[j,j+1] = b[j]
		A[j+1,j] = b[j]
	}
	symeigensystem(A,vec,nodes)
	weights = (vec[1,]:^2:*muzero)'
	weights = weights[order(nodes',1)]
	nodes = nodes'[order(nodes',1)']
	st_matrix(weightsname,weights)
	st_matrix(nodesname,nodes)
}
		
end
