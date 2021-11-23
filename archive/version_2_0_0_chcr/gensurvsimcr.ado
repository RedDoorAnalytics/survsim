*! version 1.0.0 ?????2012 MJC

/*
History
MJC ?????2012: version 1.0.0
*/

program define gensurvsimcr
	version 11.2
	
	syntax newvarname(min=1 max=2), 										///
																			///
									[										/// -Options-
																			///
										MAXTime(string)						///	-Maximum simulated time-
										MINTime(real 1E-05)					///	-minimum time for use in user-defined generation-
										NODES(int 15)						///	-Quadrature nodes-
																			///
										COVariates(string)					///	-Baseline covariates, e.g, (sex 0.5 race -0.4)-
										TDE(string)							///	-Time dependent effects to interact with tdefunc()-
										TDEFUNCtion(string)					///	-function of time to interact with time-dependent effects-
																			///
										ITerations(int 1000)				///
										CENTOL(real 0.00001)				///
										ENTER(varname)						///
										*									///
									]
		
	//=================================================================================================================================================//
	//basic error checks
		
		local nvars : word count `varlist'
		local nvar1 : word 1 of `varlist'
		if `nvars'==2 local nvar2 : word 2 of `varlist'
		
		cap which lmoremata.mlib
		if _rc {
			display in yellow "You need to install the moremata package. This can be installed using,"
			display in yellow ". {stata ssc install moremata}"
			exit 198
		}
		
		if "`maxtime'"=="" {
			di as error "maxtime() must be specified"
			exit 198
		}
		
		cap confirm number `maxtime'
		if _rc {
			di as error "maxtime() must be a number >0"
			exit 198
		}
		if `maxtime'<=0 {
			di as error "maxtime() must be a number >0"
			exit 198
		}
				
		if "`enter'"!="" {
			local enter enter(`enter')
		}
				
		qui ds							//keep here
		local allvars `r(varlist)'

	//=================================================================================================================================================//
	//Competing risks syntax
	
		local i = 1
		local 0 ,`options'
		syntax ,[						///
					H`i'(string) 		///
					LOGH`i'(string)		///
					CUMH`i'(string)		///
					LOGCUMH`i'(string)	///
					*					///
				]
				
		while wordcount("`h`i'' `logh`i'' `cumh`i'' `logcumh`i''")>0 {
			local `++i'
			local 0 ,`options'
			syntax ,[						///
						H`i'(string) 		///
						LOGH`i'(string)		///
						CUMH`i'(string)		///
						LOGCUMH`i'(string)	///
						*					///
					]
			
		}
		
		local Ncr = `i' - 1
		di "Ncr: `Ncr'"

	//===============================================================================================================================================================//
	// Baseline covariates
		
		if "`covariates'"!="" {
			tokenize `covariates'
			local ncovlist : word count `covariates'	
			local ncovvars = `ncovlist'/`=`Ncr'+1'
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
				forvalues j=1/`Ncr' {
					cap confirm num ``=`ind'+`j'''
					if _rc {
						local errortxt "invalid covariates(... ``ind'' ``=`ind'+`j''' ...)"
						local error = 1
					}
					/* Create effect for ith variable and jth risk */
					tempvar vareffect_`i'_`j'
					gen double `vareffect_`i'_`j'' = ``ind''*``=`ind'+`j''' 
				}
				local ind = `ind' + `Ncr' + 1
			}
			if "`error'"=="1" {
				di as error "`errortxt'"
				exit 198
			}
			forvalues k=1/`Ncr' {
				local cov_linpred_`k' "`vareffect_1_`k''"
			}
			if `ncovvars'>1 {
				forvalues p=2/`ncovvars' {
					forvalues m=1/`Ncr' {
						local cov_linpred_`m' "`cov_linpred_`m'' + `vareffect_`p'_`m''"
					}
				}
			}
			forvalues k=1/`Ncr' {
				tempvar xbcovs`k'
				qui gen double `xbcovs`k'' = exp(`cov_linpred_`k'')
				local xbvars `xbvars' `xbcovs`k''
			}			
			
		}
		else {
			forvalues k=1/`Ncr' {
				tempvar xbcovs`k'
				qui gen byte `xbcovs`k'' = 1
				local xbvars `xbvars' `xbcovs`k''
			}			
		}
	
	//===============================================================================================================================================================//
	// Time-dependent effects
	
		if "`tde'"!="" {
			tokenize `tde'
			local ntdelist : word count `tde'	
			local ntdevars = `ntdelist'/`=`Ncr'+1'
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
				forvalues j=1/`Ncr' {
					cap confirm num ``=`ind'+`j'''
					if _rc {
						local errortxt "invalid tde(... ``ind'' ``=`ind'+`j''' ...)"
						local error = 1
					}
					/* Create effect for ith variable and jth risk */
					tempvar tdeeffect_`i'_`j'
					gen double `tdeeffect_`i'_`j'' = ``ind''*``=`ind'+`j''' 
				}
				local ind = `ind' + `Ncr' + 1
			}
			if "`error'"=="1" {
				di as error "`errortxt'"
				exit 198
			}
			forvalues k=1/`Ncr' {
				local tde_linpred_`k' "`tdeeffect_1_`k''"
			}
			if `ntdevars'>1 {
				forvalues p=2/`ntdevars' {
					forvalues m=1/`Ncr' {
						local tde_linpred_`m' "`tde_linpred_`m'' + `tdeeffect_`p'_`m''"
					}
				}
			}

			forvalues k=1/`Ncr' {
				tempvar tdexb`k'
				qui gen double `tdexb`k'' = `tde_linpred_`k''
				local tdexbvars `tdexbvars' `tdexb`k''
			}
			local tdexbvars tdevars(`tdexbvars')
			local tdexbvars2 `tdexbvars'	//needed for event indicator
		}
		
	//===============================================================================================================================================================//
	// Core stuff

		qui gen double `nvar1' = .
		qui gen double `nvar2' = .
		tempvar tempu
		qui gen double `tempu' = runiform()
		cap drop _survsim_rc
		qui gen _survsim_rc = .
		
		// hazard or log hazard
		if "`h1'"!="" | "`logh1'"!="" {
		
			//nodes and weights
			gaussquad_ss, n(`nodes') leg
								
			// handle loghazard() or hazard()
			mata: st_local("tempcumhaz1",subinstr("`logh1'`h1'","#t","tnodes"))
			if "`logh1'"!="" local tempcumhaz1 exp(`tempcumhaz1')
			local tempcumhaz1 (`tempcumhaz1'):*xb[,1]
			forvalues i=2/`Ncr' {
				mata: st_local("tempcumhaz`i'",subinstr("`logh`i''`h`i''","#t","tnodes"))
				if "`logh`i''"!="" local tempcumhaz`i' exp(`tempcumhaz`i'')
				local tempcumhaz`i' (`tempcumhaz`i''):*xb[,`i']
			}
			
			//tde's
			if "`tde'"!="" {
				if "`tdefunction'"=="" local tdefunction tnodes
				else mata: st_local("tdefunction",subinstr(st_local("tdefunction"),"#t","tnodes"))
				local tempcumhaz1 (`tempcumhaz1') :* exp(tdexb[,1] :* (`tdefunction'))
				forvalues i=2/`Ncr' {
					local tempcumhaz`i' (`tempcumhaz`i'') :* exp(tdexb[,`i'] :* (`tdefunction'))
				}	
			}
			
			//chuck variables found in loghazard()/hazard() into Mata		//!! redo into matrix
			local nhazvars = 0
			macro drop overallsyntax1 overallsyntax2 mmrootsyntax1 mmrootsyntax2
			foreach varname in `allvars' {
				forvalues i=1/`Ncr' {
					mata: st_local("test",strofreal(regexm("`h`i''`logh`i''","`varname'"))) 
					if `test' {
						mata: st_view(`varname'=.,.,"`varname'")
						local nhazvars = `nhazvars' + 1
						if `nhazvars'==1 {
							global overallsyntax1 numeric matrix `varname'
							global overallsyntax2 `varname'
						}
						else {
							global overallsyntax1 $overallsyntax1 , numeric matrix `varname'
							global overallsyntax2 $overallsyntax2 ,`varname'
						}
						global mmrootsyntax1 $mmrootsyntax1 , `varname'[i]
						global mmrootsyntax2 $mmrootsyntax2 , `varname'
					}
				}
			}
			
			//overall hazard function
			local allhaz `tempcumhaz1' 
			forvalues i=2/`Ncr' {
				local allhaz `allhaz' :+ `tempcumhaz`i''
			}
			global cumhaz `allhaz'	
			global cumhaz0 0 		//ignore
						
			cap pr drop gensurvsimcr_core
			gensurvsimcr_core, 	maxtime(`maxtime') 				///
								tempnewvarname(`nvar1') 		///
								tempu(`tempu') 					///
								rc(_survsim_rc) 				///
								tempcovs(`xbvars')	 			///
								`tdexbvars'						///
								iterations(`iterations') 		///
								centol(`centol') 				///
								mintime(`mintime') 				///
								`enter' `mixture'
								
			cap macro drop cumhaz cumhaz0
			if `nhazvars'!=0 {
				cap macro drop overallsyntax1 overallsyntax2 mmrootsyntax1 mmrootsyntax2
			}
		}
				
		if ("`cumh1'"!="" | "`logcumh1'"!="") {

			// handle logcumh() or cumh()
			mata: st_local("tempcumhaz1",subinstr("`logcumh1'`cumh1'","#t","t"))
			mata: st_local("tempcumhaz1_0",subinstr("`logcumh1'`cumh1'","#t","enter"))
			if "`logcumh1'"!="" {
				local tempcumhaz1 exp(`tempcumhaz1')
				local tempcumhaz1_0 exp(`tempcumhaz1_0')
			}
			local tempcumhaz1 (`tempcumhaz1'):*xb[,1]
			local tempcumhaz1_0 (`tempcumhaz1_0'):*xb[,1]
			forvalues i=2/`Ncr' {
				mata: st_local("tempcumhaz`i'",subinstr("`logcumh`i''`cumh`i''","#t","t"))
				mata: st_local("tempcumhaz`i'_0",subinstr("`logcumh`i''`cumh`i''","#t","enter"))
				if "`logcumh`i''"!="" {
					local tempcumhaz`i' exp(`tempcumhaz`i'')
					local tempcumhaz`i'_0 exp(`tempcumhaz`i'_0')
				}
				local tempcumhaz`i' (`tempcumhaz`i''):*xb[,`i']
				local tempcumhaz`i'_0 (`tempcumhaz`i'_0'):*xb[,`i']
			}
			
			//tde's
			if "`tde'"!="" {
				if "`tdefunction'"=="" {
					local tdefunction1 t
					local tdefunction0 enter
				}
				else {
					mata: st_local("tdefunction1",subinstr(st_local("tdefunction"),"#t","t"))
					mata: st_local("tdefunction0",subinstr(st_local("tdefunction"),"#t","enter"))
				}
				local tempcumhaz1 (`tempcumhaz1') :* exp(tdexb[,1] :* (`tdefunction1'))
				local tempcumhaz1_0 (`tempcumhaz1_0') :* exp(tdexb[,1] :* (`tdefunction0'))
				forvalues i=2/`Ncr' {
					local tempcumhaz`i' (`tempcumhaz`i'') :* exp(tdexb[,`i'] :* (`tdefunction1'))
					local tempcumhaz`i'_0 (`tempcumhaz`i'_0') :* exp(tdexb[,`i'] :* (`tdefunction0'))
				}	
			}
			
			//chuck variables found in loghazard()/hazard() into Mata		//!! redo into matrix
			local nhazvars = 0
			macro drop overallsyntax1 overallsyntax2 mmrootsyntax1 mmrootsyntax2
			foreach varname in `allvars' {
				forvalues i=1/`Ncr' {
					mata: st_local("test",strofreal(regexm("`cumh`i''`logcumh`i''","`varname'"))) 
					if `test' {
						mata: st_view(`varname'=.,.,"`varname'")
						local nhazvars = `nhazvars' + 1
						if `nhazvars'==1 {
							global overallsyntax1 numeric matrix `varname'
							global overallsyntax2 `varname'
						}
						else {
							global overallsyntax1 $overallsyntax1 , numeric matrix `varname'
							global overallsyntax2 $overallsyntax2 ,`varname'
						}
						global mmrootsyntax1 $mmrootsyntax1 , `varname'[i]
						global mmrootsyntax2 $mmrootsyntax2 , `varname'
					}
				}
			}
			
			//overall hazard function
			local allhaz `tempcumhaz1' 
			local allhaz0 `tempcumhaz1_0' 
			forvalues i=2/`Ncr' {
				local allhaz `allhaz' :+ `tempcumhaz`i''
				local allhaz0 `allhaz0' :+ `tempcumhaz`i'_0'
			}
			global cumhaz `allhaz'	
			global cumhaz0 `allhaz0'

			cap pr drop gensurvsimcr_core
			gensurvsimcr_core, 	maxtime(`maxtime') 				///
								tempnewvarname(`nvar1') 		///
								tempu(`tempu') 					///
								rc(_survsim_rc) 				///
								tempcovs(`xbvars')	 			///
								`tdexbvars'						///
								iterations(`iterations') 		///
								centol(`centol') 				///
								mintime(`mintime') 				///
								`enter' `mixture' ch
								
			cap macro drop cumhaz cumhaz0
			cap macro drop overallsyntax1 overallsyntax2 mmrootsyntax1 mmrootsyntax2

		}
		
		
	//===============================================================================================================================================================//
	//summarise _survsim_rc
	
		qui su _survsim_rc if _survsim_rc==1, meanonly
		if r(N)>0 {
			di in yellow "Warning: `r(N)' survival times did not converge"
			di in yellow "         They have been set to final iteration"
			di in yellow "         You can identify them by _survsim_rc = 1"
		}
		qui su _survsim_rc if _survsim_rc==2, meanonly
		if r(N)>0 {
			di in yellow "Warning: `r(N)' survival times were below the lower limit of `mintime'"
			di in yellow "         They have been set to `mintime'"
			di in yellow "         You can identify them by _survsim_rc = 2"
		}
		qui su _survsim_rc if _survsim_rc==3, meanonly
		if r(N)>0 {
			di in yellow "Warning: `r(N)' survival times were above the upper limit of `maxtime'"
			di in yellow "         They have been set to `maxtime' and can be considered censored"
			di in yellow "         You can identify them by _survsim_rc = 3"
		}
		
		//event indicator
		//need total hazard and cause-specific hazard
		mata: st_view(tnodes=.,.,st_local("nvar1"))
		mata: st_view(xb=.,.,tokens(st_local("xbvars")))
		if "`tde'"!="" {
			mata: st_view(tdexb=.,.,tokens(st_local("tdexbvars2")))
		}
		mata: tothaz = `allhaz'
		
		forvalues i=1/`Ncr' {
			tempvar p`i'
			qui gen double `p`i'' = .
			mata: st_store(.,st_local("p`i'"),(`tempcumhaz`i''):/tothaz)
			local pvars `pvars' `p`i''
		}
		mata: gencrstatus()
		qui replace `nvar2' = 0 if `nvar1'==`maxtime'

		
		
		
end

mata:
mata set matastrict off
void gencrstatus()
{
	st_view(final=.,.,st_local("nvar2"))
	st_view(pmat=.,.,tokens(st_local("pvars")))
	N = st_nobs()
	for (i=1; i<=N; i++) {
		final[i,] = rdiscrete(1,1,pmat[i,])
	}
}
end

program define gaussquad_ss, rclass
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
	mata ss_gq("`weights'","`nodes'")
	return matrix weights = `weights'
	return matrix nodes = `nodes'
end

mata:
	void ss_gq(string scalar weightsname, string scalar nodesname)
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
