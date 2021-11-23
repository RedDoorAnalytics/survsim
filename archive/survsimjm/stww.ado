program stww
	version 11.0
	
 	Estimate `0'
	
end

program Estimate, eclass sortpreserve
		syntax [varlist(default=empty)] [if] [in] [,									///
										*									///
										]
	
				mlopts mlopts , `options'									/* Extract any ml options to pass to ml model */
				marksample touse
				
				/* constraints for covariates */
				cap constraint drop _all
				constraint 1 [xb][_cons] = 0
				local constopts "constraints(1)"
				
				ml model d0 stww_d0											///
							/logit_p_mix			 						///
							(ln_lambda1: _t _d=)							/// 
							/ln_gamma1 										///
							/ln_lambda2										///
							/ln_gamma2 										///
							(xb: = `varlist')								///
							if `touse',										///
							`constopts'										///
							`mlopts'										///
							maximize
	
				ml display
		
end








