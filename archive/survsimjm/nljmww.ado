program define nljmww
	syntax varlist(min=1 max=1) [if], at(name) unif(real) lambda1(real) gamma1(real) alpha(real) beta0(real) beta1(real) scov(real) 
	
	tempname A B
	scalar `A' = `at'[1, 1]
	scalar `B' = `at'[1, 2]
	
	tempvar yh
	gen double 	`yh' = 		`unif' + exp((-1) * `lambda1' * `A'^(`gamma1') * exp(`alpha'*(`beta0' + `beta1' * `A') + `scov')) in 1
	replace 	`yh' = -1 + `unif' + exp((-1) * `lambda1' * `B'^(`gamma1') * exp(`alpha'*(`beta0' + `beta1' * `B') + `scov')) in 2
	replace `varlist' = `yh'
end


