program stww_d0
	version 11

	args todo b lnf

	tempvar p lambda1 gamma1 lambda2 gamma2 coef
	mleval `p' 	= `b', eq(1) scalar
	scalar `p' 	= invlogit(`p')
	
	mleval `lambda1' 	= `b', eq(2) scalar
	scalar `lambda1' 	= exp(`lambda1')	

	mleval `gamma1' 	= `b', eq(3) scalar
	scalar `gamma1' 	= exp(`gamma1')	
	
	mleval `lambda2' 	= `b', eq(4) scalar
	scalar `lambda2' 	= exp(`lambda2')	
	
	mleval `gamma2' 	= `b', eq(5) scalar
	scalar `gamma2' 	= exp(`gamma2')	
	
	mleval `coef' = `b', eq(6)
		
	local t "$ML_y1"
	local d "$ML_y2"

	tempvar surv haz
	
	qui {
		gen double `surv' = (`p'*exp(-`lambda1'*`t'^`gamma1') + (1-`p')*exp(-`lambda2'*`t'^`gamma2') )^(exp(`coef'))
		gen double `haz' = 	(exp(`coef')*(`lambda1'*`gamma1'*`p'*`t'^(`gamma1'-1)*exp(-`lambda1'*`t'^`gamma1') + `lambda2'*`gamma2'*(1-`p')*`t'^(`gamma2'-1)*exp(-`lambda2'*`t'^`gamma2') )/(`p'*exp(-`lambda1'*`t'^`gamma1') + (1-`p')*exp(-`lambda2'*`t'^`gamma2')))^`d'
		mlsum `lnf' = log(`surv'*`haz') 
		
	}
end

