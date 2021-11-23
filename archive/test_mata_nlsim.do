clear
pr drop _all
mata
void nltest(numeric matrix knodes, numeric matrix kweights, real scalar u) 
{
222
	st_view(x=.,.,st_local("yh"))
	A = strtoreal(st_local("A"))
	B = strtoreal(st_local("B"))	
	l = strtoreal(st_local("lambda"))
	g = strtoreal(st_local("gamma"))
	x[1,] = u + rowsum((A:/2):*knodes :* ( l:*g:*(A:/2:*knodes :+ A:/2):^(g:-1) ))
	x[2,] = -1 :+ u :+ rowsum((B:/2):*knodes :* ( l:*g:*(B:/2:*knodes :+ B:/2):^(g:-1) ))	
	x
}
end

program define nlww
	syntax varlist(min=1 max=1) [if], at(name) unif(real) lambda(real) gamma(real)
	
	tempname A B
	scalar `A' = `at'[1, 1]
	scalar `B' = `at'[1, 2]
	n di "test"
	tempvar yh
	 
	mata: nltest(st_matrix("knodes"),st_matrix("kweights"),`unif')
	
	replace `varlist' = `yh'
end





local obs = 250
local lambda = 0.1
local gamma = 1.5

set obs `obs'
local count = 0
gen id_var = .
gen stime = .

mat knodes      = J(1,15,.)
mat kweights    = J(1,15,.)

local i=1
foreach n of numlist 0.991455371120813 -0.991455371120813 0.949107912342759 -0.949107912342759 0.864864423359769 -0.864864423359769 0.741531185599394 -0.741531185599394 0.586087235467691 -0.586087235467691 0.405845151377397 -0.405845151377397 0.207784955007898 -0.207784955007898 0 {
		mat knodes[1,`i'] = `n'
		local `++i'
}
local i=1
foreach n of numlist 0.022935322010529 0.022935322010529 0.063092092629979 0.063092092629979 0.104790010322250 0.104790010322250 0.140653259715525 0.140653259715525 0.169004726639267 0.169004726639267 0.190350578064785 0.190350578064785 0.204432940075298 0.204432940075298 0.209482141084728  {
		mat kweights[1,`i'] = `n'
		local `++i'
}



forvalues i=1/`obs' {
	generate y=0
	qui replace y=1 in 1

	local u = runiform()
	
	nl ww @ y in 1/2, parameters(A B) initial(A 1 B 1) 				///
			unif(`u') lambda(`lambda') gamma(`gamma') 

	while (_b[_cons]<0) {
		local u = runiform()
		local count = `count' +1
		nl ww @ y in 1/2, parameters(A B) initial(A 1 B 1) 			///
				unif(`u') lambda(`lambda') gamma(`gamma') 
	}
	qui replace id_var = `i' in `i'
	qui replace stime = _b[_cons] in `i'
	drop y
}
gen died = stime<5
replace stime = 5 if died==0
stset stime, f(died)
streg , dist(w) nohr
