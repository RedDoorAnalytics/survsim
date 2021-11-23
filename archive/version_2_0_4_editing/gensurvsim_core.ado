*! version 1.1.0 17oct2012 MJC

/*
History
MJC 17oct2012: version 1.1.0 - added globals to allow for extra variables within logh()/h() for tvc's
							 - synched for cumhazard() and logcumhazard()
MJC 27aug2012: version 1.0.0
*/

program define gensurvsim_core
	version 11.2
	syntax , maxtime(string) mintime(real) tempnewvarname(varname) tempu(varname) rc(varname) tempcovs(varname) centol(string) ///
			[temptde(varname) iterations(int 1000) enter(varname) CH]
	mata: gensurvsim_mata($overallsyntax2)
end

mata
function gensurvsim_func(t,U,nodes,weights,i,xb,enter $mmrootsyntax2)
{
	tnodes = (t:-enter):*0.5:*nodes :+ (t:+enter):/2
	tweights = (t:-enter):*weights:/2	
	cumhaz_nodes = $cumhaz
	lhs = exp(-cumhaz_nodes * tweights) - U
	return(lhs)
}

function gensurvsim_func_tde(t,U,nodes,weights,i,xb,enter,tdexb $mmrootsyntax2)
{
	tnodes = (t:-enter):*0.5:*nodes :+ (t:+enter):/2
	tweights = (t:-enter):*weights:/2
	cumhaz_nodes = $cumhaz
	lhs = exp(-cumhaz_nodes * tweights) - U
	return(lhs)
}

function gensurvsim_func_ch(t_cumhaz_temp,U,i,xb $mmrootsyntax2)
{
	return(exp(-($cumhaz)) :- U)
}

function gensurvsim_func_ch_tde(t_cumhaz_temp,U,i,xb,tdexb $mmrootsyntax2)
{
	return(exp(-($cumhaz)) :- U)
}

void gensurvsim_mata($overallsyntax1)
{
	N = st_nobs()
	st_view(time=.,.,st_local("tempnewvarname"))
	st_view(xb=.,.,st_local("tempcovs"))
	st_view(rc=.,.,st_local("rc"))
	st_view(U=.,.,st_local("tempu"))
	if (st_local("ch")=="") {
		nodes 	= st_matrix("r(nodes)")'
		weights = st_matrix("r(weights)")
	}
	maxt 	= strtoreal(st_local("maxtime"))
	mint 	= strtoreal(st_local("mintime"))
	maxit 	= strtoreal(st_local("iterations"))
	tol 	= strtoreal(st_local("centol"))
	
	//overide minimum limit of quadrature
	if (st_local("enter")!="") {
		st_view(enter=.,.,st_local("enter"))
	}
	else enter = J(N,1,0)
	
	if (st_local("ch")!="") {
		if (st_local("temptde")=="") {
			for (i=1;i<=N;i++){
				rc[i] = mm_root(t_cumhaz_temp=1,&gensurvsim_func_ch(),mint,maxt,tol,maxit,U[i],i,xb[i] $mmrootsyntax1) 	//lower,upper,tol,maxit
				time[i] = t_cumhaz_temp
			}
		}
		else {
			st_view(tde=.,.,st_local("temptde"))
			for (i=1;i<=N;i++){
				rc[i] = mm_root(t_cumhaz_temp=1,&gensurvsim_func_ch_tde(),mint,maxt,tol,maxit,U[i],i,xb[i],tde[i] $mmrootsyntax1) 	//lower,upper,tol,maxit
				time[i] = t_cumhaz_temp
			}		
		}
	}
	else {
		if (st_local("temptde")=="") {
			for (i=1;i<=N;i++){
				rc[i] = mm_root(t=1,&gensurvsim_func(),mint,maxt,tol,maxit,U[i],nodes,weights,i,xb[i],enter[i] $mmrootsyntax1) 	//lower,upper,tol,maxit
				time[i] = t
			}
		}
		else {
			st_view(tde=.,.,st_local("temptde"))
			for (i=1;i<=N;i++){
				rc[i] = mm_root(t=1,&gensurvsim_func_tde(),mint,maxt,tol,maxit,U[i],nodes,weights,i,xb[i],enter[i],tde[i]$mmrootsyntax1)
				time[i] = t
			}
		}
	}
}
end
