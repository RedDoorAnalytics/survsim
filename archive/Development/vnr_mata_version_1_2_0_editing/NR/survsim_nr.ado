*! version 1.0.0 Feb2012

program define survsim_nr
 
	mata: nr("$newvarname",n,nr_time_old,nr_time,nodes,weights,u,centol,linpred,tdelinpred)
	
end


mata:
mata set matastrict off		
	void nr(	string scalar newvar,
				real scalar n,
				numeric matrix nr_time_old,
				numeric matrix nr_time,
				numeric matrix nodes, 
				numeric matrix weights, 
				numeric matrix u,
				real scalar centol,
				numeric matrix linpred,
				numeric matrix tdelinpred)
{	

	done = 0
	while (done==0) {
	
		/* Adjust nodes and weights */
			x = nodes:*(nr_time_old:/2) :+ (nr_time_old:/2)		//x needs to contain transformed nodes
			newweights = (nr_time_old:/2):*weights
		
		/* Survival function using quadrature and xb */
			haz_nodes = $newhaz :* exp(linpred)
			cumhaz = quadrowsum(newweights:*haz_nodes)		
			surv = exp((-1):*cumhaz)
			nr_xb = surv:-u	
		
		/* User specified hazard at current time values */
			x = nr_time_old										//x now needs to contain current time values
			haz = ($newhaz) :* exp(linpred)
			nr_dxb = (-1):*surv:*haz
		
		/* NR step */	
			nr_time = nr_time_old :- (nr_xb:/nr_dxb)
		
		/* Check */
			nr_time = (nr_time:<0):*0.00000000001 :+ (nr_time:>0):*nr_time
		
		/* Minimum difference */		
			error = abs(nr_time:-nr_time_old)
			max = max(error)
			if (max<centol) {
				done = 1
			}
			else nr_time_old = nr_time

	}
	
	st_view(final=.,.,newvar)
	final[.,.] = nr_time
	
}		
end
