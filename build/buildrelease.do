//build new version of survsim
//either for website or SSC
// --> run whole do file

//!!
// No compiled code, so can be built in any Stata >=14.2
//!!


//local drive Z:/
local drive /Users/Michael/Documents
cd `drive'/survsim/

local sscbuild = 1

//=======================================================================================================================//

if `sscbuild' {											
	
	//build for SSC -> current version up is 4.0.3
	
	local sscversion 4_0_6
	cap mkdir ./ssc/version_`sscversion'
	local fdir /Users/Michael/Documents/survsim/ssc/version_`sscversion'/
}
else {													
	
	//build for website -> current version up is 4.0.5
	
	local fdir /Users/Michael/Documents/website/static/code/survsim/
}

//=======================================================================================================================//

//pkg files
if `sscbuild' {
	copy ./build/survsim_details.txt `fdir', replace
}
else {
	copy ./build/survsim.pkg `fdir', replace
	copy ./build/stata.toc `fdir', replace
}
	
//=======================================================================================================================//

//survsim

	copy ./survsim/survsim.ado `fdir', replace
	copy ./survsim/survsim_user.ado `fdir', replace
	copy ./survsim/survsim_user_core.ado `fdir', replace
	copy ./survsim/survsim_model.ado `fdir', replace
	copy ./survsim/survsim_msm.ado `fdir', replace
	
	//help files
	copy ./survsim/survsim.sthlp `fdir', replace
	copy ./survsim/survsim_parametric.sthlp `fdir', replace
	copy ./survsim/survsim_user.sthlp `fdir', replace
	copy ./survsim/survsim_model.sthlp `fdir', replace
	copy ./survsim/survsim_msm.sthlp `fdir', replace
	
