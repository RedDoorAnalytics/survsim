{smcl}
{* Michael Crowther 08aug2011 }{...}
{hline}
{cmd:help survsim}{right: }
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:survsim} {hline 2}}Simulate complex survival data{p_end}
{p2colreset}{...}


{title:Syntax}

{phang2}
{cmd: survsim} {it:varname1} [{it:varname2}] [{cmd:,} {it:options}]


{synoptset 34}{...}
{synopthdr}
{synoptline}
{synopt:{opt n(int)}} specifies the number of survival times to generate.{p_end}
{synopt:{opt lambdas(numlist)}} scale parameters for Weibull distributions.{p_end}
{synopt:{opt gammas(numlist)}} shape parameters for Weibull distributions.{p_end}
{synopt:{opt scov:ariates(string)}} baseline covariates and coefficient values to include in the survival model.{p_end}
{synopt:{opt centol(real)}} set the tolerance for the Newton-Raphson scheme. Default is 0.0001.{p_end}

{syntab:2-component mixture}
{synopt:{opt mixture}} simulate survival times from a 2-component mixture Weibull model.{p_end}
{synopt:{opt pm:ix(real)}} mixture parameter.{p_end}

{syntab:Competing risks}
{synopt:{opt cr}} simulate times from cause-specific hazards.{p_end}
{synopt:{opt ncr(int)}} specifies the number of competeing risks.{p_end}

{syntab:Joint model}
{synopt:{opt jm}} simulate survival times from a joint model of longitudinal and survival data.{p_end}
{synopt:{opt betas(numlist)}} specifies the means of the nultivariate normal distribution of the intercept and slope.{p_end}
{synopt:{opt sds(numlist)}} specifies the std. dev.'s of the nultivariate normal distribution of the intercept and slope.{p_end}
{synopt:{opt corr:elation(real)}} specifies the correlation between intercept and slope.{p_end}
{synopt:{opt alpha(real)}} specifies the current value association between longitudinal and survival submodels.{p_end}
{synopt:{opt lcov:ariates(string)}} baseline covariates and coefficient values to include in the longitudinal submodel.{p_end}

{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}{cmd:survsim} generates survival times from a variety of Weibull based distributions. Newton-Raphson iterations are used to generate survival times. {it:varname1} specifies the variable name
of which survival times are generated. {it:varname2} is required when generating competing risks data to create the status indicator.{p_end}

{title:Options}

{phang}{opt n} specifies the number of survival times to generate.{p_end}

{phang}{opt lambdas} defines the scale parameters in the Weibull distributions.The number of values needed depends on the model choice. Default is a single number. Under a 
{cmd:mixture} model 2 values are required, under a competing risks, {cmd:cr}, model the number of values are defined by {cmd:ncr}. A joint model {cmd:jm} requires a single value.{p_end}

{phang}{opt gammas} defines the shape parameters of the Weibull distributions. Number of entries must be equal to that of {cmd:lambdas}{p_end}

{phang}{opt scovariates} defines baseline covariates to be included in the linear predictor of the survival model, along with the value of the corresponding coefficient. For example
 a treatent variable coded 0/1 can be included, with a log hazard ratio of 0.5, by {cmd:scovarites}(treat 0.5). Variable treat must be in the dataset before {cmd:survsim} is run.{p_end}

{phang}{opt centol} specifies the tolerance of the Newton-Raphson scheme. Default is 0.0001.{p_end}

{dlgtab:Mixture model}

{phang}{opt mixture} specifies that survival times are simulated from a 2-component mixture Weibull distribution. {cmd:lambdas} and {cmd:gammas} must be of length 2.{p_end}

{phang}{opt pmix} defines the mixture parameter. Default is 0.5.{p_end}

{dlgtab:Competing risks}

{phang}{opt cr} specifies that survival times are simulated from the all-cause distribution from {cmd:ncr} cause-specific hazards.{p_end}

{phang}{opt ncr} defines the number of competing risks. {cmd:lambdas} and {cmd:gammas} must be of length {cmd:ncr}.{p_end}

{dlgtab:Joint model}

{title:Remarks}

{title:Examples}

{pstd}Generate times from a Weibull model{p_end}
{phang}{cmd:. survsim stime, n(1000) lambdas(0.1) gammas(1.5)}{p_end}

{pstd}Generate times from a 2-component mixture Weibull model{p_end}
{phang}{cmd:. survsim stime, n(1000) mixture lambdas(0.1 0.05) gammas(0.5 1.5) pmix(0.5)}{p_end}

{pstd}Generate times from a competing risks model with 4 cause-specific hazards{p_end}
{phang}{cmd:. survsim stime status, n(1000) cr ncr(4) lambdas(0.1 0.05 0.1 0.05) gammas(0.5 1.5 1 1.2)}{p_end}

{title:References}

{phang}Bender, R.; Augustin, T. and Blettner, M. Generating survival times to simulate Cox proportional hazards models. Stat Med, 2005, 24, 1713-1723{p_end}



{title:Authors}

{pstd}Michael J. Crowther, University of Leicester, United Kingdom. {browse "mailto:mjc76@le.ac.uk":mjc76@le.ac.uk}.{p_end}

{phang}Please report any errors you may find.{p_end}

