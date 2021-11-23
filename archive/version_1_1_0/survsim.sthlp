{smcl}
{* *! version 1.1.0 10sep2011 }{...}
{hline}
{cmd:help survsim}{right: }
{hline}

{title:Title}

{p2colset 5 16 20 2}{...}
{p2col :{cmd:survsim} {hline 2}}Simulate survival data{p_end}
{p2colreset}{...}


{title:Syntax}

{phang2}
{cmd: survsim} {it:newvarname1} [{it:newvarname2}] [{cmd:,} {it:options}]


{synoptset 34 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt n(int)}}specifies the number of survival times to generate, default is _N{p_end}
{synopt:{opt l:ambdas(numlist)}}scale parameters{p_end}
{synopt:{opt g:ammas(numlist)}}shape parameters{p_end}
{synopt:{cmdab:d:istribution(}{cmdab:gom:pertz)}}Gompertz survival distribution{p_end}
{synopt:{cmdab:d:istribution(}{cmdab:w:eibull)}}Weibull survival distribution(default){p_end}
{synopt:{opt cov:ariates(vn # [# ...] ...)}}baseline covariates{p_end}
{synopt:{opt tde(vn # [# ...] ...)}}time-dependent effects{p_end}

{syntab:2-component mixture}
{synopt:{opt mix:ture}}simulate survival times from a 2-component mixture model{p_end}
{synopt:{opt pm:ix(real)}}mixture parameter, default is 0.5{p_end}

{syntab:Competing risks}
{synopt:{opt cr}}simulate survival times from the all-cause distribution of cause-specific hazards{p_end}
{synopt:{opt ncr(int)}}specifies the number of competing risks{p_end}

{syntab:Newton-Raphson scheme}
{synopt:{opt centol(real)}}set the tolerance for the Newton-Raphson scheme, default is 0.0001{p_end}
{synopt:{opt show:diff}}display the maximum difference in estimates between iterations{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}{cmd:survsim} simulates survival times from Weibull or Gompertz parametric distributions. Newton-Raphson iterations are used to generate 
survival times under a mixture or competing risks model. Non-proportional hazards can be included with all models except a mixture model; under 
a Weibull model covariates are interacted with log time, under a Gompertz model covariates are interacted with time. {it:newvarname1} 
specifies the new variable name to contain the generated survival times. {it:newvarname2} is required when generating competing risks data to 
create the status indicator.{p_end}

{title:Options}

{phang}{opt n} specifies the number of survival times to generate. Default is _N.{p_end}

{phang}{opt lambdas}({it:numlist}) defines the scale parameters in the Weibull/Gompertz distribution(s). The number of values required depends on the model choice. 
Default is a single number corresponding to a standard parametric distribution. Under a {cmd:mixture} model 2 values are required. Under a 
competing risks model, {cmd:cr}, the number of values are defined by {cmd:ncr}.{p_end}

{phang}{opt gammas}({it:numlist}) defines the shape parameters of the parametric distribution(s). Number of entries must be equal to that of {cmd:lambdas}.{p_end}

{phang}{opt distribution}({it:string}) specifies the parametric survival distrubution to use. {cmd:gompertz} or {cmd:weibull} can be used, with {cmd:weibull}
the default.{p_end}

{phang}{opt covariates}({it:vn # [# ...] ...}) defines baseline covariates to be included in the linear predictor of the survival model, along with the value of the 
corresponding coefficient. For example a treatent variable coded 0/1 can be included, with a log hazard ratio of 0.5, by {cmd:covariates}(treat 0.5). 
Variable treat must be in the dataset before {cmd:survsim} is run. If {cmd:cr} is used with {cmd:ncr(4)}, then a value of each covariate must be 
inputted for each competing risk, e.g. {cmd:covariates}(treat 0.5 -0.2 0.1 0.25).{p_end}

{phang}{opt tde}({it:vn # [# ...] ...}) creates non-proportional hazards by interacting covariates with either log time or time under a Weibull or Gompertz model, respectively. 
This option is not available under mixture models. Values should be entered as {cmd:tde}(trt 0.5), for example.{p_end}
 
{dlgtab:Mixture model}

{phang}{opt mixture} specifies that survival times are simulated from a 2-component mixture model. {cmd:lambdas} and {cmd:gammas} must be of length 2.{p_end}

{phang}{opt pmix}({it:real}) defines the mixture parameter. Default is 0.5.{p_end}

{dlgtab:Competing risks}

{phang}{opt cr} specifies that survival times are simulated from the all-cause distribution from {cmd:ncr} cause-specific hazards.{p_end}

{phang}{opt ncr}({it:int}) defines the number of competing risks. {cmd:lambdas} and {cmd:gammas} must be of length {cmd:ncr}.{p_end}

{dlgtab:Newton-Raphson scheme}

{phang}{opt centol}({it:real}) specifies the tolerance of the Newton-Raphson scheme. Default is 0.0001.{p_end}

{phang}{opt showdiff} display the maximum difference in estimates between iterations. This can be used to monitor convergence.{p_end}

{title:Remarks}

{pstd}On rare occasions the Newton-Raphson scheme may not converge. The user should experiment with appropriate parameter values and tolerance 
levels.{p_end}

{title:Examples}

{pstd}Generate times from a Gompertz model:{p_end}
{phang}{cmd:. survsim stime, n(1000) lambdas(0.1) gammas(0.05) dist(gompertz)}{p_end}

{pstd}Generate times from a 2-component mixture Weibull model:{p_end}
{phang}{cmd:. survsim stime, n(1000) mixture lambdas(0.1 0.05) gammas(1 1.5) pmix(0.5)}{p_end}

{pstd}Generate times from a competing risks model with 4 cause-specific Weibull hazards:{p_end}
{phang}{cmd:. survsim stime status, n(1000) cr ncr(4) lambdas(0.1 0.05 0.1 0.05) gammas(0.5 1.5 1 1.2)}{p_end}

{pstd}Generate times from a Weibull model with diminishing treatment effect:{p_end}
{phang}{cmd:. survsim stime, n(1000) lambdas(0.1) gammas(1.5) cov(trt -0.5) tde(trt 0.05)}{p_end}

{pstd}Generate times from a Weibull model including a binary treatment variable, with log(hazard ratio) = -0.5, and censoring after 5 years:{p_end}
{phang}{cmd:. clear}{p_end}
{phang}{cmd:. set obs 1000}{p_end}
{phang}{cmd:. gen trt = rbinomial(1,0.5)}{p_end}
{phang}{cmd:. survsim stime, lambdas(0.1) gammas(1.5) cov(trt -0.5)}{p_end}
{phang}{cmd:. gen died = stime<=5}{p_end}
{phang}{cmd:. replace stime = 5 if died==0}{p_end}
{phang}{cmd:. stset stime, f(died=1)}{p_end}
{phang}{cmd:. streg trt, dist(weibull) nohr}{p_end}

{title:References}

{phang}Bender, R.; Augustin, T. and Blettner, M. Generating survival times to simulate Cox proportional hazards models. Stat Med, 2005, 24, 1713-1723{p_end}

{phang}Beyersmann, J.; Latouche, A.; Buchholz, A. & Schumacher, M. Simulating competing risks data in survival analysis. Stat Med, 2009, 28, 956-971{p_end}

{title:Authors}

{pstd}Michael J. Crowther, University of Leicester, United Kingdom. {browse "mailto:mjc76@le.ac.uk":mjc76@le.ac.uk}.{p_end}

{phang}Please report any errors you may find.{p_end}

