{smcl}
{* *! version 1.0 14 Nov 2020}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install command2" "ssc install command2"}{...}
{vieweralsosee "Help command2 (if installed)" "help command2"}{...}
{viewerjumpto "Syntax" "spacereg##syntax"}{...}
{viewerjumpto "Description" "spacereg##description"}{...}
{viewerjumpto "Options" "spacereg##options"}{...}
{viewerjumpto "Remarks" "spacereg##remarks"}{...}
{viewerjumpto "Examples" "spacereg##examples"}{...}
{title:Title}
{phang}
{bf:spacereg} {hline 2} a command to compute spatial standard errors for M-estimators

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:spacereg}
depvar regressorlist
[{help if}]
[{help in}]
[{cmd:,}
{it:options}]

{pstd}
{it:Recommended syntax:}

{p 8 17 2}
{cmdab:spacereg} {it:depvar} {it:regressorlist}{cmd:,} {opt coords(varlist)} {opt cutoffs(numlist)} {opt model(str)} [{it:other_options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt coords(varlist)}} coordinate variables used to compute spatial distance {p_end}
{synopt:{opt cutoffs(numlist)}} cutoff(s) for the spatial kernel; either one number (replicated across coordinates) or one per coordinate {p_end}
{synopt:{opt mod:el(str)}} model to estimate; one of: ols, logit, probit, poisson, nb, areg, reghdfe {p_end}
{synopt:{opt ker:nel(str)}} spatial kernel; {opt bartlett} (default) or {opt uniform} {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}
{pstd}

{pstd}
 {cmd:spacereg} calculates spatial standard errors for several commonly
 used non-linear M-estimators. 

{pstd}
In the recommended syntax, spatial correction is defined by:
{p_end}
{phang}
{cmd:coords()} the coordinate variables (e.g., latitude/longitude)
{p_end}
{phang}
{cmd:cutoffs()} the Bartlett window cutoffs as numeric values (either one value, replicated across dimensions, or one value per coordinate)
{p_end}

{pstd}
{cmd:depvar} is the dependent variable. May be continuous, binary, or count.

{pstd}
{cmd:regressorlist} is the regressor list. Must include a constant unless 
{cmd:areg} or {cmd:reghdfe} are used. Fixed effects are included within {cmd:model}. 
See examples for reference.

{pstd}
Important: for {cmd:ols}, {cmd:logit}, {cmd:probit}, {cmd:poisson}, and {cmd:nb}, the underlying estimator is run with {cmd:noconstant}. If you want an intercept, include a constant regressor (e.g., {cmd:gen byte const=1} and include {cmd:const} in {it:regressorlist}).

{marker options}{...}
{title:Options}
{dlgtab:Main}
{phang}
{opt coords(varlist)} coordinate variables used to compute economic distance

{pstd}
{p_end}
{phang}
{opt cutoffs(numlist)} cutoff(s) for the Bartlett window. Supply either one number or one per coordinate.

{pstd}
{p_end}
{phang}
{opt mod:el(str)}     takes 1 of 7 acceptable models:
(1) ols (2) logit (3) probit (4) poisson (5) nb (6) areg (7) reghdfe
{p_end}
{phang}
{opt ker:nel(str)} kernel/weighting scheme for spatial dependence. Supported:
{p_end}
{phang2}{opt bartlett} (default) Bartlett weights within cutoffs
{p_end}
{phang2}{opt uniform} uniform weights within cutoffs (rectangle kernel)
{p_end}
{phang}
Other options are currently ignored.
{p_end}


{marker examples}{...}
{title:Examples}
{pstd}

{pstd}
Recommended syntax (numeric cutoffs):

{pstd}
{stata spacereg dep indep1 const, coords(C1 C2) cutoffs(100 100) model(ols)}

{pstd}
Fixed effects (model option contains the FE spec):

{pstd}
{stata spacereg dep indep1, coords(C1 C2) cutoffs(100 100) model(reghdfe, fe1 fe2)}

{pstd}
Uniform (rectangle) kernel:

{pstd}
{stata spacereg dep indep1 const, coords(C1 C2) cutoffs(100 100) model(ols) kernel(uniform)}

{pstd}

{title:Stored results}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}} coefficient row vector (posted results) {p_end}
{synopt:{cmd:e(V)}} Conley variance-covariance matrix (posted results) {p_end}
{synopt:{cmd:e(V_uncorrected)}} uncorrected model variance-covariance matrix {p_end}
{synopt:{cmd:e(se)}} Conley standard errors (vector) {p_end}
{synopt:{cmd:e(hessian)}} Hessian matrix used in the sandwich {p_end}
{synopt:{cmd:e(opg)}} Outer product of gradients ("meat") {p_end}

{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(vce)}} {cmd:conley} {p_end}
{synopt:{cmd:e(vcetype)}} {cmd:Conley} {p_end}
{synopt:{cmd:e(kernel)}} kernel used ({cmd:bartlett} or {cmd:uniform}) {p_end}
{synopt:{cmd:e(conley_model)}} model string passed to {cmd:model()} {p_end}


{title:References}
{pstd}

{pstd}
Jenish, N., & Prucha, I. R. (2009). Central limit theorems and uniform laws of 
large numbers for arrays of random fields. Journal of econometrics, 150(1), 86-98.

{pstd}
Conley, T. G. (1999). GMM estimation with cross sectional dependence. 
Journal of econometrics, 92(1), 1-45.


{title:Author}
{p}

Luis Calderon and Leander Heldring, University of Bonn and Briq.

Email {browse "mailto:luis.calderon@uni-bonn.de":luis.calderon@uni-bonn.de}



{title:See Also}
Related commands:



