{smcl}
{* *! version 1.0 14 Nov 2020}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install command2" "ssc install command2"}{...}
{vieweralsosee "Help command2 (if installed)" "help command2"}{...}
{viewerjumpto "Syntax" "space_reg##syntax"}{...}
{viewerjumpto "Description" "space_reg##description"}{...}
{viewerjumpto "Options" "space_reg##options"}{...}
{viewerjumpto "Remarks" "space_reg##remarks"}{...}
{viewerjumpto "Examples" "space_reg##examples"}{...}
{title:Title}
{phang}
{bf:space_reg} {hline 2} a command to compute spatial standard errors for M-estimators

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:space_reg}
varlist(min=4)
[{help if}]
[{help in}]
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Optional}
{synopt:{opt x:reg(#)}}  denotes # of regressors

{pstd}
{p_end}
{synopt:{opt co:ord(#)}}  denotes # of coordinates

{pstd}
{p_end}
{synopt:{opt mod:el(str)}}  takes 1 of 7 acceptable models:
(1) ols (2) logit (3) probit (4) poisson (5) nb (6) areg (7) reghdfe
{p_end}
{synopt:{opt *}} {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}
{pstd}

{pstd}
 {cmd:space_reg} calculates spatial standard errors for several commonly
 used non-linear M-estimators. 

{pstd}
 To use, 4 arguments are required. 

{pstd}
 {cmd:coordlist} is the list of coordinates. This could be latitude and longitude.  

{pstd}
 {cmd:cutofflist} is the list of cutoff variables. There needs to be as many cutoff
 variables as coordinates. 

{pstd}
 {cmd:depvar} is the dependent variable. May be continuous, binary, or count.

{pstd}
 {cmd:regressorlist} is the regressor list. Must include a constant unless 
 {cmd:areg} or {cmd:reghdfe} are used. Fixed effects are included within {cmd:model}. 
 See example for reference.

{marker options}{...}
{title:Options}
{dlgtab:Main}
{phang}
{opt x:reg(#)}     denotes # of regressors

{pstd}
{p_end}
{phang}
{opt co:ord(#)}     denotes # of coordinates

{pstd}
{p_end}
{phang}
{opt mod:el(str)}     takes 1 of 7 acceptable models:
(1) ols (2) logit (3) probit (4) poisson (5) nb (6) areg (7) reghdfe
{p_end}
{phang}
{opt *}  {p_end}


{marker examples}{...}
{title:Examples}
{pstd}

{pstd}
 {stata space_reg C1 C2 cutoff1 cutoff2 dep indep1 const, xreg(2) coord(2) model(ols)}:

{pstd}
{stata space_reg C1 C2 cutoff1 cutoff2 dep indep1, xreg(1) coord(2) model(reghdfe, fe1 fe2)}

{pstd}

{title:Stored results}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(se)}}  Standard error {p_end}
{synopt:{cmd:e(hessian)}}  Hessian matrix {p_end}
{synopt:{cmd:e(opg)}}  The outer product of gradients {p_end}
{synopt:{cmd:e(conley_V)}}  The variance-covariance matrix {p_end}


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

Email {browse "mailto:leander.heldring@briq-institute.org":leander.heldring@briq-institute.org}



{title:See Also}
Related commands:



