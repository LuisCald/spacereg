/******************************************************************************/
/* spacereg.ado                                                              */
/* for STATA 15.0                                                             */
/*				Spatial standard errors for M-estimators			          */
/*				by	Luis Calderon (1) and Leander Heldring (2)                */
/*				(1) University of Bonn  (2) Northwestern University 	      */
/*				June 03, 2020                                                 */
/*			                                    		                      */
/*											                                  */
/******************************************************************************/
/******************************************************************************/


* WARNING: This program is offered without any guarantee. 
* If you find errors, please contact s6lucald@uni-bonn.de - or - (Leander's email)



/*  To use, type:                                        */
/*	>> spacereg depvar regressorlist, coords(varlist) cutoffs(numlist) model() */
/*																					 */
/*  NOTE: (1) If you want a constant in the regression, specify one of               */
/*	your input variables as a 1. (ie. include it in list of					         */
/*	regressors). *areg* and *reghdfe* do not need a constant.						 */
/*																					 */
/*	(2) Use coords() to provide coordinate variables.				     */
/*	(3) cutoffs() is a numeric list with one cutoff per coordinate.	     */
/*	(5)	model() takes 1 of 7 acceptable models: 								 	 */
/*   			- ols, logit, probit, poisson, nb, areg, reghdfe 				     */
/*				- nb stands for negative binomial							         */
/*																					 */
/*	(6) Your cutofflist must correspond to coordlist (same order)			         */
/* cutofflist is the values of LM. Zero weight is put on terms                       */
/* more than LM units apart in that dimmension in the calculation                    */
/* of the covariance  term															 */
/*																					 */
/*																					 */
/******************************************************************************/

program define spacereg, eclass
version 13.0
/*
Syntax:
	spacereg depvar regressorlist, coords(varlist) cutoffs(numlist) model(str)
*/
#delimit ;

/* New syntax only (legacy syntax removed) */
syntax varlist(min=2) [if] [in] , COORDS(varlist min=1) CUTOffs(numlist min=1) MODel(str) [KERnel(str) * ] ;

local kernel = lower("`kernel'");
if "`kernel'"=="" {;
	local kernel "bartlett";
};
if !("`kernel'"=="bartlett" | "`kernel'"=="uniform") {;
	di in red "kernel(`kernel') is invalid. Supported values: bartlett (default), uniform";
	exit 198;
};

/* raise errors if model option incorrectly specified */
if strpos("`model'", "areg"){;
	if "`model'"=="areg" {;
		di in red "areg requires a comma followed by your fixed effect";
		exit 198;
	};
	else if "`model'"!="areg"{;
			if strpos("`model'", "areg,"){;
				gettoken model opts : model, parse(","); 
				gettoken comma opts : opts, parse(","); 
				local byby `"model(`model', `opts')"'; 
				if "`opts'"==""{;
					di "`opts'";
					di in red "You are missing the fixed effect for areg";
					exit 198;
				};
			};
	};
	
};

if strpos("`model'", "reghdfe"){;
	if "`model'"=="reghdfe" {;
		di in red "areg requires a comma followed by your fixed effect";
		exit 198;
	};
	else if "`model'"!="reghdfe"{;
			if strpos("`model'", "reghdfe,"){;
				gettoken model opts : model, parse(","); 
				gettoken comma opts : opts, parse(","); 
				local byby `"model(`model', `opts')"'; 
				if "`opts'"==""{;
					di in red "You are missing the fixed effect for reghdfe";
					exit 198;
				};
			};
	};
	
};

/* raise errors if syntax is incorrectly specified */
preserve;

/* Separate input variables: coordinates, cutoffs, dependent var, regressors */
local coordvars "";
local depend "";
local regvars "";

local coordvars "`coords'";
local coord : word count `coordvars';
if `coord'<1 {;
	di in red "coords() must contain at least 1 variable";
	exit 198;
};

gettoken depend regvars : varlist;
local xreg : word count `regvars';
if `xreg'<1 {;
	di in red "You must specify at least 1 regressor after the dependent variable.";
	exit 198;
};

/* Handle numeric cutoffs: allow single scalar replicated across dimensions */
local _use_numeric_cutoffs 0;
local _cutoff_numlist "`cutoffs'";
if "`_cutoff_numlist'"!="" {;
	local _use_numeric_cutoffs 1;
	local _ncut : word count `_cutoff_numlist';
	if `_ncut'==1 & `coord'>1 {;
		local _c1 : word 1 of `_cutoff_numlist';
		local _cutoff_numlist "";
		forval j=1/`coord' {;
			local _cutoff_numlist "`_cutoff_numlist' `_c1'";
		};
	};
	local _ncut2 : word count `_cutoff_numlist';
	if `_ncut2'!=`coord' {;
		di in red "cutoffs() must be either a single number or have the same length as coords() / coord().";
		exit 198;
	};
};

/* Create working copies of coordinates */
tokenize "`coordvars'";
forval j=1/`coord' {;
	tempvar coord`j';
	gen double `coord`j'' = ``j'';
};

/* Create working cutoffs: either constant numeric cutoffs or copied cutoff vars */
if `_use_numeric_cutoffs'==1 {;
	tokenize "`_cutoff_numlist'";
	forval j=1/`coord' {;
		tempvar cut`j';
		gen double `cut`j'' = ``j'';
	};
};
else {;
	/* With new syntax, numeric cutoffs are always provided via cutoffs() */
	di in red "cutoffs() must be provided as a numeric list (e.g., cutoffs(100 100)).";
	exit 198;
};

/* Dependent variable */
tempvar Y;
gen double `Y' = `depend';

/* Regressors */
tokenize "`regvars'";
local Xlist "";
forval k=1/`xreg' {;
	tempvar X`k';
	gen double `X`k'' = ``k'';
	local Xlist "`Xlist' `X`k''";
	local ind`k' "``k''";
};


if "`if'"!=""{;
	keep `if';
};


/*Estimate the model*/
quietly{;
tempname esthold;
tempname bread;
if "`model'"=="reghdfe" | "`model'"=="areg"{;
	tempname d_list ry;
	local partialed_varlist;
	/* Estimate on original variables so stored stripes match user input */
	if "`model'"=="reghdfe" {;
		reghdfe `depend' `regvars', absorb(`opts') noconstant;
	};
	else {;
		/* areg supports one absorb variable */
		areg `depend' `regvars', absorb(`opts');
	};
	tempvar _esamp;
	gen byte `_esamp' = e(sample);
	keep if `_esamp'==1;
	estimates store `esthold';
	matrix betalgt = e(b);
	matrix cov_nd = e(V);

	/* Frisch–Waugh–Lovell residualization to compute Conley variance */
	if "`model'"=="reghdfe" {;
		reghdfe `Y', absorb(`opts') resid;
		predict `ry', resid;
	};
	else {;
		areg `Y', absorb(`opts');
		predict `ry', resid;
	};
	replace `ry'=0 if `ry'==.;

	if `xreg'==1 {;
		tempvar rX1;
		if "`model'"=="reghdfe" {;
			reghdfe `X1', absorb(`opts') resid;
			predict `rX1', resid;
		};
		else {;
			areg `X1', absorb(`opts');
			predict `rX1', resid;
		};
		replace `rX1'=0 if `rX1'==.;
		reg `ry' `rX1', noconstant;
		tempvar uhat;
		predict `uhat', resid;
		mat accum `bread'= `rX1', noconstant;
		local X1 `rX1';
	};
	else {;
		tempvar uhat;
		/* Residualize y for FE model to obtain uhat */
		if "`model'"=="reghdfe" {;
			reghdfe `Y' `Xlist', absorb(`opts') resid;
			predict `uhat', residuals;
		};
		else {;
			areg `Y' `Xlist', absorb(`opts');
			predict `uhat', resid;
		};
		forvalues k=1(1)`xreg' {;
			tempvar rX`k';
			if "`model'"=="reghdfe" {;
				reghdfe `X`k'', absorb(`opts') resid;
				predict `rX`k'', resid;
			};
			else {;
				areg `X`k'', absorb(`opts');
				predict `rX`k'', resid;
			};
			replace `rX`k''=0 if `rX`k''==.;
			local partialed_varlist `partialed_varlist' `rX`k'';
			replace `X`k''=`rX`k'';
		};
		mat accum `bread'=`partialed_varlist', noconstant;
		capture drop _reghdfe_resid;
	};
};


if "`model'"=="ols"{; 
	if `xreg'==1 {;
		reg `depend' `regvars', noconstant robust;
		tempvar _esamp3;
		gen byte `_esamp3' = e(sample);
		keep if `_esamp3'==1;
		estimates store `esthold';
		mat accum `bread'=`X1', noconstant; /*Creates X'X matrix*/
		};	
	else{;
		reg `depend' `regvars', noconstant;
		tempvar _esamp4;
		gen byte `_esamp4' = e(sample);
		keep if `_esamp4'==1;
		estimates store `esthold';
		mat accum `bread'=`Xlist' `if',noconstant;
		};	
	tempvar uhat;
	matrix betalgt = e(b); /*save parameter estimates and VCV matrix*/
	matrix cov_nd = e(V);
	predict `uhat', residuals;
};

if "`model'"=="logit"{;
	if `xreg'==1{;
			logit `depend' `regvars',noconstant;
			tempvar _esamp5;
			gen byte `_esamp5' = e(sample);
			keep if `_esamp5'==1;
			estimates store `esthold';
			mat `bread' = -1*invsym(e(V)); /*Hessian*/
			};
	else{;
			logit `depend' `regvars' , noconstant;
			tempvar _esamp6;
			gen byte `_esamp6' = e(sample);
			keep if `_esamp6'==1;
			estimates store `esthold';
			mat `bread' = -1*invsym(e(V));
			};
	tempvar phat uhat;
	matrix betalgt = e(b);
	matrix cov_nd = e(V);
	predict `phat'; 			
	gen `uhat' = `Y' - `phat'; 
};

else if "`model'"=="probit"{;
	if `xreg'==1{;
			probit `depend' `regvars',noconstant;
			tempvar _esamp7;
			gen byte `_esamp7' = e(sample);
			keep if `_esamp7'==1;
			estimates store `esthold';
			mat `bread' = -1*invsym(e(V)); /*Hessian*/
			};
	else{;
			probit `depend' `regvars' , noconstant;
			tempvar _esamp8;
			gen byte `_esamp8' = e(sample);
			keep if `_esamp8'==1;
			estimates store `esthold';
			mat `bread' = -1*invsym(e(V));
			};
	/* For probit, the per-observation score contribution rescales (y-p) by
	   phi(xb)/(p(1-p)). We build uhat as this generalized residual so that
	   X*uhat corresponds to the score vector used in the sandwich meat. */
	tempvar phat uhat xb denom;
	matrix betalgt = e(b);
	matrix cov_nd = e(V);
	predict double `xb', xb;
	predict double `phat', pr;
	gen double `denom' = `phat'*(1-`phat');
	replace `denom' = 1e-12 if `denom' < 1e-12;
	gen double `uhat' = (`Y' - `phat')*normalden(`xb')/`denom';
    replace `uhat' = 0 if `uhat'==.;
};

else if "`model'"=="poisson"{;
	if `xreg'==1{;
			poisson `depend' `regvars',noconstant;
			tempvar _esamp9;
			gen byte `_esamp9' = e(sample);
			keep if `_esamp9'==1;
			estimates store `esthold';
			mat `bread' = -1*invsym(e(V)); /*Hessian*/
			};
	else{;
			poisson `depend' `regvars' , noconstant;
			tempvar _esamp10;
			gen byte `_esamp10' = e(sample);
			keep if `_esamp10'==1;
			estimates store `esthold';
			mat `bread' = -1*invsym(e(V));
			};
	tempvar phat uhat;
	matrix betalgt = e(b);
	matrix cov_nd = e(V);
	predict `phat'; 			
	gen `uhat' = `Y' - `phat'; 
};

else if "`model'"=="nb"{;
tempvar temp_hessian;
	if `xreg'==1{;
			nbreg `depend' `regvars',noconstant;
			tempvar _esamp11;
			gen byte `_esamp11' = e(sample);
			keep if `_esamp11'==1;
			estimates store `esthold';
			mat `temp_hessian' = -1*invsym(e(V)); /*Hessian*/
			mat `bread' = J(`xreg',`xreg',0);
			forval i = 1/`xreg' {; /*we recreate the hessian, removing the alpha parameter*/
				forval j = 1/`xreg' {;
					mat `bread'[`i', `j'] = `temp_hessian'[`i', `j'];
				};
			};
		};
	else{;
			nbreg `depend' `regvars' , noconstant;
			tempvar _esamp12;
			gen byte `_esamp12' = e(sample);
			keep if `_esamp12'==1;
			estimates store `esthold';
			mat `temp_hessian' = -1*invsym(e(V));
			mat `bread' = J(`xreg',`xreg',0);
			forval i = 1/`xreg' {; 
				forval j = 1/`xreg' {;
					mat `bread'[`i', `j'] = `temp_hessian'[`i', `j'];
				};
			};
		};
	/* For nbreg, the score contribution for beta is proportional to
	   x_i * (y_i - mu_i)/(1 + alpha*mu_i). Build uhat as this generalized
	   residual so X*uhat corresponds to the score vector used in the meat. */
	tempvar phat uhat;
	tempname alpha;
	matrix betalgt = e(b);
	matrix cov_nd = e(V);
	/* Extract alpha (overdispersion). Prefer e(alpha); fall back to exp(_b[/lnalpha]). */
	capture scalar `alpha' = e(alpha);
	if _rc {;
		capture scalar `alpha' = exp(_b[/lnalpha]);
	};
	if _rc {;
		di in red "Could not extract alpha from nbreg results (expected e(alpha) or _b[/lnalpha]).";
		exit 498;
	};
	/* nbreg does not support the mu option; use n for predicted mean count */
	predict double `phat', n;
	gen double `uhat' = (`Y' - `phat')/(1 + `alpha'*`phat');
    replace `uhat' = 0 if `uhat'==.;
};

};


/*Bartlett window estimator for spatial correlation correction  */
quietly{;
	tempname filling complete_filling_for_ik filling_for_i; 	/* Declare a set of matrices */
	tempvar half_filling_for_i window;
	matrix `filling' = J(`xreg',`xreg',0);  /* Intializes the matrix */
	gen `half_filling_for_i'=0;
	gen `window'=1;		/*initializes mat.s/var.s to be used*/
	local i=1;

	while `i'<=_N {;		/*loop through all observations*/
			local j=1;
			replace `window'=1;
			while `j'<=`coord' {;	/*loop through coordinates*/
					if `i'==1{;
							tempvar dis`j';
							gen `dis`j''=0;
							};

					replace `dis`j''=abs(`coord`j''-`coord`j''[`i']);
					if "`kernel'"=="bartlett" {;
						replace `window'=`window'*(1-`dis`j''/`cut`j'');
					};
					replace `window'=0 if `dis`j''>=`cut`j'';
					local j=`j'+1;
					};			

			/* End of j loop */

			capture mat drop `filling_for_i';
			local k=1;
			while `k'<=`xreg' {; /*Iterate over each regressor, correcting for spatial dependence*/
					replace `half_filling_for_i'=`X`k''[`i']*`uhat'*`uhat'[`i']*`window'; /*part of the filling for i and regressor k*/
					if `xreg'==1{;
							mat vecaccum `complete_filling_for_ik'=`half_filling_for_i' `X1', noconstant; /*generates an a'x matrix, (1 x n)(n x 1)*/
							};
					else{; /*Each iteration creates complete_filling_for_ik, which is filling for this regressor*/
							mat vecaccum `complete_filling_for_ik'=`half_filling_for_i' `Xlist', noconstant;  /*generates an a'x row vector, (1 x n)(n x k)*/
							};
					mat `filling_for_i'=nullmat(`filling_for_i') \ `complete_filling_for_ik'; /* the result is a compilation of all the rows, where each row pertains to a regressor (k x k) */
					local k=`k'+1;
					};
			
			mat `filling'=`filling' +`filling_for_i'; /*Adding variance contribution to filling*/
			local i=`i'+1;
			};
};
/*Generate corrected VCV matrix*/
tempname inv_hessian;
matrix `inv_hessian' = inv(`bread');
mat cov_dep = (`inv_hessian'*`filling'*`inv_hessian'); 


mat se = J(`xreg',1,0);
forval i = 1(1)`xreg' {;
		local j=1;
		mat se[`i', `j'] = sqrt(cov_dep[`i', `i']);
	};

	/* Restore original estimation results so replay output (LL, R2, chi2, etc.)
	   matches Stata's native output for the chosen model, then swap in the
	   Conley covariance via ereturn repost. */
	tempname V_conley V_uncorr V_full;
	matrix `V_conley' = cov_dep;
	quietly estimates restore `esthold';
	matrix `V_uncorr' = e(V);
	matrix `V_full' = `V_uncorr';
	/* Overwrite the coefficient block with Conley covariance.
	   For nbreg, this preserves the ancillary-parameter variance/covariances. */
	forval ii=1/`xreg' {;
		forval jj=1/`xreg' {;
			matrix `V_full'[`ii',`jj'] = `V_conley'[`ii',`jj'];
		};
	};
	/* Post Conley V while keeping the original model results active */
	ereturn repost V=`V_full';

	ereturn local vcetype "Conley";
	ereturn local vce "conley";
	ereturn local kernel "`kernel'";
	ereturn local conley_cmd "spacereg";
	ereturn local conley_model "`model'";

	ereturn matrix se = se;
	ereturn matrix hessian = `bread';
	ereturn matrix V_uncorrected = `V_uncorr';
	ereturn matrix opg = `filling';

	if "`model'"=="ols" {;
		regress;
	};
	else if "`model'"=="logit" {;
		logit;
	};
	else if "`model'"=="probit" {;
		probit;
	};
	else if "`model'"=="poisson" {;
		poisson;
	};
	else if "`model'"=="nb" {;
		nbreg;
	};
	else if "`model'"=="areg" {;
		areg;
	};
	else if "`model'"=="reghdfe" {;
		reghdfe;
	};
	else {;
		ereturn display;
	};

restore;
exit;

end;

/* START HELP FILE
title[a command to compute spatial standard errors for M-estimators]

desc[
 {cmd:spacereg} calculates spatial standard errors for several commonly
 used non-linear M-estimators. 
 
 In the recommended syntax, spatial correction is defined by:
 {p_end}
 {phang}
 {cmd:coords()} the coordinate variables (e.g., latitude/longitude)
 {p_end}
 {phang}
 {cmd:cutoffs()} the kernel cutoffs as numeric values (either one value, replicated across dimensions, or one value per coordinate)
 {p_end}

 {cmd:depvar} is the dependent variable. May be continuous, binary, or count.

 {cmd:regressorlist} is the regressor list.

 Important: for {cmd:ols}, {cmd:logit}, {cmd:probit}, {cmd:poisson}, and {cmd:nb}, the underlying estimator is run with {cmd:noconstant}. If you want an intercept, include a constant regressor (e.g., {cmd:gen byte const=1} and include {cmd:const} in {cmd:regressorlist}).
]
Syntax[{cmd:spacereg} depvar regressorlist]
opt[model() takes 1 of 7 acceptable models:
(1) ols (2) logit (3) probit (4) poisson (5) nb (6) areg (7) reghdfe
]

opt[kernel() chooses the spatial kernel/weighting scheme:
(1) bartlett (default) (2) uniform
]

example[
 {stata spacereg dep indep1 const, coords(C1 C2) cutoffs(100 100) model(ols)}
 {stata spacereg dep indep1, coords(C1 C2) cutoffs(100 100) model(reghdfe, fe1 fe2)}
 {stata spacereg dep indep1 const, coords(C1 C2) cutoffs(100 100) model(ols) kernel(uniform)}
]
author[Luis Calderon and Leander Heldring]
institute[University of Bonn and Briq]
email[leander.heldring@briq-institute.org]

return[se Standard error]
return[hessian Hessian matrix]
return[opg The outer product of gradients]
return[V_uncorrected Uncorrected model variance-covariance matrix]


references[
Jenish, N., & Prucha, I. R. (2009). Central limit theorems and uniform laws of 
large numbers for arrays of random fields. Journal of econometrics, 150(1), 86-98.

Conley, T. G. (1999). GMM estimation with cross sectional dependence. 
Journal of econometrics, 92(1), 1-45.
]

seealso[

]

END HELP FILE */
