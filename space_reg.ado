/******************************************************************************/
/* space_reg.ado                                                              */
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
/*	>> space_reg coordlist cutofflist depvar regressorlist,  xreg()  coord() model() */
/*																					 */
/*  NOTE: (1) If you want a constant in the regression, specify one of               */
/*	your input variables as a 1. (ie. include it in list of					         */
/*	regressors). *areg* and *reghdfe* do not need a constant.						 */
/*																					 */
/*	(2) MUST specify positive integer for xreg() and coord() options.			     */
/*	(3)	xreg() denotes # of regressors												 */
/*	(4)	coord() denotes # of coordinates 											 */
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

program define space_reg, eclass
version 13.0
syntax varlist(min=4) [if] [in] [, Xreg(int-1) COord(int-1) MODel(str) * ]
/* everything is required. */
#delimit ;				

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
if `xreg'<1{;
	if `xreg'==-1{;
		di in red "xreg() is required. Takes an integer equal to the amount of regressors in the varlist.";
		exit 198;
		};
	di in red "xreg(`xreg') is invalid";
	exit 198;
	};	


if `coord'<1{;
	if `coord'==-1{;
		di in red "coord() is required. Takes an integer equal to the amount of coordinates in the varlist.";
		exit 198;
		};
	di in red "coord(`coord') is invalid";
	exit 198;
	};

/*Separate input variables: coordinates, cutoffs, dependent var, regressors */

/*get coordinates*/
parse "`varlist'", parse(" ");	
local a=1;
while `a'<=`coord' {;
	tempvar coord`a';
	gen `coord`a''=``a'';	
local a=`a'+1;
};

/*get cutoffs*/
local aa=1;
while `aa'<=`coord' {;
	tempvar cut`aa';
	gen `cut`aa''=``a'';	
	local a=`a'+1;
local aa=`aa'+1;
};

/*get dep variable*/
tempvar Y;
gen `Y'=``a'';			
local depend : word `a' of `varlist';

local a=`a'+1;

local b=1;
while `b'<=`xreg' {;
	tempvar X`b';
	local ind`b'="`b'";
	gen `X`b''= ``a'';
	local ind`b' : word `a' of `varlist';
	local a=`a'+1;
local b=`b'+1;
};

tempfile all;
save `all';
if "`if'"==""{;
	disp("No if conditions specified");
};
else if "`if'"!=""{;
	keep `if';
};


/*Estimate the model*/
quietly{;
tempname bread;
if "`model'"=="reghdfe" | "`model'"=="areg"{;
tempname d_list ry;
	reghdfe `Y' `X1', absorb(`opts');
	keep if e(sample)==1;
	matrix betalgt = e(b);
	matrix cov_nd = e(V); /*Degrees of freedom corrected VCV matrix*/
	reghdfe `Y', absorb(`opts') resid; /*Frisch–Waugh–Lovell (FWL) theorem*/
	predict `ry', resid; 
	replace `ry'=0 if `ry'==.;
	local partialed_varlist;
	if `xreg'==1{;
		tempvar rX1;
		reghdfe `X1', absorb(`opts') resid;
		predict `rX1', resid;
		replace `rX1'=0 if `rX1'==.; /*leads to same standard errors in unclustered models*/
		reg `ry' `rX1', noconstant;
		tempvar uhat;
		predict `uhat', resid;
		mat accum `bread'= `rX1', noconstant; /*hessian matrix*/
		local X1 `rX1'; /*replace regressor with partialed regressor*/
	};
	else{;
	tempvar uhat;
	reghdfe `Y' `X1'-`X`xreg'', absorb(`opts') resid; /*Same process as above, but for more regressors*/
	keep if e(sample)==1;
	predict `uhat', residuals;
	matrix betalgt = e(b);
	matrix cov_nd = e(V);
	forvalues k=1(1)`xreg'{;
		tempvar rX`k';
		reghdfe `X`k'', absorb(`opts') resid;
		predict `rX`k'', resid;
		replace `rX`k''=0 if `rX`k''==.;
		local partialed_varlist `partialed_varlist' `rX`k'';
		replace `X`k''=`rX`k'';
	};
	mat accum `bread'=`partialed_varlist', noconstant;
	drop _reghdfe_resid;
	};
};


if "`model'"=="ols"{; 
	if `xreg'==1 {;
		reg `Y' `X1', noconstant robust;
		keep if e(sample)==1;
		mat accum `bread'=`X1', noconstant; /*Creates X'X matrix*/
		};	
	else{;
		reg `Y' `X1'-`X`xreg'', noconstant;
		keep if e(sample)==1;
		mat accum `bread'=`X1'-`X`xreg'' `if',noconstant;
		};	
	tempvar uhat;
	matrix betalgt = e(b); /*save parameter estimates and VCV matrix*/
	matrix cov_nd = e(V);
	predict `uhat', residuals;
};

if "`model'"=="logit"{;
	if `xreg'==1{;
			logit `Y' `X1',noconstant;
			keep if e(sample)==1;
			mat `bread' = -1*invsym(e(V)); /*Hessian*/
			};
	else{;
			logit `Y' `X1'-`X`xreg'' , noconstant;
			keep if e(sample)==1;
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
			probit `Y' `X1',noconstant;
			keep if e(sample)==1;
			mat `bread' = -1*invsym(e(V)); /*Hessian*/
			};
	else{;
			probit `Y' `X1'-`X`xreg'' , noconstant;
			keep if e(sample)==1;
			mat `bread' = -1*invsym(e(V));
			};
	tempvar phat uhat;
	matrix betalgt = e(b);
	matrix cov_nd = e(V);
	predict `phat'; 			
	gen `uhat' = `Y' - `phat'; 
};

else if "`model'"=="poisson"{;
	if `xreg'==1{;
			poisson `Y' `X1',noconstant;
			keep if e(sample)==1;
			mat `bread' = -1*invsym(e(V)); /*Hessian*/
			};
	else{;
			poisson `Y' `X1'-`X`xreg'' , noconstant;
			keep if e(sample)==1;
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
			nbreg `Y' `X1',noconstant;
			keep if e(sample)==1;
			mat `temp_hessian' = -1*invsym(e(V)); /*Hessian*/
			mat `bread' = J(`xreg',`xreg',0);
			forval i = 1/`xreg' {; /*we recreate the hessian, removing the alpha parameter*/
				forval j = 1/`xreg' {;
					mat `bread'[`i', `j'] = `temp_hessian'[`i', `j'];
				};
			};
		};
	else{;
			nbreg `Y' `X1'-`X`xreg'' , noconstant;
			keep if e(sample)==1;
			mat `temp_hessian' = -1*invsym(e(V));
			mat `bread' = J(`xreg',`xreg',0);
			forval i = 1/`xreg' {; 
				forval j = 1/`xreg' {;
					mat `bread'[`i', `j'] = `temp_hessian'[`i', `j'];
				};
			};
		};
	tempvar phat uhat;
	matrix betalgt = e(b);
	matrix cov_nd = e(V);
	predict `phat'; 			
	gen `uhat' = `Y' - `phat'; 
};

};


/*Bartlett window estimator for spatial correlation correction  */
quietly{;
	tempname filling half_filling_for_i window complete_filling_for_ik filling_for_i; 	/* Declare a set of variables */
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

					replace `dis`j''=abs(`coord`j''-`coord`j''[`i']); // replace a = `dis`j'';
					replace `window'=`window'*(1-`dis`j''/`cut`j'');
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
							mat vecaccum `complete_filling_for_ik'=`half_filling_for_i' `X1'-`X`xreg'', noconstant;  /*generates an a'x row vector, (1 x n)(n x k)*/
							};
					mat `filling_for_i'=nullmat(`filling_for_i') \ `complete_filling_for_ik'; /* the result is a compilation of all the rows, where each row pertains to a regressor (k x k) */
					local k=`k'+1;
					};
			
			mat `filling'=`filling' +`filling_for_i'; /*Adding variance contribution to filling*/
			local i=`i'+1;
			};
};
/*Generate corrected VCV matrix*/
tempvar inv_hessian;
matrix `inv_hessian' = inv(`bread');
mat cov_dep = (`inv_hessian'*`filling'*`inv_hessian'); 

/*Add postestimation matrices*/
mat se = J(`xreg',1,0);
forval i = 1(1)`xreg' {;
		local j=1;
		mat se[`i', `j'] = sqrt(cov_dep[`i', `i']);
	};

ereturn matrix se = se;
ereturn matrix hessian = `bread';
ereturn matrix opg = `filling';





/*Create output. Includes variable name, parameter estimate, uncorrected standard errors, standard t-test, corrected standard errors */

local v=`a';
di _newline(2)
"Results for Spatial Errors for M-estimators";
di _newline	_col(60)	" Number of observations =  "_N;
di _col(60) "Dependent variable = `depend'";
di _newline
"Variable" _col(13) "Coef Est." _col(29) "Standard SE" _col(45) "Standard t-stat" _col(65) "Spatial SE";
di 
"---------" _col(13) "----------" _col(29) "------------" _col(45) "---------------" _col(65) "------------";
local z=1;
while `z'<=`xreg' {;
	tempvar beta`z' beta1`z' se`z' see`z' se1`z' se2`z' t1`z';
	gen `beta`z''=betalgt[1,`z']; /*parameter estimate*/
	gen `se`z''=cov_nd[`z',`z'];
	gen `see`z''=sqrt(`se`z''); /*uncorrected standard error*/
	gen `se1`z''=cov_dep[`z',`z'];
	gen `se2`z''=sqrt(`se1`z''); /*corrected standard error*/
	gen `t1`z'' = `beta`z'' / `see`z'';
	di "`ind`z''" _col(13)  `beta`z''  _col(29)   `see`z'' _col(45)   `t1`z''  _col(65)  `se2`z'';
local z=`z'+1;
};
ereturn matrix conley_V = cov_dep; /*placed here because for some reason, placing it before the output throws an error*/
use `all', clear;
exit;

end;

/* START HELP FILE
title[a command to compute spatial standard errors for M-estimators]

desc[
 {cmd:space_reg} calculates spatial standard errors for several commonly
 used non-linear M-estimators. 
 
 To use, 4 arguments are required. 
 
 {cmd:coordlist} is the list of coordinates. This could be latitude and longitude.  
 
 {cmd:cutofflist} is the list of cutoff variables. There needs to be as many cutoff
 variables as coordinates. 
 
 {cmd:depvar} is the dependent variable. May be continuous, binary, or count.
 
 {cmd:regressorlist} is the regressor list. Must include a constant unless 
 {cmd:areg} or {cmd:reghdfe} are used. Fixed effects are included within {cmd:model}. 
 See example for reference.
]
Syntax[{cmd:space_reg} coordlist cutofflist depvar regressorlist]
opt[xreg() denotes # of regressors]
opt[coord() denotes # of coordinates]
opt[model() takes 1 of 7 acceptable models:
(1) ols (2) logit (3) probit (4) poisson (5) nb (6) areg (7) reghdfe
]

example[
 {stata space_reg C1 C2 cutoff1 cutoff2 dep indep1 const, xreg(2) coord(2) model(ols)}:
 
{stata space_reg C1 C2 cutoff1 cutoff2 dep indep1, xreg(1) coord(2) model(reghdfe, fe1 fe2)}

]
author[Luis Calderon and Leander Heldring]
institute[University of Bonn and Briq]
email[leander.heldring@briq-institute.org]

return[se Standard error]
return[hessian Hessian matrix]
return[opg The outer product of gradients]
return[conley_V The variance-covariance matrix]


references[
Jenish, N., & Prucha, I. R. (2009). Central limit theorems and uniform laws of 
large numbers for arrays of random fields. Journal of econometrics, 150(1), 86-98.

Conley, T. G. (1999). GMM estimation with cross sectional dependence. 
Journal of econometrics, 92(1), 1-45.
]

seealso[

]

END HELP FILE */
