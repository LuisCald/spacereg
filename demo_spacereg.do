						* DEMO ON spacereg *
						* ----------------- *
		* Welcome to our demo. Be sure to fully read all comments as well as the 
		* Preamble in the spacereg.ado
		
* Change PATH to the directory containing the project files
local path "add_your/path/here"
local ado_path "add_your/ado/path/here"
cd `path'

* Import data from project files folder
clear all
use "spatial_data.dta"

* To use our package, you must add the .ado file to one of the .ado locations
* To see where your .ado files are located, run:
sysdir

* Store the spacereg.ado file in the ado folder or in the project files folder.
* Set the PERSONAL directory to the project files path. 
sysdir set PERSONAL `ado_path'

* To use our package, you need the following dependencies
// ssc install reghdfe 


						* VALIDATION *
						* ---------- *
			* We validate our results by comparing OLS and logit results from 
			* Conley's packages. Results are equivalent.
			
* OLS model ran in the paper
reg dep indep
x_ols C1 C2 cutoff1 cutoff2 dep indep1 const, xreg(2) coord(2) /* For comparison. This is Tim Conley's implementation. x_ols is in the supporting .ado's */
// drop epsilon window dis1 dis2  // Conley generates these. 
* If cutoffs are stored as variables, convert to numeric values once.
* (Works whether cutoffs are constant or vary slightly; uses the sample mean.)
quietly summarize cutoff1
global cutoff1_val = r(mean) // can also be a local

quietly summarize cutoff2
global cutoff2_val = r(mean)

cap gen const_v = 1

spacereg dep indep1 const_v, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(ols)
spacereg dep indep1 const_v, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(ols) kernel(uniform)
* Logit
logit binarydep indep1
xgmlt C1 C2 cutoff1 cutoff2 binarydep indep1 const, xreg(2) coord(2) /*check supporting .ado folder */
spacereg binarydep indep1 const, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(logit)

						* ADDITIONAL FUNCTIONS *
						* -------------------- *
* Probit
spacereg binarydep indep1 const, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(probit)

* Poisson
spacereg poissondep indep1 const, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(poisson)

* Negative Binomial
spacereg poissondep indep1 const, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(nb) debug

* Fixed Effects
spacereg dep indep1, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(areg, fe1)
spacereg dep indep1, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(reghdfe, fe1 fe2)


						* SPEED *
						* ----- *
* 1000 obs
expand 10
spacereg dep indep1, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(reghdfe, fe1 fe2)

* 3000 obs
expand 3
spacereg dep indep1, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(reghdfe, fe1 fe2)

* 6000 obs
expand 2
spacereg dep indep1, coords(C1 C2) cutoffs($cutoff1_val $cutoff2_val) model(reghdfe, fe1 fe2)
