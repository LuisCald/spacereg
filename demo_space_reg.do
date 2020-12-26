						* DEMO ON SPACE_REG *
						* ----------------- *
		* Welcome to our demo. Be sure to fully read all comments as well as the 
		* Preamble in the space_reg.ado
		
* Change PATH to the directory containing the project files
local path ""
cd `path'

* Import data from project files folder
clear all
use "spatial_data.dta"

* To use our package, you must add the .ado file to one of the .ado locations
* To see where your .ado files are located, run:
sysdir

* The space_reg.ado file can be found in the project files folder.
* Set the PERSONAL directory to the project files path. 
sysdir set PERSONAL `path'

* To use our package, you need the following dependencies
ssc install reghdfe 


						* VALIDATION *
						* ---------- *
			* We validate our results by comparing OLS and logit results from 
			* Conley's packages. Results are equivalent.
			
* OLS model ran in the paper
reg dep indep
x_ols C1 C2 cutoff1 cutoff2 dep indep1 const, xreg(2) coord(2) /* For comparison. This is Tim Conley's implementation. Download x_ols from his website if you want to use it */
drop epsilon window dis1 dis2  // Conley generates these. 

space_reg C1 C2 cutoff1 cutoff2 dep indep1 const, xreg(2) coord(2) model(ols)

* Logit
logit binary_dep indep1
xgmlt C1 C2 cutoff1 cutoff2 binary_dep indep1 const, xreg(2) coord(2) /*For comparison. this implementation is available from Tim Conley's website */
space_reg C1 C2 cutoff1 cutoff2 binary_dep indep1 const, xreg(2) coord(2) model(logit)

						* ADDITIONAL FUNCTIONS *
						* -------------------- *
* Probit
space_reg C1 C2 cutoff1 cutoff2 binary_dep indep1 const, xreg(2) coord(2) model(probit)

* Poisson
space_reg C1 C2 cutoff1 cutoff2 poisson_dep indep1 const, xreg(2) coord(2) model(poisson)

* Negative Binomial
space_reg C1 C2 cutoff1 cutoff2 poisson_dep indep1 const, xreg(2) coord(2) model(nb)

* Fixed Effects
space_reg C1 C2 cutoff1 cutoff2 dep indep1, xreg(1) coord(2) model(areg, fe1)
space_reg C1 C2 cutoff1 cutoff2 dep indep1, xreg(1) coord(2) model(reghdfe, fe1 fe2)


						* SPEED *
						* ----- *
* 1000 obs
expand 10
space_reg C1 C2 cutoff1 cutoff2 dep indep1, xreg(1) coord(2) model(reghdfe, fe1 fe2)

* 3000 obs
expand 3
space_reg C1 C2 cutoff1 cutoff2 dep indep1, xreg(1) coord(2) model(reghdfe, fe1 fe2)

* 6000 obs
expand 2
space_reg C1 C2 cutoff1 cutoff2 dep indep1, xreg(1) coord(2) model(reghdfe, fe1 fe2)
