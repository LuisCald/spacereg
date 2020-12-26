README
------

Authors: Luis Calderon, Leander Heldring
Version: 2020-06-10

This software is accompanied by our paper 'Spatial standard errors for several commonly used M-estimators'.
If you find any errors, please go to www.leanderheldring.com and reach out

Files:

1) demo_space_reg.do: this is a do-file with examples on how to use our package space_reg for STATA
2) demo_space_reg.ipynb: this is a Jupyter notebook with examples on how to use our package space_reg for python
3) space_reg.ado: this is the ado file which implements our spatial standard errors for STATA. Please place it in 
your ado folder
Dependencies: reghdfe
4) space_reg.py: this is the ado file which implements our spatial standard errors for python. Place it in 
the same folder as the notebook for easy demo'ing
Dependencies: Pandas, StatsModels, NumPy
5) spatial_data.dta: example data used in our demo's. This data set is taken from Tim Conley's website. We added
two limited dependent variables
6) space_reg.sthlp: this is the help file for our package. Contains information on syntax, options, references and some examples. Contact information can be found below.
-------

Requirements:

To run the python code, adhere to the following requirements
These packages are necessary to run the code.

* Using pip
pip install pandas 
pip install numpy 
pip install statsmodels

* Using anaconda
conda install -c anaconda pandas
conda install -c anaconda numpy
conda install -c conda-forge statsmodels