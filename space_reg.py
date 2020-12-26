"""Implement Conley Standard Errors for M-estimators."""
import pandas as pd
import numpy as np 
import statsmodels.discrete.discrete_model as sm
import statsmodels.regression.linear_model as lm


class SpatialStandardErrorsComputer(object):
    """Compute Conley standard errors for supported M-estimators"""

    def __init__(self, data, coordinates, cutoffs):
        """Initialize, with data, and constants
        Args:
            data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
            coordinates:  (list of str) list of two column names for the coordinates of each obs.
            cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                    distance weights
        """
        # inputs
        self._data = data
        self._coordinates = coordinates
        self._cutoffs = cutoffs

        # other attributes
        self._y_data = None
        self._x_data = None
        self._filling = None
        self._bread = None
        self._se = None

        # inputs to populate later
        self._model = None
        self._outcome_var = None
        self._rhs_vars = None

    @property
    def show_data(self):
        """Show dataset"""
        return self._data

    @property
    def show_coordinates(self):
        """Show variable names for coordinates"""
        return self._coordinates

    @property
    def show_cutoffs(self):
        """Show variable names for cutoffs"""
        return self._cutoffs

    @property
    def show_sandwich_bread(self):
        """Show the bread matrix for se computation"""
        return self._bread

    @property
    def show_sandwich_filling(self):
        """Show the filling matrix for se computation"""
        return self._filling

    @property
    def show_outcome_rhs_vars(self):
        """Return outcome and rhs vars"""
        return self._outcome_var, self._rhs_vars

    @property
    def show_model_estimated(self):
        """Return estimated model"""
        return self._model

    def compute_conley_standard_errors_all_models(self, model, y, x):
        """Compute Conley standard errors for supported M-estimators

        Supported models:
            "OLS": linear model, estimated using Ordinary Least Squares
            "logit": logistic model for binary outcome data, estimated using maximum likelihood
            "probit": probit model for binary outcome data, estimated using maximum likelihood
            "poisson": poisson GLM for count data,  estimated using maximum likelihood
            "nb": negative binomial model GLM for overdispersed count data, estimated using
                   maximum likelihood

        Args:
            model: (str) the type of model to estimate:"OLS", "logit", "probit", "poisson", or "nb"
            y: (str) name of dependent variable
            x: (list of str) list of regressors *with* constant

        Returns:
            se_df: (pd.DataFrame) column of k + 1 standard errors, for k regressors and the constant
    
        """
        # collect inputs
        self._model = model
        self._rhs_vars = self._data[x]
        self._outcome_var = self._data[y]

        # get inputs for sandwich estimator
        self._filling, self._bread = _compute_sandwich_elements(
            self._model, self._data, self._outcome_var, self._rhs_vars,
            self._coordinates, self._cutoffs)

        # sandwich estimator for standard errors of M-estimators
        self._se = np.sqrt(np.diagonal(np.dot(np.dot(self._bread, self._filling), self._bread)))

        # construct return df
        se_df = _construct_inference_table(x, self._se)

        return se_df


def compute_conley_standard_errors_all_models(model, data, y, x, coordinates, cutoffs):
    """Compute Conley standard errors for supported M-estimators

    Supported models:
            "OLS": linear model, estimated using Ordinary Least Squares
            "logit": logistic model for binary outcome data, estimated using maximum likelihood
             "probit": probit model for binary outcome data, estimated using maximum likelihood
             "poisson": poisson GLM for count data,  estimated using maximum likelihood
             "nb": negative binomial model GLM for overdispersed count data, estimated using
                   maximum likelihood

    Args:
        model: (str) the type of model to estimate:"OLS", "logit", "probit", "poisson", or "nb"
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        y: (str) name of dependent variable
        x: (list of str) list of regressors *with* constant
        coordinates:  (list of str) list of two column names for the coordinates of each obs.
        cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        se_df: (pd.DataFrame) column of k + 1 standard errors, for k regressors and the constant
        
    """
    # collect data
    x_mat = data[x]
    y_vec = data[y]

    # get inputs for sandwich estimator
    filling, bread = _compute_sandwich_elements(model, data, y_vec, x_mat, coordinates, cutoffs)

    # sandwich estimator for standard errors of M-estimators
    se = np.sqrt(np.diagonal(np.dot(np.dot(bread, filling), bread)))

    # construct return df
    se_df = _construct_inference_table(x, se)
    
    return se_df


def _construct_inference_table(column_names, se):
    return pd.Series(data=se, name="Conley s.e.", index=column_names)


def _compute_sandwich_elements(model, data, y, x, coordinates, cutoffs):
    """Compute 'bread' and 'filling' parts of sandwich estimator of variance covariance matrix

    Args:
        model: (str) the type of model to estimate:"OLS", "logit", "probit", "poisson", or "nb"
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        y: (str) name of dependent variable
        x: (list of str) list of regressors *with* constant
        coordinates:  (list of str) list of two column names for the coordinates of each obs.
        cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        filling: (np.array) K x K matrix. Computed using the Bartlett kernel, see paper
        bread: (np.array) K x K matrix, matrix of second derivatives w.r.t.
                          parameters, the 'Hessian'
                          
    """
    # select models and collect sandwich estimator inputs
    if model == "OLS":
        filling, bread = _compute_ols_conley_se(data, y, x, coordinates, cutoffs)

    elif model == "logit":
        filling, bread = _compute_logit_conley_se(data, y, x, coordinates, cutoffs)
        
    elif model == "probit":
        filling, bread = _compute_probit_conley_se(data, y, x, coordinates, cutoffs)
    
    elif model == "poisson": 
        filling, bread = _compute_poisson_conley_se(data, y, x, coordinates, cutoffs)
    
    elif model == "nb":
        filling, bread = _compute_nb_conley_se(data, y, x, coordinates, cutoffs)
    
    else:
        raise ValueError("Specified model is in list of acceptable models. "
                         "Acceptable values: OLS, logit, probit, poisson, or nb")
    
    return filling, bread


def _bartlett_window_estimator(data, x, residuals, coordinates, cutoffs):
    """Construct weights using Bartlett method, and reweight residuals

    Args:
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        x: (list of str) list of regressors *with* constant
        residuals: (np.array) residuals from estimated model, (y - predicted values)
        coordinates: (list of str) list of two column names for the coordinates of each obs.
        cutoffs: (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        bartlett_matrix: (np.array) K x K matrix. Computed using the Bartlett kernel, see paper
        
    """
    # initialize dicts and matrices
    dist_dict = {}
    bartlett_matrix = np.zeros([len(x.columns), len(x.columns)])
    for c, cutoff in enumerate(cutoffs):
        data["cutoff{}".format(c+1)] = data[cutoff]

    # implement weighting
    for i, row in enumerate(data.values):
        window = 1

        for w, coordinate in enumerate(data[coordinates].columns):
            dist_dict["dist{}".format(w+1)] = abs(data[coordinate] - data[coordinate][i])
            window *= abs((1 - dist_dict["dist{}".format(w+1)] / data["cutoff{}".format(w+1)]))
            window[dist_dict["dist{}".format(w+1)] >= data["cutoff{}".format(w+1)]] = 0
        
        window = window[:, None]
        xi_ri = residuals[i:i+1] * x[i:i+1].T 
        bartlett_for_i = np.dot(np.dot(xi_ri, residuals.T * window.T), x)

        bartlett_matrix += bartlett_for_i
        
    return bartlett_matrix


def _compute_ols_conley_se(data, y, x, coordinates, cutoffs):
    """Linear implementation of Conley Standard Errors.
            - Outer product of score vectors calculation

    Args:
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        y: (str) name of dependent variable
        x: (list of str) list of regressors *with* constant
        coordinates:  (list of str) list of two column names for the coordinates of each obs.
        cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        filling: (np.array) K x K matrix. Computed using the Bartlett kernel, see paper
        bread: (np.array) K x K matrix, matrix of second derivatives w.r.t. 
                          parameters, the 'Hessian'
                          
    """
    # Parameter estimation
    ols_model = lm.OLS(y, x)
    params = ols_model.fit().params
    # Predicted values and residuals
    pv = ols_model.predict(params)[:, None]
    resid = y[:, None] - pv
    # Outer product of score vectors with spatial dependence correction
    filling = _bartlett_window_estimator(data, x, resid, coordinates, cutoffs)
    bread = np.linalg.inv(np.dot(x.T, x))
    
    return filling, bread


def _compute_logit_conley_se(data, y, x, coordinates, cutoffs):
    """Compute Conley Standard Errors for logit estimation

    Args:
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        y: (str) name of dependent variable
        x: (list of str) list of regressors *with* constant
        coordinates:  (list of str) list of two column names for the coordinates of each obs.
        cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        filling: (np.array) K x K matrix. Computed using the Bartlett kernel, see paper
        bread: (np.array) K x K matrix, matrix of second derivatives w.r.t. 
                          parameters, the 'Hessian'
    """
    # Parameter estimation
    logit_model = sm.Logit(y, x)
    params = logit_model.fit().params
    # Predicted values, residuals, hessian
    pv = logit_model.predict(params)[:, None]
    resid = y[:, None] - pv
    hes = logit_model.hessian(params)
    # Outer product of score vectors with spatial dependence correction
    filling = _bartlett_window_estimator(data, x, resid, coordinates, cutoffs)
    bread = np.linalg.inv(hes)
    
    return filling, bread


def _compute_probit_conley_se(data, y, x, coordinates, cutoffs):
    """Compute Conley Standard Errors for Probit regressions
            - Outer product of score vectors calculation

    Args:
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        y: (str) name of dependent variable
        x: (list of str) list of regressors *with* constant
        coordinates:  (list of str) list of two column names for the coordinates of each obs.
        cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        filling: (np.array) K x K matrix. Computed using the Bartlett kernel, see paper
        bread: (np.array) K x K matrix, matrix of second derivatives w.r.t. 
                          parameters, the 'Hessian'    
                          
    """
    # Parameter estimation
    probit_model = sm.Probit(y, x)
    params = probit_model.fit().params
    # Predicted values, residuals, hessian
    pv = probit_model.predict(params)[:, None]
    resid = y[:, None] - pv
    hes = probit_model.hessian(params)
    # Outer product of score vectors with spatial dependence correction
    filling = _bartlett_window_estimator(data, x, resid, coordinates, cutoffs)
    bread = np.linalg.inv(hes)
    
    return filling, bread


def _compute_poisson_conley_se(data, y, x, coordinates, cutoffs):
    """Compute Conley Standard Errors for Poisson regression

    Args:
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        y: (str) name of dependent variable
        x: (list of str) list of regressors *with* constant
        coordinates:  (list of str) list of two column names for the coordinates of each obs.
        cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        filling: (np.array) K x K matrix. Computed using the Bartlett kernel, see paper
        bread: (np.array) K x K matrix, matrix of second derivatives w.r.t. 
                          parameters, the 'Hessian'    
    """
    # Parameter estimation
    poisson_model = sm.Poisson(y, x)
    params = poisson_model.fit().params
    # Predicted values, residuals, hessian
    pv = poisson_model.predict(params)[:, None]
    resid = y[:, None] - pv
    hes = poisson_model.hessian(params)
    # Outer product of score vectors with spatial dependence correction
    filling = _bartlett_window_estimator(data, x, resid, coordinates, cutoffs)
    bread = np.linalg.inv(hes)
    
    return filling, bread


def _compute_nb_conley_se(data, y, x, coordinates, cutoffs):
    """Compute Conley Standard Errors for negative binomial regression.

            Note: we compute standard errors for the the 'NB2' model

    Args:
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        y: (str) name of dependent variable
        x: (list of str) list of regressors *with* constant
        coordinates:  (list of str) list of two column names for the coordinates of each obs.
        cutoffs:  (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        filling: (np.array) K x K matrix. Computed using the Bartlett kernel, see paper
        bread: (np.array) K x K matrix, matrix of second derivatives w.r.t. 
                          parameters, the 'Hessian'   
                          
    """
    # Parameter estimation
    nb_model = sm.NegativeBinomial(y, x)
    params = nb_model.fit().params
    # Predicted values, residuals, hessian
    pv = nb_model.predict(params)[:, None]
    resid = y[:, None] - pv
    hes = nb_model.hessian(params)[:2, :2]
    # Outer product of score vectors with spatial dependence correction
    filling = _bartlett_window_estimator(data, x, resid, coordinates, cutoffs)
    bread = np.linalg.inv(hes)
    
    return filling, bread

