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
        self._coordinates = coordinates
        self._data = data
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

    def compute_conley_standard_errors_all_models(self, model, y, x, kernel="bartlett"):
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
            kernel: (str) weighting scheme for spatial dependence. Supported: "bartlett" (default),
                    "uniform" (uniform weights within cutoffs)

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
            self._coordinates, self._cutoffs, kernel=kernel)

        # Numerical hygiene: enforce symmetry before forming the sandwich.
        # The Conley/HAC meat should be symmetric under symmetric weights.
        self._filling = 0.5 * (self._filling + self._filling.T)
        self._bread = 0.5 * (self._bread + self._bread.T)

        # sandwich estimator for standard errors of M-estimators
        # Finite-sample degrees-of-freedom correction: n / (n - k)
        n = len(self._outcome_var)
        k = self._rhs_vars.shape[1]
        var = (n / (n - k)) * self._bread @ self._filling @ self._bread.T
        self._se = np.sqrt(np.clip(np.diagonal(var), 0.0, None))

        # construct return df
        se_df = _construct_inference_table(x, self._se)

        return se_df


def compute_conley_standard_errors_all_models(model, data, y, x, coordinates, cutoffs, kernel="bartlett"):
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
        kernel: (str) weighting scheme for spatial dependence. Supported: "bartlett" (default),
            "uniform" (uniform weights within cutoffs)

    Returns:
        se_df: (pd.DataFrame) column of k + 1 standard errors, for k regressors and the constant
        
    """
    # collect data
    x_mat = data[x]
    y_vec = data[y]

    # get inputs for sandwich estimator
    filling, bread = _compute_sandwich_elements(
        model, data, y_vec, x_mat, coordinates, cutoffs, kernel=kernel)

    filling = 0.5 * (filling + filling.T)
    bread = 0.5 * (bread + bread.T)

    # sandwich estimator for standard errors of M-estimators
    # Finite-sample degrees-of-freedom correction: n / (n - k)
    n = len(y_vec)
    k = x_mat.shape[1]
    var = (n / (n - k)) * bread @ filling @ bread.T
    se = np.sqrt(np.clip(np.diagonal(var), 0.0, None))

    # construct return df
    se_df = _construct_inference_table(x, se)
    
    return se_df


def _construct_inference_table(column_names, se):
    return pd.Series(data=se, name="Conley s.e.", index=column_names)


def _compute_sandwich_elements(model, data, y, x, coordinates, cutoffs, kernel="bartlett"):
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
        filling, bread = _compute_ols_conley_se(data, y, x, coordinates, cutoffs, kernel=kernel)

    elif model == "logit":
        filling, bread = _compute_logit_conley_se(data, y, x, coordinates, cutoffs, kernel=kernel)
        
    elif model == "probit":
        filling, bread = _compute_probit_conley_se(data, y, x, coordinates, cutoffs, kernel=kernel)
    
    elif model == "poisson": 
        filling, bread = _compute_poisson_conley_se(data, y, x, coordinates, cutoffs, kernel=kernel)
    
    elif model == "nb":
        filling, bread = _compute_nb_conley_se(data, y, x, coordinates, cutoffs, kernel=kernel)
    
    else:
        raise ValueError("Specified model is in list of acceptable models. "
                         "Acceptable values: OLS, logit, probit, poisson, or nb")
    
    return filling, bread


def _bartlett_window_estimator(data, x, residuals, coordinates, cutoffs, kernel="bartlett"):
    """Construct weights and reweight residuals

    Args:
        data: (pd.DataFrame) dataframe containing y, x, coordinates, and cutoffs
        x: (list of str) list of regressors *with* constant
        residuals: (np.array) residuals from estimated model, (y - predicted values)
        coordinates: (list of str) list of two column names for the coordinates of each obs.
        cutoffs: (list of str) list of two column names for the cutoffs to use when computing
                distance weights

    Returns:
        weight_matrix: (np.array) K x K matrix. Computed using the selected kernel
        
    """
    kernel = str(kernel).strip().lower()
    if kernel not in {"bartlett", "uniform"}:
        raise ValueError(
            f"Unsupported kernel '{kernel}'. Supported values: 'bartlett', 'uniform'.")

    x_values = x.to_numpy(dtype=float, copy=False)
    residuals = np.asarray(residuals, dtype=float)
    if residuals.ndim == 1:
        residuals = residuals.reshape(-1, 1)

    coord_values = data[coordinates].to_numpy(dtype=float, copy=False)
    cutoff_values = np.column_stack([data[c].to_numpy(dtype=float, copy=False) for c in cutoffs])

    # initialize matrices
    weight_matrix = np.zeros((x_values.shape[1], x_values.shape[1]), dtype=float)

    # implement weighting
    n_obs = x_values.shape[0]
    for i in range(n_obs):
        window = np.ones(n_obs, dtype=float)

        for w in range(coord_values.shape[1]):
            dist = np.abs(coord_values[:, w] - coord_values[i, w])
            cutoff_w = cutoff_values[:, w]
            if kernel == "bartlett":
                window *= np.abs(1 - dist / cutoff_w)
            window[dist >= cutoff_w] = 0

        xi_ri = residuals[i:i+1, :] * x_values[i:i+1, :].T
        bartlett_for_i = xi_ri @ ((residuals.T * window[None, :]) @ x_values)

        weight_matrix += bartlett_for_i
        
    return weight_matrix


def _window_estimator_from_scores(data, scores, coordinates, cutoffs, kernel="bartlett"):
    """Construct Conley/HAC meat matrix from per-observation score vectors.

    For an M-estimator with score contributions g_i (Kx1), the Conley/HAC meat is
    \sum_i \sum_j w_ij g_i g_j'. This routine implements the same rectangular-
    cutoff + kernel weighting used elsewhere in this module.

    Args:
        data: (pd.DataFrame) dataframe containing coordinates and cutoffs
        scores: (np.ndarray) n_obs x K matrix of per-observation scores
        coordinates: (list[str]) two coordinate column names
        cutoffs: (list[str]) two cutoff column names
        kernel: (str) 'bartlett' or 'uniform'

    Returns:
        meat: (np.ndarray) K x K weighted outer product of scores
    """
    kernel = str(kernel).strip().lower()
    if kernel not in {"bartlett", "uniform"}:
        raise ValueError(
            f"Unsupported kernel '{kernel}'. Supported values: 'bartlett', 'uniform'.")

    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 2:
        raise ValueError("scores must be a 2D array (n_obs x K)")

    coord_values = data[coordinates].to_numpy(dtype=float, copy=False)
    cutoff_values = np.column_stack([data[c].to_numpy(dtype=float, copy=False) for c in cutoffs])

    n_obs, k = scores.shape
    meat = np.zeros((k, k), dtype=float)

    for i in range(n_obs):
        window = np.ones(n_obs, dtype=float)

        for w in range(coord_values.shape[1]):
            dist = np.abs(coord_values[:, w] - coord_values[i, w])
            cutoff_w = cutoff_values[:, w]
            if kernel == "bartlett":
                window *= np.abs(1 - dist / cutoff_w)
            window[dist >= cutoff_w] = 0

        weighted_sum = (window[:, None] * scores).sum(axis=0)  # (k,)
        meat += np.outer(scores[i, :], weighted_sum)

    return meat


def _compute_ols_conley_se(data, y, x, coordinates, cutoffs, kernel="bartlett"):
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
    pv = np.asarray(ols_model.predict(params), dtype=float).reshape(-1, 1)
    y_col = np.asarray(y, dtype=float).reshape(-1, 1)
    resid = y_col - pv
    # Outer product of score vectors with spatial dependence correction
    filling = _bartlett_window_estimator(data, x, resid, coordinates, cutoffs, kernel=kernel)
    x_values = x.to_numpy(dtype=float, copy=False)
    bread = np.linalg.inv(x_values.T @ x_values)
    
    return filling, bread


def _compute_logit_conley_se(data, y, x, coordinates, cutoffs, kernel="bartlett"):
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
    params = logit_model.fit(disp=0).params
    # Score vectors and hessian
    scores = logit_model.score_obs(params)
    hes = logit_model.hessian(params)
    # Outer product of score vectors with spatial dependence correction
    filling = _window_estimator_from_scores(data, scores, coordinates, cutoffs, kernel=kernel)
    bread = np.linalg.inv(-hes)
    
    return filling, bread


def _compute_probit_conley_se(data, y, x, coordinates, cutoffs, kernel="bartlett"):
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
    params = probit_model.fit(disp=0).params
    # Score vectors and hessian
    scores = probit_model.score_obs(params)
    hes = probit_model.hessian(params)
    # Outer product of score vectors with spatial dependence correction
    filling = _window_estimator_from_scores(data, scores, coordinates, cutoffs, kernel=kernel)
    bread = np.linalg.inv(-hes)
    
    return filling, bread


def _compute_poisson_conley_se(data, y, x, coordinates, cutoffs, kernel="bartlett"):
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
    params = poisson_model.fit(disp=0).params
    # Score vectors and hessian
    scores = poisson_model.score_obs(params)
    hes = poisson_model.hessian(params)
    # Outer product of score vectors with spatial dependence correction
    filling = _window_estimator_from_scores(data, scores, coordinates, cutoffs, kernel=kernel)
    bread = np.linalg.inv(-hes)
    
    return filling, bread


def _compute_nb_conley_se(data, y, x, coordinates, cutoffs, kernel="bartlett"):
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
    res = nb_model.fit(disp=0)
    params = res.params
    # Score vectors and hessian
    scores_full = nb_model.score_obs(params)
    hes_full = nb_model.hessian(params)

    # Compute standard errors for the regression coefficients only (exclude alpha).
    k = x.shape[1]
    scores = scores_full[:, :k]
    hes = hes_full[:k, :k]

    filling = _window_estimator_from_scores(data, scores, coordinates, cutoffs, kernel=kernel)
    bread = np.linalg.inv(-hes)
    
    return filling, bread

