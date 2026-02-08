"""
HC2/HC3 Leverage Corrections for Spatial Conley Standard Errors
===============================================================

RESEARCH QUESTION: Does adjusting scores by leverage (HC2/HC3) improve
finite-sample coverage of Conley spatial standard errors?

WHAT'S THEORETICALLY CLEAN vs. AD HOC:

In the NON-SPATIAL setting (diagonal meat, i=j only):
  HC0:  Var = (X'X)^-1 (sum_i  e_i^2       x_i x_i') (X'X)^-1
  HC1:  Var = n/(n-k) * HC0
  HC2:  Var = (X'X)^-1 (sum_i  e_i^2/(1-h_ii)   x_i x_i') (X'X)^-1
  HC3:  Var = (X'X)^-1 (sum_i  e_i^2/(1-h_ii)^2 x_i x_i') (X'X)^-1

  HC2 is unbiased for sigma^2 because E[e_i^2] = (1-h_ii)*sigma^2.
  HC3 is the jackknife estimator. Both are well-established.

In the SPATIAL Conley setting (cross-products, i != j):
  Omega = sum_i sum_j  w_ij  g_i  g_j'

  For i=j (diagonal): the HC2/HC3 logic applies cleanly.
  For i!=j (off-diagonal): E[e_i * e_j] = -h_ij * sigma^2,
    where h_ij is the OFF-DIAGONAL element of the hat matrix.
    The HC2/HC3 correction using h_ii and h_jj does NOT correctly
    address this bias. It is an AD HOC extension.

  The "naive HC2" approach used here replaces:
    g_i  -->  g_i / sqrt(1 - h_ii)     (HC2)
    g_i  -->  g_i / (1 - h_ii)         (HC3)
  for ALL pairs (i,j), not just i=j.

  This is NOT theoretically justified for the off-diagonal terms.
  Whether it helps coverage is an EMPIRICAL question, which this
  Monte Carlo aims to answer.

For MLE models (logit, probit, poisson, NB):
  The leverage analog is the generalized leverage:
    h_ii = w_i * x_i' (X'WX)^-1 x_i
  where W = diag(Var(y_i|x_i)) is the GLM weight matrix.
  This is even more ad hoc than for OLS.

EXPERIMENT: Run the same Monte Carlo as simulate_finite_sample.ipynb
(200 reps, 10x10 grid, Gaussian copula) and compare coverage of:
  - HC0 (no correction, original code)
  - HC1 (n/(n-k) DoF correction)
  - HC2 (leverage: 1/sqrt(1-h))
  - HC3 (leverage: 1/(1-h))
for Bartlett kernel only (the main kernel of interest).
"""

import math
import warnings
import time
import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm_api
import statsmodels.discrete.discrete_model as sm
import statsmodels.regression.linear_model as lm

warnings.filterwarnings('ignore')

# ── DGP (copied from simulate_finite_sample.ipynb) ──────────────────────

SIDE = 10
REPS = 200
SEED = 12345
BETA_TRUE = 0.5
PHI = 1.0
CUTOFF = 3.0
NB_ALPHA = 0.8


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def make_grid(side):
    c1 = np.repeat(np.arange(1, side + 1), side).astype(float)
    c2 = np.tile(np.arange(1, side + 1), side).astype(float)
    return c1, c2


def pairwise_dist(c1, c2):
    dx = c1[:, None] - c1[None, :]
    dy = c2[:, None] - c2[None, :]
    return np.sqrt(dx * dx + dy * dy)


def gaussian_copula_uniform(c1, c2, phi, rng):
    d = pairwise_dist(c1, c2)
    R = np.exp(-d / phi) + np.eye(len(c1)) * 1e-10
    L = np.linalg.cholesky(R)
    z = L @ rng.standard_normal(len(c1))
    u = st.norm.cdf(z)
    return np.clip(u, 1e-12, 1.0 - 1e-12)


def build_df(y, x, c1, c2):
    return pd.DataFrame({
        'y': y, 'indep1': x, 'const': 1.0,
        'C1': c1, 'C2': c2,
        'cutoff1': CUTOFF, 'cutoff2': CUTOFF,
    })


# ── Spatial weight matrix (Bartlett kernel, vectorized) ──────────────────

def compute_W_bartlett(c1, c2, cutoff):
    """Full (n,n) Bartlett weight matrix."""
    dist1 = np.abs(c1[:, None] - c1[None, :])
    dist2 = np.abs(c2[:, None] - c2[None, :])
    W = np.maximum(0.0, 1.0 - dist1 / cutoff) * np.maximum(0.0, 1.0 - dist2 / cutoff)
    return W


# ── Generalized leverage for each model ─────────────────────────────────

def glm_weights(model_name, X, params):
    """GLM working weights: w_i = Var(y_i | x_i) for each model."""
    eta = X @ params
    if model_name == 'logit':
        p = sigmoid(eta)
        return p * (1 - p)
    elif model_name == 'probit':
        p = st.norm.cdf(eta)
        phi = st.norm.pdf(eta)
        # Weight for probit IRLS: (phi^2) / (p*(1-p))
        p = np.clip(p, 1e-12, 1 - 1e-12)
        return phi**2 / (p * (1 - p))
    elif model_name == 'poisson':
        return np.exp(eta)
    elif model_name == 'nb':
        mu = np.exp(eta)
        return mu / (1 + NB_ALPHA * mu)
    else:
        raise ValueError(model_name)


def compute_leverage(X, weights):
    """Generalized leverage: h_ii = w_i * x_i' (X'WX)^-1 x_i."""
    W = np.diag(weights)
    XtWX_inv = np.linalg.inv(X.T @ W @ X)
    # Vectorized: h_i = w_i * sum_j (X[i,j] * (XtWX_inv @ X[i])_j)
    h = weights * np.sum(X @ XtWX_inv * X, axis=1)
    return np.clip(h, 0, 0.999)  # clip to avoid division by zero


# ── Compute Conley SE with HC0/HC1/HC2/HC3 ──────────────────────────────

def conley_se_all_variants(model_name, df):
    """Compute Conley SEs (Bartlett) with HC0, HC1, HC2, HC3.

    Returns dict with keys: se_hc0, se_hc1, se_hc2, se_hc3
    (each is the SE for indep1, the first regressor).
    """
    y = df['y'].to_numpy(dtype=float)
    X = df[['indep1', 'const']].to_numpy(dtype=float)
    c1 = df['C1'].to_numpy(dtype=float)
    c2 = df['C2'].to_numpy(dtype=float)
    n, k = X.shape

    # 1. Fit model and get scores + bread
    if model_name == 'logit':
        mod = sm.Logit(y, X)
        res = mod.fit(disp=0)
        scores = mod.score_obs(res.params)
        bread = np.linalg.inv(-mod.hessian(res.params))
    elif model_name == 'probit':
        mod = sm.Probit(y, X)
        res = mod.fit(disp=0)
        scores = mod.score_obs(res.params)
        bread = np.linalg.inv(-mod.hessian(res.params))
    elif model_name == 'poisson':
        mod = sm.Poisson(y, X)
        res = mod.fit(disp=0)
        scores = mod.score_obs(res.params)
        bread = np.linalg.inv(-mod.hessian(res.params))
    elif model_name == 'nb':
        mod = sm.NegativeBinomial(y, X)
        res = mod.fit(disp=0)
        scores_full = mod.score_obs(res.params)
        hes_full = mod.hessian(res.params)
        scores = scores_full[:, :k]
        bread = np.linalg.inv(-hes_full[:k, :k])
    else:
        raise ValueError(model_name)

    # 2. Compute Bartlett weight matrix
    W = compute_W_bartlett(c1, c2, CUTOFF)

    # 3. Compute generalized leverage
    w_glm = glm_weights(model_name, X, res.params[:k])
    h = compute_leverage(X, w_glm)

    # 4. Compute meat for each HC variant
    # HC0: raw scores
    meat_hc0 = scores.T @ W @ scores

    # HC1: same meat, just scale variance by n/(n-k) later
    # (meat is the same as HC0)

    # HC2: scores_adj = scores / sqrt(1 - h_ii)
    adj2 = 1.0 / np.sqrt(1.0 - h)
    scores_hc2 = scores * adj2[:, None]
    meat_hc2 = scores_hc2.T @ W @ scores_hc2

    # HC3: scores_adj = scores / (1 - h_ii)
    adj3 = 1.0 / (1.0 - h)
    scores_hc3 = scores * adj3[:, None]
    meat_hc3 = scores_hc3.T @ W @ scores_hc3

    # 5. Symmetrize and sandwich
    results = {}
    for label, meat, dof_factor in [
        ('hc0', meat_hc0, 1.0),
        ('hc1', meat_hc0, n / (n - k)),     # same meat, scaled
        ('hc2', meat_hc2, 1.0),
        ('hc3', meat_hc3, 1.0),
        ('hc2_plus_dof', meat_hc2, n / (n - k)),  # HC2 + DoF combined
    ]:
        meat_sym = 0.5 * (meat + meat.T)
        bread_sym = 0.5 * (bread + bread.T)
        var = dof_factor * bread_sym @ meat_sym @ bread_sym.T
        se = np.sqrt(np.clip(np.diagonal(var), 0.0, None))
        results[f'se_{label}'] = float(se[0])  # SE for indep1

    return results


# ── Monte Carlo ──────────────────────────────────────────────────────────

def generate_data(model_name, rng):
    c1, c2 = make_grid(SIDE)
    x = rng.standard_normal(len(c1))
    u = gaussian_copula_uniform(c1, c2, PHI, rng)
    eta = BETA_TRUE * x

    if model_name == 'logit':
        p = np.clip(sigmoid(eta), 1e-12, 1 - 1e-12)
        y = (u < p).astype(float)
    elif model_name == 'probit':
        p = np.clip(st.norm.cdf(eta), 1e-12, 1 - 1e-12)
        y = (u < p).astype(float)
    elif model_name == 'poisson':
        mu = np.clip(np.exp(eta), 1e-12, 1e6)
        y = st.poisson(mu).ppf(u).astype(float)
    elif model_name == 'nb':
        mu = np.clip(np.exp(eta), 1e-12, 1e6)
        n_param = 1.0 / NB_ALPHA
        p_param = n_param / (n_param + mu)
        y = st.nbinom(n_param, p_param).ppf(u).astype(float)

    return build_df(y, x, c1, c2)


def simulate_one(model_name, rng):
    df = generate_data(model_name, rng)
    ses = conley_se_all_variants(model_name, df)

    # Also get beta_hat
    y = df['y'].to_numpy()
    X = df[['indep1', 'const']].to_numpy()
    if model_name == 'logit':
        beta_hat = float(sm.Logit(y, X).fit(disp=0).params[0])
    elif model_name == 'probit':
        beta_hat = float(sm.Probit(y, X).fit(disp=0).params[0])
    elif model_name == 'poisson':
        beta_hat = float(sm.Poisson(y, X).fit(disp=0).params[0])
    elif model_name == 'nb':
        beta_hat = float(sm.NegativeBinomial(y, X).fit(disp=0).params[0])

    # Coverage
    z = 1.959963984540054
    row = {'model': model_name, 'beta_hat': beta_hat}
    for key, se in ses.items():
        row[key] = se
        covers = (beta_hat - z * se <= BETA_TRUE) and (BETA_TRUE <= beta_hat + z * se)
        row[key.replace('se_', 'cover_')] = float(covers)

    return row


def run_experiment():
    models = ['logit', 'probit', 'poisson', 'nb']

    print("=" * 80)
    print("HC2/HC3 LEVERAGE CORRECTIONS FOR CONLEY SEs — MONTE CARLO EXPERIMENT")
    print("=" * 80)
    print()
    print("Setup: %d reps, %dx%d grid (n=%d), beta_true=%.1f, phi=%.1f, cutoff=%.1f"
          % (REPS, SIDE, SIDE, SIDE*SIDE, BETA_TRUE, PHI, CUTOFF))
    print("Kernel: Bartlett only")
    print()

    all_results = []

    for model_name in models:
        rng = np.random.default_rng(SEED)
        rows = []
        failures = 0
        t0 = time.perf_counter()

        for rep in range(REPS):
            try:
                rows.append(simulate_one(model_name, rng))
            except Exception as e:
                failures += 1
                continue

        elapsed = time.perf_counter() - t0
        df = pd.DataFrame(rows)
        print("%s: %d/%d reps (%.1fs, %d failures)" %
              (model_name, len(df), REPS, elapsed, failures))

        summary = {
            'model': model_name,
            'mean_beta': df['beta_hat'].mean(),
            'emp_sd': df['beta_hat'].std(ddof=1),
        }
        for variant in ['hc0', 'hc1', 'hc2', 'hc3', 'hc2_plus_dof']:
            summary['mean_se_%s' % variant] = df['se_%s' % variant].mean()
            summary['coverage_%s' % variant] = df['cover_%s' % variant].mean()

        all_results.append(summary)

    # Print results table
    print()
    print()
    print("RESULTS: 95%% COVERAGE RATES (Bartlett kernel)")
    print("=" * 80)
    print()
    print("%-8s | %5s | %5s | %5s | %5s | %5s | %7s" %
          ("Model", "HC0", "HC1", "HC2", "HC3", "HC2+df", "Emp.SD"))
    print("-" * 60)
    for r in all_results:
        print("%-8s | %.3f | %.3f | %.3f | %.3f | %.3f  | %.4f" % (
            r['model'],
            r['coverage_hc0'],
            r['coverage_hc1'],
            r['coverage_hc2'],
            r['coverage_hc3'],
            r['coverage_hc2_plus_dof'],
            r['emp_sd'],
        ))
    print()

    print()
    print("RESULTS: MEAN STANDARD ERRORS")
    print("=" * 80)
    print()
    print("%-8s | %8s | %8s | %8s | %8s | %8s | %8s" %
          ("Model", "Emp.SD", "HC0", "HC1", "HC2", "HC3", "HC2+df"))
    print("-" * 75)
    for r in all_results:
        print("%-8s | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f" % (
            r['model'],
            r['emp_sd'],
            r['mean_se_hc0'],
            r['mean_se_hc1'],
            r['mean_se_hc2'],
            r['mean_se_hc3'],
            r['mean_se_hc2_plus_dof'],
        ))
    print()

    print()
    print("INTERPRETATION")
    print("=" * 80)
    print()
    print("Target: 95%% nominal coverage. Closer to 0.950 is better.")
    print()
    print("HC0 = raw Conley (no correction)")
    print("HC1 = n/(n-k) DoF correction (the scalar fix)")
    print("HC2 = leverage: scores / sqrt(1 - h_ii)  [ad hoc for spatial]")
    print("HC3 = leverage: scores / (1 - h_ii)       [ad hoc for spatial]")
    print("HC2+df = HC2 combined with n/(n-k) DoF correction")
    print()
    print("IMPORTANT CAVEAT: HC2/HC3 are theoretically justified only for")
    print("the diagonal (i=j) of the meat matrix. For off-diagonal terms")
    print("(i!=j, which dominate in the spatial setting), the correction")
    print("is AD HOC. The bias of E[e_i * e_j] involves h_ij (off-diagonal")
    print("hat matrix), not h_ii or h_jj individually.")
    print()
    print("For MLE models, the 'leverage' is the generalized leverage from")
    print("the GLM working weights, which is itself an approximation.")

    return all_results


if __name__ == '__main__':
    run_experiment()
