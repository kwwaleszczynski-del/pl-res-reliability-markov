"""Validation of the Markov models against the series they were fitted on and held-out years.

* state residence times (sojourns): empirical vs geometric 1/(1 - P_ii) vs simulated;
* autocorrelation of the normalised output, monthly means (seasonality);
* chronological cross-validation (train 2015-2022, test 2023-2024) and leave-one-year-out:
  out-of-sample log-likelihood per transition of i.i.d., homogeneous and month-dependent
  models, and the Anderson-Goodman likelihood-ratio test that the test-period transition
  probabilities equal the training ones. Unlike the chi-square / KS tests on state
  frequencies used in the thesis, these work on transition counts, so serial dependence of
  hourly data does not invalidate them.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from res_reliability import markov as mk


def sojourns(states: np.ndarray, n_states: int) -> list[np.ndarray]:
    """Lengths of uninterrupted runs in each state (runs touching an invalid value are dropped)."""
    s = np.asarray(states)
    change = np.flatnonzero(np.diff(s) != 0) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(s)]
    out = [[] for _ in range(n_states)]
    for a, b in zip(starts, ends):
        k = s[a]
        if k < 0 or a == 0 or b == len(s):
            continue
        out[k].append(b - a)
    return [np.asarray(v) for v in out]


def acf(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Autocorrelation of a (possibly 2-D, one row per year) series, pooled over rows."""
    x = np.atleast_2d(np.asarray(x, dtype=float))
    x = x - np.nanmean(x)
    var = np.nanmean(x * x)
    out = np.empty(max_lag + 1)
    for k in range(max_lag + 1):
        a, b = x[:, : x.shape[1] - k], x[:, k:]
        out[k] = np.nanmean(a * b) / var
    return out


def log_likelihood(states: np.ndarray, P: np.ndarray, months: np.ndarray | None = None,
                   valid: np.ndarray | None = None, floor: float = 1e-6) -> tuple[float, int]:
    """Sum of log P(s_t+1 | s_t) over consecutive valid observations; P may be (12, S, S)."""
    s = np.asarray(states)
    keep = (s >= 0) if valid is None else (np.asarray(valid, bool) & (s >= 0))
    idx = np.flatnonzero(keep)
    a, b = s[idx[:-1]], s[idx[1:]]
    if P.ndim == 3:
        m = np.asarray(months)[idx[:-1]] - 1
        p = P[m, a, b]
    else:
        p = P[a, b]
    return float(np.log(np.maximum(p, floor)).sum()), int(len(a))


def anderson_goodman(test_counts: np.ndarray, P_train: np.ndarray) -> dict:
    """LR test H0: test transitions follow P_train. Returns statistic, dof, p-value."""
    n = np.asarray(test_counts, dtype=float)
    rows = n.sum(axis=1)
    P_test = mk.normalise_rows(n)
    mask = (n > 0) & (P_train > 0)
    g2 = 2.0 * float((n[mask] * np.log(P_test[mask] / P_train[mask])).sum())
    # dof: free parameters of the rows that were observed, over cells allowed by P_train
    dof = int(sum(max((P_train[i] > 0).sum() - 1, 0) for i in range(len(rows)) if rows[i] > 0))
    impossible = int(((n > 0) & (P_train <= 0)).sum())
    return {"g2": g2, "dof": dof, "p_value": float(stats.chi2.sf(g2, dof)) if dof else float("nan"),
            "transitions": int(rows.sum()), "impossible_transitions": impossible,
            "max_abs_diff": float(np.abs(P_test - P_train)[rows > 0].max())}


def cross_validate(
    x: np.ndarray,
    months_year: np.ndarray,
    train_rows: np.ndarray,
    test_rows: np.ndarray,
    n_states: int = 5,
    valid: np.ndarray | None = None,
) -> dict:
    """Fit on ``train_rows`` (years) of a (years, hours) matrix, evaluate on ``test_rows``."""
    x = np.asarray(x, dtype=float)
    n_years = x.shape[0]
    months = np.tile(months_year, n_years).reshape(n_years, -1)
    v = np.ones_like(x, dtype=bool) if valid is None else np.asarray(valid, bool)

    def flat(rows, arr):
        return arr[rows].reshape(-1)

    hom = mk.fit(flat(train_rows, x), n_states, valid=flat(train_rows, v))
    sea = mk.fit(flat(train_rows, x), n_states, months=flat(train_rows, months), valid=flat(train_rows, v))
    s_test = mk.discretize(np.where(flat(test_rows, v), flat(test_rows, x), np.nan), n_states)
    m_test = flat(test_rows, months)
    v_test = flat(test_rows, v)
    pi = mk.stationary(hom.P)
    iid = np.tile(pi, (n_states, 1))
    out = {}
    for name, P, mm in [("iid", iid, None), ("homogeneous", hom.P, None), ("seasonal", sea.P, m_test)]:
        ll, n = log_likelihood(s_test, P, mm, v_test)
        out[f"ll_{name}"] = ll / n
    counts = mk.transition_counts(s_test, n_states, v_test & (s_test >= 0))
    test_fit = mk.normalise_rows(counts)
    ll_oracle, n = log_likelihood(s_test, test_fit, None, v_test)
    out["ll_test_fit"] = ll_oracle / n
    out["n_transitions"] = n
    out.update({f"ag_{k}": val for k, val in anderson_goodman(counts, hom.P).items()})
    freq = np.bincount(s_test[v_test & (s_test >= 0)], minlength=n_states)
    freq = freq / freq.sum()
    out["pi_train"] = pi.tolist()
    out["freq_test"] = freq.tolist()
    out["pi_max_abs_diff"] = float(np.abs(pi - freq).max())
    out["mean_train_model"] = float(hom.mean())
    out["mean_test_data"] = float(np.nanmean(np.where(flat(test_rows, v), flat(test_rows, x), np.nan)))
    return out


def low_runs(daily: np.ndarray, threshold: float) -> np.ndarray:
    """Lengths (days) of runs of consecutive days below ``threshold`` in every row of a (years, days) matrix.

    Unlike ``sojourns``, runs that touch the start or the end of a year are kept (a year is a closed window).
    """
    out = []
    for row in np.atleast_2d(np.asarray(daily, dtype=float)):
        low = np.r_[False, row < threshold, False].astype(np.int8)
        d = np.diff(low)
        out.append(np.flatnonzero(d == -1) - np.flatnonzero(d == 1))
    return np.concatenate(out) if out else np.array([], dtype=np.int64)


def run_length_counts(lengths: np.ndarray, n_years: float, max_len: int = 10) -> np.ndarray:
    """Runs per year of length 1, 2, ..., max_len - 1 and >= max_len days."""
    lengths = np.asarray(lengths, dtype=np.int64)
    c = np.bincount(np.minimum(lengths, max_len), minlength=max_len + 1)[1:]
    return c / n_years


def rolling_min_mean(daily: np.ndarray, window: int) -> np.ndarray:
    """Lowest ``window``-day mean in every row (year) of a (years, days) matrix."""
    x = np.atleast_2d(np.asarray(daily, dtype=float))
    c = np.cumsum(np.c_[np.zeros(x.shape[0]), x], axis=1)
    return ((c[:, window:] - c[:, :-window]) / window).min(axis=1)
