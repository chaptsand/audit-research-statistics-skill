#!/usr/bin/env python3
"""Numeric paired-statistics helper. This does not certify pairing or provenance.

python paired_stats.py --pairs validated_pairs.json
python paired_stats.py --bh validated_family.json
Family JSON: {"family_id": ..., "expected_members": [ids], "p_by_id": {id: p}}
Only stdout is written. Requires numpy and scipy; versions are recorded.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy import stats


def finite_vector(values):
    a = np.asarray(values, dtype=float)
    if a.ndim != 1 or not a.size or not np.isfinite(a).all():
        raise ValueError("Expected a nonempty, finite 1D array")
    return a


def ranked_difference(diff, rank_decimals=None, zero_tolerance=0.0):
    d = finite_vector(diff).copy()
    if rank_decimals is not None:
        if isinstance(rank_decimals, bool) or not isinstance(rank_decimals, int) or rank_decimals < 0:
            raise ValueError("rank_decimals must be a nonnegative integer or null")
        d = np.round(d, rank_decimals)
    if not np.isfinite(zero_tolerance) or zero_tolerance < 0:
        raise ValueError("zero_tolerance must be finite and nonnegative")
    return d[np.abs(d) > zero_tolerance]


def rank_sums(diff, rank_decimals=None, zero_tolerance=0.0):
    d = ranked_difference(diff, rank_decimals, zero_tolerance)
    ranks = stats.rankdata(np.abs(d), method="average")
    wp = float(ranks[d > 0].sum())
    wm = float(ranks[d < 0].sum())
    return d, ranks, wp, wm


def exact_signed_rank(diff, alternative="two-sided", rank_decimals=None, zero_tolerance=0.0):
    if alternative not in ("two-sided", "greater", "less"):
        raise ValueError("Unsupported alternative")
    d, ranks, wp, wm = rank_sums(diff, rank_decimals, zero_tolerance)
    if len(d) > 20:
        raise ValueError("Full enumeration limited to 20 nonzero pairs; choose an explicit alternative method")
    if not len(d):
        return 0.0, 1.0
    # Average ranks are half-integers: integer doubling avoids float-tail errors.
    r2 = np.rint(2 * ranks).astype(np.int64)
    null = np.array([0], dtype=np.int64)
    for r in r2:
        null = np.concatenate((null, null + r))
    obs = int(round(2 * wp))
    if alternative == "greater":
        p = np.mean(null >= obs)
    elif alternative == "less":
        p = np.mean(null <= obs)
    else:
        total = int(r2.sum())
        p = np.mean(np.abs(2 * null - total) >= abs(2 * obs - total))
    return (min(wp, wm) if alternative == "two-sided" else wp), float(p)


def rank_biserial(diff, rank_decimals=None, zero_tolerance=0.0):
    _, _, wp, wm = rank_sums(diff, rank_decimals, zero_tolerance)
    return (wp - wm) / (wp + wm) if wp + wm else 0.0


def hodges_lehmann(diff):
    d = finite_vector(diff)
    i, j = np.triu_indices(len(d))
    return float(np.median((d[i] + d[j]) / 2))


def bootstrap_hl(diff, n_resamples=10000, seed=20260916):
    a = np.asarray(diff, dtype=float)
    if a.ndim not in (1, 2) or not a.size or not np.isfinite(a).all():
        raise ValueError("Expected finite paired differences, 1D units or 2D Run-by-Fold blocks")
    if n_resamples < 10000:
        raise ValueError("At least 10000 bootstrap resamples required")
    rng = np.random.default_rng(seed)
    n_blocks = a.shape[0]
    n_values = a.size
    if n_values > 200:
        raise ValueError("Helper limited to 200 values for Walsh bootstrap; use an approved scalable implementation")
    i, j = np.triu_indices(n_values)
    vals = np.empty(n_resamples)
    for start in range(0, n_resamples, 128):
        stop = min(start + 128, n_resamples)
        ids = rng.integers(0, n_blocks, (stop - start, n_blocks))
        samples = a[ids].reshape(stop - start, n_values)
        vals[start:stop] = np.median((samples[:, i] + samples[:, j]) / 2, axis=1)
    return tuple(float(x) for x in np.quantile(vals, [0.025, 0.975], method="linear"))


def bh_adjust(p_values):
    p = finite_vector(p_values)
    if ((p < 0) | (p > 1)).any():
        raise ValueError("p-values must be in [0,1]")
    order = np.argsort(p, kind="stable")
    scaled = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.clip(adjusted, 0, 1)
    return q


def analyse_pairs(spec):
    left = np.asarray(spec["left"], dtype=float)
    right = np.asarray(spec["right"], dtype=float)
    if left.shape != right.shape or not left.size or not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("Finite, nonempty left/right arrays of identical shape required")
    unit = spec["unit"]
    if unit not in ("run_mean", "fixed_run", "task_mean", "fold_supplement"):
        raise ValueError("Unsupported unit; CV aggregation must be performed and audited before helper input")
    if unit == "fold_supplement":
        if left.ndim != 2:
            raise ValueError("Fold supplement requires real Run-by-Fold 2D blocks")
    elif left.ndim != 1:
        raise ValueError("Primary analyses require 1D paired unit values")
    ids = spec["pair_ids"]
    if len(ids) != left.shape[0] or len(set(ids)) != len(ids) or any(not str(x) for x in ids):
        raise ValueError("Supply unique, nonempty pair_ids for units or fold-supplement Run blocks")
    evidence = spec.get("pairing_evidence", "")
    if not evidence or "REPLACE" in evidence:
        raise ValueError("Supply a pairing evidence reference; actual verification remains external")
    diff = left - right
    d = diff.reshape(-1)
    alt = spec["alternative"]
    if alt not in ("two-sided", "greater", "less"):
        raise ValueError("Unsupported alternative")
    decimals = spec.get("rank_decimals")
    tolerance = spec.get("zero_tolerance", 0.0)
    ranked, _, wp, wm = rank_sums(d, decimals, tolerance)
    p_method = spec["p_method"]
    if unit == "fold_supplement" and p_method != "asymptotic":
        raise ValueError("Fold-supplement protocol requires explicit asymptotic p")
    if p_method == "exact_signflip":
        w, p = exact_signed_rank(d, alt, decimals, tolerance)
        method_label = "conditional exact signflip"
    elif p_method == "asymptotic":
        if len(ranked):
            result = stats.wilcoxon(ranked, alternative=alt, method="approx", zero_method="wilcox", correction=False)
            w, p = float(result.statistic), float(result.pvalue)
        else:
            w, p = 0.0, 1.0
        method_label = "asymptotic; correction=False"
    else:
        raise ValueError("Explicit exact_signflip or asymptotic p_method required")
    ci = bootstrap_hl(diff, spec.get("bootstrap_resamples", 10000), spec.get("seed", 20260916))
    return {
        "comparison_id": spec["comparison_id"],
        "left_name": spec["left_name"], "right_name": spec["right_name"],
        "unit": unit, "pair_ids": ids, "difference_direction": "Left minus Right",
        "n_pairs": int(d.size), "n_nonzero": int(len(ranked)),
        "n_resampling_units": int(left.shape[0]), "alternative": alt,
        "test_method": method_label, "zero_method": "wilcox",
        "rank_decimals": decimals, "zero_tolerance": tolerance,
        "W_plus": wp, "W_minus": wm, "statistic": w, "p": p,
        "r_rb": rank_biserial(d, decimals, tolerance),
        "mean_diff": float(d.mean()), "hl_diff": hodges_lehmann(d),
        "ci_low": ci[0], "ci_high": ci[1],
        "ci_method": "paired percentile bootstrap" if unit != "fold_supplement" else "Run-block percentile bootstrap",
        "ci_level": 0.95, "ci_side": "two-sided", "ci_multiplicity_adjusted": False,
        "bootstrap_resamples": spec.get("bootstrap_resamples", 10000),
        "seed": spec.get("seed", 20260916),
        "pairing_evidence": evidence, "pairing_verified_by_helper": False,
        "differences": diff.tolist(),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "scope": "Numeric calculation only; provenance and inference assumptions require external audit",
    }


def analyse_family(spec):
    ids = spec["expected_members"]
    p_by_id = spec["p_by_id"]
    if not ids or len(ids) != len(set(ids)) or set(ids) != set(p_by_id):
        raise ValueError("Exactly the full declared family, with unique IDs, is required")
    q = bh_adjust([p_by_id[i] for i in ids])
    return {"family_id": spec["family_id"], "m": len(ids), "method": "BH", "q_by_id": dict(zip(ids, q.tolist()))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pairs", type=Path)
    group.add_argument("--bh", type=Path)
    args = parser.parse_args()
    path = args.pairs or args.bh
    spec = json.loads(path.read_text(encoding="utf-8"))
    result = analyse_pairs(spec) if args.pairs else analyse_family(spec)
    result["input_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
