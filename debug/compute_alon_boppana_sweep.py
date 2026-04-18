"""
Alon-Boppana sweep: extend Sprint 1 Ihara zeta / Ramanujan computation to
larger graph sizes, then report the growth of the sub-Ramanujan deviation.

Sprint 2 driver (Track RH-D), answers Paper 29 §6.1 open question.

Graphs computed here (new sizes — Sprint 1 had max_n=2,3 / N_max=2,3):

  - S^3 Coulomb at max_n = 4 (V = 30) and max_n = 5 (V = 55).
  - S^5 Bargmann-Segal at N_max = 4 (V = 35) and N_max = 5 (V = 56).
  - Dirac-S^3 rule A at n_max = 4 (V = 60, 2E = 144).
  - Dirac-S^3 rule B at n_max = 4 (V = 60, 2E = 624) — large but tractable.

Reuses geovac/ihara_zeta.py and geovac/ihara_zeta_dirac.py without
modification; the driver merely iterates.  The previous Sprint 1 rows
at smaller max_n / N_max are copied verbatim from their source JSONs so
that the Sprint 2 memo and fit can operate on the full deviation series
(max_n = 3, 4, 5 for S^3; N_max = 2, 3, 4, 5 for S^5; Dirac n_max
= 2, 3, 4 for both rules).

Outputs:
  debug/data/alon_boppana_sweep.json
  debug/alon_boppana_memo.md (written separately)
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import sympy as sp

from geovac.ihara_zeta import (
    _count_components,
    functional_equation_report,
    hashimoto_matrix,
    ihara_zeta_bass,
    is_ramanujan,
    zeta_zeros_from_hashimoto,
)
from geovac.lattice import GeometricLattice
from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.ihara_zeta_dirac import (
    build_dirac_s3_graph,
    describe_adjacency_rule,
)


ROOT = Path(__file__).parent
DATA_DIR = ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------------
# Core per-graph computation.  Only numerical invariants — no symbolic
# closed forms (factoring the Bass determinant at V > 30 is slow and not
# needed for the Alon-Boppana fit).
# ----------------------------------------------------------------------------

def _compute_numerical_record(label: str,
                              A_pattern: np.ndarray,
                              n_label: int) -> Dict:
    """Compute numerical Ihara invariants (no symbolic closed form).

    Parameters
    ----------
    label
        Graph family label: "S3_Coulomb", "S5_Bargmann", "Dirac_A", "Dirac_B".
    A_pattern
        V x V, 0/1 symmetric adjacency.
    n_label
        Size parameter for this graph family (max_n or N_max).

    Returns
    -------
    dict of JSON-serializable numerical data.
    """
    A = (A_pattern != 0).astype(int)
    V = int(A.shape[0])
    E = int(A.sum() // 2)
    c = _count_components(A)
    r = E - V + c

    deg = A.sum(axis=1).astype(int)
    q_max = int(deg.max()) - 1 if V > 0 else 0
    q_min = int(deg.min()) - 1 if V > 0 else 0
    regular = bool(deg.min() == deg.max())

    rec = {
        "label": label,
        "n_label": n_label,
        "V": V,
        "E": E,
        "c": c,
        "beta_1": r,
        "q_max": q_max,
        "q_min": q_min,
        "regular": regular,
        "max_degree": int(deg.max()) if V > 0 else 0,
    }

    if E == 0:
        rec["rho_T"] = 0.0
        rec["max_abs_nontrivial"] = 0.0
        rec["sqrt_q_max"] = 0.0
        rec["deviation"] = -1.0  # trivial forest — use same sentinel as Sprint 1
        rec["is_ramanujan"] = True
        return rec

    # Hashimoto spectrum
    t0 = time.time()
    T = hashimoto_matrix(A).astype(float)
    # numpy numeric eigensolve
    ev = np.linalg.eigvals(T)
    t_hashi = time.time() - t0
    mags = np.abs(ev)
    rho = float(mags.max())
    nontriv_mask = mags < rho - 1e-9
    mu_nt = float(mags[nontriv_mask].max()) if nontriv_mask.any() else 0.0
    sqrt_q_max = float(np.sqrt(max(q_max, 0)))
    dev = mu_nt - sqrt_q_max
    is_ram = dev <= 1e-9

    rec.update({
        "rho_T": rho,
        "max_abs_nontrivial": mu_nt,
        "sqrt_q_max": sqrt_q_max,
        "deviation": float(dev),
        "is_ramanujan": bool(is_ram),
        "hashimoto_dim": int(T.shape[0]),
        "hashimoto_eigensolve_seconds": float(t_hashi),
    })
    return rec


# ----------------------------------------------------------------------------
# Sprint 1 data copy-in (so the Alon-Boppana fit operates on the full series)
# ----------------------------------------------------------------------------

SPRINT1_S3_JSON = DATA_DIR / "ihara_zeta_geovac_hopf.json"
SPRINT1_DIRAC_JSON = DATA_DIR / "ihara_zeta_dirac_s3.json"


def _copy_sprint1_S3_S5() -> List[Dict]:
    """Read Sprint 1 data for S^3 (max_n=2,3) and S^5 (N_max=2,3).

    Returns a list of normalized records matching the Sprint 2 schema.
    """
    with SPRINT1_S3_JSON.open("r") as f:
        sp1 = json.load(f)

    rows: List[Dict] = []

    for mn in (2, 3):
        key = f"S3_Coulomb_max_n_{mn}"
        if key not in sp1:
            continue
        g = sp1[key]
        rec = {
            "label": "S3_Coulomb",
            "n_label": mn,
            "V": int(g["V"]),
            "E": int(g["E"]),
            "c": int(g["c"]),
            "beta_1": int(g["beta_1"]),
            "q_max": int(g.get("critical_radii", {}).get("q_max", 0)),
            "q_min": int(g.get("critical_radii", {}).get("q_min", 0)),
            "regular": bool(g["regular"]),
            "max_degree": int(g["degree_sequence"]["max"]),
            "rho_T": float(g.get("hashimoto_spectral_radius", 0.0)),
            "max_abs_nontrivial": float(g["ramanujan"]["deviation"]) +
                                  _parse_sqrt_q_max(g["ramanujan"]["explanation"]),
            "sqrt_q_max": _parse_sqrt_q_max(g["ramanujan"]["explanation"]),
            "deviation": float(g["ramanujan"]["deviation"]),
            "is_ramanujan": bool(g["ramanujan"]["is_ramanujan"]),
            "source": "sprint1",
        }
        rows.append(rec)

    for nm in (2, 3):
        key = f"S5_Bargmann_Segal_N_max_{nm}"
        if key not in sp1:
            continue
        g = sp1[key]
        rec = {
            "label": "S5_Bargmann",
            "n_label": nm,
            "V": int(g["V"]),
            "E": int(g["E"]),
            "c": int(g["c"]),
            "beta_1": int(g["beta_1"]),
            "q_max": int(g.get("critical_radii", {}).get("q_max", 0)),
            "q_min": int(g.get("critical_radii", {}).get("q_min", 0)),
            "regular": bool(g["regular"]),
            "max_degree": int(g["degree_sequence"]["max"]),
            "rho_T": float(g.get("hashimoto_spectral_radius", 0.0)),
            "max_abs_nontrivial": float(g["ramanujan"]["deviation"]) +
                                  _parse_sqrt_q_max(g["ramanujan"]["explanation"]),
            "sqrt_q_max": _parse_sqrt_q_max(g["ramanujan"]["explanation"]),
            "deviation": float(g["ramanujan"]["deviation"]),
            "is_ramanujan": bool(g["ramanujan"]["is_ramanujan"]),
            "source": "sprint1",
        }
        rows.append(rec)

    return rows


def _parse_sqrt_q_max(explanation: str) -> float:
    """Extract sqrt(q_max)=X from the is_ramanujan explanation string."""
    marker = "sqrt(q_max)="
    i = explanation.find(marker)
    if i < 0:
        return 0.0
    j = explanation.find(",", i)
    return float(explanation[i + len(marker): j])


def _copy_sprint1_dirac() -> List[Dict]:
    """Read Sprint 1 Dirac data for rule A/B at n_max=2,3."""
    with SPRINT1_DIRAC_JSON.open("r") as f:
        sp1 = json.load(f)

    rows: List[Dict] = []
    for row in sp1.get("rows", []):
        if "error" in row:
            continue
        if row["n_max"] not in (2, 3):
            # n_max=1 is trivial (2 nodes, no cycles) and already in Sprint 1
            # if we want the full history.  Keep it for completeness.
            if row["n_max"] != 1:
                continue

        rec = {
            "label": f"Dirac_{row['adjacency_rule']}",
            "n_label": int(row["n_max"]),
            "V": int(row["V"]),
            "E": int(row["E"]),
            "c": int(row["c"]),
            "beta_1": int(row["beta1"]),
            "q_max": int(row["q_max"]),
            "q_min": int(min(row["degree_sequence"])) - 1 if row["degree_sequence"] else 0,
            "regular": bool(min(row["degree_sequence"]) == max(row["degree_sequence"]))
                       if row["degree_sequence"] else False,
            "max_degree": int(max(row["degree_sequence"])) if row["degree_sequence"] else 0,
            "rho_T": float(row["spectral_radius_rho"]),
            "max_abs_nontrivial": float(row["max_abs_nontrivial"]),
            "sqrt_q_max": float(row["sqrt_q_max"]),
            "deviation": float(row["ramanujan_deviation"]),
            "is_ramanujan": bool(row["ramanujan_verdict"]),
            "source": "sprint1",
        }
        rows.append(rec)
    return rows


# ----------------------------------------------------------------------------
# Alon-Boppana fit
# ----------------------------------------------------------------------------

def _fit_deviation_models(sizes: List[int],
                          devs: List[float]) -> Dict[str, Dict]:
    """Fit deviation against candidate functional forms and report R^2.

    Two fit targets:
      (A) |deviation|  (magnitude; classical Alon-Boppana framing — the
          sub-Ramanujan gap should approach zero from below for a
          Ramanujan-sequence).
      (B) deviation    (signed; tests whether the graphs stay Ramanujan
          asymptotically OR cross the bound and become non-Ramanujan).

    Six candidate features tested for each target:
      - 1/V
      - 1/sqrt(V)
      - 1/log(V)
      - log(V)/sqrt(V)
      - sqrt(V)        (grows — captures non-Ramanujan departure)
      - V              (grows linearly — captures strong departure)

    Plus power-law log-log fits |dev| = a * V^b for the magnitude.

    Returns {model_name: {...}} nested per target.
    """
    if len(sizes) < 3:
        return {"insufficient_data": True, "n_points": len(sizes)}

    sizes_arr = np.asarray(sizes, dtype=float)
    devs_arr = np.asarray(devs, dtype=float)

    def _lin_fit(x: np.ndarray, y: np.ndarray) -> Dict:
        # Fit y = a + b*x via least squares; return a, b, R^2.
        xm = x.mean()
        ym = y.mean()
        sx = ((x - xm) ** 2).sum()
        sxy = ((x - xm) * (y - ym)).sum()
        b = sxy / sx if sx != 0 else 0.0
        a = ym - b * xm
        y_pred = a + b * x
        ss_res = ((y - y_pred) ** 2).sum()
        ss_tot = ((y - ym) ** 2).sum()
        r2 = 1.0 - ss_res / ss_tot if ss_tot != 0 else None
        return {
            "a": float(a),
            "b": float(b),
            "r_squared": float(r2) if r2 is not None else None,
            "predicted": [float(z) for z in y_pred],
            "feature": [float(z) for z in x],
        }

    features = {
        "inv_V":        1.0 / sizes_arr,
        "inv_sqrt_V":   1.0 / np.sqrt(sizes_arr),
        "inv_log_V":    1.0 / np.log(sizes_arr),
        "logV_over_sqrtV": np.log(sizes_arr) / np.sqrt(sizes_arr),
        "sqrt_V":       np.sqrt(sizes_arr),
        "V":            sizes_arr,
    }

    models_abs = {name: _lin_fit(x, np.abs(devs_arr)) for name, x in features.items()}
    models_signed = {name: _lin_fit(x, devs_arr) for name, x in features.items()}

    # Power-law fit on magnitude, |dev| = a * V^b -> log|dev| = log(a) + b*logV
    pos_mask = np.abs(devs_arr) > 1e-12
    if pos_mask.sum() >= 3:
        logV = np.log(sizes_arr[pos_mask])
        log_abs_dev = np.log(np.abs(devs_arr[pos_mask]))
        power_fit = _lin_fit(logV, log_abs_dev)
        power_fit["exponent_b"] = power_fit["b"]
        power_fit["prefactor_a"] = float(np.exp(power_fit["a"]))
    else:
        power_fit = {"insufficient_data": True}

    best_abs = max(models_abs.items(),
                   key=lambda kv: (kv[1]["r_squared"] or -1.0))
    best_signed = max(models_signed.items(),
                      key=lambda kv: (kv[1]["r_squared"] or -1.0))

    return {
        "n_points": len(sizes),
        "sizes": [float(z) for z in sizes_arr],
        "deviations": [float(z) for z in devs_arr],
        "abs_deviations": [float(z) for z in np.abs(devs_arr)],
        "models_abs_deviation": models_abs,
        "models_signed_deviation": models_signed,
        "power_law_log_log": power_fit,
        "best_model_abs": best_abs[0],
        "best_model_signed": best_signed[0],
    }


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("Sprint 2 (Track RH-D): Alon-Boppana sweep")
    print("=" * 78)

    rows: List[Dict] = []

    # --- Sprint 2 new computations ---

    # S^3 Coulomb at max_n = 4, 5
    for mn in (4, 5):
        t0 = time.time()
        L = GeometricLattice(max_n=mn)
        A = L.adjacency.toarray()
        rec = _compute_numerical_record("S3_Coulomb", A, mn)
        rec["source"] = "sprint2"
        rec["wall_seconds"] = time.time() - t0
        print(f"S3_Coulomb max_n={mn}: V={rec['V']} 2E={rec['hashimoto_dim']} "
              f"dev={rec['deviation']:+.4f} Ram={rec['is_ramanujan']} "
              f"({rec['wall_seconds']:.1f}s)")
        rows.append(rec)

    # S^5 Bargmann-Segal at N_max = 4, 5
    for nm in (4, 5):
        t0 = time.time()
        g = build_bargmann_graph(nm)
        A = g.adjacency_dense()
        rec = _compute_numerical_record("S5_Bargmann", A, nm)
        rec["source"] = "sprint2"
        rec["wall_seconds"] = time.time() - t0
        print(f"S5_Bargmann N_max={nm}: V={rec['V']} 2E={rec['hashimoto_dim']} "
              f"dev={rec['deviation']:+.4f} Ram={rec['is_ramanujan']} "
              f"({rec['wall_seconds']:.1f}s)")
        rows.append(rec)

    # Dirac-S^3 rule A at n_max = 4
    try:
        t0 = time.time()
        A, labels, deg, desc = build_dirac_s3_graph(4, "A")
        rec = _compute_numerical_record("Dirac_A", A, 4)
        rec["source"] = "sprint2"
        rec["wall_seconds"] = time.time() - t0
        print(f"Dirac_A n_max=4: V={rec['V']} 2E={rec['hashimoto_dim']} "
              f"dev={rec['deviation']:+.4f} Ram={rec['is_ramanujan']} "
              f"({rec['wall_seconds']:.1f}s)")
        rows.append(rec)
    except Exception as e:
        print(f"Dirac_A n_max=4: ERROR {type(e).__name__}: {e}")

    # Dirac-S^3 rule B at n_max = 4 (large: 2E = 624)
    try:
        t0 = time.time()
        A, labels, deg, desc = build_dirac_s3_graph(4, "B")
        rec = _compute_numerical_record("Dirac_B", A, 4)
        rec["source"] = "sprint2"
        rec["wall_seconds"] = time.time() - t0
        print(f"Dirac_B n_max=4: V={rec['V']} 2E={rec['hashimoto_dim']} "
              f"dev={rec['deviation']:+.4f} Ram={rec['is_ramanujan']} "
              f"({rec['wall_seconds']:.1f}s)")
        rows.append(rec)
    except Exception as e:
        print(f"Dirac_B n_max=4: ERROR {type(e).__name__}: {e}")

    # --- Copy in Sprint 1 rows for the full deviation series ---
    sp1_scalar = _copy_sprint1_S3_S5()
    sp1_dirac = _copy_sprint1_dirac()
    rows.extend(sp1_scalar)
    rows.extend(sp1_dirac)

    # --- Alon-Boppana fit per graph family ---

    fits: Dict[str, Dict] = {}
    for family in ("S3_Coulomb", "S5_Bargmann", "Dirac_A", "Dirac_B"):
        # Include ALL points (Ramanujan AND non-Ramanujan), but exclude
        # trivial forests where beta_1 = 0 (no cycles, deviation is
        # meaningless) and exclude the sentinel q_max = 0 rows.
        fam_rows = [r for r in rows
                    if r["label"] == family
                    and r.get("q_max", 0) > 0
                    and r.get("beta_1", 0) > 0
                    and abs(r.get("deviation", 0.0)) < 5.0]  # sanity
        fam_rows.sort(key=lambda r: r["V"])
        sizes = [r["V"] for r in fam_rows]
        devs = [r["deviation"] for r in fam_rows]

        fits[family] = {
            "family": family,
            "points": [{"n_label": r["n_label"], "V": r["V"],
                        "E": r["E"], "deviation": r["deviation"]}
                       for r in fam_rows],
            "fit": _fit_deviation_models(sizes, devs),
        }

    # --- Write JSON ---

    out = {
        "sprint": "RH-D",
        "description": (
            "Alon-Boppana sweep: Sprint 1 (max_n,N_max = 2..3) + Sprint 2 "
            "(max_n,N_max = 4..5; Dirac n_max = 4), fitted to four candidate "
            "asymptotic forms."
        ),
        "rows": rows,
        "fits": fits,
    }
    out_path = DATA_DIR / "alon_boppana_sweep.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {out_path}")

    # --- Terse fit summary ---
    print("\nFit summary:")
    for family, fit_obj in fits.items():
        fit = fit_obj["fit"]
        if "insufficient_data" in fit:
            print(f"  {family}: only {fit['n_points']} points, no fit")
            continue
        best_abs = fit["best_model_abs"]
        r2_abs = fit["models_abs_deviation"][best_abs]["r_squared"]
        best_sig = fit["best_model_signed"]
        r2_sig = fit["models_signed_deviation"][best_sig]["r_squared"]
        plaw = fit.get("power_law_log_log", {})
        plaw_str = ""
        if "r_squared" in plaw:
            plaw_str = (f", power-law |dev|~V^{plaw['exponent_b']:+.3f}, "
                        f"R^2={plaw['r_squared']:.3f}")
        print(f"  {family}: {fit['n_points']} points, "
              f"best(|dev|) = {best_abs} R^2={r2_abs:.4f}, "
              f"best(dev) = {best_sig} R^2={r2_sig:.4f}{plaw_str}")

    return out


if __name__ == "__main__":
    main()
