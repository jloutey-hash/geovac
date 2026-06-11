"""Calc Track He-Oscillator v3: Phase D extension — extended angular CI.

The Phase B+C result (debug/data/he_oscillator_v2.json) showed that:
  * Path A (Slater per-orbital lambda, n_max=4): f = 0.380 (+37% vs Drake)
    with E(1S), E(2P) accurate to 0.1% and omega accurate to 0.06%.
  * Path C5 (saturated, 12 orbitals at hand-picked lambdas): f = 0.278
    (+0.6%) at cond(S) ~ 10^10.

The Phase A→C5 jump is mostly at the dipole matrix element level
(omega only changed 3%, dipole^2 changed 25%). The Phase B+C subblock
builder restricts to (s, s) configurations in 1S and (s, p) in 1P.

Phase D extends the angular CI to include ALL singlet-allowed (l_a, l_b)
pairs — for 1S that adds (p, p), (d, d); for 1P that adds (p, d).

This is option (3) from §5.3 of the implementation memo: extend the
1S subblock to include 2p^2 configurations.

The hypothesis:
  * Path C5's f closure was partly due to its `2p_quad` orbital
    enrichment compensating for the missing (p, p) admixture in 1S.
  * Adding (p, p) admixture to the 1S CI explicitly should close the
    f gap at well-conditioned modest basis (cond < 10^4) without
    needing the C5 saturation.

Phase D test plan:
  D1: Diagnose Path A vs C5 mechanism (compare config dominance).
  D2: Run extended angular CI on Slater-rules basis with full m
      coverage (he_extended_spec). Compare to Path A.
  D3: Sweep n_max for the extended CI, track conditioning vs accuracy.

Date: 2026-05-09
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.internal_multifocal import (
    MultifocalOrbital,
    MultifocalSpec,
    _enumerate_singlet_l_pairs,
    build_singlet_LM_subblock_multifocal,
    build_singlet_LM_subblock_multifocal_extended,
    compute_he_oscillator_strength_multifocal,
    compute_he_oscillator_strength_multifocal_extended,
    he_extended_spec,
    he_slater_spec,
    overlap_matrix,
    solve_generalized_singlet,
    transition_dipole_multifocal,
)

# Reference value
F_DRAKE = 0.27616
DRAKE_E_1S = -2.903724  # Hartree, NR infinite-mass
DRAKE_E_2P = -2.123843
DRAKE_OMEGA = DRAKE_E_2P - DRAKE_E_1S


def _summarize(r: Dict[str, Any], label: str) -> None:
    f = r["f_length"]
    err = (f - F_DRAKE) / F_DRAKE * 100.0
    omega_err = (r["omega_Ha"] - DRAKE_OMEGA) / DRAKE_OMEGA * 100.0
    e1s_err = (r["E_1S_Ha"] - DRAKE_E_1S) * 1000.0
    e2p_err = (r["E_2P_Ha"] - DRAKE_E_2P) * 1000.0
    print(f"  {label}")
    print(f"    f = {f:.6f}  ({err:+.2f}% vs Drake)")
    print(f"    omega = {r['omega_Ha']:.4f} Ha  ({omega_err:+.3f}% vs Drake)")
    print(f"    E_1S = {r['E_1S_Ha']:.4f} ({e1s_err:+.1f} mHa), "
          f"E_2P = {r['E_2P_Ha']:.4f} ({e2p_err:+.1f} mHa)")
    print(f"    n_orb={r['n_orbitals']}, "
          f"configs_1S={r['n_configs_1S']}, configs_2P={r['n_configs_2P']}, "
          f"cond_S_1S={r['cond_S_1S']:.2e}, cond_S_2P={r['cond_S_2P']:.2e}")
    if "config_l_pairs_1S" in r:
        print(f"    l-pairs 1S: {r['config_l_pairs_1S']}, "
              f"l-pairs 2P: {r['config_l_pairs_2P']}")


# -------------------------------------------------------------------------
# D1 — Diagnose Path C5 mechanism
# -------------------------------------------------------------------------

def diagnose_path_c5(spec: MultifocalSpec, label: str) -> Dict[str, Any]:
    """Run the standard (non-extended) driver on a custom spec and dump
    detailed config-level dipole decomposition."""
    r_basic = compute_he_oscillator_strength_multifocal(spec, n_quad=100)

    # Re-run to get detailed breakdown
    H_1S, S_1S, configs_1S = build_singlet_LM_subblock_multifocal(
        spec, L_target=0, M_L_target=0, n_quad=100,
    )
    eigvals_1S, eigvecs_1S = solve_generalized_singlet(H_1S, S_1S)
    psi_1S = eigvecs_1S[:, 0]

    H_2P, S_2P, configs_2P = build_singlet_LM_subblock_multifocal(
        spec, L_target=1, M_L_target=0, n_quad=100,
    )
    eigvals_2P, eigvecs_2P = solve_generalized_singlet(H_2P, S_2P)
    psi_2P = eigvecs_2P[:, 0]

    # Compute per-(I,J) dipole contributions (un-summed)
    # |<Psi_init | z | Psi_final>|^2 = (sum_{I,J} c_I c_J <I|z|J>)^2
    # Decompose:  <I|z|J> matrix
    n_1S, n_2P = len(configs_1S), len(configs_2P)
    me_matrix = np.zeros((n_1S, n_2P))
    for I in range(n_1S):
        for J in range(n_2P):
            psi_I = np.zeros(n_1S); psi_I[I] = 1.0
            psi_J = np.zeros(n_2P); psi_J[J] = 1.0
            me_matrix[I, J] = transition_dipole_multifocal(
                psi_I, configs_1S, psi_J, configs_2P,
                spec_init=spec, spec_final=spec,
            )

    # Total dipole
    total_dipole = float(psi_1S @ me_matrix @ psi_2P)

    # Per-config contributions (which I, J dominate?)
    contrib = np.outer(psi_1S, psi_2P) * me_matrix
    flat = contrib.flatten()
    abs_idx = np.argsort(np.abs(flat))[::-1]

    print(f"\n  C5 dipole decomposition ({label}):")
    print(f"    Total dipole = {total_dipole:.4f}")
    print(f"    Top-5 (I, J, c_I, c_J, <I|z|J>, contribution):")
    for k in range(min(5, len(abs_idx))):
        I = abs_idx[k] // n_2P
        J = abs_idx[k] % n_2P
        i, j = configs_1S[I]
        p, q = configs_2P[J]
        oi, oj = spec.orbitals[i], spec.orbitals[j]
        op, oq = spec.orbitals[p], spec.orbitals[q]
        print(f"      I=({oi.label}, {oj.label})  "
              f"J=({op.label}, {oq.label})  "
              f"c_I={psi_1S[I]:+.3f}  c_J={psi_2P[J]:+.3f}  "
              f"<I|z|J>={me_matrix[I, J]:+.4f}  "
              f"contrib={contrib[I, J]:+.4f}")

    # Top-coefficient configs in psi_1S, psi_2P
    print(f"    Top-3 1S CI coefficients:")
    for k in np.argsort(np.abs(psi_1S))[::-1][:3]:
        i, j = configs_1S[k]
        oi, oj = spec.orbitals[i], spec.orbitals[j]
        print(f"      c={psi_1S[k]:+.4f}  config=({oi.label}, {oj.label})")
    print(f"    Top-3 2P CI coefficients:")
    for k in np.argsort(np.abs(psi_2P))[::-1][:3]:
        i, j = configs_2P[k]
        oi, oj = spec.orbitals[i], spec.orbitals[j]
        print(f"      c={psi_2P[k]:+.4f}  config=({oi.label}, {oj.label})")

    return {
        "label": label,
        "f_length": r_basic["f_length"],
        "total_dipole": total_dipole,
        "top_contributions": [
            {
                "I_label": (spec.orbitals[configs_1S[abs_idx[k]//n_2P][0]].label,
                            spec.orbitals[configs_1S[abs_idx[k]//n_2P][1]].label),
                "J_label": (spec.orbitals[configs_2P[abs_idx[k]%n_2P][0]].label,
                            spec.orbitals[configs_2P[abs_idx[k]%n_2P][1]].label),
                "c_I": float(psi_1S[abs_idx[k]//n_2P]),
                "c_J": float(psi_2P[abs_idx[k]%n_2P]),
                "me_IJ": float(me_matrix[abs_idx[k]//n_2P, abs_idx[k]%n_2P]),
                "contrib": float(contrib[abs_idx[k]//n_2P, abs_idx[k]%n_2P]),
            }
            for k in range(min(5, len(abs_idx)))
        ],
    }


def main() -> Dict[str, Any]:
    print("=" * 72)
    print("He 2^1P -> 1^1S oscillator strength: Phase D (extended angular CI)")
    print("=" * 72)
    print(f"Drake reference: f = {F_DRAKE:.5f}, omega = {DRAKE_OMEGA:.5f} Ha")
    print(f"Track 4 v1 (single-focal):    f = 0.444 (+61%)")
    print(f"Phase A (Slater n_max=4):     f = 0.380 (+37%)")
    print(f"Phase C2 (well-conditioned):  f = 0.314 (+14%)")
    print(f"Phase C5 (saturated, cond=1e10): f = 0.278 (+0.6%)")
    print()

    results: Dict[str, Any] = {
        "reference": {
            "f_drake": F_DRAKE,
            "E_1S_drake": DRAKE_E_1S,
            "E_2P_drake": DRAKE_E_2P,
            "omega_drake": DRAKE_OMEGA,
        },
        "phase_d1_diagnosis": [],
        "phase_d2_extended_slater": [],
        "phase_d3_extended_enriched": [],
    }

    # ---------------------------------------------------------------
    # D1 — Diagnose Path A and Path C5 mechanism
    # ---------------------------------------------------------------
    print("=" * 72)
    print("D1: DIAGNOSIS — Path A and Path C5 dipole decomposition")
    print("=" * 72)

    # Path A n_max=4 (m_p=0 only)
    spec_A = he_slater_spec(n_max=4)
    print(f"\nPath A: Slater n_max=4 (m_p=0 only):")
    diag_A = diagnose_path_c5(spec_A, "PathA_n4")
    results["phase_d1_diagnosis"].append(diag_A)

    # Path C5 (reproduce from v2)
    Z = 2.0
    spec_C5_orbs = []
    for lam in [1.4, 27.0/16.0, 1.95]:
        spec_C5_orbs.append(MultifocalOrbital(n=1, l=0, m=0, lam=lam,
                                              label=f"1s_lam{lam:.3f}"))
    for lam in [0.4, 0.7]:
        spec_C5_orbs.append(MultifocalOrbital(n=2, l=0, m=0, lam=lam,
                                              label=f"2s_lam{lam:.3f}"))
    for lam in [0.3, 0.5]:
        spec_C5_orbs.append(MultifocalOrbital(n=3, l=0, m=0, lam=lam,
                                              label=f"3s_lam{lam:.3f}"))
    for lam in [0.3, 0.5, 0.75, 1.0]:
        spec_C5_orbs.append(MultifocalOrbital(n=2, l=1, m=0, lam=lam,
                                              label=f"2p_lam{lam:.3f}"))
    spec_C5_orbs.append(MultifocalOrbital(n=3, l=1, m=0, lam=0.333,
                                          label="3p_lam0.333"))
    spec_C5 = MultifocalSpec(orbitals=spec_C5_orbs, Z_nuc=Z, label="C5")
    print(f"\nPath C5 (saturated, 12 orbitals):")
    diag_C5 = diagnose_path_c5(spec_C5, "PathC5_saturated")
    results["phase_d1_diagnosis"].append(diag_C5)

    # ---------------------------------------------------------------
    # D2 — Extended angular CI on Slater-rules basis with full m
    # ---------------------------------------------------------------
    print()
    print("=" * 72)
    print("D2: EXTENDED ANGULAR CI on Slater-rules basis (all m sublevels)")
    print("=" * 72)

    for n_max in [2, 3, 4]:
        spec_ext = he_extended_spec(n_max=n_max)
        print(f"\n  n_max = {n_max} (he_extended_spec, all m for l>=1):")

        # First: standard subblocks (Phase B style) on this spec
        try:
            t0 = time.time()
            r_basic = compute_he_oscillator_strength_multifocal(
                spec_ext, n_quad=80,
            )
            elapsed = time.time() - t0
            r_basic["mode"] = "phase_b_subblocks"
            r_basic["n_max"] = n_max
            r_basic["wall_seconds"] = elapsed
            results["phase_d2_extended_slater"].append(r_basic)
            _summarize(r_basic, f"Phase-B-style (s, s)/(s, p) only, n_max={n_max}")
        except Exception as e:
            print(f"    Phase-B-style FAILED: {e}")

        # Second: extended subblocks (Phase D style) with all l-pairs
        try:
            t0 = time.time()
            r_ext = compute_he_oscillator_strength_multifocal_extended(
                spec_ext, n_quad=80,
            )
            elapsed = time.time() - t0
            r_ext["mode"] = "phase_d_extended"
            r_ext["n_max"] = n_max
            r_ext["wall_seconds"] = elapsed
            results["phase_d2_extended_slater"].append(r_ext)
            _summarize(r_ext, f"Phase-D extended (all (l_a, l_b)), n_max={n_max}")
        except Exception as e:
            print(f"    Phase-D extended FAILED: {e}")
            results["phase_d2_extended_slater"].append({
                "n_max": n_max, "mode": "phase_d_extended", "error": str(e),
            })

    # ---------------------------------------------------------------
    # D3 — Extended angular CI with custom enrichment
    # ---------------------------------------------------------------
    print()
    print("=" * 72)
    print("D3: ENRICHMENT SWEEP on extended angular CI")
    print("=" * 72)

    Z = 2.0

    enrichment_specs = [
        # E1: Slater n_max=2 with multiple p-exponents (capture 2p shape)
        ("E1: n_max=2 + 2p_double",
         {"s_lams": [27.0/16.0, 0.575],
          "p_lams": [0.4, 0.7],
          "d_lams": []}),

        # E2: n_max=3 + p-exponent enrichment
        ("E2: n_max=3 + 2p_triple",
         {"s_lams": [27.0/16.0, 0.575, 0.333],
          "p_lams": [0.4, 0.7, 1.0],
          "d_lams": [0.333]}),

        # E3: n_max=3 + 1s_double + 2p_triple (mirrors C3 but with extended CI)
        ("E3: 1s_double + 2p_triple + d_block",
         {"s_lams": [1.5, 27.0/16.0, 0.575, 0.333],
          "p_lams": [0.4, 0.7, 1.0],
          "d_lams": [0.333]}),

        # E4: equivalent of C5 but with extended CI machinery
        ("E4: 1s_triple + 2s_double + 2p_quad (C5-equivalent extended)",
         {"s_lams": [1.4, 27.0/16.0, 1.95, 0.4, 0.7, 0.3, 0.5],
          "p_lams": [0.3, 0.5, 0.75, 1.0, 0.333],
          "d_lams": []}),

        # E5: BIG extended basis — 2 1s + 2p_quad + d_double, all m
        ("E5: 1s_double + 2p_quad + d_double",
         {"s_lams": [1.5, 27.0/16.0, 0.575, 0.333],
          "p_lams": [0.3, 0.5, 0.75, 1.0],
          "d_lams": [0.333, 0.5]}),
    ]

    for label, kwargs in enrichment_specs:
        spec = he_extended_spec(n_max=3, **kwargs)
        print(f"\n  {label}:")
        try:
            t0 = time.time()
            r = compute_he_oscillator_strength_multifocal_extended(
                spec, n_quad=80,
            )
            elapsed = time.time() - t0
            r["label"] = label
            r["wall_seconds"] = elapsed
            results["phase_d3_extended_enriched"].append(r)
            _summarize(r, label)
        except Exception as e:
            print(f"    FAILED ({type(e).__name__}: {e})")
            results["phase_d3_extended_enriched"].append({
                "label": label, "error": str(e),
            })

    # ---------------------------------------------------------------
    # VERDICT
    # ---------------------------------------------------------------
    print()
    print("=" * 72)
    print("PHASE D VERDICT")
    print("=" * 72)

    # Find best well-conditioned (cond < 1e4) result
    all_results = results["phase_d2_extended_slater"] + results["phase_d3_extended_enriched"]
    well_cond = [r for r in all_results if "f_length" in r and
                 r.get("cond_S_1S", 1e20) < 1e4 and r.get("cond_S_2P", 1e20) < 1e4]
    if well_cond:
        best_wc = min(well_cond,
                      key=lambda r: abs(r["f_length"] - F_DRAKE))
        f_wc = best_wc["f_length"]
        err_wc = (f_wc - F_DRAKE) / F_DRAKE * 100.0
        print(f"\nBest WELL-CONDITIONED (cond_S < 1e4) result:")
        label = best_wc.get("label", best_wc.get("mode", f"n_max={best_wc.get('n_max')}"))
        print(f"  {label}")
        print(f"  f = {f_wc:.6f}  (err = {err_wc:+.2f}%)")
        print(f"  cond_S_1S = {best_wc['cond_S_1S']:.2e}, "
              f"cond_S_2P = {best_wc['cond_S_2P']:.2e}")
        print(f"  configs_1S = {best_wc['n_configs_1S']}, "
              f"configs_2P = {best_wc['n_configs_2P']}")

        if abs(err_wc) < 5:
            verdict = "WIN — extended angular CI closes residual to <5% at well-conditioned basis"
        elif abs(err_wc) < 10:
            verdict = "STRONG PARTIAL — extended angular CI closes to <10%"
        elif abs(err_wc) < 14:
            verdict = "PARTIAL — extended angular CI improves on Phase C2 (+14%) plateau"
        else:
            verdict = "WEAK — extended angular CI does not improve on Phase C2 plateau"
        print(f"\n  Verdict: {verdict}")
        results["verdict"] = verdict

    # Best overall (any conditioning)
    all_with_f = [r for r in all_results if "f_length" in r]
    if all_with_f:
        best_all = min(all_with_f,
                       key=lambda r: abs(r["f_length"] - F_DRAKE))
        f_all = best_all["f_length"]
        err_all = (f_all - F_DRAKE) / F_DRAKE * 100.0
        print(f"\nBest OVERALL (any conditioning):")
        label = best_all.get("label", best_all.get("mode", f"n_max={best_all.get('n_max')}"))
        print(f"  {label}")
        print(f"  f = {f_all:.6f}  (err = {err_all:+.2f}%)")
        print(f"  cond_S_1S = {best_all['cond_S_1S']:.2e}")

    # Save
    out_file = PROJECT_ROOT / "debug" / "data" / "he_oscillator_v3.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as fp:
        json.dump(results, fp, indent=2, default=lambda x: float(x) if hasattr(x, '__float__') else str(x))
    print(f"\nSaved: {out_file}")
    return results


if __name__ == "__main__":
    main()
