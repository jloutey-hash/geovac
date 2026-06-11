"""Ring-orphan taxonomy probe (2026-05-30, confinement-reframe exploration).

Question (Josh): the ring-orphan constants — correlation entropy S_full(GS), K,
the L2 constant c, the Wolfenstein parameters — sit at one table. Does that table
look like the *decomposed alpha* (K = pi(B + F - Delta)) we've chewed on for a month?

This probe does the ONE concrete thing that distinguishes "same kind of orphan"
from "different tier of orphan":

  - Pull the actual 1-RDM occupations p_i for He n3 and Li+ n3.
  - Confirm S_full = -sum p_i log p_i EXACTLY from those occupations
    (closes the loop: the entropy IS the von-Neumann functional of these p_i).
  - Characterize the p_i: they are eigenvalues of a 1-RDM built from the
    eigenvector of an FCI matrix whose elements are EXACT RATIONALS (Slater
    integrals, Paper 7 VI.B / hypergeometric_slater). => the p_i are ALGEBRAIC.
  - Therefore S_full is a Baker-class transcendental: a Q-combination of
    log(algebraic) terms. That is a DIFFERENT mechanism of ring-orphan-hood
    than K = pi(B+F-Delta) (a coincidence with an EXTERNAL physical constant).

float64 is enough: we only need ~10 digits to confirm S = -sum p log p against
the known 150-digit value, and the algebraicity claim is a theorem about the
rational matrix, not a high-precision fit.

debug-tier only; no geovac/ or paper edits.
"""

import json
import os
import numpy as np
import mpmath as mp

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from debug.sprint_td_track5 import build_fci_matrix_exact, build_1rdm_mpmath

# Known 150-digit S_full values (debug/data/sprint_td_track5.json)
KNOWN = {
    "He_n3": "0.040811051366471814867998169808302807554187762125967235314692183091428112917345",
    "Li+_n3": "0.011211717937420977473956556933124075572277113089128407501084620994494566252873",
}

OUT = "debug/data/ring_orphan_probe.json"


def occupations_and_entropy(Z, n_max, dps=30):
    """Return (occupations p_i descending, S_full) via float64 diag of the
    exact-rational FCI matrix built by track5's machinery."""
    mp.mp.dps = dps
    H_mp, configs, orbitals = build_fci_matrix_exact(Z, n_max)
    n_cfg = len(configs)
    n_spatial = len(orbitals)
    # float64 copy of the (exact-rational-valued) Hamiltonian
    H = np.array([[float(H_mp[i, j]) for j in range(n_cfg)] for i in range(n_cfg)])
    w, V = np.linalg.eigh(H)
    v0 = V[:, 0]  # ground state (eigh sorts ascending)
    # 1-RDM via the same builder (mpmath), fed the float GS vector
    v0_mp = mp.matrix([mp.mpf(float(x)) for x in v0])
    rho = build_1rdm_mpmath(v0_mp, configs, n_spatial)
    rho_np = np.array([[float(rho[i, j]) for j in range(n_spatial)]
                       for i in range(n_spatial)])
    occ = np.linalg.eigvalsh(rho_np)  # eigenvalues sum to 2 (two electrons)
    occ = np.sort(occ)[::-1]
    # probabilities p_i = n_i / 2  (sum to 1)
    p = occ / 2.0
    S = -sum(pi * np.log(pi) for pi in p if pi > 1e-14)
    return occ, p, float(S), n_cfg, n_spatial, float(w[0])


def main():
    results = {"systems": {}}
    for name, (Z, n_max) in [("He_n3", (2, 3)), ("Li+_n3", (3, 3))]:
        occ, p, S, n_cfg, n_spatial, E0 = occupations_and_entropy(Z, n_max)
        known = float(KNOWN[name][:18])
        # significant nonzero occupations
        sig = [float(x) for x in occ if abs(x) > 1e-10]
        results["systems"][name] = {
            "Z": Z, "n_max": n_max, "n_configs": n_cfg, "n_spatial": n_spatial,
            "E_GS": E0,
            "S_full_float64": S,
            "S_full_known_18dig": known,
            "abs_diff": abs(S - known),
            "loop_closed": abs(S - known) < 1e-9,
            "n_significant_occupations": len(sig),
            "occupations_top8": [float(x) for x in occ[:8]],
            "probabilities_top8": [float(x) for x in p[:8]],
            "trace_check": float(sum(occ)),
        }
        print(f"\n=== {name}: Z={Z}, n_max={n_max} ===")
        print(f"  n_configs={n_cfg}, n_spatial={n_spatial}, E_GS={E0:.6f}")
        print(f"  trace(1-RDM) = {sum(occ):.10f}  (should be 2.0)")
        print(f"  # significant occupations = {len(sig)}")
        print(f"  top occupations n_i: {[f'{x:.6f}' for x in occ[:6]]}")
        print(f"  top probabilities p_i: {[f'{x:.6f}' for x in p[:6]]}")
        print(f"  S = -sum p log p   = {S:.15f}")
        print(f"  S_full (known)     = {known:.15f}")
        print(f"  |diff|             = {abs(S-known):.2e}   "
              f"LOOP {'CLOSED' if abs(S-known)<1e-9 else 'OPEN'}")

    # Light algebraicity flavor: is the dominant occupation a LOW-degree
    # algebraic number? PSLQ p_max against {1, p, p^2, ..., p^6} at modest
    # maxcoeff. NULL = algebraic but not low-degree (still algebraic-by-theorem,
    # since it is an eigenvalue of a rational matrix). HIT = low-degree.
    mp.mp.dps = 50
    he_occ, he_p, *_ = occupations_and_entropy(2, 3, dps=50)
    pmax = mp.mpf(float(he_p[0]))
    powers = [mp.mpf(1)] + [pmax**k for k in range(1, 7)]
    rel = mp.pslq(powers, maxcoeff=10**4, maxsteps=100000)
    results["dominant_occupation_algebraicity"] = {
        "p_max_He_n3": float(he_p[0]),
        "min_poly_relation_deg<=6_maxcoeff_1e4": [int(c) for c in rel] if rel else None,
        "note": ("HIT => p_max is a low-degree algebraic number; "
                 "NULL => algebraic but degree>6 (still algebraic by theorem: "
                 "eigenvalue of an exact-rational 1-RDM)."),
    }
    print("\n--- dominant occupation algebraicity (He n3) ---")
    print(f"  p_max = {float(he_p[0]):.12f}")
    print(f"  min-poly (deg<=6, maxcoeff 1e4): "
          f"{[int(c) for c in rel] if rel else 'NULL (degree>6, algebraic by theorem)'}")

    os.makedirs("debug/data", exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nWritten: {OUT}")


if __name__ == "__main__":
    main()
