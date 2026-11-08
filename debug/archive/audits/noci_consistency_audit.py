"""NOCI sandbox step 3a: the surrogate-integral consistency audit (LiH lab).

Question (the step-3 bridge question, v4.73.1): NOCI needs a consistent (S, h, g)
triple.  GeoVac's builders carry surrogate h/g (no cross-center one-electron
off-diagonals = W1d; no cross-center overlap-density ERIs = the Gaunt-sparsity
exclusion).  Does bolting genuine overlaps S onto surrogate h/g spoil NOCI's
variational behavior -- and WHICH surrogate sector does the damage?

Method: in the LiH probe basis (where the fully genuine answer is known), degrade
the genuine integrals layer by layer, GeoVac-style, and rerun the ladder:

  A  genuine S, h, g                      (baseline, = step-2 result)
  B  genuine S, g;  h cross-center off-diagonals zeroed          (W1d mimic)
  C  genuine S, h;  g cross-center overlap-density ERIs zeroed   (sparsity mimic)
  D  genuine S;     both surrogates                              (builder-style H)
  E  S = I;         both surrogates      (orthonormal-pretense analog of the
                                          current native treatment)

Per run: 3-det NOCI + 15-det in-model FCI, well parameters against the same run's
R=12 energy (each model is its own dissociation reference), plus a non-variational
collapse check (in-model energies dipping below the genuine-FCI minimum signals the
Loewdin-retrofit failure mode).

Caveat (disclosed): the surrogates MIMIC builder structure on a textbook basis;
this localizes which zeroed sector breaks what, it does not reproduce the builders.
"""

from __future__ import annotations

import json
import os
from itertools import combinations
from typing import Dict, List, Tuple

import numpy as np

from noci_h2_probe import integral_set
from noci_lih_probe import (
    DN,
    H1S,
    LI1S,
    LI2S,
    UP,
    build_orbitals,
    noci_ground,
)

HARTREE_TO_EV = 27.211386
Z1S, Z2S, ZH = 2.69, 0.65, 1.0

CENTER = {LI1S: 0, LI2S: 0, H1S: 1}


def same_center(i: int, j: int) -> bool:
    return CENTER[i] == CENTER[j]


def make_surrogates(s: np.ndarray, h: np.ndarray, g: np.ndarray):
    h_surr = h.copy()
    n = h.shape[0]
    for i in range(n):
        for j in range(n):
            if i != j and not same_center(i, j):
                h_surr[i, j] = 0.0
    g_surr = g.copy()
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    if not (same_center(i, j) and same_center(k, l)):
                        g_surr[i, j, k, l] = 0.0
    return h_surr, g_surr


def main() -> None:
    core = [(LI1S, UP), (LI1S, DN)]
    det_cov_a = core + [(LI2S, UP), (H1S, DN)]
    det_cov_b = core + [(H1S, UP), (LI2S, DN)]
    det_ion_h = core + [(H1S, UP), (H1S, DN)]
    dets3 = [det_cov_a, det_cov_b, det_ion_h]
    spin_orbitals = [(p, sp) for p in (LI1S, LI2S, H1S) for sp in (UP, DN)]
    fci_dets = [list(c) for c in combinations(spin_orbitals, 4)]

    grid = [round(r, 2) for r in np.concatenate([np.arange(2.0, 6.01, 0.25), [8.0, 10.0, 12.0]])]

    runs = ["A genuine", "B h-surr", "C g-surr", "D h+g surr", "E S=I + h+g surr"]
    rows: List[Dict[str, float]] = []
    for r in grid:
        orbs, nuclei = build_orbitals(r, Z1S, Z2S, ZH)
        s, h, g = integral_set(orbs, nuclei)
        h_surr, g_surr = make_surrogates(s, h, g)
        eye = np.eye(3)
        vnn = 3.0 / r
        cases = {
            "A genuine": (s, h, g),
            "B h-surr": (s, h_surr, g),
            "C g-surr": (s, h, g_surr),
            "D h+g surr": (s, h_surr, g_surr),
            "E S=I + h+g surr": (eye, h_surr, g_surr),
        }
        row: Dict[str, float] = {"R": r}
        for label, (ss, hh, gg) in cases.items():
            e3, _, _ = noci_ground(dets3, ss, hh, gg)
            ef, _, _ = noci_ground(fci_dets, ss, hh, gg)
            row[f"{label} | noci3"] = e3 + vnn
            row[f"{label} | fci"] = ef + vnn
        rows.append(row)

    print("R-dependence (3-det NOCI columns):")
    print("{:>5}".format("R"), *[f"{lab:>18}" for lab in runs])
    for row in rows:
        print("{:>5.2f}".format(row["R"]), *[f"{row[f'{lab} | noci3']:>18.6f}" for lab in runs])

    def well(key: str) -> Tuple[float, float, float]:
        es = np.array([row[key] for row in rows])
        idx = int(np.argmin(es))
        r_eq, e_min = rows[idx]["R"], float(es[idx])
        if 0 < idx < len(rows) - 1:
            x = np.array([rows[idx - 1]["R"], rows[idx]["R"], rows[idx + 1]["R"]])
            y = np.array([es[idx - 1], es[idx], es[idx + 1]])
            c = np.polyfit(x, y, 2)
            r_eq = float(-c[1] / (2.0 * c[0]))
            e_min = float(np.polyval(c, r_eq))
        e_ref = float(es[-1])  # R = 12 within the SAME model
        return r_eq, e_min, e_ref - e_min

    e_fci_genuine_min = well("A genuine | fci")[1]

    print("\n[wells]  (D_e vs the SAME model's R=12 energy; genuine baseline: "
          "R_eq 3.26 a0, D_e 1.33 eV FCI / 1.30 eV noci3)")
    summary = {}
    for lab in runs:
        for kind in ("noci3", "fci"):
            key = f"{lab} | {kind}"
            r_eq, e_min, d_e = well(key)
            interior = r_eq < grid[-4]
            binds = d_e > 1e-4 and interior
            collapse = e_min < e_fci_genuine_min - 5e-3
            print(
                f"  {key:28s} R_eq = {r_eq:7.3f}  E_min = {e_min:11.6f}  "
                f"D_e = {d_e * HARTREE_TO_EV:7.3f} eV  "
                f"{'BINDS ' if binds else 'UNBOUND'}"
                f"{'  ** NON-VARIATIONAL COLLAPSE **' if collapse else ''}"
            )
            summary[key] = {
                "R_eq": r_eq, "E_min": e_min, "D_e_eV": d_e * HARTREE_TO_EV,
                "binds": bool(binds), "collapse": bool(collapse),
            }

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                            "noci_consistency_audit_results.json")
    with open(out_path, "w") as fh:
        json.dump({"rows": rows, "wells": summary}, fh, indent=1)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()
