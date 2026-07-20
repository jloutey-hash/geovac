"""NOCI sandbox probe, step N4-as-diagnosis: NaH all-electron, genuine integrals.

Sandbox exploration (branch sandbox/noci, repo frozen at v4.76.0 -- no release
path).  PI-authorized 2026-07-19 ("Avery style, soup-to-nuts") after the N3b
gate STOP: this run makes NO sparsity claim.  The question is purely
diagnostic: NaH is the system the corpus framework could NEVER bind (W1e --
balanced/composed PES monotonically descends into small R; Sprint B.1
explicit-core HF, R3-B DMRG, F1-F6 all negative).  If a compact NOCI over
GENUINE end-to-end integrals binds it with roughly the right geometry, the
W1e localization (wall lives at Hamiltonian-specification, i.e. the missing/
surrogate cross-center physics) is confirmed by construction.

Design (mirrors N2, scaled up):
  - All-electron, 12 electrons.  Basis: Na {1s, 2s, 2px, 2py, 2pz, 3s},
    H {1s} -- M = 7 spatial orbitals (14 spin-orbitals).
  - Integrals: McMurchie-Davidson s/p engine (noci_md_engine, validated at
    machine precision against N1 closed forms + FD derivatives).
  - Shapes: hardcoded STO-6G 1s and STO-3G 2s (N1/N2 lineage); 2p and 3s
    least-squares 6-Gaussian fits of the zeta=1 Slater shapes (<fit|STO> =
    1.000000).  Exponents variational ON THE ISOLATED Na ATOM ONLY (in-basis
    FCI, dim 12); H zeta = 1.  The molecular curve is a fragment prediction.
  - NOCI: Loewdin cofactor rules (validated vs permutation machinery and vs
    the stored N2 LiH ladder at 1e-15).  Ladder: covalent (2 dets) ->
    + ionic Na+H- (3) -> + ionic Na-H+ (4) -> in-basis FCI (91 dets,
    bitstring, Loewdin basis).  Loewdin-orbital truncations for comparison.
  - Span identity check at one R: non-orthogonal complete space (91 dets via
    cofactor machinery) must equal the Loewdin bitstring FCI.
  - Rotational invariance check at one R: molecule along z vs along
    (1,1,1)/sqrt(3) -- FCI energies must agree (p-block covariance).

Pre-registered gates (values verified 2026-07-19):
  G1 (headline): interior minimum in the NOCI-3 PES.
  G2: R_eq near experiment 3.566 a0 (r_e = 1.887 A, NIST CCCBDB).
  G3: D_e > 0 and order-1 eV (exp D_e = 15,815 cm^-1 = 1.961 eV, Huang et al.
      JCP 133, 044301 (2010)); minimal basis expected to underbind (N2's LiH
      captured 53% of experimental D_e).
  G4: NOCI-3 captures >~90% of in-basis FCI binding (N2 pattern: 97.6%).

Known limitations (disclosed): no Na 3p (no sp polarization), no diffuse H
(anion character underdescribed), minimal basis, no BSSE correction (atom
references computed in atom-only bases).
"""

from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from noci_md_engine import (
    HARTREE_TO_EV,
    STO3G_2S,
    STO6G_1S,
    BasisFn,
    det_pair_gensc,
    fci_ground,
    fit_sto_shape,
    integral_set_md,
    lowdin_orbitals,
    noci_ground_gensc,
    transform_integrals,
)

UP, DN = 0, 1
R_EQ_EXP = 3.566          # a0  (1.887 A, NIST CCCBDB)
D_E_EXP_EV = 1.961        # eV  (15,815 cm^-1, Huang et al. JCP 133 044301)


# ------------------------------------------------------------------ basis

def build_shapes() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    shapes = {"1s": STO6G_1S, "2s": STO3G_2S}
    for kind, (l, n_r) in (("2p", (1, 2)), ("3s", (0, 3))):
        a, dco, q = fit_sto_shape(l, n_r)
        print(f"[shapes] fitted {kind}: <fit|STO> = {q:.6f}")
        shapes[kind] = (a, dco)
    return shapes


def na_orbitals(center: np.ndarray, z: Dict[str, float], shapes) -> List[BasisFn]:
    a1, d1 = shapes["1s"]
    a2, d2 = shapes["2s"]
    ap, dp = shapes["2p"]
    a3, d3 = shapes["3s"]
    return [
        BasisFn(center, (0, 0, 0), a1 * z["1s"] ** 2, d1),
        BasisFn(center, (0, 0, 0), a2 * z["2s"] ** 2, d2),
        BasisFn(center, (1, 0, 0), ap * z["2p"] ** 2, dp),
        BasisFn(center, (0, 1, 0), ap * z["2p"] ** 2, dp),
        BasisFn(center, (0, 0, 1), ap * z["2p"] ** 2, dp),
        BasisFn(center, (0, 0, 0), a3 * z["3s"] ** 2, d3),
    ]


def h_orbital(center: np.ndarray, zeta: float, shapes) -> BasisFn:
    a1, d1 = shapes["1s"]
    return BasisFn(center, (0, 0, 0), a1 * zeta ** 2, d1)


def na_atom_fci(z: Dict[str, float], shapes) -> float:
    orig = np.array([0.0, 0.0, 0.0])
    orbs = na_orbitals(orig, z, shapes)
    s, h, g = integral_set_md(orbs, [(orig, 11.0)])
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    return fci_ground(ht, gt, 11)


def select_zetas(shapes) -> Tuple[Dict[str, float], float]:
    """Variational in-basis FCI on the isolated Na atom (atoms only)."""
    grids = {
        "1s": [10.3, 10.63, 11.0],
        "2s": [3.1, 3.3, 3.5],
        "2p": [3.2, 3.44, 3.7],
        "3s": [0.75, 0.836, 0.92, 1.0],
    }
    best = None
    t0 = time.time()
    for z1 in grids["1s"]:
        for z2 in grids["2s"]:
            for zp in grids["2p"]:
                for z3 in grids["3s"]:
                    z = {"1s": z1, "2s": z2, "2p": zp, "3s": z3}
                    e = na_atom_fci(z, shapes)
                    if best is None or e < best[1]:
                        best = (z, e)
    print(f"[atoms] zeta grid ({time.time() - t0:.0f} s): "
          f"best Na zetas = {best[0]}  E(Na, in-basis FCI) = {best[1]:.6f} Ha")
    print("[atoms] anchor: minimal-STO Na RHF ~ -161.12 Ha (Clementi); "
          "in-basis value should sit in that region")
    return best


# ------------------------------------------------------------------ dets

NA1S, NA2S, NA2PX, NA2PY, NA2PZ, NA3S, H1S = range(7)
CORE = [(o, sp) for o in (NA1S, NA2S, NA2PX, NA2PY, NA2PZ) for sp in (UP, DN)]

DET_COV_A = CORE + [(NA3S, UP), (H1S, DN)]
DET_COV_B = CORE + [(H1S, UP), (NA3S, DN)]
DET_ION_H = CORE + [(H1S, UP), (H1S, DN)]     # Na+ H-
DET_ION_NA = CORE + [(NA3S, UP), (NA3S, DN)]  # Na- H+

LADDERS = {
    "cov (2 dets)": [DET_COV_A, DET_COV_B],
    "cov+ionH (3 dets)": [DET_COV_A, DET_COV_B, DET_ION_H],
    "all 4 dets": [DET_COV_A, DET_COV_B, DET_ION_H, DET_ION_NA],
}


def molecule_integrals(r: float, z_na: Dict[str, float], z_h: float, shapes,
                       axis: np.ndarray = None):
    axis = np.array([0.0, 0.0, 1.0]) if axis is None else axis / np.linalg.norm(axis)
    pos_na = np.array([0.0, 0.0, 0.0])
    pos_h = r * axis
    orbs = na_orbitals(pos_na, z_na, shapes) + [h_orbital(pos_h, z_h, shapes)]
    nuclei = [(pos_na, 11.0), (pos_h, 1.0)]
    return integral_set_md(orbs, nuclei), 11.0 / r


# ------------------------------------------------------------------ main

def main() -> None:
    print("=== N4-as-diagnosis: NaH all-electron NOCI, genuine integrals ===\n")
    shapes = build_shapes()

    z_na, e_na = select_zetas(shapes)
    z_h = 1.0
    orig = np.array([0.0, 0.0, 0.0])
    hchi = h_orbital(orig, z_h, shapes)
    s1, h1, _ = integral_set_md([hchi], [(orig, 1.0)])
    e_h = h1[0, 0] / s1[0, 0]
    e_diss = e_na + e_h
    print(f"[atoms] E(H) = {e_h:.6f} Ha  ->  dissociation reference "
          f"{e_diss:.6f} Ha\n")

    # --- consistency checks at R = 3.5 ---
    (s, h, g), vnn = molecule_integrals(3.5, z_na, z_h, shapes)
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    e_fci_z = fci_ground(ht, gt, 12) + vnn

    (s2, h2, g2), _ = molecule_integrals(3.5, z_na, z_h, shapes,
                                         axis=np.array([1.0, 1.0, 1.0]))
    x2 = lowdin_orbitals(s2)
    ht2, gt2 = transform_integrals(x2, h2, g2)
    e_fci_rot = fci_ground(ht2, gt2, 12) + vnn
    print(f"[checks] rotational invariance |dE_FCI| (z vs 111 axis) = "
          f"{abs(e_fci_z - e_fci_rot):.2e}")

    t0 = time.time()
    sos = [(p, sp) for p in range(7) for sp in (UP, DN)]
    all_dets = [[sos[k] for k in c] for c in combinations(range(14), 12)]
    e_span, _ = noci_ground_gensc(all_dets, s, h, g)
    print(f"[checks] span identity |FCI_nonorth - FCI_Loewdin| = "
          f"{abs(e_span + vnn - e_fci_z):.2e}  "
          f"({len(all_dets)} dets, {time.time() - t0:.0f} s)\n")

    # --- PES sweep ---
    grid = [round(r, 2) for r in np.concatenate(
        [np.arange(2.0, 6.01, 0.25), [7.0, 8.0, 10.0]])]
    rows: List[Dict[str, float]] = []
    for r in grid:
        (s, h, g), vnn = molecule_integrals(r, z_na, z_h, shapes)
        row: Dict[str, float] = {"R": r, "S_3s_H": s[NA3S, H1S],
                                 "S_2pz_H": s[NA2PZ, H1S]}
        for label, dets in LADDERS.items():
            e, cond = noci_ground_gensc(dets, s, h, g)
            row[label] = e + vnn
            if label == "cov+ionH (3 dets)":
                row["cond_S_config"] = cond
        x = lowdin_orbitals(s)
        ht, gt = transform_integrals(x, h, g)
        row["FCI (91 dets)"] = fci_ground(ht, gt, 12) + vnn
        st = np.eye(7)
        for label, key in (("cov (2 dets)", "Loewdin cov (2 dets)"),
                           ("cov+ionH (3 dets)", "Loewdin cov+ionH (3 dets)")):
            e, _ = noci_ground_gensc(LADDERS[label], st, ht, gt)
            row[key] = e + vnn
        rows.append(row)
        print(f"  R = {r:5.2f} done  (FCI = {row['FCI (91 dets)']:.6f})")

    keys = ["cov (2 dets)", "cov+ionH (3 dets)", "all 4 dets", "FCI (91 dets)",
            "Loewdin cov (2 dets)", "Loewdin cov+ionH (3 dets)"]
    print("\n{:>5} {:>8}".format("R", "S_3s,H"), *[f"{k:>24}" for k in keys])
    for row in rows:
        print("{:>5.2f} {:>8.4f}".format(row["R"], row["S_3s_H"]),
              *[f"{row[k]:>24.6f}" for k in keys])

    def well(key: str) -> Tuple[float, float, float]:
        es = np.array([row[key] for row in rows])
        idx = int(np.argmin(es))
        r_eq, e_min = rows[idx]["R"], float(es[idx])
        if 0 < idx < len(rows) - 1:
            xx = np.array([rows[idx - 1]["R"], rows[idx]["R"], rows[idx + 1]["R"]])
            yy = np.array([es[idx - 1], es[idx], es[idx + 1]])
            c = np.polyfit(xx, yy, 2)
            r_eq = float(-c[1] / (2.0 * c[0]))
            e_min = float(np.polyval(c, r_eq))
        return r_eq, e_min, e_diss - e_min

    print(f"\n[wells]  (vs in-basis atoms; experiment: R_eq = {R_EQ_EXP} a0, "
          f"D_e = {D_E_EXP_EV} eV)")
    summary = {}
    d_e_fci = None
    for key in keys:
        r_eq, e_min, d_e = well(key)
        interior = r_eq < grid[-4]
        binds = d_e > 1e-4 and interior
        print(f"  {key:30s} R_eq = {r_eq:6.3f} a0  E_min = {e_min:11.6f} Ha  "
              f"D_e = {d_e:9.5f} Ha = {d_e * HARTREE_TO_EV:6.3f} eV  "
              f"{'BINDS' if binds else 'UNBOUND'}")
        if key == "FCI (91 dets)":
            d_e_fci = d_e
        summary[key] = {"R_eq": r_eq, "E_min": e_min, "D_e_Ha": d_e,
                        "binds": bool(binds)}

    if d_e_fci and d_e_fci > 0:
        print("\n[compactness]  fraction of in-basis FCI binding captured:")
        for key in keys:
            if key == "FCI (91 dets)":
                continue
            print(f"  {key:30s} {summary[key]['D_e_Ha'] / d_e_fci:8.1%}")

    # --- pre-registered gates ---
    noci3 = summary["cov+ionH (3 dets)"]
    frac3 = noci3["D_e_Ha"] / d_e_fci if d_e_fci and d_e_fci > 0 else 0.0
    g1 = noci3["binds"]
    g2 = abs(noci3["R_eq"] - R_EQ_EXP) / R_EQ_EXP <= 0.15
    g3 = 0.3 <= noci3["D_e_Ha"] * HARTREE_TO_EV <= 3.0
    g4 = frac3 >= 0.85
    print("\n[gates]  (pre-registered, debug/noci_sandbox_notes.md N4)")
    print(f"  G1 interior minimum (NOCI-3):        {'PASS' if g1 else 'FAIL'}")
    print(f"  G2 R_eq within 15% of {R_EQ_EXP} a0:     "
          f"{'PASS' if g2 else 'FAIL'}  ({noci3['R_eq']:.3f} a0, "
          f"{100 * (noci3['R_eq'] - R_EQ_EXP) / R_EQ_EXP:+.1f}%)")
    print(f"  G3 D_e sane (0.3-3 eV; exp 1.96):    "
          f"{'PASS' if g3 else 'FAIL'}  ({noci3['D_e_Ha'] * HARTREE_TO_EV:.3f} eV)")
    print(f"  G4 NOCI-3 >= 85% of in-basis FCI:    "
          f"{'PASS' if g4 else 'FAIL'}  ({frac3:.1%})")
    verdict = "GO (all gates pass)" if all((g1, g2, g3, g4)) else \
        "PARTIAL/FAIL -- see gate lines"
    print(f"  VERDICT: {verdict}")
    print("\n[diagnosis] corpus anchors: balanced/composed NaH PES monotonically")
    print("  descends, no interior minimum at any tested configuration (W1e;")
    print("  Sprint B.1 explicit-core HF, R3-B DMRG, F1-F6, kwarg sweep).")
    print("  A binding verdict here = the wall is the Hamiltonian specification")
    print("  (surrogate/missing cross-center physics), not the correlation")
    print("  treatment -- confirmed by construction with genuine integrals.")

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "data", "noci_nah_probe_results.json")
    with open(out_path, "w") as fh:
        json.dump({"zetas_na": z_na, "zeta_h": z_h, "E_Na": e_na, "E_H": e_h,
                   "rows": rows, "wells": summary,
                   "gates": {"G1": bool(g1), "G2": bool(g2), "G3": bool(g3),
                             "G4": bool(g4)},
                   "exp": {"R_eq_a0": R_EQ_EXP, "D_e_eV": D_E_EXP_EV},
                   "checks": {"rot_invariance": abs(e_fci_z - e_fci_rot),
                              "span_identity": abs(e_span + vnn - e_fci_z)}},
                  fh, indent=1)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()
