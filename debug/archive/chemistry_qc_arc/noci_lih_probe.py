"""NOCI sandbox probe, step 2: LiH with genuine integrals and N-electron machinery.

Sandbox exploration (branch sandbox/noci, repo frozen at v4.76.0 -- no release path).
Follow-on to debug/noci_h2_probe.py (step 1, machinery validated on H2).

Design:
  - Basis: three contracted s-type orbitals -- Li 1s (STO-6G shape, zeta_1s), Li 2s
    (STO-3G 2s shape with its negative inner coefficient, zeta_2s), H 1s (STO-6G,
    zeta_H).  All integrals are exact closed forms over s Gaussians.  Exponents are
    chosen VARIATIONALLY ON THE ISOLATED ATOMS ONLY (no molecular tuning), so the
    molecular curve is a prediction of the atomic fragments.
  - States: non-orthogonal Slater determinants over spin-orbitals built directly on
    the atom-centered orbitals (no orthogonalization anywhere).  Matrix elements by
    explicit permutation expansion -- exact for any overlap pattern; S_IJ is
    cross-checked against det(M) every call.
  - Ladder: covalent (2 dets) -> + ionic Li+H- (3) -> + ionic Li-H+ (4) -> in-basis
    FCI (all C(6,4)=15 dets).  Loewdin-orbital truncations for comparison; the
    complete 15-det space must agree between representations (global span check).

Known limitations (disclosed): no Li 2p (no sp hybridization -> underbinds), no
diffuse function on H (Li+H- anion character underdescribed), minimal basis.
"""

from __future__ import annotations

import json
import math
from itertools import combinations, permutations
from typing import Dict, List, Sequence, Tuple

import numpy as np

from noci_h2_probe import STO3G_1S, STO6G_1S, Basis1s, integral_set, gen_eig_ground

HARTREE_TO_EV = 27.211386

# STO-3G 2s shape (zeta=1): standard 2sp-shell exponents, 2s contraction coefficients.
STO3G_2S = (
    np.array([0.9942030, 0.2310310, 0.0751386]),
    np.array([-0.09996723, 0.39951283, 0.70011547]),
)


def perm_sign(perm: Sequence[int]) -> int:
    inv = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inv += 1
    return -1 if inv % 2 else 1


def det_pair_elements(
    so_i: Sequence[Tuple[int, int]],
    so_j: Sequence[Tuple[int, int]],
    s: np.ndarray,
    h: np.ndarray,
    g: np.ndarray,
) -> Tuple[float, float]:
    """<D_I|D_J> and <D_I|H|D_J> for non-orthogonal determinants via permutation sums.

    Spin-orbitals are (spatial_index, spin) with spin in {0,1}.  Chemist eri
    g[p,q,r,s] = (pq|rs): electron 1 carries (p,q), electron 2 carries (r,s).
    """
    n = len(so_i)
    s_tot = 0.0
    h_tot = 0.0
    for perm in permutations(range(n)):
        sgn = perm_sign(perm)
        ov = []
        for k in range(n):
            p, sp = so_i[k]
            q, sq = so_j[perm[k]]
            ov.append(s[p, q] if sp == sq else 0.0)

        prod_all = 1.0
        for v in ov:
            prod_all *= v
        s_tot += sgn * prod_all

        for a in range(n):
            p, sp = so_i[a]
            q, sq = so_j[perm[a]]
            if sp != sq:
                continue
            rest = 1.0
            for k in range(n):
                if k != a:
                    rest *= ov[k]
            if rest != 0.0:
                h_tot += sgn * h[p, q] * rest

        for a in range(n):
            pa, spa = so_i[a]
            qa, sqa = so_j[perm[a]]
            if spa != sqa:
                continue
            for b in range(a + 1, n):
                pb, spb = so_i[b]
                qb, sqb = so_j[perm[b]]
                if spb != sqb:
                    continue
                rest = 1.0
                for k in range(n):
                    if k != a and k != b:
                        rest *= ov[k]
                if rest != 0.0:
                    h_tot += sgn * g[pa, qa, pb, qb] * rest
    return s_tot, h_tot


def det_overlap_check(
    so_i: Sequence[Tuple[int, int]], so_j: Sequence[Tuple[int, int]], s: np.ndarray
) -> float:
    m = np.zeros((len(so_i), len(so_j)))
    for a, (p, sp) in enumerate(so_i):
        for b, (q, sq) in enumerate(so_j):
            m[a, b] = s[p, q] if sp == sq else 0.0
    return float(np.linalg.det(m))


def noci_ground(
    dets: List[Sequence[Tuple[int, int]]], s: np.ndarray, h: np.ndarray, g: np.ndarray
) -> Tuple[float, float, float]:
    """Ground energy over a determinant list; returns (E, condS, max |S_perm - det(M)|)."""
    nd = len(dets)
    smat = np.zeros((nd, nd))
    hmat = np.zeros((nd, nd))
    worst = 0.0
    for i in range(nd):
        for j in range(i, nd):
            sij, hij = det_pair_elements(dets[i], dets[j], s, h, g)
            worst = max(worst, abs(sij - det_overlap_check(dets[i], dets[j], s)))
            smat[i, j] = smat[j, i] = sij
            hmat[i, j] = hmat[j, i] = hij
    e, cond = gen_eig_ground(hmat, smat)
    return e, cond, worst


UP, DN = 0, 1
LI1S, LI2S, H1S = 0, 1, 2


def build_orbitals(r: float, z1s: float, z2s: float, zh: float) -> Tuple[list, list]:
    pos_li = np.array([0.0, 0.0, 0.0])
    pos_h = np.array([0.0, 0.0, r])
    orbs = [
        Basis1s(pos_li, z1s, *STO6G_1S),
        Basis1s(pos_li, z2s, *STO3G_2S),
        Basis1s(pos_h, zh, *STO6G_1S),
    ]
    nuclei = [(pos_li, 3.0), (pos_h, 1.0)]
    return orbs, nuclei


def atom_energies(z1s: float, z2s: float, zh: float) -> Tuple[float, float]:
    """In-basis E(Li) (1s^2 2s doublet, single det) and E(H)."""
    pos = np.array([0.0, 0.0, 0.0])
    li_orbs = [Basis1s(pos, z1s, *STO6G_1S), Basis1s(pos, z2s, *STO3G_2S)]
    s, h, g = integral_set(li_orbs, [(pos, 3.0)])
    det_li = [(0, UP), (0, DN), (1, UP)]
    sij, hij = det_pair_elements(det_li, det_li, s, h, g)
    e_li = hij / sij

    h_orb = [Basis1s(pos, zh, *STO6G_1S)]
    s1, h1, _ = integral_set(h_orb, [(pos, 1.0)])
    e_h = h1[0, 0] / s1[0, 0]
    return e_li, e_h


def lowdin_transform(s: np.ndarray, h: np.ndarray, g: np.ndarray):
    w, v = np.linalg.eigh(s)
    x = v @ np.diag(w**-0.5) @ v.T
    ht = x @ h @ x
    gt = np.einsum("pi,qj,rk,sl,ijkl->pqrs", x, x, x, x, g)
    return np.eye(s.shape[0]), ht, gt


def main() -> None:
    # --- atomic exponent selection (variational, atoms only) ---
    zh = 1.0
    z1s_grid = [2.60, 2.69, 2.80]
    z2s_grid = [0.60, 0.65, 0.70, 0.75, 0.80, 0.90]
    best = None
    for z1 in z1s_grid:
        for z2 in z2s_grid:
            e_li, e_h = atom_energies(z1, z2, zh)
            if best is None or e_li < best[2]:
                best = (z1, z2, e_li, e_h)
    z1s, z2s, e_li, e_h = best
    e_diss = e_li + e_h
    print(f"[atoms] zeta_Li1s = {z1s}  zeta_Li2s = {z2s}  zeta_H = {zh}  (variational on atoms only)")
    print(f"[atoms] E(Li) = {e_li:.6f} Ha  (ROHF/STO-3G anchor ~ -7.3155)")
    print(f"[atoms] E(H)  = {e_h:.6f} Ha   ->  dissociation reference {e_diss:.6f} Ha\n")

    # --- determinant sets ---
    core = [(LI1S, UP), (LI1S, DN)]
    det_cov_a = core + [(LI2S, UP), (H1S, DN)]
    det_cov_b = core + [(H1S, UP), (LI2S, DN)]
    det_ion_h = core + [(H1S, UP), (H1S, DN)]      # Li+ H-
    det_ion_li = core + [(LI2S, UP), (LI2S, DN)]   # Li- H+

    ladders = {
        "cov (2 dets)": [det_cov_a, det_cov_b],
        "cov+ionH (3 dets)": [det_cov_a, det_cov_b, det_ion_h],
        "all 4 dets": [det_cov_a, det_cov_b, det_ion_h, det_ion_li],
    }

    spin_orbitals = [(p, sp) for p in (LI1S, LI2S, H1S) for sp in (UP, DN)]
    fci_dets = [list(c) for c in combinations(spin_orbitals, 4)]

    grid = [round(r, 2) for r in np.concatenate([np.arange(2.0, 6.01, 0.25), [8.0, 10.0, 12.0]])]
    rows: List[Dict[str, float]] = []
    worst_check = 0.0
    worst_span = 0.0
    for r in grid:
        orbs, nuclei = build_orbitals(r, z1s, z2s, zh)
        s, h, g = integral_set(orbs, nuclei)
        vnn = 3.0 / r
        row: Dict[str, float] = {"R": r, "S_2s_H": s[LI2S, H1S], "S_1s_H": s[LI1S, H1S]}

        for label, dets in ladders.items():
            e, cond, chk = noci_ground(dets, s, h, g)
            worst_check = max(worst_check, chk)
            row[label] = e + vnn
        e_fci, cond_fci, chk = noci_ground(fci_dets, s, h, g)
        worst_check = max(worst_check, chk)
        row["FCI (15 dets)"] = e_fci + vnn
        row["cond_fci"] = cond_fci

        st, ht, gt = lowdin_transform(s, h, g)
        e_lc, _, _ = noci_ground(ladders["cov (2 dets)"], st, ht, gt)
        row["Loewdin cov (2 dets)"] = e_lc + vnn
        e_l3, _, _ = noci_ground(ladders["cov+ionH (3 dets)"], st, ht, gt)
        row["Loewdin cov+ionH (3 dets)"] = e_l3 + vnn
        e_lf, _, _ = noci_ground(fci_dets, st, ht, gt)
        worst_span = max(worst_span, abs(e_lf - e_fci))
        rows.append(row)

    keys = ["cov (2 dets)", "cov+ionH (3 dets)", "all 4 dets", "FCI (15 dets)",
            "Loewdin cov (2 dets)", "Loewdin cov+ionH (3 dets)"]
    print("{:>5} {:>8}".format("R", "S_2s,H"), *[f"{k:>22}" for k in keys])
    for row in rows:
        print(
            "{:>5.2f} {:>8.4f}".format(row["R"], row["S_2s_H"]),
            *[f"{row[k]:>22.6f}" for k in keys],
        )

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
        return r_eq, e_min, e_diss - e_min

    print("\n[wells]  (vs in-basis atoms; experiment: R_eq = 3.015 a0, D_e = 0.0924 Ha = 2.515 eV)")
    summary = {}
    d_e_fci = None
    for key in keys:
        r_eq, e_min, d_e = well(key)
        interior = r_eq < grid[-4]
        binds = d_e > 1e-4 and interior
        line = (
            f"  {key:28s} R_eq = {r_eq:6.3f} a0  E_min = {e_min:10.6f} Ha  "
            f"D_e = {d_e:9.5f} Ha = {d_e * HARTREE_TO_EV:6.3f} eV  {'BINDS' if binds else 'UNBOUND'}"
        )
        if key == "FCI (15 dets)":
            d_e_fci = d_e
        print(line)
        summary[key] = {"R_eq": r_eq, "E_min": e_min, "D_e_Ha": d_e, "binds": bool(binds)}

    if d_e_fci and d_e_fci > 0:
        print("\n[compactness]  fraction of in-basis FCI binding captured:")
        for key in keys:
            if key == "FCI (15 dets)":
                continue
            frac = summary[key]["D_e_Ha"] / d_e_fci
            print(f"  {key:28s} {frac:8.1%}")

    print(f"\n[checks] permutation-sum vs det(M) overlap, worst   = {worst_check:.2e}")
    print(f"[checks] complete-space identity |FCI - FCI_Loewdin| = {worst_span:.2e}")
    print(f"[checks] worst FCI det-overlap condition number      = {max(r['cond_fci'] for r in rows):.2e}")

    print("\n[measurement ledger, NOQE-style, LiH]")
    for nd in (2, 3, 4):
        print(f"  {nd} dets -> S_IJ: {nd * (nd + 1) // 2}  H_IJ: {nd * (nd + 1) // 2}  "
              f"(each H_IJ costed against the NATIVE Pauli decomposition)")
    print("  corpus anchors: composed LiH native = 334 Pauli @ Q30; Loewdin retrofit")
    print("  inflated Pauli 17.9x (v4.73.1).  NOCI route: H stays at native sparsity;")
    print("  overlap cost -> n_det(n_det+1)/2 state-overlap circuits per PES point.")

    import os

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "noci_lih_probe_results.json")
    with open(out_path, "w") as fh:
        json.dump({"zetas": {"Li1s": z1s, "Li2s": z2s, "H": zh},
                   "E_Li": e_li, "E_H": e_h, "rows": rows, "wells": summary,
                   "worst_perm_check": worst_check, "worst_span": worst_span}, fh, indent=1)
    print("\n[saved] debug/data/noci_lih_probe_results.json")


if __name__ == "__main__":
    main()
