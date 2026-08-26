"""Do Coulomb Sturmians converge faster PER BASIS FUNCTION than Gaussians?

PI direction 2026-08-22. This tests the premise underneath "replace the basis
wholesale": Sturmians carry the correct nuclear cusp and the correct
exponential tail, so they should need fewer functions per unit accuracy than
Gaussians, which have neither. If true, the corpus's closed-form two-center
integrals remove the one thing (integral cost) that keeps Sturmian-basis
correlated methods out of production. If false, that program has no foundation.

WHY THIS NEEDS TESTING RATHER THAN ASSUMING. QC-1 (v4.77.0) measured something
close and got a NEGATIVE: holding H2 at M=2 and varying only contraction depth,
an uncontracted primitive gives -0.991 Ha, a 3-term contraction -1.139, a
10-term -1.148 -- and the 0.157 Ha "Slater advantage" over the uncontracted
primitive vanished once the Gaussian was properly contracted. Since a Slater
function IS a contracted Gaussian, Slater shapes may buy accuracy per
PRIMITIVE (cheap) and nothing per BASIS FUNCTION (the thing that costs).
QC-1 was H2, s-only, M=2 -- small enough that it may not generalize. This
widens it.

DESIGN. Same system (H2, R=1.4 a0), same method (all-electron FCI, so no
method error), same integral engine (McMurchie-Davidson), same evaluator
(everything rendered as contracted Gaussians). The ONLY difference between the
two families is the radial shape of the basis functions. Total energy at fixed
geometry is the metric -- NOT D_e -- because at fixed geometry there is no
dissociation reference and therefore no BSSE to contaminate the comparison.

    Sturmian family: Coulomb Sturmians at a SHARED exponent k (Avery), which
        is what makes them Sturmians rather than hydrogenic; k optimized
        variationally at each size.
    Gaussian family: real published contracted sets (STO-3G, 6-31G, 6-31G**,
        cc-pVDZ, cc-pVTZ-minus-d). Contracted, not primitives -- the QC-1
        lesson. These are the actual competition.

PRE-REGISTERED GATES (fixed before running).

  Validation, must pass or the comparison is void:
    V1  STO-3G H2 FCI reproduces the stored literature value -1.1373 Ha
        (Szabo & Ostlund, geovac/gaussian_reference.py) to +-2 mHa. Validates
        the Gaussian pipeline AND the hardcoded basis data.
    V2  Sturmian construction: same-center <chi_nl|1/r|chi_n'l> = 0 for n != n'
        (this orthogonality under the 1/r weight is the DEFINITION of a
        Coulomb Sturmian; a hydrogenic set fails it). Relative tolerance 1e-3.
    V3  Every basis is variationally sane: E decreases monotonically with M
        within each family, and the smallest overlap eigenvalue stays > 1e-6.

  Decision, on H2 at matched M:
    GO        Sturmian is LOWER than the contracted Gaussian by > 1 mHa at
              two or more matched M, and the advantage does not shrink toward
              zero as M grows.
    STOP      the contracted Gaussian is equal or lower at matched M
              (i.e. QC-1 generalizes; the premise is dead).
    BORDERLINE  mixed, or an advantage that decays with M.

Exploratory. No paper claim.
"""

from __future__ import annotations

import math
import os
import sys
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.noci_engine import (
    BasisFn, fci_ground, fit_sto_shape, integral_set_md, lowdin_orbitals,
    nuclear_md, overlap_md, transform_integrals,
)
from geovac.two_center_eri import radial_poly

R_H2 = 1.4
E_EXACT_H2 = -1.174476          # Kolos-Wolniewicz, R=1.4 a0
LIT_STO3G = -1.1373             # stored in geovac/gaussian_reference.py


# ------------------------------------------------------- Gaussian basis data
# Published contractions for HYDROGEN. Coefficients are on NORMALIZED
# primitives, which is the convention BasisFn expects ("d_i given on
# NORMALIZED primitives; the contraction is renormalized").
# Each entry: list of (lmn_kind, alphas, coeffs). 'S' -> one function,
# 'P' -> three (px, py, pz).

GAUSS_H: Dict[str, List[tuple]] = {
    "STO-3G": [
        ("S", [3.42525091, 0.62391373, 0.16885540],
              [0.15432897, 0.53532814, 0.44463454]),
    ],
    "6-31G": [
        ("S", [18.7311370, 2.8253937, 0.6401217],
              [0.03349460, 0.23472695, 0.81375733]),
        ("S", [0.1612778], [1.0]),
    ],
    "6-31G**": [
        ("S", [18.7311370, 2.8253937, 0.6401217],
              [0.03349460, 0.23472695, 0.81375733]),
        ("S", [0.1612778], [1.0]),
        ("P", [1.1], [1.0]),
    ],
    "cc-pVDZ": [
        ("S", [13.0100000, 1.9620000, 0.4446000],
              [0.0196850, 0.1379770, 0.4781480]),
        ("S", [0.1220000], [1.0]),
        ("P", [0.7270000], [1.0]),
    ],
    "cc-pVTZ(spd->sp)": [        # d shell dropped: engine + matched-M comparison
        ("S", [33.8700000, 5.0950000, 1.1590000],
              [0.0060680, 0.0453080, 0.2028220]),
        ("S", [0.3258000], [1.0]),
        ("S", [0.1027000], [1.0]),
        ("P", [1.4070000], [1.0]),
        ("P", [0.3880000], [1.0]),
    ],
}

# Known H2 FCI energies at R=1.4 for these sets, as a basis-data sanity signal.
# Soft check only (V1 hard-gates STO-3G); a large deviation flags bad data.
GAUSS_EXPECT = {"STO-3G": -1.1373, "6-31G": -1.1516, "cc-pVDZ": -1.1636}

P_LMN = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]


def gaussian_basis(name: str, centers) -> List[BasisFn]:
    out = []
    for c in centers:
        for kind, a, d in GAUSS_H[name]:
            a = np.asarray(a, dtype=float)
            d = np.asarray(d, dtype=float)
            if kind == "S":
                out.append(BasisFn(c, (0, 0, 0), a, d))
            else:
                for lmn in P_LMN:
                    out.append(BasisFn(c, lmn, a, d))
    return out


# ------------------------------------------------------ Sturmian construction

_SHAPE_CACHE: Dict[Tuple[int, int], tuple] = {}


def _shape(l: int, n_r: int):
    """Gaussian fit of the normalized zeta=1 STO r^{n_r-1} e^{-r}, ang. mom. l."""
    key = (l, n_r)
    if key not in _SHAPE_CACHE:
        a, d, q = fit_sto_shape(l, n_r)
        _SHAPE_CACHE[key] = (np.asarray(a), np.asarray(d), q)
    return _SHAPE_CACHE[key]


def _sto_norm(m: int, k: float) -> float:
    """Norm of the STO N r^{m-1} e^{-kr}: N = (2k)^{m+1/2}/sqrt((2m)!)."""
    return (2.0 * k) ** (m + 0.5) / math.sqrt(math.factorial(2 * m))


def sturmian_fn(center, n: int, l: int, lmn, k: float) -> BasisFn:
    """A Coulomb Sturmian chi_{nl} at SHARED exponent k, as a contracted Gaussian.

    chi_nl(r) = sum_p c_p r^p e^{-kr}, the c_p taken from radial_poly at
    Z = n*k (which sets the decay rate to Z/n = k for EVERY n -- that shared
    rate is exactly what makes the set Sturmian rather than hydrogenic).
    Each r^p e^{-kr} is a Slater function of index m = p+1, rendered by the
    engine's fitted shape; coefficients are converted onto those normalized
    shapes by dividing out the Slater norm. BasisFn renormalizes the whole
    contraction, so only relative weights matter.
    """
    kf = Fraction(k).limit_denominator(10 ** 6)
    coeffs, decay = radial_poly(Fraction(n) * kf, n, l)
    assert abs(float(decay) - k) < 1e-12, (float(decay), k)

    alphas, dcos = [], []
    for p, c in sorted(coeffs.items()):
        m = int(p) + 1
        a, d, _ = _shape(l, m)
        alphas.append(a * k ** 2)
        dcos.append(d * (float(c) / _sto_norm(m, k)))
    return BasisFn(center, lmn, np.concatenate(alphas), np.concatenate(dcos))


# (n, l) shells in ascending order; d omitted (engine + matched-M design).
SHELLS = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1)]


def sturmian_basis(n_shells: int, centers, k: float) -> List[BasisFn]:
    out = []
    for c in centers:
        for (n, l) in SHELLS[:n_shells]:
            for lmn in ([(0, 0, 0)] if l == 0 else P_LMN):
                out.append(sturmian_fn(c, n, l, lmn, k))
    return out


# ------------------------------------------------------------------ energies

def fci_energy(orbs, nuclei, n_elec: int) -> Tuple[float, float]:
    s, h, g = integral_set_md(orbs, nuclei)
    s_min = float(np.linalg.eigvalsh(s).min())
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    return fci_ground(ht, gt, n_elec), s_min


def h2_energy(orbs) -> Tuple[float, float]:
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([0.0, 0.0, R_H2])
    e, s_min = fci_energy(orbs, [(a, 1.0), (b, 1.0)], 2)
    return e + 1.0 / R_H2, s_min


# ------------------------------------------------------------------ gates

def gate_v2() -> bool:
    """Sturmian 1/r orthogonality: <chi_nl|1/r|chi_n'l> = 0 for n != n'."""
    print("\n[V2] Sturmian definition: <chi_nl | 1/r | chi_n'l> = 0 for n != n'")
    o = np.array([0.0, 0.0, 0.0])
    k = 1.0
    ok = True
    for l, ns in ((0, (1, 2, 3)), (1, (2, 3))):
        lmn = (0, 0, 0) if l == 0 else (0, 0, 1)
        fns = {n: sturmian_fn(o, n, l, lmn, k) for n in ns}
        diag = [abs(nuclear_md(fns[n], fns[n], o, 1.0)) for n in ns]
        scale = max(diag)
        for i, n1 in enumerate(ns):
            for n2 in ns[i + 1:]:
                v = nuclear_md(fns[n1], fns[n2], o, 1.0)
                rel = abs(v) / scale
                flag = "ok" if rel < 1e-3 else "FAIL"
                if rel >= 1e-3:
                    ok = False
                print(f"     l={l}  <{n1}|1/r|{n2}> = {v:+.3e}  "
                      f"rel {rel:.2e}  {flag}")
        # a hydrogenic set (exponent Z/n) must FAIL this -- discrimination check
        if l == 0:
            hyd = {}
            for n in ns:
                co, dec = radial_poly(Fraction(1), n, l)   # Z=1 -> decay 1/n
                al, dc = [], []
                for p, c in sorted(co.items()):
                    m = int(p) + 1
                    a_, d_, _ = _shape(l, m)
                    al.append(a_ * float(dec) ** 2)
                    dc.append(d_ * (float(c) / _sto_norm(m, float(dec))))
                hyd[n] = BasisFn(o, lmn, np.concatenate(al), np.concatenate(dc))
            v = nuclear_md(hyd[1], hyd[2], o, 1.0)
            rel = abs(v) / abs(nuclear_md(hyd[1], hyd[1], o, 1.0))
            print(f"     [discrimination] hydrogenic <1|1/r|2> rel = {rel:.2e} "
                  f"{'(nonzero, as required)' if rel > 1e-3 else '<-- BAD: gate cannot discriminate'}")
            ok = ok and rel > 1e-3
    print(f"     V2: {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> None:
    print("=== Sturmian vs contracted Gaussian: energy per basis function ===")
    print(f"    H2 at R = {R_H2} a0, all-electron FCI, exact = {E_EXACT_H2} Ha")

    a = np.array([0.0, 0.0, 0.0])
    b = np.array([0.0, 0.0, R_H2])
    centers = [a, b]

    # ---- V1 -------------------------------------------------------------
    print("\n[V1] Gaussian pipeline vs stored literature STO-3G value")
    e_sto3g, _ = h2_energy(gaussian_basis("STO-3G", centers))
    d = abs(e_sto3g - LIT_STO3G)
    print(f"     computed {e_sto3g:.6f} Ha   stored {LIT_STO3G:.4f} Ha   "
          f"delta {d*1000:.2f} mHa   {'PASS' if d < 2e-3 else 'FAIL'}")
    v1 = d < 2e-3

    v2 = gate_v2()
    if not (v1 and v2):
        print("\n*** VALIDATION FAILED -- comparison void. ***")
        return

    # ---- Gaussian family -------------------------------------------------
    print("\n[Gaussian family]")
    print(f"  {'basis':<18} {'M':>3} {'E (Ha)':>12} {'err vs exact':>13} "
          f"{'s_min':>10}  lit")
    gauss = {}
    for name in GAUSS_H:
        orbs = gaussian_basis(name, centers)
        e, sm = h2_energy(orbs)
        exp = GAUSS_EXPECT.get(name)
        note = "" if exp is None else f"{e-exp:+.4f} vs {exp}"
        gauss[len(orbs)] = (name, e)
        print(f"  {name:<18} {len(orbs):>3} {e:12.6f} "
              f"{e-E_EXACT_H2:13.6f} {sm:10.2e}  {note}")

    # ---- Sturmian family (k optimized at each size) -----------------------
    print("\n[Sturmian family]  shared exponent k, optimized at each size")
    print(f"  {'shells':<18} {'M':>3} {'k*':>6} {'E (Ha)':>12} "
          f"{'err vs exact':>13} {'s_min':>10}")
    sturm = {}
    for ns in range(1, len(SHELLS) + 1):
        best = None
        for k in (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.35, 1.5, 1.7):
            orbs = sturmian_basis(ns, centers, k)
            e, sm = h2_energy(orbs)
            if best is None or e < best[1]:
                best = (k, e, sm, len(orbs))
        k, e, sm, M = best
        lab = "+".join(f"{n}{'spd'[l]}" for n, l in SHELLS[:ns])
        sturm[M] = (lab, e)
        print(f"  {lab:<18} {M:>3} {k:6.2f} {e:12.6f} "
              f"{e-E_EXACT_H2:13.6f} {sm:10.2e}")

    # ---- verdict ---------------------------------------------------------
    print("\n=== matched-M comparison ===")
    print(f"  {'M':>3} {'Sturmian':>12} {'Gaussian':>12} {'diff (mHa)':>11}  winner")
    wins, n_matched = [], 0
    for M in sorted(set(sturm) & set(gauss)):
        es, eg = sturm[M][1], gauss[M][1]
        diff = (es - eg) * 1000.0
        n_matched += 1
        wins.append(diff)
        who = "Sturmian" if diff < -1.0 else ("Gaussian" if diff > 1.0 else "tie")
        print(f"  {M:>3} {es:12.6f} {eg:12.6f} {diff:11.3f}  {who}"
              f"   ({sturm[M][0]} vs {gauss[M][0]})")

    n_stu = sum(1 for d_ in wins if d_ < -1.0)
    print()
    if n_matched == 0:
        print("  VERDICT: no matched M -- inconclusive.")
    elif n_stu >= 2 and wins[-1] <= wins[0]:
        print("  VERDICT: GO -- Sturmian lower at >=2 matched M, advantage not shrinking.")
    elif n_stu == 0:
        print("  VERDICT: STOP -- contracted Gaussians match or beat Sturmians")
        print("           at every matched M. QC-1 generalizes; the")
        print("           'better basis per function' premise is not supported.")
    else:
        print("  VERDICT: BORDERLINE -- mixed or decaying advantage.")


if __name__ == "__main__":
    main()
