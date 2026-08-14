"""Step 1: the first end-to-end molecule computed with the native engine.

Everything validated so far has been PER-INTEGRAL, against another integral
routine. Nobody has built a full tensor with this machinery, run a solver, and
compared an energy. Until that happens, "the engine works" is a claim about
integrals, not about chemistry.

WHAT IS NATIVE HERE, and what is not. All three ingredients come from exact
GeoVac machinery -- no Gaussian fits anywhere in the native path:

  S   two-center overlap   exact Mulliken/Ruedenberg (N3b `two_center_UV`),
                           value = e^{-p}(U e^{q} + V e^{-q}), U/V exact rationals
  h   one-electron         cross block by the hydrogenic eigen-trick
                           h_ab = E_b S_ab - Z_A <a|1/r_A|b>  (N3b);
                           diagonal by E_a - Z_B <a|1/r_B|a>, the nuclear
                           attraction taken from this arc's own V_L_radial
  g   two-electron         THIS ARC. All four classes:
                             (AA|AA),(BB|BB)  one-center, analytic for 1s
                             (AA|BB)          aabb_closed_form      [1c]
                             (AA|AB) etc.     hybrid_closed_form    [2]
                             (AB|AB)          exchange_value,
                                              exact_xi=True         [3b+3e]

N3b's own docstring named the two-center ERI engine as the piece it lacked
("that engine is N4's build"). This is that engine, so the two halves compose
into a complete native calculation for the first time.

REFERENCE. The identical basis evaluated through the Gaussian-fit pipeline
(10-primitive fits, <fit|STO> = 0.999999998). Any disagreement beyond the fit
floor is a real discrepancy in one of the two paths.

Run from repo root:  python debug/step1_native_molecule.py
"""

from __future__ import annotations

import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "debug"))

from geovac import noci_engine as GE  # noqa: E402
from geovac.two_center_eri import (  # noqa: E402
    V_L_radial, aabb_closed_form, exchange_value, hybrid_closed_form,
    multipole_decomposition, r_s,
)
from noci_n3b_census import norm_const, two_center_UV, uv_float, uv_is_zero  # noqa: E402

Z = Fraction(1)                     # H2: both centres Z = 1
RF = Fraction(7, 5)                 # R = 1.4 bohr
RV = float(RF)
ORB = (1, 0, 0)                     # 1s on each centre -> M = 2
CEN = ["A", "B"]


# ----------------------------------------------------------------- S and h

def _uv(kernel):
    """Normalised <1s_A | K | 1s_B> from the exact rational route."""
    p, q, U, V = two_center_UV(Z, 1, 0, Z, 1, 0, 0, RF, kernel=kernel)
    if uv_is_zero(U, V):
        return 0.0
    NN = norm_const(Z, 1, 0, 0) * norm_const(Z, 1, 0, 0)
    return uv_float(p, q, U, V) / float(sp.N(NN, 30))


def nuclear_attraction_same_centre():
    """<1s_A | 1/r_B | 1s_A> from this arc's own V_L_radial.

    The density |1s_A|^2 is one-centre, so its potential at distance R is the
    L = 0 multipole evaluated there. Cross-checked against the closed form
    (1/R)[1 - e^{-2R}(1 + R)] for a unit-exponent 1s.
    """
    terms = multipole_decomposition(Z, 1, 0, 0, Z, 1, 0, 0)
    L, M, g, rad, b = terms[0]
    val = float(sp.re(sp.N(
        g * V_L_radial(rad, b, L).subs(r_s, sp.nsimplify(RV))
        / sp.sqrt(4 * sp.pi), 30)))
    closed = (1.0 / RV) * (1 - np.exp(-2 * RV) * (1 + RV))
    assert abs(val - closed) < 1e-12, f"V_ne: {val} vs closed {closed}"
    return val


def build_S_h():
    S = np.eye(2)
    S[0, 1] = S[1, 0] = _uv(None)
    E1s = -float(Z) ** 2 / 2
    v_other = nuclear_attraction_same_centre()
    h = np.zeros((2, 2))
    h[0, 0] = h[1, 1] = E1s - float(Z) * v_other
    # cross block, ket trick: h_ab = E_b S_ab - Z_A <a|1/r_A|b>
    h[0, 1] = h[1, 0] = E1s * S[0, 1] - float(Z) * _uv("inv_ra")
    return S, h


# --------------------------------------------------------------------- g

def one_centre_1s():
    """(1s 1s | 1s 1s) on one centre = 5Z/8, exact."""
    return 5.0 * float(Z) / 8.0


def build_g(verbose=True):
    g = np.zeros((2, 2, 2, 2))
    seen = {}
    for p in range(2):
        for q in range(2):
            for r in range(2):
                for s in range(2):
                    cs = (CEN[p], CEN[q], CEN[r], CEN[s])
                    key = cs
                    if key in seen:
                        g[p, q, r, s] = seen[key]
                        continue
                    lhs, rhs = {cs[0], cs[1]}, {cs[2], cs[3]}
                    if len(set(cs)) == 1:                       # one-centre
                        v, tag = one_centre_1s(), "one-centre"
                    elif len(lhs) == 1 and len(rhs) == 1:       # (AA|BB)
                        ZA_, ZB_ = Z, Z
                        v = float(sp.re(sp.N(aabb_closed_form(
                            ZA_, ORB, ORB, ZB_, ORB, ORB,
                            sp.nsimplify(RV)), 30)))
                        tag = "(AA|BB) 1c"
                    elif len(lhs) == 1 or len(rhs) == 1:        # hybrid
                        v = float(sp.re(sp.N(hybrid_closed_form(
                            Z, ORB, ORB, ORB, Z, ORB, sp.nsimplify(RV)), 30)))
                        tag = "hybrid  2"
                    else:                                       # exchange
                        v = exchange_value(Z, ORB, ORB, Z, ORB, ORB, RV,
                                           tau_max=8, exact_xi=True)
                        tag = "exchange 3b+3e"
                    seen[key] = v
                    g[p, q, r, s] = v
                    if verbose:
                        print(f"    ({cs[0]}{cs[1]}|{cs[2]}{cs[3]})  {v: .12f}"
                              f"   [{tag}]")
    return g


# --------------------------------------------------------------- reference

def gaussian_reference():
    shapes = {}
    arr, dco, q = GE.fit_sto_shape(0, 1, n_gauss=10)
    shapes["1s"] = (arr, dco)
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., RV])
    orbs = [GE.sto_shape_basis(pa, "1s", float(Z), shapes, (0, 0, 0)),
            GE.sto_shape_basis(pb, "1s", float(Z), shapes, (0, 0, 0))]
    nuc = [(pa, float(Z)), (pb, float(Z))]
    return GE.integral_set_md(orbs, nuc) + (q,)


def solve(S, h, g):
    X = GE.lowdin_orbitals(S)
    ht, gt = GE.transform_integrals(X, h, g)
    return GE.fci_ground(ht, gt, 2) + float(Z) ** 2 / RV


def main() -> None:
    print("Step 1 -- first end-to-end molecule on the native engine\n")
    print(f"H2, R = {RV} bohr, 1s on each centre (M = 2), 2 electrons, FCI\n")

    t0 = time.time()
    S, h = build_S_h()
    print("  NATIVE S and h (exact rational / eigen-trick):")
    print(f"    S_AB = {S[0,1]: .12f}")
    print(f"    h_AA = {h[0,0]: .12f}    h_AB = {h[0,1]: .12f}")
    print("\n  NATIVE g (this arc, all four classes):")
    g = build_g()
    t_native = time.time() - t0

    Sg, hg, gg, fitq = gaussian_reference()
    print(f"\n  GAUSSIAN reference (10-primitive fits, <fit|STO> = {fitq:.9f}):")
    print(f"    S_AB = {Sg[0,1]: .12f}   (native - ref = {S[0,1]-Sg[0,1]:+.2e})")
    print(f"    h_AA = {hg[0,0]: .12f}   (native - ref = {h[0,0]-hg[0,0]:+.2e})")
    print(f"    h_AB = {hg[0,1]: .12f}   (native - ref = {h[0,1]-hg[0,1]:+.2e})")
    dg = np.max(np.abs(g - gg))
    print(f"    max |g_native - g_ref| = {dg:.2e}")

    e_nat, e_ref = solve(S, h, g), solve(Sg, hg, gg)
    print(f"\n  TOTAL ENERGY (FCI + V_NN)")
    print(f"    native  : {e_nat:.12f} Ha   [{t_native:.1f}s to build]")
    print(f"    gaussian: {e_ref:.12f} Ha")
    print(f"    diff    : {abs(e_nat-e_ref):.2e} Ha")
    print(f"\n  chemical accuracy is 1.6e-03 Ha; the fit floor here is ~1e-09.")
    verdict = ("agree within the fit floor" if abs(e_nat - e_ref) < 1e-6
               else "DISAGREE beyond the fit floor -- investigate")
    print(f"  => the two paths {verdict}.")


if __name__ == "__main__":
    main()
