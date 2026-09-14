"""Probe: the contraction seam behind Paper 60's two-centre metric.

QUESTION (PI, 2026-09-12).  Where does overcompleteness enter *mathematically*,
and where does the transcendental enter along the continuous (decompactified)
direction?

HYPOTHESIS under test.  The one-centre content of the Fock-sphere basis is
SO(4) zonal (Gegenbauer lambda=1); the cross-centre symbol j0(pR) is the E(3)
zonal spherical function.  The two are joined by the Inonu-Wigner contraction
SO(4) -> E(3), realised analytically as the Mehler-Heine limit.  The two-centre
degeneracy sits at p = 0 -- the trivial character of the translation group --
which is also the pole at which S^3 looks locally flat.  So the transcendental
(Bessel) enters the molecular metric exactly at the contraction.

PARAMETRISATION (matches geovac/sturmian_sigma_law.py).  chi in [0,pi] with
p = k cot(chi/2), so chi = pi is p = 0 (IR, the conditioning pole) and chi = 0
is p = infinity (UV, the chirp/locality pole).  theta = pi - chi is the standard
Fock polar angle measured from the p = 0 pole.

Five tests, A-E.  A-C are identities (should be exact/classical); D locates the
degeneracy in the contraction window; E is the new falsifiable prediction --
that the law is independent of a radial weighting potential.

Run:  python debug/p60_contraction_seam_probe.py
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy.special import eval_gegenbauer, roots_legendre

from geovac.sturmian_sigma_law import sw_cross_block, COLLAPSE_CONSTANT

OUT: dict = {}


def _j0(x: np.ndarray) -> np.ndarray:
    """Spherical Bessel j0 = sin(x)/x, safe at 0."""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    nz = x != 0.0
    out[nz] = np.sin(x[nz]) / x[nz]
    return out


# ---------------------------------------------------------------- A
def test_A_sine_basis_is_zonal_gegenbauer() -> None:
    """sin(a chi) == sin(chi) * C_{a-1}^{(1)}(cos chi), exactly.

    The sine basis Paper 60 works in IS the S^3 zonal (Gegenbauer lambda=1)
    basis in half-density form -- the sin^2(chi) measure of S^3 absorbed as
    sin(chi) into each function.  This is what makes the one-centre block the
    identity.
    """
    chi = np.linspace(1e-6, np.pi - 1e-6, 4001)
    worst = 0.0
    for a in range(1, 25):
        lhs = np.sin(a * chi)
        rhs = np.sin(chi) * eval_gegenbauer(a - 1, 1.0, np.cos(chi))
        worst = max(worst, float(np.max(np.abs(lhs - rhs))))
    OUT["A_sine_eq_zonal_gegenbauer_maxabs"] = worst
    print(f"[A] max |sin(a chi) - sin(chi) C^(1)_(a-1)(cos chi)|  = {worst:.3e}"
          f"   (a = 1..24)")


# ---------------------------------------------------------------- B
def test_B_mehler_heine_contraction() -> None:
    """n^-1 C_{n-1}^{(1)}(cos(z/n)) -> j0(z):  SO(4) zonal -> E(3) zonal.

    This IS the Inonu-Wigner contraction of the sphere to the plane, in the
    zonal sector: hold the geodesic arc z = n*theta fixed while the effective
    radius n -> infinity.  Measure the rate.
    """
    z = np.linspace(0.05, 10.0, 400)
    target = _j0(z)
    rows = []
    for n in (10, 20, 40, 80, 160, 320, 640):
        approx = eval_gegenbauer(n - 1, 1.0, np.cos(z / n)) / n
        err = float(np.max(np.abs(approx - target)))
        rows.append((n, err))
    ns = np.array([r[0] for r in rows], float)
    es = np.array([r[1] for r in rows], float)
    slope = float(np.polyfit(np.log(ns), np.log(es), 1)[0])
    OUT["B_mehler_heine"] = {"rows": rows, "log_slope": slope}
    print("[B] Mehler-Heine  n^-1 C^(1)_(n-1)(cos(z/n)) -> j0(z),  z in (0,10]")
    for n, e in rows:
        print(f"       n = {n:4d}   max err = {e:.3e}")
    print(f"     fitted exponent = {slope:.3f}   (expect -2)")


# ---------------------------------------------------------------- C
def test_C_plane_wave_average_is_j0() -> None:
    """(1/4pi) \\int e^{i p.R} dOmega_p = j0(pR):  the E(3) zonal function.

    The Shibuya-Wulfman operator is multiplication by the translation phase on
    the Fock sphere; its angular average over the momentum direction is the
    symbol.  Gauss-Legendre in cos(angle); the phi integral is trivial.
    """
    x, w = roots_legendre(200)  # nodes in cos(angle) on [-1,1]
    worst = 0.0
    vals = []
    for pr in (0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0):
        avg = 0.5 * np.sum(w * np.exp(1j * pr * x))
        exact = float(_j0(np.array([pr]))[0])
        vals.append((pr, float(avg.real), exact, float(abs(avg.imag))))
        worst = max(worst, abs(avg.real - exact), abs(avg.imag))
    OUT["C_plane_wave_average"] = {"rows": vals, "max_abs_err": worst}
    print(f"[C] (1/4pi) int exp(i p.R) dOmega == j0(pR):  max err = {worst:.3e}")


# ---------------------------------------------------------------- D
def test_D_degeneracy_lives_in_the_contraction_window() -> None:
    """The near-null direction concentrates at theta ~ 1/n, the M-H window.

    Top singular vector v of the cross block C; its profile in the Fock polar
    angle theta = pi - chi is f(theta) = sum_a v_a sin(a chi) evaluated on the
    theta grid.  Report the L2 centroid <theta> and width sqrt(<theta^2>).
    If the degeneracy is a contraction-window phenomenon both scale as 1/n.
    """
    s = 2.0
    chi = np.linspace(1e-9, np.pi, 200001)
    theta = np.pi - chi
    rows = []
    for n in (20, 40, 80, 160, 320):
        C = sw_cross_block(s, n, M=300001)
        U, sig, _ = np.linalg.svd(C)
        v = U[:, 0]
        a = np.arange(1, n + 1)
        f = v @ np.sin(np.outer(a, chi))
        w = f ** 2
        norm = np.trapezoid(w, chi)
        mean = float(np.trapezoid(w * theta, chi) / norm)
        rms = float(np.sqrt(np.trapezoid(w * theta ** 2, chi) / norm))
        rows.append((n, float(1.0 - sig[0]), mean, rms, mean * n, rms * n))
    OUT["D_window"] = {
        "s": s,
        "cols": ["n", "1-sigma_max", "<theta>", "rms_theta", "n<theta>", "n*rms"],
        "rows": rows,
    }
    print(f"[D] near-null direction, profile in the Fock polar angle (s = {s})")
    print("       n   1-sigma_max     <theta>    rms(theta)   n<theta>   n*rms")
    for n, d, m, r, mn, rn in rows:
        print(f"    {n:4d}   {d:.3e}    {m:.5f}    {r:.5f}    {mn:7.3f}  {rn:7.3f}")


# ---------------------------------------------------------------- E
def _weighted_blocks(s: float, n: int, W, M: int = 300001):
    """A_ab = (2/pi) int sin(a chi) sin(b chi) W dchi;  B_ab the same with W*j0.

    A is the one-centre (intra) block under a radial weight W on the Fock
    sphere; B is the cross-centre block.  W = 1 must reproduce (I, C).
    """
    chi = np.linspace(1e-9, np.pi, M)
    wq = np.full(M, chi[1] - chi[0])
    wq[0] *= 0.5
    wq[-1] *= 0.5
    Wv = W(chi)
    sym = _j0(s / np.tan(chi / 2.0))
    a = np.arange(1, n + 1)
    smat = np.sin(np.outer(a, chi))
    A = (2.0 / np.pi) * (smat * (wq * Wv)) @ smat.T
    B = (2.0 / np.pi) * (smat * (wq * Wv * sym)) @ smat.T
    return A, B


def _gen_sigma_max(A: np.ndarray, B: np.ndarray) -> float:
    """Largest generalized singular value: spectrum of A^-1/2 B A^-1/2."""
    ev, Q = np.linalg.eigh(A)
    ev = np.maximum(ev, 1e-300)
    Ainv_half = Q @ np.diag(ev ** -0.5) @ Q.T
    return float(np.max(np.abs(np.linalg.eigvalsh(Ainv_half @ B @ Ainv_half))))


def test_E_V0_independence() -> None:
    """PREDICTION: the law is independent of a smooth positive radial weight.

    The generalized symbol is (W * j0)/W = j0, so sigma(p) does not see W at
    all.  Predict: exponent 2 AND the constant pi^2/24 survive for any W smooth
    and nonzero at the degeneracy (chi = pi).  CONTROL: a W that VANISHES there
    should move the constant -- if it does not, the test is not sensitive and
    the agreement above proves nothing.
    """
    s = 2.0
    weights = {
        "W = 1 (SW, the reference)": lambda c: np.ones_like(c),
        "W = 1 + 0.8 cos(chi)": lambda c: 1.0 + 0.8 * np.cos(c),
        "W = 2 + sin(chi)": lambda c: 2.0 + np.sin(c),
        "W = exp(-chi)": lambda c: np.exp(-c),
        "CONTROL W = 1 + cos(chi)  [vanishes at chi=pi]":
            lambda c: 1.0 + np.cos(c),
    }
    print(f"[E] weighting-potential independence  (s = {s}); "
          f"target collapse constant pi^2/24 = {COLLAPSE_CONSTANT:.6f}")
    print("                                              "
          "n=40      n=80     n=160   exponent")
    table = {}
    for label, W in weights.items():
        vals, ns = [], (40, 80, 160)
        for n in ns:
            A, B = _weighted_blocks(s, n, W)
            smax = _gen_sigma_max(A, B)
            vals.append((1.0 - smax) * (n / s) ** 2)
        d = [(1.0 - _gen_sigma_max(*_weighted_blocks(s, n, W))) for n in ns]
        expo = float(np.polyfit(np.log(np.array(ns, float)),
                                np.log(np.array(d)), 1)[0])
        table[label] = {"collapse": vals, "exponent": expo}
        print(f"    {label:46s} {vals[0]:.5f}  {vals[1]:.5f}  {vals[2]:.5f}"
              f"   {expo:+.3f}")
    OUT["E_weight_independence"] = {"s": s, "target": COLLAPSE_CONSTANT,
                                    "table": table}


def main() -> None:
    test_A_sine_basis_is_zonal_gegenbauer()
    print()
    test_B_mehler_heine_contraction()
    print()
    test_C_plane_wave_average_is_j0()
    print()
    test_D_degeneracy_lives_in_the_contraction_window()
    print()
    test_E_V0_independence()
    os.makedirs("debug/data", exist_ok=True)
    with open("debug/data/p60_contraction_seam.json", "w") as fh:
        json.dump(OUT, fh, indent=2, default=float)
    print("\nwrote debug/data/p60_contraction_seam.json")


if __name__ == "__main__":
    main()
