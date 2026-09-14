"""Follow-up to p60_contraction_seam_probe.py test D.

D measured the near-null direction's rms spread in the Fock polar angle theta
and found n * rms -> ~3.13, rising.  If that limit is pi, then the whole
conditioning law is an identity in the contraction window:

    1 - sigma_max = <1 - j0> = (s^2/24) <theta^2> = (s^2/24) (pi^2/n^2),

i.e. the pi^2 of the Kac-Murdock-Szego constant is exactly the minimal
mean-square spread of a band-limited function, a TRUNCATION quantity with no
Bessel content.  Three tests:

  F1  Richardson-extrapolate n*rms for the SW near-null direction.
  F2  Minimise <theta^2> over the band directly -- no Bessel anywhere -- and
      check n^2 * min <theta^2> -> pi^2.  This is the KMS problem in position
      space; if F1 and F2 agree, the pi^2 is provably truncation-side.
  F3  Check the factorisation 1 - sigma_max = (s^2/24) <theta^2> numerically
      at several s, i.e. that the near-null vector lives where the quadratic
      approximation to 1 - j0 is exact.

Run:  python debug/p60_window_constant_probe.py
"""
from __future__ import annotations

import json
import os

import numpy as np

from geovac.sturmian_sigma_law import sw_cross_block

OUT: dict = {}
M_GRID = 400001


def _theta_moments(v: np.ndarray, n: int):
    """(<theta>, rms theta) of f(theta) = sum_a v_a sin(a chi), chi = pi-theta."""
    chi = np.linspace(1e-9, np.pi, M_GRID)
    theta = np.pi - chi
    a = np.arange(1, n + 1)
    f = v @ np.sin(np.outer(a, chi))
    w = f ** 2
    norm = np.trapezoid(w, chi)
    mean = float(np.trapezoid(w * theta, chi) / norm)
    ms = float(np.trapezoid(w * theta ** 2, chi) / norm)
    return mean, float(np.sqrt(ms)), ms


def _theta2_matrix(n: int) -> np.ndarray:
    """T_ab = (2/pi) int_0^pi sin(a chi) sin(b chi) (pi-chi)^2 dchi."""
    chi = np.linspace(0.0, np.pi, M_GRID)
    wq = np.full(M_GRID, chi[1] - chi[0])
    wq[0] *= 0.5
    wq[-1] *= 0.5
    th2 = (np.pi - chi) ** 2
    a = np.arange(1, n + 1)
    smat = np.sin(np.outer(a, chi))
    return (2.0 / np.pi) * (smat * (wq * th2)) @ smat.T


def test_F1_richardson() -> None:
    s = 2.0
    ns = (20, 40, 80, 160, 320, 640)
    rows = []
    for n in ns:
        U, sig, _ = np.linalg.svd(sw_cross_block(s, n, M=300001))
        mean, rms, ms = _theta_moments(U[:, 0], n)
        rows.append((n, float(1.0 - sig[0]), rms, rms * n, ms * n * n))
    print("[F1] SW near-null direction: does n*rms -> pi?")
    print("        n    1-sigma_max      rms(theta)    n*rms     n^2<theta^2>")
    for n, d, r, rn, mn in rows:
        print(f"     {n:5d}   {d:.4e}    {r:.6f}   {rn:.5f}    {mn:.5f}")
    a = np.array([r[3] for r in rows])
    rich = 2.0 * a[1:] - a[:-1]          # first-order Richardson in 1/n
    print(f"     Richardson (n*rms): {np.array2string(rich, precision=5)}")
    print(f"     pi = {np.pi:.5f}   pi^2 = {np.pi ** 2:.5f}")
    OUT["F1"] = {"rows": rows, "richardson_n_rms": rich.tolist()}


def test_F2_band_limited_minimum() -> None:
    """min <theta^2> over the band, with NO Bessel function present at all."""
    print("\n[F2] min <theta^2> over span{sin(a chi)}_{a<=n}  (no Bessel anywhere)")
    print("        n    min<theta^2>    n^2 * min      ratio to pi^2")
    rows = []
    for n in (20, 40, 80, 160, 320, 640):
        T = _theta2_matrix(n)
        lam = float(np.min(np.linalg.eigvalsh(T)))
        rows.append((n, lam, lam * n * n, lam * n * n / np.pi ** 2))
        print(f"     {n:5d}   {lam:.6e}   {lam * n * n:.6f}     "
              f"{lam * n * n / np.pi ** 2:.6f}")
    a = np.array([r[2] for r in rows])
    rich = 2.0 * a[1:] - a[:-1]
    print(f"     Richardson: {np.array2string(rich, precision=5)}    "
          f"pi^2 = {np.pi ** 2:.5f}")
    OUT["F2"] = {"rows": rows, "richardson": rich.tolist()}


def test_F3_factorisation() -> None:
    """1 - sigma_max  ==  (s^2/24) <theta^2>  on the near-null direction."""
    print("\n[F3] factorisation 1-sigma_max = (s^2/24)<theta^2>, near-null vector")
    print("        s      n    1-sigma_max     (s^2/24)<theta^2>     ratio")
    rows = []
    for s in (0.5, 1.0, 2.0, 4.0):
        for n in (80, 320):
            U, sig, _ = np.linalg.svd(sw_cross_block(s, n, M=300001))
            _, _, ms = _theta_moments(U[:, 0], n)
            lhs = float(1.0 - sig[0])
            rhs = (s ** 2 / 24.0) * ms
            rows.append((s, n, lhs, rhs, lhs / rhs))
            print(f"     {s:4.1f}  {n:5d}   {lhs:.6e}    {rhs:.6e}    "
                  f"{lhs / rhs:.6f}")
    OUT["F3"] = rows


def main() -> None:
    test_F1_richardson()
    test_F2_band_limited_minimum()
    test_F3_factorisation()
    os.makedirs("debug/data", exist_ok=True)
    with open("debug/data/p60_window_constant.json", "w") as fh:
        json.dump(OUT, fh, indent=2, default=float)
    print("\nwrote debug/data/p60_window_constant.json")


if __name__ == "__main__":
    main()
