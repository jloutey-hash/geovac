"""Two-center shared-scale Coulomb--Sturmian integral engine (Paper 60).

The molecular-Sturmian / Shibuya--Wulfman integral layer, promoted from the Paper-60
sprint drivers to tracked code so the Paper-60 resource claims (notably the molecular
block-encoding 1-norm scaling, ``lambda ~ n_orb^2.2``) are backed by a regression suite
rather than transient ``debug/`` scripts.

An orbital is a triple ``(n, center, a)`` with principal number ``n``, ``center in
{'A','B'}`` (A at the origin, B at ``(0, 0, R)``), and radial decay ``a``.  The radial
part is the hydrogenic s-function

    R_n(r) = norm * exp(-a r) * L^1_{n-1}(2 a r),   int R^2 r^2 dr = 1,

so putting ``a = Q/n`` places the Goscinskian effective charge ``Q`` on principal number
``n``.  All integrals are evaluated numerically by a multipole (Legendre) expansion about
nucleus A; this is a *numerical* evaluator that complements the *exact closed forms* in
:mod:`geovac.two_center_eri` (against which it is validated to ~1e-4 --- see
:func:`validate`).

Provides, for orbitals ``i, j, k, l`` on any centers/charges:

  * :meth:`overlap`  ``S_ij   = <chi_i|chi_j>``
  * :meth:`nuclear`  ``V_ij   = <chi_i| (-1/r_A - 1/r_B) |chi_j>``
  * :meth:`kinetic`  ``T_ij   = <chi_i| -1/2 grad^2 |chi_j>``  (shared-scale Sturmian ODE)
  * :meth:`eri`      ``(ij|kl)= <chi_i chi_j | 1/r12 | chi_k chi_l>``  (chemist notation)
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.special import genlaguerre, eval_legendre
from scipy.integrate import cumulative_trapezoid

# numpy renamed trapz -> trapezoid; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# An orbital: (principal number n, center 'A'/'B', radial decay a).
Orb = Tuple[int, str, float]


class GoscinskianIntegrals:
    """Numerical two-center integral engine for shared-scale Coulomb--Sturmian s-orbitals.

    Parameters
    ----------
    R : float
        Internuclear separation (bohr); A at the origin, B at ``(0, 0, R)``.
    Lmax : int
        Highest Legendre multipole retained in the expansion about A.
    nr, nth : int
        Radial and angular (``cos theta``) grid sizes.
    rmax : float
        Radial grid extent (bohr).
    """

    def __init__(self, R: float, Lmax: int = 24, nr: int = 3000,
                 nth: int = 200, rmax: float = 60.0) -> None:
        self.R: float = R
        self.Lmax: int = Lmax
        self.r: np.ndarray = np.linspace(1e-5, rmax, nr)
        self.dr: float = self.r[1] - self.r[0]
        self.u: np.ndarray = np.sort(np.cos(np.linspace(0.0, np.pi, nth)))  # cos(theta) about A
        self.PL: List[np.ndarray] = [eval_legendre(L, self.u) for L in range(Lmax + 1)]
        self.RR, self.UU = np.meshgrid(self.r, self.u, indexing="ij")

    def _radial_norm(self, n: int, a: float) -> float:
        f = np.exp(-a * self.r) * genlaguerre(n - 1, 1)(2 * a * self.r)
        return 1.0 / np.sqrt(_trapz(f * f * self.r * self.r, self.r))

    def phi(self, orb: Orb) -> np.ndarray:
        """Full s-orbital ``/ sqrt(4 pi)`` sampled on the ``(r, u)`` grid about A."""
        n, c, a = orb
        norm = self._radial_norm(n, a)
        d = self.RR if c == "A" else np.sqrt(
            self.RR ** 2 + self.R ** 2 - 2 * self.R * self.RR * self.UU)
        return norm * np.exp(-a * d) * genlaguerre(n - 1, 1)(2 * a * d) / np.sqrt(4 * np.pi)

    def A_L(self, oi: Orb, oj: Orb) -> np.ndarray:
        """Multipole radial moments ``A^L(r) = 2 pi int phi_i phi_j P_L d(cos theta)``."""
        prod = self.phi(oi) * self.phi(oj)
        return np.array([2 * np.pi * _trapz(prod * self.PL[L][None, :], self.u, axis=1)
                         for L in range(self.Lmax + 1)])

    def overlap(self, oi: Orb, oj: Orb) -> float:
        return float(_trapz(self.A_L(oi, oj)[0] * self.r ** 2, self.r))

    def coulomb_center(self, oi: Orb, oj: Orb, C: str,
                       AL: Optional[np.ndarray] = None) -> float:
        """``<i| 1/r_C |j>`` for point center ``C in {'A','B'}`` (positive)."""
        if AL is None:
            AL = self.A_L(oi, oj)
        r, R = self.r, self.R
        if C == "A":
            return float(_trapz(AL[0] * r, r))                       # 1/rA spherical about A
        rlt = np.minimum(r, R)
        rgt = np.maximum(r, R)
        v = 0.0
        for L in range(self.Lmax + 1):
            if np.max(np.abs(AL[L])) < 1e-15:
                continue
            v += _trapz(AL[L] * (rlt ** L / rgt ** (L + 1)) * r ** 2, r)
        return float(v)                                             # 1/rB Legendre-expanded about A

    def nuclear(self, oi: Orb, oj: Orb) -> float:
        """``<i|(-1/r_A - 1/r_B)|j>``."""
        AL = self.A_L(oi, oj)
        return -(self.coulomb_center(oi, oj, "A", AL) + self.coulomb_center(oi, oj, "B", AL))

    def kinetic(self, oi: Orb, oj: Orb, kscale: float) -> float:
        """``<i| -1/2 grad^2 |j>`` via the shared-scale Sturmian ODE.

        ``-1/2 grad^2 chi_n = (n k / r_c) chi_n - 1/2 k^2 chi_n`` (``c`` = center of
        ``chi_n``), symmetrized over bra/ket.
        """
        AL = self.A_L(oi, oj)
        Sij = _trapz(AL[0] * self.r ** 2, self.r)
        Tj = oj[0] * kscale * self.coulomb_center(oi, oj, oj[1], AL) - 0.5 * kscale ** 2 * Sij
        Ti = oi[0] * kscale * self.coulomb_center(oi, oj, oi[1], AL) - 0.5 * kscale ** 2 * Sij
        return float(0.5 * (Ti + Tj))

    def eri(self, oi: Orb, oj: Orb, ok: Orb, ol: Orb) -> float:
        """``(ij|kl)`` chemist ``= <rho_ij | 1/r12 | rho_kl>``, ``rho_ij = phi_i phi_j``."""
        Aij = self.A_L(oi, oj)
        Akl = self.A_L(ok, ol)
        r, dr = self.r, self.dr
        total = 0.0
        for L in range(self.Lmax + 1):
            a, b = Aij[L], Akl[L]
            if np.max(np.abs(a)) < 1e-15 or np.max(np.abs(b)) < 1e-15:
                continue
            g = b * r * r
            inner = np.concatenate(([0.0], cumulative_trapezoid(g * r ** L, dx=dr))) * r ** (-(L + 1))
            outer = np.concatenate(
                ([0.0], cumulative_trapezoid((g * r ** (-(L + 1)))[::-1], dx=dr)))[::-1] * r ** L
            total += _trapz(a * (inner + outer) * r * r, r)
        return float(total)


def validate(R: float = 1.5, a: float = 1.3,
             Lmax: int = 28, nr: int = 5000, nth: int = 220, rmax: float = 60.0
             ) -> Dict[str, Tuple[float, float]]:
    """Validate the engine against known closed forms; return ``{name: (numeric, exact)}``.

    Checks: one-center ``<1s|1/r_A|1s> = a``; two-center ``<1s_A|1/r_B|1s_A> =
    (1/R)(1-(1+aR)e^{-2aR})``; one-center ``(1s1s|1s1s) = 5a/8``; two-center ``(AA|BB)``
    against :func:`geovac.two_center_eri.aabb_value`; and the mixed-scale one-center
    overlap ``(2 sqrt(ab)/(a+b))^3``.
    """
    from fractions import Fraction

    g = GoscinskianIntegrals(R=R, Lmax=Lmax, nr=nr, nth=nth, rmax=rmax)
    gfar = GoscinskianIntegrals(R=80.0, Lmax=20, nr=nr, nth=200, rmax=rmax)
    A1: Orb = (1, "A", a)
    B1: Orb = (1, "B", a)
    out: Dict[str, Tuple[float, float]] = {}

    # one-center <1s|1/rA|1s> = a  (isolate 1/rA by pushing B far away)
    out["self_A"] = (-gfar.nuclear(A1, A1) - 1.0 / 80.0, a)
    # two-center <1s_A|1/rB|1s_A>
    v_AB = (-g.nuclear(A1, A1)) - a
    out["cross_B"] = (v_AB, (1.0 / R) * (1 - (1 + a * R) * np.exp(-2 * a * R)))
    # one-center ERI 5a/8
    out["eri_1c"] = (g.eri(A1, A1, A1, A1), 5 * a / 8)
    # mixed-scale one-center overlap
    b = 0.9
    out["overlap_mixed"] = (g.overlap((1, "A", a), (1, "A", b)), (2 * np.sqrt(a * b) / (a + b)) ** 3)
    # two-center (AA|BB) vs exact closed form
    try:
        from geovac.two_center_eri import aabb_value
        af = Fraction(a).limit_denominator(1000)
        v_ex = float(aabb_value(af, (1, 0, 0), (1, 0, 0), af, (1, 0, 0), (1, 0, 0), R, prec=25))
        out["eri_aabb"] = (g.eri(A1, A1, B1, B1), v_ex)
    except Exception:  # pragma: no cover - closed-form engine optional at validate time
        pass
    return out


if __name__ == "__main__":
    for name, (num, ex) in validate().items():
        print(f"  {name:<16} numeric={num:.6f}  exact={ex:.6f}  err={abs(num - ex):.2e}")
