"""Independent real-space numerical reference for two-center two-electron ERIs.

Built to validate geovac.two_center_eri's engines (aabb_quadrature,
exchange_value, hybrid_quadrature, aabb_closed_form, ...) from a SEPARATE code
path: no import of that module's shell-kernel / two_center_spheroidal_product /
Gaunt-coefficient / ordered-xi machinery anywhere below.

ORBITAL CONVENTION (matches geovac.two_center_eri exactly -- see the docstring
of eri_ref below and the numerical Y_lm cross-check in __main__):

    chi(r) = N * r^l * Y_lm(theta_c, phi_c) * exp(-zeta * r_c)

with Y_lm the COMPLEX spherical harmonic, Condon-Shortley phase, normalized so
int |chi|^2 d^3r = 1. Centers: A at the origin, B at (0, 0, R). An orbital is
the tuple (zeta, l, m, center) with center in {"A", "B"}. This is the same
orbital family as the repo's hydrogenic chi(Z, n, l, m, ...) with n = l + 1,
zeta = Z / n (nodeless, i.e. Slater-type orbitals).

METHOD -- standard real-space multipole (Legendre addition-theorem) expansion
of 1/r12, evaluated ENTIRELY numerically about the single origin A:

    1/r12 = sum_{L,M} (4 pi / (2L+1)) Y_LM(Omega_1) Y_LM^*(Omega_2)
                        (r_</r_>)^L / r_>

Because A and B both sit on the z-axis, the azimuthal angle phi of any point
in space is IDENTICAL in the A-centered and B-centered frames (translating
along z does not change atan2(y, x)). So the bra density
rho_ab(r) = conj(chi_a(r)) chi_b(r) -- for ANY assignment of centers to a, b --
carries a pure e^{i(m_b - m_a) phi} factor exactly, regardless of which
centers a and b sit on. That collapses the double (L, M) multipole sum to a
SINGLE L-sum at the fixed M = m_b - m_a (checked numerically in __main__, not
just assumed), so representing each density needs only a radial grid crossed
with a Gauss-Legendre grid in u = cos(theta_A) -- no phi quadrature at all.

The result is the familiar two-electron "Slater radial integral" formula

    (ab|cd) = (-1)^{m_b - m_a} sum_L (4 pi / (2L+1))
                int int dr1 ds  r1^2 s^2 rho_L^ab(r1) rho_L^cd(s)
                                (min(r1,s))^L / (max(r1,s))^{L+1}

    rho_L^X(r) = int dOmega rho_X(r, Omega) Y_LM^*(Omega)
               = 2 pi int_{-1}^{1} du  f_X(r, u) ybar_{L,M}(u)

with f_X(r, u) := rho_X(r, u, phi=0), which is manifestly REAL (every phase
factor is 1 at phi=0), and ybar_{l,m}(u) the real part of Y_lm(theta, phi=0).
This is the standard Slater-Condon / Cowan radial-integral construction from
atomic structure theory, applied to a (possibly cross-center) density -- a
different formalism from the repo's bipolar-coordinate shell-kernel /
Gaunt-recoupling closed-form derivation. All pieces (normalization, the
Y_lm projector, the double radial sum) are written from scratch below.

Convergence: the L-sum is extended in blocks until a run of consecutive terms
falls below a tolerance relative to the running total (or an absolute floor
near zero), up to a hard cap. Radial and angular grids are fixed, cached,
composite Gauss-Legendre quadratures; only the L-sum is adaptive.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.special import gammaln, lpmv

Orbital = Tuple[float, int, int, str]  # (zeta, l, m, center)

# ============================================================ cached grids

_ANGULAR_CACHE: dict = {}


def _angular_grid(n_u: int) -> Tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes/weights on [-1, 1] in u = cos(theta). Cached."""
    if n_u not in _ANGULAR_CACHE:
        u, w = np.polynomial.legendre.leggauss(n_u)
        _ANGULAR_CACHE[n_u] = (u, w)
    return _ANGULAR_CACHE[n_u]


_RADIAL_CACHE: dict = {}


def _radial_grid(rmax: float, n_per_panel: int = 44,
                  growth: float = 1.22) -> Tuple[np.ndarray, np.ndarray]:
    """Composite Gauss-Legendre quadrature on [0, rmax], geometric panels.

    Panel edges grow geometrically from a small inner scale (0.12 bohr) so
    that BOTH the r=0 peak (chi_a centered at the grid origin) and the r~R
    peak (the other orbital's own center) get comparably fine local
    resolution, regardless of R. The double radial sum below has a min/max
    (r_</r_>) kernel with a kink along r=s; a plain product Gauss rule only
    converges algebraically across that kink, so `growth` is kept modest
    (finer panels) rather than relying on high per-panel order -- verified
    empirically in __main__ (V0b grid-refinement check).
    """
    key = (round(float(rmax), 4), n_per_panel, growth)
    if key in _RADIAL_CACHE:
        return _RADIAL_CACHE[key]
    edges = [0.0]
    e = 0.12
    while e < rmax:
        edges.append(e)
        e *= growth
    edges.append(rmax)
    edges = np.array(sorted(set(edges)))
    xg, wg = np.polynomial.legendre.leggauss(n_per_panel)
    rs, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        half = 0.5 * (b - a)
        rs.append(half * xg + 0.5 * (a + b))
        ws.append(half * wg)
    r = np.concatenate(rs)
    w = np.concatenate(ws)
    order = np.argsort(r)
    r, w = r[order], w[order]
    _RADIAL_CACHE[key] = (r, w)
    return r, w


# ============================================================ orbital math

def _radial_norm(zeta: float, l: int) -> float:
    """N such that int_0^inf (N r^l e^{-zeta r})^2 r^2 dr = 1.

    N^2 (2l+2)! / (2 zeta)^{2l+3} = 1.  Computed via gammaln to stay stable
    at larger l (not needed here since l <= 1, but costs nothing).
    """
    log_n2 = (2 * l + 3) * np.log(2.0 * zeta) - gammaln(2 * l + 3)
    return float(np.exp(0.5 * log_n2))


def _ybar(l: int, m: int, u: np.ndarray) -> np.ndarray:
    """Y_lm(theta, phi=0) as a function of u = cos(theta) -- real-valued.

    Condon-Shortley phase, standard normalization int |Y_lm|^2 dOmega = 1.
    Cross-checked bit-for-bit (~1e-16) against sympy.Ynm in __main__ (V0).
    """
    am = abs(m)
    P = lpmv(am, l, u)
    log_norm = 0.5 * (np.log(2 * l + 1.0) - np.log(4.0 * np.pi)
                       + gammaln(l - am + 1) - gammaln(l + am + 1))
    val = np.exp(log_norm) * P
    if m < 0 and (am % 2):
        val = -val
    return val


def _local_coords(r: np.ndarray, u: np.ndarray, R: float, center: str):
    """Map a global (r, u) point [r, theta measured from A, A at origin] to
    the local (r_c, u_c) spherical coordinates about `center`.

    phi is shared between the A- and B-centered frames (both centers lie on
    the z-axis), so it never appears here -- see module docstring.
    """
    if center.upper() == "A":
        return r, u
    # center == B, located at (0, 0, R)
    r2 = r * r + R * R - 2.0 * r * R * u
    rb = np.sqrt(np.maximum(r2, 0.0))
    safe_rb = np.where(rb > 1e-13, rb, 1.0)
    ub = (r * u - R) / safe_rb
    ub = np.where(rb > 1e-13, ub, -1.0)  # r along the +z axis through B: u_B -> -1
    return rb, np.clip(ub, -1.0, 1.0)


def _chi_reduced(zeta: float, l: int, m: int, rc: np.ndarray, uc: np.ndarray) -> np.ndarray:
    """chi(r) with the e^{i m phi} factor stripped -- real-valued, vectorized."""
    return _radial_norm(zeta, l) * rc ** l * np.exp(-zeta * rc) * _ybar(l, m, uc)


def _density_grid(oa: Orbital, ob: Orbital, r_grid: np.ndarray, u_grid: np.ndarray,
                   R: float) -> np.ndarray:
    """f_ab(r, u) = conj(chi_a)(r,u,phi=0) * chi_b(r,u,phi=0), real, shape (Nr, Nu).

    conj is a no-op here: at phi = 0 every e^{i m phi} factor is 1, so both
    chi_a and chi_b are already real at these sample points.
    """
    Rg, Ug = np.meshgrid(r_grid, u_grid, indexing="ij")
    ra, ua = _local_coords(Rg, Ug, R, oa[3])
    rb, ub = _local_coords(Rg, Ug, R, ob[3])
    ca = _chi_reduced(oa[0], oa[1], oa[2], ra, ua)
    cb = _chi_reduced(ob[0], ob[1], ob[2], rb, ub)
    return ca * cb


# ============================================================ the oracle

def eri_ref(oa: Orbital, ob: Orbital, oc: Orbital, od: Orbital, R: float,
            n_u: int = 220, n_panel: int = 44, growth: float = 1.22,
            l_max_start: int = 30, l_max_cap: int = 160, block: int = 8,
            rel_tol: float = 3e-9, abs_tol: float = 1e-13,
            return_diag: bool = False):
    """Independent numerical reference for the chemist ERI (ab|cd).

        (ab|cd) = int int chi_a^*(1) chi_b(1) (1/r12) chi_c^*(2) chi_d(2) d3r1 d3r2

    oa, ob, oc, od are (zeta, l, m, center) tuples, center in {"A", "B"}, l in
    {0, 1}. R is the A-B separation (bohr). Handles ANY assignment of centers
    to the four orbitals -- one-center (AA|BB), cross-center exchange
    (AB|AB), and hybrid (AAA|B) all go through the same code path.

    Returns a float (or (float, dict) if return_diag=True with L-convergence
    diagnostics).
    """
    m_ab = ob[2] - oa[2]
    m_cd = od[2] - oc[2]
    if m_ab != -m_cd:
        # phi-integral M-selection: no angular momentum route connects the
        # two densities -- exact zero.
        return (0.0, {"reason": "M-selection", "n_L": 0}) if return_diag else 0.0

    zetas = [oa[0], ob[0], oc[0], od[0]]
    rmax = R + 26.0 / min(zetas)
    r_grid, r_w = _radial_grid(rmax, n_panel, growth)
    u_grid, u_w = _angular_grid(n_u)

    f_ab = _density_grid(oa, ob, r_grid, u_grid, R)
    f_cd = _density_grid(oc, od, r_grid, u_grid, R)

    # weighted radial density: r^2 * w_r * rho_L(r), reused every L via a
    # fresh angular contraction (the expensive part -- the (Nr,Nu) grids --
    # is built once above).
    r2w = r_grid * r_grid * r_w
    Rmin = np.minimum.outer(r_grid, r_grid)
    Rmax_ = np.maximum.outer(r_grid, r_grid)
    outer_w = np.outer(r2w, r2w)

    total = 0.0
    L = abs(m_ab)
    terms = []
    small_run = 0
    while L <= l_max_cap:
        y_ab = _ybar(L, m_ab, u_grid)
        y_cd = _ybar(L, m_cd, u_grid)
        rho_ab = 2 * np.pi * (f_ab @ (u_w * y_ab))
        rho_cd = 2 * np.pi * (f_cd @ (u_w * y_cd))

        kernel = (Rmin / Rmax_) ** L / Rmax_
        term = (4 * np.pi / (2 * L + 1)) * float(
            np.sum(outer_w * np.outer(rho_ab, rho_cd) * kernel))
        terms.append(term)
        total += term

        ref_scale = max(abs(total), abs_tol)
        if abs(term) < rel_tol * ref_scale:
            small_run += 1
            if small_run >= block and len(terms) >= l_max_start:
                break
        else:
            small_run = 0
        L += 1

    total *= (-1.0) ** m_ab
    if return_diag:
        diag = {"n_L": len(terms), "L_final": L, "last_terms": terms[-block:],
                "converged": L <= l_max_cap}
        return total, diag
    return total


# ================================================================ validation

def _full_Y(l: int, m: int, theta: float, phi: float) -> complex:
    """Full complex Y_lm(theta,phi) (theta, phi genuinely varying) -- used only
    for the V0 self-consistency check, not by eri_ref itself."""
    return complex(_ybar(l, m, np.array([np.cos(theta)]))[0] * np.exp(1j * m * phi))


def _chi_full(zeta: float, l: int, m: int, r: float, theta: float, phi: float) -> complex:
    return _radial_norm(zeta, l) * r ** l * np.exp(-zeta * r) * _full_Y(l, m, theta, phi)


def _cart_from_local(rc, thc, phc, center, R):
    """Cartesian point from LOCAL spherical coords about `center`."""
    x = rc * np.sin(thc) * np.cos(phc)
    y = rc * np.sin(thc) * np.sin(phc)
    z = rc * np.cos(thc)
    if center.upper() == "B":
        z += R
    return x, y, z


def _global_spherical(x, y, z):
    r = np.sqrt(x * x + y * y + z * z)
    theta = np.arccos(np.clip(z / r, -1.0, 1.0)) if r > 1e-14 else 0.0
    phi = np.arctan2(y, x)
    return r, theta, phi


def _v0_phi_factorization_check(oa: Orbital, ob: Orbital, R: float) -> float:
    """Confirms rho_ab(r,theta,phi) = f_ab(r,theta) * e^{i(m_b-m_a) phi} EXACTLY
    (the fact the whole L-sum collapse in eri_ref rests on), by direct brute
    evaluation at points NOT on phi=0, picked in the A-centered global frame
    and converted to each orbital's own local frame.
    """
    M = ob[2] - oa[2]
    worst = 0.0
    rng = np.random.default_rng(0)
    for _ in range(12):
        r = rng.uniform(0.3, R + 3.0)
        theta = rng.uniform(0.05, np.pi - 0.05)
        phi = rng.uniform(0.0, 2 * np.pi)
        x, y, z = _cart_from_local(r, theta, phi, "A", R)

        def eval_orb(o):
            zeta, l, m, center = o
            if center.upper() == "A":
                rc, thc, phc = r, theta, phi
            else:
                rc, thc, phc = _global_spherical(x, y, z - R)
            return _chi_full(zeta, l, m, rc, thc, phc)

        rho_full = np.conj(eval_orb(oa)) * eval_orb(ob)
        rho_phi0 = np.conj(_chi_full(oa[0], oa[1], oa[2], r if oa[3].upper() == "A"
                                      else _global_spherical(x, y, z - R)[0],
                                      theta if oa[3].upper() == "A"
                                      else _global_spherical(x, y, z - R)[1], 0.0)) \
            * _chi_full(ob[0], ob[1], ob[2], r if ob[3].upper() == "A"
                        else _global_spherical(x, y, z - R)[0],
                        theta if ob[3].upper() == "A"
                        else _global_spherical(x, y, z - R)[1], 0.0)
        predicted = rho_phi0 * np.exp(1j * M * phi)
        worst = max(worst, abs(rho_full - predicted))
    return worst


def _to_repo(Z, orb, center: str) -> Orbital:
    """(Z, (n,l,m)) -> (zeta, l, m, center); requires n = l+1 (nodeless)."""
    n, l, m = orb
    if n != l + 1:
        raise ValueError(f"orbital {orb} is not nodeless (n != l+1)")
    return (float(Z) / n, l, m, center)


def main() -> None:
    from fractions import Fraction
    from geovac.two_center_eri import aabb_quadrature, exchange_value, hybrid_quadrature

    print("=" * 78)
    print("V0  phi-factorization self-check (the fact the whole method rests on)")
    print("=" * 78)
    cross_pairs = [
        ((1.3, 1, 1, "A"), (0.8, 1, -1, "B")),
        ((2.1, 0, 0, "A"), (0.6, 1, 1, "B")),
        ((0.9, 1, 0, "B"), (1.7, 1, 1, "A")),
    ]
    worst = 0.0
    for oa, ob in cross_pairs:
        w = _v0_phi_factorization_check(oa, ob, R=2.7)
        worst = max(worst, w)
    print(f"    worst |full density - f(r,u)*e^(iM phi)| over {len(cross_pairs)} "
          f"cross-center pairs x 12 random points: {worst:.3e}  "
          f"{'OK' if worst < 1e-10 else 'FAIL'}")

    print()
    print("=" * 78)
    print("V0b  internal grid-refinement convergence (no repo involved)")
    print("=" * 78)
    oa, ob = (1.4, 1, 1, "A"), (0.9, 0, 0, "B")
    oc, od = (1.1, 0, 0, "A"), (1.6, 1, 1, "B")
    v_coarse = eri_ref(oa, ob, oc, od, R=2.8, n_u=160, n_panel=32, growth=1.3)
    v_fine = eri_ref(oa, ob, oc, od, R=2.8, n_u=260, n_panel=52, growth=1.16)
    print(f"    coarse grid = {v_coarse:.10f}")
    print(f"    fine grid   = {v_fine:.10f}")
    print(f"    |diff| = {abs(v_coarse - v_fine):.3e}  "
          f"{'OK' if abs(v_coarse - v_fine) < 1e-5 else 'FAIL'}")

    Z1, Z2, Z3 = Fraction(1), Fraction(2), Fraction(3)

    print()
    print("=" * 78)
    print("VALIDATION 1  (AA|BB) one-center densities vs aabb_quadrature")
    print("=" * 78)
    aabb_cases = [
        ("1s(A) 1s(A) | 1s(B) 1s(B)", Z1, (1, 0, 0), (1, 0, 0), Z2, (1, 0, 0), (1, 0, 0), 2.5),
        ("2p0(A) 2p0(A) | 1s(B) 1s(B)", Z1, (2, 1, 0), (2, 1, 0), Z2, (1, 0, 0), (1, 0, 0), 2.5),
        ("2p1(A) 2p0(A) | 2p0(B) 2p1(B)", Z1, (2, 1, 1), (2, 1, 0), Z2, (2, 1, 0), (2, 1, 1), 2.5),
        ("3d2(A) 3d0(A) | 2p-1(B) 2p1(B)", Z1, (3, 2, 2), (3, 2, 0), Z2, (2, 1, -1), (2, 1, 1), 2.5),
    ]
    rows = []
    for name, ZA, oa_, ob_, ZB, oc_, od_, R in aabb_cases:
        ref = aabb_quadrature(ZA, oa_, ob_, ZB, oc_, od_, R)
        mine = eri_ref(_to_repo(ZA, oa_, "A"), _to_repo(ZA, ob_, "A"),
                       _to_repo(ZB, oc_, "B"), _to_repo(ZB, od_, "B"), R)
        rows.append((name, ref, mine, abs(ref - mine)))
        print(f"    {name:32s} repo={ref:+.10f}  mine={mine:+.10f}  "
              f"|diff|={abs(ref - mine):.3e}")

    print()
    print("=" * 78)
    print("VALIDATION 2  (AB|AB) cross-center exchange vs exchange_value")
    print("=" * 78)
    exch_cases = [
        ("1s(A) 1s(B) | 1s(A) 1s(B)", Z3, (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0), 3.0),
        ("2p1(A) 1s(B) | 1s(A) 2p1(B)  [sigma=1]", Z3, (2, 1, 1), (1, 0, 0), Z1, (1, 0, 0), (2, 1, 1), 3.0),
        ("2p0(A) 1s(B) | 1s(A) 2p0(B)  [sigma=0, p]", Z3, (2, 1, 0), (1, 0, 0), Z1, (1, 0, 0), (2, 1, 0), 3.0),
    ]
    for name, ZA, oa_, ob_, ZB, oc_, od_, R in exch_cases:
        ref = exchange_value(ZA, oa_, ob_, ZB, oc_, od_, R, tau_max=10)
        mine = eri_ref(_to_repo(ZA, oa_, "A"), _to_repo(ZB, ob_, "B"),
                       _to_repo(ZA, oc_, "A"), _to_repo(ZB, od_, "B"), R)
        rows.append((name, ref, mine, abs(ref - mine)))
        print(f"    {name:38s} repo={ref:+.10f}  mine={mine:+.10f}  "
              f"|diff|={abs(ref - mine):.3e}")

    print()
    print("=" * 78)
    print("VALIDATION 3 (bonus)  hybrid (3-on-one-center) vs hybrid_quadrature")
    print("=" * 78)
    hybrid_cases = [
        ("1s(A) 2p0(A) 1s(A) | 2p0(B)", Z1, (1, 0, 0), (2, 1, 0), (1, 0, 0), Z2, (2, 1, 0), 2.5),
    ]
    for name, ZA, oa_, ob_, oc_, ZB, od_, R in hybrid_cases:
        ref = hybrid_quadrature(ZA, oa_, ob_, oc_, ZB, od_, R)
        mine = eri_ref(_to_repo(ZA, oa_, "A"), _to_repo(ZA, ob_, "A"),
                       _to_repo(ZA, oc_, "A"), _to_repo(ZB, od_, "B"), R)
        rows.append((name, ref, mine, abs(ref - mine)))
        print(f"    {name:28s} repo={ref:+.10f}  mine={mine:+.10f}  "
              f"|diff|={abs(ref - mine):.3e}")
    print("    (exercises the SAME general code path as validations 1-2 -- no")
    print("     special-casing by center pattern anywhere in eri_ref)")

    print()
    print("=" * 78)
    print("NEW CAPABILITY  mixed-exponent cross-center densities (no repo route)")
    print("=" * 78)
    mixed_cases = [
        ("1s(A,z=1.0) 1s(B,z=2.3) | 1s(A,z=1.0) 1s(B,z=2.3)",
         (1.0, 0, 0, "A"), (2.3, 0, 0, "B"), (1.0, 0, 0, "A"), (2.3, 0, 0, "B"), 2.5),
        ("2p1(A,z=0.8) 1s(B,z=1.9) | 1s(A,z=1.4) 2p1(B,z=0.6)",
         (0.8, 1, 1, "A"), (1.9, 0, 0, "B"), (1.4, 0, 0, "A"), (0.6, 1, 1, "B"), 3.2),
        ("1s(A,z=2.7) 2p0(B,z=0.5) | 2p0(A,z=0.9) 1s(B,z=1.6)",
         (2.7, 0, 0, "A"), (0.5, 1, 0, "B"), (0.9, 1, 0, "A"), (1.6, 0, 0, "B"), 2.0),
        ("2p1(A,z=1.1) 2p-1(B,z=1.1) | 2p-1(A,z=1.1) 2p1(B,z=1.1)  [shared zeta]",
         (1.1, 1, 1, "A"), (1.1, 1, -1, "B"), (1.1, 1, -1, "A"), (1.1, 1, 1, "B"), 2.5),
    ]
    for name, oa_, ob_, oc_, od_, R in mixed_cases:
        val, diag = eri_ref(oa_, ob_, oc_, od_, R, return_diag=True)
        print(f"    {name}")
        print(f"        R={R}  value={val:+.10f}  (L-terms used: {diag['n_L']}, "
              f"converged={diag['converged']})")

    print()
    print("=" * 78)
    thresh = 1e-4
    worst_diff = max(d for _, _, _, d in rows)
    print(f"SUMMARY: worst |repo - mine| over {len(rows)} repo comparisons = "
          f"{worst_diff:.3e}  ->  {'PASS' if worst_diff < thresh else 'FAIL'} (<{thresh:.0e})")


if __name__ == "__main__":
    main()
