r"""Route C / C2 (the ONE hard new primitive): analytic, R-accurate two-center
prolate core-valence ERIs for all-electron LiH.

The C2 diagnostic (debug/track_logs/prolate_native_lih.md; debug/eri_core_grid_diagnostic.py)
showed grid ERIs are R-INACCURATE for tight-core densities: (cc|cc) swings 26 mHa
across the bond range and sits +132 mHa off exact, (cc|vv)/(cv|cv) ~4 mHa / +21 mHa.
That kills the all-electron increment-7 route (0.36 Ha core wall). C1 already showed
the ONE-BODY core integrals go analytic and R-flat to ~1e-50. This module is C2:
the TWO-BODY (ERI) analog.

WHAT IS NEW (Task 1).  The analytic quadrature-free V_ee (geovac/prolate_recondition.py
vee_mp / _build_Xtab_mp; float sibling geovac/neumann_vee_general_m.py build_Xtab)
is a Neumann partial-wave sum over l of X_l(xi-side) . Y1(eta) . Y2(eta), with the
xi-side X-table

    X_l^{m,s}(P1,P2; c) = int int xi1^P1 xi2^P2 (xi^2-1)^s e^{-c(xi1+xi2)}
                                   P_l^m(xi<) Q_l^m(xi>) dxi1 dxi2.

Every existing routine HARD-WIRES a SINGLE per-electron rate c = 2*alpha (line 334
of prolate_recondition, and the P1<->P2 symmetrization).  Core+valence ERIs need
DISTINCT per-electron rates c1 = alpha_p + alpha_q, c2 = alpha_r + alpha_s (the two
electrons carry different combined orbital exponents).  ``build_Xtab_pair`` below
generalizes _build_Xtab_mp to (c1, c2): the two orderings xi1<xi2 and xi2<xi1 are
now DISTINCT (no P1<->P2 symmetry when c1 != c2), each assembled from A-moments /
B-table at its OWN rate with the IBP-tail correction B-table read at ccorr = c1+c2.
At c1 = c2 it reduces to _build_Xtab_mp bit-for-bit (the frozen falsifier, G1).

WHAT IS ASSEMBLED (Task 2).  ``eri_sigma`` is the SIGMA-ONLY (all mu=0 => m=0,s=0)
ERI assembler for two-center prolate orbitals that are either
  - VALENCE: an ordinary prolate ProductFn  xi^j eta^l e^{-alpha_v xi}, or
  - CORE (Route A): a center-A 1s STO e^{-zeta r_A} = e^{-alpha xi} . e^{-alpha eta},
    alpha = zeta*R/2, with the eta-exponential expanded in Legendre polynomials
    e^{-alpha eta} = sum_l (2l+1)(-1)^l i_l(alpha) P_l(eta)  (i_l = modified spherical
    Bessel first kind), i.e. represented as e^{-alpha xi} . (eta-polynomial).
The two electron densities (phi_p.phi_q and phi_r.phi_s) become prolate product
densities with per-electron xi-rate c1/c2 and eta-polynomial = product of the two
orbital eta-polynomials; the ERI is the Neumann contraction with the (xi^2-eta^2)
Jacobian, mirroring vee_mp's assembly kernel.  (cc|cc) is special-cased to 5*zeta/8.

Run from root:  python debug/prolate_mixed_eri.py
"""
from __future__ import annotations

import os
import sys
import time
from functools import lru_cache
from typing import Dict, List, Sequence, Tuple

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geovac import neumann_vee_general_m as ngm    # noqa: E402
from geovac import prolate_recondition as pr        # noqa: E402

# High working precision: the m=0 Neumann prefactor is (2l+1) (NO factorial
# suppression, unlike m>0), so high-l B-recurrence errors are not damped in the
# energy -- give the forward Q_l recurrence headroom.
mp.mp.dps = 60

ZC_LI: float = 3.0 - 5.0 / 16.0          # 2.6875 (He-like Li^2+ 1s exponent)

_pm = ngm._polymul
_pa = ngm._polyadd
_shift = ngm._shift


# ============================================================
# Task 1: build_Xtab_pair -- two per-electron rates (c1, c2)
# ============================================================
def build_Xtab_pair(ms_pairs: List[Tuple[int, int]], l_neumann: int, p_max: int,
                    alpha1: float, alpha2: float,
                    l_caps: Dict[Tuple[int, int], int]
                    ) -> Dict[Tuple[int, int, int], List[List[mp.mpf]]]:
    r"""Two-rate X table:  X_l^{m,s}(P1,P2; c1,c2), c1 = 2*alpha1, c2 = 2*alpha2.

    Generalizes prolate_recondition._build_Xtab_mp (single rate c = 2*alpha) by
    splitting the double integral into its two orderings, which are NO LONGER
    related by P1<->P2 when c1 != c2:

        X_l(P1,P2; c1,c2) = J1 + J2,
        J1 = A_l(P1; c1) . B_l(P2; c2) - corr(W[P1], P2; inner rate c1, B at c1+c2),
        J2 = A_l(P2; c2) . B_l(P1; c1) - corr(W[P2], P1; inner rate c2, B at c1+c2),

    where J1 is the ordering xi1<xi2 (electron 1 = P_l inner at rate c1, electron 2
    = Q_l outer at rate c2) and J2 the reverse.  A_l = <xi^P (xi^2-1)^s d^mP_l>_c,
    B_l = <xi^P (xi^2-1)^s d^mQ_l>_c (ngm._A_moment / ngm._B_table), W = _W_poly
    (rate-independent).  The IBP-tail correction integrates
    e^{-c_inner xi_o} . [tail of the inner integral] against the outer Q_l e^{-c_o xi_o},
    so its B-table is at the SUM rate c1+c2; the 1/c_inner powers use the inner rate.

    Returns dict (l,m,s) -> (p_max+1)x(p_max+1) NON-SYMMETRIC list-of-lists of mpf.
    At alpha1 == alpha2 it reproduces _build_Xtab_mp bit-for-bit.
    """
    c1 = mp.mpf(2.0 * alpha1)
    c2 = mp.mpf(2.0 * alpha2)
    csum = c1 + c2
    Xtab: Dict[Tuple[int, int, int], List[List[mp.mpf]]] = {}
    for (m, s) in ms_pairs:
        l_hi = min(l_neumann, l_caps[(m, s)])
        if l_hi < m:
            continue
        deg_extra = p_max + 2 * s + (l_hi - m)
        p_corr_max = p_max + deg_extra
        n_mono = p_max + 2 * s + (l_hi - m) + 2
        Amono1 = ngm._mono_moments(c1, n_mono)
        Amono2 = ngm._mono_moments(c2, n_mono)
        Bc1 = ngm._B_table(m, s, l_hi, p_max, c1)
        Bc2 = ngm._B_table(m, s, l_hi, p_max, c2)
        Bsum = ngm._B_table(m, s, l_hi, p_corr_max, csum)
        for l in range(m, l_hi + 1):
            mat = [[mp.mpf(0)] * (p_max + 1) for _ in range(p_max + 1)]
            Av1 = [ngm._A_moment(l, m, s, P, Amono1) for P in range(p_max + 1)]
            Av2 = [ngm._A_moment(l, m, s, P, Amono2) for P in range(p_max + 1)]
            Wf = {P: [mp.mpf(w) for w in ngm._W_poly(l, m, s, P)]
                  for P in range(p_max + 1)}
            for P1 in range(p_max + 1):
                for P2 in range(p_max + 1):
                    J1 = Av1[P1] * Bc2[(l, P2)] - pr._corr_mp(Wf[P1], P2, l, c1, Bsum)
                    J2 = Av2[P2] * Bc1[(l, P1)] - pr._corr_mp(Wf[P2], P1, l, c2, Bsum)
                    mat[P1][P2] = J1 + J2
            Xtab[(l, m, s)] = mat
    return Xtab


# ============================================================
# Orbital representation (sigma-only, two-center prolate)
# ============================================================
def _i_sph(l: int, x) -> mp.mpf:
    """Modified spherical Bessel first kind i_l(x) = sqrt(pi/2x) I_{l+1/2}(x)."""
    x = mp.mpf(x)
    return mp.sqrt(mp.pi / (2 * x)) * mp.besseli(l + mp.mpf(1) / 2, x)


def sto_eta_poly(alpha, L: int = 24) -> List[mp.mpf]:
    """Monomial coeffs (low->high) of the Legendre truncation of e^{-alpha*eta}:
    e^{-a eta} = sum_{l<L} (2l+1)(-1)^l i_l(a) P_l(eta)."""
    alpha = mp.mpf(alpha)
    coeffs: List[mp.mpf] = [mp.mpf(0)]
    for l in range(L):
        a_l = (2 * l + 1) * mp.mpf(-1) ** l * _i_sph(l, alpha)
        coeffs = _pa(coeffs, [a_l * c for c in ngm._leg_coeffs(l)])
    return coeffs


class Orbital:
    """One sigma prolate orbital: xi_power (int j) . eta_poly(eta) . e^{-alpha xi},
    times ``norm``.  ``zeta``/``is_core`` tag a 1s STO (for the (cc|cc) closed form)."""

    __slots__ = ("xi_power", "eta_poly", "alpha", "norm", "zeta", "is_core")

    def __init__(self, xi_power: int, eta_poly: List[mp.mpf], alpha,
                 norm, zeta=None, is_core: bool = False):
        self.xi_power = xi_power
        self.eta_poly = eta_poly
        self.alpha = mp.mpf(alpha)
        self.norm = mp.mpf(norm)
        self.zeta = None if zeta is None else mp.mpf(zeta)
        self.is_core = is_core


def sto_orbital(zeta, R, L: int = 24, is_core: bool = False) -> Orbital:
    """Unit-normalized center-A 1s STO  sqrt(zeta^3/pi) e^{-zeta r_A}, r_A=(R/2)(xi+eta),
    as a prolate ProductFn:  norm . e^{-alpha xi} . [Legendre(e^{-alpha eta})], alpha=zeta*R/2."""
    zeta = mp.mpf(zeta)
    R = mp.mpf(R)
    alpha = zeta * R / 2
    norm = mp.sqrt(zeta ** 3 / mp.pi)
    return Orbital(0, sto_eta_poly(alpha, L), alpha, norm, zeta=zeta, is_core=is_core)


def valence_prolate_orbital(j: int, l: int, alpha) -> Orbital:
    """General (unnormalized) prolate ProductFn xi^j eta^l e^{-alpha xi} (sigma)."""
    return Orbital(j, _shift([mp.mpf(1)], l), alpha, mp.mpf(1))


# ============================================================
# eta moment  int_{-1}^1 eta^Q P_l(eta) deta  (sigma, m=0)
# ============================================================
@lru_cache(maxsize=None)
def _Ymom(l: int, Q: int) -> mp.mpf:
    return pr._mom_eta(_pm(_shift([mp.mpf(1)], Q), list(ngm._leg_coeffs(l))))


def _sum_Y(epoly: Sequence[mp.mpf], l: int, dQ: int) -> mp.mpf:
    """sum_Q epoly[Q] . int eta^{Q+dQ} P_l(eta) deta."""
    s = mp.mpf(0)
    for Q, co in enumerate(epoly):
        if co == 0:
            continue
        ym = _Ymom(l, Q + dQ)
        if ym != 0:
            s += co * ym
    return s


# ============================================================
# Task 2: sigma-only mixed ERI assembler
# ============================================================
_JAC = [(1, 2, 0, 2, 0), (-1, 2, 0, 0, 2), (-1, 0, 2, 2, 0), (1, 0, 2, 0, 2)]


def _leff(beta, tol=mp.mpf("1e-22"), lmax: int = 70) -> int:
    """Smallest l past which (2l+1)|i_l(beta)| < tol (the eta Legendre content cutoff)."""
    beta = mp.mpf(beta)
    for l in range(lmax):
        if (2 * l + 1) * abs(_i_sph(l, beta)) < tol:
            return l
    return lmax


def _auto_lneu(op: Orbital, oq: Orbital, orr: Orbital, os_: Orbital,
               e1: Sequence, e2: Sequence) -> int:
    """l cutoff: min over the two electrons.  The l-sum TERMINATES at the eta
    polynomial degree (Ymom(l,Q)=0 for l>Q), so min(deg e1, deg e2)+2 is EXACT;
    for STO densities the Legendre content dies earlier (Leff), so take the min."""
    l_finite = min(len(e1) - 1, len(e2) - 1) + 2
    all_sto = all(o.zeta is not None for o in (op, oq, orr, os_))
    if all_sto:
        c1eta = op.alpha + oq.alpha          # eta-exponent of density 1
        c2eta = orr.alpha + os_.alpha
        l_adapt = min(_leff(c1eta), _leff(c2eta)) + 3
        return min(l_finite, l_adapt)
    return l_finite


def _eri_neumann_sigma(op: Orbital, oq: Orbital, orr: Orbital, os_: Orbital,
                       R, l_neumann: int) -> mp.mpf:
    """Full Neumann sigma ERI (pq|rs); no (cc|cc) shortcut (used for diagnostics)."""
    R = mp.mpf(R)
    e1 = _pm(op.eta_poly, oq.eta_poly)       # electron-1 density eta polynomial
    e2 = _pm(orr.eta_poly, os_.eta_poly)
    p1 = op.xi_power + oq.xi_power            # electron-1 density xi power
    p2 = orr.xi_power + os_.xi_power
    c1 = op.alpha + oq.alpha                  # per-electron xi rate = 2*alpha1
    c2 = orr.alpha + os_.alpha
    N = op.norm * oq.norm * orr.norm * os_.norm
    alpha1 = c1 / 2
    alpha2 = c2 / 2
    p_max = max(p1, p2) + 2
    l_caps = {(0, 0): l_neumann}
    Xtab = build_Xtab_pair([(0, 0)], l_neumann, p_max, alpha1, alpha2, l_caps)

    h6 = (R / 2) ** 6
    pref = (2 / R) * h6 * (2 * mp.pi) ** 2
    tot = mp.mpf(0)
    for l in range(0, l_neumann + 1):
        X = Xtab.get((l, 0, 0))
        if X is None:
            continue
        npre = mp.mpf(2 * l + 1)
        for (sgn, dP1, dQ1, dP2, dQ2) in _JAC:
            P1, P2 = p1 + dP1, p2 + dP2
            if P1 > p_max or P2 > p_max:
                continue
            Y1 = _sum_Y(e1, l, dQ1)
            if Y1 == 0:
                continue
            Y2 = _sum_Y(e2, l, dQ2)
            if Y2 == 0:
                continue
            tot += sgn * npre * X[P1][P2] * Y1 * Y2
    return N * pref * tot


def eri_sigma(op: Orbital, oq: Orbital, orr: Orbital, os_: Orbital,
              R, l_neumann: int = None) -> mp.mpf:
    """Sigma-only chemist ERI (pq|rs) = int phi_p phi_q(1) 1/r12 phi_r phi_s(2).

    (cc|cc) -- all four the SAME core 1s STO -- is returned as the closed form
    5*zeta/8 (the grid was +132 mHa off on this one).  Everything else goes
    through the two-rate Neumann sum.
    """
    if (op.is_core and oq.is_core and orr.is_core and os_.is_core
            and op.zeta == oq.zeta == orr.zeta == os_.zeta):
        return 5 * mp.mpf(op.zeta) / 8
    if l_neumann is None:
        e1 = _pm(op.eta_poly, oq.eta_poly)
        e2 = _pm(orr.eta_poly, os_.eta_poly)
        l_neumann = _auto_lneu(op, oq, orr, os_, e1, e2)
    return _eri_neumann_sigma(op, oq, orr, os_, R, l_neumann)


# ============================================================
# Exact single-center references (radial quadrature, unit-normalized densities)
# ============================================================
def _hartree_pot_spherical(r, amp, beta) -> mp.mpf:
    """V(r) = int amp e^{-beta r'} / |r-r'| d3r'  for a spherical density (closed form)."""
    r, amp, beta = mp.mpf(r), mp.mpf(amp), mp.mpf(beta)
    lower = (2 / beta ** 3) - mp.e ** (-beta * r) * (r ** 2 / beta + 2 * r / beta ** 2
                                                     + 2 / beta ** 3)
    upper = mp.e ** (-beta * r) * (r / beta + 1 / beta ** 2)
    return 4 * mp.pi * amp * (lower / r + upper)


def ref_J(amp_a, beta_a, amp_b, beta_b) -> mp.mpf:
    """J = int int rho_a(1) rho_b(2)/r12,  rho_x(r)=amp_x e^{-beta_x r}, both spherical."""
    amp_b, beta_b = mp.mpf(amp_b), mp.mpf(beta_b)

    def f(r):
        return (4 * mp.pi * amp_b * mp.e ** (-beta_b * r)
                * _hartree_pot_spherical(r, amp_a, beta_a) * r ** 2)

    return mp.quad(f, [0, mp.mpf("0.5"), 1, 2, 5, 10, mp.inf])


# ============================================================
# Gates
# ============================================================
def _relerr(a, b) -> float:
    a, b = mp.mpf(a), mp.mpf(b)
    if b == 0:
        return float(abs(a))
    return float(abs(a - b) / abs(b))


def gate_G1() -> bool:
    """Reduction FALSIFIER: build_Xtab_pair(alpha,alpha) == _build_Xtab_mp bit-for-bit."""
    print("G1  REDUCTION FALSIFIER  (build_Xtab_pair|_{c1=c2} vs _build_Xtab_mp)", flush=True)
    alpha = 0.8
    ms_pairs = [(0, 0), (1, 1), (2, 1), (0, 1)]     # exercise s>0, m>0 too
    l_neumann, p_max = 8, 4
    l_caps = {ms: l_neumann for ms in ms_pairs}
    Xref = pr._build_Xtab_mp(ms_pairs, l_neumann, p_max, alpha, l_caps)
    Xpair = build_Xtab_pair(ms_pairs, l_neumann, p_max, alpha, alpha, l_caps)
    worst = 0.0
    nz = 0
    assert set(Xref.keys()) == set(Xpair.keys()), "key mismatch"
    for key, mref in Xref.items():
        mpair = Xpair[key]
        for i in range(len(mref)):
            for j in range(len(mref[i])):
                a, b = mpair[i][j], mref[i][j]
                if b != 0:
                    nz += 1
                    d = _relerr(a, b)
                    worst = max(worst, d)
                elif a != 0:
                    worst = max(worst, float(abs(a)))
    ok = worst < 1e-25
    print(f"    {len(Xref)} blocks, {nz} nonzero entries; max rel diff = {worst:.2e}", flush=True)
    print(f"    G1: {'PASS (bit-for-bit)' if ok else 'FAIL'}", flush=True)
    return ok


def gate_G2_G5() -> bool:
    zc, zv = mp.mpf(ZC_LI), mp.mpf("0.8")
    Rs = [2.70, 2.85, 3.015, 3.20, 3.45]
    ok_all = True

    # exact references (R-independent, single-center, unit-normalized)
    amp_c, beta_c = zc ** 3 / mp.pi, 2 * zc          # cc density = |1s_c|^2
    amp_v, beta_v = zv ** 3 / mp.pi, 2 * zv          # vv density = |1s_v|^2
    amp_cv, beta_cv = mp.sqrt(zc ** 3 * zv ** 3) / mp.pi, zc + zv   # cross density
    ref_cccc = 5 * zc / 8
    ref_vvvv = 5 * zv / 8
    ref_ccvv = ref_J(amp_c, beta_c, amp_v, beta_v)
    ref_cvcv = ref_J(amp_cv, beta_cv, amp_cv, beta_cv)
    # cross-check the reference machine itself against the 5z/8 closed forms
    chk_cc = _relerr(ref_J(amp_c, beta_c, amp_c, beta_c), ref_cccc)
    chk_vv = _relerr(ref_J(amp_v, beta_v, amp_v, beta_v), ref_vvvv)
    print(f"\n  reference self-check (radial-quad J vs 5z/8): "
          f"cc {chk_cc:.1e}, vv {chk_vv:.1e}", flush=True)
    print(f"  refs:  (cc|cc)=5zc/8={float(ref_cccc):.9f}  (vv|vv)=5zv/8={float(ref_vvvv):.9f}",
          flush=True)
    print(f"         (cc|vv)={float(ref_ccvv):.9f}  (cv|cv)={float(ref_cvcv):.9f}", flush=True)
    # honest note: the task's V_Hcore=(2/r)[...] is a 2-electron-core convention;
    # unit-normalized orbitals (forced by G2/G3 = 5z/8) give the (1/r)[...] potential,
    # i.e. exactly ref_ccvv above.  Show the factor-2 sibling so the choice is visible.
    print(f"         [task's (2/r) 2e-core sibling of (cc|vv) would be "
          f"{float(2*ref_ccvv):.6f}]", flush=True)

    # G2 (cc|cc) special case
    print("\nG2  (cc|cc) special case = 5 zc/8", flush=True)
    R0 = 3.015
    oc = sto_orbital(zc, R0, is_core=True)
    val = eri_sigma(oc, oc, oc, oc, R0)
    g2 = _relerr(val, ref_cccc)
    print(f"    eri_sigma(cc|cc) = {float(val):.9f}   rel = {g2:.2e}   "
          f"{'PASS' if g2 < 1e-12 else 'FAIL'}", flush=True)
    # and the Neumann path on this hardest (tightest) density, as a diagnostic
    ncc = _eri_neumann_sigma(oc, oc, oc, oc, R0,
                             _auto_lneu(oc, oc, oc, oc,
                                        _pm(oc.eta_poly, oc.eta_poly),
                                        _pm(oc.eta_poly, oc.eta_poly)))
    print(f"    [(cc|cc) via full Neumann path = {float(ncc):.9f}  rel {_relerr(ncc, ref_cccc):.2e}]",
          flush=True)
    ok_all &= g2 < 1e-12

    # G3 (vv|vv) via Neumann
    print("\nG3  (vv|vv) via Neumann = 5 zv/8", flush=True)
    ov = sto_orbital(zv, R0, is_core=False)
    val = eri_sigma(ov, ov, ov, ov, R0)
    g3 = _relerr(val, ref_vvvv)
    print(f"    eri_sigma(vv|vv) = {float(val):.9f}   rel = {g3:.2e}   "
          f"{'PASS' if g3 < 1e-10 else 'FAIL'}", flush=True)
    ok_all &= g3 < 1e-10

    # G4 (cc|vv) and (cv|cv) via Neumann at R0
    print("\nG4  (cc|vv), (cv|cv) vs exact single-center reference", flush=True)
    v_ccvv = eri_sigma(oc, oc, ov, ov, R0)
    v_cvcv = eri_sigma(oc, ov, oc, ov, R0)
    g4a = _relerr(v_ccvv, ref_ccvv)
    g4b = _relerr(v_cvcv, ref_cvcv)
    print(f"    (cc|vv) = {float(v_ccvv):.9f}  ref {float(ref_ccvv):.9f}  rel {g4a:.2e}  "
          f"{'PASS' if g4a < 1e-9 else 'FAIL'}", flush=True)
    print(f"    (cv|cv) = {float(v_cvcv):.9f}  ref {float(ref_cvcv):.9f}  rel {g4b:.2e}  "
          f"{'PASS' if g4b < 1e-9 else 'FAIL'}", flush=True)
    ok_all &= (g4a < 1e-9 and g4b < 1e-9)

    # G5 R-independence: every core-touching ERI flat across R
    print("\nG5  R-INDEPENDENCE (core-touching ERIs flat across R)", flush=True)
    print(f"    {'R':>6} {'(cc|cc)':>14} {'(cc|vv)':>14} {'(cv|cv)':>14} "
          f"{'ccvv relerr':>13} {'cvcv relerr':>13}", flush=True)
    cccc_vals, ccvv_vals, cvcv_vals = [], [], []
    for R in Rs:
        oc_R = sto_orbital(zc, R, is_core=True)
        ov_R = sto_orbital(zv, R, is_core=False)
        e_cccc = eri_sigma(oc_R, oc_R, oc_R, oc_R, R)         # special-cased
        e_ccvv = eri_sigma(oc_R, oc_R, ov_R, ov_R, R)
        e_cvcv = eri_sigma(oc_R, ov_R, oc_R, ov_R, R)
        cccc_vals.append(e_cccc)
        ccvv_vals.append(e_ccvv)
        cvcv_vals.append(e_cvcv)
        print(f"    {R:6.3f} {float(e_cccc):14.9f} {float(e_ccvv):14.9f} "
              f"{float(e_cvcv):14.9f} {_relerr(e_ccvv, ref_ccvv):13.2e} "
              f"{_relerr(e_cvcv, ref_cvcv):13.2e}", flush=True)

    def spread(vals):
        fv = [float(v) for v in vals]
        return max(fv) - min(fv)

    s_cccc = spread(cccc_vals)
    s_ccvv = spread(ccvv_vals)
    s_cvcv = spread(cvcv_vals)
    print(f"    R-SPREAD (Ha):  (cc|cc) {s_cccc:.2e}   (cc|vv) {s_ccvv:.2e}   "
          f"(cv|cv) {s_cvcv:.2e}", flush=True)
    print(f"    (grid wall was: (cc|cc) 26 mHa, (cc|vv)/(cv|cv) ~4.2 mHa spread)", flush=True)
    g5 = s_cccc < 1e-3 and s_ccvv < 1e-3 and s_cvcv < 1e-3
    print(f"    G5: {'PASS (<< 1 mHa)' if g5 else 'FAIL'}", flush=True)
    ok_all &= g5
    return ok_all


if __name__ == "__main__":
    t0 = time.time()
    print(f"Route C / C2: analytic two-center prolate core-valence ERIs "
          f"(mp.dps={mp.mp.dps}, zc={ZC_LI})\n")
    r1 = gate_G1()
    r2 = gate_G2_G5()
    print(f"\n=== ALL GATES: {'PASS' if (r1 and r2) else 'FAIL'} "
          f"({time.time() - t0:.0f}s) ===", flush=True)
