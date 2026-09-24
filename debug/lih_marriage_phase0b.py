r"""Phase 0b of the LiH "marriage" build (2026-09-22; plan debug/lih_marriage_build_plan.md):
the EXACT ordered-integral prolate-Neumann operator (geovac/lih_r12ci/neumann_exact.py, switch
geovac.lih_r12ci.kernels.USE_EXACT_NEUMANN) -- gates G-a .. G-e plus the legacy bit-identity check.

Run from root:
    python debug/lih_marriage_phase0b.py [--skip-g0] [--skip-ge | --ge-log FILE] > debug/data/lih_marriage_phase0b.log 2>&1
    python debug/lih_marriage_phase0b.py --ge-worker > debug/data/lih_marriage_phase0b_ge.log   (G-e worker)

Sections
  0    LEGACY BIT-IDENTITY: the default (legacy) path vs the PRISTINE pre-edit capture
       debug/data/lih_marriage_phase0b_legacy_pristine.npz (24 arrays) -> max |Delta| must be 0.0.
       The capture's provenance is established against `git HEAD` by debug/lih_marriage_phase0b_headcheck.py
       (git-archive copy of the committed package vs the working tree, 27 arrays, 0.0 exactly, 2026-09-23).
  diag INTERPOLANT: the design assumption that the eta-moments g_l(xi) are represented by the degree-71
       Lagrange interpolant through the GL nodes, measured (isotropic STOs zeta = 1.0 .. 4.5, l = 0..2,
       400 off-node points) together with the n = 71 Legendre coefficient of e^{-2 zeta a (xi-1)}.
  G-a  Neumann self-Coulomb of isotropic 1s STOs vs 5 zeta/8, zeta in {1.0 (B), 1.6, 2.6875, 4.5 (A)}:
       exact rel <= 1e-8 (legacy values reported alongside: the wrong answers).
  G-b  Two-centre Coulomb closed forms: (aa|bb), (aa|ab), (bb|ab) -- both dressing directions of the
       exact operator vs the Hartree-dressed closed-form anchor (energy.py _hartree_1s, the C_aabb /
       C_aaab / C_bbab controls of stage2_E0): rel <= 1e-8.  NOTE energy.py's V_aabb / V_aaab / V_bbab
       are the f-GEMINAL integrals (Stage 1), not Coulomb; the Coulomb anchors are the Hartree ones.
       (ab|ab): exact vs legacy-72 vs the 400x160 Neumann of stage2_E0 vs importance-MC (report only).
  G-c  General-m anchors.  (i) CLOSED FORM: the mode-m density B_m = rho_cyl^m e^{-2 zeta r} is the
       cos(m phi) component of the solid-harmonic density r^m sin^m(theta) cos(m phi) e^{-2 zeta r},
       whose potential is (4 pi/(2m+1)) r^m Y [ r^{-2m-1} INT_0^r r'^{2m+2} e^{-2 zeta r'} dr'
       + INT_r^inf r' e^{-2 zeta r'} dr' ] (incomplete gamma); triangle.coul_mode_potential returns
       2 pi a^3 times the cos(m phi) coefficient.  m = 0..MMAX, centres A (zeta 2.6875) and B (1.0):
       rel <= 1e-8.  (ii) INDEPENDENT ROUTES for the real non-isotropic rho_(NO0,NO1) x rho_cyl^m
       (m = 0, 1, 2): (A) the integrated self-value E_m[B,B] from the LEGACY formula on scratch grids
       NXI = 144, 288, 576, 1152 Richardson-extrapolated in h^2 vs the exact operator at NXI = 72;
       (B) pointwise V^(m)(xi_i, eta_k) at 12 grid nodes vs a BRUTE-FORCE ordered radial integral in
       which the density is evaluated analytically at every sub-quadrature point (no interpolant;
       composite Gauss, panels graded geometrically toward xi = 1 on the P side and toward xi_i on
       the Q side; converged by comparing 40 vs 60 points per panel).  Both rel <= 1e-6.
       A pointwise Richardson of the legacy formula at a FIXED off-node target is reported as a
       diagnostic only: its leading error depends on where the target falls inside a scratch cell,
       which is not a smooth function of h, so the h^2 expansion (and the extrapolation) fails there
       -- the first run of this driver gated on it and read a 7e-4 residual as a FAIL.
  G-d  Gate G0 re-run on the exact path (debug/lih_marriage_phase0.py --exact): <V_ee> grid vs engine
       <= 1 mHa; <T>, <V_ne> unchanged; total vs E_trunc.
  G-e  (time-capped 15 min) geovac.lih_r12ci energy('linexp') / energy('exp') with the switch ON vs
       the banked legacy -7.9168 / -7.9420 and the same-geminal VMC references (report only).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
import warnings
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
from numpy.polynomial.legendre import leggauss

warnings.filterwarnings("ignore", category=DeprecationWarning)
from scipy.special import gammainc, lpmn, lqmn  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")
PRISTINE = os.path.join(DATA, "lih_marriage_phase0b_legacy_pristine.npz")
REF = os.path.join(DATA, "lih_marriage_phase0_ref.npz")
REF_EXACT = os.path.join(DATA, "lih_marriage_phase0_ref_exact.npz")
G0_LOG = os.path.join(DATA, "lih_marriage_phase0_exact.log")
GE_LOG = os.path.join(DATA, "lih_marriage_phase0b_ge.log")
GE_TIMEOUT_S = 15 * 60

GATE_A = 1e-8
GATE_B = 1e-8
GATE_C1 = 1e-8
GATE_C2 = 1e-6
GATE_D_MHA = 1.0


def _t(t0: float) -> str:
    return f"[{time.time() - t0:6.1f}s]"


def rel(x: float, y: float) -> float:
    return abs(x - y) / max(abs(y), 1e-300)


# =========================================================================== #
# closed-form anchor for the mode-m potential of B_m = rho_cyl^m e^{-2 zeta r}
# =========================================================================== #
def phi_solid_harmonic(r: np.ndarray, rho_cyl: np.ndarray, zeta: float, m: int) -> np.ndarray:
    """cos(m phi) coefficient of the Coulomb potential of rho(r) = r^m sin^m(theta) cos(m phi) e^{-2 zeta r}
    (= rho_cyl^m cos(m phi) e^{-2 zeta r}), a degree-m solid harmonic times a 1s radial factor."""
    s = 2 * m + 3
    lower = gammainc(s, 2.0 * zeta * r) * math.factorial(s - 1) / (2.0 * zeta) ** s   # INT_0^r r'^{2m+2} e^{-2z r'}
    upper = np.exp(-2.0 * zeta * r) * (r / (2.0 * zeta) + 1.0 / (2.0 * zeta) ** 2)  # INT_r^inf r' e^{-2z r'}
    rr = np.maximum(r, 1e-300)
    return (4.0 * np.pi / (2 * m + 1)) * rho_cyl ** m * (lower / rr ** (2 * m + 1) + upper)


# =========================================================================== #
# rho_(NO0,NO1) analytically on an arbitrary (XI, ETA) grid (from the Phase 0 artifact)
# =========================================================================== #
def make_rho01(z) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    prims = list(z['prims'])
    T_no = z['T_no']

    def rho01(XI: np.ndarray, ETA: np.ndarray) -> np.ndarray:
        pv = np.empty((len(prims),) + XI.shape)
        for i, p in enumerate(prims):
            g = np.polynomial.polynomial.polyval(ETA, p['eta'])
            pv[i] = p['N'] * XI ** p['j'] * g * np.exp(-p['alpha'] * XI)
        n0 = np.tensordot(T_no[:, 0], pv, axes=(0, 0))
        n1 = np.tensordot(T_no[:, 1], pv, axes=(0, 0))
        return n0 * n1
    return rho01


# =========================================================================== #
# LEGACY general-m mode potential on a scratch GL grid, evaluated at arbitrary target points
# (the same cumulative-GL formula as triangle.coul_mode_potential, kernel P_l^m(xi<)Q_l^m(xi>))
# =========================================================================== #
def legacy_mode_scratch(Bfun: Callable, m: int, NXI_s: int, NETA_s: int, LMAX: int, xi_max: float,
                        R: float, a: float, targets: Sequence[Tuple[float, float]]) -> Tuple[np.ndarray, float]:
    """Returns (V at the targets, integrated self-value E_m[B,B] on the scratch grid)."""
    xg, wxg = leggauss(NXI_s)
    xi = 1.0 + 0.5 * (xg + 1.0) * (xi_max - 1.0)
    wxi = 0.5 * (xi_max - 1.0) * wxg
    eta, weta = leggauss(NETA_s)
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    JAC = XI ** 2 - ETA ** 2
    B = Bfun(XI, ETA)
    WB = JAC * B
    Peta = np.array([lpmn(m, LMAX, e)[0][m] for e in eta]).T          # (LMAX+1, NETA_s), Ferrers
    Pxi = np.array([lpmn(m, LMAX, x)[0][m] for x in xi]).T            # (LMAX+1, NXI_s)
    Qxi = np.array([lqmn(m, LMAX, x)[0][m] for x in xi]).T
    norm = {l: (math.factorial(l - m) / math.factorial(l + m)) ** 2 for l in range(m, LMAX + 1)}
    Wm = (2 * np.pi) ** 2 if m == 0 else 2 * np.pi ** 2
    pref = (2.0 - (1.0 if m == 0 else 0.0)) * (2.0 / R) * a ** 6 * Wm
    mi = np.minimum.outer(np.arange(NXI_s), np.arange(NXI_s))
    ma = np.maximum.outer(np.arange(NXI_s), np.arange(NXI_s))
    Vt = np.zeros(len(targets))
    Vgrid = np.zeros((NXI_s, NETA_s))
    tx = np.array([t[0] for t in targets]); te = np.array([t[1] for t in targets])
    Pt_xi = np.array([lpmn(m, LMAX, x)[0][m] for x in tx]).T            # (LMAX+1, nt)
    Qt_xi = np.array([lqmn(m, LMAX, x)[0][m] for x in tx]).T
    Pt_eta = np.array([lpmn(m, LMAX, e)[0][m] for e in te]).T
    for l in range(m, LMAX + 1):
        gB = (WB * Peta[l][None, :]) @ weta                              # (NXI_s,)
        c = (-1) ** m * (2 * l + 1) * norm[l]
        # targets: kernel P(min) Q(max) with the off-node target xi_t
        lo = tx[:, None] <= xi[None, :]                                  # xi_t <= xi_j
        K_t = np.where(lo, Pt_xi[l][:, None] * Qxi[l][None, :], Pxi[l][None, :] * Qt_xi[l][:, None])
        rad_t = K_t @ (wxi * gB)
        Vt += c * rad_t * Pt_eta[l]
        Kl = Pxi[l][mi] * Qxi[l][ma]
        Vgrid += c * np.outer(Kl @ (wxi * gB), Peta[l])
    Vt *= pref; Vgrid *= pref
    Eself = float(np.sum(np.outer(wxi, weta) * JAC * B * Vgrid))
    return Vt, Eself


def brute_mode_potential_at(Bfun: Callable, m: int, targets: Sequence[Tuple[float, float]], LMAX: int,
                            xi_max: float, R: float, a: float, Peta_m: np.ndarray, ETA: np.ndarray,
                            WETA: np.ndarray, NORM: Dict[Tuple[int, int], float], npts: int = 60,
                            ratio: float = 2.0) -> np.ndarray:
    """V^(m)(xi_t, eta_t) of the mode density B (analytic callable on (XI, ETA)) by BRUTE-FORCE ordered
    radial integrals: for each target, radial_l = Q_l^m(xi_t) INT_1^{xi_t} P_l^m g_l + P_l^m(xi_t)
    INT_{xi_t}^{xi_max} Q_l^m g_l with g_l(x) = sum_eta w JAC B P_l^m(eta) evaluated at every Gauss
    sub-point (no Lagrange interpolant anywhere).  Panels graded geometrically in (x-1): toward 1 on
    the P side (the core density's scale), toward xi_t on the Q side (the log of Q_l^m).  Same
    lpmn/lqmn tables, prefactor and (-1)^m (2l+1) norm assembly as triangle.coul_mode_potential."""
    xg, wg = leggauss(npts)
    Wm = (2 * np.pi) ** 2 if m == 0 else 2 * np.pi ** 2
    pref = (2.0 - (1.0 if m == 0 else 0.0)) * (2.0 / R) * a ** 6 * Wm

    def gmom(xs: np.ndarray) -> np.ndarray:                     # (LMAX+1, nx)
        XI_, ETA_ = np.meshgrid(xs, ETA, indexing='ij')
        WB = (XI_ ** 2 - ETA_ ** 2) * Bfun(XI_, ETA_)
        return Peta_m @ (WB * WETA[None, :]).T

    def rule(bounds: List[float]) -> Tuple[np.ndarray, np.ndarray]:
        xs, ws = [], []
        for b0, b1 in zip(bounds[:-1], bounds[1:]):
            xs.append(0.5 * (b1 - b0) * (xg + 1.0) + b0); ws.append(0.5 * (b1 - b0) * wg)
        return np.concatenate(xs), np.concatenate(ws)

    out = np.zeros(len(targets))
    for it, (xt, et) in enumerate(targets):
        d = xt - 1.0
        bP = [xt]
        while bP[-1] - 1.0 > 1e-7 * d:
            bP.append(1.0 + (bP[-1] - 1.0) / ratio)
        bP.append(1.0); bP = bP[::-1]
        bQ = [xt]
        while True:
            nxt = 1.0 + (bQ[-1] - 1.0) * ratio
            if nxt >= xi_max:
                bQ.append(xi_max); break
            bQ.append(nxt)
        xP, wP = rule(bP); xQ, wQ = rule(bQ)
        Pp = np.array([lpmn(m, LMAX, x)[0][m] for x in xP]).T   # (LMAX+1, nP)
        Qq = np.array([lqmn(m, LMAX, x)[0][m] for x in xQ]).T
        IP = (Pp * gmom(xP)) @ wP; IQ = (Qq * gmom(xQ)) @ wQ
        radial = lqmn(m, LMAX, xt)[0][m] * IP + lpmn(m, LMAX, xt)[0][m] * IQ
        Pe = lpmn(m, LMAX, et)[0][m]
        out[it] = pref * sum((-1) ** m * (2 * l + 1) * NORM[(m, l)] * Pe[l] * radial[l] for l in range(m, LMAX + 1))
    return out


def richardson(vals: Dict[int, float]) -> Dict[str, float]:
    """h^2 Richardson on a doubling sequence {N: value}; returns the R1 chain and R2 (h^3, h^4 forms)."""
    Ns = sorted(vals)
    out: Dict[str, float] = {}
    r1 = {}
    for n1, n2 in zip(Ns[:-1], Ns[1:]):
        r1[n2] = (4.0 * vals[n2] - vals[n1]) / 3.0
        out[f"R1({n1},{n2})"] = r1[n2]
    ks = sorted(r1)
    if len(ks) >= 2:
        out["R2_h3"] = (8.0 * r1[ks[-1]] - r1[ks[-2]]) / 7.0
        out["R2_h4"] = (16.0 * r1[ks[-1]] - r1[ks[-2]]) / 15.0
    return out


# =========================================================================== #
# G-e worker: the 2-MO engine energies with the exact switch ON (separate process)
# =========================================================================== #
def ge_worker() -> int:
    import geovac.lih_r12ci.kernels as KG
    KG.USE_EXACT_NEUMANN = True
    # NOT `from geovac.lih_r12ci import energy`: once a submodule import (kernels above) has loaded
    # geovac.lih_r12ci.energy the package attribute is the MODULE, not the function __init__ exports.
    from geovac.lih_r12ci.assembly import energy
    t0 = time.time()
    for gem in ("linexp", "exp"):
        r = energy(gem)
        pieces = {k: (None if v is None else float(v)) for k, v in r.pieces.items()}
        rec = dict(geminal=gem, E_R12=float(r.E_R12), E0=float(r.E0), dE_mHa=float(r.dE_mHa),
                   sigma2=float(r.sigma2), h=float(r.h), g=float(r.g), variational=bool(r.variational),
                   pieces=pieces, wall_s=time.time() - t0)
        print("GE_RESULT " + json.dumps(rec), flush=True)
    print(f"GE_DONE wall {time.time() - t0:.0f} s", flush=True)
    return 0


def parse_ge_log(path: str) -> List[dict]:
    out = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith("GE_RESULT "):
                out.append(json.loads(line[len("GE_RESULT "):]))
    return out


# =========================================================================== #
def main(args: argparse.Namespace) -> str:
    T0 = time.time()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    print("=" * 100)
    print("LiH MARRIAGE -- PHASE 0b: exact ordered-integral prolate-Neumann operator on the 72x44 grid")
    print("=" * 100, flush=True)

    import geovac.lih_r12ci.kernels as KG
    from geovac.lih_r12ci.hVee import (neumann_potential, make_kernel_coul, d_aa, d_ab, d_bb, P00, rho_g,
                                       NXI, NETA)
    from geovac.lih_r12ci.triangle import coul_mode_potential, MMAX, LMAX, grid2_int
    from geovac.lih_r12ci.gVee import psi_coul, psi_yuk
    from geovac.lih_r12ci.hT import yukawa_pot_iso
    from geovac.lih_r12ci.energy import _hartree_1s, R, a, ZA, ZB, stage2_E0, _mc_abab_coulomb
    print(f"  package imported {_t(T0)}  (gVee builds the phi-kernels at import)", flush=True)
    assert KG.USE_EXACT_NEUMANN is False, "default switch must be False"
    grid_int = KG.grid_int
    rA = KG.rA.ravel(); rB = KG.rB.ravel(); RHO = KG.RHO_CYL.ravel()
    shape = (NXI, NETA)
    gates: Dict[str, bool] = {}

    # ------------------------------------------------------------------ 0. bit-identity
    print("\n(0) LEGACY BIT-IDENTITY vs the pristine pre-edit capture")
    z0 = np.load(REF, allow_pickle=True)
    pairs = [tuple(p) for p in z0['pairs']]
    rho01_grid = z0['rho_no'][pairs.index((0, 1))]
    prist = np.load(PRISTINE, allow_pickle=True)
    sto = {}
    for zeta, cen in ((1.0, 'B'), (1.6, 'A'), (2.6875, 'A'), (4.5, 'A')):
        r = rA if cen == 'A' else rB
        sto[f"rho1s_{zeta}"] = (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * r)
    now = {}
    now['NP_daa'] = neumann_potential(d_aa.reshape(shape))
    now['NP_dab'] = neumann_potential(d_ab.reshape(shape))
    now['NP_dbb'] = neumann_potential(d_bb.reshape(shape))
    for k, v in sto.items():
        now['NP_' + k] = neumann_potential(v.reshape(shape))
    now['NP_rho01'] = neumann_potential(rho01_grid.reshape(shape))
    for m in range(MMAX + 1):
        Bm = rho01_grid * (np.maximum((KG.Xg.ravel() ** 2 - 1.0) * (1.0 - KG.Eg.ravel() ** 2), 0.0)) ** (0.5 * m)
        now[f'CMP_rho01_m{m}'] = coul_mode_potential(Bm.reshape(shape), m)
    for m in range(3):
        now[f'CMP_P00_m{m}'] = coul_mode_potential(P00.reshape(shape), m)
        now[f'CMP_rho_m{m}'] = coul_mode_potential(rho_g.reshape(shape), m)
    now['PSI_coul_dab'] = psi_coul(d_ab)
    now['PSI_yuk_dab'] = psi_yuk(d_ab, KG.GAM)
    now['PSI_yuk_daa'] = psi_yuk(d_aa, KG.GAM)
    Kc_leg = make_kernel_coul()
    now['W_coul'] = Kc_leg['W']
    now['Psi_coul_ab'] = Kc_leg['Psi']['ab']
    worst = 0.0
    for k in prist.files:
        d = float(np.max(np.abs(now[k] - prist[k])))
        worst = max(worst, d)
        print(f"    {k:18s} max|Delta| = {d:.1e}")
    gates['bit-identity'] = (worst == 0.0)
    print(f"  => legacy path vs pristine: max|Delta| over {len(prist.files)} arrays = {worst:.1e}  "
          f"-> {'BIT-IDENTICAL' if worst == 0.0 else 'DIFFERS'}   {_t(T0)}", flush=True)

    # operator build info
    t1 = time.time()
    op = KG.exact_neumann(LMAX, MMAX)
    print(f"  exact operator: A shape {op.A.shape} ({op.A.nbytes / 1e6:.1f} MB), n_sub={op.n_sub}, "
          f"n_gauss={op.n_gauss}, ratio={op.ratio}, build {time.time() - t1:.2f} s")

    # ------------------------------------------------------------------ diag: interpolant accuracy
    print("\n(diag) degree-71 GL-node interpolant of the eta-moments g_l(xi) (the design assumption), isotropic STOs")
    from geovac.lih_r12ci.neumann_exact import barycentric_weights, lagrange_matrix
    from scipy.special import eval_legendre
    scale = 2.0 / (KG.xi_max - 1.0)
    t_nodes = scale * (KG.XI - 1.0) - 1.0
    xi_off = 1.0 + 0.5 * (leggauss(400)[0] + 1.0) * (KG.xi_max - 1.0)          # never coincide with the 72 nodes
    Lam = lagrange_matrix(t_nodes, barycentric_weights(t_nodes), scale * (xi_off - 1.0) - 1.0)
    xg2, wg2 = leggauss(2000)
    print(f"  {'zeta':>7s} {'centre':>6s} {'interp g0':>10s} {'interp g1':>10s} {'interp g2':>10s} {'|c71/c0| of e^(-2 zeta a (xi-1))':>32s}   (rel to max|g_l|; 400 off-node points)")
    for zeta, cen in ((1.0, 'B'), (1.6, 'A'), (2.6875, 'A'), (4.5, 'A')):
        sgn = 1.0 if cen == 'A' else -1.0

        def gl(xi1d, l, zeta=zeta, sgn=sgn):
            XI_, ETA_ = np.meshgrid(xi1d, KG.ETA, indexing='ij')
            W_ = (XI_ ** 2 - ETA_ ** 2) * (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * a * (XI_ + sgn * ETA_))
            return (W_ * eval_legendre(l, KG.ETA)[None, :]) @ KG.WETA
        errs = [np.max(np.abs(Lam @ gl(KG.XI, l) - gl(xi_off, l))) / np.max(np.abs(gl(xi_off, l))) for l in range(3)]
        f2 = np.exp(-2.0 * zeta * a * 0.5 * (xg2 + 1.0) * (KG.xi_max - 1.0))
        c0 = 0.5 * np.sum(wg2 * f2); c71 = 71.5 * np.sum(wg2 * f2 * eval_legendre(71, xg2))
        print(f"  {zeta:7.4f} {cen:>6s} {errs[0]:10.1e} {errs[1]:10.1e} {errs[2]:10.1e} {abs(c71 / c0):32.1e}")
    print("  (the design memo's '~1e-13 by n=72' is NOT what is measured at zeta=4.5: the n=71 coefficient is ~2e-9 and the"
          " interpolant error ~2e-11 relative to the peak; the integrals below are exact to 1e-13 regardless because the"
          " error is oscillatory and the moments are weighted toward xi ~ 1)", flush=True)

    # ------------------------------------------------------------------ G-a
    print("\n(G-a) Neumann self-Coulomb of isotropic 1s STOs vs 5 zeta/8 (exact operator vs legacy)")
    print(f"  {'zeta':>8s} {'centre':>6s} {'J_exact':>14s} {'rel_exact':>11s} {'J_legacy':>14s} {'rel_legacy':>11s} {'Hartree-dress rel':>18s}")
    worst_a = 0.0
    for zeta, cen in ((1.0, 'B'), (1.6, 'A'), (2.6875, 'A'), (4.5, 'A')):
        r = rA if cen == 'A' else rB
        rho = sto[f"rho1s_{zeta}"]
        Ve = neumann_potential(rho.reshape(shape), exact=True).reshape(-1)
        Vl = neumann_potential(rho.reshape(shape), exact=False).reshape(-1)
        Je = grid_int(rho * Ve); Jl = grid_int(rho * Vl); Jh = grid_int(rho * _hartree_1s(r, zeta))
        ref = 5 * zeta / 8
        worst_a = max(worst_a, rel(Je, ref))
        print(f"  {zeta:8.4f} {cen:>6s} {Je:14.10f} {(Je - ref) / ref:+11.2e} {Jl:14.10f} {(Jl - ref) / ref:+11.2e} "
              f"{(Jh - ref) / ref:+18.2e}")
    gates['G-a'] = worst_a <= GATE_A
    print(f"  => G-a max exact rel = {worst_a:.2e} (gate {GATE_A:.0e}) -> {'PASS' if gates['G-a'] else 'FAIL'}", flush=True)

    # ------------------------------------------------------------------ G-b
    print("\n(G-b) two-centre Coulomb closed forms (Hartree-dressed anchors) vs the exact operator, both dressing directions")
    VA_h = _hartree_1s(rA, ZA); VB_h = _hartree_1s(rB, ZB)
    NPe = lambda dd: neumann_potential(dd.reshape(shape), exact=True).reshape(-1)
    NPl = lambda dd: neumann_potential(dd.reshape(shape), exact=False).reshape(-1)
    Ve_aa, Ve_bb, Ve_ab = NPe(d_aa), NPe(d_bb), NPe(d_ab)
    Vl_aa, Vl_bb, Vl_ab = NPl(d_aa), NPl(d_bb), NPl(d_ab)
    rows_b = [
        ("(aa|bb)", grid_int(d_bb * VA_h), grid_int(d_aa * VB_h), grid_int(d_bb * Ve_aa), grid_int(d_aa * Ve_bb),
         grid_int(d_bb * Vl_aa), grid_int(d_aa * Vl_bb)),
        ("(aa|ab)", grid_int(d_ab * VA_h), None, grid_int(d_ab * Ve_aa), grid_int(d_aa * Ve_ab),
         grid_int(d_ab * Vl_aa), grid_int(d_aa * Vl_ab)),
        ("(bb|ab)", grid_int(d_ab * VB_h), None, grid_int(d_ab * Ve_bb), grid_int(d_bb * Ve_ab),
         grid_int(d_ab * Vl_bb), grid_int(d_bb * Vl_ab)),
    ]
    worst_b = 0.0
    print(f"  {'integral':10s} {'anchor(Hartree)':>16s} {'anchor alt':>12s} {'exact dress-iso':>17s} {'rel':>9s} "
          f"{'exact dress-other':>18s} {'rel':>9s} {'legacy iso rel':>15s} {'legacy other rel':>17s}")
    for name, anc, anc2, e1, e2, l1, l2 in rows_b:
        r1, r2 = rel(e1, anc), rel(e2, anc)
        worst_b = max(worst_b, r1, r2)
        alt = f"{anc2:12.9f}" if anc2 is not None else f"{'':12s}"
        print(f"  {name:10s} {anc:16.10f} {alt} {e1:17.10f} {r1:9.1e} {e2:18.10f} {r2:9.1e} "
              f"{rel(l1, anc):15.1e} {rel(l2, anc):17.1e}")
        if anc2 is not None:
            print(f"    (anchor self-consistency: dress-A vs dress-B Hartree {rel(anc, anc2):.1e})")
    gates['G-b'] = worst_b <= GATE_B
    print(f"  => G-b max rel (6 anchored values) = {worst_b:.2e} (gate {GATE_B:.0e}) -> {'PASS' if gates['G-b'] else 'FAIL'}")
    # (ab|ab): report only
    abab_e = grid_int(d_ab * Ve_ab); abab_l = grid_int(d_ab * Vl_ab)
    t1 = time.time()
    _, info2 = stage2_E0(mc_check=False)
    abab_400 = info2['C_abab']
    mc_v, mc_e = _mc_abab_coulomb()
    print(f"  (ab|ab) Coulomb: exact-72 {abab_e:.9f}  legacy-72 {abab_l:.9f} (rel to exact {rel(abab_l, abab_e):+.1e})  "
          f"stage2 400x160 legacy Neumann {abab_400:.9f} (rel {rel(abab_400, abab_e):+.1e})  "
          f"importance-MC {mc_v:.6f}+/-{mc_e:.1e} ({(mc_v - abab_e) / mc_e:+.1f} sigma)   [{time.time() - t1:.0f}s]")
    print("  (kernels.VABAB_REF = 0.004350 is the (ab|f|ab) GEMINAL integral, not Coulomb -- not an anchor here)")
    KG.USE_EXACT_NEUMANN = True
    Kc_ex = make_kernel_coul()
    KG.USE_EXACT_NEUMANN = False
    We, Wl = Kc_ex['W'], Kc_leg['W']
    print(f"  make_kernel_coul W^coul symmetry max|W-W^T|: legacy {np.max(np.abs(Wl - Wl.T)):.1e}  exact {np.max(np.abs(We - We.T)):.1e};  "
          f"max|W_exact - W_legacy| = {np.max(np.abs(We - Wl)):.2e}", flush=True)

    # ------------------------------------------------------------------ G-c (i): closed-form solid-harmonic anchors
    print("\n(G-c i) general-m CLOSED-FORM anchors: B_m = rho_cyl^m e^{-2 zeta r} (solid harmonic x 1s), m = 0..MMAX")
    # anchor scaffolding check: m=0 vs Hartree (pi/zeta^3 normalisation)
    ph0 = phi_solid_harmonic(rA, RHO, ZA, 0); hh = (np.pi / ZA ** 3) * _hartree_1s(rA, ZA)
    print(f"  anchor formula m=0 vs (pi/zeta^3) Hartree: max rel {np.max(np.abs(ph0 - hh) / np.abs(hh)):.1e}")
    print(f"  {'m':>2s} {'centre':>6s} {'zeta':>7s} {'sup|Ve-Va|/sup|Va|':>20s} {'E_m self exact':>16s} {'E_m self anchor':>16s} "
          f"{'rel':>9s} {'legacy sup-rel':>15s} {'legacy E_m rel':>15s}")
    worst_c1 = 0.0
    for zeta, cen in ((2.6875, 'A'), (1.0, 'B')):
        r = rA if cen == 'A' else rB
        for m in range(MMAX + 1):
            Bm = (RHO ** m) * np.exp(-2.0 * zeta * r)
            Va = (2 * np.pi * a ** 3) * phi_solid_harmonic(r, RHO, zeta, m)
            Ve = coul_mode_potential(Bm.reshape(shape), m, exact=True).reshape(-1)
            Vl = coul_mode_potential(Bm.reshape(shape), m, exact=False).reshape(-1)
            sup_e = np.max(np.abs(Ve - Va)) / np.max(np.abs(Va))
            sup_l = np.max(np.abs(Vl - Va)) / np.max(np.abs(Va))
            Ee = grid2_int((Bm * Ve).reshape(shape)); Ea = grid2_int((Bm * Va).reshape(shape))
            El = grid2_int((Bm * Vl).reshape(shape))
            worst_c1 = max(worst_c1, sup_e, rel(Ee, Ea))
            print(f"  {m:2d} {cen:>6s} {zeta:7.4f} {sup_e:20.2e} {Ee:16.10e} {Ea:16.10e} {rel(Ee, Ea):9.1e} {sup_l:15.1e} {rel(El, Ea):15.1e}")
    gates['G-c(i)'] = worst_c1 <= GATE_C1
    print(f"  => G-c(i) max rel (sup-norm and integrated, m=0..{MMAX}, both centres) = {worst_c1:.2e} (gate {GATE_C1:.0e}) "
          f"-> {'PASS' if gates['G-c(i)'] else 'FAIL'}   {_t(T0)}", flush=True)

    # ------------------------------------------------------------------ G-c (ii): independent routes, rho_(NO0,NO1) x rho_cyl^m
    print("\n(G-c ii) INDEPENDENT ROUTES for the real non-isotropic rho_(NO0,NO1) x rho_cyl^m (m = 0, 1, 2)")
    print("  (A) integrated E_m[B,B]: exact-72 vs the LEGACY formula on scratch grids NXI = 144..1152 (NETA 88), Richardson in h^2")
    print("  (B) pointwise V^(m) at 12 nodes: exact-72 vs a BRUTE-FORCE ordered integral with the density evaluated analytically")
    print("      at every sub-point (no interpolant; composite Gauss graded toward xi=1 / xi_i; 60 pts/panel, checked vs 40)")
    from geovac.lih_r12ci.triangle import Peta as _Peta, _NORM as _TNORM
    rho01 = make_rho01(z0)
    chk = np.max(np.abs(rho01(KG.Xg, KG.Eg).ravel() - rho01_grid)) / np.max(np.abs(rho01_grid))
    print(f"  analytic rho01 on the 72 grid vs the artifact: max rel {chk:.1e}")
    ix = [0, 4, 12, 24, 40, 60]; ie = [8, 30]
    targets = [(float(KG.XI[i]), float(KG.ETA[k])) for i in ix for k in ie]
    NS = [144, 288, 576, 1152]
    worst_c2 = 0.0
    for m in (0, 1, 2):
        def Bfun(XI, ETA, m=m):
            rc = a * np.sqrt(np.maximum((XI ** 2 - 1.0) * (1.0 - ETA ** 2), 0.0))
            return rho01(XI, ETA) * rc ** m
        Bm72 = Bfun(KG.Xg, KG.Eg)
        Ve = coul_mode_potential(Bm72, m, exact=True)
        Vl72 = coul_mode_potential(Bm72, m, exact=False)
        Ee = grid2_int(Bm72 * Ve)
        Vt_e = np.array([Ve[i, k] for i in ix for k in ie])
        Vt_l72 = np.array([Vl72[i, k] for i in ix for k in ie])
        # (A) integrated, legacy on scratch grids + Richardson
        vals_t: Dict[int, np.ndarray] = {}; vals_E: Dict[int, float] = {}
        for N in NS:
            Vt, Es = legacy_mode_scratch(Bfun, m, N, 88, LMAX, float(KG.xi_max), R, a, targets)
            vals_t[N] = Vt; vals_E[N] = Es
        rE = richardson(vals_E)
        keyR1 = f"R1({NS[-2]},{NS[-1]})"
        rA_int = rel(rE[keyR1], Ee)
        print(f"  m={m} (A) exact-72 E_m self = {Ee:.10e};  legacy " +
              "  ".join(f"N{N}={vals_E[N]:.8e}" for N in NS))
        print(f"        Richardson " + "  ".join(f"{k}={v:.10e}" for k, v in rE.items()))
        print(f"        |{keyR1} - exact|/|exact| = {rA_int:.2e}   (R2_h3 {rel(rE['R2_h3'], Ee):.2e}, R2_h4 {rel(rE['R2_h4'], Ee):.2e})")
        # (B) pointwise, brute-force ordered integral with the analytic density
        t1 = time.time()
        Vb60 = brute_mode_potential_at(Bfun, m, targets, LMAX, float(KG.xi_max), R, a, _Peta[m], KG.ETA, KG.WETA, _TNORM, npts=60)
        Vb40 = brute_mode_potential_at(Bfun, m, targets, LMAX, float(KG.xi_max), R, a, _Peta[m], KG.ETA, KG.WETA, _TNORM, npts=40)
        conv = np.max(np.abs(Vb60 - Vb40) / np.abs(Vb60))
        pwB = np.abs(Vt_e - Vb60) / np.abs(Vb60)
        pwL = np.abs(Vt_l72 - Vb60) / np.abs(Vb60)
        print(f"  m={m} (B) brute-force convergence (60 vs 40 pts/panel) max rel {conv:.1e};  exact-72 vs brute max rel {pwB.max():.2e};  "
              f"legacy-72 vs brute max rel {pwL.max():.2e}   [{time.time() - t1:.0f}s]")
        print(f"        per node (xi_i, eta_k) [exact rel | legacy rel]: " +
              " ".join(f"({KG.XI[i]:.2f},{KG.ETA[k]:+.2f}):{pwB[j]:.0e}|{pwL[j]:.0e}" for j, (i, k) in enumerate((i, k) for i in ix for k in ie)))
        # diagnostic: the fixed-off-node-target pointwise Richardson of the legacy formula is NOT an h^2 sequence
        seq = [np.max(np.abs(vals_t[N] - Vb60) / np.abs(Vb60)) for N in NS]
        ratios = [seq[j] / seq[j + 1] for j in range(len(seq) - 1)]
        r1_t = (4.0 * vals_t[NS[-1]] - vals_t[NS[-2]]) / 3.0
        print(f"        (diag) legacy scratch pointwise max rel vs brute: " + " ".join(f"N{N}:{s:.1e}" for N, s in zip(NS, seq)) +
              f";  ratios per doubling {', '.join(f'{r:.2f}' for r in ratios)} (h^2 would be 4.00) -> R1({NS[-2]},{NS[-1]}) "
              f"max rel {np.max(np.abs(r1_t - Vb60) / np.abs(Vb60)):.1e}, not gated")
        worst_c2 = max(worst_c2, rA_int, float(pwB.max()))
    gates['G-c(ii)'] = worst_c2 <= GATE_C2
    print(f"  => G-c(ii) max rel over (A) integrated Richardson and (B) pointwise brute-force, m=0,1,2 = {worst_c2:.2e} "
          f"(gate {GATE_C2:.0e}) -> {'PASS' if gates['G-c(ii)'] else 'FAIL'}   {_t(T0)}", flush=True)

    # ------------------------------------------------------------------ psi_yuk diagnostic (for G-e / Phase 1)
    print("\n(diag) psi_yuk(rho_aa) = Neumann - smooth vs the closed-radial yukawa_pot_iso (the gT.py:92-102 ~0.5% bias)")
    ykc = yukawa_pot_iso(rA, ZA, KG.GAM)
    msk = d_aa > 1e-6 * d_aa.max()
    for flag in (False, True):
        KG.USE_EXACT_NEUMANN = flag
        yk = psi_yuk(d_aa, KG.GAM)
        wrel = np.abs(yk - ykc)[msk].mean() / np.abs(ykc[msk]).mean()
        Iy = grid_int(d_aa * yk); Iyc = grid_int(d_aa * ykc)
        print(f"    {'exact ' if flag else 'legacy'}: weighted rel {wrel:.2e};  I_Y[aa,aa] {Iy:.8f} vs closed {Iyc:.8f} (rel {rel(Iy, Iyc):.1e})")
    KG.USE_EXACT_NEUMANN = False

    # ------------------------------------------------------------------ G-d
    print("\n(G-d) gate G0 re-run on the exact path (debug/lih_marriage_phase0.py --exact)")
    if not args.skip_g0:
        t1 = time.time()
        with open(G0_LOG, 'w', encoding='utf-8') as fh:
            rc = subprocess.call([sys.executable, os.path.join(HERE, "lih_marriage_phase0.py"), "--exact"],
                                 stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)
        print(f"  phase0 --exact exit={rc}  [{time.time() - t1:.0f}s]  log {G0_LOG}")
    else:
        print(f"  (--skip-g0: reading the existing artifact {REF_EXACT} / log {G0_LOG})")
    if os.path.exists(REF_EXACT):
        ze = np.load(REF_EXACT, allow_pickle=True)
        tab = ze['g0_table']; rows = list(ze['g0_rows'])
        print(f"  exact_neumann flag in artifact: {bool(ze['exact_neumann'])};  E_trunc = {float(ze['E_trunc']):.6f}")
        print(f"  {'term':22s} {'grid':>14s} {'engine':>14s} {'Delta (mHa)':>13s}")
        for name, (gv, ev, dm) in zip(rows, tab):
            print(f"  {name:22s} {gv:14.8f} {ev:14.8f} {dm:+13.6f}")
        dVee = abs(float(tab[3, 2])); dT = abs(float(tab[0, 2])); dV = abs(float(tab[2, 2]))
        gates['G-d'] = (dVee <= GATE_D_MHA) and (dT <= GATE_D_MHA) and (dV <= GATE_D_MHA)
        print(f"  => G-d dT={dT:.6f} dVne={dV:.6f} dVee={dVee:.6f} mHa (gate {GATE_D_MHA}); total grid {tab[4, 0]:.8f} vs "
              f"E_trunc {float(ze['E_trunc']):.8f} -> {'PASS' if gates['G-d'] else 'FAIL'}")
        if os.path.exists(REF):
            zl = np.load(REF, allow_pickle=True)
            print(f"  (legacy-path Phase 0 for comparison: dVee = {float(zl['g0_table'][3, 2]):+.4f} mHa)")
        # the self-Coulomb control lines + the per-NO diagnostics from the log
        try:
            with open(G0_LOG, 'r', encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    if ('J_Neumann-5z/8' in line and ('coreLi ' in line or 'coreLi2' in line)) or \
                       'max|J_grid-J_engine|' in line or line.startswith('    NO0:'):
                        print("  " + line.rstrip())
        except OSError:
            pass
    else:
        gates['G-d'] = False
        print("  artifact missing -> FAIL")

    # ------------------------------------------------------------------ G-e
    print("\n(G-e) 2-MO engine energies with the exact switch ON (time-capped 15 min)")
    from geovac.lih_r12ci.assembly import ANALYTIC_REF, VMC_TARGETS
    ge_recs: List[dict] = []
    if args.ge_log:
        ge_recs = parse_ge_log(args.ge_log)
        print(f"  (parsed worker log {args.ge_log}: {len(ge_recs)} results)")
    elif args.skip_ge:
        print("  (--skip-ge)")
    else:
        t1 = time.time()
        try:
            with open(GE_LOG, 'w', encoding='utf-8') as fh:
                subprocess.run([sys.executable, os.path.abspath(__file__), "--ge-worker"], stdout=fh,
                               stderr=subprocess.STDOUT, cwd=ROOT, timeout=GE_TIMEOUT_S, check=False)
            ge_recs = parse_ge_log(GE_LOG)
            print(f"  worker done [{time.time() - t1:.0f}s], {len(ge_recs)} results (log {GE_LOG})")
        except subprocess.TimeoutExpired:
            ge_recs = parse_ge_log(GE_LOG) if os.path.exists(GE_LOG) else []
            print(f"  worker TIMED OUT after {GE_TIMEOUT_S} s; partial results: {len(ge_recs)}")
    for rec in ge_recs:
        gem = rec['geminal']; ref = ANALYTIC_REF[gem]; vmc = VMC_TARGETS[gem]
        print(f"  {gem}: E_R12 exact-switch = {rec['E_R12']:.6f}  (legacy banked {ref['E_R12']:.6f}, shift {(rec['E_R12'] - ref['E_R12']) * 1e3:+.3f} mHa; "
              f"VMC {vmc['E_R12']:.5f}, dev {(rec['E_R12'] - vmc['E_R12']) * 1e3:+.3f} mHa)  dE={rec['dE_mHa']:+.2f} mHa (VMC dE {vmc['dE_mHa']})  "
              f"variational={rec['variational']}  [{rec['wall_s']:.0f}s]")
        print(f"      sigma2 {rec['sigma2']:.6f} (legacy {ref['sigma2']:.6f}, VMC {vmc['sigma2']})  h {rec['h']:+.6f} (legacy {ref['h']:+.6f}, VMC {vmc['h']})  "
              f"g {rec['g']:+.6f} (legacy {ref['g']:+.6f}, VMC {vmc['g']})")
        pcs = rec['pieces']
        print("      pieces (exact - legacy): " + "  ".join(f"{k}:{pcs[k]:+.5f}({(pcs[k] - ref[k]) * 1e3:+.2f}m)"
                                                        for k in ('h_T', 'h_Vne', 'h_Vee', 'g_T', 'g_Vne', 'g_Vee')))
        if pcs.get('CovFY') is not None:
            print(f"      CovFY {pcs['CovFY']:+.5f}" + (f"  CovFE {pcs['CovFE']:+.5f}" if pcs.get('CovFE') is not None else ""))
    if not ge_recs:
        print("  (no G-e results)")

    # ------------------------------------------------------------------ verdict
    print("\n" + "=" * 100)
    for k, v in gates.items():
        print(f"  {k:14s} {'PASS' if v else 'FAIL'}")
    verdict = "GO" if all(gates.values()) else "STOP"
    print(f"VERDICT (Phase 0b, exact Neumann operator on the 72x44 grid): {verdict}")
    print(f"wall time {time.time() - T0:.0f} s")
    print("=" * 100, flush=True)
    return verdict


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--skip-g0', action='store_true', help='read the existing --exact artifact instead of re-running phase0')
    ap.add_argument('--skip-ge', action='store_true', help='skip the 2-MO engine energies')
    ap.add_argument('--ge-log', default=None, help='parse an existing --ge-worker log instead of running it')
    ap.add_argument('--ge-worker', action='store_true', help='(internal) run the energies with the exact switch ON')
    a_ = ap.parse_args()
    if a_.ge_worker:
        sys.exit(ge_worker())
    main(a_)
