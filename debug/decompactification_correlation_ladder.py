"""Decompactification correlation ladder (DIAGNOSTIC, measurement only).

Question (PI, 2026-09-06): does electron-electron interaction move the two-center
"decompactification front" ONLY through screening (Z -> Z_eff), or in some other way?

Background (debug/sprint_decompactification_R_sweep_memo.md): for one-electron H2+
the per-shell front -- the R at which the cross-center overlap of the ns_A/ns_B pair
reaches S = 1/sqrt2 (projector principal angle 45 deg, max ||[P_A,P_B]||) -- is set by
the orbital EXPONENTIAL DECAY LENGTH ell_n = n/Z, not the mean radius n^2/Z.
Paper 60 sec:manyelectron: the k-electron configuration overlap is the k-th compound
matrix of the one-electron overlap, so at the overlap level the front is one-electron
by identity.  Hypothesis under test: the front moves by Z -> Z_eff and by NOTHING else.

THE LADDER
 RUNG 0 (a) one-electron heteronuclear law: R*(nA,ZA,nB,ZB) where <nA s_A|nB s_B>=1/sqrt2
            (exact Mulliken/Ruedenberg overlap, aha_t1_core.two_center_s_overlap).  Fit
            R* = c (ell_A + ell_B), ell = n/Z, against max / geometric / harmonic / free-power
            alternatives.  Plus a dense 1s-1s exponent-ratio scan (scale-free form).
        (b) HeH2+ (Z_A=2, Z_B=1) eta-equation l-mixing vs R (heteronuclear term b = R(Z_B-Z_A)),
            c^2(R) from the exact prolate solver; compared to the homonuclear H2+ curve.
 RUNG 1  RHF on H2 and HeH+ (screening only).
 RUNG 2  singlet FCI on the same (screening + correlation).
         Substrate: per-center EVEN-TEMPERED 1s-STO basis (5 exponents/center; doubled = 10
         interleaved over the same span) so the per-center decay length is a variational
         OUTPUT, not an input (the HARD CONTROL).  One-electron integrals: closed form
         (same center) / prolate Gauss-Legendre (cross center, ~1e-9); ERIs: multipole
         expansion about A on a graded radial grid (the Paper 60 route, validated ~1e-4).
         At each R, from the one-body density:
          (i)  ell_eff per center: route A = log-slope of the dominant natural orbital's
               A-centered component over its [r50, r95] norm window; route B = moment
               zeta = 3/(2<r>) (exact for a pure 1s); route C = moment of the AA-block density.
          (ii) measured front: M1 = occupation-weighted |cos| of the principal angle between
               the A- and B-centered components of each natural orbital (= the compound-matrix
               object); M2 = signed normalized cross-center coherence
               Tr(P_AB S_BA)/sqrt(Tr(P_AA S_AA) Tr(P_BB S_BB)) (= M1 for a single MO).
          (iii) predicted front: the exact 1s(zeta_A)-1s(zeta_B) two-center overlap at the
               empirical zeta_eff(R) -- the Rung-0 one-electron law at Z_eff.
          (iv) residual (ii)-(iii) vs R, per rung, per molecule; R* = 1/sqrt2 crossings.

Guardrail (CLAUDE.md 3.5, Papers 8-9) acknowledged: genuine two-center machinery, R-dependence
measured; no single-center/shared-p0 encoding proposed, no binding claim.
Clean room: driver in debug/, data in debug/data/.  No paper / CLAUDE.md / test edits.

Usage:  python debug/decompactification_correlation_ladder.py [--quick] [--no-double]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from math import factorial, pi, sqrt

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import eval_legendre

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)


def _load(mod_name: str, rel_path: str):
    spec = importlib.util.spec_from_file_location(mod_name, os.path.join(REPO, rel_path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


core = _load("aha_t1_core", "debug/aha_t1_core.py")
two_center_s_overlap = core.two_center_s_overlap        # exact hydrogenic s-s (any a_A, a_B)

SQRT_HALF = 1.0 / sqrt(2.0)
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
DATA_DIR = os.path.join(REPO, "debug", "data")


def log(msg: str) -> None:
    print(msg, flush=True)


# =============================================================== RUNG 0 (a)
def S_1e(nA: int, ZA: float, nB: int, ZB: float, R: float) -> float:
    """<nA s (Z_A) on A | nB s (Z_B) on B> exact, hydrogenic decay a = Z/n."""
    return two_center_s_overlap(nA, ZA / nA, nB, ZB / nB, R)


def crossing_from_above(fun, R_lo: float = 1e-3, R_hi: float = 400.0):
    """Two crossings of |fun(R)| (fun = a two-center overlap):
      R_abs : largest R where |S| crosses 1/sqrt2 from above (the 45-deg principal angle;
              None when |S| never exceeds 1/sqrt2 -- no decompactification point exists);
      R_rel : R past the peak where |S| falls to |S|_max/sqrt2 (tail-reach scale; equals
              R_abs when S(0)=1).
    Returns (R_abs, S_max, R_at_max, R_rel)."""
    Rs = np.geomspace(max(R_lo, 1e-2), R_hi, 600)
    v = np.abs(np.array([fun(R) for R in Rs], float))
    v[~np.isfinite(v)] = -1.0                      # Mulliken B-auxiliary can cancel catastrophically at tiny R
    imax = int(np.argmax(v))
    Smax = float(v[imax])

    def last_cross(level):
        above = v > level
        if not above.any() or above[-1]:
            return None
        i_last = int(np.where(above)[0].max())
        return float(brentq(lambda R: abs(fun(R)) - level, Rs[i_last], Rs[i_last + 1], xtol=1e-11))

    R_abs = last_cross(SQRT_HALF)
    R_rel = last_cross(Smax * SQRT_HALF) if Smax > 1e-6 else None
    return R_abs, Smax, float(Rs[imax]), R_rel


def fit_laws(lA: np.ndarray, lB: np.ndarray, Rs: np.ndarray) -> dict:
    """Fit R* against candidate combinations of the two decay lengths (log-space)."""
    cands = {
        "additive  c(lA+lB)": lA + lB,
        "max       c max(lA,lB)": np.maximum(lA, lB),
        "geometric c 2sqrt(lA lB)": 2.0 * np.sqrt(lA * lB),
        "harmonic  c 2lAlB/(lA+lB)": 2.0 * lA * lB / (lA + lB),
    }
    out = {}
    for name, X in cands.items():
        lr = np.log(Rs) - np.log(X)
        c = float(np.exp(lr.mean()))
        res = lr - lr.mean()
        out[name] = {"c": c, "rms_log_resid": float(np.sqrt((res ** 2).mean())),
                     "max_abs_pct_dev": float(100.0 * (np.exp(np.abs(res)).max() - 1.0))}
    # free power law, symmetrized in (A,B): log R* = log c + p log lA + q log lB
    Xd = np.vstack([np.column_stack([np.ones_like(lA), np.log(lA), np.log(lB)]),
                    np.column_stack([np.ones_like(lA), np.log(lB), np.log(lA)])])
    yd = np.concatenate([np.log(Rs), np.log(Rs)])
    coef, *_ = np.linalg.lstsq(Xd, yd, rcond=None)
    res = yd - Xd @ coef
    out["free power  c lA^p lB^q (sym)"] = {"c": float(np.exp(coef[0])), "p": float(coef[1]),
                                          "q": float(coef[2]),
                                          "rms_log_resid": float(np.sqrt((res ** 2).mean()))}
    return out


def rung0a() -> dict:
    log("\n=== RUNG 0(a): one-electron heteronuclear front law ===")
    charge_pairs = [(1, 1), (2, 1), (3, 1), (2, 2)]
    ns = [1, 2, 3, 4]
    rows = []
    for ZA, ZB in charge_pairs:
        for nA in ns:
            for nB in ns:
                Rc, Smax, Rmax, Rrel = crossing_from_above(lambda R: S_1e(nA, ZA, nB, ZB, R))
                rows.append({"ZA": ZA, "ZB": ZB, "nA": nA, "nB": nB,
                             "ellA": nA / ZA, "ellB": nB / ZB,
                             "R_star": Rc, "S_max": Smax, "R_at_Smax": Rmax, "R_rel": Rrel})

    def fits_of(sub, key):
        sub = [r for r in sub if r[key] is not None]
        if len(sub) < 3:
            return {}
        return fit_laws(np.array([r["ellA"] for r in sub]), np.array([r["ellB"] for r in sub]),
                        np.array([r[key] for r in sub]))

    def show(title, fits):
        log(f"  {title}")
        for k, v in fits.items():
            log(f"    {k:32s} c={v['c']:.4f}  rms_log={v['rms_log_resid']:.4f}"
                + (f"  p={v['p']:.3f} q={v['q']:.3f}" if 'p' in v else f"  max_dev={v['max_abs_pct_dev']:.1f}%"))

    have = [r for r in rows if r["R_star"] is not None]
    same_n = [r for r in rows if r["nA"] == r["nB"]]
    fits = {
        "abs_all": fits_of(rows, "R_star"),
        "abs_same_n": fits_of(same_n, "R_star"),
        "rel_all": fits_of(rows, "R_rel"),
        "rel_same_n": fits_of(same_n, "R_rel"),
        "rel_same_n_nge2": fits_of([r for r in same_n if r["nA"] >= 2], "R_rel"),
        "rel_1s1s": fits_of([r for r in rows if r["nA"] == 1 and r["nB"] == 1], "R_rel"),
    }
    # dense 1s-1s exponent-ratio scan (scale-free): a=1, b=t; ell_A=1, ell_B=1/t.
    ts = np.geomspace(1.0, 8.0, 36)
    scan = []
    for t in ts:
        Rc, Smax, _, Rrel = crossing_from_above(lambda R: two_center_s_overlap(1, 1.0, 1, float(t), R))
        scan.append({"t": float(t), "S0": float(S_1c(1.0, float(t))), "R_star": Rc, "S_max": Smax, "R_rel": Rrel})
    have_t = [s for s in scan if s["R_star"] is not None]
    fits["scan_abs"] = fit_laws(np.ones(len(have_t)), 1.0 / np.array([s["t"] for s in have_t]),
                                np.array([s["R_star"] for s in have_t])) if len(have_t) >= 3 else {}
    fits["scan_rel"] = fit_laws(np.ones(len(scan)), 1.0 / np.array([s["t"] for s in scan]),
                                np.array([s["R_rel"] for s in scan]))
    t_c = None
    for s in scan:
        if s["R_star"] is None:
            t_c = s["t"]; break

    log(f"  grid: {len(rows)} (ZA,ZB,nA,nB) entries, {len(have)} with an absolute 1/sqrt2 crossing")
    show("fits ABSOLUTE crossing (|S|=1/sqrt2), all crossing entries:", fits["abs_all"])
    show("fits RELATIVE crossing (|S|=|S|max/sqrt2), same-n entries:", fits["rel_same_n"])
    show("fits RELATIVE, same-n with n>=2:", fits["rel_same_n_nge2"])
    show("fits RELATIVE, 1s-1s entries:", fits["rel_1s1s"])
    show("fits RELATIVE, all 64 entries (mixed n included):", fits["rel_all"])
    show(f"dense 1s-1s ratio scan a=1,b=t (ABSOLUTE; crossing lost at t_c~{t_c}):", fits["scan_abs"])
    show("dense 1s-1s ratio scan a=1,b=t (RELATIVE, all t):", fits["scan_rel"])
    log("  table (same-n entries):")
    log("   ZA ZB  n   ellA  ellB   S(0)    R*abs   R*abs/(lA+lB)   R_rel   R_rel/(lA+lB)  R_rel/max")
    for r in same_n:
        ra = "  none" if r["R_star"] is None else f"{r['R_star']:6.3f}"
        rr = "  none" if r["R_rel"] is None else f"{r['R_rel']:6.3f}"
        ra2 = "     -" if r["R_star"] is None else f"{r['R_star']/(r['ellA']+r['ellB']):6.3f}"
        rr2 = "     -" if r["R_rel"] is None else f"{r['R_rel']/(r['ellA']+r['ellB']):6.3f}"
        rr3 = "     -" if r["R_rel"] is None else f"{r['R_rel']/max(r['ellA'],r['ellB']):6.3f}"
        log(f"   {r['ZA']:2d} {r['ZB']:2d} {r['nA']:2d}  {r['ellA']:5.2f} {r['ellB']:5.2f}  {r['S_max']:.3f}  "
            f"{ra}    {ra2}       {rr}    {rr2}      {rr3}")
    log("  ratio scan (a=1, b=t):   t     S(0)    R*abs   R_rel   R_rel*t/(1+t) [=c if additive]  R_rel [=c if max-law]")
    for s in scan[::5]:
        ra = "  none" if s["R_star"] is None else f"{s['R_star']:6.3f}"
        log(f"     {s['t']:5.2f}  {s['S0']:.3f}  {ra}  {s['R_rel']:6.3f}   {s['R_rel']*s['t']/(1+s['t']):6.3f}"
            f"                     {s['R_rel']:6.3f}")
    return {"rows": rows, "fits": fits, "ratio_scan": scan, "t_c_abs_crossing_lost": t_c}


# =============================================================== RUNG 0 (b)
def _eta_matrix(c2: float, b: float, N: int = 50) -> np.ndarray:
    """H_eta = -l(l+1) + c^2 eta^2 + b eta in the orthonormal Legendre basis (m=0),
    identical construction to geovac.molecular_sturmian._angular_sep_const."""
    r_vals = np.arange(0, N, dtype=float)
    norms = np.array([2.0 / (2 * r + 1) for r in r_vals])
    nu = np.zeros((N, N))
    for i in range(N - 1):
        r = r_vals[i]
        val = (r + 1) * np.sqrt(norms[i + 1]) / ((2 * r + 1) * np.sqrt(norms[i]))
        nu[i + 1, i] = val
        nu[i, i + 1] = val
    H = np.diag(-r_vals * (r_vals + 1)) + c2 * (nu @ nu) + b * nu
    return H, nu


def eta_ground(c2: float, b: float, N: int = 50):
    H, nu = _eta_matrix(c2, b, N)
    w, V = np.linalg.eigh(H)
    k = int(np.argmax(w))
    v = V[:, k]
    p = v ** 2 / np.sum(v ** 2)
    nz = p[p > 1e-300]
    return {"A": float(w[k]), "participation_deficit": float(1.0 - p.max()),
            "entropy_bits": float(-(nz * np.log2(nz)).sum()),
            "eta_mean": float(v @ nu @ v), "p_l": p.tolist()}


def rung0b(R_grid: np.ndarray) -> dict:
    log("\n=== RUNG 0(b): HeH2+ eta-equation l-mixing vs R (b = R(Z_B-Z_A)) ===")
    from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice
    from geovac.molecular_sturmian import _angular_sep_const
    out = {"R": R_grid.tolist(), "HeH2+": [], "H2+": [], "max_A_mismatch": 0.0}
    for R in R_grid:
        for key, (ZA, ZB) in (("HeH2+", (2, 1)), ("H2+", (1, 1))):
            try:
                lat = ProlateSpheroidalLattice(float(R), ZA, ZB, N_xi=2500, xi_max=30.0,
                                               radial_method='spectral', n_basis=30)
                E, c2, A_solver = lat.solve()
                b = float(R) * (ZB - ZA)
                g = eta_ground(c2, b)
                A_ref = _angular_sep_const(0, 0, np.sqrt(max(c2, 1e-15)), b=b, n_basis=50)
                out["max_A_mismatch"] = max(out["max_A_mismatch"], abs(g["A"] - A_ref), abs(g["A"] - A_solver))
                rec = {"R": float(R), "E_elec": float(E), "c2": float(c2), "b": b, **g}
            except Exception as ex:      # pragma: no cover
                rec = {"R": float(R), "error": str(ex)}
            rec.pop("p_l", None) if key == "H2+" else None
            out[key].append(rec)
    # yesterday's H2+ curve for cross-check
    try:
        with open(os.path.join(DATA_DIR, "decompactification_R_sweep.json")) as fh:
            prev = json.load(fh)
        out["H2+_yesterday"] = {"R": prev["R_main"], **prev["bonus_eta_lmixing"]}
    except Exception:
        out["H2+_yesterday"] = None
    log("     R     HeH2+: c2    part.def  entropy  <eta>   |  H2+: part.def  entropy")
    for a, h in zip(out["HeH2+"], out["H2+"]):
        if "error" in a or "error" in h:
            continue
        log(f"  {a['R']:6.2f}   {a['c2']:8.3f}   {a['participation_deficit']:.4f}   {a['entropy_bits']:.3f}  "
            f"{a['eta_mean']:+.3f}   |   {h['participation_deficit']:.4f}   {h['entropy_bits']:.3f}")
    log(f"  eta-eigensolver A vs solver / _angular_sep_const: max mismatch {out['max_A_mismatch']:.2e}")
    return out


# =============================================================== RUNG 1/2 integrals
def sto_norm(z: float) -> float:
    return sqrt(z ** 3 / pi)


def S_1c(a: float, b: float) -> float:
    return (2.0 * sqrt(a * b) / (a + b)) ** 3


def V_1c_own(a: float, b: float) -> float:
    """<a|1/r|b>, both 1s STOs on the same center."""
    return 4.0 * (a * b) ** 1.5 / (a + b) ** 2


def V_1c_other(a: float, b: float, R: float) -> float:
    """<a_A|1/r_B|b_A>: Coulomb potential at distance R of the one-center product."""
    c = a + b
    return S_1c(a, b) * (1.0 / R - np.exp(-c * R) * (1.0 / R + 0.5 * c))


def T_1c(a: float, b: float) -> float:
    return 4.0 * (a * b) ** 2.5 / (a + b) ** 3


_GL = {}


def _gl(n: int):
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    return _GL[n]


def ab_integrals(a: float, b: float, R: float, n_xi: int = 400, n_eta: int = 120):
    """(S, <a_A|1/r_A|b_B>, <a_A|1/r_B|b_B>) for 1s STOs, prolate-spheroidal Gauss-Legendre."""
    half = R / 2.0
    rate = half * (a + b)
    xi_max = 1.0 + 40.0 / rate
    xe, we = _gl(n_eta)
    xx, wx = _gl(n_xi)
    xi = 0.5 * (xi_max - 1.0) * xx + 0.5 * (xi_max + 1.0)
    wxi = wx * 0.5 * (xi_max - 1.0)
    XI, ETA = np.meshgrid(xi, xe, indexing="ij")
    W = np.outer(wxi, we)
    r1 = half * (XI + ETA)
    r2 = half * (XI - ETA)
    f = sto_norm(a) * sto_norm(b) * np.exp(-a * r1 - b * r2) * half ** 3 * (XI ** 2 - ETA ** 2) * W
    f *= 2.0 * pi                                   # azimuth
    S = float(f.sum())
    VA = float((f / r1).sum())
    VB = float((f / r2).sum())
    return S, VA, VB


class TwoCenterSTO:
    """Even-tempered 1s-STO basis on two centers (A at origin, B at R zhat).

    One-electron integrals: closed form / prolate GL (exact to ~1e-9).
    ERIs: multipole expansion about A on a graded radial grid (Paper 60 route)."""

    def __init__(self, R: float, ZA: float, ZB: float, zA, zB,
                 nr: int = 3000, nth: int = 300, Lmax: int = 24, rmax: float | None = None):
        self.R, self.ZA, self.ZB = float(R), float(ZA), float(ZB)
        self.zetas = np.array(list(zA) + list(zB), float)
        self.cen = np.array([0] * len(zA) + [1] * len(zB))
        self.idxA = np.where(self.cen == 0)[0]
        self.idxB = np.where(self.cen == 1)[0]
        N = self.N = len(self.zetas)
        # ---- one-electron
        S = np.zeros((N, N)); VA = np.zeros((N, N)); VB = np.zeros((N, N)); T = np.zeros((N, N))
        for i in range(N):
            for j in range(i, N):
                a, b = self.zetas[i], self.zetas[j]
                if self.cen[i] == self.cen[j]:
                    s = S_1c(a, b); own = V_1c_own(a, b); oth = V_1c_other(a, b, R); t = T_1c(a, b)
                    if self.cen[i] == 0:
                        va, vb = own, oth
                    else:
                        va, vb = oth, own
                else:                                   # i on A, j on B (ordering guarantees)
                    s, va, vb = ab_integrals(a, b, R)
                    # kinetic via the STO identity -1/2 lap chi_z = (z/r_c) chi_z - z^2/2 chi_z, symmetrized
                    t = 0.5 * ((b * vb - 0.5 * b * b * s) + (a * va - 0.5 * a * a * s))
                S[i, j] = S[j, i] = s; VA[i, j] = VA[j, i] = va
                VB[i, j] = VB[j, i] = vb; T[i, j] = T[j, i] = t
        self.S, self.VA, self.VB, self.T = S, VA, VB, T
        self.h = T - self.ZA * VA - self.ZB * VB
        # gates on the one-electron layer: cross-center S vs exact Mulliken form
        errS = 0.0
        for i in self.idxA:
            for j in self.idxB:
                errS = max(errS, abs(S[i, j] - two_center_s_overlap(1, self.zetas[i], 1, self.zetas[j], R)))
        self.gate_S_cross_vs_exact = errS
        # ---- ERIs (multipole about A, graded radial grid)
        if rmax is None:
            rmax = max(60.0, 2.0 * R + 30.0)
        t = np.linspace(0.0, 1.0, nr)
        r = rmax * t ** 2 + 1e-4
        u = np.sort(np.cos(np.linspace(0.0, pi, nth)))
        PL = np.array([eval_legendre(L, u) for L in range(Lmax + 1)])          # (L, nth)
        RR, UU = np.meshgrid(r, u, indexing="ij")
        dB = np.sqrt(RR * RR + R * R - 2.0 * R * RR * UU)
        phi = []
        for k in range(N):
            d = RR if self.cen[k] == 0 else dB
            phi.append(sto_norm(self.zetas[k]) * np.exp(-self.zetas[k] * d))     # 3D-normalized chi
        pairs = [(i, j) for i in range(N) for j in range(i, N)]
        self.pair_index = {}
        for p, (i, j) in enumerate(pairs):
            self.pair_index[(i, j)] = p; self.pair_index[(j, i)] = p
        P = len(pairs)
        # multipole moments A^L_ij(r) = 2 pi int phi_i phi_j P_L du  (trapezoid in u, non-uniform)
        du = np.diff(u); wu = np.zeros(nth); wu[:-1] += 0.5 * du; wu[1:] += 0.5 * du
        PLw = PL * wu[None, :]                                                   # (L, nth)
        A = np.zeros((P, Lmax + 1, nr))
        for p, (i, j) in enumerate(pairs):
            A[p] = (2.0 * pi * (phi[i] * phi[j]) @ PLw.T).T                      # (nr,L) -> (L,nr)
        del phi
        Ls = np.arange(Lmax + 1, dtype=float)
        rL = r[None, :] ** Ls[:, None]                                           # (L, nr)
        rmL1 = r[None, :] ** (-(Ls[:, None] + 1.0))
        g = A * (r * r)[None, None, :]
        inner = cumulative_trapezoid(g * rL[None], x=r, axis=-1, initial=0.0) * rmL1[None]
        outer = cumulative_trapezoid((g * rmL1[None])[..., ::-1], x=r[::-1], axis=-1, initial=0.0)[..., ::-1]
        outer = -outer * rL[None]                                                # reversed x -> sign
        Wpot = inner + outer
        del inner, outer, g
        dr = np.diff(r); wr = np.zeros(nr); wr[:-1] += 0.5 * dr; wr[1:] += 0.5 * dr
        G = np.einsum('pLr,qLr,r->pq', A, Wpot, wr * r * r, optimize=True)
        self.gate_G_asym = float(np.abs(G - G.T).max())
        G = 0.5 * (G + G.T)
        eri = np.zeros((N, N, N, N))
        for (i, j), p in self.pair_index.items():
            for (k, l), q in self.pair_index.items():
                eri[i, j, k, l] = G[p, q]
        self.eri = eri
        # gates: grid overlap (L=0 moment) vs exact S ; one-center (aa|aa)=5a/8
        Sg = np.array([_trapz(A[self.pair_index[(i, j)], 0] * r * r, r) for i in range(N) for j in range(N)]).reshape(N, N)
        self.gate_S_grid_vs_exact = float(np.abs(Sg - S).max())
        self.gate_eri_1c = float(max(abs(eri[i, i, i, i] - 5.0 * self.zetas[i] / 8.0) for i in range(N)))
        # grid potential of each A-centered density at nucleus B (L=0 term of W at r=R) vs closed form
        errV = 0.0
        for i in self.idxA:
            p = self.pair_index[(i, i)]
            v_grid = float(np.interp(R, r, Wpot[p, 0]))
            errV = max(errV, abs(v_grid - V_1c_other(self.zetas[i], self.zetas[i], R)))
        self.gate_V_AA_at_B = errV
        del A, Wpot

    # ---------------- orthonormal frame
    def canonical_X(self, thresh: float = 1e-7):
        w, U = np.linalg.eigh(self.S)
        keep = w > thresh
        X = U[:, keep] / np.sqrt(w[keep])[None, :]
        return X, int((~keep).sum()), float(w.min())

    # ---------------- RHF (2 electrons)
    def rhf(self, X: np.ndarray, maxit: int = 300, tol: float = 1e-10):
        h, eri = self.h, self.eri
        Fo = X.T @ h @ X
        _, v = np.linalg.eigh(Fo)
        C = X @ v[:, :1]
        P = 2.0 * C @ C.T
        E_old = 0.0
        for it in range(maxit):
            J = np.einsum('ijkl,kl->ij', eri, P)
            K = np.einsum('ikjl,kl->ij', eri, P)
            F = h + J - 0.5 * K
            E = 0.5 * np.sum(P * (h + F))
            Fo = X.T @ F @ X
            _, v = np.linalg.eigh(Fo)
            C = X @ v[:, :1]
            P_new = 2.0 * C @ C.T
            if abs(E - E_old) < tol and np.abs(P_new - P).max() < 1e-8:
                P = P_new
                break
            P = 0.5 * P + 0.5 * P_new if it > 3 else P_new
            E_old = E
        J = np.einsum('ijkl,kl->ij', eri, P); K = np.einsum('ikjl,kl->ij', eri, P)
        E = 0.5 * np.sum(P * (2 * h + J - 0.5 * K))
        return float(E), C[:, 0], P, it

    # ---------------- singlet FCI (2 electrons)
    def fci(self, X: np.ndarray):
        hm = X.T @ self.h @ X
        em = np.einsum('ip,jq,kr,ls,ijkl->pqrs', X, X, X, X, self.eri, optimize=True)
        M = hm.shape[0]
        I = np.eye(M)
        H = np.kron(hm, I) + np.kron(I, hm) + em.transpose(0, 2, 1, 3).reshape(M * M, M * M)
        H = 0.5 * (H + H.T)
        w, V = np.linalg.eigh(H)
        C = V[:, 0].reshape(M, M)
        asym = float(np.abs(C - C.T).max())
        C = 0.5 * (C + C.T); C /= np.linalg.norm(C)
        gamma_mo = 2.0 * C @ C.T
        occ, U = np.linalg.eigh(gamma_mo)
        order = np.argsort(occ)[::-1]
        occ, U = occ[order], U[:, order]
        NOs = X @ U                                   # raw coefficients, columns
        P_raw = X @ gamma_mo @ X.T
        return float(w[0]), occ, NOs, P_raw, asym


# =============================================================== density analysis
R_FINE = np.linspace(0.0, 60.0, 12001)


def radial_component(zetas: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    """phi(r) = sum_mu c_mu N_mu e^{-zeta_mu r} on R_FINE."""
    return (coeffs * np.array([sto_norm(z) for z in zetas])) @ np.exp(-np.outer(zetas, R_FINE))


def decay_length_routes(zetas: np.ndarray, coeffs: np.ndarray) -> dict:
    """Effective decay exponent of a one-center s-function: route A (window log-slope),
    route B (moment 3/(2<r>)), with the window stated and the local exponent drift."""
    phi = radial_component(zetas, coeffs)
    r = R_FINE
    w = phi * phi * r * r
    norm = _trapz(w, r)
    if norm <= 0:
        return {"zeta_A": np.nan, "zeta_B": np.nan}
    rmean = _trapz(w * r, r) / norm
    zeta_B = 1.5 / rmean
    cum = np.concatenate(([0.0], cumulative_trapezoid(w, r))) / norm
    r50 = float(np.interp(0.50, cum, r)); r95 = float(np.interp(0.95, cum, r))
    mask = (r >= r50) & (r <= r95) & (np.abs(phi) > 1e-300)
    sgn = np.sign(phi[mask][0])
    y = np.log(np.abs(phi[mask])); x = r[mask]
    p = np.polyfit(x, y, 1)
    yhat = np.polyval(p, x)
    r2 = 1.0 - np.sum((y - yhat) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-300)
    dphi = np.gradient(phi, r)
    def zloc(rr):
        k = int(np.searchsorted(r, rr))
        k = min(max(k, 1), len(r) - 2)
        return float(-dphi[k] / phi[k]) if abs(phi[k]) > 1e-300 else np.nan
    return {"zeta_A": float(-p[0]), "zeta_B": float(zeta_B), "fit_R2": float(r2),
            "window": [r50, r95], "zeta_local_r50": zloc(r50), "zeta_local_r95": zloc(r95),
            "zeta_local_2r95": zloc(min(2 * r95, 59.0)), "sign": float(sgn),
            "min_basis_zeta": float(zetas.min()), "max_basis_zeta": float(zetas.max())}


def density_block_zeta(zetas: np.ndarray, Pblock: np.ndarray) -> float:
    """Route C: moment exponent of the same-center block density rho(r)=sum P_mn chi_m chi_n."""
    Nn = np.array([sto_norm(z) for z in zetas])
    E = np.exp(-np.outer(zetas, R_FINE)) * Nn[:, None]           # (n, r)
    rho = np.einsum('mr,mn,nr->r', E, Pblock, E)
    r = R_FINE
    n0 = _trapz(rho * r * r, r)
    if n0 <= 0:
        return np.nan
    return float(1.5 / (_trapz(rho * r ** 3, r) / n0))


def analyze(bas: TwoCenterSTO, NOs: np.ndarray, occ: np.ndarray, P: np.ndarray) -> dict:
    S = bas.S; iA, iB = bas.idxA, bas.idxB
    SAA, SBB, SAB = S[np.ix_(iA, iA)], S[np.ix_(iB, iB)], S[np.ix_(iA, iB)]
    per_no = []
    for k in range(NOs.shape[1]):
        if occ[k] < 1e-6:
            continue
        cA, cB = NOs[iA, k], NOs[iB, k]
        nA2 = float(cA @ SAA @ cA); nB2 = float(cB @ SBB @ cB); x = float(cA @ SAB @ cB)
        cos_k = x / sqrt(max(nA2 * nB2, 1e-300))
        per_no.append({"occ": float(occ[k]), "normA2": nA2, "normB2": nB2, "cos": cos_k})
    ntot = sum(d["occ"] for d in per_no)
    M1 = sum(d["occ"] * abs(d["cos"]) for d in per_no) / ntot
    M0 = abs(per_no[0]["cos"])
    M1_signed = sum(d["occ"] * d["cos"] for d in per_no) / ntot
    # M3: subspace route -- cos of the smallest principal angle between span{phi_k^A} and
    # span{phi_k^B} over NOs with n_k > 0.01 (SVD of the cross overlap in orthonormal frames)
    ks = [k for k in range(NOs.shape[1]) if occ[k] > 0.01]
    CA, CB = NOs[np.ix_(iA, ks)], NOs[np.ix_(iB, ks)]
    def _orth(C, Sblk):
        M = C.T @ Sblk @ C
        w, U = np.linalg.eigh(M)
        keep = w > 1e-10 * w.max()
        return C @ (U[:, keep] / np.sqrt(w[keep])[None, :])
    QA, QB = _orth(CA, SAA), _orth(CB, SBB)
    sig = np.linalg.svd(QA.T @ SAB @ QB, compute_uv=False)
    M3 = float(np.clip(sig.max(), 0.0, 1.0))
    M3_angles = np.degrees(np.arccos(np.clip(sig, 0.0, 1.0))).tolist()
    PAA, PBB, PAB = P[np.ix_(iA, iA)], P[np.ix_(iB, iB)], P[np.ix_(iA, iB)]
    tAA = float(np.sum(PAA * SAA)); tBB = float(np.sum(PBB * SBB)); tAB = float(np.sum(PAB * SAB))
    M2 = tAB / sqrt(max(tAA * tBB, 1e-300))
    # Mulliken gross populations
    qA = tAA + tAB; qB = tBB + tAB
    # decay lengths from the dominant NO's per-center components
    c1 = NOs[:, 0]
    dA = decay_length_routes(bas.zetas[iA], c1[iA])
    dB = decay_length_routes(bas.zetas[iB], c1[iB])
    zC_A = density_block_zeta(bas.zetas[iA], PAA); zC_B = density_block_zeta(bas.zetas[iB], PBB)
    return {"per_NO": per_no, "M0": M0, "M1": M1, "M1_signed": M1_signed, "M2": M2,
            "M3": M3, "M3_principal_angles_deg": M3_angles,
            "popA": qA, "popB": qB, "netA": tAA, "netB": tBB, "overlap_pop": tAB,
            "zetaA": dA, "zetaB": dB, "zetaC_A": zC_A, "zetaC_B": zC_B}


def cos_pred(zA: float, zB: float, R: float) -> float:
    """Rung-0 law at Z_eff: exact 1s(zA)_A - 1s(zB)_B two-center overlap."""
    if not (np.isfinite(zA) and np.isfinite(zB)) or zA <= 0 or zB <= 0:
        return np.nan
    return float(two_center_s_overlap(1, zA, 1, zB, R))


def find_crossing(Rs: np.ndarray, vals: np.ndarray, level: float = SQRT_HALF):
    """Largest R at which vals crosses `level` from above (PCHIP + brentq). None if absent."""
    m = np.isfinite(vals)
    Rs, vals = Rs[m], vals[m]
    if len(Rs) < 3:
        return None
    f = PchipInterpolator(Rs, vals - level)
    idx = [i for i in range(len(Rs) - 1) if (vals[i] - level) > 0 >= (vals[i + 1] - level)]
    if not idx:
        return None
    i = idx[-1]
    return float(brentq(f, Rs[i], Rs[i + 1], xtol=1e-9))


# =============================================================== molecule sweep
def et_set(z0: float, z1: float, n: int) -> list:
    return [float(z) for z in np.geomspace(z0, z1, n)]


BASES = {
    "H":  {"base": et_set(0.5, 2.53125, 5), "double": et_set(0.5, 2.53125, 10),
           "wide": et_set(0.35, 3.62, 8)},
    "He": {"base": et_set(0.9, 4.55625, 5), "double": et_set(0.9, 4.55625, 10),
           "wide": et_set(0.63, 6.52, 8)},
}
MOLECULES = {
    "H2":   {"ZA": 1.0, "ZB": 1.0, "setA": "H", "setB": "H",
             "R": [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0,
                   2.4, 2.8, 3.2, 4.0, 5.0, 6.0],
             "E_exact": {"1.4": -1.17447}, "E_atoms": -1.0},
    "HeH+": {"ZA": 2.0, "ZB": 1.0, "setA": "He", "setB": "H",
             "R": [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.46, 1.8, 2.2, 2.6,
                   3.0, 3.5, 4.0, 5.0],
             "E_exact": {"1.46": -2.9787}, "E_atoms": -2.9037},
}


def sweep_molecule(name: str, which: str, R_list, thresh: float = 1e-7,
                   nr: int = 3000, nth: int = 300) -> dict:
    spec = MOLECULES[name]
    zA = BASES[spec["setA"]][which]; zB = BASES[spec["setB"]][which]
    log(f"\n=== {name} [{which} basis: {len(zA)}+{len(zB)} 1s-STOs]  ZA={spec['ZA']} ZB={spec['ZB']}"
        f"  (ERI grid nr={nr}, nth={nth}) ===")
    log(f"    zetas A: {np.round(zA, 3).tolist()}\n    zetas B: {np.round(zB, 3).tolist()}")
    recs = []
    for R in R_list:
        t0 = time.time()
        bas = TwoCenterSTO(R, spec["ZA"], spec["ZB"], zA, zB, nr=nr, nth=nth)
        X, ndrop, wmin = bas.canonical_X(thresh)
        E_hf, C_hf, P_hf, nit = bas.rhf(X)
        E_ci, occ, NOs, P_ci, asym = bas.fci(X)
        V_nn = spec["ZA"] * spec["ZB"] / R
        a_hf = analyze(bas, C_hf[:, None], np.array([2.0]), P_hf)
        a_ci = analyze(bas, NOs, occ, P_ci)
        rec = {"R": float(R), "V_nn": V_nn, "E_HF": E_hf, "E_FCI": E_ci,
               "E_HF_tot": E_hf + V_nn, "E_FCI_tot": E_ci + V_nn,
               "n_dropped": ndrop, "S_min_eig": wmin, "hf_iters": nit, "fci_C_asym": asym,
               "occ_top4": occ[:4].tolist(),
               "gates": {"S_cross_vs_exact": bas.gate_S_cross_vs_exact,
                         "S_grid_vs_exact": bas.gate_S_grid_vs_exact,
                         "eri_1c_5a8": bas.gate_eri_1c, "G_asym": bas.gate_G_asym,
                         "V_AA_at_B_vs_exact": bas.gate_V_AA_at_B},
               "HF": a_hf, "FCI": a_ci, "secs": time.time() - t0}
        for rung in ("HF", "FCI"):
            a = rec[rung]
            a["cos_pred_B"] = cos_pred(a["zetaA"]["zeta_B"], a["zetaB"]["zeta_B"], R)
            a["cos_pred_A"] = cos_pred(a["zetaA"]["zeta_A"], a["zetaB"]["zeta_A"], R)
            a["cos_pred_C"] = cos_pred(a["zetaC_A"], a["zetaC_B"], R)
        recs.append(rec)
        log(f"  R={R:5.2f}  E_HF={rec['E_HF_tot']:.5f} E_FCI={rec['E_FCI_tot']:.5f}  "
            f"occ={np.round(occ[:2],4).tolist()}  drop={ndrop}  "
            f"| HF: M1={a_hf['M1']:.4f} zA={a_hf['zetaA']['zeta_B']:.3f} zB={a_hf['zetaB']['zeta_B']:.3f} "
            f"pred={a_hf['cos_pred_B']:.4f} "
            f"| FCI: M1={a_ci['M1']:.4f} M2={a_ci['M2']:.4f} zA={a_ci['zetaA']['zeta_B']:.3f} "
            f"zB={a_ci['zetaB']['zeta_B']:.3f} pred={a_ci['cos_pred_B']:.4f}  "
            f"gates S={rec['gates']['S_grid_vs_exact']:.1e} eri={rec['gates']['eri_1c_5a8']:.1e} "
            f"({rec['secs']:.0f}s)")
    # fronts
    Rs = np.array([r["R"] for r in recs])
    fronts = {}
    for rung in ("HF", "FCI"):
        f = {}
        for key in ("M0", "M1", "M2", "M3", "cos_pred_B", "cos_pred_A", "cos_pred_C"):
            f[key] = find_crossing(Rs, np.array([r[rung][key] for r in recs]))
        fronts[rung] = f
    return {"molecule": name, "basis": which, "zetasA": zA, "zetasB": zB, "records": recs,
            "fronts": fronts}


def summarize(res: dict) -> None:
    name, which = res["molecule"], res["basis"]
    log(f"\n  --- {name} [{which}] fronts (R where measure crosses 1/sqrt2) ---")
    for rung in ("HF", "FCI"):
        f = res["fronts"][rung]
        def fmt(x):
            return "  none " if x is None else f"{x:6.3f}"
        log(f"   {rung:3s}: R*(M1)={fmt(f['M1'])}  R*(M0)={fmt(f['M0'])}  R*(M3 subspace)={fmt(f['M3'])}  "
            f"R*(M2 coherence)={fmt(f['M2'])}  "
            f"| predicted R* at Z_eff: routeB={fmt(f['cos_pred_B'])} routeA={fmt(f['cos_pred_A'])} "
            f"routeC={fmt(f['cos_pred_C'])}")
        if f["M1"] and f["cos_pred_B"]:
            log(f"        residual R*(M1)-R*pred(B) = {f['M1']-f['cos_pred_B']:+.3f}  "
                f"({100*(f['M1']-f['cos_pred_B'])/f['cos_pred_B']:+.1f}%)")
        if f["M2"] and f["cos_pred_B"]:
            log(f"        residual R*(M2)-R*pred(B) = {f['M2']-f['cos_pred_B']:+.3f}  "
                f"({100*(f['M2']-f['cos_pred_B'])/f['cos_pred_B']:+.1f}%)")
    log("   R     | HF: zA(B) zA(A) zB(B)  M1    pred   dcos  | FCI: zA(B) zB(B) n_u    M1    pred   dcos   M2    predM2^2")
    for r in res["records"]:
        h, c = r["HF"], r["FCI"]
        log(f"  {r['R']:5.2f}  | {h['zetaA']['zeta_B']:.3f} {h['zetaA']['zeta_A']:.3f} {h['zetaB']['zeta_B']:.3f}  "
            f"{h['M1']:.4f} {h['cos_pred_B']:.4f} {h['M1']-h['cos_pred_B']:+.4f} | "
            f"{c['zetaA']['zeta_B']:.3f} {c['zetaB']['zeta_B']:.3f} {r['occ_top4'][1]:.3f} "
            f"{c['M1']:.4f} {c['cos_pred_B']:.4f} {c['M1']-c['cos_pred_B']:+.4f} {c['M2']:+.4f} "
            f"{c['cos_pred_B']**2:.4f}")


# =============================================================== validation
def validate() -> dict:
    log("\n=== VALIDATION of the Rung-1/2 substrate ===")
    out = {}
    # (1) prolate GL cross-center integrals vs exact Mulliken overlap, a != b
    errs = []
    for (a, b, R) in [(1.0, 1.0, 1.4), (0.5, 2.5, 1.0), (2.0, 0.7, 3.0), (4.5, 0.5, 0.6)]:
        S, VA, VB = ab_integrals(a, b, R)
        errs.append(abs(S - two_center_s_overlap(1, a, 1, b, R)))
    out["ab_overlap_vs_exact_max_err"] = float(max(errs))
    # (2) H atom in the H set: E -> -0.5 (needs at least one zeta near 1; ET set has 1.125)
    bas = TwoCenterSTO(40.0, 1.0, 0.0, BASES["H"]["base"], [1.0])   # B is a dummy, ZB=0
    X, nd, _ = bas.canonical_X(1e-7)
    hm = X.T @ bas.h @ X
    E_H = float(np.linalg.eigvalsh(hm)[0])
    out["H_atom_E_in_H_set"] = E_H
    # (3) He atom FCI in the He set: R large, ZB=0 dummy on B; s-limit ~ -2.879
    bas = TwoCenterSTO(40.0, 2.0, 0.0, BASES["He"]["base"], [1.0])
    X, nd, _ = bas.canonical_X(1e-7)
    # restrict to A functions only (drop the dummy) -- simplest: rebuild with tiny B contribution ignored
    E_He_hf = bas.rhf(X)[0]; E_He_ci = bas.fci(X)[0]
    out["He_atom_HF_in_He_set"] = E_He_hf; out["He_atom_FCI_in_He_set"] = E_He_ci
    out["He_ref"] = {"HF_limit": -2.8617, "FCI_s_limit_approx": -2.8790, "exact": -2.9037}
    # (4) two-center ERI vs exact aabb_value, one point (1s_A a | 1s_B b), a=b=1.3, R=1.5
    try:
        from fractions import Fraction
        from geovac.two_center_eri import aabb_value
        bas = TwoCenterSTO(1.5, 1.0, 1.0, [1.3], [1.3])
        v_grid = bas.eri[0, 0, 1, 1]
        v_ex = aabb_value(Fraction(13, 10), (1, 0, 0), (1, 0, 0), Fraction(13, 10), (1, 0, 0), (1, 0, 0), 1.5, prec=25)
        out["aabb_grid_vs_exact"] = {"grid": float(v_grid), "exact": float(v_ex), "err": float(abs(v_grid - v_ex))}
    except Exception as ex:
        out["aabb_grid_vs_exact"] = {"error": str(ex)}
    # (5) FCI vs the validated Paper-60 driver (shared-scale minimal basis, zeta=1.2, R=1.4)
    try:
        sys.path.insert(0, os.path.join(REPO, "debug"))
        from sturmian_goscinskian_integrals import GoscinskianIntegrals
        from sturmian_h2_ci_1norm import build_raw, loewdin, fci2e
        gi = GoscinskianIntegrals(R=1.4, Lmax=20, nr=2600, nth=160, rmax=60.0)
        h, S, eri = build_raw(gi, 1, 1.2)
        hm, em = loewdin(h, S, eri)
        E_ref = fci2e(hm, em)
        bas = TwoCenterSTO(1.4, 1.0, 1.0, [1.2], [1.2])
        X, _, _ = bas.canonical_X(1e-7)
        E_mine = bas.fci(X)[0]
        out["fci_vs_paper60_driver_min_basis"] = {"paper60": float(E_ref), "mine": float(E_mine),
                                                   "diff": float(abs(E_ref - E_mine))}
    except Exception as ex:
        out["fci_vs_paper60_driver_min_basis"] = {"error": str(ex)}
    for k, v in out.items():
        log(f"  {k}: {v}")
    return out


# =============================================================== main
def run(quick: bool = False, do_double: bool = True) -> dict:
    t_all = time.time()
    out = {"date": "2026-09-06", "type": "DIAGNOSTIC", "conventions": {
        "front": "R where measure crosses 1/sqrt2 from above (principal angle 45 deg, max ||[P_A,P_B]||)",
        "M0": "|cos| between A- and B-centered components of the dominant natural orbital",
        "M1": "occupation-weighted sum_k n_k |cos_k| / sum_k n_k over natural orbitals (compound-matrix object)",
        "M2": "signed coherence Tr(P_AB S_BA)/sqrt(Tr(P_AA S_AA) Tr(P_BB S_BB)); = M1 for one MO",
        "zeta routes": "A: log-slope of dominant-NO per-center component over [r50,r95]; "
                       "B: 3/(2<r>) moment (exact for pure 1s); C: moment of the AA-block density",
        "cos_pred_X": "Rung-0 law at Z_eff: exact 1s(zeta_A^X)-1s(zeta_B^X) two-center overlap at R",
        "basis": "even-tempered 1s STOs per center; double = 10 interleaved over the same span",
        "integrals": "1e: closed form / prolate GL; ERI: multipole about A, graded radial grid (Paper 60 route)",
    }}
    out["validation"] = validate()
    out["rung0a"] = rung0a()
    R0b = np.round(np.geomspace(0.2, 12.0, 45), 5) if not quick else np.array([0.5, 1.0, 2.0, 4.0, 8.0])
    out["rung0b"] = rung0b(R0b)
    out["ladder"] = {}
    for name in ("H2", "HeH+"):
        R_list = MOLECULES[name]["R"] if not quick else [MOLECULES[name]["R"][4]]
        res = sweep_molecule(name, "base", R_list)
        summarize(res)
        out["ladder"][f"{name}|base"] = res
        for which in (("double", "wide") if do_double else ()):
            res2 = sweep_molecule(name, which, R_list)
            summarize(res2)
            out["ladder"][f"{name}|{which}"] = res2
            # basis-enlargement deltas
            log(f"\n  --- {name}: basis control ({which} - base) ---")
            deltas = {}
            for rung in ("HF", "FCI"):
                f1, f2 = res["fronts"][rung], res2["fronts"][rung]
                d = {k: (None if (f1[k] is None or f2[k] is None) else f2[k] - f1[k]) for k in f1}
                deltas[rung] = d
                zl1 = [r[rung]["zetaA"]["zeta_B"] for r in res["records"]]
                zl2 = [r[rung]["zetaA"]["zeta_B"] for r in res2["records"]]
                zr1 = [r[rung]["zetaB"]["zeta_B"] for r in res["records"]]
                zr2 = [r[rung]["zetaB"]["zeta_B"] for r in res2["records"]]
                m1a = [r[rung]["M1"] for r in res["records"]]; m1b = [r[rung]["M1"] for r in res2["records"]]
                ea = [r["E_FCI_tot" if rung == "FCI" else "E_HF_tot"] for r in res["records"]]
                eb = [r["E_FCI_tot" if rung == "FCI" else "E_HF_tot"] for r in res2["records"]]
                deltas[rung]["max_abs_dzetaA_pct"] = float(100 * np.max(np.abs(np.array(zl2) / np.array(zl1) - 1)))
                deltas[rung]["max_abs_dzetaB_pct"] = float(100 * np.nanmax(np.abs(np.array(zr2) / np.array(zr1) - 1)))
                deltas[rung]["max_abs_dM1"] = float(np.max(np.abs(np.array(m1b) - np.array(m1a))))
                deltas[rung]["max_dE_mHa"] = float(1000 * np.max(np.abs(np.array(eb) - np.array(ea))))
                log(f"   {rung}: dR*(M1)={d['M1']}  dR*(M2)={d['M2']}  dR*pred(B)={d['cos_pred_B']}  "
                    f"max|dzetaA|={deltas[rung]['max_abs_dzetaA_pct']:.2f}%  "
                    f"max|dzetaB|={deltas[rung]['max_abs_dzetaB_pct']:.2f}%  "
                    f"max|dM1|={deltas[rung]['max_abs_dM1']:.4f}  max|dE|={deltas[rung]['max_dE_mHa']:.2f} mHa")
            out["ladder"][f"{name}|{which}_deltas"] = deltas
    # ---- ERI-grid resolution control: base basis at 2x angular / 2x radial resolution
    if not quick:
        out["resolution_control"] = {}
        for name in ("H2", "HeH+"):
            R_sub = [R for R in MOLECULES[name]["R"] if R in (0.6, 0.8, 1.0, 1.3, 1.46, 2.0, 4.0, 6.0)]
            hi = sweep_molecule(name, "base", R_sub, nr=6000, nth=600)
            lo = {r["R"]: r for r in out["ladder"][f"{name}|base"]["records"]}
            rows = []
            for r in hi["records"]:
                l = lo[r["R"]]
                rows.append({"R": r["R"],
                             "dE_FCI_mHa": 1000 * (r["E_FCI_tot"] - l["E_FCI_tot"]),
                             "dM1_HF": r["HF"]["M1"] - l["HF"]["M1"], "dM1_FCI": r["FCI"]["M1"] - l["FCI"]["M1"],
                             "dM2_FCI": r["FCI"]["M2"] - l["FCI"]["M2"],
                             "dzetaA_FCI_pct": 100 * (r["FCI"]["zetaA"]["zeta_B"] / l["FCI"]["zetaA"]["zeta_B"] - 1),
                             "dzetaB_FCI_pct": 100 * (r["FCI"]["zetaB"]["zeta_B"] / l["FCI"]["zetaB"]["zeta_B"] - 1),
                             "gates_hi": r["gates"], "gates_lo": l["gates"]})
            out["resolution_control"][name] = {"rows": rows, "fronts_hi": hi["fronts"],
                                               "fronts_lo": out["ladder"][f"{name}|base"]["fronts"]}
            log(f"\n  --- {name}: ERI-grid resolution control (nr 3000->6000, nth 300->600), base basis ---")
            for w in rows:
                log(f"   R={w['R']:4.2f}  dE_FCI={w['dE_FCI_mHa']:+.3f} mHa  dM1(HF)={w['dM1_HF']:+.5f}  "
                    f"dM1(FCI)={w['dM1_FCI']:+.5f}  dM2(FCI)={w['dM2_FCI']:+.5f}  dzetaA={w['dzetaA_FCI_pct']:+.3f}%  "
                    f"dzetaB={w['dzetaB_FCI_pct']:+.3f}%  gate eri1c lo/hi={w['gates_lo']['eri_1c_5a8']:.1e}/"
                    f"{w['gates_hi']['eri_1c_5a8']:.1e}  S_grid lo/hi={w['gates_lo']['S_grid_vs_exact']:.1e}/"
                    f"{w['gates_hi']['S_grid_vs_exact']:.1e}")
            fh, fl = hi["fronts"], out["ladder"][f"{name}|base"]["fronts"]
            log(f"   fronts on the sub-grid, hi-res: HF M1={fh['HF']['M1']}  FCI M1={fh['FCI']['M1']}  "
                f"FCI M2={fh['FCI']['M2']}  pred(B) HF={fh['HF']['cos_pred_B']} FCI={fh['FCI']['cos_pred_B']}")
        # ---- multipole-truncation control: L_max 24 -> 40 (the large-R (aa|aa) gate is L_max-limited)
        out["lmax_control"] = []
        for name, Rs in (("H2", (1.3, 4.0, 6.0)), ("HeH+", (0.8, 4.0))):
            sp = MOLECULES[name]; zA = BASES[sp["setA"]]["base"]; zB = BASES[sp["setB"]]["base"]
            for R in Rs:
                got = {}
                for Lmax in (24, 40):
                    bas = TwoCenterSTO(R, sp["ZA"], sp["ZB"], zA, zB, Lmax=Lmax)
                    X, _, _ = bas.canonical_X(1e-7)
                    E_hf, C, P, _ = bas.rhf(X); E_ci, occ, NOs, Pci, _ = bas.fci(X)
                    a_hf = analyze(bas, C[:, None], np.array([2.0]), P); a_ci = analyze(bas, NOs, occ, Pci)
                    got[Lmax] = {"E_FCI_tot": E_ci + sp["ZA"] * sp["ZB"] / R, "M1_HF": a_hf["M1"], "M1_FCI": a_ci["M1"],
                                 "M2_FCI": a_ci["M2"], "zetaA_FCI": a_ci["zetaA"]["zeta_B"],
                                 "zetaB_FCI": a_ci["zetaB"]["zeta_B"], "gate_eri_1c": bas.gate_eri_1c}
                a, b = got[24], got[40]
                row = {"molecule": name, "R": R, "L24": a, "L40": b,
                       "dE_FCI_mHa": 1000 * (b["E_FCI_tot"] - a["E_FCI_tot"]),
                       "dM1_HF": b["M1_HF"] - a["M1_HF"], "dM1_FCI": b["M1_FCI"] - a["M1_FCI"],
                       "dM2_FCI": b["M2_FCI"] - a["M2_FCI"],
                       "dzetaA_pct": 100 * (b["zetaA_FCI"] / a["zetaA_FCI"] - 1),
                       "dzetaB_pct": 100 * (b["zetaB_FCI"] / a["zetaB_FCI"] - 1)}
                out["lmax_control"].append(row)
                log(f"  L_max control {name} R={R}: dE_FCI={row['dE_FCI_mHa']:+.3f} mHa  dM1(HF)={row['dM1_HF']:+.5f}  "
                    f"dM1(FCI)={row['dM1_FCI']:+.5f}  dM2={row['dM2_FCI']:+.5f}  dzetaA={row['dzetaA_pct']:+.3f}%  "
                    f"dzetaB={row['dzetaB_pct']:+.3f}%  gate(aa|aa) {a['gate_eri_1c']:.1e} -> {b['gate_eri_1c']:.1e}")
    out["total_secs"] = time.time() - t_all
    os.makedirs(DATA_DIR, exist_ok=True)
    dst = os.path.join(DATA_DIR, "decompactification_correlation_ladder" + ("_quick" if quick else "") + ".json")

    def _clean(o):
        if isinstance(o, dict):
            return {str(k): _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_clean(v) for v in o]
        if isinstance(o, np.ndarray):
            return _clean(o.tolist())
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        return o
    with open(dst, "w") as fh:
        json.dump(_clean(out), fh, indent=1)
    log(f"\nwrote {dst}  ({out['total_secs']:.0f}s)")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-double", action="store_true")
    args = ap.parse_args()
    run(quick=args.quick, do_double=not args.no_double)
