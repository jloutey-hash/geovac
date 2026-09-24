r"""Phase 1 of the LiH "marriage" build (2026-09-23; plan debug/lih_marriage_build_plan.md, memo
debug/sprint_lih_marriage_memo.md) -- GATE G-LEAF.

The only genuinely 4-body-connected term of the explicit-r12 energy,

    I = INT rho1(r1) rho2(r2) rho3(r3) rho4(r4) f(r12) f(r34) / r13  d^12 r ,   f = exp(-GAM r),

evaluated with a REAL non-isotropic pair density from the Phase-0 CI reference on the leaf
electrons (2, 4) by the RI-free reduction and by brute-force Monte Carlo.

Configurations (bridge = electrons 1, 3 sharing the Coulomb; leaves = 2, 4, each f-attached to one
bridge electron):
  C0  control = the sigma gate's model (debug/lih_r12_4body_integral.py): bridge |1s_B|^2 (alpha=1.0),
      leaves |1s_A|^2 (beta=2.70).  Isotropic leaf -> a closed-form dressing exists (below).
  C1  THE G-leaf test: leaves rho_(NO0,NO1) (signed core x bond transition density, p-like about A),
      bridge rho_(NO1,NO1).
  C2  leaves rho_(NO1,NO1) (positive bond density), bridge rho_(NO0,NO0) (core).  Optional.

THE REDUCTION (general form, no isotropic shortcut):  D1(r1) = rho1(r1) Psi2(r1), Psi2(r1) =
INT f(r12) rho2(r2) d^3r2 (the axially-averaged f-kernel dressing, gVee.dress(Kf, .)), D3 likewise;
then I = INT D1 D3 / r13 by the EXACT prolate-Neumann operator (Phase 0b, exact=True).  Two versions:
  (i)  CONVERGED -- Psi evaluated by a direct high-accuracy 3-D quadrature (target-centred spherical
       rule for nodes near the leaf, leaf-centred rule far away; two rule refinements R1 < R2 must
       agree to <= 1e-4 relative), bridge on the 72x44 grid; plus a scratch-grid Neumann check.
  (ii) PRODUCTION -- Psi = gVee.dress(Kf, rho) on the 72x44 kernels.py grid (28-point phi kernel).
The gate is judged on (i); (i) - (ii) is the coarse-grid f-dressing error (a finding for Phase 3).

BRUTE FORCE: 12-D importance-sampled MC.  Every electron is drawn from a positive mixture of 1s
STO densities (sampled analytically) and the integrand weighted by rho/w.  Because the C1 leaf is
signed with INT rho_01 = 0 (orthonormal NOs), the leaf factor carries an exact-mean CONTROL VARIATE:
   Phi(1;2) = (rho_l/w_l)(2) [f(r12) - C(r1, r2)] + M(r1),
   C = Taylor polynomial of f(|r1 - r2|) in u = r2 - c to 2nd order (c = leaf centre; the
       quadrupole term switched off smoothly at s = |r1 - c| -> 0 by chi(s) = s^2/(s^2 + s0^2)),
   M(r1) = INT rho_l C d^3r2 = f M0 + gam f (s^.d) + chi (f/2)[gam^2 s^.Q.s^ - (gam/s)(Tr Q - s^.Q.s^)]
with the moments M0, d, Q of the leaf about c computed on the grid (smooth integrands, exact to
1e-10) and cross-checked by an independent spherical quadrature.  E[Phi(1;2) | r1] = Psi(r1) exactly.
Batch means (heavy 1/r13 tail) + a 3-way split scatter; the same samples also give the plain
(no-control) estimator for C0 (positive leaf) as a check.

Run from root (Git Bash):
    python debug/lih_marriage_phase1.py > debug/data/lih_marriage_phase1.log 2>&1; echo $? > debug/data/lih_marriage_phase1.exit
Options: --configs C0,C1,C2  --mc-minutes 25  --workers 12  --skip-mc  --no-scratch  --quick (tiny rules; smoke test)
Writes debug/data/lih_marriage_phase1.npz (frozen values + MC bars + grid fields for tests/test_lih_r12ci_bridge.py)
and debug/data/lih_marriage_phase1_cache.npz (dressings; reused on rerun).
Nothing in geovac/lih_r12ci/ or debug/lih_vmc.py is modified.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

# Pin BLAS to one thread per process BEFORE numpy imports.  With the spawn Pool (workers re-import this
# module), 12 workers x multi-threaded BLAS oversubscribes 16 cores and hard-crashes on Windows with no
# Python traceback -- observed 2026-09-23 in the C1 R1 dressing (the coarser RC dressing stayed below
# BLAS's multi-threading size threshold and survived; R1's larger matmuls crossed it).  Embarrassingly
# parallel across processes, so single-threaded BLAS per worker is also the fast choice (no slowdown).
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.polynomial import polyval
from scipy.special import eval_legendre

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
REF_EXACT = os.path.join(DATA, "lih_marriage_phase0_ref_exact.npz")
OUT = os.path.join(DATA, "lih_marriage_phase1.npz")
CACHE = os.path.join(DATA, "lih_marriage_phase1_cache.npz")

GAM = 0.50                 # geminal f = exp(-GAM r)   (the exp geminal; its kernels exist)
BETA = 2.70                # sigma-gate leaf 1s_A exponent
ALPHA = 1.00               # sigma-gate bridge 1s_B exponent
LMAX = 34                  # production Neumann L (hVee.neumann_potential default)
LMAX_HI = 50               # convergence report
GATE_REL = 1e-3            # the G-leaf gate
CONV_REL = 1e-4            # two refinements of the direct dressing must agree to this
SIGMA_GATE_MC_BAR = 1.5e-4  # the sigma gate's reduced-vs-MC agreement (CHANGELOG v5.15.10), MC-scatter-limited
SIGMA_GATE_REDUCED_LIMIT = 2.85625e-2
SIGMA_GATE_BRUTE = 2.85582e-2

ZHAT = np.array([0.0, 0.0, 1.0])
YHAT = np.array([0.0, 1.0, 0.0])
PART_W = 1.5               # width (bohr) of the quartic partition bump at B, chi_B = exp(-(rB/W)^4)
MC_MINUTES = {'C0': 3.0, 'C1': 22.0, 'C2': 8.0}   # per-config MC wall budget at --mc-minutes 22 (scaled)

# direct-quadrature rules: R1 (coarse) and R2 (fine); the gate uses R2, R1 is the convergence check.
#   n_gl points per panel; x panels graded toward the pole (eps_x, ratio_x) AND capped at width x_cap
#   (the H lobe seen from a target near A is an analytic-but-sharp feature at generic angles: its
#   trapezoid/GL convergence is set by a narrow strip of analyticity, so the angular rules must be
#   dense everywhere, not only at the poles); s panels graded toward the pole distance, capped at s_cap;
#   n_phi trapezoid points on [0, 2pi) (tc) / n_phi_ac (ac); ac r panels graded toward 0 and R, capped r_cap.
TC_S_BLOCK = 200           # ~thousand quadrature points per density-evaluation block in dress_tc (memory cap; 400->200 2026-09-23 to lower peak per-worker RSS)
RULES: Dict[str, dict] = {
    'R1': dict(n_gl=7, ratio_x=4.0, eps_x=1e-8, x_cap=0.25, n_phi=32, ratio_s=2.0, eps_s=0.01, s_cap=0.8,
               ratio_r=2.0, eps_r=0.01, r_cap=0.8, eps_x_ac=1e-6, n_phi_ac=32, s_extra=20.0, r_max=24.0),
    'R2': dict(n_gl=9, ratio_x=3.0, eps_x=1e-9, x_cap=0.25, n_phi=40, ratio_s=1.7, eps_s=0.005, s_cap=0.7,
               ratio_r=1.7, eps_r=0.005, r_cap=0.7, eps_x_ac=1e-7, n_phi_ac=48, s_extra=20.0, r_max=24.0),
    # RC (the retune, 2026-09-23): the GATE uses R1 (few x 1e-7 on the C0 closed form, 100x inside the
    #   1e-3 gate) and RC is the CHEAPER convergence lower-bracket -- NOT R2, which for the 12-primitive
    #   NO density projected to ~40x R1 (~5.6 h full run).  RC keeps R1's cusp resolution (n_gl=7,
    #   eps_x=1e-8) and only thins the azimuthal / far-panel counts, so RC-vs-R1 stays << the gate.
    'RC': dict(n_gl=7, ratio_x=4.5, eps_x=1e-8, x_cap=0.30, n_phi=24, ratio_s=2.2, eps_s=0.012, s_cap=1.0,
               ratio_r=2.2, eps_r=0.012, r_cap=1.0, eps_x_ac=1e-6, n_phi_ac=24, s_extra=20.0, r_max=24.0),
    'Q': dict(n_gl=4, ratio_x=6.0, eps_x=1e-5, x_cap=0.5, n_phi=8, ratio_s=3.0, eps_s=0.05, s_cap=3.0,
              ratio_r=3.0, eps_r=0.05, r_cap=3.0, eps_x_ac=1e-4, n_phi_ac=8, s_extra=20.0, r_max=25.0),
}


def _t(t0: float) -> str:
    return f"[{time.time() - t0:7.1f}s]"


def rel(x: float, y: float) -> float:
    return abs(x - y) / max(abs(y), 1e-300)


# =========================================================================== #
# closed forms: dressings of a normalised 1s(zeta) density by exp(-gam r) and exp(-gam r)/r
# =========================================================================== #
def psi_f_iso(s: np.ndarray, zeta: float, gam: float) -> np.ndarray:
    """INT (zeta^3/pi) e^{-2 zeta r} e^{-gam |s - r|} d^3r, s = distance from the 1s centre (closed form).
    Angular average of e^{-gam d} over the sphere |r| = r2:  [(1+gam|s-r2|)e^{-gam|s-r2|}
    - (1+gam(s+r2))e^{-gam(s+r2)}] / (2 s r2 gam^2); the radial integral is elementary."""
    s = np.maximum(np.asarray(s, dtype=float), 1e-14)
    A = 2.0 * zeta
    al = A - gam
    be = A + gam
    if al <= 0.0:
        raise ValueError("closed form needs 2 zeta > gam")
    eas = np.exp(-al * s)
    egs = np.exp(-gam * s)
    eAs = np.exp(-A * s)
    I1 = egs * ((1.0 + gam * s) * (1.0 - (1.0 + al * s) * eas) / al ** 2
                - gam * (2.0 - (2.0 + 2.0 * al * s + al ** 2 * s ** 2) * eas) / al ** 3)
    I2 = eAs * ((1.0 - gam * s) * (s / be + 1.0 / be ** 2)
                + gam * (s ** 2 / be + 2.0 * s / be ** 2 + 2.0 / be ** 3))
    I3 = egs * ((1.0 + gam * s) / be ** 2 + 2.0 * gam / be ** 3)
    return (2.0 * zeta ** 3 / (gam ** 2 * s)) * (I1 + I2 - I3)


def psi_yuk_iso(s: np.ndarray, zeta: float, gam: float) -> np.ndarray:
    """INT (zeta^3/pi) e^{-2 zeta r} e^{-gam |s - r|} / |s - r| d^3r (closed form).  V(0) = 4 zeta^3/(2 zeta+gam)^2."""
    s = np.maximum(np.asarray(s, dtype=float), 1e-14)
    A = 2.0 * zeta
    al = A - gam
    be = A + gam
    if al <= 0.0:
        raise ValueError("closed form needs 2 zeta > gam")
    T1 = np.exp(-gam * s) * (1.0 - (1.0 + al * s) * np.exp(-al * s)) / al ** 2
    T2 = np.exp(-A * s) * (1.0 + be * s) / be ** 2
    T3 = np.exp(-gam * s) / be ** 2
    return (2.0 * zeta ** 3 / (gam * s)) * (T1 + T2 - T3)


# =========================================================================== #
# orbital / density evaluation at arbitrary Cartesian points (value only; the prim dicts of the
# Phase-0 artifact = lih_vmc.extract_prim: phi = N xi^j e^{-alpha xi} g(eta), sigma only)
# =========================================================================== #
def prim_values(prims: Sequence[dict], pts: np.ndarray, R: float) -> np.ndarray:
    """(M, N) values of the M primitives at the N Cartesian points."""
    x = pts[:, 0]
    y = pts[:, 1]
    z = pts[:, 2]
    rA = np.sqrt(x * x + y * y + (z + 0.5 * R) ** 2)
    rB = np.sqrt(x * x + y * y + (z - 0.5 * R) ** 2)
    xi = (rA + rB) / R
    eta = (rA - rB) / R
    out = np.empty((len(prims), pts.shape[0]))
    exp_cache: Dict[float, np.ndarray] = {}
    pow_cache: Dict[int, np.ndarray] = {}
    for i, p in enumerate(prims):
        if int(p['mu']) != 0:
            raise ValueError("sigma orbitals only")
        al = float(p['alpha'])
        j = int(p['j'])
        if al not in exp_cache:
            exp_cache[al] = np.exp(-al * xi)
        if j not in pow_cache:
            pow_cache[j] = xi ** j if j > 0 else np.ones_like(xi)
        out[i] = float(p['N']) * pow_cache[j] * exp_cache[al] * polyval(eta, np.asarray(p['eta'], dtype=float))
    return out


def eval_density(spec: dict, pts: np.ndarray, ctx: dict) -> np.ndarray:
    """Density values at pts (N, 3).  spec kinds: 'iso' (1s STO density), 'no_pair' (NO_p NO_q);
    optional 'part' in {'A', 'B'}: multiply by the smooth partition chi_A = rB^2/(rA^2+rB^2), chi_B = 1 - chi_A."""
    R = ctx['R']
    kind = spec['kind']
    if kind == 'iso':
        c = np.asarray(spec['center'], dtype=float)
        zeta = float(spec['zeta'])
        r = np.sqrt(np.sum((pts - c) ** 2, axis=1))
        rho = (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * r)
    elif kind == 'no_pair':
        pv = prim_values(ctx['prims'], pts, R)
        T_no = ctx['T_no']
        p, q = int(spec['p']), int(spec['q'])
        rho = (T_no[:, p] @ pv) * (T_no[:, q] @ pv)
    else:
        raise ValueError(kind)
    part = spec.get('part')
    if part:
        # analytic partition of unity localised at B: chi_B = exp(-(rB/W)^4) (entire in rB^2), chi_A = 1 - chi_B.
        # The B-part is concentrated within ~W of B (its own pole) and carries only e^{-(R/W)^4} ~ 1e-7 of the
        # A-cusp; in the A-part the B-kink is multiplied by ~rB^4, i.e. smoothed to |r|^5 (C^4), while the
        # A-cusp is carried intact.
        x = pts[:, 0]
        y = pts[:, 1]
        z = pts[:, 2]
        rB2 = x * x + y * y + (z - 0.5 * R) ** 2
        chiB = np.exp(-(rB2 / PART_W ** 2) ** 2)
        rho = rho * ((1.0 - chiB) if part == 'A' else chiB)
    return rho


def kernel_fn(kernel: Tuple[str, float]) -> Callable[[np.ndarray], np.ndarray]:
    kind, gam = kernel
    if kind == 'f':
        return lambda s: np.exp(-gam * s)
    if kind == 'yuk':
        return lambda s: np.exp(-gam * s) / np.maximum(s, 1e-300)
    raise ValueError(kind)


# =========================================================================== #
# composite Gauss rules graded toward foci
# =========================================================================== #
def graded_bounds(lo: float, hi: float, foci: Sequence[float], eps: float, ratio: float,
                  max_width: Optional[float] = None) -> np.ndarray:
    """Panel boundaries on [lo, hi], geometrically graded (smallest width eps, ratio) toward every focus;
    with max_width, every panel wider than that is subdivided uniformly."""
    b = {float(lo), float(hi)}
    for fc in foci:
        fc = float(fc)
        if lo <= fc <= hi:
            b.add(fc)
        w = eps
        while True:
            lo_pt = fc - w
            hi_pt = fc + w
            added = False
            if lo_pt > lo + 0.5 * eps:
                b.add(lo_pt)
                added = True
            if hi_pt < hi - 0.5 * eps:
                b.add(hi_pt)
                added = True
            if not added:
                break
            w *= ratio
    bounds = np.array(sorted(b))
    if max_width is not None:
        out = [bounds[0]]
        for b0, b1 in zip(bounds[:-1], bounds[1:]):
            n = int(math.ceil((b1 - b0) / max_width))
            out.extend(list(np.linspace(b0, b1, n + 1)[1:]))
        bounds = np.array(out)
    return bounds


def gl_on_panels(bounds: np.ndarray, n: int) -> Tuple[np.ndarray, np.ndarray]:
    xg, wg = leggauss(n)
    xs, ws = [], []
    for b0, b1 in zip(bounds[:-1], bounds[1:]):
        xs.append(0.5 * (b1 - b0) * (xg + 1.0) + b0)
        ws.append(0.5 * (b1 - b0) * wg)
    return np.concatenate(xs), np.concatenate(ws)


def phi_rule(n_phi: int) -> Tuple[np.ndarray, np.ndarray]:
    """Trapezoid on [0, 2 pi) for an integrand even in phi: k = 0..n/2 with weights (1,2,..,2,1) 2pi/n."""
    k = np.arange(n_phi // 2 + 1)
    phis = 2.0 * np.pi * k / n_phi
    w = np.full(k.size, 2.0 * (2.0 * np.pi / n_phi))
    w[0] = 2.0 * np.pi / n_phi
    w[-1] = 2.0 * np.pi / n_phi
    return phis, w


# =========================================================================== #
# direct 3-D dressing quadrature
#   tc: target-centred spherical rule, polar axis toward `pole` (the cusp-carrying nucleus); the
#       kernel cusp at the target sits at s = 0 (times s^2), the nuclear cusp at the pole x = 1 and
#       s = |pole - t| -> both on coordinate lines, graded panels.
#   ac: centre-anchored spherical rule about `center` with the pole toward `other` (the second
#       nucleus at r = |other - center|, x = 1): the density is evaluated once on the (r, x) grid and
#       only the smooth kernel K(|t - r'|) is target-dependent.  Valid when the density at the target
#       is negligible (the kernel cusp then carries no weight).
# =========================================================================== #
def dress_tc(spec: dict, kernel: Tuple[str, float], targets: np.ndarray, rule: dict, ctx: dict,
             pole: np.ndarray) -> np.ndarray:
    K = kernel_fn(kernel)
    xb = graded_bounds(-1.0, 1.0, [1.0], rule['eps_x'], rule['ratio_x'], rule['x_cap'])
    xs, wx = gl_on_panels(xb, rule['n_gl'])
    phis, wphi = phi_rule(rule['n_phi'])
    sinth = np.sqrt(np.maximum(1.0 - xs ** 2, 0.0))
    out = np.empty(targets.shape[0])
    for it in range(targets.shape[0]):
        t = targets[it]
        v = pole - t
        sC = float(np.linalg.norm(v))
        nhat = v / sC if sC > 1e-12 else ZHAT.copy()
        e2 = YHAT
        e1 = np.cross(e2, nhat)
        e1 /= np.linalg.norm(e1)
        sb = graded_bounds(0.0, sC + rule['s_extra'], [sC], rule['eps_s'], rule['ratio_s'], rule['s_cap'])
        ss, ws = gl_on_panels(sb, rule['n_gl'])
        # directions (nx, nphi, 3)
        dirs = (xs[:, None, None] * nhat[None, None, :]
                + sinth[:, None, None] * (np.cos(phis)[None, :, None] * e1[None, None, :]
                                          + np.sin(phis)[None, :, None] * e2[None, None, :]))
        wang = wx[:, None] * wphi[None, :]                                     # (nx, nphi)
        wrad = ws * ss ** 2 * K(ss)                                            # (ns,)
        acc = 0.0
        blk = max(1, TC_S_BLOCK // max(1, xs.size * phis.size // 1000))        # ~TC_S_BLOCK k points per block
        for i0 in range(0, ss.size, blk):
            sblk = ss[i0:i0 + blk]
            pts = t[None, None, None, :] + sblk[:, None, None, None] * dirs[None, :, :, :]
            rho = eval_density(spec, pts.reshape(-1, 3), ctx).reshape(sblk.size, xs.size, phis.size)
            acc += float(np.sum(wrad[i0:i0 + blk][:, None, None] * wang[None, :, :] * rho))
        out[it] = acc
    return out


def dress_ac(spec: dict, kernel: Tuple[str, float], targets: np.ndarray, rule: dict, ctx: dict,
             center: np.ndarray, other: np.ndarray) -> np.ndarray:
    K = kernel_fn(kernel)
    Rd = float(np.linalg.norm(other - center))
    pz = float(np.sign((other - center)[2])) if Rd > 0 else 1.0     # pole direction along +/- z
    rb = graded_bounds(0.0, rule['r_max'], [0.0, Rd], rule['eps_r'], rule['ratio_r'], rule['r_cap'])
    rs, wr = gl_on_panels(rb, rule['n_gl'])
    xb = graded_bounds(-1.0, 1.0, [1.0], rule['eps_x_ac'], rule['ratio_x'], rule['x_cap'])
    xs, wx = gl_on_panels(xb, rule['n_gl'])
    phis, wphi = phi_rule(rule['n_phi_ac'])
    sinth = np.sqrt(np.maximum(1.0 - xs ** 2, 0.0))
    RC = rs[:, None] * sinth[None, :]                       # cylindrical radius of the source point
    ZZ = center[2] + pz * rs[:, None] * xs[None, :]         # z of the source point
    pts = np.stack([RC.ravel(), np.zeros(RC.size), ZZ.ravel()], axis=1)
    rho = eval_density(spec, pts, ctx).reshape(rs.size, xs.size)
    W = (wr * rs ** 2)[:, None] * wx[None, :] * rho         # (nr, nx)
    cphi = np.cos(phis)
    out = np.empty(targets.shape[0])
    for it in range(targets.shape[0]):
        t = targets[it]
        rc_t = float(np.hypot(t[0], t[1]))
        dz2 = (t[2] - ZZ) ** 2
        base = rc_t ** 2 + RC ** 2 + dz2                    # (nr, nx)
        d = np.sqrt(np.maximum(base[:, :, None] - 2.0 * rc_t * RC[:, :, None] * cphi[None, None, :], 0.0))
        out[it] = float(np.sum(W[:, :, None] * wphi[None, None, :] * K(d)))
    return out


def dress_direct_chunk(args: tuple) -> np.ndarray:
    """Worker: (kernel, targets, rule_name, ctx, parts, mode) -> Psi at targets.
    parts = list of dict(spec, center, other, r_near): each part of the (partitioned) density is anchored at
    its own centre; a target within r_near of that centre uses the tc rule with the pole at the centre,
    otherwise the ac rule about the centre (pole toward `other`).  mode: 'auto' | 'tc' | 'ac' (forced)."""
    kernel, targets, rule_name, ctx, parts, mode = args
    rule = RULES[rule_name]
    targets = np.asarray(targets, dtype=float)
    out = np.zeros(targets.shape[0])
    for part in parts:
        c = np.asarray(part['center'], dtype=float)
        o = np.asarray(part['other'], dtype=float)
        dist = np.linalg.norm(targets - c[None, :], axis=1)
        if mode == 'tc':
            near = np.ones(targets.shape[0], dtype=bool)
        elif mode == 'ac':
            near = np.zeros(targets.shape[0], dtype=bool)
        else:
            near = dist < part['r_near']
        if np.any(near):
            out[near] += dress_tc(part['spec'], kernel, targets[near], rule, ctx, c)
        if np.any(~near):
            out[~near] += dress_ac(part['spec'], kernel, targets[~near], rule, ctx, c, o)
    return out


def dress_direct(kernel: Tuple[str, float], targets: np.ndarray, rule_name: str, ctx: dict, parts: list,
                 pool=None, n_chunks: int = 48, mode: str = 'auto') -> np.ndarray:
    idx = np.array_split(np.arange(targets.shape[0]), min(n_chunks, targets.shape[0]))
    jobs = [(kernel, targets[i], rule_name, ctx, parts, mode) for i in idx]
    if pool is None:
        res = [dress_direct_chunk(j) for j in jobs]
    else:
        res = pool.map(dress_direct_chunk, jobs)
    out = np.empty(targets.shape[0])
    for i, r in zip(idx, res):
        out[i] = r
    return out


# =========================================================================== #
# the exact prolate-Neumann bridge, per-l, on any GL (xi, eta) grid
# =========================================================================== #
def make_grid(NXI: int, NETA: int, xi_max: float, a: float) -> dict:
    xg, wxg = leggauss(NXI)
    XI = 1.0 + 0.5 * (xg + 1.0) * (xi_max - 1.0)
    WXI = 0.5 * (xi_max - 1.0) * wxg
    ETA, WETA = leggauss(NETA)
    Xg, Eg = np.meshgrid(XI, ETA, indexing='ij')
    JAC = Xg ** 2 - Eg ** 2
    W2D = np.outer(WXI, WETA)
    RHO_CYL = a * np.sqrt(np.maximum((Xg ** 2 - 1.0) * (1.0 - Eg ** 2), 0.0))
    ZC = a * Xg * Eg
    pts = np.stack([RHO_CYL.ravel(), np.zeros(Xg.size), ZC.ravel()], axis=1)
    return dict(NXI=NXI, NETA=NETA, XI=XI, WXI=WXI, ETA=ETA, WETA=WETA, Xg=Xg, Eg=Eg, JAC=JAC, W2D=W2D,
                geo=(W2D * JAC).ravel(), pts=pts, a=a, xi_max=xi_max, rA=(a * (Xg + Eg)).ravel(),
                rB=(a * (Xg - Eg)).ravel(), RHO_CYL=RHO_CYL.ravel(), ZC=ZC.ravel())


def grid_int_g(G: dict, F: np.ndarray) -> float:
    return float(2.0 * np.pi * G['a'] ** 3 * np.sum(G['geo'] * F.ravel()))


def bridge_per_l(D1: np.ndarray, D3: np.ndarray, G: dict, op, lmax: int, R: float) -> np.ndarray:
    """c_l, l = 0..lmax:  INT D1 V_l[D3],  V_l = (2/R)(2 pi) a^3 (2l+1) radial_l(xi) P_l(eta),
    radial_l = exact ordered integral of the eta-moment g_l = INT (xi^2-eta^2) D3 P_l(eta) deta."""
    NXI, NETA = G['NXI'], G['NETA']
    a = G['a']
    pref = (2.0 / R) * (2.0 * np.pi) * a ** 3
    W3 = G['JAC'] * D3.reshape(NXI, NETA)
    W1 = (G['W2D'] * G['JAC'] * D1.reshape(NXI, NETA)) * (2.0 * np.pi * a ** 3)
    out = np.empty(lmax + 1)
    for l in range(lmax + 1):
        Pl = eval_legendre(l, G['ETA'])
        g_l = (W3 * Pl[None, :]) @ G['WETA']
        radial = op.radial(g_l, l, 0)
        out[l] = pref * (2 * l + 1) * float(np.sum(W1 * np.outer(radial, Pl)))
    return out


def legendre_content(F: np.ndarray, G: dict, lmax: int = 4) -> np.ndarray:
    """a_l = |INT F r_A^l P_l(cos theta_A)| / INT |F| r_A^l  about centre A (Phase-0 definition)."""
    cA = ((1.0 + G['Xg'] * G['Eg']) / (G['Xg'] + G['Eg'])).ravel()
    rA = G['rA']
    return np.array([abs(grid_int_g(G, F * rA ** l * eval_legendre(l, cA))) / grid_int_g(G, np.abs(F) * rA ** l)
                     for l in range(lmax + 1)])


# =========================================================================== #
# Monte Carlo (12-D importance-sampled, exact-mean control variate on the leaf factor)
# =========================================================================== #
def sample_1s(rng: np.random.Generator, n: int, zeta: float, center: np.ndarray) -> np.ndarray:
    r = rng.gamma(3.0, 1.0 / (2.0 * zeta), size=n)
    ct = rng.uniform(-1.0, 1.0, size=n)
    ph = rng.uniform(0.0, 2.0 * np.pi, size=n)
    st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
    return np.stack([r * st * np.cos(ph), r * st * np.sin(ph), r * ct], axis=1) + center[None, :]


def sample_mixture(rng: np.random.Generator, mix: List[Tuple[np.ndarray, float, float]], n: int) -> np.ndarray:
    p = np.array([m[2] for m in mix])
    p = p / p.sum()
    comp = rng.choice(len(mix), size=n, p=p)
    out = np.empty((n, 3))
    for k, (c, zeta, _) in enumerate(mix):
        idx = np.nonzero(comp == k)[0]
        if idx.size:
            out[idx] = sample_1s(rng, idx.size, zeta, np.asarray(c, dtype=float))
    return out


def mixture_density(mix: List[Tuple[np.ndarray, float, float]], pts: np.ndarray) -> np.ndarray:
    p = np.array([m[2] for m in mix])
    p = p / p.sum()
    out = np.zeros(pts.shape[0])
    for pk, (c, zeta, _) in zip(p, mix):
        r = np.sqrt(np.sum((pts - np.asarray(c, dtype=float)[None, :]) ** 2, axis=1))
        out += pk * (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * r)
    return out


def control_terms(r1: np.ndarray, r2: np.ndarray, cv: dict, gam: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Returns (C_dip, C_quad, M_dip, M_quad): the per-sample control polynomial and its exact mean, at
    dipole order and at dipole+quadrupole order (quadrupole switched off near s = 0)."""
    c = np.asarray(cv['center'], dtype=float)
    M0, dz, Qxx, Qzz, s0 = cv['M0'], cv['dz'], cv['Qxx'], cv['Qzz'], cv['s0']
    sv = r1 - c[None, :]
    s = np.sqrt(np.sum(sv ** 2, axis=1))
    s = np.maximum(s, 1e-12)
    sh = sv / s[:, None]
    f = np.exp(-gam * s)
    u = r2 - c[None, :]
    su = np.sum(sh * u, axis=1)
    uu = np.sum(u * u, axis=1)
    C_dip = f * (1.0 + gam * su)
    chi = s ** 2 / (s ** 2 + s0 ** 2)
    quad = 0.5 * f * (gam ** 2 * su ** 2 - (gam / s) * (uu - su ** 2))
    C_quad = C_dip + chi * quad
    sQs = Qxx * (sh[:, 0] ** 2 + sh[:, 1] ** 2) + Qzz * sh[:, 2] ** 2
    trQ = 2.0 * Qxx + Qzz
    M_dip = f * (M0 + gam * sh[:, 2] * dz)
    M_quad = M_dip + chi * 0.5 * f * (gam ** 2 * sQs - (gam / s) * (trQ - sQs))
    return C_dip, C_quad, M_dip, M_quad


def mc_worker(args: tuple) -> np.ndarray:
    """(cfg, seed, n_batches, batch) -> (n_batches, 3) batch means: [plain, dipole-control, dip+quad-control]."""
    cfg, seed, n_batches, batch = args
    rng = np.random.default_rng(seed)
    ctx = cfg['ctx']
    gam = cfg['gam']
    means = np.empty((n_batches, 3))
    for b in range(n_batches):
        r1 = sample_mixture(rng, cfg['mix_b'], batch)
        r3 = sample_mixture(rng, cfg['mix_b'], batch)
        r2 = sample_mixture(rng, cfg['mix_l'], batch)
        r4 = sample_mixture(rng, cfg['mix_l'], batch)
        wb1 = eval_density(cfg['spec_b'], r1, ctx) / mixture_density(cfg['mix_b'], r1)
        wb3 = eval_density(cfg['spec_b'], r3, ctx) / mixture_density(cfg['mix_b'], r3)
        wl2 = eval_density(cfg['spec_l'], r2, ctx) / mixture_density(cfg['mix_l'], r2)
        wl4 = eval_density(cfg['spec_l'], r4, ctx) / mixture_density(cfg['mix_l'], r4)
        inv13 = 1.0 / np.maximum(np.linalg.norm(r1 - r3, axis=1), 1e-12)
        f12 = np.exp(-gam * np.linalg.norm(r1 - r2, axis=1))
        f34 = np.exp(-gam * np.linalg.norm(r3 - r4, axis=1))
        Cd12, Cq12, Md1, Mq1 = control_terms(r1, r2, cfg['cv'], gam)
        Cd34, Cq34, Md3, Mq3 = control_terms(r3, r4, cfg['cv'], gam)
        base = wb1 * wb3 * inv13
        e0 = base * (wl2 * f12) * (wl4 * f34)
        e1 = base * (wl2 * (f12 - Cd12) + Md1) * (wl4 * (f34 - Cd34) + Md3)
        e2 = base * (wl2 * (f12 - Cq12) + Mq1) * (wl4 * (f34 - Cq34) + Mq3)
        means[b] = [e0.mean(), e1.mean(), e2.mean()]
    return means


def mc_run(cfg: dict, n_batches: int, batch: int, pool, n_workers: int, seed0: int) -> np.ndarray:
    per = [n_batches // n_workers + (1 if k < n_batches % n_workers else 0) for k in range(n_workers)]
    jobs = [(cfg, seed0 + 1000 * k + 17, nb, batch) for k, nb in enumerate(per) if nb > 0]
    if pool is None:
        res = [mc_worker(j) for j in jobs]
    else:
        res = pool.map(mc_worker, jobs)
    return np.concatenate(res, axis=0)


def mc_summary(bm: np.ndarray) -> dict:
    """bm (n_batches, n_est): mean, batch-means error, and the 3-way split scatter per estimator."""
    n = bm.shape[0]
    mean = bm.mean(axis=0)
    err = bm.std(axis=0, ddof=1) / np.sqrt(n)
    parts = np.array_split(np.arange(n), 3)
    split = np.array([bm[p].mean(axis=0) for p in parts])
    split_sd = split.std(axis=0, ddof=1) / np.sqrt(3)
    return dict(mean=mean, err=err, split_sd=split_sd, n_batches=n)


# =========================================================================== #
# 6-D semi-brute cross-check: MC over the BRIDGE electrons (1, 3) only, with the leaf electrons
# integrated by the trusted direct-quadrature dressing (grid Psi[rule_hi], few x 1e-7 -- section B)
# barycentrically interpolated onto the sample points.  This isolates the Coulomb-bridge (Neumann)
# step from the dressing step: the dressing IS the direct 3-D quadrature; only 1/r13 is sampled, and
# its analytic Neumann value is (i).  Cheap (interpolation, no per-sample dressing), so N reaches the
# 1e-3 bar; it is the fallback the plan asks for when the 12-D estimator cannot.
# =========================================================================== #
def interp_prolate(Fg2d: np.ndarray, pts: np.ndarray, XI: np.ndarray, ETA: np.ndarray,
                   wXI: np.ndarray, wETA: np.ndarray, R: float, lag: Callable) -> np.ndarray:
    """2-D barycentric interpolation of a grid field Fg2d (NXI x NETA) onto Cartesian points via
    (xi, eta) = ((rA+rB)/R, (rA-rB)/R).  Points beyond the grid are clamped (F ~ 0 there by decay)."""
    a = 0.5 * R
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    rA = np.sqrt(x * x + y * y + (z + a) ** 2)
    rB = np.sqrt(x * x + y * y + (z - a) ** 2)
    xi = np.clip((rA + rB) / R, XI[0], XI[-1])
    eta = np.clip((rA - rB) / R, ETA[0], ETA[-1])
    Lx = lag(XI, wXI, xi)                              # (N, NXI)
    Le = lag(ETA, wETA, eta)                           # (N, NETA)
    return np.einsum('ni,ij,nj->n', Lx, Fg2d, Le, optimize=True)


def mc6_run(D_hi_2d: np.ndarray, mix_b: list, G: dict, wXI: np.ndarray, wETA: np.ndarray, R: float,
            lag: Callable, rng: np.random.Generator, budget_s: float, batch: int = 100_000,
            min_batches: int = 40, max_batches: int = 200_000) -> dict:
    """I = INT D_hi(r1) D_hi(r3) / r13, D_hi = bridge*Psi[hi], by MC over r1, r3 ~ mix_b (batch means)."""
    XI, ETA = G['XI'], G['ETA']
    means: List[float] = []
    t0 = time.time()
    while True:
        r1 = sample_mixture(rng, mix_b, batch)
        r3 = sample_mixture(rng, mix_b, batch)
        D1 = interp_prolate(D_hi_2d, r1, XI, ETA, wXI, wETA, R, lag) / mixture_density(mix_b, r1)
        D3 = interp_prolate(D_hi_2d, r3, XI, ETA, wXI, wETA, R, lag) / mixture_density(mix_b, r3)
        inv13 = 1.0 / np.maximum(np.linalg.norm(r1 - r3, axis=1), 1e-12)
        means.append(float((D1 * D3 * inv13).mean()))
        nb = len(means)
        if nb >= max_batches or (nb >= min_batches and time.time() - t0 > budget_s):
            break
    bm = np.array(means)
    return dict(mean=float(bm.mean()), err=float(bm.std(ddof=1) / np.sqrt(bm.size)),
                n_batches=int(bm.size), N=int(bm.size * batch), wall=time.time() - t0)


# =========================================================================== #
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--configs', default='C0,C1,C2')
    ap.add_argument('--mc-minutes', type=float, default=22.0, help='MC wall budget for C1; C0/C2 scaled per MC_MINUTES')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--skip-mc', action='store_true')
    ap.add_argument('--no-scratch', action='store_true')
    ap.add_argument('--quick', action='store_true', help='tiny rules / tiny MC: smoke test only')
    ap.add_argument('--no-cache', action='store_true')
    ap.add_argument('--merge-out', action='store_true', help='seed the artifact from an existing OUT (C2 follow-on)')
    ap.add_argument('--mc6-minutes', type=float, default=2.5, help='6-D bridge cross-check wall budget per real config')
    args = ap.parse_args(argv)
    T0 = time.time()
    np.set_printoptions(linewidth=150, precision=6, suppress=True)
    configs = [c.strip() for c in args.configs.split(',') if c.strip()]
    rule_lo, rule_hi = ('Q', 'Q') if args.quick else ('RC', 'R1')  # RETUNE: gate = R1, check = RC (see RULES)

    print("=" * 100)
    print("LiH MARRIAGE -- PHASE 1: gate G-leaf, the 4-body bridge I = <rho1 rho2 rho3 rho4 f12 f34 / r13> on REAL leaves")
    print(f"  configs {configs}; rules {rule_lo} (check) / {rule_hi} (gate); GAM={GAM}; gate {GATE_REL:.0e} rel; "
          f"conv {CONV_REL:.0e}; workers {args.workers}; MC budget {args.mc_minutes} min" + ("  [QUICK]" if args.quick else ""))
    print("=" * 100, flush=True)

    # ------------------------------------------------------------------ package + artifact
    import geovac.lih_r12ci.kernels as KG
    from geovac.lih_r12ci.hVee import neumann_potential
    from geovac.lih_r12ci.energy import R as R_ENG, f_dress_iso
    from geovac.lih_r12ci.hT import yukawa_pot_iso
    from geovac.lih_r12ci.neumann_exact import ExactNeumann, barycentric_weights, lagrange_matrix
    assert KG.USE_EXACT_NEUMANN is False
    z = np.load(REF_EXACT, allow_pickle=True)
    R = float(z['R'])
    a = R / 2.0
    assert abs(R - R_ENG) < 1e-12 and abs(a - float(KG.a)) < 1e-12
    A_C = np.array([0.0, 0.0, -a])
    B_C = np.array([0.0, 0.0, +a])
    prims = [dict(N=float(p['N']), j=int(p['j']), alpha=float(p['alpha']), eta=np.asarray(p['eta'], dtype=float),
                  mu=int(p['mu']), s=int(p['s'])) for p in z['prims']]
    T_no = np.asarray(z['T_no'], dtype=float)
    ctx = dict(prims=prims, T_no=T_no, R=R)
    pairs = [tuple(int(x) for x in p) for p in z['pairs']]
    rho_no = z['rho_no']
    i00, i01, i11 = pairs.index((0, 0)), pairs.index((0, 1)), pairs.index((1, 1))
    G = make_grid(KG.NXI, KG.NETA, float(KG.xi_max), a)
    assert np.max(np.abs(G['XI'] - KG.XI)) == 0.0 and np.max(np.abs(G['geo'] - KG.geo_f)) < 1e-15
    pts = G['pts']
    NG = pts.shape[0]
    print(f"  artifact {os.path.basename(REF_EXACT)}: R={R}, exact_neumann={bool(z['exact_neumann'])}, E_trunc={float(z['E_trunc']):.6f}; "
          f"grid {KG.NXI}x{KG.NETA}, xi_max={G['xi_max']:.4f}, NG={NG}   {_t(T0)}")
    # evaluator check vs the artifact's NO values on the grid
    pv = prim_values(prims, pts, R)
    nv = T_no.T @ pv
    ev_err = np.max(np.abs(nv - z['no_val'])) / np.max(np.abs(z['no_val']))
    rho_chk = max(np.max(np.abs(nv[0] * nv[1] - rho_no[i01])), np.max(np.abs(nv[1] ** 2 - rho_no[i11])))
    print(f"  value-only evaluator vs artifact no_val: max rel {ev_err:.1e}; rho_01/rho_11 rebuilt max|d| {rho_chk:.1e}")
    print(f"  INT rho_01 = {grid_int_g(G, rho_no[i01]):+.2e}   INT rho_11 - 1 = {grid_int_g(G, rho_no[i11]) - 1:+.2e}   "
          f"INT rho_00 - 1 = {grid_int_g(G, rho_no[i00]) - 1:+.2e}", flush=True)

    # exact operator on the production grid (Phase 0b) + the LMAX_HI one
    op34 = KG.exact_neumann(LMAX, 0)
    op50 = KG.exact_neumann(LMAX_HI, 0)
    # consistency: my per-l bridge == grid_int(D * neumann_potential(D, 34, exact=True))
    Dt = rho_no[i00]
    c_l = bridge_per_l(Dt, Dt, G, op34, LMAX, R)
    ref_np = grid_int_g(G, Dt * neumann_potential(Dt.reshape(KG.NXI, KG.NETA), LMAX, exact=True).ravel())
    print(f"  bridge_per_l consistency: sum_l c_l = {c_l.sum():.12e} vs grid_int(D*neumann_potential(exact)) = {ref_np:.12e} "
          f"(rel {rel(c_l.sum(), ref_np):.1e});  5 zeta/8 control for rho_00 not applicable (NO0 is a mixture)", flush=True)

    # ------------------------------------------------------------------ closed-form anchors
    print("\n(A) CLOSED FORMS for the isotropic dressings (anchors for both the production kernels and the direct quadrature)")
    s_test = np.array([0.02, 0.3, 1.0, 2.5, 6.0, 12.0])
    for zeta in (1.0, 2.6875, 4.5):
        cf = psi_f_iso(s_test, zeta, GAM)
        fd = f_dress_iso(s_test, zeta)
        cf2 = psi_f_iso(s_test, zeta, 2 * GAM)
        cy = psi_yuk_iso(s_test, zeta, GAM)
        yk = yukawa_pot_iso(s_test, zeta, GAM)
        print(f"  zeta={zeta:6.4f}: Psi_f closed vs energy.f_dress_iso max rel {np.max(np.abs(cf - fd) / np.abs(cf)):.1e};  "
              f"Psi_Y closed vs hT.yukawa_pot_iso max rel {np.max(np.abs(cy - yk) / np.abs(cy)):.1e};  "
              f"Psi_f(0)->8z^3/(2z+g)^3: {psi_f_iso(np.array([1e-6]), zeta, GAM)[0]:.10f} vs {8 * zeta ** 3 / (2 * zeta + GAM) ** 3:.10f};  "
              f"Psi_f2(0): {psi_f_iso(np.array([1e-6]), zeta, 2 * GAM)[0]:.8f} vs {8 * zeta ** 3 / (2 * zeta + 2 * GAM) ** 3:.8f}")
    print(flush=True)

    # ------------------------------------------------------------------ pool
    import multiprocessing as mp
    pool = mp.get_context('spawn').Pool(processes=args.workers) if args.workers > 1 else None

    # ------------------------------------------------------------------ direct-quadrature validation on C0 (closed form)
    print("(B) DIRECT 3-D QUADRATURE validated on the isotropic C0 leaf |1s_A|^2 (beta=2.70) vs the closed form, all 3168 nodes")
    specA = dict(kind='iso', zeta=BETA, center=A_C)
    cf_all = psi_f_iso(G['rA'], BETA, GAM)
    R_NEAR = 5.0
    partsA = [dict(spec=specA, center=A_C, other=B_C, r_near=R_NEAR)]
    val_B = {}
    for rn in (rule_lo, rule_hi):
        t1 = time.time()
        dd = dress_direct(('f', GAM), pts, rn, ctx, partsA, pool)
        near = G['rA'] < R_NEAR
        e = np.abs(dd - cf_all) / np.abs(cf_all)
        val_B[rn] = dd
        print(f"  rule {rn}: max rel err {e.max():.2e} (near/tc {e[near].max():.2e}, far/ac {e[~near].max():.2e}), "
              f"sup-rel {np.max(np.abs(dd - cf_all)) / np.max(cf_all):.2e}   [{time.time() - t1:.0f}s]", flush=True)
    # ------------------------------------------------------------------ configurations
    spec01 = dict(kind='no_pair', p=0, q=1)
    spec11 = dict(kind='no_pair', p=1, q=1)
    spec00 = dict(kind='no_pair', p=0, q=0)

    def two_parts(spec: dict, r_near_A: float, r_near_B: float) -> list:
        return [dict(spec=dict(spec, part='A'), center=A_C, other=B_C, r_near=r_near_A),
                dict(spec=dict(spec, part='B'), center=B_C, other=A_C, r_near=r_near_B)]

    CFG = {
        'C0': dict(leaf=specA, bridge=dict(kind='iso', zeta=ALPHA, center=B_C),
                   leaf_grid=(BETA ** 3 / np.pi) * np.exp(-2.0 * BETA * G['rA']),
                   bridge_grid=(ALPHA ** 3 / np.pi) * np.exp(-2.0 * ALPHA * G['rB']),
                   parts=partsA, center=A_C, closed=cf_all,
                   tag='bridge |1s_B|^2 (1.0), leaves |1s_A|^2 (2.70)',
                   mix_b=[(B_C, 1.0, 1.0)], mix_l=[(A_C, BETA, 1.0)]),
        'C1': dict(leaf=spec01, bridge=spec11, leaf_grid=rho_no[i01], bridge_grid=rho_no[i11],
                   parts=two_parts(spec01, 5.0, 2.5), center=A_C, closed=None,
                   tag='bridge rho_(NO1,NO1), leaves rho_(NO0,NO1) [signed, p-like]',
                   mix_b=[(B_C, 1.0, 0.35), (A_C, 0.65, 0.2), ((A_C + B_C) / 2, 0.55, 0.2), (A_C, 2.7, 0.15), (B_C, 0.6, 0.1)],
                   mix_l=[(A_C, 1.1, 0.15), (A_C, 1.5, 0.25), (A_C, 2.7, 0.35), (A_C, 4.5, 0.2), (B_C, 1.0, 0.05)]),
        'C2': dict(leaf=spec11, bridge=spec00, leaf_grid=rho_no[i11], bridge_grid=rho_no[i00],
                   parts=two_parts(spec11, 8.0, 2.5), center=A_C, closed=None,
                   tag='bridge rho_(NO0,NO0) [core], leaves rho_(NO1,NO1) [bond, positive]',
                   mix_b=[(A_C, 2.7, 0.5), (A_C, 4.5, 0.2), (A_C, 1.6, 0.2), (A_C, 1.0, 0.1)],
                   mix_l=[(B_C, 1.0, 0.35), (A_C, 0.65, 0.2), ((A_C + B_C) / 2, 0.55, 0.2), (A_C, 2.7, 0.15), (B_C, 0.6, 0.1)]),
    }

    # forced-tc vs forced-ac on an overlap band for the REAL leaves (both parts): the two routes must agree
    band = np.nonzero((G['rA'] > 3.0) & (G['rA'] < 8.0) & (G['rB'] > 3.0))[0]
    band = band[np.linspace(0, band.size - 1, 12).astype(int)]
    for nm in ('C1', 'C2'):
        if nm not in configs:
            continue
        d_tc = dress_direct(('f', GAM), pts[band], rule_hi, ctx, CFG[nm]['parts'], pool, mode='tc')
        d_ac = dress_direct(('f', GAM), pts[band], rule_hi, ctx, CFG[nm]['parts'], pool, mode='ac')
        print(f"  {nm} leaf, forced-tc vs forced-ac ({rule_hi}) on 12 nodes with rA in [3,8], rB > 3: max rel "
              f"{np.max(np.abs(d_tc - d_ac) / np.abs(d_tc)):.2e}, sup-rel {np.max(np.abs(d_tc - d_ac)) / np.max(np.abs(d_tc)):.2e} "
              f"(values {np.abs(d_tc).min():.2e}..{np.abs(d_tc).max():.2e})   {_t(T0)}", flush=True)

    # ------------------------------------------------------------------ production kernels (gVee import: ~80 s)
    print("\n(C) production kernels: importing geovac.lih_r12ci.gVee (builds Kf, Kf2, smooth kernels) ...", flush=True)
    t1 = time.time()
    from geovac.lih_r12ci import gVee
    Kf_local = KG.build_kernel(GAM)
    same = np.array_equal(gVee.Kf, Kf_local)
    print(f"  gVee imported [{time.time() - t1:.0f}s]; kernels.build_kernel(GAM) == gVee.Kf bit-for-bit: {same}", flush=True)
    del Kf_local

    def prod_dress(h: np.ndarray, K: np.ndarray) -> np.ndarray:
        return gVee.dress(K, h)

    def prod_psi_yuk(h: np.ndarray) -> np.ndarray:
        KG.USE_EXACT_NEUMANN = True
        try:
            return gVee.psi_yuk(h, GAM)
        finally:
            KG.USE_EXACT_NEUMANN = False

    cache: Dict[str, np.ndarray] = {}
    if os.path.exists(CACHE) and not args.no_cache and not args.quick:
        zc = np.load(CACHE)
        cache = {k: zc[k] for k in zc.files}
        print(f"  cache {os.path.basename(CACHE)} loaded: {len(cache)} arrays")

    results: Dict[str, dict] = {}
    save: Dict[str, np.ndarray] = {}
    if args.merge_out and os.path.exists(OUT) and not args.quick:
        zo = np.load(OUT, allow_pickle=True)
        save = {k: zo[k] for k in zo.files}
        print(f"  merge-out: seeded {len(save)} arrays from {os.path.basename(OUT)}")
    # global metadata (written up-front so a per-config bank produces a self-contained artifact)
    save.update(R=R, a=a, GAM=GAM, BETA=BETA, ALPHA=ALPHA, LMAX=LMAX, LMAX_HI=LMAX_HI, rule_hi=rule_hi, rule_lo=rule_lo,
                configs=np.array(configs), grid_NXI=KG.NXI, grid_NETA=KG.NETA, grid_xi_max=G['xi_max'],
                rho_01=rho_no[i01], rho_11=rho_no[i11], rho_00=rho_no[i00], sigma_gate_mc_bar=SIGMA_GATE_MC_BAR,
                sigma_gate_reduced_limit=SIGMA_GATE_REDUCED_LIMIT, sigma_gate_brute=SIGMA_GATE_BRUTE, quick=args.quick)

    def bank(nm: str) -> None:
        """Serialize one completed config into `save` and rewrite OUT (retune (c): nothing lost on a kill)."""
        r_ = results[nm]
        for k in ('I_conv', 'I_conv_lo', 'I_conv50', 'I_prod', 'I_L2', 'I_L0', 'dPsi', 'I_mc', 's_mc'):
            save[f"{nm}_{k}"] = np.array(r_[k], dtype=float)
        save[f"{nm}_I_scratch"] = np.array(np.nan if r_['I_scratch'] is None else r_['I_scratch'])
        save[f"{nm}_verdict"] = np.array(r_['verdict'])
        for k in ('cl_hi', 'cl_prod', 'Psi_hi', 'Psi_lo', 'Psi_prod', 'aL_leaf', 'aL_bridge', 'aL_Psi', 'aL_D'):
            save[f"{nm}_{k}"] = np.asarray(r_[k])
        if 'I_closed' in r_:
            save[f"{nm}_I_closed"] = np.array(r_['I_closed'])
        if r_.get('mc') is not None:
            mc_ = r_['mc']
            save[f"{nm}_mc_means"] = np.asarray(mc_['mean'])
            save[f"{nm}_mc_errs"] = np.asarray(mc_['err'])
            save[f"{nm}_mc_split"] = np.asarray(mc_['split_sd'])
            save[f"{nm}_mc_N"] = np.array(mc_['N'])
            save[f"{nm}_mc_best"] = np.array(mc_['best'])
            save[f"{nm}_mc_wall"] = np.array(mc_['wall'])
        if r_.get('mc6') is not None:
            m6 = r_['mc6']
            save[f"{nm}_mc6_mean"] = np.array(m6['mean'])
            save[f"{nm}_mc6_err"] = np.array(m6['err'])
            save[f"{nm}_mc6_N"] = np.array(m6['N'])
        if not args.quick:
            np.savez(OUT, **save)

    wXI = barycentric_weights(G['XI'])            # for the 6-D bridge cross-check's grid interpolation
    wETA = barycentric_weights(G['ETA'])
    for name in configs:
        C = CFG[name]
        print("\n" + "=" * 100)
        print(f"({name}) {C['tag']}")
        print("=" * 100, flush=True)
        leaf, bridge = C['leaf'], C['bridge']
        parts = C['parts']
        res: Dict[str, object] = {}
        # ---------------- converged dressing (R1, R2) at the production nodes
        Psi: Dict[str, np.ndarray] = {}
        for rn in (rule_lo, rule_hi):
            key = f"{name}_Psi_{rn}"
            if key in cache:
                Psi[rn] = cache[key]
                print(f"  Psi[{rn}] from cache")
                continue
            t1 = time.time()
            Psi[rn] = dress_direct(('f', GAM), pts, rn, ctx, parts, pool)
            cache[key] = Psi[rn]
            if not args.quick and not args.no_cache:
                np.savez(CACHE, **cache)   # bank the expensive dressing IMMEDIATELY (survives an MC-stage kill)
            print(f"  direct dressing Psi[{rn}] at {NG} nodes [{time.time() - t1:.0f}s]", flush=True)
        dPsi = np.max(np.abs(Psi[rule_hi] - Psi[rule_lo])) / np.max(np.abs(Psi[rule_hi]))
        if C['closed'] is not None:
            e_cf = np.max(np.abs(Psi[rule_hi] - C['closed'])) / np.max(np.abs(C['closed']))
            print(f"  Psi[{rule_hi}] vs closed form: sup-rel {e_cf:.2e}")
        print(f"  Psi[{rule_hi}] vs Psi[{rule_lo}]: sup-rel {dPsi:.2e};  max|Psi| = {np.max(np.abs(Psi[rule_hi])):.6e}")
        # ---------------- bridge on the production grid: converged (i), production (ii)
        D_hi = C['bridge_grid'] * Psi[rule_hi]
        D_lo = C['bridge_grid'] * Psi[rule_lo]
        cl_hi = bridge_per_l(D_hi, D_hi, G, op50, LMAX_HI, R)
        cl_lo = bridge_per_l(D_lo, D_lo, G, op34, LMAX, R)
        I_conv = float(cl_hi[:LMAX + 1].sum())
        I_conv_lo = float(cl_lo.sum())
        I_conv50 = float(cl_hi.sum())
        if C['closed'] is not None:
            D_cf = C['bridge_grid'] * C['closed']
            I_cf = float(bridge_per_l(D_cf, D_cf, G, op34, LMAX, R).sum())
            print(f"  (i) converged: I[closed-form Psi] = {I_cf:.10e}  I[{rule_hi}] = {I_conv:.10e} (rel {rel(I_conv, I_cf):.1e})  "
                  f"I[{rule_lo}] = {I_conv_lo:.10e} (rel {rel(I_conv_lo, I_cf):.1e})")
            res['I_closed'] = I_cf
        else:
            print(f"  (i) converged: I[{rule_hi}] = {I_conv:.10e}   I[{rule_lo}] = {I_conv_lo:.10e}   "
                  f"rel({rule_hi},{rule_lo}) = {rel(I_conv, I_conv_lo):.2e}  ({'OK' if rel(I_conv, I_conv_lo) <= CONV_REL else 'NOT CONVERGED'} at {CONV_REL:.0e})")
        print(f"      L-series: LMAX 34 -> 50 shift = {(I_conv50 - I_conv) / I_conv:+.2e} rel;  per-l fractions l=0..7: "
              + " ".join(f"{cl_hi[l] / I_conv:+.3e}" for l in range(8)))
        print(f"      tail fractions: l>=10 {cl_hi[10:].sum() / I_conv:.2e}, l>=20 {cl_hi[20:].sum() / I_conv:.2e}, l>=35 {cl_hi[35:].sum() / I_conv:.2e}")
        t1 = time.time()
        Psi_prod = prod_dress(C['leaf_grid'], gVee.Kf)
        D_prod = C['bridge_grid'] * Psi_prod
        V_prod = neumann_potential(D_prod.reshape(KG.NXI, KG.NETA), LMAX, exact=True).ravel()
        I_prod = grid_int_g(G, D_prod * V_prod)
        cl_prod = bridge_per_l(D_prod, D_prod, G, op34, LMAX, R)
        print(f"  (ii) production: I = {I_prod:.10e}  (per-l route {cl_prod.sum():.10e}, rel {rel(cl_prod.sum(), I_prod):.1e})   [{time.time() - t1:.1f}s]")
        print(f"      (i)-(ii): {(I_conv - I_prod):+.3e} = {(I_conv - I_prod) / I_conv:+.2e} rel  <- the coarse-grid f-dressing error on I")
        print(f"      Psi_prod vs Psi[{rule_hi}]: sup-rel {np.max(np.abs(Psi_prod - Psi[rule_hi])) / np.max(np.abs(Psi[rule_hi])):.2e};  "
              f"INT bridge*Psi: conv {grid_int_g(G, D_hi):.8e} prod {grid_int_g(G, D_prod):.8e} (rel {rel(grid_int_g(G, D_prod), grid_int_g(G, D_hi)):.1e})")
        # L-truncation and the isotropic-shortcut wrong answers (for the guard's built-in fire tests)
        I_L2 = float(cl_hi[:3].sum())
        I_L0 = float(cl_hi[:1].sum())
        # Legendre content about A
        aL_leaf = legendre_content(C['leaf_grid'], G)
        aL_bridge = legendre_content(C['bridge_grid'], G)
        aL_Psi = legendre_content(Psi[rule_hi], G)
        aL_D = legendre_content(D_hi, G)
        print(f"  Legendre content about A (a0..a4): leaf " + " ".join(f"{x:.3f}" for x in aL_leaf)
              + " | dressed leaf Psi " + " ".join(f"{x:.3f}" for x in aL_Psi)
              + " | bridge " + " ".join(f"{x:.3f}" for x in aL_bridge)
              + " | D=bridge*Psi " + " ".join(f"{x:.3f}" for x in aL_D), flush=True)
        # ---------------- scratch-grid Neumann check of the dressed density
        I_scr = None
        if not args.no_scratch:
            NXs, NEs = (48, 30) if args.quick else (90, 54)
            t1 = time.time()
            Gs = make_grid(NXs, NEs, G['xi_max'], a)
            key = f"{name}_Psi_scratch"
            if key in cache and cache[key].shape[0] == Gs['pts'].shape[0]:
                Psi_s = cache[key]
            else:
                Psi_s = dress_direct(('f', GAM), Gs['pts'], rule_lo, ctx, parts, pool)
                cache[key] = Psi_s
            rho_b_s = eval_density(bridge, Gs['pts'], ctx)
            op_s = ExactNeumann(Gs['XI'], 1.0, Gs['xi_max'], LMAX, 0)
            D_s = rho_b_s * Psi_s
            I_scr = float(bridge_per_l(D_s, D_s, Gs, op_s, LMAX, R).sum())
            print(f"  scratch grid {NXs}x{NEs} (rule {rule_lo} dressing, own exact operator): I = {I_scr:.10e}  "
                  f"vs (i)[{rule_lo}] on 72x44 {I_conv_lo:.10e}: rel {rel(I_scr, I_conv_lo):.2e}   [{time.time() - t1:.0f}s]", flush=True)
        # ---------------- Monte Carlo
        mc = None
        if not args.skip_mc:
            # leaf moments about the control centre (grid) + independent A-centred spherical check
            c0 = C['center']
            rho_l = C['leaf_grid']
            uz = G['ZC'] - c0[2]
            cv = dict(center=c0, M0=grid_int_g(G, rho_l), dz=grid_int_g(G, rho_l * uz),
                      Qxx=0.5 * grid_int_g(G, rho_l * G['RHO_CYL'] ** 2), Qzz=grid_int_g(G, rho_l * uz ** 2), s0=1.2)
            # independent moment check: leaf-centred spherical quadrature (r, x) with the ac-rule panels
            rb = graded_bounds(0.0, 36.0, [0.0, R], 0.004, 1.7)
            rs, wr = gl_on_panels(rb, 10)
            xb = graded_bounds(-1.0, 1.0, [1.0], 1e-7, 3.0)
            xs, wx = gl_on_panels(xb, 10)
            RCm = rs[:, None] * np.sqrt(1 - xs ** 2)[None, :]
            ZZm = c0[2] + rs[:, None] * xs[None, :]
            ptsm = np.stack([RCm.ravel(), np.zeros(RCm.size), ZZm.ravel()], axis=1)
            rhom = eval_density(leaf, ptsm, ctx).reshape(rs.size, xs.size)
            Wm = 2 * np.pi * (wr * rs ** 2)[:, None] * wx[None, :] * rhom
            chk = dict(M0=float(Wm.sum()), dz=float(np.sum(Wm * (rs[:, None] * xs[None, :]))),
                       Qxx=float(0.5 * np.sum(Wm * RCm ** 2)), Qzz=float(np.sum(Wm * (rs[:, None] * xs[None, :]) ** 2)))
            print(f"  MC control moments about c={c0}: M0={cv['M0']:+.3e} dz={cv['dz']:+.6f} Qxx={cv['Qxx']:+.6f} Qzz={cv['Qzz']:+.6f}; "
                  f"independent spherical quadrature: M0 {chk['M0']:+.3e} dz {chk['dz']:+.6f} Qxx {chk['Qxx']:+.6f} Qzz {chk['Qzz']:+.6f}")
            mcfg = dict(ctx=ctx, gam=GAM, spec_b=bridge, spec_l=leaf, mix_b=C['mix_b'], mix_l=C['mix_l'], cv=cv)
            batch = 100_000 if args.quick else 1_000_000
            n_pilot = 2 * args.workers
            t1 = time.time()
            bm = mc_run(mcfg, n_pilot, batch, pool, args.workers, seed0=20260923)
            dt = time.time() - t1
            sm = mc_summary(bm)
            sig1 = bm.std(axis=0, ddof=1) * np.sqrt(batch)       # per-sample std
            I_ref = I_conv
            print(f"  MC pilot: {n_pilot} batches x {batch} in {dt:.0f}s ({n_pilot * batch / dt:.2e} samples/s); "
                  f"estimators [plain, dipole, dip+quad]: means {sm['mean']}  batch-err {sm['err']}")
            print(f"    per-sample rel std: {sig1 / abs(I_ref)};  N needed for {GATE_REL:.0e}: {(sig1 / (GATE_REL * abs(I_ref))) ** 2}")
            best = int(np.argmin(sig1))
            N_need = float((sig1[best] / (GATE_REL * abs(I_ref))) ** 2)
            rate = n_pilot * batch / dt
            N_cap = rate * 60.0 * MC_MINUTES[name] * (args.mc_minutes / 22.0)
            # aim at sigma_MC <= 0.5 x gate (margin for the heavy-tail optimism of batch means), floor 120 batches
            N_run = int(min(max(4.0 * N_need, (12 if args.quick else 120) * batch), N_cap))
            n_batches = max(int(N_run // batch), 3 * args.workers)
            print(f"    chosen estimator #{best}; N_needed(1e-3)={N_need:.2e}, cap (budget) {N_cap:.2e} -> running {n_batches} batches "
                  f"({n_batches * batch:.2e} samples, est. {n_batches * batch / rate / 60:.1f} min)", flush=True)
            t1 = time.time()
            bm2 = mc_run(mcfg, n_batches, batch, pool, args.workers, seed0=20260924)
            bm_all = np.concatenate([bm, bm2], axis=0)
            sm2 = mc_summary(bm_all)
            wall_mc = time.time() - t1 + dt
            N_tot = bm_all.shape[0] * batch
            mc = dict(mean=sm2['mean'], err=sm2['err'], split_sd=sm2['split_sd'], best=best, N=N_tot, wall=wall_mc,
                      n_batches=bm_all.shape[0], sig1=sig1)
            for k, lab in enumerate(('plain', 'dipole', 'dip+quad')):
                print(f"    {lab:9s}: I_MC = {sm2['mean'][k]:.8e} +/- {sm2['err'][k]:.1e} (batch means; 3-way split scatter {sm2['split_sd'][k]:.1e})"
                      f"  rel bar {sm2['err'][k] / abs(sm2['mean'][k]):.1e};  vs (i) {(sm2['mean'][k] - I_conv) / I_conv:+.2e} rel "
                      f"= {(sm2['mean'][k] - I_conv) / max(sm2['err'][k], 1e-300):+.1f} sigma;  vs (ii) {(sm2['mean'][k] - I_prod) / I_prod:+.2e} rel")
            print(f"    N = {N_tot:.3e} samples, MC wall {wall_mc / 60:.1f} min", flush=True)
        # ---------------- 6-D bridge cross-check (isolates the Coulomb/Neumann step; carries the gate if the 12-D can't)
        mc6 = None
        if not args.skip_mc:
            rng6 = np.random.default_rng(20260925 + sum(ord(ch) for ch in name))
            m6 = mc6_run(D_hi.reshape(KG.NXI, KG.NETA), C['mix_b'], G, wXI, wETA, R, lagrange_matrix, rng6,
                         budget_s=60.0 * args.mc6_minutes)
            mc6 = m6
            print(f"    6-D bridge (grid-Psi[{rule_hi}] interp, MC over r1,r3 ~ mix_b): I_MC6 = {m6['mean']:.8e} "
                  f"+/- {m6['err']:.1e}  rel bar {m6['err'] / abs(m6['mean']):.1e};  vs (i) {(m6['mean'] - I_conv) / I_conv:+.2e} rel "
                  f"= {(m6['mean'] - I_conv) / max(m6['err'], 1e-300):+.1f} sigma  (N={m6['N']:.2e}, {m6['wall'] / 60:.1f} min)", flush=True)
        # ---------------- verdict per config
        carrier = 'none'
        if mc is not None:
            k = mc['best']
            I12, s12 = float(mc['mean'][k]), float(mc['err'][k])
            bar = GATE_REL * abs(I_conv)
            conv_ok = (rel(I_conv, I_conv_lo) <= CONV_REL) if C['closed'] is None else (rel(I_conv, res['I_closed']) <= CONV_REL)
            # pick the MC route that reaches the 1e-3|I| bar; prefer the 12-D full-brute, fall back to the 6-D bridge
            s6 = float(mc6['err']) if mc6 is not None else float('inf')
            I6 = float(mc6['mean']) if mc6 is not None else float('nan')
            if s12 <= bar:
                carrier, I_mc, s_mc = '12D', I12, s12
            elif s6 <= bar:
                carrier, I_mc, s_mc = '6D', I6, s6
            else:
                carrier = '12D' if s12 <= s6 else '6D'
                I_mc, s_mc = (I12, s12) if s12 <= s6 else (I6, s6)
            tol = max(bar, 2.0 * s_mc)
            dev = abs(I_conv - I_mc)
            if s_mc > bar:
                verdict = 'INCONCLUSIVE'
            elif dev <= tol and conv_ok:
                verdict = 'PASS'
            elif not conv_ok:
                verdict = 'INCONCLUSIVE'
            else:
                verdict = 'FAIL'
            loc = ''
            if verdict != 'PASS' and mc6 is not None and s6 <= bar:
                # the 6-D uses the trusted grid dressing, so its (dis)agreement localises the residual
                loc = f";  6-D bridge vs (i) {(I6 - I_conv) / I_conv:+.2e} (dressing-independent -> residual is in the {'bridge' if abs(I6 - I_conv) > tol else 'dressing/MC-variance'})"
            print(f"  VERDICT {name}: {verdict}  carrier={carrier}  |(i) - MC| = {dev:.2e} vs tol {tol:.2e} (= max({GATE_REL:.0e}|I|, 2 sigma_MC)), "
                  f"sigma_MC/|I| = {s_mc / abs(I_conv):.1e}, converged={conv_ok}{loc}", flush=True)
        else:
            verdict = 'NO-MC'
            I_mc, s_mc = float('nan'), float('nan')
        res.update(I_conv=I_conv, I_conv_lo=I_conv_lo, I_conv50=I_conv50, I_prod=I_prod, I_L2=I_L2, I_L0=I_L0, I_scratch=I_scr,
                   cl_hi=cl_hi, cl_prod=cl_prod, Psi_hi=Psi[rule_hi], Psi_lo=Psi[rule_lo], Psi_prod=Psi_prod, dPsi=dPsi,
                   aL_leaf=aL_leaf, aL_bridge=aL_bridge, aL_Psi=aL_Psi, aL_D=aL_D, verdict=verdict, mc=mc, mc6=mc6,
                   carrier=carrier, I_mc=I_mc, s_mc=s_mc)
        results[name] = res
        # persist the cache and BANK this config's result (retune (c): C1 is safe before C2 starts)
        if not args.quick:
            np.savez(CACHE, **cache)
        bank(name)

    # ------------------------------------------------------------------ (D) dressing-error table
    print("\n" + "=" * 100)
    print("(D) DRESSING ERRORS on the 72x44 production grid (Kf, Kf2 [28-pt phi kernels], psi_yuk [exact Neumann - smooth])")
    print("=" * 100)
    print("  metrics: sup = max|Delta|/max|exact|;  ptw = max |Delta|/|exact| over nodes with |exact| >= 1e-3 max;  "
          "int = rel error of INT rho_bridge * Psi (the quantity entering I; bridge = the config's bridge density)")
    rows_D = []
    # 'int' partner: a 1s_A(2.70) density for the dressed 1s_B, a 1s_B(1.0) density for the dressed 1s_A's
    rho_A_iso = (BETA ** 3 / np.pi) * np.exp(-2.0 * BETA * G['rA'])
    rho_B_iso = (ALPHA ** 3 / np.pi) * np.exp(-2.0 * ALPHA * G['rB'])
    iso_rows = [(1.0, 'B', G['rB'], rho_A_iso), (2.6875, 'A', G['rA'], rho_B_iso), (4.5, 'A', G['rA'], rho_B_iso)]

    def metrics(got: np.ndarray, ex: np.ndarray, wb: np.ndarray) -> Tuple[float, float, float]:
        msk = np.abs(ex) >= 1e-3 * np.max(np.abs(ex))
        return (float(np.max(np.abs(got - ex)) / np.max(np.abs(ex))),
                float(np.max(np.abs(got - ex)[msk] / np.abs(ex)[msk])),
                float(rel(grid_int_g(G, wb * got), grid_int_g(G, wb * ex))))

    header = f"  {'density':22s} | {'Kf vs exact':>30s} | {'Kf2 vs exact(2g)':>30s} | {'psi_yuk vs exact':>30s}"
    print(header)
    print(f"  {'':22s} | {'sup':>9s} {'ptw':>9s} {'int':>9s} | {'sup':>9s} {'ptw':>9s} {'int':>9s} | {'sup':>9s} {'ptw':>9s} {'int':>9s}")
    for zeta, cen, r, wb in iso_rows:
        rho = (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * r)
        ex_f = psi_f_iso(r, zeta, GAM)
        ex_f2 = psi_f_iso(r, zeta, 2 * GAM)
        ex_y = psi_yuk_iso(r, zeta, GAM)
        m1 = metrics(prod_dress(rho, gVee.Kf), ex_f, wb)
        m2 = metrics(prod_dress(rho, gVee.Kf2), ex_f2, wb)
        m3 = metrics(prod_psi_yuk(rho), ex_y, wb)
        rows_D.append((f"1s zeta={zeta} on {cen}", m1, m2, m3))
        # the corpus's own 1-D "exact" functions vs the closed forms (bonus)
        e_fd = np.max(np.abs(f_dress_iso(r, zeta) - ex_f) / np.abs(ex_f))
        e_yk = np.max(np.abs(yukawa_pot_iso(r, zeta, GAM) - ex_y) / np.abs(ex_y))
        print(f"  {rows_D[-1][0]:22s} | " + " | ".join(" ".join(f"{v:9.2e}" for v in m) for m in (m1, m2, m3))
              + f"   [energy.f_dress_iso vs closed: {e_fd:.1e}; hT.yukawa_pot_iso vs closed: {e_yk:.1e}]")
    # the real leaf rho_01: Kf vs the converged direct dressing at ALL nodes; Kf2 and psi_yuk at 20 nodes (tc rule)
    if 'C1' in results:
        C = CFG['C1']
        rho = C['leaf_grid']
        wb = C['bridge_grid']
        Psi_hi = results['C1']['Psi_hi']
        m1 = metrics(prod_dress(rho, gVee.Kf), Psi_hi, wb)
        # 20 nodes: stratified in r_A over the region where the bridge is significant
        cand = np.nonzero(wb > 1e-3 * wb.max())[0]
        cand = cand[np.argsort(G['rA'][cand])]
        sel = cand[np.linspace(0, cand.size - 1, 20).astype(int)]
        t1 = time.time()
        d_f2 = dress_direct(('f', 2 * GAM), pts[sel], rule_hi, ctx, C['parts'], pool, mode='tc')
        d_y = dress_direct(('yuk', GAM), pts[sel], rule_hi, ctx, C['parts'], pool, mode='tc')
        d_f2_lo = dress_direct(('f', 2 * GAM), pts[sel], rule_lo, ctx, C['parts'], pool, mode='tc')
        d_y_lo = dress_direct(('yuk', GAM), pts[sel], rule_lo, ctx, C['parts'], pool, mode='tc')
        p_f2 = prod_dress(rho, gVee.Kf2)[sel]
        p_y = prod_psi_yuk(rho)[sel]
        sup2 = np.max(np.abs(p_f2 - d_f2)) / np.max(np.abs(d_f2))
        supy = np.max(np.abs(p_y - d_y)) / np.max(np.abs(d_y))
        ptw2 = np.max(np.abs(p_f2 - d_f2) / np.maximum(np.abs(d_f2), 1e-3 * np.max(np.abs(d_f2))))
        ptwy = np.max(np.abs(p_y - d_y) / np.maximum(np.abs(d_y), 1e-3 * np.max(np.abs(d_y))))
        int2 = rel(float(np.sum(wb[sel] * p_f2)), float(np.sum(wb[sel] * d_f2)))
        inty = rel(float(np.sum(wb[sel] * p_y)), float(np.sum(wb[sel] * d_y)))
        m2 = (float(sup2), float(ptw2), float(int2))
        m3 = (float(supy), float(ptwy), float(inty))
        rows_D.append(("real leaf rho_(NO0,NO1)", m1, m2, m3))
        print(f"  {'real leaf rho_(NO0,NO1)':22s} | " + " | ".join(" ".join(f"{v:9.2e}" for v in m) for m in (m1, m2, m3))
              + f"   [Kf: all {NG} nodes vs direct {rule_hi}; Kf2/psi_yuk: 20 nodes, rA {G['rA'][sel].min():.2f}..{G['rA'][sel].max():.2f}; "
              f"direct {rule_hi} vs {rule_lo} at those nodes: f2 {np.max(np.abs(d_f2 - d_f2_lo)) / np.max(np.abs(d_f2)):.1e}, yuk {np.max(np.abs(d_y - d_y_lo)) / np.max(np.abs(d_y)):.1e}; "
              f"'int' here = weighted sum over the 20 nodes]   [{time.time() - t1:.0f}s]", flush=True)
        save.update(D_sel_nodes=sel, D_sel_f2_direct=d_f2, D_sel_yuk_direct=d_y, D_sel_f2_prod=p_f2, D_sel_yuk_prod=p_y)
    print(flush=True)

    # ------------------------------------------------------------------ summary + save
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    for name in configs:
        r_ = results[name]
        mc = r_['mc']
        line = (f"  {name}: (i) {r_['I_conv']:.8e}  (ii) {r_['I_prod']:.8e}  (i)-(ii) {(r_['I_conv'] - r_['I_prod']) / r_['I_conv']:+.2e}"
                f"  {rule_hi}/{rule_lo} {rel(r_['I_conv'], r_['I_conv_lo']):.1e}  L34->50 {(r_['I_conv50'] - r_['I_conv']) / r_['I_conv']:+.1e}")
        if r_['I_scratch'] is not None:
            line += f"  scratch {rel(r_['I_scratch'], r_['I_conv_lo']):.1e}"
        if mc is not None:
            k = mc['best']
            line += (f"  MC12 {mc['mean'][k]:.8e} +/- {mc['err'][k]:.1e} (N={mc['N']:.2e})"
                     f"  (i)-MC12 {(r_['I_conv'] - mc['mean'][k]) / r_['I_conv']:+.2e}")
        if r_.get('mc6') is not None:
            m6 = r_['mc6']
            line += f"  MC6 {m6['mean']:.8e} +/- {m6['err']:.1e}  (i)-MC6 {(r_['I_conv'] - m6['mean']) / r_['I_conv']:+.2e}"
        line += f"  [carrier {r_.get('carrier', 'none')}] -> {r_['verdict']}"
        print(line)
    if 'C0' in results:
        r0 = results['C0']
        print(f"  C0 vs the sigma-gate record: reduced(grid limit) {SIGMA_GATE_REDUCED_LIMIT:.6e} (rel to (i) {rel(SIGMA_GATE_REDUCED_LIMIT, r0['I_conv']):+.1e}), "
              f"brute 3x120M {SIGMA_GATE_BRUTE:.6e} (rel {rel(SIGMA_GATE_BRUTE, r0['I_conv']):+.1e})")
    print(f"wall {time.time() - T0:.0f} s", flush=True)

    # frozen artifact for the guard (per-config keys + metadata already banked incrementally by bank();
    # here we add only the section-(D) dressing table and the C1 spherical-average leaf for the fire test)
    save['dressing_table'] = np.array([[*m1_, *m2_, *m3_] for (_, m1_, m2_, m3_) in rows_D])
    save['dressing_rows'] = np.array([r_[0] for r_ in rows_D])
    for name in configs:
        if results[name].get('mc6') is not None:
            save[f"{name}_mc6_mean"] = np.array(results[name]['mc6']['mean'])
            save[f"{name}_mc6_err"] = np.array(results[name]['mc6']['err'])
    # spherical average of the C1 leaf about A (for the isotropic-shortcut fire test in the guard)
    if 'C1' in results:
        rb = graded_bounds(0.0, 30.0, [0.0, R], 0.004, 1.7)
        rs, _ = gl_on_panels(rb, 10)
        xg, wg = leggauss(160)
        RCm = rs[:, None] * np.sqrt(1 - xg ** 2)[None, :]
        ZZm = A_C[2] + rs[:, None] * xg[None, :]
        ptsm = np.stack([RCm.ravel(), np.zeros(RCm.size), ZZm.ravel()], axis=1)
        rhom = eval_density(spec01, ptsm, ctx).reshape(rs.size, xg.size)
        save['C1_leaf_sph_r'] = rs
        save['C1_leaf_sph_rho'] = 0.5 * (rhom @ wg)
    if not args.quick:
        np.savez(OUT, **save)
        print(f"saved {OUT} ({os.path.getsize(OUT) / 1e6:.2f} MB)")
    else:
        print("(quick mode: artifact NOT written)")
    if pool is not None:
        pool.close()
        pool.join()
    return 0


if __name__ == '__main__':
    sys.exit(main())
