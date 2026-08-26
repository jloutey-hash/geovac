"""BeH2 conical-intersection landscape with EXACT (non-interpolated) two-center overlaps.

Goal (2026-08-25 follow-up): the interpolated-overlap Berry phase FLIPS with loop radius
around the central CI (pi at r=0.15,0.22; 0 at r=0.3; pi at 0.38,0.5), suggesting additional
CIs off the symmetric line.  But that used cubic interpolation of a coarse (0.6-bohr-spaced)
1D overlap grid.  Here we recompute with the ACTUAL overlap at every geometry and decide:
  GO       = >=1 additional CI confirmed with exact overlaps.
  ARTIFACT = the outer flips vanish (only the central CI is real).
  PARTIAL  = inconclusive within budget.

The mpmath overlap (compute_topos3_two_center_meet.overlap_two_center) is ~14s/eval (55s per
2x2 S_block) -- far too slow for a grid + loops.  The integrand is polynomial x exp(-a*xi) x
exp(-b*eta) x (xi^2-eta^2) in prolate spheroidal coords (all s,p, so the r^l * P_lm rationals
cancel to polynomials).  We evaluate it EXACTLY per geometry with Gauss-Laguerre in xi (exact
for poly x exp) + Gauss-Legendre in eta, in float64 -- validated below to ~1e-9 vs mpmath.
This is a true overlap at each (d1,d2), NOT an interpolation across geometries.
"""
from __future__ import annotations
import json, time, importlib.util
import numpy as np
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_genlaguerre, lpmv, factorial

STATES = [(1, 0), (2, 1)]              # 1s, 2p0  (sigma pair, m=0)
PAR = np.diag([1.0, -1.0])            # 2p0 odd under z->-z
NS = 2

# ---- fast float64 hydrogenic radial (modern convention, matches topos3 R_nl_modern) ----
def R_nl(Z, n, l, r):
    Z = float(Z); rho = 2 * Z * r / n
    norm = np.sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * eval_genlaguerre(n - l - 1, 2 * l + 1, rho)

def angular(l, m, ct):
    m = abs(m)
    tn = np.sqrt((2 * l + 1) / 2.0 * factorial(l - m) / factorial(l + m))
    return tn * lpmv(m, l, ct)        # scipy lpmv carries the Condon-Shortley (-1)^m

# ---- Gauss nodes (fixed, high order) ----
_LAG_X, _LAG_W = laggauss(64)         # int_0^inf f(t) e^{-t} dt = sum w f(t)
_LEG_X, _LEG_W = leggauss(96)         # eta in [-1,1]

def overlap_fast(Z1, n1, l1, Z2, n2, l2, m, R):
    """<phi_{n1 l1 m}(0,Z1) | phi_{n2 l2 m}(R zhat, Z2)>, exact per-geometry quadrature."""
    R = float(R); half = R / 2.0
    a = half * (Z1 / n1 + Z2 / n2)     # xi-decay rate: integrand ~ e^{-a*xi}
    # xi = 1 + s/a  (s = Gauss-Laguerre variable);  int_1^inf g(xi) dxi
    #   = (e^{-a}/a) * sum_k w_k [ g(1+s_k/a) e^{+a*(xi-1)} ]  with the e^{-a*xi} folded in.
    xi = 1.0 + _LAG_X / a                                   # (Nxi,)
    eta = _LEG_X                                            # (Neta,)
    XI, ETA = np.meshgrid(xi, eta, indexing="ij")          # (Nxi,Neta)
    r1 = half * (XI + ETA); r2 = half * (XI - ETA)
    with np.errstate(divide="ignore", invalid="ignore"):
        ct1 = (1 + XI * ETA) / (XI + ETA)
        ct2 = (XI * ETA - 1) / (XI - ETA)
    # exp(-b*eta) part is inside R_nl(r2) already; strip the e^{-a*xi} that GL weight supplies.
    val = (R_nl(Z1, n1, l1, r1) * angular(l1, m, ct1)
           * R_nl(Z2, n2, l2, r2) * angular(l2, m, ct2)
           * half ** 3 * (XI ** 2 - ETA ** 2))
    bad = (r1 <= 0) | (r2 <= 0) | ~np.isfinite(val)
    val = np.where(bad, 0.0, val)
    # restore the exponential weight: integrand = val, GL gives int val dxi via
    #   int_1^inf val dxi = (1/a) int_0^inf val(1+s/a) ds = (1/a) sum w_k e^{s_k} val_k
    # but val already contains e^{-a*xi}=e^{-a} e^{-s_k}; e^{-s_k} cancels GL's e^{-s_k} -> keep e^{s} factor:
    wl = (_LAG_W * np.exp(_LAG_X) / a)[:, None]            # (Nxi,1) restore full measure
    return float(np.sum(wl * (_LEG_W[None, :] * val)))

def S_block_fast(Rsep, Za, Zb):
    S = np.zeros((NS, NS))
    for i, (na, la) in enumerate(STATES):
        for j, (nb, lb) in enumerate(STATES):
            S[i, j] = overlap_fast(Za, na, la, Zb, nb, lb, 0, Rsep)
    return S

# ---- exact projectors for asymmetric (d1,d2) using cached S_block_fast ----
_CACHE_BEH = {}; _CACHE_HH = {}
def _beh(d):
    k = round(d, 6)
    if k not in _CACHE_BEH: _CACHE_BEH[k] = S_block_fast(d, 2, 1)
    return _CACHE_BEH[k]
def _hh(s):
    k = round(s, 6)
    if k not in _CACHE_HH: _CACHE_HH[k] = S_block_fast(s, 1, 1)
    return _CACHE_HH[k]

def projectors_exact(d1, d2):
    S_BeH1 = _beh(d1)
    S_BeH2 = PAR @ _beh(d2) @ PAR
    S_HH = PAR @ _hh(d1 + d2) @ PAR
    I = np.eye(NS)
    G = np.block([[I, S_BeH1, S_BeH2], [S_BeH1.T, I, S_HH], [S_BeH2.T, S_HH.T, I]])
    if np.linalg.eigvalsh(G).min() < 1e-9:
        return None
    Xh = np.linalg.cholesky(G).T
    return [Xh[:, 2 * k:2 * k + 2] @ np.linalg.pinv(Xh[:, 2 * k:2 * k + 2]) for k in range(3)]

def fgap(d1, d2):
    Ps = projectors_exact(d1, d2)
    if Ps is None: return np.nan
    w = np.linalg.eigvalsh(Ps[0] + Ps[1] + Ps[2])
    return float(np.min(np.diff(w)))

def berry(cx, cy, r, N=48):
    """Z2 Berry phase: parallel-transport lower-band-of-min-gap F-eigenvector; sign(v_end.v0)."""
    v0 = vp = None
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        Ps = projectors_exact(cx + r * np.cos(th), cy + r * np.sin(th))
        if Ps is None: return None
        w, V = np.linalg.eigh(Ps[0] + Ps[1] + Ps[2])
        if i == 0:
            k = int(np.argmin(np.diff(w))); v = V[:, k]; v0 = v.copy(); vp = v
        else:
            ov = V.T @ vp; j = int(np.argmax(np.abs(ov)))
            v = V[:, j] * np.sign(ov[j]); vp = v
    return float(np.sign(vp @ v0))

# ---------------------------------------------------------------- validation
def _validate():
    spec = importlib.util.spec_from_file_location("t3", "debug/compute_topos3_two_center_meet.py")
    t3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t3)
    import mpmath as mp; mp.mp.dps = 25
    cases = [(1,1,0,1,1,0,0,2.445), (2,1,0,1,1,0,0,2.445), (2,2,1,1,1,0,0,2.445),
             (2,2,1,1,2,1,0,2.445), (1,1,0,1,1,0,0,4.89), (1,2,1,1,2,1,0,4.89)]
    print("--- validation: fast float64 vs mpmath adaptive quad ---")
    worst = 0.0
    for c in cases:
        ref = float(t3.overlap_two_center(*c))
        got = overlap_fast(*c)
        err = abs(got - ref)
        worst = max(worst, err)
        print(f"  {c}: mpmath={ref:+.10f} fast={got:+.10f} |err|={err:.2e}")
    print(f"  worst abs err = {worst:.2e}")
    return worst

# ---------------------------------------------------------------- robust per-band Berry
def berry_band(cx, cy, r, N=360):
    """Full-frame continuity-tracked Z2 phase; returns per-band sign(<v_end|v_0>) for all 6
    bands.  Robust to which-vector-in-a-degenerate-pair (unlike single-vector transport),
    but still corrupted when the loop GRAZES a CI (passes within ~loop-step of it)."""
    Vprev = None; frames = []
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        Ps = projectors_exact(cx + r * np.cos(th), cy + r * np.sin(th))
        w, V = np.linalg.eigh(Ps[0] + Ps[1] + Ps[2])
        if Vprev is None:
            Vc = V.copy()
        else:
            M = Vprev.T @ V; perm = np.argmax(np.abs(M), axis=0); Vc = np.zeros_like(V)
            for j in range(6):
                Vc[:, perm[j]] = V[:, j] * np.sign(M[perm[j], j])
        frames.append(Vc); Vprev = Vc
    return [float(np.sign(frames[-1][:, k] @ frames[0][:, k])) for k in range(6)]


def run_landscape():
    """Full EXACT-overlap CI-landscape analysis; writes debug/data/beh2_ci_exact_landscape.json."""
    from scipy.optimize import minimize
    t0 = time.time(); out = {}

    # [1] central CI on the symmetric line
    ds = np.linspace(1.5, 3.5, 201); g = np.array([fgap(d, d) for d in ds])
    i0 = int(np.argmin(g))
    r1 = minimize(lambda d: fgap(d[0], d[0]), [ds[i0]], method="Nelder-Mead",
                  options=dict(xatol=1e-7, fatol=1e-14))
    dstar = float(r1.x[0]); gstar = float(r1.fun)
    out["central_CI"] = dict(dstar=dstar, gap=gstar)
    print(f"[1] central CI: d*={dstar:.6f}  exact F-gap={gstar:.2e}", flush=True)

    # [2] complete 2D gap census over the full loop domain [1.8,3.1]^2
    box = np.linspace(1.8, 3.1, 66); GAP = np.full((66, 66), np.nan)
    for ia, d1 in enumerate(box):
        for ib, d2 in enumerate(box):
            GAP[ia, ib] = fgap(d1, d2)

    def islm(A, i, j):
        v = A[i, j]
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == dj == 0: continue
                ii, jj = i + di, j + dj
                if 0 <= ii < 66 and 0 <= jj < 66 and A[ii, jj] < v: return False
        return True
    cands = sorted([(box[i], box[j], GAP[i, j]) for i in range(66) for j in range(66)
                    if islm(GAP, i, j)], key=lambda t: t[2])
    true = []
    for (d1c, d2c, gc) in cands[:12]:
        rr = minimize(lambda x: fgap(x[0], x[1]), [d1c, d2c], method="Nelder-Mead",
                      options=dict(xatol=1e-7, fatol=1e-14, maxiter=600))
        a, b, c = float(rr.x[0]), float(rr.x[1]), float(rr.fun)
        if c < 5e-4 and 1.8 < a < 3.1 and 1.8 < b < 3.1 and \
           not any(abs(a - x) < 3e-3 and abs(b - y) < 3e-3 for (x, y, _) in true):
            true.append((a, b, c))
    true.sort()
    out["census_box"] = [1.8, 3.1]
    out["true_CIs"] = [dict(d1=a, d2=b, gap=c, on_axis=bool(abs(a - b) < 3e-3),
                            dist_from_center=float(np.hypot(a - dstar, b - dstar)))
                       for (a, b, c) in true]
    out["true_CI_count"] = len(true)
    print(f"[2] complete census [1.8,3.1]^2: {len(true)} conical intersections", flush=True)
    for (a, b, c) in true:
        print(f"    ({a:.5f},{b:.5f}) gap={c:.1e} "
              f"{'ON-AXIS' if abs(a-b)<3e-3 else 'off-axis'} "
              f"dist={np.hypot(a-dstar,b-dstar):.4f}", flush=True)

    # [3] per-CI confirmation: tight loop (each must carry PI) + conical (linear) dispersion
    perci = []
    for (a, b, c) in true:
        s = berry(a, b, 0.05, N=72)
        signs = berry_band(a, b, 0.04, N=200)
        eps = [0.01, 0.02, 0.04, 0.08]
        gl = [fgap(a + e / np.sqrt(2), b + e / np.sqrt(2)) for e in eps]
        slopes = [gl[k] / eps[k] for k in range(4)]
        perci.append(dict(d1=a, d2=b, tight_berry=(None if s is None else int(s)),
                          crossing_bands=[k for k in range(6) if signs[k] < 0],
                          conical_slope_mean=float(np.mean(slopes)), gap_ray=gl))
    out["per_ci"] = perci
    print("[3] per-CI: tight-loop Berry + crossing-band + conical slope", flush=True)
    for p in perci:
        print(f"    ({p['d1']:.5f},{p['d2']:.5f}) tight-Berry="
              f"{'PI' if p['tight_berry']<0 else '0'} crossing-bands={p['crossing_bands']} "
              f"conical-slope~{p['conical_slope_mean']:.3f}", flush=True)

    # [4] concentric loops about the central CI: fragile single-vec vs robust per-band,
    #     documenting the grazing artifact at r ~ off-axis-CI radius (0.31)
    radii = [0.12, 0.15, 0.22, 0.30, 0.38, 0.50, 0.60]
    conc = {}
    for r in radii:
        sv = berry(dstar, dstar, r, N=64)
        sb = berry_band(dstar, dstar, r, N=360)
        conc[str(r)] = dict(single_vec=(None if sv is None else int(sv)),
                            band2=int(sb[2]), band3=int(sb[3]))
    out["concentric_loops"] = conc
    out["interp_reference"] = {"0.15": "PI", "0.22": "PI", "0.30": "0", "0.38": "PI", "0.50": "PI"}
    print("[4] concentric loops about central CI (crossing band 2 is the physical phase):", flush=True)
    for r in radii:
        c = conc[str(r)]
        print(f"    r={r:<5} band2={c['band2']:+d} (={'PI' if c['band2']<0 else '0'})"
              f"  [single-vec {c['single_vec']:+d}]", flush=True)

    import json, os
    os.makedirs("debug/data", exist_ok=True)
    json.dump(out, open("debug/data/beh2_ci_exact_landscape.json", "w"), indent=2)
    print(f"\ntotal {time.time()-t0:.0f}s; wrote debug/data/beh2_ci_exact_landscape.json")
    return out


if __name__ == "__main__":
    import sys
    t0 = time.time()
    if "--novalidate" not in sys.argv:
        worst = _validate()
        print(f"validation done in {time.time()-t0:.1f}s (worst |err| vs mpmath = {worst:.1e})\n")
        assert worst < 1e-7, "fast evaluator not accurate enough"
        print("FAST EVALUATOR VALIDATED.  Proceeding.\n")
    run_landscape()
