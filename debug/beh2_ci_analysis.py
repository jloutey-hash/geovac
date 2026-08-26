"""Stage 2: EXACT-overlap CI landscape analysis (run after evaluator validated).
Imports the validated fast evaluator from beh2_ci_exact_landscape.
"""
from __future__ import annotations
import json, time
import numpy as np
from scipy.optimize import minimize
import beh2_ci_exact_landscape as E   # validated fast exact evaluator + projectors_exact/fgap/berry

t0 = time.time()
out = {}

# --- speed check ---
tb = time.time(); _ = E.S_block_fast(2.445, 2, 1); dt = time.time() - tb
print(f"S_block_fast timing: {dt*1000:.2f} ms/call", flush=True)

# ============================================================ 1. central CI on symmetric line
print("\n[1] symmetric-line gap scan (EXACT) over d in [1.5,3.5] ...", flush=True)
ds = np.linspace(1.5, 3.5, 201)
g = np.array([E.fgap(d, d) for d in ds])
# find ALL local minima on the line (to catch a 2nd on-axis CI if it exists)
loc = [i for i in range(1, len(ds)-1) if g[i] < g[i-1] and g[i] < g[i+1]]
loc_sorted = sorted(loc, key=lambda i: g[i])
print("  local minima on d1=d2 (gap, d):", [(round(float(g[i]),5), round(float(ds[i]),4)) for i in loc_sorted[:6]])
# refine the global one
i0 = int(np.argmin(g)); d0 = ds[i0]
r1 = minimize(lambda d: E.fgap(d[0], d[0]), [d0], method="Nelder-Mead",
              options=dict(xatol=1e-6, fatol=1e-12))
dstar = float(r1.x[0]); gstar = float(r1.fun)
print(f"  refined central CI: d* = {dstar:.6f}, exact F-gap = {gstar:.3e}", flush=True)
out["central_CI"] = dict(dstar=dstar, gap=gstar,
                         line_minima=[(float(g[i]), float(ds[i])) for i in loc_sorted[:6]])

# ============================================================ 2. EXACT Berry phase vs radius
print("\n[2] EXACT Berry phase around central CI at the flip radii ...", flush=True)
radii = [0.12, 0.15, 0.22, 0.30, 0.38, 0.50, 0.60]
berry_exact = {}
for r in radii:
    s = E.berry(dstar, dstar, r, N=48)
    tag = "PI" if (s is not None and s < 0) else ("0" if s is not None else "singular")
    berry_exact[r] = (None if s is None else int(s), tag)
    print(f"  r={r:<5}: sign={s}  => {tag}", flush=True)
# far control
sf = E.berry(4.8, 4.8, 0.4, N=48)
print(f"  FAR control @(4.8,4.8) r=0.4: sign={sf} => {'PI' if sf and sf<0 else '0' if sf else 'singular'}")
out["berry_exact_vs_radius"] = {str(k): v for k, v in berry_exact.items()}
out["berry_far_control"] = (None if sf is None else int(sf))
out["interp_reference"] = {"0.15":"PI","0.22":"PI","0.30":"0","0.38":"PI","0.50":"PI"}

# ============================================================ 3. EXACT 2D gap grid -> locate all CIs
print("\n[3] EXACT 2D F-gap grid over [2.0,2.9]^2 ...", flush=True)
box = np.linspace(2.0, 2.9, 46)   # step 0.02
GAP = np.full((46, 46), np.nan)
for a, d1 in enumerate(box):
    for b, d2 in enumerate(box):
        GAP[a, b] = E.fgap(d1, d2)
    if a % 10 == 0:
        print(f"  row {a}/46 (min so far {np.nanmin(GAP):.2e}) t={time.time()-t0:.0f}s", flush=True)
gmin = float(np.nanmin(GAP)); im = np.unravel_index(np.nanargmin(GAP), GAP.shape)
print(f"  grid min gap {gmin:.3e} at (d1,d2)=({box[im[0]]:.3f},{box[im[1]]:.3f})", flush=True)

# local minima of the 2D field (candidate CIs), then refine each with Nelder-Mead
def is_locmin(A, i, j):
    v = A[i, j]
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0: continue
            ii, jj = i+di, j+dj
            if 0 <= ii < A.shape[0] and 0 <= jj < A.shape[1] and A[ii, jj] < v:
                return False
    return True
cands = [(box[i], box[j], GAP[i, j]) for i in range(46) for j in range(46) if is_locmin(GAP, i, j)]
cands.sort(key=lambda t: t[2])
print(f"  {len(cands)} grid-local-minima; smallest gaps:")
refined = []
for (d1c, d2c, gc) in cands[:8]:
    rr = minimize(lambda x: E.fgap(x[0], x[1]), [d1c, d2c], method="Nelder-Mead",
                  options=dict(xatol=1e-6, fatol=1e-13, maxiter=400))
    d1r, d2r = float(rr.x[0]), float(rr.x[1]); gr = float(rr.fun)
    refined.append((d1r, d2r, gr, float(gc)))
    print(f"    grid({d1c:.3f},{d2c:.3f}) gap {gc:.2e} -> refined ({d1r:.5f},{d2r:.5f}) gap {gr:.3e}", flush=True)
out["grid_min"] = dict(gap=gmin, at=[float(box[im[0]]), float(box[im[1]])])
out["refined_CI_candidates"] = [dict(d1=a, d2=b, gap_refined=c, gap_grid=d) for (a,b,c,d) in refined]

# how many DISTINCT true CIs (gap below a degeneracy threshold, dedup by location)
THRESH = 5e-4
true_cis = []
for (a, b, c, d) in refined:
    if c < THRESH:
        if not any(abs(a-x)<1e-3 and abs(b-y)<1e-3 for (x,y,_) in true_cis):
            true_cis.append((a, b, c))
print(f"\n  TRUE CIs (refined gap < {THRESH}): {len(true_cis)}")
for (a, b, c) in true_cis:
    print(f"    ({a:.5f},{b:.5f})  gap={c:.2e}  {'[ON-AXIS]' if abs(a-b)<2e-3 else '[off-axis]'}")
out["true_CIs"] = [dict(d1=a, d2=b, gap=c, on_axis=bool(abs(a-b)<2e-3)) for (a,b,c) in true_cis]
out["true_CI_count"] = len(true_cis)

with open("debug/data/beh2_ci_exact_landscape.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"\ntotal {time.time()-t0:.0f}s; wrote debug/data/beh2_ci_exact_landscape.json")


# ============================================================ 4. per-CI confirmation
# (appended stage: tight loops around each located CI + conical linear-dispersion check)
def stage4():
    import numpy as np, json
    cis = [(2.445082, 2.445082, "central/on-axis"),
           (2.697930, 2.266390, "off-axis +"),
           (2.266390, 2.697930, "off-axis - (mirror)")]
    res = []
    print("\n[4] per-CI confirmation: tight Berry loop (r=0.05) + conical dispersion")
    for (a, b, lab) in cis:
        # refine locally first
        rr = minimize(lambda x: E.fgap(x[0], x[1]), [a, b], method="Nelder-Mead",
                      options=dict(xatol=1e-7, fatol=1e-14, maxiter=600))
        a, b = float(rr.x[0]), float(rr.x[1]); g0 = float(rr.fun)
        s = E.berry(a, b, 0.05, N=64)
        # conical (linear) dispersion: gap along a ray from the CI should be ~linear in distance
        eps = [0.01, 0.02, 0.04, 0.08]
        gline = [E.fgap(a + e/np.sqrt(2), b + e/np.sqrt(2)) for e in eps]
        # slope estimate (gap/eps) should be ~constant (conical), not ->0
        slopes = [gline[k]/eps[k] for k in range(len(eps))]
        tag = "PI" if (s is not None and s < 0) else ("0" if s is not None else "singular")
        print(f"  {lab:22s} ({a:.5f},{b:.5f}) gap0={g0:.2e}  tight-loop Berry={tag}")
        print(f"      gap along 45deg ray at eps={eps}: {[f'{x:.4f}' for x in gline]}  slope~{np.mean(slopes):.3f}")
        res.append(dict(label=lab, d1=a, d2=b, gap0=g0, tight_berry=(None if s is None else int(s)),
                        conical_slopes=slopes, gap_ray=gline))
    return res

if __name__ == "__main__" and False:
    pass
