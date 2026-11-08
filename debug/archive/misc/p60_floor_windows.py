"""JOB A: window-sensitivity of the fitted convergence floor.

gap(K) = c + b*K**-q  (3 free params) fitted on sliding and cumulative windows of
the (K, gap) ladder, for BOTH the ground root and the 2^1S root.  A 3-parameter
free-floor fit can manufacture a floor out of slow convergence, so the deliverable
is the DRIFT of c across windows, not a single number.

Cross-checks per window:
  - pure power law  gap = b*K**-q  (c == 0), for RMS-residual comparison
  - 3-point Richardson/Shanks limit on the last three gaps (independent route)
"""
import json, os, sys
import numpy as np

CHEM = 1.5936014616
lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
rows = sorted(json.load(open("debug/data/p60_excited_l%d.json" % lmax)), key=lambda r: r["K"])
K = np.array([r["K"] for r in rows], float)
G = {0: np.array([r["gap0"] for r in rows]), 1: np.array([r["gap1"] for r in rows])}


def fit_floor(k, g, qgrid=np.linspace(0.05, 8.0, 15991)):
    """min over q of ||g - (c + b*K^-q)||, (c,b) by linear LS at each q."""
    best = None
    for q in qgrid:
        A = np.column_stack([np.ones_like(k), k ** (-q)])
        try:
            sol, *_ = np.linalg.lstsq(A, g, rcond=None)
        except np.linalg.LinAlgError:
            continue
        res = g - A @ sol
        r = float(np.sqrt(np.mean(res ** 2)))
        if best is None or r < best[0]:
            best = (r, float(sol[0]), float(sol[1]), float(q))
    return dict(rms=best[0], c=best[1], b=best[2], q=best[3])


def fit_pure(k, g, qgrid=np.linspace(0.05, 8.0, 15991)):
    """c forced to 0: g = b*K^-q."""
    best = None
    for q in qgrid:
        a = k ** (-q)
        b = float(a @ g / (a @ a))
        r = float(np.sqrt(np.mean((g - b * a) ** 2)))
        if best is None or r < best[0]:
            best = (r, b, float(q))
    return dict(rms=best[0], c=0.0, b=best[1], q=best[2])


def shanks(g):
    """Aitken delta^2 on the last three values (independent limit estimate)."""
    a, b, c = g[-3], g[-2], g[-1]
    d = c - 2 * b + a
    return float(c - (c - b) ** 2 / d) if abs(d) > 1e-30 else float("nan")


out = {"lmax": lmax, "K": K.tolist(), "roots": {}}
for root in (0, 1):
    g = G[root]
    label = "ground (1^1S)" if root == 0 else "2^1S"
    print("=" * 78)
    print("ROOT %d  %s   gaps (mHa): %s" % (root, label, np.array2string(g, precision=4)))
    print("=" * 78)
    entries = []
    for wname, sl in (
        [("slide4[%d:%d]" % (i, i + 4), slice(i, i + 4)) for i in range(len(K) - 3)]
        + [("slide5[%d:%d]" % (i, i + 5), slice(i, i + 5)) for i in range(len(K) - 4)]
        + [("cumul[%d:]" % i, slice(i, len(K))) for i in range(len(K) - 3)]
    ):
        kk, gg = K[sl], g[sl]
        f = fit_floor(kk, gg)
        p = fit_pure(kk, gg)
        sh = shanks(gg) if len(gg) >= 3 else float("nan")
        e = dict(window=wname, Kmin=float(kk[0]), Kmax=float(kk[-1]), n=len(kk),
                 c=f["c"], b=f["b"], q=f["q"], rms=f["rms"],
                 pure_q=p["q"], pure_rms=p["rms"], shanks=sh,
                 c_over_chem=f["c"] / CHEM)
        entries.append(e)
        print("  %-14s K=%4d..%-4d n=%d | c=%+9.4f mHa (%.2fx chem) q=%5.3f b=%11.3e rms=%.2e "
              "| pure-law rms=%.2e (q=%.3f) | Shanks=%+8.4f"
              % (wname, kk[0], kk[-1], len(kk), e["c"], e["c_over_chem"], e["q"], e["b"],
                 e["rms"], e["pure_rms"], e["pure_q"], sh))
    out["roots"][str(root)] = entries

    for fam, tag in (("slide4", "4-pt sliding"), ("slide5", "5-pt sliding"), ("cumul", "cumulative")):
        cs = [e["c"] for e in entries if e["window"].startswith(fam)]
        if len(cs) < 2:
            continue
        lo, hi = min(cs), max(cs)
        span = hi - lo
        ref = max(abs(np.mean(cs)), 1e-12)
        trend = "DOWN" if cs[-1] < cs[0] else ("UP" if cs[-1] > cs[0] else "flat")
        print("  --> %-13s c: first=%+.4f last=%+.4f  range [%+.4f, %+.4f]  span=%.4f mHa "
              "= %.1f%% of mean  trend %s"
              % (tag, cs[0], cs[-1], lo, hi, span, 100 * span / ref, trend))
    print()

json.dump(out, open("debug/data/p60_floor_windows_l%d.json" % lmax, "w"), indent=1)
print("wrote debug/data/p60_floor_windows_l%d.json" % lmax)
