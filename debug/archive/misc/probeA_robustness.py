"""Probe A -- leg 4 (robustness) + the subnormalisation check on the shift arm.

(a) Does the spectrum-aware exploit survive on the ASYMPTOTIC LAW alone?
    Build the degree-d polynomial from the sigma-law-PREDICTED eigenvalues
    (1 - sigma_k = (s^2/24)(k pi)^2/n^2, Paper 60 eq:sigma_law) and then measure its
    relative error on the TRUE (SVD) spectrum.

(b) Is the "shift" convention's sqrt(kappa) advantage robust to the block-encoding
    subnormalisation?  Shifting to [-1,1] needs a block-encoding of
    (2S - (lam_max+lam_min)I)/(lam_max-lam_min) whose LCU 1-norm is ~3, i.e. the
    spectrum actually lands in [-1/3, 1/3] unless it is re-amplified.  Arm "shift3"
    repeats the shift arm with the spectrum compressed by 1/3.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from geovac.sturmian_sigma_law import sw_cross_block, sigma_spectrum  # noqa: E402
import probeA_specaware_lp as core                                     # noqa: E402
from probeA_specaware_scan import bisect_degree, HEADROOM, N_GRID      # noqa: E402

core.BND_PER_DEG = 6
OUT = os.path.join(HERE, "data", "probeA_robust.json")


def spectrum(s, n):
    sig = sigma_spectrum(sw_cross_block(s, n))
    return np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))


def law_spectrum(s, n):
    """1 - sigma_k = (s^2/24) (k pi)^2 / n^2, k = 1..n, clipped to (0, 1]."""
    k = np.arange(1, n + 1)
    gap = np.clip((s ** 2 / 24.0) * (k * np.pi) ** 2 / n ** 2, 1e-12, 1.0)
    sig = 1.0 - gap
    return np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))


def eval_rel_err(coef, lam_eval, lam_ref, conv, headroom=HEADROOM, squeeze=1.0):
    lo, hi = float(lam_ref[0]), float(lam_ref[-1])
    if conv == "p60":
        x = lam_eval / hi
    else:
        x = (2 * lam_eval - (lo + hi)) / (hi - lo) / squeeze
    y = headroom * np.sqrt(lo / lam_eval)
    p = core.cheb_design(np.clip(x, -1, 1), len(coef) - 1) @ coef
    return float(np.max(np.abs(p - y) / y)), float(np.max(np.abs(x)))


def bisect_degree_b(x, y, eps, seed, dmax, bound):
    """bisect_degree with a non-unit QSP bound (head-room test)."""
    old = core.min_rel_error
    core.min_rel_error = lambda a, b, d, bound=bound: old(a, b, d, bound=bound)
    try:
        return bisect_degree(x, y, eps, seed, dmax)
    finally:
        core.min_rel_error = old


def main():
    out = []
    # ---------------- (a) law-vs-true robustness, aware/p60 arm ------------------
    for s in (1.4, 2.0, 3.0):
        for n in (6, 10, 14, 18):
            lt, ll = spectrum(s, n), law_spectrum(s, n)
            for eps in (1e-3, 1e-5):
                xt, yt, _ = core.make_problem(lt, "aware", "p60", headroom=HEADROOM)
                d, r, coef, _ = bisect_degree(xt, yt, eps, max(4, int(0.35 * lt[-1] / lt[0])), 900)
                if d < 0:
                    continue
                # same degree, but the polynomial is DESIGNED on the law spectrum
                # (accuracy imposed at law nodes, boundedness identical)
                xl = ll / lt[-1]
                yl = HEADROOM * np.sqrt(lt[0] / ll)
                r_law, _ = core.min_rel_error(xl, yl, d)
                err_true, _ = eval_rel_err(coef, lt, lt, "p60")
                # law-designed polynomial evaluated on the TRUE spectrum
                _, cl = core.min_rel_error(xl, yl, d)
                err_cross, _ = eval_rel_err(cl, lt, lt, "p60")
                out.append(dict(leg="law", s=s, n=n, eps=eps, d=d,
                                kappa=float(lt[-1] / lt[0]),
                                lam_min_true=float(lt[0]), lam_min_law=float(ll[0]),
                                err_true_design=err_true, err_law_design_on_true=err_cross,
                                r_law=float(r_law)))
                print(f"(a) s={s} n={n:2d} eps={eps:g} d={d:3d} "
                      f"lam_min true/law={lt[0]:.4e}/{ll[0]:.4e} "
                      f"err(true-designed)={err_true:.2e} "
                      f"err(law-designed on true)={err_cross:.2e}", flush=True)
                json.dump(out, open(OUT, "w"))

    # ---------------- (a2) BAND design: law-predicted eigenvalues + a safety band ---
    # The deployable version: require accuracy on a band [lam/(1+w), lam*(1+w)] around
    # each LAW-predicted eigenvalue (w covers the law's own error), never on the true
    # eigenvalues.  If d_band ~ d_aware the exploit is usable from the closed-form law.
    for s in (1.4, 2.0, 3.0):
        for n in (6, 10, 14, 18):
            lt, ll = spectrum(s, n), law_spectrum(s, n)
            lo, hi = lt[0], lt[-1]
            for w in (0.10, 0.25):
                band = []
                for lk in ll:
                    band.append(np.linspace(lk / (1 + w), min(lk * (1 + w), hi), 7))
                lb = np.clip(np.unique(np.concatenate(band)), lo * 0.999, hi)
                x = lb / hi
                y = HEADROOM * np.sqrt(lo / lb)
                d, r, coef, _ = bisect_degree(x, y, 1e-3, max(4 * n, int(0.3 * hi / lo)), 900)
                if d < 0:
                    continue
                err, _ = eval_rel_err(coef, lt, lt, "p60")
                out.append(dict(leg="band", s=s, n=n, w=w, d=int(d), eps=1e-3,
                                kappa=float(hi / lo), err_on_true=err))
                print(f"(a2) s={s} n={n:2d} w={w:.2f} d={d:3d} "
                      f"err_on_TRUE_spectrum={err:.2e}", flush=True)
                json.dump(out, open(OUT, "w"))

    # ---------------- (c) QSP head-room: bound 1.0 vs 0.9 ------------------------
    for s in (2.0,):
        for n in (6, 10, 14, 18):
            lt = spectrum(s, n)
            x, y, _ = core.make_problem(lt, "aware", "p60", headroom=HEADROOM)
            for b in (1.0, 0.9):
                d, r, coef, _ = bisect_degree_b(x, y, 1e-3, 4 * n, 900, b)
                out.append(dict(leg="qspbound", s=s, n=n, bound=b, d=int(d)))
                print(f"(c) s={s} n={n:2d} QSP bound={b} d={d}", flush=True)
                json.dump(out, open(OUT, "w"))

    # ---------------- (b) shift arm under LCU subnormalisation -------------------
    for s in (2.0,):
        for n in (4, 6, 8, 10):
            lam = spectrum(s, n)
            kap = float(lam[-1] / lam[0])
            for sq in (1.0, 3.0):
                lo, hi = lam[0], lam[-1]
                for mode in ("aware", "generic"):
                    if mode == "aware":
                        lm = lam
                    else:
                        lm = np.unique(np.concatenate([
                            lo * (hi / lo) ** np.linspace(0, 1, N_GRID),
                            0.5 * (lo + hi) - 0.5 * (hi - lo) *
                            np.cos(np.pi * np.arange(N_GRID) / (N_GRID - 1)), lam]))
                    x = (2 * lm - (lo + hi)) / (hi - lo) / sq
                    y = HEADROOM * np.sqrt(lo / lm)
                    d, r, coef, _ = bisect_degree(x, y, 1e-3,
                                                  max(4, int(2.0 * np.sqrt(kap))), 900)
                    out.append(dict(leg="squeeze", s=s, n=n, kappa=kap, squeeze=sq,
                                    mode=mode, eps=1e-3, d=int(d)))
                    print(f"(b) s={s} n={n:2d} kap={kap:7.2f} squeeze={sq} {mode:7s} "
                          f"d={d}", flush=True)
                    json.dump(out, open(OUT, "w"))


if __name__ == "__main__":
    main()
