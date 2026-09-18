"""Step 2 of the #3 wiring sprint: validate `build_one_body_direct` ABOVE the
truncation ceiling its debug predecessor was validated at.

WHY.  `debug/direct_onebody_engine.build_direct_full` was validated at
(2,2,1)/(2,2,2)/(3,3,2) only -- but the headline point is (5,5)+delta.  That is
exactly the v5.13.6 failure shape: a validation ceiling sitting below the regime
that breaks (there, high (m,s) + large l_neumann in the closed-form B-seeds).  So
before the direct engine is wired in as the default, it is checked at (4,4,2) and
(5,5,2) against the definition of correctness, `one_body_mp` + `_factored_cob`.

SECOND reason this is not a re-run of the old check: the promoted engine takes a
`basis` argument.  The debug engine hardcoded the mu-adapted (gegenbauer) family;
`recondition_energy`'s DEFAULT is `laguerre_legendre`, which the direct build has
therefore NEVER been validated on.  Every truncation is checked in both families.

Usage:
    python debug/direct_onebody_validate.py                    # fast tier
    python debug/direct_onebody_validate.py 5 5 2 gegenbauer   # one point
"""
from __future__ import annotations

import sys
import time
from typing import Tuple

import mpmath as mp
import numpy as np

from geovac import prolate_recondition as pr

FAST_TIER = [(2, 2, 1), (2, 2, 2), (3, 3, 2)]
BIG_TIER = [(4, 4, 2), (5, 5, 2)]
BASES = ("laguerre_legendre", "gegenbauer")


def ground_truth(j_max: int, l_max: int, mu_max: int, alpha: float,
                 basis: str, R: float = pr.R_DEFAULT
                 ) -> Tuple[np.ndarray, np.ndarray]:
    """The definition: mpf monomial one-body, re-based by the factored change of
    basis, then downcast.  This is what `recondition_energy` uses today."""
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        n_mom = 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20
        A = pr.ngm._mono_moments(2.0 * alpha, n_mom)
        S, H1 = pr.one_body_mp(fns, alpha, R, A)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        Tr, Ta = pr._transforms_per_mu(basis, j_max, l_max, mu_max, alpha)
        S_o = pr._factored_cob(S, mu_max + 1, Nr, Na, Tr, Ta)
        H_o = pr._factored_cob(H1, mu_max + 1, Nr, Na, Tr, Ta)
        return pr._to_f64(S_o), pr._to_f64(H_o)


def metrics(A: np.ndarray, B: np.ndarray, floor: float = 1e-12) -> dict:
    """Both deviation metrics, plus enough to adjudicate an entry-relative spike.

    WHY TWO.  Entry-relative deviation (max |dX| / |X| over entries above
    `floor`) is NOT a soundness criterion for these matrices: H1's entries span
    many orders of magnitude and some pass near zero by cancellation between T
    and V_ne, so a tiny ABSOLUTE deviation divided by a near-zero denominator
    produces a large relative number with nothing wrong.  Measured directly:
    sweeping R at (2,2,1) the entry-relative deviation spikes to 6.7e-13 at
    R = 2.4 -- where the worst entry is |H| = 1.6e-3 against max|H| = 3.4e4,
    with an absolute deviation of 1.1e-15 -- and returns to 8.5e-16 at R = 3.2,
    while SCALE-relative (max |dX| / max|X|) stays flat at ~2e-16 at every R.

    So scale-relative is the primary criterion and entry-relative is a
    secondary diagnostic.  The `worst_*` fields let a spike be adjudicated
    rather than guessed at: a spike whose worst entry is orders below the matrix
    scale, with scale-relative in band, is cancellation and not a defect.
    """
    d = np.abs(A - B)
    m = np.abs(B) > floor
    out = {"scale": float(d.max() / np.abs(B).max()), "maxabs": float(np.abs(B).max())}
    if not m.any():
        out.update(entry=0.0, worst_idx=None, worst_val=0.0, worst_dev=0.0)
        return out
    rel = d[m] / np.abs(B[m])
    k = int(np.argmax(rel))
    flat = np.argmax(np.where(m, d / np.where(m, np.abs(B), 1.0), -1.0))
    out.update(entry=float(rel.max()),
               worst_idx=tuple(int(v) for v in np.unravel_index(flat, B.shape)),
               worst_val=float(np.abs(B[m])[k]), worst_dev=float(d[m][k]))
    return out


def check(j_max: int, l_max: int, mu_max: int, basis: str,
          alpha: float = 1.0) -> None:
    t0 = time.time()
    Sd, Hd = pr.build_one_body_direct(j_max, l_max, mu_max, alpha, basis=basis)
    t_direct = time.time() - t0
    t1 = time.time()
    Sg, Hg = ground_truth(j_max, l_max, mu_max, alpha, basis)
    t_mpf = time.time() - t1
    mS, mH = metrics(Sd, Sg), metrics(Hd, Hg)
    speed = t_mpf / t_direct if t_direct > 0 else float("inf")
    print(f"  ({j_max},{l_max}) mu<={mu_max:d} {basis:18s} N={Sd.shape[0]:5d}  "
          f"S: entry={mS['entry']:.2e} scale={mS['scale']:.2e}  |  "
          f"H: entry={mH['entry']:.2e} scale={mH['scale']:.2e}  "
          f"[direct {t_direct:6.2f}s vs mpf {t_mpf:7.1f}s = {speed:7.1f}x]",
          flush=True)
    # Adjudication data for the entry-relative number, so an identical max
    # across truncations can be checked rather than assumed: the orthogonal 1D
    # blocks for indices <= min(j_max) are the SAME polynomials at every larger
    # truncation, so a worst entry inside the shared sub-block makes the max
    # saturate and become truncation-independent.
    print(f"      worst H entry at {mH['worst_idx']}: |H|={mH['worst_val']:.3e}  "
          f"|dH|={mH['worst_dev']:.2e}  max|H|={mH['maxabs']:.3e}  "
          f"(ratio to scale {mH['worst_val'] / mH['maxabs']:.2e})", flush=True)
    # PRIMARY bar: scale-relative.  Reference band measured at (2,2,1)/(2,2,2)/
    # (3,3,2) in both families: S 2.5e-17..5.2e-16, H 4.7e-17..4.2e-16.
    if max(mS["scale"], mH["scale"]) > 1e-13:
        print("    ^^ SCALE-relative ABOVE band -- do NOT wire this in; "
              "diagnose first.", flush=True)
    elif max(mS["entry"], mH["entry"]) > 1e-10:
        print("    ^^ entry-relative high but scale-relative in band: check the "
              "worst-entry line above for cancellation before alarming.",
              flush=True)


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        b = sys.argv[4] if len(sys.argv) >= 5 else "gegenbauer"
        print(f"=== direct one-body vs mpf re-based ground truth ===", flush=True)
        check(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), b)
    else:
        print("=== FAST TIER: re-validate the port, both families ===", flush=True)
        print("(gegenbauer reproduces the debug engine's recorded 7e-16 / 9e-13;\n"
              " laguerre_legendre is a NEW code path -- recondition's default)\n",
              flush=True)
        for (j, l, mu) in FAST_TIER:
            for basis in BASES:
                check(j, l, mu, basis)
        print("\n=== BIG TIER runs separately (mpf ground truth is slow):",
              flush=True)
        for (j, l, mu) in BIG_TIER:
            for basis in BASES:
                print(f"    python debug/direct_onebody_validate.py "
                      f"{j} {l} {mu} {basis}", flush=True)
