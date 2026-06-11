"""Deeper Dirac-S^3 probe: look at the γ_n vs Dirac sequence without multiplicity stacking.

The Dirac-on-S^3 has a *degenerate* spectrum ({n+3/2} with multiplicity 2(n+1)(n+2)).
The compare_to_dirac routine stacks multiplicities and truncates at dim H, which
means the "Dirac reference spectrum" is dominated by the first 1-2 shells and
therefore compressed. Let's also test comparison against the *distinct* Dirac
eigenvalue sequence {3/2, 5/2, 7/2, ...} without multiplicity, and against
the Riemann zero imaginary parts t_n = Im(rho_n_zeta) as a sanity check.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from geovac.hp_operator import load_rhm_zeros


def main():
    data_path = ROOT / "debug" / "data" / "spectral_zero_stats.json"
    for which in ("D_full", "D_even", "D_odd"):
        ims = load_rhm_zeros(str(data_path), which=which)
        n = len(ims)
        print(f"\n=== {which}: n={n} ===")

        # Sequence A: Dirac distinct eigenvalues {n+3/2}
        dirac_distinct = np.array([k + 1.5 for k in range(n)], dtype=float)

        # Sequence B: Dirac stacked with multiplicity
        dirac_stacked = []
        nn = 0
        while len(dirac_stacked) < n:
            lam = nn + 1.5
            g = 2 * (nn + 1) * (nn + 2)
            dirac_stacked.extend([lam] * g)
            nn += 1
        dirac_stacked = np.array(sorted(dirac_stacked)[:n], dtype=float)

        # Sequence C: integers (dumb reference)
        integers = np.arange(1, n + 1, dtype=float)

        # Sequence D: sqrt(integers) (a common RMT heuristic)
        sqrt_n = np.sqrt(np.arange(1, n + 1, dtype=float)) * np.mean(ims) / np.mean(np.sqrt(np.arange(1, n + 1, dtype=float)))

        # Sequence E: known Riemann zeta zero heights (hard-coded first ~22)
        riemann_first = [
            14.134725141734693790,
            21.022039638771554993,
            25.010857580145688764,
            30.424876125859513210,
            32.935061587739189691,
            37.586178158825671257,
            40.918719012147495187,
            43.327073280914999519,
            48.005150881167159728,
            49.773832477672302182,
            52.970321477714460644,
            56.446247697063394803,
            59.347044002602353079,
            60.831778524609809844,
            65.112544048081606660,
            67.079810529494173714,
            69.546401711173979253,
            72.067157674481907582,
            75.704690699083933168,
            77.144840068874800482,
            79.337375020249367922,
            82.910380854086030183,
        ]
        riemann_first = np.array(riemann_first[:n], dtype=float)

        for name, ref in (
            ("Dirac distinct {n+3/2}", dirac_distinct),
            ("Dirac stacked w/ mult", dirac_stacked),
            ("n (integer)", integers),
            ("sqrt(n) scaled", sqrt_n),
            ("Riemann zeros (first n)", riemann_first),
        ):
            # rescale ims to ref range, report residual rms after linear fit
            if ref.std() == 0 or ims.std() == 0:
                continue
            # sort both
            a = np.sort(ims)
            b = np.sort(ref)[:len(a)]
            # pearson r
            pr = float(np.corrcoef(a, b)[0, 1])
            # linear fit deviation
            slope, intercept = np.polyfit(b, a, 1)
            resid = a - (slope * b + intercept)
            rms = float(np.sqrt(np.mean(resid ** 2)))
            nlrms = rms / a.std()
            # ratio a / b
            print(f"  {name:30s}  pearson r = {pr: .4f}  resid rms = {rms:7.3f}  rel nonlin = {nlrms:.4f}")

        # Also: what does ims[n]/n look like? Is it linear in n (Weyl), sqrt(n), n log n?
        ns = np.arange(1, n + 1)
        ratios = np.sort(ims) / ns
        print(f"  ims/n: first = {ratios[0]:.2f}, last = {ratios[-1]:.2f}, mean = {ratios.mean():.2f}")
        # ims/sqrt(n)
        sqrt_ratios = np.sort(ims) / np.sqrt(ns)
        print(f"  ims/sqrt(n): first = {sqrt_ratios[0]:.2f}, last = {sqrt_ratios[-1]:.2f}, mean = {sqrt_ratios.mean():.2f}")
        # ims/(n log n)  — Weyl asymptotic for Riemann zeros
        nlogn = np.where(ns >= 2, ns * np.log(np.maximum(ns, 2)), ns)
        wlog_ratios = np.sort(ims) / nlogn
        print(f"  ims/(n log n): first = {wlog_ratios[0]:.2f}, last = {wlog_ratios[-1]:.2f}, mean = {wlog_ratios.mean():.2f}")


if __name__ == "__main__":
    main()
