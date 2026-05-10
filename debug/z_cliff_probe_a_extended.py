"""
Probe (a) extended: compare CR67 fits ACROSS shells of Cs (Z=55).

Test whether the CR67 cliff is restricted to outer shells (5p, 5s, 4d) or
extends to inner shells too. If inner shells (1s, 2s, 2p) are well-fit but
outer shells systematically fail, this confirms (a) and pins it cleanly to
the outer-shell regime where the screening transition is steepest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from geovac.neon_core import (
    _hydrogenic_radial,
    _CLEMENTI_ZETA_XE,
)


def radial_pdf_stats(R, r):
    """<r>, r_peak, r_rms from R(r) on grid r."""
    pdf = R * R * r * r
    norm = np.trapezoid(pdf, r)
    if norm <= 0:
        return float('nan'), float('nan'), float('nan')
    pdf = pdf / norm
    r_peak = float(r[np.argmax(pdf)])
    r_mean = float(np.trapezoid(pdf * r, r))
    r2 = float(np.trapezoid(pdf * r * r, r))
    return r_peak, r_mean, float(np.sqrt(r2))


# Physical reference <r> values for Cs INNER and OUTER shells.
# Source: Sobelman "Atomic Spectra and Radiative Transitions" + RHF tables;
# uncertainties widen for outer shells.
PHYS_REF_CS = {
    # (n, l, label): (<r>_phys_bohr, uncertainty_pct)
    (1, 0, "Cs 1s"): (0.0185, 5.0),     # ~ 1/Z = 1/55 ~ 0.018
    (2, 0, "Cs 2s"): (0.0762, 8.0),     # ~ 4/(Z-2) ~ 0.076
    (2, 1, "Cs 2p"): (0.0688, 8.0),     # ~ 0.07
    (3, 0, "Cs 3s"): (0.190, 10.0),     # screened by ~10
    (3, 1, "Cs 3p"): (0.184, 10.0),
    (3, 2, "Cs 3d"): (0.230, 12.0),     # 3d sees Z_eff ~ 25-30
    (4, 0, "Cs 4s"): (0.420, 15.0),
    (4, 1, "Cs 4p"): (0.443, 15.0),
    (4, 2, "Cs 4d"): (1.05, 25.0),      # 4d is the kicker — penetrates inner
    (5, 0, "Cs 5s"): (1.27, 20.0),
    (5, 1, "Cs 5p"): (1.95, 20.0),
}


def main():
    Z = 55
    zetas = _CLEMENTI_ZETA_XE[Z]
    # (z1s, z2s, z2p, z3s, z3p, z3d, z4s, z4p, z4d, z5s, z5p)
    shell_zeta_map = {
        (1, 0): zetas[0], (2, 0): zetas[1], (2, 1): zetas[2],
        (3, 0): zetas[3], (3, 1): zetas[4], (3, 2): zetas[5],
        (4, 0): zetas[6], (4, 1): zetas[7], (4, 2): zetas[8],
        (5, 0): zetas[9], (5, 1): zetas[10],
    }

    r = np.geomspace(1e-5, 100.0, 10000)

    print("=" * 86)
    print("Probe (a) extended — Cs Z=55 inner-vs-outer shell CR67 fit faithfulness")
    print("=" * 86)
    print(f"{'shell':10s} {'zeta_CR':>9s} {'Zeff_n*z':>10s} "
          f"{'r_peak_CR':>11s} {'<r>_CR':>10s} {'<r>_phys':>10s} "
          f"{'rel_err%':>10s} {'overshoot':>11s}")
    print("-" * 86)

    rows = []
    for (n, l, label), (r_phys, unc) in PHYS_REF_CS.items():
        zeta_cr = shell_zeta_map[(n, l)]
        Z_eff = n * zeta_cr
        R = _hydrogenic_radial(n, l, Z_eff, r)
        r_peak, r_mean, _ = radial_pdf_stats(R, r)

        rel_err = 100.0 * (r_mean - r_phys) / r_phys
        Z_eff_back = (3.0 * n * n - l * (l + 1)) / (2.0 * r_phys)
        overshoot = Z_eff / Z_eff_back

        print(f"{label:10s} {zeta_cr:9.3f} {Z_eff:10.2f} "
              f"{r_peak:11.4f} {r_mean:10.4f} {r_phys:10.4f} "
              f"{rel_err:10.1f} {overshoot:11.3f}")

        rows.append({
            "n": n, "l": l, "label": label,
            "zeta_CR": zeta_cr,
            "Z_eff_CR": Z_eff,
            "r_peak_CR": r_peak,
            "r_mean_CR": r_mean,
            "r_phys": r_phys,
            "rel_err_pct": rel_err,
            "Z_eff_back": Z_eff_back,
            "overshoot": overshoot,
            "uncertainty_pct": unc,
        })

    print()
    print("Group analysis:")
    inner = [r for r in rows if r["n"] <= 2]
    middle = [r for r in rows if r["n"] == 3]
    outer = [r for r in rows if r["n"] >= 4]

    def _stats(group, label):
        if not group:
            return
        oss = [r["overshoot"] for r in group]
        errs = [r["rel_err_pct"] for r in group]
        print(f"  {label:20s} N={len(group):2d}  "
              f"mean_overshoot={np.mean(oss):.3f}  "
              f"mean_rel_err_<r>%={np.mean(errs):+.1f}  "
              f"range=[{min(oss):.2f}, {max(oss):.2f}]")

    _stats(inner, "Inner (n=1,2):")
    _stats(middle, "Middle (n=3):")
    _stats(outer, "Outer (n=4,5):")

    print()

    # Verdict: cliff or gradual?
    inner_overshoot = np.mean([r["overshoot"] for r in inner]) if inner else float('nan')
    middle_overshoot = np.mean([r["overshoot"] for r in middle]) if middle else float('nan')
    outer_overshoot = np.mean([r["overshoot"] for r in outer]) if outer else float('nan')

    print("Cliff or gradual?")
    print(f"  Inner-to-outer overshoot ratio: {outer_overshoot/inner_overshoot:.2f}x")
    print(f"  Middle-to-outer overshoot ratio: {outer_overshoot/middle_overshoot:.2f}x")
    print()

    if outer_overshoot > 1.5 * inner_overshoot:
        verdict = (
            f"CONFIRMED CLIFF AT n>=4 IN Cs CR67. "
            f"Inner shells overshoot mean={inner_overshoot:.2f}x; "
            f"outer shells overshoot mean={outer_overshoot:.2f}x — "
            f"ratio {outer_overshoot/inner_overshoot:.2f}x. "
            f"CR67 single-zeta is non-faithful at the OUTER-SHELL level "
            f"of heavy atoms; inner shells are reasonably well-fit. "
            f"This sharpens the closeout-sprint diagnosis."
        )
    elif outer_overshoot > 1.2 * inner_overshoot:
        verdict = (
            f"PARTIAL outer-shell cliff. "
            f"Inner overshoot {inner_overshoot:.2f}x, "
            f"outer overshoot {outer_overshoot:.2f}x. "
            f"Effect present but gradual."
        )
    else:
        verdict = (
            f"NO CLIFF at outer-shell. CR67 overshoots roughly uniformly "
            f"({inner_overshoot:.2f}x inner, {outer_overshoot:.2f}x outer)."
        )

    print(f"VERDICT (Probe a-extended): {verdict}")

    out = {
        "probe": "a-extended",
        "atom": "Cs (Z=55)",
        "rows_per_shell": rows,
        "inner_overshoot_mean": float(inner_overshoot) if inner else None,
        "middle_overshoot_mean": float(middle_overshoot) if middle else None,
        "outer_overshoot_mean": float(outer_overshoot) if outer else None,
        "verdict": verdict,
    }

    out_path = Path(ROOT) / "debug" / "data" / "z_cliff_probe_a_extended.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
