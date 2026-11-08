"""
Z>20 cliff diagnostic — Probe (a): CR67 single-zeta fit faithfulness.

Hypothesis: The Clementi-Raimondi 1967 single-zeta exponents themselves are
non-faithful for heavy-atom outer shells, regardless of zeta count. The
hydrogenic R_nl(r; Z_eff = n*zeta_CR) has the wrong radial peak position
and wrong <r> compared to the physically-known orbital extent.

Test: For Cs 5p, Cs 6s (where the V.C.6 sprint surfaced the cliff), and
Sr 4d / 5s (intermediate Z) and Na 3s (anchor Z=11 where framework works),
compute the analytical hydrogenic R_nl(r), find the radial peak r_peak and
expectation <r>, and compare to:
  - Roughly-known physical <r> from atomic-physics tables (NIST,
    Roberts-Ginges-style effective Z reasoning).
  - The implied physical Z_eff (back-solve from <r> = n^2/Z_eff for
    hydrogenic).

Output a quantitative table of the fit faithfulness.

Read-only: uses geovac.neon_core._hydrogenic_radial as the production
implementation. NO modifications.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

# Add project root to path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from geovac.neon_core import (
    _hydrogenic_radial,
    _CLEMENTI_ZETA_NE,
    _CLEMENTI_ZETA_AR,
    _CLEMENTI_ZETA_KR,
    _CLEMENTI_ZETA_XE,
)


def radial_peak_and_expectation_value(R, r):
    """Find radial peak position and <r> from a hydrogenic R(r) on grid r.

    Returns
    -------
    r_peak : float    -- argmax of |R(r)|^2 r^2 (the true radial PDF)
    r_mean : float    -- <r> = integral |R|^2 r^3 dr (the expectation value)
    rms_r  : float    -- sqrt(<r^2>) for an additional peak-vs-extent check
    """
    pdf = R * R * r * r  # |R|^2 r^2 is the radial probability density
    # Normalize (in case of small numerical drift)
    norm = np.trapezoid(pdf, r)
    if norm <= 0:
        return float('nan'), float('nan'), float('nan')
    pdf = pdf / norm
    r_peak = float(r[np.argmax(pdf)])
    r_mean = float(np.trapezoid(pdf * r, r))
    r2_mean = float(np.trapezoid(pdf * r * r, r))
    return r_peak, r_mean, float(np.sqrt(r2_mean))


def hydrogenic_expected_r_mean(n: int, l: int, Z_eff: float) -> float:
    """Analytic hydrogenic <r>_{n,l} = (a_0 / 2 Z_eff) [3 n^2 - l(l+1)] in
    bohr (a_0 = 1)."""
    return (3.0 * n * n - l * (l + 1)) / (2.0 * Z_eff)


def hydrogenic_expected_peak(n: int, l: int, Z_eff: float) -> float:
    """For the n_r=0 case (l = n-1, e.g. 1s, 2p, 3d), the peak of |R|^2 r^2
    is at r = n^2 / Z_eff exactly. For higher l (multiple lobes), the peak
    is the OUTERMOST node — approximately n^2 / Z_eff."""
    return n * n / Z_eff


# Physical reference values for outer-shell <r> (sources cited in memo)
# Taken from atomic-physics handbooks / NIST / quasi-classical estimates;
# these are not multi-decimal accurate, but order-of-magnitude correct.
PHYSICAL_REF = {
    # (Z, n, l, label): (<r>_physical_bohr, source_note, approx_uncertainty_pct)
    (11, 3, 0, "Na 3s"): (3.43, "Slater rules + NIST: r_max(3s)~3.5", 15.0),
    (19, 4, 0, "K 4s"): (4.21, "Slater rules: r_max(4s)~4-5", 15.0),
    (38, 5, 0, "Sr 5s"): (4.56, "Slater rules + RHF tables: r_max(5s)~4.5-5",
                          20.0),
    (38, 4, 2, "Sr 4d"): (1.7, "RHF: 4d compact, r_max~1.5-2", 25.0),
    (55, 6, 0, "Cs 6s"): (5.64, "Sobelman/RHF: r_max(6s)~5.5-6", 15.0),
    (55, 5, 1, "Cs 5p"): (1.95, "Sobelman/RHF: r_max(5p)~1.8-2.2", 20.0),
    (55, 5, 0, "Cs 5s"): (1.27, "Sobelman/RHF: r_max(5s)~1.2-1.4", 20.0),
    (55, 4, 2, "Cs 4d"): (1.05, "RHF: 4d compact, r_max~1.0-1.1", 25.0),
}


def get_cr_zeta(Z: int, n: int, l: int) -> float | None:
    """Look up Clementi-Raimondi single-zeta for shell (n,l) at atom Z."""
    # Map Z to which table holds the data.
    if Z in _CLEMENTI_ZETA_NE:
        zetas = _CLEMENTI_ZETA_NE[Z]
        # (zeta_1s, zeta_2s, zeta_2p)
        if (n, l) == (1, 0): return zetas[0]
        if (n, l) == (2, 0): return zetas[1]
        if (n, l) == (2, 1): return zetas[2]
        # Outer valence at Z=11..18: 3s, 3p — NOT in [Ne] table; this is
        # the valence orbital itself, not a core. Skip.
        return None
    if Z in _CLEMENTI_ZETA_AR:
        zetas = _CLEMENTI_ZETA_AR[Z]
        # (zeta_1s, zeta_2s, zeta_2p, zeta_3s, zeta_3p)
        labels = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1)]
        for k, lab in enumerate(labels):
            if lab == (n, l): return zetas[k]
        return None
    if Z in _CLEMENTI_ZETA_KR:
        zetas = _CLEMENTI_ZETA_KR[Z]
        # (z1s, z2s, z2p, z3s, z3p, z3d, z4s, z4p)
        labels = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1),
                  (3, 2), (4, 0), (4, 1)]
        for k, lab in enumerate(labels):
            if lab == (n, l): return zetas[k]
        # Z=37, 38: 5s is the valence (4d is empty for Sr ground state)
        return None
    if Z in _CLEMENTI_ZETA_XE:
        zetas = _CLEMENTI_ZETA_XE[Z]
        # (z1s, z2s, z2p, z3s, z3p, z3d, z4s, z4p, z4d, z5s, z5p)
        labels = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2),
                  (4, 0), (4, 1), (4, 2), (5, 0), (5, 1)]
        for k, lab in enumerate(labels):
            if lab == (n, l): return zetas[k]
        return None
    return None


def main():
    results = []

    print("=" * 78)
    print("Z>20 CLIFF DIAGNOSTIC — PROBE (a): CR67 fit faithfulness")
    print("=" * 78)
    print()
    print(f"{'System':18s} {'Z_eff_CR':>10s} {'r_peak_CR':>11s} "
          f"{'<r>_CR':>10s} {'<r>_phys':>10s} "
          f"{'rel_err_<r>%':>14s} {'Z_eff_back':>11s}")
    print("-" * 78)

    # Build a fine geometric grid that captures both inner and outer behavior.
    r = np.geomspace(1e-4, 100.0, 8000)

    for (Z, n, l, label), (r_phys, source, unc) in PHYSICAL_REF.items():
        zeta_cr = get_cr_zeta(Z, n, l)
        if zeta_cr is None:
            print(f"{label:18s}  CR67 zeta NOT in core table (valence orbital)")
            continue

        # CR67 hydrogenic Z_eff = n * zeta (the convention in neon_core.py
        # production code).
        Z_eff_cr = n * zeta_cr

        # Compute the analytical hydrogenic R_nl(r; Z_eff_CR).
        R = _hydrogenic_radial(n, l, Z_eff_cr, r)

        r_peak_cr, r_mean_cr, _ = radial_peak_and_expectation_value(R, r)

        # Analytical hydrogenic <r> for cross-check:
        r_mean_analytical = hydrogenic_expected_r_mean(n, l, Z_eff_cr)
        r_peak_analytical = hydrogenic_expected_peak(n, l, Z_eff_cr)
        # (Numerical <r>_cr should agree with analytical to ~1%.)

        # Back-solve for the implied PHYSICAL Z_eff from the physical <r>:
        Z_eff_back = (3.0 * n * n - l * (l + 1)) / (2.0 * r_phys)

        # Relative error on <r>:
        rel_err_pct = 100.0 * (r_mean_cr - r_phys) / r_phys

        print(f"{label:18s} {Z_eff_cr:10.3f} {r_peak_cr:11.3f} "
              f"{r_mean_cr:10.3f} {r_phys:10.3f} "
              f"{rel_err_pct:14.1f} {Z_eff_back:11.3f}")

        results.append({
            "Z": Z, "n": n, "l": l, "label": label,
            "zeta_CR": zeta_cr,
            "Z_eff_CR_n_times_zeta": Z_eff_cr,
            "r_peak_CR_numerical": r_peak_cr,
            "r_mean_CR_numerical": r_mean_cr,
            "r_peak_CR_analytical": r_peak_analytical,
            "r_mean_CR_analytical": r_mean_analytical,
            "r_phys_reference_bohr": r_phys,
            "phys_ref_uncertainty_pct": unc,
            "rel_err_<r>_pct": rel_err_pct,
            "Z_eff_implied_by_phys_<r>": Z_eff_back,
            "Z_eff_overshoot_factor": Z_eff_cr / Z_eff_back,
            "source_note": source,
        })

    print()
    print("Z_eff_CR = n * zeta_CR (the convention production code uses).")
    print("Z_eff_back = back-solved physical Z_eff from <r>_phys.")
    print("Overshoot factor = Z_eff_CR / Z_eff_back (>1 means CR67 too compact).")
    print()

    # Compute summary statistics
    overshoots = [r["Z_eff_overshoot_factor"] for r in results]
    avg_overshoot = np.mean(overshoots) if overshoots else float('nan')
    print(f"Mean overshoot factor: {avg_overshoot:.3f}x")
    print(f"Max overshoot factor:  {max(overshoots):.3f}x")
    print(f"Min overshoot factor:  {min(overshoots):.3f}x")
    print()

    # Categorize by Z range to expose the cliff:
    rows_low_z = [r for r in results if r["Z"] <= 20]
    rows_mid_z = [r for r in results if 20 < r["Z"] <= 38]
    rows_hi_z = [r for r in results if r["Z"] > 38]

    def _avg_over(rows):
        if not rows:
            return float('nan')
        return float(np.mean([r["Z_eff_overshoot_factor"] for r in rows]))

    print(f"Z<=20  mean overshoot:  {_avg_over(rows_low_z):.3f}x  "
          f"({len(rows_low_z)} rows)")
    print(f"20<Z<=38 mean overshoot: {_avg_over(rows_mid_z):.3f}x  "
          f"({len(rows_mid_z)} rows)")
    print(f"Z>38   mean overshoot:  {_avg_over(rows_hi_z):.3f}x  "
          f"({len(rows_hi_z)} rows)")
    print()

    # Interpretive verdict
    if rows_hi_z and _avg_over(rows_hi_z) > 1.5 * _avg_over(rows_low_z + rows_mid_z):
        verdict = ("CONFIRMED: CR67 fits systematically more compact at high Z. "
                   "CR67 single-zeta is non-faithful for outer shells of "
                   "heavy atoms.")
    elif _avg_over(rows_low_z + rows_mid_z + rows_hi_z) > 1.3:
        verdict = ("PARTIAL: CR67 fits are moderately compact across Z range; "
                   "high-Z is somewhat worse but the effect is gradual, not "
                   "cliff-shaped.")
    else:
        verdict = ("NOT CONFIRMED: CR67 fits are reasonably faithful across "
                   "the Z range tested; the cliff is NOT in the CR67 fit.")

    print(f"VERDICT (Probe a): {verdict}")

    # Save data
    out = {
        "probe": "a",
        "hypothesis": "CR67 single-zeta exponents themselves are non-faithful for outer shells at high Z",
        "results_per_orbital": results,
        "summary": {
            "mean_overshoot": avg_overshoot,
            "max_overshoot": max(overshoots) if overshoots else None,
            "min_overshoot": min(overshoots) if overshoots else None,
            "low_Z_mean_overshoot": _avg_over(rows_low_z),
            "mid_Z_mean_overshoot": _avg_over(rows_mid_z),
            "high_Z_mean_overshoot": _avg_over(rows_hi_z),
        },
        "verdict": verdict,
    }

    out_path = Path(ROOT) / "debug" / "data" / "z_cliff_probe_a.json"
    out_path.parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
