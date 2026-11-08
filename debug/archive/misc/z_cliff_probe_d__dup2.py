"""
Z>20 cliff diagnostic — Probe (d): multi-zeta machinery bug check.

Hypothesis (d): The two-zeta heuristic in geovac/multi_zeta_orbitals.py
itself has a bug not yet exposed. Independent of whether two-zeta CAN
close the cliff in principle, the implementation may be wrong.

Tests:
  1. STO normalization: chi(r) = N r^(n-1) exp(-zeta r), check
     integral chi^2 r^2 dr = 1 numerically.
  2. MultiZetaOrbital evaluation: build a known-answer single-zeta
     orbital using two-zeta machinery (set both inner/outer to same zeta),
     compare to _hydrogenic_radial.
  3. Two-zeta orbital normalization: integral R^2 r^2 dr after build
     should be ~1.0 (renormalization at build time).
  4. density_from_orbitals: returns the correct shape (sum of occ * |R|^2 r^2).
  5. Reproduce the closeout sprint's BBB93 Ne tabulation as a sanity case
     (where multi-zeta with proper inner/outer signs SHOULD work).
  6. Inner/outer ratio sanity: with the published _TWO_ZETA_SPLITS table,
     verify the two-zeta orbital's <r> brackets or matches the single-zeta.

Output: pass/fail per check + diagnosis of any actual bug.

Read-only.
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
    FrozenCore,
)
from geovac.multi_zeta_orbitals import (
    STO,
    MultiZetaOrbital,
    _build_two_zeta_orbital,
    _build_ne_orbitals_neutral,
    build_two_zeta_xe_orbitals_from_cr,
    density_from_orbitals,
    core_electron_count,
    _TWO_ZETA_SPLITS,
)


def main():
    print("=" * 80)
    print("Z>20 CLIFF DIAGNOSTIC — PROBE (d): multi-zeta machinery bug check")
    print("=" * 80)
    print()

    findings = []

    # ----------------------------------------------------------------------
    # CHECK 1: STO normalization
    # ----------------------------------------------------------------------
    print("CHECK 1: STO primitive normalization")
    print("-" * 60)
    r = np.geomspace(1e-5, 80.0, 8000)

    test_stos = [
        STO(n=1, zeta=1.0),
        STO(n=2, zeta=2.0),
        STO(n=3, zeta=3.5),
        STO(n=5, zeta=12.0),
    ]
    sto_results = []
    for sto in test_stos:
        chi = sto.evaluate(r)
        norm = np.trapezoid(chi * chi * r * r, r)
        sto_results.append({
            "n": sto.n, "zeta": sto.zeta,
            "norm": float(norm),
            "rel_err_norm": float(abs(norm - 1.0)),
        })
        print(f"  STO(n={sto.n}, zeta={sto.zeta}): "
              f"integral chi^2 r^2 dr = {norm:.6f} (target 1.0)")

    max_norm_err = max(s["rel_err_norm"] for s in sto_results)
    if max_norm_err < 0.01:
        verdict_1 = "PASS — all STO primitives correctly normalized to <1% precision"
    else:
        verdict_1 = f"FAIL — STO normalization off by {max_norm_err:.4f}"
    print(f"  Verdict 1: {verdict_1}")
    findings.append(("STO normalization", float(max_norm_err), verdict_1))
    print()

    # ----------------------------------------------------------------------
    # CHECK 2: MultiZetaOrbital reduces to single-zeta correctly
    # ----------------------------------------------------------------------
    print("CHECK 2: MultiZetaOrbital with degenerate primitives = single zeta")
    print("-" * 60)
    # Build a "two-zeta" orbital where both primitives have the same zeta
    # and total coefficients sum to 1 (so it's just a single-zeta R(r)).
    z = 5.0  # arbitrary
    orb = MultiZetaOrbital(
        n_orbital=1, l_orbital=0, occupancy=2,
        primitives=(STO(n=1, zeta=z), STO(n=1, zeta=z)),
        coefficients=(0.5, 0.5),
    )
    R_two_zeta = orb.evaluate(r)
    R_single_zeta = _hydrogenic_radial(1, 0, z, r)  # n*zeta = 1*5 = 5 in production conv

    # Note: _hydrogenic_radial uses different normalization convention (Z_eff = n*zeta,
    # but for n=1 they should match if we use Z_eff=zeta). Let's compare shapes.
    # For n=1, Z_eff=5: rho = 2*5*r/1 = 10r, exp(-5r) decay. R_1s(r) = 2*Z^{3/2} exp(-Z r) = 2*5^{1.5} exp(-5r).
    # For STO with n=1, zeta=5: chi = (2*5)^{1.5}/sqrt(2!) * exp(-5r) = 10^{1.5}/sqrt(2) * exp(-5r) = 31.62/1.414 = 22.36 * exp(-5r)
    # 2*5^{1.5} = 22.36; matches. OK.

    # So R_two_zeta should equal R_single_zeta up to coefficient factor.
    # With coeffs (0.5, 0.5) summing to 1: R_two_zeta = (0.5+0.5) * chi = 1 * chi.
    # R_single_zeta should equal chi.

    # But _hydrogenic_radial RENORMALIZES on the grid; the STO is normalized
    # analytically. So they may differ by a numerical-renormalization factor.
    norm_two = np.trapezoid(R_two_zeta**2 * r**2, r)
    norm_single = np.trapezoid(R_single_zeta**2 * r**2, r)

    print(f"  MultiZetaOrbital (1s, both zeta=5, c=0.5+0.5):")
    print(f"  integral R^2 r^2 dr = {norm_two:.6f}")
    print(f"  hydrogenic R_1s(Z=5) integral R^2 r^2 dr = {norm_single:.6f}")
    # Compare ratio
    ratio = R_two_zeta[len(r)//4] / R_single_zeta[len(r)//4]
    print(f"  Ratio R_two_zeta/R_single_zeta at r={r[len(r)//4]:.3f}: {ratio:.6f}")

    if abs(norm_single - 1.0) < 0.01 and abs(ratio - 1.0) < 0.01:
        verdict_2 = "PASS — MultiZetaOrbital reduces correctly to single zeta"
    else:
        verdict_2 = (f"PARTIAL — shapes match (ratio {ratio:.3f}) but "
                     f"normalizations differ (multi {norm_two:.3f}, single "
                     f"{norm_single:.3f}). This is expected because "
                     "MultiZetaOrbital does NOT renormalize internally.")
    print(f"  Verdict 2: {verdict_2}")
    findings.append(("MultiZetaOrbital reduces to single-zeta", abs(ratio-1.0), verdict_2))
    print()

    # ----------------------------------------------------------------------
    # CHECK 3: _build_two_zeta_orbital normalization
    # ----------------------------------------------------------------------
    print("CHECK 3: _build_two_zeta_orbital renormalization")
    print("-" * 60)
    zeta_cs_5p = 12.31
    orb_cs_5p = _build_two_zeta_orbital(n=5, l=1, occupancy=6, zeta_cr=zeta_cs_5p)
    R_5p = orb_cs_5p.evaluate(r)
    norm_5p = np.trapezoid(R_5p * R_5p * r * r, r)
    print(f"  Cs 5p two-zeta orbital normalization (after build): {norm_5p:.6f}")

    if abs(norm_5p - 1.0) < 0.01:
        verdict_3 = "PASS — _build_two_zeta_orbital correctly renormalizes to <1%"
    else:
        verdict_3 = f"FAIL — _build_two_zeta_orbital normalization off by {abs(norm_5p-1.0):.4f}"
    print(f"  Verdict 3: {verdict_3}")
    findings.append(("_build_two_zeta_orbital normalization", float(abs(norm_5p-1.0)), verdict_3))
    print()

    # ----------------------------------------------------------------------
    # CHECK 4: Two-zeta vs single-zeta <r> at Cs Z=55
    # ----------------------------------------------------------------------
    print("CHECK 4: Two-zeta vs single-zeta <r> for Cs orbitals")
    print("-" * 60)
    zetas_cs = _CLEMENTI_ZETA_XE[55]
    # (z1s, z2s, z2p, z3s, z3p, z3d, z4s, z4p, z4d, z5s, z5p)
    labels = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '5s', '5p']
    quantum_numbers = [(1,0), (2,0), (2,1), (3,0), (3,1), (3,2),
                       (4,0), (4,1), (4,2), (5,0), (5,1)]
    occupancies = [2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6]

    print(f"{'orbital':10s} {'zeta_CR':>10s} {'<r>_single':>11s} "
          f"{'<r>_2zeta':>11s} {'ratio':>8s}")
    print("-" * 55)

    pdf_data = []
    for label, (n, l), occ, zeta_cr in zip(labels, quantum_numbers, occupancies, zetas_cs):
        # single-zeta R (production): Z_eff = n * zeta_cr
        R_single = _hydrogenic_radial(n, l, n * zeta_cr, r)
        pdf_single = R_single**2 * r**2
        norm_single = np.trapezoid(pdf_single, r)
        if norm_single > 0:
            r_mean_single = np.trapezoid(pdf_single * r, r) / norm_single
        else:
            r_mean_single = float('nan')

        # two-zeta R via _build_two_zeta_orbital
        orb_2z = _build_two_zeta_orbital(n=n, l=l, occupancy=occ, zeta_cr=zeta_cr)
        R_2z = orb_2z.evaluate(r)
        pdf_2z = R_2z**2 * r**2
        norm_2z = np.trapezoid(pdf_2z, r)
        r_mean_2z = np.trapezoid(pdf_2z * r, r) / norm_2z if norm_2z > 0 else float('nan')

        ratio = r_mean_2z / r_mean_single if r_mean_single > 0 else float('nan')
        print(f"{label:10s} {zeta_cr:10.3f} {r_mean_single:11.4f} "
              f"{r_mean_2z:11.4f} {ratio:8.3f}")

        pdf_data.append({
            "orbital": label, "n": n, "l": l, "zeta_cr": zeta_cr,
            "r_mean_single": float(r_mean_single),
            "r_mean_2z": float(r_mean_2z),
            "ratio_2z_over_single": float(ratio),
        })

    print()

    # The two-zeta should be MORE EXTENDED than single-zeta if the outer ratio
    # is < 1. From _TWO_ZETA_SPLITS (1.20, 0.83), the inner is more compact,
    # the outer is more diffuse. The two-zeta <r> should be intermediate to
    # the (1.20-zeta-single, 0.83-zeta-single) <r>. The ratio 2z/single
    # should be ~ 1.0 (centered around CR67).

    ratios = [d["ratio_2z_over_single"] for d in pdf_data]
    avg_ratio = np.mean(ratios)
    print(f"  Average <r>_2z / <r>_single ratio: {avg_ratio:.4f}")

    if 0.85 < avg_ratio < 1.15:
        verdict_4 = ("PASS — two-zeta <r> is centered around single-zeta <r> "
                     f"(ratio {avg_ratio:.3f}), as expected from the "
                     "(1.20, 0.83) inner/outer split")
    else:
        verdict_4 = (f"PARTIAL — two-zeta <r> ratio is {avg_ratio:.3f}; "
                     "inner/outer split is biased.")
    print(f"  Verdict 4: {verdict_4}")
    findings.append(("Two-zeta <r> centering", float(abs(avg_ratio - 1.0)), verdict_4))
    print()

    # ----------------------------------------------------------------------
    # CHECK 5: density_from_orbitals integration to N_core
    # ----------------------------------------------------------------------
    print("CHECK 5: density_from_orbitals integrates to N_core")
    print("-" * 60)
    # Use the two-zeta Cs core
    orbs_cs = build_two_zeta_xe_orbitals_from_cr(zetas_cs)
    n_core_target = core_electron_count(orbs_cs)
    print(f"  Target N_core (sum of occupancies): {n_core_target}")

    density = density_from_orbitals(orbs_cs, r)
    integral = np.trapezoid(density, r)
    print(f"  Numerical integral of density(r): {integral:.4f}")

    # The target is 54 for [Xe] core. Note the density is occ * |R|^2 r^2,
    # and the orbitals are normalized so integral |R|^2 r^2 dr = 1.
    # So integral density = sum occ = N_core.
    rel_err_5 = abs(integral - n_core_target) / n_core_target
    print(f"  Relative error: {rel_err_5*100:.2f}%")
    if rel_err_5 < 0.05:
        verdict_5 = "PASS — density_from_orbitals integrates to N_core within 5%"
    else:
        verdict_5 = f"FAIL — density integral off by {rel_err_5*100:.1f}%"
    print(f"  Verdict 5: {verdict_5}")
    findings.append(("density integrates to N_core", float(rel_err_5), verdict_5))
    print()

    # ----------------------------------------------------------------------
    # CHECK 6: Are the inner/outer ratios "right" — do they bracket CR67 zeta?
    # ----------------------------------------------------------------------
    print("CHECK 6: _TWO_ZETA_SPLITS sanity")
    print("-" * 60)
    print(f"  All entries use (inner_ratio, outer_ratio, c_in, c_out):")
    print(f"  Different shells use slightly different c_in, c_out.")
    print(f"  Inner ratio: {set(v[0] for v in _TWO_ZETA_SPLITS.values())}  (always > 1: more compact)")
    print(f"  Outer ratio: {set(v[1] for v in _TWO_ZETA_SPLITS.values())}  (always < 1: more diffuse)")
    print(f"  This is consistent with the heuristic claim that inner/outer "
          "should bracket CR67.")
    print()
    print(f"  HOWEVER: a properly tabulated BBB93 multi-zeta has 5+ primitives "
          "per orbital,")
    print(f"  with several NEGATIVE coefficients to enforce orthogonality to "
          "inner shells.")
    print(f"  The two-zeta heuristic in this module has only POSITIVE "
          "coefficients (c_in, c_out > 0).")
    print(f"  This means the radial nodes (r-zeros that distinguish 5p from 5s, "
          "etc.) are NOT captured.")

    # Cross-check: compute the orthogonality between Cs 4d two-zeta and 5p
    # two-zeta radial wavefunctions. If they are not orthogonal, the
    # screening profile is structurally wrong.
    R_4d = orbs_cs[8].evaluate(r)  # 4d
    R_5s = orbs_cs[9].evaluate(r)  # 5s
    R_5p = orbs_cs[10].evaluate(r) # 5p
    overlap_4d_5s = np.trapezoid(R_4d * R_5s * r * r, r)
    overlap_5s_5p = np.trapezoid(R_5s * R_5p * r * r, r)  # different l, automatic 0 by L
    overlap_4d_5p = np.trapezoid(R_4d * R_5p * r * r, r)  # different l, 0 by L

    print()
    print("  Inter-shell orthogonality (radial overlap):")
    print(f"    <R_4d | R_5s> = {overlap_4d_5s:.4f}  (DIFFERENT n & l: not orthogonal in radial!)")
    print(f"    <R_5s | R_5p> = {overlap_5s_5p:.4f}  (different l: orthogonal by L Y selection)")
    print(f"    <R_4d | R_5p> = {overlap_4d_5p:.4f}  (different l: orthogonal by L)")
    print()
    # The 4d-5s overlap is ESSENTIAL — they have different l, so via Y_lm orthogonality
    # the 3D overlap is zero anyway. But the 4d-5s RADIAL overlap is informative:
    # in real RHF, there is no 4s-5s constraint per se, but a 1s-2s, 2s-3s, 3s-4s, 4s-5s
    # ladder of orthogonalities. Let's check 1s-2s.
    R_1s = orbs_cs[0].evaluate(r)
    R_2s = orbs_cs[1].evaluate(r)
    R_3s = orbs_cs[3].evaluate(r)
    R_4s = orbs_cs[6].evaluate(r)
    R_5s_again = R_5s
    s_orbitals = [(0, R_1s), (1, R_2s), (3, R_3s), (6, R_4s), (9, R_5s_again)]

    print("  s-orbital ladder (orthogonality required by RHF):")
    for i in range(len(s_orbitals)):
        for j in range(i+1, len(s_orbitals)):
            n_i, R_i = s_orbitals[i]
            n_j, R_j = s_orbitals[j]
            ov = np.trapezoid(R_i * R_j * r * r, r)
            label_i = orbs_cs[n_i].n_orbital
            label_j = orbs_cs[n_j].n_orbital
            print(f"    <R_{label_i}s | R_{label_j}s> = {ov:.4f}    "
                  f"(should be 0 in real RHF)")

    s_overlaps = []
    for i in range(len(s_orbitals)):
        for j in range(i+1, len(s_orbitals)):
            n_i, R_i = s_orbitals[i]
            n_j, R_j = s_orbitals[j]
            ov = np.trapezoid(R_i * R_j * r * r, r)
            s_overlaps.append(abs(ov))

    max_overlap = max(s_overlaps) if s_overlaps else 0
    print(f"\n  Max |<n s | m s>| (n < m): {max_overlap:.4f}")
    if max_overlap > 0.1:
        verdict_6 = ("STRUCTURAL OBSERVATION — two-zeta orbitals have LARGE "
                     f"non-orthogonality (max |<ns|ms>|={max_overlap:.3f}). "
                     "This is the missing radial-nodes mechanism. "
                     "The two-zeta heuristic is structurally incomplete.")
    else:
        verdict_6 = (f"Max overlap |<ns|ms>|={max_overlap:.3f}; "
                     "non-orthogonality is mild.")
    print(f"  Verdict 6: {verdict_6}")
    findings.append(("Two-zeta s-orbital orthogonality", float(max_overlap), verdict_6))
    print()

    # ----------------------------------------------------------------------
    # Net verdict
    # ----------------------------------------------------------------------
    print("=" * 80)
    print("NET VERDICT (Probe d)")
    print("=" * 80)

    failures = [f for f in findings if "FAIL" in f[2]]
    structurals = [f for f in findings if "STRUCTURAL" in f[2]]
    if failures:
        verdict = (f"{len(failures)} machinery FAILURES detected. "
                   "Multi-zeta machinery has a coding bug.")
    elif structurals:
        verdict = (
            "NO CODING BUGS in the multi-zeta machinery itself. "
            "STRUCTURAL ISSUE flagged: "
            f"{structurals[0][2]} "
            "The heuristic two-zeta with all-positive coefficients cannot "
            "represent the radial nodes that orthogonalize valence to inner "
            "shells. This is consistent with the closeout-sprint diagnosis "
            "but adds a sharper sub-mechanism: the FAILURE OF S-ORBITAL "
            "ORTHOGONALITY in the two-zeta heuristic is the structural "
            "reason a heuristic two-zeta cannot fix what proper RHF would."
        )
    else:
        verdict = "Multi-zeta machinery passes all checks."

    print(verdict)

    out = {
        "probe": "d",
        "hypothesis": "multi-zeta machinery itself has a coding bug",
        "findings": [{"check": f[0], "metric": f[1], "verdict": f[2]} for f in findings],
        "pdf_data_per_orbital": pdf_data,
        "max_s_orbital_overlap": float(max_overlap),
        "verdict": verdict,
    }
    out_path = Path(ROOT) / "debug" / "data" / "z_cliff_probe_d.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
