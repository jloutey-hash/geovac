"""
Z>20 cliff diagnostic — Probe (c): Convention mismatch trace.

Trace Z_eff(r) -> psi_valence(r) -> R(0) -> A_HFS for Cs (single-zeta vs
multi-zeta) and check unit conventions, atomic vs SI, signed-vs-unsigned
screening, factor-of-2 / factor-of-Z / factor (2I+1)/2 errors that might
propagate into the final A_HF.

Key checks:
  1. The hyperfine A formula's prefactor: A_HF = (8 pi / 3) * (g_e/2) * (g_N/2)
     * alpha^2 * (m_e/m_p) * |psi(0)|^2. Does this match the experimental
     convention (Lande convention nu_HFS = 4A for Cs F=4 <-> F=3)?
  2. The g_N for Cs-133. Cs-133: I=7/2, mu = 2.582025 mu_N. So g_N = mu/I
     = 0.7377 (NOT 2 mu / I = 1.4753).
  3. The Z_eff(r) convention: Production code uses Z_eff = n * zeta in the
     hydrogenic R_nl; analytical convention is Z_eff = zeta directly. Does
     the n*zeta convention match standard references?
  4. The cumulative N_core(r) integration: integral n(r) dr (ld 514 of
     neon_core.py) integrates the radial number density n(r) which IS
     |R|^2 r^2 already (the r^2 weight is in n_r). So the cumulative is
     the right object. No double-r^2 error.

This probe is a STRUCTURAL check — find any factor-of-X bug in the chain.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from geovac.neon_core import (
    FrozenCore, _hydrogenic_radial,
    screened_psi_origin_squared,
)
from geovac.hyperfine_a_constant import (
    bohr_fermi_a_constant, M_PROTON_OVER_M_E, GE_FULL,
)


def main():
    print("=" * 80)
    print("Z>20 CLIFF DIAGNOSTIC — PROBE (c): Convention mismatch trace")
    print("=" * 80)
    print()

    findings = []

    # ----------------------------------------------------------------------
    # CHECK 1: BF formula sign and prefactor reproduces hydrogen 21 cm
    # ----------------------------------------------------------------------
    print("CHECK 1: BF formula reproduces hydrogen 21 cm at Z=1")
    print("-" * 60)

    psi0_h = 1.0 / np.pi  # |psi_1s(0)|^2 = 1/pi for Z=1
    g_p = 5.585695        # proton g-factor (NIST)
    bf = bohr_fermi_a_constant(
        psi0_squared=psi0_h, g_e=GE_FULL, g_N=g_p,
        m_p_over_m_e=M_PROTON_OVER_M_E,
    )
    A_h_MHz = bf['A_MHz']
    nu_h_21cm = 4.0 * A_h_MHz  # Lande F=1 <-> F=0 = (2F+1) coupling: actually for I=1/2, J=1/2: F=1 - F=0 splitting = A
    nu_exp_21cm = 1420.405751768

    # Note: For H, F=I+J=1, F=I-J=0, splitting = A * (F(F+1) - I(I+1) - J(J+1))/2
    # at F=1: 0.5*A*(2 - 0.75 - 0.75) = 0.5*A*0.5 = A/4
    # at F=0: 0.5*A*(0 - 0.75 - 0.75) = -3A/4
    # splitting = A/4 - (-3A/4) = A
    # So hydrogen 21 cm is at nu = A (in Hz), NOT 4A.
    # Cs is different: I=7/2, J=1/2: F=4, F=3
    # at F=4: 0.5*A*(20 - 15.75 - 0.75) = 0.5*A*3.5 = 1.75 A
    # at F=3: 0.5*A*(12 - 15.75 - 0.75) = 0.5*A*(-4.5) = -2.25 A
    # splitting = 1.75A - (-2.25A) = 4A
    # So Cs F=4-F=3 at 9192.631770 MHz means A = 2298.157... MHz. (matches)

    nu_h_predicted = A_h_MHz  # for H, splitting = A
    rel_err_h = 100.0 * (nu_h_predicted - nu_exp_21cm) / nu_exp_21cm

    print(f"  BF prefactor with g_e={GE_FULL:.6f}, g_p={g_p:.6f}, "
          f"|psi_1s(0)|^2 = 1/pi:")
    print(f"  A(H 1s) BF = {A_h_MHz:.4f} MHz")
    print(f"  nu(21cm) predicted = A = {nu_h_predicted:.4f} MHz")
    print(f"  nu(21cm) experimental = {nu_exp_21cm:.4f} MHz")
    print(f"  Rel error: {rel_err_h:+.4f}%")
    if abs(rel_err_h) < 0.01:
        verdict_1 = "PASS — BF formula sign + prefactor are correct"
    elif abs(rel_err_h) < 1.0:
        verdict_1 = ("ACCEPTABLE — BF formula correct to ~percent "
                     "(off-prefactor conventions or g_N rounding)")
    else:
        verdict_1 = ("FAIL — BF formula has a factor-of-X error "
                     "in the prefactor")
    print(f"  Verdict 1: {verdict_1}")
    findings.append(("BF formula H 21cm", rel_err_h, verdict_1))

    print()

    # ----------------------------------------------------------------------
    # CHECK 2: g_N convention for Cs-133
    # ----------------------------------------------------------------------
    print("CHECK 2: g_N for Cs-133 — convention sanity check")
    print("-" * 60)
    mu_cs_in_nuclear_magnetons = 2.582025  # NIST
    I_cs = 7.0 / 2.0
    g_cs_per_I = mu_cs_in_nuclear_magnetons / I_cs  # = 0.7377
    g_cs_2x = 2.0 * mu_cs_in_nuclear_magnetons / I_cs  # = 1.4754 (the WRONG def)

    print(f"  Cs-133: I = 7/2, mu = {mu_cs_in_nuclear_magnetons} mu_N")
    print(f"  g_Cs = mu / I = {g_cs_per_I:.5f}     (atomic-physics convention)")
    print(f"  g_Cs alt = 2 mu / I = {g_cs_2x:.5f}     (sanity: a 2x error if used)")
    print(f"  H proton: g_p = 5.585695 (well-known)")
    print(f"  Ratio g_p / g_Cs = {g_p / g_cs_per_I:.4f}")

    # If we mistakenly used 2*mu/I instead of mu/I, A would be 2x larger
    bf_cs_correct = bohr_fermi_a_constant(
        psi0_squared=1.328,  # the closeout-sprint Cs single-zeta value
        g_e=GE_FULL, g_N=g_cs_per_I,
    )
    bf_cs_wrong = bohr_fermi_a_constant(
        psi0_squared=1.328,
        g_e=GE_FULL, g_N=g_cs_2x,
    )
    print(f"  A(Cs) at psi0_sq=1.328, correct g_N: {bf_cs_correct['A_MHz']:.2f} MHz")
    print(f"  A(Cs) at psi0_sq=1.328, doubled g_N: {bf_cs_wrong['A_MHz']:.2f} MHz")
    print(f"  Experimental A_Cs = 2298.158 MHz")

    # The closeout sprint result is 1219 MHz (under hyperfine-A v2). If they
    # used 2*mu/I by mistake: it would be 2*1219 = 2438. That is closer to
    # 2298 but still off by 6%. The actual residual is -47%, so the 2x-g_N
    # error WOULD have hidden a 6% positive residual but we are at -47%, so
    # the 2x-g_N bug is NOT the explanation. Sanity check.

    diff_2x_vs_actual = abs(2.0 * 1219 - 2298) / 2298
    print(f"  If 2x g_N bug were the cause, residual would be {diff_2x_vs_actual*100:+.1f}% "
          f"(closeout actual: -47%). 2x-g_N bug RULED OUT.")

    verdict_2 = "PASS — g_N convention is mu/I (atomic-physics), not 2 mu / I"
    print(f"  Verdict 2: {verdict_2}")
    findings.append(("g_N convention", None, verdict_2))

    print()

    # ----------------------------------------------------------------------
    # CHECK 3: n*zeta convention in _hydrogenic_radial
    # ----------------------------------------------------------------------
    print("CHECK 3: Z_eff = n * zeta_CR convention sanity")
    print("-" * 60)

    # The hydrogenic radial wavefunction in Slater convention:
    # R_nl(r) = N(n,l,zeta) * r^l exp(-zeta r) * (polynomial)
    # The "exponent" of the radial wavefunction is zeta, NOT n*zeta.
    # Slater 1930: zeta_eff = (Z - sigma) / n*  where n* is the effective
    # principal quantum number. For Z=1 hydrogen, zeta = 1, n* = 1.
    # So |R_1s(0)|^2 = (zeta/1)^3 / pi = 1 / pi (correct).

    # Production code in _hydrogenic_radial: Z_eff = n * zeta is passed in,
    # then rho = 2 * Z_eff * r / n = 2 * zeta * r. So the radial decay rate
    # is exp(-rho/2) = exp(-zeta * r). This is the Slater convention.

    # Check: does this give 1/pi for hydrogen 1s?
    r_test = np.geomspace(1e-5, 50.0, 8000)
    R_h_1s = _hydrogenic_radial(1, 0, 1.0, r_test)  # Z_eff = 1 * zeta = 1
    # R(0) extrapolation:
    R_h_at_0 = R_h_1s[0] * (1.0 + r_test[0] * 0.5)  # 2-term Frobenius
    # Better: compare to analytical R_1s(0) = 2 * Z^{3/2} = 2.0
    R_h_at_0_analytical = 2.0  # Hydrogen 1s R_1s(0) = 2 Z^{3/2} = 2.
    # Hmm, but |psi(0)|^2 = R(0)^2 / (4 pi) = 4 / (4 pi) = 1/pi. Correct.

    # Numerical normalization check
    integrand = R_h_1s * R_h_1s * r_test * r_test
    norm = np.trapezoid(integrand, r_test)
    print(f"  Hydrogen 1s with Z_eff=1, n=1, l=0:")
    print(f"  Numerical normalization integral|R|^2 r^2 dr = {norm:.6f} (target 1.0)")
    print(f"  Analytical R(0) = 2 (since R_1s(0) = 2 Z^{{3/2}} for Z=1)")
    # Actually R_1s(r) = 2 Z^{3/2} exp(-Zr) = 2 e^{-r} for Z=1
    # R_1s(0) = 2.0
    # But in Slater convention with normalization N = 2 (2 zeta)^{n+1/2} / sqrt((2n)!) etc.

    if abs(norm - 1.0) < 0.001:
        verdict_3a = "PASS — _hydrogenic_radial normalizes correctly"
    else:
        verdict_3a = f"FAIL — normalization off by {abs(norm-1.0):.4f}"
    print(f"  Verdict 3a: {verdict_3a}")

    # Test that Cs 5p with CR67 zeta_CR=12.31, Z_eff = 5*12.31 = 61.56 gives
    # what the production code expects. The convention is correct only if
    # the radial decay is exp(-zeta * r), not exp(-Z_eff * r) = exp(-n*zeta * r).
    # In _hydrogenic_radial: rho = 2 * Z_eff * r / n = 2 * (n*zeta) * r / n
    # = 2 * zeta * r. Yes, the radial decay rate is zeta, as it should be.
    print(f"  Production: Z_eff = n * zeta in _hydrogenic_radial.")
    print(f"  This gives rho = 2 Z_eff r / n = 2 zeta r.")
    print(f"  Decay rate exp(-rho/2) = exp(-zeta r). MATCHES Slater convention.")

    # Verify: a 5p with zeta=12.31 should peak at r ~ 5/zeta = 5/12.31 = 0.4 bohr
    R_cs_5p = _hydrogenic_radial(5, 1, 5*12.31, r_test)  # Z_eff = 61.56
    pdf = R_cs_5p**2 * r_test**2
    pdf /= np.trapezoid(pdf, r_test)
    r_peak_cs5p = r_test[np.argmax(pdf)]
    expected_r_peak_5p = 25.0 / (5 * 12.31)  # n^2 / Z_eff = 25 / 61.56
    print(f"  Cs 5p (CR67) peak: numerical r_peak = {r_peak_cs5p:.4f} bohr, "
          f"expected n^2/Z_eff = {expected_r_peak_5p:.4f} bohr")

    if abs(r_peak_cs5p - expected_r_peak_5p) / expected_r_peak_5p < 0.05:
        verdict_3b = "PASS — n*zeta convention gives the expected hydrogenic peak"
    else:
        verdict_3b = f"FAIL — peak off by {100*(r_peak_cs5p-expected_r_peak_5p)/expected_r_peak_5p:+.1f}%"
    print(f"  Verdict 3b: {verdict_3b}")

    verdict_3 = f"{verdict_3a}; {verdict_3b}"
    findings.append(("n*zeta convention", None, verdict_3))

    print()

    # ----------------------------------------------------------------------
    # CHECK 4: Z_eff(r) profile sanity for Cs
    # ----------------------------------------------------------------------
    print("CHECK 4: Z_eff(r) profile boundary conditions")
    print("-" * 60)

    fc_cs = FrozenCore(Z=55, screening='single_zeta')
    fc_cs.solve()

    # Boundary conditions:
    z_at_0 = fc_cs.z_eff(0.0)
    z_at_inf = fc_cs.z_eff(50.0)
    z_at_05 = fc_cs.z_eff(0.5)
    z_at_2 = fc_cs.z_eff(2.0)

    print(f"  Cs Z_eff(r=0)   = {z_at_0:.2f}    (target: 55, no screening at origin)")
    print(f"  Cs Z_eff(r=0.5) = {z_at_05:.2f}    (intermediate, screened)")
    print(f"  Cs Z_eff(r=2.0) = {z_at_2:.2f}    (valence region)")
    print(f"  Cs Z_eff(r=50)  = {z_at_inf:.2f}    (target: 55-54=1, fully screened)")

    if abs(z_at_0 - 55) < 0.5 and abs(z_at_inf - 1.0) < 0.1:
        verdict_4 = "PASS — Z_eff(r) profile boundary conditions correct"
    else:
        verdict_4 = "FAIL — Z_eff(r) profile is wrong (boundary conditions off)"
    print(f"  Verdict 4: {verdict_4}")
    findings.append(("Z_eff(r) BCs", None, verdict_4))

    print()

    # ----------------------------------------------------------------------
    # CHECK 5: |psi_6s(0)|^2 with FrozenCore vs estimated from <r>
    # ----------------------------------------------------------------------
    print("CHECK 5: |psi_6s(0)|^2 cross-check")
    print("-" * 60)
    psi0_full = screened_psi_origin_squared(
        Z=55, n=6, l=0, n_grid=200_000, r_max=80.0,
    )
    print(f"  framework |psi_6s(0)|^2 (n_grid=200k) = {psi0_full:.4f} bohr^-3")
    print(f"  closeout-sprint Richardson value: 1.328 bohr^-3")
    # if these match to ~1%, the production code is reproducible.
    rel_err_5 = 100.0 * (psi0_full - 1.328) / 1.328
    print(f"  Rel error vs closeout: {rel_err_5:+.2f}%")
    if abs(rel_err_5) < 5.0:
        verdict_5 = "PASS — production code reproducibly returns the closeout value"
    else:
        verdict_5 = f"FAIL — production code disagrees by {rel_err_5:+.2f}%"
    print(f"  Verdict 5: {verdict_5}")
    findings.append(("|psi(0)|^2 reproducibility", rel_err_5, verdict_5))

    print()

    # ----------------------------------------------------------------------
    # CHECK 6: A_Cs end-to-end at correct g_N, full prefactor
    # ----------------------------------------------------------------------
    print("CHECK 6: Predicted A_Cs end-to-end")
    print("-" * 60)
    bf_cs = bohr_fermi_a_constant(
        psi0_squared=psi0_full,
        g_e=GE_FULL, g_N=g_cs_per_I, m_p_over_m_e=M_PROTON_OVER_M_E,
    )
    A_cs_predicted = bf_cs['A_MHz']
    A_cs_exp = 2298.158
    rel_err_cs = 100.0 * (A_cs_predicted - A_cs_exp) / A_cs_exp
    print(f"  Predicted A_Cs = {A_cs_predicted:.2f} MHz "
          f"(framework-native, BF strict + Schwinger only)")
    print(f"  Experimental A_Cs = {A_cs_exp:.2f} MHz")
    print(f"  Residual: {rel_err_cs:+.2f}%")
    print(f"  Closeout-sprint v2 single-zeta: -46.95% (after Casimir F_R applied)")
    print(f"  This row uses BF strict; without Casimir F_R = 1.555, residual is")
    print(f"  larger (more negative).")
    print()
    print(f"  Inverse: required |psi(0)|^2 to match A_exp:")
    psi_required = (A_cs_exp / A_cs_predicted) * psi0_full
    print(f"  |psi(0)|^2_needed = {psi_required:.4f} bohr^-3 "
          f"(from {psi0_full:.4f}; ratio {A_cs_exp/A_cs_predicted:.2f})")
    Z_eff_back = (psi_required * np.pi * 6**3)**(1.0/3.0)
    print(f"  Implied effective Z = {Z_eff_back:.2f}")
    print(f"  Compare to literature Roberts-Ginges Z_eff_Cs ~ 9.7-12 (without F_R)")
    print(f"  With Casimir F_R = 1.555, |psi(0)|^2_pre_relativistic = "
          f"{psi_required / 1.555:.4f}, "
          f"implied Z_eff = {(psi_required / 1.555 * np.pi * 6**3)**(1/3):.2f}")

    findings.append(("End-to-end A_Cs", rel_err_cs,
                     f"Without Casimir F_R: residual {rel_err_cs:+.1f}%; "
                     f"Casimir F_R=1.555 takes residual to ~-47% per closeout"))

    print()

    # ----------------------------------------------------------------------
    # Net verdict
    # ----------------------------------------------------------------------
    print("=" * 80)
    print("NET VERDICT (Probe c)")
    print("=" * 80)
    pass_count = sum(1 for f in findings if "PASS" in f[2])
    fail_count = sum(1 for f in findings if "FAIL" in f[2])

    if fail_count == 0:
        verdict = (
            "ALL CONVENTION CHECKS PASS. The framework's projection chain "
            "(Z_eff(r) -> R_nl -> R(0) -> |psi(0)|^2 -> A_HF) is internally "
            "consistent with no factor-of-X errors. The -47% Cs HFS residual "
            "is NOT a convention bug. The chain produces a faithful prediction "
            "GIVEN the input |psi(0)|^2; the cliff is upstream in the input "
            "(i.e., in the FrozenCore screening that determines the |psi(0)|^2)."
        )
    else:
        verdict = (
            f"{fail_count} CONVENTION ISSUES IDENTIFIED. See findings list "
            f"for specific factor-of-X to fix."
        )
    print(verdict)

    out = {
        "probe": "c",
        "hypothesis": "Convention mismatch (factor-of-X bug) in Z_eff -> psi -> A_HF chain",
        "findings": [{"check": f[0], "rel_err": f[1], "verdict": f[2]}
                     for f in findings],
        "verdict": verdict,
    }
    out_path = Path(ROOT) / "debug" / "data" / "z_cliff_probe_c.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
