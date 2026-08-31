"""Backing tests for the Paper 34 SecV.C Roothaan-autopsy baselines.

WHY THIS FILE EXISTS (2026-08-28).  SecV.C -- ~19 autopsies, ~2,500 lines -- had
NO test file.  That is precisely why a 496.5 ppm error in the deuterium
Bohr-Fermi baseline survived three certifying /qa runs: no code reviewer enters a
region with no tests, and the defect was column arithmetic, which prose reviewers
walk past.  The region was a logged "honest ceiling" for three certifications.

The defect: the drivers used an atomic-Hamiltonian g that is TWICE the CODATA
g_CODATA = mu_I/(I mu_N) the A-constant requires, together with a mass factor
m_e/m_N where it must be m_e/m_p for every nucleus (mu_N = e.hbar/(2 m_p) is
defined with the proton mass).  Net error: 2*(m_p/m_N).  At I=1 the deuteron's
m_d/m_p = 1.999 nearly cancels the doubling, which is why D looked right to three
digits.  H and T escaped because THEIR drivers used g_CODATA and m_e/m_p -- a
property of those drivers, not of I=1/2 (He-3 is I=1/2 and IS affected).

NOTE 2026-08-28: an earlier version of this docstring gave the divisor as 2I.
It is 2; the two agree only at I=1, the case the rule was derived from.

These tests encode the check that caught it -- cross-anchoring the hyperfine
baselines against each other -- so the class cannot recur silently.
"""
from __future__ import annotations

import math

import pytest

# CODATA nuclear magnetic moments, in units of the nuclear magneton, and spins.
MOMENTS = {
    "H": (2.79284734463, 0.5),    # proton
    "D": (0.85743823,    1.0),    # deuteron
    "T": (2.97896246,    0.5),    # triton
}

# Paper 34 SecV.C framework-native Bohr-Fermi baselines (MHz), point nucleus,
# g_e = 2, no recoil.  D corrected 2026-08-28 from 327.397464.
BASELINES = {
    "H": 1421.1595,
    "D": 327.234993,
    "T": 1515.865482,
}

ALPHA = 7.2973525693e-3
SCHWINGER = 1.0 + ALPHA / (2.0 * math.pi)

# reduced-mass factors and leading-Zemach shifts as the autopsy tables list them
CHAIN = {                      # (reduced-mass factor, leading Zemach in ppm, experiment MHz)
    "H": (0.998366, 39.495, 1420.405751768),
    "D": (0.999183, 98.001,  327.384352522),
    "T": (0.999454, 66.594, 1516.701470773),   # r_Z(t) = 1.762 fm
}


def _weight(species: str) -> float:
    """A_hf * (I + 1/2) up to constants shared by all 1s, Z=1 hyperfine splittings.

    The A-constant carries mu_I / I; the splitting between F = I +/- 1/2 carries a
    further (I + 1/2).  The mass factor m_e/m_p and everything else is common, so
    it cancels in any ratio taken at fixed Z and n.
    """
    mu, spin = MOMENTS[species]
    return (mu / spin) * (spin + 0.5)


@pytest.mark.parametrize("species", ["D", "T"])
def test_paper34_hfs_baselines_cross_anchor_to_hydrogen(species):
    """Every 1s Z=1 hyperfine baseline must be consistent with H 21cm.

    H, D and T are all 1s, Z = 1, and every autopsy table lists reduced mass as a
    SEPARATE multiplicative row, so the ratio of baselines is fixed purely by the
    nuclear moment and the (I + 1/2) multiplicity.  This is the check that caught
    the D defect; tritium is the control that proves the relation (it reproduces
    to 0.02 ppm).
    """
    predicted = BASELINES["H"] * _weight(species) / _weight("H")
    got = BASELINES[species]
    dev_ppm = (got - predicted) / predicted * 1e6
    assert abs(dev_ppm) < 30.0, (
        f"{species} Bohr-Fermi baseline {got} deviates {dev_ppm:+.1f} ppm from the "
        f"H-anchored prediction {predicted:.6f}. This is the 2026-08-28 defect class "
        f"(2I vs m_N/m_p in the g_atomic convention, and m_e/m_N vs m_e/m_p)."
    )


def test_paper34_hfs_cross_anchor_would_catch_the_2I_bug():
    """Non-tautology guard: the retired D baseline must FAIL the cross-anchor.

    Without this the test above could pass by construction if BASELINES were ever
    regenerated from the same buggy routine it is meant to police.
    """
    retired = 327.397464                      # the value live until 2026-08-28
    predicted = BASELINES["H"] * _weight("D") / _weight("H")
    dev_ppm = (retired - predicted) / predicted * 1e6
    assert dev_ppm > 400.0, "the cross-anchor no longer detects the retired D baseline"
    # and the corrected value must sit where the bad one did not
    assert abs((BASELINES["D"] - predicted) / predicted * 1e6) < 30.0
    # the bug is exactly the ratio (m_d/m_p) / 2
    m_d_over_m_p = 1.99900750139
    assert math.isclose(retired / BASELINES["D"], 2.0 / m_d_over_m_p, rel_tol=1e-5)


@pytest.mark.parametrize("species", ["H", "D", "T"])
def test_paper34_autopsy_chain_reproduces_its_stated_residual(species):
    """Each autopsy's printed chain must reproduce its own printed final value."""
    red, zem_ppm, exp = CHAIN[species]
    chain = BASELINES[species] * SCHWINGER * red * (1.0 - zem_ppm / 1e6)
    # values as PRINTED in the SecV.C autopsy tables (read, not inferred)
    expected_final = {"H": 1420.4318, "D": 327.315305, "T": 1516.696899}[species]
    assert abs(chain - expected_final) < 5e-3, (
        f"{species}: chain {chain:.6f} does not reproduce the stated {expected_final}"
    )


def test_paper34_deuteron_structure_correction_is_positive():
    """SecV.C D autopsy: the deuteron's net nuclear-structure correction is POSITIVE,
    opposite in sign to the proton's and triton's.

    This is the structural finding that replaced the withdrawn '+285.6 ppm matches
    the -286 ppm PY budget' closure claim.  Dividing out nu_F and the Schwinger
    factor isolates the structure correction each observable requires; for the
    weakly bound (2.2 MeV) deuteron, polarizability outweighs the negative Zemach
    term and flips the sign, so a leading-Zemach-only projection points the wrong
    way at I = 1.
    """
    needed = {}
    for sp, (red, zem_ppm, exp) in CHAIN.items():
        nu_F = BASELINES[sp] * red
        needed[sp] = ((exp / nu_F) / SCHWINGER - 1.0) * 1e6

    # the sign contrast is the claim
    assert needed["H"] < 0, f"H structure correction not negative: {needed['H']:+.1f}"
    assert needed["T"] < 0, f"T structure correction not negative: {needed['T']:+.1f}"
    assert needed["D"] > 0, f"D structure correction not positive: {needed['D']:+.1f}"

    # magnitudes as tabulated in the paper
    assert abs(needed["H"] - (-55.9)) < 2.0
    assert abs(needed["T"] - (-63.2)) < 2.0
    assert abs(needed["D"] - (+112.9)) < 2.0

    # hydrogen's leading Zemach captures most of its structure correction ...
    assert 0.6 < CHAIN["H"][1] / abs(needed["H"]) < 0.8
    # Tritium: the Zemach term is the right SIGN, which is the load-bearing
    # half of the contrast.  Its MAGNITUDE is deliberately not asserted:
    # the 2026-08-28 citation audit found the autopsy feeds r_Z(t) = 1.762 fm,
    # which is the triton CHARGE radius (1.7591 fm), not its Zemach radius.
    # The proton ratio r_Z/r_charge = 1.243 implies r_Z(t) ~ 2.2 fm, which
    # would move the applied shift from -66.6 to ~-86 ppm and the gap from
    # +3.4 to ~+22.6 ppm.  Re-sourcing that radius is an open PI call, so
    # this leg pins only the sign until it is settled.
    gap_T = needed["T"] + CHAIN["T"][1]
    assert gap_T > -20.0, (
        f"T Zemach term over-corrects: gap {gap_T:+.1f} ppm. The applied term "
        f"should not exceed the required structure correction in magnitude by "
        f"more than the r_Z uncertainty."
    )

    # ... while for D the applied Zemach has the WRONG SIGN, leaving a +210.9 ppm gap
    gap_D = needed["D"] + CHAIN["D"][1]
    assert abs(gap_D - 210.9) < 2.0, f"D non-Zemach structure gap {gap_D:+.1f} != +210.9"

    # and that gap is far larger than the +44 ppm the catalogue attributes to
    # deuteron polarizability -- the inconsistency flagged for re-sourcing
    assert gap_D > 4.0 * 44.0


# ---------------------------------------------------------------------------
# Widened coverage (2026-08-28, second pass).  The first version of this file
# covered H/D/T only -- 3 of 19 SecV.C autopsies -- and the Li-7 baseline and the
# alkali-cliff table sat in the uncovered 16, each carrying a half of the SAME
# convention defect that had just been fixed for deuterium.  Cover the rest of
# the Fermi-contact family and pin the root-cause formula itself.
# ---------------------------------------------------------------------------

MU_LI7, I_LI7 = 3.256427, 1.5          # CODATA moment, nuclear spin
Z_EFF_LI7, N_LI7 = 1.279, 2            # Clementi-Roetti 1967 screening

# nuclear masses in electron-mass units (atomic mass minus Z electron masses)
_U = 5.4857990907e-4
M_NUC_ME = {
    "D":    (2.014102 - 1 * _U) / _U,
    "T":    (3.016049 - 1 * _U) / _U,
    "He-3": (3.016029 - 2 * _U) / _U,
    "Li-7": (7.016003 - 3 * _U) / _U,
}
M_P_ME = 1836.15267343


def test_paper34_li7_baseline_cross_anchors_to_hydrogen():
    """Li-7 2s HFS splitting must cross-anchor to H 21cm at its stated Z_eff.

    Same relation as H/D/T, extended by the hydrogenic contact-density factor
    Z_eff^3 / n^3.  The retired baseline (82.977 MHz) fails this by 3.48x.
    """
    w_li = (MU_LI7 / I_LI7) * (I_LI7 + 0.5)
    predicted = BASELINES["H"] * (w_li / _weight("H")) * Z_EFF_LI7 ** 3 / N_LI7 ** 3
    assert abs(predicted - 288.913) < 0.3, f"Li-7 prediction {predicted:.3f} != 288.913"

    retired = 82.977
    assert abs(retired / predicted - 1) > 0.5, "cross-anchor no longer rejects the retired Li-7 value"


def test_paper34_track5_convention_error_formula():
    """The unified root cause: the buggy convention gives correct * 2*(m_p/m_N).

    g_atomic = 2*mu_I/mu_N must be divided by 2I to recover mu_I/I, and the mass
    factor must be m_e/m_p (mu_N is defined with the proton mass).  Doing neither
    scales the result by 2*(m_p/m_N).  This single expression reproduces every
    known instance, which is what makes it a root cause rather than a pattern.
    """
    def factor(species):
        return 2.0 * M_P_ME / M_NUC_ME[species]

    # deuterium: high by 496.5 ppm (retracted 2026-08-28)
    assert abs((factor("D") - 1.0) * 1e6 - 496.5) < 2.0

    # Li-7: low by 3.4818x -- matches the observed 82.977/288.913 to 6 digits
    assert abs(1.0 / factor("Li-7") - 3.4818) < 1e-3
    assert abs(82.977 / 288.913 - factor("Li-7")) < 1e-5, (
        "the Li-7 deficit no longer matches 2*(m_p/m_N); the root-cause "
        "identification in sec:autopsy_he3_hfs would need revisiting"
    )

    # He-3: the 'off by 3/2' the He-3 autopsy reports
    assert abs(1.0 / factor("He-3") - 1.5) < 0.01

    # tritium escaped the bug entirely -- its driver used the standard convention,
    # which is why it serves as the cross-anchor control above
    predicted_T = BASELINES["H"] * _weight("T") / _weight("H")
    assert abs((BASELINES["T"] - predicted_T) / predicted_T) < 3e-5


def test_paper34_li7_autopsy_and_alkali_table_agree():
    """The Li-7 autopsy reports a SPLITTING (2*A_hf); the alkali table reports A_hf.

    Before the 2026-08-28 correction these two loci disagreed (autopsy -89.7%,
    table -94.8%) because each carried a different half of the same defect.  They
    must now be consistent, which is the cheapest guard against a partial re-fix.
    """
    autopsy_splitting = 288.913
    table_A = 144.590
    assert abs(autopsy_splitting / table_A - 2.0) < 5e-3, (
        f"autopsy/table ratio {autopsy_splitting / table_A:.4f} != 2 "
        f"(the autopsy is a splitting, the table an A-constant)"
    )


def test_paper34_alkali_cliff_increases_with_Z():
    """Corrected alkali cliff is monotonically deeper with Z, not shallower.

    The retired table showed K/Rb/Cs all at '0.3 MHz' -- three different atoms
    coinciding to one decimal, the tell for the suppressed mass slot -- and
    supported a '~Z^2.5' growth and an 'inverse-Z' reading.  Both are withdrawn.
    """
    Z = [3, 11, 19, 37, 55]
    A_fw = [144.590, 219.756, 10.497, 23.465, 36.506]
    A_exp = [401.752, 885.813, 230.860, 1011.911, 2298.158]
    cliff = [(f - e) / e * 100.0 for f, e in zip(A_fw, A_exp)]

    for a, b in zip(cliff, cliff[1:]):
        assert b < a, f"cliff not monotonically deeper with Z: {cliff}"
    assert abs(cliff[0] - (-64.0)) < 0.5 and abs(cliff[-1] - (-98.4)) < 0.5

    # no two atoms may share a suppressed A_fw to one decimal (the original tell)
    rounded = [round(a, 1) for a in A_fw]
    assert len(set(rounded)) == len(rounded), f"degenerate A_fw entries: {rounded}"

    # The fitted log-log slope of the ratio column is 1.16.  This is pinned
    # as a DATA fact about the corrected table (it discriminates the retired
    # ~Z^2.5), NOT as a growth law: the fit misses Na by 126% because the
    # CR67 denominator is non-monotone in Z.  The structural account is
    # eq:alkali_cliff, tested in ..._cliff_closed_form below.
    import math
    n = len(Z)
    lx = [math.log(z) for z in Z]
    ly = [math.log(e / f) for e, f in zip(A_exp, A_fw)]
    mx, my = sum(lx) / n, sum(ly) / n
    slope = sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)
    assert abs(slope - 1.16) < 0.15, f"density-ratio exponent {slope:.2f} != 1.16"


# ---------------------------------------------------------------------------
# Deuteron polarizability channel (follow-on sprint, 2026-08-29)
#
# Primary source: Bonilla et al., arXiv:2508.18776 -- electronic D 1S
# two-photon-exchange decomposition.  Backs Paper 34 sec:autopsy_d_hfs
# "What the channel has to deliver" and the corrected V.B / V.D.5 magnitudes.
# ---------------------------------------------------------------------------

TPE_KHZ = {"elastic": -41.50, "polarizability": +110.16,      # as published
           "single_p": -33.14, "single_n": +8.95}


def _ppm_of_nu(khz, species="D"):
    """A kHz shift, as ppm of that species' measured splitting."""
    nu_mhz = CHAIN[species][2]
    return khz * 1e3 / (nu_mhz * 1e6) * 1e6


def _chain_no_structure(species):
    """BF x Schwinger x reduced mass -- the chain before any nuclear structure."""
    rm = CHAIN[species][0]
    return BASELINES[species] * SCHWINGER * rm


def test_paper34_d_hfs_tpe_decomposition_and_sign_flip():
    """The sourced TPE totals +44.5 kHz = +135.9 ppm, and the inelastic
    channel outweighs ALL elastic-class pieces -- the sign flip the autopsy's
    I=1 argument predicts on structural grounds alone."""
    total = sum(TPE_KHZ.values())
    assert abs(total - 44.5) < 0.1, f"TPE total {total:.2f} != 44.5 kHz"
    assert abs(_ppm_of_nu(total) - 135.9) < 0.5

    inelastic = _ppm_of_nu(TPE_KHZ["polarizability"])
    elastic = _ppm_of_nu(TPE_KHZ["elastic"] + TPE_KHZ["single_p"]
                         + TPE_KHZ["single_n"])
    assert abs(inelastic - 336.5) < 0.5, f"inelastic {inelastic:.1f}"
    assert abs(elastic - (-200.6)) < 0.5, f"elastic-class {elastic:.1f}"
    # load-bearing: polarizability wins, so the total flips positive
    assert inelastic > -elastic > 0
    assert _ppm_of_nu(total) > 0


def test_paper34_d_hfs_channel_target_and_h_parallel():
    """What the channel must deliver, and why the leftover is credible:
    consuming the sourced TPE puts D at +23.0 ppm -- same sign and class as
    the architecturally identical H 21cm chain at +18.4 ppm."""
    chain_noZ = _chain_no_structure("D")
    assert abs(chain_noZ - 327.347385) < 2e-5

    exp_d = CHAIN["D"][2]
    required = (exp_d - chain_noZ) / chain_noZ * 1e6
    assert abs(required - 112.9) < 0.5, f"required structure {required:.1f} ppm"
    # the paper's quoted gap is required minus the leading-Zemach shift
    assert abs((required + CHAIN["D"][1]) - 210.9) < 0.5

    with_tpe = chain_noZ * (1.0 + _ppm_of_nu(sum(TPE_KHZ.values())) * 1e-6)
    resid = (with_tpe - exp_d) / exp_d * 1e6
    assert abs(resid - 23.0) < 0.5, f"D-with-TPE residual {resid:.1f} ppm"

    # H parallel: identical chain shape with its own structure input consumed
    h_chain = _chain_no_structure("H") * (1.0 - CHAIN["H"][1] * 1e-6)
    h_resid = (h_chain - CHAIN["H"][2]) / CHAIN["H"][2] * 1e6
    assert 15.0 < h_resid < 22.0, f"H 21cm residual {h_resid:.1f} ppm"
    # both positive, both tens of ppm -- the same missing-QED class
    assert 0 < h_resid < resid < 100.0


def test_paper34_d_polarizability_unit_slip_guard():
    """The retired '+44 ppm deuteron polarizability' was the kHz total read as
    ppm.  Pin both halves so neither can return: 44.5 kHz IS 136 ppm, and the
    total is not the polarizability."""
    assert abs(_ppm_of_nu(44.5) - 135.9) < 0.5
    exp_d = CHAIN["D"][2]
    chain = 327.315305 * (1.0 + 44e-6)          # the retired reading
    assert (chain - exp_d) / exp_d * 1e6 < -150.0
    ratio = TPE_KHZ["polarizability"] / sum(TPE_KHZ.values())
    assert ratio > 2.0, f"polarizability/total = {ratio:.2f}"


def test_paper34_alkali_cliff_closed_form():
    """Paper 34 eq:alkali_cliff (2026-08-29 diagnosis).  The cliff is the ratio
    of the Fermi-Segre contact density to the hydrogenic single-zeta estimate:

        cliff = (Z / Z_eff^3) * (n / nu)^3 * F_rel(Z alpha)

    with nu = sqrt(Ry / E_ion) the effective quantum number from the MEASURED
    ionization energy and F_rel = 1/[g(2g-1)], g = sqrt(1-(Z alpha)^2).  This
    reproduces every tabulated cliff to <=10%, where the retired ~Z^1.2 power
    law misses Na by 126%.  It also pins the structural reading: the
    quantum-defect factor dominates, and nu is nearly CONSTANT along the
    series while n runs 2..6.
    """
    Z = [3, 11, 19, 37, 55]
    n = [2, 3, 4, 5, 6]
    Zeff = [1.279, 2.507, 2.162, 2.771, 3.475]
    cliff = [2.8, 4.0, 22.0, 43.1, 63.0]
    E_ion = [5.3917, 5.1391, 4.3407, 4.1771, 3.8939]     # eV, NIST
    RY = 13.605693

    nu = [math.sqrt(RY / e) for e in E_ion]
    gam = [math.sqrt(1.0 - (z * ALPHA) ** 2) for z in Z]
    F = [1.0 / (g * (2.0 * g - 1.0)) for g in gam]

    pred = [(Z[i] / Zeff[i] ** 3) * (n[i] / nu[i]) ** 3 * F[i] for i in range(5)]
    for i in range(5):
        rel = abs(pred[i] - cliff[i]) / cliff[i]
        assert rel < 0.11, (
            f"Z={Z[i]}: closed form {pred[i]:.1f} vs tabulated {cliff[i]} "
            f"({rel*100:.0f}% off)")

    # the retired power law is genuinely worse -- this is why it was retired
    lx = [math.log(z) for z in Z]
    ly = [math.log(c) for c in cliff]
    mx, my = sum(lx) / 5, sum(ly) / 5
    slope = (sum((a - mx) * (b - my) for a, b in zip(lx, ly))
             / sum((a - mx) ** 2 for a in lx))
    pl = [math.exp(my + slope * (a - mx)) for a in lx]
    worst_pl = max(abs(pl[i] - cliff[i]) / cliff[i] for i in range(5))
    worst_cf = max(abs(pred[i] - cliff[i]) / cliff[i] for i in range(5))
    assert worst_pl > 1.0, f"power-law worst error {worst_pl:.2f} (expected >100%)"
    assert worst_cf < worst_pl / 5.0, "closed form must beat the power law decisively"

    # STRUCTURAL: nu is nearly constant while n runs 2..6 -- the cliff's driver
    assert max(nu) - min(nu) < 0.35, f"nu should be near-constant: {nu}"
    assert all(1.5 < v < 2.0 for v in nu)
    qd_factor = [(n[i] / nu[i]) ** 3 for i in range(5)]
    charge_factor = [Z[i] / Zeff[i] ** 3 for i in range(5)]
    # quantum-defect factor spans ~17x and is monotone; charge factor spans
    # <3x and is NOT monotone -- so the defect term is the driver
    assert qd_factor[-1] / qd_factor[0] > 15.0
    assert all(b > a for a, b in zip(qd_factor, qd_factor[1:]))
    assert max(charge_factor) / min(charge_factor) < 3.0
    assert not all(b > a for a, b in zip(charge_factor, charge_factor[1:]))

    # and the CR67 estimate itself is non-monotone (the K "jump" artifact)
    psi_cr67 = [Zeff[i] ** 3 / (math.pi * n[i] ** 3) for i in range(5)]
    assert psi_cr67[2] < psi_cr67[1], "K should dip below Na in the CR67 estimate"
    psi_need = [cliff[i] * psi_cr67[i] for i in range(5)]
    assert all(b > a for a, b in zip(psi_need, psi_need[1:])), (
        f"the NEEDED density must be smooth/monotone: {psi_need}")
