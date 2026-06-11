"""
Sprint HF Track 4 — Zemach correction diagnostic
================================================

Tests whether GeoVac's existing form-factor / finite-size machinery
(geovac/nuclear/form_factor.py + nuclear_electronic.py::finite_size_coupling_pauli)
produces the Zemach correction to A_hf at the right size, OR whether Zemach
requires native machinery the framework doesn't have.

Track NI's finite_size_coupling_pauli is a correction to the 1s *binding energy*
(Foldy correction): δE_1s = (2/5) Z^4 R_p^2 / n^3 in Hartree, parameterized by
R_PROTON_BOHR (charge radius, 0.8414 fm).

The Zemach correction is to the *hyperfine coupling*:
    Δν_Zemach / ν_F = −2 Z α (r_Z / a_0)
where r_Z = ∫ d^3r d^3r' ρ(r) m(r') |r − r'| is the convolution of the proton
charge distribution ρ(r) with its magnetization distribution m(r).
For hydrogen: r_Z ≈ 1.045 fm; Δν/ν ≈ −33 ppm (Eides Ch. 7, Karshenboim 2005).

The diagnostic:
1. Read form_factor.py and nuclear_electronic.py end-to-end. Verify there is no
   magnetization distribution operator, no separate r_Z input, no size-dependent
   modulation of A_hf in hyperfine_coupling_pauli.
2. Predict Δν/ν under three scenarios:
   (A) charge-radius-as-Zemach-proxy: r_substitute = R_PROTON_BOHR
   (B) r_Z external: r_substitute = 1.045 fm
   (C) native Zemach: zero (no machinery exists)
3. Compute updated A_hf and residual against experimental 1420.4058 MHz.
4. Place Zemach radius in Paper 34 vocabulary.

Author: GeoVac Development Team
Date: 2026-05-07 (HF Track 4)
"""

from __future__ import annotations

import inspect
import json
import re
from pathlib import Path
from typing import Any, Dict


# ---------------------------------------------------------------------------
# Constants (CODATA 2018 + Eides Ch. 7)
# ---------------------------------------------------------------------------

ALPHA = 7.2973525693e-3                    # CODATA 2018 fine-structure
A0_FM = 52917.72108                        # Bohr radius in fm
R_PROTON_FM = 0.8414                       # CODATA 2018 charge radius
R_PROTON_BOHR = R_PROTON_FM / A0_FM        # 1.5901e-5 bohr
R_ZEMACH_FM = 1.045                        # Eides §7.2 (extracted from ep
                                           # scattering and muonic-H spectroscopy;
                                           # quoted ~1.045(16) fm)
R_ZEMACH_BOHR = R_ZEMACH_FM / A0_FM        # 1.974e-5 bohr
Z = 1                                      # hydrogen

# Track HF-2 closing baseline:  GeoVac BF + recoil + Schwinger a_e
A_HF_HF2_BASELINE_MHZ = 1420.4879
A_HF_EXPERIMENTAL_MHZ = 1420.4058          # NIST/CODATA 21cm line


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def ppm(x_predicted_mhz: float, x_actual_mhz: float = A_HF_EXPERIMENTAL_MHZ) -> float:
    return 1e6 * (x_predicted_mhz - x_actual_mhz) / x_actual_mhz


# ---------------------------------------------------------------------------
# Step 1 — Diagnostic of Track NI's proton-structure machinery
# ---------------------------------------------------------------------------

def diagnose_track_ni_proton_structure() -> Dict[str, Any]:
    """
    Read form_factor.py and nuclear_electronic.py.  Identify every reference
    to a proton-size parameter and every mechanism that could couple proton
    structure to A_hf.
    """
    repo = Path(__file__).resolve().parents[1]
    ff_path = repo / "geovac" / "nuclear" / "form_factor.py"
    ne_path = repo / "geovac" / "nuclear" / "nuclear_electronic.py"
    ff_src = ff_path.read_text(encoding='utf-8')
    ne_src = ne_path.read_text(encoding='utf-8')

    # Search for the canonical Zemach-relevant terms (word-boundary, no
    # case-folding for short identifiers to avoid spurious substring matches
    # such as 'R_M' inside 'r_max' or 'HA_PER_MEV').
    keywords_zemach = ["Zemach", "magnetization", "r_Z", "R_M", "magnetic moment distribution"]
    keywords_charge = ["R_PROTON", "R_nuc", "charge_radius", "charge radius"]

    src_combined = ff_src + ne_src

    def _word_count(needle: str, src: str) -> int:
        # Use word boundaries (avoids 'R_M' matching 'r_max').
        return len(re.findall(r"\b" + re.escape(needle) + r"\b", src))

    def _substr_count(needle: str, src: str) -> int:
        # Plain substring match for prefixes like 'R_PROTON' that are always
        # followed by suffixes ('_BOHR', '_FM').
        return len(re.findall(re.escape(needle), src))

    # For Zemach: short tokens 'r_Z', 'R_M' need word boundaries; multi-word
    # phrases use simple substring; canonical names ('Zemach', 'magnetization')
    # use substring (case-insensitive ok because no false positives).
    zemach_hits = {
        "Zemach": len(re.findall("zemach", src_combined, re.IGNORECASE)),
        "magnetization": len(re.findall("magnetization", src_combined, re.IGNORECASE)),
        "r_Z": _word_count("r_Z", src_combined),
        "R_M": _word_count("R_M", src_combined),
        "magnetic moment distribution": len(re.findall(
            "magnetic moment distribution", src_combined, re.IGNORECASE)),
    }
    # For charge-radius: prefixes need substring; phrases work either way.
    charge_hits = {
        "R_PROTON": _substr_count("R_PROTON", src_combined),
        "R_nuc": _word_count("R_nuc", src_combined),
        "charge_radius": _substr_count("charge_radius", src_combined),
        "charge radius": len(re.findall("charge radius", src_combined, re.IGNORECASE)),
    }

    # Inspect the hyperfine coupling signature
    from geovac.nuclear.nuclear_electronic import (
        hyperfine_coupling_pauli,
        finite_size_coupling_pauli,
        finite_size_correction,
    )
    hf_sig = inspect.signature(hyperfine_coupling_pauli)
    fs_sig = inspect.signature(finite_size_coupling_pauli)
    fsc_sig = inspect.signature(finite_size_correction)

    hf_takes_R_nuc = "R_nuc" in hf_sig.parameters
    hf_takes_r_Z = "r_Z" in hf_sig.parameters or "r_zemach" in hf_sig.parameters
    fs_takes_R_nuc = "R_nuc" in fs_sig.parameters

    # Inspect the hyperfine coupling source for any size dependence
    hf_source = inspect.getsource(hyperfine_coupling_pauli)
    has_size_in_hf = bool(re.search(r"R_nuc|R_PROTON|r_Z|magnetization",
                                    hf_source, re.IGNORECASE))

    fs_source = inspect.getsource(finite_size_coupling_pauli)
    fs_targets_hyperfine = bool(re.search(
        r"hyperfine|HF|A_hf|I\.S|I_dot_S", fs_source, re.IGNORECASE))

    # The closed-form finite_size_correction returns a perturbation to E_ns,
    # not to A_hf.  Verify by inspecting the docstring + the math.
    fsc_doc = (finite_size_correction.__doc__ or "")
    fsc_targets_binding_energy = "ns" in fsc_doc.lower() or "energy correction" in fsc_doc.lower()

    return {
        "zemach_keyword_counts_in_geovac_nuclear": zemach_hits,
        "charge_radius_keyword_counts_in_geovac_nuclear": charge_hits,
        "hyperfine_coupling_pauli_signature": str(hf_sig),
        "hyperfine_takes_R_nuc": hf_takes_R_nuc,
        "hyperfine_takes_r_Z": hf_takes_r_Z,
        "hyperfine_source_has_size_dependence": has_size_in_hf,
        "finite_size_coupling_pauli_signature": str(fs_sig),
        "finite_size_takes_R_nuc": fs_takes_R_nuc,
        "finite_size_targets_hyperfine_energy": fs_targets_hyperfine,
        "finite_size_correction_signature": str(fsc_sig),
        "finite_size_correction_targets_binding_energy": fsc_targets_binding_energy,
        "finite_size_correction_returns_scalar": True,  # by construction (returns float Ha)
        "summary_no_magnetization_operator": all(v == 0 for v in zemach_hits.values()),
    }


# ---------------------------------------------------------------------------
# Step 2 — Three Zemach scenarios
# ---------------------------------------------------------------------------

def zemach_shift_ppm(r_substitute_bohr: float) -> float:
    """
    Standard Eides Ch. 7 Zemach correction (Eides Eq. (7.4)):

        Delta nu_Z / nu_F = -2 Z alpha m_e r_Z

    in natural units.  In atomic units (m_e = 1, a_0 = 1, c = 1/alpha,
    Compton wavelength lambda_C = hbar/(m_e c) = alpha bohr),

        m_e r_Z = r_Z / lambda_C = r_Z[bohr] / alpha

    so

        Delta nu_Z / nu_F = -2 Z alpha (r_Z[bohr] / alpha) = -2 Z r_Z[bohr]

    This is the canonical formula and reproduces the textbook -41 ppm
    (Eides; with corrections -33 ppm after polarizability subtractions).

    Note: the briefing's expression "-2 Z alpha (r_Z / a_0)" was off by a
    factor of alpha.  The correct factor combines a_0 and lambda_C such
    that the alpha cancels.

    Returns the shift in ppm.
    """
    return 1e6 * (-2.0 * Z * r_substitute_bohr)


def predict_three_scenarios() -> Dict[str, Any]:
    A_baseline = A_HF_HF2_BASELINE_MHZ
    A_exp = A_HF_EXPERIMENTAL_MHZ
    baseline_residual_ppm = ppm(A_baseline, A_exp)

    # (A) Charge-radius-as-Zemach proxy: use R_PROTON_BOHR
    delta_A_ppm = zemach_shift_ppm(R_PROTON_BOHR)
    A_pred_A = A_baseline * (1.0 + delta_A_ppm * 1e-6)
    res_A_ppm = ppm(A_pred_A, A_exp)

    # (B) r_Z as external focal length: 1.045 fm
    delta_B_ppm = zemach_shift_ppm(R_ZEMACH_BOHR)
    A_pred_B = A_baseline * (1.0 + delta_B_ppm * 1e-6)
    res_B_ppm = ppm(A_pred_B, A_exp)

    # (C) Native Zemach in framework: zero (no machinery)
    A_pred_C = A_baseline
    res_C_ppm = baseline_residual_ppm

    # Cross-check: the "−33 ppm" canonical Eides value
    canonical_eides_ppm = zemach_shift_ppm(R_ZEMACH_BOHR)

    return {
        "baseline": {
            "name": "HF-2 closing: BF + recoil + Schwinger",
            "A_hf_MHz": A_baseline,
            "residual_ppm": baseline_residual_ppm,
        },
        "scenario_A_charge_radius_substitution": {
            "description": (
                "Use R_PROTON_BOHR (charge radius 0.8414 fm) as a stand-in for "
                "the Zemach radius. This is what the framework would natively "
                "produce IF its R_PROTON_BOHR parameter were re-purposed as a "
                "magnetization radius — but the framework has no quantum operator "
                "for magnetization, so this remains a hand-substitution."
            ),
            "r_substitute_fm": R_PROTON_FM,
            "r_substitute_bohr": R_PROTON_BOHR,
            "delta_nu_over_nu_ppm": delta_B_ppm,  # for diagnostic; uses 1.045 fm
            "delta_nu_using_R_PROTON_ppm": delta_A_ppm,
            "A_hf_predicted_MHz": A_pred_A,
            "residual_ppm": res_A_ppm,
        },
        "scenario_B_r_Z_external": {
            "description": (
                "Inject the experimentally extracted Zemach radius r_Z ≈ 1.045 fm "
                "as an external Layer-2 calibration focal length. Framework has no "
                "magnetization quantum operator, so this is calibrated supply, "
                "structurally identical to the recoil correction in HF-1."
            ),
            "r_Z_fm": R_ZEMACH_FM,
            "r_Z_bohr": R_ZEMACH_BOHR,
            "delta_nu_over_nu_ppm": delta_B_ppm,
            "A_hf_predicted_MHz": A_pred_B,
            "residual_ppm": res_B_ppm,
        },
        "scenario_C_native_zemach_zero": {
            "description": (
                "Framework has no Zemach machinery. Without external supply, "
                "the prediction is unchanged from the HF-2 baseline (+58 ppm)."
            ),
            "A_hf_predicted_MHz": A_pred_C,
            "residual_ppm": res_C_ppm,
        },
        "canonical_eides_zemach_correction_ppm": canonical_eides_ppm,
        "ratio_charge_to_zemach_radius": R_PROTON_FM / R_ZEMACH_FM,
        "ratio_zemach_to_charge_radius": R_ZEMACH_FM / R_PROTON_FM,
    }


# ---------------------------------------------------------------------------
# Step 3 — Verdict assembly
# ---------------------------------------------------------------------------

def assemble_verdict(
    diagnostic: Dict[str, Any],
    scenarios: Dict[str, Any],
) -> Dict[str, Any]:
    """Assemble the three-way verdict required by the briefing."""

    # The framework has zero Zemach-relevant machinery (no magnetization
    # operator, no separate r_Z input, no size dependence in HF coupling).
    # So:
    #   - There is NO native machinery for Zemach (rules out positive).
    #   - The single proton-structure parameter R_PROTON_BOHR exists but is
    #     wired to E_1s binding-energy correction, NOT to A_hf.
    #   - However, Eides §7.2 Zemach has the exact same functional form as
    #     a "−2 Z α (r/a_0)" correction; if R_PROTON_BOHR were used as a proxy,
    #     the framework would give the right *form* with the wrong *size*.
    # That is the partial-with-charge-radius-substitution outcome.

    no_zemach_machinery = (
        not diagnostic["hyperfine_takes_R_nuc"]
        and not diagnostic["hyperfine_takes_r_Z"]
        and not diagnostic["hyperfine_source_has_size_dependence"]
        and not diagnostic["finite_size_targets_hyperfine_energy"]
        and diagnostic["summary_no_magnetization_operator"]
    )

    # The "partial-with-substitution" lens: even charge-radius substitution
    # requires the user to declare "treat R_PROTON_BOHR as a magnetization
    # radius" — the framework itself does not perform any size-dependent
    # modulation of A_hf.  So technically the framework supplies neither
    # outcome (A) nor outcome (B); it supplies outcome (C), with two
    # available external substitutions of different fidelity.
    if no_zemach_machinery:
        verdict = "NEGATIVE"
        verdict_long = (
            "Track NI has no Zemach-relevant machinery. The hyperfine coupling "
            "is point-like in the proton spatial coordinate; there is no "
            "magnetization distribution operator, no separate r_Z input, and "
            "no size dependence in hyperfine_coupling_pauli. R_PROTON_BOHR is "
            "the framework's single proton-structure parameter and is wired "
            "into a one-body electronic operator (Foldy binding-energy "
            "correction), not the hyperfine coupling. Zemach is structurally "
            "external — like recoil — and must be supplied as a calibration "
            "focal length with r_Z as input."
        )
    else:
        verdict = "PARTIAL"
        verdict_long = "Framework has at least some Zemach-relevant size dependence."

    # Headline updated A_hf prediction = HF-2 baseline (no native Zemach)
    A_hf_native = scenarios["scenario_C_native_zemach_zero"]["A_hf_predicted_MHz"]
    res_native_ppm = scenarios["scenario_C_native_zemach_zero"]["residual_ppm"]

    # If the PI accepts r_Z as external Layer-2 calibration:
    A_hf_with_rZ = scenarios["scenario_B_r_Z_external"]["A_hf_predicted_MHz"]
    res_with_rZ_ppm = scenarios["scenario_B_r_Z_external"]["residual_ppm"]

    return {
        "verdict": verdict,
        "verdict_long": verdict_long,
        "headline_native_A_hf_MHz": A_hf_native,
        "headline_native_residual_ppm": res_native_ppm,
        "with_r_Z_external_A_hf_MHz": A_hf_with_rZ,
        "with_r_Z_external_residual_ppm": res_with_rZ_ppm,
    }


# ---------------------------------------------------------------------------
# Step 4 — Paper 34 projection vocabulary placement
# ---------------------------------------------------------------------------

def paper_34_vocabulary() -> Dict[str, Any]:
    """
    Place the proton charge radius, Zemach radius, and recoil correction
    in the Paper 34 projection-taxonomy vocabulary.

    Paper 34's three-axis tagging is (variable, dimension, transcendental).
    Each "focal length" external to the dimensionless graph is a projection
    that adds one or more of these axes.
    """
    return {
        "charge_radius_R_p": {
            "physical_meaning": "RMS charge radius of the proton",
            "value_fm": R_PROTON_FM,
            "what_it_modulates_in_GeoVac": (
                "Foldy finite-size correction to the 1s binding energy: "
                "δE_1s = (2/5) Z^4 R_p^2 / n^3 (Ha)"
            ),
            "paper_34_axes_added": "(variable: length scale; dimension: [L]; transcendental: none)",
            "paper_34_projection_class": "rest-mass / size / Layer-2 calibration",
            "category_per_HF_taxonomy": "C — external CODATA",
            "wired_into_GeoVac": True,
            "wiring": "form_factor.finite_size_correction → finite_size_coupling_pauli",
        },
        "zemach_radius_r_Z": {
            "physical_meaning": (
                "First moment of the convolution of the proton charge density "
                "ρ(r) with its magnetization density m(r): "
                "r_Z = ∫ d^3r d^3r' ρ(r) m(r') |r − r'|"
            ),
            "value_fm": R_ZEMACH_FM,
            "what_it_modulates_in_GeoVac": (
                "NOTHING (native): hyperfine_coupling_pauli is point-like in "
                "proton spatial coordinate. With external substitution, would "
                "modulate A_hf via Δν/ν = −2 Z α (r_Z / a_0)."
            ),
            "paper_34_axes_added": (
                "(variable: length scale (DIFFERENT from R_p; magnetization-weighted, "
                "not charge-weighted); dimension: [L]; transcendental: none)"
            ),
            "paper_34_projection_class": (
                "QCD-internal hadronic Layer-2 input. Categorically the same axis "
                "as r_charge but a STRUCTURALLY DIFFERENT focal length: charge "
                "distribution is ρ(r), magnetization distribution is m(r), and "
                "the framework distinguishes neither because both live downstream "
                "of an SU(3) gauge structure that GeoVac admits as a Wilson "
                "construction (Paper 30 + Sprint ST-SU3) but cannot autonomously "
                "compute matter elements within."
            ),
            "category_per_HF_taxonomy": "C — external (Eides §7, Karshenboim 2005)",
            "wired_into_GeoVac": False,
            "could_be_wired": (
                "Yes, as a Layer-2 multiplicative calibration on A_hf. Doing so "
                "would NOT be a derivation; it would be a calibration input "
                "structurally identical to the recoil prescription in HF-1."
            ),
        },
        "recoil_factor": {
            "physical_meaning": "Reduced-mass replacement m_e → μ in |ψ(0)|²",
            "value": "(1 + m_e/m_p)^{-3} = 0.998368",
            "what_it_modulates_in_GeoVac": (
                "|ψ_{1s}(0)|² in the Bohr-Fermi formula. Applied externally; "
                "Track NI cannot produce two-body coordinate coupling natively "
                "(HF-3 verdict)."
            ),
            "paper_34_axes_added": "(variable: mass ratio; dimension: dimensionless; transcendental: none)",
            "paper_34_projection_class": "two-body / center-of-mass / Layer-2 calibration",
            "category_per_HF_taxonomy": "C — external (textbook NRQM)",
            "wired_into_GeoVac": False,
        },
        "Schwinger_a_e": {
            "physical_meaning": "Electron anomalous magnetic moment α/(2π)",
            "value": ALPHA / (6.283185307179586),
            "what_it_modulates_in_GeoVac": (
                "g_e → 2(1 + a_e) in the Bohr-Fermi formula. GeoVac DOES produce "
                "this autonomously (Paper 28 §curvature_coefficients, Track HF-2)."
            ),
            "paper_34_projection_class": "graph-native + Hopf-base measure 1/(2π)",
            "category_per_HF_taxonomy": "A — derivable from existing GeoVac machinery",
            "wired_into_GeoVac": True,
        },
        "synthesis": (
            "The proton structure radii r_charge and r_Zemach are TWO DIFFERENT "
            "Layer-2 focal lengths, both QCD-internal. Track NI implements "
            "r_charge (one-body electronic operator on E_1s) but not r_Zemach. "
            "Both are out of the structural-skeleton scope of GeoVac (the "
            "framework has Wilson SU(3) at the gauge level — Paper 30, ST-SU3 "
            "— but cannot autonomously compute matter-distribution elements "
            "within the (N,0) representations because of the CG obstruction "
            "found in Sprint 5 Track S5). Both require external supply. The "
            "fact that R_PROTON_BOHR is wired into the binding-energy "
            "correction does not transfer to the hyperfine coupling, because "
            "the two corrections live on different focal lengths (charge "
            "distribution vs magnetization distribution) and the framework "
            "carries no operator that distinguishes them."
        ),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> Dict[str, Any]:
    diagnostic = diagnose_track_ni_proton_structure()
    scenarios = predict_three_scenarios()
    verdict = assemble_verdict(diagnostic, scenarios)
    paper34 = paper_34_vocabulary()

    output = {
        "constants": {
            "alpha": ALPHA,
            "a0_fm": A0_FM,
            "R_proton_fm": R_PROTON_FM,
            "R_proton_bohr": R_PROTON_BOHR,
            "R_zemach_fm": R_ZEMACH_FM,
            "R_zemach_bohr": R_ZEMACH_BOHR,
            "Z": Z,
            "A_hf_HF2_baseline_MHz": A_HF_HF2_BASELINE_MHZ,
            "A_hf_experimental_MHz": A_HF_EXPERIMENTAL_MHZ,
            "ratio_R_charge_to_R_zemach": R_PROTON_FM / R_ZEMACH_FM,
        },
        "step1_diagnostic": diagnostic,
        "step2_scenarios": scenarios,
        "step3_verdict": verdict,
        "step4_paper_34_vocabulary": paper34,
    }

    print("=" * 70)
    print("Sprint HF Track 4 -- Zemach correction diagnostic")
    print("=" * 70)
    print()
    print(f"R_proton (charge): {R_PROTON_FM} fm = {R_PROTON_BOHR:.4e} bohr")
    print(f"r_Z (Zemach):      {R_ZEMACH_FM} fm = {R_ZEMACH_BOHR:.4e} bohr")
    print(f"Ratio r_Z / R_charge: {R_ZEMACH_FM / R_PROTON_FM:.4f}")
    print()
    print("Step 1 -- Track NI proton-structure machinery")
    print("-" * 70)
    print(f"  Zemach keyword counts in geovac/nuclear:")
    for kw, ct in diagnostic["zemach_keyword_counts_in_geovac_nuclear"].items():
        print(f"     {kw!r}: {ct}")
    print(f"  Charge-radius keyword counts:")
    for kw, ct in diagnostic["charge_radius_keyword_counts_in_geovac_nuclear"].items():
        print(f"     {kw!r}: {ct}")
    print(f"  hyperfine_coupling_pauli signature: {diagnostic['hyperfine_coupling_pauli_signature']}")
    print(f"  hyperfine takes R_nuc?           {diagnostic['hyperfine_takes_R_nuc']}")
    print(f"  hyperfine takes r_Z?             {diagnostic['hyperfine_takes_r_Z']}")
    print(f"  hyperfine source has size dep.?  {diagnostic['hyperfine_source_has_size_dependence']}")
    print(f"  finite_size_coupling targets HF? {diagnostic['finite_size_targets_hyperfine_energy']}")
    print(f"  -> No magnetization operator in framework: {diagnostic['summary_no_magnetization_operator']}")
    print()
    print("Step 2 -- Three Zemach scenarios")
    print("-" * 70)
    print(f"  Baseline (HF-2 closing):  A_hf = {scenarios['baseline']['A_hf_MHz']:.4f} MHz, residual = {scenarios['baseline']['residual_ppm']:+.2f} ppm")
    print(f"  (A) charge-radius proxy:  dnu/nu = {scenarios['scenario_A_charge_radius_substitution']['delta_nu_using_R_PROTON_ppm']:+.3f} ppm")
    print(f"      A_hf = {scenarios['scenario_A_charge_radius_substitution']['A_hf_predicted_MHz']:.4f} MHz, residual = {scenarios['scenario_A_charge_radius_substitution']['residual_ppm']:+.2f} ppm")
    print(f"  (B) r_Z external 1.045fm: dnu/nu = {scenarios['scenario_B_r_Z_external']['delta_nu_over_nu_ppm']:+.3f} ppm")
    print(f"      A_hf = {scenarios['scenario_B_r_Z_external']['A_hf_predicted_MHz']:.4f} MHz, residual = {scenarios['scenario_B_r_Z_external']['residual_ppm']:+.2f} ppm")
    print(f"  (C) native Zemach (none): A_hf = {scenarios['scenario_C_native_zemach_zero']['A_hf_predicted_MHz']:.4f} MHz, residual = {scenarios['scenario_C_native_zemach_zero']['residual_ppm']:+.2f} ppm")
    print()
    print("Step 3 -- Verdict")
    print("-" * 70)
    print(f"  {verdict['verdict']}")
    print(f"  Headline native A_hf:           {verdict['headline_native_A_hf_MHz']:.4f} MHz, residual = {verdict['headline_native_residual_ppm']:+.2f} ppm")
    print(f"  With r_Z=1.045fm calibrated in: {verdict['with_r_Z_external_A_hf_MHz']:.4f} MHz, residual = {verdict['with_r_Z_external_residual_ppm']:+.2f} ppm")
    print()

    # Write JSON
    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "sprint_hf_track4.json"
    with open(out_path, "w", encoding='utf-8') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"Output written to: {out_path}")

    return output


if __name__ == "__main__":
    main()
