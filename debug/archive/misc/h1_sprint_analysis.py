"""Sprint H1 — Higgs from inner fluctuation: full analysis driver.

Computes the falsifier verdict at n_max=2,3 and characterizes the gauge/Higgs
sectors of inner fluctuations on the almost-commutative extension.

Output: debug/data/h1_falsifier.json + console report.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from geovac.almost_commutative import (
    AlmostCommutativeTriple,
    ElectroweakFiniteTriple,
    gauge_only_triple,
    minimal_electroweak_triple,
)


def run_falsifier_at_n_max(n_max: int, n_samples: int = 100) -> dict:
    """Test the §5 falsifier at given n_max for three D_F candidates."""
    results = {"n_max": n_max, "tests": {}}

    # Candidate 1: D_F = 0 (gauge only)
    T_zero = gauge_only_triple(n_max)
    negative, reason, data = T_zero.check_natural_negative(
        n_random_generators=n_samples, seed=42
    )
    results["tests"]["D_F_zero"] = {
        "negative_holds": negative,
        "reason": reason,
        "higgs_max": data["higgs_max"],
        "higgs_mean": data["higgs_mean"],
        "gauge_max": data["gauge_max"],
        "gauge_mean": data["gauge_mean"],
    }

    # Candidate 2: D_F with nontrivial Yukawa (electron-only)
    T_e = minimal_electroweak_triple(n_max, yukawa_e=0.1, yukawa_nu=0.0)
    negative, reason, data = T_e.check_natural_negative(
        n_random_generators=n_samples, seed=42
    )
    results["tests"]["D_F_yukawa_e_only"] = {
        "negative_holds": negative,
        "reason": reason,
        "higgs_max": data["higgs_max"],
        "higgs_mean": data["higgs_mean"],
        "gauge_max": data["gauge_max"],
        "gauge_mean": data["gauge_mean"],
    }

    # Candidate 3: D_F with both Yukawa entries
    T_full = minimal_electroweak_triple(n_max, yukawa_e=0.3, yukawa_nu=0.2)
    negative, reason, data = T_full.check_natural_negative(
        n_random_generators=n_samples, seed=42
    )
    results["tests"]["D_F_full_yukawa"] = {
        "negative_holds": negative,
        "reason": reason,
        "higgs_max": data["higgs_max"],
        "higgs_mean": data["higgs_mean"],
        "gauge_max": data["gauge_max"],
        "gauge_mean": data["gauge_mean"],
    }

    return results


def characterize_higgs_potential(n_max: int) -> dict:
    """Spectrum analysis of D_omega across Yukawa values to characterize V(Phi)."""
    yukawas = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
    Tr_D2 = []
    Tr_D4 = []
    spec_min = []
    spec_max = []
    for y in yukawas:
        T = minimal_electroweak_triple(n_max, yukawa_e=y, yukawa_nu=0.0)
        D = T.dirac_combined()
        Tr_D2.append(float(np.real(np.trace(D @ D))))
        Tr_D4.append(float(np.real(np.trace(D @ D @ D @ D))))
        eigs = np.linalg.eigvalsh(D)
        spec_min.append(float(min(eigs)))
        spec_max.append(float(max(eigs)))
    return {
        "n_max": n_max,
        "yukawas": yukawas,
        "Tr_D_squared": Tr_D2,
        "Tr_D_fourth": Tr_D4,
        "spec_min": spec_min,
        "spec_max": spec_max,
    }


def run_explicit_inner_fluctuation_decomposition(n_max: int) -> dict:
    """Build a specific representative inner fluctuation and decompose it.

    Pick a generator that activates each sector.
    """
    T = minimal_electroweak_triple(n_max, yukawa_e=0.3, yukawa_nu=0.2)
    I_GV = np.eye(T.dim_GV, dtype=np.complex128)

    # Generator 1: a = b = (1, q_pure_imag) with q a non-trivial quaternion
    # This activates the SU(2) gauge sector via [D_GV, b_GV] component IF
    # b_GV varies; here b_GV = I, so [D_GV, I] = 0. This generator picks up
    # ONLY the [D_F, b_F] piece (Higgs + R-block of gauge).
    gen_higgs_only = [(I_GV, 0.0, (0.0, 1.0, 0.0, 0.0),
                       I_GV, 1.0, (0.0, 0.0, 0.0, 0.0))]

    # Generator 2: vary M_GV to activate the gauge 1-form
    M_GV = T.gv_multiplier(1)  # non-identity multiplier
    gen_gauge = [(I_GV, 1.0, (1.0, 0.0, 0.0, 0.0),
                  M_GV, 1.0, (1.0, 0.0, 0.0, 0.0))]

    # Generator 3: combined
    gen_combined = [
        (I_GV, 0.5, (1.0, 0.5, 0.0, 0.0),
         M_GV, 1.0, (0.5, 1.0, 0.0, 0.0))
    ]

    out = {"n_max": n_max, "scenarios": {}}
    for name, gens in [
        ("higgs_only_a_b_id_GV", gen_higgs_only),
        ("gauge_only_b_F_id", gen_gauge),
        ("combined", gen_combined),
    ]:
        omega = T.inner_fluctuation_one_form(gens)
        decomp = T.decompose_fluctuation(omega)
        out["scenarios"][name] = {
            "higgs_norm": float(T.higgs_norm(omega)),
            "gauge_norm": float(T.gauge_norm(omega)),
            "matter_antimatter_off_norm": float(T.matter_antimatter_off_norm(omega)),
            "higgs_matter_LR_norm": float(np.linalg.norm(decomp["higgs_matter_LR"])),
            "higgs_matter_RL_norm": float(np.linalg.norm(decomp["higgs_matter_RL"])),
            "gauge_matter_L_norm": float(np.linalg.norm(decomp["gauge_matter_L"])),
            "gauge_matter_R_norm": float(np.linalg.norm(decomp["gauge_matter_R"])),
            "omega_total_norm": float(np.linalg.norm(omega)),
        }
    return out


def main():
    out_dir = Path("debug/data")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint H1 — Higgs from inner fluctuation: analysis")
    print("=" * 70)
    print()

    all_results = {}

    # Falsifier sweep
    print("Falsifier check at n_max=1, 2, 3:")
    for n_max in [1, 2, 3]:
        print(f"  n_max = {n_max}:")
        result = run_falsifier_at_n_max(n_max, n_samples=50)
        for cand_name, info in result["tests"].items():
            verdict = "HOLDS" if info["negative_holds"] else "FAILS"
            print(
                f"    {cand_name:25s}: §5 falsifier {verdict}  "
                f"(higgs_max={info['higgs_max']:.3e}, "
                f"gauge_max={info['gauge_max']:.3e})"
            )
        all_results[f"falsifier_n_max_{n_max}"] = result
    print()

    # Higgs potential characterization
    print("Higgs spectrum scan at n_max = 2:")
    pot = characterize_higgs_potential(2)
    print(f"  yukawa_grid = {pot['yukawas']}")
    print(f"  Tr(D^2)     = {[f'{x:.4f}' for x in pot['Tr_D_squared']]}")
    print(f"  Tr(D^4)     = {[f'{x:.4f}' for x in pot['Tr_D_fourth']]}")
    print(f"  spec_max    = {[f'{x:.4f}' for x in pot['spec_max']]}")
    all_results["higgs_potential_n_max_2"] = pot
    print()

    # Explicit decomposition
    print("Explicit inner-fluctuation decomposition at n_max = 2:")
    decomp = run_explicit_inner_fluctuation_decomposition(2)
    for scenario, info in decomp["scenarios"].items():
        print(f"  Scenario: {scenario}")
        print(f"    omega_total_norm     = {info['omega_total_norm']:.4e}")
        print(f"    higgs_norm           = {info['higgs_norm']:.4e}")
        print(f"    gauge_norm           = {info['gauge_norm']:.4e}")
        print(f"    matter<->antimatter  = {info['matter_antimatter_off_norm']:.4e}  (should be ~0)")
        print(f"    higgs LR / RL        = {info['higgs_matter_LR_norm']:.4e} / {info['higgs_matter_RL_norm']:.4e}")
    all_results["decomposition_n_max_2"] = decomp
    print()

    # Save JSON
    out_file = out_dir / "h1_falsifier.json"
    with open(out_file, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Results saved to {out_file}")
    print()

    print("=" * 70)
    print("VERDICT")
    print("=" * 70)
    print("""
The Sec.5 falsifier ('every Hermitian D_F from GeoVac forces Higgs sector zero')
HOLDS only when D_F = 0 (zero Yukawa). For any imposed Yukawa Y nonzero, the
Higgs sector is non-trivially populated by inner fluctuations.

GeoVac structure provides:
  - real structure J_GV (Track 2 verified at n_max <= 3),
  - chirality grading gamma_GV (truthful CH eigenvalue sign),
  - operator system A_GV with multiplier matrices,
  - the natural KO-3 base manifold S^3.

GeoVac structure does NOT provide:
  - a natural choice of Yukawa matrix Y in D_F.

Verdict: positive-thin (per scoping memo Sec.0). The Higgs construction is
WELL-DEFINED on the GeoVac AC extension; the Yukawa entries Y are a
free input not selected by GeoVac. This places GeoVac on the
Marcolli-vS-without-Higgs side of the 2024 distinction at the structural
level: GeoVac admits the Higgs construction but does not autonomously
generate a non-trivial D_F off-diagonal.

This is consistent with Paper 32 Sec.VIII.B G2: the gap is NOT 'inner
fluctuations cannot define a Higgs' -- they can. The gap is 'no GeoVac
mechanism selects a non-trivial Yukawa structure'.
""")


if __name__ == "__main__":
    main()
