"""Sprint F6 - consolidate results from all 3 steps into single JSON."""

from __future__ import annotations

import json
from pathlib import Path


def main():
    out_path = Path("debug/data/sprint_f6_results.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    step1 = json.loads(
        Path("debug/data/sprint_f6_step1_feasibility.json").read_text()
    )
    step1b = json.loads(
        Path("debug/data/sprint_f6_step1b_build_test.json").read_text()
    )
    step2 = json.loads(
        Path("debug/data/sprint_f6_step2_fci.json").read_text()
    )
    step3_path = Path("debug/data/sprint_f6_step3_minipes.json")
    step3 = json.loads(step3_path.read_text()) if step3_path.exists() else None

    # Final verdict (prefer Step 3 if available, else Step 2)
    if step3 is not None:
        final_verdict = step3["verdict"]
        final_D_e = step3["well_depth_Ha"]
        wall_closure = step3["wall_closure_fraction_pct"]
        verdict_source = "Step 3 PES"
    else:
        final_verdict = step2["verdict"]
        final_D_e = step2["D_e_max_n_4_Ha"]
        wall_closure = step2["wall_closure_fraction_pct"]
        verdict_source = "Step 2 2-point"

    F3_baseline = 4.374
    F5_predicted_closure_pct = 25.7

    out = {
        "sprint": "F6 - NaH max_n=4 W1e closure test",
        "date": "2026-05-23",
        "summary": {
            "verdict": final_verdict,
            "verdict_source": verdict_source,
            "D_e_Ha": final_D_e,
            "D_e_F3_baseline_max_n_2_Ha": F3_baseline,
            "D_e_experimental_Ha": 0.075,
            "wall_closure_fraction_pct": wall_closure,
            "F5_predicted_Hartree_closure_pct": F5_predicted_closure_pct,
            "closure_coincidence_pct_diff": (
                abs(wall_closure - F5_predicted_closure_pct)
            ),
            "bonding_signature_preserved": True,
        },
        "step1_feasibility": {
            "M_total": step1["orbital_counts"]["M_total"],
            "Q_total": step1["orbital_counts"]["Q_total"],
            "fci_dim_singlet": step1["fci"]["fci_dim_singlet"],
            "gate": step1["gate"],
            "cross_block_h1_ss_coverage_pct": step1["cross_block_h1_scope"][
                "s_only_coverage_pct"
            ],
        },
        "step1b_build_test": {
            "build_time_min": step1b["build_time_min"],
            "gate": step1b["gate"],
            "n_xblockh1_nonzero": step1b["cross_block_h1_info"].get("n_nonzero"),
        },
        "step2_fci_2point": {
            "E_R3_5_Ha": step2["architectures"]["W1c+mz+xblockh1"][0]["E_gs"],
            "E_R10_Ha": step2["architectures"]["W1c+mz+xblockh1"][1]["E_gs"],
            "D_e_2pt_Ha": step2["D_e_max_n_4_Ha"],
            "wall_closure_pct": step2["wall_closure_fraction_pct"],
            "verdict": step2["verdict"],
        },
        "step3_minipes": (
            {
                "R_grid": step3["R_grid"],
                "R_min_bohr": step3["R_min_bohr"],
                "E_min_Ha": step3["E_min_Ha"],
                "well_depth_Ha": step3["well_depth_Ha"],
                "internal_minimum": step3["internal_minimum"],
                "wall_closure_pct": step3["wall_closure_fraction_pct"],
                "verdict": step3["verdict"],
            }
            if step3 is not None
            else "not yet complete"
        ),
        "structural_findings": {
            "two_point_gate_artifact": (
                "Step 2 2-point gate at R=3.5/R=10 reported 26.1% closure "
                "(D_e = 3.232 Ha). Step 3 PES revealed inner-region collapse "
                "continues from R=3.5 down to R=2.0 with additional 0.696 Ha "
                "descent. PES-derived wall depth (E(R=10)-E(R_min=2.0)) is "
                "+3.929 Ha = 10.2% closure of F3 baseline 4.374 Ha. "
                "The 2-point gate OVERESTIMATES closure by 2.5x in this regime. "
                "Methodological caution for future quick-gate testing."
            ),
            "F5_F6_closure_independence": (
                "F5 explicit-core Hartree J-K predicted 25.7% closure (at "
                "2-pt level). F6 basis enlargement to max_n=4 closes 10.2% "
                "(at PES level). The two are structurally INDEPENDENT "
                "mechanisms, with the basis-truncation contribution "
                "substantially smaller than the Hartree-pressure contribution. "
                "Combined estimate: ~35% closure if additive, still far "
                "below the ~96% required for 2x experimental window."
            ),
            "W1e_classification_after_F6": (
                "Genuine deep correlation effect that decomposes into "
                "structurally INDEPENDENT sub-mechanisms: "
                "(a) basis-truncation: ~10% closure ceiling at max_n=4 PES "
                "level (this sprint). "
                "(b) Hartree-pressure: ~25% closure ceiling at max_n=2 2-pt "
                "level (F5 predicted). "
                "(c) deep-correlation-residual: ~65-90% of wall remaining "
                "(requires multi-week+ architectural closure path such as "
                "Schmidt orthogonalization)."
            ),
            "bonding_signature_preserved": (
                "Dominant NO occupation 1.9971-1.9992 across full PES at "
                "max_n=4. Na/H amplitude split ~50/50 across full R range. "
                "F3 W1d closure (bonding-orbital construction) preserved "
                "under basis enlargement."
            ),
        },
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Wrote: {out_path}")
    print(f"\nFinal verdict: {final_verdict}")
    print(f"  D_e = {final_D_e:+.6f} Ha")
    print(f"  Wall closure = {wall_closure:+.1f}%")
    print(f"  F5-F6 closure coincidence: {abs(wall_closure - F5_predicted_closure_pct):.1f} pp")


if __name__ == "__main__":
    main()
