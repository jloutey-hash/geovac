"""Consolidate Step 1+2+3 results into sprint_f1_maxn3_results.json."""
import json
from pathlib import Path

step1 = json.load(open('debug/data/sprint_f1_maxn3_step1_kernel.json'))
step2 = json.load(open('debug/data/sprint_f1_maxn3_step2_fci.json'))
step3 = json.load(open('debug/data/sprint_f1_maxn3_step3_pes.json'))

consolidated = {
    'sprint': 'F1 max_n=3 -- combined W1c x multi-zeta NaH FCI predictions test',
    'date': '2026-05-23',
    'verdict': 'CLEAN NEGATIVE',
    'verdict_summary': (
        'P1 FAILED (R_min at 2.0 bohr, smallest tested = falsification criterion). '
        'P2 FAILED (D_e=+0.71 Ha = ~10x experimental, spurious-binding signature). '
        'P3 PASSED (mz differential 0.21 Ha within predicted [0.02, 0.30] range). '
        'Substantive new content: dominant natural orbital IS now a true bonding '
        'combination (Na 3s + H 1s 50/50 at R=3.5) but with energy higher than '
        'separated configuration at every R. Bonding orbital is constructible '
        '(structural success) but cannot lower energy (energetic failure). '
        'Sub-layer 3 reclassifies from basis-closable cross-shift to genuine endomorphism.'
    ),
    'predictions': {
        'P1_internal_minimum': {
            'predicted': 'R_min in [3.0, 4.5] bohr (experimental NaH R_eq=3.566)',
            'actual': f"R_min = {step3['R_min']} bohr (smallest tested)",
            'verdict': 'FAIL',
        },
        'P2_binding_energy': {
            'predicted': 'D_e in [0.0375, 0.150] Ha (within 2x of experimental 0.075)',
            'actual': f"D_e = +{step3['D_e_pes_Ha']:.4f} Ha (PES descent depth, spurious)",
            'verdict': 'FAIL (conditional on P1)',
        },
        'P3_multizeta_differential': {
            'predicted': '|E_W1c - E_W1c+mz| in [0.02, 0.30] Ha at R_eq',
            'actual': f"|mz_diff at R=2.0| = {step3['predictions']['P3_mz_diff_in_0p02_to_0p30']['actual_mz_diff_Ha']:.4f} Ha",
            'verdict': 'PASS',
        },
    },
    'system': {
        'name': 'NaH',
        'max_n': 3,
        'Q_total': 56,
        'M_spatial': 28,
        'n_electrons': 2,
        'fci_dim': 784,
        'experimental_R_eq_bohr': 3.566,
        'experimental_D_e_Ha': 0.075,
    },
    'step1_kernel_differential': {
        'description': 'Cross-V_ne diagonal on Na 3s at R=3.566 (mz vs no_mz)',
        'h1_cv_na3s_bare_no_mz_Ha': -0.279403,
        'h1_cv_na3s_bare_mz_Ha': -0.226047,
        'differential_bare_Ha': +0.053356,
        'differential_w1c_Ha': +0.053356,
        'comparison_to_alpha_pes_step1': 'Same +0.05 Ha order as F1-P1+P2 at max_n=2',
        'sanity_check': 'PASS — multi-zeta machinery works correctly at max_n=3',
    },
    'step2_single_point_fci': {
        'R_eq_results': step2['arch_results']['R_eq'],
        'R_diss_results': step2['arch_results']['R_diss'],
        'D_e_2point': step2['D_e_2point'],
        'multi_zeta_differential': step2['multi_zeta_differential'],
        'bonding_signature': step2['bonding_signature'],
        'decision_gate': step2['decision_gate'],
    },
    'step3_mini_pes': {
        'R_values_bohr': step3['all_R_sorted'],
        'E_combined_by_R': step3['E_combined_by_R'],
        'E_w1c_alone_by_R': step3['E_w1c_alone_by_R'],
        'top_natural_occupation_by_R': step3['top_no_combined_by_R'],
        'dom_NO_NaH_split_by_R': step3['dom_no_NaH_split_combined'],
        'R_min': step3['R_min'],
        'E_min': step3['E_min'],
        'E_diss': step3['E_diss_R_10'],
        'D_e_pes': step3['D_e_pes_Ha'],
        'is_internal_min': step3['is_internal_min'],
        'verdict': step3['verdict'],
    },
    'structural_findings': [
        ('STRUCTURAL SUCCESS at max_n=3: the dominant natural orbital of W1c+mz '
         'at R=3.5 is a true bonding combination Na 3s + H 1s with 50/50 '
         'amplitude split (mixing coefficients -0.698 Na 3s + -0.687 H 1s). '
         'At max_n=2 the dominant NO was H-dominant (Na amp2 = 0). The basis '
         'enlargement DID provide the orbital mixing that max_n=2 lacked.'),
        ('ENERGETIC FAILURE: the constructed bonding orbital has higher energy '
         'than the separated configuration at every tested R. PES is '
         'monotonically descending across R in {2.0, 2.5, 3.0, 3.5, 4.0, 5.0, '
         '7.0, 10.0}. R_min = 2.0 bohr (smallest tested) = exactly the '
         'falsification criterion.'),
        ('NATURAL OCCUPATIONS remain [1.0, 1.0] across all R (open-shell '
         'singlet of bonding + antibonding singly occupied), NOT a closed-'
         'shell bond [2.0, 0.0]. The 2nd natural orbital is the antibonding '
         'combination -0.698 Na 3s + +0.687 H 1s.'),
        ('Multi-zeta is LOAD-BEARING at max_n=3 (0.21 Ha differential at '
         'R=2.0, 0.06 Ha at R=3.5, bit-zero at R=10). Mz differential '
         'magnitude similar to max_n=2 (F1-P1+P2 reported 0.21 at R=2.0, '
         '0.056 at R=3.5). Endomorphism classification preserved across '
         'basis sizes.'),
        ('W1c-alone PES descent depth at max_n=3: 0.916 Ha. Combined '
         'W1c+mz: 0.710 Ha (~22% reduction). At max_n=2 the reduction was '
         '23% (F1-P1+P2). Multi-zeta contribution is stable across basis.'),
        ('SUB-LAYER 3 RECLASSIFICATION: not basis-closable cross-shift, but '
         'genuine endomorphism class. The bonding orbital is structurally '
         'constructible (cross-shift activates as predicted), but the '
         'energetic gap that would bring the bonding combination below the '
         'separated configuration requires content the framework cannot '
         'generate from inside (likely cross-V_ne kernel SHAPE on the '
         'partner side, per Track 3 named target P2).'),
    ],
    'next_sprint_recommendations': [
        ('Paper edits (LATER, gated on PI review): Paper 19 sec:w1c_residual '
         'extension documenting sub-layer 3 reclassification; CLAUDE.md §3 '
         'dead ends row for basis-closable hypothesis falsification at max_n=3; '
         'CLAUDE.md §1.7 M-Z partition update (keep 2-bucket framing, drop '
         'proposed 3-bucket refinement).'),
        ('Pivot to F2 named target: cross-V_ne kernel-shape substitution on '
         'the partner sub-block (replace bare cross-V_ne integration with '
         'physical-n hydrogenic shape). Track 3 diagnostic (2026-05-09) '
         'identified this as the load-bearing residual after the diagonal '
         'h1 substitution leg of W1c was exhausted. Estimated 1-2 weeks.'),
        ('Alternative: pivot to a different system class (LiH+, HCl, MgH2). '
         'The W1c-residual wall is empirically Z-decreasing; at higher Z the '
         'residual is smaller fraction of bond energy and partial closure '
         'may be sufficient for physical binding.'),
        ('Do NOT pursue further max_n enlargement: at max_n=3 the structural '
         'gap (bonding orbital constructibility) is closed; further basis '
         'enlargement would not address the energetic gap that the '
         'CLEAN-NEGATIVE PES reveals.'),
    ],
}

with open('debug/data/sprint_f1_maxn3_results.json', 'w') as f:
    json.dump(consolidated, f, indent=2)

print(f"Consolidated results saved to debug/data/sprint_f1_maxn3_results.json")
print(f"Verdict: {consolidated['verdict']}")
print(f"  P1: {consolidated['predictions']['P1_internal_minimum']['verdict']}")
print(f"  P2: {consolidated['predictions']['P2_binding_energy']['verdict']}")
print(f"  P3: {consolidated['predictions']['P3_multizeta_differential']['verdict']}")
