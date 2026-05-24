"""Append F1 max_n=3 verification to the predictions JSON."""
import json
from pathlib import Path

# Load predictions
preds_path = Path('debug/data/sprint_w1c_mz_partition_predictions.json')
with open(preds_path) as f:
    preds = json.load(f)

# Load Step 3 results
step3_path = Path('debug/data/sprint_f1_maxn3_step3_pes.json')
with open(step3_path) as f:
    step3 = json.load(f)

# Append verification section
preds['verification'] = {
    'sprint': 'F1 max_n=3 -- Sprint of 2026-05-23 evening',
    'date': '2026-05-23',
    'verdict': 'CLEAN NEGATIVE',
    'verdict_summary': (
        'P1 FAILED (R_min at smallest tested R=2.0 bohr, exactly the '
        'falsification criterion). P2 FAILED (conditional on P1). '
        'P3 PASSED (multi-zeta differential within predicted [0.02, 0.30] Ha '
        'range at R=2.0). The basis-closable cross-shift hypothesis '
        'for sub-layer 3 is REFUTED at max_n=3; sub-layer 3 reclassifies '
        'as a genuine endomorphism class (cannot be closed by basis '
        'enlargement at this scale).'
    ),
    'prediction_1': {
        'predicted_range_bohr': preds['predictions_F1_max_n_3']['prediction_1']['statement'],
        'actual_R_min_bohr': step3['R_min'],
        'falsifier_criterion': 'R_min at smallest tested R',
        'verdict': 'FAIL',
        'detail': (
            'R_min = 2.0 bohr (smallest tested) at max_n=3 under combined '
            'W1c+multi-zeta architecture. PES monotonically descending '
            'across R in {2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0}. '
            'NO internal equilibrium minimum.'
        ),
    },
    'prediction_2': {
        'predicted_range_Ha': [0.0375, 0.150],
        'actual_D_e_Ha': step3['D_e_pes_Ha'],
        'verdict': 'FAIL (conditional on P1)',
        'detail': (
            'D_e (PES, smallest-R to dissociation) = +0.7097 Ha -- about '
            '10x larger than the experimental NaH 0.075 Ha, indicating '
            'continued over-attraction at small R. This is the '
            'spurious-binding signature, not a real D_e.'
        ),
    },
    'prediction_3': {
        'predicted_range_Ha': [0.02, 0.30],
        'actual_mz_diff_Ha': abs(step3['predictions']['P3_mz_diff_in_0p02_to_0p30']['actual_mz_diff_Ha']),
        'verdict': 'PASS',
        'detail': (
            '|E_W1c - E_W1c+mz| at R=R_min = 2.0 bohr is 0.2066 Ha, '
            'within the predicted [0.02, 0.30] range. Multi-zeta is '
            'load-bearing at max_n=3 (consistent with F1-P1+P2 finding '
            'at max_n=2). The endomorphism classification of multi-zeta '
            'is preserved across basis sizes.'
        ),
    },
    'structural_findings': [
        ('CRITICAL POSITIVE: at max_n=3 the dominant natural orbital of '
         'combined W1c+mz at R=3.5 shows true bonding combination '
         '(Na 3s + H 1s, ~50/50 amplitude split with mixing '
         'coefficients -0.698 Na 3s + -0.687 H 1s). At max_n=2 it was '
         'H-dominant (1.0 H amplitude). Larger basis DID enable '
         'bonding-combination construction. The 2nd natural orbital is '
         'the antibonding combination (-0.698 Na 3s + +0.687 H 1s).'),
        ('BUT: natural occupations remain [1.0, 1.0] (not [2.0, 0.0]), '
         'so the FCI is still an open-shell singlet of two singly-'
         'occupied orbitals (bonding and antibonding), NOT a closed-'
         'shell bonded pair.'),
        ('PES is monotonically descending across full range '
         'R in {2.0, ..., 10.0} bohr; no internal minimum at any '
         'tested R. PES descent depth = 0.71 Ha at R=2.0 vs R=10.0 '
         '(combined W1c+mz).'),
        ('Multi-zeta differential R-dependence: smoothly decays from '
         '-0.207 Ha at R=2.0 to bit-zero at R=10.0 (same R-dependent '
         'pattern as F1-P1+P2 max_n=2).'),
        ('W1c-alone PES descent depth at max_n=3: 0.916 Ha. Combined '
         'W1c+mz: 0.710 Ha. Multi-zeta reduces descent depth by ~22% '
         'at max_n=3, similar to ~23% at max_n=2.'),
        ('Sub-layer 3 RECLASSIFIES: NOT basis-closable cross-shift, '
         'but GENUINE ENDOMORPHISM CLASS that the framework cannot '
         'close from inside its bimodule machinery.'),
    ],
    'sub_layer_3_reclassification': {
        'previous_classification': 'BIMODULE CROSS-SHIFT (basis-closable)',
        'new_classification': 'MODULE ENDOMORPHISM (basis-irreducible at tested scales)',
        'justification': (
            'Sprint F1 max_n=3 tested whether the basis enlargement from '
            'max_n=2 to max_n=3 (Q=20 to Q=56) provides the dimensional '
            'richness for the FCI to construct a true bonding combination '
            'with energy lower than the separated configuration. The result '
            'is mixed: the FCI DID construct a true bonding combination as '
            'the dominant natural orbital (Na 3s + H 1s, 50/50 amplitude '
            'split at R=3.5), but the bonding combination still has higher '
            'energy than the separated orbital configuration at every '
            'tested R. The PES is still monotonically descending; no internal '
            'minimum emerges. The basis enlargement closed the STRUCTURAL '
            'gap (bonding orbital is now constructible) but the ENERGETIC '
            'gap remains. Sub-layer 3 is therefore not basis-closable in '
            'the sense of producing binding-with-equilibrium; it requires '
            'a deeper closure mechanism beyond what the framework can '
            'generate from within its bimodule machinery at the tested scales.'
        ),
        'updated_taxonomy': (
            '2-bucket M-Z partition unchanged: cross-shift (handled by '
            'framework, sub-layers 1 + sub-layer 3 substructure) vs '
            'endomorphism (external input, sub-layers 2 + reclassified '
            'sub-layer 3). Sub-layer 3 was originally hypothesized as '
            'basis-closable cross-shift; F1 max_n=3 falsifies this and '
            'confirms sub-layer 3 lives in the same endomorphism class as '
            'sub-layer 2 in the sense of being external to framework-side '
            'closure.'
        ),
        'consistency_with_M_Z_original': (
            'STRONG -- the M-Z 2-bucket partition (cross-shift / '
            'endomorphism) holds; the proposed 3-bucket refinement '
            '(cross-shift / basis-closable cross-shift / endomorphism) is '
            'not supported. The original W3 and LS-8a wall framing extends '
            'cleanly: there is calibration data the framework needs that '
            'it cannot generate from inside, even at modest basis enlargement.'
        ),
        'orbital_structure_at_max_n_3': (
            'At max_n=3 the framework constructs the right bonding orbital '
            '(structural success) but cannot lower its energy below the '
            'separated configuration (energetic failure). The two-singly-'
            'occupied natural occupation pattern [1.0, 1.0] persists from '
            'max_n=2 even though the orbital shapes change. This is a deeper '
            'finding than just "basis too small" -- it indicates the '
            'framework needs additional content (Track 3 P2: cross-V_ne '
            'kernel shape, or a structurally different solver) that is not '
            'accessible via basis enlargement on the current architecture.'
        ),
    },
    'next_sprint_recommendations': [
        ('Apply paper edits per CLEAN NEGATIVE branch: Paper 19 '
         'sec:w1c_residual extension documenting the sub-layer 3 '
         'reclassification; CLAUDE.md §3 dead ends row for '
         'basis-closable hypothesis falsification.'),
        ('Pivot to Track 3 named target P2 (cross-V_ne kernel-shape '
         'substitution): replace the bare cross-V_ne kernel integration '
         'on Na valence with physical-n hydrogenic shape (from FrozenCore '
         'solver) rather than just diagonal eigenvalue substitution. This '
         'addresses a structurally distinct mechanism from the orbital-pair '
         'flexibility tested here.'),
        ('Alternative: pivot to a different system class (LiH+, HCl, '
         'MgH2) where the W1c-residual wall is empirically smaller and a '
         'partial closure may give physical binding even without complete '
         'W1c removal.'),
        ('Do NOT pursue further max_n enlargement on NaH: at max_n=3 the '
         'dominant bonding orbital is already constructible, so the '
         'structural gap is closed at this basis. Further enlargement '
         'would not change the wall character.'),
    ],
}

with open(preds_path, 'w') as f:
    json.dump(preds, f, indent=2)

print('Verification appended to predictions JSON')
print(f"verdict: {preds['verification']['verdict']}")
print(f"P1: {preds['verification']['prediction_1']['verdict']}")
print(f"P2: {preds['verification']['prediction_2']['verdict']}")
print(f"P3: {preds['verification']['prediction_3']['verdict']}")
