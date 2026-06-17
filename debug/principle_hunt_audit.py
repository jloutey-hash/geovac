"""Sprint C2 principle hunt audit.

Tests five candidate discriminator principles against the 60-entry
forced/free seam catalogue (docs/forced_free_seam.md):

  P1 — multi-focal depth                (Paper 57 §5.1, MF axis)
  P2 — period-content classifiability   (Paper 57 §5.2, P axis)
  P3 — dimensional character            (Paper 57 §5.3, Dim axis)
  P4 — compactness inheritance          (Paper 57 §5.4, refined: Peter-Weyl
                                          derivable AND no continuum input)
  P5 — packing-reachability             (NEW; the candidate family-2 axis)

  Plus conjunctions: P1 ∧ P5, P1 ∨ P5, etc.

For each principle, count true positives, false positives, true negatives,
false negatives, and per-family breakdowns. Output JSON to
debug/data/principle_hunt_audit.json.
"""

import json
from dataclasses import dataclass, field, asdict
from typing import Callable
from pathlib import Path


# =====================================================================
# CATALOGUE
# =====================================================================

@dataclass
class Entry:
    id: str
    name: str
    status: str            # 'F' / 'A' / 'C'
    mf: int                # 1 or 2 (2 means MF > 1)
    p: str                 # 'M1' / 'M2' / 'M3' / 'M1+M2' / 'M1/M2/M2' / 'outside' / 'none' / 'meta' / 'N/A'
    dim: str               # 'dimensionless' / 'dimensionful' / 'N/A'
    packing_reachable: str # 'yes' / 'no' / 'conditional'
    family: str            # 'none' / 'multi_focal_composition' / 'inner_factor_input_data' / 'admitted_only'
    domain: str
    note: str = ""


# Each entry tagged on all axes. packing_reachable values are justified
# by the witness-derivation rule (yes), explicit no-go theorems (no), or
# conditional adoption of Upgrade B / DAS (conditional).

ENTRIES = [
    # A. Foundations / skeleton (all F, all packing-reachable from Paper 0)
    Entry('A1', 'Graph Laplacian spectrum n^2-1', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'skeleton', 'direct from packing axiom'),
    Entry('A2', 'S^3 as forced outer manifold', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'skeleton', 'Bertrand + SO(4) closure'),
    Entry('A3', 'kappa = -1/16 universal kinetic scale', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'skeleton', 'Fock projection exact rational'),
    Entry('A4', 'Shell degeneracy g_n = n^2', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'skeleton', 'representation theory'),
    Entry('A5', 'Discreteness <-> compactness', 'F', 1, 'N/A', 'N/A',
          'yes', 'none', 'skeleton', 'Peter-Weyl organizing observation'),

    # B. Gauge structure
    Entry('B1', 'Gauge group U(1)xSU(2)xSU(3)', 'F', 1, 'M1', 'dimensionless',
          'yes', 'none', 'gauge', 'Bertrand x Hopf-tower n<=3 (Paper 32 VIII.B)'),
    Entry('B2', 'Inner algebra factor count = 3', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'inner', 'Hurwitz + Hopf-tower (Door 4b/4d)'),
    Entry('B3', 'C factor at n=1', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'inner', 'unique div algebra of S^1'),
    Entry('B4', 'M_3(C) factor at n=3', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'inner', 'Hurwitz fallback at S^5'),
    Entry('B5', 'H vs M_2(C) at n=2', 'A', 1, 'none', 'dimensionless',
          'conditional', 'admitted_only', 'inner',
          'Door 4c NEGATIVE on J-sign-table; closed by Upgrade B (Door 4e)'),
    Entry('B6', 'dim M(D_F) = 128 per generation', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'inner', 'Forced-Count Theorem (Paper 32 VIII)'),
    Entry('B7', 'Higgs admission (Mexican-hat allowed)', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'inner', 'Sprint H1 POSITIVE-THIN'),
    Entry('B8', 'Gauge lower-bound forcing (saturated)', 'C', 1, 'none', 'dimensionless',
          'no', 'inner_factor_input_data', 'gauge',
          'Bertrand x Hopf gives upper bound only'),

    # C. Spectral / period content (all F, all packing-reachable via theorems)
    Entry('C1', 'Master Mellin engine M[Tr(D^k e^{-tD^2})] k in {0,1,2}', 'F', 1, 'meta', 'dimensionless',
          'yes', 'none', 'periods', 'case-exhaustion theorem'),
    Entry('C2', 'M1 prefactor Vol(S^2)/4 = pi', 'F', 1, 'M1', 'dimensionless',
          'yes', 'none', 'periods', 'Hopf-base measure signature'),
    Entry('C3', 'Seeley-DeWitt a_0 = a_1 = sqrt(pi) on S^3', 'F', 1, 'M2', 'dimensionless',
          'yes', 'none', 'periods', 'Paper 28 T9 closed form'),
    Entry('C4', 'zeta_{D^2}(s) = 2^{2s-1}[lambda(2s-2)-lambda(2s)]', 'F', 1, 'M2', 'dimensionless',
          'yes', 'none', 'periods', 'Paper 28 T9'),
    Entry('C5', 'chi_{-4} identity D_even - D_odd = 2^{s-1}(beta(s)-beta(s-2))', 'F', 1, 'M3', 'dimensionless',
          'yes', 'none', 'periods', 'Paper 28 IV vertex-parity Hurwitz'),
    Entry('C6', 'M3 subset MT(Z[i,1/2]) level <= 4', 'F', 1, 'M3', 'dimensionless',
          'yes', 'none', 'periods', 'Sprint Mixed-Tate; Paper 55'),
    Entry('C7', 'Universal propinquity rate 4/pi', 'F', 1, 'M1', 'dimensionless',
          'yes', 'none', 'oa', 'WH1 PROVEN; Paper 38 + Paper 40 rank-invariant'),
    Entry('C8', 'Pythagorean <H_local, D_W^L>_HS = 0', 'F', 1, 'M1', 'dimensionless',
          'yes', 'none', 'oa', 'Paper 43 10.2 bit-exact'),
    Entry('C9', 'S_min irreducible 200+ dps', 'F', 2, 'outside', 'dimensionless',
          'yes', 'none', 'periods',
          'NOT in {1,pi^2,pi^4,G,beta(4),zeta(3),zeta(5)}; derivable from Paper 0 but lives outside Mellin engine'),

    # D. Gravity
    Entry('D1', 'Spectral action two-term-exact on S^3', 'F', 1, 'M2', 'dimensionless',
          'yes', 'none', 'gravity', 'Bernoulli identity zeta_unit(-k)=0'),
    Entry('D2', 'S_BH = A Lambda^2/(12 pi); G_N = 3 pi/Lambda^2', 'F', 1, 'M1', 'dimensionful',
          'yes', 'none', 'gravity', 'Paper 51 G4 closed form (dimensionful via external Lambda)'),
    Entry('D3', 'Cone coefficient -1/12 (Sommerfeld-Cheeger)', 'F', 1, 'M2', 'dimensionless',
          'yes', 'none', 'gravity', 'Paper 51 G4-4c bit-exact'),
    Entry('D4', 'Replica derivative d Delta_K/d alpha |_{alpha=1} = +1/6', 'F', 1, 'M2', 'dimensionless',
          'yes', 'none', 'gravity', 'Paper 51 G4-4f closed form'),
    Entry('D5', 'Cutoff function f moments phi(0), phi(1), phi(2)', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'gravity',
          'no autonomous selection; CC fine-tuning phi(2)/phi(1)^2 ~ 10^-124 required'),
    Entry('D6', 'Cosmological constant value Lambda_cc', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'gravity', 'inherits D5'),
    Entry('D7', 'Relations among G_N, Lambda_cc, S_BH (Wald)', 'F', 1, 'M1', 'dimensionless',
          'yes', 'none', 'gravity',
          'two-term-exactness => pure Einstein => action-G = entropy-G'),

    # E. QED and alpha
    Entry('E1', 'Self-energy structural zero Sigma(n_ext=0) = 0', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'qed', 'Paper 28 theorem'),
    Entry('E2', 'Selection rules 8/8 (1+6+1 partition)', 'F', 1, 'M1', 'dimensionless',
          'yes', 'none', 'qed', 'Paper 33 vector-photon promotion'),
    Entry('E3', 'F_2 = 5 sqrt(2)/3 graph-native (pi-free)', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'qed', 'Paper 33 GN arc'),
    Entry('E4', 'F_2/(alpha/2 pi) = 1.084 Parker-Toms (R/12)', 'F', 1, 'M2', 'dimensionless',
          'yes', 'none', 'qed', 'Paper 33 bit-exact 0.5% match'),
    Entry('E5', 'B = 42, F = pi^2/6, Delta = 1/40 spectral homes', 'F', 1, 'M1/M2/M2', 'dimensionless',
          'yes', 'none', 'alpha', 'three independent spectral homes proven'),
    Entry('E6', 'Combination rule K = pi(B + F - Delta) for alpha^-1', 'C', 2, 'M1+M2', 'dimensionless',
          'no', 'multi_focal_composition', 'alpha',
          'numerical coincidence 8.8e-8; 12 mechanisms eliminated; Observation per 13.5'),
    Entry('E7', 'Two-loop SE counterterms Z_2, delta m', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'qed', 'LS-8a wall'),
    Entry('E8', 'Multi-loop QED corrections beyond LS-7', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'qed', 'HF-5 wall'),

    # F. SM inner factor (all C, all packing-unreachable)
    Entry('F1', 'Yukawa values (9 charged-lepton + quark)', 'C', 2, 'outside', 'dimensionless',
          'no', 'multi_focal_composition', 'inner',
          'Sprint H1 non-selection + 162-cell Yukawa-PSLQ clean negative'),
    Entry('F2', 'Generation count N_gen = 3', 'C', 1, 'none', 'dimensionless',
          'no', 'inner_factor_input_data', 'inner',
          'Direction 2 + Read 2 NO-GO (packing-unreachable)'),
    Entry('F3', 'Inner KO-dim signature', 'C', 1, 'none', 'dimensionless',
          'no', 'inner_factor_input_data', 'inner', 'Direction 2 NO-GO (packing-unreachable)'),
    Entry('F4', 'Higgs VEV v', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'inner', 'mass scale calibration input'),
    Entry('F5', 'CKM mixing angles', 'C', 2, 'outside', 'dimensionless',
          'no', 'multi_focal_composition', 'inner',
          'Boyle-Farnsworth gives shape not values'),
    Entry('F6', 'PMNS mixing angles', 'C', 2, 'outside', 'dimensionless',
          'no', 'multi_focal_composition', 'inner', 'same as F5'),
    Entry('F7', 'Neutrino mass values', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'inner', 'inherits F1 + F4'),

    # G. Multi-focal-composition walls
    Entry('G1', 'HF-3 recoil cross-register V_eN(r_e, R_n)', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'qed', 'multi-focal wall pattern instance 1'),
    Entry('G2', 'HF-4 Zemach magnetization density', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'qed', 'multi-focal wall instance 2 (W1b extended)'),
    Entry('G3', 'HF-5 multi-loop QED on hyperfine', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'qed', 'multi-focal wall instance 3'),
    Entry('G4', 'LS-8a two-loop SE renormalization (= E7)', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'qed', 'multi-focal wall instance 4'),
    Entry('G5', 'W1e chemistry inner-region overattraction', 'C', 2, 'outside', 'dimensionful',
          'no', 'multi_focal_composition', 'chemistry',
          'multi-focal wall instance 5; chemistry-side analog of H1'),
    Entry('G6', 'H1 Yukawa non-selection (= F1)', 'C', 2, 'outside', 'dimensionless',
          'no', 'multi_focal_composition', 'inner', 'multi-focal wall instance 6'),

    # H. Chemistry / QC scaling
    Entry('H1', 'Angular sparsity 1.44% depends only on l_max', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'chemistry', 'Paper 22 potential-independent theorem'),
    Entry('H2', 'N_Pauli = 11.10 x Q across 38 molecules', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'qc', 'Paper 14 isostructural invariance'),
    Entry('H3', 'Pauli ratio balanced/composed ~ O(B^2)', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'qc', 'proven structural'),
    Entry('H4', 'Nuclear magic numbers 2, 8, 20, 40, 70, 112', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'nuclear', 'HO shell closures from graph counting; Paper 23'),
    Entry('H5', 'Fock projection rigidity (S^3 unique to -Z/r)', 'F', 1, 'none', 'dimensionless',
          'yes', 'none', 'nuclear', 'Paper 23 NJ theorem'),
    Entry('H6', 'Clementi-Raimondi Z_eff exponents', 'C', 1, 'outside', 'dimensionless',
          'no', 'inner_factor_input_data', 'chemistry',
          'external atomic-physics fits; W1e period-class memo'),
    Entry('H7', 'Multi-zeta coefficients in second-row valence basis', 'C', 1, 'outside', 'dimensionless',
          'no', 'inner_factor_input_data', 'chemistry', 'inherits H6'),

    # I. Foundational calibration
    Entry('I1', 'Fine-structure constant alpha (value)', 'C', 2, 'M1+M2', 'dimensionless',
          'no', 'multi_focal_composition', 'alpha', 'inherits E6 (combination rule is an Observation)'),
    Entry('I2', 'Born rule p = |<a|psi>|^2', 'C', 1, 'N/A', 'N/A',
          'no', 'inner_factor_input_data', 'foundational',
          'Gleason via Hilbert inheritance; framework does not improve on Gleason'),
    Entry('I3', 'Higgs direction n_hat in S^2', 'C', 1, 'M1', 'dimensionless',
          'conditional', 'inner_factor_input_data', 'inner',
          'Boyle-Farnsworth input; possible Hopf-base identification (flagged)'),
]


# =====================================================================
# PREDICATES
# =====================================================================

def true_status(e: Entry) -> str:
    """Reduce F/A/C to F/C for the binary test.
    Treat 'admitted' as 'forced' (the framework structurally allows it).
    """
    return 'F' if e.status in ('F', 'A') else 'C'


def P1(e: Entry) -> str:
    """Multi-focal depth: forced iff MF = 1."""
    return 'F' if e.mf == 1 else 'C'


def P2(e: Entry) -> str:
    """Period-content classifiability: forced iff P != 'outside'."""
    return 'C' if e.p == 'outside' else 'F'


def P3(e: Entry) -> str:
    """Dimensional character: forced iff dimensionless."""
    return 'F' if e.dim == 'dimensionless' else 'C'


def P4(e: Entry) -> str:
    """Compactness inheritance: forced iff packing-reachable AND
    (not multi-focal OR period not 'outside'). This refines Paper 57 §5.4
    so the criterion is testable cell-by-cell.
    """
    if e.packing_reachable == 'no':
        return 'C'
    if e.packing_reachable == 'conditional':
        return 'F'  # admit conditional reachability
    if e.mf == 2 and e.p == 'outside':
        return 'C'
    return 'F'


def P5(e: Entry) -> str:
    """Packing-reachability: forced iff a witness derivation exists from
    {Paper 0 axiom + standard NCG axioms + Upgrade B}.
    """
    return 'F' if e.packing_reachable in ('yes', 'conditional') else 'C'


def P1_and_P5(e: Entry) -> str:
    """Conjunction: predict F iff both P1 and P5 predict F."""
    return 'F' if P1(e) == 'F' and P5(e) == 'F' else 'C'


def P1_or_P5(e: Entry) -> str:
    """Disjunction: predict F iff either P1 or P5 predicts F."""
    return 'F' if P1(e) == 'F' or P5(e) == 'F' else 'C'


PREDICATES = {
    'P1_multifocal_depth': P1,
    'P2_period_classifiable': P2,
    'P3_dimensional': P3,
    'P4_compactness': P4,
    'P5_packing_reachable': P5,
    'P1_and_P5': P1_and_P5,
    'P1_or_P5': P1_or_P5,
}


# =====================================================================
# EVALUATION
# =====================================================================

def evaluate(name: str, predicate: Callable[[Entry], str], entries):
    tp = fp = tn = fn = 0
    misclass = []
    per_family_fail = {
        'none': 0,
        'multi_focal_composition': 0,
        'inner_factor_input_data': 0,
        'admitted_only': 0,
    }
    for e in entries:
        pred = predicate(e)
        truth = true_status(e)
        if pred == 'F' and truth == 'F':
            tp += 1
        elif pred == 'F' and truth == 'C':
            fp += 1
            misclass.append({
                'id': e.id, 'name': e.name,
                'predicted': 'F', 'actual': 'C',
                'family': e.family, 'reason': e.note,
            })
            per_family_fail[e.family] += 1
        elif pred == 'C' and truth == 'C':
            tn += 1
        else:
            fn += 1
            misclass.append({
                'id': e.id, 'name': e.name,
                'predicted': 'C', 'actual': 'F',
                'family': e.family, 'reason': e.note,
            })
            per_family_fail[e.family] += 1

    n = len(entries)
    return {
        'predicate': name,
        'n': n,
        'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn,
        'accuracy': (tp + tn) / n,
        'misclassification_count': fp + fn,
        'misclassified_entries': misclass,
        'per_family_failures': per_family_fail,
    }


def summarize_catalogue(entries):
    status_counts = {}
    family_counts = {}
    packing_counts = {}
    for e in entries:
        status_counts[e.status] = status_counts.get(e.status, 0) + 1
        family_counts[e.family] = family_counts.get(e.family, 0) + 1
        packing_counts[e.packing_reachable] = packing_counts.get(e.packing_reachable, 0) + 1
    return {
        'total': len(entries),
        'by_status': status_counts,
        'by_family': family_counts,
        'by_packing_reachable': packing_counts,
    }


def main():
    summary = summarize_catalogue(ENTRIES)
    print("=" * 70)
    print("FORCED/FREE SEAM CATALOGUE — Sprint C2 principle hunt audit")
    print("=" * 70)
    print(f"\nCatalogue summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    results = {}
    for name, pred in PREDICATES.items():
        results[name] = evaluate(name, pred, ENTRIES)
        r = results[name]
        print(f"\n--- {name} ---")
        print(f"  accuracy: {r['accuracy']:.3f}  ({r['tp']+r['tn']}/{r['n']})")
        print(f"  TP={r['tp']}  FP={r['fp']}  TN={r['tn']}  FN={r['fn']}")
        print(f"  misclassifications: {r['misclassification_count']}")
        if r['misclassified_entries']:
            print(f"  failing entries by family:")
            for fam, cnt in r['per_family_failures'].items():
                if cnt > 0:
                    print(f"    {fam}: {cnt}")
            print(f"  failing entry IDs: {[m['id'] for m in r['misclassified_entries']]}")

    # Save results.
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "principle_hunt_audit.json"
    out_path.write_text(json.dumps({
        'catalogue_summary': summary,
        'predicate_results': results,
        'entries': [asdict(e) for e in ENTRIES],
    }, indent=2))
    print(f"\nResults written to {out_path}")


if __name__ == "__main__":
    main()
