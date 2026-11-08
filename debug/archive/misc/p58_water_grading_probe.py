"""Paper 58 Prediction 1 probe: does the molecular point group's abelian
grading predict the surviving symmetry-based qubit reduction?

Prediction as drafted (paper_58 pred:water): H2O has point group C2v, order 4,
abelian = Z2 x Z2 => a clean TWO-bit spatial grading => 2 qubits of SPATIAL
reduction, where linear molecules carry a continuous U(1) => Z grading =>
a multiplicative factor instead.

--- ROUND 1 (superseded, kept as a record of the artifact) ---
First run called extended_tapered_from_spec WITHOUT `nuclei` and measured
spatial dQ = 0 for all three of LiH/BeH2/H2O.  That was a PROBE BUG, not a
finding: inside the pipeline `nuclei_list = nuclei or []` and the branch is
guarded by `if use_atom_swap and nuclei_list:`, so the entire spatial branch is
silently skipped when nuclei are absent.  MolecularSpec.nuclei is None for all
hydride_spec-derived systems.

--- ROUND 2 (this file) ---
Two corrections:

(a) Supply nuclei.  Two sources are legitimate:
      * spec-supplied: _diatomic_multi_center_spec populates nuclei, so N2/F2
        (homonuclear, D_inf_h, centrosymmetric, equivalent atoms) are the
        well-posed linear-symmetric test cases.
      * angle-free hand-built: BeH2 is linear, so H at (0,0,+/-R) with Be at
        the origin is unambiguous and involves no invented geometry.

(b) DO NOT hand-build bent H2O geometry.  molecular_spec.py contains no bond
    angle anywhere (grep: angle|104|bent|theta -> zero hits); h2o_spec is
    hydride_spec(8), parameterized by Z and a single O-H distance across 5
    blocks (O_core, 2 bond pairs, 2 lone pairs).  The composed H2O Hamiltonian
    therefore carries NO angular geometry, so C2v is not represented in it.
    Feeding invented bent nuclei to the stabilizer builder would test a
    symmetry the Hamiltonian does not encode -- a confound, not a measurement.
    Prediction 1 is NOT WELL-POSED on this builder; that is the finding.

Run from repo root:  python debug/p58_water_grading_probe.py
"""

from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.extended_tapering import (  # noqa: E402
    extended_tapered_from_spec,
    find_equivalent_atom_pairs,
    is_centrosymmetric,
)
from geovac import molecular_spec  # noqa: E402


def linear_symmetric_nuclei(Z_center: float, Z_out: float, R: float,
                            sym_c: str, sym_o: str):
    """X-A-X on the z-axis.  Angle-free, so no geometry is invented."""
    return [
        {'Z': float(Z_center), 'position': (0.0, 0.0, 0.0), 'label': sym_c},
        {'Z': float(Z_out), 'position': (0.0, 0.0, float(R)), 'label': sym_o + '1'},
        {'Z': float(Z_out), 'position': (0.0, 0.0, -float(R)), 'label': sym_o + '2'},
    ]


def build_cases():
    """(label, spec, nuclei, point group, nuclei provenance)."""
    cases = []

    # Spec-supplied nuclei: homonuclear diatomics are D_inf_h.
    for name, pg in (("N2", "D_inf_h"), ("F2", "D_inf_h")):
        spec = getattr(molecular_spec, name.lower() + "_spec")()
        cases.append((name, spec, spec.nuclei, pg, "spec-supplied"))

    # Heteronuclear diatomic control: no equivalent atoms, not centrosymmetric.
    lih = molecular_spec.lih_spec()
    cases.append(("LiH", lih, [
        {'Z': 3.0, 'position': (0.0, 0.0, 0.0), 'label': 'Li'},
        {'Z': 1.0, 'position': (0.0, 0.0, float(lih.R)), 'label': 'H'},
    ], "C_inf_v", "hand-built (linear, angle-free)"))

    # Linear symmetric hydride: angle-free, unambiguous.
    beh2 = molecular_spec.beh2_spec()
    cases.append(("BeH2", beh2,
                  linear_symmetric_nuclei(4.0, 1.0, beh2.R, 'Be', 'H'),
                  "D_inf_h", "hand-built (linear, angle-free)"))

    return cases


CONFIGS = [
    ("hopf+ell",        True, True, False, False),
    ("+swap",           True, True, True,  False),
    ("+swap+inversion",  True, True, True,  True),
]


def n_pauli(qubit_op) -> int:
    if qubit_op is None:
        return -1
    return sum(1 for term in qubit_op.terms if term)


def spatial_kinds(kinds):
    return [k for k in kinds if ('swap' in k.lower() or 'inv' in k.lower())]


def main() -> None:
    results = {"_h2o_not_well_posed": {
        "reason": "composed builder carries no bond angle; C2v unrepresented",
        "evidence": "molecular_spec.py has zero hits for angle|104|bent|theta; "
                    "h2o_spec == hydride_spec(8), single O-H distance, 5 blocks",
    }}

    for label, spec, nuclei, pg, provenance in build_cases():
        print(f"\n{'=' * 70}\n{label}  ({pg})   nuclei: {provenance}\n{'=' * 70}")
        entry = {"point_group": pg, "nuclei_provenance": provenance,
                 "configs": {}}
        try:
            pairs = find_equivalent_atom_pairs(spec, nuclei)
            entry["n_equivalent_atom_pairs"] = len(pairs) if pairs else 0
            entry["centrosymmetric"] = bool(is_centrosymmetric(nuclei))
            print(f"  equivalent atom pairs : {entry['n_equivalent_atom_pairs']}")
            print(f"  centrosymmetric       : {entry['centrosymmetric']}")
        except Exception as exc:
            entry["geometry_probe_error"] = f"{type(exc).__name__}: {exc}"
            print(f"  [geometry probe: {type(exc).__name__}: {exc}]")

        for cfg, hopf, ell, swap, inv in CONFIGS:
            try:
                out = extended_tapered_from_spec(
                    spec, use_hopf=hopf, use_ell_parity=ell,
                    use_atom_swap=swap, use_inversion=inv,
                    nuclei=nuclei,
                )
                kinds = list(out.get("kinds_kept") or [])
                row = {
                    "Q_naive": out.get("Q_naive"),
                    "Q_tapered": out.get("Q_tapered"),
                    "delta_Q": out.get("delta_Q"),
                    "n_stabs_kept": out.get("n_stabs_kept"),
                    "n_kinds": len(kinds),
                    "spatial_kinds": spatial_kinds(kinds),
                    "n_pauli_tapered": n_pauli(out.get("qubit_op_tapered")),
                    "n_pauli_naive": n_pauli(out.get("qubit_op_naive")),
                }
                entry["configs"][cfg] = row
                print(f"  {cfg:<17} Q {row['Q_naive']}->{row['Q_tapered']} "
                      f"(dQ={row['delta_Q']})  "
                      f"Pauli {row['n_pauli_naive']}->{row['n_pauli_tapered']}  "
                      f"spatial_kinds={row['spatial_kinds']}")
            except Exception as exc:
                entry["configs"][cfg] = {"error": f"{type(exc).__name__}: {exc}",
                                         "traceback": traceback.format_exc(limit=3)}
                print(f"  {cfg:<17} FAILED: {type(exc).__name__}: {exc}")

        results[label] = entry

    print(f"\n{'=' * 70}\nSPATIAL CONTRIBUTION\n{'=' * 70}")
    print(f"{'system':<7}{'group':<10}{'base dQ':>9}{'+swap':>7}{'+sw+inv':>9}"
          f"{'spatial':>9}{'Pauli base':>11}{'Pauli sp.':>10}")
    for label in [k for k in results if not k.startswith('_')]:
        c = results[label]["configs"]
        base = c.get("hopf+ell", {}).get("delta_Q")
        sw = c.get("+swap", {}).get("delta_Q")
        si = c.get("+swap+inversion", {}).get("delta_Q")
        pb = c.get("hopf+ell", {}).get("n_pauli_tapered")
        ps = c.get("+swap+inversion", {}).get("n_pauli_tapered")
        spatial = (si - base) if isinstance(si, int) and isinstance(base, int) else None
        results[label]["spatial_delta_Q"] = spatial
        pg = results[label]["point_group"]
        print(f"{label:<7}{pg:<10}{str(base):>9}{str(sw):>7}{str(si):>9}"
              f"{str(spatial):>9}{str(pb):>11}{str(ps):>10}")

    print("\nH2O: Prediction 1 is NOT WELL-POSED on the composed builder "
          "(no bond angle => C2v unrepresented).")

    out_path = REPO / "debug" / "data" / "p58_water_grading_probe.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
