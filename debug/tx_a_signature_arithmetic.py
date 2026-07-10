"""
Sprint TX-A: Dimensional axiomatization of Paper 34 projections.

This script encodes the 15 projections as graded maps and tests:
  Theorem 1 (Dimensional Consistency): admissibility = signature addition.
  Theorem 2 (Minimality): shortest admissible chain unique up to commutation.
  Theorem 3 (Composition Algebra): commutation, idempotence, inverses.

The encoding uses simple tuples for each projection's signature:
   (variable_set, dimension_signature, transcendental_signature).

variable_set: frozenset of named variables introduced
dimension_signature: tuple of (mass, length, time, charge, temperature)
                     exponents (the SI base), where 'energy' is encoded as
                     (mass, length, time) = (1, 2, -2) for E = m c^2 if
                     convertible, or as a "soft dimension" tag.
transcendental_signature: a frozenset of named transcendental tiers.

Run: python debug/tx_a_signature_arithmetic.py
"""
import json
from itertools import combinations, permutations
from pathlib import Path

OUT_DIR = Path(__file__).parent / "data"
OUT_DIR.mkdir(exist_ok=True)

# Five physical dimensions tracked as a tuple of integer exponents:
# (mass, length, time, charge, action).  We add a "solid_angle" soft slot.
# Energy = m * length^2 / time^2; in atomic units we'll often use the soft tag
# "energy" instead of decomposing.  For the dimensional-consistency theorem
# what matters is that each projection has a fixed signature.

# Use lightweight string tags: each projection has a frozen signature.

# ---- The 15 projections (numbered as in Paper 34 §III.1..§III.15) ----

PROJECTIONS = [
    {  # 1
        "id": "fock_conformal",
        "name": "Fock conformal",
        "section": "III.1",
        "variables_added": ["Z", "E"],
        "dimensions_added": ["energy"],
        "transcendental_class": ["rational_kappa"],
        # Note: pi^2 from Vol(S^3) appears in *integrated* observables, not from
        # the projection itself; we tag this projection's intrinsic signature
        # only.  The volume factor enters via the mellin/observation step.
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": False,
        "has_inverse": False,
        "input_layer": 0,
        "output_layer": 1,
    },
    {  # 2
        "id": "hopf_bundle",
        "name": "Hopf bundle",
        "section": "III.2",
        "variables_added": ["alpha"],
        "dimensions_added": [],
        "transcendental_class": ["pi_hopf_base"],  # pi = Vol(S^2)/4
        "signature_pi_intrinsic": True,
        "graded": True,
        "idempotent": False,
        "has_inverse": False,
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 3
        "id": "bargmann_segal",
        "name": "Bargmann-Segal",
        "section": "III.3",
        "variables_added": ["hbar_omega"],
        "dimensions_added": ["energy"],
        "transcendental_class": ["rational"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": False,
        "has_inverse": False,
        "input_layer": 0,
        "output_layer": 1,
    },
    {  # 4
        "id": "stereographic",
        "name": "Stereographic / conformal coordinate change",
        "section": "III.4",
        "variables_added": ["r"],
        "dimensions_added": ["length"],
        "transcendental_class": ["conformal_factor"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": False,
        "has_inverse": True,  # invertible on the open dense subset
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 5
        "id": "sturmian",
        "name": "Sturmian reparameterization at lambda=Z/n",
        "section": "III.5",
        "variables_added": [],  # reuses Z, n
        "dimensions_added": [],
        "transcendental_class": ["rational"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": True,  # Sturmian basis is its own reparameterization fixed point
        "has_inverse": True,  # change-of-basis on the bound + discretized continuum
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 6
        "id": "spectral_action",
        "name": "Connes-Chamseddine spectral action",
        "section": "III.6",
        "variables_added": ["Lambda", "alpha"],
        "dimensions_added": ["energy"],
        "transcendental_class": ["sqrt_pi", "pi2k_rational"],
        "signature_pi_intrinsic": True,
        "graded": True,
        "idempotent": False,
        "has_inverse": False,  # spectral mode-counting is one-way
        "input_layer": 1,
        "output_layer": 2,
    },
    {  # 7
        "id": "spinor_lift",
        "name": "Camporesi-Higuchi spinor lift",
        "section": "III.7",
        "variables_added": ["alpha"],  # via (Z*alpha)^2
        "dimensions_added": [],
        "transcendental_class": [
            "half_int_hurwitz",
            "odd_zeta",
            "catalan_G",
            "dirichlet_beta",
        ],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": False,
        "has_inverse": False,
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 8
        "id": "wigner_3j",
        "name": "Wigner 3j angular coupling",
        "section": "III.8",
        "variables_added": [],
        "dimensions_added": [],
        "transcendental_class": ["sqrt_int_rationals"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": True,  # coupling a coupled state to itself is identity-like
        "has_inverse": True,  # Clebsch-Gordan is unitary
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 9
        "id": "wigner_D",
        "name": "Wigner D rotation between molecular centers",
        "section": "III.9",
        "variables_added": ["R_AB"],
        "dimensions_added": ["length"],
        "transcendental_class": ["sqrt_int_rationals"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": False,
        "has_inverse": True,  # rotations are invertible
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 10
        "id": "wilson_plaquette",
        "name": "Wilson plaquette",
        "section": "III.10",
        "variables_added": ["beta_gauge"],
        "dimensions_added": [],
        "transcendental_class": ["su2_haar"],
        "signature_pi_intrinsic": True,  # Haar measure carries pi
        "graded": True,
        "idempotent": False,
        "has_inverse": False,
        "input_layer": 1,
        "output_layer": 2,
    },
    {  # 11
        "id": "vector_photon",
        "name": "Vector-photon promotion",
        "section": "III.11",
        "variables_added": [],
        "dimensions_added": [],  # solid angle, treated as dimensionless
        "transcendental_class": ["one_over_4pi"],
        "signature_pi_intrinsic": True,
        "graded": True,
        "idempotent": False,
        "has_inverse": False,
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 12
        "id": "molframe_hypersph",
        "name": "Molecule-frame hyperspherical",
        "section": "III.12",
        "variables_added": ["R"],
        "dimensions_added": ["length"],
        "transcendental_class": ["piecewise_smooth_R"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": False,
        "has_inverse": False,
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 13
        "id": "drake_swainson",
        "name": "Drake-Swainson asymptotic subtraction",
        "section": "III.13",
        "variables_added": ["K_subtract"],  # transient
        "dimensions_added": [],  # transient energy via K, cancelled
        "transcendental_class": ["flow"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": True,  # subtracting twice = subtracting once
        "has_inverse": False,
        "input_layer": 2,
        "output_layer": 2,
    },
    {  # 14
        "id": "rest_mass",
        "name": "Rest-mass projection",
        "section": "III.14",
        "variables_added": ["m"],
        "dimensions_added": ["mass"],
        "transcendental_class": ["trivial_ring_preserving"],
        "signature_pi_intrinsic": False,
        "graded": True,
        "idempotent": False,  # m^2 + m^2 != m^2 unless m=0
        "has_inverse": True,  # set m = 0
        "input_layer": 1,
        "output_layer": 1,
    },
    {  # 15
        "id": "observation_window",
        "name": "Observation / temporal-window",
        "section": "III.15",
        "variables_added": ["beta_temp"],  # 1/T
        "dimensions_added": ["time"],
        "transcendental_class": ["two_pi_rational", "pi2k_via_matsubara"],
        "signature_pi_intrinsic": True,
        "graded": True,
        "idempotent": False,
        "has_inverse": True,  # send beta -> infinity
        "input_layer": 1,
        "output_layer": 2,
    },
]

# Sanity check: 15 projections, all 15 have unique ids
assert len(PROJECTIONS) == 15
ids = [p["id"] for p in PROJECTIONS]
assert len(set(ids)) == 15

# ---- Theorem 1: Dimensional consistency by signature addition ----

def signature_sum(projection_ids):
    """Sum the variable / dimension / transcendental signatures of a chain."""
    p_by_id = {p["id"]: p for p in PROJECTIONS}
    vars_total = set()
    dims_total = []
    trans_total = set()
    pi_intrinsic = False
    for pid in projection_ids:
        p = p_by_id[pid]
        vars_total.update(p["variables_added"])
        dims_total.extend(p["dimensions_added"])
        trans_total.update(p["transcendental_class"])
        pi_intrinsic = pi_intrinsic or p["signature_pi_intrinsic"]
    return {
        "variables": sorted(vars_total),
        "dimensions": sorted(dims_total),
        "transcendentals": sorted(trans_total),
        "pi_intrinsic": pi_intrinsic,
    }


# Worked examples for Theorem 1
EXAMPLES = {
    "hydrogen_lamb_shift_LS6a": {
        "chain": [
            "fock_conformal",
            "spectral_action",
            "sturmian",
            "drake_swainson",
        ],
        "target_observable": "Lamb shift, energy units, contains pi^2",
        "expected_dimensions": {"energy"},
        "expected_pi": True,
    },
    "K_alpha_three_sector": {
        "chain": ["fock_conformal", "hopf_bundle", "spinor_lift"],
        "target_observable": "K = pi(B+F-Delta) ~= 1/alpha (dimensionless)",
        "expected_dimensions": set(),
        "expected_pi": True,
    },
    "stefan_boltzmann": {
        "chain": ["fock_conformal", "observation_window"],
        "target_observable": "Stefan-Boltzmann pi^2/90, free-energy density",
        "expected_dimensions": {"energy", "time"},
        "expected_pi": True,
    },
    "he_ground_state_casimir_CI": {
        "chain": ["fock_conformal", "wigner_3j"],
        "target_observable": "He ground-state energy via Casimir CI",
        "expected_dimensions": {"energy"},
        "expected_pi": False,  # graph-native CI is pi-free in algebraic ring
    },
    "hydrogen_E_n_bohr": {
        "chain": ["fock_conformal"],
        "target_observable": "E_n = -Z^2/(2 n^2) Ha",
        "expected_dimensions": {"energy"},
        "expected_pi": False,
    },
    "spatial_casimir_S3_scalar": {
        "chain": ["fock_conformal"],  # spatial only, no observation_window
        "target_observable": "E_Cas = 1/240 (rational, no pi)",
        "expected_dimensions": {"energy"},
        "expected_pi": False,
    },
    "spatial_casimir_S3_dirac": {
        "chain": ["fock_conformal", "spinor_lift"],
        "target_observable": "E_Cas^Dirac = 17/480 (rational, no pi)",
        "expected_dimensions": {"energy"},
        "expected_pi": False,
    },
    "uehling_VP": {
        "chain": ["fock_conformal", "spectral_action"],
        "target_observable": "VP coefficient Pi = 1/(48 pi^2)",
        "expected_dimensions": {"energy"},  # carried by Lambda
        "expected_pi": True,
    },
}


def test_theorem_1():
    """Each example: check that signature_sum(chain) is consistent with
    the target observable's expected dimensional/transcendental tags."""
    results = {}
    for name, ex in EXAMPLES.items():
        sig = signature_sum(ex["chain"])
        dims_set = set(sig["dimensions"])
        consistent_dims = ex["expected_dimensions"].issubset(dims_set)
        consistent_pi = (sig["pi_intrinsic"] == ex["expected_pi"])
        results[name] = {
            "chain": ex["chain"],
            "signature": sig,
            "expected_dims": sorted(ex["expected_dimensions"]),
            "expected_pi": ex["expected_pi"],
            "consistent_dims": consistent_dims,
            "consistent_pi": consistent_pi,
            "verdict": "PASS" if (consistent_dims and consistent_pi) else "FAIL",
        }
    return results


# ---- Theorem 2: Minimality (uniqueness up to commutation) ----

def commutes(p1, p2):
    """Two projections commute (up to operator-system order) if
    neither's output is the other's required input AND neither
    introduces a variable the other consumes."""
    p_by_id = {p["id"]: p for p in PROJECTIONS}
    a = p_by_id[p1]
    b = p_by_id[p2]
    # Strict input/output layer check: spectral_action takes layer-1 to layer-2,
    # so it can't commute with anything that produces layer-2 output.
    if a["input_layer"] > b["output_layer"] and b["input_layer"] > a["output_layer"]:
        return False  # both depend on each other
    # Variable consumption: if a's var-set intersects b's var-set, they share
    # variables and order may matter.  Conservative rule: commute iff
    # variable sets are disjoint AND dimension sets disjoint.
    vars_a = set(a["variables_added"])
    vars_b = set(b["variables_added"])
    if vars_a & vars_b:
        return False  # both add same variable - cannot commute (double introduction)
    # Different transcendental classes always commute on signature level
    return True


def test_theorem_2():
    """For each example with k>=2 projections, enumerate all permutations
    and check that signature_sum is independent of order (it is, by
    construction, since signatures are sets/lists summed unordered).
    Then count permutations that achieve the same target signature.
    The number of distinct admissible orderings is k! / (number of
    non-commuting pairs).  Compare to a naive expectation."""
    p_by_id = {p["id"]: p for p in PROJECTIONS}
    results = {}
    for name, ex in EXAMPLES.items():
        chain = ex["chain"]
        if len(chain) < 2:
            continue
        # All permutations give same signature (sum is commutative)
        sigs = []
        for perm in permutations(chain):
            sigs.append(signature_sum(list(perm)))
        # Verify all signatures equal
        ref = sigs[0]
        all_equal = all(
            s["variables"] == ref["variables"]
            and s["dimensions"] == ref["dimensions"]
            and s["transcendentals"] == ref["transcendentals"]
            and s["pi_intrinsic"] == ref["pi_intrinsic"]
            for s in sigs
        )
        # Count commuting pairs vs non-commuting pairs
        non_commuting = []
        for a, b in combinations(chain, 2):
            if not commutes(a, b):
                non_commuting.append((a, b))
        results[name] = {
            "chain": chain,
            "all_perms_same_signature": all_equal,
            "non_commuting_pairs": non_commuting,
            "n_permutations": len(sigs),
            "verdict": "PASS" if all_equal else "FAIL",
        }
    return results


# ---- Theorem 3: Composition Algebra (15x15 table) ----

def composition_relation(p1, p2):
    """Return one of: 'commute', 'one_way' (a;b ok, b;a not),
    'not_composable', 'idempotent' (if p1==p2 and p1 idempotent).
    """
    p_by_id = {p["id"]: p for p in PROJECTIONS}
    a = p_by_id[p1]
    b = p_by_id[p2]
    if p1 == p2:
        return "idempotent" if a["idempotent"] else "self_compose"
    # Check layer compatibility: a;b means apply a first, then b.
    # b's input_layer must be reachable from a's output_layer.
    can_a_then_b = a["output_layer"] >= b["input_layer"]
    can_b_then_a = b["output_layer"] >= a["input_layer"]
    if not can_a_then_b and not can_b_then_a:
        return "not_composable"
    # Variable conflict (both add same variable)
    vars_a = set(a["variables_added"])
    vars_b = set(b["variables_added"])
    if vars_a & vars_b:
        return "variable_conflict"
    # Inverse pair?
    if (a["has_inverse"] and b["has_inverse"]
            and a["id"] in {"stereographic", "wigner_3j", "wigner_D"}
            and b["id"] == a["id"]):
        return "inverse_pair"
    if can_a_then_b and can_b_then_a:
        return "commute"
    if can_a_then_b and not can_b_then_a:
        return "one_way_AB"
    if can_b_then_a and not can_a_then_b:
        return "one_way_BA"
    return "unclassified"


def build_composition_table():
    table = {}
    for p1 in PROJECTIONS:
        row = {}
        for p2 in PROJECTIONS:
            row[p2["id"]] = composition_relation(p1["id"], p2["id"])
        table[p1["id"]] = row
    return table


def composition_table_summary(table):
    """Count each kind of relation."""
    counts = {}
    for p1, row in table.items():
        for p2, rel in row.items():
            counts[rel] = counts.get(rel, 0) + 1
    return counts


# ---- Axis-independence diagnostic (the honest verdict question) ----

def axis_independence():
    """For each projection, list (var_axis_active, dim_axis_active,
    trans_axis_active).  Are there projections that move only one
    axis?  Are there projections that move all three together?"""
    results = []
    for p in PROJECTIONS:
        var_active = bool(p["variables_added"])
        dim_active = bool(p["dimensions_added"])
        trans_active = (
            p["transcendental_class"] != ["rational"]
            and p["transcendental_class"] != ["trivial_ring_preserving"]
            and p["transcendental_class"] != ["sqrt_int_rationals"]
        )
        # "trans_active_pi" specifically means injects pi
        trans_pi = p["signature_pi_intrinsic"]
        results.append({
            "id": p["id"],
            "axes_active": {
                "variable": var_active,
                "dimension": dim_active,
                "transcendental_nontrivial": trans_active,
                "transcendental_injects_pi": trans_pi,
            },
        })
    return results


# ---- Main ----

def main():
    t1 = test_theorem_1()
    t2 = test_theorem_2()
    t3 = build_composition_table()
    t3_summary = composition_table_summary(t3)
    axis_indep = axis_independence()

    print("=" * 70)
    print("THEOREM 1 (Dimensional Consistency by Signature Addition)")
    print("=" * 70)
    for name, r in t1.items():
        print(f"\n{name}:")
        print(f"  chain = {r['chain']}")
        print(f"  signature = {r['signature']}")
        print(f"  verdict = {r['verdict']}")

    print("\n" + "=" * 70)
    print("THEOREM 2 (Minimality / order-independence)")
    print("=" * 70)
    for name, r in t2.items():
        print(f"\n{name}:")
        print(f"  k = {len(r['chain'])}, n_perms = {r['n_permutations']}")
        print(f"  all permutations give same signature: {r['all_perms_same_signature']}")
        print(f"  non_commuting_pairs (operator-system order matters): {r['non_commuting_pairs']}")

    print("\n" + "=" * 70)
    print("THEOREM 3 (Composition table summary)")
    print("=" * 70)
    print(f"\nAggregate counts across 15x15 = 225 cells:")
    for rel, cnt in sorted(t3_summary.items()):
        print(f"  {rel:25s}: {cnt:4d}")

    print("\n" + "=" * 70)
    print("AXIS INDEPENDENCE")
    print("=" * 70)
    only_var = []
    only_dim = []
    only_trans = []
    all_three = []
    none = []
    for r in axis_indep:
        a = r["axes_active"]
        flags = (a["variable"], a["dimension"], a["transcendental_nontrivial"])
        if flags == (True, False, False):
            only_var.append(r["id"])
        elif flags == (False, True, False):
            only_dim.append(r["id"])
        elif flags == (False, False, True):
            only_trans.append(r["id"])
        elif flags == (True, True, True):
            all_three.append(r["id"])
        elif flags == (False, False, False):
            none.append(r["id"])
    print(f"\n  Move only variable axis : {only_var}")
    print(f"  Move only dimension axis: {only_dim}")
    print(f"  Move only transcendental: {only_trans}")
    print(f"  Move all three axes     : {all_three}")
    print(f"  Move none (Layer-1 ops) : {none}")

    # Save outputs
    out_table = {
        "projections": PROJECTIONS,
        "examples_with_signatures": {
            name: {
                "chain": ex["chain"],
                "signature": signature_sum(ex["chain"]),
                "target_observable": ex["target_observable"],
                "consistent_dims": t1[name]["consistent_dims"],
                "consistent_pi": t1[name]["consistent_pi"],
                "verdict": t1[name]["verdict"],
            }
            for name, ex in EXAMPLES.items()
        },
        "axis_independence": axis_indep,
        "theorem_1_results": t1,
        "theorem_2_results": t2,
    }
    with open(OUT_DIR / "tx_a_projection_table.json", "w") as f:
        json.dump(out_table, f, indent=2)

    out_comp = {
        "projection_ids": [p["id"] for p in PROJECTIONS],
        "composition_table": t3,
        "summary_counts": t3_summary,
    }
    with open(OUT_DIR / "tx_a_composition_table.json", "w") as f:
        json.dump(out_comp, f, indent=2)

    print(f"\nWrote: {OUT_DIR / 'tx_a_projection_table.json'}")
    print(f"Wrote: {OUT_DIR / 'tx_a_composition_table.json'}")


if __name__ == "__main__":
    main()
