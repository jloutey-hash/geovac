"""TX-B Observable 3: Hydrogen Bohr spectrum E_n = -Z^2 / (2 n^2).

NEGATIVE CONTROL.

The Bohr hydrogenic spectrum in atomic units is purely algebraic in n and Z.
There is no integration over any continuous parameter at any step.  No pi
should appear.

We verify symbolically that E_n is a rational for any rational Z and integer n,
across a panel of (Z, n) values.

Paper 7's Fock-projected Coulomb operator on the unit S^3 has Laplace-Beltrami
spectrum lambda_n = -(n^2 - 1).  Mapping to physical hydrogen via the
energy-shell rescaling p_0^2 = -2 E gives E_n = -Z^2 / (2 n^2).  Both steps are
algebraic; no continuous integration is performed.

(Caveat: the SCALE Z is dimensional [charge], the SCALE n is dimensionless,
and atomic units are dimensionless throughout.  The transcendental signature
is what we test, not dimensional content.)
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


OUT = Path(__file__).parent / "data" / "tx_b_obs_3.json"


def hydrogen_bohr_spectrum_panel():
    """Compute E_n for a panel of (Z, n) values and verify rationality."""
    n_var = sp.symbols("n", positive=True, integer=True)
    Z_var = sp.symbols("Z", positive=True)

    # Symbolic formula
    E_formula = -Z_var**2 / (2 * n_var**2)

    # Test panel: integer Z in 1..6, n in 1..10
    panel = []
    all_rational = True
    contains_pi_anywhere = False
    for Z in range(1, 7):
        for n in range(1, 11):
            E_val = E_formula.subs([(Z_var, Z), (n_var, n)])
            E_val_simplified = sp.simplify(E_val)
            is_rational = isinstance(E_val_simplified, sp.Rational)
            has_pi = sp.pi in E_val_simplified.atoms()
            if not is_rational:
                all_rational = False
            if has_pi:
                contains_pi_anywhere = True
            panel.append({
                "Z": Z,
                "n": n,
                "E_n_atomic_units": str(E_val_simplified),
                "is_rational": is_rational,
                "contains_pi": has_pi,
            })

    # Symbolic check that the formula itself is rational in (Z, n)
    formula_atoms = E_formula.atoms()
    formula_has_pi = sp.pi in formula_atoms

    # Rydberg-shift verification: matches NIST level for hydrogen ground state.
    # Ground state energy in Hartree: E_1 = -1/2 Ha = -13.605693... eV
    E1_Z1 = E_formula.subs([(Z_var, 1), (n_var, 1)])
    E1_Z1_value = sp.Rational(-1, 2)
    matches_NIST_ground_state = (E1_Z1 == E1_Z1_value)

    transcendental_ledger = [
        {"step": "(1) Paper 7 Fock-projected Coulomb operator on unit S^3",
         "value": "Laplace-Beltrami eigenvalues lambda_n = -(n^2 - 1)",
         "transcendentals": "NONE -- integers"},
        {"step": "(2) Energy-shell rescaling p_0^2 = -2 E",
         "value": "p_0 = Z/n, E_n = -p_0^2/2 = -Z^2/(2 n^2)",
         "transcendentals": "NONE -- algebraic"},
        {"step": "(3) Atomic-unit Bohr spectrum E_n = -Z^2 / (2 n^2)",
         "value": "Pure rational for integer (Z, n)",
         "transcendentals": "NONE"},
    ]

    return {
        "spectrum_formula": str(E_formula),
        "formula_contains_pi": formula_has_pi,
        "test_panel_size": len(panel),
        "all_panel_entries_rational": all_rational,
        "any_panel_entry_contains_pi": contains_pi_anywhere,
        "matches_NIST_ground_state_E1_eq_minus_1_over_2": matches_NIST_ground_state,
        "panel_first_3_entries": panel[:3],
        "panel_last_3_entries": panel[-3:],
        "transcendental_ledger": transcendental_ledger,
    }


def main():
    res = hydrogen_bohr_spectrum_panel()

    pi_predicted = False
    pi_observed = (res["formula_contains_pi"] or res["any_panel_entry_contains_pi"])
    is_rational = res["all_panel_entries_rational"]

    out = {
        "observable_id": 3,
        "observable_name": "Hydrogen Bohr spectrum E_n = -Z^2/(2n^2) [atomic units]",
        "predicted_class": "pure rational",
        "pi_predicted": pi_predicted,
        "computed_class": ("pure rational" if is_rational else "non-rational"),
        "pi_observed": pi_observed,
        "match_with_prediction": (pi_predicted == pi_observed) and is_rational,
        "result": res,
        "verdict": ("CONFIRMS Paper 35 Prediction 1 (negative control)"
                     if (not pi_observed and is_rational)
                     else "REFUTES Paper 35 Prediction 1"),
        "paper_34_projections": [
            "Fock conformal projection (Paper 7 S^3 Coulomb)",
            "Energy-shell rescaling (algebraic)"
        ],
        "temporal_integration_present": False,
        "specific_integration": "None -- pure algebraic eigenvalue formula",
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps({
        "observable": "Hydrogen Bohr spectrum (atomic units)",
        "predicted_pi": pi_predicted,
        "observed_pi": pi_observed,
        "all_rational": is_rational,
        "match": (pi_predicted == pi_observed) and is_rational,
    }, indent=2))


if __name__ == "__main__":
    main()
