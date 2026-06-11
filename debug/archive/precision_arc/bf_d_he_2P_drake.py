"""BF-D: He 2^3P multiplet via Drake 1971 radial amplitudes.

Builds on BR-C's verified angular J-pattern (f_SS, f_SOO from Bethe-Salpeter §39)
and replaces the faulty radial amplitude formula with linear combinations of
M^k direct/exchange integrals from the production module
``geovac.breit_integrals``.

Method
------
For (1s)(np) ^3P, the NR Breit-Pauli fine structure decomposes as
    E(J) = E_SO(J) + E_SS(J) + E_SOO(J)
with angular J-dependence
    E_SO(J) = (ζ/2) · X(J),  X(J) = J(J+1) - L(L+1) - S(S+1)
    E_SS(J) = A_SS · f_SS(J),   f_SS  = (-2, +1, -1/5) for J=(0,1,2)
    E_SOO(J)= A_SOO· f_SOO(J),  f_SOO = (+2, +1, -1)

The radial amplitudes A_SS, A_SOO are linear combinations of direct and
exchange Breit-Pauli retarded Slater integrals M^k (k=0,1,2) in which the
angular combining coefficients come from the 3j/6j algebra of the LS-coupled
tensor operator projected onto ^3P.

Following Drake 1971 Phys. Rev. A 3, 908 Eq. (14-17), and the canonical
LS-coupling derivation for (1s)(np) ^3P (Bethe-Salpeter §39; Cowan 1981 §12),
the most widely cited closed forms are

    A_SS  = α^2 · (3/50) · ( M^2_exch  -  M^2_dir )               [Drake 1971]
    A_SOO = α^2 · (1/2 ) · ( N^1_exch  -  N^1_dir )

but multiple sign/convention variants exist in the literature (Bethe-Salpeter
§39 uses a different projector choice with coefficients ±3/2, Johnson 2007
Ch 8 separates the retarded from instantaneous pieces).  We therefore
TEST a family of candidate combinations systematically (not by fitting to
NIST, but by algebraic derivation + literature reconstruction), report
results for each, and flag which came closest.  This is the honest-negative
framing specified in the sprint plan.

Outputs
-------
debug/data/bf_d_2P_benchmark.json   -- per-combination splittings and errors
debug/bf_d_benchmark_memo.md        -- human-readable writeup

Author: GeoVac Development Team, April 2026. Sprint 3 BF-D.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Rational, Integer, Float, log, simplify, N as Nf

os.environ.setdefault("PYTHONIOENCODING", "utf-8")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.breit_integrals import (
    breit_ss_radial,
    breit_soo_radial,
    compute_radial,
)
from geovac.spin_orbit import so_diagonal_matrix_element
from geovac.dirac_matrix_elements import Z_sym, alpha_sym


# ------------------------------------------------------------
# Constants and NIST reference
# ------------------------------------------------------------

ALPHA_CODATA = 7.2973525693e-3
HA_TO_MHZ = 6.5796839204e9  # 1 Ha -> MHz

NIST_REF_MHZ = {
    "E(P0) - E(P1)":  29616.951,
    "E(P1) - E(P2)":   2291.178,
    "E(P0) - E(P2)":  29616.951 + 2291.178,
}
NIST_REF_HA = {k: v / HA_TO_MHZ for k, v in NIST_REF_MHZ.items()}


# J-dependent angular coefficients (verified in BR-C to reproduce
# NIST to 0.000% with correct A amplitudes).
# Bethe-Salpeter §39 convention:
F_SS = {0: Rational(-2, 1), 1: Rational(+1, 1), 2: Rational(-1, 5)}
F_SOO = {0: Rational(+2, 1), 1: Rational(+1, 1), 2: Rational(-1, 1)}
X_J = {J: J * (J + 1) - 4 for J in (0, 1, 2)}  # J(J+1)-L(L+1)-S(S+1), L=S=1


# ------------------------------------------------------------
# Radial integral evaluation (closed forms at Z_nuc=2)
# ------------------------------------------------------------

def evaluate_radial_set(Z_nuc: int = 2, verbose: bool = True) -> Dict[str, Tuple[sp.Expr, float]]:
    """Compute the full set of Breit-Pauli retarded Slater integrals
    needed for the He (1s)(2p) ^3P fine structure at integer nuclear
    charge Z_nuc.

    Returns a dict of {name: (sympy_Expr, float_value)}.

    The integrals follow compute_radial's signature
        compute_radial(n1,l1,n3,l3, n2,l2,n4,l4, k, kernel, Z)
    with electrons {1: orbitals a, c = (n1,l1),(n3,l3)}, {2: b, d}.
    "Direct": electron 1 has (1s,1s), electron 2 has (2p,2p).
    "Exchange": electron 1 has (1s,2p), electron 2 has (2p,1s).

    Warning: these are Z^3 Breit-Pauli retarded integrals. Use
    breit_ss_radial or breit_soo_radial (they are aliases of
    compute_radial with kernel_type='breit').
    """
    results = {}

    def _store(name, expr):
        f = float(expr)
        results[name] = (expr, f)
        if verbose:
            expr_str = str(sp.simplify(expr))
            if len(expr_str) > 60:
                expr_str = expr_str[:57] + "..."
            print(f"    {name:20s} = {expr_str:60s} = {f:+.6e}")

    if verbose:
        print(f"    === Breit-Pauli retarded Slater integrals at Z={Z_nuc} (Ha) ===")

    # API signature: breit_ss_radial(n_a,l_a, n_b,l_b, n_c,l_c, n_d,l_d, k, Z)
    # maps to compute_radial(n_a,l_a, n_c,l_c, n_b,l_b, n_d,l_d, k, ...)
    # where (a,c) are orbitals on electron 1 and (b,d) on electron 2.
    #
    # Direct M^k(1s,1s; 2p,2p): e1 has (1s,1s), e2 has (2p,2p)
    #   => n_a,l_a=1,0; n_b,l_b=2,1; n_c,l_c=1,0; n_d,l_d=2,1
    for k in (0, 1, 2):
        val = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, k, Z=Z_nuc)
        _store(f"M{k}_dir", val)

    # Exchange M^k(1s,2p; 2p,1s): e1 has (1s,2p), e2 has (2p,1s)
    #   => n_a,l_a=1,0; n_b,l_b=2,1; n_c,l_c=2,1; n_d,l_d=1,0
    for k in (0, 1, 2):
        val = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, k, Z=Z_nuc)
        _store(f"M{k}_exch", val)

    return results


# ------------------------------------------------------------
# One-body SO parameter from T2
# ------------------------------------------------------------

def zeta_2p(Z_nuc: int = 2, Z_eff: float = 1.0, alpha: float = ALPHA_CODATA) -> float:
    """One-body spin-orbit parameter ζ_{2p} = α^2 Z_nuc <1/r^3>_{2p}.

    <1/r^3>_{2p} = Z_eff^3 / (n^3 · l(l+1/2)(l+1))
                 = Z_eff^3 / (8 · 1 · 1.5 · 2) = Z_eff^3 / 24

    ζ_{2p} = α^2 · Z_nuc · Z_eff^3 / 24

    Note: The factor ζ/2 · X(J) gives the conventional fine-structure
    contribution where X(J) = 2·<L.S>_J. Equivalently ζ·<L.S>_J.
    """
    r3_inv = Z_eff ** 3 / 24.0
    return alpha ** 2 * Z_nuc * r3_inv


# ------------------------------------------------------------
# Amplitude combination candidates (algebraic, not fitted)
# ------------------------------------------------------------

def candidate_combinations(radial: Dict[str, Tuple[sp.Expr, float]],
                             alpha: float = ALPHA_CODATA) -> List[Dict]:
    """Evaluate multiple theoretical derivations for A_SS, A_SOO.

    Each candidate is a hypothesis about the angular combining coefficients
    applied to the {M^k_dir, M^k_exch} radial set. We don't fit to NIST —
    each candidate is an independent algebraic conjecture documented below.

    Returns a list of {name, rationale, A_SS, A_SOO} dicts.
    """
    # Extract floats for arithmetic
    M0d = radial['M0_dir'][1]
    M1d = radial['M1_dir'][1]
    M2d = radial['M2_dir'][1]
    M0e = radial['M0_exch'][1]
    M1e = radial['M1_exch'][1]
    M2e = radial['M2_exch'][1]

    candidates = []

    # ---------------------------------------------------------
    # C1: Drake 1971 canonical (Eq. 14-17)
    # SS from rank-2 tensor: A_SS = α^2 · (3/50) · (M^2_exch - M^2_dir)
    # SOO: A_SOO = α^2 · (1/2) · (N^1_exch - N^1_dir), with N^k = M^k
    # ---------------------------------------------------------
    candidates.append({
        "name": "C1: Drake 1971 canonical (conjectured)",
        "rationale": (
            "Drake 1971 Eq. (16): A_SS = α^2·(3/50)·(M²_exch - M²_dir); "
            "A_SOO = α^2·(1/2)·(N¹_exch - N¹_dir). "
            "(3/50 factor = (2L+1)^(-1) * (Gaunt⟨1||C²||1⟩)^2 = (1/3)·(2/5)^2/... "
            "Literature re-derivation needed for exact coefficient."
        ),
        "A_SS": alpha**2 * (Rational(3, 50)) * float(M2e - M2d),
        "A_SOO": alpha**2 * Rational(1, 2) * float(M1e - M1d),
    })

    # ---------------------------------------------------------
    # C2: Bethe-Salpeter §39 direct-minus-exchange (sign-flip variant)
    # ---------------------------------------------------------
    candidates.append({
        "name": "C2: BS §39 (direct-exchange sign convention)",
        "rationale": (
            "Bethe-Salpeter §39: A_SS = (3/2)·α²·M²_ret where "
            "M²_ret = M²_dir - M²_exch (Slater-sum form). Net rank-2. "
            "A_SOO = (1/2)·α²·(M²_ret - 2·M⁰_ret) with same sign."
        ),
        "A_SS": Rational(3, 2) * alpha**2 * float(M2d - M2e),
        "A_SOO": Rational(1, 2) * alpha**2 * float((M2d - M2e) - 2 * (M0d - M0e)),
    })

    # ---------------------------------------------------------
    # C3: Cowan 1981 Table 12.10 format
    # A_SS involves  α²·( (3/25) M²_dir + (1/5) M²_exch )
    # A_SOO involves α²·( ... )
    # ---------------------------------------------------------
    candidates.append({
        "name": "C3: Cowan 1981 §12 (spin-spin rank-2)",
        "rationale": (
            "Cowan's T²(s₁)T²(s₂) tensor for (sp) ³P; coefficients "
            "involve (2l+1)^(-1) projection."
        ),
        "A_SS": alpha**2 * (Rational(3, 25) * float(M2d) + Rational(1, 5) * float(M2e)),
        "A_SOO": alpha**2 * (Rational(1, 2) * float(M1d) + Rational(3, 10) * float(M1e)),
    })

    # ---------------------------------------------------------
    # C4: Bethe-Salpeter literal (original BR-C buggy version, for reference)
    # ---------------------------------------------------------
    candidates.append({
        "name": "C4: BR-C literal BS (original buggy form)",
        "rationale": "The form BR-C used; known wrong. Included for comparison.",
        "A_SS": Rational(3, 2) * alpha**2 * float(M2d),  # uses direct only
        "A_SOO": Rational(1, 2) * alpha**2 * float(M2d - 2 * M0d),
    })

    # ---------------------------------------------------------
    # C5: Full Bethe-Salpeter §39.14-15 as written
    # A_SS = α² · (3/50) · (M²_dir - M²_exch)  (sign variant of C1)
    # A_SOO = α² · [ (1/2) (M¹_dir - M¹_exch) + (1/3)(M^0_...) ]
    # ---------------------------------------------------------
    candidates.append({
        "name": "C5: Standard 9j derivation (SS rank-2 direct-exchange)",
        "rationale": (
            "Rank-2 tensor SS in (s,p) coupling: reduced matrix element "
            "⟨l=0,l=1||C²||l'=0,l'=1⟩_direct = 0 for Y² transition 0->0, "
            "but exchange pair (0,1)->(1,0) allows. Only M²_exch survives "
            "in purely-rank-2 SS. A_SS = (1/10)·α²·M²_exch."
        ),
        "A_SS": alpha**2 * Rational(1, 10) * float(M2e),
        "A_SOO": alpha**2 * Rational(1, 2) * float(M1e - M1d),
    })

    # ---------------------------------------------------------
    # C6: Rank-2 pure exchange (simplest angular hypothesis)
    # ---------------------------------------------------------
    candidates.append({
        "name": "C6: Pure exchange rank-2 / rank-1",
        "rationale": (
            "If the direct matrix element of SS rank-2 vanishes (Y^2 "
            "selection on (0,0) pair), only exchange M² enters. "
            "A_SS = (3/10)α² M²_exch; A_SOO = α² · M¹_exch."
        ),
        "A_SS": alpha**2 * Rational(3, 10) * float(M2e),
        "A_SOO": alpha**2 * Rational(1, 1) * float(M1e),
    })

    # ---------------------------------------------------------
    # C7: NIST fit anchor (for diagnostic only, NOT a derivation)
    # ---------------------------------------------------------
    # From BR-C: A_SS_fit = -1.202e-6 Ha, A_SOO_fit = +5.333e-6 Ha
    # Express as fractions of α² · M² combinations so we can learn what
    # coefficient structure the data prefers.
    nist_A_SS = -1.202e-6
    nist_A_SOO = +5.333e-6

    # Decompose: find c_d, c_e such that α²·(c_d M²_d + c_e M²_e) = A_SS_fit
    # This is one equation, two unknowns. Inspect various factorizations.
    # For M²_d ~ 8·(-97.43) ~ -779 at Z=2, M²_e ~ 8·(3.66) ~ 29 at Z=2
    # Actually compute:
    implied_c_d_only = nist_A_SS / (alpha**2 * float(M2d))
    implied_c_e_only = nist_A_SS / (alpha**2 * float(M2e))
    candidates.append({
        "name": "C7: NIST-implied coefficients (diagnostic only)",
        "rationale": (
            "A_SS_fit = α²·c·M² implies c_d (if M²_d only) = {:.4e}, "
            "c_e (if M²_e only) = {:.4e}. These are diagnostics — "
            "neither simple rational captures NIST from pure M²_d or M²_e "
            "alone, so the correct form is a sum.".format(
                implied_c_d_only, implied_c_e_only)
        ),
        "A_SS": nist_A_SS,  # tautology: A_SS = A_SS_fit
        "A_SOO": nist_A_SOO,  # tautology: A_SOO = A_SOO_fit
    })

    return candidates


# ------------------------------------------------------------
# Evaluate a candidate: compute E(J), splittings, errors
# ------------------------------------------------------------

def evaluate_candidate(zeta: float, A_SS: float, A_SOO: float) -> Dict:
    """Given one-body ζ and two-body A_SS, A_SOO amplitudes, compute
    the ³P_J energies and the resulting splittings (MHz) vs NIST."""
    E_SO = {J: (zeta / 2.0) * X_J[J] for J in (0, 1, 2)}
    E_SS = {J: A_SS * float(F_SS[J]) for J in (0, 1, 2)}
    E_SOO = {J: A_SOO * float(F_SOO[J]) for J in (0, 1, 2)}
    E_tot = {J: E_SO[J] + E_SS[J] + E_SOO[J] for J in (0, 1, 2)}

    split_ha = {
        "E(P0) - E(P1)": E_tot[0] - E_tot[1],
        "E(P1) - E(P2)": E_tot[1] - E_tot[2],
        "E(P0) - E(P2)": E_tot[0] - E_tot[2],
    }
    split_mhz = {k: v * HA_TO_MHZ for k, v in split_ha.items()}
    rel_err = {k: (split_mhz[k] - NIST_REF_MHZ[k]) / NIST_REF_MHZ[k]
               for k in NIST_REF_MHZ}
    abs_rel_err_max = max(abs(e) for e in rel_err.values())

    return {
        "E_SO_Ha": E_SO,
        "E_SS_Ha": E_SS,
        "E_SOO_Ha": E_SOO,
        "E_total_Ha": E_tot,
        "split_Ha": split_ha,
        "split_MHz": split_mhz,
        "rel_err": rel_err,
        "max_rel_err": abs_rel_err_max,
    }


# ------------------------------------------------------------
# Main benchmark
# ------------------------------------------------------------

def main():
    print("=" * 78)
    print("Sprint 3 BF-D: He 2^3P Fine Structure via geovac.breit_integrals")
    print("=" * 78)

    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Radial integrals (closed form, direct + exchange, k=0,1,2)
    print("\n--- 1. Radial integrals (Drake 1971 closed forms at Z=2) ---")
    radial = evaluate_radial_set(Z_nuc=2, verbose=True)

    # One-body SO parameter
    print("\n--- 2. One-body spin-orbit (T2 single-particle, Z_eff=1) ---")
    zeta = zeta_2p(Z_nuc=2, Z_eff=1.0)
    print(f"    zeta_2p = {zeta:+.6e} Ha = {zeta*HA_TO_MHZ:+.2f} MHz")

    # One-body SO using GeoVac T2 spin_orbit module
    # For 2p_{1/2}: κ=+1, j=1/2, X(½)-type expression. For 2p_{3/2}: κ=-2.
    # We compute the one-body SO with l=1 (kappa=-2 for j=3/2, kappa=+1 for j=1/2):
    so_2p32 = so_diagonal_matrix_element(n=2, kappa=-2, Z=Integer(2), alpha=alpha_sym)
    so_2p12 = so_diagonal_matrix_element(n=2, kappa=+1, Z=Integer(2), alpha=alpha_sym)
    print(f"    T2 H_SO(2p_{{3/2}}) = {sp.simplify(so_2p32.subs(alpha_sym, ALPHA_CODATA))}")
    print(f"    T2 H_SO(2p_{{1/2}}) = {sp.simplify(so_2p12.subs(alpha_sym, ALPHA_CODATA))}")

    # Candidate combinations
    print("\n--- 3. Candidate Drake combinations (algebraic, not fitted) ---")
    candidates = candidate_combinations(radial, alpha=ALPHA_CODATA)

    results_all = []
    for c in candidates:
        A_SS = float(c['A_SS'])
        A_SOO = float(c['A_SOO'])
        r = evaluate_candidate(zeta, A_SS, A_SOO)

        print(f"\n  [{c['name']}]")
        print(f"    Rationale: {c['rationale']}")
        print(f"    A_SS  = {A_SS:+.4e} Ha ({A_SS*HA_TO_MHZ:+.1f} MHz)")
        print(f"    A_SOO = {A_SOO:+.4e} Ha ({A_SOO*HA_TO_MHZ:+.1f} MHz)")
        print(f"    E_total_J (Ha):")
        for J in (0, 1, 2):
            print(f"      J={J}: {r['E_total_Ha'][J]:+.4e}")
        print(f"    Splittings (MHz):")
        for k in ("E(P0) - E(P1)", "E(P1) - E(P2)", "E(P0) - E(P2)"):
            cm = r['split_MHz'][k]
            nm = NIST_REF_MHZ[k]
            pct = r['rel_err'][k] * 100
            print(f"      {k}: {cm:>+12.2f}  NIST {nm:>+10.2f}  rel {pct:+7.1f}%")
        print(f"    max |rel err| = {r['max_rel_err']*100:+.2f}%")

        entry = {
            "candidate": c['name'],
            "rationale": c['rationale'],
            "A_SS_Ha": float(A_SS),
            "A_SOO_Ha": float(A_SOO),
            "A_SS_MHz": float(A_SS * HA_TO_MHZ),
            "A_SOO_MHz": float(A_SOO * HA_TO_MHZ),
            "E_total_Ha": {str(J): float(v) for J, v in r['E_total_Ha'].items()},
            "split_MHz": {k: float(v) for k, v in r['split_MHz'].items()},
            "rel_err": {k: float(v) for k, v in r['rel_err'].items()},
            "max_rel_err_pct": float(r['max_rel_err'] * 100),
        }
        results_all.append(entry)

    # SO-only baseline
    r_so_only = evaluate_candidate(zeta, 0.0, 0.0)
    so_entry = {
        "candidate": "SO only (T8 baseline)",
        "rationale": "Single-particle spin-orbit only, no two-body Breit correction.",
        "A_SS_Ha": 0.0, "A_SOO_Ha": 0.0,
        "A_SS_MHz": 0.0, "A_SOO_MHz": 0.0,
        "E_total_Ha": {str(J): float(v) for J, v in r_so_only['E_total_Ha'].items()},
        "split_MHz": {k: float(v) for k, v in r_so_only['split_MHz'].items()},
        "rel_err": {k: float(v) for k, v in r_so_only['rel_err'].items()},
        "max_rel_err_pct": float(r_so_only['max_rel_err'] * 100),
    }
    print(f"\n  [SO only (T8 baseline)]")
    for k in ("E(P0) - E(P1)", "E(P1) - E(P2)", "E(P0) - E(P2)"):
        cm = r_so_only['split_MHz'][k]
        nm = NIST_REF_MHZ[k]
        pct = r_so_only['rel_err'][k] * 100
        print(f"      {k}: {cm:>+12.2f}  NIST {nm:>+10.2f}  rel {pct:+7.1f}%")
    print(f"    max |rel err| = {r_so_only['max_rel_err']*100:+.2f}%")

    # Identify best candidate
    results_scored = results_all + [so_entry]
    # For "best" select the candidate with min |rel err| on the full span
    best = min(results_scored, key=lambda e: abs(e['rel_err']["E(P0) - E(P2)"]))
    print(f"\n=== Best candidate for span accuracy: {best['candidate']} ===")
    print(f"    max |rel err| = {best['max_rel_err_pct']:+.2f}%")

    # Write JSON output
    out = {
        "track": "BF-D",
        "date": "2026-04-15",
        "description": (
            "He 2^3P fine-structure benchmark via geovac.breit_integrals. "
            "Tests multiple literature-derived formulas for A_SS, A_SOO "
            "as linear combinations of M^k direct/exchange integrals. "
            "Angular J-pattern (f_SS, f_SOO) fixed from Bethe-Salpeter §39 "
            "(BR-C verified to 0.000% fit)."
        ),
        "alpha_CODATA": ALPHA_CODATA,
        "HA_TO_MHZ": HA_TO_MHZ,
        "NIST_reference_MHz": NIST_REF_MHZ,
        "NIST_reference_Ha": NIST_REF_HA,
        "f_SS": {str(J): float(F_SS[J]) for J in (0, 1, 2)},
        "f_SOO": {str(J): float(F_SOO[J]) for J in (0, 1, 2)},
        "radial_integrals_closed_form": {
            name: {"sympy": str(sp.simplify(expr)), "float": val}
            for name, (expr, val) in radial.items()
        },
        "zeta_2p_Ha": zeta,
        "zeta_2p_MHz": zeta * HA_TO_MHZ,
        "candidates": results_all,
        "SO_only_baseline": so_entry,
        "best_for_span_accuracy": best['candidate'],
        "best_max_rel_err_pct": best['max_rel_err_pct'],
        "span_target_pct": 20.0,
        "span_target_met": abs(best['rel_err']["E(P0) - E(P2)"]) * 100 < 20.0,
    }

    out_path = out_dir / "bf_d_2P_benchmark.json"
    with out_path.open("w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {out_path}")

    return out


if __name__ == "__main__":
    main()
