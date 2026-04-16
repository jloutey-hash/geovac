"""Complete α⁴ (Breit-Pauli) fine-structure corrections for hydrogenic atoms.

Track T8 of the Dirac-on-S³ Tier 3 sprint. Extends T2's spin-orbit
Hamiltonian (``geovac.spin_orbit``) with the Darwin and mass-velocity
corrections to form the full α⁴ fine-structure ladder.

Physics
-------
The full α⁴ fine-structure correction has three one-body diagonal terms:

1. **Spin-orbit** (T2, already in ``geovac.spin_orbit``):
       E_SO = −Z⁴α²(κ+1) / [4n³ l(l+½)(l+1)]    (l ≥ 1; zero for l = 0)

2. **Darwin** (Zitterbewegung / contact interaction):
       E_D = Z⁴α² / (2n³)  · δ_{l,0}
   Nonzero ONLY for s-states (l = 0, κ = −1). For l ≥ 1 it is exactly zero.
   This is the complement of SO: SO lives on l ≥ 1, Darwin on l = 0.

3. **Mass-velocity** (relativistic kinetic energy correction):
       E_MV = −(Z⁴α²) / (2n⁴) · [n/(l+½) − ¾]
   Applies to ALL (n, l) states.

Combined Dirac fine-structure formula:
    E_FS = E_SO + E_D + E_MV = −(Z⁴α²) / (2n⁴) · [n/(j+½) − ¾]

The combined result depends on (n, j) only — not on l separately.
This is the "accidental" Dirac degeneracy: 2s_{1/2} and 2p_{1/2}
have the same fine-structure energy shift despite different l values.

Paper 18 taxonomy
-----------------
All three terms are algebraic in (Z, n, l, α) — no π, no
transcendentals beyond α itself. The Darwin and mass-velocity terms
are α²-order corrections to the kinetic energy operator and the
contact interaction respectively. They are spinor-intrinsic content
in the Paper 18 classification.

See ``docs/spin_orbit_design_memo.md`` and ``docs/dirac_t8_memo.md``.
"""

from __future__ import annotations

from typing import Dict, Iterable, Optional

import sympy as sp
from sympy import Expr, Integer, Rational, Symbol, simplify

from geovac.dirac_matrix_elements import (
    DiracLabel,
    alpha_sym,
    Z_sym,
    iter_dirac_labels,
    kappa_to_j,
    kappa_to_l,
)
from geovac.spin_orbit import so_diagonal_matrix_element


__all__ = [
    "darwin_diagonal",
    "mass_velocity_diagonal",
    "fine_structure_total",
    "verify_dirac_formula",
]


# ---------------------------------------------------------------------------
# Darwin term
# ---------------------------------------------------------------------------


def darwin_diagonal(
    n: int,
    kappa: int,
    Z=Z_sym,
    alpha: Expr = None,
) -> Expr:
    """Closed-form Darwin correction ⟨n, κ | H_D | n, κ⟩.

    E_D = Z⁴α² / (2n³)   for l = 0 (κ = −1)
    E_D = 0                for l ≥ 1

    The Darwin term is the complement of the spin-orbit term: it is
    nonzero precisely where SO is zero (s-states), and zero where SO
    is nonzero (l ≥ 1 states).

    Parameters
    ----------
    n : int
        Principal quantum number, ≥ 1.
    kappa : int
        Dirac κ, nonzero integer.
    Z : sympy Expr or int, default ``Z_sym``
        Nuclear charge.
    alpha : sympy Expr or None
        Fine-structure constant. Default is ``alpha_sym`` (symbolic).

    Returns
    -------
    sympy Expr
        Darwin energy shift in Hartree. Zero for l ≥ 1.
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1, got {n}")
    if kappa == 0:
        raise ValueError("κ = 0 is not allowed")

    l = kappa_to_l(kappa)
    if l >= n:
        raise ValueError(f"κ={kappa} implies l={l}, but l must be < n={n}")

    if alpha is None:
        alpha = alpha_sym

    # Darwin is nonzero ONLY for s-states (l = 0).
    if l != 0:
        return Integer(0)

    # E_D = Z⁴ α² / (2 n³)
    Z_expr = Z if isinstance(Z, Expr) else Integer(Z)
    return Z_expr**4 * alpha**2 / (Integer(2) * Integer(n)**3)


# ---------------------------------------------------------------------------
# Mass-velocity term
# ---------------------------------------------------------------------------


def mass_velocity_diagonal(
    n: int,
    kappa: int,
    Z=Z_sym,
    alpha: Expr = None,
) -> Expr:
    """Closed-form mass-velocity correction ⟨n, κ | H_MV | n, κ⟩.

    E_MV = −(Z⁴α²) / (2n⁴) · [n/(l+½) − ¾]

    This is the relativistic kinetic energy correction −p⁴/(8m³c²).
    It applies to ALL (n, l) states, including s-states.

    Parameters
    ----------
    n : int
        Principal quantum number, ≥ 1.
    kappa : int
        Dirac κ, nonzero integer.
    Z : sympy Expr or int, default ``Z_sym``
        Nuclear charge.
    alpha : sympy Expr or None
        Fine-structure constant. Default is ``alpha_sym`` (symbolic).

    Returns
    -------
    sympy Expr
        Mass-velocity energy shift in Hartree. Always negative (lowers energy).
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1, got {n}")
    if kappa == 0:
        raise ValueError("κ = 0 is not allowed")

    l = kappa_to_l(kappa)
    if l >= n:
        raise ValueError(f"κ={kappa} implies l={l}, but l must be < n={n}")

    if alpha is None:
        alpha = alpha_sym

    Z_expr = Z if isinstance(Z, Expr) else Integer(Z)
    n_sym = Integer(n)
    l_plus_half = Rational(2 * l + 1, 2)

    # E_MV = −(Z⁴α²) / (2n⁴) · [n/(l+½) − ¾]
    bracket = n_sym / l_plus_half - Rational(3, 4)
    return -Z_expr**4 * alpha**2 * bracket / (Integer(2) * n_sym**4)


# ---------------------------------------------------------------------------
# Total fine-structure (SO + Darwin + MV)
# ---------------------------------------------------------------------------


def fine_structure_total(
    n: int,
    kappa: int,
    Z=Z_sym,
    alpha: Expr = None,
) -> Expr:
    """Total α⁴ fine-structure correction E_SO + E_D + E_MV.

    Should equal the Dirac formula:
        E_FS = −(Z⁴α²) / (2n⁴) · [n/(j+½) − ¾]

    Parameters
    ----------
    n : int
        Principal quantum number, ≥ 1.
    kappa : int
        Dirac κ, nonzero integer.
    Z : sympy Expr or int, default ``Z_sym``
        Nuclear charge.
    alpha : sympy Expr or None
        Fine-structure constant. Default is ``alpha_sym`` (symbolic).

    Returns
    -------
    sympy Expr
        Total fine-structure energy shift in Hartree.
    """
    if alpha is None:
        alpha = alpha_sym

    e_so = so_diagonal_matrix_element(n, kappa, Z=Z, alpha=alpha)
    e_d = darwin_diagonal(n, kappa, Z=Z, alpha=alpha)
    e_mv = mass_velocity_diagonal(n, kappa, Z=Z, alpha=alpha)

    return sp.nsimplify(sp.expand(e_so + e_d + e_mv), rational=False)


# ---------------------------------------------------------------------------
# Dirac formula verification
# ---------------------------------------------------------------------------


def _dirac_formula(
    n: int,
    j: Rational,
    Z=Z_sym,
    alpha: Expr = None,
) -> Expr:
    """Standard Dirac fine-structure formula.

    E_FS = −(Z⁴α²) / (2n⁴) · [n/(j+½) − ¾]

    Parameters
    ----------
    n : int
        Principal quantum number.
    j : sympy Rational
        Total angular momentum (half-integer).
    Z : sympy Expr or int
        Nuclear charge.
    alpha : sympy Expr or None
        Fine-structure constant.

    Returns
    -------
    sympy Expr
    """
    if alpha is None:
        alpha = alpha_sym

    Z_expr = Z if isinstance(Z, Expr) else Integer(Z)
    n_sym = Integer(n)
    j_plus_half = j + Rational(1, 2)

    bracket = n_sym / j_plus_half - Rational(3, 4)
    return -Z_expr**4 * alpha**2 * bracket / (Integer(2) * n_sym**4)


def verify_dirac_formula(
    n_max: int,
    Z=Integer(1),
    alpha: Expr = None,
) -> dict:
    """Verify E_SO + E_D + E_MV = Dirac formula for all (n, l, j) up to n_max.

    For each state (n, κ), computes the three-term sum and the Dirac
    formula, and checks they agree symbolically (via sympy simplify).

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.
    Z : sympy Expr or int, default 1
        Nuclear charge.
    alpha : sympy Expr or None
        Fine-structure constant (default symbolic).

    Returns
    -------
    dict
        Keys: 'all_match' (bool), 'n_states' (int), 'failures' (list),
        'details' (list of dicts with n, kappa, l, j, e_total, e_dirac, match).
    """
    if alpha is None:
        alpha = alpha_sym

    details = []
    failures = []

    for n in range(1, n_max + 1):
        for l in range(n):
            # Enumerate kappa values for this l
            kappa_list = [-(l + 1)]
            if l >= 1:
                kappa_list.append(l)

            for kappa in kappa_list:
                j = kappa_to_j(kappa)

                e_so = so_diagonal_matrix_element(n, kappa, Z=Z, alpha=alpha)
                e_d = darwin_diagonal(n, kappa, Z=Z, alpha=alpha)
                e_mv = mass_velocity_diagonal(n, kappa, Z=Z, alpha=alpha)

                e_total = sp.expand(e_so + e_d + e_mv)
                e_dirac = sp.expand(_dirac_formula(n, j, Z=Z, alpha=alpha))

                diff = simplify(e_total - e_dirac)
                match = (diff == 0)

                entry = {
                    'n': n, 'kappa': kappa, 'l': l,
                    'j': str(j),
                    'e_so': str(e_so),
                    'e_darwin': str(e_d),
                    'e_mv': str(e_mv),
                    'e_total': str(e_total),
                    'e_dirac': str(e_dirac),
                    'match': match,
                }
                details.append(entry)

                if not match:
                    failures.append(entry)

    return {
        'all_match': len(failures) == 0,
        'n_states': len(details),
        'failures': failures,
        'details': details,
    }
