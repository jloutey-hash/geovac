"""Breit-Pauli spin-orbit coupling Hamiltonian on the Dirac-on-S³ basis.

Track T2 of the Dirac-on-S³ Tier 2 sprint. Consumes T1 closed-form
matrix elements (``geovac.dirac_matrix_elements``) and assembles the
one-electron spin-orbit operator H_SO as a *diagonal* block in the
(n, κ, m_j) basis.

Physics
-------
Breit-Pauli leading-order spin-orbit operator for a one-electron
hydrogenic system in atomic units:

    H_SO = ξ(r) · L·S ,    ξ(r) = (1/(2 m² c²)) · (1/r) · dV/dr .

For V(r) = −Z/r (Coulomb) and atomic units with α = 1/c,

    ξ(r) = Z · α² / (2 · r³) ,

hence

    H_SO = (Z α² / 2) · (1/r³) · L·S .

In the κ-native basis both (1/r³) [on the radial factor] and L·S
[by Szmytkowski Eq. 2.7 reduction] are diagonal:

    ⟨n, κ, m_j | L·S | n, κ, m_j⟩    = −(κ + 1)/2               (T1)
    ⟨n, l      | 1/r³ | n, l     ⟩   = Z³ / [n³ · l(l+½)(l+1)]  (T1)

Therefore H_SO is diagonal in (n, κ, m_j) with eigenvalue

    H_SO(n, κ) = (Z α² / 2) · [ −(κ+1)/2 ] · Z³ / [n³ · l(l+½)(l+1)]
              =  −Z⁴ α² (κ+1) / [4 · n³ · l(l+½)(l+1)] ,

where l = kappa_to_l(κ). For s-states (l=0, κ=−1 only) the factor
(κ+1) = 0 vanishes and H_SO = 0 identically — the Kramers cancellation.
The ⟨1/r³⟩ divergence at l=0 is never evaluated.

Paper 18 taxonomy
-----------------
The α² prefactor is spinor-intrinsic content: first-order in α², no
transcendentals, no π. The radial factor is a closed rational in
(Z, n, l). T5 formalizes the classification; see the design memo.

Scope
-----
This is the *non-relativistic-limit* Breit-Pauli SO. The full Dirac
operator carries γ = √(1 − (Zα)²) relativistic corrections in the
radial factor (via Martínez-y-Romero recursions); those are deferred.
At Z α ≪ 1 (Z ≲ 30) the error is < 1%.

See ``docs/spin_orbit_design_memo.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Union

import sympy as sp
from sympy import Expr, Integer, Rational, Symbol

from geovac.dirac_matrix_elements import (
    DiracLabel,
    alpha_sym,
    Z_sym,
    angular_matrix_L_dot_S,
    inverse_r_cubed_hydrogenic,
    iter_dirac_labels,
    kappa_to_l,
)


__all__ = [
    "so_diagonal_matrix_element",
    "build_so_hamiltonian_block",
    "verify_z4_scaling",
]


# ---------------------------------------------------------------------------
# Core closed-form matrix element
# ---------------------------------------------------------------------------


def so_diagonal_matrix_element(
    n: int,
    kappa: int,
    Z=Z_sym,
    alpha: Expr = None,
) -> Expr:
    """Closed-form Breit-Pauli spin-orbit matrix element ⟨n, κ, m_j | H_SO | n, κ, m_j⟩.

    H_SO is diagonal in (n, κ, m_j) and independent of m_j (L·S commutes
    with J_z; ⟨1/r³⟩ is m-independent).

    Parameters
    ----------
    n : int
        Principal quantum number, ≥ 1.
    kappa : int
        Dirac κ, nonzero integer. κ<0 for j=l+1/2; κ>0 for j=l−1/2.
    Z : sympy Expr or int, default ``Z_sym``
        Nuclear charge. Pass an integer for a concrete numeric answer
        (still symbolic in α unless α is also given).
    alpha : sympy Expr or None
        Fine-structure constant. Default is ``alpha_sym`` (symbolic).

    Returns
    -------
    sympy Expr
        Closed-form eigenvalue in Hartree. Zero for l=0 (Kramers).

    Notes
    -----
    For l=0 states (κ=−1 only) this returns exactly 0 — the κ+1
    factor from ⟨L·S⟩ = −(κ+1)/2 vanishes. The divergent ⟨1/r³⟩
    is *not* evaluated; the algebra short-circuits correctly.
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

    # Kramers: l=0 ⇒ κ=−1 ⇒ (κ+1)=0. Return 0 without touching ⟨1/r³⟩.
    if l == 0:
        return Integer(0)

    # L·S diagonal eigenvalue, from T1:
    ls_val = Rational(-(kappa + 1), 2)  # = −(κ+1)/2

    # ⟨1/r³⟩ diagonal, from T1 (closed Bethe-Salpeter form):
    r3_inv = inverse_r_cubed_hydrogenic(n, l, Z=Z)

    # Prefactor: Z·α²/2 times L·S and ⟨1/r³⟩.
    # Note: r3_inv already carries a Z³ factor, so the full result is Z⁴.
    return sp.simplify((Z * alpha**2 / Integer(2)) * ls_val * r3_inv)


# ---------------------------------------------------------------------------
# Hamiltonian block assembly
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SOHamiltonianBlock:
    """Diagonal Breit-Pauli spin-orbit Hamiltonian over a finite (n ≤ n_max)
    set of Dirac labels.

    The operator is diagonal in (n, κ, m_j); the dict ``diag`` maps each
    DiracLabel to its eigenvalue (sympy Expr, symbolic in α and Z unless
    bound).
    """
    n_max: int
    Z: Expr
    alpha: Expr
    diag: Dict[DiracLabel, Expr]

    def __len__(self) -> int:
        return len(self.diag)

    def eigenvalues(self) -> list:
        return list(self.diag.values())

    def as_vector(self, labels: Optional[Iterable[DiracLabel]] = None) -> list:
        if labels is None:
            return list(self.diag.values())
        return [self.diag[lab] for lab in labels]


def build_so_hamiltonian_block(
    n_max: int,
    Z=Z_sym,
    alpha: Expr = None,
) -> SOHamiltonianBlock:
    """Build the Breit-Pauli SO Hamiltonian block for all Dirac labels at n ≤ n_max.

    Every matrix element is computed in closed form. No quadrature, no
    diagonalization — H_SO is already diagonal in the κ-basis.

    Parameters
    ----------
    n_max : int
        Upper bound on the principal quantum number.
    Z : sympy Expr or int, default Z_sym
        Nuclear charge.
    alpha : sympy Expr or None
        Fine-structure constant; default symbolic ``alpha_sym``.

    Returns
    -------
    SOHamiltonianBlock
        Mapping from every DiracLabel (n ≤ n_max) to its diagonal
        eigenvalue. l=0 entries are exactly 0.
    """
    if n_max < 1:
        raise ValueError(f"n_max must be ≥ 1, got {n_max}")
    if alpha is None:
        alpha = alpha_sym

    diag: Dict[DiracLabel, Expr] = {}
    for lab in iter_dirac_labels(n_max):
        diag[lab] = so_diagonal_matrix_element(lab.n_fock, lab.kappa, Z=Z, alpha=alpha)

    return SOHamiltonianBlock(n_max=n_max, Z=Z, alpha=alpha, diag=diag)


# ---------------------------------------------------------------------------
# Scaling verification helper
# ---------------------------------------------------------------------------


def verify_z4_scaling(
    n: int, kappa: int, Z_values: Iterable[int], alpha: Expr = None,
) -> Dict[int, Expr]:
    """Evaluate H_SO(n, κ) at each Z ∈ ``Z_values`` and return the mapping.

    For any (n, κ) with l ≥ 1 the ratio H_SO(Z) / H_SO(Z_ref) must equal
    (Z/Z_ref)⁴ exactly (symbolic). Test helper for T2 regression tests.

    Parameters
    ----------
    n, kappa : int
        Quantum numbers.
    Z_values : iterable of int
        Nuclear charges to evaluate at. Intended use: Z ∈ {1, 3, 4, 38}.
    alpha : sympy Expr or None
        Default symbolic α.

    Returns
    -------
    dict[int, sympy Expr]
        {Z_i: H_SO(n, κ, Z_i)}.
    """
    if alpha is None:
        alpha = alpha_sym
    return {Z: so_diagonal_matrix_element(n, kappa, Z=Integer(Z), alpha=alpha)
            for Z in Z_values}
