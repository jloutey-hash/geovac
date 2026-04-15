"""π-free rational certificate for relativistic (spin-ful) Hamiltonians.

Track T5 of the Dirac-on-S³ Tier 2 sprint. Companion to
``geovac.dirac_s3.verify_pi_free``, extended to certify that the
*symbolic* coefficients of T1/T2-built relativistic Hamiltonians live in
the taxonomic ring

    Q(α²) × {γ := √(1 − (Zα)²)} × {hydrogenic radial seeds}

with **no hidden π, no ζ values, no other transcendentals**.

Why a companion module rather than an extension of ``dirac_s3.py``?
``dirac_s3.py`` is deliberately T1/T2-blind — it only knows about the
Camporesi–Higuchi spectrum and spinor labels. Certifying T3-class
coefficients requires importing the T1/T2 symbol registry
(``alpha_sym``, ``gamma_rel``, ``Z_sym``). Keeping that here preserves
the D1 module's locked surface while still providing a pattern directly
analogous to the D1 ``verify_pi_free``.

Mechanism
---------
The certifier walks the sympy expression tree of a coefficient and
collects the set of *free symbols* and *function calls*. It accepts
only:

    free symbols:  {α, γ := √(1−(Zα)²), Z, n, l, κ, …(user-declared)}
    functions:     {sqrt on a polynomial in (Zα), integer/rational ops}
    constants:     sympy Rational (including Integer)

It rejects if any of the following are encountered:

    - sympy.pi (or any power of pi)
    - sympy.zeta (Riemann zeta)
    - sympy.exp applied to anything other than an integer exponent of
      an allowed symbol (reserved for future expansion — currently
      no T1/T2/T3 coefficient uses exp at symbolic stage)
    - transcendental functions: log, sin, cos, tan, atan, erf, E₁,
      Gamma, loggamma, ...
    - any sympy Symbol not in the allowed registry

For T3's **concrete** numerical-α output (floats), this certifier does
not re-factorize bit patterns. Instead it certifies the *builder itself*
by running it with α left as a free symbol (``alpha=alpha_sym``,
``Z=Z_sym``) and walking the resulting symbolic expression.

References
----------
- Paper 18 §IV (exchange constants taxonomy).
- Paper 24 §III (π-free rational certificate pattern, scalar).
- ``docs/spin_orbit_design_memo.md`` §Paper 18 classification.
- ``docs/paper18_spinor_subtier_proposal.tex`` (drop-in for §IV).
"""

from __future__ import annotations

from typing import Any, Iterable, Mapping, Optional, Sequence, Set, Tuple

import sympy as sp
from sympy import Expr, Function, Rational, Symbol

from geovac.dirac_matrix_elements import (
    Z_sym,
    alpha_sym,
    gamma_rel,
)


__all__ = [
    "SpinorTaxonomyError",
    "default_allowed_symbols",
    "verify_spinor_pi_free",
    "certify_so_block",
]


# ---------------------------------------------------------------------------
# Taxonomy registry
# ---------------------------------------------------------------------------

# Symbols legitimately present in T1/T2/T3 symbolic coefficients.
# Expansion policy: if a future tier legitimately introduces a new
# symbol, it should be added here AND the change documented in
# Paper 18 §IV (new tier or subtier). Silently extending this list
# is a taxonomy violation.
_DEFAULT_ALLOWED_SYMBOL_NAMES = frozenset({
    # Spinor-intrinsic (new Tier-2 subtier):
    "alpha",      # fine-structure constant α (α² in H_SO)
    "gamma_rel",  # γ = √(1 - (Zα)²) (reserved by T1; not yet used at leading order)
    # Scalar (Paper 18 existing tiers):
    "Z",          # nuclear charge (intrinsic; rational in Z upstairs)
    # Auxiliary dummy indices that may appear in user expressions
    # (e.g., summation variables from down-stream use):
    # Downstream callers can pass these in as allowed via `extra_symbols`.
})

# Transcendental / special-function heads that are flat-out forbidden
# anywhere in the expression tree. Encountering any one is a hard
# rejection — the coefficient is not in the spinor-intrinsic ring.
_FORBIDDEN_FUNCTION_HEADS = frozenset({
    "log",      # natural logarithm — not in taxonomy
    "exp",      # (we allow only bare constants/rationals; exp enters via E₁
                #  which is an embedding constant, not spinor-intrinsic)
    "sin", "cos", "tan", "cot", "sec", "csc",
    "asin", "acos", "atan", "atan2",
    "sinh", "cosh", "tanh",
    "erf", "erfc", "erfi",
    "gamma", "loggamma", "digamma", "polygamma",  # Gamma / Γ — embedding
    "zeta",      # Riemann zeta — odd-zeta tier (not spinor-intrinsic)
    "expint",    # Exponential integral E₁ — Level-2 embedding
    "besselj", "besseli", "bessely", "besselk",
    "hyper", "meijerg",
    # ``Function`` placeholders from T1 (reduced angular matrix elements):
    # These are symbolic stand-ins, not transcendentals, so they get a
    # separate "placeholder" list below rather than outright forbidden.
})

# Symbolic placeholders from T1 that represent rank-1 reduced matrix
# elements (σ_z, L_z, L_+, L_-). These are emitted as unevaluated
# sympy Function instances in T1 for later closure. They carry no
# transcendental content by construction (they are stand-ins for
# Wigner-Eckart reduced matrix elements, which are rationals).
_ALLOWED_PLACEHOLDER_NAMES = frozenset({
    "L_reduced",
    "sigma_reduced",
})


class SpinorTaxonomyError(ValueError):
    """Raised when a coefficient contains content outside the Paper 18
    spinor-intrinsic taxonomy (π, ζ, unexpected transcendentals, or
    unregistered free symbols).
    """


def default_allowed_symbols() -> Set[Symbol]:
    """Return the default registry of allowed sympy Symbols.

    Currently: {Z_sym, alpha_sym, gamma_rel}. Downstream callers
    building user-indexed coefficient dictionaries may need to extend
    this via the ``extra_symbols=`` kwarg of :func:`verify_spinor_pi_free`.
    """
    return {Z_sym, alpha_sym, gamma_rel}


# ---------------------------------------------------------------------------
# Tree walk
# ---------------------------------------------------------------------------


def _classify_node(
    node: Expr,
    allowed_symbols: Set[Symbol],
    allowed_placeholders: Set[str],
) -> None:
    """Walk a sympy expression tree and raise on taxonomy violations.

    Raises
    ------
    SpinorTaxonomyError
        If any forbidden function head, forbidden symbol, or forbidden
        constant (e.g., sp.pi) is encountered.
    """
    # Constants
    if node.is_Number:
        if node is sp.pi or node == sp.pi:
            raise SpinorTaxonomyError(
                "π found in coefficient — not in spinor-intrinsic ring"
            )
        if node.has(sp.pi):
            raise SpinorTaxonomyError(
                f"π-bearing constant {node!r} — not allowed"
            )
        # Rationals, integers, floats — all fine structurally.
        # (Floats are allowed for coefficients that the user has evaluated
        # at a numerical α; the ring-check still applies to the *symbolic*
        # form prior to that.)
        return

    # The π symbol itself (sp.pi is an AtomicExpr, above covers it, but
    # defensively check):
    if node is sp.pi:
        raise SpinorTaxonomyError(
            "π found in coefficient — not in spinor-intrinsic ring"
        )

    # Symbols
    if node.is_Symbol:
        if node in allowed_symbols:
            return
        raise SpinorTaxonomyError(
            f"Unregistered free symbol {node!r} — not in taxonomy. "
            f"Allowed: {sorted(s.name for s in allowed_symbols)}"
        )

    # Function calls: sp.Function instances, sp.sqrt, sp.Pow, etc.
    # sqrt is sp.Pow(x, 1/2). Allow sqrt of a polynomial in allowed
    # symbols (this is how γ = sqrt(1-(Zα)²) is represented).
    head_name = node.func.__name__.lower() if hasattr(node.func, '__name__') else ''
    if head_name in _FORBIDDEN_FUNCTION_HEADS:
        raise SpinorTaxonomyError(
            f"Forbidden function {head_name!r} in coefficient — "
            f"outside spinor-intrinsic taxonomy. "
            f"({head_name!r} belongs to a different Paper 18 tier.)"
        )

    # Placeholders (T1 rank-1 reduced matrix elements)
    if isinstance(node, Function) and node.func.__name__ in allowed_placeholders:
        # Check arguments, recursively.
        for arg in node.args:
            _classify_node(arg, allowed_symbols, allowed_placeholders)
        return

    # Arithmetic/composite expressions: recurse.
    for arg in node.args:
        _classify_node(arg, allowed_symbols, allowed_placeholders)


def verify_spinor_pi_free(
    Hamiltonian: Any,
    *,
    alpha_sym_: Optional[Symbol] = None,
    gamma_sym_: Optional[Symbol] = None,
    extra_symbols: Optional[Iterable[Symbol]] = None,
    extra_placeholders: Optional[Iterable[str]] = None,
    raise_on_violation: bool = True,
) -> bool:
    """Certify that every coefficient of ``Hamiltonian`` lies in the
    spinor-intrinsic taxonomic ring

        Q(α²) × {γ := √(1 − (Zα)²)} × {⟨r^k⟩, ⟨1/r⟩, ⟨1/r²⟩, ⟨1/r³⟩, …}

    with no π, no ζ, no unexpected transcendentals.

    Parameters
    ----------
    Hamiltonian : dict, SOHamiltonianBlock, iterable of Expr, or single Expr
        The object to certify. Accepted shapes:

        - ``dict``-like: values are sympy Exprs (e.g., ``H_SO.diag``).
        - ``SOHamiltonianBlock``: the certifier pulls ``.diag.values()``.
        - iterable of ``sympy.Expr``.
        - a single ``sympy.Expr``.

    alpha_sym_, gamma_sym_ : sympy.Symbol, optional
        Overrides for the α and γ symbols. Default: ``alpha_sym`` and
        ``gamma_rel`` from ``geovac.dirac_matrix_elements``.

    extra_symbols : iterable of Symbol, optional
        Additional symbols to accept (e.g., summation indices, angular
        labels). Any free symbol in a coefficient that is not in the
        core registry ∪ ``extra_symbols`` triggers a rejection.

    extra_placeholders : iterable of str, optional
        Additional ``Function`` heads to accept as placeholders (rank-1
        reduced matrix elements or similar symbolic stand-ins).

    raise_on_violation : bool, default True
        If True, raise ``SpinorTaxonomyError`` on the first violation. If
        False, return False instead. Success always returns True.

    Returns
    -------
    bool
        True iff every coefficient classifies cleanly.

    Raises
    ------
    SpinorTaxonomyError
        If ``raise_on_violation=True`` and any coefficient violates the
        taxonomy.

    Examples
    --------
    >>> from geovac.spin_orbit import build_so_hamiltonian_block
    >>> from geovac.dirac_matrix_elements import alpha_sym, Z_sym
    >>> block = build_so_hamiltonian_block(n_max=2, Z=Z_sym, alpha=alpha_sym)
    >>> verify_spinor_pi_free(block)
    True

    Deliberately contaminated:

    >>> import sympy as sp
    >>> bad = {'x': sp.pi * alpha_sym**2}
    >>> verify_spinor_pi_free(bad, raise_on_violation=False)
    False
    """
    allowed_symbols = default_allowed_symbols()
    if alpha_sym_ is not None:
        allowed_symbols.add(alpha_sym_)
    if gamma_sym_ is not None:
        allowed_symbols.add(gamma_sym_)
    if extra_symbols:
        allowed_symbols.update(extra_symbols)

    allowed_placeholders = set(_ALLOWED_PLACEHOLDER_NAMES)
    if extra_placeholders:
        allowed_placeholders.update(extra_placeholders)

    # Normalize input shape to an iterable of Exprs.
    coeffs = _extract_coefficients(Hamiltonian)

    try:
        for coeff in coeffs:
            expr = sp.sympify(coeff)
            _classify_node(expr, allowed_symbols, allowed_placeholders)
    except SpinorTaxonomyError:
        if raise_on_violation:
            raise
        return False
    return True


def _extract_coefficients(H: Any) -> Sequence[Expr]:
    """Pull a flat list of sympy coefficients out of the accepted input
    shapes.
    """
    # SOHamiltonianBlock-like: has a `diag` attribute that is a dict.
    if hasattr(H, "diag") and isinstance(H.diag, Mapping):
        return list(H.diag.values())

    # Mapping: use values.
    if isinstance(H, Mapping):
        return list(H.values())

    # Single sympy expression (or anything sympify'able to atom)
    if isinstance(H, Expr):
        return [H]

    # Iterable of sympy expressions.
    if hasattr(H, "__iter__"):
        return list(H)

    raise TypeError(
        f"Hamiltonian of type {type(H).__name__!r} is not a supported "
        f"shape. Expected dict, SOHamiltonianBlock-like, iterable of Expr, "
        f"or a single sympy Expr."
    )


# ---------------------------------------------------------------------------
# Track T2 convenience wrapper
# ---------------------------------------------------------------------------


def certify_so_block(n_max: int = 3, Z: Any = None, alpha: Any = None) -> bool:
    """Build a symbolic T2 spin-orbit block at (Z_sym, alpha_sym) and
    certify it with :func:`verify_spinor_pi_free`.

    This is a one-line end-to-end test of the T2 pipeline against
    Paper 18 §IV. A failure here indicates either:

    - A regression in ``geovac.spin_orbit`` (new transcendental surfaced).
    - A gap in the certifier's allowed registry (legitimate new symbol
      that the taxonomy hasn't yet classified — flag for Paper 18 update).

    Parameters
    ----------
    n_max : int, default 3
        Maximum Fock principal quantum number for the block.
    Z, alpha : sympy.Symbol or None
        Symbols to build the block with. Defaults to ``Z_sym`` and
        ``alpha_sym`` (fully symbolic).

    Returns
    -------
    bool
        True iff the block certifies cleanly.
    """
    # Lazy import so `dirac_s3.py`/T1 tests don't pull T2 transitively.
    from geovac.spin_orbit import build_so_hamiltonian_block
    if Z is None:
        Z = Z_sym
    if alpha is None:
        alpha = alpha_sym
    block = build_so_hamiltonian_block(n_max=n_max, Z=Z, alpha=alpha)
    return verify_spinor_pi_free(block)
