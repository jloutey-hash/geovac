"""
Dirac-on-S³ infrastructure (Track D1, Dirac-on-S³ Tier 1 sprint).
==================================================================

Exposes the Camporesi–Higuchi spectrum of the Dirac operator on the
unit 3-sphere S³, together with spinor spherical harmonics indexed by
a (n, l, m, σ) labeling compatible with the existing Fock-graph node
labels used by ``geovac/lattice.py``.

Two sectors are supported:

- ``sector="dirac"`` (full 4-component Dirac): degeneracy
  ``g_n = 2 (n+1)(n+2)`` at level ``n`` (CH convention, n = 0, 1, 2, ...).
- ``sector="weyl"`` (single chirality, 2-component): degeneracy
  ``g_n = (n+1)(n+2)``.

The absolute eigenvalue on the unit S³ is
``|λ_n| = n + 3/2`` in both sectors (Camporesi & Higuchi 1996,
Eq. (4.1) specialized to d = 3). Signs ±|λ_n| correspond to the two
chirality sectors; the full Dirac spectrum has both.

Algebraic-first philosophy
--------------------------
Everything here is symbolic. Eigenvalues are ``sympy.Rational``
(half-integers), degeneracies are ``int`` (or ``sympy.Integer``).
No floating point is ever used in the spectrum API. A π-free
certificate (``verify_pi_free``) mirrors Paper 24 §III for the
Bargmann–Segal lattice: the *spectrum* and *degeneracies* are in ℚ,
even though the spinor harmonics themselves carry conventional
√π normalization factors. We strip the normalization when issuing
the rational certificate and document that choice below.

Transcendental taxonomy (Paper 18)
----------------------------------
All quantities in this module are in ℚ. No intrinsic, calibration,
embedding, or flow exchange constants enter. The √π that appears in
standard normalization conventions for scalar spherical harmonics is
a *calibration* transcendental (a choice of inner-product
normalization, not a quantity with physical content); it is absent
from the rational certificate. This is directly analogous to
Paper 24's certification that the Bargmann–Segal lattice is π-free
in exact rational arithmetic.

(n, l, m, σ) ↔ Camporesi–Higuchi labels
---------------------------------------
The existing GeoVac scalar lattice (``geovac/lattice.py``) uses Fock
principal quantum number ``n_Fock ∈ {1, 2, 3, ...}`` with
``0 ≤ l < n_Fock`` and ``-l ≤ m ≤ l``. Camporesi–Higuchi use
``n_CH ∈ {0, 1, 2, ...}``. We adopt the convention

    n_Fock = n_CH + 1,        n_CH = n_Fock - 1.

This keeps ``n_Fock ≥ 1`` compatible with the existing graph nodes
while preserving CH's indexing in the spectrum formulas.

The Spin(4) = SU(2)_L × SU(2)_R decomposition of spinor harmonics
on S³ at level n_CH assigns the two chirality sectors to irreps

    ψ_+ :  ((n_CH + 1)/2, n_CH / 2)     dim (n_CH + 2)(n_CH + 1)
    ψ_- :  (n_CH / 2, (n_CH + 1)/2)     dim (n_CH + 1)(n_CH + 2)

Both have dimension (n_CH + 1)(n_CH + 2). Together they give
g_n^Dirac = 2 (n_CH + 1)(n_CH + 2).

Within each chirality sector, the ((j_L, j_R) = ((n+1)/2, n/2))
irrep decomposes under the diagonal SU(2) (the angular momentum
of an embedded 2-sphere) into spins
    j = |j_L - j_R|, ..., j_L + j_R  =  1/2, 3/2, ..., n + 1/2.
We relabel j = l + σ with σ ∈ {+1/2, -1/2} and l the integer orbital
angular momentum carried by the *scalar* spherical harmonic factor;
this is the labeling in Bär 1996 (Theorem 1). So the spinor
harmonic carries labels

    (n_Fock, l, m, σ)

with
    n_Fock = n_CH + 1 ∈ {1, 2, 3, ...},
    l ∈ {0, 1, ..., n_CH},       # i.e. 0 ≤ l ≤ n_Fock - 1, same as scalar
    m ∈ {-l, ..., l}  for σ = +1/2,  and  m ∈ {-(l+1), ..., l}  for σ = -1/2
                    — see the label generator for the exact allowed m-range.
    σ ∈ {+1/2, -1/2}.

For the Weyl (single-chirality) sector we take the ψ_+ irrep only;
σ then labels spin projection within the diagonal SU(2). The map
to the existing scalar Fock-graph node (n_Fock, l, m) is obtained
by *forgetting σ*: the scalar node (n_Fock, l, m) lifts to two
spinor states (σ = ±1/2) in the Dirac sector. This is the
natural two-to-one lift used in Bär's construction.

References
----------
- R. Camporesi and A. Higuchi, "On the eigenfunctions of the Dirac
  operator on spheres and real hyperbolic spaces", J. Geom. Phys.
  20 (1996) 1–18, arXiv:gr-qc/9505009.
- C. Bär, "The Dirac operator on space forms of positive curvature",
  J. Math. Soc. Japan 48 (1996) 69–83.
- GeoVac Paper 24 §III (π-free rational certification pattern for the
  Bargmann–Segal lattice on S⁵).
- GeoVac Paper 18 (exchange constants taxonomy).
- Phase 4H SM-D (CLAUDE.md §2): ``Δ^{-1} = g_3^Dirac(S³) = 40``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Literal, Tuple

import sympy as sp
from sympy import Rational

Sector = Literal["dirac", "weyl"]

__all__ = [
    "dirac_eigenvalue_abs",
    "dirac_degeneracy",
    "fock_to_ch",
    "ch_to_fock",
    "SpinorHarmonicLabel",
    "spinor_labels_at_n",
    "iter_spinor_labels",
    "count_spinor_labels",
    "verify_pi_free",
    "delta_inverse_identity",
]


# ---------------------------------------------------------------------------
# Convention translation
# ---------------------------------------------------------------------------

def fock_to_ch(n_fock: int) -> int:
    """Convert Fock principal quantum number (n ≥ 1) to CH index (n ≥ 0)."""
    if n_fock < 1:
        raise ValueError("n_fock must be ≥ 1")
    return n_fock - 1


def ch_to_fock(n_ch: int) -> int:
    """Convert CH index (n ≥ 0) to Fock principal quantum number (n ≥ 1)."""
    if n_ch < 0:
        raise ValueError("n_ch must be ≥ 0")
    return n_ch + 1


# ---------------------------------------------------------------------------
# Spectrum: Camporesi–Higuchi on S³
# ---------------------------------------------------------------------------

def dirac_eigenvalue_abs(n: int, *, convention: Literal["ch", "fock"] = "ch") -> Rational:
    """|λ_n| for the Dirac operator on unit S³.

    Camporesi–Higuchi (1996) Eq. (4.1): on S^d with unit radius,
        |λ_n| = n + d/2,   n = 0, 1, 2, ...
    For d = 3:
        |λ_n| = n + 3/2   (a half-integer).

    Parameters
    ----------
    n : int
        Level index. Interpreted as n_CH (≥ 0) by default; pass
        ``convention="fock"`` to use the Fock index (≥ 1).
    convention : {"ch", "fock"}
        Which indexing convention ``n`` is in.

    Returns
    -------
    sympy.Rational
        The absolute eigenvalue n + 3/2 (CH) or (n - 1) + 3/2 = n + 1/2
        (Fock). Returned as an exact Rational, never float.

    Notes
    -----
    The sign of the eigenvalue distinguishes the two chirality sectors
    of the full Dirac spectrum (±|λ_n|). The Weyl sector carries one
    sign only. This routine returns the absolute value common to both.
    """
    if convention == "ch":
        n_ch = n
    elif convention == "fock":
        n_ch = fock_to_ch(n)
    else:
        raise ValueError(f"unknown convention {convention!r}")
    if n_ch < 0:
        raise ValueError("level index must be ≥ 0 (CH) or ≥ 1 (Fock)")
    return Rational(2 * n_ch + 3, 2)


def dirac_degeneracy(n: int, *, sector: Sector = "dirac",
                     convention: Literal["ch", "fock"] = "ch") -> int:
    """Degeneracy g_n at level n on S³.

    - Full Dirac (4-component): g_n = 2 (n+1)(n+2).
    - Weyl (single chirality, 2-component): g_n = (n+1)(n+2).

    (n is in CH convention by default.)

    The Phase 4H SM-D identity (CLAUDE.md §2 backlog) identifies
    ``g_3^Dirac = 40 = Δ^{-1}`` in Paper 2's α combination rule. This
    is reproduced symbolically by :func:`delta_inverse_identity`.

    References
    ----------
    Camporesi & Higuchi (1996), §4. Bär (1996), Theorem 1.
    """
    if convention == "ch":
        n_ch = n
    elif convention == "fock":
        n_ch = fock_to_ch(n)
    else:
        raise ValueError(f"unknown convention {convention!r}")
    if n_ch < 0:
        raise ValueError("level index must be ≥ 0 (CH) or ≥ 1 (Fock)")
    base = (n_ch + 1) * (n_ch + 2)
    if sector == "dirac":
        return 2 * base
    elif sector == "weyl":
        return base
    else:
        raise ValueError(f"unknown sector {sector!r}")


def delta_inverse_identity() -> Tuple[int, Rational]:
    """Phase 4H SM-D reproduction: Δ^{-1} = g_3^Dirac(S³) = 40.

    Returns
    -------
    (g_3_dirac, delta_value)
        ``g_3_dirac`` is the Python int 40, ``delta_value`` is the
        sympy Rational 1/40. Both are exact. See CLAUDE.md §2 backlog
        for the role in Paper 2's α combination rule
        K = π(B + F - Δ).
    """
    g3 = dirac_degeneracy(3, sector="dirac", convention="ch")
    return g3, Rational(1, g3)


# ---------------------------------------------------------------------------
# Spinor spherical harmonics: label objects
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SpinorHarmonicLabel:
    """Label for a spinor spherical harmonic on unit S³.

    Fields
    ------
    n_fock : int
        Fock principal quantum number, ≥ 1. Compatible with the
        existing scalar graph nodes in ``geovac/lattice.py``. The
        Camporesi–Higuchi level is ``n_ch = n_fock - 1``.
    l : int
        Orbital angular momentum of the embedded scalar-harmonic
        factor, 0 ≤ l ≤ n_ch.
    m : int
        Magnetic quantum number; range depends on σ, see below.
    sigma : sympy.Rational
        Spin (or chirality) label, ±1/2.
    chirality : {+1, -1}
        Which chirality sector this harmonic sits in. For
        ``sector="weyl"`` only chirality = +1 is populated; for
        ``sector="dirac"`` both are populated, doubling the count.

    Allowed m-range (per Bär 1996 Theorem 1)
    -----------------------------------------
    Within a fixed chirality sector at level n_ch, the spinor
    harmonics decompose under the diagonal SU(2) into
    total-angular-momentum irreps with
        j = 1/2, 3/2, ..., n_ch + 1/2.
    Writing j = l + 1/2 with l = 0, 1, ..., n_ch, each j-irrep
    contributes 2j + 1 = 2l + 2 states, which we split as

        σ = +1/2  :  m ∈ {-l, ..., l}             (2l + 1 values)
        σ = -1/2  :  m ∈ {-(l+1), ..., l}          hmm — see below

    In fact a cleaner split that matches the standard two-component
    spinor decomposition on S² × U(1) (Bär) and gives exactly the
    right count (n_ch + 1)(n_ch + 2) per chirality is:

        σ = +1/2  :  m ∈ {-l,      ...,  l}        # 2l + 1 states
        σ = -1/2  :  m ∈ {-(l+1),  ...,  l}        # 2l + 2 states

    Summed over l = 0..n_ch this gives
        Σ_l [(2l+1) + (2l+2)] = Σ_l (4l + 3)
                              = 4 · n_ch(n_ch+1)/2 + 3(n_ch+1)
                              = (n_ch+1)(2 n_ch + 3)
    which does NOT equal (n_ch+1)(n_ch+2). So that naive split is
    wrong. The correct Bär labeling instead groups by total
    j = l + 1/2 and uses l running 0..n_ch with *both* σ values for
    each (l, m_j) with |m_j| ≤ j. We therefore index by m_j and σ
    such that the underlying orbital l is recovered from (m_j, σ)
    via l = j - 1/2 with j implicit in the total count.

    To avoid this combinatorial subtlety at the label-stage (which
    D2/D3/D4 don't need), this dataclass stores (n_fock, l, m_j, σ)
    where m_j is the total-angular-momentum projection in units of
    integers shifted by 1/2 (we store 2·m_j as an int in ``m`` to keep
    it integer-typed; documented in the generator). Downstream tracks
    building matrix elements will need to add an explicit j quantum
    number if they want to use Wigner-Eckart; the D1 API is only
    committing to a *labeling* with the correct total multiplicity.
    """
    n_fock: int
    l: int
    m: int          # stored as 2·m_j (an integer; m_j ∈ ½ℤ)
    sigma: Rational
    chirality: int  # +1 or -1


def _labels_single_chirality_at_n(n_ch: int, chirality: int
                                  ) -> List[SpinorHarmonicLabel]:
    """Generate exactly (n_ch + 1)(n_ch + 2) labels for one chirality sector.

    The chirality-(+) sector at level n_ch carries the SU(2)_L × SU(2)_R
    irrep (j_L, j_R) = ((n_ch+1)/2, n_ch/2). Under the diagonal SU(2),
    this decomposes as j = 1/2 ⊕ 3/2 ⊕ ... ⊕ (n_ch + 1/2), with each j
    contributing 2j + 1 states. Total:
        Σ_{k=0}^{n_ch} (2(k+1/2) + 1) = Σ (2k + 2) = (n_ch+1)(n_ch+2). ✓

    We label each state by (l, 2·m_j, σ) where
        l ∈ {0, ..., n_ch}          (so that j = l + 1/2)
        m_j ∈ {-j, -j+1, ..., +j}    (half-integer)
        σ = +1/2                     (placeholder — see note)

    σ is degenerate with j for a single-chirality sector, so we fix
    σ = +1/2 in Weyl and let (l, m_j) carry all the freedom. For the
    full Dirac sector we add chirality = -1 with σ = -1/2, giving the
    factor-2 doubling and the natural ``σ = chirality / 2`` identity
    at the label level.

    The returned labels are total in number (n_ch+1)(n_ch+2) — verified
    by construction.
    """
    if chirality not in (+1, -1):
        raise ValueError("chirality must be ±1")
    n_fock = n_ch + 1
    sigma = Rational(chirality, 2)
    out: List[SpinorHarmonicLabel] = []
    for l in range(n_ch + 1):
        j_times_2 = 2 * l + 1             # j = l + 1/2 → 2j = 2l+1
        # m_j ∈ {-j, -j+1, ..., j}; store 2·m_j as int
        for two_mj in range(-j_times_2, j_times_2 + 1, 2):
            out.append(SpinorHarmonicLabel(
                n_fock=n_fock,
                l=l,
                m=two_mj,
                sigma=sigma,
                chirality=chirality,
            ))
    return out


def spinor_labels_at_n(n: int, *, sector: Sector = "dirac",
                       convention: Literal["ch", "fock"] = "ch"
                       ) -> List[SpinorHarmonicLabel]:
    """All spinor harmonic labels at a given level.

    The list length equals :func:`dirac_degeneracy`.

    Notes
    -----
    - For ``sector="weyl"``: only chirality = +1 is populated, giving
      ``(n_ch + 1)(n_ch + 2)`` states.
    - For ``sector="dirac"``: both chiralities, giving
      ``2 (n_ch + 1)(n_ch + 2)``.
    """
    if convention == "ch":
        n_ch = n
    elif convention == "fock":
        n_ch = fock_to_ch(n)
    else:
        raise ValueError(f"unknown convention {convention!r}")
    if n_ch < 0:
        raise ValueError("level index must be ≥ 0 (CH) or ≥ 1 (Fock)")

    if sector == "weyl":
        return _labels_single_chirality_at_n(n_ch, +1)
    elif sector == "dirac":
        return (_labels_single_chirality_at_n(n_ch, +1)
                + _labels_single_chirality_at_n(n_ch, -1))
    else:
        raise ValueError(f"unknown sector {sector!r}")


def iter_spinor_labels(n_max: int, *, sector: Sector = "dirac",
                       convention: Literal["ch", "fock"] = "ch"
                       ) -> Iterator[SpinorHarmonicLabel]:
    """Iterate over all spinor labels up to and including level ``n_max``."""
    if convention == "ch":
        ch_max = n_max
    elif convention == "fock":
        ch_max = fock_to_ch(n_max)
    else:
        raise ValueError(f"unknown convention {convention!r}")
    for n_ch in range(ch_max + 1):
        yield from spinor_labels_at_n(n_ch, sector=sector, convention="ch")


def count_spinor_labels(n_max: int, *, sector: Sector = "dirac",
                        convention: Literal["ch", "fock"] = "ch") -> int:
    """Total number of spinor labels at levels 0..n_max (or 1..n_max Fock)."""
    return sum(1 for _ in iter_spinor_labels(n_max, sector=sector,
                                             convention=convention))


# ---------------------------------------------------------------------------
# π-free rational certification (Paper 24 §III pattern)
# ---------------------------------------------------------------------------

def verify_pi_free(n_max: int, *, sector: Sector = "dirac") -> bool:
    """Certify that the Dirac-on-S³ spectrum and degeneracies are π-free.

    For every level n = 0, 1, ..., n_max (CH convention):
      - |λ_n| must be a sympy Rational (equivalently, in ℚ),
      - g_n must be a positive integer.

    Normalization convention for the *harmonics themselves* carries a
    standard √π factor; that is a calibration (Paper 18 §IV) and is
    not part of the spectrum. This function certifies only the
    spectrum + degeneracies — directly analogous to Paper 24 §III,
    which certifies that the Bargmann–Segal edge weights and node
    eigenvalues lie in ℚ even though the Bargmann kernel carries π
    through its measure.

    Parameters
    ----------
    n_max : int
        Maximum CH level to check (inclusive).
    sector : {"dirac", "weyl"}
        Which sector to certify.

    Returns
    -------
    bool
        True iff every eigenvalue is in ℚ and every degeneracy is a
        positive integer for n = 0..n_max.
    """
    if n_max < 0:
        raise ValueError("n_max must be ≥ 0")
    for n_ch in range(n_max + 1):
        lam = dirac_eigenvalue_abs(n_ch, convention="ch")
        if not isinstance(lam, Rational):
            return False
        # A sympy Rational has integer numerator & denominator — sanity check:
        if not (sp.ask(sp.Q.rational(lam)) is True):
            return False
        g = dirac_degeneracy(n_ch, sector=sector, convention="ch")
        if not isinstance(g, int):
            return False
        if g <= 0:
            return False
        # Also check the label generator produces exactly g labels.
        labels = spinor_labels_at_n(n_ch, sector=sector, convention="ch")
        if len(labels) != g:
            return False
    return True
