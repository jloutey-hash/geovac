"""Backing tests for Paper 58 (Angular Sparsity Is an Atomic-Sector Property).

Covers the two observations MEASURED for the first time in the Paper-58 work,
which exist in no other backing:

  test_paper58_no_angle
      Obs. "The composed builder carries no angular geometry".  The paper's
      H2O discussion, and the claim that Corollary 1 is not well posed in the
      bent case, both rest on this.  If someone later adds angular geometry to
      the builder, this test SHOULD fail -- that is the point.

  test_paper58_swap_cost
      Obs. "Permutation symmetry costs qubits on this builder".  Enabling
      equivalent-atom-swap tapering REDUCES the qubit saving and inflates the
      Pauli count.  This also pins the documentation discrepancy: the
      docstring of extended_tapered_from_spec advertises "BeH2 +1"; the
      measured value is -1 (Hopf-only baseline) / -3 (Hopf + ell baseline).

Backfill still owed for Paper 58 (census, bra/ket certificate, m-rule theorem,
NaH ladder): those live in exploratory drivers and must be migrated
self-contained into tests/ before the paper is a finished artifact, per the
cite-permanent-records policy.  See the paper's backing table.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

pytest.importorskip("openfermion",
                    reason="tapering pipeline requires openfermion")


# ---------------------------------------------------------------------------
# Obs: the composed builder carries no angular geometry
# ---------------------------------------------------------------------------

def test_paper58_no_angle():
    """molecular_spec encodes no bond angle; H2O is radial-only.

    Three independent legs, so a single cosmetic edit cannot silently
    flip the result:
      (a) no angular vocabulary anywhere in the module source;
      (b) h2o_spec delegates to the generic hydride factory;
      (c) the built spec exposes a single R and no nuclei.
    """
    from geovac import molecular_spec as ms

    src = Path(ms.__file__).read_text(encoding="utf-8", errors="replace")

    # (a) No angular vocabulary. 'angle' as a whole word, degree symbols,
    #     or the water bond angle in either convention.
    angular_hits = re.findall(
        r"\bangle\b|\btheta\b|\bbent\b|104\.[0-9]|\bdeg\b", src, re.IGNORECASE)
    assert angular_hits == [], (
        f"molecular_spec.py now contains angular vocabulary {angular_hits!r}; "
        "Paper 58 Obs. 'no angular geometry' and the not-well-posed status of "
        "Corollary 1 in the bent case both need revisiting."
    )

    # (b) h2o_spec is the generic hydride factory, i.e. Z + one distance.
    h2o_src = src[src.index("def h2o_spec"):]
    h2o_body = h2o_src[:h2o_src.index("def hf_spec")]
    assert "hydride_spec(8" in h2o_body, (
        "h2o_spec no longer delegates to hydride_spec(8); its geometry "
        "parameterization has changed."
    )

    # (c) The built spec has one radial parameter and no nuclear positions.
    spec = ms.h2o_spec()
    assert spec.nuclei is None, (
        "h2o_spec now supplies nuclei; C2v may now be representable and "
        "Paper 58 Corollary 1 may have become testable."
    )
    assert isinstance(spec.R, float) and spec.R > 0
    # Five blocks (core, two bond pairs, two lone pairs), all sharing one R.
    assert len(spec.blocks) == 5, f"expected 5 H2O blocks, got {len(spec.blocks)}"


# ---------------------------------------------------------------------------
# Obs: permutation symmetry costs qubits
# ---------------------------------------------------------------------------

def _linear_symmetric_nuclei(Z_c: float, Z_o: float, R: float):
    """X-A-X on the z-axis.  Linear, so angle-free and unambiguous."""
    return [
        {"Z": float(Z_c), "position": (0.0, 0.0, 0.0), "label": "A"},
        {"Z": float(Z_o), "position": (0.0, 0.0, float(R)), "label": "X1"},
        {"Z": float(Z_o), "position": (0.0, 0.0, -float(R)), "label": "X2"},
    ]


def _n_pauli(qubit_op) -> int:
    return sum(1 for term in qubit_op.terms if term)


def _taper(spec, nuclei, *, ell: bool, swap: bool):
    from geovac.extended_tapering import extended_tapered_from_spec
    out = extended_tapered_from_spec(
        spec, use_hopf=True, use_ell_parity=ell,
        use_atom_swap=swap, use_inversion=False, nuclei=nuclei,
    )
    return out["delta_Q"], _n_pauli(out["qubit_op_tapered"])


@pytest.mark.slow
def test_paper58_swap_cost():
    """Atom-swap tapering reduces delta_Q and inflates Pauli count (BeH2).

    Guards the sign, not a fitted magnitude: swap must make delta_Q strictly
    WORSE and the Pauli count strictly larger, under both baselines.
    """
    from geovac import molecular_spec as ms

    spec = ms.beh2_spec()
    nuclei = _linear_symmetric_nuclei(4.0, 1.0, spec.R)

    # Sanity: the geometry really is swap-eligible, else the test is vacuous.
    from geovac.extended_tapering import (
        find_equivalent_atom_pairs, is_centrosymmetric,
    )
    assert is_centrosymmetric(nuclei), "BeH2 test geometry is not centrosymmetric"
    assert len(find_equivalent_atom_pairs(spec, nuclei)) >= 1, (
        "no equivalent atom pairs found; the swap leg would be vacuous"
    )

    for ell in (False, True):
        dq_base, pauli_base = _taper(spec, nuclei, ell=ell, swap=False)
        dq_swap, pauli_swap = _taper(spec, nuclei, ell=ell, swap=True)

        label = "hopf+ell" if ell else "hopf-only"
        assert dq_swap < dq_base, (
            f"[{label}] atom-swap did NOT cost qubits: "
            f"delta_Q {dq_base} -> {dq_swap}.  Paper 58 Obs. "
            "'permutation symmetry costs qubits' and the documentation "
            "discrepancy it reports would both need revisiting."
        )
        assert pauli_swap > pauli_base, (
            f"[{label}] atom-swap did NOT inflate Pauli count: "
            f"{pauli_base} -> {pauli_swap}"
        )

    # Pin the documentation discrepancy explicitly: the docstring claims +1.
    dq_base, _ = _taper(spec, nuclei, ell=False, swap=False)
    dq_swap, _ = _taper(spec, nuclei, ell=False, swap=True)
    assert dq_swap - dq_base < 0, (
        "extended_tapered_from_spec docstring advertises BeH2 atom-swap as "
        f"+1 qubit; measured {dq_swap - dq_base:+d} on the Hopf-only baseline."
    )


@pytest.mark.slow
def test_paper58_swap_null_control():
    """LiH has no equivalent atoms, so the swap flag must be a no-op.

    Without this, test_paper58_swap_cost could be passing because the swap
    path is broken in general rather than because permutation costs qubits.
    """
    from geovac import molecular_spec as ms

    spec = ms.lih_spec()
    nuclei = [
        {"Z": 3.0, "position": (0.0, 0.0, 0.0), "label": "Li"},
        {"Z": 1.0, "position": (0.0, 0.0, float(spec.R)), "label": "H"},
    ]
    from geovac.extended_tapering import find_equivalent_atom_pairs
    assert len(find_equivalent_atom_pairs(spec, nuclei)) == 0

    dq_base, pauli_base = _taper(spec, nuclei, ell=True, swap=False)
    dq_swap, pauli_swap = _taper(spec, nuclei, ell=True, swap=True)
    assert (dq_base, pauli_base) == (dq_swap, pauli_swap), (
        "swap flag changed the result for a molecule with no equivalent "
        f"atoms: delta_Q {dq_base}->{dq_swap}, Pauli {pauli_base}->{pauli_swap}"
    )


# ===========================================================================
# Thm. 1 (angular selection survives molecularization only in the abelian
# sector).  Two independent legs, matching the two halves of the proof.
#
# Leg (i) -- why m survives: L_z is the SAME operator about every point on the
# internuclear axis, so m is a label the two centers share.  Proved by
# symbolic translation covariance, together with the two contrasts that make
# the hypothesis load-bearing (L_x does NOT survive a shift along z; L_z does
# NOT survive a shift OFF the axis -- the latter is exactly why the bent case
# of Cor. 1 is not merely untested but ill-posed, cf. test_paper58_no_angle).
#
# Leg (ii) -- why l dies: the SAME (l, l') pair that is exactly zero on one
# center is nonzero across two.  Exact rational arithmetic, so vanishing is
# decidable rather than inferred from a tolerance.
#
# Two-center machinery below mirrors debug/compute_topos3_exact_meet.py and is
# gated by its own exact-identity certification (test_paper58_machinery_cert)
# so these results never rest on an unvalidated copy.
# ===========================================================================

from fractions import Fraction  # noqa: E402
from math import comb, factorial  # noqa: E402

sp = pytest.importorskip("sympy", reason="exact rational two-center machinery")

_xi, _eta = sp.symbols("xi eta")


def _radial_poly_coeffs(Z: Fraction, n: int, l: int):
    lam = sp.Rational(2 * Z.numerator, Z.denominator * n)
    lag = [sp.Rational((-1) ** j * comb(n + l, n - l - 1 - j), factorial(j))
           for j in range(n - l)]
    return ({l + j: c * lam ** j for j, c in enumerate(lag)},
            sp.Rational(Z.numerator, Z.denominator * n))


def _radial_norm_sq(Z: Fraction, n: int, l: int):
    coeffs, a = _radial_poly_coeffs(Z, n, l)
    return sp.nsimplify(sum(
        c1 * c2 * sp.factorial(k1 + k2 + 2) / (2 * a) ** (k1 + k2 + 3)
        for k1, c1 in coeffs.items() for k2, c2 in coeffs.items()))


def _theta_norm_sq(l: int, m: int):
    return sp.Rational(2 * factorial(l + m), (2 * l + 1) * factorial(l - m))


def _G_lm(l: int, m: int, x):
    t = sp.Symbol("_t")
    return sp.diff(sp.legendre(l, t), t, m).subs(t, x)


def _two_center_UV(Z1, n1, l1, Z2, n2, l2, m, R, kernel="S"):
    """Two-center matrix element as e^{-p}(U e^{q} + V e^{-q}), U/V exact.

    Element == 0 <=> U == V == 0, by Lindemann independence of the
    exponentials -- so vanishing is DECIDABLE.

    kernel selects the operator between the two orbitals:
      'S'    : identity        -> overlap
      '1/rA' : 1/r_A           -> nuclear attraction to center A
      '1/rB' : 1/r_B           -> nuclear attraction to center B
    In prolate spheroidal coordinates r_A = (R/2)(xi + eta) and
    r_B = (R/2)(xi - eta), while the volume element carries
    (xi^2 - eta^2) = (xi + eta)(xi - eta).  So each Coulomb kernel cancels
    one factor of the volume element exactly, leaving a polynomial -- which
    is why the `den == 1` assertion below is a real check and not a formality.
    """
    m = abs(m)
    Rs = sp.Rational(R.numerator, R.denominator)
    half = Rs / 2
    c1, a = _radial_poly_coeffs(Z1, n1, l1)
    c2, b = _radial_poly_coeffs(Z2, n2, l2)
    ra, rb = half * (_xi + _eta), half * (_xi - _eta)
    ct_a = (1 + _xi * _eta) / (_xi + _eta)
    ct_b = (_xi * _eta - 1) / (_xi - _eta)
    s2 = (_xi ** 2 - 1) * (1 - _eta ** 2)
    geom = {
        "S": half ** 3 * (_xi ** 2 - _eta ** 2),
        "1/rA": half ** 2 * (_xi - _eta),
        "1/rB": half ** 2 * (_xi + _eta),
    }[kernel]
    expr = (sum(c * ra ** k for k, c in c1.items())
            * sum(c * rb ** k for k, c in c2.items())
            * s2 ** m / ((_xi + _eta) * (_xi - _eta)) ** m
            * _G_lm(l1, m, ct_a) * _G_lm(l2, m, ct_b)
            * geom)
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    assert den == 1, f"prolate cancellation failed for kernel {kernel!r}"
    poly = sp.Poly(sp.expand(num), _xi, _eta)
    p = (a + b) * Rs / 2
    q = (a - b) * Rs / 2
    max_i = max(mm[0] for mm in poly.monoms())
    max_j = max(mm[1] for mm in poly.monoms())
    a_c = {0: 1 / p}
    for i in range(1, max_i + 1):
        a_c[i] = (1 + i * a_c[i - 1]) / p
    if q != 0:
        u, v = {0: 1 / q}, {0: -1 / q}
        for j in range(1, max_j + 1):
            u[j] = sp.Rational((-1) ** j) / q + j * u[j - 1] / q
            v[j] = sp.Rational(-1) / q + j * v[j - 1] / q
        U = sum(c * a_c[i] * u[j]
                for (i, j), c in zip(poly.monoms(), poly.coeffs()))
        V = sum(c * a_c[i] * v[j]
                for (i, j), c in zip(poly.monoms(), poly.coeffs()))
        return sp.nsimplify(U), sp.nsimplify(V)
    U = sum(c * a_c[i] * sp.Rational(2, j + 1)
            for (i, j), c in zip(poly.monoms(), poly.coeffs()) if j % 2 == 0)
    return sp.nsimplify(U), None


def _overlap_UV(Z1, n1, l1, Z2, n2, l2, m, R):
    """Overlap: the identity-kernel case of _two_center_UV."""
    return _two_center_UV(Z1, n1, l1, Z2, n2, l2, m, R, kernel="S")


def _hydrogenic_E(Z: Fraction, n: int):
    z = sp.Rational(Z.numerator, Z.denominator)
    return -z ** 2 / (2 * n ** 2)


def _certificate_residuals(Z1, n1, l1, Z2, n2, l2, m, R,
                           *, bra_n=None, swap_kernels=False):
    """Residuals of  E_a S + Z_A I_A  ==  E_b S + Z_B I_B, coefficient-wise.

    The identity is the hydrogenic eigen-trick: the kinetic energy obtained by
    letting (-1/2 grad^2 - Z_A/r_A) act LEFT on the bra must equal the value
    obtained by letting (-1/2 grad^2 - Z_B/r_B) act RIGHT on the ket.  It is
    linear in the bra/ket pair, so any common normalization cancels and
    unnormalized radial functions are legitimate.

    S, I_A and I_B share the same (p, q), so e^{-p}e^{+q} and e^{-p}e^{-q} are
    linearly independent functions of R and the identity must hold on EACH
    coefficient separately -- a strictly stronger check than agreement at one R.

    bra_n / swap_kernels are deliberate corruptions used by the mutation guard.
    """
    S = _two_center_UV(Z1, n1, l1, Z2, n2, l2, m, R, "S")
    IA = _two_center_UV(Z1, n1, l1, Z2, n2, l2, m, R, "1/rA")
    IB = _two_center_UV(Z1, n1, l1, Z2, n2, l2, m, R, "1/rB")
    if swap_kernels:
        IA, IB = IB, IA
    E1 = _hydrogenic_E(Z1, bra_n if bra_n is not None else n1)
    E2 = _hydrogenic_E(Z2, n2)
    z1 = sp.Rational(Z1.numerator, Z1.denominator)
    z2 = sp.Rational(Z2.numerator, Z2.denominator)
    out = []
    for k in (0, 1):
        if S[k] is None:
            continue
        out.append(sp.simplify((E1 * S[k] + z1 * IA[k])
                               - (E2 * S[k] + z2 * IB[k])))
    return out


_CERT_PANEL = [
    (1, 1, 0, 1, 1, 0, 0),
    (1, 1, 0, 1, 2, 1, 0),
    (1, 2, 1, 3, 2, 1, 1),
    (1, 1, 0, 3, 2, 0, 0),
    (1, 2, 0, 3, 3, 2, 0),
    (2, 2, 1, 1, 3, 1, 1),
]


def test_paper58_braket_certificate():
    """Eq. (1): E_a S - Z_B I_B == E_b S - Z_A I_A, exactly, entrywise.

    This is the machinery certificate for the cross-center one-body integrals:
    it fails if S, I_A or I_B is computed wrong, so it is the check that lets
    the census claim genuine (rather than merely self-consistent) integrals.
    """
    R = Fraction(1)
    for Z1, n1, l1, Z2, n2, l2, m in _CERT_PANEL:
        res = _certificate_residuals(Fraction(Z1), n1, l1,
                                     Fraction(Z2), n2, l2, m, R)
        assert res, f"no coefficient legs produced for {(Z1, n1, l1, Z2, n2, l2, m)}"
        assert all(r == 0 for r in res), (
            f"bra/ket certificate violated at "
            f"Z1={Z1}({n1},{l1}) Z2={Z2}({n2},{l2}) m={m}: residuals {res}"
        )


def test_paper58_braket_certificate_is_sensitive():
    """Mutation guard: the certificate must DETECT a corrupted integral.

    Without this, the certificate could be an algebraic tautology that holds
    regardless of whether the integrals are right.  Two independent
    corruptions must both break it.
    """
    R = Fraction(1)
    # Use a case with Z1 != Z2 and n1 != n2 so both corruptions are visible.
    Z1, n1, l1, Z2, n2, l2, m = 1, 1, 0, 3, 2, 0, 0

    # (a) wrong bra eigenvalue (pretend the 1s orbital is a 2s).
    bad = _certificate_residuals(Fraction(Z1), n1, l1, Fraction(Z2), n2, l2,
                                 m, R, bra_n=n1 + 1)
    assert any(r != 0 for r in bad), (
        "certificate did not notice a wrong hydrogenic eigenvalue; it is "
        "tautological and cannot certify the integrals"
    )

    # (b) I_A and I_B interchanged (i.e. attraction to the wrong nucleus).
    bad = _certificate_residuals(Fraction(Z1), n1, l1, Fraction(Z2), n2, l2,
                                 m, R, swap_kernels=True)
    assert any(r != 0 for r in bad), (
        "certificate did not notice I_A and I_B being interchanged"
    )


# ===========================================================================
# Census, Table 1 -- the S and h rows (the COMPUTED rows).
#
# Configuration: Z_A = 3 at the origin, Z_B = 1 at R = 3 zhat, hydrogenic
# n_max = 2 per center, so 5 orbitals per center and M = 10.
#
# The paper's counts decompose exactly, which is what makes them checkable:
#   S = 32 = 10 same-center (orthonormality of same-Z hydrogenics) + 22 cross
#   h = 44 = 22 same-center (off-center nuclear attraction) + 22 cross
# with 22 = 2 blocks x 11 matching-m pairs, since the m multiset {0,0,-1,0,1}
# gives 3*3 + 1*1 + 1*1 = 11.
#
# Mixed exactness, stated: the cross blocks are exact rationals (this file's
# machinery); the same-center h block uses the production float routine
# geovac.shibuya_wulfman.compute_cross_center_vne_element, so its zeros are
# threshold decisions, not decided ones.  The g row is NOT covered here -- it
# is symmetry-rule counting in the paper and no two-center ERI engine exists.
# ===========================================================================

_NMAX2_STATES = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]
_CENSUS_ZA, _CENSUS_ZB, _CENSUS_R = Fraction(3), Fraction(1), Fraction(3)


def _matching_m_pairs():
    return [(a, b) for a in _NMAX2_STATES for b in _NMAX2_STATES
            if a[2] == b[2]]


def test_paper58_census_matching_m_count():
    """The combinatorial backbone of Table 1: 11 matching-m pairs per block."""
    assert len(_NMAX2_STATES) == 5
    assert len(_matching_m_pairs()) == 11, (
        "matching-m pair count changed; Table 1's 22 = 2 x 11 cross-block "
        "figure no longer decomposes as stated"
    )


def test_paper58_census_cross_blocks_exact():
    """Cross-center S and h are nonzero on every matching-m pair, exactly.

    Establishes the 22-entry cross-block contribution to BOTH rows of
    Table 1, in exact rational arithmetic.  h is assembled by the certified
    eigen-trick h = E_a S - Z_B I_B (Eq. 1), so this leg also exercises the
    certificate on the census configuration rather than only on its own panel.
    """
    zb = sp.Rational(_CENSUS_ZB.numerator, _CENSUS_ZB.denominator)
    dead_S, dead_h = [], []
    for (n1, l1, m1), (n2, l2, m2) in _matching_m_pairs():
        S = _two_center_UV(_CENSUS_ZA, n1, l1, _CENSUS_ZB, n2, l2,
                           m1, _CENSUS_R, "S")
        IB = _two_center_UV(_CENSUS_ZA, n1, l1, _CENSUS_ZB, n2, l2,
                            m1, _CENSUS_R, "1/rB")
        if not _nonzero(*S):
            dead_S.append(((n1, l1, m1), (n2, l2, m2)))
        E1 = _hydrogenic_E(_CENSUS_ZA, n1)
        h = tuple(None if S[k] is None else sp.simplify(E1 * S[k] - zb * IB[k])
                  for k in (0, 1))
        if not _nonzero(*h):
            dead_h.append(((n1, l1, m1), (n2, l2, m2)))

    assert dead_S == [], f"cross-center S vanished on matching-m pairs: {dead_S}"
    assert dead_h == [], f"cross-center h vanished on matching-m pairs: {dead_h}"


def test_paper58_census_same_center_h_block():
    """Same-center h block: off-center nuclear attraction, 11 nonzeros.

    Uses the PRODUCTION routine, which independently corroborates
    Thm. 1(i): it returns 0.0 for m1 != m2 by an explicit early return, and
    imposes no condition on l whatsoever.
    """
    sw = pytest.importorskip("geovac.shibuya_wulfman")

    Z_orb, Z_nuc, R = 3.0, 1.0, 3.0
    nonzero_pairs, m_mismatch_nonzero = 0, []
    for (n1, l1, m1) in _NMAX2_STATES:
        for (n2, l2, m2) in _NMAX2_STATES:
            val = sw.compute_cross_center_vne_element(
                Z_orb, n1, l1, m1, n2, l2, m2, Z_nuc, R, L_max=l1 + l2)
            if abs(val) > 1e-10:
                nonzero_pairs += 1
                if m1 != m2:
                    m_mismatch_nonzero.append(((n1, l1, m1), (n2, l2, m2), val))

    assert m_mismatch_nonzero == [], (
        "production cross-center V_ne is nonzero for m1 != m2, contradicting "
        f"Thm. 1(i): {m_mismatch_nonzero}"
    )
    assert nonzero_pairs == 11, (
        f"expected 11 nonzero same-center h entries per center, got "
        f"{nonzero_pairs}; Table 1's h = 44 no longer decomposes as 22 + 22"
    )


def test_paper58_census_totals():
    """Table 1 totals: S = 32 and h = 44, from the verified decomposition.

    Same-center S is 5 per center by orthonormality of same-Z hydrogenic
    orbitals (distinct (n,l,m) are orthogonal; the diagonal is nonzero), which
    is asserted rather than recomputed.  Everything else is checked above.
    """
    cross_per_block = len(_matching_m_pairs())
    same_center_S_per_center = len(_NMAX2_STATES)      # orthonormal diagonal
    same_center_h_per_center = cross_per_block          # matching-m coupling

    S_total = 2 * same_center_S_per_center + 2 * cross_per_block
    h_total = 2 * same_center_h_per_center + 2 * cross_per_block

    assert S_total == 32, f"Table 1 S total is {S_total}, paper says 32"
    assert h_total == 44, f"Table 1 h total is {h_total}, paper says 44"
    assert 2 * len(_NMAX2_STATES) == 10, "M should be 10 for this configuration"


def _nonzero(U, V) -> bool:
    return (U != 0) if V is None else (U != 0 or V != 0)


def test_paper58_nonzero_predicate_can_report_zero():
    """Guard against the vacuity mode: a predicate that always says True.

    The cross-center legs below assert `_nonzero(...)`.  If `_nonzero` could
    never return False those assertions would be vacuous, so pin both
    branches (V is None, and V present) directly.

    Note on scope: exact vanishing of a two-center overlap happens at isolated
    ROOTS in R (the Topos-3/4 vanishing lemma: S = prefactor * P(t), with
    sporadic rational residual roots), not for whole (n, l, m) pairs.  A scan
    over same-Z and rate-coincident cross-Z pairs at R = 1, including
    |dn| up to 4, found no identically-vanishing pair.  The cross-center
    assertions are therefore GENERIC-SUPPORT claims at a stated R, and the
    discriminating control for Thm. 1(ii) is the same-center integral that
    vanishes for the same (l, l') pair -- not a hunted zero.
    """
    assert _nonzero(sp.Integer(0), None) is False
    assert _nonzero(sp.Integer(0), sp.Integer(0)) is False
    assert _nonzero(sp.Integer(1), None) is True
    assert _nonzero(sp.Integer(0), sp.Integer(1)) is True


def test_paper58_machinery_cert():
    """Gate: normalized <1s|1s>(R) equals (1 + R + R^2/3) e^{-R} exactly.

    Without this the two-center legs below would rest on copied, unvalidated
    code.  Exact rational identity at three separations.
    """
    nsq = _radial_norm_sq(Fraction(1), 1, 0) * _theta_norm_sq(0, 0)
    for R, want in {Fraction(1): sp.Rational(7, 3),
                    Fraction(2): sp.Rational(13, 3),
                    Fraction(7, 2): sp.Rational(103, 12)}.items():
        U, V = _overlap_UV(Fraction(1), 1, 0, Fraction(1), 1, 0, 0, R)
        assert V is None, f"equal decay rates should give q=0 at R={R}"
        assert sp.simplify(U / nsq - want) == 0, f"machinery wrong at R={R}"


def test_paper58_m_rule_lz_shared_along_axis():
    """Thm. 1(i): L_z is origin-independent ALONG the axis, and only there.

    Three legs.  The first is the theorem's mechanism; the other two show the
    hypothesis is load-bearing rather than decorative.
    """
    x, y, z, c = sp.symbols("x y z c", real=True)
    f = x ** 2 * y + y * z ** 2 + x * y * z + z ** 3

    def L_z(g):
        return -sp.I * (x * sp.diff(g, y) - y * sp.diff(g, x))

    def L_x(g):
        return -sp.I * (y * sp.diff(g, z) - z * sp.diff(g, y))

    def shift(g, var):
        return g.subs(var, var - c)

    # (1) L_z commutes with translation along z: shifting the origin to the
    #     other nucleus leaves the generator alone, so m is a SHARED label.
    assert sp.simplify(shift(L_z(f), z) - L_z(shift(f, z))) == 0, (
        "L_z failed to commute with translation along z; Thm. 1(i) mechanism "
        "is wrong."
    )

    # (2) L_x does NOT: the non-abelian generators are per-center, which is
    #     why l (needing all of SO(3)) cannot be a shared label.
    assert sp.simplify(shift(L_x(f), z) - L_x(shift(f, z))) != 0, (
        "L_x unexpectedly commuted with a z-translation"
    )

    # (3) L_z does NOT survive a shift OFF the axis.  This is why bent
    #     geometry has no axial m label at all -- cf. test_paper58_no_angle.
    assert sp.simplify(shift(L_z(f), x) - L_z(shift(f, x))) != 0, (
        "L_z unexpectedly commuted with an off-axis translation; the "
        "collinearity hypothesis of Thm. 1 would be vacuous"
    )


def test_paper58_l_selection_is_same_center_only():
    """Thm. 1(ii): the same (l, l') pair is zero on one center, nonzero on two.

    This is the decisive contrast.  <1s|2p_0> vanishes identically at a single
    center by orthogonality of the angular factors, for ANY radial functions.
    Across two centers on the axis, at the same m = 0, it is nonzero.
    """
    Z = Fraction(1)
    R = Fraction(1)

    # Same center: the theta factor alone kills it, independent of radial part.
    ct = sp.Symbol("ct")
    ang = sp.integrate(_G_lm(0, 0, ct) * _G_lm(1, 0, ct), (ct, -1, 1))
    assert sp.simplify(ang) == 0, (
        "same-center 1s/2p_0 angular overlap is not zero; the contrast this "
        "test rests on has evaporated"
    )

    # Two centers: same (l, l') = (0, 1), same m = 0 -> nonzero, exactly.
    U, V = _overlap_UV(Z, 1, 0, Z, 2, 1, 0, R)
    assert _nonzero(U, V), (
        "cross-center <1s_A|2p0_B> vanished; Thm. 1(ii) (no cross-center "
        "l-selection) would be false"
    )

    # Support spreads across ALL available l at fixed m: the re-expansion
    # does not terminate, which is the l-density the census measured.
    states = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
    dead = [s for s in states
            if not _nonzero(*_overlap_UV(Z, 1, 0, Z, s[0], s[1], 0, R))]
    assert dead == [], (
        f"cross-center overlap of 1s_A vanished against {dead}; l-support "
        "is expected to be full at fixed m"
    )
