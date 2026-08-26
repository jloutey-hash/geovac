"""Backing tests for Paper 14 (sec:tc_atomic_sparsity) -- the atomic xTC operator
INHERITS GeoVac's angular Gaunt sparsity, exactly for s-referenced systems and
near-exactly (small, structured, sub-% fill-in) for open-p references.

These are the tracked artifacts for the v5.0.3--v5.0.5 result (memos
``debug/sprint_xtc_pinclusive_memo.md``, ``debug/sprint_xtc_pblock_memo.md``).  Each
count is EXACT (angular support is radial-independent), so the assertions are exact
integers, not tolerances.

Every claim is mapped to an assertion that would FAIL if the claim were false, and the
suite carries its own non-tautology discriminator: the SAME engine that returns 0
fill-in for an s-reference returns 4 (s+p+d) and 96 (s+p+d+f) for an open-p reference,
so "0 fill-in" is a measured structural fact, not a hard-coded zero.

  (1) s-reference: 0 fill-in AND 0 fill-out at every basis -> support identical ->
      the Jordan-Wigner Pauli term set is preserved (the qubit sparsity survives xTC);
  (2) open-p reference: 0 fill-in through s+p (l<=1 protected), then a small structured
      fill-in once d orbitals enter (4 at s+p+d, 96 at s+p+d+f), fill-out always 0, and
      the density stays single-digit % and tracks Coulomb within <=2.5% relative;
  (3) mechanism: the reference-density multipole set is {(0,0)} for a spherical
      (s / closed-shell) reference and {(0,0),(2,0)} for an open-p reference -- the L'=2
      leg is exactly what makes fill-in possible;
  (4) m-conservation: the contracted operator conserves total m (no fill-in ever
      violates it), inherited from the vertex m-rule.
"""
import pytest

from geovac.xtc_angular_sparsity import (
    REF_ANG, fast_angular, four_Y, is_spherical, open_subshell,
    reference_density_multipoles,
)

_LABEL = {1: "s+p", 2: "s+p+d", 3: "s+p+d+f"}


# ===========================================================================
# (1) s-reference: 0 fill-in AND 0 fill-out at EVERY basis (support preserved)
# ===========================================================================
@pytest.mark.parametrize("lmax", [1, 2, 3])
def test_s_reference_zero_fillin_all_bases(lmax):
    """Be 1s^2 2s^2 (closed/spherical): xTC support == Coulomb support, exactly."""
    r = fast_angular(lmax, REF_ANG["Be_s2"])
    assert r["n_fill_in"] == 0, f"{_LABEL[lmax]}: s-ref fill-in {r['n_fill_in']} != 0"
    assert r["n_fill_out"] == 0
    # identical support => identical nonzero-block set => Pauli term set preserved
    assert r["n_l3"] == r["n_coulomb"]
    # the reference density is a pure monopole (this is *why* it is bit-exact)
    assert r["refmult"] == [(0, 0)]


# ===========================================================================
# (2) open-p reference: protected through s+p, small structured fill-in at d/f
# ===========================================================================
@pytest.mark.parametrize("name", ["C_2p2", "O_2p4"])
@pytest.mark.parametrize("lmax,expected_fill", [(1, 0), (2, 4), (3, 96)])
def test_pblock_reference_structured_fillin(name, lmax, expected_fill):
    """Open-p reference: 0 fill-in at s+p; 4 at s+p+d; 96 at s+p+d+f (fill-out 0)."""
    r = fast_angular(lmax, REF_ANG[name])
    assert r["n_fill_in"] == expected_fill, (
        f"{name} {_LABEL[lmax]}: fill-in {r['n_fill_in']} != {expected_fill}")
    assert r["n_fill_out"] == 0
    # the operator stays sparse: single-digit % density, tracking Coulomb closely
    assert r["l3_density"] < 0.15
    assert r["l3_density"] >= r["coul_density"]                 # fill-in can only add
    assert (r["l3_density"] - r["coul_density"]) / r["coul_density"] <= 0.025  # <=2.5% rel


def test_pblock_fillin_is_a_real_measurement_not_a_zero():
    """Non-tautology guard: the engine that returns 0 for s-ref DOES fill in for p-ref.

    Without this, "0 fill-in for the s-reference" could be a machinery artifact.  The
    s-ref and p-ref differ ONLY in the reference occupation fed to the same counter.
    """
    # exercised at BOTH lmax where the p-reference fills in (2 and 3), so the
    # s-reference 0 has a nonzero sibling at each basis -- an "if s_ref: return 0"
    # shortcut cannot pass this.
    for lmax, expect_p in [(2, 4), (3, 96)]:
        s_ref = fast_angular(lmax, REF_ANG["Be_s2"])["n_fill_in"]
        p_ref = fast_angular(lmax, REF_ANG["C_2p2"])["n_fill_in"]
        assert s_ref == 0 and p_ref == expect_p and p_ref > s_ref


# ===========================================================================
# (3) mechanism: the reference-density multipole set gates the fill-in
# ===========================================================================
def test_reference_density_multipole_sets():
    """Spherical reference -> {(0,0)} only; open-p reference -> {(0,0),(2,0)}."""
    s_mult = sorted(reference_density_multipoles(REF_ANG["Be_s2"], Lmax=4))
    p_mult = sorted(reference_density_multipoles(REF_ANG["C_2p2"], Lmax=4))
    assert s_mult == [(0, 0)]
    assert p_mult == [(0, 0), (2, 0)]
    # the density multipole is even-L only and M=0 (diagonal, real-density) in both
    assert all(Lp % 2 == 0 and Mp == 0 for (Lp, Mp) in p_mult)


# ===========================================================================
# (4) m-conservation of the contracted vertex
# ===========================================================================
def test_vertex_m_conservation():
    """The four-harmonic vertex is nonzero only when M+M' = m_a - m_c (total m)."""
    # a genuine m-conserving channel exists (p0 <-> p0 through the reference monopole)
    val, active = four_Y(1, 0, 0, 0, 0, 0, 1, 0)
    assert abs(val) > 1e-12 and active
    # an m-NON-conserving request returns exactly zero with no active channel
    val_bad, active_bad = four_Y(1, 1, 0, 0, 0, 0, 1, -1)
    assert val_bad == 0.0 and active_bad == []


def test_vertex_is_genuinely_multichannel():
    """The vertex is non-abelian: some (a,c) couple to MORE THAN ONE resultant Lambda.

    This is the property that makes the 3-body operator not collapse to 2-body (so xTC
    is needed at all); if every vertex had a single resultant the sparsity question
    would be trivial.
    """
    # d<->d external pair with a quadrupole correlator leg reaches several Lambda
    _, active = four_Y(2, 0, 2, 0, 2, 0, 2, 0)
    assert len(active) >= 2


# ===========================================================================
# (5) monoatomic sweep: the Unsoeld spherical law governs 0 fill-in across the
#     whole Z=1..56 library (closed / half-filled-high-spin / s-open -> exact 0)
# ===========================================================================
@pytest.mark.parametrize("l,k", [
    (0, 1), (0, 2),                       # s-open: always spherical
    (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),   # p^1..p^6
    (2, 1), (2, 2), (2, 5), (2, 8), (2, 10),          # d^1,d^2,d^5,d^8,d^10
])
def test_unsold_spherical_law(l, k):
    """xTC fill-in at s+p+d+f is ZERO iff the reference density is spherical (Unsoeld).

    General, exact test of the whole monoatomic result: no atom-specific hard-coding --
    the occupation is generated from (l,k) and the predicted 0/nonzero is from
    ``is_spherical`` alone.  Catches both directions (a spurious fill-in on a spherical
    shell, or a missed fill-in on an open one)."""
    fill = fast_angular(3, open_subshell(l, k), Lmax=6)["n_fill_in"]
    if is_spherical(l, k):
        assert fill == 0, f"(l={l},k={k}) spherical but fill-in {fill} != 0"
    else:
        assert fill > 0, f"(l={l},k={k}) non-spherical but fill-in {fill} == 0"


@pytest.mark.parametrize("l,k", [(1, 3), (2, 5)])
def test_half_filled_high_spin_is_spherical(l, k):
    """Half-filled high-spin p^3 (N,P) and d^5 (Cr,Mn): spherical density -> 0 fill-in
    at EVERY basis, and a pure-monopole reference -- the s-reference guarantee extends
    to half-filled shells."""
    # the reference density collapses to the L'=0 monopole (Unsoeld) -- the mechanism
    assert sorted(reference_density_multipoles(open_subshell(l, k), Lmax=6)) == [(0, 0)]
    for lmax in (1, 2, 3):
        r = fast_angular(lmax, open_subshell(l, k), Lmax=2 * lmax)
        assert r["n_fill_in"] == 0 and r["n_fill_out"] == 0 and r["refmult"] == [(0, 0)]


def test_open_d_reference_multipole_and_fillin():
    """Open-d reference (transition metals) carries the extra L'=4 density multipole,
    but through s+p+d+f its fill-in matches open-p exactly (0/4/96)."""
    d2 = open_subshell(2, 2)
    assert sorted(reference_density_multipoles(d2, Lmax=6)) == [(0, 0), (2, 0), (4, 0)]
    for lmax, expect in [(1, 0), (2, 4), (3, 96)]:
        assert fast_angular(lmax, d2, Lmax=2 * lmax)["n_fill_in"] == expect


def test_open_p_and_open_d_diverge_at_g():
    """The extra L'=4 multipole of an open-d reference first BITES once g orbitals enter:
    at s+p+d+f+g the open-p fill-in is 264 and the open-d fill-in is 268 (+4).  Below g
    they are identical -- so the divergence is a genuine, localized L'=4 effect, not a
    generic densification (both remain ~0.07% of the 390625 blocks)."""
    p_fill = fast_angular(4, open_subshell(1, 2), Lmax=8)["n_fill_in"]
    d_fill = fast_angular(4, open_subshell(2, 2), Lmax=8)["n_fill_in"]
    assert p_fill == 264
    assert d_fill == 268
    assert d_fill > p_fill  # the L'=4 multipole adds fill-in the L'=2 case cannot reach
