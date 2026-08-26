"""Backing tests for Paper 14 (sec:tc) -- the atomic transcorrelated (TC/xTC)
non-Hermitian effective operator and its quantum-algorithm cost.

Each load-bearing claim of the v5.0.6 sprint (memos
``debug/sprint_tc_quantum_cost_measure_memo.md`` + its accuracy-validation section)
is mapped to a test that GENUINELY constrains it -- an assertion that would FAIL if
the claim were false, not a tautology:

  (1) the xTC 3-body -> 2-body contraction ``v2`` is Hermitian to machine precision
      (falsifier: the sibling convective ``K`` term, built by the same machinery, has
      deviation ~0.2, so a <1e-14 bound is a real discriminator);
  (2) the ONLY non-Hermitian term is the convective ``K``: ``D`` and ``xTC-L3`` are
      Hermitian, ``K`` deviation ~0.2, and deleting ``K`` restores a Hermitian operator;
  (3) the full non-Hermitian operator (He, Li) has a real spectrum, a real ground
      state, benign eigenvector conditioning ``kappa_V in [1.0, 2.2]``, and an LCU
      1-norm ``<= 1.13x`` plain;
  (4) ``K`` roughly DOUBLES the Jordan-Wigner Pauli count (117 -> 249); symmetrizing
      restores 117 with 1-norm ~0.98x plain;
  (5) the symmetrization "escape hatch" is an ACCURACY false economy: at the best
      tested basis ``E_sym`` overshoots BELOW the exact energy while ``E_TC`` stays
      ABOVE it -- the two land on OPPOSITE sides of exact.

The engine is the tracked port :mod:`geovac.transcorrelated_sturmian`, validated
bit-for-bit against the historical ``debug/xtc_poc_li.py`` driver and the pinned
``debug/data/xtc_quantum_cost.json`` / ``xtc_sym_accuracy.json`` numbers.
"""
import numpy as np
import pytest
import scipy.linalg as sla

from geovac import transcorrelated_sturmian as TC

# Cost metrics are grid-converged (kappa_V, lambda-ratio, Pauli count are invariant to
# ~1e-4 from Ng=400..800); use a light grid for the default-collected fixtures.
_NG, _NX = 500, 96
# Accuracy runs need the k-optimized grid of the accuracy-validation driver.
_NG_ACC, _NX_ACC = 600, 96


# ---------------------------------------------------------------------------
# module-scoped systems (built once)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def he_cost():
    """He 1s^2, s-only ns=3, k=1.7: 2-body TC ``D + K`` (no genuine 3-body)."""
    return TC.build_atomic_system(3, 1.7, 1.0, Z=2, n_elec=2,
                                  Ng=_NG, nx=_NX, with_L3=False)


@pytest.fixture(scope="module")
def li_cost():
    """Li 1s^2 2s, s-only ns=3, k=1.5: genuine 3-body -> ``D + K + xTC-L3``."""
    return TC.build_atomic_system(3, 1.5, 1.0, Z=3, n_elec=3,
                                  Ng=_NG, nx=_NX, with_L3=True)


def _full_variant(sys):
    """The full non-Hermitian TC operator variant for this system."""
    return sys.xtc_full() if sys.v2 is not None else sys.tc2()


def _ground_evector_max_imag(H):
    """max|Im| of the phase-fixed ground right-eigenvector (0 == real ground state)."""
    w, VR = sla.eig(H, right=True)
    i0 = int(np.argmin(w.real))
    v = VR[:, i0]
    v = v / v[int(np.argmax(np.abs(v)))]
    return float(np.max(np.abs(v.imag)))


# ===========================================================================
# (1) xTC 3-body -> 2-body contraction is Hermitian to machine precision
# ===========================================================================
def test_xtc_contraction_is_hermitian(li_cost):
    """The xTC contraction v2 is Hermitian for the Li s-reference (<1e-14).

    Discriminator: the convective K term, produced by the same TC machinery, has a
    Hermiticity deviation ~0.2 (see test_only_K_is_nonhermitian), so a <1e-14 bound
    genuinely rules out non-Hermitian fill-in from the contraction itself.
    """
    dev_v2 = TC.herm_dev(li_cost.v2)
    dev_K = TC.herm_dev(li_cost.asym_K)
    assert dev_v2 < 1e-14, f"xTC contraction not Hermitian: {dev_v2:.2e}"
    # sanity that the bound is meaningful: the sibling K term is nowhere near it
    assert dev_K > 1e6 * dev_v2


# ===========================================================================
# (2) the ONLY non-Hermitian term is the convective K
# ===========================================================================
def test_only_K_is_nonhermitian(li_cost):
    """D and xTC-L3 are Hermitian; K is non-Hermitian (~0.2); deleting K restores it."""
    assert TC.herm_dev(li_cost.asym_w) < 1e-12          # D (finite TC kernel)
    assert TC.herm_dev(li_cost.v2) < 1e-14              # xTC-contracted L3
    assert 0.15 < TC.herm_dev(li_cost.asym_K) < 0.35    # convective K (~0.226)

    # the FULL operator is non-Hermitian, but D + xTC-L3 (K removed) is Hermitian
    H_full = TC.build_fci_matrix(li_cost, li_cost.xtc_full())
    H_noK = TC.build_fci_matrix(li_cost, li_cost.xtc_noK())
    assert TC.herm_dev_matrix(H_full) > 1e-2            # K makes it non-Hermitian
    assert TC.herm_dev_matrix(H_noK) < 1e-10            # remove K -> Hermitian again


def test_only_K_is_nonhermitian_he(he_cost):
    """He 2-body control: D Hermitian, K non-Hermitian (~0.29)."""
    assert TC.herm_dev(he_cost.asym_w) < 1e-12
    assert 0.15 < TC.herm_dev(he_cost.asym_K) < 0.35    # ~0.287


# ===========================================================================
# (3) full operator: real spectrum, real ground state, benign kappa_V, mild lambda
# ===========================================================================
# kV_ref pinned ONLY where grid-stable.  Li's low-lying block conditioning is stable
# (1.208 at Ng=400..800); He's bounces (1.28-1.98 across grids) because its full-CI
# low-lying block sits near a degeneracy -- so for He only the robust O(1) range is
# asserted, not a tight value.  lam-ratio and non-normality ARE grid-stable and pinned.
@pytest.mark.parametrize("sysname,kV_tol,lamr_ref,nonnorm_ref", [
    ("he", None,  1.096, 0.0694),   # He D+K   (kV grid-sensitive -> range only)
    ("li", 0.02,  1.122, 0.0147),   # Li D+K+xTC-L3 (kV stable)
])
def test_full_operator_real_and_conditioned(sysname, kV_tol, lamr_ref, nonnorm_ref,
                                            he_cost, li_cost):
    sys = he_cost if sysname == "he" else li_cost
    hso, asym, v0 = _full_variant(sys)
    H = TC.build_H(sys.dets, sys.didx, hso, asym, sys.nso, v0=v0)
    mm = TC.measure_matrix(H)

    # real spectrum + real ground state
    assert mm["all_real"] is True
    assert mm["im_spread"] < 1e-9
    assert mm["im0"] < 1e-9
    assert _ground_evector_max_imag(H) < 1e-8

    # eigenvector conditioning is benign O(1); 2.2 is the STOP guard (robust for both)
    assert 1.0 <= mm["kV_low"] <= 2.2
    if kV_tol is not None:
        assert abs(mm["kV_low"] - 1.208) < kV_tol

    # non-normality (grid-stable): small but nonzero for the genuine non-Hermitian op
    assert abs(mm["nonnormality"] - nonnorm_ref) < 1e-3
    assert mm["nonnormality"] > 1e-3          # genuinely non-normal (not Hermitian)

    # LCU 1-norm inflation stays well under the 1.5x STOP (measured <= 1.13x)
    lam_full = TC.lcu_lambda(hso, asym, sys.nso)["lam"]
    lam_plain = TC.lcu_lambda(sys.hso, sys.asym_coul, sys.nso)["lam"]
    ratio = lam_full / lam_plain
    assert ratio <= 1.13
    assert abs(ratio - lamr_ref) < 0.01


def test_li_noK_is_well_conditioned_and_sparse(li_cost):
    """The Hermitian D+xTC-L3 (no K) is kappa_V=1 and lowers lambda to ~0.85x plain."""
    hso, asym, v0 = li_cost.xtc_noK()
    H = TC.build_H(li_cost.dets, li_cost.didx, hso, asym, li_cost.nso, v0=v0)
    mm = TC.measure_matrix(H)
    assert mm["nonnormality"] < 1e-12          # Hermitian
    assert abs(mm["kV_low"] - 1.0) < 1e-6
    lam = TC.lcu_lambda(hso, asym, li_cost.nso)
    lam_plain = TC.lcu_lambda(li_cost.hso, li_cost.asym_coul, li_cost.nso)["lam"]
    assert lam["n_pauli"] == 117               # no K -> no Pauli doubling
    assert abs(lam["lam"] / lam_plain - 0.852) < 0.01


# ===========================================================================
# (4) K doubles the Pauli count; symmetrization restores it (lambda ~0.98x plain)
# ===========================================================================
@pytest.mark.parametrize("sysname", ["he", "li"])
def test_K_doubles_pauli_symmetrization_restores(sysname, he_cost, li_cost):
    sys = he_cost if sysname == "he" else li_cost

    n_plain = TC.lcu_lambda(sys.hso, sys.asym_coul, sys.nso)["n_pauli"]
    lam_plain = TC.lcu_lambda(sys.hso, sys.asym_coul, sys.nso)["lam"]
    assert n_plain == 117

    hso, asym, v0 = _full_variant(sys)
    full = TC.lcu_lambda(hso, asym, sys.nso)
    assert full["n_pauli"] == 249                        # K doubles the count
    # K contributes genuine imaginary Pauli coefficients (non-Hermitian signature)
    assert full["max_imag_coeff"] > 1e-3

    asym_sym = 0.5 * (asym + TC.dagger_asym(asym))
    symm = TC.lcu_lambda(hso, asym_sym, sys.nso)
    assert symm["n_pauli"] == 117                        # symmetrization restores plain
    assert symm["max_imag_coeff"] < 1e-12               # Hermitian again
    ratio_sym = symm["lam"] / lam_plain
    assert 0.90 <= ratio_sym <= 1.0                      # ~0.977 / 0.978
    assert abs(ratio_sym - 0.977) < 0.01


# ===========================================================================
# (5) symmetrization is an ACCURACY false economy (opposite sides of exact)
# ===========================================================================
def test_symmetrization_false_economy_he():
    """He (best-basis fast): E_TC is ABOVE exact, E_sym overshoots BELOW exact."""
    sys = TC.build_atomic_system(3, 2.0, 1.0, Z=2, n_elec=2,
                                 Ng=_NG_ACC, nx=_NX_ACC, with_L3=False)
    E_TC, im, E_sym = TC.tc_energies(sys)
    exact = TC.EXACT_NR["He"]
    assert im < 1e-9
    assert (E_TC - exact) > 0.0, f"E_TC not above exact: {E_TC - exact:+.5f}"
    assert (E_sym - exact) < 0.0, f"E_sym not below exact: {E_sym - exact:+.5f}"


@pytest.mark.slow
def test_symmetrization_false_economy_li_basis_dependence():
    """Li: E_sym only overshoots below exact at the BEST basis (ns=5); at ns=3 both
    E_TC and E_sym are above exact.  The sign flip with basis is the load-bearing
    CAUTION finding (symmetrization degrades as the basis improves)."""
    exact = TC.EXACT_NR["Li"]

    # small basis: both above exact (no overshoot yet)
    s3 = TC.build_atomic_system(3, 1.6, 1.0, Z=3, n_elec=3,
                                Ng=_NG_ACC, nx=_NX_ACC, with_L3=True)
    E_TC3, _, E_sym3 = TC.tc_energies(s3)
    assert (E_TC3 - exact) > 0.0
    assert (E_sym3 - exact) > 0.0                     # ns=3: sym still above exact

    # best basis: E_TC above, E_sym overshoots below (opposite sides of exact)
    s5 = TC.build_atomic_system(5, 1.6, 1.0, Z=3, n_elec=3,
                                Ng=_NG_ACC, nx=_NX_ACC, with_L3=True)
    E_TC5, im5, E_sym5 = TC.tc_energies(s5)
    assert im5 < 1e-9
    assert (E_TC5 - exact) > 0.0, f"E_TC(ns5) not above exact: {E_TC5 - exact:+.5f}"
    assert (E_sym5 - exact) < 0.0, f"E_sym(ns5) not below exact: {E_sym5 - exact:+.5f}"


# ===========================================================================
# (1', secondary/slow) open-p reference: contraction still Hermitian
# ===========================================================================
@pytest.mark.slow
def test_xtc_contraction_hermitian_open_p():
    """Open-p reference (C 3P) generalizes the s-reference Hermiticity result.

    Reuses the debug p-block engine READ-ONLY; skipped if unavailable.
    """
    import os
    import sys as _sys
    dbg = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "debug")
    if dbg not in _sys.path:
        _sys.path.insert(0, dbg)
    try:
        import xtc_poc_li_pinclusive as P  # noqa: F401
        from xtc_pblock_engine import REFERENCES
    except Exception as exc:  # pragma: no cover - environment-dependent
        pytest.skip(f"p-block engine unavailable: {exc}")

    spec = REFERENCES["C_3P"]
    o = P.assemble(max_n=2, k=2.0, gamma=1.0, Z=spec["Z"], n_elec=spec["n_elec"],
                   Ng=600, nx=120, want_exact3=True, ref_occ=spec["ref_occ"])
    v2 = o["_arr"]["v2"]
    assert TC.herm_dev(v2) < 1e-14
