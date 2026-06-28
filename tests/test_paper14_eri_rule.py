"""Paper 14/20 ERI-rule characterization tests — the A/B dual-rule framing.

Pins the documented framing (criteria.md "Dual-rule ERI framing"; CF-1;
debug/sprint_group4_prework_memo.md):

- the **QC product** (atomic `lattice_index` + composed `composed_qubit`) realizes the
  **pair-diagonal rule A** (drops the m-swap ERIs) — this is what makes the reported
  sparsity numbers what they are;
- the **precision-physics** path (`casimir_ci._gaunt_ck`) realizes the **exact global-M_L
  rule B** (keeps the m-swap ERIs);
- switching the product to B re-prices N_Pauli by a **constant factor** — 2.51× (main-group)
  / 3.25× (d-block) — leaving the scaling unchanged.

These lock the A-state: if the product is ever switched to global-M_L (option B) or the rule
drifts, they fail and force the papers (which disclose A) back in sync with the code.
"""
from __future__ import annotations
import pytest


def test_product_realizes_pair_diagonal_composed():
    """composed_qubit._ck_coefficient (rule A) drops the m-swap; casimir_ci._gaunt_ck (B) keeps it."""
    import geovac.composed_qubit as cq
    from geovac.casimir_ci import _gaunt_ck
    # c^2(p_{+1}, p_{-1}): the m-swap coupling. Rule A (q=mc-ma) -> 3j vanishes; rule B keeps it.
    assert abs(cq._ck_coefficient(1, 1, 1, -1, 2)) < 1e-12, \
        "composed product must realize pair-diagonal A (drop the m-swap ERI)"
    assert abs(_gaunt_ck(1, 1, 1, -1, 2)) > 1e-12, \
        "global-M_L B must keep the m-swap ERI (the exact Coulomb rule)"


def test_product_realizes_pair_diagonal_atomic():
    """The atomic lattice_index ERI table (rule A) contains ZERO m-transfer entries."""
    from geovac.lattice_index import LatticeIndex
    li = LatticeIndex(n_electrons=2, max_n=2, nuclear_charge=2,
                      vee_method='slater_full', h1_method='hybrid')
    states = li.lattice.states  # (n, l, m)
    n_mswap = sum(1 for (a, b, c, d) in li._eri if states[a][2] != states[c][2])
    assert n_mswap == 0, \
        f"atomic product must be pair-diagonal A; found {n_mswap} m-transfer ERIs (would be rule B)"


@pytest.mark.slow
def test_global_rule_reprices_constant_factor():
    """Switching the product to exact global-M_L re-prices N_Pauli 2.51× (main) / 3.25× (d-block)."""
    import geovac.composed_qubit as cq
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.casimir_ci import _gaunt_ck
    from geovac.ecosystem_export import _rebuild_spec

    def npauli(name, use_global):
        orig = cq._ck_coefficient
        if use_global:
            cq._ck_coefficient = lambda la, ma, lc, mc, k: _gaunt_ck(la, ma, lc, mc, k)
        try:
            spec = _rebuild_spec(name, R=None, max_n=2, core_method='pk')
            res = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
            return int(res['N_pauli'])
        finally:
            cq._ck_coefficient = orig

    lih_ratio = npauli('LiH', True) / npauli('LiH', False)    # main-group s/p
    sch_ratio = npauli('ScH', True) / npauli('ScH', False)    # d-block
    assert abs(lih_ratio - 2.51) < 0.05, f"main-group re-pricing {lih_ratio:.3f} != 2.51"
    assert abs(sch_ratio - 3.25) < 0.05, f"d-block re-pricing {sch_ratio:.3f} != 3.25"
    # the d-block re-prices HARDER than main-group (the "d-block sparser" reversal)
    assert sch_ratio > lih_ratio, "d-block must re-price harder than main-group (the reversal)"
