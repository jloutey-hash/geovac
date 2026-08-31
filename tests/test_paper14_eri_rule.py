"""Paper 14/20 ERI-rule tests — the EXACT global-M_L rule, everywhere.

HISTORY.  Until 2026-08-29 this file locked a documented "dual-rule framing"
(CF-1): the QC product realized a *pair-diagonal rule A* (each electron
conserves its own m), the precision paths the *exact global-M_L rule B*
(m_a + m_b = m_c + m_d), and switching A->B was claimed to re-price N_Pauli by
a constant factor (2.51x main-group / 3.25x d-block) "leaving the scaling
unchanged".

Both halves of that framing fell:

1. **Rule A was a sign error, not a convention.**  It was produced by
   `q = mc - ma` in `_ck_coefficient`, which makes the Wigner-3j bottom row
   sum to 2(mc - ma) instead of 0, so the coefficient vanishes unless
   ma == mc — silently deleting 88 of 126 nonzero c^k at l<=2 (69.8%).
   Arbitrated against 4-D quadrature of the angular factor from its
   definition (Condon–Shortley 4, the shipped order 0); the corrected
   evaluators agree with each other and improve the He variational energy
   (0.349% -> 0.255%, still bounded).

2. **The "constant factor, scaling unchanged" claim was a single-point
   artifact.**  At max_n=2 the B/A ratio is uniformly 2.51x across
   main-group molecules, but along the BASIS axis it grows
   (1.00x / 2.51x / 5.40x at max_n = 1/2/3 on LiH), so the within-molecule
   Pauli exponent moves 2.54 -> 3.17.  What survives exactly is the
   cross-molecule linearity N_Pauli = c*Q at fixed max_n, with c = 27.90
   (was 11.10 under rule A).

PI direction 2026-08-29: exact rule everywhere; production now realizes B on
every path.  See debug/sprint_eri_evaluator_defects_memo.md.
"""
from __future__ import annotations

import pytest


def test_all_paths_realize_exact_global_ml():
    """Every production c^k implementation keeps the m-swap coupling.

    The witness is c^2(p_{+1}, p_{-1}) — the coupling the sign error zeroed.
    Its exact value is -sqrt(6)/5.
    """
    import numpy as np

    expected = -np.sqrt(6.0) / 5.0

    import geovac.composed_qubit as cq
    from geovac.casimir_ci import _gaunt_ck
    from geovac.lattice_index import LatticeIndex
    from geovac.sturmian_solver import _ck_coefficient as ck_sturm
    from geovac.tc_integrals import _ck_coefficient as ck_tc
    from geovac.nuclear.harmonic_shell import _ck_coefficient as ck_ho
    from geovac.nuclear.potential_sparsity import ck_coefficient as ck_ps

    li = LatticeIndex.__new__(LatticeIndex)
    impls = {
        "composed_qubit": cq._ck_coefficient,
        "casimir_ci": _gaunt_ck,
        "lattice_index": lambda *a: li._ck_coefficient(*a),
        "sturmian_solver": ck_sturm,
        "tc_integrals": ck_tc,
        "nuclear.harmonic_shell": ck_ho,
        "nuclear.potential_sparsity": ck_ps,
    }
    for name, fn in impls.items():
        val = fn(1, 1, 1, -1, 2)
        assert abs(val - expected) < 1e-12, (
            f"{name}: c^2(1,+1;1,-1) = {val} != -sqrt(6)/5 -- the m-swap "
            "coupling is being dropped again (wrong-sign q regression)"
        )


def test_atomic_table_carries_m_transfer_and_conserves_ml():
    """The atomic ERI table keeps the 42 m-transfer entries and every entry
    obeys the exact Coulomb selection rule."""
    from geovac.lattice_index import LatticeIndex

    li = LatticeIndex(n_electrons=2, max_n=2, nuclear_charge=2,
                      vee_method='slater_full', h1_method='hybrid')
    states = li.lattice.states  # (n, l, m)
    n_mswap = sum(1 for (a, b, c, d) in li._eri if states[a][2] != states[c][2])
    assert n_mswap == 42, (
        f"atomic table should carry the 42 m-transfer ERIs of the exact rule; "
        f"got {n_mswap}"
    )
    for (a, b, c, d) in li._eri:
        assert states[a][2] + states[b][2] == states[c][2] + states[d][2], (
            f"M_L-violating ERI ({a},{b},{c},{d}) in the atomic table"
        )


def test_composed_repricing_pins():
    """Pin the measured exact-rule re-pricing so neither a silent regression
    to rule A nor a drift of the corrected numbers can pass.

    LiH composed at max_n=2: 837 non-identity Pauli terms (rule A gave 333).
    Cross-molecule universality: N_Pauli = 27.900 x Q exactly, three rows of
    the periodic table (rule A gave 11.10).
    """
    import warnings

    from geovac.ecosystem_export import _rebuild_spec
    from geovac.composed_qubit import build_composed_hamiltonian

    for name, want_np in (("LiH", 837), ("NaH", 558), ("KH", 558)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spec = _rebuild_spec(name, R=None, max_n=2, core_method="pk")
            res = build_composed_hamiltonian(spec, pk_in_hamiltonian=True,
                                             verbose=False)
        Q, NP = int(res["Q"]), int(res["N_pauli"])
        assert NP == want_np, f"{name}: N_pauli {NP} != {want_np}"
        assert abs(NP / Q - 27.900) < 1e-6, (
            f"{name}: N/Q = {NP/Q} != 27.900 -- the exact-rule universality "
            "coefficient moved"
        )
        # discrimination: the retired rule-A count must NOT reappear
        assert NP != 333 and abs(NP / Q - 11.10) > 1.0, (
            f"{name}: rule-A counts have re-surfaced (wrong-sign q regression)"
        )


@pytest.mark.slow
def test_composed_exponent_no_longer_2p5():
    """The within-molecule Pauli exponent under the exact rule is ~3.17
    (LiH, max_n=1..3) -- NOT the retired 2.5, and still below the Gaussian
    3.9--4.3.  Guards both directions: a regression to the sign error would
    restore ~2.5; a further defect inflating terms would push past 3.9."""
    import warnings

    import numpy as np

    import geovac.composed_qubit as cq

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sw = cq.composed_lih_scaling_sweep(max_n_values=[1, 2, 3],
                                           verbose=False)
    Q = np.array([d["Q"] for d in sw["sweep_data"]], float)
    P = np.array([d["N_pauli"] for d in sw["sweep_data"]], float)
    alpha = float(np.polyfit(np.log(Q), np.log(P), 1)[0])
    assert 3.0 <= alpha <= 3.35, (
        f"exact-rule composed exponent {alpha:.3f} outside [3.0, 3.35] "
        "(measured 3.172; the retired rule-A fit gave 2.54)"
    )
