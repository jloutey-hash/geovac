"""Backing test for Paper 19's sharpened R_eq-drift localization.

Paper 19 (Sec. "Fixed-geometry energy versus well shape") now states that the
balanced-coupled LiH well is a STRICT TWO-TERM balance: only the nuclear
repulsion V_NN and the cross-center V_ne carry any R-dependence, and the 8.8%
outward R_eq drift is specifically the (too-weak) R-slope of the cross-center
V_ne. Every OTHER electronic term -- within-block one-body, within-block ERI, and
cross-block ERI -- is bit-exactly R-independent.

This test backs the structural core of that claim (the cheap, decisive half): the
R-dependence partition. The force-magnitude half (F[V_NN]=-0.331, F[cross]=+0.302,
sum=-0.029 = tilt; 8.8% deficit) is an FCI-PES measurement chronicled in
debug/sprint_balanced_coupled_reqdrift_memo.md + debug/data/balanced_reqdrift_termdecomp.json.

It is NOT vacuous by construction: it asserts both that the three within/cross-block
terms do NOT change with R AND that V_ne + V_NN DO change. A build that returned
identical or all-zero matrices would fail the second set; a build that let a
within-block term drift with R would fail the first (the wrong answer the
sharpened claim excludes).
"""
import numpy as np
import pytest

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec

MAX_N = 2
R1, R2 = 3.015, 3.30          # true min and an outward point
N_GRID = 2000                 # R-independence is exact in R regardless of grid
L_MAX = 4


def _components(R):
    spec = lih_spec(R=R, max_n=MAX_N)
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=N_GRID, L_max=L_MAX)
    comp = build_composed_hamiltonian(spec, pk_in_hamiltonian=False, verbose=False)
    eri_within = comp['eri']
    return {
        'h1_no_pk': ham['h1_no_pk'],
        'h1_cross_vne': ham['h1_cross_vne'],
        'eri_balanced': ham['eri'],
        'eri_within': eri_within,
        'cross_block_eri': ham['eri'] - eri_within,
        'vnn': float(ham['nuclear_repulsion']),
    }


@pytest.mark.slow
def test_balanced_lih_well_is_two_term_vnn_crossvne_balance():
    c1, c2 = _components(R1), _components(R2)

    # --- the R-INDEPENDENT terms (must be bit-exactly identical across R) ------
    R_INDEP_TOL = 1e-10
    for key in ('h1_no_pk', 'eri_within', 'cross_block_eri'):
        d = float(np.max(np.abs(c1[key] - c2[key])))
        assert d < R_INDEP_TOL, (
            f"{key} is supposed to be bit-exactly R-independent (the paper's "
            f"two-term-balance claim), but max|Delta(R)| = {d:.2e} between "
            f"R={R1} and R={R2}"
        )

    # --- the R-DEPENDENT terms (must actually move; guards against vacuity) ----
    d_vne = float(np.max(np.abs(c1['h1_cross_vne'] - c2['h1_cross_vne'])))
    d_vnn = abs(c1['vnn'] - c2['vnn'])
    assert d_vne > 1e-3, (
        f"cross-center V_ne did not change with R (max|Delta| = {d_vne:.2e}); "
        f"if it is R-independent the whole two-term-balance story is wrong"
    )
    assert d_vnn > 1e-3, (
        f"V_NN did not change with R (|Delta| = {d_vnn:.2e}); Z_A Z_B/R must move"
    )
    # V_NN change must be the analytic 3/R difference (Z_Li Z_H = 3)
    expected_dvnn = 3.0 * abs(1.0 / R1 - 1.0 / R2)
    assert abs(d_vnn - expected_dvnn) < 1e-9, (
        f"V_NN R-change {d_vnn:.6f} != analytic 3(1/R1-1/R2) {expected_dvnn:.6f}"
    )


if __name__ == '__main__':
    import sys
    test_balanced_lih_well_is_two_term_vnn_crossvne_balance()
    print("PASS")
    sys.exit(0)
