"""Paper 20 tab:molecules balanced lambda_ni column pins (7th cert, 2026-07-02).

The 7th group4 cert's completeness-critic asked for the never-verified
balanced-table cells; a live sweep found the entire second/third-row
lambda_ni column stale at drafting vintage (0.3-7% drift; Pauli COUNTS all
exact) -- the same stale-table-not-regression class as the v4.56.0 first-row
rel-lambda finding, which produced tests/test_paper14_rel_lambda.py. This is
the balanced-table analog: pin every lambda_ni cell so the column can never
drift unguarded again.

lambda_ni = sum |coeff| over non-identity Pauli terms of the balanced build
(max_n=2, spec-factory defaults; NaH/HCl at their table R where the factory
does not carry it).
"""
from __future__ import annotations

import warnings

import pytest

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import (
    nah_spec, mgh2_spec, hcl_spec, h2s_spec, ph3_spec, sih4_spec,
    kh_spec, cah2_spec, hbr_spec, h2se_spec, ash3_spec, geh4_spec,
)


def _lambda_ni(res: dict) -> float:
    op = res.get('qubit_op') or res.get('hamiltonian')
    return sum(abs(complex(c)) for t, c in op.terms.items() if t)


def _n_pauli_nonid(res: dict) -> int:
    op = res.get('qubit_op') or res.get('hamiltonian')
    return sum(1 for t in op.terms if t)


# (name, spec factory, R or None, table N_pauli (incl. identity), lambda_ni)
# lambda values recomputed live 2026-07-02 (paper cells synced same day);
# ALL cells re-measured 2026-08-29 under the exact global-M_L rule (the
# wrong-sign-q fix; see debug/sprint_eri_evaluator_defects_memo.md).
# Isostructural invariance survives exactly (NaH == KH, HCl == HBr, ...).
_BALANCED_TABLE = [
    ('NaH',  nah_spec,  3.566,   575,  19.6),
    ('MgH2', mgh2_spec, None,   4861, 110.5),
    ('HCl',  hcl_spec,  2.409,  9824, 866.1),
    ('H2S',  h2s_spec,  None,  13863, 879.6),
    ('PH3',  ph3_spec,  None,  18854, 895.3),
    ('SiH4', sih4_spec, None,  24745, 914.1),
    ('KH',   kh_spec,   None,    575,  32.2),  # factory R=4.243
    ('CaH2', cah2_spec, None,   4861, 128.0),
    ('HBr',  hbr_spec,  None,   9824, 877.0),
    ('H2Se', h2se_spec, None,  13863, 883.4),
    ('AsH3', ash3_spec, None,  18854, 885.2),
    ('GeH4', geh4_spec, None,  24745, 877.3),
]

# Retired rule-A cells (must never re-surface as live values): NaH/KH 239,
# MgH2/CaH2 1501, HCl/HBr 2936, H2S/H2Se 4119, PH3/AsH3 5582, SiH4/GeH4 7273.


def test_rule_a_counts_do_not_resurface():
    """A wrong-sign-q regression would restore the retired counts exactly;
    catch it on the cheapest cell (NaH)."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = build_balanced_hamiltonian(nah_spec(max_n=2), R=3.566)
    n = _n_pauli_nonid(res) + 1
    assert n != 239, "retired pair-diagonal NaH count has re-surfaced"
    assert n == 575, f"NaH balanced count drifted: {n}"



@pytest.mark.slow
@pytest.mark.parametrize(
    "name,factory,R,n_pauli,lam", _BALANCED_TABLE,
    ids=[r[0] for r in _BALANCED_TABLE],
)
def test_balanced_table_cell(name, factory, R, n_pauli, lam):
    """Pin the tab:molecules row: Pauli count (incl. identity) exact,
    lambda_ni within 0.5."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        spec = factory(max_n=2)
        res = (build_balanced_hamiltonian(spec, R=R) if R
               else build_balanced_hamiltonian(spec))
    assert _n_pauli_nonid(res) + 1 == n_pauli, name
    got = _lambda_ni(res)
    assert abs(got - lam) < 0.5, f"{name}: lambda_ni {got:.1f} vs table {lam}"
