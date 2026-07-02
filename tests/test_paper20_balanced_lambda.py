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
# lambda values recomputed live 2026-07-02 (paper cells synced same day).
_BALANCED_TABLE = [
    ('NaH',  nah_spec,  3.566, 239,  18.9),
    ('MgH2', mgh2_spec, None, 1501, 101.9),
    ('HCl',  hcl_spec,  2.409, 2936, 798.8),
    ('H2S',  h2s_spec,  None, 4119, 816.9),
    ('PH3',  ph3_spec,  None, 5582, 834.9),
    ('SiH4', sih4_spec, None, 7273, 853.6),
    ('KH',   kh_spec,   None,  239,  31.6),   # factory R=4.243; the v4.56-era
    # 28.15 probe value was unreproducible at any candidate geometry (superseded)
    ('CaH2', cah2_spec, None, 1501, 118.9),
    ('HBr',  hbr_spec,  None, 2936, 809.4),
    ('H2Se', h2se_spec, None, 4119, 820.2),
    ('AsH3', ash3_spec, None, 5582, 824.2),
    ('GeH4', geh4_spec, None, 7273, 815.8),
]


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
