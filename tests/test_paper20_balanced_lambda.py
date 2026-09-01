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


# (name, spec factory, R or None, table N_pauli (incl. identity),
#  lambda_ni, QWC groups)
#
# GEOMETRY, corrected 2026-08-31.  R is None for every row now, which
# means 'use the molecule's own bond length' -- the builder resolves it
# from spec.R.  It previously meant 'use 3.015', LiH's bond length, for
# every molecule, because build_balanced_hamiltonian defaulted to that and
# ignored the spec.  Ten of these twelve rows therefore carried
# LiH-geometry values, and the paper columns were re-synced from them.
# Only NaH and HCl were right, because this file passed their real R
# explicitly -- and those two are exactly the cells a 2026-08-30 pass
# 'corrected' in the papers, away from this test.  Retired (LiH-geometry)
# values, which must never re-surface: 110.5 / 879.6 / 895.3 / 914.1 /
# 32.2 / 128.0 / 877.0 / 877.3, QWC 1282 / 1477 / 1551 / 1175.
# See debug/qa/balanced_lambda_geometry_finding.md.
#
# QWC is pinned here too, but it is NOT a structural quantity: the greedy
# grouper's count moves with geometry (unlike the Pauli count, which does
# not).  See tests/test_paper20_geometry_independence.py.
_BALANCED_TABLE = [
    ('NaH',  nah_spec,  None,    575,  19.6,   69),
    ('MgH2', mgh2_spec, None,   4861, 111.8,  903),
    ('HCl',  hcl_spec,  None,   9824, 866.1, 1173),
    ('H2S',  h2s_spec,  None,  13863, 873.6, 1264),
    ('PH3',  ph3_spec,  None,  18854, 889.3, 1427),
    ('SiH4', sih4_spec, None,  24745, 909.2, 1566),
    ('KH',   kh_spec,   None,    575,  28.1,   69),
    ('CaH2', cah2_spec, None,   4861, 123.8,  898),
    ('HBr',  hbr_spec,  None,   9824, 875.3, 1176),
    ('H2Se', h2se_spec, None,  13863, 883.4, 1264),
    ('AsH3', ash3_spec, None,  18854, 885.2, 1468),
    ('GeH4', geh4_spec, None,  24745, 874.8, 1562),
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



# The two 575-term rows build in well under a second, so they run in the
# DEFAULT suite.  The 2026-08-31 geometry defect survived a month partly
# because every cell of this file sat behind --slow: 12 of its 13 cases
# were skipped on a normal run, so the guard that held the correct values
# never spoke up while the papers drifted away from it.
_FAST = {'NaH', 'KH'}


@pytest.mark.parametrize(
    "name,factory,R,n_pauli,lam,qwc",
    [pytest.param(*row, marks=() if row[0] in _FAST
                  else pytest.mark.slow) for row in _BALANCED_TABLE],
    ids=[r[0] for r in _BALANCED_TABLE],
)
def test_balanced_table_cell(name, factory, R, n_pauli, lam, qwc):
    """Pin the tab:molecules row at the molecule's OWN geometry:
    Pauli count (incl. identity) exact, lambda_ni within 0.5, QWC exact.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        spec = factory(max_n=2)
        res = (build_balanced_hamiltonian(spec, R=R) if R
               else build_balanced_hamiltonian(spec))
    assert _n_pauli_nonid(res) + 1 == n_pauli, name
    got = _lambda_ni(res)
    assert abs(got - lam) < 0.5, f"{name}: lambda_ni {got:.1f} vs table {lam}"
    assert res['n_qwc_balanced'] == qwc, (
        f"{name}: QWC {res['n_qwc_balanced']} vs table {qwc}")
