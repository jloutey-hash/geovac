"""Paper 14 / Paper 20 — relativistic non-identity 1-norm (lambda_ni) regression pins.

These values back the lambda_ni^rel column of tab:spinor_resource (Paper 14) and
tab:paper20_tier2 (Paper 20).  They had NO test until 2026-06-29 and silently drifted:
the published first-row (LiH/BeH) cells were stale (pre-chemistry-re-entry) values that
contradicted the live code by up to 3.6x (BeH n=2 table 40.26 vs code 143.96), while the
frozen-core CaH path matched exactly.  A QA cert surfaced the drift; diagnosis (CHANGELOG
v4.56.0) showed the code is authoritative (the table's BeH rel lambda 40.26 was physically
impossible: 3.5x BELOW its own scalar lambda 139.12).  This test pins the live values so
the column can never drift unguarded again.

lambda_ni = sum of |coeff| over the non-identity Pauli terms of the qubit operator
(matches the paper's lambda_ni^rel definition).
"""
import warnings

import pytest

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import (
    lih_spec_relativistic,
    beh_spec_relativistic,
    cah_spec_relativistic,
)


def _lambda_ni(res: dict) -> float:
    op = res["qubit_op"] if "qubit_op" in res else res["hamiltonian"]
    return sum(abs(complex(c)) for term, c in op.terms.items() if term)


# (name, spec factory, expected lambda_ni^rel at n_max=2)  -- the default (include_breit=False)
# relativistic build, i.e. the encoding the papers describe.
_REL_LAMBDA_NMAX2 = [
    ("LiH", lih_spec_relativistic, 40.59),
    ("BeH", beh_spec_relativistic, 143.96),
    ("CaH", cah_spec_relativistic, 18.68),
]


@pytest.mark.slow
@pytest.mark.parametrize("name,factory,expected", _REL_LAMBDA_NMAX2)
def test_rel_lambda_ni_nmax2(name: str, factory, expected: float) -> None:
    """Pin lambda_ni^rel at n_max=2 (the native-Q sweet spot reported in the papers)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = build_composed_hamiltonian(factory(max_n=2))
    got = _lambda_ni(res)
    assert abs(got - expected) < 0.05, (
        f"{name} n_max=2 rel lambda_ni = {got:.2f}, expected {expected:.2f} "
        f"(tab:spinor_resource / tab:paper20_tier2 cell)"
    )


@pytest.mark.slow
def test_rel_lambda_ni_exceeds_pauli_count_sanity() -> None:
    """Physical sanity that caught the stale table: the relativistic 1-norm must NOT fall
    far below the scalar 1-norm (the stale BeH cell claimed rel lambda 3.5x BELOW scalar)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rel = _lambda_ni(build_composed_hamiltonian(beh_spec_relativistic(max_n=2)))
        sc_spec = beh_spec_relativistic(max_n=2)
        sc_spec.relativistic = False
        sca = _lambda_ni(build_composed_hamiltonian(sc_spec))
    # rel lambda is a modest increase over scalar (~+3.5% here), never a 3.5x DECREASE.
    assert rel >= 0.9 * sca, (
        f"BeH rel lambda_ni {rel:.2f} fell below 0.9x scalar {sca:.2f} "
        f"-- the stale-table failure mode (40.26 vs scalar 139.12)"
    )
