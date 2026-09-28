"""Library He qubit Hamiltonian: physical h1 vs hop-augmented ('hybrid') h1.

Re-measures Pauli count (non-identity), identity-inclusive count, non-identity
1-norm lambda, and greedy QWC groups for He at a given n_max, for:
  mode='exact'    LatticeIndex(..., h1_method='exact')    -- the PHYSICAL h1
  mode='hybrid'   LatticeIndex(..., h1_method='hybrid')   -- legacy hop model
  mode='composed' build_composed_hamiltonian(He atomic-block spec) -- an
                  independent second route (composed_qubit ERI + h1 code)
  mode='library'  geovac.ecosystem_export.hamiltonian('He', max_n=n)

Usage:  python debug/library_he_fix.py <n_max> <mode> [--qwc]
Appends one JSON line to debug/data/library_he_fix.jsonl.
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

OUT = os.path.join(os.path.dirname(__file__), "data", "library_he_fix.jsonl")


def build(n: int, mode: str):
    warnings.simplefilter("ignore")
    if mode in ("exact", "hybrid"):
        from geovac.lattice_index import LatticeIndex
        from geovac.qubit_encoding import JordanWignerEncoder
        li = LatticeIndex(n_electrons=2, max_n=n, nuclear_charge=2,
                          vee_method="slater_full", h1_method=mode)
        H1 = li._H1_spatial.toarray()
        import numpy as np
        offdiag = float(np.abs(H1 - np.diag(np.diag(H1))).max())
        return JordanWignerEncoder(li).build_qubit_operator(), offdiag
    if mode == "composed":
        from geovac.ecosystem_export import _rebuild_spec
        from geovac.composed_qubit import build_composed_hamiltonian
        import numpy as np
        spec = _rebuild_spec("He", R=None, max_n=n, core_method="pk")
        res = build_composed_hamiltonian(spec, verbose=False)
        h1 = np.asarray(res["h1"])
        offdiag = float(np.abs(h1 - np.diag(np.diag(h1))).max())
        return res["qubit_op"], offdiag
    if mode == "library":
        from geovac.ecosystem_export import hamiltonian
        return hamiltonian("He", max_n=n)._qubit_op, None
    raise ValueError(mode)


def main() -> None:
    n = int(sys.argv[1])
    mode = sys.argv[2]
    do_qwc = "--qwc" in sys.argv
    t0 = time.time()
    qop, offdiag = build(n, mode)
    t_build = time.time() - t0
    terms = qop.terms
    n_nonid = sum(1 for k in terms if k)
    lam = float(sum(abs(c) for k, c in terms.items() if k))
    ident = complex(terms.get((), 0.0)).real
    nq = 1 + max((q for k in terms for q, _ in k), default=-1)
    row = dict(n_max=n, mode=mode, Q=nq, n_pauli_nonid=n_nonid,
               n_pauli_incl_id=len(terms), lambda_nonid=lam,
               identity_coeff=ident, h1_max_offdiag=offdiag,
               t_build_s=round(t_build, 1))
    if do_qwc:
        from geovac.measurement_grouping import count_qwc_groups
        t1 = time.time()
        row["qwc_groups"] = count_qwc_groups(qop)
        row["t_qwc_s"] = round(time.time() - t1, 1)
    print(json.dumps(row), flush=True)
    with open(OUT, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
