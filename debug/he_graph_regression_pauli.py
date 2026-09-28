"""He qubit Hamiltonian (ecosystem 'He' = LatticeIndex slater_full, h1_method='hybrid')
vs the physical h1_method='exact': Pauli count and 1-norm (report-only for the PM;
Paper 14 / group4 registry keys he_n*_pauli / he_n*_lambda were measured on 'hybrid')."""
import os
import sys
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
warnings.filterwarnings("ignore")

from geovac.lattice_index import LatticeIndex  # noqa: E402
from geovac.qubit_encoding import JordanWignerEncoder  # noqa: E402


def stats(op):
    terms = {k: v for k, v in op.terms.items() if k != ()}
    return len(terms), sum(abs(v) for v in terms.values())


for n in [int(a) for a in sys.argv[1:]] or [2, 3]:
    for mode in ("hybrid", "exact"):
        li = LatticeIndex(n_electrons=2, max_n=n, nuclear_charge=2,
                          vee_method="slater_full", h1_method=mode)
        op = JordanWignerEncoder(li).build_qubit_operator()
        npauli, lam = stats(op)
        print(f"He n_max={n} h1={mode:6s}: Q={li.n_sp} Pauli(non-id)={npauli} lambda={lam:.3f}",
              flush=True)
