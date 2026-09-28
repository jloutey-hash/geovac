"""Paper 14 H-atom 1-norm control: physical h1 vs retired hop model.

Physical h1 (diagonal -1/(2n^2)) gives non-identity lambda = n_max/2 exactly
(each of the 2n^2 spin-orbitals of shell n contributes |eps_n|/2 via
n_p = (1 - Z_p)/2), i.e. lambda ~ Q^{1/2}.  The retired 1.57 exponent is
checked here against the hybrid (hop) model.
"""
import os
import sys
import warnings

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
warnings.simplefilter("ignore")

from geovac.lattice_index import LatticeIndex  # noqa: E402
from geovac.qubit_encoding import JordanWignerEncoder  # noqa: E402

for mode in ("exact", "hybrid"):
    Qs, lams = [], []
    for n in (2, 3):
        li = LatticeIndex(n_electrons=1, max_n=n, nuclear_charge=1,
                          vee_method="slater_full", h1_method=mode)
        qop = JordanWignerEncoder(li).build_qubit_operator()
        lam = float(sum(abs(c) for k, c in qop.terms.items() if k))
        lam_all = float(sum(abs(c) for c in qop.terms.values()))
        Qs.append(2 * n * n)
        lams.append(lam)
        print(mode, n, 2 * n * n, "lambda_nonid", round(lam, 6),
              "lambda_incl_id", round(lam_all, 6), "n/2", n / 2)
    a = np.polyfit(np.log(Qs), np.log(lams), 1)[0]
    print(mode, "exponent (4 pts, non-identity)", round(a, 4))
