"""He FCI energy of the library (physical-h1) Hamiltonian at n_max=2,3."""
import os
import sys
import warnings
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
warnings.simplefilter("ignore")
from geovac.lattice_index import LatticeIndex  # noqa: E402

for mode in ("exact", "hybrid"):
    for n in (2, 3):
        li = LatticeIndex(n_electrons=2, max_n=n, nuclear_charge=2,
                          vee_method="slater_full", h1_method=mode)
        E, _ = li.compute_ground_state(n_states=1)
        print("RESULT", mode, n, round(float(E[0]), 6), flush=True)
