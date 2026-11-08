"""He regression, independent route: LatticeIndex (determinant FCI, slater_full)
with h1_method='hybrid' (paper_fci_atoms Table I) vs 'exact' (physical).

Second algorithmic route for the physical-Hamiltonian He numbers produced by
debug/he_graph_regression.py (casimir_ci singlet M_L=0 route).
"""
import os
import sys
import time
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
warnings.filterwarnings("ignore")

from geovac.lattice_index import LatticeIndex  # noqa: E402

E_EX = -2.903724377034


def main() -> None:
    nmaxes = [int(a) for a in sys.argv[1:]] or [2, 3, 4]
    for n in nmaxes:
        for mode in ("hybrid", "exact"):
            t0 = time.time()
            li = LatticeIndex(n_electrons=2, max_n=n, nuclear_charge=2,
                              vee_method="slater_full", h1_method=mode)
            res = li.compute_ground_state(n_states=1)
            E = float(res[0][0]) if isinstance(res, tuple) else float(res[0])
            print(f"n_max={n} h1={mode:6s} N_SD={li.n_sd if hasattr(li,'n_sd') else '?'} "
                  f"E={E:.6f} err={100*(E-E_EX)/abs(E_EX):+.3f}%  [{time.time()-t0:.0f}s]",
                  flush=True)


if __name__ == "__main__":
    main()
