"""He n_max=5, route 2 for the OpenFermion drop audit: dict-built fermion operator, plain (unscaled) jordan_wigner
with the deletion threshold lowered to 1e-13. Counts non-identity terms above 1e-11 and prints the decade histogram."""
import sys, warnings, numpy as np
warnings.simplefilter("ignore"); sys.path.insert(0, "debug")
import openfermion_tolerance_audit as A
from openfermion import jordan_wigner
from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import JordanWignerEncoder
n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
li = LatticeIndex(n_electrons=2, max_n=n, nuclear_charge=2, vee_method='slater_full', h1_method='exact')
with A.tolerance(1e-13):
    q = jordan_wigner(JordanWignerEncoder(li).build_fermion_operator())
m = np.array([abs(v) for k, v in q.terms.items() if k])
print(f"He n={n}: tight-route non-identity terms above 1e-11: {(m > 1e-11).sum()}, "
      f"lambda_nonid {m[m > 1e-11].sum():.6f}")
print("decades 1e-16..1e-5:", np.histogram(np.log10(m[m > 0]), bins=range(-16, -4))[0].tolist(), flush=True)
