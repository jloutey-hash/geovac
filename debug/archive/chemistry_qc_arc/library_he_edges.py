"""Edge count of GeometricLattice(max_n) -> predicted 1-norm shift of the
retired hop (each edge contributes 2 spins x (1/16) to the JW 1-norm when
its XX+YY strings do not cancel against two-body strings)."""
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from geovac.lattice import GeometricLattice  # noqa: E402

old = {2: 11.175482, 3: 74.207417, 4: 275.718, 5: 790.007}
for n in (2, 3, 4, 5):
    A = GeometricLattice(max_n=n, nuclear_charge=2).adjacency
    E = int((A != 0).sum() // 2)
    print(n, "edges", E, "shift", E / 8, "pred lambda_exact", old[n] - E / 8)
