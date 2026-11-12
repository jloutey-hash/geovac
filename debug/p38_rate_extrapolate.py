import numpy as np
from scipy.linalg import eigh
import importlib.util, os
spec = importlib.util.spec_from_file_location("r", os.path.join(os.path.dirname(__file__), "p38_rate_reduction.py"))
r = importlib.util.module_from_spec(spec); spec.loader.exec_module(r)
ns = [800, 1600, 3200, 6400]; v = []
for n in ns:
    lam = eigh(r.A_matrix(np.arange(1, n + 1)), eigvals_only=True, subset_by_index=[0, 0])[0]
    v.append(n * lam / 2); print(n, "%.10f" % v[-1], flush=True)
v = np.array(v)
r1 = 2 * v[1:] - v[:-1]; print("Richardson1:", ["%.10f" % x for x in r1])
r2 = (4 * r1[1:] - r1[:-1]) / 3; print("Richardson2:", ["%.10f" % x for x in r2])
print("Galerkin lambda_2: 2.7547547422")
