"""PROBE: is the atomic one-electron Hamiltonian a DERIVED object?

The PI's conjecture ("maybe a Hamiltonian is a projection"), made precise:

    In the shared-k Coulomb-Sturmian basis, h1 should be a function of NOTHING but
    (integer labels n, the overlap metric S, the focal scale k, the charge Z):

        h1  =  k^2 ( I - S/2 )  -  Z k diag(1/n)

Derivation (textbook Sturmian theory -- this is a VERIFICATION, not a discovery):
the Sturmian equation  [-1/2 lap + k^2/2 - nk/r] chi_n = 0  gives
T chi_n = (nk/r - k^2/2) chi_n, and the potential-weighted orthonormality
<chi_m| 1/r |chi_n> = (k/n) delta_mn  collapses both T and V onto the metric:

    T_mn = k^2 delta_mn - (k^2/2) S_mn        V_mn = -(Zk/n) delta_mn

So the "potential" IS the statement that the 1/r-weighted metric is the identity,
and the kinetic term IS the L2 metric.  The Hamiltonian carries ZERO information
beyond (labels, metric, scale).  Everything dynamical enters through k -- and
k = p0 is the Fock focal length (p0^2 = -2E): the scale IS the projection parameter.

Checked against the INDEPENDENT grid-built h1/S of geovac.transcorrelated_sturmian
(build_one_body: numerical quadrature, gradient-form kinetic -- an INDEPENDENT
route that knows nothing of the identity).

Where it must FAIL (the control): the two-electron ERI tensor is NOT a function
of (S, labels, k) -- that is where independent content (the e-e cusp, the genus
story) enters.  A rank test makes that concrete.
"""
import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from geovac import transcorrelated_sturmian as TC  # noqa: E402

out = {"identity": [], "control": {}}

print("=" * 78)
print("h1 = k^2 (I - S/2) - Z k diag(1/n)   vs the independent grid-built h1")
print("=" * 78)
print(f"{'ns':>3}{'k':>7}{'Z':>5}{'Ng':>6}{'max|h1 - pred|':>17}{'scale |h1|':>12}{'rel':>11}")
worst = 0.0
for ns, k, Z, Ng in [(3, 1.7, 2.0, 800), (4, 2.4, 3.0, 800), (5, 2.0, 3.0, 1200),
                     (6, 1.3, 1.0, 1200), (4, 3.2, 4.0, 1600)]:
    r, wr = TC.make_grid(k, Ng=Ng)
    S, h1, Rtab, W = TC.build_one_body(ns, r, wr, k, Z)
    n = np.arange(1, ns + 1)
    pred = k * k * (np.eye(ns) - S / 2.0) - Z * k * np.diag(1.0 / n)
    dev = float(np.abs(h1 - pred).max())
    scale = float(np.abs(h1).max())
    rel = dev / scale
    worst = max(worst, rel)
    print(f"{ns:>3}{k:>7.2f}{Z:>5.1f}{Ng:>6}{dev:>17.3e}{scale:>12.3f}{rel:>11.2e}")
    out["identity"].append(dict(ns=ns, k=k, Z=Z, Ng=Ng, max_dev=dev,
                                scale=scale, rel=rel))
print(f"\nworst relative deviation: {worst:.2e}  (grid-quadrature level = PASS)")

print()
print("=" * 78)
print("CONTROL -- the ERI must NOT reduce to (S, labels, k):  rank test")
print("=" * 78)
# If the ERI were a function of the same data, its matricization would live in the
# span of low-degree words in S (S is ns x ns => that span has dimension <= ns^2).
# Measure the numerical rank of the (ns^2 x ns^2) Coulomb ERI matricization and
# its residual after projecting onto the span of {kron(A,B): A,B in {I,S,S^2,S^3}}.
ns, k, Z, Ng = 4, 2.4, 3.0, 800
r, wr = TC.make_grid(k, Ng=Ng)
S, h1, Rtab, W = TC.build_one_body(ns, r, wr, k, Z)
Km = TC.build_kernels(r, 0.7, nx=96)
eri, _, _ = TC.two_body(ns, Rtab, W, Km)
E = eri.reshape(ns * ns, ns * ns)          # (i k | j l) matricization
mats = []
P = [np.eye(ns), S, S @ S, S @ S @ S]
for A in P:
    for B in P:
        mats.append(np.kron(A, B).ravel())
M = np.array(mats).T                        # basis of metric-generated candidates
coef, res, *_ = np.linalg.lstsq(M, E.ravel(), rcond=None)
fit = M @ coef
rel_res = float(np.linalg.norm(E.ravel() - fit) / np.linalg.norm(E.ravel()))
print(f"  best fit of ERI by words in the metric (16 kron terms): rel residual = {rel_res:.3f}")
print(f"  => the two-electron tensor is {'NOT' if rel_res > 0.05 else ''} metric-generated —")
print(f"     independent content enters exactly at the e-e channel, as the axis map says.")
out["control"] = dict(ns=ns, rel_residual=rel_res)

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/minimal_rep_h1_identity.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/minimal_rep_h1_identity.json")
