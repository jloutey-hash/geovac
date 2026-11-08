"""DLMF 30.3 cross-check of the GeoVac angular separation constant A(c^2).

The Paper 58 / Paper 11 eta-equation (m=0, b=0, homonuclear) is
    d/deta[(1-eta^2) G'] + (-A + c^2 eta^2) G = 0
DLMF 30.2.1 spheroidal equation is
    d/dz[(1-z^2) w'] + (lambda + gamma^2 (1-z^2) - mu^2/(1-z^2)) w = 0
Matching term by term gives the dictionary
    gamma^2 = -c^2,     lambda_DLMF = c^2 - A.
Hence DLMF 30.3.3  lambda(0) = n(n+1)      ->  A(0) = -n(n+1)
and  DLMF 30.3.4  -1 < dlambda/dgamma^2 < 0 ->  0 < dA/dc^2 < 1.

Independent route: scipy.special.obl_cv, per the docstring claim in
geovac/molecular_sturmian.py::_angular_sep_const ("equals -obl_cv").
"""
import numpy as np
from scipy.special import obl_cv
from geovac.molecular_sturmian import _angular_sep_const

def A_of_c2(c2, m=0, n_sph=0, n_basis=60):
    return _angular_sep_const(m, n_sph, np.sqrt(max(c2, 0.0)), b=0.0, n_basis=n_basis)

print("=" * 78)
print("DLMF 30.3.3  ->  A(c^2 = 0) = -n(n+1)")
print("=" * 78)
ok33 = True
for m in (0, 1, 2):
    for n_sph in (0, 1, 2):
        n = m + n_sph
        A0 = A_of_c2(0.0, m, n_sph)
        pred = -n * (n + 1)
        good = abs(A0 - pred) < 1e-9
        ok33 &= good
        print(f"  m={m} n_sph={n_sph} (n={n}): A(0)={A0:+.12f}  -n(n+1)={pred:+d}  {'OK' if good else 'FAIL'}")

print()
print("=" * 78)
print("Independent route: A(c^2) vs -obl_cv(m, n, c)")
print("=" * 78)
maxdev = 0.0
for m in (0, 1):
    for n_sph in (0, 1, 2):
        n = m + n_sph
        for c2 in (0.5, 2.0, 8.0, 20.0, 50.0):
            c = np.sqrt(c2)
            ours = A_of_c2(c2, m, n_sph, n_basis=80)
            ref = -obl_cv(m, n, c)
            dev = abs(ours - ref)
            maxdev = max(maxdev, dev)
            if dev > 1e-6:
                print(f"  m={m} n={n} c2={c2:5.1f}: ours={ours:+.9f} scipy={ref:+.9f} dev={dev:.2e}  <-- ")
print(f"  max |A_geovac - (-obl_cv)| over grid = {maxdev:.3e}")

print()
print("=" * 78)
print("DLMF 30.3.4  ->  0 < dA/dc^2 < 1   (central differences)")
print("=" * 78)
h = 1e-4
viol = []
rows = []
for m in (0, 1, 2):
    for n_sph in (0, 1, 2):
        n = m + n_sph
        slopes = []
        for c2 in np.linspace(0.25, 60.0, 40):
            d = (A_of_c2(c2 + h, m, n_sph, 80) - A_of_c2(c2 - h, m, n_sph, 80)) / (2 * h)
            slopes.append(d)
            if not (0.0 < d < 1.0):
                viol.append((m, n, c2, d))
        rows.append((m, n, min(slopes), max(slopes)))
        print(f"  m={m} n={n}: dA/dc^2 in [{min(slopes):.6f}, {max(slopes):.6f}]")
print(f"  violations of 0 < dA/dc^2 < 1 : {len(viol)}")
for v in viol[:5]:
    print("   ", v)

print()
print("=" * 78)
print("Monotonicity of A across the Paper 58 front sweep (H2+, Z=1)")
print("=" * 78)
# c^2 = -R^2 E/2 ; sweep R over the front range used in Paper 58 (R* ~ 1.5-17)
for R in (1.0, 2.0, 4.0, 8.0, 16.0):
    E = -1.1  # representative bound electronic energy scale (Ha), sign only matters via c2>0
    c2 = -R**2 * E / 2.0
    A = A_of_c2(c2, 0, 0, 80)
    lam = c2 - A
    print(f"  R={R:5.1f}  c^2={c2:8.3f}  A={A:+10.5f}  lambda_DLMF=c^2-A={lam:+10.5f}")

print()
print("SUMMARY:", "30.3.3 PASS" if ok33 else "30.3.3 FAIL",
      "| independent-route max dev %.1e" % maxdev,
      "| 30.3.4 %s" % ("PASS" if not viol else f"FAIL ({len(viol)})"))
