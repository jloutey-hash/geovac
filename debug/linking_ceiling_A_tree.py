"""PL1 check (docs/predictions/2026-10-16_linking_ceiling.md): a four-density tree's
fibre equals a one-dimensional four-Bessel moment, for star and chain alike.
Uses the fast-decaying variant g_v = (c_v k^2 + 1)^{-3/2}, whose 1-D Fourier transform is
|y| K1(a|y|)/(pi c), so the identity reads
  int_{R^3} prod_v g_v(L_v(k)) d^3k = 2^4/(pi prod c_v) int_0^inf y^4 prod_v K1(a_v y) dy.
The delta-function mechanism is the same as for 1/Delta; this variant converges fast in 3-D."""
import numpy as np, mpmath as mp, json, time

def tree3d(c, forms, N, U):
    u = np.linspace(-U, U, N); h = u[1] - u[0]
    k = np.sinh(u); w = np.cosh(u) * h
    P, Q, R = np.meshgrid(k, k, k, indexing="ij", sparse=True)
    W = w[:, None, None] * w[None, :, None] * w[None, None, :]
    tot = W.copy()
    for cv, L in zip(c, forms(P, Q, R)):
        tot = tot * (cv * L**2 + 1.0) ** -1.5
    return tot.sum()

star  = lambda p, q, r: (p + q + r, p, q, r)
chain = lambda p, q, r: (p, q - p, r - q, r)
mp.mp.dps = 30
def moment(c):
    a = [1 / mp.sqrt(cv) for cv in c]
    I = mp.quad(lambda y: y**4 * mp.besselk(1, a[0]*y) * mp.besselk(1, a[1]*y)
                * mp.besselk(1, a[2]*y) * mp.besselk(1, a[3]*y), [0, 1, mp.inf])
    return 16 / (mp.pi * mp.fprod(c)) * I

out = []
for c in [(0.2, 0.1, 0.15, 0.05), (0.25, 0.08, 0.12, 0.2), (0.1, 0.1, 0.1, 0.1)]:
    ref = moment(c)
    row = {"c": c, "bessel": float(ref)}
    for name, f in [("star", star), ("chain", chain)]:
        v1 = tree3d(c, f, 241, 9.0); v2 = tree3d(c, f, 321, 9.0)
        row[name] = v2; row[name + "_conv"] = abs(v2 - v1)
        row[name + "_rel"] = abs(v2 - float(ref)) / float(ref)
    out.append(row); print(row)
json.dump(out, open("debug/data/linking_ceiling_A_tree.json", "w"), indent=1)
