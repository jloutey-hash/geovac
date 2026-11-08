"""EXPERIMENT B: does a per-shell lambda destroy the angular (Gaunt) sparsity?

The corpus's objection to mixed exponents is Track DF Sprint 5: heterogeneous per-pair
Z_eff needed Loewdin orthogonalization, which inflated the Pauli count 14x (1711 vs 120).
But that result is MOLECULAR -- the orthogonalization ran ACROSS centers.

Claim under test: on a SINGLE center the overlap matrix is block-diagonal in (l, m) by
angular orthogonality REGARDLESS of the radial exponents, so S^{-1/2} is block-diagonal
too, Loewdin mixes only same-(l,m) functions, and the Gaunt support of the ERI tensor --
which depends only on angular labels -- is untouched.  Mixed exponents would then be free
of the sparsity cost FOR ATOMS.

The load-bearing measurement is NOT "block-diagonal transforms don't mix blocks" (trivial).
It is whether the PHYSICAL situation produces a block-diagonal transform.  So the
different-l overlaps are computed by explicit 3D quadrature with genuinely different
exponents -- measured, not assumed -- on one center and on two, using the validated
mixed-exponent engine debug/two_center_grid_lm.py.
"""
import importlib.util
import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

_spec = importlib.util.spec_from_file_location(
    "tcg", os.path.join(HERE, "two_center_grid_lm.py"))
tcg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tcg)

out = {}

# Deliberately MISMATCHED exponents: a tight core-like and a diffuse valence-like value.
BASIS = [(2.70, 0, 0), (0.65, 0, 0), (1.90, 1, 0), (0.55, 1, 0),
         (1.40, 1, 1), (0.80, 2, 0)]

print("=" * 80)
print("B1  same-centre overlaps between DIFFERENT-l functions at DIFFERENT exponents")
print("    (3D quadrature; if these are not zero the whole argument fails)")
print("=" * 80)
eng = tcg.TwoCenterLM(R=2.0, nr=900, nu=40, nphi=40, rmax=45.0, Lmax=10, real=True)
print(f"{'zeta_i':>8}{'l_i':>4}{'m_i':>4}{'zeta_j':>8}{'l_j':>4}{'m_j':>4}{'overlap':>14}")
same_diff_l = []
same_same_l = []
for i, (za, la, ma) in enumerate(BASIS):
    for j, (zb, lb, mb) in enumerate(BASIS):
        if j <= i:
            continue
        s = eng.overlap((za, la, ma, "A"), (zb, lb, mb, "A"))
        tag = "same l,m" if (la, ma) == (lb, mb) else "diff l,m"
        (same_same_l if (la, ma) == (lb, mb) else same_diff_l).append(abs(s))
        if (la, ma) != (lb, mb) or True:
            print(f"{za:>8.2f}{la:>4}{ma:>4}{zb:>8.2f}{lb:>4}{mb:>4}{s:>14.3e}   {tag}")
print()
print(f"  max |overlap| across DIFFERENT (l,m), mixed exponents : {max(same_diff_l):.3e}")
print(f"  typical |overlap| within the SAME (l,m)               : {max(same_same_l):.3e}")
print(f"  => single-centre S is (l,m)-block-diagonal to quadrature precision: "
      f"{max(same_diff_l) < 1e-10}")
out["B1"] = dict(max_diff_lm=float(max(same_diff_l)),
                 max_same_lm=float(max(same_same_l)),
                 block_diagonal=bool(max(same_diff_l) < 1e-10))

print()
print("=" * 80)
print("B2  CONTROL -- the same functions across TWO centres")
print("    (this is the situation Track DF Sprint 5 measured; l-mixing must appear)")
print("=" * 80)
print(f"{'zeta_i':>8}{'l_i':>4}{'zeta_j':>8}{'l_j':>4}{'overlap A-B':>14}")
cross = []
for (za, la, ma) in BASIS[:4]:
    for (zb, lb, mb) in BASIS[:4]:
        if (la, ma) == (lb, mb):
            continue
        s = eng.overlap((za, la, ma, "A"), (zb, lb, mb, "B"))
        cross.append(abs(s))
        print(f"{za:>8.2f}{la:>4}{zb:>8.2f}{lb:>4}{s:>14.3e}")
print()
print(f"  max |cross-centre overlap| between DIFFERENT l : {max(cross):.3e}")
ratio = max(cross) / max(max(same_diff_l), 1e-300)
print(f"  ratio to the single-centre value              : {ratio:.2e}x")
out["B2"] = dict(max_cross_diff_l=float(max(cross)), ratio_to_single_centre=float(ratio))

print()
print("=" * 80)
print("B3  consequence: Gaunt support of the ERI tensor under Loewdin")
print("=" * 80)
n = len(BASIS)
S = np.zeros((n, n))
for i, (za, la, ma) in enumerate(BASIS):
    for j, (zb, lb, mb) in enumerate(BASIS):
        S[i, j] = eng.overlap((za, la, ma, "A"), (zb, lb, mb, "A"))
ev, U = np.linalg.eigh(S)
X = U @ np.diag(ev ** -0.5) @ U.T
lm = [(l, m) for (_, l, m) in BASIS]
offblock = max(abs(X[i, j]) for i in range(n) for j in range(n) if lm[i] != lm[j])
print(f"  max |X_ij| connecting different (l,m)  : {offblock:.3e}   "
      f"(X = S^-1/2, single centre)")

# Gaunt support mask on the (l,m) labels, then transform a support-respecting tensor
def gaunt_allowed(a, b, c, d):
    la, ma = lm[a]; lb, mb = lm[b]; lc, mc = lm[c]; ld, md = lm[d]
    ok = False
    for L in range(0, 7):
        if abs(la - lc) <= L <= la + lc and abs(lb - ld) <= L <= lb + ld \
           and (la + lc + L) % 2 == 0 and (lb + ld + L) % 2 == 0 \
           and (ma - mc) == -(mb - md) and abs(ma - mc) <= L:
            ok = True
    return ok

mask = np.zeros((n, n, n, n), dtype=bool)
for a in range(n):
    for b in range(n):
        for c in range(n):
            for d in range(n):
                mask[a, b, c, d] = gaunt_allowed(a, b, c, d)
rng = np.random.default_rng(0)
T = rng.normal(size=(n,) * 4) * mask
def transform(T, M):
    T = np.einsum("pqrs,pa->aqrs", T, M)
    T = np.einsum("aqrs,qb->abrs", T, M)
    T = np.einsum("abrs,rc->abcs", T, M)
    return np.einsum("abcs,sd->abcd", T, M)
Tl = transform(T, X)
nz0 = int(mask.sum())
nz1 = int((np.abs(Tl) > 1e-12).sum())
# CONTROL: a transform that DOES mix l (what two centres force)
Xmix = np.linalg.qr(rng.normal(size=(n, n)))[0]
Tm = transform(T, Xmix)
nz2 = int((np.abs(Tm) > 1e-12).sum())
print(f"  nonzero ERI entries, Gaunt support        : {nz0:5d} / {n**4}")
print(f"  after single-centre Loewdin (X)           : {nz1:5d}   "
      f"fill-in {nz1 - nz0:+d}")
print(f"  after an l-MIXING transform (control)     : {nz2:5d}   "
      f"fill-in {nz2 - nz0:+d}  ({nz2 / max(nz0,1):.2f}x)")
out["B3"] = dict(max_offblock_X=float(offblock), nz_gaunt=nz0,
                 nz_after_single_centre=nz1, nz_after_mixing_control=nz2)

print()
print("=" * 80)
print("VERDICT")
print("=" * 80)
print(f"  single centre, mixed exponents: different-l overlap {max(same_diff_l):.1e} "
      f"-> S^-1/2 off-block {offblock:.1e} -> ERI fill-in {nz1 - nz0:+d}")
print(f"  two centres  (control)        : different-l overlap {max(cross):.1e}  "
      f"-> l-mixing transform gives fill-in {nz2 - nz0:+d}")

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/r12ci_lambda_sparsity.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_lambda_sparsity.json")
