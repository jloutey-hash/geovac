"""Step 2: verify the two angular rules the 3-electron term inventory rests on.

Writing out the dot products in <Phi| F H F |Phi> (F = sum_{i<j} f_ij, multiplicative)
for an s-only basis gives only four shapes:

  (i)   same pair                      -> two-body machinery
  (ii)  2 legs, shared vertex, scalar   -> RULE A below
  (iii) 2 legs, shared vertex, vec-vec  -> RULE B below  (= the module's L3/kA object)
  (iv)  3 legs, closed triangle, scalar -> the delta_abc/(2a+1)^2 rule (step 1, validated)

There is NO vector-leg closed triangle: a triangle needs three factors, and the only
three-factor terms come from F^2 x (1/r_kl), whose legs are all scalar.

RULE A:  <P_a(r1.r2) P_b(r1.r3)>  =  delta_a0 delta_b0
         (two legs at a shared vertex see independent orientations -> L=0 only.
          This is *why* the module's two-body kernels are L=0 kernels.)

RULE B:  < (r12hat . r13hat) h(r12) g(r13) >  =  <h (r12hat.r1hat)> <g (r13hat.r1hat)>
         (exactness of the kA-factorization used by geovac.transcorrelated_sturmian.
          three_body -- an independent check of promoted code.)
"""
import numpy as np
from numpy.polynomial.legendre import legval

rng = np.random.default_rng(20260825)


def rand_dirs(m):
    v = rng.normal(size=(m, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def P(L, x):
    cf = np.zeros(L + 1); cf[L] = 1.0
    return legval(x, cf)


print("=" * 74)
print("RULE A:  <P_a(r1.r2) P_b(r1.r3)> = delta_a0 delta_b0")
print("=" * 74)
n = 600000
u1, u2, u3 = rand_dirs(n), rand_dirs(n), rand_dirs(n)
x12 = np.einsum("ij,ij->i", u1, u2)
x13 = np.einsum("ij,ij->i", u1, u3)
print(f"{'a':>2}{'b':>3}{'rule':>10}{'sampled':>12}{'mc_err':>10}{'ok':>6}")
okA = True
for a, b in [(0,0), (1,1), (2,2), (1,0), (0,1), (2,1), (3,3), (1,2)]:
    rule = 1.0 if (a == 0 and b == 0) else 0.0
    v = P(a, x12) * P(b, x13)
    got, err = float(v.mean()), float(v.std() / np.sqrt(n))
    ok = abs(got - rule) < max(5 * err, 2e-4); okA &= ok
    print(f"{a:>2}{b:>3}{rule:>10.4f}{got:>12.6f}{err:>10.1e}{'OK' if ok else 'FAIL':>6}")
print(f"RULE A verified: {okA}")

print()
print("=" * 74)
print("RULE B:  <(r12hat.r13hat) h g> = <h (r12hat.r1hat)> <g (r13hat.r1hat)>")
print("         (the kA factorization behind geovac ... three_body, checked exactly)")
print("=" * 74)
# fixed radii; sample orientations.  h, g are the u'-style radial factors.
print(f"{'r1':>5}{'r2':>5}{'r3':>5}{'gamma':>7}{'factorized':>14}{'sampled':>14}"
      f"{'mc_err':>10}{'ok':>6}")
okB = True
for (r1, r2, r3, gam) in [(1.0, 1.0, 1.0, 0.7), (0.7, 1.9, 1.3, 0.7),
                          (2.1, 0.5, 1.6, 1.4), (1.2, 1.2, 0.4, 0.3)]:
    v1 = r1 * u1; v2 = r2 * u2; v3 = r3 * u3
    d12 = v1 - v2; d13 = v1 - v3
    n12 = np.linalg.norm(d12, axis=1); n13 = np.linalg.norm(d13, axis=1)
    h = 0.5 * np.exp(-gam * n12)              # u'(r12)
    g = 0.5 * np.exp(-gam * n13)
    cos_1213 = np.einsum("ij,ij->i", d12, d13) / (n12 * n13)
    full = h * g * cos_1213
    # factorized: each leg projected on r1hat, averaged independently
    projA = h * np.einsum("ij,ij->i", d12, u1) / n12
    projB = g * np.einsum("ij,ij->i", d13, u1) / n13
    fac = float(projA.mean()) * float(projB.mean())
    got, err = float(full.mean()), float(full.std() / np.sqrt(n))
    ok = abs(got - fac) < max(5 * err, 1e-5); okB &= ok
    print(f"{r1:>5.1f}{r2:>5.1f}{r3:>5.1f}{gam:>7.2f}{fac:>14.8f}{got:>14.8f}"
          f"{err:>10.1e}{'OK' if ok else 'FAIL':>6}")
print(f"RULE B verified: {okB}")

print()
print("=" * 74)
print("TERM INVENTORY for <Phi| F H F |Phi>, s-only, N=3")
print("=" * 74)
rows = [
    ("F^2 * v(r_i)            ", "<=2 legs, shared vertex/same pair", "RULE A (L=0) / two_body"),
    ("F^2 * 1/r_kl            ", "3 legs, CLOSED TRIANGLE, scalar  ", "delta_abc/(2a+1)^2  [step 1]"),
    ("F^2 |grad_i Phi|^2      ", "<=2 legs, shared vertex          ", "RULE A (L=0)"),
    ("F Phi grad_i F . grad_i Phi", "2 legs; dot -> scalar pair kernel", "RULE A (L=0)"),
    ("Phi^2 |grad_i F|^2  j=k ", "same pair                        ", "two_body"),
    ("Phi^2 |grad_i F|^2  j!=k", "2 legs, shared vertex, VEC-VEC   ", "RULE B = kA / three_body"),
]
print(f"{'term':<28}{'shape':<35}{'handled by':<30}")
for t, s, h in rows:
    print(f"{t:<28}{s:<35}{h:<30}")
print()
print("=> no vector-leg closed triangle exists; remaining work is ASSEMBLY, not")
print("   new angular derivation.  L>0 multipoles are needed ONLY in the triangle.")
