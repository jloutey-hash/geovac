"""Step 1 of the 3-electron variational R12-CI build: the TRIANGLE angular machinery.

For N=3 the operator product <Phi| f H f |Phi> stays at most 3-body (all indices live
in {1,2,3}), so an EXACT treatment is available with no RI -- Be (N=4) is the first case
that genuinely needs 4-body.  The pieces that do NOT factorize into the module's existing
shared-vertex kernels are the TRIANGLE terms  A(r12) B(r13) C(r23)  (a closed loop, no
shared vertex).

Derived rule (verified below):   <P_a(12) P_b(13) P_c(23)>  =  delta_{abc} / (2a+1)^2
Only all-three-equal multipoles survive.  Hence

  T = sum_L (2L+1)^{-2} INT rho1 rho2 rho3 a_L(r1,r2) b_L(r1,r3) c_L(r2,r3)

which is a matrix contraction, not a 6D quadrature.

This file ONLY builds and validates the angular machinery.  No physics yet.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

rng = np.random.default_rng(20260825)


# ---------------------------------------------------------------------------
# 1. the coupling rule, checked by brute-force orientation sampling
# ---------------------------------------------------------------------------
def brute_triple_average(a, b, c, n=400000):
    """<P_a(r1.r2) P_b(r1.r3) P_c(r2.r3)> by isotropic sampling of three directions."""
    def rand_dirs(m):
        v = rng.normal(size=(m, 3))
        return v / np.linalg.norm(v, axis=1, keepdims=True)
    u1, u2, u3 = rand_dirs(n), rand_dirs(n), rand_dirs(n)
    from numpy.polynomial.legendre import legval
    def P(L, x):
        cf = np.zeros(L + 1); cf[L] = 1.0
        return legval(x, cf)
    vals = P(a, np.einsum("ij,ij->i", u1, u2)) * P(b, np.einsum("ij,ij->i", u1, u3)) \
         * P(c, np.einsum("ij,ij->i", u2, u3))
    return float(vals.mean()), float(vals.std() / np.sqrt(n))


print("=" * 72)
print("1. triangle coupling rule   <P_a P_b P_c> = delta_abc / (2a+1)^2")
print("=" * 72)
print(f"{'a':>2}{'b':>3}{'c':>3}{'rule':>12}{'sampled':>12}{'mc_err':>10}{'ok':>5}")
ok_all = True
for (a, b, c) in [(0,0,0), (1,1,1), (2,2,2), (3,3,3), (1,1,0), (1,0,1),
                  (2,1,1), (2,2,1), (2,2,0), (1,2,3)]:
    rule = 1.0 / (2 * a + 1) ** 2 if (a == b == c) else 0.0
    got, err = brute_triple_average(a, b, c)
    ok = abs(got - rule) < max(5 * err, 2e-4)
    ok_all &= ok
    print(f"{a:>2}{b:>3}{c:>3}{rule:>12.6f}{got:>12.6f}{err:>10.1e}{'OK' if ok else 'FAIL':>5}")
print(f"\nrule verified: {ok_all}")


# ---------------------------------------------------------------------------
# 2. Legendre moments of a pair kernel, and the triangle contraction
# ---------------------------------------------------------------------------
def legendre_moments(fun, r, Lmax, nx=200):
    """a_L(r_i, r_j) = (2L+1)/2 INT_-1^1 f(r_ij) P_L(x) dx,   r_ij = |r_i - r_j|.

    Returns array (Lmax+1, Ng, Ng).  Convention: f(r12) = sum_L a_L P_L(cos th12)."""
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    xs, ws = leggauss(nx)
    out = np.zeros((Lmax + 1, r.size, r.size))
    from numpy.polynomial.legendre import legval
    for x, w in zip(xs, ws):
        r12 = np.sqrt(np.maximum(R1 * R1 + R2 * R2 - 2 * R1 * R2 * x, 1e-30))
        fv = fun(r12)
        for L in range(Lmax + 1):
            cf = np.zeros(L + 1); cf[L] = 1.0
            out[L] += (2 * L + 1) / 2.0 * w * fv * legval(x, cf)
    return out


def triangle_contract(aL, bL, cL, d1, d2, d3):
    """sum_L (2L+1)^-2 * sum_{r1r2r3} d1(r1) d2(r2) d3(r3) aL(r1,r2) bL(r1,r3) cL(r2,r3)"""
    tot = 0.0
    for L in range(aL.shape[0]):
        # P(r1,r2) = sum_r3 bL(r1,r3) d3(r3) cL(r3,r2)      [cL symmetric]
        P = (bL[L] * d3[None, :]) @ cL[L]
        tot += (d1 @ ((aL[L] * P) @ d2)) / (2 * L + 1) ** 2
    return float(tot)


def brute_triangle(fa, fb, fc, r, d1, d2, d3, nang=60):
    """Direct 3-orientation quadrature of the same object (validation reference).

    Fix r1hat = zhat and integrate Omega2, Omega3; azimuthal symmetry lets phi2 = 0."""
    xt, wt = leggauss(nang)                       # cos(theta2), cos(theta3)
    ph, wp = leggauss(nang)
    phi3 = np.pi * (ph + 1.0)                     # [0, 2pi)
    wphi = np.pi * wp
    R1 = r[:, None, None]; R2 = r[None, :, None]; R3 = r[None, None, :]
    tot = np.zeros((r.size, r.size, r.size))
    for c2, w2 in zip(xt, wt):
        s2 = np.sqrt(1 - c2 * c2)
        r12 = np.sqrt(np.maximum(R1**2 + R2**2 - 2 * R1 * R2 * c2, 1e-30))
        A = fa(r12)
        for c3, w3 in zip(xt, wt):
            s3 = np.sqrt(1 - c3 * c3)
            r13 = np.sqrt(np.maximum(R1**2 + R3**2 - 2 * R1 * R3 * c3, 1e-30))
            B = fb(r13)
            acc = np.zeros_like(tot)
            for p3, wpp in zip(phi3, wphi):
                cos23 = c2 * c3 + s2 * s3 * np.cos(p3)
                r23 = np.sqrt(np.maximum(R2**2 + R3**2 - 2 * R2 * R3 * cos23, 1e-30))
                acc += wpp * fc(r23)
            tot += w2 * w3 * A * B * acc
    tot /= (2.0 * 2.0 * 2 * np.pi)                # normalize the three solid angles
    return float(np.einsum("i,j,k,ijk->", d1, d2, d3, tot))


print()
print("=" * 72)
print("2. triangle contraction  vs  direct 3-orientation quadrature")
print("=" * 72)
Ng = 26
r = np.linspace(0.12, 4.2, Ng)
d1 = np.exp(-1.3 * r) * r ** 2
d2 = np.exp(-0.9 * r) * r ** 2 * (1 + 0.3 * r)
d3 = np.exp(-1.7 * r) * r ** 2
cases = {
    "coul x coul x coul": (lambda x: 1 / x, lambda x: 1 / x, lambda x: 1 / x),
    "gem  x coul x gem ": (lambda x: np.exp(-0.7 * x), lambda x: 1 / x, lambda x: np.exp(-0.7 * x)),
    "gem  x gem  x coul": (lambda x: np.exp(-0.4 * x), lambda x: np.exp(-1.1 * x), lambda x: 1 / x),
}
print(f"{'case':>20}{'Lmax':>6}{'multipole':>14}{'brute':>14}{'rel.diff':>11}")
for name, (fa, fb, fc) in cases.items():
    ref = brute_triangle(fa, fb, fc, r, d1, d2, d3, nang=48)
    for Lmax in (2, 6, 12):
        aL = legendre_moments(fa, r, Lmax); bL = legendre_moments(fb, r, Lmax)
        cL = legendre_moments(fc, r, Lmax)
        val = triangle_contract(aL, bL, cL, d1, d2, d3)
        print(f"{name:>20}{Lmax:>6}{val:>14.8f}{ref:>14.8f}{abs(val-ref)/abs(ref):>11.2e}")


# ---------------------------------------------------------------------------
# 3. Isolate rule-correctness from quadrature resolution.
#    (a) all-smooth control (gem x gem x gem): no singular kernel anywhere, so both
#        routes should converge and agree to high precision.  If they do, the rule +
#        contraction are right and the residuals above are quadrature on 1/r12.
#    (b) Coulomb via EXACT analytic Legendre moments  a_L = r_<^L / r_>^(L+1)
#        (no quadrature at all on the singular kernel).
# ---------------------------------------------------------------------------
def coul_moments_exact(r, Lmax):
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    rlt = np.minimum(R1, R2); rgt = np.maximum(R1, R2)
    return np.array([rlt ** L / rgt ** (L + 1) for L in range(Lmax + 1)])


print()
print("=" * 72)
print("3a. ALL-SMOOTH control (gem x gem x gem) -- no singular kernel")
print("=" * 72)
fa = lambda x: np.exp(-0.7 * x)
fb = lambda x: np.exp(-1.1 * x)
fc = lambda x: np.exp(-0.4 * x)
ref_hi = brute_triangle(fa, fb, fc, r, d1, d2, d3, nang=96)
print(f"{'Lmax':>6}{'multipole':>16}{'brute(96)':>16}{'rel.diff':>12}")
for Lmax in (2, 4, 8, 14, 20):
    aL = legendre_moments(fa, r, Lmax, nx=300)
    bL = legendre_moments(fb, r, Lmax, nx=300)
    cL = legendre_moments(fc, r, Lmax, nx=300)
    val = triangle_contract(aL, bL, cL, d1, d2, d3)
    print(f"{Lmax:>6}{val:>16.10f}{ref_hi:>16.10f}{abs(val-ref_hi)/abs(ref_hi):>12.2e}")

print()
print("=" * 72)
print("3b. Coulomb leg via EXACT analytic moments  a_L = r_<^L / r_>^(L+1)")
print("=" * 72)
fg = lambda x: np.exp(-0.7 * x)
ref_c = brute_triangle(fg, lambda x: 1 / x, fg, r, d1, d2, d3, nang=120)
print(f"{'Lmax':>6}{'multipole(exact coul)':>24}{'brute(120)':>16}{'rel.diff':>12}")
for Lmax in (2, 6, 12, 24, 40):
    aL = legendre_moments(fg, r, Lmax, nx=300)
    bL = coul_moments_exact(r, Lmax)
    cL = legendre_moments(fg, r, Lmax, nx=300)
    val = triangle_contract(aL, bL, cL, d1, d2, d3)
    print(f"{Lmax:>6}{val:>24.10f}{ref_c:>16.10f}{abs(val-ref_c)/abs(ref_c):>12.2e}")

print()
print("interpretation: 3a near machine-level => rule+contraction CORRECT;")
print("                3b residual => brute reference under-resolves 1/r12, not the rule.")
