"""Be R12-CI, Part 1: the reference determinant Phi_0 (Be 1s^2 2s^2) and its energy.

Minimal Slater basis: 1s = e^{-z1 r}, 2s = Schmidt-orthogonalized (r e^{-z2 r}) against 1s.
Closed-shell RHF energy for two doubly-occupied s-orbitals a=1s, b=2s:
    E0 = 2 h_a + 2 h_b + J_aa + J_bb + 4 J_ab - 2 K_ab
All integrals are RADIAL (s-orbitals): computed by Gauss-Legendre radial quadrature.
Exponents (z1,z2) optimized to minimise E0 -> Phi_0 ~ HF (so the r12 lowering downstream
is CORRELATION, not reference error).  Validate vs the known Be HF limit (-14.573 Ha).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

Z = 4.0  # Be nucleus

# radial grid on [0, Rmax]
NR = 600
Rmax = 25.0
xg, wg = leggauss(NR)
rr = 0.5 * Rmax * (xg + 1.0)
wr = 0.5 * Rmax * wg


def norm_radial(R):
    """normalize a radial function so INT R^2 r^2 dr = 1."""
    n2 = np.sum(R * R * rr ** 2 * wr)
    return R / np.sqrt(n2)


def orbitals(z1, z2):
    """return normalized radial 1s, 2s (2s Schmidt-orthogonalized to 1s)."""
    s1 = norm_radial(np.exp(-z1 * rr))
    chi = rr * np.exp(-z2 * rr)                       # nodeless 2s primitive
    ov = np.sum(s1 * chi * rr ** 2 * wr)              # <1s|chi>
    s2 = norm_radial(chi - ov * s1)
    return s1, s2


def kinetic(R):
    """T = 1/2 INT (dR/dr)^2 r^2 dr  (s-orbital)."""
    dR = np.gradient(R, rr)
    return 0.5 * np.sum(dR * dR * rr ** 2 * wr)


def v_ne(R):
    """V_ne = -Z INT R^2 r dr."""
    return -Z * np.sum(R * R * rr * wr)


def slater_R0(Ra, Rb):
    """R^0[a,b] = INT INT Ra(r1)^2 Rb(r2)^2 (1/r>) r1^2 r2^2 dr1 dr2  (L=0 Coulomb, s-orbs).
       = Coulomb J between densities |a|^2 and |b|^2."""
    a2 = Ra * Ra * rr ** 2 * wr
    b2 = Rb * Rb * rr ** 2 * wr
    r_gt = np.maximum.outer(rr, rr)
    return a2 @ (1.0 / r_gt) @ b2


def slater_R0_exchange(Ra, Rb):
    """K exchange for s-orbitals = R^0 with the cross density Ra*Rb on each electron:
       K = INT INT [Ra(r1)Rb(r1)][Ra(r2)Rb(r2)] (1/r>) r1^2 r2^2 dr1 dr2."""
    ab1 = Ra * Rb * rr ** 2 * wr
    r_gt = np.maximum.outer(rr, rr)
    return ab1 @ (1.0 / r_gt) @ ab1


def E0(params):
    z1, z2 = params
    if z1 <= 0 or z2 <= 0:
        return 1e6
    a, b = orbitals(z1, z2)
    ha = kinetic(a) + v_ne(a)
    hb = kinetic(b) + v_ne(b)
    Jaa = slater_R0(a, a)
    Jbb = slater_R0(b, b)
    Jab = slater_R0(a, b)
    Kab = slater_R0_exchange(a, b)
    return 2 * ha + 2 * hb + Jaa + Jbb + 4 * Jab - 2 * Kab


if __name__ == "__main__":
    print("=" * 70)
    print("Be reference Phi_0 = 1s^2 2s^2 (minimal Slater, single-zeta each)")
    print("=" * 70)
    # optimize exponents
    res = minimize(E0, x0=[3.7, 0.7], method="Nelder-Mead",
                   options=dict(xatol=1e-4, fatol=1e-6))
    z1, z2 = res.x
    a, b = orbitals(z1, z2)
    ha = kinetic(a) + v_ne(a); hb = kinetic(b) + v_ne(b)
    Jaa = slater_R0(a, a); Jbb = slater_R0(b, b)
    Jab = slater_R0(a, b); Kab = slater_R0_exchange(a, b)
    Eref = 2 * ha + 2 * hb + Jaa + Jbb + 4 * Jab - 2 * Kab
    print(f"optimized exponents: z(1s)={z1:.4f}  z(2s)={z2:.4f}")
    print(f"  <1s|1s>={np.sum(a*a*rr**2*wr):.6f}  <2s|2s>={np.sum(b*b*rr**2*wr):.6f}"
          f"  <1s|2s>={np.sum(a*b*rr**2*wr):.2e}")
    print(f"  h_1s={ha:.5f}  h_2s={hb:.5f}")
    print(f"  J_1s1s={Jaa:.5f}  J_2s2s={Jbb:.5f}  J_1s2s={Jab:.5f}  K_1s2s={Kab:.5f}")
    print(f"\n  E0 (Phi_0, this minimal basis) = {Eref:.5f} Ha")
    print(f"  Be HF limit (reference)        = -14.5730 Ha")
    print(f"  Be exact (nonrel)              = -14.6674 Ha")
    print(f"  => minimal-basis HF gap to HF limit: {Eref-(-14.5730):+.4f} Ha")
    print(f"     correlation energy to capture (exact-thisPhi0): {(-14.6674)-Eref:+.4f} Ha")
    np.savez("debug/data/be_r12ci_ref.npz", z1=z1, z2=z2, a=a, b=b, rr=rr, wr=wr, Eref=Eref)
    print("\nsaved orbitals -> debug/data/be_r12ci_ref.npz")
