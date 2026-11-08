"""N=4 WALL -- Be-relevant 4-electron integral, reduced closed form vs brute force.

Confirms (with a number) the diagnostic verdict of r12ci_4e_wall_diagnostic.py /
r12ci_4e_wall_q2_reducibility.py: the genuinely-4-body bridging term of the variational
explicit-r12 CI is EXACT, RI-FREE, and REDUCIBLE -- not an RI/decidability wall.

The object (the only genuinely 4-body-connected term in <Phi|F H F|Phi>, F=sum f_pq
multiplicative -- kinetic gradients are pair-local and V_ne is one-body, so neither can
bridge two disjoint correlation edges; only the two-body Coulomb can):

    I = INT rho1(r1) rho2(r2) rho3(r3) rho4(r4)  f(r12) f(r34) (1/r13)  d3r1..d3r4      (*)

a CHAIN 2-1-3-4.  Be-relevant orbitals: electrons 1,3 (the BRIDGE) are 2p (Be's angular
correlating space); electrons 2,4 (the LEAVES) are 1s (the core).  This exercises the
bridge multipoles L in {0,2} (Q1: bridge terminates at L<=2*l_bridge=2).

Two independent evaluations of (*):
  BRUTE   -- the full 12-D integral by importance-sampled Monte-Carlo (no reduction used).
  REDUCED -- the closed-form path the diagnostic predicts:
             (1) each 1s LEAF integrates out into a 1D radial DRESSING of its bridge partner
                 Phi_leaf(r1) = INT |phi_1s(r2)|^2 f(r12) d3r2   (spherical average of f),
             (2) leaving a standard TWO-electron Slater-Condon Coulomb integral between the
                 two dressed 2p densities:  I = sum_{L=0,2} Theta^L * R^L[Drad,Drad].
             No resolution-of-identity; L-sum finite; radial = 1D leaf conv + 2D bridge.

If REDUCED == BRUTE (to MC precision) the 4-electron integral is confirmed exact & RI-free.
Pure validation: one integral, two methods; no energies.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

rng = np.random.default_rng(20260920)

# --------------------------------------------------------------------------- #
# system: Be-relevant Slater orbitals + a simple geminal
# --------------------------------------------------------------------------- #
Z1S = 3.70     # Be 1s Slater exponent (leaves = core)
Z2P = 1.00     # diffuse 2p (bridge = angular correlating space)
GAM = 0.50     # geminal  f(r) = exp(-GAM r)


def f_gem(r):
    return np.exp(-GAM * r)


# normalized densities (|phi|^2):
#   1s :  rho = (Z1S^3/pi) e^{-2 Z1S r}                      (isotropic)
#   2pz:  rho = (Z2P^5/pi) r^2 e^{-2 Z2P r} cos^2(theta)     (angular cos^2)
def rho_1s(r):
    return (Z1S ** 3 / np.pi) * np.exp(-2 * Z1S * r)


def rho_2pz(r, ct):
    return (Z2P ** 5 / np.pi) * r ** 2 * np.exp(-2 * Z2P * r) * ct ** 2


# --------------------------------------------------------------------------- #
# BRUTE : 12-D importance-sampled Monte-Carlo of (*)
#   sample ri from each density, estimator = < f(r12) f(r34) / r13 >
# --------------------------------------------------------------------------- #
def sample_1s(n):
    r = rng.gamma(3.0, 1.0 / (2 * Z1S), size=n)         # r^2 e^{-2Z r}
    ct = rng.uniform(-1, 1, size=n)                     # isotropic
    phi = rng.uniform(0, 2 * np.pi, size=n)
    st = np.sqrt(1 - ct * ct)
    return np.stack([r * st * np.cos(phi), r * st * np.sin(phi), r * ct], axis=1)


def sample_2pz(n):
    r = rng.gamma(5.0, 1.0 / (2 * Z2P), size=n)         # r^4 e^{-2Z r}
    u = rng.uniform(0, 1, size=n)
    ct = np.cbrt(2 * u - 1.0)                           # pdf ~ cos^2 : CDF^{-1}
    phi = rng.uniform(0, 2 * np.pi, size=n)
    st = np.sqrt(np.maximum(1 - ct * ct, 0.0))
    return np.stack([r * st * np.cos(phi), r * st * np.sin(phi), r * ct], axis=1)


def brute(n_tot, batch=4_000_000):
    acc = 0.0; acc2 = 0.0; ntot = 0
    while ntot < n_tot:
        n = min(batch, n_tot - ntot)
        r1 = sample_2pz(n); r2 = sample_1s(n); r3 = sample_2pz(n); r4 = sample_1s(n)
        d12 = np.linalg.norm(r1 - r2, axis=1)
        d34 = np.linalg.norm(r3 - r4, axis=1)
        d13 = np.linalg.norm(r1 - r3, axis=1)
        g = f_gem(d12) * f_gem(d34) / np.maximum(d13, 1e-12)
        acc += g.sum(); acc2 += (g * g).sum(); ntot += n
    mean = acc / ntot
    var = acc2 / ntot - mean * mean
    return mean, np.sqrt(var / ntot)


# --------------------------------------------------------------------------- #
# REDUCED : leaf dressing (1D) + 2-electron Slater-Condon Coulomb (deterministic)
# --------------------------------------------------------------------------- #
# radial grid (Gauss-Legendre mapped to [0, Rmax])
NR = 400
Rmax = 30.0
xg, wg = leggauss(NR)
rr = 0.5 * Rmax * (xg + 1.0)
wr = 0.5 * Rmax * wg

# angle grid for the spherical average of f (leaf monopole)
NX = 120
xx, wx = leggauss(NX)


def f_monopole(r1, r2):
    """f_0(r1,r2) = (1/2) INT_-1^1 f(|r1-r2|) dx  (spherical avg of f over the leaf angle)."""
    r1 = np.asarray(r1)[:, None]; r2 = np.asarray(r2)[None, :]        # (Nr1, Nr2)
    out = np.zeros((r1.shape[0], r2.shape[1]))
    for x, w in zip(xx, wx):
        r12 = np.sqrt(np.maximum(r1 * r1 + r2 * r2 - 2 * r1 * r2 * x, 1e-30))
        out += 0.5 * w * f_gem(r12)
    return out


def leaf_dressing(r_bridge):
    """Phi_leaf(r1) = INT |phi_1s(r2)|^2 f(r12) d3r2
       = INT [4 Z1S^3 r2^2 e^{-2Z1S r2}] f_0(r1,r2) dr2   (dOmega2 gives 4pi; rho_1s*4pi=4Z^3..)."""
    f0 = f_monopole(r_bridge, rr)                                    # (Nb, NR)
    radial_1s = 4 * Z1S ** 3 * rr ** 2 * np.exp(-2 * Z1S * rr)       # |phi_1s|^2 * 4pi * r2^2
    return f0 @ (radial_1s * wr)                                     # (Nb,)


def slater_RL(Drad, L):
    """R^L[D,D] = INT INT Drad(r1) Drad(r3) r<^L/r>^{L+1} r1^2 r3^2 dr1 dr3."""
    r_lt = np.minimum.outer(rr, rr)
    r_gt = np.maximum.outer(rr, rr)
    kernel = r_lt ** L / r_gt ** (L + 1)
    wD = Drad * rr ** 2 * wr
    return wD @ kernel @ wD


def theta_L(L, na=64):
    """Theta^L = INT INT (cos^2 th1)(cos^2 th3) P_L(u1.u3) dOmega1 dOmega3.

    Computed by DETERMINISTIC quadrature (no hand-derived constant): fix phi1=0 (x 2pi by
    azimuthal symmetry), integrate cos(th1), cos(th3), phi3 by Gauss-Legendre."""
    from numpy.polynomial.legendre import legval
    cf = np.zeros(L + 1); cf[L] = 1.0
    ct, wct = leggauss(na)                     # cos(theta1), cos(theta3)
    pp, wpp = leggauss(na)
    phi3 = np.pi * (pp + 1.0); wphi = np.pi * wpp
    tot = 0.0
    for c1, w1 in zip(ct, wct):
        s1 = np.sqrt(1 - c1 * c1)
        for c3, w3 in zip(ct, wct):
            s3 = np.sqrt(1 - c3 * c3)
            cosg = s1 * s3 * np.cos(phi3) + c1 * c3          # u1.u3
            integ = (c1 * c1) * (c3 * c3) * legval(cosg, cf)
            tot += w1 * w3 * np.sum(wphi * integ)
    return 2 * np.pi * tot                     # the phi1 integral (x 2pi)


def reduced():
    Phi_leaf = leaf_dressing(rr)                         # dressing on the bridge radial grid
    # dressed 2p RADIAL density (angular cos^2 pulled into Theta^L):
    #   rho_2pz = (Z2P^5/pi) r^2 e^{-2Z2P r} cos^2 ;  dress radial by Phi_leaf(r)
    Drad = (Z2P ** 5 / np.pi) * rr ** 2 * np.exp(-2 * Z2P * rr) * Phi_leaf
    total = 0.0
    parts = {}
    for L in (0, 2):
        RL = slater_RL(Drad, L)
        TL = theta_L(L)
        parts[L] = RL * TL
        total += RL * TL
    return total, parts


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #
print("=" * 78)
print("Be-relevant 4-electron bridging integral   I = <rho1 rho2 rho3 rho4  f12 f34 / r13>")
print(f"  bridge (e1,e3) = 2p (Z={Z2P}) ; leaves (e2,e4) = 1s (Z={Z1S}) ; f=exp(-{GAM} r)")
print("=" * 78)

red_total, red_parts = reduced()
print("\nREDUCED (deterministic closed form: leaf dressing + 2e Slater-Condon):")
for L in (0, 2):
    print(f"   L={L}:  Theta^L * R^L = {red_parts[L]:+.8e}")
print(f"   I_reduced = {red_total:+.8e}   (L-sum terminates at L=2, no RI)")

print("\nBRUTE (12-D importance-sampled Monte-Carlo, no reduction):")
for N in (20_000_000, 80_000_000):
    m, e = brute(N)
    rel = abs(m - red_total) / abs(red_total)
    print(f"   N={N:>11,}:  I_brute = {m:+.8e} +/- {e:.1e}   "
          f"rel.diff vs reduced = {rel:.2e}  ({(m-red_total)/e:+.1f} sigma)")

print("\n" + "-" * 78)
print("If I_brute -> I_reduced within MC error: the genuinely-4-body bridging integral is")
print("EXACT and RI-FREE -- it reduces to a 1D leaf dressing + a standard 2-electron")
print("Slater-Condon Coulomb integral, L-sum terminating at 2*l_bridge. The 'N=4 needs")
print("4-body operators / no <=3-body reduction' wall is SOFT for this (scalar Coulomb) term.")
print("Caveats unchanged: (i) exchange terms are the same chain class (permuted labels);")
print("(ii) the QUANTUM-ENCODING wall (4-body Pauli) is separate and real.")
