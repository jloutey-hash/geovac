"""Two-center LiH R12-CI energy -- assembling the Be {Phi0, F Phi0} 2x2 (be_r12ci_full.py)
in the prolate two-center geometry, with the validated two-center 4-body reduction
(lih_r12_4body_integral.py sigma + lih_r12_4body_pi_channel.py pi/delta) living in g.

Ansatz (Be-style, ansatz B of debug/lih_r12_build_plan.md):
    Psi = Phi0 + c (F - Fbar) Phi0 ,   G = (F - Fbar) Phi0 ,   F = sum_{i<j} f(r_ij).
2x2 generalized eigenproblem
    H = [[E0, h],[h, g]] ,  S = [[1, 0],[0, sigma2]] ,
    E0 = <Phi0|H|Phi0>,  Fbar = <Phi0|F|Phi0>,
    sigma2 = <Phi0|(F-Fbar)^2|Phi0>,  h = <Phi0|H(F-Fbar)|Phi0>,  g = <G|H|G>.
The genuinely-4-body content is the chain f12 (1/r13) f34 inside g's V_ee part -- exactly the
integral validated this session.

Minimal reference (the two-center analog of Be's minimal 1s^2 2s^2): a 4-electron determinant
    Phi0 = |1s_A^2 1s_B^2|      (Li core on focus A, H-side pair on focus B)
with 1s_A (exponent za, tight, Li-core-like) and 1s_B (exponent zb, diffuse). A PoC of the
machinery, not a spectroscopic LiH (as the Be R12-CI was a 19%-of-correlation PoC).

STAGE 1 (this file): the f-integral primitives + Fbar (well-conditioned), validated vs MC.
Clean decomposition: every AO 2-body f-integral <p q|f|r s> = INT rho_{pr}(1) rho_{qs}(2) f
reduces to prolate QUADRATURE by dressing the ISOTROPIC member of a pair (Psi_aa, Psi_bb are
radial), EXCEPT (ab|f|ab) where the two-center product rho_ab sits on both sides -> one clean
6-D importance-MC.

STAGE 2 (this file): E0 = <Phi0|H|Phi0>, the two-center closed-shell determinant energy
(well-conditioned), assembled from one-electron h = T + V_ne and the two-electron Coulomb
tensor, each integral cross-checked against a closed form or an independent route (all controls
< 5e-3, most < 1e-11).  E0 = -7.8878 Ha, above exact -8.070 (variational); the 182 mHa gap is
basis + correlation, recovered by Stages 3-4.  See stage2_E0().

The ill-conditioned sigma2 / h (quadrature) and g (the 4-body reduction) are the next stages.

Geminal f(r) = exp(-GAM r).  Run from root:  python debug/lih_r12ci_energy.py
"""
import os
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lih_r12_4body_integral import (   # noqa: E402
    R, a, CENTER_A, CENTER_B, f_gem, sample_1s)

rng = np.random.default_rng(20260921)

# ---- minimal two-center model --------------------------------------------- #
Z_A, Z_B = 3.0, 1.0        # Li, H nuclear charges (used later for E0)
ZA, ZB = 2.70, 1.00        # 1s_A (Li-core-like), 1s_B (valence) orbital exponents
GAM = 0.50                 # geminal f(r) = exp(-GAM r)
N_A = np.sqrt(ZA ** 3 / np.pi)   # 1s Slater norm: |1s|^2 = (z^3/pi) e^{-2 z r}, INT=1
N_B = np.sqrt(ZB ** 3 / np.pi)

# ---- prolate grid --------------------------------------------------------- #
NXI, NETA = 240, 96
_xg, _wxg = leggauss(NXI)
_xi_max = 1.0 + 44.0 / (2 * min(ZA, ZB) * a)
XI1D = 1.0 + 0.5 * (_xg + 1.0) * (_xi_max - 1.0)
WXI = 0.5 * (_xi_max - 1.0) * _wxg
_eg, _weg = leggauss(NETA)
ETA1D = _eg.copy(); WETA = _weg.copy()
XI, ETA = np.meshgrid(XI1D, ETA1D, indexing='ij')
RA = a * (XI + ETA); RB = a * (XI - ETA); JAC = (XI ** 2 - ETA ** 2)


def _grid_int(fg):
    """INT fg dtau over the prolate grid (fg on (XI,ETA); phi-independent)."""
    return 2 * np.pi * a ** 3 * np.einsum('i,j,ij,ij->', WXI, WETA, JAC, fg)


rho_A = (ZA ** 3 / np.pi) * np.exp(-2 * ZA * RA)      # |1s_A|^2
rho_B = (ZB ** 3 / np.pi) * np.exp(-2 * ZB * RB)      # |1s_B|^2
orb_A = N_A * np.exp(-ZA * RA)                        # 1s_A
orb_B = N_B * np.exp(-ZB * RB)                        # 1s_B
rho_ab = orb_A * orb_B                                # 2-center transition density
S_AB = _grid_int(rho_ab)                              # <1s_A|1s_B>

# ---- isotropic f-dressing (Be-style spherical average; exact for a 1s) ----- #
_NR = 400; _Rmax = 30.0
_xr, _wr = leggauss(_NR); _rr = 0.5 * _Rmax * (_xr + 1.0); _wrr = 0.5 * _Rmax * _wr
_NX = 160; _xx, _wx = leggauss(_NX)


def f_dress_iso(s_vals, zeta):
    """Psi(s) = INT |1s(zeta)|^2(r2) f(|s-r2|) d3r2  (radial; 1s isotropic about its center)."""
    s = np.asarray(s_vals)[:, None]; r2 = _rr[None, :]
    f0 = np.zeros((s.shape[0], _rr.shape[0]))
    for x, w in zip(_xx, _wx):
        r12 = np.sqrt(np.maximum(s * s + r2 * r2 - 2 * s * r2 * x, 1e-30))
        f0 += 0.5 * w * f_gem(r12)
    radial = 4 * zeta ** 3 * _rr ** 2 * np.exp(-2 * zeta * _rr)   # |1s|^2 * 4pi r2^2
    return f0 @ (radial * _wrr)


Psi_aa = f_dress_iso(RA.ravel(), ZA).reshape(RA.shape)   # f-dressing of |1s_A|^2 (fn of r_A)
Psi_bb = f_dress_iso(RB.ravel(), ZB).reshape(RB.shape)   # f-dressing of |1s_B|^2 (fn of r_B)

# ---- AO 2-body f-integrals ------------------------------------------------- #
# <p q|f|r s> = INT rho_{pr}(1) rho_{qs}(2) f ; pair densities aa,bb,ab. Dress the isotropic one.
V_aaaa = _grid_int(rho_A * Psi_aa)                        # (aa|f|aa)
V_bbbb = _grid_int(rho_B * Psi_bb)                        # (bb|f|bb)
V_aabb = _grid_int(rho_A * Psi_bb)                        # (aa|f|bb) = INT rho_A Psi_bb
V_aaab = _grid_int(rho_ab * Psi_aa)                       # (aa|f|ab) = INT rho_ab Psi_aa
V_bbab = _grid_int(rho_ab * Psi_bb)                       # (bb|f|ab) = INT rho_ab Psi_bb


def mc_ff(sample1, sample2, n=24_000_000, batch=3_000_000):
    """<f> for e1~sample1(n), e2~sample2(n); bounded f -> low variance, batch-means error."""
    means = []; ntot = 0
    while ntot < n:
        r1 = sample1(batch); r2 = sample2(batch)
        means.append(f_gem(np.linalg.norm(r1 - r2, axis=1)).mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


def s_A(n): return sample_1s(n, ZA, CENTER_A)
def s_B(n): return sample_1s(n, ZB, CENTER_B)


def mc_abab(n=24_000_000, batch=3_000_000):
    """(ab|f|ab) = INT rho_ab(1) rho_ab(2) f, by importance from q=|1s_A|^2:
    weight u = rho_ab/q = (N_B/N_A) e^{ZA r_A - ZB r_B} (finite variance since ZA>ZB)."""
    means = []; ntot = 0
    while ntot < n:
        r1 = s_A(batch); r2 = s_A(batch)
        rA1 = np.linalg.norm(r1 - CENTER_A, axis=1); rB1 = np.linalg.norm(r1 - CENTER_B, axis=1)
        rA2 = np.linalg.norm(r2 - CENTER_A, axis=1); rB2 = np.linalg.norm(r2 - CENTER_B, axis=1)
        u1 = (N_B / N_A) * np.exp(ZA * rA1 - ZB * rB1)
        u2 = (N_B / N_A) * np.exp(ZA * rA2 - ZB * rB2)
        d = np.linalg.norm(r1 - r2, axis=1)
        means.append((u1 * u2 * f_gem(d)).mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


# =========================================================================== #
# STAGE 2 -- E0 = <Phi0|H|Phi0>, the two-center closed-shell determinant energy.
#   Phi0 = |m0^2 m1^2|,  m0,m1 = Loewdin(1s_A, 1s_B)  (spans the ionic Li+(1s_A^2)
#   H-(1s_B^2) reference).  A determinant is invariant under a nonsingular mix of
#   its occupied orbitals, so <Phi0|H|Phi0> is the closed-shell RHF energy
#     E0 = 2(h00 + h11) + J00 + J11 + 4 J01 - 2 K01 + Z_A Z_B / R .
#   h = T + V_ne assembled from closed forms + prolate-grid attraction integrals
#   (each with a closed-form control).  Two-electron Coulomb: the 5 dressable
#   integrals by exact Hartree-dressing (closed-form 1s potential), (ab|ab) by the
#   VALIDATED prolate Neumann machinery (_prolate_neumann_coulomb) cross-checked vs MC.
#   Well-conditioned (no cancellation of large numbers), unlike the coming sigma2/h.
# =========================================================================== #
from lih_r12_4body_integral import (   # noqa: E402
    _prolate_neumann_coulomb, build_grid, rho_1s)

_PAIR = np.array([[0, 2], [2, 1]])     # AO index pair (0=a,1=b) -> W index (0=aa,1=bb,2=ab)
_G2_CACHE = {}


def _grid2(NXI=400, NETA=160, LMAX=44):
    key = (NXI, NETA, LMAX)
    if key not in _G2_CACHE:
        _G2_CACHE[key] = build_grid(NXI, NETA, LMAX)
    return _G2_CACHE[key]


def _gint(G, fg):
    """INT fg dtau over prolate grid G (phi-independent integrand)."""
    return 2 * np.pi * a ** 3 * np.einsum('i,j,ij,ij->', G['wxi'], G['weta'], G['JAC'], fg)


def _pt_charge_attr(Z):
    """<1s(Z)|1/r_other|1s(Z)> = potential of a unit 1s(Z) density at distance R (closed form)."""
    return (1.0 / R) * (1.0 - (1.0 + Z * R) * np.exp(-2.0 * Z * R))


def _hartree_1s(r, Z):
    """Coulomb potential of a normalized 1s(Z) density at radius r: (1/r)[1-(1+Zr)e^{-2Zr}]."""
    r = np.maximum(np.asarray(r, dtype=float), 1e-30)
    return (1.0 / r) * (1.0 - (1.0 + Z * r) * np.exp(-2.0 * Z * r))


def _mc_abab_coulomb(n=12_000_000, batch=3_000_000):
    """(ab|1/r12|ab)=INT rho_ab(1)rho_ab(2)/r12, importance from |1s_A|^2 (weight u=rho_ab/rho_A
    =(N_B/N_A)e^{ZA rA - ZB rB}; variance = INT rho_B rho_B /r12^2, finite).  1/r12 heavy tail
    -> batch-means error.  A cross-check on the Neumann value, not the value used in E0."""
    means = []; ntot = 0
    while ntot < n:
        r1 = sample_1s(batch, ZA, CENTER_A); r2 = sample_1s(batch, ZA, CENTER_A)
        rA1 = np.linalg.norm(r1 - CENTER_A, axis=1); rB1 = np.linalg.norm(r1 - CENTER_B, axis=1)
        rA2 = np.linalg.norm(r2 - CENTER_A, axis=1); rB2 = np.linalg.norm(r2 - CENTER_B, axis=1)
        u1 = (N_B / N_A) * np.exp(ZA * rA1 - ZB * rB1)
        u2 = (N_B / N_A) * np.exp(ZA * rA2 - ZB * rB2)
        d = np.maximum(np.linalg.norm(r1 - r2, axis=1), 1e-12)
        means.append((u1 * u2 / d).mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


def stage2_E0(NXI=400, NETA=160, LMAX=44, mc_check=True):
    """E0 = <Phi0|H|Phi0> for Phi0=|m0^2 m1^2|, m=Loewdin(1s_A,1s_B).  Returns (E0, info)."""
    G = _grid2(NXI, NETA, LMAX)
    rA_, rB_ = G['rA'], G['rB']
    rhoA = rho_1s(rA_, ZA); rhoB = rho_1s(rB_, ZB)
    orbA = N_A * np.exp(-ZA * rA_); orbB = N_B * np.exp(-ZB * rB_)
    rho_ab = orbA * orbB
    S = _gint(G, rho_ab)                                    # <1s_A|1s_B> on G2

    # -- attraction integrals <p|1/r_c|q> (Jacobian cancels the 1/r focus singularity) --
    def attr(dens, rc):
        return _gint(G, dens / np.maximum(rc, 1e-30))
    A_ArA = attr(rhoA, rA_);   A_BrB = attr(rhoB, rB_)      # ctrl ZA, ZB
    A_ArB = attr(rhoA, rB_);   A_BrA = attr(rhoB, rA_)      # ctrl pt(ZA), pt(ZB)
    A_abrA = attr(rho_ab, rA_); A_abrB = attr(rho_ab, rB_)  # hybrid (transition density)

    # -- kinetic: -1/2 del^2 1s(z) = (-z^2/2 + z/r) 1s(z); T_AB two ways (Hermiticity) --
    T_AA = 0.5 * ZA ** 2; T_BB = 0.5 * ZB ** 2
    T_AB_qB = -0.5 * ZB ** 2 * S + ZB * A_abrB
    T_AB_qA = -0.5 * ZA ** 2 * S + ZA * A_abrA
    T_AB = 0.5 * (T_AB_qB + T_AB_qA)

    # -- one-electron h = T + V_ne (V_ne = -Z_A/rA - Z_B/rB, nuclear charges) --
    h_AA = T_AA - Z_A * A_ArA - Z_B * A_ArB
    h_BB = T_BB - Z_A * A_BrA - Z_B * A_BrB
    h_AB = T_AB - Z_A * A_abrA - Z_B * A_abrB
    h_AO = np.array([[h_AA, h_AB], [h_AB, h_BB]])

    # -- two-electron Coulomb: 5 dressable via exact Hartree potential, (ab|ab) via Neumann --
    VA = _hartree_1s(rA_, ZA); VB = _hartree_1s(rB_, ZB)
    C_aaaa = _gint(G, rhoA * VA)                            # ctrl 5 ZA/8
    C_bbbb = _gint(G, rhoB * VB)                            # ctrl 5 ZB/8
    C_aabb = _gint(G, rhoB * VA)                            # INT rho_B V_A
    C_aabb_alt = _gint(G, rhoA * VB)                        # INT rho_A V_B (symmetry check)
    C_aaab = _gint(G, rho_ab * VA)                          # dress isotropic rho_A
    C_bbab = _gint(G, rho_ab * VB)                          # dress isotropic rho_B
    C_abab = _prolate_neumann_coulomb(rho_ab, rho_ab, G)[0]
    C_abab_mc = _mc_abab_coulomb() if mc_check else None
    W = np.array([[C_aaaa, C_aabb, C_aaab],
                  [C_aabb, C_bbbb, C_bbab],
                  [C_aaab, C_bbab, C_abab]])

    # -- Loewdin AO->MO, transform h and the 2-electron tensor --
    Smat = np.array([[1.0, S], [S, 1.0]])
    sval, svec = np.linalg.eigh(Smat)
    X = svec @ np.diag(1.0 / np.sqrt(sval)) @ svec.T        # m_p = sum_mu X[mu,p] AO_mu
    V4 = np.empty((2, 2, 2, 2))                             # V4[mu,nu,la,si] = (mu la|nu si)
    for mu in range(2):
        for nu in range(2):
            for la in range(2):
                for si in range(2):
                    V4[mu, nu, la, si] = W[_PAIR[mu, la], _PAIR[nu, si]]
    g = np.einsum('mp,nq,lr,ks,mnlk->pqrs', X, X, X, X, V4, optimize=True)   # g[p,q,r,s]=(pr|qs)
    h_MO = X.T @ h_AO @ X
    T_AO = np.array([[T_AA, T_AB], [T_AB, T_BB]])
    T_MO = X.T @ T_AO @ X
    Vne_MO = h_MO - T_MO                                    # h = T + V_ne
    T_exp = 2.0 * (T_MO[0, 0] + T_MO[1, 1])                # <Phi0|sum_i T_i|Phi0>
    Vne_exp = 2.0 * (Vne_MO[0, 0] + Vne_MO[1, 1])          # <Phi0|sum_i V_ne(i)|Phi0>

    J00, J11 = g[0, 0, 0, 0], g[1, 1, 1, 1]                 # (00|00), (11|11)
    J01 = g[0, 1, 0, 1]                                     # (00|11)
    K01 = g[0, 1, 1, 0]                                     # (01|10)
    h00, h11 = h_MO[0, 0], h_MO[1, 1]
    V_NN = Z_A * Z_B / R
    E1 = 2.0 * (h00 + h11)
    E2 = J00 + J11 + 4.0 * J01 - 2.0 * K01
    E0 = E1 + E2 + V_NN

    rel = lambda x, y: (abs(x - y) / abs(y)) if y else abs(x - y)
    info = dict(
        S=S, A_ArA=A_ArA, e_ArA=rel(A_ArA, ZA), A_BrB=A_BrB, e_BrB=rel(A_BrB, ZB),
        A_ArB=A_ArB, pt_A=_pt_charge_attr(ZA), e_ArB=rel(A_ArB, _pt_charge_attr(ZA)),
        A_BrA=A_BrA, pt_B=_pt_charge_attr(ZB), e_BrA=rel(A_BrA, _pt_charge_attr(ZB)),
        T_AB_qB=T_AB_qB, T_AB_qA=T_AB_qA, e_TAB=rel(T_AB_qB, T_AB_qA),
        C_aaaa=C_aaaa, e_aaaa=rel(C_aaaa, 5 * ZA / 8), C_bbbb=C_bbbb, e_bbbb=rel(C_bbbb, 5 * ZB / 8),
        C_aabb=C_aabb, C_aabb_alt=C_aabb_alt, e_aabb=rel(C_aabb, C_aabb_alt),
        C_abab=C_abab, C_abab_mc=C_abab_mc,
        h00=h00, h11=h11, J00=J00, J11=J11, J01=J01, K01=K01, E1=E1, E2=E2, V_NN=V_NN,
        T_exp=T_exp, Vne_exp=Vne_exp)
    return E0, info


if __name__ == "__main__":
    print(f"[model] Z_A={Z_A} Z_B={Z_B} R={R}  za={ZA} zb={ZB}  <a|b>={S_AB:.6f}")
    print(f"[grid ] norm(rho_A)={_grid_int(rho_A):.6f}  norm(rho_B)={_grid_int(rho_B):.6f} (exact 1)")

    print("\n--- AO 2-body f-integrals (f=exp(-GAM r)): quadrature vs MC ---")
    checks = [("(aa|f|aa)", V_aaaa, mc_ff(s_A, s_A)),
              ("(bb|f|bb)", V_bbbb, mc_ff(s_B, s_B)),
              ("(aa|f|bb)", V_aabb, mc_ff(s_A, s_B))]
    for nm, q, (m, e) in checks:
        print(f"  {nm}: quad={q:.6f}  MC={m:.6f}+/-{e:.1e}  rel={abs(q - m) / abs(m):.2e}")
    print(f"  (aa|f|ab): quad={V_aaab:.6f}   (bb|f|ab): quad={V_bbab:.6f}   [1 iso dressing each]")
    Vabab, eab = mc_abab()
    print(f"  (ab|f|ab): importance-MC={Vabab:.6f}+/-{eab:.1e}   [2-center on both sides]")

    # ---- AO f-tensor over pair densities (0=aa,1=bb,2=ab), Loewdin, Fbar -----
    W = np.array([[V_aaaa, V_aabb, V_aaab],
                  [V_aabb, V_bbbb, V_bbab],
                  [V_aaab, V_bbab, Vabab]])
    dens = {('a', 'a'): 0, ('b', 'b'): 1, ('a', 'b'): 2, ('b', 'a'): 2}
    AOs = ['a', 'b']

    def Vf(mu, nu, la, si):     # <mu nu|f|la si> = W[pair(mu,la), pair(nu,si)]
        return W[dens[(mu, la)], dens[(nu, si)]]

    Smat = np.array([[1.0, S_AB], [S_AB, 1.0]])
    sval, svec = np.linalg.eigh(Smat)
    Xlow = svec @ np.diag(1 / np.sqrt(sval)) @ svec.T     # MO_p = sum_mu Xlow[mu,p] AO_mu

    def Vf_MO(p, q, r, s):
        tot = 0.0
        for mi, mu in enumerate(AOs):
            for ni, nu in enumerate(AOs):
                for li, la in enumerate(AOs):
                    for si, sg in enumerate(AOs):
                        tot += (Xlow[mi, p] * Xlow[ni, q] * Xlow[li, r] * Xlow[si, s]
                                * Vf(mu, nu, la, sg))
        return tot

    Jf = lambda p, q: Vf_MO(p, q, p, q)
    Kf = lambda p, q: Vf_MO(p, q, q, p)
    Fbar = Jf(0, 0) + Jf(1, 1) + 4 * Jf(0, 1) - 2 * Kf(0, 1)
    print(f"\n--- Fbar over the Loewdin determinant |m0^2 m1^2| ---")
    print(f"  J^f_00={Jf(0,0):.5f}  J^f_11={Jf(1,1):.5f}  J^f_01={Jf(0,1):.5f}  K^f_01={Kf(0,1):.5f}")
    print(f"  Fbar = <Phi0| sum_ij f_ij |Phi0> = {Fbar:.6f}")

    print("\n" + "-" * 76)
    print("STAGE 1 done: AO f-integral primitives (quad==MC on the 3 isotropic ones) + Fbar.")

    print("\n" + "=" * 76)
    print("STAGE 2 -- E0 = <Phi0|H|Phi0>, two-center closed-shell determinant energy")
    print("  Phi0 = |m0^2 m1^2|, m=Loewdin(1s_A,1s_B)  [ionic Li+(1s_A^2) H-(1s_B^2) ref]")
    print("=" * 76)
    E0, info = stage2_E0(mc_check=True)

    print("\n  controls (each AO integral vs its closed form / independent route):")
    print(f"    S_AB (G2)   = {info['S']:.6f}  vs Stage-1 grid {S_AB:.6f}   (rel {abs(info['S']-S_AB)/abs(S_AB):.1e})")
    print(f"    <A|1/rA|A>  = {info['A_ArA']:.6f}  vs ZA={ZA:.4f}            (rel {info['e_ArA']:.1e})")
    print(f"    <B|1/rB|B>  = {info['A_BrB']:.6f}  vs ZB={ZB:.4f}            (rel {info['e_BrB']:.1e})")
    print(f"    <A|1/rB|A>  = {info['A_ArB']:.6f}  vs pt(ZA)={info['pt_A']:.6f}     (rel {info['e_ArB']:.1e})")
    print(f"    <B|1/rA|B>  = {info['A_BrA']:.6f}  vs pt(ZB)={info['pt_B']:.6f}     (rel {info['e_BrA']:.1e})")
    print(f"    T_AB        : expand@B {info['T_AB_qB']:.6f} = expand@A {info['T_AB_qA']:.6f}  "
          f"(Hermiticity rel {info['e_TAB']:.1e})")
    print(f"    (aa|aa)     = {info['C_aaaa']:.6f}  vs 5ZA/8={5*ZA/8:.6f}     (rel {info['e_aaaa']:.1e})")
    print(f"    (bb|bb)     = {info['C_bbbb']:.6f}  vs 5ZB/8={5*ZB/8:.6f}     (rel {info['e_bbbb']:.1e})")
    print(f"    (aa|bb)     = {info['C_aabb']:.6f}  vs Hartree-dress(rho_A V_B)={info['C_aabb_alt']:.6f}"
          f"  (rel {info['e_aabb']:.1e})")
    if info['C_abab_mc'] is not None:
        mcv, mce = info['C_abab_mc']
        print(f"    (ab|ab)     = {info['C_abab']:.6f}  [Neumann]  vs importance-MC {mcv:.6f}+/-{mce:.1e}"
              f"  (rel {abs(info['C_abab']-mcv)/abs(mcv):.1e})")

    print(f"\n  energy pieces (analytic targets for the VMC cross-check):")
    print(f"    <T>   = {info['T_exp']:+.6f}   <V_ne> = {info['Vne_exp']:+.6f}   "
          f"<V_ee> = {info['E2']:+.6f}   V_NN = {info['V_NN']:+.6f}")
    print(f"\n  one-electron (MO):  h00 = {info['h00']:+.6f}   h11 = {info['h11']:+.6f}")
    print(f"  two-electron (MO):  J00={info['J00']:.5f}  J11={info['J11']:.5f}  "
          f"J01={info['J01']:.5f}  K01={info['K01']:.5f}")
    print(f"  E1 = 2(h00+h11) = {info['E1']:+.6f}   E2 = {info['E2']:+.6f}   V_NN = {info['V_NN']:+.6f}")
    print(f"\n  >>> E0 = <Phi0|H|Phi0> = {E0:+.6f} Ha")
    ok = E0 > -8.070
    print(f"  variational check: E0 {'>' if ok else '<= (BUG!!)'} LiH exact -8.070   "
          f"(gap to exact = {E0 - (-8.070):+.4f} Ha; this is what Stages 3-4 recover)")

    print("\n" + "-" * 76)
    print("STAGE 1+2 done.  E0 (well-conditioned) is assembled and control-validated.")
    print("NEXT: sigma2 = <F^2>-Fbar^2 and h = <Phi0|HF|Phi0>-Fbar*E0 -- the ILL-CONDITIONED")
    print("  pieces (small residuals of large numbers -> quadrature/near-exact, no cancellation,")
    print("  per the Be lesson; do with fresh care).  Then g = <G|H|G> with the 4-body chain")
    print("  f12(1/r13)f34 via the validated sigma+pi/delta reduction; then the 2x2 -> E_R12.")
