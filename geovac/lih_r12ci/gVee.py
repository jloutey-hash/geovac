"""LiH R12-CI Stage 4b (part 2c): ANALYTIC (RI-free) g_Vee = <Phi0|(F-Fbar)^2 V_ee|Phi0>.

  g_Vee = <F^2 V_ee> - 2 Fbar <F V_ee> + Fbar^2 <V_ee>
        = <F^2 V_ee> - 2 Fbar h_Vee - Fbar^2 E2         (h_Vee = Cov[F,V_ee], E2 = <V_ee>).
Only <F^2 V_ee> = <(sum_p f_p)^2 (sum_r coul_r)> is new -- a sum over 216 (p,q,r) triples of
three pairwise interactions on the 4-electron block density |Phi0|^2 = D_p(1,2) D_p(3,4)/4.

UNIFIED reducer.  For each triple, sum over the D_p AO-pair components on each block
(base density per electron), then integrate the two non-Coulomb electrons by message-passing:
  * a spectator (no f-edge) -> multiply by its base-density integral (INT P00=INT P11=1, INT P01=0);
  * an f-edge to a kept (Coulomb) electron -> DRESS that electron ( *= Psi^f or Psi^{f^2} );
  * an f-edge between two integrated electrons -> a scalar f-interaction;
  * two f-edges meeting at an integrated vertex whose BOTH neighbours are the Coulomb pair
    -> the non-separable TRIANGLE (validated in lih_r12ci_triangle_gate.py).
The Coulomb kernel on the kept pair is 1/r (no f on r), Yukawa e^{-gam r}/r (one f on r), or
e^{-2gam r}/r (f^2 on r).  All separable Coulomb/Yukawa are m=0 (axially symmetric one-electron
densities); only the triangle needs the general-m prolate-Neumann machinery.

Normalisation cross-checked analytically: this scheme reproduces E2 = <V_ee> exactly
(GATE E2 below), and subsumes the triangle gate.

Run from debug/:  python lih_r12ci_gVee_analytic.py
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from . import kernels as _KG                                   # USE_EXACT_NEUMANN keys the cache
from .kernels import (
    a, GAM, grid_int, dens, build_kernel, rho_cyl_f, zc_f, geo_f)
from .hVee import (
    neumann_potential, P00, P11, P01, rho_g, NXI, NETA, cov_FA_FB, make_kernel_f,
    make_kernel_coul, make_kernel_Y, E2_REF)
from .energy import R
# reuse the VALIDATED triangle machinery
from .triangle import (
    build_kernel_m, triangle_raw, MMAX, NG, _MODE_CACHE)

NG_ = NXI * NETA
Kf = build_kernel(GAM)
Kf2 = build_kernel(2 * GAM)

# --------------------------------------------------------------------------- #
# Yukawa (screened-Coulomb) m=0 potential of an axially-symmetric density:
#   Psi^Y(h;gam) = Psi^coul(h) - Psi^{smooth}(h),  smooth kernel (1-e^{-gam d})/d (no singularity)
# --------------------------------------------------------------------------- #
def build_kernel_smooth(gam, nphi=48):
    xp, wp = leggauss(nphi); phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = rho_cyl_f; z = zc_f; rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    K = np.zeros((NG_, NG_))
    for ph, w in zip(phi, wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * np.cos(ph), 1e-24))
        K += 2.0 * w * (1.0 - np.exp(-gam * d)) / d
    return K


_KS = {GAM: build_kernel_smooth(GAM), 2 * GAM: build_kernel_smooth(2 * GAM)}
_NEU_CACHE = {}


def psi_coul(h):
    key = (h.tobytes(), bool(_KG.USE_EXACT_NEUMANN))          # legacy / exact cached separately
    if key not in _NEU_CACHE:
        _NEU_CACHE[key] = neumann_potential(h.reshape(NXI, NETA)).reshape(-1)
    return _NEU_CACHE[key]


def psi_yuk(h, gam):
    return psi_coul(h) - a ** 3 * (_KS[gam] @ (geo_f * h))


def dress(K, h):
    return a ** 3 * (K @ (geo_f * h))


def I_coul(L, Rt):
    return grid_int(L * psi_coul(Rt))


def I_yuk(L, Rt, gam):
    return grid_int(L * psi_yuk(Rt, gam))


# --------------------------------------------------------------------------- #
# the triangle Fourier kernels (built once) + per-mode middle-density modes
# --------------------------------------------------------------------------- #
Fm = [build_kernel_m(GAM, m) for m in range(MMAX + 1)]


def triangle_value(mu, mu_name, left, right):
    v, _ = triangle_raw(Fm, mu, left, right, mu_name)
    return v


# --------------------------------------------------------------------------- #
# the unified per-triple reducer
# --------------------------------------------------------------------------- #
ELECS = (1, 2, 3, 4)
BLOCKU, BLOCKD = (1, 2), (3, 4)
PAIRS = [(1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)]
CINTRA = {frozenset((1, 2)), frozenset((3, 4))}
# D_p components: electron (lower index)->g, (higher index)->h, with sign kappa
COMPS = [(P00, 'P00', P11, 'P11', 1.0), (P11, 'P11', P00, 'P00', 1.0),
         (P01, 'P01', P01, 'P01', -2.0)]


def _reduce(base, bname, e, gg, rem_f, kernel_kind, gam_k):
    """message-pass integrate the non-kept electrons; return the triple's value (one kappa,kappa')."""
    D = dict(base); NM = dict(bname); scalar = 1.0
    kept = {e, gg}
    edges = [set(ed) for ed in rem_f]                       # f-edges NOT on the coul pair
    alive = set(ELECS)
    tri_mu = None; tri_name = None

    def dnb(u):
        s = set()
        for ed in edges:
            if u in ed:
                s |= (ed - {u})
        return s

    while True:
        nk = [u for u in alive if u not in kept]
        if not nk:
            break
        u = next((c for c in nk if len(dnb(c)) <= 1), None)
        if u is None:                                        # TRIANGLE: 2 distinct nbrs, both kept
            u = next(c for c in nk if dnb(c) == kept)
            tri_mu, tri_name = D[u], NM[u]
            edges = [ed for ed in edges if u not in ed]; alive.discard(u)
            continue
        inc = [ed for ed in edges if u in ed]
        nb = dnb(u)
        if len(nb) == 0:                                     # spectator
            scalar *= grid_int(D[u])
        else:                                                # dress the single neighbour w
            w = next(iter(nb)); mult = len(inc)
            D[w] = D[w] * dress(Kf if mult == 1 else Kf2, D[u])
            NM[w] = None                                     # w now dressed -> not a clean base name
        alive.discard(u); edges = [ed for ed in edges if u not in ed]

    if scalar == 0.0:
        return 0.0
    L, Rt = D[e], D[gg]
    if tri_mu is not None:                                   # kernel is coul (triangle -> 1/r)
        return scalar * triangle_value(tri_mu, tri_name, L, Rt)
    if kernel_kind == 'coul':
        return scalar * I_coul(L, Rt)
    return scalar * I_yuk(L, Rt, gam_k)


def eval_triple(p, q, r):
    """<f_p f_q coul_r> over |Phi0|^2 (=D_p D_p/4)."""
    e, gg = r
    rset = frozenset(r)
    fedges = [frozenset(p), frozenset(q)]
    n_on_r = sum(1 for ed in fedges if ed == rset)
    rem_f = [tuple(ed) for ed in fedges if ed != rset]
    kernel_kind = 'coul' if n_on_r == 0 else 'yuk'
    gam_k = None if n_on_r == 0 else (GAM if n_on_r == 1 else 2 * GAM)
    tot = 0.0
    for gU, gUn, hU, hUn, kU in COMPS:
        for gD, gDn, hD, hDn, kD in COMPS:
            base = {1: gU, 2: hU, 3: gD, 4: hD}
            bname = {1: gUn, 2: hUn, 3: gDn, 4: hDn}
            tot += kU * kD * _reduce(base, bname, e, gg, rem_f, kernel_kind, gam_k)
    return 0.25 * tot


def group_of(p, q, r):
    """(f-type SS/SI/II, coul-type CS/CI) for the 6-product MC breakdown."""
    def ptype(x):
        return 'S' if frozenset(x) in CINTRA else 'I'
    a_, b_ = ptype(p), ptype(q)
    ft = 'SS' if a_ == b_ == 'S' else ('II' if a_ == b_ == 'I' else 'SI')
    ct = 'CS' if frozenset(r) in CINTRA else 'CI'
    return (ft, ct)


def F2Vee():
    """<F^2 V_ee> = sum over 216 ordered (p,q,r) triples, grouped."""
    groups = {}
    for p in PAIRS:
        for q in PAIRS:
            for r in PAIRS:
                v = eval_triple(p, q, r)
                key = group_of(p, q, r)
                groups[key] = groups.get(key, 0.0) + v
    return groups


def eval_pair_fc(p, r):
    """<f_p coul_r> over |Phi0|^2 (one f-edge, one Coulomb) -- for grid-consistent <F V_ee>."""
    e, gg = r
    rset = frozenset(r); pf = frozenset(p)
    if pf == rset:
        kernel_kind, gam_k, rem_f = 'yuk', GAM, []
    else:
        kernel_kind, gam_k, rem_f = 'coul', None, [tuple(p)]
    tot = 0.0
    for gU, gUn, hU, hUn, kU in COMPS:
        for gD, gDn, hD, hDn, kD in COMPS:
            base = {1: gU, 2: hU, 3: gD, 4: hD}; bname = {1: gUn, 2: hUn, 3: gDn, 4: hDn}
            tot += kU * kD * _reduce(base, bname, e, gg, rem_f, kernel_kind, gam_k)
    return 0.25 * tot


def eval_coul_only(r):
    """<coul_r> over |Phi0|^2 (no f-edge) -- grid-consistent <V_ee> piece."""
    e, gg = r
    tot = 0.0
    for gU, gUn, hU, hUn, kU in COMPS:
        for gD, gDn, hD, hDn, kD in COMPS:
            base = {1: gU, 2: hU, 3: gD, 4: hD}; bname = {1: gUn, 2: hUn, 3: gDn, 4: hDn}
            tot += kU * kD * _reduce(base, bname, e, gg, [], 'coul', None)
    return 0.25 * tot


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI g_Vee (analytic, RI-free): <(F-Fbar)^2 V_ee>")
    print(f"  sigma2 grid {NXI}x{NETA}, MMAX={MMAX}, f=exp(-{GAM} r)")
    print("=" * 78)

    # ---- GATE Y: Yukawa potential vs the h_T closed form (isotropic 1s_A) ----
    from .hT import yukawa_pot_iso, rA_f, d_aa, ZA
    yk = psi_yuk(d_aa, GAM); ykc = yukawa_pot_iso(rA_f, ZA, GAM)
    msk = d_aa > 1e-6 * d_aa.max()
    print(f"\n  GATE Y -- psi_yuk(rho_aa) vs closed yukawa_pot_iso: rel "
          f"{np.abs(yk-ykc)[msk].mean()/np.abs(ykc[msk]).mean():.1e}")

    # ---- GATE E2: the reducer reproduces <V_ee> with NO f-edges ----
    E2_enum = 2 * eval_coul_only((1, 2)) + 4 * eval_coul_only((1, 3))
    print(f"  GATE E2 -- reducer <V_ee> = {E2_enum:.6f}  vs Stage-2 {E2_REF:.6f}"
          f"   (rel {abs(E2_enum-E2_REF)/E2_REF:.1e})")

    # ---- <F^2 V_ee> by the enumerator, grouped into the 6 MC products ----
    print("\n  --- <F^2 V_ee> 6-product breakdown (analytic vs banked MC) ---")
    MC = {('SS', 'CS'): 0.141, ('SS', 'CI'): 0.578, ('SI', 'CS'): 0.895,
          ('SI', 'CI'): 3.994, ('II', 'CS'): 1.524, ('II', 'CI'): 7.388}
    grp = F2Vee()
    labels = {('SS', 'CS'): '<S^2 C_S>', ('SS', 'CI'): '<S^2 C_I>', ('SI', 'CS'): '2<S I C_S>',
              ('SI', 'CI'): '2<S I C_I>', ('II', 'CS'): '<I^2 C_S>', ('II', 'CI'): '<I^2 C_I>'}
    F2V = 0.0
    for key in [('SS', 'CS'), ('SS', 'CI'), ('SI', 'CS'), ('SI', 'CI'), ('II', 'CS'), ('II', 'CI')]:
        val = grp.get(key, 0.0); F2V += val
        print(f"    {labels[key]:12s} = {val:+.5f}   (MC ~ {MC[key]:+.3f})")
    print(f"    {'<F^2 V_ee>':12s} = {F2V:+.5f}   (MC ~ +14.519)")

    # ---- assemble g_Vee ---- (two ways: mixed cov-vs-enum, and grid-consistent all-enum) ----
    Af = make_kernel_f(GAM); Kc = make_kernel_coul(); Ky = make_kernel_Y(GAM)
    h_Vee, info = cov_FA_FB(Af, Kc, Ky)
    Fbar = info['FbarA']; E2 = info['FbarB']
    g_Vee_mixed = F2V - 2 * Fbar * h_Vee - Fbar ** 2 * E2

    # grid-consistent: <F V_ee> and <V_ee> from the SAME enumerator machinery -> clean cancellation
    FVee = sum(eval_pair_fc(p, r) for p in PAIRS for r in PAIRS)
    hVee_enum = FVee - Fbar * E2_enum                      # = Cov[F,V_ee] on the grid
    g_Vee = F2V - 2 * Fbar * FVee + Fbar ** 2 * E2_enum    # = <(F-Fbar)^2 V_ee>, all grid-consistent
    print(f"\n  Fbar = {Fbar:.5f}   <F V_ee>_enum = {FVee:.5f}   <V_ee>_enum = {E2_enum:.5f}")
    print(f"  h_Vee: cov={h_Vee:+.5f}  enum={hVee_enum:+.5f}   (VMC +0.3082)")
    print(f"  g_Vee (mixed cov/enum)      = {g_Vee_mixed:+.6f}")
    print(f"  >>> g_Vee (grid-consistent) = {g_Vee:+.6f}   (VMC target +0.53605)")
