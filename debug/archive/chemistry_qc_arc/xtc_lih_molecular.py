"""
Molecular xTC angular-sparsity test on LiH (4 electrons).
=========================================================

Question (extends the atomic p-inclusive GO, sprint_xtc_pinclusive_memo.md):
The atomic run showed the xTC-contracted effective 2-body operator inherits
GeoVac's angular Gaunt sparsity EXACTLY (0 fill-in) -- but ONLY because
contracting a correlator leg over a single-center s reference orbital (l=0)
forces its leg multipole L'=0 (monopole), collapsing the four_Y vertex to a
single multipole (== Coulomb's support). A MOLECULE breaks that premise: the
reference occupied MO is a TWO-CENTER object (dense in the atomic basis about
any single center), so contracting a leg over it need NOT collapse to L'=0.

We build a minimal two-center LiH:
  Li (center A): core-s(2.7), val-s(0.65), p_{-1,0,+1}(0.65)
  H  (center B): s(1.0)
=> 6 spatial / 12 spin-orbitals, 4 electrons (Li 1s^2 2s + H 1s), C(12,4)=495 dets.

All two-center one-/two-electron integrals come from the VALIDATED TwoCenterLM
engine (two_center_grid_lm.py; s/p vs oracles ~1e-6). Complex spherical harmonics
=> every basis function has a definite m about the molecular z-axis (both centers
on z), so m-conservation is the exact angular selection rule (l is NOT a good
number for a two-center density -- the point of the test).

The geminal-TC Hamiltonian:
  * Hermitian TC 2-body kernel  w(r12) = 1/r + D  (Slater geminal, Kato cusp),
    two-center ERI via the same multipole assembly as Coulomb with the radial
    kernel replaced by the Legendre multipole of w (validated vs direct r12
    double-grid quadrature).
  * 3-body L3 (grid-moment reformulation of the atomic four_Y machinery: correlator
    expanded about the common origin A, orbital angular content taken as numeric
    grid moments; radial = scalar multipole of u' -- angular SELECTION exact, radial
    magnitudes model-level, exactly as the atomic run).

xTC = replace Coulomb by w + contract L3 against the LiH reference 1-RDM (the
occupied RHF sigma MOs). We measure fill-in (xTC-nonzero where plain-zero) and
1-norm, plain vs xTC, in the orthonormal AO basis -- exactly as the atomic run.

Reused verbatim from xtc_poc_li.py: lowdin, transform_1/2/3, spatial, spin,
make_dets, apply_ops, h_spin, asym_from_phys, build_H, build_H3, ground,
xtc_contract, op_metrics.
"""
import os, sys, json, time
import numpy as np
from math import pi

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
_ROOT = os.path.abspath(os.path.join(_HERE, '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from two_center_grid_lm import TwoCenterLM
from tc_threebody_collapse_angular import g_pair, gaunt
from xtc_poc_li import (lowdin, transform_1, transform_2, transform_3,
                        spatial, spin, make_dets, h_spin, asym_from_phys,
                        build_H, build_H3, ground, xtc_contract, op_metrics)

FOURPI = 4.0 * pi


# ----------------------------------------------------------------------------
# geminal kernels (same as ctf12 / xtc_poc_li)
# ----------------------------------------------------------------------------
def w_kernel(r, g):
    r = np.asarray(r, float)
    small = r < 1e-8
    rr = np.where(small, 1.0, r)
    out = -np.expm1(-g * rr) / rr + (g / 2) * np.exp(-g * rr) - 0.25 * np.exp(-2 * g * rr)
    return np.where(small, 1.5 * g - 0.25, out)

def uprime(r, g):    # u'(r) = 1/2 e^{-g r}
    return 0.5 * np.exp(-g * r)


# ----------------------------------------------------------------------------
# radial trapezoid weights on the engine's graded grid
# ----------------------------------------------------------------------------
def radial_weights(r):
    w = np.zeros_like(r)
    w[1:-1] = (r[2:] - r[:-2]) / 2.0
    w[0] = (r[1] - r[0]) / 2.0
    w[-1] = (r[-1] - r[-2]) / 2.0
    return w


def multipole_kernels(r, kernel, Lmax, nx=80):
    """FL[L][i,j] = (2L+1)/2 int_-1^1 kernel(r12(r_i,r_j,x)) P_L(x) dx  (nr x nr)."""
    R1, R2 = np.meshgrid(r, r, indexing='ij')
    xs, ws = np.polynomial.legendre.leggauss(nx)
    Pl = {L: np.polynomial.legendre.Legendre.basis(L) for L in range(Lmax + 1)}
    acc = {L: np.zeros_like(R1) for L in range(Lmax + 1)}
    for x, wx in zip(xs, ws):
        r12 = np.sqrt(np.maximum(R1 * R1 + R2 * R2 - 2 * R1 * R2 * x, 1e-30))
        f = kernel(r12)
        for L in range(Lmax + 1):
            acc[L] += wx * f * Pl[L](x)
    return {L: (2 * L + 1) / 2.0 * acc[L] for L in range(Lmax + 1)}


# ----------------------------------------------------------------------------
# One-/two-electron matrices in the two-center STO basis (chemist -> physicist)
# ----------------------------------------------------------------------------
def build_S_h1(eng, orbs, Z_A=3.0, Z_B=1.0):
    """S and h1 = kinetic - Z_A/r_A - Z_B/r_B (correct LiH nuclear charges)."""
    ns = len(orbs)
    S = np.zeros((ns, ns)); h1 = np.zeros((ns, ns))
    for a in range(ns):
        for b in range(a, ns):
            S[a, b] = S[b, a] = eng.overlap(orbs[a], orbs[b])
            T = eng.kinetic(orbs[a], orbs[b])
            Vne = -(Z_A * eng.coulomb_center(orbs[a], orbs[b], 'A')
                    + Z_B * eng.coulomb_center(orbs[a], orbs[b], 'B'))
            h1[a, b] = h1[b, a] = T + Vne
    return S, h1


def eri_kernel_chemist(eng, oi, oj, ok, ol, FL, Lmax):
    """(oi oj | ok ol) chemist, general radial kernel FL[L]."""
    A = eng._moment(oi, oj, "b"); B = eng._moment(ok, ol, "f")
    r = eng.r; wt = r * r * radial_weights(r)
    total = 0.0
    for L in range(Lmax + 1):
        pref = FOURPI / (2 * L + 1)
        for M in range(-L, L + 1):
            a, b = A[(L, M)], B[(L, M)]
            if np.max(np.abs(a)) < 1e-14 or np.max(np.abs(b)) < 1e-14:
                continue
            total += pref * np.real((a * wt) @ FL[L] @ (b * wt))
    return float(np.real(total))


def build_eri_phys(eng, orbs, FL, Lmax):
    """Physicist <ij|kl> = chemist (ik|jl); e1 pair (i,k), e2 pair (j,l)."""
    ns = len(orbs)
    eri = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        mi = orbs[i][2]
        for k in range(ns):
            mk = orbs[k][2]
            M1 = mk - mi                       # e1 density multipole M
            for j in range(ns):
                mj = orbs[j][2]
                for l in range(ns):
                    ml = orbs[l][2]
                    if (mi + mj) != (mk + ml):  # total-m conservation about z
                        continue
                    eri[i, j, k, l] = eri_kernel_chemist(eng, orbs[i], orbs[k],
                                                         orbs[j], orbs[l], FL, Lmax)
    return eri


# ----------------------------------------------------------------------------
# Two-center 3-body L3 (grid-moment reformulation; model radial = multipole of u')
#   V3[a,b,c,d,e,f]: vertex particle-1 (a,d); legs particle-2 (b,e) line L,
#   particle-3 (c,f) line L'.  Matches xtc_poc_li_pinclusive.build_V3 convention.
# ----------------------------------------------------------------------------
def build_V3(eng, orbs, uML, Lcorr, tol=1e-13):
    ns = len(orbs)
    r = eng.r; wt = r * r * radial_weights(r)
    pairs = [(i, j) for i in range(ns) for j in range(ns)]
    npair = len(pairs)
    pidx = {p: n for n, p in enumerate(pairs)}
    corr = [(L, M) for L in range(Lcorr + 1) for M in range(-L, L + 1)]

    # leg moments Mb_pair(L,M)(r) = int conj(chi_b) chi_e Y*_LM dOmega, per pair
    Mb = eng  # engine caches moments; grab "b" moment dicts per pair
    JL = {}    # JL[(L,M)] : (npair, nr)  J_L(r1;pair)=sum_r2 uML_L[r1,r2] Mb(L,M)(r2) r2^2 wr2
    Mb_arr = {}
    for (L, M) in corr:
        arr = np.zeros((npair, len(r)), dtype=complex)
        for n, (b, e) in enumerate(pairs):
            arr[n] = eng._moment(orbs[b], orbs[e], "b")[(L, M)]
        Mb_arr[(L, M)] = arr
        JL[(L, M)] = (arr * wt[None, :]) @ uML[L].T   # (npair, nr)

    # vertex moments Mf_pair(Lam,N)(r) = int conj(chi_a) chi_d Y_{Lam N} dOmega (f-moment)
    LamMax = 2 * Lcorr
    Mf_arr = {}
    for Lam in range(LamMax + 1):
        for N in range(-Lam, Lam + 1):
            arr = np.zeros((npair, len(r)), dtype=complex)
            for n, (a, d) in enumerate(pairs):
                arr[n] = eng._moment(orbs[a], orbs[d], "f")[(Lam, N)]
            Mf_arr[(Lam, N)] = arr

    V3 = np.zeros((ns, ns, ns, ns, ns, ns))
    # enumerate correlator channels (L,M,L',M') with nonzero g_pair -> Lambda
    channels = []
    for (L, M) in corr:
        for (Lp, Mp) in corr:
            N = M + Mp
            gp = {}
            for Lam in range(abs(L - Lp), L + Lp + 1):
                if Lam > LamMax or abs(N) > Lam:
                    continue
                c = g_pair(L, M, Lp, Mp, Lam)
                if c != 0.0:
                    gp[Lam] = c
            if gp:
                channels.append((L, M, Lp, Mp, N, gp))

    for pv, (a, d) in enumerate(pairs):
        for (L, M, Lp, Mp, N, gp) in channels:
            # vertex radial factor (nr,)
            vf = np.zeros(len(r), dtype=complex)
            for Lam, c in gp.items():
                vf += c * Mf_arr[(Lam, N)][pv]
            base = vf * wt                                  # r1^2 wr1 folded
            if np.max(np.abs(base)) < tol:
                continue
            jL = JL[(L, M)]; jLp = JL[(Lp, Mp)]
            pref = (-0.5) * (FOURPI / (2 * L + 1)) * (FOURPI / (2 * Lp + 1))
            # contrib[be,cf] = pref sum_r1 base(r1) jL[be,r1] jLp[cf,r1]
            contrib = pref * ((jL * base[None, :]) @ jLp.T)   # (npair,npair) complex
            cr = np.real(contrib)
            if np.max(np.abs(cr)) < tol:
                continue
            c2 = cr.reshape(ns, ns, ns, ns)                  # (b,e,c,f)
            V3[a, :, :, d, :, :] += c2.transpose(0, 2, 1, 3) # (b,c,e,f)
    return V3


# ----------------------------------------------------------------------------
# closed-shell RHF in the orthonormal AO basis, physicist eri_o[p,q,r,s]=<pq|rs>
#   e1 pair (p,r), e2 pair (q,s). F = h + sum_rs D_rs (<pr|qs> - 1/2 <pr|sq>)
#   D = 2 C_occ C_occ^T (spatial density).  Returns (C, occ_list, E_rhf, D).
# ----------------------------------------------------------------------------
def rhf(h1o, eri_o, n_occ, maxit=200, tol=1e-9):
    ns = h1o.shape[0]
    # core guess
    ev, C = np.linalg.eigh(h1o)
    E_old = 0.0
    for it in range(maxit):
        D = 2.0 * C[:, :n_occ] @ C[:, :n_occ].T
        # G[p,q] = sum_rs D_rs ( <pr|qs> - 1/2 <pr|sq> )
        J = np.einsum('rs,prqs->pq', D, eri_o)
        K = np.einsum('rs,prsq->pq', D, eri_o)
        F = h1o + J - 0.5 * K
        E = 0.5 * np.sum(D * (h1o + F))
        ev, C = np.linalg.eigh(F)
        if abs(E - E_old) < tol:
            break
        E_old = E
    occ = list(range(n_occ))
    return C, occ, float(E), D


# ----------------------------------------------------------------------------
# xTC contraction over a MOLECULAR reference: rotate V3 to MO basis, contract
# over occupied MO spin-orbitals (diagonal 1-RDM there), rotate v2 back to AO.
#   Returns v2_ao (nso^4 antisym-ready), v1_ao (nso^2), v0.
# ----------------------------------------------------------------------------
def xtc_contract_mo(V3o_ao, nso, C_spatial, occ_spatial):
    ns = V3o_ao.shape[0]
    # spin-orbital MO coefficient matrix U[ao_so, mo_so]
    U = np.zeros((nso, nso))
    for a in range(ns):
        for p in range(ns):
            U[2 * a, 2 * p] = C_spatial[a, p]
            U[2 * a + 1, 2 * p + 1] = C_spatial[a, p]
    # rotate V3 (spatial) to MO spatial
    V3mo = transform_3(V3o_ao, C_spatial.T)   # C.T @ ... = MO from AO (columns are MOs)
    # occupied MO spin-orbitals (closed shell: both spins of each occ spatial)
    occ_so = []
    for o in occ_spatial:
        occ_so += [2 * o, 2 * o + 1]
    v2_mo, v1_mo, v0 = xtc_contract(V3mo, nso, occ_so)
    # rotate v2_mo, v1_mo back to AO spin-orbitals: v_ao = U v_mo U^T (indices)
    v1_ao = U @ v1_mo @ U.T
    t = np.einsum('pi,ijkl->pjkl', U, v2_mo)
    t = np.einsum('qj,pjkl->pqkl', U, t)
    t = np.einsum('rk,pqkl->pqrl', U, t)
    v2_ao = np.einsum('sl,pqrl->pqrs', U, t)
    return v2_ao, v1_ao, v0, v2_mo, v1_mo


# ----------------------------------------------------------------------------
# sparsity metrics on a spatial physicist 2-body tensor
# ----------------------------------------------------------------------------
def support(tensor4, tol=1e-9):
    nz = np.abs(tensor4) > tol
    return set(map(tuple, np.argwhere(nz)))


def sparsity_report(base_tensor, xtc_tensor, tol=1e-9):
    sb = support(base_tensor, tol); sx = support(xtc_tensor, tol)
    fill_in = sx - sb
    ns = base_tensor.shape[0]
    return dict(
        base_nnz=len(sb), xtc_nnz=len(sx), total=ns ** 4,
        base_density=len(sb) / ns ** 4, xtc_density=len(sx) / ns ** 4,
        fill_in=len(fill_in),
        l1_base=float(np.sum(np.abs(base_tensor))),
        l1_xtc=float(np.sum(np.abs(xtc_tensor))),
        l1_ratio=float(np.sum(np.abs(xtc_tensor)) / (np.sum(np.abs(base_tensor)) + 1e-300)),
    )


# ============================================================================
# Assemble LiH and run the sparsity study
# ============================================================================
def assemble(R=3.0, gamma=1.0, zeta_core=2.7, zeta_val=0.65, zeta_H=1.0,
             nr=600, nu=28, nphi=28, rmax=45.0, Lmax=10, Lcorr=2,
             n_elec=4, include_Hp=False, want_fci=True):
    orbs = [
        (zeta_core, 0, 0, 'A'),   # Li core s
        (zeta_val, 0, 0, 'A'),    # Li val s
        (zeta_val, 1, -1, 'A'),   # Li p-1
        (zeta_val, 1, 0, 'A'),    # Li p0
        (zeta_val, 1, 1, 'A'),    # Li p+1
        (zeta_H, 0, 0, 'B'),      # H s
    ]
    if include_Hp:
        orbs.append((zeta_H, 1, 0, 'B'))   # H p0 (sigma)
    ns = len(orbs)
    eng = TwoCenterLM(R, nr=nr, nu=nu, nphi=nphi, rmax=rmax, Lmax=Lmax, real=False)

    S, h1 = build_S_h1(eng, orbs, Z_A=3.0, Z_B=1.0)
    FLc = multipole_kernels(eng.r, lambda rr: 1.0 / np.maximum(rr, 1e-30), Lmax)
    FLw = multipole_kernels(eng.r, lambda rr: w_kernel(rr, gamma), Lmax)
    uML = multipole_kernels(eng.r, lambda rr: uprime(rr, gamma), Lcorr)

    eri_coul = build_eri_phys(eng, orbs, FLc, Lmax)
    eri_w = build_eri_phys(eng, orbs, FLw, Lmax)
    V3 = build_V3(eng, orbs, uML, Lcorr)

    # Loewdin (full S; two-center couples l within m)
    X = lowdin(S)
    h1o = transform_1(h1, X)
    eri_coul_o = transform_2(eri_coul, X)
    eri_w_o = transform_2(eri_w, X)
    V3o = transform_3(V3, X)

    n_occ = n_elec // 2
    C, occ, E_rhf, D = rhf(h1o, eri_coul_o, n_occ)

    nso = 2 * ns
    out = dict(R=R, gamma=gamma, ns=ns, nso=nso, n_elec=n_elec, orbs=orbs,
               E_rhf=E_rhf, n_occ=n_occ, include_Hp=include_Hp)

    # xTC contraction over molecular reference (MO)
    v2_ao, v1_ao, v0, v2_mo, v1_mo = xtc_contract_mo(V3o, nso, C, occ)

    # effective 2-body tensors, SPATIAL physicist, for sparsity:
    #   base  = Coulomb eri_coul_o
    #   xtc   = w + (contracted L3 in spatial physicist form)
    # v2_ao is spin-orbital <pq||rs>-ready coefficient; fold to spatial physicist
    # by taking the alpha-alpha-alpha-alpha spin block WITHOUT antisym so it lines up
    # with eri_coul_o (which is spatial physicist). We compare the spatial physicist
    # 2-body operator supports directly.
    v2_spatial = spin_to_spatial_phys(v2_ao, nso, ns)
    xtc_2b = eri_w_o + v2_spatial

    rep_full = sparsity_report(eri_coul_o, xtc_2b)
    rep_w = sparsity_report(eri_coul_o, eri_w_o)          # w vs coulomb (same support?)
    rep_v2 = sparsity_report(eri_coul_o, v2_spatial)      # contracted-L3 alone vs coulomb

    out['sparsity_full'] = rep_full
    out['sparsity_w_vs_coul'] = rep_w
    out['sparsity_v2_vs_coul'] = rep_v2

    # ---- faithful spin-orbital (qubit-Hamiltonian) fill-in + 1-norm ----
    asym_coul = asym_from_phys(eri_coul_o, nso)
    asym_w = asym_from_phys(eri_w_o, nso)
    asym_x = asym_w + v2_ao
    out['sparsity_asym'] = asym_sparsity(asym_coul, asym_x, nso)
    out['sparsity_asym_w_vs_coul'] = asym_sparsity(asym_coul, asym_w, nso)

    # ---- mechanistic: reference-density multipole spectrum (L'=0 collapse?) ----
    out['ref_multipoles'] = ref_density_multipoles(eng, orbs, X, C, occ)

    # ---- control: AO-diagonal reference (atomic-style, aufbau on h1o diag) ----
    diag_e = np.diag(h1o)
    occ_ao = list(np.argsort(diag_e)[:n_occ])
    occ_so_ao = []
    for o in occ_ao:
        occ_so_ao += [2 * o, 2 * o + 1]
    v2_ao_diag, v1_ao_diag, v0_diag = xtc_contract(V3o, nso, occ_so_ao)
    asym_x_diag = asym_w + v2_ao_diag
    out['sparsity_asym_aodiag'] = asym_sparsity(asym_coul, asym_x_diag, nso)
    out['occ_ao_control'] = [int(x) for x in occ_ao]

    # m-block cross check: does anything cross total-m?
    out['m_violations'] = count_m_violations(v2_spatial, orbs)

    if want_fci:
        out.update(fci_gates(h1o, eri_coul_o, eri_w_o, V3o, nso, occ, C, n_elec))

    out['_arr'] = dict(eri_coul_o=eri_coul_o, eri_w_o=eri_w_o, V3o=V3o,
                       v2_spatial=v2_spatial, C=C, occ=occ, orbs=orbs, X=X,
                       h1o=h1o, nso=nso, ns=ns, uML=uML, eng_r=eng.r)
    return out


def spin_to_spatial_phys(v2_ao, nso, ns):
    """Extract a spatial physicist tensor from the spin-orbital coefficient v2_ao
    (coefficient of a_p^ a_q^ a_t a_s in the effective 2-body; xtc_contract's v2
    convention v2[p,q,s,t]). We want the SAME-SPIN, opposite-pair spatial content
    aligned with eri_phys[i,j,k,l]=<ij|kl>, e1=(i,k) e2=(j,l).
    v2[p,q,s,t] ~ <pq||st>-type; take the aa,bb blocks (p,q,s,t all alpha) and map
    spatial i=sp(p), j=sp(q), k=sp(s), l=sp(t)."""
    out = np.zeros((ns, ns, ns, ns))
    for p in range(0, nso, 2):        # alpha only, opposite-spin channel captured via q beta
        for q in range(1, nso, 2):    # q beta so no exchange antisymmetry collapses it
            for s in range(0, nso, 2):
                for t in range(1, nso, 2):
                    out[spatial(p), spatial(q), spatial(s), spatial(t)] += v2_ao[p, q, s, t]
    return out


def asym_sparsity(asym_base, asym_xtc, nso, tol=1e-9):
    """Fill-in + 1-norm on the spin-orbital antisymmetrized <pq||rs> tensors
    (the actual qubit-Hamiltonian 2-body object)."""
    sb = support(asym_base, tol); sx = support(asym_xtc, tol)
    return dict(base_nnz=len(sb), xtc_nnz=len(sx), total=nso ** 4,
                fill_in=len(sx - sb),
                l1_base=float(np.sum(np.abs(asym_base))),
                l1_xtc=float(np.sum(np.abs(asym_xtc))),
                l1_ratio=float(np.sum(np.abs(asym_xtc)) /
                               (np.sum(np.abs(asym_base)) + 1e-300)))


def ref_density_multipoles(eng, orbs, X, C, occ, Lmax=6):
    """Multipole spectrum (about center A) of the RHF reference density
    rho_ref = sum_{i occ} 2 |phi_i|^2 in the ORIGINAL AO basis, on the grid.
    phi_i (orthonormal-MO coeff C) in original AO: coeff_ao = X @ C[:, i].
    Returns |moment_L|^2 summed over M and the fraction of angular norm in L>0
    -- the mechanistic measure of whether the s-reference L'=0 collapse survives."""
    ns = len(orbs)
    # coefficient of each ORIGINAL AO in each occupied MO
    Cocc_ao = X @ C[:, occ]                    # (ns, n_occ)
    chis = [eng.chi(o) for o in orbs]          # complex grid arrays
    rho = np.zeros_like(chis[0], dtype=complex)
    for j, i in enumerate(occ):
        phi = sum(Cocc_ao[a, j] * chis[a] for a in range(ns))
        rho += 2.0 * np.conj(phi) * phi        # closed shell (2 per MO)
    # angular moments about A at each radius: int rho Y*_LM dOmega
    specL = {}
    for L in range(Lmax + 1):
        s = 0.0
        for M in range(-L, L + 1):
            mom = np.tensordot(rho, eng.Yb[(L, M)], axes=([1, 2], [0, 1]))  # (nr,) uses conj(Y)*wang
            # radial norm of this (L,M) moment
            s += float(np.real(np.sum(np.abs(mom) ** 2 * eng.r ** 2 * radial_weights(eng.r))))
        specL[L] = s
    tot = sum(specL.values())
    frac_L0 = specL[0] / tot if tot else 0.0
    return dict(specL={int(k): v for k, v in specL.items()},
                frac_L0=frac_L0, frac_Lpos=1.0 - frac_L0)


def count_m_violations(tensor4, orbs, tol=1e-9):
    ns = len(orbs)
    viol = 0
    for i in range(ns):
        for j in range(ns):
            for k in range(ns):
                for l in range(ns):
                    if abs(tensor4[i, j, k, l]) > tol:
                        if (orbs[i][2] + orbs[j][2]) != (orbs[k][2] + orbs[l][2]):
                            viol += 1
    return viol


def fci_gates(h1o, eri_coul_o, eri_w_o, V3o, nso, occ, C, n_elec):
    """FCI sanity: plain, xTC; geminal->0 handled by caller via gamma sweep."""
    dets, didx = make_dets(nso, n_elec)
    hso = h_spin(h1o, nso)
    asym_coul = asym_from_phys(eri_coul_o, nso)
    asym_w = asym_from_phys(eri_w_o, nso)
    H_plain = build_H(dets, didx, hso, asym_coul, nso)
    E_plain, _ = ground(H_plain, hermitian=True)
    # xTC energy: w + contracted-L3 over MO ref
    v2_ao, v1_ao, v0, _, _ = xtc_contract_mo(V3o, nso, C, occ)
    asym_x = asym_w + v2_ao
    hso_x = hso + v1_ao
    H_xTC = build_H(dets, didx, hso_x, asym_x, nso, v0=v0)
    E_xTC, imag_xTC = ground(H_xTC)
    return dict(ndet=len(dets), E_plain=E_plain, E_xTC=E_xTC, imag_xTC=imag_xTC,
                xtc_v0=v0)


if __name__ == '__main__':
    t0 = time.time()
    o = assemble(nr=600, want_fci=True)
    E_nn = 3.0 / o['R']
    print('ns=%d nso=%d ndet=%d' % (o['ns'], o['nso'], o.get('ndet', 0)))
    print('E_rhf(elec)=%.5f  E_plain(elec)=%.6f  +E_nn=%.4f -> E_tot~%.5f'
          % (o['E_rhf'], o['E_plain'], E_nn, o['E_plain'] + E_nn))
    print('E_xTC(elec)=%.6f imag=%.1e' % (o['E_xTC'], o['imag_xTC']))
    print('\n-- spatial 2-body support (Coulomb vs xTC) --')
    print(o['sparsity_full'])
    print('\n-- spin-orbital asym <pq||rs> (qubit Hamiltonian) --')
    print('xtc:', o['sparsity_asym'])
    print('w-only vs coul:', o['sparsity_asym_w_vs_coul'])
    print('AO-diagonal ref control:', o['sparsity_asym_aodiag'])
    print('\n-- mechanistic: RHF sigma-reference multipole spectrum --')
    print(o['ref_multipoles'])
    print('m_violations:', o['m_violations'])
    print('elapsed %.1fs' % (time.time() - t0))
