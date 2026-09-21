"""All-electron prolate FCI for diatomics — the variational-CORE test.

The frozen-core prolate LiH (increment 6, CHANGELOG v5.15.6) still drifts because
the Li 1s^2 core is FROZEN (static screening, can't polarize/contract with R). The
verdict: the cure is a variational CORE = an all-electron CI where ALL electrons,
core included, are active and R-adaptive.

This module does exactly that, reusing the validated grid machinery in
geovac/prolate_scf.py:
  1. Generate M one-electron eigen-MOs of h = T - Z_A/r_A - Z_B/r_B at the PHYSICAL
     charges (LiH: Z_A=3, Z_B=1) via get_orbital_on_grid.  These are mutually
     orthogonal (same Hermitian h) and R-adaptive; h is DIAGONAL = their E_elec.
  2. ERIs (pq|rs) over the MOs via compute_vee_integral (elliptic-K azimuthal avg).
  3. A determinant FCI (Slater-Condon) over the 2M spin-orbitals, N electrons.
     V_ee enters PAIRWISE (2-body) -> this sidesteps the N=4 explicit-r12 4-body
     wall (which is specific to Hylleraas explicit-r12, not to CI).

No frozen core, no PK, no explicit r12.  The exact test of whether a variational
core cures the R_eq drift.

Run from root:  python debug/prolate_allelectron_fci.py [test]
"""
import os
import sys
import time
from itertools import combinations

import numpy as np
from scipy.special import ellipk, ellipe

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac.prolate_scf import get_orbital_on_grid, compute_vee_integral   # noqa: E402


# ===========================================================================
# Generalized-m (pi/delta) ERI machinery.
# compute_vee_integral (prolate_scf) uses only the mu=0 azimuthally-averaged
# Coulomb kernel 8 pi K/sqrt(s).  For m != 0 orbitals the ERI needs the
# azimuthal-transfer-mu kernel.  With orbital Phi_p = psi_p(xi,eta) e^{i m_p phi}
# (psi normalized so INT psi^2 J dxi deta * 2pi = 1), the chemist ERI is
#   (pq|rs) = delta(m_p+m_r, m_q+m_s) * SUM_{1,2} [psi_p psi_q J w](1)
#             K_mu[1,2] [psi_r psi_s J w](2),   mu = m_p - m_q,
#   K_mu[1,2] = 2 pi F_|mu|(a,b),  a = rho1^2+rho2^2+dz^2, b = 2 rho1 rho2,
#   F_0 = 4 K(k)/sqrt(a+b),  F_1 = (4/b)[a K/sqrt(a+b) - sqrt(a+b) E],
#   F_{mu+1} = [2 mu (a/b) F_mu - (mu-1/2) F_{mu-1}]/(mu+1/2)   (toroidal recur.)
#   (k^2 = 2b/(a+b) is scipy's parameter m).  mu=0 reproduces 8 pi K/sqrt(s).
# ===========================================================================
def _azimuthal_kernels(orb, mu_max=2):
    """Precompute the ng x ng azimuthal-transfer kernels K_mu, mu=0..mu_max."""
    R = orb['R']
    XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
    rho = (R / 2) * np.sqrt(np.maximum((XI ** 2 - 1) * (1 - ETA ** 2), 0.0))
    z = (R / 2) * XI * ETA
    rf = rho.ravel(); zf = z.ravel()
    A = rf[:, None] ** 2 + rf[None, :] ** 2 + (zf[:, None] - zf[None, :]) ** 2
    B = 2.0 * rf[:, None] * rf[None, :]
    apb = A + B
    on_axis = B < 1e-14
    m_ell = np.clip(2.0 * B / np.maximum(apb, 1e-300), 0.0, 1.0 - 1e-15)
    Kk = ellipk(m_ell); Ek = ellipe(m_ell)
    sq = np.sqrt(apb)
    F = [4.0 * Kk / sq]                                   # F_0
    Bsafe = np.where(on_axis, 1.0, B)
    F1 = (4.0 / Bsafe) * (A * Kk / sq - sq * Ek)
    F1 = np.where(on_axis, 0.0, F1)
    F.append(F1)
    for mu in range(1, mu_max):
        Fn = (2.0 * mu * (A / Bsafe) * F[mu] - (mu - 0.5) * F[mu - 1]) / (mu + 0.5)
        Fn = np.where(on_axis, 0.0, Fn)
        F.append(Fn)
    F[0] = np.where(on_axis, 2.0 * np.pi / sq, F[0])      # F_0 on-axis limit
    return {mu: 2.0 * np.pi * F[mu] for mu in range(mu_max + 1)}


def _density_weighted(orb_p, orb_q):
    """[psi_p psi_q J w] flattened, for the ERI contraction."""
    R = orb_p['R']
    XI, ETA = np.meshgrid(orb_p['xi'], orb_p['eta'], indexing='ij')
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    WX, WE = np.meshgrid(orb_p['w_xi'], orb_p['w_eta'], indexing='ij')
    return (orb_p['psi'] * orb_q['psi'] * J * WX * WE).ravel()


def vee_m(orb_p, orb_q, orb_r, orb_s, m_p, m_q, m_r, m_s, Kmats):
    """chemist (pq|rs) with azimuthal quantum numbers; 0 unless m_p+m_r=m_q+m_s."""
    if m_p + m_r != m_q + m_s:
        return 0.0
    mu = abs(m_p - m_q)
    if mu not in Kmats:
        raise ValueError(f"kernel mu={mu} not precomputed")
    rpq = _density_weighted(orb_p, orb_q)
    rrs = _density_weighted(orb_r, orb_s)
    return float(rpq @ (Kmats[mu] @ rrs))


# ---------------------------------------------------------------------------
# MO generation + one/two-electron integrals in an orthonormal MO basis
# ---------------------------------------------------------------------------
def build_mo_integrals(R, M, Z_A, Z_B, m_channel=0, N_xi_solve=6000,
                       N_grid=48, xi_max=14.0, verbose=False):
    """M sigma (m=0) eigen-MOs of h at (Z_A,Z_B,R); return orthonormal (h1, eri).

    Returns (h1[M,M], eri[M,M,M,M] chemist (pq|rs), eps[M], Smax_off).
    Small grid non-orthogonality is removed by a Loewdin transform; h1 is built
    as the symmetrized 1/2(eps_p+eps_q) S_pq and rotated to the orthonormal basis.
    """
    orbs = []
    for na in range(M):
        o = get_orbital_on_grid(R=R, Z_A=Z_A, Z_B=Z_B, n_angular=na, m=m_channel,
                                N_xi_solve=N_xi_solve, N_xi_grid=N_grid,
                                N_eta_grid=N_grid, xi_max_grid=xi_max)
        o['R'] = R
        orbs.append(o)
    eps = np.array([o['E_elec'] for o in orbs])

    # grid overlap S_pq
    xi, eta = orbs[0]['xi'], orbs[0]['eta']
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    WX, WE = np.meshgrid(orbs[0]['w_xi'], orbs[0]['w_eta'], indexing='ij')
    wvol = J * WX * WE * 2 * np.pi
    P = np.array([o['psi'] for o in orbs])                     # [M, nxi, neta]
    S = np.einsum('pij,qij,ij->pq', P, P, wvol)
    Smax_off = np.max(np.abs(S - np.diag(np.diag(S))))

    # symmetrized one-electron matrix in the (non-orthonormal) MO basis
    H1 = 0.5 * (eps[:, None] + eps[None, :]) * S

    # ERIs (pq|rs)_chem = INT INT psi_p psi_q (1) 1/r12 psi_r psi_s (2)
    # compute_vee_integral(a,b,c,d) forms rho_ac=psi_a psi_c, rho_bd=psi_b psi_d
    # -> call (p, r, q, s) to get (pq|rs)_chem.
    eri = np.zeros((M, M, M, M))
    for p in range(M):
        for q in range(p, M):
            for r in range(M):
                for s in range(r, M):
                    val = compute_vee_integral(orbs[p], orbs[r], orbs[q], orbs[s])
                    for (a, b) in {(p, q), (q, p)}:
                        for (c, d) in {(r, s), (s, r)}:
                            eri[a, b, c, d] = val
    # Loewdin orthonormalization X = S^{-1/2}
    sval, svec = np.linalg.eigh(S)
    X = svec @ np.diag(1.0 / np.sqrt(np.clip(sval, 1e-12, None))) @ svec.T
    h1_o = X.T @ H1 @ X
    eri_o = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, eri, optimize=True)
    if verbose:
        print(f"    MO eps: {np.round(eps, 4)}  Smax_off={Smax_off:.1e}")
    return h1_o, eri_o, eps, Smax_off


# ---------------------------------------------------------------------------
# Multi-exponent basis: orbitals at several (Z_eff_A, Z_eff_B) scales, so the
# 1-particle space spans both the tight core (Z~3) and the diffuse valence (Z~1).
# Orbitals from different generating potentials are NOT co-eigenstates, so h1
# needs the true kinetic matrix, obtained from  T phi_q = eps_q phi_q + A(gen_q) phi_q
# (A = attraction matrix, a MULTIPLICATIVE operator computable on the grid).
# ---------------------------------------------------------------------------
def _attraction_matrix(orbs, R, Z_A, Z_B):
    """A[p,q] = <phi_p| Z_A/r_A + Z_B/r_B |phi_q>  (positive; V_ne = -A).
    Uses the ANALYTIC identity (Z_A/r_A + Z_B/r_B)*J = (R/2)^2[Z_A(xi-eta)+Z_B(xi+eta)]
    (the 1/r focus singularity cancels the Jacobian), so it is grid-robust even when
    the quadrature samples the focus (xi->1, eta->-1) -- essential for a graded grid."""
    xi, eta = orbs[0]['xi'], orbs[0]['eta']
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    hR2 = (R / 2.0) ** 2
    attrJ = hR2 * (Z_A * (XI - ETA) + Z_B * (XI + ETA))   # = (Z_A/r_A + Z_B/r_B) * J
    WX, WE = np.meshgrid(orbs[0]['w_xi'], orbs[0]['w_eta'], indexing='ij')
    wvol = attrJ * WX * WE * 2 * np.pi
    P = np.array([o['psi'] for o in orbs])
    return np.einsum('pij,qij,ij->pq', P, P, wvol)


def _overlap_matrix(orbs, R):
    xi, eta = orbs[0]['xi'], orbs[0]['eta']
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    WX, WE = np.meshgrid(orbs[0]['w_xi'], orbs[0]['w_eta'], indexing='ij')
    wvol = J * WX * WE * 2 * np.pi
    P = np.array([o['psi'] for o in orbs])
    return np.einsum('pij,qij,ij->pq', P, P, wvol)


def build_mo_integrals_multiexp(R, scales, Z_A_phys, Z_B_phys, n_ang_max=1,
                                N_xi_solve=6000, N_grid=48, xi_max=14.0,
                                verbose=False):
    """Multi-exponent orthonormal MO integrals.

    scales : list of (Z_eff_A, Z_eff_B) generating-charge pairs.  For each scale
             and each n_angular in 0..n_ang_max, one sigma orbital is generated.
    Returns (h1_o[M,M], eri_o chemist, M, Smax_off, cond_S).
    """
    orbs, gen, eps = [], [], []
    for (za, zb) in scales:
        for na in range(n_ang_max + 1):
            o = get_orbital_on_grid(R=R, Z_A=za, Z_B=zb, n_angular=na, m=0,
                                    N_xi_solve=N_xi_solve, N_xi_grid=N_grid,
                                    N_eta_grid=N_grid, xi_max_grid=xi_max)
            o['R'] = R
            orbs.append(o); gen.append((za, zb)); eps.append(o['E_elec'])
    M = len(orbs)
    eps = np.array(eps)
    S = _overlap_matrix(orbs, R)

    # attraction matrices: physical + one per distinct generating charge set
    A_phys = _attraction_matrix(orbs, R, Z_A_phys, Z_B_phys)
    Agen_cache = {}
    for (za, zb) in set(gen):
        Agen_cache[(za, zb)] = _attraction_matrix(orbs, R, za, zb)

    # kinetic: <q|T|p> = eps_p S_qp + <q|A(gen_p)|p>  (symmetrize p<->q)
    T = np.zeros((M, M))
    for p in range(M):
        Ap = Agen_cache[gen[p]]
        for q in range(M):
            Aq = Agen_cache[gen[q]]
            tqp = eps[p] * S[q, p] + Ap[q, p]
            tpq = eps[q] * S[p, q] + Aq[p, q]
            T[p, q] = 0.5 * (tqp + tpq)
    H1 = T - A_phys                                   # h = T + V_ne(phys) = T - A_phys

    # ERIs (pq|rs)_chem = vee(p, r, q, s)
    eri = np.zeros((M, M, M, M))
    for p in range(M):
        for q in range(p, M):
            for r in range(M):
                for s in range(r, M):
                    val = compute_vee_integral(orbs[p], orbs[r], orbs[q], orbs[s])
                    for (a, b) in {(p, q), (q, p)}:
                        for (c, d) in {(r, s), (s, r)}:
                            eri[a, b, c, d] = val

    sval, svec = np.linalg.eigh(S)
    cond_S = sval[-1] / max(sval[0], 1e-30)
    X = svec @ np.diag(1.0 / np.sqrt(np.clip(sval, 1e-10, None))) @ svec.T
    h1_o = X.T @ H1 @ X
    eri_o = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, eri, optimize=True)
    if verbose:
        print(f"    M={M} eps={np.round(eps,3)} cond(S)={cond_S:.1e} "
              f"Smax_off={np.max(np.abs(S-np.diag(np.diag(S)))):.1e}")
    return h1_o, eri_o, M, cond_S


# ---------------------------------------------------------------------------
# Full builder with pi (m != 0) orbitals.  One-electron S/h/V_ne are BLOCK
# DIAGONAL in m (cross-m phi-integral = 0); ERIs couple m-blocks via vee_m.
# ---------------------------------------------------------------------------
def build_mo_integrals_full(R, specs, Z_A_phys, Z_B_phys,
                            N_xi_solve=6000, N_grid=44, xi_max=15.0, verbose=False):
    """specs : list of (Z_eff_A, Z_eff_B, n_angular, m).  For m != 0 the (xi,eta)
    part is the |m| solution.  Returns (h1_o, eri_o, M, m_list, cond_S)."""
    orbs, gen, eps, mlist = [], [], [], []
    for (za, zb, na, mm) in specs:
        o = get_orbital_on_grid(R=R, Z_A=za, Z_B=zb, n_angular=na, m=abs(mm),
                                N_xi_solve=N_xi_solve, N_xi_grid=N_grid,
                                N_eta_grid=N_grid, xi_max_grid=xi_max)
        o['R'] = R
        orbs.append(o); gen.append((za, zb, abs(mm))); eps.append(o['E_elec']); mlist.append(mm)
    M = len(orbs)
    eps = np.array(eps); mlist = np.array(mlist)
    same_m = (mlist[:, None] == mlist[None, :])

    S = _overlap_matrix(orbs, R) * same_m
    A_phys = _attraction_matrix(orbs, R, Z_A_phys, Z_B_phys) * same_m
    Agen = {}
    for g in set(gen):
        Agen[g] = _attraction_matrix(orbs, R, g[0], g[1])
    T = np.zeros((M, M))
    for p in range(M):
        Ap = Agen[gen[p]]
        for q in range(M):
            if not same_m[p, q]:
                continue
            Aq = Agen[gen[q]]
            T[p, q] = 0.5 * ((eps[p] * S[q, p] + Ap[q, p]) + (eps[q] * S[p, q] + Aq[p, q]))
    H1 = T - A_phys

    mu_max = int(np.max(np.abs(mlist[:, None] - mlist[None, :])))
    Kmats = _azimuthal_kernels(orbs[0], mu_max=mu_max)
    eri = np.zeros((M, M, M, M))
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    if mlist[p] + mlist[r] != mlist[q] + mlist[s]:
                        continue
                    eri[p, q, r, s] = vee_m(orbs[p], orbs[q], orbs[r], orbs[s],
                                            mlist[p], mlist[q], mlist[r], mlist[s], Kmats)

    sval, svec = np.linalg.eigh(S)
    cond_S = sval[-1] / max(sval[0], 1e-30)
    X = svec @ np.diag(1.0 / np.sqrt(np.clip(sval, 1e-10, None))) @ svec.T
    h1_o = X.T @ H1 @ X
    eri_o = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, eri, optimize=True)
    if verbose:
        print(f"    M={M} mu_max={mu_max} cond(S)={cond_S:.1e}")
    return h1_o, eri_o, M, mlist, cond_S


# ---------------------------------------------------------------------------
# Spin-orbital determinant FCI (Slater-Condon), Sz=0 block
# ---------------------------------------------------------------------------
def _phys(eri, i, j, k, l):
    """physicist <ij|kl> (spatial) = chemist (ik|jl)."""
    return eri[i, k, j, l]


def _dets(M, na, nb):
    """all determinants with na alpha + nb beta spatial orbitals; spin-orbital
    index = 2*spatial + spin (0=alpha,1=beta). Returns list of sorted tuples."""
    out = []
    for a in combinations(range(M), na):
        for b in combinations(range(M), nb):
            so = sorted([2 * p for p in a] + [2 * p + 1 for p in b])
            out.append(tuple(so))
    return out


def _spatial(so):
    return so // 2, so % 2   # (spatial, spin)


def _h_so(h1, i, j):
    pi, si = _spatial(i); pj, sj = _spatial(j)
    return h1[pi, pj] if si == sj else 0.0


def _eri_so(eri, i, j, k, l):
    """spin-orbital physicist <ij|kl> = d(si,sk) d(sj,sl) (spatial phys)."""
    pi, si = _spatial(i); pj, sj = _spatial(j)
    pk, sk = _spatial(k); pl, sl = _spatial(l)
    if si == sk and sj == sl:
        return _phys(eri, pi, pj, pk, pl)
    return 0.0


def _antisym(eri, i, j, k, l):
    return _eri_so(eri, i, j, k, l) - _eri_so(eri, i, j, l, k)


def _diff(D1, D2):
    """spin-orbitals in D1 not in D2, in D2 not in D1 (as sorted lists)."""
    s1, s2 = set(D1), set(D2)
    return sorted(s1 - s2), sorted(s2 - s1)


def _phase(D, orbs_removed):
    """sign from annihilating the listed occupied spin-orbitals in D (in order)."""
    D = list(D)
    sign = 1
    for o in orbs_removed:
        idx = D.index(o)
        sign *= (-1) ** idx
        D.pop(idx)
    return sign


def _matel(h1, eri, D1, D2):
    """<D1|H|D2> by Slater-Condon."""
    if D1 == D2:
        e = sum(_h_so(h1, i, i) for i in D1)
        for a in range(len(D1)):
            for b in range(a + 1, len(D1)):
                e += _antisym(eri, D1[a], D1[b], D1[a], D1[b])
        return e
    o1, o2 = _diff(D1, D2)      # o1 in D1 not D2 ; o2 in D2 not D1
    ndiff = len(o1)
    if ndiff == 1:
        i, a = o1[0], o2[0]
        common = sorted(set(D1) & set(D2))
        val = _h_so(h1, i, a)
        for m in common:
            val += _antisym(eri, i, m, a, m)
        # phase: bring i and a to matching positions
        ph = _phase(D1, [i]) * _phase(D2, [a])
        return ph * val
    if ndiff == 2:
        i, j = o1
        a, b = o2
        val = _antisym(eri, i, j, a, b)
        ph = _phase(D1, [i, j]) * _phase(D2, [a, b])
        return ph * val
    return 0.0


def fci_energy(h1, eri, M, nelec, n_states=1):
    """Ground-state FCI energy in the Sz=0 block (nelec even, na=nb=nelec/2)."""
    na = nb = nelec // 2
    dets = _dets(M, na, nb)
    nd = len(dets)
    H = np.zeros((nd, nd))
    for a in range(nd):
        for b in range(a, nd):
            v = _matel(h1, eri, dets[a], dets[b])
            H[a, b] = H[b, a] = v
    evals = np.linalg.eigvalsh(H)
    return (evals[0], nd) if n_states == 1 else (evals[:n_states], nd)


# ---------------------------------------------------------------------------
# validation: 2-electron H2 (known) — validates MO-gen + integrals + FCI
# ---------------------------------------------------------------------------
def validate_h2():
    R = 1.4
    print(f"VALIDATION: prolate all-electron FCI on H2 (Z=1,1), R={R}")
    print("  known: E_HF ~ -1.128, full CI/exact ~ -1.174 Ha (D_e ~ 0.17)")
    for M in (2, 3, 4, 6):
        t = time.time()
        h1, eri, eps, soff = build_mo_integrals(R, M, 1.0, 1.0, N_xi_solve=5000,
                                                N_grid=44, xi_max=13.0)
        E_elec, nd = fci_energy(h1, eri, M, nelec=2)
        E_tot = E_elec + 1.0 / R
        print(f"  M={M}  n_det={nd}  E_tot={E_tot:.5f}  Smax_off={soff:.1e}  "
              f"[{time.time()-t:.0f}s]", flush=True)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        validate_h2()
    else:
        validate_h2()
