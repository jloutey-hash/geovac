"""
p-inclusive xTC PoC on Li (1s^2 2s, doublet) in the Coulomb-Sturmian s+p basis.
================================================================================

Extends the validated s-only engine (xtc_poc_li.py) to a MINIMAL s+p basis
(max_n=2 => 1s, 2s, 2p_-1, 2p_0, 2p_1  == 2 s-orbitals + one p-shell = 5 spatial,
10 spin-orbitals, C(10,3)=120 dets), to answer the one question the s-only run
could NOT: does the xTC-contracted effective TWO-body operator inherit the ~5%
angular Gaunt sparsity (Track 1's silver lining), i.e. does the 3-body -> 2-body
contraction PRESERVE the angular block structure or FILL IN the zero blocks?

Angular machinery is Track 1's VALIDATED four-harmonic W (four_Y) reused directly
(tc_threebody_collapse_angular.py, checked 1.9e-15 vs quadrature).  All FCI /
contraction plumbing is reused verbatim from xtc_poc_li.py.

Conventions (IMPORTANT):
  * All angular factors use the SAME correct complex-Y Gaunt convention as Track 1.
    The framework's production `_ck_coefficient` (SturmianCI) uses a q-sign
    convention that silently DROPS physical m-changing 2-body Coulomb multipoles
    (e.g. <2p+1 2p-1|1/r12|2p0 2p0> = -0.0342 -> returned 0).  We therefore build
    BOTH Coulomb and L3 with the correct gaunt so the plain-vs-xTC comparison is
    physically complete and internally consistent.  Our Coulomb matches the
    framework EXACTLY on every m-conserving block (validated).
  * Coulomb 2-body: electron-1 pair sees Y*_{LM}, electron-2 pair sees Y_{LM}
    (from 1/r12 = sum_L (4pi/(2L+1)) Y*_{LM}(1) Y_{LM}(2)).
  * L3 correlator lines: the vertex electron sees Y_{LM} (un-conjugated, via
    four_Y), each leg electron sees Y*_{LM} (conjugated, via Track 1's G_leg).
"""
import os, sys
import numpy as np
from math import pi

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
_ROOT = os.path.abspath(os.path.join(_HERE, '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# --- reuse the validated s-only FCI / contraction plumbing --------------------
from xtc_poc_li import (make_grid, lowdin, transform_1, transform_2, transform_3,
                        spatial, spin, make_dets, h_spin, asym_from_phys,
                        build_H, build_H3, xtc_contract, ground, op_metrics,
                        w_kernel)
# --- reuse Track 1's validated four-harmonic angular machinery -----------------
from tc_threebody_collapse_angular import four_Y, gaunt
# --- framework radial pieces (for validation + shared-k Sturmian radial) -------
from geovac.sturmian_solver import hydrogenic_radial, _slater_rk

FOURPI = 4.0 * pi


# ----------------------------------------------------------------------------
# Basis: (n,l,m) up to max_n.  max_n=2 -> [1s, 2s, 2p-1, 2p0, 2p1].
# ----------------------------------------------------------------------------
def build_states(max_n=2):
    st = []
    for n in range(1, max_n + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                st.append((n, l, m))
    return st


# ----------------------------------------------------------------------------
# Angular factors (correct complex-Y Gaunt convention; consistent with four_Y)
# ----------------------------------------------------------------------------
def gA(l, m, L, M, lp, mp):
    """int Y*_{lm} Y*_{LM} Y_{lp mp} dOmega   (= Track 1 G_leg; conjugated correlator)."""
    return (-1) ** (m + M) * gaunt(l, -m, L, -M, lp, mp)

def gB(l, m, L, M, lp, mp):
    """int Y*_{lm} Y_{LM} Y_{lp mp} dOmega   (un-conjugated correlator)."""
    return (-1) ** (m) * gaunt(l, -m, L, M, lp, mp)


# ----------------------------------------------------------------------------
# Radial: shared-k Coulomb-Sturmian R_{nl}(r) = hydrogenic_radial(r,n,l,Z=n*k)
#   (Z_eff = n*k gives common decay exp(-k r) for all n; the shared-p0 basis.)
# ----------------------------------------------------------------------------
def radial_table(states, r, k):
    return [hydrogenic_radial(r, n, l, n * k) for (n, l, m) in states]


def build_overlap(states, R, W):
    """S[a,b] = int R_a R_b r^2 dr, (l,m)-diagonal."""
    ns = len(states)
    S = np.zeros((ns, ns))
    for a in range(ns):
        for b in range(ns):
            if states[a][1] == states[b][1] and states[a][2] == states[b][2]:
                S[a, b] = np.sum(R[a] * R[b] * W)
    return S


def build_h1(states, S, k, Z):
    """Analytic Sturmian one-body (framework formula): h1 = -k^2/2 S + diag(k^2 - Z k / n)."""
    ns = len(states)
    h1 = -k ** 2 / 2.0 * S.copy()
    for a in range(ns):
        n = states[a][0]
        h1[a, a] += k ** 2 - Z * k / n
    return h1


# ----------------------------------------------------------------------------
# Radial multipole kernels M_L[r1,r2]
# ----------------------------------------------------------------------------
def coulomb_ML(r, Lmax):
    R1, R2 = np.meshgrid(r, r, indexing='ij')
    rlo = np.minimum(R1, R2); rhi = np.maximum(R1, R2)
    rhi = np.where(rhi < 1e-30, 1e-30, rhi)
    out = {}
    for L in range(Lmax + 1):
        out[L] = rlo ** L / rhi ** (L + 1)
    return out


def scalar_ML(r, kernel, Lmax, nx=200):
    """Legendre multipole M_L[r1,r2] = (2L+1)/2 int_-1^1 kernel(r12) P_L(x) dx of a
    scalar function kernel(r12)."""
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
# Two-body scalar-kernel ERI  <ab|f|cd>  (physicist; e1 pair (a,c), e2 pair (b,d))
#   = sum_L (4pi/(2L+1)) gA(a,L,M,c) gB(b,L,M,d) (dens_ac . M_L . dens_bd),  M=mc-ma
# ----------------------------------------------------------------------------
def build_eri_scalar(states, R, W, ML, Lmax):
    ns = len(states)
    dens = {}
    for a in range(ns):
        for c in range(ns):
            dens[(a, c)] = R[a] * R[c] * W
    # precompute radial R^L(a,c,b,d) lazily via matrix products
    eri = np.zeros((ns, ns, ns, ns))
    # angular tables
    for a in range(ns):
        la, ma = states[a][1], states[a][2]
        for c in range(ns):
            lc, mc = states[c][1], states[c][2]
            M = mc - ma
            # e1 angular per L
            aL = {}
            for L in range(Lmax + 1):
                v = gA(la, ma, L, M, lc, mc)
                if v != 0.0:
                    aL[L] = v
            if not aL:
                continue
            dac = dens[(a, c)]
            # radial vectors r1: y_L(r2) = M_L . dac  (contract electron-1 side)
            yL = {L: ML[L] @ dac for L in aL}
            for b in range(ns):
                lb, mb = states[b][1], states[b][2]
                for d in range(ns):
                    ld, md = states[d][1], states[d][2]
                    if ma + mb != mc + md:
                        continue
                    val = 0.0
                    for L, av in aL.items():
                        bv = gB(lb, mb, L, M, ld, md)
                        if bv == 0.0:
                            continue
                        rad = float(dens[(b, d)] @ yL[L])
                        val += (FOURPI / (2 * L + 1)) * av * bv * rad
                    if abs(val) > 1e-14:
                        eri[a, b, c, d] = val
    return eri


# ----------------------------------------------------------------------------
# Three-body L3 tensor  V3[a,b,c, d,e,f]  (vertex particle-1 (a,d);
#   legs particle-2 (b,e) line L, particle-3 (c,f) line L')
#
#   V3 = -1/2 sum_{L,M,L',M'} four_Y(a;L,M,L',M';d) G_leg(b;L,M;e) G_leg(c;L',M';f)
#            * (4pi/(2L+1))(4pi/(2L'+1)) * RadT
#   RadT = sum_r1 rho_ad(r1) * J_L(r1;b,e) * J_L'(r1;c,f)
#   J_L(r1;b,e) = sum_r2 uL_L[r1,r2] dens_be(r2),  uL_L = multipole of u'(r12)
#   rho_ad(r1) = R_a R_d r1^2 W
#
# Angular support is EXACTLY Track 1's four_Y (the required control).  Radial uses
# the scalar multipole of u' = 1/2 e^{-g r12} for each correlator line (a scalar-
# harmonic model of grad_i u; magnitudes are model-level, the angular SELECTION
# and hence the sparsity verdict are exact).
# ----------------------------------------------------------------------------
def build_V3(states, R, W, uML, Lmax):
    ns = len(states)
    npair = ns * ns
    pairs = [(i, j) for i in range(ns) for j in range(ns)]
    pidx = {p: n for n, p in enumerate(pairs)}

    dens_mat = np.array([R[i] * R[j] * W for (i, j) in pairs])          # (npair, Ng)
    rho_mat = dens_mat                                                   # rho_ad same form
    # J_L[L] : (npair, Ng)  J_L(r1;b,e) = uML[L] . dens_be
    JL = {L: (uML[L] @ dens_mat.T).T for L in range(Lmax + 1)}           # (npair,Ng)

    # angular precompute
    # G_leg array per (L,M): vector over pairs (leg = (l_b m_b)->(l_e m_e))
    def Gleg_arr(L, M):
        v = np.zeros(npair)
        for n, (b, e) in enumerate(pairs):
            lb, mb = states[b][1], states[b][2]
            le, me = states[e][1], states[e][2]
            v[n] = gA(lb, mb, L, M, le, me)
        return v
    Gleg_cache = {}
    for L in range(Lmax + 1):
        for M in range(-L, L + 1):
            Gleg_cache[(L, M)] = Gleg_arr(L, M)

    # four_Y per vertex pair (a,d)
    FY = {}
    for (a, d) in pairs:
        la, ma = states[a][1], states[a][2]
        ld, md = states[d][1], states[d][2]
        entries = []
        for L in range(Lmax + 1):
            for M in range(-L, L + 1):
                for Lp in range(Lmax + 1):
                    for Mp in range(-Lp, Lp + 1):
                        v, _ = four_Y(la, ma, L, M, Lp, Mp, ld, md)
                        if v != 0.0:
                            entries.append((L, M, Lp, Mp, v))
        if entries:
            FY[(a, d)] = entries

    V3 = np.zeros((ns, ns, ns, ns, ns, ns))
    for (a, d), entries in FY.items():
        ia = pidx  # unused
        rho = rho_mat[pidx[(a, d)]]                                      # (Ng,)
        for (L, M, Lp, Mp, fyv) in entries:
            g2 = Gleg_cache[(L, M)]                                      # (npair,) leg2
            g3 = Gleg_cache[(Lp, Mp)]                                    # (npair,) leg3
            if not (np.any(g2) and np.any(g3)):
                continue
            pref = fyv * (FOURPI / (2 * L + 1)) * (FOURPI / (2 * Lp + 1))
            base = (rho[None, :] * JL[L])                               # (npair,Ng) rows=be
            RadT = base @ JL[Lp].T                                       # (npair,npair) [be,cf]
            contrib = (-0.5) * pref * (g2[:, None] * g3[None, :]) * RadT  # (be, cf)
            # scatter: contrib[be,cf] -> V3[a,b,c,d,e,f]
            c2 = contrib.reshape(ns, ns, ns, ns)                        # (b,e,c,f)
            V3[a, :, :, d, :, :] += c2.transpose(0, 2, 1, 3)            # (b,c,e,f)
    return V3


# ----------------------------------------------------------------------------
# Angular density metric of a spatial 2-body tensor: fraction / support / 1-norm
# ----------------------------------------------------------------------------
def angular_density(tensor4, tol=1e-9):
    ns = tensor4.shape[0]
    nz = np.abs(tensor4) > tol
    support = set(map(tuple, np.argwhere(nz)))
    return dict(nnz=int(nz.sum()), total=ns ** 4,
                density=float(nz.sum()) / ns ** 4,
                l1=float(np.sum(np.abs(tensor4))),
                support=support)


# ----------------------------------------------------------------------------
# Pauli count via openfermion JW (Hermitian effective operator only)
#   h_so : (nso,nso) spin-orbital 1-body ; asym : (nso,nso,nso,nso) <pq||rs>
#   Reconstruct a symmetric physicist 2-body tensor V[p,q,r,s] from asym for an
#   InteractionOperator, then JW.  Returns (#pauli_terms, 1-norm_of_paulis).
# ----------------------------------------------------------------------------
def pauli_count(h_so, eri_phys_so, const=0.0):
    """eri_phys_so[p,q,r,s] = <pq|rs> physicist (spin-orbital, chemistry-free).
    Build OF InteractionOperator with one_body=h, two_body in OF convention."""
    import openfermion as of
    nso = h_so.shape[0]
    # OpenFermion InteractionOperator: H = const + sum h[p,q] a_p^ a_q
    #   + sum two_body[p,q,r,s] a_p^ a_q^ a_r a_s
    # with two_body[p,q,r,s] = 1/2 <pq|sr> (physicist) mapping to a_p^ a_q^ a_r a_s.
    two = np.zeros((nso, nso, nso, nso))
    for p in range(nso):
        for q in range(nso):
            for r in range(nso):
                for s in range(nso):
                    # <pq|sr> physicist -> coefficient of a_p^ a_q^ a_r a_s
                    two[p, q, r, s] = 0.5 * eri_phys_so[p, q, s, r]
    io = of.InteractionOperator(const, h_so, two)
    qop = of.jordan_wigner(of.get_fermion_operator(io))
    qop.compress(1e-9)
    terms = [c for t, c in qop.terms.items() if t != ()]  # drop identity
    return dict(n_pauli=len(terms), l1_pauli=float(sum(abs(c) for c in terms)))


def eri_phys_spinorbital(eri_spatial, nso):
    """<pq|rs> physicist in spin-orbitals from spatial physicist eri[i,j,k,l]=<ij|kl>."""
    out = np.zeros((nso, nso, nso, nso))
    for p in range(nso):
        for q in range(nso):
            for r in range(nso):
                for s in range(nso):
                    if spin(p) == spin(r) and spin(q) == spin(s):
                        out[p, q, r, s] = eri_spatial[spatial(p), spatial(q),
                                                      spatial(r), spatial(s)]
    return out


# ============================================================================
# Assemble the p-inclusive system at (max_n, k, gamma)
# ============================================================================
def assemble(max_n=2, k=1.5, gamma=1.0, Z=3, n_elec=3, Ng=1200, nx=200,
             Lmax=2, want_exact3=True, ref_occ=(0, 1, 2)):
    states = build_states(max_n)
    ns = len(states)
    r, wr = make_grid(k, Ng=Ng)
    W = r * r * wr
    R = radial_table(states, r, k)

    S = build_overlap(states, R, W)
    h1 = build_h1(states, S, k, Z)

    # radial multipole kernels on the SAME grid
    Mcoul = coulomb_ML(r, Lmax)
    Mw = scalar_ML(r, lambda rr: w_kernel(rr, gamma), Lmax, nx=nx)
    Muprime = scalar_ML(r, lambda rr: 0.5 * np.exp(-gamma * rr), Lmax, nx=nx)

    eri_coul = build_eri_scalar(states, R, W, Mcoul, Lmax)
    eri_w = build_eri_scalar(states, R, W, Mw, Lmax)
    V3 = build_V3(states, R, W, Muprime, Lmax) if want_exact3 else None

    # Loewdin (l,m-diagonal S -> block S^{-1/2}) then transform integrals
    X = lowdin(S)
    h1o = transform_1(h1, X)
    eri_coul_o = transform_2(eri_coul, X)
    eri_w_o = transform_2(eri_w, X)
    V3o = transform_3(V3, X) if V3 is not None else None

    nso = 2 * ns
    dets, didx = make_dets(nso, n_elec)
    hso = h_spin(h1o, nso)
    asym_coul = asym_from_phys(eri_coul_o, nso)
    asym_w = asym_from_phys(eri_w_o, nso)

    out = dict(max_n=max_n, ns=ns, nso=nso, k=k, gamma=gamma, ndet=len(dets),
               states=states, ref_occ=ref_occ)

    # plain FCI (Hermitian Coulomb)
    H_plain = build_H(dets, didx, hso, asym_coul, nso)
    out['E_plain'], _ = ground(H_plain, hermitian=True)

    # xTC: replace Coulomb by w (Hermitian TC 2-body) + contracted-L3 v2
    if V3o is not None:
        v2, v1, v0 = xtc_contract(V3o, nso, ref_occ)
        hso_x = hso + v1
        asym_x = asym_w + v2
        H_xTC = build_H(dets, didx, hso_x, asym_x, nso, v0=v0)
        out['E_xTC'], out['imag_xTC'] = ground(H_xTC)
        out['xtc_v0'] = v0

        # exact 3-body (for the contraction-fidelity sanity, cheap at 120 dets)
        if want_exact3:
            H3 = build_H3(dets, didx, V3o, nso)
            H_w = build_H(dets, didx, hso, asym_w, nso)
            out['E_TC2'], _ = ground(H_w)
            out['E_exactTC'], _ = ground(H_w + H3)

        # ---- keep arrays for the sparsity study ----
        out['_arr'] = dict(V3o=V3o, v2=v2, v1=v1, v0=v0, nso=nso,
                           asym_coul=asym_coul, asym_w=asym_w, asym_x=asym_x,
                           hso=hso, hso_x=hso_x,
                           eri_coul_o=eri_coul_o, eri_w_o=eri_w_o)
    return out


if __name__ == '__main__':
    o = assemble()
    print('plain', o['E_plain'], 'xTC', o.get('E_xTC'))
