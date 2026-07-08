"""
The sparsity-DESTROYING option (PI-directed 2026-07-07): the balanced builder
NEGLECTS the cross-center overlap (treats S=I; PK is an operator approximation of
strict Schmidt orthogonalization, balanced_coupled.py l.476). Here we RESTORE the
real overlap S and Loewdin-orthogonalize the balanced integrals (h1,eri) properly,
accept the densification, and measure:
  (1) does destroying the sparsity FIX the geometry? (omega_e +45%, R_eq +7% sparse)
  (2) how much does the Pauli count inflate? (Track DF Sprint 5 measured 14x but
      NOT the geometry-improvement -- that is the new read-out here.)

S structure (LiH, M=15): 3 sub-blocks of 5 (1s,2s,2p-1,2p0,2p+1):
  Li_core (Z=3 @ origin), LiH_bond_center (Z=1 @ origin), LiH_bond_partner (Z=1 @ H@R).
  within-block = I; core<->bond_center = same-centre different-Z radial overlap
  (delta_l delta_m); {core,bond_center}<->partner = two-centre overlap (m-conserving).
Loewdin: X=S^-1/2; h1'=X h1 X; eri'=einsum(X X X X, eri). FCI via coupled_fci_energy.
"""
from __future__ import annotations
import importlib.util, json, time
import numpy as np
from scipy.linalg import sqrtm
from scipy.special import genlaguerre, factorial
from geovac.balanced_coupled import (build_balanced_hamiltonian, _get_block_geometry,
                                      _get_nuclei_for_lih)
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec
from openfermion import jordan_wigner
from geovac.qubit_encoding import build_fermion_op_from_integrals

from scipy.special import lpmv

# --- robust two-centre overlap via FIXED Gauss-Legendre in prolate spheroidal
# coords (no adaptive mpmath hang; the topos3 quad is the documented WRONG TOOL). ---
_gl_xi_u, _gl_xi_w = np.polynomial.legendre.leggauss(80)
_gl_eta, _gl_eta_w = np.polynomial.legendre.leggauss(60)
def _ang(l, m, ct):
    m = abs(m)
    norm = np.sqrt((2*l+1)/2.0 * factorial(l-m)/factorial(l+m))
    return norm * lpmv(m, l, ct)
def overlap_two_center(Z1, n1, l1, Z2, n2, l2, m, R):
    """<phi_{n1 l1 m}(0,Z1) | phi_{n2 l2 m}(R zhat, Z2)>, m conserved."""
    R = float(R); half = R/2.0
    xi_max = 1.0 + 90.0/max((float(Z1)+float(Z2))*R, 1e-6)   # integrand ~ e^{-(Z1+Z2)(R/2)xi}
    xi = 1.0 + (xi_max-1.0)*(_gl_xi_u+1.0)/2.0
    wxi = _gl_xi_w*(xi_max-1.0)/2.0
    XI, ETA = np.meshgrid(xi, _gl_eta, indexing='ij')
    W = np.outer(wxi, _gl_eta_w)
    r1 = half*(XI+ETA); r2 = half*(XI-ETA)
    ct1 = (1+XI*ETA)/(XI+ETA); ct2 = (XI*ETA-1)/(XI-ETA)
    ct1 = np.clip(ct1, -1, 1); ct2 = np.clip(ct2, -1, 1)
    integ = (R_nl(Z1,n1,l1,r1)*_ang(l1,m,ct1)
             * R_nl(Z2,n2,l2,r2)*_ang(l2,m,ct2)
             * half**3 * (XI**2 - ETA**2))
    return float(np.sum(integ*W))

R_TRUE = 3.015
CM1 = 1.0/219474.6313702
MU = (7.016004*1.007825)/(7.016004+1.007825)*1822.888486
K_TRUE = MU*(1405.65*CM1)**2

_rgrid = np.linspace(1e-6, 80.0, 6000)
def R_nl(Z, n, l, r):
    Z = float(Z); rho = 2*Z*r/n
    norm = np.sqrt((2*Z/n)**3 * factorial(n-l-1)/(2*n*factorial(n+l)))
    return norm*np.exp(-rho/2)*rho**l*genlaguerre(n-l-1, 2*l+1)(rho)
def radial_overlap(Z1,n1,Z2,n2,l):
    f = R_nl(Z1,n1,l,_rgrid)*R_nl(Z2,n2,l,_rgrid)*_rgrid**2
    return float(np.trapezoid(f,_rgrid))

def build_S(sub, R):
    """Full non-orthogonal overlap: within-block = I; core<->bond_center (both at
    origin, Z=3 vs Z=1) = same-centre radial overlap (delta_l delta_m); {core,
    bond_center}<->partner (at R) = two-centre overlap (m-conserving). Robust GL quad."""
    orbs = []  # (Z, n, l, m, at_origin)
    for sb in sub:
        at_origin = (sb['side'] != 'partner')
        for (n,l,m) in sb['states']:
            orbs.append((float(sb['Z_orb']), n, l, m, at_origin))
    M = len(orbs); S = np.eye(M); n_tc = 0
    for i in range(M):
        Zi,ni,li,mi,oi = orbs[i]
        for j in range(i+1, M):
            Zj,nj,lj,mj,oj = orbs[j]
            if mi != mj:
                v = 0.0
            elif oi and oj:                                  # same centre (origin)
                v = radial_overlap(Zi,ni,Zj,nj,li) if li == lj else 0.0
            else:                                            # origin <-> partner: two-centre
                v = overlap_two_center(Zi,ni,li, Zj,nj,lj, mi, R); n_tc += 1
            S[i,j] = S[j,i] = v
    return S, orbs, n_tc

def lowdin_floor_X(S, eps=0.05):
    """Symmetric Loewdin with eigenvalue flooring (keeps dimension M so the FCI's
    block/qubit structure is preserved). The restored-overlap basis is near-linearly-
    dependent (redundant s-functions); the near-null eigenvalue is floored to eps
    rather than dropped. Returns X (symmetric, MxM) and n floored."""
    w, U = np.linalg.eigh(S)
    n_floored = int(np.sum(w < eps))
    wf = np.maximum(w, eps)
    return U @ np.diag(1.0/np.sqrt(wf)) @ U.T, n_floored

def transform_eri(eri, X):                    # X: (M, M) symmetric
    return np.einsum('ip,jq,kr,ls,ijkl->pqrs', X, X, X, X, eri, optimize=True)

def run(R):
    spec = lih_spec(R=R, max_n=2)
    sub = _get_block_geometry(spec)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=2000, L_max=4,
                                     screened_cross_center=False, verbose=False)
    h1 = np.array(ham['h1']); eri = np.array(ham['eri']); nuc = ham['nuclear_repulsion']
    S, orbs, n_tc = build_S(sub, R)
    smin = float(np.min(np.linalg.eigvalsh(S)))
    X, n_dropped = lowdin_floor_X(S)          # floored symmetric Loewdin (keeps M)
    h1p = X @ h1 @ X
    erip = transform_eri(eri, X)
    # sparse (S=I) baseline FCI (what the builder currently does)
    E_sparse = float(coupled_fci_energy(ham, n_electrons=n_e, verbose=False)['E_coupled'])
    # Loewdin FCI: swap in transformed integrals
    ham2 = dict(ham); ham2['h1'] = h1p; ham2['eri'] = erip
    E_lowdin = float(coupled_fci_energy(ham2, n_electrons=n_e, verbose=False)['E_coupled'])
    # Pauli counts
    n_pauli_sparse = len(ham['qubit_op'].terms)
    qop = jordan_wigner(build_fermion_op_from_integrals(h1p, erip, nuc))
    n_pauli_lowdin = len(qop.terms)
    return dict(R=R, S_min_eig=smin, n_dropped=n_dropped, E_sparse=E_sparse,
                E_lowdin=E_lowdin, n_pauli_sparse=n_pauli_sparse,
                n_pauli_lowdin=n_pauli_lowdin,
                pauli_ratio=n_pauli_lowdin/max(n_pauli_sparse,1))

if __name__ == "__main__":
    grid = [2.7, 2.9, 3.1, 3.3, 3.5]
    rows = []; t0 = time.perf_counter()
    for R in grid:
        ta = time.perf_counter(); r = run(R); rows.append(r)
        print(f"R={R}: S_min={r['S_min_eig']:+.4f} drop={r['n_dropped']}  "
              f"E_sparse={r['E_sparse']:.5f}  E_lowdin={r['E_lowdin']:.5f}  "
              f"Pauli {r['n_pauli_sparse']}->{r['n_pauli_lowdin']} "
              f"({r['pauli_ratio']:.1f}x)  [{time.perf_counter()-ta:.0f}s]", flush=True)
    json.dump(rows, open('debug/data/lowdin_tradeoff.json','w'), indent=2)

    def autopsy(key):
        R = np.array([r['R'] for r in rows]); E = np.array([r[key] for r in rows])
        c = np.polyfit(R-R_TRUE, E, 3); p = np.poly1d(c); dp=p.deriv(1); ddp=p.deriv(2)
        roots = dp.r[np.isreal(dp.r)].real
        mins = [x for x in roots if ddp(x)>0 and R.min()-R_TRUE<x<R.max()-R_TRUE]
        if not mins: return None
        xeq = min(mins, key=lambda x: p(x)); Req = xeq+R_TRUE
        curv = float(ddp(xeq)); omega = np.sqrt(max(curv,0)/MU)/CM1
        return dict(R_eq=Req, R_eq_err=abs(Req-R_TRUE)/R_TRUE*100, curv=curv,
                    curv_over_k=curv/K_TRUE, omega_e=omega, omega_err=(omega/1405.65-1)*100)
    print("\n=== GEOMETRY: does destroying sparsity fix it? ===")
    print(f"  experiment: R_eq=3.015, omega_e=1406 cm^-1, k={K_TRUE:.4f}")
    for key in ('E_sparse','E_lowdin'):
        a = autopsy(key)
        if a: print(f"  {key:9s}: R_eq={a['R_eq']:.3f} (+{a['R_eq_err']:.1f}%)  "
                    f"omega_e={a['omega_e']:.0f} ({a['omega_err']:+.0f}%)  curv {a['curv_over_k']:.2f}x k")
        else: print(f"  {key:9s}: no interior minimum")
    print(f"  Pauli inflation (mean): {np.mean([r['pauli_ratio'] for r in rows]):.1f}x")
    print(f"  total {time.perf_counter()-t0:.0f}s")
