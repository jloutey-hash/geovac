"""Validation gates for the xTC PoC (run before the study)."""
import numpy as np
import xtc_poc_li as X
from geovac.sturmian_solver import _radial_overlap, _slater_rk, hydrogenic_radial

def gate_integrals():
    print("=== GATE 1: grid integrals vs analytic (s-only Sturmian) ===")
    ns, k, Z = 4, 1.3, 3
    r, wr = X.make_grid(k, Ng=2000)
    S, h1, Rtab, W = X.build_one_body(ns, r, wr, k, Z)
    # analytic overlap
    maxS = 0.0
    for i in range(ns):
        for j in range(ns):
            Zi, Zj = (i+1)*k, (j+1)*k
            ana = _radial_overlap(i+1, 0, Zi, j+1, 0, Zj)
            maxS = max(maxS, abs(S[i, j]-ana))
    # analytic h1 diagonal:  Sturmian  h1_ii = k^2/2 - Z*k/n_i ; offdiag = -k^2/2 S_ij
    maxh = 0.0
    for i in range(ns):
        for j in range(ns):
            ni, nj = i+1, j+1
            Zi, Zj = ni*k, nj*k
            Sij = _radial_overlap(ni, 0, Zi, nj, 0, Zj)
            ana = (k**2 - Z*k/ni) if i == j else 0.0
            ana += -k**2/2.0 * Sij
            maxh = max(maxh, abs(h1[i, j]-ana))
    # coulomb eri vs slater R^0
    eri_c, _, _ = X.two_body(ns, Rtab, W, X.build_kernels(r, 1.0, nx=64))
    maxe = 0.0
    for (i, j, kk, l) in [(0,0,0,0),(0,1,0,1),(1,1,1,1),(0,1,1,0),(0,2,1,3),(2,3,2,3)]:
        Zi,Zj,Zk,Zl=(i+1)*k,(j+1)*k,(kk+1)*k,(l+1)*k
        ana = _slater_rk(i+1,0,Zi, j+1,0,Zj, kk+1,0,Zk, l+1,0,Zl, 0)
        maxe = max(maxe, abs(eri_c[i,j,kk,l]-ana))
    print(f"  max|dS|={maxS:.2e}  max|dh1|={maxh:.2e}  max|d eri_coul|={maxe:.2e}")
    return maxS < 1e-4 and maxh < 1e-4 and maxe < 2e-4

def gate_he_plain():
    print("=== GATE 2: He s-only plain FCI (should approach s-limit -2.87903) ===")
    best = 1e9
    for k in [1.6,1.7,1.8]:
        o = X.assemble(4, k, 1.0, Z=2, n_elec=2, Ng=800, nx=64, want_exact3=False)
        best = min(best, o['E_plain'])
    print(f"  best He s-only plain FCI (ns=4) = {best:.6f}  (He s-limit -2.879029)")
    return -2.90 < best < -2.86

def gate_xtc_exact(ns=3, k=1.5, gamma=1.2):
    print(f"=== GATE 3: xTC effective op vs exact 3-body on ref-row (ns={ns},k={k},g={gamma}) ===")
    o = X.assemble(ns, k, gamma, Z=3, n_elec=3, Ng=800, nx=128, want_exact3=True)
    A = o['_arrays']; H3 = A['H3']; dets = A['dets']; didx = A['didx']
    ref = A['ref_occ']; nso = A['nso']
    v2, h1b, v0b = A['v2'], A['v1'], A['v0']    # now BARE coeffs from module
    Iref = didx[tuple(sorted(ref))]
    H_eff = X.build_H(dets, didx, h1b, v2, nso, v0=v0b)  # effective-only (no TC2)
    refset = set(ref)
    max_le2 = 0.0; max_triple = 0.0
    for J, dJ in enumerate(dets):
        nexc = 3 - len(refset & set(dJ))
        d = abs(H_eff[Iref, J] - H3[Iref, J])
        if nexc <= 2:
            max_le2 = max(max_le2, d)
        else:
            max_triple = max(max_triple, abs(H3[Iref, J]))
    print(f"  ref diag: H_eff={H_eff[Iref,Iref]:.8f}  H3={H3[Iref,Iref]:.8f}")
    print(f"  ref-row max|H_eff-H3| over <=2-exc = {max_le2:.2e}  (must be ~0: exact on ref+S+D)")
    print(f"  dropped pure-3body residual (triples) max|H3| = {max_triple:.2e}  (nonzero = genuine 3-body)")
    print(f"  E_TC2={o['E_TC2']:.6f}  E_exactTC={o['E_exactTC']:.6f}  E_xTC={o['E_xTC']:.6f}")
    print(f"  L3 exact shift = {(o['E_exactTC']-o['E_TC2'])*1e3:+.3f} mHa ; "
          f"xTC shift = {(o['E_xTC']-o['E_TC2'])*1e3:+.3f} mHa")
    return max_le2 < 1e-9

def gate_geminal_zero():
    print("=== GATE 4: geminal->0 (large gamma) must reproduce plain FCI ===")
    ns, k = 3, 1.5
    Ep = None
    for gamma in [4.0, 10.0, 25.0]:
        o = X.assemble(ns, k, gamma, Z=3, n_elec=3, Ng=1000, nx=160, want_exact3=True)
        Ep = o['E_plain']
        print(f"  g={gamma:5.1f}: plain={o['E_plain']:.6f} TC2={o['E_TC2']:.6f} "
              f"exactTC={o['E_exactTC']:.6f} xTC={o['E_xTC']:.6f}  "
              f"|xTC-plain|={abs(o['E_xTC']-o['E_plain']):.2e}")
    return abs(o['E_xTC']-Ep) < 5e-3

if __name__ == '__main__':
    import time
    t0 = time.time()
    g1 = gate_integrals()
    g2 = gate_he_plain()
    g3 = gate_xtc_exact()
    g4 = gate_geminal_zero()
    print("\nGATES:", dict(integrals=g1, he_plain=g2, xtc_exact=g3, geminal_zero=g4))
    print(f"wall {time.time()-t0:.1f}s")
