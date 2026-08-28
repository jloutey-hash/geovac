"""Is the l>0 rank payoff FLAT in basis size?  (the load-bearing claim)"""
import importlib.util, io, json, os, sys
import numpy as np
from scipy.linalg import eigh
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ".")
spec = importlib.util.spec_from_file_location("sp", "debug/ee_split_sp_fci.py")
sp = importlib.util.module_from_spec(spec); spec.loader.exec_module(sp)
from geovac import transcorrelated_sturmian as TC
from geovac.xtc_angular_sparsity import gA, gB

def radmap(orbs, rp_index):
    n = len(orbs)
    idx = np.zeros((n, n), dtype=int)
    for a, (na, la, _) in enumerate(orbs):
        for c, (nc, lc, _) in enumerate(orbs):
            idx[a, c] = rp_index[(min((na, la), (nc, lc)), max((na, la), (nc, lc)))]
    return idx

def assemble_fast(orbs, idx, RLd, U, Lmax):
    n = len(orbs)
    g = np.zeros((n,)*4, dtype=complex)
    for L in range(Lmax + 1):
        pref = 4*np.pi/(2*L+1); RL = RLd[L]
        RLfull = RL[idx[:, :, None, None], idx[None, None, :, :]]   # [a,c,b,d]
        for M in range(-L, L+1):
            A = np.array([[gA(la, ma, L, M, lc, mc) for (_, lc, mc) in orbs]
                          for (_, la, ma) in orbs])
            B = np.array([[gB(lb, mb, L, M, ld, md) for (_, ld, md) in orbs]
                          for (_, lb, mb) in orbs])
            if np.abs(A).max() < 1e-14 or np.abs(B).max() < 1e-14: continue
            g += pref*np.einsum("ac,bd,acbd->abcd", A, B, RLfull, optimize=True)
    gr = np.einsum("ap,bq,cr,ds,pqrs->abcd", U.conj(), U.conj(), U, U, g, optimize=True)
    return np.real(gr), float(np.abs(gr.imag).max())

print(f"{'basis':>10}{'n_orb':>7}{'n_rp':>6}{'E_full':>15}{'r@1mHa':>9}{'r@0.01mHa':>11}")
out = {}
for ns, npp in [(3,1), (4,1), (3,2), (4,2), (5,2)]:
    orbs, rad, S, h1, rpi, RL, RLs, RLw = sp.build(ns, npp, 2.0, 2.0, Ng=700, Lmax=2)
    U = sp.real_harmonic_transform(orbs); idx = radmap(orbs, rpi)
    Sr, hr, _ = sp.one_body_real(S, h1, U)
    g, im = assemble_fast(orbs, idx, RL, U, 2)
    gs, _ = assemble_fast(orbs, idx, RLs, U, 2)
    assert im < 1e-12, f"imag {im:.1e}"
    E_full, nd = sp.fci(Sr, hr, g, 2)
    r1 = r001 = None
    for m in range(1, len(rpi)+1):
        RLt = {}
        for L in range(3):
            lam, V = np.linalg.eigh(RLw[L]); o = np.argsort(-np.abs(lam))[:m]
            RLt[L] = (V[:, o]*lam[o]) @ V[:, o].T
        gwt, _ = assemble_fast(orbs, idx, RLt, U, 2)
        d = abs(sp.fci(Sr, hr, gs-gwt, 2)[0] - E_full)*1000
        if r1 is None and d < 1.0: r1 = m
        if r001 is None and d < 0.01: r001 = m; break
    print(f"{f'{ns}s+{npp}p':>10}{len(orbs):>7}{len(rpi):>6}{E_full:>15.8f}{str(r1):>9}{str(r001):>11}")
    out[f"{ns}s+{npp}p"] = dict(n_orb=len(orbs), n_rp=len(rpi), E=E_full, r1=r1, r001=r001)
os.makedirs("debug/data", exist_ok=True)
io.open("debug/data/ee_split_sp_sweep.json","w",encoding="utf-8").write(json.dumps(out, indent=2, default=float))
print("\nwrote debug/data/ee_split_sp_sweep.json")
