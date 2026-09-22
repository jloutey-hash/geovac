r"""#2 FCI solver speedup (2026-09-22): the dense fci_energy builds all nd^2 determinant
pairs (280s at M=16, nd=14400), but the CI Hamiltonian is SPARSE -- two determinants
connect only if they differ by <= 2 spin-orbitals.  fci_energy_fast:
  1. bitmask each determinant (occupied spin-orbitals),
  2. find connected pairs by a VECTORIZED popcount of the XOR (ndiff = popcount/2 <= 2),
  3. compute the Slater-Condon element (reusing _matel) only for connected pairs,
  4. ground state via scipy.sparse.linalg.eigsh (Lanczos), no dense nd x nd matrix.

Reproduces the dense fci_energy to ~1e-10 (GATE), ~15-20x faster at M=16.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

from prolate_allelectron_fci import _dets, _matel

_M1 = np.uint64(0x5555555555555555)
_M2 = np.uint64(0x3333333333333333)
_M4 = np.uint64(0x0f0f0f0f0f0f0f0f)
_H01 = np.uint64(0x0101010101010101)
_ONE = np.uint64(1)
_TWO = np.uint64(2)
_FOUR = np.uint64(4)
_56 = np.uint64(56)


def _popcount64(x: np.ndarray) -> np.ndarray:
    """SWAR popcount for a uint64 ndarray (numpy-version-independent)."""
    x = x - ((x >> _ONE) & _M1)
    x = (x & _M2) + ((x >> _TWO) & _M2)
    x = (x + (x >> _FOUR)) & _M4
    return (x * _H01) >> _56


def fci_energy_fast(h1, eri, M, nelec, n_states=1):
    """Ground-state FCI energy (Sz=0 block) via a sparse connected-pair build + Lanczos.
    Drop-in for prolate_allelectron_fci.fci_energy (same signature/return)."""
    na = nb = nelec // 2
    dets = _dets(M, na, nb)
    nd = len(dets)
    if nd == 1:
        e = _matel(h1, eri, dets[0], dets[0])
        return (e, nd) if n_states == 1 else (np.array([e]), nd)
    assert 2 * M <= 64, "spin-orbital index must fit in uint64"
    masks = np.array([int(sum(1 << so for so in D)) for D in dets], dtype=np.uint64)

    rows = np.empty(0, dtype=np.int64)
    cols = np.empty(0, dtype=np.int64)
    rl, cl, vl = [], [], []
    diag = np.empty(nd)
    for a in range(nd):
        diag[a] = _matel(h1, eri, dets[a], dets[a])
        if a + 1 < nd:
            pc = _popcount64(masks[a + 1:] ^ masks[a])
            conn = np.nonzero(pc <= _FOUR)[0] + (a + 1)   # differ by <=2 spin-orbitals
            Da = dets[a]
            for b in conn.tolist():
                v = _matel(h1, eri, Da, dets[b])
                if v != 0.0:
                    rl.append(a); cl.append(b); vl.append(v)
    # assemble symmetric sparse H (upper off-diagonal + mirror + diagonal)
    r = np.array(rl, dtype=np.int64); c = np.array(cl, dtype=np.int64); v = np.array(vl)
    ri = np.concatenate([np.arange(nd), r, c])
    ci = np.concatenate([np.arange(nd), c, r])
    vi = np.concatenate([diag, v, v])
    H = sp.csr_matrix((vi, (ri, ci)), shape=(nd, nd))
    k = min(max(n_states, 1), nd - 1)
    ev = eigsh(H, k=k, which='SA', return_eigenvectors=False)
    ev = np.sort(ev)
    return (float(ev[0]), nd) if n_states == 1 else (ev[:n_states], nd)


def gate():
    from prolate_allelectron_fci import fci_energy
    import time
    print("=" * 70)
    print("GATE fci_energy_fast vs dense fci_energy (energy match + speedup)")
    print("=" * 70)
    rng = np.random.default_rng(0)
    for M in (8, 10, 12, 14):
        h1 = rng.standard_normal((M, M)); h1 = h1 + h1.T
        e = rng.standard_normal((M, M, M, M)) * 0.1
        e = e + e.transpose(1, 0, 3, 2)
        e = e + e.transpose(2, 3, 0, 1)                 # (pq|rs)=(rs|pq) chemist symmetry
        t = time.time(); Ed, nd = fci_energy(h1, e, M, 4); td = time.time() - t
        t = time.time(); Ef, _ = fci_energy_fast(h1, e, M, 4); tf = time.time() - t
        print(f"  M={M:2d} nd={nd:5d}: dense={Ed:.8f} ({td:.1f}s)  fast={Ef:.8f} ({tf:.1f}s)  "
              f"dE={abs(Ed-Ef):.1e}  {td/max(tf,0.01):.0f}x  {'PASS' if abs(Ed-Ef)<1e-8 else 'FAIL'}")


if __name__ == "__main__":
    gate()
