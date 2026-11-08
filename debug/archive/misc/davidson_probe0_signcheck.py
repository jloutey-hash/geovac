"""
Probe 0 -- independent-convention check of geovac.coupled_composition.coupled_fci_energy.

Builds a brute-force second-quantized FCI in the FULL Fock space (bitmask
determinants, explicit a^dag / a with standard sign counting), restricts to the
(N_up, N_down) sector, and compares the ground energy to the library builder on
a small random-but-8-fold-symmetric integral set.

This is the "independent computation" leg of CLAUDE.md 13.4a: it does not reuse
any of the library's Slater-rule code.
"""
from __future__ import annotations
import itertools
import numpy as np
from geovac.coupled_composition import coupled_fci_energy


def random_integrals(M: int, seed: int = 0):
    rng = np.random.default_rng(seed)
    h1 = rng.standard_normal((M, M))
    h1 = 0.5 * (h1 + h1.T)
    # 8-fold-symmetric ERI via a Cholesky-like factorization: (pq|rs) = sum_P B[p,q,P] B[r,s,P]
    nP = M + 2
    B = rng.standard_normal((M, M, nP))
    B = 0.5 * (B + B.transpose(1, 0, 2))
    eri = np.einsum('pqP,rsP->pqrs', B, B)
    return h1, eri


# ---------------------------------------------------------------- brute force
def _apply_op(det: int, ops):
    """Apply a sequence of (kind, orbital) from RIGHT to LEFT. kind: 'a' or 'c'.
    Sign convention: a_P |...> = (-1)^{# occupied spin-orbitals with index < P} ..."""
    sign = 1
    for kind, P in ops:
        bit = 1 << P
        if kind == 'a':
            if not (det & bit):
                return 0, 0
            sign *= (-1) ** bin(det & (bit - 1)).count('1')
            det ^= bit
        else:
            if det & bit:
                return 0, 0
            sign *= (-1) ** bin(det & (bit - 1)).count('1')
            det |= bit
    return det, sign


def brute_force_fci(h1, eri, M, n_up, n_down, e_core=0.0):
    nso = 2 * M           # spin-orbital P: alpha p -> p, beta p -> M + p
    dets = []
    for a in itertools.combinations(range(M), n_up):
        for b in itertools.combinations(range(M), n_down):
            d = 0
            for p in a:
                d |= 1 << p
            for p in b:
                d |= 1 << (M + p)
            dets.append(d)
    idx = {d: i for i, d in enumerate(dets)}
    n = len(dets)
    H = np.zeros((n, n))

    def spin(P):
        return 0 if P < M else 1

    def orb(P):
        return P if P < M else P - M

    for I, det in enumerate(dets):
        # one-body:  sum_PQ h_PQ a^dag_P a_Q
        for Q in range(nso):
            if not (det >> Q) & 1:
                continue
            for P in range(nso):
                if spin(P) != spin(Q):
                    continue
                newdet, sg = _apply_op(det, [('a', Q), ('c', P)])
                if sg == 0:
                    continue
                J = idx.get(newdet)
                if J is None:
                    continue
                H[J, I] += sg * h1[orb(P), orb(Q)]
        # two-body:  1/2 sum_PQRS <PQ|v|RS> a^dag_P a^dag_Q a_S a_R
        #            <PQ|v|RS> = (pr|qs) delta_{sP,sR} delta_{sQ,sS}
        for R in range(nso):
            if not (det >> R) & 1:
                continue
            for S in range(nso):
                if S == R or not (det >> S) & 1:
                    continue
                for P in range(nso):
                    if spin(P) != spin(R):
                        continue
                    for Q in range(nso):
                        if spin(Q) != spin(S):
                            continue
                        newdet, sg = _apply_op(det, [('a', R), ('a', S), ('c', Q), ('c', P)])
                        if sg == 0:
                            continue
                        J = idx.get(newdet)
                        if J is None:
                            continue
                        H[J, I] += 0.5 * sg * eri[orb(P), orb(R), orb(Q), orb(S)]
    H += e_core * np.eye(n)
    return H, dets


if __name__ == '__main__':
    for M, n_up, n_down in [(4, 2, 2), (5, 2, 2), (4, 1, 1), (5, 2, 1)]:
        h1, eri = random_integrals(M, seed=M * 10 + n_up)
        e_core = 0.37
        H, dets = brute_force_fci(h1, eri, M, n_up, n_down, e_core)
        asym = np.abs(H - H.T).max()
        w = np.linalg.eigvalsh(H)
        res = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': e_core}
        out = None
        if n_up == n_down:
            out = coupled_fci_energy(res, n_electrons=n_up + n_down, verbose=False)
        print(f"M={M} (n_up,n_down)=({n_up},{n_down})  n_det={len(dets)}  "
              f"|H-H^T|max={asym:.2e}  E0_brute={w[0]:+.12f}"
              + (f"  E0_lib={out['E_coupled']:+.12f}  DIFF={w[0]-out['E_coupled']:+.3e}"
                 if out else "  (lib requires n_up==n_down)"))
