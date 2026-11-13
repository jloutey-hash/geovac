"""Item #2 (WH1 almost-commutative extension): kernel checks on the true CH window tensored with M_k.

(A)  Left translations alone: kernel on O_n (x) M_k is 1 (x) M_k (dim k^2) -- the matrix factor is not metrized.
(A') Left translations + (right Spin(4) rotations (x) Ad pi_k): kernel is the scalars (dim 1) -- the combined
     SU(2)_L x SU(2)_R action is ergodic, so it gives a Lip-norm with no extra input.
(B)  Dirac route with the finite factor (M_k on C^k, D_F hermitian): ker [D (x) 1, .] cap ker [1 (x) D_F, .]
     contains 1 (x) C[D_F], so it is never the scalars for k >= 2 (the Connes metric on states is degenerate).
(C)  The compression P_n commutes with the right (+spin) generators: the truncation respects the action.
Usage: python debug/p38_ac_extension_check.py
"""
import sys, os
from fractions import Fraction
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests"))
import test_p38_true_ch_window as W  # reuse the frozen true-window construction

SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex), np.diag([1, -1]).astype(complex)]


def generators(n_max):
    labels, idx, D, ev, U = W._build(n_max)
    N = len(labels)
    Lg = [np.zeros((N, N), complex) for _ in range(3)]
    Rg = [np.zeros((N, N), complex) for _ in range(3)]
    js = sorted({lab[0] for lab in labels})
    for j in js:
        d = int(2 * j) + 1
        Jm = W._spin_mats(j)
        ms = [j - k for k in range(d)]
        ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
        I_d, I_2 = np.eye(d), np.eye(2)
        for a in range(3):
            Lg[a][np.ix_(ids, ids)] = np.kron(np.kron(Jm[a], I_d), I_2)
            Rg[a][np.ix_(ids, ids)] = np.kron(np.kron(I_d, Jm[a]), I_2) + np.kron(np.kron(I_d, I_d), SIG[a] / 2)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = U[:, keep]
    return Lg, Rg, V, D


def op_basis(n_max):
    _, vecs, Dc, dimH, _ = W._window(n_max)
    dim_O = int(np.linalg.matrix_rank(vecs, tol=1e-9))
    _, _, Vh = np.linalg.svd(vecs, full_matrices=False)
    return [b.reshape(dimH, dimH) for b in Vh[:dim_O]], Dc, dimH


def kernel_dim(basis, gens):
    """dim of {T in span(basis): [G, T] = 0 for all G in gens} (T linear in basis coefficients)."""
    rows = [np.concatenate([(G @ B - B @ G).reshape(-1) for G in gens]) for B in basis]
    M = np.array(rows).T
    return len(basis) - int(np.linalg.matrix_rank(M, tol=1e-8))


def spin_rep(k):
    return W._spin_mats(Fraction(k - 1, 2))


def main():
    rng = np.random.default_rng(1)
    for n_max in (1, 2):
        Lg, Rg, V, D = generators(n_max)
        # (C) P_n commutes with right(+spin) and left generators; D commutes with both
        P = V @ V.conj().T
        cP = max(np.abs(P @ G - G @ P).max() for G in Lg + Rg)
        cD = max(np.abs(D @ G - G @ D).max() for G in Lg + Rg)
        Lc = [V.conj().T @ G @ V for G in Lg]
        Rc = [V.conj().T @ G @ V for G in Rg]
        basisO, Dc, dimH = op_basis(n_max)
        print(f"n_max={n_max}: dim O={len(basisO)}  [P,gen]={cP:.1e}  [D,gen]={cD:.1e}")
        for k in (1, 2, 3):
            Ik = np.eye(k)
            E = [np.zeros((k, k), complex) for _ in range(k * k)]
            for i in range(k * k):
                E[i][i // k, i % k] = 1
            basis = [np.kron(B, e) for B in basisO for e in E]
            pik = spin_rep(k)
            left = [np.kron(G, Ik) for G in Lc]
            right = [np.kron(G, Ik) + np.kron(np.eye(dimH), pik[a]) for a, G in enumerate(Rc)]
            kA = kernel_dim(basis, left)
            kA2 = kernel_dim(basis, left + right)
            X = rng.standard_normal((k, k)) + 1j * rng.standard_normal((k, k))
            DF = X + X.conj().T
            kB = kernel_dim(basis, [np.kron(Dc, Ik), np.kron(np.eye(dimH), DF)])
            print(f"  k={k}: (A) left-only kernel={kA} (k^2={k*k})  (A') left+right(x)Ad kernel={kA2}  "
                  f"(B) Dirac+D_F kernel={kB} (k={k})")
            assert kA == k * k and kA2 == 1 and kB == k
        assert cP < 1e-9 and cD < 1e-9
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
