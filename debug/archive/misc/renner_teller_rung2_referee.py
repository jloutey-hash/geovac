"""Referee for the pi-state energy at R = 2 (Z_c = 1, two protons, linear):
VARIATIONAL Galerkin in a fixed-k Coulomb-Sturmian basis (a complete set as
n_max grows), optimised over k.  With functions normalised <1/rho> = 1/n:
H(k) = k^2 (I - O/2) - k T(kR),  H c = E O c,  O_ij = int S_i S_j rho^2.
At k = p0 with E = -k^2/2 this reproduces the secular equation T c = k c,
so the secular solution is a stationary point of the Galerkin family.
Restricted to the pi_x subspace {n p_x} (exact at s = 0)."""
import importlib.util
import os

import mpmath as mp
import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    'rs', os.path.join(ROOT, 'tests', 'test_renner_teller_rung2_sturmian.py'))
rs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rs)
mp.mp.dps = 40
R = 2

for nm in (3, 5, 7, 9):
    B, T = rs.t_series(nm, l_max=1)
    idx = [i for i, b in enumerate(B) if b[2] == 'x']
    O = mp.matrix(len(idx), len(idx))
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            O[a, b] = mp.mpf(sp.N(sp.integrate(rs.sturmian(B[i][0], 1) * rs.sturmian(B[j][0], 1)
                                               * rs.rho**2, (rs.rho, 0, sp.oo)), 40))
    T0 = sp.lambdify(rs.Ss, T[0], 'mpmath')
    L = mp.cholesky(O)
    Li = mp.inverse(L)
    def lowest(k):
        Tk = mp.matrix(T0(k * R))
        Ts = mp.matrix(len(idx), len(idx))
        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                Ts[a, b] = Tk[i, j]
        H = k**2 * (mp.eye(len(idx)) - O / 2) - k * Ts
        E = mp.eigsy(Li * H * Li.T)[0]
        return min(E[t] for t in range(len(idx)))
    kopt = mp.findroot(lambda k: mp.diff(lowest, k), mp.mpf('1.0'))
    print('n_max=%d  p-functions=%d  Galerkin min_k E_pi = %.10f at k = %.6f'
          % (nm, len(idx), float(lowest(kopt)), float(kopt)))
