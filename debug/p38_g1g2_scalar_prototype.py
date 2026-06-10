"""Sprint P38-G1G2 Phase B: scalar end-to-end prototype of the GE-vS
dual (lifted-state) map on SU(2), with the central Fejer measure.

Construction (all exact in the Peter-Weyl coefficient picture):
  - Window W = {j : j <= j_max} in L^2(SU(2)); orthonormal basis
    e^j_{ab} = sqrt(2j+1) D^j_{ab}.
  - Fejer amplitude h = sum_j sqrt(d_j) chi_j;  in the orthonormal
    basis h has coefficient 1 on every diagonal element e^j_{aa} and 0
    elsewhere; ||h||^2 = Z = sum d_j.  xi = h/sqrt(Z).  The induced
    measure |xi|^2 dHaar is EXACTLY the central Fejer kernel measure.
  - Lifted state: mu(T) = <xi, T xi>.
  - Dual map: upsilon(T)(g) = <U_g xi, T U_g xi>  (U_g = left
    translation on the window).

Checks:
  T1 (exact structural identity): upsilon(P M_f P)(g) = (K*f)(g)
     = sigma_J f(g) for single-band f = D^J_{AB}, with sigma_J the
     corrected multiplier symbol c_J/(d_J Z),
     c_J = sum_{j1,j2 in W, triangle(j1,j2,J)} sqrt(d_{j1} d_{j2}).
     This simultaneously re-verifies the corrected Lemma-L2 symbols
     from an independent (vector-state) direction.
  T2 (almost-inverse rate): ||K*f - f||_grid <= gamma * Lip_grid(f),
     gamma = int K(x) d(e,x) dx computed on the class measure with the
     same metric convention as Lip_grid.
  T3 (tautological contractivity): Lip_grid(upsilon(T)) <=
     L_action(T) := sup_grid ||U_g T U_g^+ - T|| / d(g,e), for random
     T in the system span.

Run from repo root: python debug/p38_g1g2_scalar_prototype.py
"""

from __future__ import annotations

import json
from fractions import Fraction

import numpy as np
from scipy.linalg import expm
from sympy import S
from sympy.physics.quantum.cg import CG


# ----------------------------------------------------------------------
# Representation machinery
# ----------------------------------------------------------------------

def jm_range(j):
    """m values for spin j, descending."""
    k = int(round(2 * j))
    return [Fraction(k, 2) - i for i in range(k + 1)]


def su2_generators(j):
    """Jz, Jy (numpy, dimension 2j+1) in the |j m> basis."""
    ms = [float(m) for m in jm_range(j)]
    d = len(ms)
    Jz = np.diag(ms).astype(complex)
    Jp = np.zeros((d, d), dtype=complex)
    jf = float(j)
    for i in range(1, d):
        m = ms[i]  # raise m -> m+1 lands at index i-1
        Jp[i - 1, i] = np.sqrt(jf * (jf + 1) - m * (m + 1))
    Jm = Jp.conj().T
    Jy = (Jp - Jm) / (2j)
    return Jz, Jy


def wigner_D(j, alpha, beta, gamma):
    """D^j(g) for zyz Euler angles, unitary (2j+1)x(2j+1)."""
    Jz, Jy = su2_generators(j)
    return (expm(-1j * alpha * Jz) @ expm(-1j * beta * Jy)
            @ expm(-1j * gamma * Jz))


def rotation_angle(alpha, beta, gamma):
    """Geodesic distance d(e, g) on SU(2): the rotation angle of the
    j=1/2 matrix, in [0, pi] doubled to SU(2) range [0, pi]... we use
    theta with cos(theta) = Re tr(D^{1/2}(g))/2, theta in [0, pi];
    consistent convention used for BOTH gamma and Lip."""
    tr = np.trace(wigner_D(Fraction(1, 2), alpha, beta, gamma))
    c = np.clip(np.real(tr) / 2.0, -1.0, 1.0)
    return float(np.arccos(c))  # in [0, pi]


_CG_CACHE = {}


def cg(j1, m1, j2, m2, j3, m3):
    key = (j1, m1, j2, m2, j3, m3)
    if key not in _CG_CACHE:
        val = CG(S(j1), S(m1), S(j2), S(m2), S(j3), S(m3)).doit()
        _CG_CACHE[key] = float(val)
    return _CG_CACHE[key]


# ----------------------------------------------------------------------
# Window, basis indexing, multiplier matrices
# ----------------------------------------------------------------------

class Window:
    def __init__(self, j_max):
        self.js = []
        j = Fraction(0)
        while j <= j_max:
            self.js.append(j)
            j += Fraction(1, 2)
        self.index = {}
        idx = 0
        for j in self.js:
            for a in jm_range(j):
                for b in jm_range(j):
                    self.index[(j, a, b)] = idx
                    idx += 1
        self.dim = idx
        self.Z = sum(int(2 * j) + 1 for j in self.js)

    def xi(self):
        v = np.zeros(self.dim, dtype=complex)
        for j in self.js:
            for a in jm_range(j):
                v[self.index[(j, a, a)]] = 1.0
        return v / np.linalg.norm(v)

    def U(self, alpha, beta, gamma):
        """Left translation lambda_g on the window, orthonormal basis.
        lambda_g D^j_{ab} = sum_c conj(D^j(g))_{?}... implemented as:
        coefficient blocks C^j (a-row, b-col) transform C -> D^j(g)^+ C,
        which on the orthonormal flattened basis is block kron."""
        blocks = []
        for j in self.js:
            Dg = wigner_D(j, alpha, beta, gamma)
            d = Dg.shape[0]
            blocks.append(np.kron(Dg.conj(), np.eye(d)))
        from scipy.linalg import block_diag
        return block_diag(*blocks)

    def multiplier(self, J, A, B):
        """Matrix of M_f, f = D^J_{AB}, in the orthonormal window basis.
        <e^{j1}_{a1b1}, f e^{j2}_{a2b2}> =
            sqrt(d2/d1) <J A j2 a2|j1 a1><J B j2 b2|j1 b1>."""
        M = np.zeros((self.dim, self.dim), dtype=complex)
        for j2 in self.js:
            d2 = int(2 * j2) + 1
            for j1 in self.js:
                d1 = int(2 * j1) + 1
                if abs(J - j2) > j1 or j1 > J + j2:
                    continue
                pref = np.sqrt(d2 / d1)
                for a2 in jm_range(j2):
                    a1 = A + a2
                    if abs(a1) > j1:
                        continue
                    c_a = cg(J, A, j2, a2, j1, a1)
                    if c_a == 0.0:
                        continue
                    for b2 in jm_range(j2):
                        b1 = B + b2
                        if abs(b1) > j1:
                            continue
                        c_b = cg(J, B, j2, b2, j1, b1)
                        if c_b == 0.0:
                            continue
                        M[self.index[(j1, a1, b1)],
                          self.index[(j2, a2, b2)]] = pref * c_a * c_b
        return M


def f_eval(J, A, B, alpha, beta, gamma):
    """f(g) = D^J_{AB}(g)."""
    D = wigner_D(J, alpha, beta, gamma)
    ms = jm_range(J)
    return D[ms.index(A), ms.index(B)]


def sigma_symbol(win, J):
    """Corrected L2 multiplier symbol sigma_J = c_J/(d_J Z)."""
    c = 0.0
    for j1 in win.js:
        for j2 in win.js:
            if abs(j1 - j2) <= J <= j1 + j2 and (j1 + j2 - J) % 1 == 0:
                c += np.sqrt((int(2 * j1) + 1) * (int(2 * j2) + 1))
    return c / ((int(2 * J) + 1) * win.Z)


# ----------------------------------------------------------------------
# Main checks
# ----------------------------------------------------------------------

def main():
    rng = np.random.default_rng(7)
    win = Window(Fraction(1))         # W = {0, 1/2, 1};  n = 3 window
    xi = win.xi()
    print(f"window dim = {win.dim}, Z = {win.Z}")

    grid = [tuple(rng.uniform(0, 2 * np.pi, 3) * np.array([1, 0.5, 1]))
            for _ in range(24)]
    Us = {g: win.U(*g) for g in grid}

    # ---- T1: upsilon(P M_f P) == sigma_J * f, both g-conventions ----
    print("\nT1: upsilon(PM_fP)(g) vs sigma_J * f(g)   [exact identity]")
    bands = []
    for J in [Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2)]:
        ms = jm_range(J)
        A, B = ms[0], ms[-1]          # one representative element
        M = win.multiplier(J, A, B)
        sig = sigma_symbol(win, J)
        r_fwd, r_inv = [], []
        for g in grid:
            U = Us[g]
            xv = U @ xi
            ups = np.vdot(xv, M @ xv)
            fg = f_eval(J, A, B, *g)
            fginv = f_eval(J, A, B, -g[2], -g[1], -g[0])  # g^{-1} euler
            r_fwd.append(abs(ups - sig * fg))
            r_inv.append(abs(ups - sig * fginv))
        res = min(max(r_fwd), max(r_inv))
        conv = "g" if max(r_fwd) <= max(r_inv) else "g^{-1}"
        bands.append((float(J), sig, res, conv))
        print(f"  J={float(J):4.1f}  sigma={sig:.10f}  "
              f"max-resid={res:.2e}  (convention: {conv})")
    t1_pass = all(b[2] < 1e-10 for b in bands)
    print(f"  sigma values vs corrected-L2 closed forms: "
          f"sigma(1)={sigma_symbol(win, Fraction(1)):.6f} "
          f"(expect (5+2*sqrt3)/18 = {(5+2*np.sqrt(3))/18:.6f}), "
          f"sigma(2)={sigma_symbol(win, Fraction(2)):.6f} (expect 0.1)")
    print(f"  T1 verdict: {'PASS' if t1_pass else 'FAIL'}")

    # ---- gamma: first moment of the Fejer measure (class integral) ----
    # K(theta) = |h|^2/Z with h = sum sqrt(d_j) chi_j,
    # chi_j(theta) = sin((2j+1)theta/2)/sin(theta/2);
    # class measure (2/pi) sin^2(theta/2) dtheta, theta in [0, 2pi]?
    # SU(2) angle theta in [0, 2pi); d(e,g) = theta/2 in [0, pi] to
    # match the j=1/2 arccos convention above.
    th = np.linspace(1e-9, 2 * np.pi, 400001)
    h = np.zeros_like(th)
    for j in win.js:
        n = int(2 * j) + 1
        h += np.sqrt(n) * np.sin(n * th / 2) / np.sin(th / 2)
    K = h ** 2 / win.Z
    w = (1 / np.pi) * np.sin(th / 2) ** 2
    norm = np.trapezoid(K * w, th)
    gamma = np.trapezoid(K * w * (th / 2), th) / norm
    print(f"\nFejer-measure normalization check: {norm:.8f} (expect 1)")
    print(f"gamma (first moment) = {gamma:.6f}")

    # ---- T2: ||K*f - f||_grid <= gamma * Lip_grid(f) ----
    print("\nT2: almost-inverse rate")
    t2_ok = True
    for J in [Fraction(1, 2), Fraction(1), Fraction(2)]:
        ms = jm_range(J)
        A, B = ms[0], ms[-1]
        sig = sigma_symbol(win, J)
        vals = np.array([f_eval(J, A, B, *g) for g in grid])
        err = np.max(np.abs((sig - 1) * vals))
        # Lipschitz estimate on grid pairs
        lip = 0.0
        for i in range(len(grid)):
            for k in range(i + 1, len(grid)):
                gi, gk = grid[i], grid[k]
                # distance d(gi, gk) = angle of D(gi) D(gk)^+
                D1 = wigner_D(Fraction(1, 2), *gi)
                D2 = wigner_D(Fraction(1, 2), *gk)
                c = np.clip(np.real(np.trace(D1 @ D2.conj().T)) / 2, -1, 1)
                dd = np.arccos(c)
                if dd > 1e-6:
                    lip = max(lip, abs(vals[i] - vals[k]) / dd)
        ok = err <= gamma * lip * 1.0000001
        t2_ok &= ok
        print(f"  J={float(J):4.1f}: ||K*f-f||={err:.4f}  "
              f"gamma*Lip>={gamma * lip:.4f}  {'OK' if ok else 'VIOLATED'}")
    print(f"  T2 verdict: {'PASS' if t2_ok else 'FAIL'}")

    # ---- T3: contractivity of upsilon wrt action seminorm ----
    print("\nT3: Lip(upsilon(T)) <= L_action(T), random T in system span")
    t3_ok = True
    for trial in range(3):
        T = np.zeros((win.dim, win.dim), dtype=complex)
        for J in [Fraction(1, 2), Fraction(1), Fraction(2)]:
            ms = jm_range(J)
            A, B = ms[0], ms[-1]
            T += rng.normal() * win.multiplier(J, A, B)
        T = T + T.conj().T
        ups = np.array([np.vdot(Us[g] @ xi, T @ (Us[g] @ xi)) for g in grid])
        lip_u, L_act = 0.0, 0.0
        for i in range(len(grid)):
            for k in range(i + 1, len(grid)):
                D1 = wigner_D(Fraction(1, 2), *grid[i])
                D2 = wigner_D(Fraction(1, 2), *grid[k])
                c = np.clip(np.real(np.trace(D1 @ D2.conj().T)) / 2, -1, 1)
                dd = np.arccos(c)
                if dd > 1e-6:
                    lip_u = max(lip_u, abs(ups[i] - ups[k]) / dd)
                    # action seminorm evaluated at the pair quotient
                    Q = Us[grid[i]] @ Us[grid[k]].conj().T
                    L_act = max(L_act,
                                np.linalg.norm(Q @ T @ Q.conj().T - T, 2) / dd)
        ok = lip_u <= L_act * 1.0000001
        t3_ok &= ok
        print(f"  trial {trial}: Lip(upsilon(T))={lip_u:.4f}  "
              f"L_action(T)={L_act:.4f}  {'OK' if ok else 'VIOLATED'}")
    print(f"  T3 verdict: {'PASS' if t3_ok else 'FAIL'}")

    out = {"window_jmax": 1.0, "dim": win.dim, "gamma": float(gamma),
           "T1_bands": bands, "T1": bool(t1_pass), "T2": bool(t2_ok),
           "T3": bool(t3_ok)}
    with open("debug/data/p38_g1g2_scalar_prototype.json", "w") as f:
        json.dump(out, f, indent=1, default=float)
    print("\nOVERALL:", "PASS" if (t1_pass and t2_ok and t3_ok) else "FAIL")


if __name__ == "__main__":
    main()
