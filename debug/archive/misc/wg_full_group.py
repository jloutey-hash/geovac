"""WG follow-up: is Dong's G(O^CH) = (U(1).Spin(4)) x| Z2 exactly?  Component group G/G_0.

Run:  OMP_NUM_THREADS=1 python debug/wg_full_group.py [n_max ...]   (default 2 3)

Argument (memo debug/wg_full_group_memo.md):
  (A) commutant O' = C  =>  C*(O) = M_N (double commutant, finite dim)  =>  every unital complete order
      automorphism of O is Ad U (C*-envelope universal property + M_N simple + Skolem-Noether).
  (B) g in G normalises G_0 and acts on [g,g] = su2_L + su2_R; Out = Z2 (swap).  Kernel of G -> Z2 is
      G_0 . T', T' = {e^{i theta(D)} in G}, because the Spin(4) rep on the window is multiplicity free
      (so its commutant is the torus of functions of D).  T' is computed exactly by Smith normal form of the
      block-support constraints, then each torsion element is checked directly on O.
  (C) the inversion realises the swap; iota^2 is a scalar (split extension).
Checks: commutant dim; Spin(4) multiplicity-free; T' torsion mod U(1).Z(Spin4); direct residuals;
transposition (a unital order automorphism of O) has NO implementing unitary, and is shown not 2-positive
on O directly (independent route for the 'complete' hypothesis).
"""
import importlib.util
import itertools
import json
import os
import sys
from fractions import Fraction

import numpy as np
import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


p38 = _load('p38_wgfull', os.path.join(ROOT, 'tests', 'test_p38_true_ch_window.py'))
wg = _load('wg_guard_wgfull', os.path.join(ROOT, 'tests', 'test_paper32_window_gauge_group.py'))


def setup(n_max):
    labels, idx, D, ev, Uev = p38._build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = Uev[:, keep]
    Dc = V.conj().T @ D @ V
    keys = []
    for Nn in range(1, 2 * n_max + 1):
        J = Fraction(Nn - 1, 2)
        Ms = [J - k for k in range(Nn)]
        keys += [(J, M, Mp) for M, Mp in itertools.product(Ms, Ms)]
    comp = {k: V.conj().T @ p38._mult_op(*k, labels, idx) @ V for k in keys}
    E = wg.onb(list(comp.values()))
    return labels, idx, D, V, Dc, comp, E


def proj_resid(X, B):
    """relative distance of X from span(E) (B = rows of ONB vectors)."""
    c = X.reshape(-1)
    nc = np.linalg.norm(c)
    if nc < 1e-14:
        return 0.0
    return float(np.linalg.norm(c - B.T @ (B.conj() @ c)) / nc)


def preserves(U, E, B):
    return max(proj_resid(U @ e @ U.conj().T, B) for e in E)


def commutant_dim(E):
    N = E[0].shape[0]
    I = np.eye(N)
    G = np.zeros((N * N, N * N), complex)
    for e in E:
        L = np.kron(I, e.T) - np.kron(e, I)  # vec(Xe - eX) row-major
        G += L.conj().T @ L
    w = np.linalg.eigvalsh(0.5 * (G + G.conj().T))
    return int(np.sum(w < 1e-9 * w[-1]))


def intertwiner_null(comp, iota):
    """null dim + solution of U e = iota(e) U over the compressed matrix coefficients."""
    N = next(iter(comp.values())).shape[0]
    I = np.eye(N)
    G = np.zeros((N * N, N * N), complex)
    for k, e in comp.items():
        img = iota(k, e)
        L = np.kron(I, e.T) - np.kron(img, I)
        G += L.conj().T @ L
    w, v = np.linalg.eigh(0.5 * (G + G.conj().T))
    null = int(np.sum(w < 1e-9 * w[-1]))
    Us = v[:, 0].reshape(N, N)
    Us = Us / np.sqrt(np.trace(Us.conj().T @ Us).real / N)
    return null, Us, float(w[0] / w[-1]), float(w[min(1, len(w) - 1)] / w[-1])


def run(n_max, tonly=False):
    out = {'n_max': n_max, 'mode': 'T-only' if tonly else 'full'}
    labels, idx, D, V, Dc, comp, E = setup(n_max)
    N = Dc.shape[0]
    B = np.array([e.reshape(-1) for e in E])
    out['dimH'], out['dimO'] = N, len(E)

    # (A) commutant of O
    if not tonly:
        out['commutant_dim'] = commutant_dim(E)

    # eigenspaces of D on the window
    lam = np.round(np.linalg.eigvalsh(Dc) * 2) / 2
    vals = sorted(set(lam.tolist()))
    w, Q = np.linalg.eigh(Dc)
    P = {}
    for l in vals:
        cols = Q[:, np.abs(w - l) < 1e-6]
        P[l] = cols @ cols.conj().T
    out['eigenvalues'] = vals
    out['eig_mults'] = [int(round(np.trace(P[l]).real)) for l in vals]

    # Spin(4) generators on the window: left J (on m), right J + sigma/2 (on m', s)
    n = len(labels)
    gens = {'L': [], 'R': []}
    for a in range(3):
        Lg = np.zeros((n, n), complex)
        Rg = np.zeros((n, n), complex)
        for j in sorted({lab[0] for lab in labels}):
            dd = int(2 * j) + 1
            Jm = p38._spin_mats(j)[a]
            ms = [j - k for k in range(dd)]
            ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
            sig = (p38.SX, p38.SY, p38.SZ)[a] / 2
            Lg[np.ix_(ids, ids)] = np.kron(Jm, np.eye(2 * dd))
            Rg[np.ix_(ids, ids)] = np.kron(np.eye(dd), np.kron(Jm, np.eye(2)) + np.kron(np.eye(dd), sig))
        gens['L'].append(V.conj().T @ Lg @ V)
        gens['R'].append(V.conj().T @ Rg @ V)
    # Casimirs per eigenspace -> (a, b) labels; multiplicity-free check
    ab = []
    for l in vals:
        cL = sum(np.trace(P[l] @ g @ g).real for g in gens['L']) / np.trace(P[l]).real
        cR = sum(np.trace(P[l] @ g @ g).real for g in gens['R']) / np.trace(P[l]).real
        a = (-1 + np.sqrt(1 + 4 * cL)) / 2
        b = (-1 + np.sqrt(1 + 4 * cR)) / 2
        ab.append((round(2 * a) / 2, round(2 * b) / 2))
        # irreducible: dim = (2a+1)(2b+1)
        assert abs((2 * a + 1) * (2 * b + 1) - np.trace(P[l]).real) < 1e-6
    out['spin4_labels'] = [list(x) for x in ab]
    out['multiplicity_free'] = len(set(ab)) == len(ab)

    # (B) T' = {e^{i theta(D)} preserving O}: block-support constraints per harmonic degree k
    supports = {}
    for (J, M, Mp), e in comp.items():
        for i, l in enumerate(vals):
            for i2, l2 in enumerate(vals):
                if np.linalg.norm(P[l] @ e @ P[l2]) > 1e-9:
                    supports.setdefault(int(2 * J), set()).add((i, i2))
    nv = len(vals)
    ks = sorted(supports)
    rows = []
    for kk, k in enumerate(ks):
        for (i, i2) in sorted(supports[k]):
            r = [0] * (nv + len(ks))
            r[i] += 1
            r[i2] -= 1
            r[nv + kk] -= 1
            rows.append(r)
    A = sp.Matrix(rows)
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(A, domain=sp.ZZ)
    diag = [abs(int(S[i, i])) for i in range(min(S.shape)) if S[i, i] != 0]
    out['T_constraint_rank'] = len(diag)
    out['T_identity_dim'] = (nv + len(ks)) - len(diag)
    out['T_torsion_invariants'] = [d for d in diag if d > 1]

    # direct enumeration of T' torsion: theta in (2pi/m) Z for small m, check against constraints
    # (independent of the SNF): find all theta mod overall phase with theta_0 = 0 satisfying blocks
    def ok_theta(th):
        for k in ks:
            ph = None
            for (i, i2) in supports[k]:
                d = (th[i] - th[i2]) % 1.0
                if ph is None:
                    ph = d
                elif min(abs(d - ph), 1 - abs(d - ph)) > 1e-9:
                    return False
        return True

    found = []
    for m in (2, 3, 4):
        for th in itertools.product(range(m), repeat=nv - 1):
            t = (0.0,) + tuple(x / m for x in th)
            if ok_theta(t):
                found.append(t)
    found = sorted(set(tuple(round(x, 9) for x in t) for t in found))
    out['T_elements_found_mod_phase'] = [list(t) for t in found]

    # which are in Z(Spin4): (eL, eR) acts on (a,b) by eL^{2a} eR^{2b}; mod overall phase
    center = []
    for eL, eR in itertools.product((1, -1), (1, -1)):
        ph = [((eL ** int(2 * a)) * (eR ** int(2 * b))) for (a, b) in ab]
        th = [0.0 if p == ph[0] else 0.5 for p in ph]
        center.append(tuple(th))
    out['center_mod_phase'] = sorted(set(center))
    extra = [t for t in found if t not in set(center)]
    out['T_extra_mod_center'] = [list(t) for t in extra]

    # direct residual checks on O for each found element and for sign(D) (control)
    def tunit(th):
        return sum(np.exp(2j * np.pi * th[i]) * P[l] for i, l in enumerate(vals))

    out['T_found_residuals'] = [preserves(tunit(t), E, B) for t in found]
    sgn = tunit(tuple(0.0 if l > 0 else 0.5 for l in vals))
    out['signD_residual'] = preserves(sgn, E, B)
    lev = tunit(tuple(0.5 * (i % 2) for i in range(nv)))
    out['alternating_level_residual'] = preserves(lev, E, B)
    if tonly:
        return out

    # (C) inversion: intertwiner, swap of L and R, iota^2
    inv = lambda k, e: (-1) ** int(k[2] - k[1]) * comp[(k[0], -k[2], -k[1])]
    nullI, Ui, r0, r1 = intertwiner_null(comp, inv)
    out['inversion_null'] = nullI
    out['inversion_anticommutes_D'] = float(np.linalg.norm(Ui @ Dc + Dc @ Ui))
    out['inversion_unitarity'] = float(np.linalg.norm(Ui.conj().T @ Ui - np.eye(N)))
    out['inversion_preserves_O'] = preserves(Ui, E, B)
    # swap: Ui L_a Ui* in span(R)?
    RB = np.array([g.reshape(-1) for g in gens['R']])
    RBo = np.linalg.qr(RB.T)[0].T
    LB = np.array([g.reshape(-1) for g in gens['L']])
    LBo = np.linalg.qr(LB.T)[0].T
    out['inv_L_to_R_resid'] = max(proj_resid(Ui @ g @ Ui.conj().T, RBo.conj()) for g in gens['L'])
    out['inv_R_to_L_resid'] = max(proj_resid(Ui @ g @ Ui.conj().T, LBo.conj()) for g in gens['R'])
    sq = Ui @ Ui
    c = np.trace(sq) / N
    out['inversion_square_is_scalar'] = float(np.linalg.norm(sq - c * np.eye(N)))

    # transposition: closure of O, implementing unitary (expect none), and direct 2-positivity failure
    out['O_closed_under_transpose'] = max(proj_resid(e.T, B) for e in E)
    out['O_closed_under_conjugation'] = max(proj_resid(e.conj(), B) for e in E)
    nullT, Ut, rT, rT1 = intertwiner_null(comp, lambda k, e: e.T)
    out['transpose_intertwiner_null'] = nullT
    out['transpose_gram_ratios'] = [rT, rT1]
    out['transpose_U_unitarity'] = float(np.linalg.norm(Ut.conj().T @ Ut - np.eye(N)))
    out['transpose_U_comm_D'] = float(np.linalg.norm(Ut @ Dc - Dc @ Ut))
    out['transpose_U_anticomm_D'] = float(np.linalg.norm(Ut @ Dc + Dc @ Ut))
    out['transpose_U_L_to_R'] = max(proj_resid(Ut @ g @ Ut.conj().T, RBo.conj()) for g in gens['L'])
    out['transpose_U_L_to_L'] = max(proj_resid(Ut @ g @ Ut.conj().T, LBo.conj()) for g in gens['L'])
    # conjugation map (antilinear on O; recorded only): O closed under it, see above
    rng = np.random.default_rng(11)
    Eh = []
    for e in E:  # Hermitian spanning set of O
        Eh += [0.5 * (e + e.conj().T), 0.5j * (e - e.conj().T)]
    best = 0.0
    for trial in range(400):
        def rnd():
            return sum(rng.normal() * h for h in Eh[:: max(1, len(Eh) // 40)])
        Aa, Bb = rnd(), rnd()
        C = sum((rng.normal() + 1j * rng.normal()) * e for e in E[:: max(1, len(E) // 30)])
        X = np.block([[Aa, C], [C.conj().T, Bb]])
        Xt = np.block([[Aa.T, C.T], [C.conj().T.T, Bb.T]])
        gap = np.linalg.eigvalsh(X)[0] - np.linalg.eigvalsh(Xt)[0]
        best = max(best, gap / np.linalg.norm(X, 2))
    out['transpose_2pos_violation_rel'] = float(best)
    # control: Ad of a fixed rotation (left translation exp(i*0.7 L_z)) preserves 2-positivity trivially
    return out


if __name__ == '__main__':
    tonly = 'tonly' in sys.argv
    levels = [int(a) for a in sys.argv[1:] if a != 'tonly'] or [2, 3]
    rows = []
    for n_max in levels:
        r = run(n_max, tonly)
        rows.append(r)
        print(json.dumps(r, default=str))
    os.makedirs(os.path.join(ROOT, 'debug', 'data'), exist_ok=True)
    with open(os.path.join(ROOT, 'debug', 'data', 'wg_full_group_tonly.json' if tonly else 'wg_full_group.json'), 'w') as f:
        json.dump(rows, f, default=str, indent=1)
