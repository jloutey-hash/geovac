"""Real structure J + order-zero / order-one on the TRUE Camporesi-Higuchi window.

Measures, on the object Paper 38's convergence theorem is about (complete Dirac
eigenspaces |lambda| <= n_max + 1/2 of the round-S^3 spinor Dirac, operator
system O = P C(S^3) P), the real-structure sign triple, J O J^{-1} = O, and the
order-zero / order-one residuals.  Prior measurements (Paper 32 item 5) were on
the Fock/Weyl MODEL objects only.

Realization (as debug/p38_true_ch_window.py): L^2(SU(2)) (x) C^2, basis
sqrt(2j+1) D^j_{m m'} (x) e_s, D = 2 sigma.J_R + 3/2.

Closed-form J (derivation):
  * spinor charge conjugation on S^3 = SU(2) in the left-invariant frame is
    J psi = (i sigma_y) conj(psi), conj = pointwise complex conjugation.
  * pointwise conjugation of the basis: conj D^j_{m m'} = (-1)^{m-m'} D^j_{-m,-m'}
    (Condon-Shortley; verified numerically below AND independently by sampling).
  * J_R = -i X (X real left-invariant fields) => K J_R K = -J_R;
    (i sigma_y) conj(sigma) (i sigma_y)^{-1} = -sigma  =>  sigma.J_R invariant,
    so J D = + D J  (epsilon' = +1);  (i sigma_y K)^2 = (i sigma_y)^2 = -1
    (epsilon = -1).  KO-dimension 3 row.
  So J = U K with  U e_{j,m,m',s} = (-1)^{m-m'} (i sigma_y)_{s's} e_{j,-m,-m',s'}.

Routes:
  A. Clebsch-Gordan multipliers + closed-form U (this file).
  B. Sampled Wigner-D (no CG, no conj formula): multipliers AND the conjugation
     unitary obtained by least-squares expansion of pointwise-sampled functions.
  N. Numerical solve: all U with U conj(X) = s_X X U for X in {D, J_L^a,
     J_R^a + sigma^a/2} (the Dirac + its isometry generators); solution space =
     one phase per irreducible block; J^2 is phase-independent.

Output: debug/data/true_window_real_structure.json
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import time
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tests"))
from test_p38_true_ch_window import _build, _mult_op, _spin_mats, SX, SY, SZ  # noqa: E402

ISY = np.array([[0, 1], [-1, 0]], complex)  # i sigma_y
TOL = 1e-10


# --------------------------------------------------------------------------
# Route A
# --------------------------------------------------------------------------
def closed_form_U(labels, idx) -> np.ndarray:
    n = len(labels)
    U = np.zeros((n, n), complex)
    for (j, m, mp, s) in labels:
        ph = (-1) ** int(m - mp)
        for s2 in (0, 1):
            c = ISY[s2, s]
            if c != 0:
                U[idx[(j, -m, -mp, s2)], idx[(j, m, mp, s)]] = ph * c
    return U


def ambient_generators(labels, idx):
    """J_L^a (on m), K^a = J_R^a + sigma^a/2 (on m', s): isometry generators."""
    n = len(labels)
    JL = [np.zeros((n, n), complex) for _ in range(3)]
    KR = [np.zeros((n, n), complex) for _ in range(3)]
    JR = [np.zeros((n, n), complex) for _ in range(3)]
    js = sorted({l[0] for l in labels})
    paul = (SX, SY, SZ)
    for j in js:
        d = int(2 * j) + 1
        S = _spin_mats(j)
        ms = [j - k for k in range(d)]
        ids = [idx[(j, m, mp, s)] for m in ms for mp in ms for s in (0, 1)]
        I_d, I_2 = np.eye(d), np.eye(2)
        for a in range(3):
            JL[a][np.ix_(ids, ids)] = np.kron(np.kron(S[a], I_d), I_2)
            jr = np.kron(np.kron(I_d, S[a]), I_2)
            JR[a][np.ix_(ids, ids)] = jr
            KR[a][np.ix_(ids, ids)] = jr + np.kron(np.kron(I_d, I_d), paul[a]) / 2
    return JL, JR, KR


def multipliers(n_max: int, labels, idx):
    out = []
    for N in range(1, 2 * n_max + 1):
        J = Fraction(N - 1, 2)
        Ms = [J - k for k in range(N)]
        for M, Mp in itertools.product(Ms, Ms):
            out.append(((N, M, Mp), _mult_op(J, M, Mp, labels, idx)))
    return out


def order_stats(Ts: np.ndarray, Bs: np.ndarray, V: np.ndarray, Dc: np.ndarray | None):
    """For all pairs (i, j): C = [X_i, B_j], X_i = T_i (order 0) or [Dc, T_i] (order 1).
    max-abs taken in the AMBIENT Peter-Weyl basis (V C V^dag), op norm basis-free."""
    n = len(Ts)
    Xs = Ts if Dc is None else (Dc @ Ts - Ts @ Dc)
    maxabs = np.zeros((n, n))
    opn = np.zeros((n, n))
    for i in range(n):
        C = Xs[i] @ Bs - Bs @ Xs[i]
        Camb = V @ C @ V.conj().T
        maxabs[i] = np.abs(Camb).reshape(n, -1).max(axis=1)
        opn[i] = np.linalg.norm(C, ord=2, axis=(1, 2))
    fail = maxabs > TOL
    return {
        "pairs": int(n * n),
        "failing": int(fail.sum()),
        "failing_fraction": float(fail.mean()),
        "max_abs_ambient": float(maxabs.max()),
        "max_op_norm": float(opn.max()),
    }, maxabs, opn


def route_a(n_max: int, fixed_only: bool = False) -> Dict:
    t0 = time.time()
    labels, idx, D, ev, U_eig = _build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = U_eig[:, keep]
    d = V.shape[1]
    Dc = V.conj().T @ D @ V
    U = closed_form_U(labels, idx)
    res: Dict = {"n_max": n_max, "dim_window": d}
    # ---- sign triple, ambient ----
    res["unitary_err"] = float(np.abs(U @ U.conj().T - np.eye(len(labels))).max())
    res["J2_plus1_err"] = float(np.abs(U @ U.conj() - np.eye(len(labels))).max())
    res["J2_minus1_err"] = float(np.abs(U @ U.conj() + np.eye(len(labels))).max())
    res["JD_plus_DJ_err"] = float(np.abs(U @ D.conj() - D @ U).max())
    res["JD_minus_DJ_err"] = float(np.abs(U @ D.conj() + D @ U).max())
    JL, JR, KR = ambient_generators(labels, idx)
    res["J_JR_antic_err"] = float(max(np.abs(U @ X.conj() + X @ U).max() for X in JR))
    res["J_JL_antic_err"] = float(max(np.abs(U @ X.conj() + X @ U).max() for X in JL))
    # ---- window preservation ----
    W = U @ V.conj()
    res["window_preserved_err"] = float(np.abs(V @ (V.conj().T @ W) - W).max())
    Jw = V.conj().T @ W  # J on window: psi -> Jw conj(psi)
    res["Jw_unitary_err"] = float(np.abs(Jw @ Jw.conj().T - np.eye(d)).max())
    res["Jw_J2_minus1_err"] = float(np.abs(Jw @ Jw.conj() + np.eye(d)).max())
    res["Jw_JD_plus_err"] = float(np.abs(Jw @ Dc.conj() - Dc @ Jw).max())
    # ---- multipliers ----
    mults = multipliers(n_max, labels, idx)
    labs = [l for l, _ in mults]
    lab_idx = {l: i for i, l in enumerate(labs)}
    Amb = np.array([A for _, A in mults])
    # ambient conjugation identity  U conj(A_{J M M'}) U^dag = (-1)^{M-M'} A_{J,-M,-M'}
    err_conj = 0.0
    for (N, M, Mp), A in mults:
        B = Amb[lab_idx[(N, -M, -Mp)]] * (-1) ** int(M - Mp)
        err_conj = max(err_conj, float(np.abs(U @ A.conj() @ U.conj().T - B).max()))
    res["ambient_conj_multiplier_err"] = err_conj
    Ts = V.conj().T @ Amb @ V  # V^dag A V
    flat = Ts.reshape(len(Ts), -1)
    dimO = int(np.linalg.matrix_rank(flat, tol=1e-9))
    res["dim_O"] = dimO
    # J O J^{-1}
    JTJ = Jw @ Ts.conj() @ Jw.conj().T
    _, _, Vh = np.linalg.svd(flat, full_matrices=False)
    Q = Vh[:dimO]
    proj = (JTJ.reshape(len(Ts), -1) @ Q.conj().T) @ Q
    res["JOJ_outside_O_err"] = float(np.abs(proj - JTJ.reshape(len(Ts), -1)).max())
    err_lab = 0.0
    for i, (N, M, Mp) in enumerate(labs):
        err_lab = max(err_lab, float(np.abs(JTJ[i] - (-1) ** int(M - Mp) * Ts[lab_idx[(N, -M, -Mp)]]).max()))
    res["JTJ_label_identity_err"] = err_lab
    # J b^* J^{-1} and the Toeplitz identity J b^* J^{-1} = b (b = P f P)
    Bs = Jw @ np.transpose(Ts, (0, 2, 1)) @ Jw.conj().T  # conj(T^dag)=T^T
    res["J_bstar_J_equals_b_err"] = float(np.abs(Bs - Ts).max())
    # fixed generators N = 2
    gi = [i for i, l in enumerate(labs) if l[0] == 2]
    fz, _, opz = order_stats(Ts[gi], Bs[gi], V, None)
    fo, _, opo = order_stats(Ts[gi], Bs[gi], V, Dc)
    gen_norm = float(max(np.linalg.norm(Ts[i], 2) for i in gi))
    res["fixed_N2"] = {"order0_op": fz["max_op_norm"], "order0_maxabs": fz["max_abs_ambient"],
                       "order1_op": fo["max_op_norm"], "order1_maxabs": fo["max_abs_ambient"],
                       "gen_op_norm": gen_norm, "commD_gen_op_norm":
                       float(max(np.linalg.norm(Dc @ Ts[i] - Ts[i] @ Dc, 2) for i in gi))}
    if not fixed_only:
        s0, m0, o0 = order_stats(Ts, Bs, V, None)
        s1, m1, o1 = order_stats(Ts, Bs, V, Dc)
        # model convention [a, J b J^{-1}]: same set of commutators up to relabeling
        BsJ = Jw @ Ts.conj() @ Jw.conj().T
        s0b, _, _ = order_stats(Ts, BsJ, V, None)
        res["order0"] = s0
        res["order1"] = s1
        res["order0_JbJ_convention"] = s0b
        # restricted-evaluation: which pairs vanish?  band content of zero pairs
        zero_pairs = [(labs[i], labs[j]) for i in range(len(labs)) for j in range(len(labs))
                      if m0[i, j] <= TOL]
        res["order0_zero_pairs_involving_nonscalar"] = int(sum(
            1 for a, b in zero_pairs if a[0] > 1 and b[0] > 1))
        res["_opnorms"] = {"o0": o0, "o1": o1, "m0": m0, "m1": m1}
    res["seconds"] = time.time() - t0
    return res


# --------------------------------------------------------------------------
# Route N: numerical solve for J from D + isometry equivariance
# --------------------------------------------------------------------------
def route_n(n_max: int) -> Dict:
    labels, idx, D, ev, U_eig = _build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = U_eig[:, keep]
    d = V.shape[1]
    JL, JR, KR = ambient_generators(labels, idx)
    comp = lambda X: V.conj().T @ X @ V
    Dc = comp(D)
    gens = [comp(X) for X in JL + KR]
    out = {"n_max": n_max}
    for sD in (+1, -1):
        # linear map U -> {U conj(X) - s X U}; vec column-major
        rows = []
        I = np.eye(d)
        for X, s in [(Dc, sD)] + [(G, -1) for G in gens]:
            # vec(U A) = (A^T kron I) vec U ; vec(B U) = (I kron B) vec U
            rows.append(np.kron(X.conj().T, I) - s * np.kron(I, X))
        M = np.vstack(rows)
        sv = np.linalg.svd(M, compute_uv=False)
        null = int((sv < 1e-9 * max(1.0, sv.max())).sum()) + max(0, d * d - len(sv))
        out[f"solution_dim_epsprime_{'+' if sD > 0 else '-'}"] = null
        if sD > 0:
            _, _, Vh = np.linalg.svd(M, full_matrices=False)
            sols = Vh[-null:].conj()
            # random combination -> generic solution; normalize to unitary blockwise via polar
            rng = np.random.default_rng(1)
            Ug = (rng.standard_normal(null) + 1j * rng.standard_normal(null)) @ sols
            Ug = Ug.reshape(d, d, order="F")
            u, _, vh = np.linalg.svd(Ug)
            Uu = u @ vh
            out["generic_solution_satisfies"] = float(max(
                np.abs(Uu @ X.conj() - s * X @ Uu).max()
                for X, s in [(Dc, 1)] + [(G, -1) for G in gens]))
            out["generic_J2_minus1_err"] = float(np.abs(Uu @ Uu.conj() + np.eye(d)).max())
            out["generic_J2_plus1_err"] = float(np.abs(Uu @ Uu.conj() - np.eye(d)).max())
            # closed form lies in the solution family
            Jw = V.conj().T @ closed_form_U(labels, idx) @ V.conj()
            v = Jw.reshape(-1, order="F")
            coef = sols.conj() @ v
            out["closed_form_in_family_err"] = float(np.abs(sols.T @ coef - v).max())
    out["n_irreducible_blocks"] = 2 * n_max
    return out


# --------------------------------------------------------------------------
# Route B: sampled Wigner-D; conjugation by sampling (no CG, no conj formula)
# --------------------------------------------------------------------------
_PAULI = (SX, SY, SZ)


def _b_spin(j2: int):
    j = j2 / 2
    ms = [j - k for k in range(j2 + 1)]
    jz = np.diag(ms).astype(complex)
    jp = np.zeros((j2 + 1, j2 + 1), complex)
    for a, m in enumerate(ms):
        if a > 0:
            jp[a - 1, a] = np.sqrt(j * (j + 1) - m * (m + 1))
    return [(jp + jp.conj().T) / 2, (jp - jp.conj().T) / (2j), jz]


def _b_wigner(j2: int, g: np.ndarray) -> np.ndarray:
    from scipy.linalg import expm
    w, Vv = np.linalg.eig(g)
    H = Vv @ np.diag(-np.angle(w)) @ np.linalg.inv(Vv)
    c = [np.trace(H @ s_).real for s_ in _PAULI]
    J = _b_spin(j2)
    return expm(-1j * sum(ci * Ji for ci, Ji in zip(c, J)))


def route_b(n_max: int, n_samples: int = 1600) -> Dict:
    rng = np.random.default_rng(11)
    gs = []
    for _ in range(n_samples):
        q = rng.standard_normal(4)
        q /= np.linalg.norm(q)
        gs.append(np.array([[q[0] + 1j * q[3], q[2] + 1j * q[1]],
                            [-q[2] + 1j * q[1], q[0] - 1j * q[3]]]))
    cols, labs, Dcache = [], [], {}
    for j2 in range(3 * n_max + 2):
        Ds = np.array([_b_wigner(j2, g) for g in gs])
        Dcache[j2] = Ds
        for a in range(j2 + 1):
            for b in range(j2 + 1):
                cols.append(np.sqrt(j2 + 1) * Ds[:, a, b])
                labs.append((j2, a, b))
    F = np.array(cols).T
    Fp = np.linalg.pinv(F)
    H_idx = [i for i, l in enumerate(labs) if l[0] <= n_max]
    Hlabs = [labs[i] for i in H_idx]
    nf = len(Hlabs)
    full = [(l, s) for l in Hlabs for s in (0, 1)]
    pos = {f: i for i, f in enumerate(full)}
    Dop = np.zeros((len(full), len(full)), complex)
    for (j2, a, b), s1 in full:
        J = _b_spin(j2)
        for (j2p, a2, b2), s2 in full:
            if j2p != j2 or a2 != a:
                continue
            val = sum(2 * J[k][b2, b] * _PAULI[k][s2, s1] for k in range(3))
            if b2 == b and s2 == s1:
                val += 1.5
            Dop[pos[((j2p, a2, b2), s2)], pos[((j2, a, b), s1)]] = val
    ev, Ue = np.linalg.eigh(Dop)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = Ue[:, keep]
    d = V.shape[1]
    Dc = V.conj().T @ Dop @ V
    # conjugation unitary on functions, by sampling: conj(e_i) expanded
    Uf = np.zeros((nf, nf), complex)
    for c, i in enumerate(H_idx):
        coef = Fp @ np.conj(F[:, i])
        Uf[:, c] = coef[H_idx]
    Ufull = np.kron(Uf, ISY)
    out = {"n_max": n_max, "dim_window": d}
    out["Uf_unitary_err"] = float(np.abs(Uf @ Uf.conj().T - np.eye(nf)).max())
    out["J2_minus1_err"] = float(np.abs(Ufull @ Ufull.conj() + np.eye(2 * nf)).max())
    out["JD_plus_DJ_err"] = float(np.abs(Ufull @ Dop.conj() - Dop @ Ufull).max())
    W = Ufull @ V.conj()
    out["window_preserved_err"] = float(np.abs(V @ (V.conj().T @ W) - W).max())
    Jw = V.conj().T @ W
    Ts = []
    for N in range(1, 2 * n_max + 1):
        Ds = Dcache[N - 1]
        for M, Mp in itertools.product(range(N), range(N)):
            fvals = Ds[:, M, Mp]
            A = np.zeros((nf, nf), complex)
            for c, i in enumerate(H_idx):
                A[:, c] = (Fp @ (fvals * F[:, i]))[H_idx]
            Ts.append(V.conj().T @ np.kron(A, np.eye(2)) @ V)
    Ts = np.array(Ts)
    flat = Ts.reshape(len(Ts), -1)
    dimO = int(np.linalg.matrix_rank(flat, tol=1e-7))
    out["dim_O"] = dimO
    JTJ = Jw @ Ts.conj() @ Jw.conj().T
    _, _, Vh = np.linalg.svd(flat, full_matrices=False)
    Q = Vh[:dimO]
    proj = (JTJ.reshape(len(Ts), -1) @ Q.conj().T) @ Q
    out["JOJ_outside_O_err"] = float(np.abs(proj - JTJ.reshape(len(Ts), -1)).max())
    Bs = Jw @ np.transpose(Ts, (0, 2, 1)) @ Jw.conj().T
    s0, m0, o0 = order_stats(Ts, Bs, V, None)
    s1, m1, o1 = order_stats(Ts, Bs, V, Dc)
    out["order0"] = s0
    out["order1"] = s1
    out["_opnorms"] = {"o0": o0, "o1": o1, "m0": m0, "m1": m1}
    return out


# --------------------------------------------------------------------------
# Mechanism: the residuals are Toeplitz boundary terms of the compression
# --------------------------------------------------------------------------
def boundary_mechanism(n_max: int) -> Dict:
    """[PfP, PgP] = P g (1-P) f P - P f (1-P) g P (fg = gf), and likewise for
    [[D, PfP], PgP] with [D, f] = c(df) pointwise commuting with g.  A band-N
    multiplier (spin J = (N-1)/2 on both SU(2)_L and the total right spin K,
    Dirac level n <-> (j_L, K) = (n/2, (n+1)/2) or ((n+1)/2, n/2)) moves the
    level n (|lambda| = n + 3/2) by at most N - 1, so the residual must vanish
    on every row/column at level n <= top - (max(N_a, N_b) - 1), top = n_max - 1.
    For the degree-1 generators (N = 2) that confines it to the TOP level.  Measured: the largest residual entry outside that top
    boundary strip, the normalized Hilbert-Schmidt norm ||C||_F / sqrt(dim H)
    of the degree-1 residual, and the strip's share of dim H."""
    labels, idx, D, ev, U_eig = _build(n_max)
    keep = np.abs(ev) <= n_max + 0.5 + 1e-9
    V = U_eig[:, keep]
    lev = np.rint(np.abs(ev[keep]) - 1.5).astype(int)
    top = n_max - 1
    Dc = V.conj().T @ D @ V
    Uc = closed_form_U(labels, idx)
    Jw = V.conj().T @ Uc @ V.conj()
    mults = multipliers(n_max, labels, idx)
    Ns = np.array([l[0] for l, _ in mults])
    Ts = V.conj().T @ np.array([A for _, A in mults]) @ V
    Bs = Jw @ np.transpose(Ts, (0, 2, 1)) @ Jw.conj().T
    # measured level reach of each multiplier on the WHOLE ambient space
    # (all Dirac levels present there, including those above the window)
    lev_all = np.rint(np.abs(ev) - 1.5).astype(int)
    reach = []
    for _, A in mults:
        T = U_eig.conj().T @ A @ U_eig
        nz = np.argwhere(np.abs(T) > 1e-12)
        reach.append(int(np.abs(lev_all[nz[:, 0]] - lev_all[nz[:, 1]]).max()) if len(nz) else 0)
    reach_ok = bool(all(r <= N - 1 for r, N in zip(reach, Ns)))
    out = {"n_max": n_max, "reach_le_band_minus_1": reach_ok,
           "max_reach_by_band": {int(N): int(max(r for r, M in zip(reach, Ns) if M == N))
                                 for N in sorted(set(Ns.tolist()))}}
    worst0 = worst1 = 0.0
    DT = Dc @ Ts - Ts @ Dc
    for i in range(len(Ts)):
        for order, X in ((0, Ts[i]), (1, DT[i])):
            C = X @ Bs - Bs @ X
            for j in range(len(Ts)):
                R = max(Ns[i], Ns[j]) - 1
                inner = lev <= top - R
                if not inner.any():
                    continue
                v = max(np.abs(C[j][inner, :]).max(), np.abs(C[j][:, inner]).max())
                if order == 0:
                    worst0 = max(worst0, float(v))
                else:
                    worst1 = max(worst1, float(v))
    out["order0_max_outside_boundary_strip"] = worst0
    out["order1_max_outside_boundary_strip"] = worst1
    gi = [i for i in range(len(Ts)) if Ns[i] == 2]
    hs0 = max(np.linalg.norm(Ts[a] @ Bs[b] - Bs[b] @ Ts[a]) for a in gi for b in gi)
    hs1 = max(np.linalg.norm(DT[a] @ Bs[b] - Bs[b] @ DT[a]) for a in gi for b in gi)
    d = len(lev)
    out["N2_order0_HS_over_sqrt_dim"] = float(hs0 / np.sqrt(d))
    out["N2_order1_HS_over_sqrt_dim"] = float(hs1 / np.sqrt(d))
    out["N2_strip_share_of_dim"] = float(np.mean(lev == top))
    return out


def main() -> None:
    results: Dict = {"route_A": {}, "route_N": {}, "route_B": {}, "cross": {}, "fixed_N2": {}}
    for n in (1, 2, 3, 4):
        r = route_a(n)
        on = r.pop("_opnorms")
        results["route_A"][n] = r
        print(f"[A] n_max={n}: {json.dumps({k: v for k, v in r.items() if k != 'fixed_N2'}, default=float)}")
        results["fixed_N2"][n] = r["fixed_N2"]
        if n <= 3:
            rn = route_n(n)
            results["route_N"][n] = rn
            print(f"[N] n_max={n}: {rn}")
        if n <= 3:
            rb = route_b(n)
            onb = rb.pop("_opnorms")
            results["route_B"][n] = rb
            print(f"[B] n_max={n}: {rb}")
            results["cross"][n] = {k: float(np.abs(on[k] - onb[k]).max()) for k in on}
            print(f"[A vs B] n_max={n}: {results['cross'][n]}")
    results["mechanism"] = {}
    for n in (1, 2, 3, 4, 5):
        mb = boundary_mechanism(n)
        results["mechanism"][n] = mb
        print(f"[mech] n_max={n}: {mb}")
    for n in (5, 6):
        r = route_a(n, fixed_only=True)
        results["fixed_N2"][n] = r["fixed_N2"]
        print(f"[A fixed] n_max={n}: {r['fixed_N2']}")
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    with open(os.path.join(HERE, "data", "true_window_real_structure.json"), "w") as fh:
        json.dump(results, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
