"""Diagnostic: does strong orthogonality (Q12) + 1e orthonormality separate the two
kappa(S) multipliers of the He R12-CI?  (debug/ only; READ-ONLY use of the engine.)

Three questions, all answerable from the ALREADY-ASSEMBLED (H, M):

  D1  Q12 exactly, no RI.  The "occupied x occupied product space" is literally the
      span of the orbital-pair basis functions, so the strong-orthogonality projection
      is a congruence transform of (H, M) -- no new integrals.  Measure kappa before/after.
      SELF-VALIDATION: the projection is a change of basis WITHIN THE SAME SPAN, so the
      variational energy must be invariant to machine precision.  If E moves, the
      implementation is wrong.

  D2  How much of the geminal is genuinely new?  Schur complement of the geminal block.

  D3  Mechanism for part 2: is the pair-block kappa inherited from the ONE-ELECTRON
      shared-k Sturmian overlap?  And is that 1e block orthonormal in the Sturmian
      (1/r-weighted) metric?  Plus: for a single center, is per-l Lowdin sparsity-free?
"""
import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

_spec = importlib.util.spec_from_file_location(
    "ctf12_r12ci_he", os.path.join(HERE, "ctf12_r12ci_he.py"))
ctf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ctf)

EXACT_HE = ctf.EXACT_HE
NG, NX = 600, 128
KMIN = 1.30
R_MAX = 42.0 / KMIN
_RGRID = ctf.make_grid(KMIN, NG, r_max=R_MAX)

# recorded accurate operating points (debug/data/r12ci_encoding_cost.json)
OPS = [dict(ns=3, n_gem=1, k=1.7, gamma=0.7, ref_kappa=176.94, ref_err=0.798),
       dict(ns=3, n_gem=2, k=1.7, gamma=0.7, ref_kappa=250.97, ref_err=0.640),
       dict(ns=4, n_gem=1, k=1.7, gamma=0.7, ref_kappa=189.0,  ref_err=0.474),
       dict(ns=4, n_gem=2, k=1.7, gamma=0.7, ref_kappa=236.0,  ref_err=0.449)]


def build_bfs(ns, n_gem, Rtab, dRtab):
    bfs = []
    for i in range(1, ns + 1):
        for j in range(i, ns + 1):
            bfs.append(ctf.orbital_pair(i, j, Rtab, dRtab))
    for ref in range(1, n_gem + 1):
        bfs.append(ctf.geminal(ref, Rtab, dRtab))
    return bfs


def overlap_and_H(ns, n_gem, k, gamma):
    r, wr = _RGRID
    Rtab, dRtab = ctf.build_tabs(max(ns, n_gem, 1), r, float(k))
    K = ctf.build_kernels(r, gamma, nx=NX)
    bfs = build_bfs(ns, n_gem, Rtab, dRtab)
    H, M = ctf.assemble(bfs, r, wr, K, Z=2)
    return H, M


def kappa_norm(M):
    """cond of the unit-diagonal correlation matrix (the NOQE-relevant object)."""
    d = np.sqrt(np.diag(M))
    S = M / np.outer(d, d)
    ev = np.linalg.eigvalsh(S)
    return float(ev.max() / ev.min()), float(ev.min())


def energy(H, M, thr=1e-10):
    ev, U = np.linalg.eigh(M)
    keep = ev > thr
    X = U[:, keep] / np.sqrt(ev[keep])
    return float(np.linalg.eigvalsh(X.T @ H @ X)[0])


def q12_project(H, M, n_pair):
    """Exact strong-orthogonality projection of the geminal block against the
    orbital-pair space.  T = [[I, -M_PP^-1 M_PG],[0, I]];  H'=T^T H T, M'=T^T M T."""
    n = M.shape[0]
    P, G = slice(0, n_pair), slice(n_pair, n)
    C = np.linalg.solve(M[P, P], M[P, G])       # pair-space coefficients of the geminals
    T = np.eye(n)
    T[P, G] = -C
    return T.T @ H @ T, T.T @ M @ T, C


out = {"system": "He", "exact": EXACT_HE, "grid": dict(Ng=NG, nx=NX, kmin=KMIN),
       "D1_q12": [], "D2_new_content": [], "D3_mechanism": {}}

print("=" * 78)
print("D1 -- exact Q12 strong orthogonality (congruence transform, no new integrals)")
print("=" * 78)
print(f"{'ns':>3}{'ngem':>5}{'E_raw':>13}{'E_Q12':>13}{'|dE| Ha':>11}"
      f"{'kap_raw':>10}{'kap_Q12':>10}{'ratio':>8}{'offdiag':>10}")
for op in OPS:
    H, M = overlap_and_H(op["ns"], op["n_gem"], op["k"], op["gamma"])
    n_pair = M.shape[0] - op["n_gem"]
    E0 = energy(H, M)
    k0, lmin0 = kappa_norm(M)
    Hp, Mp, C = q12_project(H, M, n_pair)
    E1 = energy(Hp, Mp)
    k1, lmin1 = kappa_norm(Mp)
    offd = float(np.abs(Mp[0:n_pair, n_pair:]).max())
    print(f"{op['ns']:>3}{op['n_gem']:>5}{E0:>13.9f}{E1:>13.9f}{abs(E1-E0):>11.2e}"
          f"{k0:>10.1f}{k1:>10.1f}{k0/k1:>8.2f}{offd:>10.1e}")
    out["D1_q12"].append(dict(ns=op["ns"], n_gem=op["n_gem"], E_raw=E0, E_q12=E1,
                              dE=abs(E1 - E0), kappa_raw=k0, kappa_q12=k1,
                              ratio=k0 / k1, max_offdiag_PG=offd,
                              err_mHa=(E0 - EXACT_HE) * 1000,
                              ref_kappa=op["ref_kappa"], ref_err_mHa=op["ref_err"]))

print()
print("=" * 78)
print("D2 -- how much of the geminal is genuinely OUTSIDE the pair space?")
print("=" * 78)
print(f"{'ns':>3}{'ngem':>5}{'gem':>5}{'||new||^2/||G||^2':>20}{'collinearity':>14}")
for op in OPS:
    H, M = overlap_and_H(op["ns"], op["n_gem"], op["k"], op["gamma"])
    n_pair = M.shape[0] - op["n_gem"]
    P, G = slice(0, n_pair), slice(n_pair, M.shape[0])
    C = np.linalg.solve(M[P, P], M[P, G])
    schur = M[G, G] - M[G, P] @ C
    for g in range(op["n_gem"]):
        frac = float(schur[g, g] / M[n_pair + g, n_pair + g])
        print(f"{op['ns']:>3}{op['n_gem']:>5}{g+1:>5}{frac:>20.4f}{np.sqrt(1-frac):>14.4f}")
        out["D2_new_content"].append(dict(ns=op["ns"], n_gem=op["n_gem"], gem=g + 1,
                                          new_frac=frac,
                                          collinearity=float(np.sqrt(1 - frac))))

print()
print("=" * 78)
print("D3 -- mechanism: is the pair-block kappa inherited from the 1e Sturmian block?")
print("=" * 78)
r, wr = _RGRID
W = r * r * wr
k = 1.7
print(f"{'ns':>3} | {'kap(S_1e) L2':>13}{'kap(S_1e)^2':>13}{'kap(pair) meas':>16}"
      f" | {'max|Sturm-I|':>13}{'kap per-l Lowdin':>17}")
for ns in range(2, 7):
    Rtab, dRtab = ctf.build_tabs(ns, r, k)
    Rm = np.array([Rtab[n] for n in range(1, ns + 1)])
    S1 = (Rm * W) @ Rm.T                                   # L2:      int R_n R_m r^2 dr
    d1 = np.sqrt(np.diag(S1)); S1n = S1 / np.outer(d1, d1)
    Ssturm = (Rm * (r * wr)) @ Rm.T                        # 1/r-wt:  int R_n R_m r   dr
    ds = np.sqrt(np.diag(Ssturm)); Ssn = Ssturm / np.outer(ds, ds)
    ev1 = np.linalg.eigvalsh(S1n)
    kap1 = ev1.max() / ev1.min()
    # measured pair-block kappa (n_gem = 0)
    _, Mp0 = overlap_and_H(ns, 0, k, 0.7)
    kap_pair, _ = kappa_norm(Mp0)
    sturm_dev = float(np.abs(Ssn - np.eye(ns)).max())
    # per-l Lowdin on the 1e block -> pair block becomes orthonormal by construction
    w1, V1 = np.linalg.eigh(S1n)
    Xlow = V1 @ np.diag(w1 ** -0.5) @ V1.T
    kap_after = float(np.linalg.cond(Xlow.T @ S1n @ Xlow))
    print(f"{ns:>3} | {kap1:>13.2f}{kap1**2:>13.2f}{kap_pair:>16.2f}"
          f" | {sturm_dev:>13.2e}{kap_after:>17.4f}")
    out["D3_mechanism"][str(ns)] = dict(kappa_S1e_L2=float(kap1),
                                        kappa_S1e_sq=float(kap1 ** 2),
                                        kappa_pair_measured=float(kap_pair),
                                        max_dev_sturmian_metric_from_I=sturm_dev,
                                        kappa_after_lowdin=kap_after)

os.makedirs("debug/data", exist_ok=True)
with open("debug/data/r12ci_strong_orthogonality_probe.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/r12ci_strong_orthogonality_probe.json")


# ---------------------------------------------------------------------------
# D4 -- THE DECIDER.  Q12 lowers kappa(S), but the projected geminal is a LINEAR
# COMBINATION  |G~> = |G> - sum_p c_p |pair_p>  of the original (preparable) basis
# functions.  On a device that combination costs an LCU with 1-norm ||(1,-c)||_1,
# normalized to the new function's own norm.  If the 1-norm penalty >= the kappa
# gain, this is the corpus's cost-conservation wall in a 5th currency, not a win.
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("D4 -- did the cost move, or did it go away?  (LCU 1-norm of the Q12 geminal)")
print("=" * 78)
print(f"{'ns':>3}{'ngem':>5}{'gem':>5}{'||c||_1':>10}{'norm_ratio':>12}"
       f"{'lcu_1norm':>11}{'kappa_gain':>12}{'net':>9}")
out["D4_preparation_cost"] = []
for op in OPS:
    H, M = overlap_and_H(op["ns"], op["n_gem"], op["k"], op["gamma"])
    n = M.shape[0]; n_pair = n - op["n_gem"]
    P, G = slice(0, n_pair), slice(n_pair, n)
    C = np.linalg.solve(M[P, P], M[P, G])
    schur = M[G, G] - M[G, P] @ C
    k0, _ = kappa_norm(M)
    Hp, Mp, _ = q12_project(H, M, n_pair)
    k1, _ = kappa_norm(Mp)
    # normalize each ORIGINAL basis function to unit norm, then express |G~>
    d = np.sqrt(np.diag(M))
    for g in range(op["n_gem"]):
        gi = n_pair + g
        # |G~> = |G> - sum_p C[p,g] |pair_p>   in UNIT-NORMALIZED function amplitudes
        amps = np.concatenate(([d[gi]], -C[:, g] * d[:n_pair]))
        l1 = float(np.abs(amps).sum())
        new_norm = float(np.sqrt(schur[g, g]))
        lcu = l1 / new_norm                       # 1-norm per unit of prepared state
        gain = k0 / k1
        print(f"{op['ns']:>3}{op['n_gem']:>5}{g+1:>5}"
              f"{float(np.abs(C[:, g]).sum()):>10.3f}{new_norm/d[gi]:>12.4f}"
              f"{lcu:>11.2f}{gain:>12.2f}{gain/lcu:>9.2f}")
        out["D4_preparation_cost"].append(
            dict(ns=op["ns"], n_gem=op["n_gem"], gem=g + 1,
                 c_1norm=float(np.abs(C[:, g]).sum()),
                 norm_ratio=new_norm / d[gi], lcu_1norm=lcu,
                 kappa_gain=gain, net=gain / lcu))

with open("debug/data/r12ci_strong_orthogonality_probe.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nnet > 1  =>  genuine win;  net ~ 1  =>  cost conserved (wall in a new currency)")
