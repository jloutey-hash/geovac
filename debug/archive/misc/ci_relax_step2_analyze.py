r"""Combine LM block matrices from debug/data/cirelax_<ref>_lm_s*.npz; solve the generalized
eigenproblem for ansatz (C) frozen 2x2, (A) CI-free + F Psi_CI, (B) CI-free + F Phi_I; jackknife
over blocks for errors on E and on the correlated differences.  Writes the optimum vectors (in the
2N {Phi_I, F Phi_I} basis) to debug/data/cirelax_<ref>_vecs.npz for the re-sampling VMC step.
Usage: python debug/ci_relax_step2_analyze.py 9det [tol]
"""
from __future__ import annotations

import os
import sys
import glob
import numpy as np
import scipy.linalg as sl

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def projectors(cvec):
    N = len(cvec); Z = np.zeros(N, complex)
    QC = np.column_stack([np.concatenate([cvec, Z]), np.concatenate([Z, cvec])])
    QA = np.zeros((2 * N, N + 1), complex); QA[:N, :N] = np.eye(N); QA[N:, N] = cvec
    QB = np.eye(2 * N, dtype=complex)
    Q0 = np.concatenate([cvec, Z])[:, None]
    return {"0": Q0, "C": QC, "A": QA, "B": QB}


def solve(S, H, Q, tol):
    Ss = Q.conj().T @ S @ Q; Hs = Q.conj().T @ H @ Q
    Ss = 0.5 * (Ss + Ss.conj().T)
    lam, U = np.linalg.eigh(Ss)
    keep = lam > tol * lam[-1]
    Uk = U[:, keep] / np.sqrt(lam[keep])
    Hp = Uk.conj().T @ Hs @ Uk
    ev, Y = sl.eig(Hp)
    ok = np.abs(ev.imag) < 1e-6 + 1e-3 * np.abs(ev.real - ev.real.min())
    idx = np.argmin(np.where(ok, ev.real, np.inf))
    x = Uk @ Y[:, idx]
    return float(ev[idx].real), Q @ x, int(keep.sum())


def main():
    ref = sys.argv[1] if len(sys.argv) > 1 else "9det"
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-10
    files = sorted(glob.glob(os.path.join(DATA, f"cirelax_{ref}_lm_s*.npz")))
    Sb, Hb = [], []
    for f in files:
        z = np.load(f, allow_pickle=True)
        Sb.append(z["S"]); Hb.append(z["H"]); cvec = z["cvec"]; E_trunc = float(z["E_trunc"])
    Sb = np.concatenate(Sb); Hb = np.concatenate(Hb); nb = len(Sb)
    # re-block to ~50 jackknife groups
    ng = min(50, nb); per = nb // ng
    Sg = Sb[:ng * per].reshape(ng, per, *Sb.shape[1:]).mean(1)
    Hg = Hb[:ng * per].reshape(ng, per, *Hb.shape[1:]).mean(1)
    Q = projectors(cvec)
    S = Sg.mean(0); H = Hg.mean(0)
    res = {}; vecs = {}
    for k in ("0", "C", "A", "B"):
        E, v, nk = solve(S, H, Q[k], tol)
        res[k] = E; vecs[k] = v
        print(f"  {k}: E = {E:.6f}   (kept {nk}/{Q[k].shape[1]})")
    jk = {k: [] for k in res}
    for g in range(ng):
        Sj = (S * ng - Sg[g]) / (ng - 1); Hj = (H * ng - Hg[g]) / (ng - 1)
        for k in res:
            jk[k].append(solve(Sj, Hj, Q[k], tol)[0])
    jk = {k: np.array(v) for k, v in jk.items()}
    def jerr(a):
        return np.sqrt((ng - 1) / ng * np.sum((a - a.mean()) ** 2))
    print(f"\nref={ref} files={len(files)} blocks={nb} groups={ng}  E_trunc={E_trunc:.6f}")
    for k in res:
        print(f"  E({k}) = {res[k]:.6f} +/- {jerr(jk[k]):.6f}")
    for a, b in (("C", "0"), ("A", "C"), ("B", "C"), ("B", "A")):
        d = (res[a] - res[b]) * 1e3; e = jerr(jk[a] - jk[b]) * 1e3
        print(f"  E({a})-E({b}) = {d:+.3f} +/- {e:.3f} mHa")
    out = {f"v{k}": vecs[k] for k in vecs}
    out.update({f"E{k}": res[k] for k in res}); out.update({f"err{k}": jerr(jk[k]) for k in res})
    np.savez(os.path.join(DATA, f"cirelax_{ref}_vecs.npz"), **out)
    # report CI-coefficient change in (A): compare normalized d_I to the reference
    vA = vecs["A"]; N = len(cvec)
    dA = vA[:N]; cA = vA[N:]  # F-part = c * cvec
    lam = (cvec.conj() @ vA[N:]) / (cvec.conj() @ cvec)
    ph = dA @ cvec.conj() / abs(dA @ cvec.conj())
    dn = dA / ph; dn = dn / np.linalg.norm(dn)
    print("\n  (A) relaxed det coefficients vs frozen (normalized):")
    for i in range(min(N, 12)):
        print(f"    {i}: frozen {cvec[i].real:+.5f}  relaxed {dn[i].real:+.5f}{dn[i].imag:+.1e}j")
    cfr = (vecs["C"][N:] @ cvec.conj()) / (vecs["C"][:N] @ cvec.conj())
    cre = (cvec.conj() @ vA[N:]) / (cvec.conj() @ cvec) / (dA @ cvec.conj() / (cvec.conj() @ cvec))
    ov = abs(dn @ cvec.conj()) / np.linalg.norm(cvec)
    print(f"  geminal coef c (frozen C) = {cfr.real:+.5f}; (A) c/d_overlap = {cre.real:+.5f}; |<d_A|c_ref>| = {ov:.6f}")


if __name__ == "__main__":
    main()
