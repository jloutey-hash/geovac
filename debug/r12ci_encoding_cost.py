"""Encoding-cost probe for the R12-CI -> quantum-algorithm pivot (He).

QUESTION.  R12-CI (variational explicitly-correlated CI) reaches He 0.45-0.80 mHa
at a qubit-proxy of ~6-8 (debug/ctf12_r12ci_he.py, sprint_ctf12_poc_memo.md).  BUT
its geminal G = e^{-g r12} R_ref(r1) R_ref(r2) is a NON-ORTHOGONAL correlated basis
function, so encoding the R12-CI eigenproblem on a quantum computer is a
non-orthogonal quantum eigensolver (NOQE) / generalized-eigenvalue problem whose
cost is governed by the overlap-matrix conditioning kappa(S).

This driver reuses the existing R12-CI integral machinery READ-ONLY (imported from
debug/ctf12_r12ci_he.py) and measures kappa2(S) at the *accurate operating point*,
plus its growth as geminals/orbitals are added and as the geminal exponent g is
varied.  It then maps kappa(S) to the NOQE resource implication:

  * fault-tolerant GEVP route (Liang et al. 2112.02554): runtime ~ kappa_B directly,
    with kappa_B the condition number of the metric (overlap) block.
  * near-term NOQE route (Baek et al. 2205.09039; improved 2608.12830): eigenvalue
    error / measurement cost controlled by kappa of the RETAINED (thresholded)
    overlap block.

No production/paper/CLAUDE.md edits; debug/ only.  The imported module is not
modified (used purely as an integral library).
"""
import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

# ---- import the existing R12-CI engine as a READ-ONLY integral library ----
_spec = importlib.util.spec_from_file_location(
    "ctf12_r12ci_he", os.path.join(HERE, "ctf12_r12ci_he.py"))
ctf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ctf)

EXACT_HE = ctf.EXACT_HE

# Fixed diffuse grid matching run_sweep() in ctf12_r12ci_he.py so energies match
# the recorded best-(k,gamma) operating points in debug/data/ctf12_r12ci_he.json.
NG, NX = 600, 128
KMIN = 1.30
R_MAX = 42.0 / KMIN
_RGRID = ctf.make_grid(KMIN, NG, r_max=R_MAX)


def build_bfs(ns, n_gem, Rtab, dRtab):
    bfs = []
    for i in range(1, ns + 1):
        for j in range(i, ns + 1):
            bfs.append(ctf.orbital_pair(i, j, Rtab, dRtab))
    for ref in range(1, n_gem + 1):
        bfs.append(ctf.geminal(ref, Rtab, dRtab))
    return bfs


def overlap_and_H(ns, n_gem, k, gamma):
    """Assemble H, M for the (ns orbital-pairs + n_gem geminals) basis at (k,gamma)."""
    r, wr = _RGRID
    Rtab, dRtab = ctf.build_tabs(max(ns, n_gem, 1), r, float(k))
    K = ctf.build_kernels(r, gamma, nx=NX)
    bfs = build_bfs(ns, n_gem, Rtab, dRtab)
    H, M = ctf.assemble(bfs, r, wr, K, Z=2)
    return H, M, len(bfs)


def cond_report(M):
    """Condition numbers of the overlap matrix.

    Returns raw kappa (unnormalized basis) and the physically-meaningful
    normalized kappa (unit-diagonal correlation matrix -- the quantum states
    prepared in NOQE are normalized, so S_ii = 1)."""
    ev_raw = np.linalg.eigvalsh(M)
    D = np.sqrt(np.diag(M))
    Sn = M / np.outer(D, D)
    Sn = 0.5 * (Sn + Sn.T)
    ev = np.linalg.eigvalsh(Sn)
    ev = np.clip(ev, 1e-300, None)
    return {
        "kappa_raw": float(ev_raw.max() / ev_raw.min()),
        "kappa_norm": float(ev.max() / ev.min()),
        "lam_min_norm": float(ev.min()),
        "lam_max_norm": float(ev.max()),
        "lam_min_raw": float(ev_raw.min()),
        "eigs_norm": [float(x) for x in ev],
    }, Sn


def retained_kappa(Sn, thresholds):
    """kappa of the retained (thresholded) overlap block = the NOQE cost driver
    (Baek et al. 2608.12830: eigenvalue error controlled by kappa of the
    canonical-orthogonalization-retained block after discarding eigenvalues < tau)."""
    ev = np.linalg.eigvalsh(Sn)
    out = {}
    lam_max = ev.max()
    for tau in thresholds:
        kept = ev[ev > tau]
        if kept.size == 0:
            out[f"{tau:.0e}"] = {"n_kept": 0, "kappa_retained": None}
        else:
            out[f"{tau:.0e}"] = {"n_kept": int(kept.size),
                                 "kappa_retained": float(lam_max / kept.min())}
    return out


def energy(H, M, thr=1e-8):
    return ctf.solve_gen(H, M, thr=thr)


# --------------------------------------------------------------------------
# Accurate operating points recovered from the R12-CI sweep (best (k,gamma)).
OPPOINTS = [
    dict(ns=3, n_gem=1, k=1.7, gamma=0.7, ref_err_mHa=0.798),
    dict(ns=3, n_gem=2, k=1.7, gamma=0.7, ref_err_mHa=0.640),
    dict(ns=4, n_gem=1, k=1.6, gamma=0.7, ref_err_mHa=0.474),
    dict(ns=4, n_gem=2, k=1.6, gamma=0.7, ref_err_mHa=0.449),
]

THRESH = [1e-4, 1e-6, 1e-8, 1e-10]


def run():
    out = {"system": "He", "exact": EXACT_HE,
           "grid": {"Ng": NG, "nx": NX, "r_max": R_MAX},
           "note": "kappa_norm = cond of unit-diagonal overlap (NOQE-relevant); "
                   "geminal G=e^{-g r12}R_ref R_ref is the non-orthogonal correlated function.",
           "operating_points": [], "orbital_only_trend": [],
           "geminal_trend_fixed_ns3": [], "gamma_scan_ns3_ngem1": [],
           "collinearity_probe": None}

    # ---- (1) the accurate operating points ----
    print("=== Accurate operating points (R12-CI, best k,gamma) ===")
    print(f"{'ns':>3}{'ngem':>5}{'nbf':>5}{'E':>13}{'err_mHa':>9}"
          f"{'kappa_norm':>12}{'kappa_raw':>12}{'lam_min':>11}")
    for op in OPPOINTS:
        H, M, nbf = overlap_and_H(op["ns"], op["n_gem"], op["k"], op["gamma"])
        E, kept = energy(H, M)
        cr, Sn = cond_report(M)
        rk = retained_kappa(Sn, THRESH)
        rec = dict(op)
        rec.update(nbf=nbf, qubit_proxy=2 * op["ns"], E=float(E),
                   err_mHa=float((E - EXACT_HE) * 1000), kept=int(kept),
                   retained_kappa=rk, **{kk: cr[kk] for kk in
                   ("kappa_raw", "kappa_norm", "lam_min_norm", "lam_max_norm", "lam_min_raw")})
        rec["eigs_norm"] = cr["eigs_norm"]
        out["operating_points"].append(rec)
        print(f"{op['ns']:>3}{op['n_gem']:>5}{nbf:>5}{E:>13.6f}"
              f"{(E-EXACT_HE)*1000:>9.3f}{cr['kappa_norm']:>12.1f}"
              f"{cr['kappa_raw']:>12.1f}{cr['lam_min_norm']:>11.2e}")

    # ---- (2) orbital-only trend (n_gem=0): conditioning of the Sturmian pair block ----
    print("\n=== Orbital-only overlap conditioning (no geminal) at k=1.7 ===")
    print(f"{'ns':>3}{'nbf':>5}{'kappa_norm':>12}{'lam_min':>11}")
    for ns in [1, 2, 3, 4, 5, 6]:
        H, M, nbf = overlap_and_H(ns, 0, 1.7, 1.0)  # gamma irrelevant for c=1 block
        cr, _ = cond_report(M)
        out["orbital_only_trend"].append(
            dict(ns=ns, nbf=nbf, kappa_norm=cr["kappa_norm"],
                 lam_min_norm=cr["lam_min_norm"], kappa_raw=cr["kappa_raw"]))
        print(f"{ns:>3}{nbf:>5}{cr['kappa_norm']:>12.1f}{cr['lam_min_norm']:>11.2e}")

    # ---- (3) geminal trend at fixed ns=3 (add geminals) ----
    print("\n=== Adding geminals at fixed ns=3, k=1.7, gamma=0.7 ===")
    print(f"{'ngem':>5}{'nbf':>5}{'E':>13}{'err_mHa':>9}{'kappa_norm':>12}{'lam_min':>11}")
    for n_gem in [0, 1, 2, 3]:
        H, M, nbf = overlap_and_H(3, n_gem, 1.7, 0.7)
        E, kept = energy(H, M)
        cr, _ = cond_report(M)
        out["geminal_trend_fixed_ns3"].append(
            dict(n_gem=n_gem, nbf=nbf, E=float(E),
                 err_mHa=float((E - EXACT_HE) * 1000), kappa_norm=cr["kappa_norm"],
                 lam_min_norm=cr["lam_min_norm"], kappa_raw=cr["kappa_raw"]))
        print(f"{n_gem:>5}{nbf:>5}{E:>13.6f}{(E-EXACT_HE)*1000:>9.3f}"
              f"{cr['kappa_norm']:>12.1f}{cr['lam_min_norm']:>11.2e}")

    # ---- (4) gamma scan at ns=3, n_gem=1: the accuracy<->conditioning tension ----
    print("\n=== gamma scan at ns=3, n_gem=1, k=1.7 (accuracy vs conditioning) ===")
    print(f"{'gamma':>7}{'E':>13}{'err_mHa':>9}{'kappa_norm':>12}{'lam_min':>11}")
    for g in [0.4, 0.5, 0.7, 0.9, 1.1, 1.3, 1.6, 2.0, 2.5, 3.0]:
        H, M, nbf = overlap_and_H(3, 1, 1.7, g)
        E, kept = energy(H, M)
        cr, _ = cond_report(M)
        out["gamma_scan_ns3_ngem1"].append(
            dict(gamma=g, E=float(E), err_mHa=float((E - EXACT_HE) * 1000),
                 kappa_norm=cr["kappa_norm"], lam_min_norm=cr["lam_min_norm"]))
        print(f"{g:>7.2f}{E:>13.6f}{(E-EXACT_HE)*1000:>9.3f}"
              f"{cr['kappa_norm']:>12.1f}{cr['lam_min_norm']:>11.2e}")

    # ---- (5) which pair is (near-)collinear?  <G|(1,1) orbital pair> normalized ----
    # Overlap of the geminal with each orbital pair at the operating point.
    r, wr = _RGRID
    Rtab, dRtab = ctf.build_tabs(3, r, 1.7)
    K = ctf.build_kernels(r, 0.7, nx=NX)
    bfs = build_bfs(3, 1, Rtab, dRtab)
    _, M = ctf.assemble(bfs, r, wr, K, Z=2)
    D = np.sqrt(np.diag(M))
    Sn = M / np.outer(D, D)
    labels = [f"({i},{j})" for i in range(1, 4) for j in range(i, 4)] + ["G(ref1)"]
    gidx = len(bfs) - 1
    coll = {labels[a]: float(Sn[gidx, a]) for a in range(len(bfs)) if a != gidx}
    out["collinearity_probe"] = {
        "description": "normalized overlap <G|orbital-pair> at ns=3,n_gem=1,k=1.7,gamma=0.7",
        "geminal_vs_orbital_pairs": coll}
    print("\n=== Geminal collinearity: normalized <G|orbital-pair> ===")
    for lab, v in coll.items():
        print(f"  <G|{lab:>7}> = {v:+.4f}")

    return out


if __name__ == "__main__":
    t0 = time.time()
    res = run()
    res["walltime_s"] = time.time() - t0
    os.makedirs("debug/data", exist_ok=True)
    with open("debug/data/r12ci_encoding_cost.json", "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nwall {time.time()-t0:.1f}s -> debug/data/r12ci_encoding_cost.json")
