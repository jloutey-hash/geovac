r"""Multi-geminal VMC LINEAR-METHOD diagnostic on the M=16 pi trunc=1e-2 reference.

Question: does a linear expansion in several geminals, Psi = (c0 + sum_a c_a F_a) Psi_CI,
close the gap between the single-geminal marriage (-8.042181, analytic) and the same-reference
2-body Pade exp-Jastrow VMC?

Mode 'lm': Metropolis on |g_guide Psi_CI|^2 (g_guide = the k=1 marriage multiplier), store per
sample (walker-averaged) the overlap/Hamiltonian matrices in the basis
    phi_0 = Psi_CI,  phi_a = A_a Psi_CI,
    A_a = F_a = sum_{i<j} f_a(r_ij) for 8 geminals {linexp,exp} x g in {0.25,0.5,1,2},
    plus A_9 = F_{linexp,0.5}^2 (the quadratic-term probe).
Estimators (non-symmetric, standard):  S_ab = < A_a A_b / g^2 >,
    H_ab = < A_a / g^2 * Re[ H(A_b Psi_CI) / Psi_CI ] >,
    H(A Psi)/Psi = A E_L0 - 1/2 sum_i lap_i A - sum_i grad_i A . grad_i ln Psi_CI  (standard form).
Mode 'pade': plain lih_vmc.vmc with PadeJastrow(gamma) on the same reference.

--nsweep is TOTAL sweeps including --nburn (as in lih_vmc).  Writes debug/data/mg_<tag>.npz.
"""
from __future__ import annotations

import os
import sys
import time
import argparse

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
DATA = os.path.join(HERE, "data")

# Namespace stub: import geovac submodules without running geovac/__init__.py (the working tree
# currently has an in-flight archive move that breaks the package __init__; we touch nothing there).
import types as _types
if "geovac" not in sys.modules:
    _pkg = _types.ModuleType("geovac")
    _pkg.__path__ = [os.path.join(os.path.dirname(HERE), "geovac")]
    sys.modules["geovac"] = _pkg

import lih_vmc
from lih_pi_m16_vmc import _build_wf_from_no

REF_NPZ = os.path.join(DATA, "lih_pi_m16_energy_trunc1e-2.npz")

GEMINALS = [(kind, g) for kind in ("linexp", "exp") for g in (0.25, 0.5, 1.0, 2.0)]
LABELS = ["Psi0"] + [f"{k}{g}" for k, g in GEMINALS] + ["F2(linexp0.5)"]
IDX_LX05 = 1 + GEMINALS.index(("linexp", 0.5))


def _f(kind: str, g: float, d: np.ndarray):
    """f, f', lap f = f'' + 2 f'/d."""
    e = np.exp(-g * d)
    if kind == "exp":
        f = e; fp = -g * e; fpp = g * g * e
    else:
        f = d * e; fp = (1.0 - g * d) * e; fpp = (-2.0 * g + g * g * d) * e
    return f, fp, fpp + 2.0 * fp / d


def basis_terms(r: np.ndarray):
    """A_a, grad_i A_a, sum_i lap_i A_a for the 8 geminals + F^2.  Returns
    A (nb,Nw), G (nb,Nw,Ne,3), L (nb,Nw) for nb = 9 non-trivial functions."""
    Nw, Ne, _ = r.shape
    ng = len(GEMINALS)
    A = np.zeros((ng + 1, Nw)); G = np.zeros((ng + 1, Nw, Ne, 3)); L = np.zeros((ng + 1, Nw))
    for i in range(Ne):
        for j in range(i + 1, Ne):
            dvec = r[:, i] - r[:, j]
            d = np.linalg.norm(dvec, axis=-1)
            d = np.where(d < 1e-12, 1e-12, d)
            u = dvec / d[:, None]
            for a, (kind, g) in enumerate(GEMINALS):
                f, fp, fl = _f(kind, g, d)
                A[a] += f
                G[a, :, i] += fp[:, None] * u
                G[a, :, j] -= fp[:, None] * u
                L[a] += 2.0 * fl
    k = IDX_LX05 - 1
    F = A[k]; gF = G[k]; lF = L[k]
    A[ng] = F * F
    G[ng] = 2.0 * F[:, None, None] * gF
    L[ng] = 2.0 * np.sum(gF * gF, axis=(1, 2)) + 2.0 * F * lF
    return A, G, L


def sample_matrices(wf, r, cg: float, Fbar_g: float):
    """Walker-averaged S (nb,nb), H (nb,nb) for one configuration set r."""
    psi, gpsi, lpsi = lih_vmc.psi_ci_full(wf, r)
    glnpsi = gpsi / psi[:, None, None]
    EL0 = lih_vmc.potential(wf, r) - 0.5 * np.sum(lpsi, axis=1) / psi      # complex
    A, G, L = basis_terms(r)
    Nw = r.shape[0]
    Aall = np.vstack([np.ones((1, Nw)), A])                                # (nb,Nw)
    HA = np.empty_like(Aall)
    HA[0] = EL0.real
    cross = np.einsum('aWed,Wed->aW', G, glnpsi)                           # complex
    HA[1:] = (A * EL0[None, :] - 0.5 * L - cross).real
    gg = 1.0 + cg * (A[IDX_LX05 - 1] - Fbar_g)
    w = 1.0 / gg ** 2
    S = np.einsum('aW,bW,W->ab', Aall, Aall, w) / Nw
    H = np.einsum('aW,bW,W->ab', Aall, HA, w) / Nw
    return S, H


def run_lm(wf, cg, Fbar_g, nwalk, nsweep, nburn, seed, thin, step=0.45):
    rng = np.random.default_rng(seed)
    Ne = wf.na + wf.nb
    r = lih_vmc._init_walkers(wf, nwalk, rng)
    kind_g = ("linexp", 0.5)

    def logp(rr):
        psi = lih_vmc.psi_ci_val(wf, rr)
        Nw = rr.shape[0]; F = np.zeros(Nw)
        for i in range(Ne):
            for j in range(i + 1, Ne):
                d = np.linalg.norm(rr[:, i] - rr[:, j], axis=-1)
                F += _f(kind_g[0], kind_g[1], np.maximum(d, 1e-12))[0]
        g = 1.0 + cg * (F - Fbar_g)
        return 2.0 * np.log(np.abs(psi) + 1e-300) + 2.0 * np.log(np.abs(g) + 1e-300)

    lp = logp(r)
    acc = ntry = 0; sd = step
    Ss, Hs = [], []
    t0 = time.time()
    for sweep in range(nsweep):
        for i in range(Ne):
            rprop = r.copy()
            rprop[:, i] += rng.normal(0, sd, size=(nwalk, 3))
            lpp = logp(rprop)
            a = np.log(rng.uniform(size=nwalk)) < (lpp - lp)
            r[a] = rprop[a]; lp[a] = lpp[a]
            acc += int(a.sum()); ntry += nwalk
        if sweep < nburn and sweep % 50 == 49:
            ar = acc / max(ntry, 1)
            if ar > 0.6: sd *= 1.1
            elif ar < 0.4: sd /= 1.1
            acc = ntry = 0
        if sweep >= nburn and (sweep - nburn) % thin == 0:
            S, H = sample_matrices(wf, r, cg, Fbar_g)
            Ss.append(S); Hs.append(H)
        if sweep % 500 == 0:
            print(f"   sweep {sweep}/{nsweep}  step={sd:.3f}  t={time.time()-t0:.0f}s", flush=True)
    return np.array(Ss), np.array(Hs), acc / max(ntry, 1), sd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='lm', choices=['lm', 'pade'])
    ap.add_argument('--gamma', type=float, default=1.5)
    ap.add_argument('--nwalk', type=int, default=1000)
    ap.add_argument('--nsweep', type=int, default=9000, help='TOTAL sweeps incl. burn-in')
    ap.add_argument('--nburn', type=int, default=4000)
    ap.add_argument('--thin', type=int, default=5)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tag', default=None)
    args = ap.parse_args()

    z = np.load(REF_NPZ, allow_pickle=True)
    E_trunc = float(z['E_trunc']); Mk = int(z['Mk'])
    cg = float(z['linexp_c']); Fbar_g = float(z['linexp_Fbar'])
    wf = _build_wf_from_no(z['ct'], z['T_no'], Mk, E_trunc)
    t0 = time.time()
    if args.mode == 'lm':
        tag = args.tag or f"lm_s{args.seed}"
        print(f"[lm] seed={args.seed} nwalk={args.nwalk} nsweep={args.nsweep} nburn={args.nburn} "
              f"guide c={cg:+.5f} Fbar={Fbar_g:.5f}", flush=True)
        Ss, Hs, acc, sd = run_lm(wf, cg, Fbar_g, args.nwalk, args.nsweep, args.nburn,
                                 args.seed, args.thin)
        np.savez(os.path.join(DATA, f"mg_{tag}.npz"), S=Ss, H=Hs, labels=np.array(LABELS),
                 seed=args.seed, nwalk=args.nwalk, nsweep=args.nsweep, nburn=args.nburn,
                 acc=acc, step=sd, E_trunc=E_trunc, c_guide=cg, Fbar_guide=Fbar_g)
        S = Ss.mean(0); H = Hs.mean(0)
        import scipy.linalg as sl
        sub = [0, IDX_LX05]
        ev = sl.eigvals(H[np.ix_(sub, sub)], S[np.ix_(sub, sub)]).real.min()
        print(f"[lm] done {time.time()-t0:.0f}s acc={acc:.2f}  E(Psi0)={H[0,0]/S[0,0]:.6f}  "
              f"E_k1(linexp0.5)={ev:.6f}", flush=True)
    else:
        tag = args.tag or f"pade_g{args.gamma}_s{args.seed}"
        jas = lih_vmc.PadeJastrow(args.gamma, wf.na, wf.nb)
        E, err, acc = lih_vmc.vmc(wf, nwalk=args.nwalk, nsweep=args.nsweep, nburn=args.nburn,
                                  jas=jas, seed=args.seed, thin=args.thin)
        np.savez(os.path.join(DATA, f"mg_{tag}.npz"), E=E, err=err, gamma=args.gamma,
                 seed=args.seed, nwalk=args.nwalk, nsweep=args.nsweep, nburn=args.nburn)
        print(f"[pade] gamma={args.gamma} seed={args.seed} E={E:.6f} +/- {err:.6f} "
              f"({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
