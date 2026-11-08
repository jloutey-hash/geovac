r"""Gate G2 analysis: analytic transition columns (debug/data/ci_relaxed_<t>.npz) vs the real-space MC
(debug/data/ci_relaxed_g2_<t>_s*.npz).  Per determinant I:
    SF = <chi_I|F|psi>, KF = <chi_I|T F|psi> (kinetic), VF = <chi_I|V F|psi>, HF = KF + VF,
    s_I = SF - Fbar d_I,   t_I = HF - Fbar <chi_I|H|psi>.
Control variates (variance reduction only; their means are EXACTLY known and independent of the
code under test): c1 = Z<chi_I^* psi/P> - d_I, c2 = Z<chi_I^* H psi/P> - (H_engine d)_I (exact-integral
Slater-Condon engine, debug/data/ci_relaxed_engine_<t>.npz).  Coefficients by least squares over
blocks; errors from blocks re-binned into 60 groups (all seeds).
Usage: python debug/ci_relaxed_g2_analyze.py 1e-2
"""
import glob
import os
import sys

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def cv_estimate(x: np.ndarray, C: np.ndarray, ngroup: int = 60):
    """mean of x - C beta (C zero-mean control columns), error from re-binned groups."""
    Cc = C - C.mean(0)
    beta = np.linalg.lstsq(Cc, x - x.mean(), rcond=None)[0]
    y = x - C @ beta
    m = len(y) // ngroup
    g = y[:m * ngroup].reshape(ngroup, m).mean(1)
    return float(y.mean()), float(g.std(ddof=1) / np.sqrt(ngroup))


def main() -> int:
    tag = sys.argv[1] if len(sys.argv) > 1 else "1e-2"
    t = float(tag)
    an = np.load(os.path.join(DATA, f"ci_relaxed_{tag}.npz"), allow_pickle=True)
    He = np.load(os.path.join(DATA, f"ci_relaxed_engine_{tag}.npz"))['H_engine']
    files = sorted(glob.glob(os.path.join(DATA, f"ci_relaxed_g2_{t:.0e}_s*.npz")))
    B = {}
    for f in files:
        z = np.load(f)
        for k in ('S0', 'SF', 'KF', 'VF', 'Hpsi'):
            B.setdefault(k, []).append(z[k])
    B = {k: np.concatenate(v) for k, v in B.items()}
    nb = len(B['SF'])
    d = an['d_ref']; N = len(d)
    Fbar = float(an['SF'] @ d)
    Hd_e = He @ d
    ana = dict(SF=an['SF'], KF=an['KF'], VF=an['V1'] + an['FVee'])
    ana['HF'] = ana['KF'] + ana['VF']
    ana['s'] = an['s']; ana['t'] = an['tT'] + an['tVne'] + an['tVee']
    B['HF'] = B['KF'] + B['VF']
    B['s'] = B['SF'] - Fbar * B['S0']
    B['t'] = B['HF'] - Fbar * B['Hpsi']
    print(f"files={len(files)} blocks={nb}  N={N}  Fbar(analytic)={Fbar:.6f}")
    print(f"exact-known checks (no CV): max |MC S0 - d|/sigma and |MC Hpsi - H_engine d|/sigma:")
    for k, ex in (('S0', d), ('Hpsi', Hd_e)):
        zs = []
        for I in range(N):
            m, e = cv_estimate(B[k][:, I], np.zeros((nb, 1)))
            zs.append((ex[I] - m) / e)
        print(f"  {k}: z = " + " ".join(f"{z:+.2f}" for z in zs))
    worst = {}
    for k in ('SF', 'KF', 'VF', 'HF', 's', 't'):
        print(f"\n  {k}:  I   analytic        MC(CV)          sigma      z      rel.sigma")
        for I in range(N):
            C = np.stack([B['S0'][:, I] - d[I], B['Hpsi'][:, I] - Hd_e[I]], axis=1)
            m, e = cv_estimate(B[k][:, I], C)
            z = (ana[k][I] - m) / e
            rel = e / max(abs(ana[k][I]), 1e-300)
            print(f"     {I:3d} {ana[k][I]:+.7f}  {m:+.7f}  {e:.2e}  {z:+6.2f}   {rel:.1e}")
            worst.setdefault(k, []).append((abs(z), rel))
    print("\n  summary: max|z|, min/max relative sigma per column")
    for k, v in worst.items():
        zz = [a for a, _ in v]; rr = [b for _, b in v]
        print(f"   {k:3s} max|z|={max(zz):.2f}  rel.sigma {min(rr):.1e} .. {max(rr):.1e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
