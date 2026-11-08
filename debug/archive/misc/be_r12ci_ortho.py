"""Be R12-CI, trustworthy energy via the ORTHOGONALIZED correlation basis.

The linear ansatz (1+cF)Phi_0 is ill-conditioned: f=1-e^{-r}->1, so F~6 (nearly parallel to
Phi_0) and the correlation energy is a difference of large near-equal numbers.  Fix: use
G = (F - <F>) Phi_0  (orthogonal to Phi_0) and compute every small quantity DIRECTLY as a
variance/covariance under |Phi_0|^2 -- no catastrophic cancellation.

2x2 in {Phi_0, G}, metric diag(1, sigma2):
   [[E0 , h ],   with  sigma2 = <G|G>   = E[(F-Fbar)^2]           (|Phi|^2)
    [ h , g ]]         h      = <Phi|H|G> = hT + E[(F-Fbar) V]     (V part |Phi|^2)
                       g      = <G|H|G>   = gT + E[(F-Fbar)^2 V]   (V part |Phi|^2; has the 4-body)
E0 = exact analytic reference (-14.539, noise-free).  Correlation shift = lowest eigenvalue - E0.
Kinetic pieces hT,gT: gradient form, smooth importance sampling (no node divergence):
   hT = 1/2 INT [ Phi0 gradPhi0.gradF + (F-Fbar) |gradPhi0|^2 ] / S00
   gT = 1/2 INT [ |gradF|^2 Phi0^2 + 2(F-Fbar) Phi0 gradF.gradPhi0 + (F-Fbar)^2 |gradPhi0|^2 ] / S00
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import be_r12ci_matelem as MM
from scipy.linalg import eigh

Phi0g = MM.phi0_and_grad
Fg = MM.F_and_grad
potential = MM.potential
E0 = MM.EREF


# ---- |Phi_0|^2 Metropolis: Fbar, sigma2, hV=E[(F-Fbar)V], gV=E[(F-Fbar)^2 V] ---- #
def phi2_stats(nw=4000, nsteps=12000, nequil=2000, step=0.4, seed=3):
    rng = np.random.default_rng(seed)
    R = rng.normal(scale=0.8, size=(nw, 4, 3))
    p = Phi0g(R)[0] ** 2
    Fs = []; Vs = []
    for it in range(nsteps):
        Rn = R + rng.normal(scale=step, size=R.shape)
        pn = Phi0g(Rn)[0] ** 2
        acc = rng.uniform(size=nw) < (pn / np.maximum(p, 1e-300))
        R[acc] = Rn[acc]; p[acc] = pn[acc]
        if it >= nequil and (it - nequil) % 5 == 0:
            Fs.append(Fg(R)[0]); Vs.append(potential(R))
    F = np.concatenate(Fs); V = np.concatenate(Vs)
    Fbar = F.mean()
    dF = F - Fbar
    sigma2 = np.mean(dF * dF)
    hV = np.mean(dF * V)
    gV = np.mean(dF * dF * V)
    Ev = V.mean()
    # simple block error bars
    nb = 20
    def berr(x):
        b = np.array_split(x, nb); return np.std([bi.mean() for bi in b], ddof=1) / np.sqrt(nb)
    return dict(Fbar=Fbar, sigma2=sigma2, hV=hV, gV=gV, Ev=Ev,
                e_sig=berr(dF * dF), e_hV=berr(dF * V), e_gV=berr(dF * dF * V))


# ---- smooth importance sampling: S00, and kinetic hT, gT (need Fbar) ---- #
def kinetic_parts(Fbar, nw=6_000_000, batch=1_000_000, seed=5):
    rng = np.random.default_rng(seed)
    S00 = 0.0; HT = 0.0; GT = 0.0; KE00 = 0.0; n = 0
    while n < nw:
        m = min(batch, nw - n)
        R, w = MM.sample(m, rng)
        phi, gph = Phi0g(R); F, gF = Fg(R)
        dF = F - Fbar
        gg_pp = np.sum(gph * gph, axis=(1, 2))          # |gradPhi0|^2
        gg_pF = np.sum(gph * gF, axis=(1, 2))           # gradPhi0.gradF
        gg_FF = np.sum(gF * gF, axis=(1, 2))            # |gradF|^2
        S00 += np.sum(w * phi * phi)
        KE00 += np.sum(w * 0.5 * gg_pp)
        HT += np.sum(w * 0.5 * (phi * gg_pF + dF * gg_pp))
        GT += np.sum(w * 0.5 * (gg_FF * phi * phi + 2 * dF * phi * gg_pF + dF * dF * gg_pp))
        n += m
    return S00 / nw, HT / nw, GT / nw, KE00 / nw


if __name__ == "__main__":
    print("=" * 66)
    print("Be R12-CI, orthogonalized basis  {Phi_0, (F-<F>)Phi_0}")
    print("=" * 66)
    st = phi2_stats()
    print(f"|Phi|^2:  Fbar={st['Fbar']:.4f}  sigma2={st['sigma2']:.5f}+/-{st['e_sig']:.5f}")
    print(f"          hV=E[(F-Fbar)V]={st['hV']:.5f}+/-{st['e_hV']:.5f}")
    print(f"          gV=E[(F-Fbar)^2 V]={st['gV']:.5f}+/-{st['e_gV']:.5f}   <- carries the 4-body Coulomb")
    S00, HT, GT, KE00 = kinetic_parts(st["Fbar"])
    hT = HT / S00; gT = GT / S00; ke00 = KE00 / S00
    print(f"smooth:   KE_00/S_00={ke00:.5f}  hT={hT:.5f}  gT={gT:.5f}")
    print(f"GATE  KE_00/S_00 + E[V] = {ke00 + st['Ev']:.5f}  (E0={E0:.5f}; diff {ke00+st['Ev']-E0:+.5f})")

    sigma2 = st["sigma2"]; h = hT + st["hV"]; g = gT + st["gV"]
    Sm = np.array([[1.0, 0.0], [0.0, sigma2]])
    Hm = np.array([[E0, h], [h, g]])
    w, v = eigh(Hm, Sm)
    print(f"\nh=<Phi|H|G>={h:.5f}   g=<G|H|G>={g:.5f}   sigma2={sigma2:.5f}")
    print(f"E_R12 = {w[0]:.5f} Ha   correlation captured {w[0]-E0:+.5f} Ha "
          f"({100*(w[0]-E0)/(-0.0944):.0f}% of the -94.4 mHa true corr)")
    print(f"variational (E>=exact -14.6674): {w[0] >= -14.6674}")
