"""Be R12-CI, trustworthy energy.  2x2 {Phi_0, F Phi_0}.

Variance control (the fix for the earlier noisy energy):
  * SINGULAR Coulomb parts  E[V F], E[F^2 V]  -- sampled from |Phi_0|^2 (Metropolis).
    Under |Phi_0|^2 the same-spin Fermi hole tames 1/r12 and opposite-spin densities are
    smooth, so the Coulomb estimators have low variance (this is why VMC samples |Psi|^2).
  * SMOOTH parts  S_01=E[F], S_11=E[F^2], and the KINETIC 1/2 INT grad.grad  -- importance
    sampled from a smooth product density (gradient form, no 1/Psi node divergence).
Both use the SAME normalization S_00 (computed in each sampler; cross-checked).

Matrix elements (V = V_ne + V_ee multiplicative; normalized so S_00 = 1):
   S = [[1, sF],[sF, sF2]] ,  sF=E|Phi|2[F]? -- NO: overlaps are SMOOTH, use smooth MC.
   H_00 = KE_00 + E|Phi|2[V]           (gate: must equal E0 = -14.539)
   H_01 = KE_01 + E|Phi|2[V F]
   H_11 = KE_11 + E|Phi|2[F^2 V]
Solve H c = E S c with NOISE-CONTROLLED (not eigensolver-biased) elements.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import be_r12ci_matelem as MM   # orbitals, phi0_and_grad, F_and_grad, potential, sample

Phi0g = MM.phi0_and_grad
Fg = MM.F_and_grad
potential = MM.potential
EREF = MM.EREF


def Phi0(R):
    return Phi0g(R)[0]


def Fsum(R):
    return Fg(R)[0]


# --------- smooth importance sampling: S_00,S_01,S_11 + kinetic KE_ij -------- #
def smooth_part(nw=4_000_000, batch=1_000_000, seed=1):
    rng = np.random.default_rng(seed)
    s = np.zeros(3); ke = np.zeros(3); n = 0
    while n < nw:
        m = min(batch, nw - n)
        R, w = MM.sample(m, rng)
        phi, gph = Phi0g(R)
        F, gF = Fg(R)
        p0 = phi; p1 = F * phi
        gp0 = gph; gp1 = gF * phi[:, None, None] + F[:, None, None] * gph
        s[0] += np.sum(w * p0 * p0); s[1] += np.sum(w * p0 * p1); s[2] += np.sum(w * p1 * p1)
        ke[0] += np.sum(w * 0.5 * np.sum(gp0 * gp0, axis=(1, 2)))
        ke[1] += np.sum(w * 0.5 * np.sum(gp0 * gp1, axis=(1, 2)))
        ke[2] += np.sum(w * 0.5 * np.sum(gp1 * gp1, axis=(1, 2)))
        n += m
    return s / nw, ke / nw          # S_00,S_01,S_11 (unnormalized) ; KE_00,KE_01,KE_11


# --------- |Phi_0|^2 Metropolis: the singular Coulomb expectations ----------- #
def phi2_part(nw=3000, nsteps=4000, nequil=1000, step=0.4, seed=2):
    rng = np.random.default_rng(seed)
    R = rng.normal(scale=0.8, size=(nw, 4, 3))
    p = Phi0(R) ** 2
    accV = []  # E[V], E[VF], E[F^2 V]
    for it in range(nsteps):
        Rn = R + rng.normal(scale=step, size=R.shape)
        pn = Phi0(Rn) ** 2
        acc = rng.uniform(size=nw) < (pn / np.maximum(p, 1e-300))
        R[acc] = Rn[acc]; p[acc] = pn[acc]
        if it >= nequil and (it - nequil) % 5 == 0:
            V = potential(R); F = Fsum(R)
            accV.append(np.stack([V, V * F, F * F * V], axis=1))
    A = np.concatenate(accV, axis=0)          # (Nsamp, 3)
    return A.mean(axis=0), A.std(axis=0) / np.sqrt(A.shape[0])


if __name__ == "__main__":
    print("=" * 68)
    print("Be R12-CI trustworthy energy  (|Phi|^2 Coulomb + smooth kinetic)")
    print("=" * 68)
    (S00, S01u, S11u), (KE00, KE01, KE11) = smooth_part()
    # normalize by S00
    S01 = S01u / S00; S11 = S11u / S00
    ke00, ke01, ke11 = KE00 / S00, KE01 / S00, KE11 / S00
    print(f"smooth: S_00={S00:.1f}  S_01/S_00={S01:.5f}  S_11/S_00={S11:.5f}")
    print(f"        KE_00/S_00={ke00:.5f}  KE_01/S_00={ke01:.5f}  KE_11/S_00={ke11:.5f}")

    (EV, EVF, EF2V), (eV, eVF, eF2V) = phi2_part()
    print(f"|Phi|^2: E[V]={EV:.5f}+/-{eV:.5f}  E[VF]={EVF:.5f}+/-{eVF:.5f}  "
          f"E[F^2 V]={EF2V:.5f}+/-{eF2V:.5f}")

    gate = ke00 + EV
    print(f"\nGATE  H_00/S_00 = KE_00/S_00 + E[V] = {gate:.5f}   (E0={EREF:.5f}; diff {gate-EREF:+.5f})")

    S = np.array([[1.0, S01], [S01, S11]])
    H = np.array([[ke00 + EV, ke01 + EVF], [ke01 + EVF, ke11 + EF2V]])
    from scipy.linalg import eigh
    w, v = eigh(H, S)
    print(f"\nS = {S.ravel()}")
    print(f"H = {H.ravel()}")
    print(f"E_R12 = {w[0]:.5f} Ha   correlation captured {w[0]-EREF:+.5f} Ha")
    print(f"Be exact (nonrel) -14.6674 ; variational: {w[0] >= -14.6674}")
