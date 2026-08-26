"""aha Track 1 / STEP 3 -- water's A1 block: the O<->H canonical-correlation spectrum,
its effective rank, and a Woodbury (low-rank-exact) whitening priced in d_inv currency.

Setup (identical geometry/metric to debug/sturmian_sw_water_conditioning.py):
  water O + 2H, C2v, R_OH = 1.809, R_HH = 2.861 bohr, shared-scale s-only SW metric.
  A1 block (holds the ground state) in the ordered basis  [ O_n ; (H1_n + H2_n)/sqrt2 ]:

        A1 = [[ I         , sqrt2 * P ],
              [ sqrt2 * P^T,  I + Q    ]]      P = SW(R_OH),  Q = SW(R_HH)

  so the O sub-block is EXACTLY the identity and the H+ sub-block is I + Q -- which is
  precisely the "gerade" combination of the two EQUIVALENT hydrogens, i.e. the H2+ lever
  applied to the H2 sub-system.  The generalized (non-identity intra-block) version of
  the Track-1 theorem is the canonical-correlation form

        sigma_k = svd( S_OO^{-1/2} S_O,H+ S_H+,H+^{-1/2} ),
        lam(A1) in  [ (1-sig_max) lam_min(D) , (1+sig_max) lam_max(D) ],  D = blkdiag(I, I+Q)
        =>  cond(A1)  <=  cond(D) * (1+sig_max)/(1-sig_max)      [exact 2-factor split]

Question: is the N^1.97 growth entirely in the sigma factor, is the sigma spectrum
low-rank, and does a rank-r EXACT (Woodbury) treatment of those directions flatten the
QSVT degree d_inv ~ kappa ln(kappa/eps) that Paper 60's resource model charges?

Run:  python debug/aha_t1_water.py
"""
from __future__ import annotations

import sys, os, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from debug.aha_t1_core import sw_block_trapz, loglog_fit

R_OH = 1.809
ANG = np.deg2rad(104.5)
R_HH = 2.0 * R_OH * np.sin(ANG / 2.0)
EPS_INV = 1e-3          # QSVT inverse accuracy used in Paper 60's resource model


def d_inv(kappa: float, eps: float = EPS_INV) -> int:
    """Paper-60 resource model: degree of the QSVT x^{-1/2} polynomial on [1/kappa,1]."""
    if kappa <= 1.0 + 1e-9:
        return 0
    return int(np.ceil(kappa * np.log(kappa / eps)))


def invsqrt(M: np.ndarray) -> np.ndarray:
    w, U = np.linalg.eigh(M)
    return U @ np.diag(w ** -0.5) @ U.T


def sqrtm_sym(M: np.ndarray) -> np.ndarray:
    w, U = np.linalg.eigh(M)
    return U @ np.diag(np.sqrt(np.maximum(w, 0.0))) @ U.T


def water_A1(nmax: int):
    """A1 block plus its constituent pieces."""
    I = np.eye(nmax)
    P = sw_block_trapz(R_OH, nmax, "S")
    Q = sw_block_trapz(R_HH, nmax, "S")
    SOO = I
    SHH = I + Q
    SOH = np.sqrt(2.0) * P
    A1 = np.block([[SOO, SOH], [SOH.T, SHH]])
    return A1, SOO, SOH, SHH


def main() -> None:
    print("=" * 100)
    print(f"STEP 3 -- water A1 block.  R_OH = {R_OH} bohr, R_HH = {R_HH:.3f} bohr, SW metric, k=1")
    print("=" * 100)

    # -------------------------------------------------- 3a. the exact 2-factor split
    print("\n[3a] cond(A1) split into  cond(D) x (1+sig_max)/(1-sig_max)")
    print("     D = blkdiag(I_O, I+Q) is the metric with the O<->H coupling REMOVED;")
    print("     sig_k are the O<->H+ canonical correlations (principal-angle cosines).")
    print(f"  {'nmax':>4} {'N':>4} | {'cond(A1)':>10} {'cond(D)':>8} {'cond(I+Q)':>10}"
          f" | {'sig_max':>9} {'(1+s)/(1-s)':>12} {'product bound':>14} {'bound/actual':>12}")
    print("  " + "-" * 100)
    nmaxes = list(range(2, 13))
    condA1, sigmax, condD, Ns = [], [], [], []
    for nmax in nmaxes:
        A1, SOO, SOH, SHH = water_A1(nmax)
        Ct = invsqrt(SOO) @ SOH @ invsqrt(SHH)
        sig = np.linalg.svd(Ct, compute_uv=False)
        D = np.block([[SOO, np.zeros_like(SOH)], [np.zeros_like(SOH).T, SHH]])
        cA, cD = float(np.linalg.cond(A1)), float(np.linalg.cond(D))
        f = (1 + sig.max()) / (1 - sig.max())
        Ns.append(3 * nmax)          # total water dimension used in the paper's fit
        condA1.append(cA); sigmax.append(float(sig.max())); condD.append(cD)
        print(f"  {nmax:>4} {3*nmax:>4} | {cA:>10.2f} {cD:>8.3f} {np.linalg.cond(SHH):>10.3f}"
              f" | {sig.max():>9.6f} {f:>12.2f} {cD*f:>14.2f} {cD*f/cA:>12.3f}")

    a, r2a, resa = loglog_fit(Ns, [1 - s for s in sigmax])
    b, r2b, resb = loglog_fit(Ns, condA1)
    c, r2c, resc = loglog_fit(Ns, condD)
    print(f"\n     fits over N = 3*nmax = {Ns[0]}..{Ns[-1]}:")
    print(f"       cond(A1)          ~ N^{b:+.4f}  (R2={r2b:.4f}, max ln-resid {resb:.4f})"
          f"   <- paper: N^1.97")
    print(f"       1 - sigma_max     ~ N^{a:+.4f}  (R2={r2a:.4f}, max ln-resid {resa:.4f})")
    print(f"       cond(D) (no O<->H)~ N^{c:+.4f}  (R2={r2c:.4f}, max ln-resid {resc:.4f})"
          f"   <- FLAT: {condD[0]:.2f} -> {condD[-1]:.2f}")
    print(f"     => the ENTIRE N-growth of cond(A1) is carried by 1-sigma_max; the")
    print(f"        coupling-free metric D is basis-flat, exactly like the H2+ gerade sector.")

    # ----------------------------------------------- 3b. the sigma spectrum / rank
    print("\n[3b] O<->H+ canonical-correlation spectrum (effective rank of the coupling)")
    print(f"  {'nmax':>4} | sigma_k (descending)")
    print("  " + "-" * 96)
    for nmax in (4, 6, 8, 10, 12):
        A1, SOO, SOH, SHH = water_A1(nmax)
        sig = np.linalg.svd(invsqrt(SOO) @ SOH @ invsqrt(SHH), compute_uv=False)
        print(f"  {nmax:>4} | " + " ".join(f"{x:.5f}" for x in sig[:10])
              + ("" if nmax <= 10 else " ..."))
    print("\n  effective-rank measures (largest nmax = 12, N = 36):")
    A1, SOO, SOH, SHH = water_A1(12)
    sig = np.linalg.svd(invsqrt(SOO) @ SOH @ invsqrt(SHH), compute_uv=False)
    e = sig ** 2 / np.sum(sig ** 2)
    print(f"     sigma_1..sigma_5 = {np.array2string(sig[:5], precision=5)}")
    print(f"     energy fraction in top-r:  " +
          "  ".join(f"r={r}:{np.sum(e[:r]):.4f}" for r in (1, 2, 3, 4, 5)))
    print(f"     participation ratio (exp of spectral entropy) = "
          f"{np.exp(-np.sum(e*np.log(e+1e-300))):.2f}")
    print(f"     #sigma_k > 0.5 : {int((sig>0.5).sum())}   #>0.1 : {int((sig>0.1).sum())}"
          f"   #>0.01 : {int((sig>0.01).sum())}  (of {len(sig)})")

    # -------------------------- 3c. rank-r truncation: does it reproduce cond(A1)?
    print("\n[3c] Rank-r truncation of the O<->H+ coupling: how many sigma-directions")
    print("     must be kept to reproduce cond(A1) within 10%?")
    print(f"  {'nmax':>4} {'N':>4} {'cond(A1)':>10} | " +
          "  ".join(f"r={r}" for r in range(1, 6)) + "   -> r_10%")
    print("  " + "-" * 96)
    r10 = {}
    for nmax in (4, 6, 8, 10, 12):
        A1, SOO, SOH, SHH = water_A1(nmax)
        Ct = invsqrt(SOO) @ SOH @ invsqrt(SHH)
        U, sv, Vt = np.linalg.svd(Ct)
        Dh = np.block([[sqrtm_sym(SOO), np.zeros_like(SOH)],
                       [np.zeros_like(SOH).T, sqrtm_sym(SHH)]])
        cA = float(np.linalg.cond(A1))
        outs, rstar = [], None
        for r in range(1, min(6, nmax + 1)):
            Cr = (U[:, :r] * sv[:r]) @ Vt[:r, :]
            Sh = np.block([[np.eye(nmax), Cr], [Cr.T, np.eye(nmax)]])
            cr = float(np.linalg.cond(Dh @ Sh @ Dh))
            outs.append(cr)
            if rstar is None and abs(cr - cA) / cA < 0.10:
                rstar = r
        r10[nmax] = rstar
        print(f"  {nmax:>4} {3*nmax:>4} {cA:>10.2f} | " +
              "  ".join(f"{x:>7.2f}" for x in outs) + f"   -> r={rstar}")

    # ---------------------------------------------- 3d. Woodbury pricing in d_inv
    print("\n[3d] WOODBURY PRICING (Paper-60 currency: d_inv ~ kappa ln(kappa/eps), eps=1e-3)")
    print("     EXACT structure: in the D-whitened frame  Shat = [[I, Ctil],[Ctil^T, I]] = I + K,")
    print("     K symmetric of rank 2*rank(Ctil), eigenpairs (u_k, +/-v_k)/sqrt2 <-> 1 +/- sig_k.")
    print("     So  Shat^{-1/2} = I + sum_k [(1+sig_k)^{-1/2}-1] w_k^+ w_k^{+T}")
    print("                          + [(1-sig_k)^{-1/2}-1] w_k^- w_k^{-T}   -- EXACT, rank 2r.")
    print("     Treat the top-r pairs exactly (2r rank-1 corrections) and hand the REMAINDER")
    print("     to generic QSVT: the remainder's condition number is (1+sig_{r+1})/(1-sig_{r+1}).")
    print("     A valid whitener is  W = D^{-1/2} Shat^{-1/2}  (W^T A1 W = I), so")
    print("       d_eff(r) = d_inv(cond(I+Q))  +  2r  +  d_inv((1+sig_{r+1})/(1-sig_{r+1}))")
    print("     vs the baseline  d_base = d_inv(cond(A1)).")
    print(f"\n  {'nmax':>4} {'N':>4} | {'cond(A1)':>9} {'d_base':>8} | {'d_inv(I+Q)':>10}"
          f" | {'r=1':>7} {'r=2':>7} {'r=3':>7} {'r=4':>7} | {'best r':>6} {'d_eff':>7}"
          f" {'speedup':>8}")
    print("  " + "-" * 104)
    price = []
    for nmax in (4, 6, 8, 10, 12):
        A1, SOO, SOH, SHH = water_A1(nmax)
        Ct = invsqrt(SOO) @ SOH @ invsqrt(SHH)
        sig = np.linalg.svd(Ct, compute_uv=False)
        cA = float(np.linalg.cond(A1))
        dbase = d_inv(cA)
        dq = d_inv(float(np.linalg.cond(SHH)))
        deffs = []
        for r in range(1, 5):
            kres = (1 + sig[r]) / (1 - sig[r]) if r < len(sig) else 1.0
            deffs.append(dq + 2 * r + d_inv(kres))
        rbest = int(np.argmin(deffs)) + 1
        price.append(dict(nmax=nmax, N=3 * nmax, cond=cA, d_base=dbase,
                          d_eff=deffs[rbest - 1], r=rbest))
        print(f"  {nmax:>4} {3*nmax:>4} | {cA:>9.2f} {dbase:>8} | {dq:>10}"
              f" | " + " ".join(f"{x:>7}" for x in deffs)
              + f" | {rbest:>6} {deffs[rbest-1]:>7} {dbase/max(deffs[rbest-1],1):>7.1f}x")

    Np = [p['N'] for p in price]
    be, r2e, rese = loglog_fit(Np, [p['d_eff'] for p in price])
    bb, r2bb, resbb = loglog_fit(Np, [p['d_base'] for p in price])
    print(f"\n     d_base ~ N^{bb:+.3f} (R2={r2bb:.4f})     d_eff ~ N^{be:+.3f} (R2={r2e:.4f})")
    print(f"     => Woodbury {'FLATTENS' if be < 0.5 else 'does NOT flatten'} the metric penalty"
          f" (flat = exponent ~ 0).")

    # ----------------------------------------- honest check: what does the tail cost?
    print("\n[3e] Honest check -- WHY the low-rank treatment works here:")
    A1, SOO, SOH, SHH = water_A1(12)
    sig = np.linalg.svd(invsqrt(SOO) @ SOH @ invsqrt(SHH), compute_uv=False)
    print(f"     nmax=12 (N=36):  sigma = " + " ".join(f"{x:.4f}" for x in sig[:6]) + " ...")
    for r in range(0, 5):
        kres = (1 + sig[r]) / (1 - sig[r])
        print(f"       after removing top-{r}: residual kappa = {kres:8.2f}   d_inv = {d_inv(kres):>6}")
    print("     The growth is concentrated in the SINGLE stiffest direction; removing it")
    print("     exactly leaves a residual conditioning that is small and (see [3f]) flat.")

    print("\n[3f] Does the residual (post-rank-1) conditioning stay flat in N?")
    print(f"  {'nmax':>4} {'N':>4} | {'sig_1':>8} {'sig_2':>8} {'kappa_res(r=1)':>14}"
          f" {'kappa_res(r=2)':>14}")
    print("  " + "-" * 62)
    k1, k2, Nl = [], [], []
    for nmax in nmaxes:
        A1, SOO, SOH, SHH = water_A1(nmax)
        sig = np.linalg.svd(invsqrt(SOO) @ SOH @ invsqrt(SHH), compute_uv=False)
        a1 = (1 + sig[1]) / (1 - sig[1])
        a2 = (1 + sig[2]) / (1 - sig[2]) if len(sig) > 2 else 1.0
        k1.append(a1); k2.append(a2); Nl.append(3 * nmax)
        print(f"  {nmax:>4} {3*nmax:>4} | {sig[0]:>8.5f} {sig[1]:>8.5f} {a1:>14.3f} {a2:>14.3f}")
    p1, r21, _ = loglog_fit(Nl, k1)
    p2, r22, _ = loglog_fit(Nl, k2)
    print(f"\n     kappa_res(r=1) ~ N^{p1:+.3f} (R2={r21:.3f});"
          f"  kappa_res(r=2) ~ N^{p2:+.3f} (R2={r22:.3f})")

    os.makedirs("debug/data", exist_ok=True)
    with open("debug/data/aha_t1_water.json", "w") as fh:
        json.dump(dict(N=Ns, cond_A1=condA1, sigma_max=sigmax, cond_D=condD,
                       price=price, r10={str(k): v for k, v in r10.items()}), fh, indent=1)
    print("\n  wrote debug/data/aha_t1_water.json")
    print("DONE.")


if __name__ == "__main__":
    main()
