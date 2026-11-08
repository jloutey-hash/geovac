"""aha Track 1 / STEPS 2 + 4 -- the sigma_max(N, R, Z) growth law, and the
non-monotone crossover between the two walls.

Track (i)  SW  : shared-scale Shibuya-Wulfman metric.  Structurally a ONE-parameter
                 family in s = k R (the nuclear charges do not appear in S at all),
                 so the (R, Z) collapse question has an exact answer on this track.
Track (ii) GOS : Goscinskian / hydrogenic L^2 metric, a = Z_c/n per center -- the
                 mixed-charge layer (v4.94.0).  Intra blocks are the identity
                 (hydrogenic radial functions at one Z are orthonormal), so the
                 same theorem applies and Z-dependence is genuine.  Exact
                 Mulliken/Ruedenberg A_m(p)/B_n(q) closed form (machine precision).

Reported per (R, Z):  the law  1 - sigma_max ~ N^(-alpha)  with R^2 and max log
residual, alongside the directly-fitted cond(S) ~ N^beta, and the commutator norm.

Run:  python debug/aha_t1_sweep.py
"""
from __future__ import annotations

import sys, os, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from debug.aha_t1_core import (sw_block_trapz, goscinskian_cross_block, hyd_s_coeffs,
                               assemble_two_center, sigma_spectrum, cond_from_sigma,
                               commutator_from_sigma, loglog_fit)

NMAXES = list(range(2, 11))          # per-center basis size; N = 2 * nmax
RGRID = [1.4, 2.0, 3.0, 4.0, 6.0, 10.0]
ZPAIRS = [(1.0, 1.0), (1.0, 3.0), (1.0, 8.0)]


def _one_center_orthonormality(nmax: int, Z: float) -> float:
    """max |<chi_n|chi_m> - delta| for the Goscinskian a = Z/n set, by 1D quadrature."""
    r = np.linspace(1e-9, 400.0, 400001)
    F = []
    for n in range(1, nmax + 1):
        d = hyd_s_coeffs(n, Z / n).astype(float)
        F.append(np.polyval(d[::-1], r) * np.exp(-(Z / n) * r))
    F = np.array(F)
    G = (F * r ** 2) @ F.T * (r[1] - r[0])
    return float(np.abs(G - np.eye(nmax)).max())


def analyse(Cfull: np.ndarray, label: str, rows: list) -> dict:
    """Sweep leading submatrices of Cfull (nmax = 2..) and fit the laws."""
    Ns, one_m_sig, conds, comms, smaxs = [], [], [], [], []
    for nmax in NMAXES:
        C = Cfull[:nmax, :nmax]
        sig = sigma_spectrum(C)
        smax = float(sig.max())
        S = assemble_two_center(C)
        cd = float(np.linalg.cond(S))
        cm = commutator_from_sigma(sig)
        Ns.append(2 * nmax); smaxs.append(smax); one_m_sig.append(1 - smax)
        conds.append(cd); comms.append(cm)
        rows.append(dict(case=label, nmax=nmax, N=2 * nmax, sigma_max=smax,
                         one_minus_sigma=1 - smax, cond=cd, comm=cm,
                         sigma_top3=[float(x) for x in sig[:3]]))
    a, r2a, resa = loglog_fit(Ns, one_m_sig)
    b, r2b, resb = loglog_fit(Ns, conds)
    return dict(label=label, N=Ns, smax=smaxs, one_m=one_m_sig, cond=conds, comm=comms,
                alpha=-a, r2_alpha=r2a, res_alpha=resa, beta=b, r2_beta=r2b, res_beta=resb)


def main() -> None:
    rows: list = []
    print("=" * 100)
    print("STEP 2 -- sigma_max growth law.  N = 2*nmax (total two-center dimension).")
    print("          Fitted:  1 - sigma_max ~ N^(-alpha)   and   cond(S) ~ N^(+beta).")
    print("          Theorem => beta = alpha asymptotically (cond = (1+s)/(1-s) ~ 2/(1-s)).")
    print("=" * 100)

    # ------------------------------------------------------------------ track (i) SW
    print("\n### TRACK (i): Shibuya-Wulfman shared-scale metric.  s = k*R is the ONLY")
    print("### parameter -- the SW metric contains no nuclear charge, so the (R,Z)")
    print("### collapse is exact and one-dimensional by construction.")
    print(f"\n  {'R (=s at k=1)':>13} | {'alpha (1-sig_max)':>17} {'R2':>7} {'maxres':>8}"
          f" | {'beta (cond)':>11} {'R2':>7} {'maxres':>8} | {'1-sig_max @N=20':>15}"
          f" {'cond @N=20':>10} {'||[PA,PB]|| @N=20':>17}")
    print("  " + "-" * 118)
    sw_res = {}
    for R in RGRID:
        Cf = sw_block_trapz(R, max(NMAXES), "S")     # k = 1  => s = R
        d = analyse(Cf, f"SW s={R}", rows)
        sw_res[R] = d
        print(f"  {R:>13.1f} | {d['alpha']:>17.4f} {d['r2_alpha']:>7.4f} {d['res_alpha']:>8.4f}"
              f" | {d['beta']:>11.4f} {d['r2_beta']:>7.4f} {d['res_beta']:>8.4f}"
              f" | {d['one_m'][-1]:>15.6f} {d['cond'][-1]:>10.2f} {d['comm'][-1]:>17.6f}")

    print("\n  Paper-60 reference exponents:  cond(SW full) ~ N^1.85 (paper) / N^1.81 (memo).")
    print(f"  This run at R=2.0 (the paper's H2+ geometry): beta = {sw_res[2.0]['beta']:.4f},"
          f" alpha = {sw_res[2.0]['alpha']:.4f}.")

    # ------------------------------------------------- collapse test for the SW family
    print("\n  [collapse] Is  1 - sigma_max(N, s) = f(N * g(s))?  A single-variable collapse")
    print("  requires alpha to be s-INDEPENDENT (parallel log-log lines).  Measured:")
    al = np.array([sw_res[R]['alpha'] for R in RGRID])
    print(f"     alpha(s) over s = {RGRID}:")
    print(f"        {np.array2string(al, precision=4)}")
    print(f"     spread: min={al.min():.4f} max={al.max():.4f} ratio={al.max()/al.min():.3f}")
    # prefactor A(s): 1-sig = A N^-alpha  ->  collapse variable x = N * A^(1/alpha)
    print("\n     Collapse attempt with a COMMON exponent alpha_bar (weighted mean):")
    ab = float(al.mean())
    print(f"        alpha_bar = {ab:.4f}")
    print(f"     {'s':>6} | {'A(s) = (1-sig)*N^alpha_bar (should be N-independent)':>56}")
    print("     " + "-" * 68)
    collapse_spread = []
    for R in RGRID:
        d = sw_res[R]
        A = np.array(d['one_m']) * np.array(d['N'], float) ** ab
        collapse_spread.append(A.max() / A.min())
        print(f"     {R:>6.1f} | " + "  ".join(f"{x:.4f}" for x in A)
              + f"   [max/min = {A.max()/A.min():.2f}]")
    print(f"\n     => worst within-s drift of A over N = {max(collapse_spread):.2f}x."
          f"  A single-parameter collapse")
    print(f"        1-sigma_max = F(N * g(s)) with a COMMON exponent is "
          f"{'SUPPORTED' if max(collapse_spread) < 1.3 else 'NOT clean'}"
          f" (>1.3x drift = the exponent genuinely moves with s).")

    # --------------------------------------------------------------- track (ii) GOS
    print("\n\n### TRACK (ii): Goscinskian / hydrogenic L^2 metric, a = Z_c/n per center.")
    print("### Intra-center orthonormality check (must be I for the theorem to apply):")
    for Z in (1.0, 3.0, 8.0):
        print(f"     Z={Z:>4.1f}, nmax=10:  max|<chi_n|chi_m> - delta| = "
              f"{_one_center_orthonormality(10, Z):.2e}   [1D quadrature]")
    print(f"\n  {'ZA,ZB':>7} {'R':>5} | {'alpha':>8} {'R2':>7} {'maxres':>8}"
          f" | {'beta(cond)':>10} {'R2':>7} {'maxres':>8}"
          f" | {'1-sig @N=20':>11} {'cond @N=20':>10} {'comm @N=20':>10}")
    print("  " + "-" * 106)
    gos_res = {}
    for (ZA, ZB) in ZPAIRS:
        for R in RGRID:
            Cf = goscinskian_cross_block(max(NMAXES), ZA, ZB, R)
            d = analyse(Cf, f"GOS Z=({ZA:g},{ZB:g}) R={R}", rows)
            gos_res[(ZA, ZB, R)] = d
            print(f"  {ZA:g},{ZB:g}".rjust(9) + f" {R:>5.1f} | {d['alpha']:>8.4f}"
                  f" {d['r2_alpha']:>7.4f} {d['res_alpha']:>8.4f} | {d['beta']:>10.4f}"
                  f" {d['r2_beta']:>7.4f} {d['res_beta']:>8.4f} | {d['one_m'][-1]:>11.6f}"
                  f" {d['cond'][-1]:>10.2f} {d['comm'][-1]:>10.6f}")
        print()

    # --------------------------------------------------------------- STEP 4 crossover
    print("=" * 100)
    print("STEP 4 -- the two walls are NON-MONOTONE in each other.")
    print("  cond  = (1+sig_max)/(1-sig_max)      -> DIVERGES as sig_max -> 1")
    print("  comm  = max_k sig_k sqrt(1-sig_k^2)  -> <= 1/2 always, and -> 0 for the")
    print("          very direction (sig -> 1) that makes cond explode.")
    print("  The commutator therefore SATURATES at 1/2 (attained when some sig_k ~ 1/sqrt2)")
    print("  and carries NO scaling information: all the scaling is in 1 - sig_max.")
    print("=" * 100)
    print(f"\n  {'case':>22} {'N':>4} | {'sig_max':>9} {'1-sig_max':>10} {'cond':>10}"
          f" | {'comm':>8} {'sig* (argmax)':>13} {'|sig*-1/sqrt2|':>14}")
    print("  " + "-" * 100)
    inv = 1 / np.sqrt(2)
    for case_C, lab in [(sw_block_trapz(2.0, 10, "S"), "SW s=2.0"),
                        (sw_block_trapz(1.4, 10, "S"), "SW s=1.4"),
                        (goscinskian_cross_block(10, 1.0, 8.0, 1.4), "GOS Z=(1,8) R=1.4")]:
        for nmax in (2, 4, 6, 8, 10):
            C = case_C[:nmax, :nmax]
            sig = sigma_spectrum(C)
            f = sig * np.sqrt(np.maximum(0, 1 - sig ** 2))
            ka = int(np.argmax(f))
            print(f"  {lab:>22} {2*nmax:>4} | {sig.max():>9.6f} {1-sig.max():>10.6f}"
                  f" {cond_from_sigma(sig):>10.2f} | {f[ka]:>8.6f} {sig[ka]:>13.6f}"
                  f" {abs(sig[ka]-inv):>14.6f}")
        print()

    # ---------------------------------------------------------- correlation summary
    allc = np.array([r['cond'] for r in rows])
    allm = np.array([r['comm'] for r in rows])
    print(f"  Over ALL {len(rows)} (case, N) points in this sweep:")
    print(f"     cond(S) spans {allc.min():.2f} .. {allc.max():.2f}  ({allc.max()/allc.min():.0f}x)")
    print(f"     ||[P_A,P_B]|| spans {allm.min():.4f} .. {allm.max():.4f}"
          f"  (theoretical ceiling 0.5)")
    print(f"     Pearson r(log cond, comm) = {np.corrcoef(np.log(allc), allm)[0,1]:+.3f}")
    hi = allc > 20
    print(f"     restricted to cond > 20 ({hi.sum()} pts): comm = "
          f"{allm[hi].mean():.4f} +/- {allm[hi].std():.4f}   (pinned at the ceiling)")

    os.makedirs("debug/data", exist_ok=True)
    with open("debug/data/aha_t1_sweep.json", "w") as fh:
        json.dump(rows, fh, indent=1)
    print("\n  wrote debug/data/aha_t1_sweep.json")
    print("DONE.")


if __name__ == "__main__":
    main()
