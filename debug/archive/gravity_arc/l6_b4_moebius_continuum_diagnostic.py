"""B4 -- Mobius alpha>1 mechanism: continuum-extrapolation diagnostic.

The alpha>1 conical slope was validated (v3.19.0 Track 5; task #25) as
    slope(alpha) = Delta_K / (1/alpha - alpha),  Delta_K = K_wedge - alpha*K_disk,
    recovery = slope / (-1/12) = alpha / (2*alpha - 1)   [Mobius]
using the FD azimuthal substrate at fixed a=0.05 and single t=1.0.

This sprint (B1) showed the FD lattice carries an O(a) Dirichlet-boundary error,
and GD-2 flagged the alpha>1 coefficient as t~1-tuned. DIAGNOSTIC: re-measure the
alpha>1 recovery with the CLEAN machinery -- SPECTRAL azimuthal (exact
m_eff=(k+1/2)/alpha) + fine-a linear (O(a)) extrapolation to continuum + multiple
t -- and ask whether the Mobius form alpha/(2alpha-1) survives, or whether it was
an FD + single-t artifact.

Outcomes:
  - recovery(continuum) ~ alpha/(2alpha-1), t-stable  => Mobius is CONTINUUM-REAL
    (mechanism still open, but it is physics not a substrate artifact)
  - recovery -> 1 (SC scalar continuation)            => Mobius was a substrate artifact
  - recovery -> something else, t-stable              => new closed form
  - recovery t-UNSTABLE                               => not a clean continuum object
"""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh_tridiagonal

OUT = Path(__file__).parent / "data" / "l6_b4_moebius_continuum_diagnostic.json"
OUT.parent.mkdir(exist_ok=True)


def radial_eigs(m, R, N_rho):
    a = R / N_rho
    k = np.arange(1, N_rho + 1)
    rho = k * a
    diag = 2.0 / a**2 + (m**2 - 0.25) / rho**2
    off = -np.ones(N_rho - 1) / a**2
    return eigvalsh_tridiagonal(diag, off)


def h(m, t, R, N_rho):
    return float(np.sum(np.exp(-t * radial_eigs(abs(m), R, N_rho))))


def K_wedge_spectral(alpha, t, R, N_rho, K_az):
    """4 sum_{j>=0} h((j+1/2)/alpha, t); alpha=1 is the disk."""
    return 4.0 * sum(h((j + 0.5) / alpha, t, R, N_rho) for j in range(K_az))


def slope_at(alpha, t, R, N_rho):
    # azimuthal cutoff scales with alpha (m_eff smaller => more modes)
    K_az = int(np.ceil(alpha * 70))
    K_disk = K_wedge_spectral(1.0, t, R, N_rho, max(80, K_az))
    K_w = K_wedge_spectral(alpha, t, R, N_rho, K_az)
    Delta_K = K_w - alpha * K_disk
    return Delta_K / (1.0 / alpha - alpha)


def main():
    res = {}
    R = 10.0
    SC = -1.0 / 12.0
    print("=" * 76)
    print("B4 -- Mobius alpha>1 continuum-extrapolation diagnostic")
    print("=" * 76)
    alphas = [2.0, 3.0, 5.0]
    ts = [0.5, 1.0, 2.0]
    Nrs = [200, 400, 800]   # for O(a) linear extrapolation a -> 0

    rows = []
    for alpha in alphas:
        moebius = alpha / (2 * alpha - 1)
        print(f"\nalpha = {alpha}   Mobius recovery alpha/(2alpha-1) = {moebius:.5f}")
        print(f"  {'t':>5} {'recov(a=R/200)':>15} {'recov(400)':>12} {'recov(800)':>12} "
              f"{'recov(a->0)':>13} {'vs Mobius':>11} {'vs SC=1':>9}")
        for t in ts:
            recs = []
            for Nr in Nrs:
                s = slope_at(alpha, t, R, Nr)
                recs.append(s / SC)
            # linear-in-a extrapolation (O(a)); a = R/Nr, Nr doubles => a halves
            # use 2 finest: rec0 = 2*rec(fine) - rec(coarse)
            rec0 = 2 * recs[-1] - recs[-2]
            dev_moeb = (rec0 - moebius) / moebius
            dev_sc = rec0 - 1.0
            print(f"  {t:>5} {recs[0]:>15.5f} {recs[1]:>12.5f} {recs[2]:>12.5f} "
                  f"{rec0:>13.5f} {dev_moeb:>+11.4f} {dev_sc:>+9.4f}")
            rows.append(dict(alpha=alpha, t=t, moebius=moebius,
                             rec_coarse=recs[0], rec_mid=recs[1], rec_fine=recs[2],
                             rec_continuum=rec0, dev_vs_moebius=dev_moeb,
                             dev_vs_SC=dev_sc))
    res["rows"] = rows

    # ---- verdict ----
    print("\n" + "=" * 76)
    # t-stability of the continuum recovery per alpha, and match to Mobius
    by_alpha = {}
    for alpha in alphas:
        r = [x for x in rows if x["alpha"] == alpha]
        rec0s = [x["rec_continuum"] for x in r]
        spread = max(rec0s) - min(rec0s)
        mean_dev_moeb = np.mean([abs(x["dev_vs_moebius"]) for x in r])
        mean_dev_sc = np.mean([abs(x["dev_vs_SC"]) for x in r])
        by_alpha[alpha] = dict(rec0_mean=float(np.mean(rec0s)),
                               t_spread=float(spread),
                               mean_dev_moebius=float(mean_dev_moeb),
                               mean_dev_SC=float(mean_dev_sc))
        print(f"  alpha={alpha}: continuum recovery {np.mean(rec0s):.4f} "
              f"(t-spread {spread:.4f}); |dev| vs Mobius {mean_dev_moeb:.3f}, "
              f"vs SC=1 {mean_dev_sc:.3f}")
    res["by_alpha"] = by_alpha

    # classify
    avg_moeb = np.mean([v["mean_dev_moebius"] for v in by_alpha.values()])
    avg_sc = np.mean([v["mean_dev_SC"] for v in by_alpha.values()])
    avg_spread = np.mean([v["t_spread"] for v in by_alpha.values()])
    if avg_spread > 0.10:
        verdict = "T-UNSTABLE (not a clean continuum object)"
    elif avg_moeb < 0.05 and avg_moeb < avg_sc:
        verdict = "MOBIUS-CONTINUUM-REAL (survives clean extrapolation; mechanism open)"
    elif avg_sc < 0.05 and avg_sc < avg_moeb:
        verdict = "SC-CONTINUATION (Mobius was a substrate/t artifact)"
    else:
        verdict = "NEITHER (new continuum form; |dev_moeb|=%.3f |dev_sc|=%.3f)" % (avg_moeb, avg_sc)
    print(f"\n[Verdict] {verdict}")
    print(f"  avg |dev| vs Mobius = {avg_moeb:.3f}, vs SC=1 = {avg_sc:.3f}, "
          f"t-spread = {avg_spread:.3f}")
    res["verdict"] = verdict
    res["avg_dev_moebius"] = float(avg_moeb)
    res["avg_dev_SC"] = float(avg_sc)
    res["avg_t_spread"] = float(avg_spread)

    with OUT.open("w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"\nsaved {OUT}")
    print("=" * 76)


if __name__ == "__main__":
    main()
