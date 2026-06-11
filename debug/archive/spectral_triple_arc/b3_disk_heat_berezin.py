"""B3 -- disk propinquity backbone via the heat-kernel Berezin map.

The cigar-tip backbone is D^2_alpha (x) S^2 (constant warp). The outer tensor is
PURE_TENSOR (Paper 45/39); S^2 is Hawkins Berezin-Toeplitz (Paper 38 L4);
azimuthal S^1 is Toeplitz (Paper 44). The genuinely-new factor is the 2D
disk-with-cone D^2_alpha itself -- NOT a pure tensor of interval (x) S^1, because
the Laplacian's 1/rho^2 term couples radial and azimuthal.

Claim (B3 crux): the disk propinquity backbone transports via the HEAT-KERNEL
Berezin map B_Lambda(f) = P_Lambda e^{-t Delta/2} M_f e^{-t Delta/2} P_Lambda
(t ~ 1/Lambda), whose four L4 properties all follow from heat-semigroup facts
that hold on any compact manifold-with-boundary (cone included):
  (a) positivity   : e^{-tD/2} positivity-preserving, congruence preserves PSD
  (b) contractivity: ||e^{-tD/2}|| <= 1, ||P M_f P|| <= ||f||_inf
  (c) approx-id    : ||B_Lambda(f) - P_Lambda M_f P_Lambda|| -> 0 at rate O(1/Lambda)
                     (heat smoothing at scale t=1/Lambda); Lambda ~ n (Weyl) => O(1/n)
  (d) Lipschitz    : ||grad(e^{-tDelta}f)||_inf <= ||grad f||_inf (heat contraction)

This driver verifies (a)(b)(c) numerically on the flat disk (alpha=1) eigenbasis
and measures the (c) rate.  alpha != 1 (cone) extends because the cone heat
semigroup (Cheeger) is also positivity-preserving and contractive.

Disk Laplacian eigenbasis (Dirichlet at rho=R):
  psi_{k,j}(rho,phi) = N_{k,j} J_{|k|}(a_{|k|,j} rho/R) e^{i k phi},
  lambda_{k,j} = (a_{|k|,j}/R)^2,  a_{m,j} = j-th positive zero of J_m,
  N_{k,j} = 1/sqrt(pi R^2 J_{|k|+1}(a_{|k|,j})^2).
Radial test function f(rho) = (rho/R)^2 (smooth, 0<=f<=1, block-diagonal in k).
"""

import json
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros, jv
from numpy.polynomial.legendre import leggauss

OUT = Path(__file__).parent / "data" / "b3_disk_heat_berezin.json"
OUT.parent.mkdir(exist_ok=True)

R = 1.0  # disk radius (scale out)


def build_modes(K, Jmax, Lambda):
    """All (k,j) modes with |k|<=K, j<=Jmax, lambda<=Lambda. Returns list of
    dicts with k, j, lam, alpha(zero), norm."""
    modes = []
    for k in range(-K, K + 1):
        m = abs(k)
        zeros = jn_zeros(m, Jmax)
        for j, a in enumerate(zeros, start=1):
            lam = (a / R) ** 2
            if lam <= Lambda:
                norm = 1.0 / np.sqrt(np.pi * R**2 * jv(m + 1, a) ** 2)
                modes.append(dict(k=k, m=m, j=j, alpha=a, lam=lam, norm=norm))
    modes.sort(key=lambda d: d["lam"])
    return modes


def radial_Mf_block(modes_k, f_rho, ng=400):
    """Matrix of M_f restricted to a fixed-k radial block (f is radial).
    <psi_{k,j}|f|psi_{k,j'}> = 2 pi N_j N_j' int_0^R J_m(a_j r/R) f(r) J_m(a_j' r/R) r dr.
    """
    x, w = leggauss(ng)
    r = 0.5 * R * (x + 1.0)          # nodes on (0,R)
    wr = 0.5 * R * w
    n = len(modes_k)
    M = np.zeros((n, n))
    fr = f_rho(r)
    Jvals = []
    for d in modes_k:
        Jvals.append(jv(d["m"], d["alpha"] * r / R))
    for i in range(n):
        for j in range(n):
            integ = np.sum(Jvals[i] * fr * Jvals[j] * r * wr)
            M[i, j] = 2 * np.pi * modes_k[i]["norm"] * modes_k[j]["norm"] * integ
    return M


def main():
    res = {}
    print("=" * 74)
    print("B3 -- disk heat-kernel Berezin: L4 properties + approximate-id rate")
    print("=" * 74)

    f_rho = lambda r: (r / R) ** 2          # 0<=f<=1, ||f||_inf = 1
    f_inf = 1.0
    # ||grad f||_inf = max |2 rho/R^2| = 2/R
    gradf_inf = 2.0 / R

    print(f"\ntest f(rho) = (rho/R)^2,  ||f||_inf = {f_inf},  ||grad f||_inf = {gradf_inf}")

    # sweep the cutoff Lambda; t = c/Lambda (heat smoothing at mode scale)
    c_heat = 1.0
    rows = []
    Lambdas = [50.0, 100.0, 200.0, 400.0, 800.0, 1600.0]
    K, Jmax = 12, 40
    print(f"\n  {'Lambda':>8} {'n_modes':>8} {'t=c/L':>9} {'min eig B':>12} "
          f"{'||B||':>9} {'||B - PMfP||':>13}")
    for Lam in Lambdas:
        modes = build_modes(K, Jmax, Lam)
        n = len(modes)
        t = c_heat / Lam
        # block-diagonal in k: assemble full B and P M_f P blockwise, take norms
        ks = sorted(set(d["k"] for d in modes))
        min_eig = np.inf
        normB = 0.0
        approx_id = 0.0
        Bnorm_f = f_inf  # for contractivity ratio
        for k in ks:
            mk = [d for d in modes if d["k"] == k]
            Mf = radial_Mf_block(mk, f_rho)
            lam = np.array([d["lam"] for d in mk])
            Dw = np.diag(np.exp(-t * lam / 2.0))       # e^{-t Delta/2}
            B = Dw @ Mf @ Dw                            # heat Berezin (this block)
            # properties
            ev = np.linalg.eigvalsh(0.5 * (B + B.T))
            min_eig = min(min_eig, float(ev.min()))
            normB = max(normB, float(np.max(np.abs(np.linalg.eigvalsh(0.5*(B+B.T))))))
            approx_id = max(approx_id, float(np.linalg.norm(B - Mf, 2)))
        rows.append(dict(Lambda=Lam, n_modes=n, t=t,
                         min_eig_B=min_eig, normB=normB, approx_id=approx_id))
        print(f"  {Lam:>8.0f} {n:>8} {t:>9.5f} {min_eig:>12.3e} "
              f"{normB:>9.5f} {approx_id:>13.5f}")

    res["rows"] = rows

    # ---- rate fit: approx_id ~ C * Lambda^{-p} ----
    Ls = np.array([r["Lambda"] for r in rows])
    aid = np.array([r["approx_id"] for r in rows])
    p = float(np.polyfit(np.log(Ls), np.log(aid), 1)[0])
    print(f"\n  approximate-identity rate: ||B-PMfP|| ~ Lambda^{p:.3f}")
    print(f"  (Lambda ~ n by Weyl on the disk, so gamma_n -> 0 at this power)")

    # ---- verdict on the four L4 properties ----
    pos_ok = all(r["min_eig_B"] > -1e-9 for r in rows)        # (a) positivity
    contr_ok = all(r["normB"] <= f_inf + 1e-9 for r in rows)  # (b) contractivity
    aid_ok = aid[-1] < aid[0] and p < -0.3                    # (c) approx id -> 0
    print("\n" + "=" * 74)
    print(f"[L4 properties on the disk]")
    print(f"  (a) positivity  (f>=0 => B PSD):    {'PASS' if pos_ok else 'FAIL'} "
          f"(min eig >= {min(r['min_eig_B'] for r in rows):.2e})")
    print(f"  (b) contractivity (||B|| <= ||f||): {'PASS' if contr_ok else 'FAIL'} "
          f"(max ||B|| = {max(r['normB'] for r in rows):.5f} <= {f_inf})")
    print(f"  (c) approx identity -> 0, rate:     {'PASS' if aid_ok else 'FAIL'} "
          f"(Lambda^{p:.2f})")
    print(f"  (d) Lipschitz: ||grad e^{{-tD}}f|| <= ||grad f|| -- heat contraction "
          f"(analytic, not computed here)")
    verdict = ("B3-CRUX-COMPRESSES (heat-kernel Berezin gives all L4 properties; "
               "rate O(Lambda^p))" if (pos_ok and contr_ok and aid_ok)
               else "B3-CRUX-CHECK")
    print(f"\n[Verdict] {verdict}")
    res["rate_p"] = p
    res["L4_positivity"] = pos_ok
    res["L4_contractivity"] = contr_ok
    res["L4_approx_id"] = aid_ok
    res["verdict"] = verdict
    with OUT.open("w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"\nsaved {OUT}")
    print("=" * 74)


if __name__ == "__main__":
    main()
