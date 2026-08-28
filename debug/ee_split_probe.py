"""PROBE: does the e-e tensor admit a partial split
          (geometric factor) x (dynamical factor) + irreducible remainder?

Follow-on to rem:minimal_presentation (PI-directed).  s-only shared-k Sturmians.

The split is EXACT, not a fit.  The L=0 Coulomb kernel is 1/max(r1,r2), and
    1/max = min(1/r1, 1/r2) = (1/r1 + 1/r2)/2  -  |1/r1 - 1/r2|/2
so
    g  =  g_sep - g_W,
    g_sep[i,j,k,l] = 1/2 [ V_ik S_jl + S_ik V_jl ],   V = <chi|1/r|chi> = k diag(1/n)
    g_W  from the kernel  K_W = |1/r1 - 1/r2|/2  >= 0.

PRE-REGISTERED PREDICTIONS
  P1  dynamical factor is TOTAL:  g(k) = k * g(1)  (machine);  g(1) rational
      (sympy exact anchors: (11|11) = 5/8 at k=1).
  P2  split identity exact;  g_sep matricization (ik),(jl) has RANK 2 exactly;
      V-formula (k/n) delta matches the grid <1/r> to grid precision.
  P3  FCI with g_sep ONLY is exactly uncorrelated: it equals the FCI of the
      dressed one-body operator h1 + (N-1)/2 * V with NO two-body term at all.
      Hence ALL correlation lives in W.
  P4  W is not small: on the (1s,1s) ground pair <W> = 3k/8 vs <g> = 5k/8 (60%).
      The split's value is structural, not perturbative.
  P5  (measured, open)  eigenvalue decay of W's matricization; correlation
      recovered vs rank of W.  Fast decay => the remainder is compressible and
      the partial split is USEFUL; flat decay => W is irreducibly high-rank.
  P6  dictionary sharpening: with the honest geometric dictionary (words in
      {I,S} PLUS V), the fit residual of g isolates to exactly W/|g|; W itself
      stays unfittable (the irreducible part is W, not an artifact of a poor
      dictionary).
"""
import io
import itertools
import json
import os
import sys

import numpy as np
from scipy.linalg import eigh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from geovac import transcorrelated_sturmian as TC  # noqa: E402

out = {}


# ---------------------------------------------------------------------------
# exact radial kernels (no angular quadrature needed at L=0)
# ---------------------------------------------------------------------------
def kernels_exact(r):
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    K_full = 1.0 / np.maximum(R1, R2)
    K_sep = 0.5 * (1.0 / R1 + 1.0 / R2)
    K_W = 0.5 * np.abs(1.0 / R1 - 1.0 / R2)
    return K_full, K_sep, K_W


def build_all(ns, k, Z, Ng=700):
    r, wr = TC.make_grid(k, Ng=Ng)
    S, h1, Rtab, W = TC.build_one_body(ns, r, wr, k, Z)
    K_full, K_sep, K_W = kernels_exact(r)
    D = {(i, kk): Rtab[i + 1][0] * Rtab[kk + 1][0] * W
         for i in range(ns) for kk in range(ns)}

    def eri_of(K):
        g = np.zeros((ns, ns, ns, ns))
        for i in range(ns):
            for j in range(ns):
                for kk in range(ns):
                    for l in range(ns):
                        g[i, j, kk, l] = D[(i, kk)] @ K @ D[(j, l)]
        return g

    g = eri_of(K_full)
    g_sep = eri_of(K_sep)
    g_W = eri_of(K_W)
    Vr = np.array([[np.sum(Rtab[i + 1][0] * Rtab[j + 1][0] * r * wr)
                    for j in range(ns)] for i in range(ns)])   # <chi|1/r|chi>
    return r, wr, S, h1, Vr, g, g_sep, g_W


def pairmat(g, ns):
    """matricization with electron-1 pair (i,k) as row, electron-2 pair (j,l) as col."""
    return g.transpose(0, 2, 1, 3).reshape(ns * ns, ns * ns)


# ---------------------------------------------------------------------------
ns, k, Z, Ng = 5, 2.0, 2.0, 700
r, wr, S, h1, Vr, g, g_sep, g_W = build_all(ns, k, Z, Ng)

print("=" * 78)
print("P1 -- the dynamical factor is TOTAL:  g(k) = k g(1);  g(1) is rational")
print("=" * 78)
_, _, _, _, _, g_k1, _, _ = build_all(ns, 1.0, Z, Ng)
ratio_dev = float(np.abs(g / k - g_k1).max() / np.abs(g_k1).max())
print(f"  max |g(k)/k - g(1)| / |g(1)|  =  {ratio_dev:.2e}   (k = {k})")
import sympy as sp                                     # exact anchors at k=1
rr1, rr2 = sp.symbols("r1 r2", positive=True)
R1s = 2 * sp.exp(-rr1)
R1s2 = 2 * sp.exp(-rr2)
inner = sp.integrate(R1s2**2 * rr2**2 / rr1, (rr2, 0, rr1)) \
      + sp.integrate(R1s2**2 * rr2**2 / rr2, (rr2, rr1, sp.oo))
exact_1111 = sp.integrate(R1s**2 * rr1**2 * inner, (rr1, 0, sp.oo))
print(f"  sympy exact (11|11) at k=1  =  {exact_1111}   "
      f"(numeric {g_k1[0,0,0,0]:.10f}, dev {abs(float(exact_1111)-g_k1[0,0,0,0]):.1e})")
out["P1"] = dict(scale_dev=ratio_dev, exact_1111=str(exact_1111),
                 numeric_1111=float(g_k1[0, 0, 0, 0]))

print()
print("=" * 78)
print("P2 -- the exact split;  rank(g_sep) = 2;  V = k diag(1/n)")
print("=" * 78)
split_dev = float(np.abs(g - (g_sep - g_W)).max() / np.abs(g).max())
n = np.arange(1, ns + 1)
V_label = k * np.diag(1.0 / n)
v_dev = float(np.abs(Vr - V_label).max())
g_sep_label = 0.5 * (np.einsum("ik,jl->ijkl", V_label, S)
                     + np.einsum("ik,jl->ijkl", S, V_label))
label_dev = float(np.abs(g_sep - g_sep_label).max() / np.abs(g_sep).max())
sv = np.linalg.svd(pairmat(g_sep, ns), compute_uv=False)
print(f"  max |g - (g_sep - g_W)| / |g|      = {split_dev:.2e}")
print(f"  max |<1/r> - (k/n) delta|          = {v_dev:.2e}")
print(f"  max |g_sep - labels-x-metric form| = {label_dev:.2e}")
print(f"  singular values of g_sep pairmat: {sv[0]:.3f}, {sv[1]:.3f}, "
      f"then {sv[2]:.1e} (rank 2: {'CONFIRMED' if sv[2] < 1e-10 * sv[0] else 'FAILED'})")
out["P2"] = dict(split_dev=split_dev, v_dev=v_dev, label_dev=label_dev,
                 sv3_over_sv1=float(sv[2] / sv[0]))

print()
print("=" * 78)
print("P3 -- FCI with g_sep only is EXACTLY uncorrelated (one-body in disguise)")
print("=" * 78)
X = TC.lowdin(S)
h1o = TC.transform_1(h1, X)
Vo = TC.transform_1(Vr, X)
res_p3 = {}
for name, N_elec in [("He-like (N=2)", 2), ("Li-like (N=3)", 3)]:
    nso = 2 * ns
    dets, didx = TC.make_dets(nso, N_elec)
    go_sep = TC.transform_2(g_sep, X)
    hso = TC.h_spin(h1o, nso)
    asym = TC.asym_from_phys(go_sep, nso)
    Hfci = TC.build_H(dets, didx, hso, asym, nso)
    E_sep = float(eigh(Hfci, eigvals_only=True)[0])
    # the dressed ONE-BODY problem: h_eff = h1 + (N-1)/2 * V, NO two-body at all
    h_eff = TC.h_spin(h1o + 0.5 * (N_elec - 1) * Vo, nso)
    H1b = TC.build_H(dets, didx, h_eff, np.zeros_like(asym), nso)
    E_1b = float(eigh(H1b, eigvals_only=True)[0])
    print(f"  {name}:  FCI(g_sep) = {E_sep:.12f}   dressed-1-body = {E_1b:.12f}"
          f"   diff = {abs(E_sep - E_1b):.2e}")
    res_p3[name] = dict(E_sep=E_sep, E_1b=E_1b, diff=abs(E_sep - E_1b))
out["P3"] = res_p3

print()
print("=" * 78)
print("P4 -- W is structural, not small")
print("=" * 78)
w_frac_1s = g_W[0, 0, 0, 0] / g[0, 0, 0, 0]
print(f"  <1s^2|W|1s^2> / <1s^2|g|1s^2> = {w_frac_1s:.6f}   "
      f"(predicted (3k/8)/(5k/8) = 0.600000)")
frob = float(np.linalg.norm(g_W) / np.linalg.norm(g))
print(f"  ||W||_F / ||g||_F             = {frob:.4f}")
out["P4"] = dict(w_frac_1s=float(w_frac_1s), frob_frac=frob)

print()
print("=" * 78)
print("P6 -- dictionary sharpening: the irreducible part IS W")
print("=" * 78)
def fit_resid(target, mats):
    M = np.array([m.ravel() for m in mats]).T
    coef, *_ = np.linalg.lstsq(M, target.ravel(), rcond=None)
    return float(np.linalg.norm(target.ravel() - M @ coef)
                 / np.linalg.norm(target.ravel()))

Gm, Sm, Wm = pairmat(g, ns), pairmat(g_sep, ns), pairmat(g_W, ns)
base_S = [np.eye(ns), S, S @ S, S @ S @ S]
base_SV = base_S + [V_label, V_label @ S, S @ V_label, V_label @ V_label]
words_S = [np.kron(A, B) for A in base_S for B in base_S]
words_SV = [np.kron(A, B) for A in base_SV for B in base_SV]
for nm, T in [("g   ", Gm), ("gsep", Sm), ("W   ", Wm)]:
    rS, rSV = fit_resid(T, words_S), fit_resid(T, words_SV)
    print(f"  {nm}: residual vs words(S) = {rS:.4f}   vs words(S,V) = {rSV:.4f}")
    out.setdefault("P6", {})[nm.strip()] = dict(words_S=rS, words_SV=rSV)

print()
print("=" * 78)
print("P5 -- spectrum of W and correlation recovered vs rank (the payoff curve)")
print("=" * 78)
lam, U = np.linalg.eigh(Wm)
order = np.argsort(-np.abs(lam))
lam_o = lam[order]
print("  leading |eigenvalues| of W pairmat: "
      + ", ".join(f"{abs(v):.4f}" for v in lam_o[:8]))
tail = np.array([np.sum(lam_o[m:] ** 2) for m in range(len(lam_o))])
print("  Frobenius tail fraction after m terms: "
      + ", ".join(f"m={m}:{np.sqrt(tail[m]/tail[0]):.3f}" for m in (1, 2, 4, 8, 12)))
out["P5_spectrum"] = [float(v) for v in lam_o]

def fci_energy(g_tensor, N_elec):
    nso = 2 * ns
    dets, didx = TC.make_dets(nso, N_elec)
    go = TC.transform_2(g_tensor, X)
    hso = TC.h_spin(h1o, nso)
    asym = TC.asym_from_phys(go, nso)
    return float(eigh(TC.build_H(dets, didx, hso, asym, nso), eigvals_only=True)[0])

res_curve = {}
for name, N_elec in [("He-like N=2", 2), ("Li-like N=3", 3)]:
    E_full = fci_energy(g, N_elec)
    E_sep_only = fci_energy(g_sep, N_elec)
    corr_total = E_full - E_sep_only          # everything W contributes to FCI
    rows = []
    for m in (0, 1, 2, 4, 8, 12, ns * ns):
        Wm_r = (U[:, order[:m]] * lam_o[:m]) @ U[:, order[:m]].T
        g_r = g_sep - Wm_r.reshape(ns, ns, ns, ns).transpose(0, 2, 1, 3)
        # enforce the real-orbital pair symmetries after truncation
        g_r = 0.25 * (g_r + g_r.transpose(2, 1, 0, 3).transpose(0, 1, 2, 3)
                      + g_r.transpose(0, 3, 2, 1) + g_r.transpose(2, 3, 0, 1))
        E_r = fci_energy(g_r, N_elec)
        rows.append((m, E_r, E_r - E_full))
    print(f"  {name}:  E_full = {E_full:.8f}   E(sep-only) = {E_sep_only:.8f}"
          f"   (W worth {1000 * (E_sep_only - E_full):+.2f} mHa)")
    for m, E_r, d in rows:
        print(f"      rank {m:>3}:  E = {E_r:.8f}   E - E_full = {1000 * d:+9.4f} mHa")
    res_curve[name] = dict(E_full=E_full, E_sep=E_sep_only,
                           curve=[(m, E, d) for m, E, d in rows])
out["P5_curve"] = res_curve

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/ee_split_probe.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, default=float)
print("\nwrote debug/data/ee_split_probe.json")
