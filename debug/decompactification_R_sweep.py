"""Two-center "decompactification" as a continuous function of R (DIAGNOSTIC).

Question (PI, 2026-09-05): is the compact->non-compact transition at two centers
-- the exact atomic-l label being lost (Paper 0 "Level 2"; Paper 8) -- a SMOOTH
function of internuclear distance R rather than a Boolean switch, and is there a
per-shell front that tracks r_n ~ R (n ~ sqrt(Z R))?

Reading under test: the LABEL is Boolean (axial m exact at all R; l stops being
conserved), but the COUPLING that l used to separate turns on continuously.

Primary system: H2+ (Z=1 per center), one electron, prolate-spheroidal EXACT.

Two deliverables from ONE R-sweep:

 D1  per-shell decompactification front.
     Measure = per-shell BICENTRIC character (option b of the prompt), made
     precise as the commutator of the two single-center shell projectors:
         char_n(R) = ||[P_A^(n), P_B^(n)]|| = S_n(R) * sqrt(1 - S_n(R)^2)
                   = (1/2)|sin 2 theta_n(R)|,   theta_n = arccos S_n(R),
     where S_n(R) = <chi^A_{n,0} | chi^B_{n,0}> is the EXACT two-center overlap of
     two proper hydrogenic ns orbitals (Mulliken/Ruedenberg auxiliaries;
     debug/aha_t1_core.two_center_s_overlap, validated to 1e-12 vs the analytic
     1s-1s form). char_n = 0 at theta=0 (S=1, united / R->0) AND at theta=90
     (S=0, atomic / R->inf); it PEAKS at theta=45 (S=1/sqrt2). The peak location
     R*(n) is the shell's decompactification point. Fit R*(n) ~ n^p and compare
     to the r_n = n^2/Z prediction (p=2, i.e. n*(R) ~ R^0.5).

 D2  four aligned curves on the SAME R-grid.
     (a) gamma(R): Paper 8 bond-sphere angle, closed form. p_R=1/R, p0=Z.
         cos gamma=(p0^2-p_R^2)/(p0^2+p_R^2), sin gamma=2 p0 p_R/(p0^2+p_R^2).
         gamma->pi at R=0 (antipodal/united), gamma->0 at R->inf (separated).
     (b) ||[P_A,P_B]|| and the principal angles between the two single-center
         orbital projectors. H2+ s-block {1s,2s,3s} sweep + the LiH sigma block
         {1s,2s,2p0} anchor (reproduces the v5.1.0 datum 0.500; 7.6/44.7/67.3).
     (c) aggregate l-mixing weight = sum_n char_n(R) (integrated from D1).
     (d) transcendental seed argument a(R) of the closed-form two-center engine
         (geovac/neumann_vee.py exp1 seed e^a E1(a)): Mulliken p = (R/2)(zA+zB).

 BONUS curve: direct spheroidal l-mixing of the H2+ ground-state eta-eigenstate
     (Paper 11 eta-equation, Legendre expansion). 1 - max_l a_l^2 and the Shannon
     entropy of {a_l^2}. This is the most literal realization of Paper 0's
     "l replaced by a separation parameter"; c^2(R) taken from the exact prolate
     solver (geovac.prolate_spheroidal_lattice).

Clean-room: driver in debug/, data in debug/data/. No paper/CLAUDE/test edits.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import exp1

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)


def _load(mod_name: str, rel_path: str):
    spec = importlib.util.spec_from_file_location(mod_name, os.path.join(REPO, rel_path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

core = _load("aha_t1_core", "debug/aha_t1_core.py")
fast = _load("fast_two_center_overlap", "debug/fast_two_center_overlap.py")

two_center_s_overlap = core.two_center_s_overlap        # exact hydrogenic s-s
commutator_from_sigma = core.commutator_from_sigma      # max_k sig_k sqrt(1-sig^2)
sigma_spectrum = core.sigma_spectrum
overlap_fast = fast.overlap_fast                        # (n,l,m) two-center, l<=1

Z = 1.0          # H2+ per-center charge
SQRT_HALF = 1.0 / np.sqrt(2.0)


# ======================================================================= helpers
def Sn(n: int, R: float, Zc: float = Z) -> float:
    """Exact two-center overlap of two hydrogenic ns orbitals (a = Zc/n)."""
    return two_center_s_overlap(n, Zc / n, n, Zc / n, R)


def char_n(n: int, R: float, Zc: float = Z) -> float:
    """Per-shell bicentric character = ||[P_A^(n),P_B^(n)]|| = S sqrt(1-S^2)."""
    s = min(max(Sn(n, R, Zc), 0.0), 1.0)
    return s * np.sqrt(max(0.0, 1.0 - s * s))


def gamma_bond(R: float, p0: float = Z) -> float:
    """Paper 8 bond-sphere polar angle gamma(R)."""
    pR = 1.0 / R
    cg = (p0 * p0 - pR * pR) / (p0 * p0 + pR * pR)
    return float(np.arccos(np.clip(cg, -1.0, 1.0)))


def block_commutator(states, R, Za=Z, Zb=Z):
    """||[P_A,P_B]|| and principal angles for a set of (n,l) states (m=0), both
    centers orthonormal within themselves -> sigma_k = svd(cross overlap)."""
    d = len(states)
    C = np.zeros((d, d))
    for i, (na, la) in enumerate(states):
        for j, (nb, lb) in enumerate(states):
            if la == 0 and lb == 0:
                C[i, j] = two_center_s_overlap(na, Za / na, nb, Zb / nb, R)
            else:
                C[i, j] = overlap_fast(Za, na, la, Zb, nb, lb, 0, R)
    sig = np.clip(np.linalg.svd(C, compute_uv=False), 0.0, 1.0)
    angles = np.degrees(np.arccos(sig))
    return float(commutator_from_sigma(sig)), sig.tolist(), angles.tolist()


# ------------------------------------------- direct spheroidal eta l-mixing (BONUS)
def _eta_ground_coeffs(c2: float, m: int = 0, lmax: int = 40, gerade: bool = True):
    """Ground (largest-A) eigenvector of the homonuclear eta-equation
    d/deta[(1-eta^2)G'] + (-A + c^2 eta^2 - m^2/(1-eta^2)) G = 0,
    expanded in NORMALIZED associated Legendre Pbar_l^m. Returns (A, coeffs, ls).
    Cross-checked bit-for-bit against prolate_spheroidal_lattice's A(c^2)."""
    c2 = float(c2)
    Lfull = lmax + 2

    def a_couple(l, mm):  # eta Pbar_l^m = a_l Pbar_{l+1}^m + a_{l-1} Pbar_{l-1}^m
        return np.sqrt((l + 1 - mm) * (l + 1 + mm) / ((2 * l + 1) * (2 * l + 3)))

    E = np.zeros((Lfull, Lfull))
    for l in range(m, Lfull - 1):
        a = a_couple(l, m)
        E[l, l + 1] = a
        E[l + 1, l] = a
    E2 = E @ E                                   # eta^2 operator
    ls = [l for l in range(m, lmax + 1) if ((l - m) % 2 == 0) == gerade]
    M = np.zeros((len(ls), len(ls)))
    for i, li in enumerate(ls):
        for j, lj in enumerate(ls):
            M[i, j] = c2 * E2[li, lj] - (li * (li + 1) if i == j else 0.0)
    w, V = np.linalg.eigh(M)
    k = int(np.argmax(w))
    return float(w[k]), V[:, k], ls


def eta_l_mixing(c2: float):
    """(participation_deficit, shannon_entropy_bits, A) of the ground eta-state."""
    A, vec, ls = _eta_ground_coeffs(c2)
    p = vec ** 2
    p = p / p.sum()
    part_deficit = float(1.0 - p.max())
    nz = p[p > 1e-300]
    entropy = float(-(nz * np.log2(nz)).sum())
    return part_deficit, entropy, A, {int(l): float(pi) for l, pi in zip(ls, p)}


# ======================================================================= main sweep
def run():
    out = {"system": "H2+", "Z_per_center": Z, "conventions": {
        "char_n": "||[P_A^n,P_B^n]|| = S_n sqrt(1-S_n^2), S_n=<ns_A|ns_B> exact",
        "gamma": "Paper 8: cos g=(p0^2-pR^2)/(p0^2+pR^2), pR=1/R, p0=Z=1",
        "a_seed": "Mulliken p=(R/2)(zA+zB), zeta=1 -> a=R ; seed e^a E1(a)",
        "front_prediction": "r_n=n^2/Z -> R*(n)~n^2 (p=2), n*(R)~sqrt(ZR) (q=0.5)",
    }}

    # ---- main physical R-grid (united-atom -> separated-atom for H2+) ----
    R_main = np.round(np.geomspace(0.2, 12.0, 45), 5)

    # (a) gamma(R)
    gamma = [gamma_bond(R) for R in R_main]

    # (b) H2+ commutator + principal angles, s-block {1s,2s,3s}
    sblock = [(1, 0), (2, 0), (3, 0)]
    comm_h2, ang_h2, sig_h2 = [], [], []
    for R in R_main:
        cval, sig, ang = block_commutator(sblock, R)
        comm_h2.append(cval)
        ang_h2.append(ang)
        sig_h2.append(sig)

    #     LiH sigma-block anchor {1s,2s,2p0}, Z=1 both (balanced bond block)
    lih_states = [(1, 0), (2, 0), (2, 1)]
    R_lih = np.round(np.geomspace(0.5, 10.0, 40), 5)
    comm_lih, ang_lih = [], []
    for R in R_lih:
        cval, sig, ang = block_commutator(lih_states, R)
        comm_lih.append(cval)
        ang_lih.append(ang)
    c_anchor, sig_anchor, ang_anchor = block_commutator(lih_states, 3.015)

    # (c) aggregate l-mixing = sum_n char_n(R), n=1..N_AGG
    N_AGG = 6
    agg = [float(sum(char_n(n, R) for n in range(1, N_AGG + 1))) for R in R_main]

    # (d) seed argument a(R) = (R/2)(zA+zB), zeta=1
    zeta = 1.0
    a_seed = [float(R / 2.0 * (zeta + zeta)) for R in R_main]      # = R
    seed_val = [float(np.exp(a) * exp1(a)) for a in a_seed]        # e^a E1(a)

    # BONUS: direct spheroidal eta l-mixing for the H2+ ground state at c^2(R)
    from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice
    eta_part, eta_entropy, eta_c2, eta_A = [], [], [], []
    for R in R_main:
        lat = ProlateSpheroidalLattice(float(R), 1, 1, N_xi=2500, xi_max=30.0,
                                       radial_method='spectral', n_basis=30)
        try:
            _, c2, _ = lat.solve()
        except Exception:
            c2 = float('nan')
        eta_c2.append(float(c2))
        if np.isfinite(c2):
            pd, ent, A, _ = eta_l_mixing(c2)
        else:
            pd = ent = A = float('nan')
        eta_part.append(pd)
        eta_entropy.append(ent)
        eta_A.append(A)

    # ---- D1: per-shell front heatmap + R*(n) + scaling fit ----
    n_shells = list(range(1, 9))
    R_front = np.round(np.geomspace(0.2, 22.0, 90), 5)
    heat = {}                                   # char_n(R) grid
    for n in n_shells:
        heat[n] = [float(char_n(n, R)) for R in R_front]

    # R*(n): S_n(R) = 1/sqrt2 (theta=45, peak character), unique (S_n decreasing)
    Rstar = {}
    for n in n_shells:
        hi = 3.0 * n * n + 5.0
        while Sn(n, hi) > SQRT_HALF:
            hi *= 1.5
        Rstar[n] = float(brentq(lambda R: Sn(n, R) - SQRT_HALF, 1e-4, hi, xtol=1e-9))

    ns = np.array(n_shells, float)
    Rs = np.array([Rstar[n] for n in n_shells])
    # fits: all shells and n>=2 (drop small-n curvature)
    def loglog(x, y):
        p = np.polyfit(np.log(x), np.log(y), 1)
        res = np.log(y) - np.polyval(p, np.log(x))
        r2 = 1 - np.sum(res ** 2) / np.sum((np.log(y) - np.log(y).mean()) ** 2)
        return float(p[0]), float(r2)
    p_all, r2_all = loglog(ns, Rs)
    p_hi, r2_hi = loglog(ns[1:], Rs[1:])

    # r_n = n^2/Z reference and the actual R*/r_n ratio (tests the premise)
    r_n = (ns ** 2) / Z
    ratio = Rs / r_n
    # decay-length reference: ell_n = n/Z ; R*/ell_n
    ell_n = ns / Z
    ratio_ell = Rs / ell_n

    front = {
        "n_shells": n_shells,
        "R_front": R_front.tolist(),
        "char_heatmap": {int(n): heat[n] for n in n_shells},
        "R_star": {int(n): Rstar[n] for n in n_shells},
        "r_n_meanradius": r_n.tolist(),
        "ratio_Rstar_over_rn": ratio.tolist(),
        "ell_n_decaylength": ell_n.tolist(),
        "ratio_Rstar_over_elln": ratio_ell.tolist(),
        "fit_Rstar_vs_n_exponent_all": p_all, "fit_R2_all": r2_all,
        "fit_Rstar_vs_n_exponent_nge2": p_hi, "fit_R2_nge2": r2_hi,
        "inverse_n_of_R_exponent_nge2": 1.0 / p_hi,
        "predicted_Rstar_exponent": 2.0,
        "predicted_n_of_R_exponent": 0.5,
    }

    out.update({
        "R_main": R_main.tolist(),
        "curve_a_gamma_rad": gamma,
        "curve_a_gamma_deg": [float(np.degrees(g)) for g in gamma],
        "curve_b_h2plus": {"block": sblock, "commutator": comm_h2,
                           "principal_angles_deg": ang_h2, "sigma": sig_h2},
        "curve_b_lih": {"block": lih_states, "R": R_lih.tolist(),
                        "commutator": comm_lih, "principal_angles_deg": ang_lih,
                        "anchor_R3.015": {"commutator": c_anchor,
                                          "sigma": sig_anchor,
                                          "principal_angles_deg": ang_anchor}},
        "curve_c_aggregate_lmixing": agg, "N_aggregate_shells": N_AGG,
        "curve_d_seed_arg_a": a_seed, "curve_d_seed_val_eaE1a": seed_val,
        "bonus_eta_lmixing": {"c2": eta_c2, "participation_deficit": eta_part,
                              "entropy_bits": eta_entropy, "A": eta_A},
        "front": front,
    })

    # -------- console summary --------
    print("=" * 74)
    print("DECOMPACTIFICATION R-SWEEP  (H2+, Z=1 per center)")
    print("=" * 74)
    print("\nD2(b) LiH anchor {1s,2s,2p0} @ R=3.015:")
    print("   sigma        =", [round(s, 4) for s in sig_anchor])
    print("   angles (deg) =", [round(a, 1) for a in ang_anchor])
    print("   ||[P_A,P_B]||=", round(c_anchor, 4), " (v5.1.0 datum: 0.500; 7.6/44.7/67.3)")

    print("\nD2 four curves -- monotonicity / continuity:")
    def dirn(v):
        v = np.array([x for x in v if np.isfinite(x)])
        d = np.diff(v)
        if np.all(d > -1e-9): return "monotone up"
        if np.all(d < 1e-9): return "monotone down"
        return "rises then falls (single peak)" if np.argmax(v) not in (0, len(v)-1) else "non-monotone"
    print("   (a) gamma(R)          :", dirn(gamma), "  [pi ->", round(np.degrees(gamma[-1]),1), "deg]")
    print("   (b) ||[P_A,P_B]|| H2+ :", dirn(comm_h2), " max", round(max(comm_h2),4))
    print("   (c) aggregate l-mix   :", dirn(agg))
    print("   (d) seed a(R)=R       :", dirn(a_seed), " (linear); e^aE1(a):", dirn(seed_val))
    print("   BONUS eta l-mixing    :", dirn(eta_part))

    print("\nD1 per-shell front  R*(n) where S_n = 1/sqrt2 (theta=45, peak char):")
    print("   n   R*(n)    r_n=n^2/Z   R*/r_n    ell_n=n/Z   R*/ell_n")
    for i, n in enumerate(n_shells):
        print(f"   {n}  {Rs[i]:7.3f}   {r_n[i]:7.1f}    {ratio[i]:6.3f}    {ell_n[i]:7.1f}   {ratio_ell[i]:6.3f}")
    print(f"\n   fit R*(n) ~ n^p :  p_all = {p_all:.3f} (R2={r2_all:.4f}) ; "
          f"p_(n>=2) = {p_hi:.3f} (R2={r2_hi:.4f})")
    print(f"   => n*(R) ~ R^{1.0/p_hi:.3f}   (predicted sqrt(ZR): 0.5 ; decay-length n/Z: 1.0)")
    print(f"   predicted R*(n)~n^2 (p=2, mean radius) vs measured p={p_hi:.2f}")

    os.makedirs(os.path.join(REPO, "debug", "data"), exist_ok=True)
    dst = os.path.join(REPO, "debug", "data", "decompactification_R_sweep.json")
    with open(dst, "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote", dst)
    return out


if __name__ == "__main__":
    run()
