"""
TC three-body operator L3: does GeoVac's Gaunt/6j angular selection rule
collapse it to two-body (as momentum conservation does in plane waves)?

Decomposer-style structural probe.  Angular-only (single center, Coulomb-Sturmian
S^3 / spherical-harmonic labels).  Uses the framework's own wigner3j.

L3 = -1/2 sum_i sum_{j!=i,k!=i,j!=k} grad_i u(r_ij) . grad_i u(r_ik)

Correlator multipole expansion (same shape as the potential, Paper 22 eq. Legendre):
    u(|r_i-r_j|) = sum_{L,M} (4 pi/(2L+1)) u_L(r_i,r_j) Y*_{LM}(rhat_j) Y_{LM}(rhat_i)

=> at the SHARED vertex i, two correlator harmonics Y_{LM}(rhat_i) and
   Y_{L'M'}(rhat_i) multiply, together with the external bra/ket harmonics of
   electron i.  Particles j and k each see ONE correlator harmonic.

The angular signature of "two correlator lines meet at i" is the four-harmonic
integral
    W(a;L,L';c) = int dOmega_i  Y*_{l_a m_a} Y_{LM} Y_{L'M'} Y_{l_c m_c}
which decomposes exactly through the SO(3) tensor product
    Y_{LM} Y_{L'M'} = sum_Lambda G_pair(L,L',Lambda) Y_{Lambda, M+M'}
as
    W(a;L,L';c) = sum_Lambda  G_pair(L,L',Lambda) * G_ext(a,Lambda,c).

Lambda is the SO(3) analog of the plane-wave resultant q+q'.  In plane waves
q+q' is a SINGLE vector (multiplicity 1) -> the shared vertex separates ->
two-body only.  Here Lambda runs over |L-L'|..L+L' intersected with the
external triangle |l_a-l_c|..l_a+l_c.  If >1 Lambda survives, the vertex does
NOT separate: irreducibly three-body.

We MEASURE:
  (1) the SO(3) channel multiplicity n_Lambda per external (a,c) pair,
  (2) the matrix rank of W as a bilinear form on (correlator L,M) x (L',M'),
  (3) full three-body angular tensor sparsity (m-conservation at 3 vertices),
  (4) an explicit factorization test: does the (i-line):(jk-line) reshape have
      Schmidt rank 1 (collapse) or > 1 (genuine 3-body)?
"""

import json
import numpy as np
from math import sqrt, pi
from itertools import product
from geovac.angular_integrals import wigner3j

FOURPI = 4.0 * pi


def gaunt(l1, m1, l2, m2, l3, m3):
    """int Y_{l1m1} Y_{l2m2} Y_{l3m3} dOmega  (complex spherical harmonics)."""
    if (m1 + m2 + m3) != 0:
        return 0.0
    pref = sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / FOURPI)
    return pref * wigner3j(l1, l2, l3, 0, 0, 0) * wigner3j(l1, l2, l3, m1, m2, m3)


def g_ext(la, ma, Lam, MLam, lc, mc):
    """int Y*_{la ma} Y_{Lam MLam} Y_{lc mc} dOmega."""
    # Y*_{la ma} = (-1)^ma Y_{la,-ma}
    return (-1) ** ma * gaunt(la, -ma, Lam, MLam, lc, mc)


def g_pair(L, M, Lp, Mp, Lam):
    """Coefficient of Y_{Lam, M+Mp} in the product Y_{LM} Y_{L'M'}.

    Y_{LM} Y_{L'M'} = sum_Lam c(L,L',Lam) Y_{Lam,M+Mp}
    c = int Y_{LM} Y_{L'M'} Y*_{Lam,M+Mp} dOmega
      = (-1)^(M+Mp) int Y_{LM} Y_{L'M'} Y_{Lam,-(M+Mp)} dOmega
    """
    MLam = M + Mp
    if abs(MLam) > Lam:
        return 0.0
    return (-1) ** (MLam) * gaunt(L, M, Lp, Mp, Lam, -MLam)


def four_Y(la, ma, L, M, Lp, Mp, lc, mc, keep_lowest_only=False):
    """Shared-vertex four-harmonic integral, exact, plus the list of active Lambda.

    keep_lowest_only=True keeps only the LOWEST allowed SO(3) resultant channel
    Lambda -- the abelian 'single resultant q+q'' plane-wave mimic.
    """
    MLam = M + Mp
    if MLam != (ma - mc):
        # m-conservation at vertex i: -ma + M + Mp + mc = 0
        return 0.0, []
    val = 0.0
    active = []
    Lam_lo = max(abs(L - Lp), abs(la - lc), abs(MLam))
    Lam_hi = min(L + Lp, la + lc)
    for Lam in range(Lam_lo, Lam_hi + 1):
        gp = g_pair(L, M, Lp, Mp, Lam)
        if gp == 0.0:
            continue
        ge = g_ext(la, ma, Lam, MLam, lc, mc)
        if ge == 0.0:
            continue
        if keep_lowest_only and active:
            continue  # discard all but the first (lowest) active channel
        val += gp * ge
        active.append(Lam)
    return val, active


def orbitals(l_max):
    return [(l, m) for l in range(l_max + 1) for m in range(-l, l + 1)]


def cross_check_quadrature(cases):
    """Validate four_Y against direct numerical spherical quadrature."""
    from scipy.special import sph_harm
    # Gauss-Legendre in cos(theta), uniform phi
    nth, nph = 40, 41
    x, w = np.polynomial.legendre.leggauss(nth)   # nodes in cos(theta)
    theta = np.arccos(x)
    phi = np.linspace(0, 2 * pi, nph, endpoint=False)
    dphi = 2 * pi / nph
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    WT = np.outer(w, np.ones(nph)) * dphi  # measure d(cos th) dphi

    def Y(l, m):
        # scipy sph_harm(m, l, phi, theta)
        return sph_harm(m, l, PH, TH)

    out = []
    for (la, ma, L, M, Lp, Mp, lc, mc) in cases:
        integrand = np.conj(Y(la, ma)) * Y(L, M) * Y(Lp, Mp) * Y(lc, mc)
        num = np.sum(integrand * WT)
        exact, _ = four_Y(la, ma, L, M, Lp, Mp, lc, mc)
        out.append((float(np.real(num)), float(np.imag(num)), exact))
    return out


def shared_vertex_analysis(l_max_ext, L_max_corr, tol=1e-12):
    """For every external (a,c), build W[(L,M),(L',M')] and measure rank + n_Lambda."""
    ext = orbitals(l_max_ext)
    corr = orbitals(L_max_corr)
    results = []
    nLambda_hist = {}
    for (la, ma) in ext:
        for (lc, mc) in ext:
            rows = []  # each row = flattened over (L',M') for fixed (L,M)
            maxLambda = 0
            any_nonzero = False
            A = np.zeros((len(corr), len(corr)), dtype=float)
            for iL, (L, M) in enumerate(corr):
                for iLp, (Lp, Mp) in enumerate(corr):
                    v, active = four_Y(la, ma, L, M, Lp, Mp, lc, mc)
                    A[iL, iLp] = v
                    if abs(v) > tol:
                        any_nonzero = True
                    if active:
                        maxLambda = max(maxLambda, len(active))
                        nLambda_hist[len(active)] = nLambda_hist.get(len(active), 0) + 1
            if not any_nonzero:
                continue
            rank = int(np.linalg.matrix_rank(A, tol=tol))
            results.append(dict(la=la, ma=ma, lc=lc, mc=mc,
                                rank=rank, max_nLambda=maxLambda,
                                nnz=int(np.sum(np.abs(A) > tol))))
    return results, nLambda_hist


def single_resultant_control(l_max, L_max_corr, tol=1e-12):
    """PW-mimic: keep ONLY the lowest SO(3) resultant channel at the shared vertex
    (the abelian 'single q+q'' analog) and measure how much of the full L3 angular
    tensor is captured vs lost to higher channels."""
    ext = orbitals(l_max)
    corr = orbitals(L_max_corr)

    def G_leg(la, ma, L, M, lc, mc):
        return (-1) ** (ma + M) * gaunt(la, -ma, L, -M, lc, mc)

    ext_pairs = [(a, c) for a in ext for c in ext]
    full = []
    trunc = []
    for (ai, ci) in ext_pairs:
        la, ma = ai; lci, mci = ci
        for (aj, cj) in ext_pairs:
            laj, maj = aj; lcj, mcj = cj
            for (ak, ck) in ext_pairs:
                lak, mak = ak; lck, mck = ck
                vf = 0.0; vt = 0.0
                for (L, M) in corr:
                    gj = G_leg(laj, maj, L, M, lcj, mcj)
                    if gj == 0.0:
                        continue
                    for (Lp, Mp) in corr:
                        gk = G_leg(lak, mak, Lp, Mp, lck, mck)
                        if gk == 0.0:
                            continue
                        wf, _ = four_Y(la, ma, L, M, Lp, Mp, lci, mci)
                        wt, _ = four_Y(la, ma, L, M, Lp, Mp, lci, mci,
                                       keep_lowest_only=True)
                        vf += wf * gj * gk
                        vt += wt * gj * gk
                full.append(vf); trunc.append(vt)
    full = np.array(full); trunc = np.array(trunc)
    fro_full = float(np.linalg.norm(full))
    fro_resid = float(np.linalg.norm(full - trunc))
    nnz_full = int(np.sum(np.abs(full) > tol))
    nnz_missed = int(np.sum((np.abs(full) > tol) & (np.abs(trunc) <= tol)))
    return dict(fro_full=fro_full,
                fro_residual_after_single_channel=fro_resid,
                frac_norm_in_higher_channels=fro_resid / fro_full if fro_full else 0.0,
                nnz_full=nnz_full,
                nnz_entries_missed_by_single_channel=nnz_missed)


def full_threebody_sparsity(l_max, L_max_corr, tol=1e-12):
    """Enumerate the full L3 angular tensor and measure density + factorization rank.

    Indices: electron i (bra a_i, ket c_i), electron j (bra a_j, ket c_j),
             electron k (bra a_k, ket c_k).  Correlator lines L (i-j) and L' (i-k).
    Angular value = sum over correlator (L,M),(L',M') of
        W(a_i;L,M,L',M';c_i) * G_j(a_j; L,M; c_j) * G_k(a_k; L',M'; c_k)
    where G_j = int Y*_{a_j} Y*_{LM} Y_{c_j} (particle j sees conj of correlator
    harmonic on rhat_j).  m-conservation is automatic through the gaunts.
    """
    ext = orbitals(l_max)
    corr = orbitals(L_max_corr)

    def G_leg(la, ma, L, M, lc, mc):
        # particle j: int Y*_{la ma} Y*_{L M} Y_{lc mc} = (-1)^(ma+M) gaunt(la,-ma,L,-M,lc,mc)
        return (-1) ** (ma + M) * gaunt(la, -ma, L, -M, lc, mc)

    # Precompute the shared-vertex tensor W[a_i,c_i][(L,M),(L',M')]
    # and legs.  Then for the FACTORIZATION rank we reshape the full tensor as
    # (electron i block) x (electrons j,k block) and take the Schmidt rank.
    n_total = 0
    n_nonzero = 0
    # sample-based factorization rank: build a matrix M[ (a_i,c_i) , (a_j,c_j,a_k,c_k) ]
    ext_pairs = [(a, c) for a in ext for c in ext]
    idx_i = {p: n for n, p in enumerate(ext_pairs)}
    # jk space is large; restrict to l_max for a rank probe
    big = np.zeros((len(ext_pairs), len(ext_pairs) * len(ext_pairs)))
    for (ai, ci) in ext_pairs:
        la, ma = ai
        lci, mci = ci
        for j_idx, (aj, cj) in enumerate(ext_pairs):
            laj, maj = aj
            lcj, mcj = cj
            for k_idx, (ak, ck) in enumerate(ext_pairs):
                lak, mak = ak
                lck, mck = ck
                val = 0.0
                for (L, M) in corr:
                    gj = G_leg(laj, maj, L, M, lcj, mcj)
                    if gj == 0.0:
                        continue
                    for (Lp, Mp) in corr:
                        gk = G_leg(lak, mak, Lp, Mp, lck, mck)
                        if gk == 0.0:
                            continue
                        w, _ = four_Y(la, ma, L, M, Lp, Mp, lci, mci)
                        if w == 0.0:
                            continue
                        val += w * gj * gk
                n_total += 1
                if abs(val) > tol:
                    n_nonzero += 1
                col = j_idx * len(ext_pairs) + k_idx
                big[idx_i[(ai, ci)], col] = val
    density = n_nonzero / n_total if n_total else 0.0
    schmidt_rank = int(np.linalg.matrix_rank(big, tol=tol))
    return dict(n_total=n_total, n_nonzero=n_nonzero, density=density,
                schmidt_rank_i_vs_jk=schmidt_rank,
                dim_i=len(ext_pairs))


def main():
    report = {}

    # --- validation gate: four_Y vs numerical quadrature -----------------
    cases = [
        (1, 0, 1, 0, 1, 0, 1, 0),
        (2, 1, 1, 0, 1, 1, 2, 0),
        (1, 1, 2, -1, 1, 0, 2, 0),
        (2, -2, 2, 1, 2, 1, 2, 0),
        (0, 0, 1, 0, 1, 0, 0, 0),
        (1, -1, 2, 2, 1, -1, 2, 0),
    ]
    qc = cross_check_quadrature(cases)
    val_err = max(abs(q[0] - q[2]) for q in qc)
    val_imag = max(abs(q[1]) for q in qc)
    report["validation"] = dict(cases=cases,
                                quad_vs_exact=qc,
                                max_abs_err=val_err,
                                max_abs_imag=val_imag)
    print(f"[validation] max|quad-exact| = {val_err:.2e}, max|imag| = {val_imag:.2e}")

    # --- shared-vertex rank / SO(3) channel multiplicity -----------------
    for (le, lc) in [(1, 1), (2, 2), (2, 3)]:
        res, hist = shared_vertex_analysis(le, lc)
        ranks = [r["rank"] for r in res]
        maxnl = [r["max_nLambda"] for r in res]
        frac_rank_ge2 = np.mean([r >= 2 for r in ranks]) if ranks else 0.0
        frac_nl_ge2 = np.mean([m >= 2 for m in maxnl]) if maxnl else 0.0
        key = f"shared_vertex_lext{le}_lcorr{lc}"
        report[key] = dict(
            n_external_pairs_nonzero=len(res),
            rank_hist={int(k): int((np.array(ranks) == k).sum()) for k in sorted(set(ranks))},
            max_nLambda_hist={int(k): int((np.array(maxnl) == k).sum()) for k in sorted(set(maxnl))},
            nLambda_perW_hist={int(k): int(v) for k, v in sorted(hist.items())},
            frac_external_rank_ge2=float(frac_rank_ge2),
            frac_external_nLambda_ge2=float(frac_nl_ge2),
            max_rank=int(max(ranks)) if ranks else 0,
        )
        print(f"[{key}] nonzero ext pairs={len(res)}, max rank={max(ranks) if ranks else 0}, "
              f"frac rank>=2={frac_rank_ge2:.3f}, frac nLambda>=2={frac_nl_ge2:.3f}")

    # --- full three-body tensor density + factorization rank -------------
    for (l_max, L_corr) in [(1, 1), (1, 2), (2, 2)]:
        fb = full_threebody_sparsity(l_max, L_corr)
        key = f"full3body_lmax{l_max}_lcorr{L_corr}"
        report[key] = fb
        print(f"[{key}] density={fb['density']:.4f} "
              f"({fb['n_nonzero']}/{fb['n_total']}), "
              f"Schmidt rank i:jk = {fb['schmidt_rank_i_vs_jk']} (dim_i={fb['dim_i']})")

    # --- PW-mimic single-resultant-channel control -----------------------
    for (l_max, L_corr) in [(1, 2), (2, 2)]:
        sc = single_resultant_control(l_max, L_corr)
        key = f"single_resultant_control_lmax{l_max}_lcorr{L_corr}"
        report[key] = sc
        print(f"[{key}] frac norm in higher SO(3) channels = "
              f"{sc['frac_norm_in_higher_channels']:.3f}; "
              f"nonzero entries missed by single channel = "
              f"{sc['nnz_entries_missed_by_single_channel']}/{sc['nnz_full']}")

    with open("debug/data/tc_threebody_collapse.json", "w") as f:
        json.dump(report, f, indent=2)
    print("\nwrote debug/data/tc_threebody_collapse.json")


if __name__ == "__main__":
    main()
