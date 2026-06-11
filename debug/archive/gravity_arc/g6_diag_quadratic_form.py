"""Sprint G6-Diag first-pass — linearized perturbation diagnostic of the
spectral action around the CH Dirac on truncated S^3 substrate.

The question
------------
Does the GeoVac discrete substrate host eigenmodes of the second-variation
quadratic form S^(2)[V, V] of the spectral action that carry spin-2 angular
momentum content (the (j_L, j_R) = (1, 1) irrep of SO(4) = SU(2)_L x SU(2)_R)?

If YES: gravitons are POSSIBLE on the substrate (necessary condition for
graviton dynamics).

If NO: gravitons are structurally blocked at the substrate level.

Setup
-----
Background: CH Dirac D_0 on truncated H_{n_max}, eigenvalues +/-(n+3/2)
with multiplicity g_n^Weyl = (n+1)(n+2) per chirality (so g_n^Dirac =
2(n+1)(n+2)).

At n_max = 1: dim H = 4 + 12 = 16. Hilbert space decomposes into 4
"sectors":
  S_1 = positive chirality n=0:  dim 2, lambda = +3/2,  SO(4) irrep (1/2, 0)
  S_2 = negative chirality n=0:  dim 2, lambda = -3/2,  SO(4) irrep (0, 1/2)
  S_3 = positive chirality n=1:  dim 6, lambda = +5/2,  SO(4) irrep (1, 1/2)
  S_4 = negative chirality n=1:  dim 6, lambda = -5/2,  SO(4) irrep (1/2, 1)

Perturb D = D_0 + eps V where V is Hermitian. Spectral action
S = Tr e^{-D^2/Lambda^2} expands as S = S_0 + eps S^(1) + eps^2/2 S^(2) + ...

S^(2)[V, V] is a quadratic form on Herm(H) (256-dim real space at n_max=1).

In the D_0 eigenbasis, V decomposes into "block matrices" connecting
each pair of sectors. The quadratic form is DIAGONAL in this basis:
each off-block has a single eigenvalue.

SO(4) classification: off-block S_i x S_j carries the tensor product
of the SO(4) irreps of S_i and S_j. The (1, 1) irrep (graviton) appears
in specific off-block tensor products.

Spin-2 candidate off-blocks at n_max = 1 (worked out by SU(2) tensor
product decomposition):
  - S_3 x S_3 (within sector pos chir n=1): (1, 1/2) ⊗ (1, 1/2)
    decomposes into ... (1, 1) ... YES, one copy
  - S_4 x S_4 (within sector neg chir n=1): (1/2, 1) ⊗ (1/2, 1)
    decomposes into ... (1, 1) ... YES, one copy
  - S_1 x S_4 (off-block pos n=0 x neg n=1): (1/2, 0) ⊗ (1/2, 1)
    decomposes into (0, 1) + (1, 1)  YES, one copy
  - S_2 x S_3 (off-block neg n=0 x pos n=1): (0, 1/2) ⊗ (1, 1/2)
    decomposes into (1, 0) + (1, 1)  YES, one copy

So four candidate off-blocks contain (1, 1) graviton content.
Are the corresponding eigenvalues of S^(2) nonzero (= propagating)?

Implementation
--------------
1. Build D_0 eigenvalue list (16 entries)
2. For each PAIR of sectors (S_i, S_j) including i=j, compute the
   eigenvalue B_{i,j} of S^(2) on the Hermitian basis modes in that
   off-block.
3. Verify the analytical formulas via finite-difference computation
   of g(eps) = Tr exp(-(D_0 + eps V)^2 / Lambda^2) for several test V.
4. Tabulate eigenvalues with multiplicities.
5. For each spin-2 candidate off-block, report the multiplicity of the
   (1, 1) irrep within it and the eigenvalue.

Verdict
-------
- POSITIVE (graviton possible): at least one (1, 1) irrep has nonzero
  eigenvalue.
- NEGATIVE (graviton blocked): all (1, 1) irreps have zero eigenvalue
  OR they don't exist in the substrate.

Honest scope of THIS FIRST-PASS
-------------------------------
- Only n_max = 1 (smallest case).
- Only checks NECESSARY condition (mode existence + nonzero eigenvalue).
- Does NOT check SUFFICIENT condition (Fierz-Pauli kinetic structure,
  gauge invariance, propagator structure).
- Does NOT verify (1, 1) irrep multiplicities via explicit CG decomp;
  uses SU(2) tensor product algebra only.

Full sufficient analysis is the 2-3 week G6-Diag-Full sprint.
"""

import json
from pathlib import Path

import numpy as np

mp_dps = 50

OUT_JSON = Path(__file__).parent / "data" / "g6_diag_quadratic_form.json"
OUT_JSON.parent.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# CH Dirac on truncated S^3
# ---------------------------------------------------------------------------

def build_D0_eigenvalues(n_max):
    """Build the eigenvalue list of CH Dirac D_0 at truncation n_max on unit S^3.

    Returns: list of (lambda, sector_index) tuples. The sectors are:
        sector 1: (n=0, pos chir),  lambda = +3/2,  dim 2
        sector 2: (n=0, neg chir),  lambda = -3/2,  dim 2
        sector 3: (n=1, pos chir),  lambda = +5/2,  dim 6
        sector 4: (n=1, neg chir),  lambda = -5/2,  dim 6
    """
    eigs = []
    sectors = []
    sector_info = []  # (sector_idx, lambda, jL, jR, dim)
    sector_idx = 0
    for n in range(n_max + 1):
        lam_pos = (n + 1.5)
        lam_neg = -(n + 1.5)
        weyl_dim = (n + 1) * (n + 2)
        # positive chirality: (j_L, j_R) = ((n+1)/2, n/2)
        jL_pos = (n + 1) / 2.0
        jR_pos = n / 2.0
        sector_idx += 1
        sector_info.append((sector_idx, lam_pos, jL_pos, jR_pos, weyl_dim))
        for _ in range(weyl_dim):
            eigs.append(lam_pos)
            sectors.append(sector_idx)
        # negative chirality: (j_L, j_R) = (n/2, (n+1)/2)
        jL_neg = n / 2.0
        jR_neg = (n + 1) / 2.0
        sector_idx += 1
        sector_info.append((sector_idx, lam_neg, jL_neg, jR_neg, weyl_dim))
        for _ in range(weyl_dim):
            eigs.append(lam_neg)
            sectors.append(sector_idx)
    return np.array(eigs), np.array(sectors), sector_info


# ---------------------------------------------------------------------------
# Analytical quadratic form coefficients
# ---------------------------------------------------------------------------

def coeff_A_diagonal_or_within_sector(lam, Lambda_sq):
    """Coefficient of S^(2) on diagonal modes V_{alpha,alpha} or within-sector
    off-diagonal modes (where the perturbation is within a degenerate
    eigenvalue subspace).

    For a Hermitian basis matrix E in the D_0 eigenbasis with support entirely
    within a degenerate eigenvalue subspace at eigenvalue lambda:
      eigenvalue of S^(2) on E  =  a_lambda * (4 lambda^2 / Lambda^4 - 2 / Lambda^2)

    where a_lambda = exp(-lambda^2 / Lambda_sq).
    """
    a_lam = np.exp(-lam * lam / Lambda_sq)
    return a_lam * (4 * lam * lam / (Lambda_sq * Lambda_sq) - 2 / Lambda_sq)


def coeff_B_cross_sector(lam_a, lam_b, Lambda_sq):
    """Coefficient of S^(2) on cross-sector off-diagonal modes E^R_{alpha,beta}
    or E^I_{alpha,beta} where alpha in sector A (eigenvalue lam_a) and beta in
    sector B (eigenvalue lam_b), with lam_a != lam_b.

    For a Hermitian basis matrix with the only nonzero entries being (alpha, beta)
    and (beta, alpha) (matched for Hermiticity), the eigenvalue of S^(2) is
    derived from second-order perturbation theory:

      eigenvalue of S^(2)  =  -(2 / Lambda^2) * (lam_a * a_lam_a - lam_b * a_lam_b)
                              / (lam_a - lam_b)

    where a_lam = exp(-lam^2 / Lambda_sq).
    """
    a_a = np.exp(-lam_a * lam_a / Lambda_sq)
    a_b = np.exp(-lam_b * lam_b / Lambda_sq)
    return -(2 / Lambda_sq) * (lam_a * a_a - lam_b * a_b) / (lam_a - lam_b)


# ---------------------------------------------------------------------------
# Numerical verification via finite differences
# ---------------------------------------------------------------------------

def spectral_action(D, Lambda_sq):
    """S[D] = Tr exp(-D^2 / Lambda^2)."""
    D_sq = D @ D
    eigs = np.linalg.eigvalsh(D_sq)
    return np.sum(np.exp(-eigs / Lambda_sq))


def compute_g_double_prime_FD(D_0, V, Lambda_sq, eps=1e-3):
    """Finite-difference computation of g''(0) where g(eps) = Tr e^{-(D_0+eps V)^2/Lambda^2}.

    Uses 5-point stencil: g''(0) = (-g(-2eps) + 16 g(-eps) - 30 g(0) + 16 g(eps) - g(2eps)) / (12 eps^2)
    """
    g_m2 = spectral_action(D_0 - 2 * eps * V, Lambda_sq)
    g_m1 = spectral_action(D_0 - eps * V, Lambda_sq)
    g_0 = spectral_action(D_0, Lambda_sq)
    g_p1 = spectral_action(D_0 + eps * V, Lambda_sq)
    g_p2 = spectral_action(D_0 + 2 * eps * V, Lambda_sq)
    return (-g_m2 + 16 * g_m1 - 30 * g_0 + 16 * g_p1 - g_p2) / (12 * eps * eps)


# ---------------------------------------------------------------------------
# Main G6-Diag computation
# ---------------------------------------------------------------------------

def main():
    results = {}
    print("=" * 72)
    print("Sprint G6-Diag first-pass:")
    print("Second-variation quadratic form on truncated CH Dirac S^3 substrate")
    print("=" * 72)

    n_max = 1
    Lambda_sq = 6.0  # Lambda^2; chosen to be O(lowest mode^2 / 1) for moderate cutoff

    eigs, sector_assignments, sector_info = build_D0_eigenvalues(n_max)
    dim_H = len(eigs)
    D_0 = np.diag(eigs)  # diagonal in eigenbasis by construction

    print(f"\n[Setup] n_max = {n_max},  Lambda^2 = {Lambda_sq}")
    print(f"  dim H = {dim_H}")
    print(f"  Number of sectors: {len(sector_info)}")
    print(f"  Sector | lambda | (j_L, j_R) | dim")
    for sec_idx, lam, jL, jR, d in sector_info:
        print(f"  {sec_idx:>6d} | {lam:>+6.2f} | ({jL:.1f}, {jR:.1f})  | {d}")

    results["setup"] = {
        "n_max": n_max,
        "Lambda_sq": Lambda_sq,
        "dim_H": dim_H,
        "sectors": [{"idx": s, "lambda": l, "jL": jL, "jR": jR, "dim": d}
                    for s, l, jL, jR, d in sector_info],
    }

    # -----------------------------------------------------------------------
    # Step 1: Compute analytical eigenvalues of S^(2) per sector pair.
    # -----------------------------------------------------------------------
    print("\n[Step 1] Analytical eigenvalues of S^(2) per sector-pair off-block:")
    print(f"  {'sector pair':>13s} | {'lambda_a':>9s} | {'lambda_b':>9s} | {'eigenvalue':>11s}")
    eigenvalue_by_block = {}
    for i, (si, lam_i, jLi, jRi, di) in enumerate(sector_info):
        for j, (sj, lam_j, jLj, jRj, dj) in enumerate(sector_info):
            if j < i:
                continue
            if lam_i == lam_j and i == j:
                # within same sector (degenerate)
                eigenvalue = coeff_A_diagonal_or_within_sector(lam_i, Lambda_sq)
                label = f"S_{si} x S_{si}"
            elif lam_i != lam_j:
                eigenvalue = coeff_B_cross_sector(lam_i, lam_j, Lambda_sq)
                label = f"S_{si} x S_{sj}"
            else:
                # different sectors but same |lambda|? shouldn't happen with our setup
                eigenvalue = float('nan')
                label = f"S_{si} x S_{sj}"
            eigenvalue_by_block[(si, sj)] = eigenvalue
            print(f"  {label:>13s} | {lam_i:>+9.4f} | {lam_j:>+9.4f} | {eigenvalue:>+11.6f}")

    results["analytical_eigenvalues"] = {
        f"({si}, {sj})": ev for (si, sj), ev in eigenvalue_by_block.items()
    }

    # -----------------------------------------------------------------------
    # Step 2: Verify analytical eigenvalues via finite difference.
    # -----------------------------------------------------------------------
    print("\n[Step 2] Verify analytical eigenvalues via finite difference (5-point stencil):")
    print(f"  {'mode type':>30s} | {'analytical':>11s} | {'finite diff':>11s} | {'rel diff':>10s}")
    fd_panel = []
    # Within-sector modes: pick first state in each sector for diagonal,
    # and first pair within sector if dim >= 2
    for sec_idx, lam, jL, jR, d in sector_info:
        # Find first state index in this sector
        states_in_sec = np.where(sector_assignments == sec_idx)[0]
        if len(states_in_sec) == 0:
            continue
        # Diagonal mode: V = |alpha><alpha|, alpha = first state
        alpha = states_in_sec[0]
        V_diag = np.zeros((dim_H, dim_H), dtype=float)
        V_diag[alpha, alpha] = 1.0
        eig_analytical = coeff_A_diagonal_or_within_sector(lam, Lambda_sq)
        eig_fd = compute_g_double_prime_FD(D_0, V_diag, Lambda_sq, eps=1e-3)
        rel = abs(eig_analytical - eig_fd) / (abs(eig_analytical) + 1e-15)
        label = f"S_{sec_idx} diag (alpha={alpha})"
        print(f"  {label:>30s} | {eig_analytical:>+11.6f} | {eig_fd:>+11.6f} | {rel:.3e}")
        fd_panel.append({
            "mode": label, "analytical": eig_analytical, "fd": eig_fd, "rel_diff": rel,
        })
        if len(states_in_sec) >= 2:
            # Within-sector off-diagonal
            beta = states_in_sec[1]
            V_R = np.zeros((dim_H, dim_H), dtype=float)
            V_R[alpha, beta] = 1.0 / np.sqrt(2)
            V_R[beta, alpha] = 1.0 / np.sqrt(2)
            eig_fd = compute_g_double_prime_FD(D_0, V_R, Lambda_sq, eps=1e-3)
            rel = abs(eig_analytical - eig_fd) / (abs(eig_analytical) + 1e-15)
            label = f"S_{sec_idx} off-diag R (alpha={alpha}, beta={beta})"
            print(f"  {label:>30s} | {eig_analytical:>+11.6f} | {eig_fd:>+11.6f} | {rel:.3e}")
            fd_panel.append({
                "mode": label, "analytical": eig_analytical, "fd": eig_fd, "rel_diff": rel,
            })

    # Cross-sector modes
    for i_idx, (si, lam_i, jLi, jRi, di) in enumerate(sector_info):
        for j_idx, (sj, lam_j, jLj, jRj, dj) in enumerate(sector_info):
            if j_idx <= i_idx:
                continue
            if lam_i == lam_j:
                continue  # degenerate, handled above
            alpha = np.where(sector_assignments == si)[0][0]
            beta = np.where(sector_assignments == sj)[0][0]
            V_R = np.zeros((dim_H, dim_H), dtype=float)
            V_R[alpha, beta] = 1.0 / np.sqrt(2)
            V_R[beta, alpha] = 1.0 / np.sqrt(2)
            eig_analytical = coeff_B_cross_sector(lam_i, lam_j, Lambda_sq)
            eig_fd = compute_g_double_prime_FD(D_0, V_R, Lambda_sq, eps=1e-3)
            rel = abs(eig_analytical - eig_fd) / (abs(eig_analytical) + 1e-15)
            label = f"S_{si}xS_{sj} cross"
            print(f"  {label:>30s} | {eig_analytical:>+11.6f} | {eig_fd:>+11.6f} | {rel:.3e}")
            fd_panel.append({
                "mode": label, "analytical": eig_analytical, "fd": eig_fd, "rel_diff": rel,
            })

    results["fd_verification"] = fd_panel

    # -----------------------------------------------------------------------
    # Step 3: SO(4) decomposition: which off-blocks carry (1, 1) graviton irrep?
    # -----------------------------------------------------------------------
    print("\n[Step 3] SO(4) = SU(2)_L x SU(2)_R tensor product decomposition:")
    print("  For each off-block S_i x S_j*, decompose into SU(2)_L x SU(2)_R irreps.")
    print("  Spin-2 graviton irrep is (1, 1) with dim 9.")
    print()

    def su2_tensor_product(j1, j2):
        """SU(2) tensor product j1 ⊗ j2 = |j1-j2| ⊕ ... ⊕ (j1+j2)."""
        result = []
        j_min = abs(j1 - j2)
        j_max = j1 + j2
        j = j_min
        while j <= j_max + 1e-9:
            result.append(j)
            j += 1.0
        return result

    print(f"  {'sector pair':>13s} | {'SU(2)_L decomp':>30s} | {'SU(2)_R decomp':>30s} | {'(1,1) mult':>10s}")
    graviton_blocks = []
    for i_idx, (si, lam_i, jLi, jRi, di) in enumerate(sector_info):
        for j_idx, (sj, lam_j, jLj, jRj, dj) in enumerate(sector_info):
            if j_idx < i_idx:
                continue
            # Compute (jLi, jRi) ⊗ (jLj, jRj)* = (jLi, jRi) ⊗ (jLj, jRj)
            # (since SU(2) irreps are self-conjugate)
            JL_decomp = su2_tensor_product(jLi, jLj)
            JR_decomp = su2_tensor_product(jRi, jRj)
            # (1, 1) appears iff 1 is in JL_decomp AND 1 is in JR_decomp
            has_1_L = any(abs(j - 1.0) < 1e-9 for j in JL_decomp)
            has_1_R = any(abs(j - 1.0) < 1e-9 for j in JR_decomp)
            grav_mult = 1 if (has_1_L and has_1_R) else 0
            label = f"S_{si} x S_{sj}"
            JL_str = ", ".join(f"{j:g}" for j in JL_decomp)
            JR_str = ", ".join(f"{j:g}" for j in JR_decomp)
            mark = "<-- GRAVITON CANDIDATE" if grav_mult > 0 else ""
            print(f"  {label:>13s} | {JL_str:>30s} | {JR_str:>30s} | {grav_mult:>10d}  {mark}")
            if grav_mult > 0:
                graviton_blocks.append({
                    "sectors": (si, sj),
                    "JL_decomp": JL_decomp,
                    "JR_decomp": JR_decomp,
                    "graviton_mult": grav_mult,
                    "eigenvalue": eigenvalue_by_block.get((si, sj), eigenvalue_by_block.get((sj, si))),
                })

    results["graviton_candidate_blocks"] = [
        {"sectors": [b["sectors"][0], b["sectors"][1]],
         "JL_decomp": b["JL_decomp"],
         "JR_decomp": b["JR_decomp"],
         "graviton_mult": b["graviton_mult"],
         "eigenvalue": b["eigenvalue"]}
        for b in graviton_blocks
    ]

    # -----------------------------------------------------------------------
    # Step 4: Verdict.
    # -----------------------------------------------------------------------
    print("\n[Step 4] Verdict:")
    print()
    print(f"  Graviton candidate off-blocks: {len(graviton_blocks)} found.")
    if not graviton_blocks:
        print("  NEGATIVE: no (1, 1) irrep content in any off-block at n_max = 1.")
        verdict = "NEGATIVE"
    else:
        all_nonzero = all(abs(b["eigenvalue"]) > 1e-9 for b in graviton_blocks)
        if all_nonzero:
            print("  POSITIVE (first-pass): all (1, 1) candidate blocks have nonzero")
            print("  S^(2) eigenvalue -> graviton-irrep modes propagate.")
            verdict = "POSITIVE-FIRST-PASS"
        else:
            zeros = [b for b in graviton_blocks if abs(b["eigenvalue"]) <= 1e-9]
            print(f"  MIXED: {len(zeros)} of {len(graviton_blocks)} (1, 1) candidate blocks")
            print(f"  have zero eigenvalue.")
            verdict = "MIXED"

    results["verdict"] = verdict
    print()
    print(f"  Verdict: {verdict}")
    print()
    print("  Honest scope: this first-pass only verifies NECESSARY condition")
    print("  (existence of nonzero (1, 1)-irrep eigenmodes). SUFFICIENT")
    print("  conditions (Fierz-Pauli kinetic, gauge invariance, propagator")
    print("  structure) require the full 2-3 week G6-Diag-Full sprint.")

    # Save JSON
    with OUT_JSON.open("w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
    print("=" * 72)
    return results


if __name__ == "__main__":
    main()
