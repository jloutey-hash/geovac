"""G6 Fierz-Pauli decomposition of the (1,1) graviton subspace.

Decomposes the (1,1) irrep of SO(4) under diagonal SU(2):
  (1,1) = J=0 (trace, 1 dim) + J=1 (antisymmetric, 3 dim) + J=2 (graviton, 5 dim)

Computes S^(2) analytically (no stencil) within each J-sector.
Compares J=2 eigenvalues to continuum Lichnerowicz TT spectrum k(k+2)-2.

The analytical S^(2) weight is:
  S^(2)[V] = sum_{k,l} w_{kl} |V_{kl}|^2
  w_{kk} = e^{-lam_k^2/L} (4 lam_k^2/L^2 - 2/L)       [within-sector]
  w_{kl} = (lam_k+lam_l)^2/L^2 * phi(lam_k^2/L, lam_l^2/L) - (a_k+a_l)/L  [cross-sector]
  phi(x,y) = (e^{-x} - e^{-y})/(y-x), phi(x,x) = e^{-x}
"""

import numpy as np
from sympy.physics.wigner import clebsch_gordan
from sympy import N as sympy_N, Rational
import json
import time
from pathlib import Path

# This module lives in tests/gravity_support/ (migrated from debug/ on
# 2026-07-03; see README.md in this directory). Standalone driver runs
# still write to the transient clean-room dir debug/data/; the mkdir is
# deferred to main() so that importing this module from the test suite
# has no filesystem side effects.
_REPO = Path(__file__).resolve().parents[2]
OUT_DIR = _REPO / "debug" / "data"

_cg_cache = {}


def cg_val(j1, m1, j2, m2, J, M):
    key = (j1, m1, j2, m2, J, M)
    if key in _cg_cache:
        return _cg_cache[key]
    if abs(m1 + m2 - M) > 1e-10:
        _cg_cache[key] = 0.0
        return 0.0
    if J < abs(j1 - j2) - 1e-10 or J > j1 + j2 + 1e-10:
        _cg_cache[key] = 0.0
        return 0.0
    if abs(m1) > j1 + 1e-10 or abs(m2) > j2 + 1e-10 or abs(M) > J + 1e-10:
        _cg_cache[key] = 0.0
        return 0.0
    val = float(sympy_N(clebsch_gordan(
        Rational(int(2 * j1), 2), Rational(int(2 * j2), 2), Rational(int(2 * J), 2),
        Rational(int(2 * m1), 2), Rational(int(2 * m2), 2), Rational(int(2 * M), 2))))
    _cg_cache[key] = val
    return val


def half_int_range(j):
    m = -j
    result = []
    while m <= j + 1e-10:
        result.append(round(2 * m) / 2)
        m += 1.0
    return result


class CHBasis:
    def __init__(self, n_max):
        self.n_max = n_max
        self.states = []
        self.sectors = {}
        for n in range(n_max + 1):
            for chi in [+1, -1]:
                j_L = (n + 1) / 2.0 if chi == +1 else n / 2.0
                j_R = n / 2.0 if chi == +1 else (n + 1) / 2.0
                key = (n, chi, j_L, j_R)
                indices = []
                for m_L in half_int_range(j_L):
                    for m_R in half_int_range(j_R):
                        indices.append(len(self.states))
                        self.states.append((n, chi, j_L, j_R, m_L, m_R))
                self.sectors[key] = indices
        self.dim = len(self.states)
        self.D0 = np.array([chi * (n + 1.5) for (n, chi, *_) in self.states])
        self.sector_keys = sorted(self.sectors.keys())


def build_fp_basis(basis):
    """Build (1,1) basis decomposed by diagonal-SU(2) J = 0, 1, 2.

    Returns: {J: list of orthonormal (dim x dim) Hermitian matrices}
    Also returns sector-pair labels for each basis vector.
    """
    dim = basis.dim
    states = basis.states
    JL, JR = 1.0, 1.0

    j_raw = {0: [], 1: [], 2: []}
    j_labels = {0: [], 1: [], 2: []}

    for a_idx, key_a in enumerate(basis.sector_keys):
        na, chia, jLa, jRa = key_a
        for b_idx, key_b in enumerate(basis.sector_keys):
            if b_idx < a_idx:
                continue
            nb, chib, jLb, jRb = key_b

            if JL < abs(jLa - jLb) - 1e-10 or JL > jLa + jLb + 1e-10:
                continue
            if JR < abs(jRa - jRb) - 1e-10 or JR > jRa + jRb + 1e-10:
                continue

            raw = {}
            for M_L in half_int_range(JL):
                for M_R in half_int_range(JR):
                    V = np.zeros((dim, dim))
                    for i_loc in basis.sectors[key_a]:
                        si = states[i_loc]
                        for j_loc in basis.sectors[key_b]:
                            sj = states[j_loc]
                            c_L = cg_val(sj[2], sj[4], JL, M_L, si[2], si[4])
                            c_R = cg_val(sj[3], sj[5], JR, M_R, si[3], si[5])
                            if abs(c_L * c_R) > 1e-15:
                                V[i_loc, j_loc] += c_L * c_R
                    if np.max(np.abs(V)) > 1e-14:
                        raw[(M_L, M_R)] = V

            if not raw:
                continue

            is_within = (a_idx == b_idx)
            lam_a = chia * (na + 1.5)
            lam_b = chib * (nb + 1.5)
            sector_type = "within" if is_within else "cross"

            for J in [0, 1, 2]:
                for M in half_int_range(float(J)):
                    VJM = np.zeros((dim, dim))
                    any_nonzero = False
                    for M_L in half_int_range(JL):
                        M_R_needed = M - M_L
                        M_R = round(2 * M_R_needed) / 2
                        if abs(M_R - M_R_needed) > 1e-10:
                            continue
                        if abs(M_R) > JR + 1e-10:
                            continue
                        if (M_L, M_R) not in raw:
                            continue
                        c = cg_val(JL, M_L, JR, M_R, float(J), M)
                        if abs(c) > 1e-15:
                            VJM += c * raw[(M_L, M_R)]
                            any_nonzero = True

                    if not any_nonzero or np.max(np.abs(VJM)) < 1e-14:
                        continue

                    lbl = f"({na},{chia})-({nb},{chib}) J={J} M={M} {sector_type}"

                    if not is_within:
                        Vh = (VJM + VJM.T) / 2
                        Va = (VJM - VJM.T) / 2
                        if np.max(np.abs(Vh)) > 1e-14:
                            j_raw[J].append(Vh)
                            j_labels[J].append(lbl + " sym")
                        if np.max(np.abs(Va)) > 1e-14:
                            j_raw[J].append(Va)
                            j_labels[J].append(lbl + " anti")
                    else:
                        Vh = (VJM + VJM.T) / 2
                        if np.max(np.abs(Vh)) > 1e-14:
                            j_raw[J].append(Vh)
                            j_labels[J].append(lbl)

    result = {}
    for J in [0, 1, 2]:
        flat = [v.flatten() for v in j_raw[J]]
        orth = []
        for v in flat:
            w = v.copy()
            for b in orth:
                w -= np.dot(w, b) * b
            n = np.linalg.norm(w)
            if n > 1e-10:
                orth.append(w / n)
        result[J] = [v.reshape((dim, dim)) for v in orth]

    return result


def s2_weight_matrix(D0, Lambda_sq):
    """Analytical S^(2) weight matrix W[i,j]."""
    dim = len(D0)
    lam = D0
    a = np.exp(-lam ** 2 / Lambda_sq)
    W = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            if abs(lam[i] ** 2 - lam[j] ** 2) < 1e-12:
                W[i, j] = a[i] * (4 * lam[i] ** 2 / Lambda_sq ** 2 - 2 / Lambda_sq)
            else:
                phi2 = Lambda_sq * (a[i] - a[j]) / (lam[j] ** 2 - lam[i] ** 2)
                W[i, j] = (lam[i] + lam[j]) ** 2 / Lambda_sq ** 2 * phi2 - (a[i] + a[j]) / Lambda_sq
    return W


def compute_s2_stencil(D0, V, Lambda_sq, eps=1e-4):
    """5-point stencil for validation."""
    D0_mat = np.diag(D0)

    def S(e):
        eigs = np.linalg.eigvalsh(D0_mat + e * V)
        return np.sum(np.exp(-eigs ** 2 / Lambda_sq))

    return (-S(2 * eps) + 16 * S(eps) - 30 * S(0) + 16 * S(-eps) - S(-2 * eps)) / (12 * eps ** 2)


def cluster_eigenvalues(eigs, tol=1e-4):
    if len(eigs) == 0:
        return []
    sorted_e = np.sort(eigs)[::-1]
    clusters = []
    i = 0
    while i < len(sorted_e):
        val = sorted_e[i]
        count = 1
        while i + count < len(sorted_e) and abs(sorted_e[i + count] - val) < max(tol, abs(val) * 1e-3):
            count += 1
        clusters.append((val, count))
        i += count
    return clusters


def run_analysis(n_max, Lambda_sq=6.0, validate=False):
    print(f"\n{'=' * 72}")
    print(f"n_max = {n_max}, Lambda^2 = {Lambda_sq}")
    print(f"{'=' * 72}")

    t0 = time.time()
    basis = CHBasis(n_max)
    print(f"  dim_H = {basis.dim}")

    t1 = time.time()
    fp = build_fp_basis(basis)
    t_basis = time.time() - t1

    total = sum(len(fp[J]) for J in [0, 1, 2])
    print(f"  FP basis built in {t_basis:.1f}s:")
    for J in [0, 1, 2]:
        name = {0: "trace/scalar", 1: "antisymmetric", 2: "sym-traceless (graviton)"}[J]
        print(f"    J={J} ({name}): dim = {len(fp[J])}")
    print(f"    Total (1,1) dim: {total}")

    W = s2_weight_matrix(basis.D0, Lambda_sq)

    if validate and len(fp[2]) > 0:
        print(f"\n  Validating analytical S^(2) vs 5-point stencil...")
        n_test = min(5, len(fp[2]))
        max_err = 0
        for idx in range(n_test):
            V = fp[2][idx]
            s2_a = np.sum(W * V * V)
            s2_n = compute_s2_stencil(basis.D0, V, Lambda_sq)
            if abs(s2_n) > 1e-15:
                err = abs(s2_a - s2_n) / abs(s2_n)
                max_err = max(max_err, err)
        print(f"    Max relative error ({n_test} J=2 vectors): {max_err:.2e}")

    results = {"n_max": n_max, "Lambda_sq": Lambda_sq, "dim_H": basis.dim}

    for J in [0, 1, 2]:
        n_b = len(fp[J])
        if n_b == 0:
            results[f"J={J}"] = {"dim": 0, "s2_clusters": [], "lap_clusters": []}
            continue

        V_stack = np.array(fp[J])
        K_s2 = np.einsum('aij,ij,bij->ab', V_stack, W, V_stack)
        eigs_s2 = np.linalg.eigvalsh(K_s2)
        cl_s2 = cluster_eigenvalues(eigs_s2)

        DiffMat = basis.D0[:, None] - basis.D0[None, :]
        DV = DiffMat[None, :, :] * V_stack
        K_lap = np.einsum('aij,bij->ab', DV, DV)
        eigs_lap = np.linalg.eigvalsh(K_lap)
        cl_lap = cluster_eigenvalues(eigs_lap, tol=0.5)

        name = {0: "J=0 (trace)", 1: "J=1 (antisym)", 2: "J=2 (graviton)"}[J]
        print(f"\n  {name}  [dim={n_b}]")
        print(f"    S^(2) eigenvalues:")
        for val, mult in cl_s2:
            tag = ""
            if abs(val) < 1e-8:
                tag = "  [zero]"
            print(f"      {val:+.8f}  x{mult}{tag}")

        print(f"    Discrete Laplacian ||[D,V]||^2:")
        for val, mult in cl_lap:
            if abs(val) > 0.5:
                k = np.sqrt(val)
                print(f"      {val:.1f} = ({k:.1f})^2  x{mult}")
            else:
                print(f"      ~0  x{mult}")

        # Within-sector vs cross-sector count
        n_within = 0
        n_cross = 0
        for v in fp[J]:
            lap_val = np.sum(DiffMat ** 2 * v ** 2)
            if lap_val < 0.5:
                n_within += 1
            else:
                n_cross += 1
        print(f"    Within-sector (||[D,V]||^2=0): {n_within}")
        print(f"    Cross-sector  (||[D,V]||^2>0): {n_cross}")

        results[f"J={J}"] = {
            "dim": n_b,
            "s2_eigenvalues": np.sort(eigs_s2)[::-1].tolist(),
            "s2_clusters": [(float(v), int(m)) for v, m in cl_s2],
            "lap_eigenvalues": np.sort(eigs_lap)[::-1].tolist(),
            "lap_clusters": [(float(v), int(m)) for v, m in cl_lap],
            "n_within_sector": n_within,
            "n_cross_sector": n_cross,
        }

    elapsed = time.time() - t0
    results["time_s"] = elapsed
    print(f"\n  Total: {elapsed:.1f}s")
    return results


def lichnerowicz_comparison(all_results):
    """Compare J=2 positive eigenvalues to Lichnerowicz TT spectrum."""
    print(f"\n\n{'=' * 72}")
    print("LICHNEROWICZ COMPARISON")
    print("Continuum TT eigenvalues on S^3: lambda_k = k(k+2)-2, k >= 2")
    print("  k=2: 6,  k=3: 13,  k=4: 22,  k=5: 33,  k=6: 46")
    print("  Ratios (vs k=2): 2.167, 3.667, 5.500, 7.667")
    print(f"{'=' * 72}")

    for key, res in sorted(all_results.items()):
        j2 = res.get("J=2", {})
        if not j2 or j2.get("dim", 0) == 0:
            continue
        pos_clusters = [(v, m) for v, m in j2["s2_clusters"] if v > 1e-6]
        if len(pos_clusters) < 2:
            continue
        min_pos = min(v for v, m in pos_clusters)
        print(f"\n  {key}:")
        print(f"    J=2 positive eigenvalue ratios (vs min={min_pos:.6f}):")
        for v, m in sorted(pos_clusters):
            ratio = v / min_pos
            print(f"      ratio={ratio:.4f}  val={v:.6f}  x{m}")


def main():
    print("=" * 72)
    print("G6 FIERZ-PAULI: (1,1) = J=0 + J=1 + J=2 under diagonal SU(2)")
    print("=" * 72)

    all_results = {}

    for n_max in [1, 2, 3]:
        r = run_analysis(n_max, Lambda_sq=6.0, validate=(n_max == 2))
        all_results[f"n_max={n_max}"] = r

    # n_max=4 if feasible
    try:
        r = run_analysis(4, Lambda_sq=6.0)
        all_results["n_max=4"] = r
    except Exception as e:
        print(f"\n  n_max=4 skipped: {e}")

    lichnerowicz_comparison(all_results)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "g6_fierz_pauli.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
