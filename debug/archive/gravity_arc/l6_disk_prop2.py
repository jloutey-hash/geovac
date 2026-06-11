"""L6 Track 1 -- prop=2 confirmation on the assembled disk operator system.

Charter step 1 (debug/gravity_campaign_R2_theorem_charter_memo.md):
    "Confirm prop=2 on the assembled disk operator system (the phase-1-check
     follow-on compute)."

The disk-with-cone backbone factor is the 2-disk D^2.  Its truncated operator
system is O_n = P_n C(disk) P_n, the compression of continuous functions on the
disk to the truncated mode space {(k,j): |k|<=K, j<J}.  We build O_n from the
compressed disk-coordinate monomials z^a z-bar^b (z = rho e^{i phi}), and compute

    prop(O_n) = smallest p with dim(O_n^p) = N^2.

Framework benchmark (Paper 32 sec III / Paper 44): prop = 2 means O^2 already
spans the full M_N(C) envelope -- the "2-step generating" Connes-vS Toeplitz
property (Prop 4.2).  The azimuthal factor IS the Toeplitz S^1 (prop=2 verbatim,
Paper 44); this driver confirms the assembled disk (radial interval (x) azimuthal
S^1) inherits prop=2.

Radial modes: eigenfunctions u_{m,j} of the Hermitian polar Laplacian
-u'' + (m^2-1/4)/rho^2 on (0,R], a = R/N_rho (same convention as warped_dirac).
In the u = sqrt(rho) f representation the radial measure is d rho, so
    <m',j'| rho^p e^{i q phi} | m,j>  (azimuthal q = k'-k)
       = delta_{k', k+q} * sum_i u_{m',j'}(i) rho_i^p u_{m,j}(i) * a
with m = |k+1/2| (fermionic), m' = |k+q+1/2|.
"""

import json
import itertools
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh_tridiagonal, eigh_tridiagonal

from geovac.operator_system import operator_system_dim, _extract_matrix_basis

OUT = Path(__file__).parent / "data" / "l6_disk_prop2.json"
OUT.parent.mkdir(exist_ok=True)


def radial_modes(m: float, R: float, N_rho: int, J: int):
    """Lowest J eigen-(value, u-vector) of the centrifugal radial operator."""
    a = R / N_rho
    k = np.arange(1, N_rho + 1)
    rho = k * a
    diag = 2.0 / a**2 + (m**2 - 0.25) / rho**2
    off = -np.ones(N_rho - 1) / a**2
    w, v = eigh_tridiagonal(diag, off, select="i", select_range=(0, J - 1))
    return w, v, rho, a


def build_disk_operator_system(K: int, J: int, R: float, N_rho: int,
                               max_deg: int):
    """Generators of O_n = P_n C(disk) P_n for |k|<=K-? , j<J.

    Azimuthal modes k in {-K, ..., K-1} (2K modes); radial j in {0..J-1}.
    Mode index n_idx = (k_index)*J + j with k_index = k+K.
    Returns list of generator matrices (compressed monomial multiplications).
    """
    ks = list(range(-K, K))          # 2K azimuthal modes
    n_modes = len(ks) * J
    # radial mode vectors per |m| = |k+1/2|
    mode_vecs = {}
    rho = None
    a = None
    for k in ks:
        m = abs(k + 0.5)
        w, v, rho, a = radial_modes(m, R, N_rho, J)
        mode_vecs[k] = v   # columns = u_{m,j}, j=0..J-1

    def idx(k, j):
        return (k + K) * J + j

    def mult_matrix(p: int, q: int):
        """Compressed multiplication by rho^p e^{i q phi}: shifts k -> k+q."""
        M = np.zeros((n_modes, n_modes), dtype=complex)
        for k in ks:
            kp = k + q
            if kp not in ks:
                continue
            Vk = mode_vecs[k]      # (N_rho, J)
            Vkp = mode_vecs[kp]    # (N_rho, J)
            # radial overlap matrix O[j',j] = sum_i Vkp[i,j'] rho_i^p Vk[i,j] * a
            wgt = (rho ** p) * a
            Ovl = (Vkp * wgt[:, None]).T @ Vk    # (J, J)
            for jp in range(J):
                for j in range(J):
                    M[idx(kp, jp), idx(k, j)] = Ovl[jp, j]
        return M

    # Monomials z^a zbar^b -> rho^{a+b} e^{i(a-b)phi}: q = a-b, p = a+b.
    gens = [np.eye(n_modes, dtype=complex)]
    seen = set()
    for a_ in range(max_deg + 1):
        for b_ in range(max_deg + 1):
            if a_ + b_ == 0:
                continue
            if a_ + b_ > max_deg:
                continue
            q = a_ - b_
            p = a_ + b_
            if abs(q) >= 2 * K:   # azimuthal shift outside truncation -> zero
                continue
            key = (p, q)
            if key in seen:
                continue
            seen.add(key)
            gens.append(mult_matrix(p, q))
    return gens, n_modes


def prop_from_generators(gens, max_k=8, tol=1e-9, verbose=False):
    """prop = smallest k with dim(O^k) = N^2 (mirrors operator_system.propagation_number)."""
    N = gens[0].shape[0]
    target = N * N
    dims = []
    cur = list(gens)
    d1 = operator_system_dim(cur, tol=tol)
    dims.append(d1)
    if verbose:
        print(f"    k=1 dim(O)={d1}/{target}")
    if d1 == target:
        return 1, dims
    for k in range(2, max_k + 1):
        basis = _extract_matrix_basis(cur, tol=tol)
        newg = [A @ B for A in basis for B in gens]
        dk = operator_system_dim(newg, tol=tol)
        dims.append(dk)
        if verbose:
            print(f"    k={k} dim(O^{k})={dk}/{target}")
        if dk == target:
            return k, dims
        if dk == dims[-2]:
            return -1, dims     # saturated below target
        cur = newg
    return -1, dims


def main():
    res = {}
    print("=" * 72)
    print("L6 Track 1 -- prop=2 on the assembled disk operator system")
    print("=" * 72)
    R, N_rho, max_deg = 8.0, 120, 8
    grid = [(2, 2), (2, 3), (3, 2), (3, 3)]
    rows = []
    for (K, J) in grid:
        gens, N = build_disk_operator_system(K, J, R, N_rho, max_deg)
        prop, dims = prop_from_generators(gens, verbose=True)
        full = N * N
        print(f"  (K={K}, J={J})  N={N}  |gens|={len(gens)}  "
              f"dim(O)={dims[0]}  dim(O^2)={dims[1] if len(dims)>1 else '-'}/{full}  "
              f"=> prop={prop}")
        rows.append(dict(K=K, J=J, N=N, n_gens=len(gens), dim_O=dims[0],
                         dim_O2=(dims[1] if len(dims) > 1 else None),
                         envelope=full, dim_sequence=dims, prop=prop))
    res["rows"] = rows
    all2 = all(r["prop"] == 2 for r in rows)
    # also report whether O^2 fills the envelope (the 2-step generating signature)
    print("\n" + "=" * 72)
    verdict = "PROP2-CONFIRMED" if all2 else "PROP2-MIXED"
    print(f"[Verdict] {verdict}")
    for r in rows:
        sig = "O^2 = full envelope" if r["dim_O2"] == r["envelope"] else \
              f"O^2 = {r['dim_O2']}/{r['envelope']}"
        print(f"  (K={r['K']},J={r['J']}) prop={r['prop']}  ({sig})")
    res["verdict"] = verdict
    res["all_prop2"] = all2
    with OUT.open("w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print(f"\nsaved {OUT}")
    print("=" * 72)


if __name__ == "__main__":
    main()
