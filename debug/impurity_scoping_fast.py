"""
Anderson impurity model scoping in the GeoVac angular momentum basis.

Fast go/no-go investigation: does Gaunt-structured sparsity help for
impurity models?

Part 1: s-orbital SIAM (no angular structure).
Part 2: d-orbital SIAM (Slater-Condon U with Gaunt selection rules).
Part 3: Compare to random U(N) rotation of the orbital basis.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from openfermion import FermionOperator, jordan_wigner
from sympy.physics.wigner import gaunt as sp_gaunt

RNG = np.random.default_rng(42)

OUT_DATA = Path(__file__).parent / "data" / "impurity_scoping_results.json"
OUT_PLOT = Path(__file__).parent / "plots" / "impurity_pauli_comparison.png"
OUT_DATA.parent.mkdir(parents=True, exist_ok=True)
OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def spin_orbital(spatial: int, spin: int) -> int:
    """Map (spatial, spin in {0,1}) -> spin-orbital index (interleaved)."""
    return 2 * spatial + spin


def pauli_1norm(qubit_op) -> float:
    return sum(abs(c) for c in qubit_op.terms.values())


def qwc_groups_greedy(qubit_op) -> int:
    """Count QWC groups via greedy coloring (qubit-wise commuting)."""
    terms = [t for t in qubit_op.terms.keys() if len(t) > 0]
    groups: list[list[tuple]] = []

    def qwc(t1, t2):
        d1 = dict(t1)
        for q, p in t2:
            if q in d1 and d1[q] != p:
                return False
        return True

    for t in terms:
        placed = False
        for g in groups:
            if all(qwc(t, s) for s in g):
                g.append(t)
                placed = True
                break
        if not placed:
            groups.append([t])
    return len(groups)


def count_metrics(fop: FermionOperator) -> dict:
    qop = jordan_wigner(fop)
    n_pauli = len([t for t in qop.terms if len(t) > 0])
    l1 = pauli_1norm(qop)
    return {
        "n_pauli": n_pauli,
        "l1_norm": float(l1),
        "qwc_groups": qwc_groups_greedy(qop),
    }


def rotate_1body(h: np.ndarray, U: np.ndarray) -> np.ndarray:
    """h'_pq = sum_ab U*_ap h_ab U_bq   (unitary similarity in MO basis)."""
    return U.conj().T @ h @ U


def rotate_2body(v: np.ndarray, U: np.ndarray) -> np.ndarray:
    """v'_{pqrs} = sum U*_ap U*_br v_{abcd} U_cq U_ds  (chemist notation)."""
    # v[a,b,c,d] -> v'[p,q,r,s]
    w = np.einsum("ap,abcd->pbcd", U.conj(), v)
    w = np.einsum("bq,pbcd->pqcd", U.conj(), w)
    w = np.einsum("cr,pqcd->pqrd", U, w)
    w = np.einsum("ds,pqrd->pqrs", U, w)
    return w


def build_fermion_op(h: np.ndarray, v: np.ndarray, thresh: float = 1e-10) -> FermionOperator:
    """Build FermionOperator from 1-body h_pq and 2-body v_{pqrs} (chemist).

    H = sum_{p,q,sigma} h_pq a†_{p,sigma} a_{q,sigma}
      + (1/2) sum_{p,q,r,s,sigma,tau} v_{pqrs} a†_{p,sigma} a†_{r,tau} a_{s,tau} a_{q,sigma}
    """
    n = h.shape[0]
    fop = FermionOperator()
    # One-body
    for p in range(n):
        for q in range(n):
            val = h[p, q]
            if abs(val) < thresh:
                continue
            for sig in (0, 1):
                fop += FermionOperator(((spin_orbital(p, sig), 1),
                                        (spin_orbital(q, sig), 0)), val)
    # Two-body
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    val = v[p, q, r, s]
                    if abs(val) < thresh:
                        continue
                    for sig in (0, 1):
                        for tau in (0, 1):
                            if sig == tau and (p == r or q == s):
                                # Pauli exclusion handled by op algebra,
                                # but skip definitely-zero fast cases:
                                pass
                            fop += FermionOperator((
                                (spin_orbital(p, sig), 1),
                                (spin_orbital(r, tau), 1),
                                (spin_orbital(s, tau), 0),
                                (spin_orbital(q, sig), 0),
                            ), 0.5 * val)
    return fop


def random_unitary(n: int, rng) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    d = np.diag(r) / np.abs(np.diag(r))
    return q * d


# ---------------------------------------------------------------------------
# Part 1: s-orbital SIAM
# ---------------------------------------------------------------------------

def build_s_siam():
    """1 impurity (spatial 0) + 4 bath (spatial 1..4). Total n=5 spatial."""
    t = 1.0
    eps_d = 0.0
    V = 0.3 * t
    U = 4.0 * t
    n_bath = 4
    n = 1 + n_bath

    h = np.zeros((n, n))
    h[0, 0] = eps_d
    for k in range(1, n_bath + 1):
        h[k, k] = -2.0 * t * np.cos(np.pi * k / (n_bath + 1))
        h[0, k] = V
        h[k, 0] = V

    # 2-body: only U on impurity. U n_{0↑} n_{0↓}.
    # In chemist notation v_{pqrs}: U * a†_{0↑}a_{0↑} * a†_{0↓}a_{0↓}
    # corresponds to v[0,0,0,0] = U (with both sigma,tau channels picking it up)
    # (1/2) * v_{0000} * sum_{sig,tau} n_sig n_tau_eff -> U n_up n_dn requires v[0,0,0,0] = U
    v = np.zeros((n, n, n, n))
    v[0, 0, 0, 0] = U
    return h, v, n


# ---------------------------------------------------------------------------
# Part 2: d-orbital SIAM
# ---------------------------------------------------------------------------

def gaunt_tensor_l2(F0=4.0, F2=0.8, F4=0.5):
    """Build U_{m1 m2 m3 m4} for l=2 using Slater-Condon:

    U_{m1 m2 m3 m4} = sum_k c^k(l,m1; l,m2) * ... * F^k
    Simplified: use Gaunt coefficients directly.

        <m1 m3 | 1/r12 | m2 m4> = sum_k (4π/(2k+1)) F^k
                                  * <Y_lm1|Y_kq|Y_lm2><Y_lm3|Y_k,-q|Y_lm4>
    with q = m1 - m2 = m4 - m3.

    We use sympy.physics.wigner.gaunt(l1,l2,l3,m1,m2,m3) which returns
    ∫ Y_{l1,m1} Y_{l2,m2} Y_{l3,m3} dΩ.

    The ERI in the l=2 subspace is:
        v[m1,m2,m3,m4] = sum_{k in {0,2,4}} F^k *
            sum_q gaunt(l,k,l, -m1, q, m2) * gaunt(l,k,l, -m3, -q, m4)
            * (-1)^{m1+m3} * (4π/(2k+1))
    (sign from Y*_{l,m} = (-1)^m Y_{l,-m}).

    We just compute numerically and rely on selection rules.
    """
    l = 2
    ms = list(range(-l, l + 1))  # -2..+2
    dim = 2 * l + 1  # 5
    Fk = {0: F0, 2: F2, 4: F4}
    v = np.zeros((dim, dim, dim, dim))

    # index: m -> m+l (0..4)
    def idx(m):
        return m + l

    for m1 in ms:
        for m2 in ms:
            for m3 in ms:
                for m4 in ms:
                    # m-conservation
                    if m1 + m3 != m2 + m4:
                        continue
                    q = m1 - m2
                    total = 0.0
                    for k in (0, 2, 4):
                        if abs(q) > k:
                            continue
                        # gaunt returns rational/symbolic; cast to float
                        g1 = float(sp_gaunt(l, k, l, -m1, q, m2))
                        g2 = float(sp_gaunt(l, k, l, -m3, -q, m4))
                        if g1 == 0.0 or g2 == 0.0:
                            continue
                        prefac = 4.0 * np.pi / (2 * k + 1)
                        sign = (-1) ** (m1 + m3)
                        total += Fk[k] * prefac * sign * g1 * g2
                    v[idx(m1), idx(m2), idx(m3), idx(m4)] = total
    return v


def build_d_siam():
    """Impurity: 5 d-orbitals (spatial 0..4).
    Bath: 4 sites (spatial 5..8), hybridize selectively to m=0,±1 (simple).
    Total n = 9 spatial -> 18 qubits.
    """
    t = 1.0
    eps_d = 0.0
    V = 0.3 * t
    n_d = 5
    n_bath = 4
    n = n_d + n_bath

    h = np.zeros((n, n))
    for m_idx in range(n_d):
        h[m_idx, m_idx] = eps_d

    for k in range(n_bath):
        site = n_d + k
        h[site, site] = -2.0 * t * np.cos(np.pi * (k + 1) / (n_bath + 1))
        # hybridize bath site k to d-orbital (k mod 5) — keeps angular character
        m_idx = k % n_d
        h[site, m_idx] = V
        h[m_idx, site] = V

    # 2-body: Slater-Condon on impurity d-block, zero elsewhere
    v_d = gaunt_tensor_l2(F0=4.0, F2=0.8, F4=0.5)
    v = np.zeros((n, n, n, n))
    v[:n_d, :n_d, :n_d, :n_d] = v_d
    return h, v, n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run():
    results = {}
    t0 = time.time()

    # --- Part 1: s-orbital SIAM ---
    print("[Part 1] s-orbital SIAM (Q=10)...")
    h_s, v_s, n_s = build_s_siam()
    fop_s = build_fermion_op(h_s, v_s)
    m_s_geo = count_metrics(fop_s)
    print(f"  GeoVac (site) basis: {m_s_geo}")

    # Rotated basis
    n_samples = 5
    rotated_s = []
    for i in range(n_samples):
        U = random_unitary(n_s, RNG)
        h_r = rotate_1body(h_s, U).real
        v_r = rotate_2body(v_s.astype(complex), U).real
        # symmetrize small imaginary residue
        fop_r = build_fermion_op(h_r, v_r, thresh=1e-8)
        rotated_s.append(count_metrics(fop_r))
        print(f"  rotated sample {i}: {rotated_s[-1]}")

    # --- Part 2: d-orbital SIAM ---
    print(f"\n[Part 2] d-orbital SIAM (Q={2*9}=18)...")
    h_d, v_d, n_d = build_d_siam()

    # Count non-zero Gaunt U elements on impurity block
    v_imp = v_d[:5, :5, :5, :5]
    nnz_gaunt = int(np.sum(np.abs(v_imp) > 1e-10))
    print(f"  non-zero impurity U elements: {nnz_gaunt} / 625")

    print("  Building GeoVac-basis Hamiltonian...")
    t1 = time.time()
    fop_d = build_fermion_op(h_d, v_d)
    m_d_geo = count_metrics(fop_d)
    print(f"  GeoVac basis: {m_d_geo}  ({time.time()-t1:.1f}s)")

    # Rotated d-orbital basis — rotate only the d-block (5 orbitals) to keep
    # bath hybridization structure but scramble the angular character.
    # This is the comparison that matters: the sparsity claim is about the
    # angular basis on the impurity.
    print("  Random U(5) rotation of d-orbital block...")
    rotated_d = []
    for i in range(n_samples):
        U5 = random_unitary(5, RNG)
        U_full = np.eye(n_d, dtype=complex)
        U_full[:5, :5] = U5
        h_r = rotate_1body(h_d.astype(complex), U_full).real
        v_r = rotate_2body(v_d.astype(complex), U_full).real
        t1 = time.time()
        fop_r = build_fermion_op(h_r, v_r, thresh=1e-8)
        m = count_metrics(fop_r)
        rotated_d.append(m)
        print(f"  rotated sample {i}: {m}  ({time.time()-t1:.1f}s)")

    results = {
        "s_siam": {
            "n_spatial": n_s,
            "n_qubits": 2 * n_s,
            "geovac": m_s_geo,
            "rotated_samples": rotated_s,
            "rotated_mean_pauli": float(np.mean([x["n_pauli"] for x in rotated_s])),
        },
        "d_siam": {
            "n_spatial": n_d,
            "n_qubits": 2 * n_d,
            "nnz_gaunt_impurity_U": nnz_gaunt,
            "nnz_ratio": nnz_gaunt / 625.0,
            "geovac": m_d_geo,
            "rotated_samples": rotated_d,
            "rotated_mean_pauli": float(np.mean([x["n_pauli"] for x in rotated_d])),
        },
        "sparsity_ratio_s": (np.mean([x["n_pauli"] for x in rotated_s])
                             / m_s_geo["n_pauli"]),
        "sparsity_ratio_d": (np.mean([x["n_pauli"] for x in rotated_d])
                             / m_d_geo["n_pauli"]),
        "runtime_seconds": time.time() - t0,
    }
    print(f"\nTotal runtime: {results['runtime_seconds']:.1f}s")
    print(f"Sparsity ratio s-SIAM (rotated/GeoVac Pauli): {results['sparsity_ratio_s']:.2f}x")
    print(f"Sparsity ratio d-SIAM (rotated/GeoVac Pauli): {results['sparsity_ratio_d']:.2f}x")

    with open(OUT_DATA, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved: {OUT_DATA}")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(7, 5))
    labels = ["s-SIAM\n(Q=10)", "d-SIAM\n(Q=18)"]
    geo_vals = [m_s_geo["n_pauli"], m_d_geo["n_pauli"]]
    rot_vals = [np.mean([x["n_pauli"] for x in rotated_s]),
                np.mean([x["n_pauli"] for x in rotated_d])]
    rot_std = [np.std([x["n_pauli"] for x in rotated_s]),
               np.std([x["n_pauli"] for x in rotated_d])]
    x = np.arange(len(labels))
    w = 0.35
    ax.bar(x - w/2, geo_vals, w, label="GeoVac angular basis", color="C0")
    ax.bar(x + w/2, rot_vals, w, yerr=rot_std, label="Random-rotated basis",
           color="C3", capsize=5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Number of Pauli terms")
    ax.set_title("Anderson impurity Pauli count: GeoVac vs rotated basis")
    ax.legend()
    for i, (g, r) in enumerate(zip(geo_vals, rot_vals)):
        ax.annotate(f"{r/g:.2f}x", xy=(i, max(g, r)), xytext=(0, 5),
                    textcoords="offset points", ha="center", fontsize=10,
                    fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=120)
    print(f"Saved: {OUT_PLOT}")
    return results


if __name__ == "__main__":
    run()
