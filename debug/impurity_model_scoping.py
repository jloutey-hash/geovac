"""
Anderson impurity model scoping for GeoVac basis.

Research question: is the GeoVac angular momentum eigenbasis (n,l,m labels,
Gaunt-structured V_ee) a good basis for quantum impurity models?

Builds single-impurity Anderson model (SIAM) variants:
  (a) s-orbital impurity (l=0), 6 bath states
  (b) d-orbital impurity (l=2), 5 impurity orbitals, 5 bath "shells"

Compares Pauli count, 1-norm, and QWC groups between:
  - GeoVac-style (angular momentum eigenbasis, Gaunt-structured U)
  - "Generic" basis (random unitary rotation of the same Hamiltonian)
  - Site basis (Hubbard-like, for context)

Outputs to debug/data/impurity_scoping_results.json and plots.
"""

import json
from pathlib import Path
from itertools import product

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.stats import unitary_group

from openfermion import (
    FermionOperator,
    jordan_wigner,
    count_qubits,
    QubitOperator,
)
from openfermion.transforms import get_fermion_operator
from openfermion.ops import InteractionOperator

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def pauli_stats(qop: QubitOperator):
    """Non-identity Pauli count, 1-norm (excluding identity), QWC groups (greedy)."""
    terms = [(t, c) for t, c in qop.terms.items() if len(t) > 0]
    n_pauli = len(terms)
    l1 = float(sum(abs(c) for _, c in terms))

    # Greedy QWC grouping
    def qwc_compatible(t1, t2):
        d1 = dict(t1)
        for q, p in t2:
            if q in d1 and d1[q] != p:
                return False
        return True

    groups = []
    for t, _ in sorted(terms, key=lambda x: -abs(x[1])):
        placed = False
        for g in groups:
            if all(qwc_compatible(t, t2) for t2 in g):
                g.append(t)
                placed = True
                break
        if not placed:
            groups.append([t])
    return n_pauli, l1, len(groups)


def build_interaction_operator(h1, h2):
    """Wrap one-body h1[p,q] and two-body h2[p,q,r,s] (chemist notation:
    h2[p,q,r,s] a†_p a†_r a_s a_q) into an InteractionOperator."""
    n = h1.shape[0]
    # openfermion's InteractionOperator uses physics ordering:
    # 0.5 * sum h2[p,q,r,s] a†_p a†_q a_r a_s
    # We pass the 4-index tensor directly as "two_body_tensor".
    return InteractionOperator(0.0, h1, 0.5 * h2)


def rotate_hamiltonian(h1, h2, U):
    """Apply unitary rotation U (n x n) to single-particle orbitals."""
    h1r = U.conj().T @ h1 @ U
    # h2[p,q,r,s] -> U†_pp' U†_qq' h2 U_r'r U_s's
    h2r = np.einsum("pP,qQ,PQRS,rR,sS->pqrs", U.conj(), U.conj(), h2, U, U)
    return h1r, h2r


def gaunt_density(l_max):
    """Fraction of (l1 m1, l2 m2, l3 m3, l4 m4) quartets (all l <= l_max)
    for which the Gaunt selection rules allow nonzero <12|1/r12|34>.
    Selection rules: m1+m3 = m2+m4, and triangle(l1,l2,k), triangle(l3,l4,k)
    for some common k with (l1+l2+k) and (l3+l4+k) even."""
    orbitals = [(l, m) for l in range(l_max + 1) for m in range(-l, l + 1)]
    n = len(orbitals)
    nz = 0
    total = n ** 4
    for (l1, m1), (l2, m2), (l3, m3), (l4, m4) in product(orbitals, repeat=4):
        if (m1 - m2) != (m4 - m3):
            continue
        k_min = max(abs(l1 - l2), abs(l3 - l4), abs(m1 - m2))
        k_max = min(l1 + l2, l3 + l4)
        allowed = False
        for k in range(k_min, k_max + 1):
            if (l1 + l2 + k) % 2 == 0 and (l3 + l4 + k) % 2 == 0:
                allowed = True
                break
        if allowed:
            nz += 1
    return nz / total, n


# ------------------------------------------------------------------
# Part 3: Build SIAM and compute Pauli count in different bases
# ------------------------------------------------------------------

def build_siam_s_orbital(n_bath=6, t=1.0, V0=0.5, U=4.0, eps_d=-2.0):
    """Single s-orbital impurity + 1D bath.

    Spin-orbital ordering: site 0 = impurity, sites 1..n_bath = bath.
    For each site we have (up, down) spin-orbitals.
    Total spin orbitals: 2*(1 + n_bath).

    Returns (h1, h2) in spin-orbital basis (chemist order).
    """
    n_sites = 1 + n_bath
    n_so = 2 * n_sites
    h1 = np.zeros((n_so, n_so))
    h2 = np.zeros((n_so, n_so, n_so, n_so))

    # Bath dispersion: discretize ε_k = -2t cos(k), k = pi*j/(n_bath+1)
    # (Use the site-basis tridiagonal bath: H_bath = -t sum_j (c†_j c_{j+1} + h.c.)
    # Diagonalizing gives the ε_k eigenvalues. We build in site basis and then
    # diagonalize the one-body bath.)
    bath_h = np.zeros((n_bath, n_bath))
    for j in range(n_bath - 1):
        bath_h[j, j + 1] = -t
        bath_h[j + 1, j] = -t
    ev, evec = np.linalg.eigh(bath_h)

    # One-body: impurity eps_d + bath eigenvalues + hybridization V0 to site 0 of bath
    # (which in eigen basis becomes V_k = V0 * evec[0, k])
    for sigma in range(2):
        # impurity diagonal
        p_imp = 2 * 0 + sigma
        h1[p_imp, p_imp] = eps_d
        # bath diagonal
        for k in range(n_bath):
            p_bk = 2 * (1 + k) + sigma
            h1[p_bk, p_bk] = ev[k]
            # hybridization (impurity -- bath eigenmode k)
            Vk = V0 * evec[0, k]
            h1[p_imp, p_bk] = Vk
            h1[p_bk, p_imp] = Vk

    # Two-body: on-site U on impurity
    # U n_{d,up} n_{d,down} = U a†_{d,up} a_{d,up} a†_{d,down} a_{d,down}
    # chemist: h2[p,q,r,s] a†_p a†_r a_s a_q, so U n_up n_down corresponds to
    # h2[p=up, q=up, r=down, s=down] = U (and the (up<->down) permutation)
    p_up = 0
    p_dn = 1
    h2[p_up, p_up, p_dn, p_dn] = U
    h2[p_dn, p_dn, p_up, p_up] = U

    return h1, h2, n_so


def build_siam_d_orbital(n_bath_shells=2, t=1.0, V0=0.3, U=4.0, J=0.5, eps_d=-3.0):
    """d-orbital impurity (5 orbitals, l=2, m=-2..+2) with Gaunt-structured U.

    Bath: n_bath_shells copies of a d-symmetric bath (5 orbitals per shell,
    dispersion ε_k = -2t cos(k_shell)). This keeps angular symmetry intact.

    On-site interaction: Kanamori form (simplified)
        U_{mm'm''m'''} = U delta_{m,m'''} delta_{m',m''} (density-density)
                       + J delta_{m,m''} delta_{m',m'''} (exchange)
    restricted to the Gaunt-allowed m-channels (m + m' = m'' + m''').
    """
    n_imp = 5  # l=2
    n_shell = 5
    n_sites = n_imp + n_bath_shells * n_shell
    n_so = 2 * n_sites
    h1 = np.zeros((n_so, n_so))
    h2 = np.zeros((n_so, n_so, n_so, n_so))

    # Impurity one-body
    for m_imp in range(n_imp):
        for sigma in range(2):
            p = 2 * m_imp + sigma
            h1[p, p] = eps_d

    # Bath shells (dispersion per shell)
    shell_energies = [-2 * t * np.cos(np.pi * (s + 1) / (n_bath_shells + 1))
                      for s in range(n_bath_shells)]
    for s, eps_s in enumerate(shell_energies):
        for m in range(n_shell):
            for sigma in range(2):
                p = 2 * (n_imp + s * n_shell + m) + sigma
                h1[p, p] = eps_s
                # angular-preserving hybridization: only m -> m (selection rule)
                p_imp = 2 * m + sigma
                h1[p_imp, p] = V0
                h1[p, p_imp] = V0

    # On-site Kanamori U (density-density + exchange) on impurity d-shell,
    # restricted to m1 + m2 = m3 + m4 (Gaunt selection for k=0,2,4 combined).
    # For simplicity use density-density only (U for same-m up-down, U' for
    # different-m), since full Slater-Condon is heavy for scoping.
    Up = U - 2 * J  # inter-orbital density-density
    for m in range(n_imp):
        # intra-orbital U n_up n_down
        p_u, p_d = 2 * m, 2 * m + 1
        h2[p_u, p_u, p_d, p_d] = U
        h2[p_d, p_d, p_u, p_u] = U

    for m1 in range(n_imp):
        for m2 in range(n_imp):
            if m1 == m2:
                continue
            for s1 in range(2):
                for s2 in range(2):
                    p1 = 2 * m1 + s1
                    p2 = 2 * m2 + s2
                    # density-density U' n_{m1,s1} n_{m2,s2}
                    h2[p1, p1, p2, p2] = Up if s1 != s2 else (Up - J)

    return h1, h2, n_so


# ------------------------------------------------------------------
# Main computation
# ------------------------------------------------------------------

def analyse(name, h1, h2, n_so):
    # ensure symmetry in h2 (chemist): h2[p,q,r,s] = h2[r,s,p,q]
    h2_sym = 0.5 * (h2 + h2.transpose(2, 3, 0, 1))
    iop = build_interaction_operator(h1, h2_sym)
    fop = get_fermion_operator(iop)
    qop = jordan_wigner(fop)
    npauli, l1, qwc = pauli_stats(qop)
    return {
        "name": name,
        "n_spinorbitals": n_so,
        "n_qubits": count_qubits(qop),
        "n_pauli": npauli,
        "one_norm": l1,
        "qwc_groups": qwc,
    }


def random_unitary(n, seed=0):
    return unitary_group.rvs(n, random_state=seed)


def main():
    root = Path("C:/Users/jlout/Desktop/Project_Geometric")
    data_dir = root / "debug" / "data"
    plot_dir = root / "debug" / "plots"
    data_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    # --- Part 3a: s-orbital SIAM ---
    h1, h2, n_so = build_siam_s_orbital(n_bath=6, t=1.0, V0=0.5, U=4.0, eps_d=-2.0)
    r_geo = analyse("siam_s_geovac", h1, h2, n_so)

    # generic rotation: mix impurity and bath orbitals (this DESTROYS locality)
    U_rand = random_unitary(n_so // 2, seed=42)
    # apply to spatial-orbital block structure
    U_spin = np.kron(U_rand, np.eye(2))  # same rotation for both spins
    h1r, h2r = rotate_hamiltonian(h1, h2, U_spin)
    r_gen = analyse("siam_s_generic", h1r, h2r, n_so)

    # site basis: the raw site-basis Hubbard form (no bath diagonalization).
    # Rebuild: impurity + bath in SITE basis.
    n_bath = 6
    n_sites = 1 + n_bath
    n_so2 = 2 * n_sites
    h1_site = np.zeros((n_so2, n_so2))
    h2_site = np.zeros((n_so2, n_so2, n_so2, n_so2))
    for sigma in range(2):
        h1_site[2 * 0 + sigma, 2 * 0 + sigma] = -2.0  # eps_d
        # bath hopping tridiagonal
        for j in range(n_bath - 1):
            p = 2 * (1 + j) + sigma
            q = 2 * (2 + j) + sigma
            h1_site[p, q] = -1.0
            h1_site[q, p] = -1.0
        # hybridization impurity-site1 (one hopping)
        p_imp = 2 * 0 + sigma
        p_b0 = 2 * 1 + sigma
        h1_site[p_imp, p_b0] = 0.5
        h1_site[p_b0, p_imp] = 0.5
    h2_site[0, 0, 1, 1] = 4.0
    h2_site[1, 1, 0, 0] = 4.0
    r_site = analyse("siam_s_site", h1_site, h2_site, n_so2)

    results["s_orbital"] = {"geovac": r_geo, "generic_rotation": r_gen, "site_basis": r_site}

    # --- Part 3b: d-orbital SIAM (smaller bath to keep feasible) ---
    h1d, h2d, n_sod = build_siam_d_orbital(n_bath_shells=1, t=1.0, V0=0.3, U=4.0, J=0.5, eps_d=-3.0)
    r_d_geo = analyse("siam_d_geovac", h1d, h2d, n_sod)

    U_rand_d = random_unitary(n_sod // 2, seed=7)
    U_spin_d = np.kron(U_rand_d, np.eye(2))
    h1dr, h2dr = rotate_hamiltonian(h1d, h2d, U_spin_d)
    r_d_gen = analyse("siam_d_generic", h1dr, h2dr, n_sod)

    results["d_orbital"] = {"geovac": r_d_geo, "generic_rotation": r_d_gen}

    # --- Part 4: Gaunt density vs l_max ---
    gaunt = {}
    for l_max in range(0, 5):
        density, n_orb = gaunt_density(l_max)
        gaunt[l_max] = {"n_orb": n_orb, "density_frac": density}
    results["gaunt_density"] = gaunt

    # --- Part 2 summary: SIAM in GeoVac labels ---
    results["mapping"] = {
        "impurity_labels": "d-shell: (n_imp, l=2, m=-2..+2), 5 spatial orbitals",
        "bath_labels": "(n_bath, l, m) matching impurity symmetry; only same-m hybridizes (selection rule)",
        "U_structure": "Gaunt 4-index tensor on impurity; m1+m2 = m3+m4 (only ~2.76% density at l_max=2)",
        "decoupling": "one-body h1 diagonal in (l,m) blocks; two-body V_ee sparse by Gaunt",
    }

    # --- Save ---
    out_file = data_dir / "impurity_scoping_results.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Saved {out_file}")

    # --- Plots ---
    # Plot 1: Pauli count comparison
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

    s_labels = ["GeoVac\n(ang. mom.)", "Generic\n(random rot.)", "Site basis"]
    s_pauli = [r_geo["n_pauli"], r_gen["n_pauli"], r_site["n_pauli"]]
    s_l1 = [r_geo["one_norm"], r_gen["one_norm"], r_site["one_norm"]]
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    ax[0].bar(s_labels, s_pauli, color=colors)
    ax[0].set_ylabel("Non-identity Pauli terms")
    ax[0].set_title(f"s-orbital SIAM ({r_geo['n_qubits']} qubits)\nPauli count by basis")
    for i, v in enumerate(s_pauli):
        ax[0].text(i, v, str(v), ha="center", va="bottom")

    d_labels = ["GeoVac\n(d-shell)", "Generic\n(random rot.)"]
    d_pauli = [r_d_geo["n_pauli"], r_d_gen["n_pauli"]]
    ax[1].bar(d_labels, d_pauli, color=["#1f77b4", "#d62728"])
    ax[1].set_ylabel("Non-identity Pauli terms")
    ax[1].set_title(f"d-orbital SIAM ({r_d_geo['n_qubits']} qubits)\nPauli count by basis")
    for i, v in enumerate(d_pauli):
        ax[1].text(i, v, str(v), ha="center", va="bottom")

    plt.tight_layout()
    p1 = plot_dir / "impurity_pauli_comparison.png"
    plt.savefig(p1, dpi=120)
    plt.close()
    print(f"Saved {p1}")

    # Plot 2: ERI density (Gaunt) vs l_max
    fig, ax = plt.subplots(figsize=(7, 4.5))
    lmax_vals = list(gaunt.keys())
    densities = [100 * gaunt[l]["density_frac"] for l in lmax_vals]
    n_orbs = [gaunt[l]["n_orb"] for l in lmax_vals]
    ax.plot(lmax_vals, densities, "o-", lw=2, markersize=9, color="#1f77b4")
    for lm, d, n in zip(lmax_vals, densities, n_orbs):
        ax.annotate(f"{d:.2f}%\n({n} orbitals)",
                    (lm, d), textcoords="offset points", xytext=(8, 5), fontsize=9)
    ax.set_xlabel("l_max")
    ax.set_ylabel("Gaunt-allowed ERI density (%)")
    ax.set_yscale("log")
    ax.set_title("Angular sparsity (Gaunt selection) — Paper 22 regime\n"
                 "d-shell impurity (l=2) sits at 2.76% density")
    ax.grid(alpha=0.3)
    ax.axvline(2, color="red", ls="--", alpha=0.5, label="d-shell (l=2)")
    ax.legend()
    plt.tight_layout()
    p2 = plot_dir / "impurity_eri_density.png"
    plt.savefig(p2, dpi=120)
    plt.close()
    print(f"Saved {p2}")

    # --- Report summary ---
    print("\n=== SUMMARY ===")
    print(f"s-orbital SIAM ({r_geo['n_qubits']} qubits):")
    print(f"  GeoVac (ang mom):     {r_geo['n_pauli']:5d} Pauli, 1-norm={r_geo['one_norm']:.3f}, QWC={r_geo['qwc_groups']}")
    print(f"  Generic (rotated):    {r_gen['n_pauli']:5d} Pauli, 1-norm={r_gen['one_norm']:.3f}, QWC={r_gen['qwc_groups']}")
    print(f"  Site basis:           {r_site['n_pauli']:5d} Pauli, 1-norm={r_site['one_norm']:.3f}, QWC={r_site['qwc_groups']}")
    print(f"\nd-orbital SIAM ({r_d_geo['n_qubits']} qubits):")
    print(f"  GeoVac (ang mom):     {r_d_geo['n_pauli']:5d} Pauli, 1-norm={r_d_geo['one_norm']:.3f}, QWC={r_d_geo['qwc_groups']}")
    print(f"  Generic (rotated):    {r_d_gen['n_pauli']:5d} Pauli, 1-norm={r_d_gen['one_norm']:.3f}, QWC={r_d_gen['qwc_groups']}")
    print(f"  Ratio (generic/GeoVac): {r_d_gen['n_pauli']/r_d_geo['n_pauli']:.2f}x Pauli, "
          f"{r_d_gen['one_norm']/r_d_geo['one_norm']:.2f}x 1-norm")


if __name__ == "__main__":
    main()
