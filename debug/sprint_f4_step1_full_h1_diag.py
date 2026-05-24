"""Sprint F4 Step 1 EXTENDED — full F3 h1 (all 10 orbitals at max_n=2).

The minimal Step 1 diagnostic restricted to {|Na 3s>, |H 1s>}. The actual
F3 build at max_n=2 has 10 spatial orbitals (5 Na + 5 H per shell n=1, 2):
Na: (1,0,0), (2,0,0), (2,1,-1), (2,1,0), (2,1,1)
H:  (1,0,0), (2,0,0), (2,1,-1), (2,1,0), (2,1,1)

The actual F3 bonding orbital from the FULL F3 h1 diagonalization may
have different composition than the 2x2 estimate above. Let's compute
its core overlaps explicitly.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from geovac.balanced_coupled import build_balanced_hamiltonian  # noqa: E402
from geovac.molecular_spec import nah_spec  # noqa: E402


def main():
    R_EQ = 3.566

    print("=" * 78)
    print("Sprint F4 Step 1 EXTENDED — full F3 h1 (10 orbital, max_n=2)")
    print(f"NaH R = {R_EQ} bohr; W1c + multi-zeta + cross-block h1")
    print("=" * 78)

    t0 = time.perf_counter()

    spec = nah_spec(max_n=2)
    nuclei = [
        {'Z': 11.0, 'position': (0.0, 0.0, 0.0), 'label': 'Na'},
        {'Z': 1.0, 'position': (0.0, 0.0, R_EQ), 'label': 'H'},
    ]

    # Build full F3 stack on NaH
    result = build_balanced_hamiltonian(
        spec=spec, R=R_EQ, nuclei=nuclei,
        screened_cross_center=True,
        multi_zeta_basis=True,
        cross_block_h1=True,
        verbose=False,
    )
    h1 = result['h1']
    M = result['M']
    print(f"M={M} spatial orbitals, h1 shape = {h1.shape}")

    # Find sub-block geometry to identify the bonding orbital structure
    from geovac.balanced_coupled import _get_block_geometry
    sub_blocks = _get_block_geometry(spec)
    print(f"\nSub-blocks: {len(sub_blocks)}")
    orbital_labels = []
    for sb in sub_blocks:
        for (n, l, m) in sb['states']:
            orbital_labels.append(f"{sb['label']}({n},{l},{m})")
    for i, lab in enumerate(orbital_labels):
        print(f"  [{i}] {lab}")

    # Diagonalize h1
    w, v = eigh(h1)
    print(f"\nh1 eigenvalues:")
    for i in range(min(M, 6)):
        # Identify dominant component
        idx_max = np.argmax(np.abs(v[:, i]))
        print(f"  E_{i} = {w[i]:+.4f} Ha, dominant: {orbital_labels[idx_max]} "
              f"(weight {v[idx_max, i]**2:.3f})")

    # Find the bonding orbital: lowest with substantial mixing on Na and H sides
    bond_idx = None
    bond_eigvec = None
    for i in range(M):
        # Compute fraction on Na vs H sides
        na_orbs = [k for k, lab in enumerate(orbital_labels) if lab.startswith('Na')]
        h_orbs = [k for k, lab in enumerate(orbital_labels) if lab.startswith('H')]
        na_frac = float(np.sum(v[na_orbs, i] ** 2))
        h_frac = float(np.sum(v[h_orbs, i] ** 2))
        if min(na_frac, h_frac) > 0.1:
            bond_idx = i
            bond_eigvec = v[:, i]
            print(f"\nFound bonding orbital at index {i}: E = {w[i]:+.4f} Ha")
            print(f"  Na-side fraction = {na_frac:.3f}, H-side fraction = {h_frac:.3f}")
            print(f"  Composition:")
            for k, lab in enumerate(orbital_labels):
                if abs(bond_eigvec[k]) > 0.05:
                    print(f"    {lab}: {bond_eigvec[k]:+.4f}")
            break

    if bond_idx is None:
        print("\nNo bonding orbital found in h1 spectrum (no eigenvector with >10% on both sides)")
        # Use lowest h1 eigenvector
        bond_idx = 0
        bond_eigvec = v[:, 0]
        print(f"Using lowest h1 eigenvector at E = {w[0]:+.4f} Ha")
        print(f"  Composition:")
        for k, lab in enumerate(orbital_labels):
            if abs(bond_eigvec[k]) > 0.05:
                print(f"    {lab}: {bond_eigvec[k]:+.4f}")

    # ------------------------------------------------------------------
    # Compute bonding-vs-core overlaps using same Step 1 machinery
    # ------------------------------------------------------------------
    print(f"\n{'-' * 78}")
    print("Overlap of bonding orbital with Na [Ne] core orbitals:")
    print(f"{'-' * 78}")

    from geovac.cross_block_h1 import (
        hydrogenic_R_nl_analytical,
        overlap_ss_axial,
    )
    from geovac.multi_zeta_orbitals import get_physical_valence_orbitals

    NA_CR_EXPONENTS = {(1, 0): 10.62608, (2, 0): 3.21766, (2, 1): 3.45200}
    NA_CORE_ENERGIES = {(1, 0): -40.4787, (2, 0): -2.7967, (2, 1): -1.5181}

    # Build per-orbital callables matching the F3 build:
    # On Na side, multi-zeta Na 3s replaces hydrogenic if (n=3, l=0); but
    # at max_n=2 the Na sub-block uses block_n in {1, 2} with n_val_offset=2,
    # mapping to physical n in {3, 4}. So block_n=1, l=0 -> Na 3s (mz),
    # block_n=2, l=0 -> Na 4s (mz), block_n=2, l=1 -> Na 4p (mz if available).
    na_valence = get_physical_valence_orbitals(11)
    mz_dict = {(orb.n_orbital, orb.l_orbital): orb for orb in na_valence}

    def na_orbital_callable(block_n, l):
        # block_n=1 -> physical_n=3 with n_val_offset=2
        physical_n = block_n + 2
        if (physical_n, l) in mz_dict:
            mz = mz_dict[(physical_n, l)]
            return lambda r, _o=mz: _o.evaluate(np.asarray(r))
        # Fall back to hydrogenic with Z_orb=1 (framework convention)
        return lambda r, _n=block_n, _l=l: hydrogenic_R_nl_analytical(
            1.0, _n, _l, np.asarray(r)
        )

    def h_orbital_callable(block_n, l):
        return lambda r, _n=block_n, _l=l: hydrogenic_R_nl_analytical(
            1.0, _n, _l, np.asarray(r)
        )

    z_Na = 0.0
    z_H = R_EQ
    QUAD = dict(n_rho=80, n_z=100, rho_max=20.0, z_max=20.0)

    # For each Na core orbital, compute its overlap with the bonding orbital
    # = sum_k c_k * <psi_k | psi_core>
    # Each <psi_k | psi_core> is a same-center or cross-center overlap.
    # Restrict to s-s overlaps (l=0 components only); 2p core for l=1 components.

    core_overlaps = {}
    pk_barrier_per_core = {}

    for (n_c, l_c), E_c in NA_CORE_ENERGIES.items():
        if l_c != 0:
            # Skip 2p for now (only s-orbital bonding components; can extend)
            # Use the same axial 2pz machinery as Step 1 above
            continue

        # Build Na core orbital R_c(r) (hydrogenic with Z_eff = n*zeta)
        zeta = NA_CR_EXPONENTS[(n_c, l_c)]
        Z_eff = n_c * zeta

        def R_core(r, _Z=Z_eff, _n=n_c, _l=l_c):
            return hydrogenic_R_nl_analytical(_Z, _n, _l, np.asarray(r))

        # Loop over bonding orbital components
        S_bond_core = 0.0
        for k, lab in enumerate(orbital_labels):
            ck = bond_eigvec[k]
            if abs(ck) < 1e-6:
                continue

            # Parse label: e.g. "Na(1,0,0)" or "H(2,1,-1)"
            side, qns = lab[:2 if lab.startswith('Na') else 1], lab.split('(')[1].rstrip(')')
            qn_parts = [int(x) for x in qns.split(',')]
            block_n, l_k, m_k = qn_parts

            # Skip non-s components (won't overlap with s-core by m-orthogonality)
            if l_k != 0 or m_k != 0:
                continue

            if side == 'Na':
                R_k = na_orbital_callable(block_n, l_k)
                S_k_core = overlap_ss_axial(
                    R_k, z_Na, R_core, z_Na, **QUAD
                )
            elif side == 'H':
                R_k = h_orbital_callable(block_n, l_k)
                S_k_core = overlap_ss_axial(
                    R_k, z_H, R_core, z_Na, **QUAD
                )
            else:
                continue

            S_bond_core += ck * S_k_core

        dH_pk = (0.0 - E_c) * S_bond_core ** 2  # absolute PK
        pk_barrier_per_core[f"({n_c},{l_c})"] = {
            "S_bond_core": S_bond_core,
            "abs_S_pct": 100.0 * abs(S_bond_core),
            "E_c_Ha": E_c,
            "dH_pk_abs": dH_pk,
        }
        print(
            f"  ({n_c},{l_c}) core (E_c = {E_c:+.4f} Ha): "
            f"S_bond_core = {S_bond_core:+.4f} (|S|={100*abs(S_bond_core):.2f}%) "
            f"dH^PK_abs = {dH_pk:+.4f} Ha"
        )
        core_overlaps[f"({n_c},{l_c})"] = S_bond_core

    abs_pk_total = sum(v["dH_pk_abs"] for v in pk_barrier_per_core.values())
    max_abs_pct = max(v["abs_S_pct"] for v in pk_barrier_per_core.values())
    print(f"\nTotal PK barrier (Ev=0, s-cores only): {abs_pk_total:+.4f} Ha")
    print(f"Max |overlap| (s-cores): {max_abs_pct:.2f}%")

    t1 = time.perf_counter()
    print(f"\nWall time: {t1 - t0:.2f} s")

    # Save
    out = {
        "sprint": "F4 Step 1 EXTENDED — full F3 h1 diag",
        "date": "2026-05-23",
        "R_AB_bohr": R_EQ,
        "M": M,
        "h1_lowest_5_eigvalues": w[:5].tolist(),
        "bonding_orbital_idx": bond_idx,
        "bonding_orbital_E_Ha": float(w[bond_idx]),
        "bonding_orbital_composition": {
            orbital_labels[k]: float(bond_eigvec[k]) for k in range(M)
            if abs(bond_eigvec[k]) > 0.05
        },
        "core_overlaps_s_cores": pk_barrier_per_core,
        "total_pk_barrier_s_cores_Ha": abs_pk_total,
        "max_abs_overlap_pct": max_abs_pct,
        "wall_time_s": t1 - t0,
    }
    out_path = os.path.join(_HERE, "data", "sprint_f4_step1_full_h1_diag.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
