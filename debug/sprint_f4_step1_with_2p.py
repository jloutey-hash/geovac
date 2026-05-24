"""Sprint F4 Step 1 EXTENDED — add Na 2p_z core overlap.

The full F3 h1 bonding orbital from `sprint_f4_step1_full_h1_diag.py`
has c_Na_3s = -0.62, c_Na_4s = +0.35, c_H_1s = -0.47, c_H_2s = +0.52.
None of these are p-orbital components, so the Na 2p_z core overlap
involves only the 2pz component of each Na valence orbital expressed
back via spherical harmonics — which for s-character orbitals is zero.

But the bonding orbital ALSO contains polarization through h1 mixing of
the Na 2p_z (block_n=2, l=1, m=0) and H 2p_z components if they couple.
Let me re-examine the full eigenvector to see if any 2p_0 component is
present.

This script just enriches the data; it does not change the verdict.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.linalg import eigh

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from geovac.balanced_coupled import build_balanced_hamiltonian, _get_block_geometry  # noqa: E402
from geovac.molecular_spec import nah_spec  # noqa: E402
from geovac.cross_block_h1 import hydrogenic_R_nl_analytical, overlap_ss_axial  # noqa: E402
from geovac.multi_zeta_orbitals import get_physical_valence_orbitals  # noqa: E402


def main():
    R_EQ = 3.566
    spec = nah_spec(max_n=2)
    nuclei = [
        {'Z': 11.0, 'position': (0.0, 0.0, 0.0), 'label': 'Na'},
        {'Z': 1.0, 'position': (0.0, 0.0, R_EQ), 'label': 'H'},
    ]
    result = build_balanced_hamiltonian(
        spec=spec, R=R_EQ, nuclei=nuclei,
        screened_cross_center=True, multi_zeta_basis=True, cross_block_h1=True,
        verbose=False,
    )
    h1 = result['h1']
    M = result['M']
    sub_blocks = _get_block_geometry(spec)
    orbital_labels = []
    for sb in sub_blocks:
        for (n, l, m) in sb['states']:
            orbital_labels.append((sb['label'], n, l, m))

    w, v = eigh(h1)
    bond_idx = 0  # lowest eigenvector (from previous diagnostic)
    bond_eigvec = v[:, bond_idx]
    E_bond = w[bond_idx]

    print(f"Full bonding orbital composition (|c| > 0.01):")
    for k, (lab, n, l, m) in enumerate(orbital_labels):
        if abs(bond_eigvec[k]) > 0.01:
            print(f"  [{k}] {lab}({n},{l},{m}): {bond_eigvec[k]:+.5f}")

    # Compute Na 2p_z overlap if any l=1, m=0 components on either side
    NA_CR_2P_ZETA = 3.45200
    NA_CORE_E_2P = -1.5181
    Z_eff_2p = 2 * NA_CR_2P_ZETA

    z_Na = 0.0
    z_H = R_EQ
    QUAD = dict(n_rho=100, n_z=120, rho_max=20.0, z_max=20.0)

    # For Na 2p_z core <-> orbital k, the overlap is nonzero only if orbital k
    # has m=0 azimuthal symmetry. Specifically:
    # - Na (block_n, l=1, m=0) maps to a 4p_z-like orbital at Na -> overlap
    #   with Na 2p_z core
    # - H (block_n, l=1, m=0) is a 1s with l=1 - this is the H "p_z" effective
    #   in the framework basis -> cross-center 2p_z overlap

    # Build axial 2D quadrature
    from scipy.special import roots_legendre
    u, wu = roots_legendre(QUAD['n_rho'])
    rho_pts = QUAD['rho_max'] * (u + 1.0) / 2.0
    rho_w = wu * QUAD['rho_max'] / 2.0
    v_z, wv_z = roots_legendre(QUAD['n_z'])
    z_pts = 0.5 * (z_Na + z_H) + QUAD['z_max'] * v_z
    z_w = wv_z * QUAD['z_max']
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')

    # Na 2p_z core: R_2p(r_Na) * cos(theta_Na) * sqrt(3/(4 pi))
    r_Na = np.sqrt(rho_grid ** 2 + (z_grid - z_Na) ** 2)
    R_2p_core = hydrogenic_R_nl_analytical(
        Z_eff_2p, 2, 1, r_Na.flatten()
    ).reshape(r_Na.shape)
    cos_theta_Na = np.where(r_Na > 1e-12, (z_grid - z_Na) / r_Na, 0.0)
    psi_2pz_core = R_2p_core * cos_theta_Na * np.sqrt(3.0 / (4.0 * np.pi))

    # Loop over all valence orbital components and accumulate
    # <bond | Na 2p_z>
    na_valence = get_physical_valence_orbitals(11)
    mz_dict = {(o.n_orbital, o.l_orbital): o for o in na_valence}

    def na_p_callable(block_n, l):
        physical_n = block_n + 2
        if (physical_n, l) in mz_dict:
            mz = mz_dict[(physical_n, l)]
            return lambda r, _o=mz: _o.evaluate(np.asarray(r))
        return lambda r, _n=block_n, _l=l: hydrogenic_R_nl_analytical(
            1.0, _n, _l, np.asarray(r)
        )

    S_bond_2pz = 0.0
    for k, (lab, n, l, m) in enumerate(orbital_labels):
        ck = bond_eigvec[k]
        if abs(ck) < 1e-6:
            continue
        if m != 0:
            continue  # only m=0 components couple to 2pz (axial)
        if l == 0:
            continue  # s-orbitals don't overlap with 2pz (different parity in cos(theta))

        side = lab.split('_')[-1]  # 'center' or 'partner'
        if side == 'center':
            # Na side l=1 m=0 orbital
            R_k = na_p_callable(n, l)
            r_k = r_Na  # centered at Na
            cos_theta_k = cos_theta_Na
        else:
            # H side l=1 m=0 orbital (hydrogenic 2p_z at z=R_EQ)
            R_k = lambda r: hydrogenic_R_nl_analytical(1.0, n, l, np.asarray(r))
            r_k = np.sqrt(rho_grid ** 2 + (z_grid - z_H) ** 2)
            cos_theta_k = np.where(r_k > 1e-12, (z_grid - z_H) / r_k, 0.0)

        R_k_val = R_k(r_k.flatten()).reshape(r_k.shape)
        # Y_{1,0} normalization for the orbital itself
        psi_k = R_k_val * cos_theta_k * np.sqrt(3.0 / (4.0 * np.pi))

        # Integral: 2 pi * sum_grid (psi_k * psi_2pz_core * rho)
        integrand = psi_k * psi_2pz_core * rho_grid
        weight_2d = np.outer(rho_w, z_w)
        raw = np.sum(weight_2d * integrand)
        S_k_core = 2.0 * np.pi * raw

        contrib = ck * S_k_core
        S_bond_2pz += contrib
        print(f"  Component {lab}({n},{l},{m}) with c={ck:+.4f}: "
              f"<psi_k|2pz>={S_k_core:+.4f} -> contrib={contrib:+.4f}")

    print(f"\nTotal S(bond, Na 2p_z core) = {S_bond_2pz:+.4f} (|S|={100*abs(S_bond_2pz):.2f}%)")
    dH_2pz = (0.0 - NA_CORE_E_2P) * S_bond_2pz ** 2
    print(f"dH^PK_(Ev=0) from 2pz core = {dH_2pz:+.4f} Ha")
    print(f"Note: 2p core has 3 m=-1,0,+1 components; here we only computed m=0.")
    print(f"By azimuthal selection, only m=0 component couples to s-symmetric bonding orbital.")

    # Combine with previous s-core result
    prev = json.load(open(os.path.join(_HERE, "data", "sprint_f4_step1_full_h1_diag.json")))
    prev_total = prev['total_pk_barrier_s_cores_Ha']
    new_total = prev_total + dH_2pz
    print(f"\nUpdated total PK barrier (s + 2p_z): {new_total:+.4f} Ha (was {prev_total:+.4f} from s-cores)")
    print(f"vs F3 wall depth ~4.37 Ha at NaH max_n=2 R=2.0 bohr -> ratio = {new_total / 4.37 * 100:.2f}%")

    out = {
        "sprint": "F4 Step 1 — add 2p_z core overlap",
        "S_bond_2pz_core": S_bond_2pz,
        "dH_pk_2pz_Ha": dH_2pz,
        "total_pk_barrier_s_plus_2pz_Ha": new_total,
        "wall_depth_F3_Ha": 4.37,
        "pk_fraction_of_wall_pct": new_total / 4.37 * 100,
    }
    out_path = os.path.join(_HERE, "data", "sprint_f4_step1_with_2p.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
