"""Josh's occupation-entropy hypothesis: does node-occupation entropy track the
which-site entanglement, across the 5-system frozen/responsive split?

Josh's word-picture (2026-05-30): graph nodes = electron occupation slots
(occupied/unoccupied). The ENTROPY of the occupation pattern should correlate
with the which-site entanglement S_mode. Reframed as the SECOND confinement
coordinate: a closed subshell = Fock-space analog of a closed manifold = one
configuration = low entropy = "confined" (NaH single pair); a partially-filled
shell = many configurations = high entropy = "open" (MgH2 two pairs). The
periodic table IS occupation-confinement.

This is the standard SINGLE-ORBITAL ENTANGLEMENT ENTROPY (Legeza/Reiher QI-of-
electronic-structure): per spatial orbital p, occupation reduced state has 4
probabilities {empty, up, dn, double}; s_p = -sum w ln w; S_occ_nodes = sum_p s_p.

PRE-REGISTERED prediction (locked before running):
  frozen single-pair (H2,NaH,KH): LOW S_occ_nodes
  responsive two-pair (LiH,MgH2): HIGH S_occ_nodes
  and S_occ_nodes should TRACK S_mode (positive correlation across systems).
FAILURE MODE: S_occ_nodes flat across all 5, OR anti-tracks S_mode.

GUARD (honest): the RAW filling fraction N/M is already ruled out
(NaH 2/20 = MgH2 4/40 = 0.1, opposite behavior). So we report filling fraction
too, to show it does NOT discriminate — only the configurational (single-orbital)
entropy should. If S_occ_nodes also fails to discriminate, hypothesis dies.

CAVEAT baked in: correlation entropy and entanglement are both multireference
measures, so positive correlation is PARTLY EXPECTED. The test confirms Josh's
geometric reframing is CONSISTENT and gives the occupation-entropy handle; it is
not an independent surprise. Value = the reframe (2nd confinement coordinate) +
the periodic-table reconnection, not a shocking correlation.

Architecture identical to discriminator. FCI bit-exact validated.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import debug.species_ii_discriminator as disc
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import nah_spec, kh_spec, lih_spec, mgh2_spec


def single_orbital_entropy(vec, a_str, b_str, nb, M):
    """Sum of per-orbital (single-node) occupation entropies, nats.

    For each spatial orbital p, occupation reduced state probs over
    {empty, up, dn, double} from the ground-state determinant amplitudes.
    s_p = -sum w ln w;  S_occ_nodes = sum_p s_p.  Also returns the mean
    natural-style occupation n_p = w_up+w_dn+2*w_double per orbital.
    """
    w = np.zeros((M, 4))  # columns: empty, up, dn, double
    for ai, al in enumerate(a_str):
        aset = set(al)
        for bi, be in enumerate(b_str):
            p2 = abs(vec[ai * nb + bi]) ** 2
            if p2 < 1e-15:
                continue
            bset = set(be)
            for p in range(M):
                ina = p in aset
                inb = p in bset
                if ina and inb:
                    w[p, 3] += p2
                elif ina:
                    w[p, 1] += p2
                elif inb:
                    w[p, 2] += p2
                else:
                    w[p, 0] += p2
    s_tot = 0.0
    occ = np.zeros(M)
    for p in range(M):
        wp = w[p]
        wp = wp[wp > 1e-15]
        s_tot += float(-np.sum(wp * np.log(wp)))
        occ[p] = w[p, 1] + w[p, 2] + 2 * w[p, 3]
    return s_tot, occ, w


def h2_build(R):
    from geovac.molecular_spec import MolecularSpec, OrbitalBlock
    spec = MolecularSpec(
        name='H2',
        blocks=[OrbitalBlock(label='H2_bond', block_type='bond', Z_center=1.0,
                             n_electrons=2, max_n=2, has_h_partner=True,
                             Z_partner=1.0, max_n_partner=2, Z_nuc_center=1.0)],
        nuclear_repulsion_constant=1.0 / R,
        description='H2 2e no core')
    nuclei = [{"Z": 1.0, "position": (0.0, 0.0, 0.0), "label": "Ha"},
              {"Z": 1.0, "position": (0.0, 0.0, float(R)), "label": "Hb"}]
    res = build_balanced_hamiltonian(spec, R=R, nuclei=nuclei,
                                     screened_cross_center=True,
                                     multi_zeta_basis=False, cross_block_h1=True)
    res['nuclear_repulsion'] = 1.0 / R
    return res


def measure(name, builder, n_el, R):
    if name == 'H2':
        result = h2_build(R)
    else:
        result = build_balanced_hamiltonian(builder(), R=R,
                                            screened_cross_center=True,
                                            multi_zeta_basis=False, cross_block_h1=True)
    M = result['M']; site = disc.site_of_spatial(result)
    vals, vecs, a_str, b_str, nb, MM = disc.build_fci(result, n_el, k=8)
    d = int(np.sum(vals - vals[0] < 1e-6))
    manifold = [vecs[:, i] for i in range(d)]
    S_mode, rank = disc.ensemble_mode_entropy(manifold, a_str, b_str, nb, site, M)
    # single-orbital entropy on the (representative) lowest eigenvector
    S_occ_nodes, occ, w = single_orbital_entropy(vecs[:, 0], a_str, b_str, nb, M)
    filling = n_el / (2 * M)  # electrons / spin-orbitals
    # number of "partially occupied" orbitals (0.05 < n_p < 1.95) = active freedom
    n_active = int(np.sum((occ > 0.05) & (occ < 1.95)))
    return {'molecule': name, 'R_bohr': float(R), 'n_electrons': n_el,
            'M_spatial': M, 'filling_fraction': filling,
            'S_mode': float(S_mode), 'S_occ_nodes': float(S_occ_nodes),
            'n_partially_occupied_orbitals': n_active,
            'degeneracy': d}


def main():
    R = 3.566  # near a common equilibrium-ish geometry
    print("=" * 72)
    print("JOSH'S OCCUPATION-ENTROPY HYPOTHESIS (single-orbital entropy)")
    print(f"All systems at R={R} bohr, common architecture")
    print("Predict: S_occ_nodes LOW for 1-pair (H2,NaH,KH), HIGH for 2-pair (LiH,MgH2)")
    print("Guard: filling fraction should NOT discriminate (NaH=MgH2=0.1)")
    print("=" * 72)
    systems = [('H2', None, 2), ('NaH', nah_spec, 2), ('KH', kh_spec, 2),
               ('LiH', lih_spec, 4), ('MgH2', mgh2_spec, 4)]
    rows = []
    for name, builder, n_el in systems:
        r = measure(name, builder, n_el, R)
        rows.append(r)
        print(f"  {name:5} {n_el}e M={r['M_spatial']:2d}  fill={r['filling_fraction']:.3f}  "
              f"S_mode={r['S_mode']:.4f}  S_occ_nodes={r['S_occ_nodes']:.4f}  "
              f"n_partial_orb={r['n_partially_occupied_orbitals']}", flush=True)

    print("\n" + "=" * 72)
    print(f"{'system':6} {'pairs':>5} {'fill':>6} {'S_mode':>8} {'S_occ_nodes':>12} {'n_partial':>10}")
    pairs = {'H2': 1, 'NaH': 1, 'KH': 1, 'LiH': 2, 'MgH2': 2}
    for r in rows:
        print(f"{r['molecule']:6} {pairs[r['molecule']]:>5} {r['filling_fraction']:>6.3f} "
              f"{r['S_mode']:>8.4f} {r['S_occ_nodes']:>12.4f} "
              f"{r['n_partially_occupied_orbitals']:>10}")

    # correlation between S_mode and S_occ_nodes
    sm = np.array([r['S_mode'] for r in rows])
    so = np.array([r['S_occ_nodes'] for r in rows])
    fil = np.array([r['filling_fraction'] for r in rows])
    pr = float(np.corrcoef(sm, so)[0, 1])
    pr_fill = float(np.corrcoef(sm, fil)[0, 1])
    # frozen (1-pair) vs responsive (2-pair) separation in S_occ_nodes
    frozen_so = [r['S_occ_nodes'] for r in rows if pairs[r['molecule']] == 1]
    resp_so = [r['S_occ_nodes'] for r in rows if pairs[r['molecule']] == 2]
    print(f"\ncorr(S_mode, S_occ_nodes) = {pr:+.3f}")
    print(f"corr(S_mode, filling)     = {pr_fill:+.3f}  (guard: should be ~0/weak)")
    print(f"S_occ_nodes: 1-pair = {[round(x,3) for x in frozen_so]}  "
          f"2-pair = {[round(x,3) for x in resp_so]}")
    sep = min(resp_so) - max(frozen_so)
    print(f"separation (min 2-pair - max 1-pair) = {sep:+.4f}  "
          f"{'CLEAN SPLIT' if sep > 0 else 'OVERLAP'}")

    out = {'R_bohr': R, 'rows': rows,
           'corr_Smode_Soccnodes': pr, 'corr_Smode_filling': pr_fill,
           'frozen_Socc': frozen_so, 'responsive_Socc': resp_so,
           'separation': sep,
           '_prediction': 'S_occ_nodes LOW 1-pair, HIGH 2-pair, tracks S_mode',
           '_guard': 'filling fraction should NOT discriminate'}
    Path('debug/data/species_ii_occupation_entropy.json').write_text(json.dumps(out, indent=2))
    print("\nWrote debug/data/species_ii_occupation_entropy.json")


if __name__ == '__main__':
    main()
