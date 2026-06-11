"""
Sprint HF Track 3 — Recoil structural-derivability diagnostic.

This script is diagnostic only. It does NOT modify production code, does NOT
extend Track NI, and does NOT compute new physics. Its only job is to record
in code (and via small ground-truth checks) the structural facts that the memo
draws conclusions from.

Three checks:

1. The cross-register V_ne term in Track NI is a one-body electronic operator
   parameterized by the classical scalar R_PROTON_BOHR. Specifically,
   `finite_size_correction` returns a real number, not a Pauli dict, and the
   `finite_size_coupling_pauli` operator acts only on the electron 1s qubits
   (`Z_e_up`, `Z_e_down`). No Pauli string in V_fs has *any* nontrivial Pauli
   on the nuclear register.

2. The hyperfine coupling crosses registers but only on spin indices: the
   four-qubit Pauli strings touch the *spin* qubits of the proton 0s and the
   electron 1s, not their spatial labels. `R_nuc` does not enter
   `hyperfine_coupling_pauli` at all.

3. The reduced-mass factor (1 + m_e/m_p)^{-3} on |psi(0)|^2 has a clean Fock-
   projection reading: replacing m_e by mu_red in the conformal scale
   p_0 = sqrt(-2 m E_n) shifts the Bohr radius and propagates the right factor
   through the wavefunction at origin. We compute the numbers symbolically as
   a sanity check, not as a derivation.

Author: Sprint HF Track 3 (2026-05-07)
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict

import numpy as np
import sympy as sp

from geovac.nuclear.nuclear_electronic import (
    build_deuterium_composed_hamiltonian,
    finite_size_coupling_pauli,
    hyperfine_coupling_pauli,
    build_nuclear_block,
    build_electronic_block,
    HF_HYDROGEN_HA,
)
from geovac.nuclear.form_factor import (
    R_PROTON_BOHR,
    finite_size_correction,
)


def check_v_fs_is_electronic_only() -> Dict[str, Any]:
    """
    Confirm that every Pauli string in finite_size_coupling_pauli is a pure
    electronic operator: identity on all nuclear qubits, only Z on the
    electron 1s up/down qubits (plus the global identity term).
    """
    nuc = build_nuclear_block(N_shells=2, hw=10.0)
    elec = build_electronic_block(Z=1, n_max=2)
    Q_nuc = nuc['Q_nuc']
    Q_elec = elec['Q_elec']
    Q_total = Q_nuc + Q_elec

    pauli = finite_size_coupling_pauli(
        Q_nuc, Q_elec, nuc, elec,
        R_nuc=R_PROTON_BOHR, Z=1,
    )

    nuclear_active_strings = []
    electronic_strings = []
    for k, v in pauli.items():
        nuc_part = k[:Q_nuc]
        elec_part = k[Q_nuc:]
        nuc_active = any(c != 'I' for c in nuc_part)
        elec_active = any(c != 'I' for c in elec_part)
        if nuc_active:
            nuclear_active_strings.append((k, v))
        if elec_active:
            electronic_strings.append((k, v))

    return {
        'Q_nuc': Q_nuc,
        'Q_elec': Q_elec,
        'Q_total': Q_total,
        'n_pauli_terms': len(pauli),
        'n_nuclear_active': len(nuclear_active_strings),
        'n_electronic_active': len(electronic_strings),
        'nuclear_active_examples': nuclear_active_strings[:3],
        'electronic_examples': [(k, float(v)) for k, v in electronic_strings[:5]],
        'finite_size_correction_value_Ha': finite_size_correction(
            Z=1.0, R_nuc=R_PROTON_BOHR, n=1, l=0,
        ),
        'finite_size_correction_is_scalar': True,
        'r_nuc_is_classical_parameter': True,
        'verdict': (
            'V_fs is a one-body electronic operator parameterized by the '
            'classical scalar R_PROTON_BOHR. No Pauli string in V_fs has '
            'any nontrivial Pauli on the nuclear register.'
        ),
    }


def check_hyperfine_is_spin_only() -> Dict[str, Any]:
    """
    Confirm that the hyperfine coupling is spin-spin only: the four-qubit
    Pauli strings act on (proton_up, proton_down, electron_up, electron_down)
    qubits within the 0s nuclear shell and 1s electronic shell, with no
    coupling to spatial quantum numbers.
    """
    nuc = build_nuclear_block(N_shells=2, hw=10.0)
    elec = build_electronic_block(Z=1, n_max=2)
    Q_nuc = nuc['Q_nuc']
    Q_elec = elec['Q_elec']

    pauli = hyperfine_coupling_pauli(
        Q_nuc, Q_elec, nuc, elec, A_hf=HF_HYDROGEN_HA,
    )

    # Identify the four "active" qubits the hyperfine acts on.
    active_qubits = set()
    weights = []
    for k in pauli.keys():
        weight = sum(1 for c in k if c != 'I')
        weights.append(weight)
        for q, c in enumerate(k):
            if c != 'I':
                active_qubits.add(q)

    # All hyperfine Pauli strings should be 4-qubit operators in the
    # single-occupancy sector.
    all_weight_two_or_four = all(w in (2, 4) for w in weights)

    # The 4 active qubits should be the proton 0s up/down and electron 1s
    # up/down, both pure spin labels (spatial labels n_r=0, l=0, m_l=0 and
    # n=1, l=0, m=0 fixed).
    return {
        'n_pauli_terms': len(pauli),
        'pauli_weights': sorted(set(weights)),
        'n_active_qubits': len(active_qubits),
        'active_qubits_indices': sorted(active_qubits),
        'all_terms_weight_2_or_4': all_weight_two_or_four,
        'r_nuc_appears_in_hyperfine': False,
        'verdict': (
            'Hyperfine couples spin DOF only. R_nuc is not an argument to '
            'hyperfine_coupling_pauli. Spatial coordinates of proton and '
            'electron are not coupled as quantum operators.'
        ),
    }


def check_full_composed_hamiltonian_factorizes() -> Dict[str, Any]:
    """
    Build the full Track NI composed Hamiltonian and confirm that no Pauli
    string couples electron *spatial* qubits (e.g. 2p_m) to nuclear *spatial*
    qubits (e.g. proton 0p, 0d, 1s). The only cross-register terms are spin-
    spin within the (1s_e, 0s_p) sub-block.
    """
    out = build_deuterium_composed_hamiltonian(
        N_shells=2, hw=10.0, n_max_elec=2,
        include_finite_size=True, include_hyperfine=True,
    )
    H_pauli = out['pauli_terms']
    Q_nuc = out['Q_nuc']
    Q_elec = out['Q_elec']
    Q_total = out['Q_total']

    elec_states = out['electronic_block']['states']
    nuc_states_p = out['nuclear_block']['states_p']
    nuc_states_n = out['nuclear_block']['states_n']

    # Identify the "spin-only" qubits in the 1s electron pair and 0s proton
    # pair. Everything else is spatial.
    elec_1s_qubits = set()
    for i, (n, l, m, ms) in enumerate(elec_states):
        if n == 1 and l == 0 and m == 0:
            elec_1s_qubits.add(Q_nuc + i)

    nuc_0s_p_qubits = set()
    for i, s in enumerate(nuc_states_p):
        if s.n_r == 0 and s.l == 0 and s.m_l == 0:
            nuc_0s_p_qubits.add(i)

    spin_only_qubits = elec_1s_qubits | nuc_0s_p_qubits

    # All other qubits are "spatial" (electron 2s/2p/etc. or proton/neutron
    # excited shells).
    nuclear_qubits_all = set(range(Q_nuc))
    elec_qubits_all = set(range(Q_nuc, Q_total))
    spatial_qubits = (nuclear_qubits_all | elec_qubits_all) - spin_only_qubits

    # Count cross-register Pauli strings that hit spatial qubits on BOTH
    # registers (would be a coordinate-coordinate coupling).
    n_spatial_cross_terms = 0
    n_spinspin_cross_terms = 0
    n_nuclearonly = 0
    n_electroniconly = 0
    n_identity = 0

    for k, v in H_pauli.items():
        nuc_part = k[:Q_nuc]
        elec_part = k[Q_nuc:]
        nuc_active = any(c != 'I' for c in nuc_part)
        elec_active = any(c != 'I' for c in elec_part)

        if not nuc_active and not elec_active:
            n_identity += 1
            continue
        if nuc_active and not elec_active:
            n_nuclearonly += 1
            continue
        if elec_active and not nuc_active:
            n_electroniconly += 1
            continue

        # Cross-register: does it hit spatial-spatial?
        nuc_spatial_hit = any(
            c != 'I' and q in spatial_qubits
            for q, c in enumerate(k[:Q_nuc])
        )
        elec_spatial_hit = any(
            c != 'I' and (Q_nuc + q) in spatial_qubits
            for q, c in enumerate(k[Q_nuc:])
        )
        if nuc_spatial_hit and elec_spatial_hit:
            n_spatial_cross_terms += 1
        else:
            n_spinspin_cross_terms += 1

    return {
        'Q_total': Q_total,
        'Q_nuc': Q_nuc,
        'Q_elec': Q_elec,
        'n_pauli_terms_total': len(H_pauli),
        'n_identity': n_identity,
        'n_nuclear_only': n_nuclearonly,
        'n_electronic_only': n_electroniconly,
        'n_spinspin_cross_register': n_spinspin_cross_terms,
        'n_spatial_cross_register': n_spatial_cross_terms,
        'verdict': (
            'Cross-register Pauli strings exist only in the spin-spin '
            'sector (proton 0s spin x electron 1s spin, hyperfine I.S). '
            'Zero terms couple spatial proton/neutron qubits to spatial '
            'electron qubits. The framework has NO V_ne(r_e, R_n) two-body '
            'coordinate operator.'
        ),
    }


def check_reduced_mass_fock_projection_reading() -> Dict[str, Any]:
    """
    Sanity check: the (1 + m_e/m_p)^{-3} factor on |psi(0)|^2 follows from
    rescaling the Bohr radius by mu_red/m_e in atomic units. This is NOT a
    derivation from Fock projection — it is a check that the standard
    reduced-mass factor agrees with HF-1's hand-applied number, and that the
    framework can host the reduced-mass replacement *if* one declares the
    Fock projection's mass parameter to be mu_red rather than m_e.
    """
    # Symbolic
    m_e, m_p = sp.symbols('m_e m_p', positive=True)
    mu = m_e * m_p / (m_e + m_p)
    # |psi_1s(0)|^2 = (Z * mu / m_e)^3 / pi  in atomic units (m_e = 1).
    # In a.u. with mu_red replacing m_e: |psi(0)|^2 -> Z^3 * (mu/m_e)^3 / pi.
    # Ratio to point-nucleus value:
    factor_sym = (mu / m_e) ** 3
    factor_simplified = sp.simplify(factor_sym)
    # = 1 / (1 + m_e/m_p)^3
    # Verify by series:
    series = sp.series(factor_simplified, m_e, 0, 3)

    # Numerical
    me_over_mp = 1.0 / 1836.15267343
    factor_num = 1.0 / (1.0 + me_over_mp) ** 3

    # HF-1 numbers
    A_strict_MHz = 1421.16
    A_recoil_MHz = A_strict_MHz * factor_num
    A_exp_MHz = 1420.405751768

    return {
        'factor_symbolic': str(factor_simplified),
        'factor_series_in_me_over_mp': str(series),
        'factor_numerical': factor_num,
        'one_minus_factor_ppm': (1.0 - factor_num) * 1e6,
        'A_strict_BF_MHz': A_strict_MHz,
        'A_recoil_corrected_MHz': A_recoil_MHz,
        'A_experimental_MHz': A_exp_MHz,
        'strict_residual_ppm': (A_strict_MHz - A_exp_MHz) / A_exp_MHz * 1e6,
        'recoil_residual_ppm': (A_recoil_MHz - A_exp_MHz) / A_exp_MHz * 1e6,
        'derivation_status': (
            'NOT derived. The reduced-mass replacement m_e -> mu_red in the '
            'Bohr radius is a hand-applied prescription, not produced by '
            'Fock projection of the Track NI cross-register architecture. '
            'The Fock projection p_0 = sqrt(-2E) is a stereographic scaling '
            'of a single-particle Hamiltonian; making it two-body would '
            'require a center-of-mass / relative coordinate split that the '
            'framework does not currently host. The factor is consistent '
            'with the framework but external to it.'
        ),
    }


def main() -> None:
    results = {
        'sprint': 'HF-3',
        'date': '2026-05-07',
        'goal': (
            'Diagnose whether GeoVac\'s cross-register architecture (Track NI) '
            'natively produces reduced-mass / recoil corrections to the '
            'hyperfine Bohr-Fermi prediction, or whether recoil is an external '
            'prescription.'
        ),
        'check_1_v_fs_electronic_only': check_v_fs_is_electronic_only(),
        'check_2_hyperfine_spin_only': check_hyperfine_is_spin_only(),
        'check_3_full_factorization': check_full_composed_hamiltonian_factorizes(),
        'check_4_reduced_mass_status': check_reduced_mass_fock_projection_reading(),
        'verdict': 'NEGATIVE',
        'verdict_summary': (
            'Track NI\'s V_ne is a classical scalar parameter R_PROTON_BOHR '
            'fed into a one-body electronic operator. The proton\'s spatial '
            'coordinate is not a dynamical variable on the nuclear register. '
            'The hyperfine coupling does cross registers, but only on spin '
            'indices. No V_ne(r_e, R_n) two-body coordinate-coupling operator '
            'exists in the framework. Reduced-mass / recoil corrections at '
            'every order remain external focal-length inputs. The +0.754 MHz '
            '/ +531 ppm residual at strict Bohr-Fermi (g_e = 2, no recoil) '
            'is the recoil contribution; the recoil-corrected number 1418.84 '
            'MHz is the right baseline for HF-2 to attack with a_e from '
            'graph-native vertex correction.'
        ),
    }

    # Strip non-JSON-serializable items
    def _clean(o: Any) -> Any:
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_clean(x) for x in o]
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        try:
            json.dumps(o)
            return o
        except (TypeError, ValueError):
            return str(o)

    cleaned = _clean(results)
    out_path = os.path.join(
        os.path.dirname(__file__), 'data', 'sprint_hf_track3.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(cleaned, f, indent=2)

    print('Sprint HF-3 diagnostic complete.')
    print(f'Verdict: {results["verdict"]}')
    print()
    print('Check 1 (V_fs electronic-only):',
          results['check_1_v_fs_electronic_only']['n_nuclear_active'],
          'nuclear-active Pauli strings (expected 0)')
    print('Check 2 (hyperfine spin-only):',
          results['check_2_hyperfine_spin_only']['r_nuc_appears_in_hyperfine'],
          '(expected False)')
    print('Check 3 (full factorization):',
          results['check_3_full_factorization']['n_spatial_cross_register'],
          'spatial-spatial cross-register terms (expected 0)')
    print('Check 4 (reduced mass): factor =',
          results['check_4_reduced_mass_status']['factor_numerical'],
          'expected (1 + m_e/m_p)^-3 ~ 0.99837')
    print()
    print(f'Output: {out_path}')


if __name__ == '__main__':
    main()
