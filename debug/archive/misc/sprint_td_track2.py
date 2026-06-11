"""Sprint TD Track 2 — Temperature decoding on atomic thermal observables.

Verifies the user's instinct: thermodynamic S(T) for an atomic system
factors as

    S_thermo(T) = k_B * S_level_mixing(T)
                = k_B * S_info(Boltzmann ensemble)

with the temperature parameter beta living entirely on the (apparatus)
S^1_beta tensor factor of T_{S^3} ⊗ T_{S^1_beta}, and the dimensionless
residual being the information-theoretic entropy of the level-population
distribution. The graph-side (S^3) part is captured at T → 0 by Paper 27's
*spatial* 1-RDM von Neumann entropy of the ground state.

Three atomic systems:
  1. Hydrogen (Z=1, 1 electron) — single-determinant eigenstates,
     Paper 27 S_full(GS) = 0 by single-determinant rigidity. Cleanest case:
     ALL of S_thermo is apparatus, residual is exactly zero.
  2. Helium (Z=2, 2 electrons) — graph-native CI ground state has
     V_ee correlation entropy (Paper 27 EP-2c). Thermal decomposition has
     both apparatus and intra-state correlation residual.
  3. Li+ (Z=3, 2 electrons) — He-like at higher Z. Paper 27 EP-2c
     predicts S_full ~ Z^{-2.6} so residual should drop sharply.

Operational decomposition (pure quantum statistical mechanics, exact):

  rho(beta) = sum_n p_n |n><n|       Gibbs state
  p_n = exp(-beta E_n) / Z(beta)     level probability (per microstate)
  Z(beta) = sum_n exp(-beta E_n)     partition function (microstate sum)

  S_vN(rho) = -sum_n p_n log p_n     # apparatus = level-mixing entropy
            = beta * U + log Z       # standard thermodynamic identity

  rho_1(beta) = Tr_2 rho(beta) = sum_n p_n rho_n^(1)    # spatial 1-RDM

  S_intra(beta) = sum_n p_n S_vN(rho_n^(1))   # weighted Paper-27 entropy

The total ensemble-on-spatial entropy is bounded by

  max_n S_vN(rho_n^(1)) <= S_vN(rho_1)
                        <= S_apparatus + S_intra (Lindblad concavity)

Limits we verify:
  T -> 0: S_apparatus -> 0,         rho_1 -> rho_GS^(1),
          residual = S_paper27(GS).
  T -> inf: S_apparatus -> log(N_microstates),
            rho_1 -> uniform/N_orb (in the truncated basis).

This is conceptual numerical verification, not new physics. The point is
to exhibit the temperature-decoding decomposition operationally on three
real atomic systems and connect the residual to Paper 27.

Author: Sprint TD Track 2 worker fork
Date: 2026-05-08
"""

from __future__ import annotations

import json
import os
import sys
import time
from typing import Dict, List, Tuple

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.archive.misc.energy_entanglement_decoupling import build_decomposed_hamiltonians
from debug.archive.misc.entanglement_geometry import (
    build_1rdm_from_singlet_ci,
    compute_entanglement_measures,
)
from geovac.casimir_ci import build_graph_native_fci


# ---------------------------------------------------------------------------
#  Generic thermal decomposition machinery
# ---------------------------------------------------------------------------

def boltzmann_weights(energies: np.ndarray,
                      degeneracies: np.ndarray,
                      beta: float) -> Tuple[np.ndarray, float, float]:
    """Compute level probabilities, log Z, and U for given (E_n, g_n, beta).

    Uses a stabilized formulation with E_min subtracted so exp(-beta dE)
    stays in [0, 1]. The returned p_n is the probability per LEVEL
    (g_n already absorbed in numerator), summing to 1.

    Returns
    -------
    p_level : array of level probabilities (sums to 1)
    logZ : log of partition function (level sum, NOT microstate sum)
    U : <E> in atomic units
    """
    E_min = float(energies.min())
    dE = energies - E_min
    # Per-level weight: g_n * exp(-beta dE_n)
    log_w = np.log(degeneracies.astype(float)) - beta * dE
    log_w_max = log_w.max()
    w = np.exp(log_w - log_w_max)
    Z_level = w.sum()
    p_level = w / Z_level
    # logZ is the "microstate" partition log Z = log sum g_n exp(-beta E_n)
    logZ = log_w_max + np.log(Z_level) - beta * E_min  # add back E_min shift
    U = float((p_level * energies).sum())
    return p_level, float(logZ), U


def thermal_decomposition(energies: np.ndarray,
                          degeneracies: np.ndarray,
                          beta: float) -> Dict[str, float]:
    """Compute the temperature-decoding decomposition at one beta.

    Returns a dict with:
        beta, T, F (free energy), U, S_thermo, S_microstate_info,
        S_level_info (treats levels as effective states with g_n weight),
        residual_info (= 0 by construction — apparatus = total at this layer).

    The "residual" against Paper-27 1-RDM entropy lives at a different
    layer: it is computed by mixing 1-RDMs across populated states,
    handled separately when state-resolved 1-RDMs are available.
    """
    p_level, logZ, U = boltzmann_weights(energies, degeneracies, beta)
    F = U - (1.0 / beta) * (logZ - 0.0)  # F = -kT log Z; logZ already includes E_min
    F_alt = -(1.0 / beta) * logZ
    S_thermo = beta * (U - F_alt)  # = beta*U + log Z, dimensionless (units of k_B)

    # Apparatus content split, treating each MICROSTATE as distinct:
    # S_micro = -sum_microstates P log P
    #         = -sum_n g_n * (p_level_n / g_n) log(p_level_n / g_n)
    #         = -sum_n p_level_n log(p_level_n / g_n)
    #         = -sum_n p_level_n log p_level_n + sum_n p_level_n log g_n
    p_clip = np.clip(p_level, 1e-300, None)
    S_level_info = float(-(p_level * np.log(p_clip)).sum())  # over levels
    log_g = np.log(degeneracies.astype(float))
    S_microstate_info = S_level_info + float((p_level * log_g).sum())

    # Verification residual: S_thermo (textbook) vs S_microstate_info
    residual = S_thermo - S_microstate_info

    return {
        'beta': float(beta),
        'T_au': 1.0 / float(beta),
        'logZ': float(logZ),
        'F': float(F_alt),
        'U': U,
        'S_thermo': float(S_thermo),
        'S_microstate_info': S_microstate_info,
        'S_level_info': S_level_info,
        'identity_residual': float(residual),  # should be ~ machine eps
        'p_level': p_level.tolist(),
    }


# ---------------------------------------------------------------------------
#  Hydrogen (Z=1, single electron)
# ---------------------------------------------------------------------------

def hydrogen_spectrum(Z: int = 1, n_max: int = 8) -> Tuple[np.ndarray, np.ndarray, List[Tuple]]:
    """Hydrogenic bound-state spectrum E_n = -Z^2/(2 n^2), g_n = n^2.

    Each principal shell contributes degeneracy n^2 (l = 0..n-1, m = -l..l).
    Spin degeneracy is NOT included — it is a constant 2 prefactor that
    multiplies every Z(beta) by 2 and shifts logZ by log(2), but does not
    change S_thermo or any temperature-decoding decomposition. We omit it
    for clarity.
    """
    energies = []
    degeneracies = []
    levels = []
    for n in range(1, n_max + 1):
        E = -float(Z) ** 2 / (2.0 * n * n)
        g = n * n
        energies.append(E)
        degeneracies.append(g)
        levels.append(('H', n, g, E))
    return np.array(energies), np.array(degeneracies), levels


# ---------------------------------------------------------------------------
#  Helium / Li+ (Z >= 2, 2-electron graph-native CI)
# ---------------------------------------------------------------------------

def helike_spectrum_with_states(Z: int, n_max: int = 3,
                                spin_label: str = 'singlet'):
    """Build He-like graph-native FCI matrix and return spectrum + eigenvectors.

    Uses build_decomposed_hamiltonians for compatibility with Paper 27's
    1-RDM machinery (configs and orbital labels match).

    Returns
    -------
    energies : sorted eigenvalues (the level structure within m_total=0,
               singlet sector at this n_max)
    degeneracies : array of 1's by default — graph-native CI in fixed
                   m_total/spin sector does not enumerate explicit
                   degeneracy. We use g_n = 1 per eigenstate.
                   (M_L = 0 sector mixes ¹S and ¹P at M_L=0; for thermal
                   decomposition the labeling does not matter as long as
                   we do not mis-attribute each eigenvalue to a separate
                   physical multiplet. We treat each eigenvalue as one
                   non-degenerate microstate of a graph-truncated basis.)
    eigenvectors : array of shape (n_configs, n_states), each column a
                   normalized CI vector
    configs : list of orbital index pairs (i, j) defining the basis
    n_spatial : int (orbital count)
    """
    decomp = build_decomposed_hamiltonians(float(Z), n_max)
    H = decomp['H_full']
    configs = decomp['configs']
    n_spatial = decomp['n_spatial']

    evals, evecs = np.linalg.eigh(H)
    return evals, np.ones(len(evals), dtype=int), evecs, configs, n_spatial


def helike_state_resolved_entropies(evals: np.ndarray,
                                    evecs: np.ndarray,
                                    configs,
                                    n_spatial: int,
                                    n_states_resolve: int = 5) -> List[float]:
    """Compute Paper 27 S_full = von Neumann entropy of spatial 1-RDM
    for the lowest `n_states_resolve` eigenstates of a He-like FCI block.

    Returns a list (length = min(n_states_resolve, len(evals))) of S_vN values.
    """
    n_take = min(n_states_resolve, evecs.shape[1])
    S_list = []
    for k in range(n_take):
        ci = evecs[:, k]
        rho = build_1rdm_from_singlet_ci(ci, configs, n_spatial)
        ent = compute_entanglement_measures(rho)
        S_list.append(float(ent['von_neumann_entropy']))
    return S_list


def helike_thermal_residual(evals: np.ndarray,
                            evecs: np.ndarray,
                            configs,
                            n_spatial: int,
                            beta: float,
                            n_states_used: int = 30) -> Dict[str, float]:
    """Compute thermal decomposition AND the spatial 1-RDM thermal mix.

    The spatial-1-RDM thermal-mix entropy is the ENSEMBLE entropy:

        rho_1(beta) = sum_n p_n * rho_n^{(1)}

        S_thermal_mix = S_vN(rho_1(beta))
                     = - Tr rho_1 log rho_1

    This is the "decoded-temperature" residual the user is asking about.

    Returns a dict combining the level-mixing decomposition (via
    thermal_decomposition) with the state-resolved residual.
    """
    n_take = min(n_states_used, evecs.shape[1])
    E_use = evals[:n_take]
    g_use = np.ones(n_take, dtype=int)

    # Level-mixing decomposition (apparatus piece)
    decomp = thermal_decomposition(E_use, g_use, beta)
    p_level = np.array(decomp['p_level'])

    # Build thermal 1-RDM as p-weighted mix of state 1-RDMs
    rho_1_thermal = np.zeros((n_spatial, n_spatial))
    S_intra_per_state = []
    for k in range(n_take):
        ci_k = evecs[:, k]
        rho_k = build_1rdm_from_singlet_ci(ci_k, configs, n_spatial)
        ent_k = compute_entanglement_measures(rho_k)
        S_intra_per_state.append(float(ent_k['von_neumann_entropy']))
        rho_1_thermal += p_level[k] * rho_k

    # Trace check
    trace_thermal = float(np.trace(rho_1_thermal))

    # Spatial-1-RDM thermal entropy (the user's "decoded temperature" residual)
    ent_thermal = compute_entanglement_measures(rho_1_thermal)
    S_spatial_thermal = float(ent_thermal['von_neumann_entropy'])

    # Concavity bounds:
    #   max_n S_n <= S_vN(rho_1) <= sum_n p_n S_n + S_apparatus
    S_intra_avg = float(np.dot(p_level, np.array(S_intra_per_state)))
    bound_lower = float(np.max(S_intra_per_state[:n_take]))
    bound_upper = S_intra_avg + decomp['S_microstate_info']

    decomp.update({
        'S_intra_per_state': S_intra_per_state,
        'S_intra_thermal_avg': S_intra_avg,
        'S_spatial_1rdm_thermal': S_spatial_thermal,
        'rho_1_thermal_trace': trace_thermal,
        'concavity_lower_bound': bound_lower,
        'concavity_upper_bound': bound_upper,
        'concavity_holds': bool(bound_lower - 1e-10 <= S_spatial_thermal <= bound_upper + 1e-10),
    })
    return decomp


# ---------------------------------------------------------------------------
#  Drivers
# ---------------------------------------------------------------------------

def run_hydrogen(n_max: int = 6, n_T: int = 30) -> Dict:
    """Hydrogen Z=1 thermal decomposition.

    Hydrogen has g_n = n^2 microstates per shell. Single-determinant
    eigenstates have spatial 1-RDM = |phi_nlm><phi_nlm|, a rank-1
    projector — Paper 27 S_full = 0 by single-determinant rigidity for
    this 1-electron pure state. Thermal mixing of pure 1-electron
    eigenstates produces a nonzero S_vN(rho_1) BUT this is purely the
    LEVEL-mixing entropy translated onto the 1-RDM (since each level's
    1-RDM is a rank-1 orthogonal projector in the chosen basis). So
    S_spatial_1rdm_thermal == S_microstate_info for H — by construction,
    when each level contributes orthogonal rank-1 1-RDMs.

    Note: in the 1-electron picture rho_1 IS the state, and it is just a
    classical Gibbs distribution over orthogonal projectors. The
    "temperature decoding" is trivially clean here.
    """
    Z = 1
    energies, degeneracies, levels = hydrogen_spectrum(Z=Z, n_max=n_max)
    E1 = energies[0]
    E2 = energies[1]
    gap = E2 - E1
    # Span beta from k_B T ~ 0.01 * gap to ~ 10 * gap
    T_au = np.logspace(np.log10(0.01 * gap), np.log10(20.0 * gap), n_T)
    betas = 1.0 / T_au

    rows = []
    for beta in betas:
        decomp = thermal_decomposition(energies, degeneracies, beta)
        # For hydrogen, the 1-electron state IS the spatial 1-RDM (no 2nd
        # electron to trace out). Within an atomic shell each (n,l,m)
        # microstate gives an orthogonal eigenstate of the 1e Hamiltonian,
        # so S_vN(rho_1) over a Gibbs mixture is exactly the level-mixing
        # entropy when we count each (n,l,m) microstate as a separate level.
        # That's the S_microstate_info we already computed.
        rows.append({
            'beta_au': decomp['beta'],
            'T_au': decomp['T_au'],
            'F_au': decomp['F'],
            'U_au': decomp['U'],
            'S_thermo_kB': decomp['S_thermo'],
            'S_microstate_info': decomp['S_microstate_info'],
            'S_level_info': decomp['S_level_info'],
            'identity_residual': decomp['identity_residual'],
            'S_paper27_residual': 0.0,    # exact: 1e GS has no V_ee correlation
            'spatial_residual_kB': 0.0,   # Paper 27 floor for 1e single det
        })

    return {
        'system': 'H',
        'Z': Z,
        'n_max_basis': n_max,
        'energies_au': energies.tolist(),
        'degeneracies': degeneracies.tolist(),
        'thermal_table': rows,
        'paper27_GS_S_full': 0.0,
        'GS_E_au': float(energies[0]),
        'first_excitation_au': float(gap),
    }


def run_helike(Z: int, n_max: int = 3, n_T: int = 30,
               n_states_used: int = 30,
               n_states_resolve: int = 5) -> Dict:
    """He-like (Z=2, 3, ...) thermal decomposition with state-resolved residual.

    Uses graph-native FCI in the m_total=0 singlet sector.
    """
    evals, _gees, evecs, configs, n_spatial = helike_spectrum_with_states(
        Z=Z, n_max=n_max
    )
    n_take = min(n_states_used, len(evals))
    evals_use = evals[:n_take]
    evecs_use = evecs[:, :n_take]

    # State-resolved Paper-27 entropies for low-lying eigenstates
    S_low = helike_state_resolved_entropies(
        evals_use, evecs_use, configs, n_spatial,
        n_states_resolve=min(n_states_resolve, n_take),
    )
    S_paper27_GS = float(S_low[0]) if S_low else 0.0

    E0 = float(evals_use[0])
    if len(evals_use) > 1:
        gap = float(evals_use[1] - evals_use[0])
    else:
        gap = 1e-3
    T_au = np.logspace(np.log10(0.005 * gap), np.log10(10.0 * gap), n_T)
    betas = 1.0 / T_au

    rows = []
    t_start = time.time()
    for beta in betas:
        decomp = helike_thermal_residual(
            evals_use, evecs_use, configs, n_spatial,
            beta=float(beta), n_states_used=n_take,
        )
        rows.append({
            'beta_au': decomp['beta'],
            'T_au': decomp['T_au'],
            'F_au': decomp['F'],
            'U_au': decomp['U'],
            'S_thermo_kB': decomp['S_thermo'],
            'S_microstate_info': decomp['S_microstate_info'],
            'identity_residual': decomp['identity_residual'],
            'S_intra_thermal_avg': decomp['S_intra_thermal_avg'],
            'S_spatial_1rdm_thermal': decomp['S_spatial_1rdm_thermal'],
            'rho1_thermal_trace': decomp['rho_1_thermal_trace'],
            'concavity_holds': decomp['concavity_holds'],
            'concavity_lower_bound': decomp['concavity_lower_bound'],
            'concavity_upper_bound': decomp['concavity_upper_bound'],
        })
    print(f"  Z={Z}: spectrum + thermal scan in {time.time() - t_start:.1f}s "
          f"(n_states_used={n_take}, n_T={n_T})")

    return {
        'system': f'He-like Z={Z}',
        'Z': Z,
        'n_max_basis': n_max,
        'spin_sector': 'singlet, m_total=0',
        'n_states_in_sector': len(evals),
        'n_states_used': n_take,
        'energies_au': evals_use.tolist(),
        'first_5_S_paper27': S_low,
        'paper27_GS_S_full': S_paper27_GS,
        'GS_E_au': E0,
        'first_excitation_au': gap,
        'thermal_table': rows,
    }


# ---------------------------------------------------------------------------
#  Top-level orchestration
# ---------------------------------------------------------------------------

def main(out_path: str):
    print("=" * 70)
    print("Sprint TD Track 2: Temperature decoding on atomic thermal observables")
    print("=" * 70)

    print("\n[1/3] Hydrogen (Z=1, 1 electron, n_max=6)")
    H_data = run_hydrogen(n_max=6, n_T=30)
    print(f"     GS energy = {H_data['GS_E_au']:.6f} Ha")
    print(f"     1->2 gap  = {H_data['first_excitation_au']:.6f} Ha")
    print(f"     Paper-27 S_full(GS) = {H_data['paper27_GS_S_full']:.6e}")

    print("\n[2/3] Helium (Z=2, 2 electrons, n_max=3)")
    He_data = run_helike(Z=2, n_max=3, n_T=30,
                         n_states_used=30, n_states_resolve=5)
    print(f"     GS energy = {He_data['GS_E_au']:.6f} Ha")
    print(f"     Paper-27 S_full(GS) = {He_data['paper27_GS_S_full']:.6e}")
    print(f"     S_full first 5 states: "
          f"{[f'{s:.4f}' for s in He_data['first_5_S_paper27']]}")

    print("\n[3/3] Li+ (Z=3, 2 electrons, n_max=3)")
    Lip_data = run_helike(Z=3, n_max=3, n_T=30,
                          n_states_used=30, n_states_resolve=5)
    print(f"     GS energy = {Lip_data['GS_E_au']:.6f} Ha")
    print(f"     Paper-27 S_full(GS) = {Lip_data['paper27_GS_S_full']:.6e}")
    print(f"     S_full first 5 states: "
          f"{[f'{s:.4f}' for s in Lip_data['first_5_S_paper27']]}")

    # Sanity verification of the textbook identity
    for label, data in [('H', H_data), ('He', He_data), ('Li+', Lip_data)]:
        residuals = [r['identity_residual'] for r in data['thermal_table']]
        print(f"\n  Identity check {label}: max |S_thermo - S_microstate_info| "
              f"= {max(abs(r) for r in residuals):.2e}")

    out = {
        'sprint': 'TD Track 2',
        'date': '2026-05-08',
        'description': ('Temperature decoding on atomic thermal observables: '
                        'verifies textbook S_thermo = k_B * S_microstate_info '
                        'and connects the T->0 thermal-mix spatial 1-RDM '
                        'entropy to Paper 27 S_full(GS).'),
        'systems': {
            'hydrogen': H_data,
            'helium': He_data,
            'lithium_plus': Lip_data,
        },
        'notes': [
            'All energies, beta, T are in atomic units (Hartree, 1/Ha, Ha).',
            'k_B = 1 throughout. Multiplying S by k_B converts to J/K.',
            'Paper-27 S_full(H, GS) = 0 by single-determinant rigidity.',
            'Paper-27 S_full(He, GS) at n_max=3 is the V_ee correlation entropy.',
            'identity_residual checks textbook S_thermo = -sum_microstates P log P; '
            'should be machine epsilon at every beta.',
        ],
    }

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote results to {out_path}")
    return out


if __name__ == '__main__':
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'sprint_td_track2.json')
    main(out_path)
