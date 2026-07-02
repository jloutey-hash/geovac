"""WH8 Step-1 probe -- does the graph's skeleton measure generate Born statistics?

Three legs (PI-directed solo sprint, 2026-07-01, group4-recert pause):

  A. COINCIDENCE: on FLAT states (uniform amplitude over their support) the
     skeleton's counting measure equals Born exactly -- the coincidence locus.
     This is why "graph-native probability" survives casual inspection.
  B. PHYSICAL DIVERGENCE: the He graph-native CI ground state's Born weights
     vs the skeleton's candidate measures (uniform over the SD basis; uniform
     over the exact support; uniform over the S_z=0 sector) -- total-variation
     distances. The skeleton offers no state-dependent measure, so these
     counting measures are its only candidates.
  C. DYNAMICAL INSTABILITY: the coincidence locus is not preserved by the
     graph's own Hamiltonian -- a flat state evolves off it (TV grows ~t^2),
     so no skeleton-side rule can track Born through dynamics.

Expected verdict: NEGATIVE (the skeleton does not generate the Born measure).
Upgrades the Born rule's Class-1 classification (external-input three-class
partition) from categorization to TESTED NEGATIVE; supports the WH8 reading
(Born measure = the exchange constant of the observation projection; unique
given the projection lattice by Gleason's theorem, dim >= 3).

Frozen falsifier: tests/test_wh8_born_probe.py.
Data: debug/data/wh8_born_measure_probe.json.
"""
from __future__ import annotations

import json
import os
import warnings

import numpy as np
from scipy.sparse.linalg import eigsh


# ------------------------------------------------------------------ helpers

def tv_distance(p: np.ndarray, q: np.ndarray) -> float:
    """Total variation distance between two probability vectors."""
    return 0.5 * float(np.sum(np.abs(np.asarray(p) - np.asarray(q))))


def born(vec: np.ndarray) -> np.ndarray:
    w = np.abs(np.asarray(vec)) ** 2
    return w / w.sum()


# --------------------------------------------------------------- Leg A ----

def leg_a() -> dict:
    """Hydrogen n=2 multiplet (2s, 2p_{-1,0,+1}): flat vs generic states."""
    d = 4
    flat = np.ones(d) / np.sqrt(d)
    generic = np.array([1.0, 2.0, 3.0, 4.0])
    generic /= np.linalg.norm(generic)
    counting = np.ones(d) / d
    return {
        'dim': d,
        'tv_flat_exact_zero': tv_distance(born(flat), counting),
        'tv_generic': tv_distance(born(generic), counting),  # exactly 1/3
        'born_generic': born(generic).tolist(),
        'counting': counting.tolist(),
    }


# --------------------------------------------------------------- Leg B ----

def leg_b():
    """He graph-native CI (zero-parameter, exact Slater integrals) ground
    state: Born weights vs the skeleton's counting-measure candidates."""
    from geovac.lattice_index import LatticeIndex

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        li = LatticeIndex(n_electrons=2, max_n=3, nuclear_charge=2,
                          vee_method='slater_full')
        H = li.assemble_hamiltonian()
    n_sd = H.shape[0]
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    p_born = born(psi)

    # Skeleton candidates (the graph's only native measures are counting-type)
    q_basis = np.ones(n_sd) / n_sd
    supp = p_born > 1e-14
    q_support = np.where(supp, 1.0, 0.0)
    q_support /= q_support.sum()
    sz = np.array([
        sum(0.5 if li.sp_states[o][3] == 0 else -0.5 for o in sd)
        for sd in li.sd_basis
    ])
    sec = np.isclose(sz, 0.0)
    q_sector = np.where(sec, 1.0, 0.0)
    q_sector /= q_sector.sum()

    order = np.argsort(p_born)[::-1]
    top = []
    for I in order[:5]:
        labels = [li.sp_states[o] for o in li.sd_basis[I]]
        top.append({
            'weight': float(p_born[I]),
            'config': ['n{}l{}m{}{}'.format(n, l, m, 'ud'[s])
                       for (n, l, m, s) in labels],
        })

    result = {
        'system': 'He graph-native CI, max_n=3, vee=slater_full (0 params)',
        'n_sd': int(n_sd),
        'E0_Ha': float(evals[0]),
        'support_size': int(supp.sum()),
        'sz0_sector_size': int(sec.sum()),
        'top_configs': top,
        'tv_born_vs_uniform_basis': tv_distance(p_born, q_basis),
        'tv_born_vs_uniform_support': tv_distance(p_born, q_support),
        'tv_born_vs_uniform_sz0_sector': tv_distance(p_born, q_sector),
    }
    return result, H, p_born


# --------------------------------------------------------------- Leg C ----

def leg_c(H, p_born: np.ndarray) -> dict:
    """Evolve a FLAT state (on the ground state's support) under the graph's
    own Hamiltonian: the coincidence locus TV(Born, counting) leaves zero."""
    supp = np.where(p_born > 1e-14)[0]
    n_sd = len(p_born)
    psi0 = np.zeros(n_sd)
    psi0[supp] = 1.0 / np.sqrt(len(supp))
    counting = np.zeros(n_sd)
    counting[supp] = 1.0 / len(supp)

    Hd = H.toarray()
    E, V = np.linalg.eigh(Hd)
    c0 = V.T @ psi0

    ts = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    tvs = []
    for t in ts:
        psi_t = V @ (np.exp(-1j * E * t) * c0)
        tvs.append(tv_distance(born(psi_t), counting))

    # small-t growth exponent from the first two nonzero points
    t1, t2 = ts[1], ts[2]
    v1, v2 = max(tvs[1], 1e-300), max(tvs[2], 1e-300)
    slope = float(np.log(v2 / v1) / np.log(t2 / t1))
    return {
        'note': 'flat state on the GS support, evolved under the CI H',
        't_grid': ts,
        'tv_vs_counting': tvs,
        'small_t_growth_exponent': slope,   # ~2 expected (perturbative t^2)
    }


# ------------------------------------------------------------------- main -

def main() -> dict:
    print('=== WH8 Step-1 probe: skeleton measure vs Born ===\n')
    a = leg_a()
    print('[A] coincidence locus (H n=2 multiplet, d=4):')
    print('    TV(flat)    = {:.3e}   (exact 0: flat states ARE the locus)'
          .format(a['tv_flat_exact_zero']))
    print('    TV(generic) = {:.6f} (= 1/3 for amplitudes prop 1,2,3,4)\n'
          .format(a['tv_generic']))

    b, H, p_born = leg_b()
    print('[B] He graph-native CI ground state (N_SD = {}, E0 = {:.6f} Ha):'
          .format(b['n_sd'], b['E0_Ha']))
    print('    top config weight  = {:.4f}  ({})'
          .format(b['top_configs'][0]['weight'],
                  ' '.join(b['top_configs'][0]['config'])))
    print('    TV vs uniform basis   = {:.4f}'.format(b['tv_born_vs_uniform_basis']))
    print('    TV vs uniform support = {:.4f}'.format(b['tv_born_vs_uniform_support']))
    print('    TV vs uniform Sz=0    = {:.4f}\n'.format(b['tv_born_vs_uniform_sz0_sector']))

    c = leg_c(H, p_born)
    print('[C] dynamical instability of the coincidence locus:')
    for t, tv in zip(c['t_grid'], c['tv_vs_counting']):
        print('    t = {:5.2f}   TV = {:.6f}'.format(t, tv))
    print('    small-t growth exponent ~ {:.2f} (t^2 expected)\n'
          .format(c['small_t_growth_exponent']))

    negative = (b['tv_born_vs_uniform_support'] > 0.5
                and c['tv_vs_counting'][-1] > 1e-3)
    verdict = ('NEGATIVE (expected): the skeleton counting measure does not '
               'generate Born statistics on physical states, and the '
               'coincidence locus is dynamically unstable.'
               if negative else
               'UNEXPECTED: review before concluding.')
    print('VERDICT:', verdict)

    out = {'leg_a': a, 'leg_b': b, 'leg_c': c, 'verdict': verdict}
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'data', 'wh8_born_measure_probe.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print('\nwrote', out_path)
    return out


if __name__ == '__main__':
    main()
