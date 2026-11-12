"""Heavy-atom panel (record HP, docs/predictions/2026-11-12_heavy_atom_panel.md).

Question: does an exponential (Slater) basis set entirely on free atoms -- Paper 62's best atom-set rung, R5H with
the 3d also set on the atom -- beat the best same-size Gaussian basis tested (ANO-R1) on molecules with a first-row
heavy atom, and does the margin grow with nuclear charge?  Paper 62 measured this rung on water only (-76.260011,
4.7 mHa behind ANO-R1).

The kit is Paper 62's, generalized from O to X = C, N, O, F with nothing else changed:
  X: Slater primitives (1s, 1s', 2s, 2p, 3p), five distinct exponents minimizing the radial average-of-configuration
     Hartree-Fock energy of the free atom (tests/test_water_ccsd.radial_atom, Nelder-Mead from 41 starts), plus one
     3d whose exponent minimizes the free atom's frozen-core CCSD(T) energy, UHF reference for the high-spin
     determinant (C, O: M_S = 1; N: M_S = 3/2; F: M_S = 1/2), grid then bounded Brent in ln k (record SL-3d rule).
  H: 1s, 2s, 2p from H2 RHF at R = 1.4 bohr (debug/data/scale_ladder_h2.json 'H3'), as on water's R5H rung.
  12 Gaussians per s/p primitive, 8 per polarization function.  14 functions on X, 5 on H: the cc-pVDZ and ANO-R1
  counts exactly.

    python debug/heavy_atom_panel.py atom-sp C|N|O|F          s/p exponents (O reproduces record SL's R5)
    python debug/heavy_atom_panel.py atom-3d C|N|O|F          3d exponent (O reproduces record SL-3d's R5)
    python debug/heavy_atom_panel.py atom-gauss C|N|O|F ccpvdz|anor1|slater   free-atom control energies
    python debug/heavy_atom_panel.py mol CH4|NH3|H2O|HF slater|ccpvdz|anor1   SCORED (after registration)
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import time
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('sl', os.path.join(HERE, 'scale_ladder.py'))
sl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sl)            # applies scale_ladder's np.diag patch; loads test_water_ccsd as sl.w
w, GE = sl.w, sl.GE
_spec2 = importlib.util.spec_from_file_location('ra', os.path.join(HERE, 'scale_ladder_r3_diag_atom.py'))
ra = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(ra)           # uhf, spin_orbital_uhf

DATA = os.path.join(HERE, 'data')
OUT = os.path.join(DATA, 'heavy_atom_panel')
BASIS = os.path.join(DATA, 'basis', '{}_H_C_N_O_F.json')
GAUSS = {'ccpvdz': 'cc-pvdz', 'anor1': 'ano-r1'}
ANG = 1.8897261246257702                # bohr per angstrom (CODATA 2018)

# Z, open-p occupancy (average of configuration), high-spin UHF (n_alpha, n_beta)
ELEM: Dict[str, Dict] = {
    'C': dict(Z=6.0, q_p=2, na=4, nb=2),
    'N': dict(Z=7.0, q_p=3, na=5, nb=2),
    'O': dict(Z=8.0, q_p=4, na=5, nb=3),
    'F': dict(Z=9.0, q_p=5, na=5, nb=4),
}
O_PRIMS = sl.O_PRIMS                    # (1,0), (1,0), (2,0), (2,1), (3,1)
H_PRIMS = sl.H_PRIMS                    # (1,0), (2,0), (2,1)
GRID_3D = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5]


def _path(name: str) -> str:
    os.makedirs(OUT, exist_ok=True)
    return os.path.join(OUT, name)


def _write_json(name: str, obj) -> None:
    p = _path(name)
    with open(p + '.tmp', 'w') as f:
        json.dump(obj, f, indent=1)
    os.replace(p + '.tmp', p)


# ---------------------------------------------------------------- geometry (fixed in record HP before scoring)
def molecule(name: str) -> Tuple[List[Tuple[np.ndarray, float]], str, int]:
    """Nuclei (bohr), heavy element, electron count."""
    if name == 'H2O':
        return w.water_nuclei(), 'O', 10                       # Paper 62: R 1.809 bohr, 104.5 deg
    if name == 'HF':
        return [(np.zeros(3), 9.0), (np.array([0.0, 0.0, 1.7325]), 1.0)], 'F', 10   # Paper 62 sec:hf
    if name == 'NH3':
        r, ang = 1.0124 * ANG, math.radians(106.67)
        # three H at equal angle ang from each other, N at the origin, C3 axis along z
        cos_t = math.sqrt((1 + 2 * math.cos(ang)) / 3.0)        # angle between N-H and -z axis
        sin_t = math.sqrt(1 - cos_t ** 2)
        Hs = [np.array([r * sin_t * math.cos(phi), r * sin_t * math.sin(phi), -r * cos_t])
              for phi in (0.0, 2 * math.pi / 3, 4 * math.pi / 3)]
        return [(np.zeros(3), 7.0)] + [(h, 1.0) for h in Hs], 'N', 10
    if name == 'CH4':
        r = 1.0870 * ANG / math.sqrt(3.0)
        Hs = [np.array(v, float) * r for v in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
        return [(np.zeros(3), 6.0)] + [(h, 1.0) for h in Hs], 'C', 10
    raise KeyError(name)


# ---------------------------------------------------------------- atom side
def atom_sp_energy(ks, el: str) -> float:
    e = ELEM[el]
    s = [(n, k) for (n, l), k in zip(O_PRIMS, ks) if l == 0]
    p = [(n, k) for (n, l), k in zip(O_PRIMS, ks) if l == 1]
    try:
        E = w.radial_atom(s, p, Z=e['Z'], q_p=e['q_p'])[0]
    except (RuntimeError, np.linalg.LinAlgError, ValueError):   # a wild start: no convergence, or a singular overlap
        return 1e3
    return float(E) if np.isfinite(E) else 1e3


def atom_sp(el: str) -> Dict:
    Z = ELEM[el]['Z']
    # distinct generic starts (two equal 1s exponents make the overlap singular); scaled by Z/8 from the
    # two-scale v5.20.12 values, not from any optimized set
    base = np.array([sl.KC, 0.7 * sl.KV, 0.8 * sl.KC, 1.4 * sl.KV, 0.9 * sl.KV]) * Z / 8.0
    cand = [base]
    rng = np.random.default_rng(0)
    cand += [base * np.exp(rng.normal(0, 0.5, 5)) for _ in range(40)]
    best = None
    for x0 in cand:
        r = minimize(lambda k: atom_sp_energy(np.abs(k), el), x0, method='Nelder-Mead',
                     options=dict(xatol=1e-8, fatol=1e-11, maxiter=6000))
        if best is None or r.fun < best.fun:
            best = r
    ks = [float(x) for x in np.abs(best.x)]
    return {'element': el, 'E_avg_config_HF': float(best.fun), 'ks_1s_1s_2s_2p_3p': ks, 'n_starts': len(cand)}


def heavy_kit(nuclei, el: str, five, k3d: float, h_ks, h_groups, n_gauss=12, ng_pol=8):
    Z = ELEM[el]['Z']
    fns, cols_of = [], []
    for center, Zc in nuclei:
        if Zc == Z:
            for (n, l), k in zip(O_PRIMS, five):
                w._shell_kit(center, [(l, {n: 1.0}, k)], n_gauss, ng_pol, fns, [], cols_of)
            w._shell_kit(center, [(2, {3: 1.0}, k3d)], n_gauss, ng_pol, fns, [], cols_of)
        elif Zc == 1.0:
            for (n, l), g in zip(H_PRIMS, h_groups):
                w._shell_kit(center, [(l, {n: 1.0}, h_ks[g])], n_gauss, ng_pol, fns, [], cols_of)
        else:
            raise ValueError(Zc)
    return fns, w._kit_columns(fns, cols_of)


def _uhf_from(S, h, g, na, nb, C0, tol=1e-10, maxit=800):
    """ra.uhf with a given starting orbital set (same damping and convergence test)."""
    X = GE.lowdin_orbitals(S)
    Da = C0[:, :na] @ C0[:, :na].T
    Db = C0[:, :nb] @ C0[:, :nb].T
    E_old = 0.0
    for it in range(maxit):
        J = np.einsum('pqrs,rs->pq', g, Da + Db)
        Ka = np.einsum('prqs,rs->pq', g, Da)
        Kb = np.einsum('prqs,rs->pq', g, Db)
        Fa, Fb = h + J - Ka, h + J - Kb
        E = 0.5 * (np.sum((Da + Db) * h) + np.sum(Da * Fa) + np.sum(Db * Fb))
        ea, Cpa = np.linalg.eigh(X.T @ Fa @ X)
        eb, Cpb = np.linalg.eigh(X.T @ Fb @ X)
        Ca, Cb = X @ Cpa, X @ Cpb
        Da_new = Ca[:, :na] @ Ca[:, :na].T
        Db_new = Cb[:, :nb] @ Cb[:, :nb].T
        if abs(E - E_old) < tol and np.abs(Da_new - Da).max() < 1e-7 and it > 2:
            return float(E), (ea, Ca), (eb, Cb)
        Da = 0.7 * Da_new + 0.3 * Da
        Db = 0.7 * Db_new + 0.3 * Db
        E_old = E
    raise RuntimeError('UHF did not converge')


def lowest_uhf(S, h, g, na, nb, n_perturbed=5, seed=1):
    """The default core-guess UHF can stop in an excited state: on the C atom it gave -37.4243 against -37.6878
    from perturbed starts (measured 2026-11-12).  Take the lowest of the default start and n_perturbed
    perturbed core-guess starts (fixed seed)."""
    best = ra.uhf(S, h, g, na, nb)
    X = GE.lowdin_orbitals(S)
    _, P = np.linalg.eigh(X.T @ h @ X)
    rng = np.random.default_rng(seed)
    for _ in range(n_perturbed):
        C0 = X @ np.linalg.qr(P + 0.3 * rng.normal(size=P.shape))[0]
        try:
            r = _uhf_from(S, h, g, na, nb, C0)
        except RuntimeError:
            continue
        if r[0] < best[0] - 1e-9:
            best = r
    return best


def atom_uhf_ccsdt(el: str, S, h, g) -> Dict:
    e = ELEM[el]
    na, nb = e['na'], e['nb']
    E, (ea, Ca), (eb, Cb) = lowest_uhf(S, h, g, na, nb)
    fc = 1
    fs, anti, nocc = ra.spin_orbital_uhf(h, g, Ca[:, fc:], ea[fc:], Cb[:, fc:], eb[fc:], na - fc, nb - fc)
    Ec, Et = w.ccsd(fs, anti, nocc)
    return {'UHF': float(E), 'CCSD_fc': float(E + Ec), 'CCSD(T)_fc': float(E + Ec + Et)}


def atom_slater(el: str, k3d: float) -> Dict:
    nuc = [(np.zeros(3), ELEM[el]['Z'])]
    five = json.load(open(_path('atom_sp_%s.json' % el)))['ks_1s_1s_2s_2p_3p']
    fns, T = heavy_kit(nuc, el, five, k3d, [1.0], [0, 0, 0])
    S, h, g = w.ao_integrals(fns, nuc, T)
    return atom_uhf_ccsdt(el, S, h, g)


def atom_3d(el: str) -> Dict:
    evals_path = _path('atom_3d_evals_%s.jsonl' % el)
    cache = {}
    if os.path.exists(evals_path):
        for line in open(evals_path):
            r = json.loads(line)
            cache[round(r['k3d'], 10)] = r

    def ev(k):
        key = round(float(k), 10)
        if key not in cache:
            t0 = time.time()
            r = {'k3d': float(k), **atom_slater(el, float(k)), 't': round(time.time() - t0, 1)}
            with open(evals_path, 'a') as f:
                f.write(json.dumps(r) + '\n')
            cache[key] = r
            print('%s k3d %.5f  UHF %.6f  fcCCSD(T) %.6f  %.0f s' % (el, k, r['UHF'], r['CCSD(T)_fc'], r['t']),
                  flush=True)
        return cache[key]

    grid = [ev(k) for k in GRID_3D]
    i = int(np.argmin([r['CCSD(T)_fc'] for r in grid]))
    if i in (0, len(GRID_3D) - 1):
        raise RuntimeError('%s: 3d grid minimum at the edge (%.2f); widen GRID_3D' % (el, GRID_3D[i]))
    lo, hi = math.log(GRID_3D[i - 1]), math.log(GRID_3D[i + 1])
    res = minimize_scalar(lambda x: ev(math.exp(x))['CCSD(T)_fc'], bounds=(lo, hi), method='bounded',
                          options=dict(xatol=1e-3))
    best = ev(math.exp(res.x))
    return {'element': el, 'k3d': best['k3d'], 'atom': best, 'n_evals': len(cache)}


def atom_gauss(el: str, basis: str) -> Dict:
    nuc = [(np.zeros(3), ELEM[el]['Z'])]
    if basis == 'slater':
        k3d = json.load(open(_path('atom_3d_%s.json' % el)))['k3d']
        r = atom_slater(el, k3d)
    else:
        fns, T = w.gaussian_basis(BASIS.format(GAUSS[basis]), nuc)
        S, h, g = w.ao_integrals(fns, nuc, T)
        r = atom_uhf_ccsdt(el, S, h, g)
    return {'element': el, 'basis': basis, **r}


# ---------------------------------------------------------------- molecule (scored)
def mol(name: str, basis: str) -> Dict:
    nuc, el, nelec = molecule(name)
    t0 = time.time()
    if basis == 'slater':
        five = json.load(open(_path('atom_sp_%s.json' % el)))['ks_1s_1s_2s_2p_3p']
        k3d = json.load(open(_path('atom_3d_%s.json' % el)))['k3d']
        h2 = json.load(open(os.path.join(DATA, 'scale_ladder_h2.json')))['H3']
        fns, T = heavy_kit(nuc, el, five, k3d, h2['ks'], h2['groups'])
    else:
        fns, T = w.gaussian_basis(BASIS.format(GAUSS[basis]), nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    enuc = w.vnn(nuc)
    t1 = time.time()
    # the default RHF start can stop in an excited state (the UHF one did on the C atom): if a perturbed start
    # (same closed-shell equations, run as UHF with na = nb) finds a lower energy, stop instead of scoring
    e_rhf = w.rhf(S, h, g, nelec, enuc)[0]
    e_alt = lowest_uhf(S, h, g, nelec // 2, nelec // 2)[0] + enuc
    if e_alt < e_rhf - 1e-6:
        raise RuntimeError('%s/%s: default RHF %.6f above a perturbed start %.6f' % (name, basis, e_rhf, e_alt))
    r = w.run(S, h, g, nelec, enuc, frozen_core=1)
    _, _, C = w.rhf(S, h, g, nelec, enuc)
    hm = C.T @ h @ C
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    hh, gg, ecore = w.freeze_core(hm, gm, 1)
    st = w.jw_pauli_stats(hh, gg, enuc + ecore)
    row = {'molecule': name, 'basis': basis, 'nbf': int(T.shape[1]), 'qubits_fc': 2 * hh.shape[0],
           'RHF': float(r['RHF']), 'CCSD_fc': float(r['CCSD']), 'CCSD(T)_fc': float(r['CCSD(T)']),
           'pauli_fc': int(st['n_pauli']), 'one_norm_fc': float(st['one_norm']),
           'rhf_stability_check': float(e_alt),
           't_integrals': round(t1 - t0, 1), 't_solve': round(time.time() - t1, 1)}
    np.savez(_path('mo_%s_%s.npz' % (name, basis)), h1=hh, eri=gg, ecore=enuc + ecore)
    return row


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'atom-sp':
        r = atom_sp(sys.argv[2])
        _write_json('atom_sp_%s.json' % sys.argv[2], r)
    elif cmd == 'atom-3d':
        r = atom_3d(sys.argv[2])
        _write_json('atom_3d_%s.json' % sys.argv[2], r)
    elif cmd == 'atom-gauss':
        r = atom_gauss(sys.argv[2], sys.argv[3])
        with open(_path('atom_controls.jsonl'), 'a') as f:
            f.write(json.dumps(r) + '\n')
    elif cmd == 'mol':
        r = mol(sys.argv[2], sys.argv[3])
        with open(_path('molecules.jsonl'), 'a') as f:
            f.write(json.dumps(r) + '\n')
    else:
        raise SystemExit(__doc__)
    print(json.dumps(r), flush=True)
