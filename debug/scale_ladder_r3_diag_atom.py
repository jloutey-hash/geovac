"""Neutral O atom (triplet, M_S = 1 determinant) UHF + fc-CCSD(T) in the ladder's Gaussian-fitted
Slater kits: the faithful one-centre probe for the R3 step (debug/scale_ladder_r3_diag_memo.md).
Reuses tests/test_water_ccsd.ccsd, which takes a canonical single determinant in spin-orbital form
(the Fock operator is diagonal per spin for canonical UHF orbitals, so fs[o, v] = 0 holds exactly).

    python debug/scale_ladder_r3_diag_atom.py
"""
import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('d', os.path.join(HERE, 'scale_ladder_r3_diag.py'))
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)
w, GE = d.w, d.sl.GE


def uhf(S, h, g, na, nb, tol=1e-10, maxit=500):
    X = GE.lowdin_orbitals(S)
    e, Cp = np.linalg.eigh(X.T @ h @ X)
    Ca = Cb = X @ Cp
    Da = Ca[:, :na] @ Ca[:, :na].T
    Db = Cb[:, :nb] @ Cb[:, :nb].T
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
            Da, Db = Da_new, Db_new
            break
        Da = 0.7 * Da_new + 0.3 * Da
        Db = 0.7 * Db_new + 0.3 * Db
        E_old = E
    else:
        raise RuntimeError('UHF did not converge')
    return E, (ea, Ca), (eb, Cb)


def spin_orbital_uhf(h, g, Ca, ea, Cb, eb, na, nb):
    """Spin orbitals ordered: occupied alpha, occupied beta, virtual alpha, virtual beta."""
    n = Ca.shape[1]
    cols = ([(0, Ca[:, i], ea[i]) for i in range(na)] + [(1, Cb[:, i], eb[i]) for i in range(nb)]
            + [(0, Ca[:, a], ea[a]) for a in range(na, n)] + [(1, Cb[:, a], eb[a]) for a in range(nb, n)])
    spin = np.array([c[0] for c in cols])
    C = np.array([c[1] for c in cols]).T
    eps = np.array([c[2] for c in cols])
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    same = (spin[:, None] == spin[None, :]).astype(float)
    phys = gm.transpose(0, 2, 1, 3) * same[:, None, :, None] * same[None, :, None, :]
    anti = phys - phys.transpose(0, 1, 3, 2)
    return np.diag(eps), anti, na + nb


def atom(oks, frozen_core=1):
    nuc = [(np.zeros(3), 8.0)]
    fns, T = d.kit(nuc, oks)
    S, h, g = w.ao_integrals(fns, nuc, T)
    na, nb = 5, 3
    E, (ea, Ca), (eb, Cb) = uhf(S, h, g, na, nb)
    fc = frozen_core
    fs, anti, nocc = spin_orbital_uhf(h, g, Ca[:, fc:], ea[fc:], Cb[:, fc:], eb[fc:], na - fc, nb - fc)
    Ec, Et = w.ccsd(fs, anti, nocc)
    s2 = na - nb  # M_S only; <S^2> not needed for a basis comparison
    return {'UHF': float(E), 'CCSD_fc': float(E + Ec), 'CCSD(T)_fc': float(E + Ec + Et),
            'corr_fc': float(Ec + Et), 'eps_a_occ': [float(x) for x in ea[:na]], 'eps_b_occ': [float(x) for x in eb[:nb]]}


if __name__ == '__main__':
    out = {}
    for name in ('R2', 'R3', 'R3dB', 'R3vB', 'R3sB'):
        t0 = time.time()
        out[name] = atom(d.SETS[name])
        e = out[name]
        print('%-5s UHF %.6f  fcCCSD %.6f  fcCCSD(T) %.6f  corr_fc %.6f  %.0f s' % (
            name, e['UHF'], e['CCSD_fc'], e['CCSD(T)_fc'], e['corr_fc'], time.time() - t0), flush=True)
    with open(os.path.join(d.DATA, 'scale_ladder_r3_diag_atom.json'), 'w') as f:
        json.dump({'sets': d.SETS, 'energies': out}, f, indent=1)
