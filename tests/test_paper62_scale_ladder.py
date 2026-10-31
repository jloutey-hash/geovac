"""Paper 62, the scale ladder (Sec. ladder, Sec. r3): regression guard.

Record: docs/predictions/2026-10-30_scale_ladder.md.  Written as a separate step after
the paper (CLAUDE.md Sec. 9, guard-writing rule), fire-tested against the wrong answers
named in each test's docstring.

What is recomputed here and what is read.  The free-oxygen-atom radial HF energies
(the ladder's objective) are RECOMPUTED from the stored exponents through
tests/test_water_ccsd.radial_atom, and the exponents are checked to be stationary
minima of that objective.  The R3 one-centre probe (neutral O atom, UHF + fc-CCSD(T) in
the Gaussian-fitted kits) is RECOMPUTED.  Only the water energies are READ from
tests/data/paper62_scale_ladder_2026_10_30.json, because each costs about twenty
minutes of integrals; the ratios the paper quotes are then derived from those stored
energies and the recomputed atom energies, never copied.
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

DATA = json.load(open(os.path.join(HERE, 'data', 'paper62_scale_ladder_2026_10_30.json')))
O_PRIMS = [tuple(p) for p in DATA['o_prims']]
H_PRIMS = [(1, 0), (2, 0), (2, 1)]
KC, KV = 7.79324757, 2.28124155          # the v5.20.12 two-scale kit (CHANGELOG v5.20.12)


# ----------------------------------------------------------------------------- objective
def atom_energy(ks, groups):
    """The ladder's objective: radial average-of-configuration HF of the free O atom."""
    s = [(n, ks[g]) for (n, l), g in zip(O_PRIMS, groups) if l == 0]
    p = [(n, ks[g]) for (n, l), g in zip(O_PRIMS, groups) if l == 1]
    return w.radial_atom(s, p)[0]


def atom_energies():
    return {r: atom_energy(DATA['rungs'][r]['o_ks'], DATA['rungs'][r]['o_groups'])
            for r in ('R2', 'R3', 'R4', 'R5')}


def test_atom_exponents_are_stationary_minima_of_the_radial_hf():
    """Every rung's stored exponents are a stationary local minimum of the objective
    the record names, and R2 is the v5.20.12 two-scale kit.

    Wrong answers rejected: exponents optimized on some other objective (a different
    occupation, a 3D fitted energy, water itself), a rung whose exponents were copied
    from another rung, and a base rung that is not the kit it claims to reproduce.
    """
    for r, E in atom_energies().items():
        ks, groups = DATA['rungs'][r]['o_ks'], DATA['rungs'][r]['o_groups']
        for g in range(len(ks)):
            h = 1e-3 * ks[g]
            up, dn = list(ks), list(ks)
            up[g] += h
            dn[g] -= h
            Eu, Ed = atom_energy(up, groups), atom_energy(dn, groups)
            grad = (Eu - Ed) / (2 * h)
            # measured 2026-10-31: |grad| = 6e-8 .. 1.7e-6 Ha/bohr^-1 on every rung
            assert abs(grad) < 5e-5, (r, g, grad)                 # stationary
            assert Eu > E - 1e-9 and Ed > E - 1e-9, (r, g, Eu - E, Ed - E)   # minimum, not saddle
    r2 = DATA['rungs']['R2']['o_ks']
    assert abs(r2[0] - KC) < 1e-6 and abs(r2[1] - KV) < 1e-6, r2


def test_the_p_split_carries_the_atom_and_water_gains_and_transfer_is_about_half():
    """SL1 and SL2 of the record, derived here rather than copied: the 2p/3p split
    (R4 -> R5) carries 0.90 of the gain on the atom (recomputed) and on water (stored),
    and the water gain is 0.52 of the atom's.  The paper's quoted 0.900 and 0.516 must
    be what these numbers give.

    Wrong answers rejected: a ladder where the s rungs carry the gain, a transfer ratio
    outside the registered band, stored water rows that do not match the paper's ratios
    (data/paper drift), and atom steps that are not the recomputed objective's.
    """
    A = atom_energies()
    W = {r: DATA['rungs'][r]['CCSD(T)_fc'] for r in ('R2', 'R3', 'R4', 'R5', 'R5H')}
    atom_total = (A['R2'] - A['R5']) * 1e3
    assert abs(atom_total - 257.3) < 0.1, atom_total
    atom_p_share = (A['R4'] - A['R5']) / (A['R2'] - A['R5'])
    assert abs(atom_p_share - 0.902) < 0.002, atom_p_share
    G = W['R2'] - W['R5']
    assert G > 0
    sl1 = (W['R4'] - W['R5']) / G
    sl2 = G * 1e3 / atom_total
    assert sl1 >= 0.6 and abs(sl1 - 0.900) < 0.001, sl1
    assert 0.3 <= sl2 <= 0.8 and abs(sl2 - 0.516) < 0.001, sl2
    # the non-monotone rung and the hydrogen rung, as reported
    assert W['R3'] > W['R2'] and A['R3'] < A['R2']
    assert W['R5H'] < W['R5']


# ----------------------------------------------------------------------- the R3 probe
def _kit(oks, h_k=1.0, n_gauss=12, ng_pol=8):
    """O-only one-centre kit: the five ladder primitives plus a 3d at oks[5]."""
    nuclei = [(np.zeros(3), 8.0)]
    fns, cols_of = [], []
    for (n, l), k in zip(O_PRIMS, oks[:5]):
        w._shell_kit(np.zeros(3), [(l, {n: 1.0}, k)], n_gauss, ng_pol, fns, [], cols_of)
    w._shell_kit(np.zeros(3), [(2, {3: 1.0}, oks[5])], n_gauss, ng_pol, fns, [], cols_of)
    return nuclei, fns, w._kit_columns(fns, cols_of)


def _uhf(S, h, g, na, nb, tol=1e-10, maxit=500):
    X = w.GE.lowdin_orbitals(S)
    _, Cp = np.linalg.eigh(X.T @ h @ X)
    Ca = Cb = X @ Cp
    Da, Db = Ca[:, :na] @ Ca[:, :na].T, Cb[:, :nb] @ Cb[:, :nb].T
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
        Da_new, Db_new = Ca[:, :na] @ Ca[:, :na].T, Cb[:, :nb] @ Cb[:, :nb].T
        if abs(E - E_old) < tol and np.abs(Da_new - Da).max() < 1e-7 and it > 2:
            break
        Da, Db = 0.7 * Da_new + 0.3 * Da, 0.7 * Db_new + 0.3 * Db
        E_old = E
    else:
        raise RuntimeError('UHF did not converge')
    return E, (ea, Ca), (eb, Cb)


def _spin_orbital_uhf(h, g, Ca, ea, Cb, eb, na, nb):
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


def _o_atom(oks, frozen_core=1):
    """Neutral O, M_S = 1 determinant (5 alpha, 3 beta), UHF + fc-CCSD(T)."""
    nuclei, fns, T = _kit(oks)
    S, h, g = w.ao_integrals(fns, nuclei, T)
    na, nb = 5, 3
    E, (ea, Ca), (eb, Cb) = _uhf(S, h, g, na, nb)
    fc = frozen_core
    fs, anti, nocc = _spin_orbital_uhf(h, g, Ca[:, fc:], ea[fc:], Cb[:, fc:], eb[fc:], na - fc, nb - fc)
    Ec, Et = w.ccsd(fs, anti, nocc)
    return float(E), float(Ec + Et)


@pytest.mark.slow
def test_r3_correlation_loss_is_a_one_centre_effect_carried_by_the_3d():
    """Sec. r3: on the free O atom the R2 -> R3 step gains at UHF and loses at
    correlation by nearly equal amounts (a wash), the correlation loss equals the
    2.72 mHa water loses, and putting the 3d alone back at R2's exponent removes
    about 80% of it.

    Wrong answers rejected: a step whose loss is mean-field rather than correlation;
    a loss that is not the 3d's (e.g. the s shift's); a correlation loss on the atom
    that does not equal the water one (which would make the water step a transfer
    failure rather than a one-centre effect).
    """
    r3 = DATA['rungs']['R3']['o_ks']                     # [core, s-valence, p family]
    KC3, KS3, KP3 = r3
    sets = {'R2': (KC, KV, KV, KV, KV, KV),
            'R3': (KC3, KS3, KS3, KP3, KP3, KP3),
            'R3dB': (KC3, KS3, KS3, KP3, KP3, KV)}
    E = {k: _o_atom(v) for k, v in sets.items()}
    uhf_gain = (E['R3'][0] - E['R2'][0]) * 1e3            # negative = gain
    corr_loss = (E['R3'][1] - E['R2'][1]) * 1e3           # positive = loss
    assert -2.95 < uhf_gain < -2.75, uhf_gain
    assert 2.6 < corr_loss < 2.85, corr_loss
    assert abs(uhf_gain + corr_loss) < 0.3                 # a wash at fc-CCSD(T)
    water_corr_loss = 1e3 * ((DATA['rungs']['R3']['CCSD(T)_fc'] - DATA['rungs']['R3']['RHF'])
                             - (DATA['rungs']['R2']['CCSD(T)_fc'] - DATA['rungs']['R2']['RHF']))
    assert abs(water_corr_loss - corr_loss) < 0.1, (water_corr_loss, corr_loss)
    corr_loss_3d_back = (E['R3dB'][1] - E['R2'][1]) * 1e3
    assert corr_loss_3d_back < 0.3 * corr_loss, (corr_loss_3d_back, corr_loss)
