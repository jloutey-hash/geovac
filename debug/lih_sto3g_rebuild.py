"""LB: rebuild of the LiH STO-3G qubit baseline (docs/predictions/2026-11-13_lih_sto3g_rebuild.md).

Control: same pipeline at 1.45 A vs OpenFermion's shipped H1-Li1_sto-3g_singlet_1.45.hdf5.
Main: R = 3.015 bohr.  Writes debug/data/lih_sto3g_rebuild.json.
"""
import json
import os

import numpy as np
import openfermion
from openfermion import MolecularData

from geovac.gaussian_reference import lih_sto3g
from geovac.measurement_grouping import count_qwc_groups
from geovac.qubit_encoding import build_fermion_op_from_integrals, exact_jordan_wigner

BOHR_A = 0.529177210903


def measure(sys_):
    fop = build_fermion_op_from_integrals(sys_['h1'], sys_['eri'], sys_['nuclear_repulsion'])
    qop = exact_jordan_wigner(fop)
    lam = sum(abs(c) for c in qop.terms.values())
    ident = qop.terms.get((), 0.0)
    return qop, dict(terms=len(qop.terms), lam=float(lam), lam_ni=float(lam - abs(ident)),
                     identity=float(np.real(ident)), qwc=count_qwc_groups(qop))


out = {}

# ---- control
d = os.path.join(os.path.dirname(openfermion.__file__), 'testing', 'data')
mol = MolecularData(filename=os.path.join(d, 'H1-Li1_sto-3g_singlet_1.45'))
mol.load()
ctrl = lih_sto3g(1.45 / BOHR_A)
_, m_ctrl = measure(ctrl)
ref_h1 = mol.one_body_integrals
ref_eri = np.asarray(mol.two_body_integrals).transpose(0, 3, 1, 2)  # OF physicist -> chemist
ref = dict(h1=ref_h1, eri=ref_eri, nuclear_repulsion=float(mol.nuclear_repulsion))
_, m_ref = measure(ref)
out['control'] = dict(
    R_bohr=1.45 / BOHR_A,
    hf_ours=ctrl['hf_energy'], hf_file=float(mol.hf_energy),
    dE=ctrl['hf_energy'] - float(mol.hf_energy),
    mo_eps_ours=[float(x) for x in ctrl['mo_energies']],
    mo_eps_file=[float(x) for x in mol.orbital_energies],
    ours=m_ctrl, file=m_ref,
    ref_nuclear_repulsion=float(mol.nuclear_repulsion),
)
print('control', json.dumps(out['control'], indent=1))

# ---- main
main = lih_sto3g(3.015)
qop, m_main = measure(main)
out['main'] = dict(R_bohr=3.015, hf=main['hf_energy'],
                   mo_eps=[float(x) for x in main['mo_energies']], **m_main)

# ---- invariance of the group count under a 1e-13 relative perturbation (3 seeds)
inv = []
for seed in range(3):
    rng = np.random.default_rng(seed)
    pert = openfermion.QubitOperator()
    for k, c in qop.terms.items():
        pert.terms[k] = c * (1.0 + 1e-13 * rng.standard_normal())
    inv.append(count_qwc_groups(pert))
out['main']['qwc_perturbed'] = inv
print('main', json.dumps(out['main'], indent=1))

with open('debug/data/lih_sto3g_rebuild.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=1)
