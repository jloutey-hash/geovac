"""Build tests/data/paper62_scale_ladder_2026_10_30.json from the stored ladder rows.

Only the quantities that cost ~20 minutes each to produce (the water energies) are
copied.  The free-O-atom energies are NOT copied as truth: the guard recomputes them
from the stored exponents through tests/test_water_ccsd.radial_atom and checks that
the exponents are stationary minima of that objective.  The R3 atom-probe numbers are
copied for reference only; the guard recomputes those too.

    python debug/build_p62_guard_data.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

atom = json.load(open(os.path.join(HERE, 'data', 'scale_ladder_atom.json')))
rows = {}
for line in open(os.path.join(HERE, 'data', 'scale_ladder.jsonl')):
    r = json.loads(line)
    rows[r['rung']] = r
probe = json.load(open(os.path.join(HERE, 'data', 'scale_ladder_r3_diag_atom.json')))

out = {
    'source': 'debug/data/scale_ladder_atom.json, scale_ladder.jsonl, scale_ladder_r3_diag_atom.json (2026-10-30)',
    'geometry': {'R_OH': 1.809, 'angle_deg': 104.5},
    'o_prims': [[1, 0], [1, 0], [2, 0], [2, 1], [3, 1]],
    'rungs': {},
    'r3_probe_reference': {k: {'UHF': v['UHF'], 'corr_fc': v['corr_fc']} for k, v in probe['energies'].items()},
}
for name, key in (('R2', 'R2rep'), ('R3', 'R3'), ('R4', 'R4'), ('R5', 'R5'), ('R5H', 'R5H')):
    r = rows[key]
    d = {'o_ks': r['o_ks'], 'o_groups': r['o_groups'], 'h_ks': r['h_ks'], 'h_groups': r['h_groups'],
         'RHF': r['RHF'], 'CCSD(T)_fc': r['CCSD(T)_fc'], 'CCSD(T)': r['CCSD(T)']}
    if name in atom:
        d['atom_E_reference'] = atom[name]['E']
    out['rungs'][name] = d

path = os.path.join(ROOT, 'tests', 'data', 'paper62_scale_ladder_2026_10_30.json')
with open(path, 'w') as f:
    json.dump(out, f, indent=1)
print('wrote', path)
