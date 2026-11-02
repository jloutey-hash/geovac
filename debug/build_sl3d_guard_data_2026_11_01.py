"""Build tests/data/paper62_sl3d_2026_11_01.json for tests/test_paper62_scale_ladder_3d.py (record SL-3d,
docs/predictions/2026-11-01_scale_ladder_3d.md, GO 2026-11-01).

A copy of stored results, no new computation.  Sources (all read, none written):
  debug/data/scale_ladder.jsonl                registered (tied-3d) water rows (record SL); R2 = the R2rep row
  debug/data/scale_ladder_atom.json            registered O s/p exponents and groups per rung
  debug/data/scale_ladder_h2.json              the H2-set H exponents (R5H)
  debug/data/scale_ladder_3d_atom_<R>.json     atom-set 3d per rung, old/new atom UHF and fc-CCSD(T)
  debug/data/scale_ladder_3d_atom_evals_R3.jsonl  every R3 atom evaluation (for the record)
  debug/data/scale_ladder_3d.jsonl             new water rows R2d (control) and R3d, R4d, R5d, R5Hd
  debug/data/scale_ladder_r3_diag_atom.json, debug/data/scale_ladder_r3_diag.jsonl   control C1 (R3dB)
Refuses to build on a duplicate row (the record scores the first run only).
    python debug/build_sl3d_guard_data_2026_11_01.py
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / 'data'
E_KEYS = ('RHF', 'CCSD(T)_fc', 'CCSD(T)')


def rows(path):
    out = {}
    for line in open(path):
        r = json.loads(line)
        assert r['rung'] not in out, ('duplicate row', r['rung'], path)
        out[r['rung']] = r
    return out


old = rows(DATA / 'scale_ladder.jsonl')
sp = json.load(open(DATA / 'scale_ladder_atom.json'))
h2 = json.load(open(DATA / 'scale_ladder_h2.json'))['H3']
new = rows(DATA / 'scale_ladder_3d.jsonl')
diag_atom = json.load(open(DATA / 'scale_ladder_r3_diag_atom.json'))
diag = rows(DATA / 'scale_ladder_r3_diag.jsonl')

out = {
    'record': 'docs/predictions/2026-11-01_scale_ladder_3d.md',
    'built_by': 'debug/build_sl3d_guard_data_2026_11_01.py',
    'o_prims': [[1, 0], [1, 0], [2, 0], [2, 1], [3, 1]],
    'h_prims': [[1, 0], [2, 0], [2, 1]],
    'h_k1': {'h_ks': [1.0], 'h_groups': [0, 0, 0]},
    'h_h2set': {'h_ks': [float(x) for x in h2['ks']], 'h_groups': h2['groups']},
    'reference_fc': {'ccpvdz': -76.241025, 'slater_dzp': -76.281739},   # registry p62_ccpvdz_fc, p62_dzp_fc
    'registered_sp': {R: {'ks': sp[R]['ks'], 'groups': sp[R]['groups']} for R in ('R2', 'R3', 'R4', 'R5')},
    'old_water': {('R2' if k == 'R2rep' else k): {e: old[k][e] for e in E_KEYS}
                  for k in ('R2rep', 'R3', 'R4', 'R5', 'R5H')},
    'atom': {},
    'new_water': {k: {'o5': r['o_k_1s_1s_2s_2p_3p'], 'k3d': r['k3d'], 'h_ks': r['h_ks'],
                      'h_groups': r['h_groups'], **{e: r[e] for e in E_KEYS},
                      't_integrals': r['t_integrals']} for k, r in new.items()},
    'control_c1': {
        'atom': {k: {e: diag_atom['energies'][k][e] for e in ('UHF', 'CCSD(T)_fc')} for k in ('R3', 'R3dB')},
        'water_R3dB': {e: diag['R3dB'][e] for e in ('RHF', 'CCSD(T)_fc')},
    },
    'atom_evals_R3': [json.loads(l) for l in open(DATA / 'scale_ladder_3d_atom_evals_R3.jsonl')],
}
for R in ('R2', 'R3', 'R4', 'R5'):
    a = json.load(open(DATA / ('scale_ladder_3d_atom_%s.json' % R)))
    out['atom'][R] = {'s_p_exponents': a['s_p_exponents'], 'k3d_old': a['k3d_old'], 'k3d_new': a['k3d_new'],
                      'old': {e: a['old'][e] for e in ('UHF', 'CCSD(T)_fc')},
                      'new': {e: a['new'][e] for e in ('UHF', 'CCSD(T)_fc')},
                      'n_evals': a['n_evals']}

dest = HERE.parent / 'tests' / 'data' / 'paper62_sl3d_2026_11_01.json'
dest.write_text(json.dumps(out, indent=1))
print('wrote', dest, sorted(out['new_water']))
