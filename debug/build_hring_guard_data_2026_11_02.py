"""Build tests/data/paper62_hring_2026_11_02.json, the stored copy the Paper 62 H6-ring guard reads
(tests/test_paper62_hydrogen_ring.py).  Copies the scored rows of debug/data/hring_study.jsonl (modes
'final' and 'labels') and the re-optimized p exponents of debug/data/hring_exponents.json, with the
sha256 of each source.  Run once; the guard never imports from debug/.
"""
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ROWS = os.path.join(HERE, 'data', 'hring_study.jsonl')
EXP = os.path.join(HERE, 'data', 'hring_exponents.json')
OUT = os.path.join(ROOT, 'tests', 'data', 'paper62_hring_2026_11_02.json')


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


final, labels = {}, {}
for line in open(ROWS):
    r = json.loads(line)
    if r.get('mode') == 'final':
        final[r['basis']] = {'nbf': r['nbf'], 'RHF': r['RHF'], 'CCSD': r['CCSD'], 'CCSD(T)': r['CCSD(T)'],
                             'qubit': r['qubit']}
    elif r.get('mode') == 'labels':
        labels[r['basis']] = {'site': r['site'], 'ring': r['ring']}
exp = json.load(open(EXP))
out = {'source': {'hring_study.jsonl': sha(ROWS), 'hring_exponents.json': sha(EXP)},
       'geometry': {'N': 6, 'spacing_ang': 1.0},
       'p_exponents': {'sto_p': exp['sto_p'], 'gto_p': exp['gto_p']},
       'final': final, 'labels': labels}
json.dump(out, open(OUT, 'w'), indent=1)
print('wrote', OUT, sorted(final), sorted(labels))
