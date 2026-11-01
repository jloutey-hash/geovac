"""Build tests/data/paper59_pp_2026_11_01.json for tests/test_paper59_hydrogen_p_functions.py (record PP,
docs/predictions/2026-10-29_hydrogen_p_functions.md, GO 2026-11-01).

A copy of stored results: the scorer's own assembly (debug/pp_score.py, everything above its scoring block, executed
unchanged) applied to the 109 exact points in debug/data/pp/.  Stores F_L, L = 0..150, to 40 significant digits for
the five scored items and the two 1s controls, and PC's parent O1 O2 q_L to L = 230.  No new computation.
    python debug/build_pp_guard_data_2026_11_01.py
"""
import json
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
src = (HERE / 'pp_score.py').read_text()
ns = {'__file__': str(HERE / 'pp_score.py')}
exec(compile(src.split('out, lines = {}, []')[0], 'pp_score_defs', 'exec'), ns)

items = {
    '1s-U': ns['q']('U', '0', '0', '0', '0'),
    '1s-E': ns['q']('E', '0', '0', '0', '0'),
    'PA-U': ns['PA']('U'),
    'PB-U': ns['PB']('U'),
    'PC-U': ns['PC']('U'),
    'PA-E': ns['PA']('E'),
    'PB-E': ns['PB']('E'),
    'O12-U-L230': ns['O12']('U', '0', '0', lm=230),
}
out = {
    'record': 'docs/predictions/2026-10-29_hydrogen_p_functions.md',
    'built_by': 'debug/build_pp_guard_data_2026_11_01.py',
    'convention': 'F_L / (4 zX^3 (z1 z2)^(3/2)) of the unnormalized integral, z1 = z2 = 6/5, zX = 13/10',
    'geometries': {'U': ['9/5', '12/5', '3/4'], 'E': ['19/10', '19/10', '1']},
    'F': {k: [mp.nstr(v, 40) for v in (seq if k.endswith('L230') else seq[:151])] for k, seq in items.items()},
}
dest = HERE.parent / 'tests' / 'data' / 'paper59_pp_2026_11_01.json'
dest.write_text(json.dumps(out, indent=0))
print('wrote', dest, {k: len(v) for k, v in out['F'].items()})
