"""ANO-family matched-accuracy ladder for water (record AM, 2026-11-11).

Paper 62's matched-accuracy ladder (record MA) is the cc-pV family only. At 24 functions
ANO-R1 already sits 17.1 mHa above the Slater DZP instead of cc-pVDZ's 40.7. This asks
which ANO rung first reaches the Slater DZP, and at what qubit cost.

ANO-R1 is a column truncation of ANO-R2 (same primitives; R1's contraction columns are
R2's first columns, checked 2026-11-11), so every rung here is ANO-R2 kept to the first
n contracted functions per shell:

    A24 = ANO-R1          O [3s2p1d]    H [2s1p]     24 functions on water
    A35                   O [3s2p2d]    H [2s2p]     35 (the size of cc-pVDZ+TZpol)
    A41                   O [4s3p2d]    H [3s2p]     41 (the size of aug-cc-pVDZ)
    A61 = ANO-R2          O [4s4p2d1f]  H [3s2p1d]   61

    python debug/ano_matched_accuracy.py prep          # write the truncated basis files
    python debug/ano_matched_accuracy.py count         # function counts on water (control)
    python debug/ano_matched_accuracy.py ne <rung>     # control: Ne atom, non-scored
    python debug/ano_matched_accuracy.py water <rung>  # scored run (after registration)
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
BASIS = os.path.join(DATA, 'basis')
DL = os.path.join(HERE, 'lit_scan', 'dl_bse_samesize')
_spec = importlib.util.spec_from_file_location(
    'ma', os.path.join(HERE, 'matched_accuracy_gaussian.py'))
ma = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ma)
w = ma.w
ma.DATA = os.path.join(DATA, 'ano_matched_accuracy')
os.makedirs(ma.DATA, exist_ok=True)

# contracted functions kept per angular momentum (s, p, d, f), per element
RUNGS = {
    'A24': {'1': (2, 1), '8': (3, 2, 1), '10': (3, 2, 1)},
    'A35': {'1': (2, 2), '8': (3, 2, 2), '10': (3, 2, 2)},
    'A41': {'1': (3, 2), '8': (4, 3, 2), '10': (4, 3, 2)},
    'A61': {'1': (3, 2, 1), '8': (4, 4, 2, 1), '10': (4, 4, 2, 1)},
}
for rung in RUNGS:
    ma.FILES[rung] = 'ano-r2_%s_H_O_Ne.json' % rung


def truncate(data: dict, keep: dict) -> dict:
    out = {'elements': {}}
    for z, el in data['elements'].items():
        shells = []
        for sh in el['electron_shells']:
            l = sh['angular_momentum'][0]
            n = keep[z][l] if l < len(keep[z]) else 0
            if n == 0:
                continue
            assert n <= len(sh['coefficients']), (z, l, n)
            shells.append({**sh, 'coefficients': sh['coefficients'][:n]})
        out['elements'][z] = {'electron_shells': shells}
    return out


def prep() -> None:
    r2 = ma._clean(json.load(open(os.path.join(DL, 'ano-r2.json'))))
    for rung, keep in RUNGS.items():
        with open(os.path.join(BASIS, ma.FILES[rung]), 'w') as f:
            json.dump(truncate(r2, keep), f, indent=1)
        print('wrote', ma.FILES[rung])


def count() -> None:
    nuclei = w.water_nuclei()
    for rung in RUNGS:
        print(rung, ma.count(rung, nuclei))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'prep':
        prep()
    elif cmd == 'count':
        count()
    elif cmd == 'ne':
        ma.ne(sys.argv[2])
    elif cmd == 'water':
        ma.water(sys.argv[2])
