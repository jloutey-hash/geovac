"""Record PP: launch every exact point the scorer needs, NWORK at a time; skips points already on disk.
    python debug/pp_launch.py [NWORK]
Point names: <geom>_<d1>_<d2>_<a>_<b>_<step>, d = offset in D1, D2 and a, b = offset in z1, z2, each in {m, 0, p}
(minus h, none, plus h), step 30 or 31 (h = 10^-step).
"""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).parent
GEOM = {'U': (Fraction(9, 5), Fraction(12, 5)), 'E': (Fraction(19, 10), Fraction(19, 10))}
ZX, Z = Fraction(13, 10), Fraction(6, 5)
OFF = {'m': -1, '0': 0, 'p': 1}
PM, ZZ = ('m', 'p'), ('m', '0', 'p')


def points():
    pts = set()
    for g in GEOM:
        for d1 in PM:                                   # PA: d/dD1, O_1 on z1
            for a in ZZ:
                pts.add((g, d1, '0', a, '0', 30))
        for d1 in PM:                                   # PB: d/dD1 d/dD2, O_1 O_2
            for d2 in PM:
                for a in ZZ:
                    for b in ZZ:
                        pts.add((g, d1, d2, a, b, 30))
    for a in ZZ:                                        # PC at U: O_1 O_2, no D derivative (L = 150 set,
        for b in ZZ:                                    # superseded for PC by Addendum 1)
            pts.add(('U', '0', '0', a, b, 30))
            pts.add(('U', '0', '0', a, b, 30, 230))     # Addendum 1: PC inputs to L = 230
    pts.add(('E', '0', '0', '0', '0', 30))              # PP-0(b) 1s control at E (U's is in the PC set)
    for d1 in PM:                                       # PP-0(a) recheck of PA-U
        for a in ZZ:
            pts.add(('U', d1, '0', a, '0', 31))
    return sorted(pts, key=lambda p: (len(p), p))


def run(p):
    g, d1, d2, a, b, st = p[:6]
    lm = p[6] if len(p) > 6 else 150
    name = '_'.join(map(str, p[:6])) + ('' if lm == 150 else f'_L{lm}')
    if (HERE / 'data' / 'pp' / f'{name}.json').exists():
        return name + ' (on disk)'
    h = Fraction(1, 10 ** st)
    D1, D2 = GEOM[g]
    args = [D1 + OFF[d1] * h, D2 + OFF[d2] * h, ZX, Z + OFF[a] * h, Z + OFF[b] * h]
    r = subprocess.run([sys.executable, str(HERE / 'pp_run.py'), name] + [str(x) for x in args] + [str(lm)],
                       cwd=HERE, capture_output=True, text=True)
    return (r.stdout.strip() or name) + ('' if r.returncode == 0 else ' FAILED ' + r.stderr[-300:])


if __name__ == '__main__':
    pts = points()
    print(len(pts), 'points', flush=True)
    with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 10) as ex:
        for msg in ex.map(run, pts):
            print(msg, flush=True)
    print('done', flush=True)
