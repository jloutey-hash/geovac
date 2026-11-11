"""Same-size Gaussian control for Paper 62's water result (record SZ, 2026-11-11).

Paper 62 compares its 24-function Slater bases only with the cc-pV family. cc-pVDZ is built
for systematic convergence, not for the lowest energy at its size, so this run scores water at
fc-CCSD(T) in three Gaussian bases of the same size built for energy: pcseg-1 (Jensen),
def2-SVP (Weigend-Ahlrichs) and ANO-R1 (Zobel-Widmark-Veryazov). Same geometry, solver and
qubit counter as debug/matched_accuracy_gaussian.py.

    python debug/samesize_gaussian_control.py prep          # re-serialize the downloads
    python debug/samesize_gaussian_control.py count         # function counts on water (control)
    python debug/samesize_gaussian_control.py ccpvdz        # control: reproduce cc-pVDZ
    python debug/samesize_gaussian_control.py ne <basis>    # control: Ne atom, non-scored
    python debug/samesize_gaussian_control.py water <basis> # scored run (after registration)

Downloaded files (untrusted) live in debug/lit_scan/dl_bse_samesize/; `prep` keeps only
exponents and coefficients and writes them to debug/data/basis/.
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
# Rows go to debug/data/samesize_control/, not into the MA record's own jsonl files.
ma.DATA = os.path.join(DATA, 'samesize_control')
os.makedirs(ma.DATA, exist_ok=True)

NEW = {'pcseg1': 'pcseg-1', 'def2svp': 'def2-svp', 'anor1': 'ano-r1'}
for key, bse in NEW.items():
    ma.FILES[key] = bse + '_H_O_Ne.json'


def prep() -> None:
    for bse in NEW.values():
        obj = ma._clean(json.load(open(os.path.join(DL, bse + '.json'))))
        with open(os.path.join(BASIS, bse + '_H_O_Ne.json'), 'w') as f:
            json.dump(obj, f, indent=1)
        print('wrote', bse)


def count() -> None:
    nuclei = w.water_nuclei()
    for key in ('ccpvdz', *NEW):
        print(key, ma.count(key, nuclei))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'prep':
        prep()
    elif cmd == 'count':
        count()
    elif cmd == 'ccpvdz':
        ma.water('ccpvdz')
    elif cmd == 'ne':
        ma.ne(sys.argv[2])
    elif cmd == 'water':
        ma.water(sys.argv[2])
