"""Pairwise three-centre driver; engine + tests in tests/test_pairwise_three_center.py.

Usage: python debug/pairwise_three_center.py
"""
import importlib.util
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('pw', os.path.join(ROOT, 'tests', 'test_pairwise_three_center.py'))
pw = importlib.util.module_from_spec(_s)
_s.loader.exec_module(pw)
channels, GROUND_TRUTH = pw.channels, pw.GROUND_TRUTH

if __name__ == '__main__':
    base, cosg = channels(40)
    fine, _ = channels(40, nk=700, k_scale=4.0, nt=90, ht=0.05)
    print('cos(angle YXZ) = %.6f' % cosg)
    for a, b in zip(base, fine):
        if a['l'] in (0, 1, 2, 3, 4, 6, 8, 10, 14, 20, 30, 40):
            print('l<=%2d  partial %.12f  (fine grid %.12f)  term %+.3e  err vs GT %+.2e'
                  % (a['l'], a['partial_sum'], b['partial_sum'], a['term'], a['partial_sum'] - GROUND_TRUTH))
    with open(os.path.join(ROOT, 'debug', 'data', 'pairwise_three_center.json'), 'w') as f:
        json.dump({'cos_gamma': cosg, 'base': base, 'fine': fine, 'ground_truth': GROUND_TRUTH}, f, indent=1)
