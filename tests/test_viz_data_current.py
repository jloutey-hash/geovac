"""
Slow regeneration check: the committed viz/public/data artifacts match a
fresh ``geovac.viz_export.export_all`` run (the single-source-of-truth rule,
docs/visualization_plan.md §2/§3 — package changes that move an exported
number must fail loudly instead of leaving the site stale).

Comparison is tolerant on floats (eigensolvers and BLAS differ in final
ulps across platforms; the committed data is generated on the PI's machine,
CI runs Linux) and exact on everything else. If CI ever flakes on ordering
rather than values (e.g. an OpenFermion version change reordering terms),
loosen the comparison deliberately — do not delete the check.
"""

import json
import math
import os

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DATA_DIR = os.path.join(REPO, 'viz', 'public', 'data')

REL_TOL = 1e-6
ABS_TOL = 1e-9


def _compare(a, b, path=''):
    """Recursively compare committed vs regenerated JSON values."""
    problems = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a:
                problems.append(f'{path}.{k}: missing in committed')
            elif k not in b:
                problems.append(f'{path}.{k}: missing in regenerated')
            else:
                problems.extend(_compare(a[k], b[k], f'{path}.{k}'))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            problems.append(f'{path}: length {len(a)} != {len(b)}')
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                problems.extend(_compare(x, y, f'{path}[{i}]'))
    elif isinstance(a, float) or isinstance(b, float):
        if not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
            problems.append(f'{path}: type mismatch {type(a)} vs {type(b)}')
        elif not math.isclose(float(a), float(b), rel_tol=REL_TOL, abs_tol=ABS_TOL):
            problems.append(f'{path}: {a} != {b}')
    elif a != b:
        problems.append(f'{path}: {a!r} != {b!r}')
    return problems


@pytest.mark.slow
def test_committed_viz_data_matches_regeneration(tmp_path):
    pytest.importorskip('openfermion')
    from geovac.viz_export import export_all

    assert os.path.isdir(DATA_DIR), (
        'viz/public/data missing — run '
        'python -m geovac.viz_export --out viz/public/data'
    )
    manifest = export_all(str(tmp_path), verbose=False)

    committed_manifest = json.load(
        open(os.path.join(DATA_DIR, 'manifest.json'), encoding='utf-8'))
    assert committed_manifest['files'] == manifest['files']

    all_problems = []
    for rel in manifest['files']:
        fresh = json.load(open(os.path.join(str(tmp_path), *rel.split('/')),
                               encoding='utf-8'))
        committed = json.load(open(os.path.join(DATA_DIR, *rel.split('/')),
                                   encoding='utf-8'))
        probs = _compare(committed, fresh, rel)
        all_problems.extend(probs[:10])
    assert not all_problems, (
        f'{len(all_problems)}+ divergences between committed viz data and '
        f'fresh regeneration (first shown):\n' + '\n'.join(all_problems[:30]))
