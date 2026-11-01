"""Guard for the numpy 2.4.0 flat-assignment memory leak (CHANGELOG v5.21.25, v5.21.29).

numpy 2.4.0 leaks about 2.8 KB per call, never freed, whenever an array is written through
``ndarray.flat`` with a slice or an index array.  np.diag, np.eye, np.identity, np.diagflat
and np.fill_diagonal all write that way, and geovac calls them at about 130 sites, so a
long optimizer run grew to 18 GB.  Measured 2026-10-31 on Python 3.14.0 in isolated
environments: 2.4.0 leaks; 2.4.1, 2.4.6 and 2.5.4 do not.

The fix is to exclude the release, not to patch every call site: setup.py requires
numpy != 2.4.0, and ``import geovac`` warns if 2.4.0 is already installed.

Written as a separate step after the fix (CLAUDE.md Sec. 9, guard-writing rule).
"""
import gc
import os
import re
import warnings

import numpy as np
import psutil

import geovac

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CALLS = 20000
# A leaking numpy grows ~2.8 KB per call: about 55 MB over CALLS.  A clean one grows < 1 MB.
LEAK_MB = 10.0


def _rss_growth_mb(f, calls=CALLS):
    proc = psutil.Process(os.getpid())
    for _ in range(200):
        f()
    gc.collect()
    a = proc.memory_info().rss
    for _ in range(calls):
        f()
    gc.collect()
    return (proc.memory_info().rss - a) / 1e6


def test_setup_excludes_the_leaking_release():
    """The install requirement for numpy excludes exactly 2.4.0 and keeps its floor.

    Wrong answers rejected: the exclusion dropped (a fresh install can pull 2.4.0), and
    an exclusion of a different release that leaves 2.4.0 allowed.
    """
    src = open(os.path.join(ROOT, 'setup.py'), encoding='utf-8').read()
    m = re.search(r"'numpy([^']*)'", src)
    assert m, "no numpy requirement in setup.py"
    spec = m.group(1).replace(' ', '')
    assert '!=2.4.0' in spec.split(','), spec
    assert '>=1.20.0' in spec.split(','), spec


def test_import_warning_fires_on_2_4_0_only(monkeypatch):
    """The import-time check warns on 2.4.0 and stays silent on the fixed releases.

    Wrong answers rejected: an empty or wrong version set (never warns), an inverted test
    (warns on every clean install), and a warning of a category users filter by default.
    """
    for ver, expect in (('2.4.0', True), ('2.4.1', False), ('2.4.6', False), ('2.5.4', False)):
        monkeypatch.setattr(np, '__version__', ver)
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter('always')
            got = geovac._warn_if_leaky_numpy()
        assert got is expect, (ver, got)
        assert len(rec) == int(expect), (ver, [str(r.message) for r in rec])
        if expect:
            assert rec[0].category is RuntimeWarning
            assert 'numpy>=2.4.1' in str(rec[0].message)


def test_installed_numpy_does_not_leak_in_flat_assignment():
    """The numpy actually installed does not leak in the five functions that write through
    ``.flat``, and the measurement can see a leak of that size (positive control).

    Fails, correctly, on numpy 2.4.0: measured 54 to 57 MB growth per form over 20,000 calls.
    Wrong answers rejected: a leaking numpy left in the environment, and a probe too
    coarse to see a 2.8 KB-per-call leak.
    """
    held = []
    control = _rss_growth_mb(lambda: held.append(bytearray(2800)))
    assert control > LEAK_MB, control          # the probe sees a leak of this size
    del held
    v = np.random.rand(50)
    z = np.zeros((50, 50))
    forms = {
        'np.diag': lambda: np.diag(v),
        'np.eye': lambda: np.eye(50),
        'np.identity': lambda: np.identity(50),
        'np.diagflat': lambda: np.diagflat(v),
        'np.fill_diagonal': lambda: np.fill_diagonal(z, v),
    }
    grown = {k: _rss_growth_mb(f) for k, f in forms.items()}
    bad = {k: round(g, 1) for k, g in grown.items() if g > LEAK_MB}
    assert not bad, f"numpy {np.__version__} leaks (MB over {CALLS} calls): {bad}"
