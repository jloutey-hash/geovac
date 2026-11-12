"""Build fc water MO integrals (RHF canonical, O 1s frozen) for one Paper 62 Gaussian basis,
route of tests/test_paper62_matched_accuracy.py water_run.  Usage: python ... ccpvdz|anor1|a35"""
import sys, os, json, time, importlib.util, tempfile
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(ROOT, 'tests')
def load(name, file):
    sp = importlib.util.spec_from_file_location(name, os.path.join(T, file)); m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m
w = load('wcc', 'test_water_ccsd.py')
which = sys.argv[1]
if which == 'ccpvdz':
    path = os.path.join(T, 'data', 'paper62_cc_pvdz_H_O_Ne_2026_10_31.json')
elif which == 'anor1':
    path = os.path.join(T, 'data', 'paper62_ano_r1_H_O_Ne_2026_11_11.json')
elif which == 'a35':
    am = load('am', 'test_paper62_ano_matched_accuracy.py')
    data = json.load(open(am.R2_FILE))
    path = os.path.join(tempfile.gettempdir(), 'p62_a35_basis.json')
    json.dump(am.truncate(data, am.RUNGS['A35']), open(path, 'w'))
t0 = time.time()
nuclei = w.water_nuclei()
fns, Tm = w.gaussian_basis(path, nuclei)
S, h, g = w.ao_integrals(fns, nuclei, Tm)
enuc = w.vnn(nuclei)
Ehf, eps, C = w.rhf(S, h, g, 10, enuc)
hm = C.T @ h @ C
gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
hh, gg, ecore = w.freeze_core(hm, gm, 1)
# dipole-free localization needs AO data: keep S and C so the driver can localize later
st = w.jw_pauli_stats(hh, gg, enuc + ecore)
out = os.path.join(ROOT, 'debug', 'data', f'df_water_{which}_mo.npz')
np.savez(out, h1=hh, eri=gg, ecore=enuc + ecore, e_rhf=Ehf, nbf=Tm.shape[1], one_norm=st['one_norm'], n_pauli=st['n_pauli'])
print(which, 'nbf', Tm.shape[1], 'RHF', Ehf, 'raw 1-norm', st['one_norm'], 'pauli', st['n_pauli'], f'{time.time()-t0:.0f}s')
