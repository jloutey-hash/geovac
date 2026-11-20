"""Validate geovac.df_lambda against Lee et al. 2021 Table XIII (Reiher FeMoco) and openfermion's factorize."""
import sys, os, time, importlib.util, json
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geovac.df_lambda import df_lambda, one_body_tprime
import openfermion.resource_estimates as ofre
D = os.path.dirname(ofre.__file__)
spec = importlib.util.spec_from_file_location('of_fdf', os.path.join(D, 'df', 'factorize_df.py'))
fdf = importlib.util.module_from_spec(spec); spec.loader.exec_module(fdf)
with h5py.File(os.path.join(D, 'integrals', 'eri_reiher.h5'), 'r') as f:
    eri = np.asarray(f['eri'][()]); h1 = np.asarray(f['h0'][()]); ec = float(f['ecore'][()])
TABLE = [(0.1, 53, 765, 253.6), (0.0125, 182, 3983, 280.8), (0.00125, 360, 13031, 294.8), (0.0001, 581, 25054, 296.7)]
lamT = np.abs(np.linalg.eigvalsh(one_body_tprime(h1, eri))).sum()
print('lambda_T prime', lamT, '(Lee: 38.6)')
res = []
for th, L, Xi, lam in TABLE:
    t = time.time(); r = df_lambda(h1, eri, threshold=th)
    _, F, rk, nv = fdf.factorize(eri, th)
    lof = 0.0
    for k in range(F.shape[2]):
        e = np.linalg.eigvalsh(F[:, :, k]); lof += 0.25 * np.abs(e).sum() ** 2
    lof += lamT
    print(th, 'own', round(r['lambda_total'], 3), r['n_leaves'], r['n_eigenvectors'], '| of', round(lof, 3), F.shape[2], nv, '| Lee', lam, L, Xi, '| rel', abs(lof - r['lambda_total']) / lof, f'{time.time()-t:.1f}s')
    res.append(dict(th=th, own=r['lambda_total'], own_L=r['n_leaves'], own_Xi=r['n_eigenvectors'], of=lof, of_L=int(F.shape[2]), of_Xi=int(nv), lee=lam, lee_L=L, lee_Xi=Xi))
json.dump(dict(lambda_Tprime=float(lamT), rows=res), open(os.path.join(os.path.dirname(__file__), 'data', 'df_lambda_reiher_check.json'), 'w'), indent=1)
