"""DF lambda (geovac.df_lambda) for Paper 62 water Hamiltonians, fc, RHF canonical orbitals,
plus the truncation's fc-CCSD(T) error and an Edmiston-Ruedenberg localized-orbital check.
Usage: python debug/df_lambda_water.py dzp|ccpvdz|anor1|a35 [--loc]"""
import sys, os, json, time, importlib.util
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from geovac.df_lambda import df_lambda
sp = importlib.util.spec_from_file_location('wcc', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(sp); sp.loader.exec_module(w)
SCR = r'C:/Users/jlout/AppData/Local/Temp/claude/71f3482566352bd0/1be4f2b5-9390-4c7c-8d69-2264cf38a863/scratchpad/water_dzp_mo.npz'
which = sys.argv[1]
d = np.load(SCR if which == 'dzp' else os.path.join(ROOT, 'debug', 'data', f'df_water_{which}_mo.npz'))
h1, eri, ecore = d['h1'], d['eri'], float(d['ecore'])
N = h1.shape[0]; NOCC = 4
THS = [1e-2, 1e-3, 1e-4, 1e-5, 0.0]


def ccsdt(h, g):
    return float(w.run(np.eye(N), h, g, 2 * NOCC, ecore)['CCSD(T)'])


def er_localize(h, g, blocks, sweeps=60, tol=1e-10):
    """Edmiston-Ruedenberg: maximize sum_i (ii|ii) by Jacobi 2x2 rotations within each block."""
    h, g = h.copy(), g.copy()
    def rot(i, j, t):
        c, s = np.cos(t), np.sin(t)
        R = np.eye(N); R[i, i] = c; R[j, j] = c; R[i, j] = -s; R[j, i] = s   # columns = new orbitals
        return R
    def apply(h, g, i, j, t):
        c, s = np.cos(t), np.sin(t)
        for ax in range(4):
            gi = np.take(g, i, axis=ax).copy(); gj = np.take(g, j, axis=ax).copy()
            idx = [slice(None)] * 4
            idx[ax] = i; g[tuple(idx)] = c * gi + s * gj
            idx[ax] = j; g[tuple(idx)] = -s * gi + c * gj
        hi, hj = h[i].copy(), h[j].copy(); h[i] = c * hi + s * hj; h[j] = -s * hi + c * hj
        hi, hj = h[:, i].copy(), h[:, j].copy(); h[:, i] = c * hi + s * hj; h[:, j] = -s * hi + c * hj
        return h, g
    def f2(sub, t):
        c, s = np.cos(t), np.sin(t)
        u = np.array([c, s]); v = np.array([-s, c])
        return sum(np.einsum('pqrs,p,q,r,s->', sub, x, x, x, x) for x in (u, v))
    for sw in range(sweeps):
        gain = 0.0
        for blk in blocks:
            for a in range(len(blk)):
                for b in range(a + 1, len(blk)):
                    i, j = blk[a], blk[b]
                    sub = g[np.ix_([i, j], [i, j], [i, j], [i, j])]
                    f0, f8, f4 = f2(sub, 0.0), f2(sub, np.pi / 8), f2(sub, np.pi / 4)
                    cc = 0.5 * (f0 + f4); A = f0 - cc; B = f8 - cc
                    t = 0.25 * np.arctan2(B, A)
                    if abs(t) > 1e-12:
                        gain += f2(sub, t) - f0
                        h, g = apply(h, g, i, j, t)
        if gain < tol:
            break
    return h, g, sw


out = {'basis': which, 'n_active': N, 'raw_one_norm_canonical': float(w.jw_pauli_stats(h1, eri, 0.0)['one_norm'])}
t0 = time.time()
out['e_exact'] = ccsdt(h1, eri)
rows = []
for th in THS:
    r = df_lambda(h1, eri, ecore, threshold=th, return_eri=th > 0)
    row = {k: r[k] for k in ('threshold', 'lambda_total', 'lambda_one_body', 'lambda_two_body', 'n_leaves', 'rank', 'n_eigenvectors', 'max_abs_eri_error')}
    if th > 0:
        row['dE_ccsdt_mEh'] = 1e3 * (ccsdt(h1, r['eri_df']) - out['e_exact'])
    rows.append(row)
    print(which, row, f'{time.time()-t0:.0f}s', flush=True)
out['rows'] = rows
if '--loc' in sys.argv:
    hl, gl, nsw = er_localize(h1, eri, [list(range(NOCC)), list(range(NOCC, N))])
    out['er_sweeps'] = nsw
    out['er_sum_iiii_canonical'] = float(np.einsum('iiii->', eri)); out['er_sum_iiii_local'] = float(np.einsum('iiii->', gl))
    out['raw_one_norm_ER'] = float(w.jw_pauli_stats(hl, gl, 0.0)['one_norm'])
    out['df_lambda_ER_1e-3'] = df_lambda(hl, gl, threshold=1e-3, reconstruct=False)['lambda_total']
    out['df_lambda_ER_0'] = df_lambda(hl, gl, threshold=0.0, reconstruct=False)['lambda_total']
    out['e_hf_check_ER'] = float(w.run(np.eye(N), hl, gl, 2 * NOCC, ecore, triples=False)['RHF'])
    print({k: v for k, v in out.items() if k != 'rows'}, flush=True)
json.dump(out, open(os.path.join(ROOT, 'debug', 'data', f'df_lambda_water_{which}.json'), 'w'), indent=1)
