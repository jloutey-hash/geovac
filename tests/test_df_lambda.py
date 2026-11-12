"""geovac.df_lambda: double-factorization 1-norm (Lee et al. 2021, arXiv:2011.03494, App. C).

Three independent checks.
1. Published reference: Lee et al. Table XIII (Reiher FeMoco, 54 orbitals), rows at
   thresholds 0.1, 0.00125 and 0.0001: lambda = 253.6, 294.8, 296.7; L = 53, 360, 581;
   eigenvectors Xi = 765, 13031, 25054; lambda_T' = 38.6 (caption).  Integrals from the
   eri_reiher.h5 file shipped with openfermion.resource_estimates (needs h5py).
2. Second route: openfermion.resource_estimates.df.factorize (loaded from its file, since
   the package __init__ hides it without pyscf/jax) with the lambda formula of its
   compute_lambda_df.py transcribed; agreement to 1e-8 relative on STO-3G water
   (frozen-core and all-electron) and an STO-3G H4 rectangle, at threshold 0 and 1e-3.
3. Spectral bound (Loaiza et al. 2023): any LCU 1-norm of H, constant shifts aside, is at
   least half the spectral range of H over the whole Fock space; checked on frozen-core
   STO-3G water (12 qubits) by sparse diagonalization.

Written as a separate step after the tool was validated (CLAUDE.md Sec. 9); each test
names the wrong answer it rejects.
"""
import importlib.util
import os

import numpy as np
import pytest

from geovac.df_lambda import df_lambda, one_body_tprime

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_df', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)


def _of_dir():
    import openfermion.resource_estimates as ofre
    return os.path.dirname(ofre.__file__)


def _sto3g_mo(nuc, nelec, fc):
    fns, T = w.gaussian_basis(w._basis('sto-3g_H_O.json'), nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    _, _, C = w.rhf(S, h, g, nelec, w.vnn(nuc))
    hm = C.T @ h @ C
    gm = np.einsum('pqrs,pi,qj,rk,sl->ijkl', g, C, C, C, C, optimize=True)
    return w.freeze_core(hm, gm, fc)


H4_NUC = [(np.array([x, y, 0.0]), 1.0) for x in (0.0, 1.8) for y in (0.0, 2.6)]


def _openfermion_lambda(h1, eri, thresh):
    spec = importlib.util.spec_from_file_location('of_factorize_df', os.path.join(_of_dir(), 'df', 'factorize_df.py'))
    fdf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fdf)
    _, F, _, _ = fdf.factorize(eri, thresh)
    # compute_lambda_df.py, transcribed (its import needs pyscf)
    T = h1 - 0.5 * np.einsum('illj->ij', eri) + np.einsum('llij->ij', eri)
    lam = np.abs(np.linalg.eigvalsh(T)).sum()
    for k in range(F.shape[2]):
        lam += 0.25 * np.abs(np.linalg.eigvalsh(F[:, :, k])).sum() ** 2
    return float(lam)


def test_reiher_femoco_reproduces_lee_table_xiii():
    """Lee et al. Table XIII to its printed precision, with L and Xi exact.

    Wrong answers rejected: the 1/4 of Eq. (C11) dropped or doubled; the two-body
    correction sum_r (pq|rr) of Eq. (A9) omitted from T' (lambda_T' 38.6 would move);
    the -1/2 sum_r (pr|rq) of Eq. (3) omitted; the C41 cut applied to |f| alone instead of
    (sum|f|)|f| (Xi would move); not stopping at the first empty leaf (L would move).
    """
    h5py = pytest.importorskip('h5py', reason='h5py needed to read the Reiher integrals')
    path = os.path.join(_of_dir(), 'integrals', 'eri_reiher.h5')
    if not os.path.exists(path):
        pytest.skip('eri_reiher.h5 not shipped with this openfermion')
    with h5py.File(path, 'r') as f:
        eri = np.asarray(f['eri'][()])
        h1 = np.asarray(f['h0'][()])
    lam_t = np.abs(np.linalg.eigvalsh(one_body_tprime(h1, eri))).sum()
    assert abs(lam_t - 38.6) < 0.05, lam_t
    for th, L, Xi, lam in [(0.1, 53, 765, 253.6), (0.00125, 360, 13031, 294.8), (0.0001, 581, 25054, 296.7)]:
        r = df_lambda(h1, eri, threshold=th, reconstruct=False)
        assert r['n_leaves'] == L, (th, r['n_leaves'])
        assert r['n_eigenvectors'] == Xi, (th, r['n_eigenvectors'])
        assert abs(r['lambda_total'] - lam) < 0.05, (th, r['lambda_total'])


@pytest.mark.parametrize('system', ['water_fc', 'water_ae', 'h4'])
@pytest.mark.parametrize('thresh', [0.0, 1e-3])
def test_two_routes_agree(system, thresh):
    """Own packed-pair-space route vs openfermion's full-supermatrix factorize, 1e-8 relative.

    Wrong answers rejected: a packing weight other than sqrt 2 on off-diagonal pairs
    (leaves rescaled), leaves taken in ascending order (truncation stops early), and any
    one-body or 1/4 error confined to one route.
    """
    if system == 'h4':
        h1, eri, _ = _sto3g_mo(H4_NUC, 4, 0)
    else:
        h1, eri, _ = _sto3g_mo(w.CRAWFORD_NUC, 10, 1 if system == 'water_fc' else 0)
    own = df_lambda(h1, eri, threshold=thresh)['lambda_total']
    ref = _openfermion_lambda(h1, eri, thresh)
    assert abs(own - ref) <= 1e-8 * ref, (own, ref)


def test_lambda_at_least_half_the_fock_space_spectral_range():
    """Frozen-core STO-3G water: lambda_DF (untruncated) = 16.50 >= (E_max - E_min)/2 = 9.78,
    over all particle numbers (12 qubits); untruncated DF reconstructs the ERIs exactly.

    Wrong answers rejected: lambda reported as the two-body part alone (7.85), the
    one-body Schatten norm dropped, and a lambda computed from a truncated tensor passed
    off as exact (reconstruction error checked).
    """
    from openfermion import InteractionOperator, get_sparse_operator
    import scipy.sparse.linalg as sla
    ha, ga, _ = _sto3g_mo(w.CRAWFORD_NUC, 10, 1)
    r = df_lambda(ha, ga, threshold=0.0)
    assert r['max_abs_eri_error'] < 1e-12
    M = ha.shape[0]
    n = 2 * M
    sp, sg = np.arange(n) // 2, np.arange(n) % 2
    same = (sg[:, None] == sg[None, :]).astype(float)
    one = ha[sp][:, sp] * same
    two = 0.5 * ga[np.ix_(sp, sp, sp, sp)].transpose(0, 2, 3, 1) * same[:, None, None, :] * same[None, :, :, None]
    H = get_sparse_operator(InteractionOperator(0.0, one, two)).tocsc()
    lo = sla.eigsh(H, k=1, which='SA')[0][0]
    hi = sla.eigsh(H, k=1, which='LA')[0][0]
    half = 0.5 * (hi - lo)
    assert abs(half - 9.776) < 0.01, half
    assert r['lambda_total'] >= half, (r['lambda_total'], half)
    assert r['lambda_two_body'] < half          # the two-body part alone would violate it


def test_df_lambda_is_invariant_under_active_orbital_rotation():
    """The supermatrix moves by O (x) O, so lambda_DF and lambda_T' are rotation invariant
    while the raw JW Pauli 1-norm is not.  Rejects an implementation whose leaves depend
    on the orbital frame (e.g. a pivoted Cholesky with frame-dependent pivots, truncated)."""
    ha, ga, _ = _sto3g_mo(w.CRAWFORD_NUC, 10, 1)
    O, _ = np.linalg.qr(np.random.default_rng(7).standard_normal(ha.shape))
    hr = O.T @ ha @ O
    gr = np.einsum('pqrs,pa,qb,rc,sd->abcd', ga, O, O, O, O, optimize=True)
    a = df_lambda(ha, ga, threshold=1e-4)
    b = df_lambda(hr, gr, threshold=1e-4)
    assert abs(a['lambda_total'] - b['lambda_total']) < 1e-9 * a['lambda_total']
    raw_a = w.jw_pauli_stats(ha, ga, 0.0)['one_norm']
    raw_b = w.jw_pauli_stats(hr, gr, 0.0)['one_norm']
    assert abs(raw_a - raw_b) > 1.0, (raw_a, raw_b)
