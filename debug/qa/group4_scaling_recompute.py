"""group4 QA backfill — recompute the four Paper-14 scaling exponents from live data.

Paper 14 claims (sec:scaling): atomic N_Pauli ~ Q^3.15 (R^2 0.9995), 1-norm ~ Q^1.69
(R^2 0.997), QWC ~ Q^3.36 (R^2 1.000), over the He sweep n_max=2..5 (Q=10,28,60,110);
composed N_Pauli ~ Q^2.5 (composed_lih_scaling_sweep, vary max_n). The existing
test_fit_scaling only checks fit_pauli_scaling on SYNTHETIC data — no published exponent
is pinned from GeoVac data. This driver computes the real exponents so the pinning tests
assert the right numbers with sensible tolerance bands.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import JordanWignerEncoder, fit_pauli_scaling
from geovac.measurement_grouping import count_qwc_groups


def _l1_nonid(qop):
    """1-norm of non-identity terms of an OpenFermion QubitOperator."""
    s = 0.0
    for term, coeff in qop.terms.items():
        if term:  # () is identity
            s += abs(coeff)
    return float(s)


def _fit(qs, ys):
    qs = np.asarray(qs, float); ys = np.asarray(ys, float)
    a, _ = fit_pauli_scaling(qs, ys)
    # R^2 in log-log
    lq, ly = np.log(qs), np.log(ys)
    exp, logc = np.polyfit(lq, ly, 1)
    pred = exp * lq + logc
    ss_res = np.sum((ly - pred) ** 2); ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return float(a), float(r2)


def atomic_sweep(nmax_values=(2, 3, 4, 5), qwc_max_nmax=4):
    """N_Pauli + 1-norm at all nmax; QWC only for nmax <= qwc_max_nmax (greedy
    partition is O(n_terms^2) and intractable at nmax=5/Q=110)."""
    rows = []
    for nmax in nmax_values:
        t0 = time.time()
        li = LatticeIndex(n_electrons=2, max_n=nmax, nuclear_charge=2,
                          vee_method='slater_full', h1_method='hybrid')
        enc = JordanWignerEncoder(li)
        an = enc.analyze()
        qop = enc.build_qubit_operator()
        l1 = _l1_nonid(qop)
        t_build = time.time() - t0
        print(f"  He nmax={nmax}: Q={an.n_qubits} N_pauli={an.n_pauli_terms} "
              f"l1={l1:.3f} (build {t_build:.0f}s)", flush=True)
        qwc = None
        if nmax <= qwc_max_nmax:
            tq = time.time()
            qwc = count_qwc_groups(qop)
            print(f"      qwc={qwc} (qwc {time.time()-tq:.0f}s)", flush=True)
        rows.append({'nmax': nmax, 'Q': an.n_qubits, 'N_pauli': an.n_pauli_terms,
                     'l1': l1, 'qwc': qwc})
    return rows


def main():
    out = Path(__file__).resolve().parent / 'group4_scaling_recompute.json'
    result = {}

    print("=== ATOMIC sweep (He, n_max=2..5) ===")
    rows = atomic_sweep()
    Q = [r['Q'] for r in rows]
    a_p, r2_p = _fit(Q, [r['N_pauli'] for r in rows])
    a_l, r2_l = _fit(Q, [r['l1'] for r in rows])
    qrows = [r for r in rows if r['qwc'] is not None]
    a_q, r2_q = (_fit([r['Q'] for r in qrows], [r['qwc'] for r in qrows])
                 if len(qrows) >= 2 else (float('nan'), float('nan')))
    print(f"  (QWC fit over Q={[r['Q'] for r in qrows]})")
    print(f"\n  N_Pauli ~ Q^{a_p:.3f}  (R^2={r2_p:.4f})   [paper 3.15, R^2 0.9995]")
    print(f"  1-norm  ~ Q^{a_l:.3f}  (R^2={r2_l:.4f})   [paper 1.69, R^2 0.997]")
    print(f"  QWC     ~ Q^{a_q:.3f}  (R^2={r2_q:.4f})   [paper 3.36, R^2 1.000]")
    result['atomic'] = {'rows': rows,
                        'exp_pauli': a_p, 'r2_pauli': r2_p,
                        'exp_l1': a_l, 'r2_l1': r2_l,
                        'exp_qwc': a_q, 'r2_qwc': r2_q}

    print("\n=== COMPOSED LiH sweep (max_n=1..3) ===")
    try:
        from geovac.composed_qubit import composed_lih_scaling_sweep
        sw = composed_lih_scaling_sweep(max_n_values=[1, 2, 3], verbose=False)
        cQ = [d['Q'] for d in sw['sweep_data']]
        cP = [d['N_pauli'] for d in sw['sweep_data']]
        a_c, r2_c = _fit(cQ, cP)
        print(f"  data: {[(d['max_n'], d['Q'], d['N_pauli']) for d in sw['sweep_data']]}")
        print(f"  composed N_Pauli ~ Q^{a_c:.3f}  (R^2={r2_c:.4f})   [paper 2.5]")
        result['composed'] = {'Q': cQ, 'N_pauli': cP, 'exp': a_c, 'r2': r2_c,
                              'fitted_exp_reported': sw.get('alpha_fit')}
    except Exception as e:  # noqa: BLE001
        print(f"  composed sweep error: {type(e).__name__}: {e}")
        result['composed'] = {'error': f'{type(e).__name__}: {e}'}

    out.write_text(json.dumps(result, indent=2, default=float), encoding='utf-8')
    print(f"\nWrote {out}")


if __name__ == '__main__':
    main()
