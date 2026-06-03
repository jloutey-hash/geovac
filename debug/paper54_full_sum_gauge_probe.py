"""
Paper 54 forward-scoping probe: FULL rotationally-invariant gauge field
========================================================================
Paper 54's eq:A_full claims the gauge field is the SUM over all multiplier
labels:  A = Σ_{NLM} (M†_{NLM} ⊗ M†_{NLM}) [D_total, (M_{NLM} ⊗ M_{NLM})].

But the actual driver tensor_product_gauged.py uses only the SINGLE
strongest commutator pair (comms[0] ⊗ comms[0]). The 75%/74% (memo) and
77%/32% (paper) connected-fraction numbers all come from single-pair runs.

Nuclear-sprint named follow-on #1: "Sum over ALL multiplier pairs in the
gauge field — should recover 100% Gaunt compatibility [and may close the
connected-fraction gap]."

THIS PROBE tests that hypothesis directly. We build the full-sum gauge
field and measure:
  (a) connected (non-factorizable) fraction of {D, A}
  (b) Gaunt compatibility (|Δl1| == |Δl2|) in the (+,+) chirality sector
  (c) m-conservation
  (d) multipole (k) decomposition

DECISION RELEVANCE:
  - If full-sum connected fraction stabilizes (does NOT crash toward 0
    with n_max) AND Gaunt compatibility reaches ~100%, the angular gap
    is closeable -> the gap is a single-pair artifact, not structural.
  - If connected fraction still collapses with n_max, the gap is
    structural (a second conformal-factor-type wall on the connected
    content itself).

NOTE: this is a DIAGNOSTIC probe (DO NOT modify production code/papers).
"""
import sys, json, time
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.operator_system import (
    build_multiplier_matrix, HyperLabel, allowed_multiplier_labels,
)
from debug.tensor_product_dirac import build_single_particle_scalar


def build_embedded_multipliers(sp, n_max):
    """Build all 3-Y multiplier matrices embedded in the chirality-doubled
    space. Returns list of dicts with label, M_full, and [D, M_full]."""
    D = sp['D']
    dim = sp['dim']
    N_scalar = sp['N_scalar']
    scalar_basis = sp['scalar_basis']

    out = []
    for (N, L, M) in allowed_multiplier_labels(n_max):
        M_sc = build_multiplier_matrix(N, L, M, scalar_basis).real
        if np.linalg.norm(M_sc) < 1e-15:
            continue
        M_full = np.zeros((dim, dim))
        M_full[:N_scalar, :N_scalar] = M_sc
        M_full[N_scalar:, N_scalar:] = M_sc
        comm = D @ M_full - M_full @ D
        out.append({
            'label': (N, L, M),
            'M_full': M_full,
            'comm': comm,
            'comm_norm': float(np.linalg.norm(comm)),
        })
    return out


def connected_fraction(DA, d1, d2):
    """Partial-trace subtraction: connected (non-factorizable) fraction."""
    I1, I2 = np.eye(d1), np.eye(d2)
    DA_t = DA.reshape(d1, d2, d1, d2)
    Tr2 = np.trace(DA_t, axis1=1, axis2=3)
    Tr1 = np.trace(DA_t, axis1=0, axis2=2)
    full_tr = np.trace(DA)
    factorizable = (np.kron(Tr2 / d2, I2)
                    + np.kron(I1, Tr1 / d1)
                    - (full_tr / (d1 * d2)) * np.eye(d1 * d2))
    connected = DA - factorizable
    den = np.linalg.norm(DA)
    return (float(np.linalg.norm(connected) / den) if den > 0 else 0.0,
            connected)


def build_full_sum_gauge(sp1, sp2, n_max):
    """Build A = Σ_lab (M†⊗M†)[D_total, M⊗M], the rotationally-invariant
    full-sum gauge field of Paper 54 eq:A_full."""
    d1, d2 = sp1['dim'], sp2['dim']
    D1, D2 = sp1['D'], sp2['D']
    gamma1 = sp1['gamma']

    D_total = np.kron(D1, np.eye(d2)) + np.kron(gamma1, D2)

    mults1 = build_embedded_multipliers(sp1, n_max)
    mults2 = build_embedded_multipliers(sp2, n_max)
    # Same operator system on both factors; pair label-by-label.
    assert len(mults1) == len(mults2)

    A = np.zeros((d1 * d2, d1 * d2))
    n_nonzero_comm = 0
    for m1, m2 in zip(mults1, mults2):
        M1, M2 = m1['M_full'], m2['M_full']
        c1, c2 = m1['comm'], m2['comm']  # [D1,M1], [D2,M2]
        if m1['comm_norm'] < 1e-12 and m2['comm_norm'] < 1e-12:
            continue
        n_nonzero_comm += 1
        # [D_total, M1⊗M2] = [D1,M1]⊗M2 + (γ1 M1)⊗[D2,M2]
        fluct = np.kron(c1, M2) + np.kron(gamma1 @ M1, c2)
        # conjugate pairing (M†⊗M†); matrices are real symmetric here
        pairing = np.kron(M1.T, M2.T)
        A += pairing @ fluct
    A = (A + A.T) / 2.0
    return D_total, A, n_nonzero_comm


def build_single_pair_gauge(sp1, sp2, n_max):
    """Reproduce the paper/memo single-strongest-pair gauge field."""
    d1, d2 = sp1['dim'], sp2['dim']
    D1, D2 = sp1['D'], sp2['D']
    gamma1 = sp1['gamma']
    D_total = np.kron(D1, np.eye(d2)) + np.kron(gamma1, D2)

    mults1 = build_embedded_multipliers(sp1, n_max)
    mults2 = build_embedded_multipliers(sp2, n_max)
    m1 = max(mults1, key=lambda x: x['comm_norm'])
    m2 = max(mults2, key=lambda x: x['comm_norm'])
    M1, M2 = m1['M_full'], m2['M_full']
    c1, c2 = m1['comm'], m2['comm']
    fluct = np.kron(c1, M2) + np.kron(gamma1 @ M1, c2)
    A = np.kron(M1, M2) @ fluct
    A = (A + A.T) / 2.0
    return D_total, A


def angular_audit(connected, scalar_basis):
    """Gaunt compatibility, m-conservation, multipole-k decomposition in
    the (+,+) chirality sector."""
    N = len(scalar_basis)
    dim_full = 2 * N
    V = np.zeros((N, N, N, N))
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    V[a, b, c, d] = connected[a * dim_full + b,
                                              c * dim_full + d]
    gaunt_ok = gaunt_bad = 0.0
    m_ok = m_bad = 0.0
    by_k = {}
    for a, ba in enumerate(scalar_basis):
        for b, bb in enumerate(scalar_basis):
            for c, bc in enumerate(scalar_basis):
                for d, bd in enumerate(scalar_basis):
                    val = V[a, b, c, d]
                    if abs(val) < 1e-14:
                        continue
                    dl1 = abs(ba.l - bc.l)
                    dl2 = abs(bb.l - bd.l)
                    if dl1 == dl2:
                        gaunt_ok += val ** 2
                        by_k[dl1] = by_k.get(dl1, 0.0) + val ** 2
                    else:
                        gaunt_bad += val ** 2
                    if ba.m + bb.m == bc.m + bd.m:
                        m_ok += val ** 2
                    else:
                        m_bad += val ** 2
    tot_g = gaunt_ok + gaunt_bad
    tot_m = m_ok + m_bad
    tot_k = sum(by_k.values())
    return {
        'gaunt_compat_pct': (gaunt_ok / tot_g * 100) if tot_g > 0 else None,
        'm_conserv_pct': (m_ok / tot_m * 100) if tot_m > 0 else None,
        'multipole_k': {int(k): (v / tot_k * 100) for k, v in sorted(by_k.items())}
                        if tot_k > 0 else {},
        'V_norm': float(np.linalg.norm(V)),
    }


def main():
    results = {}
    for n_max in [2, 3]:
        print('=' * 64)
        print(f'n_max = {n_max}')
        print('=' * 64)
        sp1 = build_single_particle_scalar(n_max)
        sp2 = build_single_particle_scalar(n_max)
        d1, d2 = sp1['dim'], sp2['dim']
        scalar_basis = sp1['scalar_basis']

        # --- single pair (reproduce paper/memo) ---
        D_total, A_sp = build_single_pair_gauge(sp1, sp2, n_max)
        DA_sp = D_total @ A_sp + A_sp @ D_total
        frac_sp, conn_sp = connected_fraction(DA_sp, d1, d2)
        ang_sp = angular_audit(conn_sp, scalar_basis)
        print(f'  [single pair]  connected fraction = {frac_sp*100:.1f}%')
        print(f'                 Gaunt = {ang_sp["gaunt_compat_pct"]}, '
              f'm-cons = {ang_sp["m_conserv_pct"]}, k = {ang_sp["multipole_k"]}')

        # --- full sum (Paper 54 eq:A_full) ---
        t0 = time.time()
        D_total, A_fs, n_nz = build_full_sum_gauge(sp1, sp2, n_max)
        DA_fs = D_total @ A_fs + A_fs @ D_total
        frac_fs, conn_fs = connected_fraction(DA_fs, d1, d2)
        ang_fs = angular_audit(conn_fs, scalar_basis)
        print(f'  [full sum]     ({n_nz} nonzero-comm generators summed)')
        print(f'                 connected fraction = {frac_fs*100:.1f}%')
        print(f'                 Gaunt = {ang_fs["gaunt_compat_pct"]}, '
              f'm-cons = {ang_fs["m_conserv_pct"]}, k = {ang_fs["multipole_k"]}')
        print(f'                 ||A_fullsum|| = {np.linalg.norm(A_fs):.4f}, '
              f'time = {time.time()-t0:.1f}s')

        results[n_max] = {
            'single_pair': {'connected_fraction': frac_sp, **ang_sp},
            'full_sum': {'connected_fraction': frac_fs,
                         'n_nonzero_generators': n_nz, **ang_fs},
        }

    with open('debug/data/paper54_full_sum_gauge_probe.json', 'w') as f:
        json.dump(results, f, indent=2)
    print('\nWrote debug/data/paper54_full_sum_gauge_probe.json')

    print('\n' + '=' * 64)
    print('GAP CHARACTERIZATION')
    print('=' * 64)
    sp_fracs = [results[n]['single_pair']['connected_fraction'] for n in [2, 3]]
    fs_fracs = [results[n]['full_sum']['connected_fraction'] for n in [2, 3]]
    print(f'  single-pair connected fraction n_max 2->3: '
          f'{sp_fracs[0]*100:.1f}% -> {sp_fracs[1]*100:.1f}%')
    print(f'  full-sum    connected fraction n_max 2->3: '
          f'{fs_fracs[0]*100:.1f}% -> {fs_fracs[1]*100:.1f}%')
    sp_gaunt = [results[n]['single_pair']['gaunt_compat_pct'] for n in [2, 3]]
    fs_gaunt = [results[n]['full_sum']['gaunt_compat_pct'] for n in [2, 3]]
    print(f'  single-pair Gaunt compat   n_max 2->3: '
          f'{sp_gaunt[0]:.1f}% -> {sp_gaunt[1]:.1f}%')
    print(f'  full-sum    Gaunt compat   n_max 2->3: '
          f'{fs_gaunt[0]:.1f}% -> {fs_gaunt[1]:.1f}%')


if __name__ == '__main__':
    main()
