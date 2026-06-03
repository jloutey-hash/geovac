"""
Paper 54 forward-scoping probe 2: does the FULL-SUM gauge field's connected
interaction match Coulomb RADIAL weights?
===========================================================================
Probe 1 showed the full-sum gauge field recovers ~97% Gaunt compatibility
(closing the angular gap) but leaves a ~15-30% non-connected residual.

The DECISION GATE is NOT the connected fraction (a partial-trace bookkeeping
quantity). It is whether the connected part's RADIAL matrix-element weights
are proportional to the Coulomb 1/r12 kernel. Paper 54 Prop 3 (radial
negative) used the single-pair gauge field: Pearson 0.58 (n=2), 0.41 (n=3),
CV > 1600%, DECREASING with n_max.

This probe recomputes the Coulomb-vs-spectral radial correlation for the
FULL-SUM gauge field. If the full sum lifts Pearson toward 1.0 and CV
toward 0, the radial gap is closeable. If Pearson stays low / decreases
with n_max, the radial gap is a second conformal-factor wall (matching the
resolvent CLEAN NEGATIVE).

Coulomb reference uses the same chordal-distance Green's function on S3
that the resolvent sprint used: the 1/(N^2-1) Laplacian-resolvent kernel is
NOT used here; instead we use the exact angular Gaunt structure with a
flat-space-like radial reference, identical in spirit to the resolvent
sprint's "exact Coulomb Slater" comparison. To stay self-contained and
avoid importing the full hypergeometric Slater machinery, we use the
NEUMANN/Green reference V_coul[a,b,c,d] = sum_NLM M_NLM[a,c] conj(M_NLM[d,b])
/ (N^2-1)  (Laplacian Green's function on S3, the chordal-distance kernel
of Paper 7) restricted to the SAME support as the spectral interaction.

This is the Paper-7 chordal-distance / Fock-projection reference. The point
of the probe is the CORRELATION trend with n_max, not the absolute kernel.
"""
import sys, json
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.operator_system import (
    build_multiplier_matrix, allowed_multiplier_labels,
)
from debug.tensor_product_dirac import build_single_particle_scalar
from debug.paper54_full_sum_gauge_probe import (
    build_embedded_multipliers, build_full_sum_gauge, connected_fraction,
)


def coulomb_reference(sp, n_max):
    """Chordal-distance Green's function kernel on S3 (Paper 7 / resolvent
    sprint reference): V[a,b,c,d] = sum_NLM M[a,c] conj(M[d,b]) / (N^2-1),
    N>=2 (zero mode excluded). Built in the SCALAR (un-doubled) basis."""
    scalar_basis = sp['scalar_basis']
    N = len(scalar_basis)
    V = np.zeros((N, N, N, N))
    for (Nq, Lq, Mq) in allowed_multiplier_labels(n_max):
        denom = Nq * Nq - 1.0
        if denom <= 0:
            continue  # zero mode N=1 excluded (Laplacian Green's fn)
        Msc = build_multiplier_matrix(Nq, Lq, Mq, scalar_basis).real
        if np.linalg.norm(Msc) < 1e-15:
            continue
        # V[a,b,c,d] += M[a,c] * M[d,b] / (N^2-1)
        # (real symmetric multipliers, conj is identity)
        V += np.einsum('ac,db->abcd', Msc, Msc) / denom
    return V


def extract_pp_tensor(connected, N):
    """(+,+) chirality sector 4-index tensor."""
    dim_full = 2 * N
    V = np.zeros((N, N, N, N))
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    V[a, b, c, d] = connected[a * dim_full + b,
                                              c * dim_full + d]
    return V


def correlate(V_spec, V_coul):
    """Pearson, CV of ratio, sign agreement, coverage — on the support
    where the spectral interaction is nonzero."""
    s = V_spec.ravel()
    c = V_coul.ravel()
    mask = np.abs(s) > 1e-12
    n_spec = int(mask.sum())
    n_coul = int((np.abs(c) > 1e-12).sum())
    # overlap support
    both = mask & (np.abs(c) > 1e-12)
    n_both = int(both.sum())
    if n_both < 2:
        return {'pearson': None, 'cv_ratio': None, 'sign_agree': None,
                'n_spec_nonzero': n_spec, 'n_coul_nonzero': n_coul,
                'n_overlap': n_both}
    ss, cc = s[both], c[both]
    pear = float(np.corrcoef(ss, cc)[0, 1])
    ratio = ss / cc
    cv = float(np.std(ratio) / np.abs(np.mean(ratio))) if np.mean(ratio) != 0 else None
    sign = float(np.mean(np.sign(ss) == np.sign(cc)))
    return {'pearson': pear, 'cv_ratio': cv, 'sign_agree': sign,
            'n_spec_nonzero': n_spec, 'n_coul_nonzero': n_coul,
            'n_overlap': n_both}


def main():
    results = {}
    for n_max in [2, 3]:
        sp1 = build_single_particle_scalar(n_max)
        sp2 = build_single_particle_scalar(n_max)
        N = sp1['N_scalar']
        d1, d2 = sp1['dim'], sp2['dim']

        # full-sum gauge connected interaction
        D_total, A_fs, n_nz = build_full_sum_gauge(sp1, sp2, n_max)
        DA = D_total @ A_fs + A_fs @ D_total
        _, conn = connected_fraction(DA, d1, d2)
        V_spec = extract_pp_tensor(conn, N)

        V_coul = coulomb_reference(sp1, n_max)

        corr = correlate(V_spec, V_coul)
        results[n_max] = {'n_nonzero_generators': n_nz, **corr}
        print(f'n_max={n_max}  [full-sum gauge vs chordal Coulomb]')
        print(f'  Pearson = {corr["pearson"]}, CV(ratio) = {corr["cv_ratio"]}, '
              f'sign-agree = {corr["sign_agree"]}')
        print(f'  support: spec_nonzero={corr["n_spec_nonzero"]}, '
              f'coul_nonzero={corr["n_coul_nonzero"]}, overlap={corr["n_overlap"]}')

    with open('debug/data/paper54_radial_under_fullsum.json', 'w') as f:
        json.dump(results, f, indent=2)
    print('\nWrote debug/data/paper54_radial_under_fullsum.json')

    print('\n' + '=' * 56)
    print('RADIAL GAP TREND (the GO/STOP-deciding metric)')
    print('=' * 56)
    p = [results[n]['pearson'] for n in [2, 3]]
    print(f'  full-sum Pearson  n_max 2->3: {p[0]:.3f} -> {p[1]:.3f}  '
          f'({"DECREASING -> structural wall" if p[1] < p[0] else "increasing"})')
    print('  (Paper 54 single-pair Prop 3: 0.58 -> 0.41, also decreasing)')
    print('  (Resolvent sprint Dirac-weight: 0.81 -> 0.75, also decreasing)')


if __name__ == '__main__':
    main()
