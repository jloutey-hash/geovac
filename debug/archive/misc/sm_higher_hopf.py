"""
Track SM-F: Higher Hopf fibrations and Standard Model couplings.

Paper 2's construction for U(1)_em uses the complex Hopf bundle
    S^1 -> S^3 -> S^2.
The spectral invariants come from
    - Base S^2:   B = sum over Fock shells (n,l) of (2l+1)*l(l+1) = 42 at n_max=3
    - Fiber S^1:  F = zeta_{S^1}(1) = 2*zeta_R(2) = pi^2/3  but paper uses F = zeta_R(2) = pi^2/6
                  (redefined as Dirichlet series D_{n^2}(d_max=4) = zeta_R(2))
    - Total S^3:  Delta = 1/(|lambda_{n_max}| * N(n_max-1)) = 1/(8*5) = 1/40

Combination:  K = pi (B + F - Delta),  alpha = 1/K.

Selection principle (paper):  B/N = dim(S^3) = 3, which is the TOTAL SPACE
dimension, not the base.  This uniquely selects n_max = 3 via a closed-form
rational ratio.

This script extends the construction to the quaternionic and octonionic
Hopf fibrations:
    Quaternionic:   S^3 -> S^7  -> S^4   (dim S^7 = 7,  dim S^4 = 4)
    Octonionic:     S^7 -> S^15 -> S^8   (dim S^15 = 15, dim S^8 = 8)

Strategy:
1. Compute Laplacian spectrum on the BASE (S^4, S^8) with its natural
   degeneracies.
2. Try B/N equal to various natural targets:
      dim(total space):  7 (quat), 15 (oct)
      dim(base):         4,         8
      dim(fiber):        3,         7
      dim division alg.: 4 (H),     8 (O)
   Look for the unique n_max at which the closed-form ratio hits.
3. For the fiber F, compute the Dirichlet series of the base degeneracy
   at the natural exponent (packing dimension), or the fiber spectral zeta.
4. Combine:  K = pi (B + F - Delta),  alpha = 1/K,
   and compare with  alpha_W ~ 1/30,  alpha_s(M_Z) ~ 0.118, etc.
5. Also try B/N targets motivated by SU(5) GUT:  sin^2 theta_W = 3/8 ->
   alpha_W/alpha_em = 8/3.

All computations in exact rational/symbolic form via sympy, mpmath at 50 digits
as a fallback.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 50

# ---------------------------------------------------------------------------
# Degeneracy formulas for S^N Laplacian eigenvalues
# ---------------------------------------------------------------------------
#   S^N Laplacian eigenvalue:   l*(l + N - 1),   l = 0, 1, 2, ...
#   Degeneracy on S^N:
#       dim V_l^N = (2l + N - 1) * (l + N - 2)! / (l! * (N - 1)!)
# For N = 2:  2l+1
# For N = 3:  (l+1)^2
# For N = 4:  (l+1)(l+2)(2l+3)/6
# For N = 5:  (l+1)(l+2)^2(l+3)/12
# For N = 8:  full expression from the formula above.

l = sp.symbols('l', integer=True, nonnegative=True)


def deg_SN(N):
    """Symbolic degeneracy of the l-th spherical harmonic space on S^N."""
    if N == 0:
        return sp.Piecewise((2, sp.Eq(l, 0)), (0, True))
    if N == 1:
        return sp.Piecewise((1, sp.Eq(l, 0)), (2, True))
    # general formula: (2l + N - 1)/(N-1) * C(l + N - 2, N - 2)
    return sp.simplify(
        (2*l + N - 1) / (N - 1) * sp.binomial(l + N - 2, N - 2)
    )


def eig_SN(N):
    """Eigenvalue l(l + N - 1) of the Laplacian on S^N."""
    return l * (l + N - 1)


# Verify low-N degeneracy formulas against known cases
def verify_degeneracies():
    checks = {}
    # S^2: expected 2l+1
    d2 = sp.simplify(deg_SN(2) - (2*l + 1))
    checks['S2_deg_ok'] = (d2 == 0)
    # S^3: expected (l+1)^2
    d3 = sp.simplify(deg_SN(3) - (l + 1)**2)
    checks['S3_deg_ok'] = (d3 == 0)
    # S^4: (l+1)(l+2)(2l+3)/6
    d4 = sp.simplify(deg_SN(4) - (l + 1)*(l + 2)*(2*l + 3)/6)
    checks['S4_deg_ok'] = (d4 == 0)
    # S^5: (l+1)(l+2)^2(l+3)/12
    d5 = sp.simplify(deg_SN(5) - (l + 1)*(l + 2)**2*(l + 3)/12)
    checks['S5_deg_ok'] = (d5 == 0)
    return checks


# ---------------------------------------------------------------------------
# Fock shell structure for the higher Hopf bundles
# ---------------------------------------------------------------------------
# Paper 2 uses the S^3 Fock shell n = 1, 2, 3 with per-shell angular momenta
# l = 0, ..., n-1. This is the S^3 Peter-Weyl decomposition of L^2(S^3) into
# Fock hydrogenic shells, and the inner sum over l reproduces the S^2 base
# Casimir per shell.
#
# For the higher Hopf bundles we do NOT have a hydrogenic Fock projection
# (that is exactly the content of Paper 24's HO rigidity theorem - only
# Coulomb/S^3 has this). So the "natural" analog is a direct cutoff on the
# base spherical harmonic level l.
#
# Two reasonable interpretations:
# (a) Direct cutoff:  truncate base spherical harmonics at l = L_max.
# (b) Shell cutoff:   use a multi-shell structure with inner angular
#                     quantum numbers, mimicking Paper 2's (n,l).
#
# Interpretation (a) is the minimal, unambiguous choice.  Interpretation
# (b) requires picking a hydrogen-like shell structure which is not
# canonical for non-Coulomb bases.  We use (a) as the primary and
# mention (b) qualitatively.

def B_base(N_base, L_max):
    """
    Degeneracy-weighted Casimir trace of the base S^{N_base}.

        B(L_max) = sum_{l=0}^{L_max} deg_SN(N_base)(l) * eig_SN(N_base)(l)

    This is the direct base analog of Paper 2's B = sum (2l+1) l(l+1) summed
    over l = 0..2 (the S^2 Casimirs contained in the S^3 Fock shells n<=3).
    """
    total = sp.Integer(0)
    d_expr = deg_SN(N_base)
    e_expr = eig_SN(N_base)
    for L in range(L_max + 1):
        total += d_expr.subs(l, L) * e_expr.subs(l, L)
    return sp.nsimplify(total)


def N_base(N_base, L_max):
    """Total state count on S^{N_base} through level L_max (sum of degeneracies)."""
    total = sp.Integer(0)
    d_expr = deg_SN(N_base)
    for L in range(L_max + 1):
        total += d_expr.subs(l, L)
    return sp.nsimplify(total)


# ---------------------------------------------------------------------------
# Selection principle search
# ---------------------------------------------------------------------------
def find_selection_cutoff(N_base_dim, targets, L_max_search=30):
    """
    For each target ratio t in `targets`, find all L_max <= L_max_search
    at which B(L_max) / N(L_max) exactly equals t (if any).

    Also tabulate the actual ratio vs L_max.
    """
    table = []
    for L in range(1, L_max_search + 1):
        B = B_base(N_base_dim, L)
        N_ = N_base(N_base_dim, L)
        ratio = sp.Rational(B, N_)
        table.append((L, int(B), int(N_), ratio, float(ratio)))
    hits = {float(t): [] for t in targets}
    for L, B, N_, ratio, rf in table:
        for t in targets:
            if ratio == sp.nsimplify(t):
                hits[float(t)].append(L)
    return table, hits


# ---------------------------------------------------------------------------
# Fiber analog F (Dirichlet series of degeneracy)
# ---------------------------------------------------------------------------
# Paper 2: F = D_{n^2}(s = d_max = 4) = sum_{n>=1} n^2 * n^{-4} = zeta_R(2).
# Here the degeneracy weight is the S^3 Fock shell degeneracy g_n = n^2
# and the exponent is the packing lattice d_max = 4.
#
# For higher Hopf:
#   - Base is S^4 with degeneracy g_l = (l+1)(l+2)(2l+3)/6  (degree 3 in l)
#   - Base is S^8 with degeneracy of degree 7 in l (see sympy formula)
# The Dirichlet series of a polynomial-degeneracy sequence in l is a
# sum of shifted Riemann zetas.  Let's compute it symbolically.

def fiber_F_dirichlet(N_base_dim, s_values):
    """
    Compute F(s) = sum_{l=1}^infty deg_SN(N_base_dim)(l) * l^{-s}
    for the base S^{N_base_dim}.

    Return a dict {s: symbolic closed form (zeta combination)}.
    """
    d = sp.expand(sp.expand_func(deg_SN(N_base_dim)))  # polynomial in l
    d = sp.expand(d)
    # write d = sum c_k * l^k
    poly = sp.Poly(d, l)
    coeffs = poly.all_coeffs()[::-1]  # c_0, c_1, ...
    results = {}
    for s in s_values:
        total = sp.Integer(0)
        for k, c in enumerate(coeffs):
            if c == 0:
                continue
            # sum_{l=1}^infty l^k * l^{-s} = zeta(s - k)
            total += c * sp.zeta(s - k)
        results[s] = sp.nsimplify(total)
    return results


# ---------------------------------------------------------------------------
# Delta analog (finite-size boundary)
# ---------------------------------------------------------------------------
# Paper 2: Delta = 1 / (|lambda_{n_max}| * N(n_max - 1)),
# where lambda_{n_max} = -(n_max^2 - 1) = -8 at n_max=3, and N(2) = 5.
#
# For the higher bases we mimic this with:
#   |lambda_{L_max}| = L_max * (L_max + N_base_dim - 1)    (top eigenvalue)
#   N(L_max - 1)     = total state count through level L_max - 1 on base
# We will compute this at whatever L_max the selection principle picks
# (if any), or at a default choice.

def Delta_analog(N_base_dim, L_max):
    """|lambda_top| * N(L_max - 1) reciprocal, analog of Paper 2 Delta."""
    if L_max <= 0:
        return sp.oo
    lam_top = sp.nsimplify(L_max * (L_max + N_base_dim - 1))
    N_below = N_base(N_base_dim, L_max - 1) if L_max >= 1 else sp.Integer(0)
    denom = lam_top * N_below
    if denom == 0:
        return sp.oo
    return sp.Rational(1, int(denom))


# ---------------------------------------------------------------------------
# Full K, alpha computation
# ---------------------------------------------------------------------------
def compute_K_alpha(B, F, Delta):
    K = sp.pi * (B + F - Delta)
    K_num = mp.mpf(str(sp.nsimplify(B + F - Delta).evalf(50))) * mp.pi
    alpha = 1 / K_num
    return K, K_num, alpha


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    results = {}

    # Sanity: verify degeneracy formulas
    results['degeneracy_checks'] = {k: bool(v) for k, v in verify_degeneracies().items()}

    # Sanity: reproduce Paper 2 B = 42
    B_paper2 = sp.Integer(0)
    for n in range(1, 4):
        for ll in range(n):
            B_paper2 += (2*ll + 1) * ll * (ll + 1)
    results['paper2_B'] = int(B_paper2)
    assert B_paper2 == 42

    # Paper 2 F and Delta
    F_paper2 = sp.pi**2 / 6
    Delta_paper2 = sp.Rational(1, 40)
    K_paper2 = sp.pi * (42 + F_paper2 - Delta_paper2)
    K_paper2_num = mp.pi * (mp.mpf(42) + mp.pi**2/6 - mp.mpf('0.025'))
    results['paper2_K'] = str(K_paper2_num)
    results['paper2_alpha_inv'] = str(K_paper2_num)

    # -----------------------------------------------------------------------
    # Direct base-truncation analog for S^4 and S^8
    # -----------------------------------------------------------------------
    # For Paper 2 itself, the "direct base truncation" equivalent on S^2
    # would cap l at some L_max and sum (2l+1)*l(l+1).  Paper 2 effectively
    # uses L_max = 2 (since Fock shell n_max=3 goes up to l=2), giving
    # B = 0 + 6 + 30 = 36 on S^2 alone; the extra 6 from shell n=2 (l=1)
    # is double-counted in Paper 2's shell structure.  Let's verify:
    B_S2_Lmax2 = 0
    for ll in range(3):
        B_S2_Lmax2 += (2*ll + 1) * ll * (ll + 1)
    results['B_S2_direct_Lmax2'] = B_S2_Lmax2  # should be 36
    # Paper 2's B=42 arises from the double count in the Fock shell scheme:
    # shell n=2 contains l=0,1 and shell n=3 contains l=0,1,2, so l=0 is
    # counted twice (contributes 0) and l=1 is counted twice (contributes
    # 2*6=12) and l=2 once (30), plus n=1 shell with l=0 (0). Total:
    #    n=1: 0
    #    n=2: 0 + 6 = 6
    #    n=3: 0 + 6 + 30 = 36
    #    sum = 42.  OK, verified.

    # -----------------------------------------------------------------------
    # Quaternionic Hopf: base S^4
    # -----------------------------------------------------------------------
    quat = {}
    table_4, hits_4 = find_selection_cutoff(
        4, targets=[3, 4, 7, sp.Rational(8, 3)], L_max_search=30
    )
    quat['ratio_table'] = [
        {'L_max': L, 'B': B, 'N': Nv, 'ratio_str': str(r), 'ratio_float': rf}
        for (L, B, Nv, r, rf) in table_4
    ]
    quat['selection_hits'] = {str(k): v for k, v in hits_4.items()}

    # Closed-form ratio B(L)/N(L) as a function of L (symbolic)
    L_sym = sp.symbols('L', integer=True, positive=True)
    # For S^4, sum symbolically
    d4 = sp.expand(sp.expand_func(deg_SN(4)))       # (l+1)(l+2)(2l+3)/6
    e4 = eig_SN(4)                  # l*(l+3)
    B4_expr = sp.simplify(sp.summation(d4 * e4, (l, 0, L_sym)))
    N4_expr = sp.simplify(sp.summation(d4,     (l, 0, L_sym)))
    ratio4_expr = sp.simplify(B4_expr / N4_expr)
    quat['B_closed_form'] = str(sp.factor(B4_expr))
    quat['N_closed_form'] = str(sp.factor(N4_expr))
    quat['ratio_closed_form'] = str(sp.factor(ratio4_expr))

    # Dirichlet-series fiber analog (natural exponent = dim of total space = 7,
    # and/or dim(S^4)+1 = 5 as 'packing' analog, plus a few near candidates)
    quat['F_dirichlet'] = {}
    for s in [5, 6, 7, 8]:
        Fs = fiber_F_dirichlet(4, [s])[s]
        Fs_num = mp.mpf(str(sp.nsimplify(Fs).evalf(50)))
        quat['F_dirichlet'][str(s)] = {'sym': str(Fs), 'num': str(Fs_num)}

    # Evaluate K at each reasonable (L_max, s) combination
    quat['K_alpha_grid'] = []
    for L_max_choice in range(1, 8):
        B_L = int(B_base(4, L_max_choice))
        N_L = int(N_base(4, L_max_choice))
        Delta_L = Delta_analog(4, L_max_choice)
        Delta_L_num = mp.mpf(str(sp.nsimplify(Delta_L).evalf(50)))
        for s in [5, 6, 7, 8]:
            Fs = fiber_F_dirichlet(4, [s])[s]
            Fs_num = mp.mpf(str(sp.nsimplify(Fs).evalf(50)))
            K_num = mp.pi * (mp.mpf(B_L) + Fs_num - Delta_L_num)
            alpha_inv = K_num / mp.pi  # for comparison with alpha-style forms
            quat['K_alpha_grid'].append({
                'L_max': L_max_choice,
                'B': B_L,
                'N': N_L,
                'BoverN': B_L / N_L if N_L else None,
                's_fiber': s,
                'F_num': str(Fs_num),
                'Delta': str(Delta_L_num),
                'K': str(K_num),
                'alpha=1/K': str(1/K_num),
                'K/pi_(for_alpha_em_style)': str(K_num/mp.pi),
            })

    # -----------------------------------------------------------------------
    # Octonionic Hopf: base S^8
    # -----------------------------------------------------------------------
    oct_ = {}
    table_8, hits_8 = find_selection_cutoff(
        8, targets=[7, 8, 15, sp.Rational(8, 3)], L_max_search=20
    )
    oct_['ratio_table'] = [
        {'L_max': L, 'B': B, 'N': Nv, 'ratio_str': str(r), 'ratio_float': rf}
        for (L, B, Nv, r, rf) in table_8
    ]
    oct_['selection_hits'] = {str(k): v for k, v in hits_8.items()}

    d8 = sp.expand(sp.expand_func(deg_SN(8)))
    e8 = eig_SN(8)
    B8_expr = sp.simplify(sp.summation(d8 * e8, (l, 0, L_sym)))
    N8_expr = sp.simplify(sp.summation(d8,     (l, 0, L_sym)))
    ratio8_expr = sp.simplify(B8_expr / N8_expr)
    oct_['B_closed_form'] = str(sp.factor(B8_expr))
    oct_['N_closed_form'] = str(sp.factor(N8_expr))
    oct_['ratio_closed_form'] = str(sp.factor(ratio8_expr))

    oct_['F_dirichlet'] = {}
    for s in [9, 10, 11, 12]:
        Fs = fiber_F_dirichlet(8, [s])[s]
        Fs_num = mp.mpf(str(sp.nsimplify(Fs).evalf(50)))
        oct_['F_dirichlet'][str(s)] = {'sym': str(Fs), 'num': str(Fs_num)}

    oct_['K_alpha_grid'] = []
    for L_max_choice in range(1, 6):
        B_L = int(B_base(8, L_max_choice))
        N_L = int(N_base(8, L_max_choice))
        Delta_L = Delta_analog(8, L_max_choice)
        Delta_L_num = mp.mpf(str(sp.nsimplify(Delta_L).evalf(50)))
        for s in [9, 10, 11, 12]:
            Fs = fiber_F_dirichlet(8, [s])[s]
            Fs_num = mp.mpf(str(sp.nsimplify(Fs).evalf(50)))
            K_num = mp.pi * (mp.mpf(B_L) + Fs_num - Delta_L_num)
            oct_['K_alpha_grid'].append({
                'L_max': L_max_choice,
                'B': B_L,
                'N': N_L,
                'BoverN': B_L / N_L if N_L else None,
                's_fiber': s,
                'F_num': str(Fs_num),
                'Delta': str(Delta_L_num),
                'K': str(K_num),
                'alpha=1/K': str(1/K_num),
                'K/pi': str(K_num/mp.pi),
            })

    # -----------------------------------------------------------------------
    # SM comparison targets
    # -----------------------------------------------------------------------
    sm_targets = {
        'alpha_em_inv': mp.mpf('137.035999084'),
        'alpha_em': 1/mp.mpf('137.035999084'),
        'alpha_W_inv_at_MZ': mp.mpf('29.5'),         # weak fine structure
        'alpha_W_at_MZ': 1/mp.mpf('29.5'),
        'alpha_s_at_MZ': mp.mpf('0.1181'),
        'alpha_s_inv_at_MZ': 1/mp.mpf('0.1181'),
        'sin2thetaW': mp.mpf('0.23122'),
        '8/3_SU5_ratio': mp.mpf(8)/3,
    }
    results['sm_targets'] = {k: str(v) for k, v in sm_targets.items()}

    # Scan the grid for near-hits
    def scan_for_hits(grid, label):
        hits = []
        for entry in grid:
            alpha_from_K = mp.mpf(entry['alpha=1/K'])
            K_val = mp.mpf(entry['K'])
            K_over_pi = mp.mpf(entry['K/pi']) if 'K/pi' in entry else mp.mpf(entry['K/pi_(for_alpha_em_style)'])
            for name, tgt in sm_targets.items():
                tgt_val = tgt
                rel_alpha = abs(alpha_from_K - tgt_val) / tgt_val
                rel_K = abs(K_val - tgt_val) / tgt_val
                rel_Kpi = abs(K_over_pi - tgt_val) / tgt_val
                for kind, rel, val in [
                    ('alpha=1/K vs ' + name, rel_alpha, alpha_from_K),
                    ('K vs ' + name, rel_K, K_val),
                    ('K/pi vs ' + name, rel_Kpi, K_over_pi),
                ]:
                    if float(rel) < 0.10:
                        hits.append({
                            'label': label,
                            'L_max': entry['L_max'],
                            's_fiber': entry['s_fiber'],
                            'kind': kind,
                            'value': str(val),
                            'target': str(tgt_val),
                            'rel_err': str(rel),
                        })
        return hits

    near_hits = scan_for_hits(quat['K_alpha_grid'], 'quat_S4') \
              + scan_for_hits(oct_['K_alpha_grid'], 'oct_S8')
    results['near_hits_within_10pct'] = near_hits

    # -----------------------------------------------------------------------
    # Selection-principle-based "canonical" results
    # -----------------------------------------------------------------------
    # Paper 2 picks L_max via B/N = dim(S^3) = 3.  Let's find the dim(base)
    # analog:  B/N = dim(base).
    # For S^4, dim = 4.  Closed-form ratio B4/N4:
    #    Already computed as ratio4_expr
    # Solve ratio4 = 4 over L_sym:
    try:
        sols4_dim_base = sp.solve(sp.Eq(ratio4_expr, 4), L_sym)
        quat['selection_B/N=4_solutions'] = [str(s) for s in sols4_dim_base]
    except Exception as e:
        quat['selection_B/N=4_solutions'] = 'solve_failed: ' + str(e)

    try:
        sols4_dim_total = sp.solve(sp.Eq(ratio4_expr, 7), L_sym)  # total space S^7
        quat['selection_B/N=7_solutions'] = [str(s) for s in sols4_dim_total]
    except Exception as e:
        quat['selection_B/N=7_solutions'] = 'solve_failed: ' + str(e)

    try:
        sols8_dim_base = sp.solve(sp.Eq(ratio8_expr, 8), L_sym)
        oct_['selection_B/N=8_solutions'] = [str(s) for s in sols8_dim_base]
    except Exception as e:
        oct_['selection_B/N=8_solutions'] = 'solve_failed: ' + str(e)

    try:
        sols8_dim_total = sp.solve(sp.Eq(ratio8_expr, 15), L_sym)
        oct_['selection_B/N=15_solutions'] = [str(s) for s in sols8_dim_total]
    except Exception as e:
        oct_['selection_B/N=15_solutions'] = 'solve_failed: ' + str(e)

    # -----------------------------------------------------------------------
    # Try the Cl(6) dimension hit
    # -----------------------------------------------------------------------
    # Furey construction uses Cl(6), dim 2^6 = 64. Does B(L_max) = 64 appear?
    for L_max_choice in range(1, 12):
        B_4 = int(B_base(4, L_max_choice))
        B_8 = int(B_base(8, L_max_choice))
        if B_4 == 64:
            quat['B=64_L_max'] = L_max_choice
        if B_8 == 64:
            oct_['B=64_L_max'] = L_max_choice

    results['quaternionic_S4_base'] = quat
    results['octonionic_S8_base'] = oct_

    # Save
    outdir = Path(__file__).resolve().parent / 'data' / 'track_alpha_sm'
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / 'sm_f_higher_hopf.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)

    # Print summary
    print('='*78)
    print('Track SM-F: Higher Hopf fibrations - summary')
    print('='*78)
    print()
    print('Paper 2 sanity: B = 42, K = {:.6f}'.format(float(K_paper2_num)))
    print()
    print('--- Quaternionic S^4 base -----------------------------------------')
    print('B(L)/N(L) closed form:', quat['ratio_closed_form'])
    for entry in table_4[:10]:
        L, B, N_, r, rf = entry
        print(f'  L={L:2d}  B={B:<10d}  N={N_:<8d}  B/N = {r}  ({rf:.4f})')
    print('  Selection B/N=4 (dim base):  ', quat.get('selection_B/N=4_solutions'))
    print('  Selection B/N=7 (dim total): ', quat.get('selection_B/N=7_solutions'))
    print()
    print('--- Octonionic S^8 base -------------------------------------------')
    print('B(L)/N(L) closed form:', oct_['ratio_closed_form'])
    for entry in table_8[:10]:
        L, B, N_, r, rf = entry
        print(f'  L={L:2d}  B={B:<12d}  N={N_:<10d}  B/N = {r}  ({rf:.4f})')
    print('  Selection B/N=8 (dim base):  ', oct_.get('selection_B/N=8_solutions'))
    print('  Selection B/N=15 (dim total):', oct_.get('selection_B/N=15_solutions'))
    print()
    print('--- SM near-hits (rel err < 10%) ----------------------------------')
    if near_hits:
        for h in near_hits[:40]:
            print(f'  {h["label"]}  L={h["L_max"]}  s={h["s_fiber"]}  '
                  f'{h["kind"]} = {h["value"][:12]}  vs {h["target"][:12]}  '
                  f'rel_err={float(mp.mpf(h["rel_err"])):.3e}')
    else:
        print('  (none)')

    return results


if __name__ == '__main__':
    main()
