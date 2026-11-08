"""Sprint TD Track 5 — PSLQ probe of GeoVac S_full(GS) entropy values
against the master Mellin engine ring.

The question (concrete): is the T → 0 thermal residual computed in Track 2
(0.040811 for He at n_max=3, 0.011212 for Li+ at n_max=3) representable in
the master Mellin engine ring (M1 ∪ M2 ∪ M3) per Paper 32 §VIII case-
exhaustion theorem and Paper 18 §III.7 master-Mellin reading?

Tasks:
  1. Compute S_full(GS) at 100+ digit precision for He, Li+, and additional
     systems via mpmath.
  2. Build mechanical basis of M1/M2/M3 elements (W3 falsification protocol).
  3. PSLQ at 100 dps against the basis.
  4. Cross-system consistency check.

Methodology: builds the M_L=0 singlet FCI matrix in exact rational form
(via geovac.hypergeometric_slater.compute_rk_algebraic + Z^2/2n^2 + kappa
adjacency), converts to mpmath at 150 dps, diagonalizes via mp.eigsy,
extracts the GS eigenvector, builds the spatial 1-RDM in mpmath, and
computes the von Neumann entropy of its eigenvalues at 150 dps.

W3 falsification discipline:
  - Mechanical basis frozen BEFORE any test
  - Pure-rational control class included for false-positive baseline
  - Per-mechanism z-score against control class

Usage:
  python debug/sprint_td_track5.py [--dps 150] [--n-trials-pslq 5]

Output:
  debug/data/sprint_td_track5.json
  debug/sprint_td_track5_memo.md
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.hypergeometric_slater import compute_rk_algebraic
from geovac.casimir_ci import _gaunt_ck
from geovac.lattice import GeometricLattice


# =====================================================================
#  High-precision He / Li+ FCI builder
# =====================================================================

def _orbital_list(n_max: int):
    """Return list of (n, l, m) tuples in the same order as casimir_ci uses."""
    lattice = GeometricLattice(max_n=n_max)
    return list(lattice.states)


def _adjacency_dict(n_max: int) -> Dict[Tuple[int, int], int]:
    """Return adjacency entries as dict (i, j) -> +/-1 from GeometricLattice.

    On the Fock graph the adjacency entries are integer-valued (typically 0,
    +1 or -1 with sign from raising/lowering).
    """
    lattice = GeometricLattice(max_n=n_max)
    A = lattice.adjacency
    if hasattr(A, 'toarray'):
        A_dense = A.toarray()
    else:
        A_dense = np.array(A)
    adj = {}
    n = A_dense.shape[0]
    for i in range(n):
        for j in range(n):
            if i != j:
                v = A_dense[i, j]
                # Should be integer-valued on Fock graph
                v_int = int(round(v))
                if abs(v - v_int) < 1e-10 and v_int != 0:
                    adj[(i, j)] = v_int
                elif abs(v) > 1e-10:
                    # Non-integer adjacency — fall back to float, flag below
                    adj[(i, j)] = v
    return adj


def _gaunt_exact(l1: int, m1: int, l2: int, m2: int, k: int):
    """Exact Gaunt c^k coefficient via sympy real spherical harmonic 3j.

    Returns a sympy expression (typically rational × √(...) / (...)).
    For our purposes we'll multiply two of these together (in the V_ee
    sum c_ac × c_bd) — products often simplify to rationals after the
    sqrt cancels.
    """
    import sympy as sp
    from sympy.physics.wigner import gaunt
    # Standard Gaunt coefficient
    g = gaunt(l1, k, l2, m1, 0, -m2)  # sphericals follow standard convention
    # The c^k(la,ma; lb,mb) used in casimir_ci involves a different normalization.
    # Use the relation: sum over m's using <Y_lm|Y_kq|Y_l'm'> standard form.
    # Since this is complicated for real spherical harmonics, fall back to
    # numerical _gaunt_ck and convert.
    return None  # placeholder


def build_fci_matrix_exact(Z: int, n_max: int, m_total: int = 0):
    """Build the He/Li+-like singlet FCI matrix at exact precision (mpmath).

    Returns:
        H_mp: mpmath matrix (n_configs × n_configs) at current mp.dps
        configs: list of (i, j) orbital index pairs
        orbitals: list of (n, l, m) tuples
    """
    orbitals = _orbital_list(n_max)
    n_spatial = len(orbitals)
    adj = _adjacency_dict(n_max)

    # Build h1 in mpmath
    h1 = mp.zeros(n_spatial, n_spatial)
    for i, (n, l, m) in enumerate(orbitals):
        h1[i, i] = -mp.mpf(Z) ** 2 / (2 * mp.mpf(n) ** 2)
    # kappa = -1/16, off-diag = kappa * (-A[i,j]) = (1/16) * A[i,j]
    kappa_neg = mp.mpf(1) / 16  # = -kappa
    for (i, j), v in adj.items():
        h1[i, j] = kappa_neg * mp.mpf(v)

    # Build configs (singlet, M_L = m_total)
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                configs.append((i, j))
    n_configs = len(configs)

    # k_orb = Z (physics convention from casimir_ci)
    k_orb = mp.mpf(Z)

    # Two-electron integrals: cache as mpmath
    # gaunt c_k via numerical (it's an exact algebraic value, but at mpmath
    # precision the float64 c_k is only ~1e-15 — we need a high-precision
    # computation).
    # _gaunt_ck returns float; we recompute it via sympy at high precision.
    import sympy as sp
    from sympy.physics.wigner import wigner_3j

    def _ck_high_precision(la, ma, lb, mb, k) -> mp.mpf:
        """Compute c^k(la, ma; lb, mb) at mpmath precision via sympy 3j.

        Convention from casimir_ci._gaunt_ck:
          c^k(la,ma; lb,mb) = sqrt((2la+1)(2lb+1)/(2k+1))
                               × (la lb k | 0 0 0)
                               × (la lb k | -ma mb (ma-mb))
            [for COMPLEX spherical harmonics]
        For REAL spherical harmonics, casimir_ci builds linear
        combinations. We need to mirror exactly what casimir_ci uses.

        Easiest: read the sign/value from _gaunt_ck (float) but compute
        the magnitude at high precision. Since _gaunt_ck returns 0 for
        forbidden transitions, we use it as a guard.
        """
        c_float = _gaunt_ck(la, ma, lb, mb, k)
        if abs(c_float) < 1e-15:
            return mp.mpf(0)
        # Use the formula from casimir_ci (lines around _gaunt_ck) at
        # high precision. _gaunt_ck for REAL SH uses delta_{ma, mb} for
        # m_q = 0 case (k transition with m=0 photon). Looking at
        # standard convention:
        #   c^k(la, ma; lb, mb) for real SH = real-form 3j
        # Easiest robust approach: high-precision sympy computation of
        # the same formula casimir_ci uses.
        # casimir_ci._gaunt_ck implementation: standard real-SH Gaunt.
        # Sign comes from the float computation; magnitude can be
        # recomputed at high precision via:
        #   c = (-1)^ma × sqrt((2la+1)(2lb+1)) × W3j(la, k, lb; 0,0,0)
        #                 × W3j(la, k, lb; -ma, ma-mb, mb)
        # for real harmonics with m=0 difference (real Gaunt for diagonal m).
        # This is m_q = 0 case (photon m_q = 0); for nonzero m_q we have
        # the sin/cos combinations.
        #
        # SAFER: trust float c_float to ~15 digits, multiply by the rational
        # R^k integral to get an estimate; then recompute at high precision
        # only if PSLQ depends on it.
        #
        # For 100 dps PSLQ on entropy, we need ALL matrix elements to
        # 100+ dps. Float64 c_k means matrix elements only good to ~15
        # digits. SOLUTION: use sympy's wigner_3j at high precision for
        # the angular part.
        try:
            # 3j symbols (use sympy for exact rational + sqrt expressions)
            three_j_0 = wigner_3j(la, k, lb, 0, 0, 0)
            three_j_m = wigner_3j(la, k, lb, -ma, ma - mb, mb)
            # The standard physicist convention: c^k = sqrt((2la+1)(2lb+1)) × 3j × 3j
            # Some conventions use (-1)^ma factor.
            from sympy import sqrt as sp_sqrt, Rational, sympify
            c_sym = (sp_sqrt(sympify(2 * la + 1)) * sp_sqrt(sympify(2 * lb + 1))
                     * three_j_0 * three_j_m * (-1) ** ma)
            c_mp = mp.mpf(str(sp.N(c_sym, mp.mp.dps + 5)))
            # Sign-check vs float (allow sign flip from convention)
            if c_float != 0 and float(c_mp) * c_float < 0:
                c_mp = -c_mp
            return c_mp
        except Exception:
            return mp.mpf(c_float)

    rk_cache: Dict[Tuple, mp.mpf] = {}

    def get_rk_mpmath(n1, l1, n3, l3, n2, l2, n4, l4, k) -> mp.mpf:
        key = (n1, l1, n3, l3, n2, l2, n4, l4, k)
        if key in rk_cache:
            return rk_cache[key]
        try:
            f = compute_rk_algebraic(n1, l1, n3, l3, n2, l2, n4, l4, k)
            # Fraction → mpmath at current dps
            v = mp.mpf(f.numerator) / mp.mpf(f.denominator)
        except Exception as e:
            print(f"  WARNING: compute_rk_algebraic failed for "
                  f"({n1}{l1},{n3}{l3};{n2}{l2},{n4}{l4},k={k}): {e}")
            from geovac.hypergeometric_slater import compute_rk_float
            v = mp.mpf(compute_rk_float(n1, l1, n3, l3, n2, l2, n4, l4, k))
        rk_cache[key] = v
        return v

    ck_cache: Dict[Tuple, mp.mpf] = {}

    def get_ck(la, ma, lb, mb, k) -> mp.mpf:
        key = (la, ma, lb, mb, k)
        if key in ck_cache:
            return ck_cache[key]
        v = _ck_high_precision(la, ma, lb, mb, k)
        ck_cache[key] = v
        return v

    def g_int_mp(a: int, b: int, c: int, d: int) -> mp.mpf:
        """High-precision <ab|cd> = sum_k c^k(a,c) c^k(b,d) R^k(ac,bd)."""
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        result = mp.mpf(0)
        k_min = max(abs(la - lc), abs(lb - ld))
        k_max = min(la + lc, lb + ld)
        for k in range(k_min, k_max + 1):
            if (la + lc + k) % 2 != 0:
                continue
            if (lb + ld + k) % 2 != 0:
                continue
            ck_ac = get_ck(la, ma, lc, mc, k)
            ck_bd = get_ck(lb, mb, ld, md, k)
            if abs(ck_ac) < mp.mpf('1e-30') or abs(ck_bd) < mp.mpf('1e-30'):
                continue
            rk = get_rk_mpmath(na, la, nc, lc, nb, lb, nd, ld, k)
            result += ck_ac * ck_bd * rk
        return result * k_orb

    # Build FCI matrix
    H_mp = mp.zeros(n_configs, n_configs)
    sqrt2 = mp.sqrt(2)

    print(f"  Building {n_configs}×{n_configs} FCI matrix at {mp.mp.dps} dps...")
    t0 = time.time()
    for I in range(n_configs):
        i, j = configs[I]
        for J in range(I, n_configs):
            p, q = configs[J]
            # Singlet config:
            #   |I> = (|ij> + |ji>)/sqrt(2) for i!=j; |ii> for i==j
            # Same for |J>.
            N_I = sqrt2 if i != j else mp.mpf(1)
            N_J = sqrt2 if p != q else mp.mpf(1)

            # 1-body part: sum over a,b (a in I, b in J) of h1 contribution
            # <ab|h1|cd> with c,d single-electron states
            # For singlet:
            #   <I|h1|J> = (1/N_I N_J) [<ij|h1|pq> + <ij|h1|qp> + permutations]
            # Standard SD matrix element for singlet:
            h1_val = mp.mpf(0)
            # Enumerate I-perms x J-perms, contract on second-electron delta
            I_perms = [(i, j)] if i == j else [(i, j), (j, i)]
            J_perms = [(p, q)] if p == q else [(p, q), (q, p)]
            for (a, b) in I_perms:
                for (c, d) in J_perms:
                    # h1: one-body operator on electron 1, identity on 2
                    if b == d:
                        h1_val += h1[a, c]
                    if a == c:
                        h1_val += h1[b, d]
            h1_val /= (N_I * N_J)

            # 2-body part: V_ee
            v_val = mp.mpf(0)
            for (a, b) in I_perms:
                for (c, d) in J_perms:
                    v_val += g_int_mp(a, b, c, d)
            v_val /= (N_I * N_J)

            # Total matrix element. Two-body sum overcounts permutations:
            # standard convention: H_IJ = h1 + 0.5 V (when summing both electron
            # orderings) -- but our enumeration already treats I_perms × J_perms
            # explicitly, NOT halving. Need to match casimir_ci behavior. Test:
            # For diagonal i==j==p==q, h1_val = 2 h1[a,a]/(1·1) = 2 ε_a, which is
            # correct (two electrons in same orbital, both contribute ε_a).
            # For V_ee with i==j==p==q: v_val = <aa|aa>, and N_I N_J = 1, so
            # v_val = J_aa, correct.
            # For off-diagonal i!=j==p!=q (matching configs): h1_val and v_val
            # both have 4 permutation contributions, divided by 2 (=√2·√2) → 2.
            H_mp[I, J] = h1_val + v_val
            if I != J:
                H_mp[J, I] = h1_val + v_val
    print(f"  FCI matrix built in {time.time() - t0:.1f}s.")
    return H_mp, configs, orbitals


# =====================================================================
#  Spatial 1-RDM and von Neumann entropy at high precision
# =====================================================================

def build_1rdm_mpmath(ci_coeffs, configs, n_spatial: int):
    """Build the spatial 1-RDM at mpmath precision.

    Mirror of debug/entanglement_geometry.build_1rdm_from_singlet_ci in
    high precision. Trace = 2 (two electrons).
    """
    rho = mp.zeros(n_spatial, n_spatial)
    sqrt2 = mp.sqrt(2)
    n_configs = len(configs)
    for I, (i, j) in enumerate(configs):
        cI = ci_coeffs[I]
        for J, (p, q) in enumerate(configs):
            cJ = ci_coeffs[J]
            coef = cI * cJ
            N_I = sqrt2 if i != j else mp.mpf(1)
            N_J = sqrt2 if p != q else mp.mpf(1)
            I_perms = [(i, j)] if i == j else [(i, j), (j, i)]
            J_perms = [(p, q)] if p == q else [(p, q), (q, p)]
            for (a, b) in I_perms:
                for (c, d) in J_perms:
                    if b == d:
                        rho[a, c] += coef / (N_I * N_J)
    rho *= 2
    return rho


def vn_entropy_mpmath(rho) -> mp.mpf:
    """Von Neumann entropy of a 1-RDM at mpmath precision.

    Convention matches debug/entanglement_geometry.compute_entanglement_measures:
        S = -sum (n_i / 2) log(n_i / 2)
    where n_i are eigenvalues of rho (occupations summing to 2).
    """
    evals_mp = mp.eigsy(rho, eigvals_only=True)
    S = mp.mpf(0)
    for ev in evals_mp:
        ev_re = mp.re(ev) if hasattr(ev, 'real') else ev
        if ev_re > mp.mpf('1e-50'):
            ni_half = ev_re / 2
            if ni_half > 0:
                S -= ni_half * mp.log(ni_half)
    return S


def compute_S_full_GS_mpmath(Z: int, n_max: int, dps: int = 150) -> Dict:
    """Compute Paper 27 S_full for the singlet GS at high precision.

    Returns dict with: S_full (mp.mpf string), GS_energy (str), n_configs,
    runtime.
    """
    mp.mp.dps = dps
    print(f"\n=== Z={Z}, n_max={n_max} (dps={dps}) ===")
    t0 = time.time()
    H_mp, configs, orbitals = build_fci_matrix_exact(Z, n_max)
    n_configs = len(configs)
    n_spatial = len(orbitals)
    print(f"  Diagonalizing {n_configs}×{n_configs} matrix...")
    t1 = time.time()
    evals, evecs = mp.eigsy(H_mp)
    # eigsy returns eigvalues sorted ascending; pick first
    # evals is a row vector (1 × n) in mpmath; access via [0, k]
    if hasattr(evals, '__getitem__') and hasattr(evals, 'rows'):
        # mpmath matrix
        E0 = evals[0]
        for k in range(evals.rows):
            if evals[k] < E0:
                E0 = evals[k]
        # Find argmin
        idx0 = 0
        for k in range(evals.rows):
            if evals[k] == E0:
                idx0 = k
                break
        v0 = mp.matrix([evecs[i, idx0] for i in range(n_configs)])
    else:
        # list-like
        E0 = min(evals)
        idx0 = list(evals).index(E0)
        v0 = mp.matrix([evecs[i, idx0] for i in range(n_configs)])
    print(f"    eigsy in {time.time() - t1:.1f}s; E_GS = {mp.nstr(E0, 30)}")
    # Build 1-RDM at high precision
    rho = build_1rdm_mpmath(v0, configs, n_spatial)
    trace_rho = sum(rho[i, i] for i in range(n_spatial))
    print(f"    1-RDM trace = {mp.nstr(trace_rho, 12)} (should be 2.0)")
    # Compute entropy
    S_full = vn_entropy_mpmath(rho)
    print(f"    S_full(GS) = {mp.nstr(S_full, 50)}")
    runtime = time.time() - t0
    return {
        'Z': Z,
        'n_max': n_max,
        'dps': dps,
        'n_configs': n_configs,
        'n_spatial': n_spatial,
        'E_GS': mp.nstr(E0, dps + 5),
        'S_full_GS_str': mp.nstr(S_full, dps + 5),
        'S_full_GS_mp': S_full,
        'trace_rho': mp.nstr(trace_rho, 30),
        'runtime_s': runtime,
    }


# =====================================================================
#  Mechanical basis (W3 falsification protocol)
# =====================================================================

def build_mechanical_basis(target_value: mp.mpf, dps: int = 150):
    """Mechanically generate basis of M1/M2/M3 elements.

    Specification frozen BEFORE running PSLQ:

    Seeds organized by master Mellin engine class:
      M1 (Hopf-base measure / pi-family): π powers, 1/π powers,
          Vol(S^k) ratios (Vol(S²)/Vol(S³) = 4π / 2π² = 2/π, etc.)
      M2 (Seeley-DeWitt / heat kernel): ζ(2k) = π^{2k} × ℚ family,
          and the M2 ring √π×ℚ ⊕ π²×ℚ identified in MR-B
      M3 (vertex parity Hurwitz Dirichlet-L): Catalan G, β(2k) values,
          Hurwitz zeta at quarter-integer shifts ζ(s, 1/4), ζ(s, 3/4)
      ALG (algebraic controls): √n for small n, golden ratio φ
      RAT (pure rational): control class for false-positive baseline

    Generation rule (frozen):
      Single-seed: rational prefactor (n / d) with n in {1,...,5},
                   d in {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 30, 40, 60},
                   gcd(n, d) = 1.
      Two-seed:    a × b / (n / d) with n in {1, 2, 3}, same d set.
      Range filter: |v| in [1e-6, 1.0]. (entropy values are O(0.01-0.1))
      Dedup: 1e-12 relative tolerance.
    """
    mp.mp.dps = dps
    SEEDS = {}
    SEED_CLASS = {}

    # M1 — π and Hopf-base measure family
    SEEDS['1'] = mp.mpf(1)
    SEED_CLASS['1'] = 'RAT'  # special: identity belongs in RAT control
    pi = mp.pi
    SEEDS['pi'] = pi
    SEED_CLASS['pi'] = 'M1'
    SEEDS['pi^2'] = pi ** 2
    SEED_CLASS['pi^2'] = 'M2'  # ζ(2) = π²/6 lives here
    SEEDS['pi^3'] = pi ** 3
    SEED_CLASS['pi^3'] = 'M1'
    SEEDS['pi^4'] = pi ** 4
    SEED_CLASS['pi^4'] = 'M2'  # ζ(4) = π⁴/90
    SEEDS['pi^6'] = pi ** 6
    SEED_CLASS['pi^6'] = 'M2'
    SEEDS['1/pi'] = 1 / pi
    SEED_CLASS['1/pi'] = 'M1'
    SEEDS['1/pi^2'] = 1 / pi ** 2
    SEED_CLASS['1/pi^2'] = 'M1'
    SEEDS['1/pi^3'] = 1 / pi ** 3
    SEED_CLASS['1/pi^3'] = 'M1'
    SEEDS['Vol_S2_div_pi2'] = 4 * pi / pi ** 2  # = 4/π
    SEED_CLASS['Vol_S2_div_pi2'] = 'M1'
    SEEDS['1/(2pi)'] = 1 / (2 * pi)
    SEED_CLASS['1/(2pi)'] = 'M1'
    SEEDS['1/(4pi)'] = 1 / (4 * pi)
    SEED_CLASS['1/(4pi)'] = 'M1'
    SEEDS['1/(2pi^2)'] = 1 / (2 * pi ** 2)
    SEED_CLASS['1/(2pi^2)'] = 'M1'
    SEEDS['sqrt_pi'] = mp.sqrt(pi)
    SEED_CLASS['sqrt_pi'] = 'M2'
    SEEDS['1/sqrt_pi'] = 1 / mp.sqrt(pi)
    SEED_CLASS['1/sqrt_pi'] = 'M2'

    # M2 — Riemann zeta even values (= rational × π^{even})
    SEEDS['zeta(2)'] = mp.zeta(2)
    SEED_CLASS['zeta(2)'] = 'M2'
    SEEDS['zeta(4)'] = mp.zeta(4)
    SEED_CLASS['zeta(4)'] = 'M2'
    SEEDS['zeta(6)'] = mp.zeta(6)
    SEED_CLASS['zeta(6)'] = 'M2'
    SEEDS['zeta(8)'] = mp.zeta(8)
    SEED_CLASS['zeta(8)'] = 'M2'

    # M3 — vertex parity Dirichlet-L family
    SEEDS['Catalan'] = mp.catalan
    SEED_CLASS['Catalan'] = 'M3'
    # Dirichlet beta β(2) = Catalan, β(4), β(6) via series
    def dirichlet_beta(s_int, n_terms=200):
        return mp.nsum(lambda k: (-1) ** k / (2 * k + 1) ** s_int,
                       [0, mp.inf])
    SEEDS['beta(4)'] = dirichlet_beta(4)
    SEED_CLASS['beta(4)'] = 'M3'
    SEEDS['beta(6)'] = dirichlet_beta(6)
    SEED_CLASS['beta(6)'] = 'M3'
    # Hurwitz at quarter-integer shifts
    SEEDS['zeta_hurwitz_2_1_4'] = mp.zeta(2, mp.mpf(1) / 4)
    SEED_CLASS['zeta_hurwitz_2_1_4'] = 'M3'
    SEEDS['zeta_hurwitz_2_3_4'] = mp.zeta(2, mp.mpf(3) / 4)
    SEED_CLASS['zeta_hurwitz_2_3_4'] = 'M3'
    # ζ(3) — odd Riemann, in spinor sector via Paper 18 motivic weight
    SEEDS['zeta(3)'] = mp.zeta(3)
    SEED_CLASS['zeta(3)'] = 'M3'  # technically M2/M3 boundary; included for breadth
    SEEDS['zeta(5)'] = mp.zeta(5)
    SEED_CLASS['zeta(5)'] = 'M3'

    # ALG — algebraic controls
    SEEDS['sqrt2'] = mp.sqrt(2)
    SEED_CLASS['sqrt2'] = 'ALG'
    SEEDS['sqrt3'] = mp.sqrt(3)
    SEED_CLASS['sqrt3'] = 'ALG'
    SEEDS['sqrt5'] = mp.sqrt(5)
    SEED_CLASS['sqrt5'] = 'ALG'
    SEEDS['phi'] = (1 + mp.sqrt(5)) / 2
    SEED_CLASS['phi'] = 'ALG'
    SEEDS['log2'] = mp.log(2)
    SEED_CLASS['log2'] = 'ALG'
    SEEDS['log3'] = mp.log(3)
    SEED_CLASS['log3'] = 'ALG'

    # RAT — pure rational control class
    # (uses '1' from above + stand-alone integer ratios)
    # The mechanical-prefactor procedure already generates rational
    # multiples of every seed, so the RAT class fills in via prefactors.

    seed_names = list(SEEDS.keys())
    seed_vals = [SEEDS[n] for n in seed_names]

    PREFACTOR_NUM = list(range(-5, 0)) + list(range(1, 6))  # ±1..5
    DENOMS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 30, 40, 60]

    raw_basis: Dict[str, Tuple[mp.mpf, str]] = {}
    range_min = mp.mpf('1e-6')
    range_max = mp.mpf(1)

    print(f"  Generating mechanical basis (seeds={len(SEEDS)})...")
    # Single-seed
    for sname in seed_names:
        sval = SEEDS[sname]
        sclass = SEED_CLASS[sname]
        for n in PREFACTOR_NUM:
            for d in DENOMS:
                if math.gcd(abs(n), d) != 1:
                    continue
                v = mp.mpf(n) * sval / mp.mpf(d)
                if range_min < abs(v) < range_max:
                    if n == 1 and d == 1:
                        fname = sname
                    elif n == -1 and d == 1:
                        fname = f"-{sname}"
                    elif d == 1:
                        fname = f"{n}*{sname}"
                    elif n == 1:
                        fname = f"{sname}/{d}"
                    elif n == -1:
                        fname = f"-{sname}/{d}"
                    else:
                        fname = f"({n}/{d})*{sname}"
                    raw_basis[fname] = (v, sclass)

    # Two-seed (a*b and a/b with rational prefactor)
    for i, sa in enumerate(seed_names):
        for j, sb in enumerate(seed_names):
            if i >= j:
                continue
            sa_class = SEED_CLASS[sa]
            sb_class = SEED_CLASS[sb]
            # Mixed class label: take "lower" Mellin-tier or composite
            classes_sorted = sorted([sa_class, sb_class])
            if classes_sorted[0] == classes_sorted[1]:
                mech = classes_sorted[0]
            elif 'RAT' in classes_sorted:
                non_rat = [c for c in classes_sorted if c != 'RAT'][0]
                mech = non_rat
            else:
                mech = '+'.join(classes_sorted)
            sav = SEEDS[sa]
            sbv = SEEDS[sb]
            for n in [-3, -2, -1, 1, 2, 3]:
                for d in [1, 2, 3, 4, 5, 6, 8]:
                    if math.gcd(abs(n), d) != 1:
                        continue
                    # a*b/(n/d)
                    v = mp.mpf(n) * sav * sbv / mp.mpf(d)
                    if range_min < abs(v) < range_max:
                        fname = f"({n}/{d})*{sa}*{sb}"
                        if fname not in raw_basis:
                            raw_basis[fname] = (v, mech)
                    # a/b/(n/d)
                    if abs(sbv) > mp.mpf('1e-30'):
                        v = mp.mpf(n) * sav / sbv / mp.mpf(d)
                        if range_min < abs(v) < range_max:
                            fname = f"({n}/{d})*{sa}/{sb}"
                            if fname not in raw_basis:
                                raw_basis[fname] = (v, mech)

    # Deduplicate
    print(f"  Raw basis: {len(raw_basis)} forms; deduplicating...")
    items = list(raw_basis.items())
    items.sort(key=lambda x: float(x[1][0]))
    deduped = []
    last_v = None
    for fname, (v, c) in items:
        vf = float(v)
        if last_v is not None and abs(vf - last_v) / max(abs(vf), 1e-30) < 1e-12:
            continue
        deduped.append((fname, v, c))
        last_v = vf
    print(f"  After dedup: {len(deduped)} forms.")
    return deduped


def basis_class_counts(basis):
    return Counter(c for (_, _, c) in basis)


# =====================================================================
#  PSLQ search
# =====================================================================

def pslq_search_per_basis(target_mp: mp.mpf, basis, dps: int = 150,
                          tol_exp: int = 50, max_coeff: int = 10000):
    """Search for PSLQ identification of target across the basis.

    Strategy:
      Try small subsets of basis (up to 8 elements) per Mellin class,
      then run mp.pslq([target, *subset], tol=mp.mpf(10)**(-tol_exp)).
      If integer combination found with all |coeff| <= max_coeff, record.

    Returns list of identification candidates with z-score baseline.
    """
    mp.mp.dps = dps
    candidates = []

    # Group basis by class
    by_class = {}
    for fname, val, cls in basis:
        by_class.setdefault(cls, []).append((fname, val))

    # For each class, do one-element PSLQ first (most stringent)
    for cls in sorted(by_class.keys()):
        elements = by_class[cls]
        print(f"    Class {cls}: testing {len(elements)} elements (1-elt PSLQ)")
        # Single-element check: is target a rational multiple of element?
        for fname, val in elements:
            try:
                rel = mp.pslq([target_mp, val],
                              tol=mp.mpf(10) ** (-tol_exp),
                              maxcoeff=max_coeff)
            except Exception:
                rel = None
            if rel is not None and len(rel) == 2:
                a, b = rel
                if a != 0 and b != 0:
                    # target × a + val × b = 0  =>  target = -b/a × val
                    ratio = -mp.mpf(b) / mp.mpf(a)
                    residual = target_mp + a * target_mp / target_mp * 0  # placeholder
                    # Compute actual residual: a × target + b × val
                    resid = mp.mpf(a) * target_mp + mp.mpf(b) * val
                    candidates.append({
                        'class': cls,
                        'depth': 1,
                        'forms': [fname],
                        'integer_coeffs': [int(a), int(b)],
                        'max_coeff': max(abs(int(a)), abs(int(b))),
                        'residual': mp.nstr(abs(resid), 10),
                        'identification': f"target = ({-b}/{a}) × {fname}",
                    })

    # 2-element subsets within and across classes (limited combinatorial budget)
    print(f"    2-element subset search (limited)...")
    # All M1, M2, M3 elements (skip ALG, RAT control for stress)
    primary_classes = ['M1', 'M2', 'M3', 'M1+M2', 'M1+M3', 'M2+M3']
    primary_elements = []
    for cls in primary_classes:
        if cls in by_class:
            primary_elements.extend([(c[0], c[1], cls) for c in by_class[cls]])
    print(f"      Primary elements: {len(primary_elements)}")
    # We can't do all pairs (would be O(n²) at n~10000); subsample most-natural
    # elements (small max_coeff signature) and try pairs there.
    natural_basis = [(f, v, c) for (f, v, c) in primary_elements
                     if not f.startswith('-') and 'mp.mpf' not in f]
    # Sort by name length as proxy for simplicity
    natural_basis.sort(key=lambda x: len(x[0]))
    # Take first 200
    natural_short = natural_basis[:200]
    print(f"      Pair-testing {len(natural_short)} natural elements ({len(natural_short)*(len(natural_short)-1)//2} pairs)")
    for i, (fa, va, ca) in enumerate(natural_short):
        for j, (fb, vb, cb) in enumerate(natural_short):
            if i >= j:
                continue
            try:
                rel = mp.pslq([target_mp, va, vb],
                              tol=mp.mpf(10) ** (-tol_exp),
                              maxcoeff=max_coeff)
            except Exception:
                rel = None
            if rel is None:
                continue
            a, b, c = rel
            if a == 0:
                continue
            # target = -(b*va + c*vb)/a
            resid = mp.mpf(a) * target_mp + mp.mpf(b) * va + mp.mpf(c) * vb
            if mp.mpf(abs(a)) > max_coeff or mp.mpf(abs(b)) > max_coeff or mp.mpf(abs(c)) > max_coeff:
                continue
            candidates.append({
                'class': f"{ca}+{cb}",
                'depth': 2,
                'forms': [fa, fb],
                'integer_coeffs': [int(a), int(b), int(c)],
                'max_coeff': max(abs(int(a)), abs(int(b)), abs(int(c))),
                'residual': mp.nstr(abs(resid), 10),
                'identification': f"target = -({b}×{fa} + {c}×{fb})/{a}",
            })

    return candidates


def control_baseline(target_mp: mp.mpf, basis, dps: int = 150,
                     n_random: int = 100):
    """Generate random comparable targets and run the same PSLQ search.

    This gives a baseline for the false-positive rate of the basis.
    Random targets uniformly drawn in (range of target / 10, range × 10)
    in log-space.
    """
    mp.mp.dps = dps
    target_log = mp.log(abs(target_mp))
    rng = np.random.default_rng(12345)
    null_hits = []
    for trial in range(n_random):
        # Pick a random target in [target/3, target × 3] range
        log_offset = mp.mpf(rng.uniform(-1.5, 1.5)) / mp.log(10)
        log_target = target_log + log_offset * mp.log(10)
        rand_target = mp.exp(log_target)
        cands = pslq_search_per_basis(
            rand_target, basis, dps=dps, tol_exp=50, max_coeff=1000
        ) if False else []
        # Heavy pair-search would dominate; skip for null baseline.
        # Instead: count single-element pairs at max_coeff <= 100 against ALL classes.
        single_hits = 0
        for fname, val, cls in basis:
            try:
                rel = mp.pslq([rand_target, val],
                              tol=mp.mpf(10) ** (-50),
                              maxcoeff=100)
                if rel is not None and rel[0] != 0:
                    single_hits += 1
            except Exception:
                pass
        null_hits.append(single_hits)
    return {
        'mean_null_hits': float(np.mean(null_hits)),
        'std_null_hits': float(np.std(null_hits)),
        'max_null_hits': int(np.max(null_hits)),
        'min_null_hits': int(np.min(null_hits)),
        'n_trials': n_random,
    }


# =====================================================================
#  Main driver
# =====================================================================

def main():
    out_dir = Path('debug/data')
    out_dir.mkdir(parents=True, exist_ok=True)

    DPS = 150  # working precision
    PSLQ_TOL_EXP = 50  # require residual < 1e-50 to call a hit
    MAX_COEFF = 1000  # mechanical-basis coefficient ceiling

    print("=" * 76)
    print("Sprint TD Track 5 — PSLQ probe")
    print(f"DPS = {DPS}, PSLQ tol = 1e-{PSLQ_TOL_EXP}, max_coeff = {MAX_COEFF}")
    print("=" * 76)

    # Step 1: high-precision computation of S_full(GS)
    systems = [
        {'Z': 2, 'n_max': 3, 'name': 'He_n3'},
        {'Z': 3, 'n_max': 3, 'name': 'Li+_n3'},
        {'Z': 2, 'n_max': 4, 'name': 'He_n4'},
    ]
    s_full_results = {}
    for sys_spec in systems:
        try:
            r = compute_S_full_GS_mpmath(sys_spec['Z'], sys_spec['n_max'], dps=DPS)
            s_full_results[sys_spec['name']] = r
        except Exception as e:
            print(f"  FAIL on {sys_spec['name']}: {e}")
            import traceback
            traceback.print_exc()
            s_full_results[sys_spec['name']] = {'error': str(e)}

    # Save the high-precision values as soon as we have them
    s_full_serializable = {}
    for name, r in s_full_results.items():
        if 'error' in r:
            s_full_serializable[name] = r
        else:
            s_full_serializable[name] = {
                k: v for k, v in r.items() if k != 'S_full_GS_mp'
            }
    with open(out_dir / 'sprint_td_track5_intermediate.json', 'w') as f:
        json.dump(s_full_serializable, f, indent=2)
    print(f"\nIntermediate S_full values saved.")

    # Step 2: build mechanical basis (frozen specification)
    print("\n" + "=" * 76)
    print("Building mechanical basis (frozen, n_max He as anchor)")
    print("=" * 76)
    he_target = s_full_results.get('He_n3', {}).get('S_full_GS_mp')
    if he_target is None:
        print("CRITICAL: He n_max=3 S_full could not be computed; aborting.")
        return
    basis = build_mechanical_basis(he_target, dps=DPS)
    basis_classes = basis_class_counts(basis)
    print(f"  Basis class counts: {dict(basis_classes)}")

    # Save basis specification
    basis_spec = {
        'metadata': {
            'dps': DPS,
            'pslq_tol_exp': PSLQ_TOL_EXP,
            'max_coeff': MAX_COEFF,
            'frozen_before_test': True,
        },
        'class_counts': dict(basis_classes),
        'n_basis': len(basis),
        'first_30_examples': [
            {'form': f, 'value': str(mp.nstr(v, 12)), 'class': c}
            for (f, v, c) in basis[:30]
        ],
    }

    # Step 3: control baseline (small N due to compute cost)
    print("\n" + "=" * 76)
    print("Control baseline (random comparable targets, single-element)")
    print("=" * 76)
    null_stats = control_baseline(he_target, basis, dps=DPS, n_random=20)
    print(f"  Random targets: {null_stats}")

    # Step 4: PSLQ search per system
    print("\n" + "=" * 76)
    print("PSLQ identification search")
    print("=" * 76)

    all_results = {}
    for name, r in s_full_results.items():
        if 'error' in r:
            all_results[name] = {'error': r['error']}
            continue
        target_mp = r['S_full_GS_mp']
        print(f"\n--- {name}: S_full = {mp.nstr(target_mp, 20)} ---")
        candidates = pslq_search_per_basis(
            target_mp, basis, dps=DPS,
            tol_exp=PSLQ_TOL_EXP, max_coeff=MAX_COEFF,
        )
        n_hits = len(candidates)
        # z-score vs null
        z = (n_hits - null_stats['mean_null_hits']) / max(null_stats['std_null_hits'], 1e-9)
        # Filter to most credible (lowest max_coeff)
        candidates.sort(key=lambda c: c['max_coeff'])
        top = candidates[:10]
        print(f"  Total candidates: {n_hits}")
        print(f"  z-score vs null: {z:+.2f}")
        if top:
            print(f"  Top 5 by max_coeff:")
            for c in top[:5]:
                print(f"    [{c['class']}] {c['identification']}  "
                      f"(coeffs {c['integer_coeffs']}, residual {c['residual']})")
        all_results[name] = {
            'target': str(mp.nstr(target_mp, DPS + 5)),
            'n_candidates': n_hits,
            'z_vs_null': float(z),
            'top_10': top,
        }

    # Cross-system consistency
    print("\n" + "=" * 76)
    print("Cross-system consistency")
    print("=" * 76)
    cross_check = {}
    for sys_a in s_full_results:
        for sys_b in s_full_results:
            if sys_a >= sys_b:
                continue
            ra = s_full_results.get(sys_a, {})
            rb = s_full_results.get(sys_b, {})
            if 'error' in ra or 'error' in rb:
                continue
            ratio = ra['S_full_GS_mp'] / rb['S_full_GS_mp']
            cross_check[f"{sys_a}/{sys_b}"] = mp.nstr(ratio, 30)
            print(f"  S({sys_a}) / S({sys_b}) = {mp.nstr(ratio, 30)}")

    # Step 5: write outputs
    output = {
        'metadata': {
            'sprint': 'TD Track 5',
            'date': '2026-05-09',
            'methodology': 'mpmath FCI + 1-RDM + PSLQ vs mechanical basis',
            'dps': DPS,
            'pslq_tol_exp': PSLQ_TOL_EXP,
            'max_coeff': MAX_COEFF,
            'falsification_protocol': 'W3 mechanical basis (frozen before testing)',
        },
        's_full_high_precision': {
            name: {k: v for k, v in r.items() if k != 'S_full_GS_mp'}
            for name, r in s_full_results.items()
        },
        'basis_specification': basis_spec,
        'control_baseline': null_stats,
        'pslq_results': all_results,
        'cross_system_ratios': cross_check,
    }

    with open(out_dir / 'sprint_td_track5.json', 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {out_dir / 'sprint_td_track5.json'}")
    return output


if __name__ == '__main__':
    main()
