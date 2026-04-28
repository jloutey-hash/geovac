"""
gn_f2_rational_function.py

Investigates F₂(t) as an explicit rational function of the hopping coupling t
in graph-native QED on the finite GeoVac Fock graph.

Four parts:
  1. n_max=2 symbolic: extract p(t)/q(t) exactly using sympy
  2. n_max=3: same if feasible; Born series polynomial approximation as fallback
  3. Convergence: C(n_max) × F₂(κ) vs α/(2π) for n_max=2,3 (and 4 if feasible)
  4. Born-series Taylor decomposition at t=0 (first 5 coefficients, n_max=2)

Scientific question:
  Is α the projection constant from graph-native scalar QED to vector QED?
  If C(n_max) × F₂(κ) → α/(2π) ≈ 0.001161, the answer is YES.

Output:
  debug/data/gn_f2_rational_function.json
  debug/gn_f2_rational_function_memo.md   (written by memo section at end)

Key references from CLAUDE.md / gn5_self_energy_memo.md:
  F₂(0) = 5√2/3 ≈ 2.357   (t=0, n_max=2, expected)
  F₂(κ) ≈ 2.3525           (t=-1/16, n_max=2, expected)
  C(n_max=3) = 50471424/1779441125  (projection exchange constant)
  α/(2π) ≈ 1.1614e-3       (Schwinger result)
"""

import json
import sys
import time
from pathlib import Path
from fractions import Fraction

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')

import sympy as sp
from sympy import Rational, sqrt, Symbol, nsimplify, cancel, simplify, factor, series
import numpy as np

# ---------------------------------------------------------------------------
# Add project root to path
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.graph_qed_propagator import DiracGraphOperator, KAPPA_SCALAR
from geovac.graph_qed_photon import build_fock_graph, compute_photon_propagator
from geovac.graph_qed_vertex import (
    build_projection_matrix, build_vertex_tensor, vertex_tensor_to_matrices
)
from geovac.graph_qed_continuum_bridge import (
    compute_projection_constant,
)

# Physical constant
ALPHA = 1.0 / 137.035999084
SCHWINGER = ALPHA / (2 * np.pi)   # α/(2π) ≈ 0.001161

# ---------------------------------------------------------------------------
# Helper: build photon propagator + vertex matrices for a given n_max
# ---------------------------------------------------------------------------
def build_qed_objects(n_max: int):
    """Return (G_gamma_sym, V_mats, V_bare, N_dirac, E_fock).

    G_gamma_sym : sympy Matrix  (photon propagator, t-independent)
    V_mats      : list of sympy Matrix, one per photon edge
    V_bare      : sympy Matrix, sum of V_mats
    N_dirac     : int
    E_fock      : int
    """
    # Photon propagator (t-independent, exact sympy)
    pp = compute_photon_propagator(n_max, exact=True)
    G_gamma_sym = pp.G_gamma if pp.G_gamma is not None else sp.Matrix(pp.G_gamma_numeric.tolist())

    # Vertex tensor — returns (entries, N_dirac, V_fock, E_fock)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(n_max)

    # vertex_tensor_to_matrices already returns sympy Matrices
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)
    V_bare = V_mats[0]
    for vm in V_mats[1:]:
        V_bare = V_bare + vm

    return G_gamma_sym, V_mats, V_bare, N_dirac, E_fock


# ---------------------------------------------------------------------------
# Part 1: n_max=2 symbolic rational function F₂(t)
# ---------------------------------------------------------------------------
def extract_f2_symbolic(n_max: int = 2, timeout_seconds: float = 300.0):
    """Extract F₂(t) = Tr(Lambda(t)) / Tr(V_bare · G_e(t)) as symbolic fraction.

    Returns dict with keys:
        f2_expr          : sympy expression (unsimplified ratio)
        f2_cancelled     : sympy expression after sp.cancel()
        f2_at_0          : F₂ evaluated at t=0
        f2_at_kappa      : F₂ evaluated at t=κ=-1/16
        numerator_poly   : coefficients of numerator (highest degree first)
        denominator_poly : coefficients of denominator
        num_degree       : degree of numerator
        den_degree       : degree of denominator
        number_field     : string description of field
        elapsed_seconds  : wall time
    """
    t = Symbol('t')
    print(f"[Part 1] n_max={n_max}: building QED objects...", flush=True)
    t0 = time.time()

    G_gamma_sym, V_mats, V_bare, N_dirac, E_fock = build_qed_objects(n_max)
    print(f"  N_dirac={N_dirac}, E_fock={E_fock}, elapsed={time.time()-t0:.1f}s", flush=True)

    # Build D(t) = Λ + t·A symbolically
    op_sym = DiracGraphOperator(n_max=n_max, t=t)
    D_sym = op_sym.matrix_sympy()
    print(f"  D(t) built, shape={D_sym.shape}, elapsed={time.time()-t0:.1f}s", flush=True)

    # Invert G_e(t) = D(t)^{-1}
    print(f"  Inverting {N_dirac}×{N_dirac} symbolic matrix...", flush=True)
    G_e_sym = D_sym.inv()
    print(f"  Inversion done, elapsed={time.time()-t0:.1f}s", flush=True)

    # Vertex correction Lambda(t) = Σ_{e',e''} G_γ[e',e''] · V_{e'} · G_e(t) · V_{e''}^T
    print(f"  Computing Lambda(t)...", flush=True)
    E = E_fock
    Lambda_total = sp.zeros(N_dirac, N_dirac)
    for e1 in range(E):
        for e2 in range(E):
            g_val = G_gamma_sym[e1, e2]
            if g_val == 0:
                continue
            contrib = g_val * V_mats[e1] * G_e_sym * V_mats[e2].T
            Lambda_total += contrib
    print(f"  Lambda(t) built, elapsed={time.time()-t0:.1f}s", flush=True)

    # Traces
    tr_lambda = Lambda_total.trace()
    norm_matrix = V_bare * G_e_sym
    tr_norm = norm_matrix.trace()

    # Simplify traces
    print(f"  Simplifying traces...", flush=True)
    tr_lambda_simplified = sp.nsimplify(sp.simplify(tr_lambda), rational=False)
    tr_norm_simplified = sp.nsimplify(sp.simplify(tr_norm), rational=False)
    print(f"  Tr(Lambda) simplified, elapsed={time.time()-t0:.1f}s", flush=True)

    # F₂(t) = Tr(Lambda) / Tr(V_bare · G_e)
    f2_expr = tr_lambda_simplified / tr_norm_simplified
    print(f"  Canceling F₂(t)...", flush=True)
    f2_cancelled = sp.cancel(f2_expr)
    print(f"  F₂(t) cancelled, elapsed={time.time()-t0:.1f}s", flush=True)

    # Evaluate at special points
    kappa_val = Rational(-1, 16)
    f2_at_0 = sp.nsimplify(f2_cancelled.subs(t, 0), rational=False)
    f2_at_kappa = sp.nsimplify(f2_cancelled.subs(t, kappa_val), rational=False)
    print(f"  F₂(0) = {f2_at_0}", flush=True)
    print(f"  F₂(κ) = {f2_at_kappa} ≈ {float(f2_at_kappa):.6f}", flush=True)

    # Extract polynomial degrees and coefficients
    # Write as ratio of polynomials in t
    numer, denom = sp.fraction(f2_cancelled)
    numer_poly = sp.Poly(numer, t)
    denom_poly = sp.Poly(denom, t)

    num_coeffs = [str(c) for c in numer_poly.all_coeffs()]
    den_coeffs = [str(c) for c in denom_poly.all_coeffs()]

    elapsed = time.time() - t0
    print(f"  Done. Elapsed={elapsed:.1f}s", flush=True)

    # Determine number field content (look for sqrt in expression)
    expr_str = str(f2_cancelled)
    has_sqrt2 = 'sqrt(2)' in expr_str or 'sqrt2' in expr_str
    has_sqrt3 = 'sqrt(3)' in expr_str
    has_sqrt6 = 'sqrt(6)' in expr_str
    fields = ['ℚ']
    if has_sqrt2:
        fields.append('√2')
    if has_sqrt3:
        fields.append('√3')
    if has_sqrt6:
        fields.append('√6')
    number_field = '(' + ', '.join(fields) + ')' if len(fields) > 1 else 'ℚ'

    return {
        'f2_expr_str': str(f2_expr),
        'f2_cancelled_str': str(f2_cancelled),
        'f2_at_0_str': str(f2_at_0),
        'f2_at_0_float': float(f2_at_0),
        'f2_at_kappa_str': str(f2_at_kappa),
        'f2_at_kappa_float': float(f2_at_kappa),
        'numerator_coefficients': num_coeffs,
        'denominator_coefficients': den_coeffs,
        'numerator_degree': numer_poly.degree(),
        'denominator_degree': denom_poly.degree(),
        'number_field': number_field,
        'N_dirac': N_dirac,
        'E_fock': E_fock,
        'elapsed_seconds': elapsed,
    }


# ---------------------------------------------------------------------------
# Part 2: n_max=3 (symbolic or Born series)
# ---------------------------------------------------------------------------
def extract_f2_nmax3_born(n_max: int = 3, born_order: int = 5):
    """Born-series polynomial approximation to F₂(t) at n_max=3.

    Uses G_e(t) ≈ Λ⁻¹ Σ_{k=0}^{born_order} (-t·Λ⁻¹A)^k

    Returns dict similar to Part 1 but with 'method': 'born_series'.
    """
    t_sym = Symbol('t')
    print(f"[Part 2] n_max={n_max}, Born order={born_order}...", flush=True)
    t0 = time.time()

    G_gamma_sym, V_mats, V_bare, N_dirac, E_fock = build_qed_objects(n_max)
    print(f"  N_dirac={N_dirac}, E_fock={E_fock}, elapsed={time.time()-t0:.1f}s", flush=True)

    # Build Λ and A (symbolic)
    op_sym = DiracGraphOperator(n_max=n_max, t=t_sym)
    D_full = op_sym.matrix_sympy()

    # Get diagonal part Λ (set t=0 in D_full to get just eigenvalues)
    D_t0 = op_sym.__class__(n_max=n_max, t=Rational(0)).matrix_sympy()
    Lambda_mat = D_t0   # diagonal

    # A = (D_full - D_t0) / t  — off-diagonal part / t
    # Since D_full = Λ + t*A, A = (D_full - Λ) / t
    A_mat = (D_full - Lambda_mat) / t_sym  # should be t-independent

    # Λ⁻¹: diagonal, invert element by element
    Lambda_inv = Lambda_mat.inv()   # diagonal → diagonal inversion, fast

    # Neumann sum: G_e(t) = Λ⁻¹ Σ_{k=0}^{K} (-t Λ⁻¹ A)^k
    # Let M = -t Λ⁻¹ A (linear in t_sym)
    M = -t_sym * Lambda_inv * A_mat
    print(f"  Built M = -t·Λ⁻¹·A, elapsed={time.time()-t0:.1f}s", flush=True)

    # Build Born series truncated polynomial in t
    # G_k = Λ⁻¹ * M^k
    G_sum = sp.eye(N_dirac)  # k=0 term (M^0 = I)
    M_power = sp.eye(N_dirac)
    for k in range(1, born_order + 1):
        M_power = M_power * M
        G_sum += M_power
        print(f"  Born k={k}, elapsed={time.time()-t0:.1f}s", flush=True)

    G_e_born = Lambda_inv * G_sum
    print(f"  G_e Born series built, elapsed={time.time()-t0:.1f}s", flush=True)

    # Vertex correction Lambda(t) with Born-approximated G_e
    Lambda_total = sp.zeros(N_dirac, N_dirac)
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_val = G_gamma_sym[e1, e2]
            if g_val == 0:
                continue
            contrib = g_val * V_mats[e1] * G_e_born * V_mats[e2].T
            Lambda_total += contrib
    print(f"  Lambda(t) built (Born), elapsed={time.time()-t0:.1f}s", flush=True)

    tr_lambda = Lambda_total.trace()
    norm_matrix = V_bare * G_e_born
    tr_norm = norm_matrix.trace()

    # Evaluate at t=0 and t=κ
    kappa_val = Rational(-1, 16)
    tr_lambda_t0 = sp.nsimplify(tr_lambda.subs(t_sym, 0), rational=False)
    tr_norm_t0 = sp.nsimplify(tr_norm.subs(t_sym, 0), rational=False)
    f2_at_0 = sp.nsimplify(tr_lambda_t0 / tr_norm_t0, rational=False)

    tr_lambda_kappa = sp.nsimplify(tr_lambda.subs(t_sym, kappa_val), rational=False)
    tr_norm_kappa = sp.nsimplify(tr_norm.subs(t_sym, kappa_val), rational=False)
    f2_at_kappa = sp.nsimplify(tr_lambda_kappa / tr_norm_kappa, rational=False)

    print(f"  F₂(0) = {f2_at_0} ≈ {float(f2_at_0):.6f}", flush=True)
    print(f"  F₂(κ) = {f2_at_kappa} ≈ {float(f2_at_kappa):.6f}", flush=True)

    elapsed = time.time() - t0
    print(f"  Done. Elapsed={elapsed:.1f}s", flush=True)

    return {
        'method': 'born_series',
        'born_order': born_order,
        'f2_at_0_str': str(f2_at_0),
        'f2_at_0_float': float(f2_at_0),
        'f2_at_kappa_str': str(f2_at_kappa),
        'f2_at_kappa_float': float(f2_at_kappa),
        'N_dirac': N_dirac,
        'E_fock': E_fock,
        'elapsed_seconds': elapsed,
    }


def extract_f2_nmax3_exact(n_max: int = 3, time_limit: float = 300.0):
    """Attempt exact symbolic inversion for n_max=3.

    If the matrix inversion + trace simplification takes > time_limit seconds,
    fall back to Born series.
    """
    t_sym = Symbol('t')
    print(f"[Part 2] n_max={n_max}: attempting exact symbolic...", flush=True)
    t0 = time.time()

    G_gamma_sym, V_mats, V_bare, N_dirac, E_fock = build_qed_objects(n_max)
    print(f"  N_dirac={N_dirac}, E_fock={E_fock}", flush=True)
    print(f"  Attempting {N_dirac}×{N_dirac} symbolic inversion...", flush=True)

    op_sym = DiracGraphOperator(n_max=n_max, t=t_sym)
    D_sym = op_sym.matrix_sympy()

    # Try inversion — this is the expensive step for n_max=3 (26×26)
    try:
        G_e_sym = D_sym.inv()
        inv_time = time.time() - t0
        print(f"  Inversion succeeded in {inv_time:.1f}s", flush=True)

        if inv_time > time_limit:
            print(f"  WARNING: inversion took {inv_time:.1f}s > {time_limit}s limit", flush=True)
    except Exception as exc:
        print(f"  Symbolic inversion failed: {exc}", flush=True)
        print(f"  Falling back to Born series...", flush=True)
        return extract_f2_nmax3_born(n_max=n_max)

    # If we get here, exact inversion succeeded
    Lambda_total = sp.zeros(N_dirac, N_dirac)
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_val = G_gamma_sym[e1, e2]
            if g_val == 0:
                continue
            Lambda_total += g_val * V_mats[e1] * G_e_sym * V_mats[e2].T

    tr_lambda = sp.simplify(Lambda_total.trace())
    tr_norm = sp.simplify((V_bare * G_e_sym).trace())

    f2_expr = tr_lambda / tr_norm
    f2_cancelled = sp.cancel(f2_expr)

    kappa_val = Rational(-1, 16)
    f2_at_0 = sp.nsimplify(f2_cancelled.subs(t_sym, 0), rational=False)
    f2_at_kappa = sp.nsimplify(f2_cancelled.subs(t_sym, kappa_val), rational=False)

    numer, denom = sp.fraction(f2_cancelled)
    numer_poly = sp.Poly(numer, t_sym)
    denom_poly = sp.Poly(denom, t_sym)

    elapsed = time.time() - t0
    print(f"  F₂(0) = {f2_at_0}", flush=True)
    print(f"  F₂(κ) ≈ {float(f2_at_kappa):.6f}", flush=True)
    print(f"  Done. Elapsed={elapsed:.1f}s", flush=True)

    return {
        'method': 'exact_symbolic',
        'f2_cancelled_str': str(f2_cancelled),
        'f2_at_0_str': str(f2_at_0),
        'f2_at_0_float': float(f2_at_0),
        'f2_at_kappa_str': str(f2_at_kappa),
        'f2_at_kappa_float': float(f2_at_kappa),
        'numerator_degree': numer_poly.degree(),
        'denominator_degree': denom_poly.degree(),
        'N_dirac': N_dirac,
        'E_fock': E_fock,
        'elapsed_seconds': elapsed,
    }


# ---------------------------------------------------------------------------
# Part 3: Convergence study C(n_max) × F₂(κ) vs α/(2π)
# ---------------------------------------------------------------------------
def convergence_study(n_max_values=None):
    """For each n_max, compute:
        F₂(κ) : anomalous moment at t=κ=-1/16
        C(n_max) : projection exchange constant
        product : C × F₂(κ)
        ratio : product / (α/2π)
    """
    if n_max_values is None:
        n_max_values = [2, 3]

    from geovac.graph_qed_self_energy import extract_anomalous_moment

    rows = []
    kappa_val = Rational(-1, 16)

    for n_max in n_max_values:
        print(f"[Part 3] n_max={n_max}...", flush=True)
        t0 = time.time()

        # F₂(κ) — use exact arithmetic
        try:
            am = extract_anomalous_moment(n_max=n_max, t=kappa_val, exact=True)
            f2_kappa_float = float(am['F2_float'])
            f2_kappa_str = str(am['F2'])
        except Exception as exc:
            print(f"  WARNING: exact F₂ failed ({exc})", flush=True)
            print(f"  Trying numeric fallback...", flush=True)
            try:
                am_num = extract_anomalous_moment(n_max=n_max,
                                                   t=float(kappa_val), exact=False)
                f2_kappa_float = float(am_num['F2_float'])
                f2_kappa_str = str(am_num['F2'])
            except Exception as exc2:
                print(f"  Both failed: {exc2}", flush=True)
                f2_kappa_float = float('nan')
                f2_kappa_str = "N/A"

        # Projection exchange constant C(n_max)
        try:
            bridge = compute_projection_constant(n_max_fock=n_max)
            C_val = bridge['projection_constant_float']
            C_str = bridge['projection_constant']
        except Exception as exc:
            print(f"  WARNING: C computation failed ({exc})", flush=True)
            C_val = None
            C_str = "N/A"

        if C_val is not None:
            product = C_val * f2_kappa_float
            ratio_to_schwinger = product / SCHWINGER
        else:
            product = None
            ratio_to_schwinger = None

        elapsed = time.time() - t0
        print(f"  F₂(κ)={f2_kappa_float:.6f}, C={C_str}, "
              f"product={product}, ratio={ratio_to_schwinger}", flush=True)
        print(f"  Elapsed={elapsed:.1f}s", flush=True)

        rows.append({
            'n_max': n_max,
            'f2_at_kappa': f2_kappa_float,
            'f2_at_kappa_str': f2_kappa_str,
            'C_projection': C_val,
            'C_projection_str': C_str,
            'product_C_times_F2': product,
            'schwinger': SCHWINGER,
            'ratio_to_schwinger': ratio_to_schwinger,
            'elapsed_seconds': elapsed,
        })

    return rows


# ---------------------------------------------------------------------------
# Part 4: Born series / Taylor series decomposition at t=0, n_max=2
# ---------------------------------------------------------------------------
def taylor_decomposition(n_max: int = 2, n_terms: int = 6):
    """Expand F₂(t) as Taylor series around t=0.

    Returns first n_terms coefficients [c_0, c_1, c_2, ...] where
    F₂(t) = c_0 + c_1*t + c_2*t^2 + ...

    For the rational function from Part 1, this is exact.
    We use the symbolic F₂(t) from Part 1 and differentiate.
    """
    t = Symbol('t')
    print(f"[Part 4] n_max={n_max}, Taylor series to order {n_terms}...", flush=True)
    t0_time = time.time()

    # We need to recompute F₂(t) (or use cached result)
    # For efficiency, build the propagator and compute F₂(t) as series directly
    G_gamma_sym, V_mats, V_bare, N_dirac, E_fock = build_qed_objects(n_max)

    op_sym = DiracGraphOperator(n_max=n_max, t=t)
    D_sym = op_sym.matrix_sympy()
    G_e_sym = D_sym.inv()

    Lambda_total = sp.zeros(N_dirac, N_dirac)
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_val = G_gamma_sym[e1, e2]
            if g_val == 0:
                continue
            Lambda_total += g_val * V_mats[e1] * G_e_sym * V_mats[e2].T

    tr_lambda = Lambda_total.trace()
    tr_norm = (V_bare * G_e_sym).trace()

    # Simplify
    tr_lambda_s = sp.simplify(tr_lambda)
    tr_norm_s = sp.simplify(tr_norm)

    f2_expr = sp.cancel(tr_lambda_s / tr_norm_s)

    # Taylor series
    f2_series = sp.series(f2_expr, t, 0, n=n_terms)
    print(f"  Series: {f2_series}", flush=True)

    # Extract coefficients
    coefficients = []
    for k in range(n_terms):
        coeff = f2_series.coeff(t, k)
        coeff_simplified = sp.nsimplify(coeff, rational=False)

        # Classify
        coeff_float = float(coeff_simplified.evalf()) if coeff_simplified != 0 else 0.0

        # Check if rational
        try:
            coeff_rational = sp.Rational(coeff_simplified)
            is_rational = True
            classification = 'rational'
        except (TypeError, ValueError):
            is_rational = False
            # Check for sqrt content
            expr_str = str(coeff_simplified)
            if 'sqrt' in expr_str:
                classification = 'algebraic (contains sqrt)'
            else:
                classification = 'algebraic'

        print(f"  c_{k} = {coeff_simplified} ≈ {coeff_float:.8g} [{classification}]",
              flush=True)
        coefficients.append({
            'order': k,
            'coefficient_str': str(coeff_simplified),
            'coefficient_float': coeff_float,
            'classification': classification,
        })

    elapsed = time.time() - t0_time
    print(f"  Done. Elapsed={elapsed:.1f}s", flush=True)

    return {
        'series_str': str(f2_series),
        'coefficients': coefficients,
        'n_terms': n_terms,
        'elapsed_seconds': elapsed,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    results = {}

    print("=" * 70, flush=True)
    print("GN F2 Rational Function Investigation", flush=True)
    print("=" * 70, flush=True)
    print(f"Schwinger α/(2π) = {SCHWINGER:.8e}", flush=True)
    print(flush=True)

    # -----------------------------------------------------------------------
    # Part 1: n_max=2 exact symbolic
    # -----------------------------------------------------------------------
    print("Part 1: n_max=2 exact symbolic F₂(t)", flush=True)
    print("-" * 50, flush=True)
    try:
        part1 = extract_f2_symbolic(n_max=2)
        results['part1_nmax2_symbolic'] = part1
        print(f"  F₂(0) = {part1['f2_at_0_str']} ≈ {part1['f2_at_0_float']:.6f}", flush=True)
        print(f"  F₂(κ) = {part1['f2_at_kappa_str']} ≈ {part1['f2_at_kappa_float']:.6f}",
              flush=True)
        print(f"  Numerator degree: {part1['numerator_degree']}", flush=True)
        print(f"  Denominator degree: {part1['denominator_degree']}", flush=True)
        print(f"  Number field: {part1['number_field']}", flush=True)
    except Exception as exc:
        print(f"  FAILED: {exc}", flush=True)
        results['part1_nmax2_symbolic'] = {'error': str(exc)}

    print(flush=True)

    # -----------------------------------------------------------------------
    # Part 4: Taylor series at n_max=2 (do before Part 2 since it reuses
    # the n_max=2 infrastructure and avoids recomputing)
    # -----------------------------------------------------------------------
    print("Part 4: Taylor series decomposition (n_max=2)", flush=True)
    print("-" * 50, flush=True)
    try:
        part4 = taylor_decomposition(n_max=2, n_terms=6)
        results['part4_taylor_nmax2'] = part4
    except Exception as exc:
        print(f"  FAILED: {exc}", flush=True)
        results['part4_taylor_nmax2'] = {'error': str(exc)}

    print(flush=True)

    # -----------------------------------------------------------------------
    # Part 2: n_max=3
    # -----------------------------------------------------------------------
    print("Part 2: n_max=3 (exact attempt, Born series fallback)", flush=True)
    print("-" * 50, flush=True)
    try:
        # Try exact for up to 120 seconds, otherwise Born order 4
        part2 = extract_f2_nmax3_exact(n_max=3, time_limit=120.0)
        results['part2_nmax3'] = part2
    except Exception as exc:
        print(f"  Exact failed ({exc}), trying Born series...", flush=True)
        try:
            part2 = extract_f2_nmax3_born(n_max=3, born_order=4)
            results['part2_nmax3'] = part2
        except Exception as exc2:
            print(f"  Born series also failed: {exc2}", flush=True)
            results['part2_nmax3'] = {'error': str(exc2)}

    print(flush=True)

    # -----------------------------------------------------------------------
    # Part 3: Convergence study
    # -----------------------------------------------------------------------
    print("Part 3: Convergence C(n_max) × F₂(κ) vs α/(2π)", flush=True)
    print("-" * 50, flush=True)
    try:
        part3 = convergence_study(n_max_values=[2, 3])
        results['part3_convergence'] = part3
        print(f"  α/(2π) = {SCHWINGER:.8e}", flush=True)
        for row in part3:
            prod = row.get('product_C_times_F2')
            ratio = row.get('ratio_to_schwinger')
            print(f"  n_max={row['n_max']}: F₂={row['f2_at_kappa']:.5f}, "
                  f"C={row['C_projection_str'][:30]}, "
                  f"C×F₂={prod}, ratio={ratio}", flush=True)
    except Exception as exc:
        print(f"  FAILED: {exc}", flush=True)
        results['part3_convergence'] = {'error': str(exc)}

    print(flush=True)

    # -----------------------------------------------------------------------
    # Save JSON
    # -----------------------------------------------------------------------
    out_path = PROJECT_ROOT / 'debug' / 'data' / 'gn_f2_rational_function.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Results saved to {out_path}", flush=True)

    return results


if __name__ == '__main__':
    results = main()
