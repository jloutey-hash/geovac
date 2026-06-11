"""
f2_rational_nmax3.py

Task 3: Extract the exact symbolic rational function F₂(t) = p(t)/q(t)
at n_max=3 using sympy. The 28×28 Dirac matrix inversion is feasible
(previously completed in ~6.4s).

Extracts:
  - Degrees of p(t) and q(t)
  - Algebraic number field of coefficients
  - First 6 Taylor coefficients c₀, c₂, c₄, c₆, c₈, c₁₀
    (odd ones zero by t→-t symmetry)
  - Comparison of n_max=2 and n_max=3 rational functions

Output: debug/data/f2_rational_nmax3.json
"""

import json
import sys
import time
from pathlib import Path

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')

import sympy as sp
from sympy import Rational, Symbol, sqrt, cancel, Poly, series, fraction

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.graph_qed_propagator import DiracGraphOperator, KAPPA_SCALAR
from geovac.graph_qed_photon import build_fock_graph, compute_photon_propagator
from geovac.graph_qed_vertex import (
    build_projection_matrix, build_vertex_tensor, vertex_tensor_to_matrices
)


# ---------------------------------------------------------------------------
# Helper: build t-independent photon propagator + vertex matrices
# ---------------------------------------------------------------------------
def build_qed_objects(n_max: int):
    """Return (G_gamma_sym, V_mats, V_bare, N_dirac, E_fock)."""
    pp = compute_photon_propagator(n_max, exact=True)
    G_gamma_sym = pp.G_gamma if pp.G_gamma is not None else sp.Matrix(pp.G_gamma_numeric.tolist())

    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(n_max)
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)
    V_bare = V_mats[0]
    for vm in V_mats[1:]:
        V_bare = V_bare + vm

    return G_gamma_sym, V_mats, V_bare, N_dirac, E_fock


# ---------------------------------------------------------------------------
# Extract F₂(t) as exact rational function
# ---------------------------------------------------------------------------
def extract_f2_rational(n_max: int):
    """Extract F₂(t) = Tr(Lambda(t)) / Tr(V_bare · G_e(t)) exactly."""
    t = Symbol('t')
    t0 = time.time()
    print(f"n_max={n_max}: building QED objects...", flush=True)

    G_gamma_sym, V_mats, V_bare, N_dirac, E_fock = build_qed_objects(n_max)
    print(f"  N_dirac={N_dirac}, E_fock={E_fock}, elapsed={time.time()-t0:.1f}s", flush=True)

    # Build D(t) = Lambda + t*A symbolically
    op_sym = DiracGraphOperator(n_max=n_max, t=t)
    D_sym = op_sym.matrix_sympy()
    print(f"  D(t) built, shape={D_sym.shape}, elapsed={time.time()-t0:.1f}s", flush=True)

    # Invert symbolically
    print(f"  Inverting {N_dirac}x{N_dirac} symbolic matrix...", flush=True)
    G_e_sym = D_sym.inv()
    inv_time = time.time() - t0
    print(f"  Inversion done in {inv_time:.1f}s", flush=True)

    # Vertex correction Lambda(t)
    print(f"  Computing Lambda(t) ({E_fock}^2 = {E_fock**2} terms)...", flush=True)
    Lambda_total = sp.zeros(N_dirac, N_dirac)
    count = 0
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_val = G_gamma_sym[e1, e2]
            if g_val == 0:
                continue
            Lambda_total += g_val * V_mats[e1] * G_e_sym * V_mats[e2].T
            count += 1
    print(f"  Lambda(t) built ({count} nonzero entries), elapsed={time.time()-t0:.1f}s", flush=True)

    # Compute traces
    print(f"  Computing traces...", flush=True)
    tr_lambda = Lambda_total.trace()
    tr_norm = (V_bare * G_e_sym).trace()
    print(f"  Traces computed, elapsed={time.time()-t0:.1f}s", flush=True)

    # Cancel to get p(t)/q(t)
    print(f"  Canceling F2(t) = Tr(Lambda)/Tr(V_bare*G_e)...", flush=True)
    f2_expr = cancel(tr_lambda / tr_norm)
    print(f"  Cancel done, elapsed={time.time()-t0:.1f}s", flush=True)

    # Extract numerator and denominator polynomials
    numer, denom = fraction(f2_expr)
    numer_poly = Poly(sp.expand(numer), t)
    denom_poly = Poly(sp.expand(denom), t)

    num_degree = numer_poly.degree()
    den_degree = denom_poly.degree()
    print(f"  F2(t) = p_{num_degree}(t) / q_{den_degree}(t)", flush=True)

    # Evaluate at special points
    kappa_val = Rational(-1, 16)
    f2_at_0 = f2_expr.subs(t, 0)
    f2_at_kappa = f2_expr.subs(t, kappa_val)
    f2_at_0_simplified = sp.nsimplify(f2_at_0, rational=False)
    f2_at_kappa_simplified = sp.nsimplify(f2_at_kappa, rational=False)

    print(f"  F2(0) = {f2_at_0_simplified} ~ {float(f2_at_0_simplified):.10f}", flush=True)
    print(f"  F2(kappa) = {f2_at_kappa_simplified} ~ {float(f2_at_kappa_simplified):.10f}", flush=True)

    # Taylor coefficients: F2(t) = c0 + c1*t + c2*t^2 + ...
    # F2(t) is EVEN so c1=c3=c5=... = 0
    print(f"  Extracting Taylor coefficients...", flush=True)
    # Use series expansion
    f2_series = series(f2_expr, t, 0, n=12)
    taylor_coeffs = {}
    for k in range(12):
        coeff = f2_series.coeff(t, k)
        if coeff != 0:
            coeff_simplified = sp.nsimplify(coeff, rational=False)
            taylor_coeffs[k] = {
                'value_str': str(coeff_simplified),
                'value_float': float(coeff_simplified),
            }
            print(f"    c_{k} = {coeff_simplified} ~ {float(coeff_simplified):.10f}", flush=True)
        else:
            taylor_coeffs[k] = {'value_str': '0', 'value_float': 0.0}

    # Number field analysis: scan all coefficients for algebraic content
    print(f"  Analyzing number field of coefficients...", flush=True)
    all_coeffs_str = str(f2_expr)
    irrationals_found = set()
    # Check for sqrt(k) for k=2..20
    for k in range(2, 21):
        if f'sqrt({k})' in all_coeffs_str:
            irrationals_found.add(f'sqrt({k})')
    # Also check in individual polynomial coefficients
    for c in numer_poly.all_coeffs():
        c_str = str(c)
        for k in range(2, 21):
            if f'sqrt({k})' in c_str:
                irrationals_found.add(f'sqrt({k})')
    for c in denom_poly.all_coeffs():
        c_str = str(c)
        for k in range(2, 21):
            if f'sqrt({k})' in c_str:
                irrationals_found.add(f'sqrt({k})')

    number_field = 'Q' if not irrationals_found else f'Q({", ".join(sorted(irrationals_found))})'
    print(f"  Number field: {number_field}", flush=True)

    # Extract polynomial coefficients as strings
    num_coeffs = [str(c) for c in numer_poly.all_coeffs()]
    den_coeffs = [str(c) for c in denom_poly.all_coeffs()]

    elapsed = time.time() - t0
    print(f"  Total elapsed: {elapsed:.1f}s", flush=True)

    return {
        'n_max': n_max,
        'N_dirac': N_dirac,
        'E_fock': E_fock,
        'method': 'exact_symbolic',
        'f2_cancelled_str': str(f2_expr),
        'numerator_degree': num_degree,
        'denominator_degree': den_degree,
        'numerator_coefficients': num_coeffs,
        'denominator_coefficients': den_coeffs,
        'f2_at_0_str': str(f2_at_0_simplified),
        'f2_at_0_float': float(f2_at_0_simplified),
        'f2_at_kappa_str': str(f2_at_kappa_simplified),
        'f2_at_kappa_float': float(f2_at_kappa_simplified),
        'taylor_coefficients': taylor_coeffs,
        'number_field': number_field,
        'irrationals_found': sorted(list(irrationals_found)),
        'elapsed_seconds': elapsed,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    results = {
        'description': (
            'Exact symbolic F2(t) rational function extraction at n_max=2 and n_max=3. '
            'Taylor coefficients, polynomial degrees, and number field analysis.'
        ),
    }

    print("=" * 70, flush=True)
    print("Task 3: Exact Symbolic F2(t) Rational Function", flush=True)
    print("=" * 70, flush=True)

    # n_max=2 (quick sanity check)
    print("\n--- n_max=2 ---", flush=True)
    r2 = extract_f2_rational(n_max=2)
    results['n_max_2'] = r2
    print(f"\nn_max=2 summary:", flush=True)
    print(f"  Degree: {r2['numerator_degree']}/{r2['denominator_degree']}", flush=True)
    print(f"  F2(0) = {r2['f2_at_0_str']} ~ {r2['f2_at_0_float']:.10f}", flush=True)
    print(f"  F2(kappa) ~ {r2['f2_at_kappa_float']:.10f}", flush=True)
    print(f"  Number field: {r2['number_field']}", flush=True)

    # n_max=3 (the main target)
    print("\n--- n_max=3 ---", flush=True)
    r3 = extract_f2_rational(n_max=3)
    results['n_max_3'] = r3
    print(f"\nn_max=3 summary:", flush=True)
    print(f"  Degree: {r3['numerator_degree']}/{r3['denominator_degree']}", flush=True)
    print(f"  F2(0) = {r3['f2_at_0_str']} ~ {r3['f2_at_0_float']:.10f}", flush=True)
    print(f"  F2(kappa) ~ {r3['f2_at_kappa_float']:.10f}", flush=True)
    print(f"  Number field: {r3['number_field']}", flush=True)

    # Compare
    print("\n" + "=" * 70, flush=True)
    print("Comparison", flush=True)
    print("-" * 70, flush=True)
    print(f"  n_max=2: deg {r2['numerator_degree']}/{r2['denominator_degree']}, "
          f"field {r2['number_field']}", flush=True)
    print(f"  n_max=3: deg {r3['numerator_degree']}/{r3['denominator_degree']}, "
          f"field {r3['number_field']}", flush=True)

    # Taylor coefficient comparison
    print("\nTaylor coefficients (even only, odd vanish):", flush=True)
    print(f"  {'k':>3} {'c_k (n=2)':>20} {'c_k (n=3)':>20}", flush=True)
    for k in range(0, 12, 2):
        c2 = r2['taylor_coefficients'].get(k, {}).get('value_float', 0.0)
        c3 = r3['taylor_coefficients'].get(k, {}).get('value_float', 0.0)
        print(f"  {k:>3} {c2:>20.10f} {c3:>20.10f}", flush=True)

    # Save
    out_path = PROJECT_ROOT / 'debug' / 'data' / 'f2_rational_nmax3.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}", flush=True)

    return results


if __name__ == '__main__':
    main()
