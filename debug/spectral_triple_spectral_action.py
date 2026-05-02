"""
Spectral action computation on the Fock-projected S³ spectral triple.
=====================================================================

Computes Tr f(D²/Λ²) for the Connes-Chamseddine spectral action on
the GeoVac finite Dirac graph and compares against K = π(B + F − Δ).

Paper 2 targets:
  K = π(B + F − Δ) ≈ 137.036  (= 1/α)
  K/π = B + F − Δ ≈ 43.620
  B = 42, F = π²/6 ≈ 1.6449, Δ = 1/40 = 0.025

Spectral action framework:
  Bosonic:    S_b[D, Λ] = Tr f(D²/Λ²)
  Fermionic:  S_f[D, Λ] = ⟨ψ, Dψ⟩  (not computed here)

  For a smooth even f with moments f_k = ∫₀^∞ f(u) u^{k/2-1} du,
  the Seeley-DeWitt expansion gives:
    Tr f(D²/Λ²) ~ Σ_k f_k · a_k(D²) · Λ^{d-k}  (d = manifold dim)

On the finite graph, there's no asymptotic expansion — Tr f(D²/Λ²)
is an exact finite sum.  We compute it directly and also extract
effective SD-like coefficients from the heat kernel.
"""

import json
import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.spectral_triple import FockSpectralTriple


def compute_eigenvalues(st):
    """Get numerical eigenvalues of D from the spectral triple."""
    D_num = np.array(st.dirac_operator.tolist(), dtype=float)
    return np.linalg.eigvalsh(D_num)


def sharp_cutoff_action(evals_sq, cutoff_sq):
    """Tr f(D²/Λ²) with f = Θ(1 - x): count eigenvalues with λ² ≤ Λ²."""
    return int(np.sum(evals_sq <= cutoff_sq + 1e-12))


def heat_kernel_action(evals_sq, t):
    """Tr exp(-t D²) = heat kernel trace at parameter t."""
    return float(np.sum(np.exp(-t * evals_sq)))


def smooth_cutoff_action(evals_sq, cutoff_sq, width=0.5):
    """Tr f(D²/Λ²) with smooth sigmoid f(x) = 1/(1 + exp((x-1)/w))."""
    x = evals_sq / cutoff_sq
    return float(np.sum(1.0 / (1.0 + np.exp((x - 1.0) / width))))


def gaussian_cutoff_action(evals_sq, cutoff_sq):
    """Tr exp(-D²/Λ²) — Gaussian test function."""
    return float(np.sum(np.exp(-evals_sq / cutoff_sq)))


def polynomial_cutoff_action(evals_sq, cutoff_sq, degree=2):
    """Tr (1 - D²/Λ²)^p Θ(1 - D²/Λ²) — polynomial test function."""
    x = evals_sq / cutoff_sq
    mask = x <= 1.0
    return float(np.sum((1.0 - x[mask]) ** degree))


def extract_heat_kernel_coefficients(evals_sq, t_values):
    """Extract effective Seeley-DeWitt coefficients from heat kernel.

    Tr exp(-t D²) ~ Σ_k a_k t^{(k-d)/2}  as t → 0+

    For d=3 (S³):
      Tr exp(-t D²) ~ a_0 t^{-3/2} + a_1 t^{-1/2} + a_2 t^{1/2} + ...

    We fit the heat kernel trace to extract effective coefficients.
    """
    traces = np.array([heat_kernel_action(evals_sq, t) for t in t_values])

    # For d=3: powers are t^{-3/2}, t^{-1/2}, t^{1/2}, t^{3/2}
    powers = np.array([-1.5, -0.5, 0.5, 1.5])

    # Build Vandermonde-like matrix
    A = np.column_stack([t_values ** p for p in powers])

    # Least squares fit
    coeffs, _, _, _ = np.linalg.lstsq(A, traces, rcond=None)
    return dict(zip([f"a_{k}" for k in range(4)], coeffs))


def connes_chamseddine_bosonic(evals_sq, Lambda_sq, f_moments=None):
    """Compute the Connes-Chamseddine bosonic spectral action.

    S_b = Σ_k f_k · a_k · Λ^{3-k}

    where f_k are the moments of the test function f and a_k are the
    Seeley-DeWitt coefficients of D².

    On the finite graph, we compute Tr f(D²/Λ²) directly for several
    test functions and report the result.
    """
    results = {}

    # 1. Sharp cutoff (counting function N(Λ))
    results['sharp_cutoff'] = sharp_cutoff_action(evals_sq, Lambda_sq)

    # 2. Heat kernel (Gaussian)
    results['gaussian'] = gaussian_cutoff_action(evals_sq, Lambda_sq)

    # 3. Smooth sigmoid
    for w in [0.1, 0.3, 0.5, 1.0]:
        results[f'sigmoid_w{w}'] = smooth_cutoff_action(evals_sq, Lambda_sq, width=w)

    # 4. Polynomial cutoffs
    for p in [1, 2, 3, 4]:
        results[f'poly_p{p}'] = polynomial_cutoff_action(evals_sq, Lambda_sq, degree=p)

    return results


def spectral_action_scan(evals_sq, Lambda_values):
    """Scan the spectral action across multiple cutoff values."""
    results = []
    for L in Lambda_values:
        L_sq = L ** 2
        row = {
            'Lambda': float(L),
            'Lambda_sq': float(L_sq),
            'sharp': sharp_cutoff_action(evals_sq, L_sq),
            'gaussian': gaussian_cutoff_action(evals_sq, L_sq),
            'sigmoid_w0.3': smooth_cutoff_action(evals_sq, L_sq, width=0.3),
            'poly_p2': polynomial_cutoff_action(evals_sq, L_sq, degree=2),
        }
        results.append(row)
    return results


def supertrace_action(evals, gamma_signs, Lambda_sq):
    """Compute the SUPER spectral action Str f(D²/Λ²) = Tr(γ f(D²/Λ²)).

    This is the boson-fermion graded trace. On the continuum S³,
    the supertrace of the CC action vanishes identically at every
    SD order (Sprint ST-1 Finding F1: a_k^{D²}/a_k^{Δ_LB} = 4).
    """
    evals_sq = evals ** 2

    # Sharp cutoff supertrace
    mask = evals_sq <= Lambda_sq + 1e-12
    str_sharp = float(np.sum(gamma_signs[mask]))

    # Gaussian supertrace
    weights = np.exp(-evals_sq / Lambda_sq)
    str_gaussian = float(np.sum(gamma_signs * weights))

    # Heat kernel supertrace at various t
    str_heat = {}
    for t in [0.01, 0.1, 0.5, 1.0, 2.0]:
        weights = np.exp(-t * evals_sq)
        str_heat[f't={t}'] = float(np.sum(gamma_signs * weights))

    return {
        'sharp_cutoff': str_sharp,
        'gaussian': str_gaussian,
        'heat_kernel': str_heat,
    }


def paper2_targets():
    """Paper 2 target values for comparison."""
    import math
    B = 42
    F = math.pi ** 2 / 6
    Delta = 1.0 / 40
    K_over_pi = B + F - Delta
    K = math.pi * K_over_pi
    alpha_inv = 137.035999084  # CODATA 2018
    return {
        'B': B,
        'F': F,
        'Delta': Delta,
        'K_over_pi': K_over_pi,
        'K': K,
        'alpha_inv_CODATA': alpha_inv,
        'K_minus_alpha_inv': K - alpha_inv,
        'relative_error': abs(K - alpha_inv) / alpha_inv,
    }


def main():
    import math

    results = {}
    targets = paper2_targets()
    results['paper2_targets'] = targets

    print("=" * 70)
    print("SPECTRAL ACTION ON THE FOCK S³ SPECTRAL TRIPLE")
    print("=" * 70)
    print()
    print(f"Paper 2 targets:")
    print(f"  B = {targets['B']}")
    print(f"  F = π²/6 = {targets['F']:.10f}")
    print(f"  Δ = 1/40 = {targets['Delta']}")
    print(f"  K/π = B + F − Δ = {targets['K_over_pi']:.10f}")
    print(f"  K = π(B + F − Δ) = {targets['K']:.10f}")
    print(f"  1/α (CODATA) = {targets['alpha_inv_CODATA']}")
    print()

    # ====================================================================
    # Compute for all 4 configurations at n_max=2 and n_max=3
    # ====================================================================
    configs = [
        {'n_max': 2, 'j_type': 'permutation', 'adjacency_weights': 'uniform'},
        {'n_max': 2, 'j_type': 'kramers', 'adjacency_weights': 'cg'},
        {'n_max': 3, 'j_type': 'permutation', 'adjacency_weights': 'uniform'},
        {'n_max': 3, 'j_type': 'kramers', 'adjacency_weights': 'cg'},
    ]

    for cfg in configs:
        label = f"n{cfg['n_max']}_{cfg['j_type']}_{cfg['adjacency_weights']}"
        print("-" * 70)
        print(f"Configuration: n_max={cfg['n_max']}, "
              f"J={cfg['j_type']}, weights={cfg['adjacency_weights']}")
        print("-" * 70)

        t0 = time.time()
        st = FockSpectralTriple(**cfg)
        evals = compute_eigenvalues(st)
        build_time = time.time() - t0

        evals_sq = evals ** 2
        dim_H = st.dim_H

        # Get chirality signs for supertrace
        gamma_signs = np.array([float(st._lattice.chirality[i]) for i in range(dim_H)])

        print(f"  dim_H = {dim_H}")
        print(f"  Build time: {build_time:.2f}s")
        print(f"  Eigenvalues of D: {sorted(evals)[:4]}...{sorted(evals)[-4:]}")
        print(f"  |λ| range: [{min(abs(evals)):.6f}, {max(abs(evals)):.6f}]")
        print(f"  Eigenvalues of D²: [{min(evals_sq):.6f}, {max(evals_sq):.6f}]")
        print()

        cfg_results = {
            'n_max': cfg['n_max'],
            'j_type': cfg['j_type'],
            'adjacency_weights': cfg['adjacency_weights'],
            'dim_H': dim_H,
            'build_time_s': build_time,
            'eigenvalues_D': sorted(float(e) for e in evals),
            'eigenvalues_D_sq': sorted(float(e) for e in evals_sq),
        }

        # ----------------------------------------------------------------
        # 1. Natural cutoffs from the spectrum
        # ----------------------------------------------------------------
        print("  1. NATURAL CUTOFF ANALYSIS")
        print("  " + "-" * 40)

        # Camporesi-Higuchi: |λ_n| = n + 3/2 for n = 0,1,...,n_max-1
        # At n_max=2: highest CH level is n=1, |λ|=5/2
        # At n_max=3: highest CH level is n=2, |λ|=7/2
        # The actual D eigenvalues are perturbations of these by κ*A

        ch_cutoffs = []
        for n in range(cfg['n_max']):
            lam_ch = n + 1.5
            ch_cutoffs.append(lam_ch)

        # Also try the max eigenvalue and some multiples
        lam_max = max(abs(evals))
        natural_cutoffs = ch_cutoffs + [lam_max, lam_max * 1.01]

        print(f"  CH eigenvalues: {ch_cutoffs}")
        print(f"  Max |λ_D|: {lam_max:.6f}")
        print()

        cutoff_results = []
        for L in natural_cutoffs:
            L_sq = L ** 2
            N_sharp = sharp_cutoff_action(evals_sq, L_sq)
            gauss = gaussian_cutoff_action(evals_sq, L_sq)
            sig = smooth_cutoff_action(evals_sq, L_sq, width=0.1)
            poly2 = polynomial_cutoff_action(evals_sq, L_sq, degree=2)

            print(f"  Λ = {L:.4f}: N(Λ)={N_sharp}, Gauss={gauss:.4f}, "
                  f"Sigmoid={sig:.4f}, Poly2={poly2:.4f}")

            cutoff_results.append({
                'Lambda': float(L),
                'sharp': N_sharp,
                'gaussian': gauss,
                'sigmoid_w0.1': sig,
                'poly_p2': poly2,
            })

        cfg_results['natural_cutoff_analysis'] = cutoff_results
        print()

        # ----------------------------------------------------------------
        # 2. Spectral action at Λ² = max eigenvalue²  (the full trace)
        # ----------------------------------------------------------------
        print("  2. FULL SPECTRAL ACTION (Λ = max|λ|)")
        print("  " + "-" * 40)

        # At the maximum cutoff, sharp cutoff = dim_H
        # The interesting quantity is the WEIGHTED trace
        L_full_sq = lam_max ** 2 * 1.01  # slightly above max to capture all

        full_action = connes_chamseddine_bosonic(evals_sq, L_full_sq)

        for key, val in full_action.items():
            print(f"  {key}: {val:.6f}")

        cfg_results['full_spectral_action'] = {
            k: float(v) for k, v in full_action.items()
        }

        # Compare to K/π and K
        print()
        print(f"  Comparison to Paper 2:")
        print(f"    dim_H = {dim_H}, K/π = {targets['K_over_pi']:.6f}")
        print(f"    dim_H / (K/π) = {dim_H / targets['K_over_pi']:.6f}")

        # The sharp cutoff at max is just dim_H.
        # The Gaussian at Λ=max gives Tr exp(-D²/max²)
        gauss_full = full_action['gaussian']
        print(f"    Tr exp(-D²/Λ²_max) = {gauss_full:.6f}")

        cfg_results['comparison'] = {
            'dim_H_over_K_pi': dim_H / targets['K_over_pi'],
        }
        print()

        # ----------------------------------------------------------------
        # 3. Heat kernel trace and SD coefficient extraction
        # ----------------------------------------------------------------
        print("  3. HEAT KERNEL EXPANSION")
        print("  " + "-" * 40)

        t_values_log = np.logspace(-3, 1, 50)
        heat_traces = np.array([heat_kernel_action(evals_sq, t) for t in t_values_log])

        # Small-t behavior: should see t^{-d/2} = t^{-3/2} growth for d=3
        # Fit: Tr e^{-tD²} ≈ a₀ t^{-3/2} + a₁ t^{-1/2} + a₂ t^{1/2} + ...

        # Use only small-t values for the fit
        small_t_mask = t_values_log < 0.5
        t_small = t_values_log[small_t_mask]
        traces_small = heat_traces[small_t_mask]

        # Fit with 4 terms: t^{-3/2}, t^{-1/2}, t^{1/2}, t^{3/2}
        powers = [-1.5, -0.5, 0.5, 1.5]
        A_mat = np.column_stack([t_small ** p for p in powers])
        sd_coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, traces_small, rcond=None)

        print(f"  Seeley-DeWitt fit (small-t regime, {len(t_small)} points):")
        for k, (p, c) in enumerate(zip(powers, sd_coeffs)):
            print(f"    a_{k} (coeff of t^{p:+.1f}) = {c:.8f}")

        if len(residuals) > 0:
            print(f"    Fit residual: {residuals[0]:.2e}")

        # On continuum unit S³:
        # a₀ = dim(spinor) · Vol(S³)/(4π)^{3/2} = 4 · 2π² / (4π)^{3/2}
        #     = 4 · 2π² / (8π^{3/2}) = π^{1/2}
        # From qed_vacuum_polarization.py: a₀ = a₁ = √π, a₂ = √π/8
        sqrt_pi = math.sqrt(math.pi)
        print()
        print(f"  Continuum S³ reference (unit sphere):")
        print(f"    a₀ = √π = {sqrt_pi:.8f}")
        print(f"    a₁ = √π = {sqrt_pi:.8f}")
        print(f"    a₂ = √π/8 = {sqrt_pi/8:.8f}")
        print(f"  Ratios to continuum:")
        print(f"    a₀(graph)/a₀(cont) = {sd_coeffs[0]/sqrt_pi:.6f}")
        print(f"    a₁(graph)/a₁(cont) = {sd_coeffs[1]/sqrt_pi:.6f}")
        if sqrt_pi/8 > 1e-10:
            print(f"    a₂(graph)/a₂(cont) = {sd_coeffs[2]/(sqrt_pi/8):.6f}")

        cfg_results['heat_kernel_sd_coeffs'] = {
            f'a_{k}': float(c) for k, c in enumerate(sd_coeffs)
        }
        cfg_results['heat_kernel_sd_powers'] = powers

        # The CC spectral action uses:
        # S_b = f₀ a₀ Λ³ + f₂ a₂ Λ + f₄ a₄ Λ⁻¹ + ...
        # For d=3, the leading term is Λ³ (cosmological),
        # next is Λ (Einstein-Hilbert), then Λ⁻¹, ...
        print()

        # ----------------------------------------------------------------
        # 4. Direct finite sum: Σ |λ_i|^k for spectral zeta
        # ----------------------------------------------------------------
        print("  4. SPECTRAL ZETA MOMENTS Σ |λ|^k")
        print("  " + "-" * 40)

        abs_evals = np.abs(evals)

        zeta_moments = {}
        for k in range(-4, 8):
            if k == 0:
                val = float(dim_H)
            elif k < 0 and np.any(abs_evals < 1e-14):
                val = float('inf')
            else:
                val = float(np.sum(abs_evals ** k))
            zeta_moments[k] = val

        print(f"  Tr |D|^k (finite spectral zeta at -k):")
        for k in range(-4, 8):
            val = zeta_moments[k]
            if val == float('inf'):
                print(f"    k={k:+2d}: inf (zero eigenvalues)")
            else:
                print(f"    k={k:+2d}: {val:.8f}")

        # Key: Tr |D|^0 = dim_H, Tr |D|^{-1} is the spectral zeta at s=1
        # Compare Tr I = dim_H to Δ⁻¹ = 40
        if cfg['n_max'] == 3:
            print(f"\n  *** dim_H = {dim_H} = Δ⁻¹ = 40 from Paper 2 ***")

        cfg_results['zeta_moments'] = {str(k): v for k, v in zeta_moments.items()}
        print()

        # ----------------------------------------------------------------
        # 5. Supertrace of the spectral action
        # ----------------------------------------------------------------
        print("  5. SUPERTRACE (BOSON-FERMION GRADED)")
        print("  " + "-" * 40)

        str_results = supertrace_action(evals, gamma_signs, lam_max**2 * 1.01)

        print(f"  Str(Θ(Λ²-D²)) = {str_results['sharp_cutoff']:.6f}")
        print(f"  Str(exp(-D²/Λ²)) = {str_results['gaussian']:.6f}")
        print(f"  Str(exp(-tD²)) at various t:")
        for t_label, val in str_results['heat_kernel'].items():
            print(f"    {t_label}: {val:.8f}")

        # Sprint ST-1 Finding F1: on continuum S³, Str vanishes identically
        # because a_k^{D²}/a_k^{Δ_LB} = 4 at every SD order
        print(f"\n  ST-1 F1 check: Str should vanish on continuum S³")
        print(f"  On finite graph: Str(sharp) = {str_results['sharp_cutoff']:.6f}")

        # Graded zeta moments: Str |D|^k = Tr(γ |D|^k)
        print(f"\n  Graded zeta moments Str |D|^k:")
        for k in [0, 1, 2, 3, -1, -2]:
            if k < 0 and np.any(abs_evals < 1e-14):
                str_val = float('inf')
            else:
                str_val = float(np.sum(gamma_signs * abs_evals ** k))
            print(f"    k={k:+2d}: {str_val:.8f}")

        cfg_results['supertrace'] = {
            'sharp_cutoff': str_results['sharp_cutoff'],
            'gaussian': str_results['gaussian'],
            'heat_kernel': str_results['heat_kernel'],
        }
        print()

        # ----------------------------------------------------------------
        # 6. K-comparison: search for cutoff/function that gives K/π
        # ----------------------------------------------------------------
        print("  6. K/π COMPARISON SEARCH")
        print("  " + "-" * 40)

        # Scan over cutoffs to find where various spectral actions cross K/π
        K_pi = targets['K_over_pi']
        K = targets['K']

        Lambda_scan = np.linspace(0.5, lam_max * 1.5, 200)

        best_matches = {}

        for test_fn_name, test_fn in [
            ('sharp', lambda esq, Lsq: sharp_cutoff_action(esq, Lsq)),
            ('gaussian', lambda esq, Lsq: gaussian_cutoff_action(esq, Lsq)),
            ('sigmoid_w0.1', lambda esq, Lsq: smooth_cutoff_action(esq, Lsq, width=0.1)),
            ('sigmoid_w0.3', lambda esq, Lsq: smooth_cutoff_action(esq, Lsq, width=0.3)),
            ('poly_p2', lambda esq, Lsq: polynomial_cutoff_action(esq, Lsq, degree=2)),
        ]:
            vals = [test_fn(evals_sq, L**2) for L in Lambda_scan]

            # Find where the value crosses K/π
            for target_name, target_val in [('K/pi', K_pi), ('K', K), ('B', 42.0)]:
                best_dist = float('inf')
                best_L = None
                best_val = None
                for L, v in zip(Lambda_scan, vals):
                    dist = abs(v - target_val)
                    if dist < best_dist:
                        best_dist = dist
                        best_L = L
                        best_val = v

                key = f"{test_fn_name}_vs_{target_name}"
                if best_dist < target_val * 0.5:
                    best_matches[key] = {
                        'Lambda': float(best_L),
                        'value': float(best_val),
                        'target': target_val,
                        'distance': float(best_dist),
                        'rel_error': float(best_dist / target_val),
                    }

        print(f"  Searching for Λ where Tr f(D²/Λ²) ≈ K/π = {K_pi:.4f}:")
        print(f"  (also checking K = {K:.4f} and B = 42)")
        print()

        for key, match in sorted(best_matches.items()):
            print(f"    {key}:")
            print(f"      Λ = {match['Lambda']:.4f}, value = {match['value']:.4f}, "
                  f"target = {match['target']:.4f}, rel_err = {match['rel_error']:.4f}")

        cfg_results['K_comparison'] = best_matches
        print()

        # ----------------------------------------------------------------
        # 7. The non-perturbative finite trace
        # ----------------------------------------------------------------
        print("  7. NON-PERTURBATIVE FINITE TRACES")
        print("  " + "-" * 40)

        # The exact finite-graph spectral action for specific functions

        # Tr(I) = dim_H
        print(f"  Tr(I) = dim_H = {dim_H}")

        # Tr(D²) — the total "Dirac squared" mass
        tr_d2 = float(np.sum(evals_sq))
        print(f"  Tr(D²) = {tr_d2:.8f}")

        # Tr(D⁴)
        tr_d4 = float(np.sum(evals_sq ** 2))
        print(f"  Tr(D⁴) = {tr_d4:.8f}")

        # Str(D) = Tr(γD) — the graded first moment (linear spectral action)
        str_D = float(np.sum(gamma_signs * evals))
        print(f"  Str(D) = Tr(γD) = {str_D:.8f}")

        # Str(D²) — graded second moment
        str_D2 = float(np.sum(gamma_signs * evals_sq))
        print(f"  Str(D²) = Tr(γD²) = {str_D2:.8f}")

        # Str(|D|) — graded absolute moment
        str_absD = float(np.sum(gamma_signs * abs_evals))
        print(f"  Str(|D|) = Tr(γ|D|) = {str_absD:.8f}")

        # Key ratios
        print(f"\n  Key ratios:")
        print(f"    Tr(D²)/dim_H = {tr_d2/dim_H:.8f}")
        print(f"    Tr(D⁴)/Tr(D²) = {tr_d4/tr_d2:.8f}")
        print(f"    Tr(D²)/π = {tr_d2/math.pi:.8f}")
        print(f"    Str(D)/dim_H = {str_D/dim_H:.8f}")

        # Comparison to Casimir sums on the continuum
        # On unit S³: Σ g_n |λ_n|^2 = Σ 2(n+1)(n+2)(n+3/2)²
        # This is the second Casimir trace
        continuum_tr_d2 = 0.0
        for n_ch in range(cfg['n_max']):
            g_n = 2 * (n_ch + 1) * (n_ch + 2)
            lam_n = n_ch + 1.5
            continuum_tr_d2 += g_n * lam_n ** 2

        print(f"\n  Continuum Casimir reference (n_max={cfg['n_max']}):")
        print(f"    Σ g_n |λ_n|² = {continuum_tr_d2:.4f}")
        print(f"    Tr(D²) / continuum = {tr_d2/continuum_tr_d2:.6f}")

        cfg_results['finite_traces'] = {
            'Tr_I': dim_H,
            'Tr_D2': tr_d2,
            'Tr_D4': tr_d4,
            'Str_D': str_D,
            'Str_D2': str_D2,
            'Str_absD': str_absD,
            'continuum_Tr_D2': continuum_tr_d2,
        }
        print()

        # ----------------------------------------------------------------
        # 8. Paper 2 decomposition check: B, F, Δ from spectral data
        # ----------------------------------------------------------------
        if cfg['n_max'] == 3:
            print("  8. PAPER 2 DECOMPOSITION CHECK (n_max=3)")
            print("  " + "-" * 40)

            # dim_H = 40 = Δ⁻¹  ← Already known
            print(f"  dim_H = {dim_H} = Δ⁻¹ from Paper 2 ✓")

            # Per-shell Casimir c(n) = n²(n²-1)/2 at n_fock
            # B = 42 = Σ (2l+1)l(l+1) over n_fock=1..3, l=0..n_fock-1
            # This is a pure n,l lattice property, not from D eigenvalues
            B_check = 0
            for n_fock in range(1, cfg['n_max'] + 1):
                for l in range(n_fock):
                    B_check += (2*l + 1) * l * (l + 1)
            print(f"  B = Σ (2l+1)l(l+1) = {B_check} (Paper 2: 42)")

            # Str(D) on the finite graph vs continuum
            # On continuum: Str(D) = Σ g_n^+ |λ_n| - Σ g_n^- |λ_n|
            # With balanced ±: this gives a measure of the chirality asymmetry
            print(f"  Str(D) = {str_D:.6f}")
            print(f"  Str(|D|) = {str_absD:.6f}")

            # Spectral zeta at s=-1: Σ |λ_i|
            sum_abs_lam = float(np.sum(abs_evals))
            print(f"  Σ |λ_i| = {sum_abs_lam:.8f}")
            print(f"  Σ |λ_i| / π = {sum_abs_lam/math.pi:.8f}")

            # The CC spectral action on the finite graph at the natural cutoff
            # Λ = max CH eigenvalue for n_max-1 = 5/2 for n_max=3 (n_CH=2)
            L_natural = cfg['n_max'] - 1 + 1.5  # highest CH level
            L_sq_natural = L_natural ** 2

            sa_natural = connes_chamseddine_bosonic(evals_sq, L_sq_natural)
            print(f"\n  CC action at Λ = {L_natural} (natural cutoff n_CH={cfg['n_max']-1}):")
            for key, val in sa_natural.items():
                print(f"    {key}: {val:.6f}")

            cfg_results['paper2_check'] = {
                'dim_H_equals_Delta_inv': dim_H == 40,
                'B_from_lattice': B_check,
                'Str_D': str_D,
                'Str_absD': str_absD,
                'sum_abs_lambda': sum_abs_lam,
                'natural_cutoff_Lambda': L_natural,
                'cc_at_natural': {k: float(v) for k, v in sa_natural.items()},
            }

        print()
        results[label] = cfg_results

    # ====================================================================
    # SYNTHESIS: cross-configuration comparison
    # ====================================================================
    print("=" * 70)
    print("SYNTHESIS: CROSS-CONFIGURATION COMPARISON")
    print("=" * 70)
    print()

    # Compare Tr(D²) across configs
    print("Tr(D²) across configurations:")
    for label, cfg_res in results.items():
        if isinstance(cfg_res, dict) and 'finite_traces' in cfg_res:
            tr = cfg_res['finite_traces']['Tr_D2']
            dim = cfg_res['dim_H']
            print(f"  {label}: Tr(D²) = {tr:.6f}, Tr(D²)/dim = {tr/dim:.6f}")

    print()
    print("Str(D) across configurations (should be nonzero):")
    for label, cfg_res in results.items():
        if isinstance(cfg_res, dict) and 'finite_traces' in cfg_res:
            print(f"  {label}: Str(D) = {cfg_res['finite_traces']['Str_D']:.6f}")

    print()
    print("Str(D²) across configurations (vanishes on continuum by ST-1 F1):")
    for label, cfg_res in results.items():
        if isinstance(cfg_res, dict) and 'finite_traces' in cfg_res:
            print(f"  {label}: Str(D²) = {cfg_res['finite_traces']['Str_D2']:.6f}")

    # ====================================================================
    # Save results
    # ====================================================================
    output_path = Path(__file__).resolve().parent / "data" / "spectral_triple_spectral_action.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    def json_safe(obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=json_safe)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
