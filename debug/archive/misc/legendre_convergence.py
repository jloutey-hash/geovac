#!/usr/bin/env python3
"""
Legendre Convergence Analysis for Nuclear Potential Expansion.

Measures how well the truncated Legendre expansion V_approx(θ; k_max)
approximates the exact nuclear potential V_exact(θ), as a function of
truncation order, charge asymmetry, and origin choice.

Key question: is convergence fast enough that l_max=6-8 suffices for
asymmetric bonds, or do we need a fundamentally different angular basis?
"""

import numpy as np
from numpy.polynomial.legendre import legval
from scipy.special import legendre
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pathlib import Path
import json

PLOT_DIR = Path(__file__).parent / "plots"
DATA_DIR = Path(__file__).parent / "data"
PLOT_DIR.mkdir(exist_ok=True)
DATA_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def f_k(s: float, rho: float, k: int) -> float:
    """f_k(s, rho) = (min/max)^k / max.  Split-region Legendre coefficient."""
    if s <= 0.0 or rho <= 0.0:
        return 0.0
    mn, mx = min(s, rho), max(s, rho)
    return (mn / mx)**k / mx


def V_exact(theta: np.ndarray, s: float,
            Z_A: float, rho_A: float,
            Z_B: float, rho_B: float) -> np.ndarray:
    """
    Exact nuclear potential for a single electron at (s, θ).
    Nucleus A at +z (distance rho_A), nucleus B at -z (distance rho_B).

    V = -Z_A / |r - rho_A z| - Z_B / |r + rho_B z|
      = -Z_A / sqrt(s² + ρ_A² - 2sρ_A cos θ)
        -Z_B / sqrt(s² + ρ_B² + 2sρ_B cos θ)
    """
    cos_th = np.cos(theta)
    r_A = np.sqrt(s**2 + rho_A**2 - 2*s*rho_A*cos_th)
    r_B = np.sqrt(s**2 + rho_B**2 + 2*s*rho_B*cos_th)
    # Avoid division by zero
    r_A = np.maximum(r_A, 1e-15)
    r_B = np.maximum(r_B, 1e-15)
    return -Z_A / r_A - Z_B / r_B


def legendre_coefficients(s: float, Z_A: float, rho_A: float,
                          Z_B: float, rho_B: float,
                          k_max: int) -> np.ndarray:
    """
    Compute Legendre expansion coefficients c_k for the nuclear potential.

    The generating function expansion of 1/|r - r'| gives:
      1/sqrt(s² + ρ² - 2sρ cos θ) = Σ_k f_k(s,ρ) P_k(cos θ)  [nucleus at +z]
      1/sqrt(s² + ρ² + 2sρ cos θ) = Σ_k (-1)^k f_k(s,ρ) P_k(cos θ)  [nucleus at -z]

    So: c_k = -Z_A * f_k(s, ρ_A) - Z_B * (-1)^k * f_k(s, ρ_B)
    """
    c = np.zeros(k_max + 1)
    for k in range(k_max + 1):
        c[k] = -Z_A * f_k(s, rho_A, k) - Z_B * ((-1)**k) * f_k(s, rho_B, k)
    return c


def V_approx(theta: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    """Evaluate truncated Legendre expansion at given angles."""
    cos_th = np.cos(theta)
    result = np.zeros_like(theta)
    for k, c_k in enumerate(coeffs):
        result += c_k * legendre(k)(cos_th)
    return result


def L2_relative_error(s: float, Z_A: float, rho_A: float,
                      Z_B: float, rho_B: float,
                      k_max: int, n_quad: int = 200) -> float:
    """
    Compute L² relative error of truncated expansion:
    ||V_exact - V_approx||² / ||V_exact||²
    integrated over θ ∈ [0, π] with sin θ dθ weight.
    Uses Gauss-Legendre quadrature on x = cos θ.
    """
    # Gauss-Legendre on [-1, 1] (x = cos θ, dx = -sin θ dθ)
    x_nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    theta_nodes = np.arccos(x_nodes)

    v_ex = V_exact(theta_nodes, s, Z_A, rho_A, Z_B, rho_B)
    coeffs = legendre_coefficients(s, Z_A, rho_A, Z_B, rho_B, k_max)
    v_ap = V_approx(theta_nodes, coeffs)

    diff_sq = np.sum(weights * (v_ex - v_ap)**2)
    norm_sq = np.sum(weights * v_ex**2)

    if norm_sq < 1e-30:
        return 0.0
    return diff_sq / norm_sq


# ---------------------------------------------------------------------------
# Verification: check that expansion matches exact at high k_max
# ---------------------------------------------------------------------------

def verify_expansion():
    """Verify that the Legendre expansion converges to V_exact."""
    print("=" * 60)
    print("VERIFICATION: Legendre expansion vs exact potential")
    print("=" * 60)

    # Simple test: symmetric H2, s=1.5, rho=0.7
    s, rho = 1.5, 0.7
    theta_test = np.array([0.0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])

    v_ex = V_exact(theta_test, s, 1.0, rho, 1.0, rho)

    for k_max in [2, 5, 10, 20, 40]:
        coeffs = legendre_coefficients(s, 1.0, rho, 1.0, rho, k_max)
        v_ap = V_approx(theta_test, coeffs)
        max_err = np.max(np.abs(v_ex - v_ap))
        print(f"  k_max={k_max:3d}: max |V_ex - V_ap| = {max_err:.2e}")

    # Also test asymmetric case
    print("\n  Asymmetric 6:1 at s=0.5:")
    s = 0.5
    v_ex = V_exact(theta_test, s, 6.0, rho, 1.0, rho)
    for k_max in [2, 5, 10, 20, 40]:
        coeffs = legendre_coefficients(s, 6.0, rho, 1.0, rho, k_max)
        v_ap = V_approx(theta_test, coeffs)
        max_err = np.max(np.abs(v_ex - v_ap))
        print(f"  k_max={k_max:3d}: max |V_ex - V_ap| = {max_err:.2e}")

    print()


# ---------------------------------------------------------------------------
# Test configurations
# ---------------------------------------------------------------------------

def build_configurations(rho: float = 0.7):
    """Build all test configurations."""
    configs = []

    for s in [0.5, 1.5]:
        # (a) Symmetric H2, midpoint
        configs.append({
            'label': f'H2 mid (s={s})',
            'Z_A': 1.0, 'Z_B': 1.0,
            'rho_A': rho, 'rho_B': rho,
            's': s, 'system': 'H2', 'origin': 'midpoint',
        })

        # (b) HeH+ midpoint
        configs.append({
            'label': f'HeH+ mid (s={s})',
            'Z_A': 2.0, 'Z_B': 1.0,
            'rho_A': rho, 'rho_B': rho,
            's': s, 'system': 'HeH+', 'origin': 'midpoint',
        })

        # (b) HeH+ charge-center
        Z_A, Z_B = 2.0, 1.0
        rho_A_cc = rho * 2*Z_B / (Z_A + Z_B)
        rho_B_cc = rho * 2*Z_A / (Z_A + Z_B)
        configs.append({
            'label': f'HeH+ CC (s={s})',
            'Z_A': Z_A, 'Z_B': Z_B,
            'rho_A': rho_A_cc, 'rho_B': rho_B_cc,
            's': s, 'system': 'HeH+', 'origin': 'charge-center',
        })

        # (c) 6:1 midpoint
        configs.append({
            'label': f'6:1 mid (s={s})',
            'Z_A': 6.0, 'Z_B': 1.0,
            'rho_A': rho, 'rho_B': rho,
            's': s, 'system': '6:1', 'origin': 'midpoint',
        })

        # (c) 6:1 charge-center
        Z_A, Z_B = 6.0, 1.0
        rho_A_cc = rho * 2*Z_B / (Z_A + Z_B)
        rho_B_cc = rho * 2*Z_A / (Z_A + Z_B)
        configs.append({
            'label': f'6:1 CC (s={s})',
            'Z_A': Z_A, 'Z_B': Z_B,
            'rho_A': rho_A_cc, 'rho_B': rho_B_cc,
            's': s, 'system': '6:1', 'origin': 'charge-center',
        })

    return configs


# ---------------------------------------------------------------------------
# Convergence sweep
# ---------------------------------------------------------------------------

def convergence_sweep(configs, k_max_range=None):
    """Compute L² error vs k_max for all configurations."""
    if k_max_range is None:
        k_max_range = np.arange(0, 21)

    results = {}
    for cfg in configs:
        label = cfg['label']
        errors = []
        for k_max in k_max_range:
            err = L2_relative_error(
                cfg['s'], cfg['Z_A'], cfg['rho_A'],
                cfg['Z_B'], cfg['rho_B'], int(k_max)
            )
            errors.append(err)
        results[label] = np.array(errors)
        print(f"  {label}: err[k=0]={errors[0]:.4f}, err[k=4]={errors[4]:.6f}, "
              f"err[k=10]={errors[10]:.2e}, err[k=20]={errors[20]:.2e}")

    return results


# ---------------------------------------------------------------------------
# Fit convergence rate
# ---------------------------------------------------------------------------

def fit_convergence(k_max_range, errors, label=""):
    """Fit algebraic and exponential decay models to the error vs k_max."""
    # Use k >= 2 to avoid initial transients
    mask = (k_max_range >= 2) & (errors > 1e-30)
    k_fit = k_max_range[mask].astype(float)
    e_fit = errors[mask]

    if len(k_fit) < 3:
        return {'type': 'insufficient', 'p': np.nan, 'beta': np.nan,
                'R2_alg': np.nan, 'R2_exp': np.nan}

    # Algebraic: log(error) = log(A) - p * log(k)
    try:
        log_k = np.log(k_fit)
        log_e = np.log(e_fit)
        p_alg = np.polyfit(log_k, log_e, 1)
        pred_alg = np.polyval(p_alg, log_k)
        ss_res = np.sum((log_e - pred_alg)**2)
        ss_tot = np.sum((log_e - np.mean(log_e))**2)
        R2_alg = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        p = -p_alg[0]
    except Exception:
        R2_alg, p = 0.0, np.nan

    # Exponential: log(error) = log(A) - beta * k
    try:
        p_exp = np.polyfit(k_fit, log_e, 1)
        pred_exp = np.polyval(p_exp, k_fit)
        ss_res = np.sum((log_e - pred_exp)**2)
        ss_tot = np.sum((log_e - np.mean(log_e))**2)
        R2_exp = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        beta = -p_exp[0]
    except Exception:
        R2_exp, beta = 0.0, np.nan

    conv_type = 'exponential' if R2_exp > R2_alg else 'algebraic'

    return {
        'type': conv_type,
        'p': p, 'beta': beta,
        'R2_alg': R2_alg, 'R2_exp': R2_exp,
    }


def find_k_for_error(errors, k_max_range, target):
    """Find minimum k_max achieving target relative L² error."""
    for i, e in enumerate(errors):
        if e <= target:
            return int(k_max_range[i])
    return None  # not achieved in range


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_A_potential_profiles(configs, k_maxes=[2, 4, 8, 16]):
    """Plot A: V_exact vs V_approx for 6:1 case at s=0.5."""
    theta = np.linspace(0.001, np.pi - 0.001, 500)

    # Find 6:1 configs at s=0.5
    cfg_mid = [c for c in configs if c['system'] == '6:1'
               and c['origin'] == 'midpoint' and c['s'] == 0.5][0]
    cfg_cc = [c for c in configs if c['system'] == '6:1'
              and c['origin'] == 'charge-center' and c['s'] == 0.5][0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, cfg, title in zip(axes, [cfg_mid, cfg_cc],
                               ['Midpoint origin', 'Charge-center origin']):
        v_ex = V_exact(theta, cfg['s'], cfg['Z_A'], cfg['rho_A'],
                       cfg['Z_B'], cfg['rho_B'])
        ax.plot(np.degrees(theta), v_ex, 'k-', lw=2, label='Exact')

        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(k_maxes)))
        for km, color in zip(k_maxes, colors):
            coeffs = legendre_coefficients(
                cfg['s'], cfg['Z_A'], cfg['rho_A'],
                cfg['Z_B'], cfg['rho_B'], km
            )
            v_ap = V_approx(theta, coeffs)
            ax.plot(np.degrees(theta), v_ap, '--', color=color, lw=1.2,
                    label=f'k_max={km}')

        ax.set_xlabel('θ (degrees)')
        ax.set_ylabel('V_nuc')
        ax.set_title(f'6:1 asymmetry, s=0.5 — {title}')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Nuclear Potential Legendre Expansion: 6:1 Charge Asymmetry',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'legendre_potential_profiles.png', dpi=150)
    plt.close()
    print(f"  Saved: {PLOT_DIR / 'legendre_potential_profiles.png'}")


def plot_B_loglog(k_max_range, results):
    """Plot B: L² relative error vs k_max (log-log)."""
    fig, ax = plt.subplots(figsize=(10, 7))

    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<']
    for i, (label, errors) in enumerate(results.items()):
        mask = errors > 1e-30
        ax.plot(k_max_range[mask], errors[mask],
                marker=markers[i % len(markers)], markersize=4,
                label=label, linewidth=1.2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('k_max')
    ax.set_ylabel('L² relative error')
    ax.set_title('Legendre Convergence (log-log)')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    ax.axhline(0.01, color='r', ls=':', alpha=0.5, label='1% error')
    ax.axhline(0.001, color='orange', ls=':', alpha=0.5, label='0.1% error')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'legendre_convergence_loglog.png', dpi=150)
    plt.close()
    print(f"  Saved: {PLOT_DIR / 'legendre_convergence_loglog.png'}")


def plot_C_loglinear(k_max_range, results):
    """Plot C: L² relative error vs k_max (log-linear)."""
    fig, ax = plt.subplots(figsize=(10, 7))

    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<']
    for i, (label, errors) in enumerate(results.items()):
        mask = errors > 1e-30
        ax.plot(k_max_range[mask], errors[mask],
                marker=markers[i % len(markers)], markersize=4,
                label=label, linewidth=1.2)

    ax.set_yscale('log')
    ax.set_xlabel('k_max')
    ax.set_ylabel('L² relative error')
    ax.set_title('Legendre Convergence (log-linear)')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    ax.axhline(0.01, color='r', ls=':', alpha=0.5)
    ax.axhline(0.001, color='orange', ls=':', alpha=0.5)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'legendre_convergence_loglinear.png', dpi=150)
    plt.close()
    print(f"  Saved: {PLOT_DIR / 'legendre_convergence_loglinear.png'}")


def plot_D_coefficients(configs, k_max: int = 20):
    """Plot D: |c_k| vs k for each configuration."""
    fig, ax = plt.subplots(figsize=(10, 7))

    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<']
    ks = np.arange(0, k_max + 1)

    for i, cfg in enumerate(configs):
        coeffs = legendre_coefficients(
            cfg['s'], cfg['Z_A'], cfg['rho_A'],
            cfg['Z_B'], cfg['rho_B'], k_max
        )
        abs_c = np.abs(coeffs)
        mask = abs_c > 1e-30
        ax.plot(ks[mask], abs_c[mask],
                marker=markers[i % len(markers)], markersize=4,
                label=cfg['label'], linewidth=1.0)

    ax.set_yscale('log')
    ax.set_xlabel('k')
    ax.set_ylabel('|c_k|')
    ax.set_title('Legendre Coefficients |c_k| vs k')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'legendre_coefficients.png', dpi=150)
    plt.close()
    print(f"  Saved: {PLOT_DIR / 'legendre_coefficients.png'}")


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary(k_max_range, results, configs):
    """Print the summary table."""
    print("\n" + "=" * 100)
    print("LEGENDRE CONVERGENCE SUMMARY")
    print("=" * 100)
    print(f"{'Configuration':<22} {'Type':<12} {'Rate':<12} "
          f"{'R2_alg':<8} {'R2_exp':<8} {'k(1%)':<7} {'k(0.1%)':<8} {'CC factor':<10}")
    print("-" * 100)

    summary_data = {}

    for cfg in configs:
        label = cfg['label']
        errors = results[label]

        fit = fit_convergence(k_max_range, errors, label)
        k_1pct = find_k_for_error(errors, k_max_range, 0.01)
        k_01pct = find_k_for_error(errors, k_max_range, 0.001)

        rate_str = (f"p={fit['p']:.2f}" if fit['type'] == 'algebraic'
                    else f"b={fit['beta']:.2f}")
        k1_str = str(k_1pct) if k_1pct is not None else '>20'
        k01_str = str(k_01pct) if k_01pct is not None else '>20'

        # Compute CC improvement factor for systems that have both origins
        cc_factor = ''
        if cfg['origin'] == 'charge-center':
            mid_label = label.replace(' CC ', ' mid ')
            if mid_label in results:
                mid_errors = results[mid_label]
                k_mid = find_k_for_error(mid_errors, k_max_range, 0.01)
                if k_mid is not None and k_1pct is not None and k_1pct > 0:
                    cc_factor = f"{k_mid/k_1pct:.2f}x"

        print(f"{label:<22} {fit['type']:<12} {rate_str:<12} "
              f"{fit['R2_alg']:<8.4f} {fit['R2_exp']:<8.4f} "
              f"{k1_str:<7} {k01_str:<8} {cc_factor:<10}")

        summary_data[label] = {
            'convergence_type': fit['type'],
            'algebraic_exponent': float(fit['p']),
            'exponential_rate': float(fit['beta']),
            'R2_algebraic': float(fit['R2_alg']),
            'R2_exponential': float(fit['R2_exp']),
            'k_max_1pct': k_1pct,
            'k_max_01pct': k_01pct,
        }

    print("=" * 100)

    # Save data
    output = {
        'rho': 0.7,
        'k_max_range': k_max_range.tolist(),
        'results': {k: v.tolist() for k, v in results.items()},
        'summary': summary_data,
    }
    out_path = DATA_DIR / 'legendre_convergence.json'
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nData saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Nuclear Potential Legendre Convergence Analysis")
    print("=" * 60)

    # Step 1: Verify expansion
    verify_expansion()

    # Step 2: Build configurations
    configs = build_configurations(rho=0.7)
    print(f"\n{len(configs)} test configurations built.\n")

    # Step 3: Convergence sweep
    print("Computing L² errors...")
    k_max_range = np.arange(0, 21)
    results = convergence_sweep(configs, k_max_range)

    # Step 4: Plots
    print("\nGenerating plots...")
    plot_A_potential_profiles(configs)
    plot_B_loglog(k_max_range, results)
    plot_C_loglinear(k_max_range, results)
    plot_D_coefficients(configs)

    # Step 5: Summary
    print_summary(k_max_range, results, configs)

    # Key conclusion
    print("\n" + "=" * 60)
    print("KEY FINDING:")
    print("=" * 60)

    # Check 6:1 midpoint vs charge-center at s=0.5
    err_mid = results.get('6:1 mid (s=0.5)', np.array([]))
    err_cc = results.get('6:1 CC (s=0.5)', np.array([]))

    if len(err_mid) > 0 and len(err_cc) > 0:
        k1_mid = find_k_for_error(err_mid, k_max_range, 0.01)
        k1_cc = find_k_for_error(err_cc, k_max_range, 0.01)

        fit_mid = fit_convergence(k_max_range, err_mid)
        fit_cc = fit_convergence(k_max_range, err_cc)

        print(f"  6:1 midpoint (s=0.5):  {fit_mid['type']}, k(1%)={k1_mid}")
        print(f"  6:1 charge-center:     {fit_cc['type']}, k(1%)={k1_cc}")

        if k1_cc is not None and k1_cc <= 8:
            print("  => Charge-center origin makes l_max=6-8 SUFFICIENT.")
            print("  => No fundamentally different angular basis needed.")
        elif k1_mid is not None and k1_mid <= 8:
            print("  => Even midpoint origin converges by l_max=6-8.")
        else:
            print("  => Convergence is SLOW. May need a different angular basis.")


if __name__ == '__main__':
    main()
