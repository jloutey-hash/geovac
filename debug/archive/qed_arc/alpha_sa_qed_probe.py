"""
alpha-SA QED Probe: QED-motivated spectral quantities on S^3.

The previous probe (alpha_sa_hopf_action.py) tested graph-topological invariants
(L0, L1 spectra, Hodge decomposition, etc.) and found a CLEAN NEGATIVE: no
graph-spectral invariant produces K(m)/pi.

This probe tests the other side: QED-on-S^3 calculations that weight the Dirac
spectrum by interaction physics. The Dirac modes are the "catalog"; QED tells
you how to weight and sum them when they interact.

Key QED quantities:
  - One-loop effective action: Gamma = -(1/2) Tr log(D^2/Lambda^2)
  - Vacuum polarization: Pi(0) ~ Sigma g_n / |lambda_n|^2
  - Photon self-energy with vertex selection rules
  - Sunset diagrams (two fermion lines + one photon line)
  - Spectral determinant ratios
  - Renormalization-group running

The target at each m is K(m)/pi = B(m) + F - Delta(m).
"""

import numpy as np
from scipy.integrate import quad
import json


ALPHA_INV = 137.035999084
K_OVER_PI = ALPHA_INV / np.pi
N_INF = 500  # truncation for "infinite" sums


# =====================================================================
# Spectral data
# =====================================================================

def dirac_ev(n):
    """Absolute Dirac eigenvalue on unit S^3: |lambda_n| = n + 3/2."""
    return n + 1.5

def dirac_deg(n):
    """Single-chirality Dirac degeneracy: g_n = 2(n+1)(n+2)."""
    return 2 * (n + 1) * (n + 2)

def scalar_ev(n):
    """Scalar Laplacian eigenvalue on unit S^3: n^2 - 1, n >= 1."""
    return n * n - 1

def scalar_deg(n):
    """Scalar degeneracy: n^2."""
    return n * n

def K_ingredients(m):
    B = m * (m - 1) * (m + 1) * (m + 2) * (2 * m + 1) / 20
    F = np.pi**2 / 6
    Delta = 1.0 / (2 * (m + 1) * (m + 2))
    return {'B': B, 'F': F, 'Delta': Delta, 'K_over_pi': B + F - Delta}


# =====================================================================
# QED-motivated coupling / vertex weights
# =====================================================================

def cg_coupling_sq(n):
    """
    Squared reduced CG coupling for E1 transitions from level n.
    T(n)^2 = 1/(2n+3)^2 -- from Paper 28 / qed_vertex.py.
    """
    return 1.0 / (2 * n + 3)**2


def vertex_allowed(n1, n2):
    """E1 dipole selection rule: |n1 - n2| = 1."""
    return abs(n1 - n2) == 1


def vertex_parity_allowed(n1, n2, n_gamma):
    """Paper 28 vertex parity: n1 + n2 + n_gamma must be odd."""
    return (n1 + n2 + n_gamma) % 2 == 1


# =====================================================================
# Probe
# =====================================================================

def compute_qed_candidates(m):
    K_data = K_ingredients(m)
    B_val = K_data['B']
    Delta_val = K_data['Delta']
    F_val = K_data['F']
    target = K_data['K_over_pi']

    c = {}  # candidates

    # Reference
    c['B + F - Delta (target)'] = B_val + F_val - Delta_val

    # -----------------------------------------------------------------
    # 1. TRUNCATED DIRAC SPECTRAL ZETAS
    # -----------------------------------------------------------------
    for s in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
        c[f'D_trunc(s={s})'] = sum(
            dirac_deg(n) * dirac_ev(n)**(-s) for n in range(m))

    # Infinite versions (well-converged)
    for s in [2.0, 3.0, 4.0]:
        c[f'D_inf(s={s})'] = sum(
            dirac_deg(n) * dirac_ev(n)**(-s) for n in range(N_INF))

    # -----------------------------------------------------------------
    # 2. TRUNCATED SPECTRAL DETERMINANTS
    # -----------------------------------------------------------------
    if m > 0:
        log_det_D2 = sum(dirac_deg(n) * np.log(dirac_ev(n)**2)
                         for n in range(m))
        c['log det(D^2)_trunc'] = log_det_D2
        c['(1/2) log det(D^2)'] = 0.5 * log_det_D2

        log_det_D = sum(dirac_deg(n) * np.log(dirac_ev(n))
                        for n in range(m))
        c['log det|D|_trunc'] = log_det_D

    # Scalar determinant (skip zero mode at n=1)
    if m >= 2:
        log_det_S = sum(scalar_deg(n) * np.log(scalar_ev(n))
                        for n in range(2, m + 1))
        c['log det(scalar)_trunc'] = log_det_S

    # -----------------------------------------------------------------
    # 3. DIRAC CASIMIR ENERGY (truncated mode sums)
    # -----------------------------------------------------------------
    E_D = sum(dirac_deg(n) * dirac_ev(n) for n in range(m))
    c['E_cas_Dirac'] = E_D
    c['E_cas_Dirac / 4'] = E_D / 4.0
    c['E_cas_Dirac / (4*pi)'] = E_D / (4 * np.pi)

    # Scalar Casimir
    E_S = sum(scalar_deg(n) * np.sqrt(max(scalar_ev(n), 0))
              for n in range(1, m + 1))
    c['E_cas_scalar'] = E_S

    # Supertrace Casimir (bosonic - fermionic/4)
    c['Str_Casimir'] = E_S - E_D / 4.0

    # -----------------------------------------------------------------
    # 4. ONE-LOOP EFFECTIVE ACTION  Gamma = -(1/2) Tr log(D^2/Lambda^2)
    # -----------------------------------------------------------------
    for lam_sq_label, lam_sq in [
        ('1', 1.0),
        ('|lam_m|^2', dirac_ev(m)**2),
        ('|lam_{m-1}|^2', dirac_ev(m - 1)**2 if m >= 1 else 1),
        ('4', 4.0),
        ('(3/2)^2', 2.25),
    ]:
        Gamma = -0.5 * sum(dirac_deg(n) * np.log(dirac_ev(n)**2 / lam_sq)
                           for n in range(m))
        c[f'Gamma_D(Lam2={lam_sq_label})'] = Gamma
        c[f'-Gamma_D(Lam2={lam_sq_label})'] = -Gamma

    # -----------------------------------------------------------------
    # 5. VACUUM POLARIZATION Pi(q=0)
    # -----------------------------------------------------------------
    # One-loop bubble: Pi ~ sum g_n / |lambda_n|^2
    Pi_trunc = sum(dirac_deg(n) / dirac_ev(n)**2 for n in range(m))
    Pi_inf = sum(dirac_deg(n) / dirac_ev(n)**2 for n in range(N_INF))
    c['Pi(0)_trunc'] = Pi_trunc
    c['Pi(0)_inf'] = Pi_inf

    # With (2/3pi) QED prefactor
    c['(2/3pi)*Pi(0)_trunc'] = (2.0 / (3 * np.pi)) * Pi_trunc
    c['(2/3pi)*Pi(0)_inf'] = (2.0 / (3 * np.pi)) * Pi_inf

    # As dressed coupling: B + Pi - Delta
    c['B + Pi(0)_trunc - Delta'] = B_val + Pi_trunc - Delta_val
    c['B + Pi(0)_inf - Delta'] = B_val + Pi_inf - Delta_val
    c['B + (2/3pi)*Pi(0)_inf - Delta'] = (
        B_val + (2.0 / (3 * np.pi)) * Pi_inf - Delta_val)

    # -----------------------------------------------------------------
    # 6. SUNSET DIAGRAMS (two fermion lines, one photon)
    # -----------------------------------------------------------------
    # Unrestricted: sum_{n1,n2} g1*g2 / (|lam1|^2 * |lam2|^2)
    sunset_unrestr = sum(
        dirac_deg(n1) * dirac_deg(n2) /
        (dirac_ev(n1)**2 * dirac_ev(n2)**2)
        for n1 in range(m) for n2 in range(m))
    c['Sunset_unrestr(2,2)'] = sunset_unrestr
    c['sqrt(Sunset_unrestr)'] = np.sqrt(sunset_unrestr)

    # Vertex-restricted (E1: |n1-n2|=1)
    sunset_E1 = sum(
        dirac_deg(n1) * dirac_deg(n2) /
        (dirac_ev(n1)**2 * dirac_ev(n2)**2)
        for n1 in range(m) for n2 in range(m)
        if vertex_allowed(n1, n2))
    c['Sunset_E1(2,2)'] = sunset_E1

    # CG-weighted vertex-restricted
    sunset_CG = sum(
        cg_coupling_sq(min(n1, n2)) *
        dirac_deg(n1) * dirac_deg(n2) /
        (dirac_ev(n1)**2 * dirac_ev(n2)**2)
        for n1 in range(m) for n2 in range(m)
        if vertex_allowed(n1, n2))
    c['Sunset_CG(2,2)'] = sunset_CG

    # Sunset with photon propagator 1/|lam_gamma|^2 (E1 photon modes)
    # Photon on S^3 has scalar eigenvalues: lam_gamma = n^2 - 1
    sunset_photon = 0.0
    for n1 in range(m):
        for n2 in range(m):
            if not vertex_allowed(n1, n2):
                continue
            n_gamma = abs(n1 - n2)  # photon angular momentum
            gamma_ev = max(scalar_ev(n_gamma + 1), 1e-30)  # +1 for Fock convention
            sunset_photon += (dirac_deg(n1) * dirac_deg(n2) /
                              (dirac_ev(n1)**2 * dirac_ev(n2)**2 * gamma_ev))
    c['Sunset_photon(2,2,1)'] = sunset_photon

    # Various exponent choices for sunset
    for s1, s2 in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1)]:
        val = sum(
            dirac_deg(n1) * dirac_deg(n2) /
            (dirac_ev(n1)**(2*s1) * dirac_ev(n2)**(2*s2))
            for n1 in range(m) for n2 in range(m)
            if vertex_allowed(n1, n2))
        c[f'Sunset_E1({s1},{s2})'] = val

    # -----------------------------------------------------------------
    # 7. RENORMALIZATION-GROUP MOTIVATED
    # -----------------------------------------------------------------
    # Each shell contributes to running: delta(1/alpha) ~ g_n * log(lambda_n/lambda_0)
    if m >= 2:
        rg_sum = sum(dirac_deg(n) * np.log(dirac_ev(n) / dirac_ev(0))
                     for n in range(1, m))
        c['RG_sum'] = rg_sum
        c['(2/3pi)*RG_sum'] = (2.0 / (3 * np.pi)) * rg_sum
        c['B + (2/3pi)*RG_sum - Delta'] = (
            B_val + (2.0 / (3 * np.pi)) * rg_sum - Delta_val)

    # -----------------------------------------------------------------
    # 8. SPECTRAL DETERMINANT RATIOS
    # -----------------------------------------------------------------
    if m >= 2:
        ratio_log = log_det_D2 - log_det_S
        c['log(det D^2 / det scalar)'] = ratio_log
        c['log(det scalar / det D^2)'] = -ratio_log

    # -----------------------------------------------------------------
    # 9. MIXED: B + (QED quantity) - Delta
    # -----------------------------------------------------------------
    for key_suffix in [
        'D_trunc(s=2.0)', 'D_trunc(s=4.0)',
        'D_inf(s=2.0)', 'D_inf(s=4.0)',
        'log det|D|_trunc', '(1/2) log det(D^2)',
        'E_cas_scalar', 'E_cas_Dirac / 4',
        'Sunset_CG(2,2)', 'Sunset_E1(2,2)',
    ]:
        if key_suffix in c:
            c[f'B + {key_suffix} - Delta'] = B_val + c[key_suffix] - Delta_val

    # -----------------------------------------------------------------
    # 10. EULER-HEISENBERG TYPE: effective Lagrangian from integrating out fermions
    # -----------------------------------------------------------------
    # L_EH ~ sum_n g_n * |lambda_n|^{-4} * (something involving field strength)
    # At zero field: just the spectral zeta at s=2
    # Already captured in D_trunc(s=4.0) and D_inf(s=4.0) (note: |lam|^{-2s})

    # Asymptotic expansion of the effective action at large Lambda:
    # S_eff ~ f_0 * a_0 + f_1 * a_1 + f_2 * a_2 + ...
    # where a_k are Seeley-DeWitt coefficients
    # a_0 = sqrt(pi), a_1 = sqrt(pi), a_2 = sqrt(pi)/8 on unit S^3
    sqrt_pi = np.sqrt(np.pi)
    # f_k = int_0^infty f(x) x^{(3-k)/2 - 1} dx  (d=3)
    # For f(x) = exp(-x): f_k = Gamma((3-k)/2)
    from scipy.special import gamma as gamma_func
    for k in range(3):
        a_k = sqrt_pi * [1.0, 1.0, 1.0/8][k]
        f_k = gamma_func((3 - k) / 2.0)
        c[f'SD_term(k={k})'] = f_k * a_k

    SD_total = sum(gamma_func((3 - k) / 2.0) * sqrt_pi * [1, 1, 1.0/8][k]
                   for k in range(3))
    c['SD_total(k=0..2)'] = SD_total

    # Ratio of SD total to Lambda^3 factor
    if m >= 1:
        Lambda = dirac_ev(m)
        c[f'SD_total * Lambda^3'] = SD_total * Lambda**3
        c[f'SD_total * Lambda^2'] = SD_total * Lambda**2
        c[f'SD_total * Lambda'] = SD_total * Lambda

    # -----------------------------------------------------------------
    # 11. NOVEL: spectral zeta derivatives and logarithmic moments
    # -----------------------------------------------------------------
    for s in [2.0, 4.0]:
        log_moment = sum(dirac_deg(n) * dirac_ev(n)**(-s) * np.log(dirac_ev(n))
                         for n in range(m))
        c[f'D_log_moment(s={s})'] = log_moment

    # -----------------------------------------------------------------
    # 12. FUNCTIONAL DETERMINANT with mass
    # -----------------------------------------------------------------
    for mass in [0.0, 0.5, 1.0, 1.5]:
        val = sum(dirac_deg(n) * np.log(dirac_ev(n)**2 + mass**2)
                  for n in range(m))
        c[f'Tr log(D^2+m^2), m={mass}'] = val
        c[f'-(1/2) Tr log(D^2+m^2), m={mass}'] = -0.5 * val

    # -----------------------------------------------------------------
    # 13. HEAT KERNEL TRACES at various t
    # -----------------------------------------------------------------
    for t in [0.1, 0.5, 1.0, np.pi/2, np.pi, 2*np.pi]:
        K_D = sum(dirac_deg(n) * np.exp(-t * dirac_ev(n)**2)
                  for n in range(m))
        K_S = sum(scalar_deg(n) * np.exp(-t * scalar_ev(n))
                  for n in range(1, m + 1))
        c[f'K_Dirac(t={t:.3f})'] = K_D
        c[f'K_scalar(t={t:.3f})'] = K_S
        if K_D > 0:
            c[f'-log K_Dirac(t={t:.3f})'] = -np.log(K_D)

    # -----------------------------------------------------------------
    # 14. SELF-ENERGY with vertex weights
    # -----------------------------------------------------------------
    # Sigma(n) = sum_{n'} |V_{nn'}|^2 / (|lam_n| + |lam_{n'}|)
    # summed over all external legs
    SE_total = 0.0
    for n in range(m):
        for np_ in range(m):
            if not vertex_allowed(n, np_):
                continue
            SE_total += (dirac_deg(n) * dirac_deg(np_) * cg_coupling_sq(min(n, np_)) /
                         (dirac_ev(n) + dirac_ev(np_)))
    c['Self_energy_CG'] = SE_total
    c['B + Self_energy_CG - Delta'] = B_val + SE_total - Delta_val

    # -----------------------------------------------------------------
    # 15. SPECTRAL ACTION with cutoff function f(x) = (1-x)*theta(1-x)
    # -----------------------------------------------------------------
    Lambda_m = dirac_ev(m)
    for cut_label, cut_fn in [
        ('sharp', lambda x: 1.0 if x < 1 else 0.0),
        ('linear', lambda x: max(0, 1 - x)),
        ('gaussian', lambda x: np.exp(-x**2)),
    ]:
        SA_D = sum(dirac_deg(n) * cut_fn(dirac_ev(n)**2 / Lambda_m**2)
                   for n in range(m + 3))
        SA_S = sum(scalar_deg(n) * cut_fn(scalar_ev(n) / Lambda_m**2)
                   for n in range(1, m + 3))
        c[f'SA_D({cut_label},Lam=|lam_m|)'] = SA_D
        c[f'SA_S({cut_label},Lam=|lam_m|)'] = SA_S
        c[f'SA_Str({cut_label})'] = SA_S - SA_D / 4.0
        c[f'B + SA_Str({cut_label}) - Delta'] = B_val + (SA_S - SA_D / 4.0) - Delta_val

    # -----------------------------------------------------------------
    # 16. PARTIAL F with QED correction
    # -----------------------------------------------------------------
    F_partial = sum(1.0 / n**2 for n in range(1, m + 1))
    c['F_partial'] = F_partial
    c['B + F_partial - Delta'] = B_val + F_partial - Delta_val

    # F_partial weighted by QED vertex factor
    F_qed = sum(cg_coupling_sq(n - 1) / n**2 for n in range(1, m + 1))
    c['F_QED_weighted'] = F_qed

    # F via Dirac zeta: D_{n^2}(4) = zeta(2) = F
    F_dirac = sum(scalar_deg(n) * float(n)**(-4) for n in range(1, N_INF + 1))
    c['F_from_D_{n^2}(4)'] = F_dirac  # should be pi^2/6

    # -----------------------------------------------------------------
    # 17. BARE + DRESSED splitting
    # -----------------------------------------------------------------
    # Idea: 1/alpha = 1/alpha_bare + (radiative corrections)
    # alpha_bare from graph topology (B), corrections from QED loops
    # Try: K/pi = B + (one-loop) - (boundary)

    # One-loop vacuum energy density on S^3
    if m >= 1:
        vac_energy = -0.5 * sum(dirac_deg(n) * dirac_ev(n) for n in range(m))
        c['Vac_energy_Dirac'] = vac_energy
        c['-Vac_energy_Dirac / (dim S^3)'] = -vac_energy / 3.0

    # -----------------------------------------------------------------
    # 18. NOVEL: Cross-sector interaction sums
    # -----------------------------------------------------------------
    # Scalar-Dirac cross zeta: Sigma_{n_S, n_D} g_S * g_D / (lam_S * lam_D^2)
    if m >= 2:
        cross = sum(
            scalar_deg(ns) * dirac_deg(nd) /
            (scalar_ev(ns) * dirac_ev(nd)**2)
            for ns in range(2, m + 1) for nd in range(m))
        c['Cross_zeta(1,2)'] = cross

    # Shell-matched: sum_n g_S(n+1) * g_D(n) * f(eigenvalues)
    shell_product = sum(
        scalar_deg(n + 1) * dirac_deg(n) /
        (scalar_ev(n + 1) * dirac_ev(n))
        for n in range(m) if n + 1 >= 2)
    c['Shell_product'] = shell_product

    # -----------------------------------------------------------------
    # 19. SEELEY-DEWITT SUPERTRACE (non-perturbative remainder)
    # -----------------------------------------------------------------
    # Already known from ST sprint: Str[f(D^2/Lam^2)] perturbative = 0 (F1).
    # Try: non-perturbative piece = exact sum - SD asymptotic
    if m >= 2:
        for lam_sq in [dirac_ev(m - 1)**2, dirac_ev(m)**2]:
            exact_D = sum(dirac_deg(n) * np.exp(-dirac_ev(n)**2 / lam_sq)
                          for n in range(N_INF))
            exact_S = sum(scalar_deg(n) * np.exp(-scalar_ev(n) / lam_sq)
                          for n in range(1, N_INF + 1))
            str_exact = exact_S - exact_D / 4.0
            c[f'Str_nonpert(Lam2={lam_sq:.2f})'] = str_exact

    return c, K_data


def main():
    print("=" * 78)
    print("alpha-SA QED PROBE: QED-motivated spectral quantities on S^3")
    print("=" * 78)
    print(f"Target: 1/alpha = {ALPHA_INV:.10f}")
    print(f"Target: K/pi    = {K_OVER_PI:.10f}")

    all_results = {}

    for m in range(1, 6):
        print(f"\n{'=' * 60}")
        print(f"  m = {m}")
        print(f"{'=' * 60}")

        candidates, K_data = compute_qed_candidates(m)
        target = K_data['K_over_pi']
        print(f"  K(m)/pi = {target:.6f}  "
              f"[B={K_data['B']:.0f}, F={K_data['F']:.6f}, "
              f"Delta={K_data['Delta']:.6f}]")

        ranked = []
        for name, val in candidates.items():
            if val is not None and np.isfinite(val) and val != 0:
                rel_err = abs(val - target) / abs(target) if target != 0 else 999
                ranked.append((name, val, rel_err))
        ranked.sort(key=lambda x: x[2])

        print(f"\n  Top 30 closest to K(m)/pi = {target:.6f}:")
        for name, val, rel_err in ranked[:30]:
            marker = "<<<" if rel_err < 1e-6 else ("<~" if rel_err < 0.05 else "")
            print(f"    {name:55s} = {val:14.6f}  "
                  f"rel_err={rel_err:.6f} {marker}")

        all_results[m] = [(name, float(val), float(rel_err))
                          for name, val, rel_err in ranked]

    # ==================================================================
    # TRACKING TEST
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("TRACKING: candidates that track K(m)/pi across m=2..5")
    print(f"{'=' * 78}")

    name_set = None
    val_dicts = {}
    for m in range(2, 6):
        d = {name: val for name, val, _ in all_results[m]}
        val_dicts[m] = d
        if name_set is None:
            name_set = set(d.keys())
        else:
            name_set &= set(d.keys())

    tracking = []
    for name in name_set:
        ratios = []
        for m in range(2, 6):
            t = K_ingredients(m)['K_over_pi']
            v = val_dicts[m][name]
            if t != 0:
                ratios.append((m, v / t))
        vals = [r for _, r in ratios]
        if len(vals) == 4 and np.mean(vals) > 0:
            cv = np.std(vals) / np.mean(vals) if np.mean(vals) != 0 else 999
            if cv < 0.15:
                tracking.append((name, cv, ratios))
    tracking.sort(key=lambda x: x[1])

    for name, cv, ratios in tracking[:20]:
        r_str = ", ".join(f"({m},{r:.4f})" for m, r in ratios)
        print(f"  {name:55s} CV={cv:.4f}  [{r_str}]")

    # ==================================================================
    # SAVE
    # ==================================================================
    out = {}
    for m in all_results:
        out[str(m)] = [{'name': n, 'value': v, 'rel_err': e}
                       for n, v, e in all_results[m]]
    with open('debug/data/alpha_sa_qed_probe.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to debug/data/alpha_sa_qed_probe.json")


if __name__ == '__main__':
    main()
