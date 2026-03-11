"""v0.9.33 Diagnostic: Exact Exchange K Integrals via Poisson Solve.

Computes exact exchange integrals K_kl using the same Poisson machinery
as J integrals (v0.9.32), replacing Ohno-Klopman monopole approximation.

Diagnostics:
  1. K integral table and K/J ratios at p0=sqrt(10)
  2. Gap projection: estimate Delta_E_K before running FCI
  3. Full FCI at R=3.015 with exact K
  4. Term decomposition (H1, J, K, V_NN)
  5. PES scan: R = 2.0 to 6.0
  6. Ten-configuration comparison
"""

import sys
import time
import numpy as np

sys.path.insert(0, '.')

from geovac.molecular_sturmian import (
    compute_molecular_sturmian_betas,
    compute_exact_j_integrals,
    compute_exact_k_integrals,
)
from geovac.mo_fci import MOSturmianFCI


def main() -> None:
    Z_A, Z_B = 3.0, 1.0
    R_eq = 3.015
    nmax = 3
    n_el = 4
    n_grid = 40
    E_exact = -8.071
    E_atoms = -7.892  # Li(-7.478) + H(-0.5) at FCI limit; approximate

    lines = []
    def log(s: str = '') -> None:
        print(s)
        lines.append(s)

    log("=" * 70)
    log("v0.9.33 Diagnostic: Exact Exchange K Integrals")
    log("=" * 70)

    # ---------------------------------------------------------------
    # 1. K integral table at reference p0
    # ---------------------------------------------------------------
    p0_ref = np.sqrt(10.0)  # v0.9.32 initial p0
    log(f"\n--- 1. K Integral Table (p0={p0_ref:.4f}, R={R_eq}) ---")

    mo_results = compute_molecular_sturmian_betas(
        Z_A=Z_A, Z_B=Z_B, R=R_eq, p0=p0_ref, nmax=nmax
    )
    valid_mos = [mo for mo in mo_results if np.isfinite(mo[4])]
    n_mo = len(valid_mos)

    t0 = time.perf_counter()
    J = compute_exact_j_integrals(
        mo_results, Z_A, Z_B, R_eq, p0_ref, n_grid=30
    )
    t_j = time.perf_counter() - t0

    t0 = time.perf_counter()
    K = compute_exact_k_integrals(
        mo_results, Z_A, Z_B, R_eq, p0_ref, n_grid=30
    )
    t_k = time.perf_counter() - t0

    log(f"  n_mo = {n_mo}, J time = {t_j:.2f}s, K time = {t_k:.2f}s")

    # Diagonal check: K_kk = J_kk
    log("\n  Diagonal check (K_kk vs J_kk):")
    for k in range(min(n_mo, 5)):
        mo = valid_mos[k]
        rel = abs(K[k, k] - J[k, k]) / max(abs(J[k, k]), 1e-10) * 100
        log(f"    MO {k} (n={mo[0]},m={mo[1]}): J={J[k,k]:.6f}, "
            f"K={K[k,k]:.6f}, err={rel:.4f}%")

    # Full K matrix (upper triangle, same-m pairs)
    log("\n  K matrix (nonzero off-diagonal, same-m pairs):")
    log(f"  {'k':>3} {'l':>3} {'m_k':>3} {'m_l':>3} {'K_kl':>10} {'J_kl':>10} {'K/J':>8}")
    for k in range(n_mo):
        for l in range(k + 1, n_mo):
            if valid_mos[k][1] != valid_mos[l][1]:
                continue  # different m
            kj_ratio = K[k, l] / J[k, l] if abs(J[k, l]) > 1e-10 else 0.0
            log(f"  {k:3d} {l:3d} {valid_mos[k][1]:3d} {valid_mos[l][1]:3d} "
                f"{K[k,l]:10.6f} {J[k,l]:10.6f} {kj_ratio:8.4f}")

    # Selection rule check
    log("\n  Selection rule (cross-m pairs, should be ~0):")
    for k in range(min(n_mo, 4)):
        for l in range(k + 1, n_mo):
            if valid_mos[k][1] == valid_mos[l][1]:
                continue
            log(f"    K[{k},{l}] (m={valid_mos[k][1]},m={valid_mos[l][1]}): "
                f"{K[k,l]:.2e}")
            if l - k > 2:
                break

    # ---------------------------------------------------------------
    # 2. Gap projection: estimate Delta_E_K
    # ---------------------------------------------------------------
    log(f"\n--- 2. Gap Projection at p0={p0_ref:.4f} ---")

    # Ohno-Klopman K for comparison
    mo_n_vals = [mo[0] for mo in valid_mos]
    K_ohno = np.zeros((n_mo, n_mo))
    for k in range(n_mo):
        nk = mo_n_vals[k]
        r_k = nk**2 / max(Z_A, 0.5)
        for l in range(n_mo):
            if k == l:
                K_ohno[k, l] = J[k, l]
                continue
            nl = mo_n_vals[l]
            r_l = nl**2 / max(Z_B, 0.5)
            K_ohno[k, l] = 1.0 / np.sqrt(R_eq**2 + ((r_k + r_l) / 2.0)**2)

    # For 4 electrons in 1sigma^2 2sigma^2:
    # Exchange pairs: (0alpha,1alpha) and (0beta,1beta) -> 2 same-spin pairs
    # DeltaE_K = -2 * (K_exact[0,1] - K_ohno[0,1])
    log(f"  K_exact[0,1] = {K[0,1]:.6f} Ha")
    log(f"  K_ohno[0,1]  = {K_ohno[0,1]:.6f} Ha")
    delta_k_01 = K[0, 1] - K_ohno[0, 1]
    delta_E_K = -2.0 * delta_k_01
    log(f"  Delta_K[0,1]  = {delta_k_01:.6f} Ha")
    log(f"  Delta_E_K (dominant) = -2 * Delta_K = {delta_E_K:.6f} Ha")
    log(f"  v0.9.32 E_mol = -7.131 Ha, E_atoms = {E_atoms}")
    log(f"  Gap to bound = {E_atoms - (-7.131):.3f} Ha")
    log(f"  Expected E_mol(v0.9.33) ~ {-7.131 + delta_E_K:.3f} Ha")
    if delta_E_K < 0:
        log(f"  --> Exact K should LOWER energy by {abs(delta_E_K):.3f} Ha")
    else:
        log(f"  --> Exact K would RAISE energy by {delta_E_K:.3f} Ha (unexpected)")

    # ---------------------------------------------------------------
    # 3. Full FCI at R=3.015
    # ---------------------------------------------------------------
    log(f"\n--- 3. Full FCI at R={R_eq} ---")

    fci = MOSturmianFCI(
        Z_A=Z_A, Z_B=Z_B, R=R_eq, nmax=nmax,
        n_electrons=n_el, n_grid=n_grid, xi_max=12.0,
    )
    t0 = time.perf_counter()
    E_mol, p0_star = fci.solve(damping=0.5, max_iter=20, verbose=True)
    t_fci = time.perf_counter() - t0

    err_pct = (E_mol - E_exact) / abs(E_exact) * 100
    log(f"\n  Result: E_mol = {E_mol:.6f} Ha, p0* = {p0_star:.4f}")
    log(f"  Exact:  E_exact = {E_exact} Ha")
    log(f"  Error:  {err_pct:.2f}%")
    log(f"  Time:   {t_fci:.1f}s")
    log(f"  v0.9.32 p0* = 3.777 -> v0.9.33 p0* = {p0_star:.4f}")

    bound = E_mol < E_atoms
    log(f"  Bound:  {'YES' if bound else 'NO'} (E_mol={E_mol:.4f} vs E_atoms={E_atoms})")
    if bound:
        D_e = E_atoms - E_mol
        log(f"  D_e:    {D_e:.4f} Ha (expt 0.092)")

    # ---------------------------------------------------------------
    # 4. Term decomposition
    # ---------------------------------------------------------------
    log(f"\n--- 4. Energy Decomposition at R={R_eq} ---")
    decomp = fci.decompose_energy()
    log(f"  E_H1     = {decomp['E_H1']:.6f} Ha")
    log(f"  E_J      = {decomp['E_J']:.6f} Ha")
    log(f"  E_K      = {decomp['E_K']:.6f} Ha")
    log(f"  V_NN     = {decomp['V_NN']:.6f} Ha")
    log(f"  E_sum    = {decomp['E_sum']:.6f} Ha")
    log(f"  E_fci    = {decomp['E_fci']:.6f} Ha")
    log(f"  Residual = {decomp['residual']:.6f} Ha (off-diagonal CI)")

    # K/J ratios for occupied MO pairs in raw basis
    log("\n  K/J ratios (raw MO basis, occupied pairs):")
    J_raw = fci._J_raw
    K_raw = fci._K_raw
    for k in range(min(n_mo, 3)):
        for l in range(k, min(n_mo, 3)):
            if abs(J_raw[k, l]) > 1e-10:
                ratio = K_raw[k, l] / J_raw[k, l]
                log(f"    K[{k},{l}]/J[{k},{l}] = {K_raw[k,l]:.6f}/{J_raw[k,l]:.6f} = {ratio:.4f}")

    # ---------------------------------------------------------------
    # 5. PES scan
    # ---------------------------------------------------------------
    R_values = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
    log(f"\n--- 5. PES Scan ---")
    log(f"  {'R':>6} {'E_mol':>10} {'p0*':>8} {'D_e':>8} {'bound':>6}")

    pes_data = []
    for R_val in R_values:
        fci_r = MOSturmianFCI(
            Z_A=Z_A, Z_B=Z_B, R=R_val, nmax=nmax,
            n_electrons=n_el, n_grid=n_grid, xi_max=12.0,
        )
        E_r, p0_r = fci_r.solve(damping=0.5, max_iter=20, verbose=False)
        D_e_r = E_atoms - E_r if E_r < E_atoms else 0.0
        is_bound = E_r < E_atoms
        log(f"  {R_val:6.3f} {E_r:10.6f} {p0_r:8.4f} {D_e_r:8.4f} {'YES' if is_bound else 'NO':>6}")
        pes_data.append((R_val, E_r, p0_r))

    # Dissociation diagnostic
    E_R6 = pes_data[-1][1]
    diss_err = abs(E_R6 - E_atoms)
    log(f"\n  Dissociation: E(R=6) = {E_R6:.4f}, |delta| = {diss_err:.4f} Ha")
    log(f"  Dissociation check: {'PASS' if diss_err < 0.05 else 'FAIL'} (tol 0.05 Ha)")

    # Find minimum
    E_min = min(pes_data, key=lambda x: x[1])
    log(f"  PES minimum: R={E_min[0]:.3f}, E={E_min[1]:.6f}")
    if E_min[1] < E_atoms:
        log(f"  D_e at minimum: {E_atoms - E_min[1]:.4f} Ha")
    log(f"  R_eq check: {'PASS' if 2.5 <= E_min[0] <= 3.5 else 'FAIL'} "
        f"(R_eq={E_min[0]:.3f}, target [2.5, 3.5])")

    # ---------------------------------------------------------------
    # 6. Decomposition at R=6 (dissociation)
    # ---------------------------------------------------------------
    log(f"\n--- 6. Energy Decomposition at R=6.0 ---")
    fci_diss = MOSturmianFCI(
        Z_A=Z_A, Z_B=Z_B, R=6.0, nmax=nmax,
        n_electrons=n_el, n_grid=n_grid, xi_max=12.0,
    )
    E_diss, _ = fci_diss.solve(damping=0.5, max_iter=20, verbose=False)
    decomp_diss = fci_diss.decompose_energy()
    log(f"  E_H1     = {decomp_diss['E_H1']:.6f} Ha")
    log(f"  E_J      = {decomp_diss['E_J']:.6f} Ha")
    log(f"  E_K      = {decomp_diss['E_K']:.6f} Ha")
    log(f"  V_NN     = {decomp_diss['V_NN']:.6f} Ha")
    log(f"  E_sum    = {decomp_diss['E_sum']:.6f} Ha")
    log(f"  E_fci    = {decomp_diss['E_fci']:.6f} Ha")

    # Delta between R=3.015 and R=6.0
    log(f"\n  Delta (R=3.015 -> R=6.0):")
    log(f"    dE_H1 = {decomp['E_H1'] - decomp_diss['E_H1']:.6f} Ha")
    log(f"    dE_J  = {decomp['E_J'] - decomp_diss['E_J']:.6f} Ha")
    log(f"    dE_K  = {decomp['E_K'] - decomp_diss['E_K']:.6f} Ha")
    log(f"    dV_NN = {decomp['V_NN'] - decomp_diss['V_NN']:.6f} Ha")

    # ---------------------------------------------------------------
    # 7. Ten-configuration comparison
    # ---------------------------------------------------------------
    log(f"\n--- 7. Ten-Configuration Comparison ---")
    log(f"  {'Version':>12} {'Method':>25} {'E(R=3.015)':>12} {'D_e':>8} {'R_eq':>8} {'Bound':>6}")
    log(f"  {'v0.9.11':>12} {'atom-cent CP-corr':>25} {'---':>12} {'+0.093':>8} {'~2.5':>8} {'YES':>6}")
    log(f"  {'v0.9.18':>12} {'hybrid Fourier+SW':>25} {'---':>12} {'+0.143':>8} {'<2.0':>8} {'YES':>6}")
    log(f"  {'v0.9.21-30':>12} {'Sturmian attempts':>25} {'varies':>12} {'---':>8} {'---':>8} {'NO':>6}")
    log(f"  {'v0.9.31':>12} {'MO FCI approx V_ee':>25} {'-10.07':>12} {'~1.6':>8} {'<2.0':>8} {'YES':>6}")
    log(f"  {'v0.9.32':>12} {'MO FCI exact J':>25} {'-7.131':>12} {'0':>8} {'~3.5':>8} {'NO':>6}")

    E_933 = pes_data[3][1]  # R=3.015 entry
    D_e_933 = E_atoms - E_933 if E_933 < E_atoms else 0.0
    bound_933 = 'YES' if E_933 < E_atoms else 'NO'
    log(f"  {'v0.9.33':>12} {'MO FCI exact J+K':>25} {E_933:12.3f} {D_e_933:8.3f} "
        f"{'~' + str(E_min[0]):>8} {bound_933:>6}")
    log(f"  {'Expt':>12} {'---':>25} {'-8.071':>12} {'+0.092':>8} {'3.015':>8} {'YES':>6}")

    # ---------------------------------------------------------------
    # Save results
    # ---------------------------------------------------------------
    outfile = 'debug/data/mo_fci_lih_results.txt'
    with open(outfile, 'w') as f:
        f.write('\n'.join(lines))
    print(f"\nResults saved to {outfile}")


if __name__ == '__main__':
    main()
