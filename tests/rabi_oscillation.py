"""
Rabi Oscillation Test — High-Precision Weak-Field Benchmark
=============================================================

v0.7.0: The Time Machine.

Proves that the GeoVac lattice produces quantitatively correct coherent
quantum dynamics by operating in the strict weak-field regime where the
rotating-wave approximation (RWA) is exact.

Lattice resolution strategy:
    The graph Laplacian spectrum becomes DENSER at higher max_n, causing
    multi-level leakage that prevents clean 2-level Rabi oscillations.
    At max_n=10 (385 states), the 85 dipole-coupled eigenstates are too
    closely spaced for any field strength to isolate a single transition.
    At max_n=4 (30 states), eigenstate 1 is spectrally isolated (81x
    the Rabi frequency from any competitor), yielding 99.98% population
    transfer — a near-perfect 2-level system.

    We use max_n=10 for norm conservation (proving the propagator at
    production scale) and max_n=4 for Rabi precision (where eigenstates
    are spectrally pure and the RWA is exact).

Perturbation:
    Weak oscillating electric field: V(t) = E0 * z * cos(omega * t)
    E0 = 0.0005 a.u. (deep in the perturbative / 2-level regime)
    Resonant frequency: omega = E_target - E_gs (from lattice eigenvalues)

Validation:
    1. Norm conservation:  ||psi(t)|| = 1 to machine precision  (max_n=10)
    2. Rabi accuracy:      Peak P_target > 0.95, period error < 1.0%  (max_n=4)
    3. Off-resonance:      Detuned drive produces suppressed transfer  (max_n=4)

Date: February 20, 2026
"""

import numpy as np
import sys
import io
import time

sys.path.insert(0, '.')

# Ensure UTF-8 output on Windows
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from scipy.linalg import eigh
from geovac import AtomicSolver
from geovac.dynamics import TimePropagator


# ======================================================================
# Shared constants
# ======================================================================
E0 = 0.0005  # Driving field amplitude (weak-field regime)


# ======================================================================
# System builder: full diag + dipole-based target selection
# ======================================================================
def build_rabi_system(max_n: int = 4):
    """
    Build the hydrogen system with full diagonalization and identify the
    best dipole-coupled eigenstate for clean 2-level Rabi oscillations.

    Uses spectral isolation criterion: picks the eigenstate whose gap
    to the nearest significantly coupled competitor is maximised relative
    to its Rabi frequency.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number for the lattice

    Returns
    -------
    info : dict
    """
    solver = AtomicSolver(max_n=max_n, Z=1)
    N = solver.n_states

    # Full diagonalization — all N eigenstates
    H_dense = solver.H.toarray()
    E_all, psi_all = eigh(H_dense)

    psi_gs = psi_all[:, 0]
    E_gs = E_all[0]

    # Build dipole operator
    V_z = TimePropagator.build_dipole_z(solver.lattice)
    V_z_dense = V_z.toarray()

    # Compute all dipole couplings to ground state
    couplings = []
    for i in range(1, N):
        d_i = abs(float(psi_gs @ V_z_dense @ psi_all[:, i]))
        dE_i = E_all[i] - E_gs
        couplings.append((i, d_i, dE_i))

    # Sorted by coupling strength for diagnostics
    couplings_sorted = sorted(couplings, key=lambda x: -x[1])
    print(f"  Full diag: {N} eigenstates")
    print(f"  Top dipole couplings to ground state:")
    for rank in range(min(5, len(couplings_sorted))):
        idx, d, dE = couplings_sorted[rank]
        print(f"    #{rank+1}: eigenstate {idx:3d}, "
              f"mu = {d:.4f}, dE = {dE:.6f} Ha")

    # Find the best-isolated dipole-coupled eigenstate.
    # Require mu > 0.5 so we pick transitions strong enough for a short
    # propagation (T_half < ~3000 a.u.) where numerical drift is negligible.
    mu_threshold = 0.5
    coupled = [(i, d, dE) for (i, d, dE) in couplings if d > mu_threshold]

    best_j = 0
    best_ratio = 0.0
    for j in range(len(coupled)):
        _, mu_j, dE_j = coupled[j]
        Omega_j = mu_j * E0
        if Omega_j <= 0:
            continue
        # Minimum gap to any other significantly coupled state
        min_gap = float('inf')
        for k in range(len(coupled)):
            if k == j:
                continue
            _, mu_k, dE_k = coupled[k]
            if mu_k < 0.1 * mu_j:
                continue
            gap = abs(dE_j - dE_k)
            if gap < min_gap:
                min_gap = gap
        if min_gap == float('inf'):
            min_gap = dE_j
        isolation = min_gap / Omega_j
        if isolation > best_ratio:
            best_ratio = isolation
            best_j = j

    idx_target = coupled[best_j][0]
    mu_12 = coupled[best_j][1]
    psi_target = psi_all[:, idx_target]
    E_target = E_all[idx_target]
    omega_res = E_target - E_gs

    print(f"\n  Selected eigenstate {idx_target} "
          f"(mu = {mu_12:.4f}, dE = {omega_res:.6f}, "
          f"isolation = {best_ratio:.1f}x)")

    return {
        'solver': solver,
        'E_all': E_all, 'psi_all': psi_all,
        'psi_gs': psi_gs, 'psi_target': psi_target,
        'E_gs': E_gs, 'E_target': E_target,
        'omega_res': omega_res,
        'V_z': V_z, 'mu_12': mu_12,
        'idx_target': idx_target,
        'isolation': best_ratio,
    }


# ======================================================================
# TEST 1: Norm Conservation (Unitarity) — max_n=10
# ======================================================================
def test_1_norm_conservation() -> dict:
    """
    Verify that Crank-Nicolson preserves ||psi|| = 1.0 under static H.

    Uses max_n=10 (385 states) to prove the propagator at production scale.
    """
    print("\n" + "#" * 70)
    print("TEST 1: NORM CONSERVATION (max_n=10, 385 states)")
    print("#" * 70)

    solver = AtomicSolver(max_n=10, Z=1)
    E, psi = solver.compute_ground_state(n_states=1)
    psi0 = psi[:, 0].astype(complex)

    print(f"\n  System: Hydrogen (Z=1), max_n=10, N_states={solver.n_states}")
    print(f"  E_gs = {E[0]:.6f} Ha")
    print(f"  Initial norm: {float(np.abs(np.vdot(psi0, psi0))):.15f}")

    dt = 0.1
    n_steps = 1000
    prop = TimePropagator(solver.H, dt=dt)

    norms = []

    def track_norm(step: int, t: float, psi_t: np.ndarray) -> None:
        norms.append(float(np.abs(np.vdot(psi_t, psi_t))))

    t0 = time.time()
    psi_final = prop.evolve(psi0, n_steps, callback=track_norm)
    t1 = time.time()

    norm_final = float(np.abs(np.vdot(psi_final, psi_final)))
    norm_max_dev = max(abs(n - 1.0) for n in norms)
    overlap = float(abs(np.vdot(psi[:, 0], psi_final)))

    print(f"  Final norm:    {norm_final:.15f}")
    print(f"  Max deviation: {norm_max_dev:.2e}")
    print(f"  Steps: {n_steps}, dt={dt}, T={n_steps*dt:.1f} a.u.")
    print(f"  Wall time: {(t1-t0)*1000:.0f} ms")
    print(f"  |<gs|psi(T)>| = {overlap:.15f}")

    passed = norm_max_dev < 1e-10
    print(f"\n  Status: {'OK PASS' if passed else 'FAIL'}"
          f" (norm deviation < 1e-10)")

    return {
        'test': 'norm_conservation',
        'norm_final': norm_final,
        'max_deviation': norm_max_dev,
        'passed': passed,
    }


# ======================================================================
# TEST 2: Weak-Field Rabi Oscillation (Precision Benchmark) — max_n=4
# ======================================================================
def test_2_rabi_oscillation(sys_info: dict) -> dict:
    """
    Drive the gs -> target transition with a weak resonant field.

    E0 = 0.0005 a.u. ensures strict 2-level dynamics (no leakage).
    max_n=4 provides spectral isolation of 81x the Rabi frequency.

    Validates:
      - Peak P_target > 0.95   (clean population transfer)
      - Period error < 1.0%    (quantitative RWA agreement)
    """
    print("\n" + "#" * 70)
    print("TEST 2: WEAK-FIELD RABI OSCILLATION (max_n=4, 30 states)")
    print("#" * 70)

    solver = sys_info['solver']
    psi_gs = sys_info['psi_gs']
    psi_target = sys_info['psi_target']
    E_gs = sys_info['E_gs']
    E_target = sys_info['E_target']
    omega_res = sys_info['omega_res']
    V_z = sys_info['V_z']
    mu_12 = sys_info['mu_12']
    idx_target = sys_info['idx_target']

    print(f"\n  System: Hydrogen (Z=1), max_n={solver.max_n}, "
          f"N_states={solver.n_states}")
    print(f"  E_gs     = {E_gs:.6f} Ha  (eigenstate 0)")
    print(f"  E_target = {E_target:.6f} Ha  (eigenstate {idx_target})")
    print(f"  Resonant frequency: omega = {omega_res:.6f} Ha")
    print(f"  Dipole moment: mu = {mu_12:.6f} a.u.")
    print(f"  Spectral isolation: {sys_info['isolation']:.1f}x")

    # Analytical predictions (RWA)
    Omega_R = mu_12 * E0
    T_half = np.pi / Omega_R
    T_rabi = 2.0 * np.pi / Omega_R

    print(f"\n  --- Analytical Predictions (RWA) ---")
    print(f"  Field strength:  E0 = {E0} a.u.")
    print(f"  Rabi frequency:  Omega_R = {Omega_R:.6e} Ha")
    print(f"  Half period:     T_half = {T_half:.2f} a.u.")
    print(f"  Full period:     T_Rabi = {T_rabi:.2f} a.u.")

    # Simulation parameters
    dt = 0.1
    T_target = T_half * 1.15  # go 15% past peak to capture it
    n_steps = int(np.ceil(T_target / dt))

    print(f"\n  --- Simulation Setup ---")
    print(f"  dt = {dt} a.u.")
    print(f"  n_steps = {n_steps}")
    print(f"  T_sim = {n_steps * dt:.2f} a.u.")
    print(f"  Steps per T_half: {T_half / dt:.0f}")

    # Evolve
    psi0 = psi_gs.astype(complex)
    prop = TimePropagator(solver.H, dt=dt)

    populations = []
    norms = []

    def track(step: int, t: float, psi_t: np.ndarray) -> None:
        p_gs = float(abs(np.vdot(psi_gs, psi_t))**2)
        p_tgt = float(abs(np.vdot(psi_target, psi_t))**2)
        norm = float(np.abs(np.vdot(psi_t, psi_t)))
        populations.append((t, p_gs, p_tgt))
        norms.append(norm)

    t0 = time.time()
    prop.evolve_driven(
        psi0, n_steps,
        H0=solver.H.tocsr(),
        V_dipole=V_z,
        E0=E0,
        omega=omega_res,
        callback=track,
    )
    t1 = time.time()
    wall_time = t1 - t0

    # Extract arrays
    times = np.array([p[0] for p in populations])
    p_gs_vals = np.array([p[1] for p in populations])
    p_tgt_vals = np.array([p[2] for p in populations])

    # Find peak
    idx_peak = int(np.argmax(p_tgt_vals))
    t_sim_peak = times[idx_peak]
    p_tgt_peak = p_tgt_vals[idx_peak]
    p_gs_at_peak = p_gs_vals[idx_peak]
    leakage = 1.0 - p_gs_at_peak - p_tgt_peak

    # Period accuracy
    period_error = abs(t_sim_peak - T_half) / T_half * 100.0
    norm_max_dev = max(abs(n - 1.0) for n in norms)

    print(f"\n  --- Results ---")
    print(f"  Peak P_target:      {p_tgt_peak:.6f}  "
          f"(target: > 0.95)")
    print(f"  P_gs at peak:       {p_gs_at_peak:.6f}")
    print(f"  Leakage (1-Pgs-Pt): {leakage:.6f}")
    print(f"  t_sim(peak):        {t_sim_peak:.2f} a.u.")
    print(f"  T_half(theory):     {T_half:.2f} a.u.")
    print(f"  Period error:       {period_error:.4f}%  "
          f"(target: < 1.0%)")
    print(f"  Norm max deviation: {norm_max_dev:.2e}")
    print(f"  Wall time:          {wall_time:.1f}s")

    # Population trace
    n_samples = min(20, len(populations))
    sample_indices = np.linspace(0, len(populations) - 1, n_samples, dtype=int)
    print(f"\n    {'t':>10}  {'P_gs':>8}  {'P_tgt':>8}  {'Sum':>8}")
    print(f"    {'----------':>10}  {'--------':>8}  {'--------':>8}  {'--------':>8}")
    for i in sample_indices:
        t_i, p1, p2 = populations[i]
        print(f"    {t_i:10.2f}  {p1:8.5f}  {p2:8.5f}  {p1+p2:8.5f}")

    # Validation
    pass_population = p_tgt_peak > 0.95
    pass_period = period_error < 1.0
    pass_norm = norm_max_dev < 1e-6
    passed = pass_population and pass_period and pass_norm

    print(f"\n  --- Validation ---")
    print(f"  [{'PASS' if pass_population else 'FAIL'}] "
          f"Peak P_target = {p_tgt_peak:.4f} > 0.95")
    print(f"  [{'PASS' if pass_period else 'FAIL'}] "
          f"Period error = {period_error:.4f}% < 1.0%")
    print(f"  [{'PASS' if pass_norm else 'FAIL'}] "
          f"Norm deviation = {norm_max_dev:.2e} < 1e-6")
    print(f"\n  Status: {'OK PASS' if passed else 'FAIL'}")

    return {
        'test': 'rabi_oscillation',
        'p_tgt_peak': p_tgt_peak,
        'p_gs_at_peak': p_gs_at_peak,
        'leakage': leakage,
        't_sim_peak': t_sim_peak,
        'T_half_theory': T_half,
        'period_error_pct': period_error,
        'norm_dev': norm_max_dev,
        'passed': passed,
    }


# ======================================================================
# TEST 3: Off-Resonance Check — max_n=4
# ======================================================================
def test_3_off_resonance(sys_info: dict) -> dict:
    """
    Drive at 2x resonant frequency. Population transfer should be
    strongly suppressed compared to the resonant case.
    """
    print("\n" + "#" * 70)
    print("TEST 3: OFF-RESONANCE CHECK (Detuned Driving)")
    print("#" * 70)

    solver = sys_info['solver']
    psi_gs = sys_info['psi_gs']
    psi_target = sys_info['psi_target']
    omega_res = sys_info['omega_res']
    V_z = sys_info['V_z']
    mu_12 = sys_info['mu_12']

    Omega_R = mu_12 * E0
    T_half = np.pi / Omega_R

    # Detune to 2x resonant frequency
    omega_detuned = 2.0 * omega_res
    detuning = omega_detuned - omega_res

    # RWA prediction for detuned max transfer
    Omega_eff = np.sqrt(Omega_R**2 + detuning**2)
    max_transfer_rwa = (Omega_R / Omega_eff)**2

    print(f"\n  Resonant omega:       {omega_res:.6f} Ha")
    print(f"  Detuned omega:        {omega_detuned:.6f} Ha")
    print(f"  Detuning delta:       {detuning:.6f} Ha")
    print(f"  Omega_R:              {Omega_R:.6e}")
    print(f"  Max transfer (RWA):   {max_transfer_rwa:.6f}")
    print(f"  E0 = {E0} a.u.")

    dt = 0.1
    n_steps = int(np.ceil(T_half * 1.15 / dt))

    psi0 = psi_gs.astype(complex)
    prop = TimePropagator(solver.H, dt=dt)

    state = {'p_tgt_max': 0.0}

    def make_track(tgt, s):
        def track(step: int, t: float, psi_t: np.ndarray) -> None:
            p = float(abs(np.vdot(tgt, psi_t))**2)
            if p > s['p_tgt_max']:
                s['p_tgt_max'] = p
        return track

    t0 = time.time()
    prop.evolve_driven(
        psi0, n_steps,
        H0=solver.H.tocsr(),
        V_dipole=V_z,
        E0=E0,
        omega=omega_detuned,
        callback=make_track(psi_target, state),
    )
    t1 = time.time()

    p_tgt_max = state['p_tgt_max']
    pass_suppressed = p_tgt_max < 0.5
    passed = pass_suppressed

    print(f"\n  Peak P_target (off-res): {p_tgt_max:.6f}")
    print(f"  Wall time: {(t1-t0)*1000:.0f} ms")
    print(f"\n  [{'PASS' if pass_suppressed else 'FAIL'}] "
          f"P_target = {p_tgt_max:.4f} < 0.50 (suppressed vs resonant)")
    print(f"\n  Status: {'OK PASS' if passed else 'FAIL'}")

    return {
        'test': 'off_resonance',
        'p_tgt_max': p_tgt_max,
        'max_transfer_rwa': max_transfer_rwa,
        'passed': passed,
    }


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    t_start = time.time()

    print("=" * 70)
    print("RABI OSCILLATION TEST: Weak-Field Precision Benchmark")
    print("=" * 70)
    print(f"\nv0.7.0: The Time Machine")
    print(f"Method: Crank-Nicolson unitary propagator")
    print(f"Regime: Weak-field (E0 = {E0} a.u.), 2-level dynamics")
    print(f"Goal:   Quantitative match to analytical Rabi theory")

    # Build Rabi system (max_n=4 for spectral purity)
    print(f"\n  Building Rabi system (max_n=4, full diag)...")
    sys_info = build_rabi_system(max_n=4)
    print(f"  Done. {sys_info['solver'].n_states} states, "
          f"mu = {sys_info['mu_12']:.4f} a.u., "
          f"isolation = {sys_info['isolation']:.1f}x")

    r1 = test_1_norm_conservation()
    r2 = test_2_rabi_oscillation(sys_info)
    r3 = test_3_off_resonance(sys_info)

    # ---- SUMMARY ----
    print(f"\n{'='*70}")
    print(f"RABI OSCILLATION SUMMARY")
    print(f"{'='*70}")

    all_passed = r1['passed'] and r2['passed'] and r3['passed']

    print(f"\n  {'Test':<30}  {'Result':>16}  {'Status':>10}")
    print(f"  {'-'*30}  {'-'*16}  {'-'*10}")

    dev_str = f"dev={r1['max_deviation']:.1e}"
    s1 = 'OK PASS' if r1['passed'] else 'FAIL'
    print(f"  {'Norm conservation (n=10)':<30}  {dev_str:>16}  {s1:>10}")

    p_str = f"P={r2['p_tgt_peak']:.4f}"
    s2 = 'OK PASS' if r2['passed'] else 'FAIL'
    print(f"  {'Rabi peak transfer (n=4)':<30}  {p_str:>16}  {s2:>10}")

    err_str = f"err={r2['period_error_pct']:.4f}%"
    print(f"  {'Rabi period accuracy':<30}  {err_str:>16}  {'':>10}")

    leak_str = f"leak={r2['leakage']:.4f}"
    print(f"  {'2-level purity':<30}  {leak_str:>16}  {'':>10}")

    off_str = f"P={r3['p_tgt_max']:.4f}"
    s3 = 'OK PASS' if r3['passed'] else 'FAIL'
    print(f"  {'Off-resonance suppression':<30}  {off_str:>16}  {s3:>10}")

    n_passed = sum(1 for r in [r1, r2, r3] if r['passed'])
    print(f"\n  Result: {n_passed}/3 tests passed")

    if all_passed:
        print(f"\n  The Lattice is not just an eigenvalue solver.")
        print(f"  It is a Time Machine.")
        print(f"  Weak-field Rabi dynamics: quantitatively confirmed.")

    t_total = time.time() - t_start
    print(f"\n  Total time: {t_total:.1f}s")
    print(f"{'='*70}")
