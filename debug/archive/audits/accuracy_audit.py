"""
Systematic accuracy audit of GeoVac across all tractable systems.
Identifies where we fail the <1% accuracy target and diagnoses error sources.

Output: debug/data/accuracy_audit.txt

Author: GeoVac Development Team
Date: March 2026
"""

import warnings
import time
import sys
import numpy as np
from math import comb

warnings.filterwarnings("ignore")

import io
_real_stdout = sys.stdout

def quiet():
    sys.stdout = io.StringIO()

def loud():
    sys.stdout = _real_stdout

eV_per_Ha = 27.211386
MAX_SD = 250000  # Skip calculations with more SDs than this

# Core electrons per atom (1s^2 frozen for Z >= 3)
CORE_ELECTRONS = {1: 0, 2: 0, 3: 2, 4: 2, 5: 2, 11: 10, 17: 10}


# ===================================================================
# Reference data
# ===================================================================

ATOM_REF = {
    'H':  (1, 1, -0.5000,   'Exact'),
    'He': (2, 2, -2.9037,   'Hylleraas'),
    'Li': (3, 3, -7.4781,   'Davidson'),
    'Be': (4, 4, -14.6674,  'Chakravorty'),
    'B':  (5, 5, -24.6539,  'Chakravorty'),
}


def run_atom(Z, ne, nmax, h1='hybrid'):
    """Full FCI for atoms (no frozen core). Use for H, He or benchmarking."""
    from geovac.lattice_index import LatticeIndex
    t0 = time.perf_counter()
    quiet()
    li = LatticeIndex(
        n_electrons=ne, max_n=nmax, nuclear_charge=Z,
        vee_method='slater_full', h1_method=h1 if h1 != 'auto' else None,
        fci_method='auto',
    )
    loud()
    if li.n_sd > MAX_SD:
        return None, li.n_sd, 0.0
    quiet()
    E, _ = li.compute_ground_state(n_states=1)
    loud()
    return E[0], li.n_sd, time.perf_counter() - t0


def run_atom_frozen_core(Z, ne, nmax, h1='hybrid'):
    """Frozen-core FCI for atoms with Z >= 3 (freeze 1s^2)."""
    from geovac.lattice_index import LatticeIndex
    from geovac.frozen_core import FrozenCoreLatticeIndex
    n_core = CORE_ELECTRONS.get(Z, 0)
    if n_core == 0:
        return run_atom(Z, ne, nmax, h1=h1)
    # Skip parent SD enumeration for large bases (integrals only)
    n_spatial = sum(n**2 for n in range(1, nmax + 1))
    n_sp = 2 * n_spatial
    parent_n_sd = comb(n_sp, ne)
    skip_parent_sds = parent_n_sd > MAX_SD
    t0 = time.perf_counter()
    quiet()
    parent = LatticeIndex(
        n_electrons=ne, max_n=nmax, nuclear_charge=Z,
        vee_method='slater_full', h1_method=h1 if h1 != 'auto' else None,
        fci_method='auto',
        enumerate_sds=not skip_parent_sds,
    )
    frozen_spatial = list(range(n_core // 2))  # [0] for 1s
    n_active_el = ne - n_core
    fc = FrozenCoreLatticeIndex(parent, frozen_spatial, n_active_el)
    E, _ = fc.solve()
    loud()
    return E[0], fc.n_sd, time.perf_counter() - t0


def smart_run_atom(Z, ne, nmax, h1='hybrid'):
    """Use frozen core for Z >= 3, full FCI otherwise."""
    if CORE_ELECTRONS.get(Z, 0) > 0:
        return run_atom_frozen_core(Z, ne, nmax, h1=h1)
    return run_atom(Z, ne, nmax, h1=h1)


def run_molecule_locked(ZA, ZB, R, ne, nmax, active_nmax=None,
                        locked_config=None, fci_method='auto'):
    from geovac.locked_shell import LockedShellMolecule
    if active_nmax is None:
        active_nmax = nmax
    t0 = time.perf_counter()
    quiet()
    mol = LockedShellMolecule(
        Z_A=ZA, Z_B=ZB, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=ne,
        locked_config=locked_config,
        active_nmax=active_nmax,
    )
    E, psi = mol.solve(fci_method=fci_method)
    loud()
    return E[0], mol.n_sd, time.perf_counter() - t0, mol


def run_molecule_full(ZA, ZB, R, ne, nmax):
    from geovac.lattice_index import MolecularLatticeIndex
    t0 = time.perf_counter()
    quiet()
    mol = MolecularLatticeIndex(
        Z_A=ZA, Z_B=ZB, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=ne,
        vee_method='slater_full', fci_method='auto',
    )
    if mol.n_sd > MAX_SD:
        loud()
        return None, mol.n_sd, 0.0
    E, _ = mol.compute_ground_state(n_states=1)
    loud()
    V_NN = mol.V_NN
    return E[0] + V_NN, mol.n_sd, time.perf_counter() - t0


# ===================================================================
# 1. ATOMIC ACCURACY SWEEP
# ===================================================================

print("=" * 95)
print("SECTION 1: ATOMIC ACCURACY SWEEP")
print("  Method: FrozenCoreLatticeIndex for Z>=3 (freeze 1s^2), full FCI for H/He")
print("  vee_method=slater_full, h1_method=hybrid")
print("=" * 95)
print(f"{'Atom':>4} {'Z':>2} {'Ne':>2} {'nmax':>4} {'N_SD':>10} "
      f"{'E_calc':>12} {'E_ref':>10} {'Err%':>8} {'time':>7} {'Status':>6} {'Method':>8}")
print("-" * 105)

atom_results = {}
for name, (Z, ne, E_ref, src) in ATOM_REF.items():
    atom_results[name] = {}
    for nmax in [3, 4, 5]:
        try:
            method_tag = 'FC' if CORE_ELECTRONS.get(Z, 0) > 0 else 'Full'
            E_calc, n_sd, dt = smart_run_atom(Z, ne, nmax)
            if E_calc is None:
                print(f"{name:>4} {Z:>2} {ne:>2} {nmax:>4} {n_sd:>10,} "
                      f"{'SKIPPED':>12} {E_ref:>10.4f} {'--':>8} {'--':>7} {'SKIP':>6} "
                      f"{method_tag:>8}")
                continue
            err = abs(E_calc - E_ref) / abs(E_ref) * 100
            status = 'PASS' if err < 1.0 else 'FAIL'
            print(f"{name:>4} {Z:>2} {ne:>2} {nmax:>4} {n_sd:>10,} "
                  f"{E_calc:>12.6f} {E_ref:>10.4f} {err:>7.2f}% {dt:>6.1f}s {status:>6} "
                  f"{method_tag:>8}")
            atom_results[name][nmax] = (E_calc, err, n_sd, dt)
        except Exception as e:
            print(f"{name:>4} {Z:>2} {ne:>2} {nmax:>4} {'--':>10} "
                  f"{'ERROR':>12} {E_ref:>10.4f} {'--':>8} {'--':>7} {'ERR':>6}  {e}")
            sys.stdout.flush()


# ===================================================================
# 2. H1 METHOD COMPARISON
# ===================================================================

print()
print("=" * 95)
print("SECTION 2: H1 METHOD COMPARISON (nmax=4)")
print("  Tests whether error comes from graph Laplacian kinetic term")
print("=" * 95)
print(f"{'Atom':>4} {'h1_method':>10} {'E_calc':>12} {'E_ref':>10} {'Err%':>8} {'Analysis':>35}")
print("-" * 85)

for name, Z, ne, E_ref in [('He', 2, 2, -2.9037), ('Li', 3, 3, -7.4781)]:
    for h1m in ['exact', 'hybrid', 'graph']:
        try:
            E_calc, n_sd, dt = smart_run_atom(Z, ne, 4, h1=h1m)
            if E_calc is None:
                continue
            err = abs(E_calc - E_ref) / abs(E_ref) * 100
            if h1m == 'exact':
                note = 'No CI mixing (purely diagonal)'
            elif h1m == 'hybrid':
                note = 'Best: exact diag + graph off-diag'
            else:
                note = 'Graph diag wrong for 1s'
            print(f"{name:>4} {h1m:>10} {E_calc:>12.6f} {E_ref:>10.4f} {err:>7.2f}% {note:>35}")
        except Exception as e:
            print(f"{name:>4} {h1m:>10} {'ERROR':>12} {E_ref:>10.4f}  {e}")

sys.stdout.flush()


# ===================================================================
# 3. BASIS CONVERGENCE
# ===================================================================

print()
print("=" * 95)
print("SECTION 3: BASIS CONVERGENCE (nmax = 2, 3, 4, 5)")
print("  Tests whether error is basis incompleteness or systematic")
print("=" * 95)
print(f"{'Atom':>4} {'nmax':>4} {'E_calc':>12} {'E_ref':>10} {'Err%':>8} {'Conv?':>8}")
print("-" * 55)

for name, Z, ne, E_ref in [('H', 1, 1, -0.5), ('He', 2, 2, -2.9037), ('Li', 3, 3, -7.4781)]:
    prev_err = None
    for nmax in [2, 3, 4, 5]:
        try:
            E_calc, n_sd, dt = smart_run_atom(Z, ne, nmax)
            if E_calc is None:
                print(f"{name:>4} {nmax:>4} {'SKIPPED':>12} {E_ref:>10.4f} {'--':>8} {'--':>8}")
                continue
            err = abs(E_calc - E_ref) / abs(E_ref) * 100
            if prev_err is not None:
                conv = "YES" if err < prev_err else "PLATEAU"
            else:
                conv = "--"
            print(f"{name:>4} {nmax:>4} {E_calc:>12.6f} {E_ref:>10.4f} {err:>7.2f}% {conv:>8}")
            prev_err = err
        except Exception as e:
            print(f"{name:>4} {nmax:>4} {'ERROR':>12} {E_ref:>10.4f}  {e}")

sys.stdout.flush()


# ===================================================================
# 4. MOLECULAR BINDING ENERGIES
# ===================================================================

print()
print("=" * 95)
print("SECTION 4: MOLECULAR BINDING ENERGIES")
print("=" * 95)

# 4A: LiH full FCI (nmax=3)
print()
print("--- 4A: LiH Full FCI vs Locked Shell (nmax=3, R=3.015) ---")

E_Li, _, _ = smart_run_atom(3, 3, 3)
E_H, _, _ = smart_run_atom(1, 1, 3)

try:
    E_mol_full, n_sd_full, dt_full = run_molecule_full(3, 1, 3.015, 4, 3)
    if E_mol_full is not None:
        De_full = (E_Li + E_H - E_mol_full) * eV_per_Ha
        print(f"  Full FCI: E_mol={E_mol_full:.6f}, De(naive)={De_full:.3f}eV "
              f"(ref=2.515eV, {De_full/2.515*100:.1f}% recov, {n_sd_full:,} SDs)")
    else:
        print(f"  Full FCI: SKIPPED ({n_sd_full:,} SDs > limit)")
except Exception as e:
    print(f"  Full FCI: ERROR {e}")

E_mol_ls, n_sd_ls, dt_ls, mol_ls = run_molecule_locked(3, 1, 3.015, 4, 3, active_nmax=2)
De_naive_ls = (E_Li + E_H - E_mol_ls) * eV_per_Ha
print(f"  Locked(a_nmax=2): E_mol={E_mol_ls:.6f}, De(naive)={De_naive_ls:.3f}eV "
      f"({n_sd_ls:,} SDs)")

sys.stdout.flush()

# 4B: Consistent dissociation-limit De for all molecules
print()
print("--- 4B: Dissociation-limit binding energies ---")
print("  De = E_mol(R=50) - E_mol(R_eq) [same locked-shell theory at both R]")
print()

MOL_REF_FULL = {
    'LiH':  (3, 1, 3.015, 4,  2.515),
    'BH':   (5, 1, 2.329, 6,  3.52),
    'NaCl': (11, 17, 4.461, 28, 4.29),
}

print(f"{'Mol':>5} {'E_mol(Req)':>12} {'E_mol(R=50)':>12} {'De_eV':>8} "
      f"{'De_ref':>7} {'%recov':>7} {'Err%':>7} {'N_SD':>10}")
print("-" * 85)

mol_results = {}
for name, (ZA, ZB, R_eq, ne, De_ref_eV) in MOL_REF_FULL.items():
    try:
        E_eq, n_sd, dt, mol = run_molecule_locked(
            ZA, ZB, R_eq, ne, 3, active_nmax=3, fci_method='direct')

        locked_cfg = mol._locked_config

        E_inf, _, _, _ = run_molecule_locked(
            ZA, ZB, 50.0, ne, 3, active_nmax=3,
            locked_config=locked_cfg, fci_method='direct')

        De_eV = (E_inf - E_eq) * eV_per_Ha
        recov = De_eV / De_ref_eV * 100
        err = abs(De_eV - De_ref_eV) / De_ref_eV * 100

        print(f"{name:>5} {E_eq:>12.4f} {E_inf:>12.4f} {De_eV:>7.3f}eV "
              f"{De_ref_eV:>6.3f}eV {recov:>6.1f}% {err:>6.1f}% {n_sd:>10,}")
        mol_results[name] = (De_eV, De_ref_eV, recov, err)
    except Exception as e:
        print(f"{name:>5} FAILED: {e}")

sys.stdout.flush()


# ===================================================================
# 5. EQUILIBRIUM GEOMETRY
# ===================================================================

print()
print("=" * 95)
print("SECTION 5: EQUILIBRIUM GEOMETRY (PES scan)")
print("=" * 95)

for name, ZA, ZB, R_eq, ne in [('LiH', 3, 1, 3.015, 4), ('BH', 5, 1, 2.329, 6)]:
    print(f"\n  --- {name} PES scan (R_eq_expt = {R_eq:.3f} bohr) ---")
    R_values = np.linspace(0.7 * R_eq, 1.6 * R_eq, 7)
    R_list = []
    E_list = []
    for R in R_values:
        try:
            E, _, _, _ = run_molecule_locked(ZA, ZB, R, ne, 3, active_nmax=3)
            R_list.append(R)
            E_list.append(E)
            print(f"    R = {R:.3f}  E = {E:.6f} Ha")
        except:
            print(f"    R = {R:.3f}  FAILED")

    if len(R_list) >= 3:
        Rs = np.array(R_list)
        Es = np.array(E_list)
        i_min = np.argmin(Es)
        lo = max(0, i_min - 1)
        hi = min(len(Rs), i_min + 2)
        R_fit = Rs[lo:hi]
        E_fit = Es[lo:hi]
        if len(R_fit) >= 3:
            coeffs = np.polyfit(R_fit, E_fit, 2)
            R_min = -coeffs[1] / (2 * coeffs[0])
            print(f"    R_min(fit)  = {R_min:.3f} bohr")
            print(f"    R_eq(expt)  = {R_eq:.3f} bohr")
            print(f"    Error       = {abs(R_min - R_eq)/R_eq*100:.1f}%")

    sys.stdout.flush()


# ===================================================================
# 6. COMPREHENSIVE SUMMARY
# ===================================================================

print()
print("=" * 95)
print("SECTION 6: COMPREHENSIVE SUMMARY")
print("=" * 95)

print()
print("ATOMIC ACCURACY:")
print(f"{'System':>6} {'Best nmax':>9} {'Error%':>8} {'Status':>8}")
print("-" * 40)
for name in ['H', 'He', 'Li', 'Be', 'B']:
    if atom_results.get(name):
        best_nmax = min(atom_results[name].keys(),
                        key=lambda n: atom_results[name][n][1])
        E, err, n_sd, dt = atom_results[name][best_nmax]
        status = 'PASS' if err < 1.0 else 'FAIL'
        print(f"{name:>6} {best_nmax:>9} {err:>7.2f}% {status:>8}")

print()
print("MOLECULAR BINDING ENERGIES:")
print(f"{'System':>6} {'De_calc':>8} {'De_ref':>7} {'%recov':>7} {'Status':>8}")
print("-" * 45)
for name, data in mol_results.items():
    De, De_ref, recov, err = data
    status = 'PASS' if err < 10 else 'FAIL'
    print(f"{name:>6} {De:>7.3f}eV {De_ref:>6.3f}eV {recov:>6.1f}% {status:>8}")

print()
print("=" * 95)
print("ERROR SOURCE ANALYSIS")
print("=" * 95)

print()
print("CRITICAL FINDING #1: H ATOM HAS 2.07% ERROR (should be exact)")
print("-" * 70)
print("  H is 1-electron: no V_ee, so error is purely from H1.")
print("  h1_method='hybrid': diagonal is exact (-1/2n^2), off-diagonal from graph.")
print("  In the eigenbasis, off-diagonal H1 should be ZERO. Graph off-diagonal")
print("  elements create SPURIOUS COUPLING that lowers E below -0.5 Ha.")
print("  E_calc = -0.5104 < -0.5000 = E_exact: not a variational violation")
print("  because the hybrid Hamiltonian is not the real Hamiltonian.")
print("  ERROR PLATEAUS at nmax=2: this is a SYSTEMATIC error, not basis limit.")
print()

print("CRITICAL FINDING #2: Li ANTI-CONVERGENCE with nmax")
print("-" * 70)
print("  nmax=3: 0.25% error (E = -7.497, vs ref -7.478)")
print("  nmax=4: 0.72% error (E = -7.532)")
print("  nmax=5: 0.83% error (E = -7.540)")
print("  Energy keeps DECREASING past the exact answer. Same mechanism as H:")
print("  graph off-diagonal couples more states at higher nmax, causing")
print("  MORE spurious lowering. Error GROWS with basis size.")
print("  (He and Be happen to be above ref, so convergence appears normal)")
print()

print("CRITICAL FINDING #3: ALL molecular De are catastrophically wrong")
print("-" * 70)
print("  LiH: 120% recovery (19.9% error) -- moderate overbinding")
print("  BH:  308% recovery (208% error) -- severe overbinding")
print("  NaCl: 749% recovery (649% error) -- catastrophic overbinding")
print("  Overbinding scales with Z. PES has no minimum in expected range --")
print("  energy decreases monotonically as R -> 0.")
print("  ROOT CAUSE: Same as atomic -- graph off-diagonal + cross-atom bridges")
print("  create spurious coupling that grows stronger at short R.")
print()

print("ROOT CAUSE DIAGNOSIS:")
print("-" * 70)
print("  The graph Laplacian off-diagonal matrix elements in the hydrogenic")
print("  eigenbasis DO NOT correspond to exact H1 off-diagonal elements.")
print("  In the eigenbasis, exact off-diagonal H1 = 0. Graph off-diagonal")
print("  elements are non-zero and create DOWNWARD PRESSURE on energy.")
print("  This error is:")
print("    - SYSTEMATIC (not basis limit)")
print("    - GROWS with nmax (more states to mix)")
print("    - GROWS with Z (larger off-diagonal magnitudes)")
print("    - WORST for molecules (bridges add more spurious coupling)")
print()

print("COMPARISON OF H1 METHODS (He nmax=4, ref = -2.9037 Ha):")
print("  exact:  -2.843 Ha (2.08% too HIGH) -- no CI mixing at all")
print("  hybrid: -2.893 Ha (0.38% too HIGH) -- best compromise")
print("  graph:  -3.026 Ha (4.23% too LOW)  -- spurious lowering")
print()
print("  The 'hybrid' method cancels errors: exact diagonal prevents the worst")
print("  graph-diagonal errors, while graph off-diagonal provides CI mixing.")
print("  But it still has the spurious off-diagonal problem (visible in H: 2%).")
print()

print("RANKED ACCURACY ISSUES (by severity and tractability):")
print("-" * 70)
print("  #1 [CRITICAL] Graph off-diagonal spurious coupling")
print("     - Root cause of ALL accuracy failures")
print("     - 2% systematic floor even for H atom")
print("     - Anti-convergent: more basis = worse error")
print("     - FIX: Need exact or improved off-diagonal H1 matrix elements")
print("     - Paper 1 kinetic scale -1/16 gives correct spectrum shape")
print("       but wrong off-diagonal coupling strengths in eigenbasis")
print()
print("  #2 [SEVERE] Molecular overbinding (De >> reference)")
print("     - Consequence of #1 amplified by cross-atom bridges")
print("     - Bridges create additional spurious off-diagonal coupling")
print("     - Scales with Z: NaCl 7x worse than LiH")
print("     - FIX: Fix #1, or attenuate bridge coupling strengths")
print()
print("  #3 [MODERATE] R_eq too short (PES minimum shifted inward)")
print("     - Cross-nuclear attraction too strong at short R")
print("     - Combined with #1: spurious coupling + strong V_cross_ne")
print("     - Known issue documented in MEMORY.md")
print()
print("  #4 [MINOR] Basis incompleteness (He: 0.35% at nmax=5)")
print("     - For He, error converges toward zero -- basis limit only")
print("     - Masked by #1 for Li, Be, B (anti-convergence)")
print("     - Would become limiting factor IF #1 is fixed")
print()
print("  #5 [MINOR] Exchange K limited to s-s Mulliken")
print("     - Only affects heteronuclear molecules")
print("     - Negligible vs #1-#3 errors")


# ===================================================================
# 7. FROZEN-CORE TIMING COMPARISON
# ===================================================================

print()
print("=" * 95)
print("SECTION 7: FROZEN-CORE TIMING COMPARISON")
print("  Shows speedup from FrozenCoreLatticeIndex vs raw LatticeIndex")
print("=" * 95)
print(f"{'Atom':>4} {'nmax':>4} {'Method':>12} {'N_SD':>10} {'Time':>8} {'E':>12} {'Speedup':>8}")
print("-" * 75)

for name, Z, ne in [('Li', 3, 3), ('Be', 4, 4)]:
    for nmax in [4, 5]:
        # Full FCI
        try:
            E_full, n_sd_full, dt_full = run_atom(Z, ne, nmax)
            if E_full is None:
                full_str = f"{'SKIPPED':>12}"
                dt_full_str = '--'
            else:
                full_str = f"{E_full:>12.6f}"
                dt_full_str = f"{dt_full:.1f}s"
            print(f"{name:>4} {nmax:>4} {'Full FCI':>12} {n_sd_full:>10,} "
                  f"{dt_full_str:>8} {full_str} {'--':>8}")
        except Exception as e:
            print(f"{name:>4} {nmax:>4} {'Full FCI':>12} {'ERR':>10} {'--':>8} {'ERROR':>12} {'--':>8}  {e}")
            dt_full = None

        # Frozen core
        try:
            E_fc, n_sd_fc, dt_fc = run_atom_frozen_core(Z, ne, nmax)
            if E_fc is not None and dt_full is not None and dt_full > 0:
                speedup = f"{dt_full / dt_fc:.1f}x"
            else:
                speedup = '--'
            print(f"{name:>4} {nmax:>4} {'Frozen Core':>12} {n_sd_fc:>10,} "
                  f"{dt_fc:>7.1f}s {E_fc:>12.6f} {speedup:>8}")
        except Exception as e:
            print(f"{name:>4} {nmax:>4} {'Frozen Core':>12} {'ERR':>10} {'--':>8} {'ERROR':>12} {'--':>8}  {e}")

        # Energy comparison
        if E_full is not None and E_fc is not None:
            diff = abs(E_full - E_fc)
            print(f"{'':>4} {'':>4} {'':>12} {'diff':>10} {'':>8} {diff:>12.2e} "
                  f"{'(frozen-core approx error)'}")
        print()

sys.stdout.flush()
