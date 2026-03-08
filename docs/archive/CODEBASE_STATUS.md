# Codebase Status Report

**Date:** February 14, 2026
**Status:** Post-Unified Architecture Refactor

## CORE MODULES (ESSENTIAL - ACTIVELY USED)

### ✓ CLEAN & REFACTORED

1. **geovac/lattice.py**
   - Status: ✓ Refactored with unified architecture
   - Contains: `GeometricLattice` class
   - Features: Computes `node_weights = -Z/n²` (potential as topology)
   - Used by: All chemistry demos, benchmark suite, MoleculeHamiltonian

2. **geovac/hamiltonian.py** ⚠️ **ESSENTIAL - DO NOT DELETE!**
   - Status: ✓ Refactored with unified architecture
   - Contains: `MoleculeHamiltonian`, `HeliumHamiltonian`, `HeliumPackingSolver`
   - Formula: `H = kinetic_scale*(D - A) + W` (unified!)
   - Used by: chemistry_lab.py, demo_h2.py, benchmark_suite.py
   - **THIS IS THE CORE OF THE SYSTEM!**

3. **geovac/__init__.py**
   - Status: ✓ Clean
   - Exports all public API
   - Defines constants (UNIVERSAL_KINETIC_SCALE, etc.)

4. **geovac/atomic_solver.py**
   - Status: ? Needs review
   - Contains: `AtomicSolver`, `solve_hydrogen`, `solve_atom`
   - Used by: tests/bulk_physics_puzzles.py
   - May duplicate functionality of MoleculeHamiltonian?

### ? SPECIALIZED MODULES (FUNCTIONAL BUT UNUSED)

These import successfully but are NOT used in any tests/demos:

5. **geovac/holographic_analysis.py**
   - Status: ? Research code, not used in benchmarks
   - Contains: Spectral dimension, holographic entropy calculations
   - Used by: NONE (only exported in __init__.py)
   - **Candidate for archive**

6. **geovac/hyperfine_contact.py**
   - Status: ? Research code, not used in benchmarks
   - Used by: NONE (only exported in __init__.py)
   - **Candidate for archive**

7. **geovac/proton_radius_improved.py**
   - Status: ? Research code, not used in benchmarks
   - Used by: NONE (only exported in __init__.py)
   - **Candidate for archive**

8. **geovac/muonic_hydrogen.py**
   - Status: ? Specialized solver
   - Contains: `MuonicHydrogenSolver`
   - Used by: NONE actively (exported in __init__.py)
   - **Keep if it's validated research, otherwise archive**

9. **geovac/dirac_hamiltonian.py**
   - Status: ? Experimental
   - Contains: `DiracHamiltonian` (relativistic)
   - Used by: demo/demo_h2_dirac.py (exists but not in benchmark suite)
   - **Keep if experimental, otherwise archive**

10. **geovac/fundamental_constants.py**
    - Status: ? Research code
    - Contains: Electromagnetic impedance, proton radius predictions
    - Used by: NONE actively
    - **Candidate for archive**

## TESTS & BENCHMARKS

### ✓ WORKING

1. **demo/chemistry_lab.py**
   - Status: ✓ Works with unified architecture
   - Tests: He (1.80%), H- (0.12%), H3+ (5.25%) - ALL PASS!
   - Uses: MoleculeHamiltonian with unified architecture

2. **demo/demo_h2.py**
   - Status: ? Needs verification with unified architecture
   - Tests: H2 molecule

### ⚠️ NEEDS UPDATE

3. **tests/benchmark_suite.py**
   - Status: ⚠️ Using OLD "hybrid mode" (not unified architecture!)
   - Tests: H, He+, Li2+, He (correlation)
   - Missing: He/H-/H3+ chemistry tests from chemistry_lab.py
   - **NEEDS REFACTOR** to use unified architecture

4. **tests/advanced_benchmarks.py**
   - Status: ? Unknown
   - **Needs review**

5. **tests/bulk_physics_puzzles.py**
   - Status: ? Using AtomicSolver
   - **Needs review**

### DIAGNOSTIC SCRIPTS (Can be archived after validation)

6. **demo/diagnose_hydride.py** - H- Z_eff diagnostic
7. **demo/diagnose_kinetic_scale.py** - Kinetic scale scan
8. **demo/find_hydride_solution.py** - H- solution finder
9. **demo/final_hydride_scan.py** - Extended H- scan

**Status:** These were used to solve H- problem. Can be moved to archive.

## DOCUMENTATION

### ✓ CURRENT

1. **docs/UNIFIED_ARCHITECTURE.md** - Complete explanation of refactor
2. **docs/H_MINUS_SOLUTION.md** - H- diagnostic results
3. **docs/DIAGNOSTIC_SUMMARY.md** - Problem-solving process
4. **docs/CHEMISTRY_EXPANSION.md** - Chemistry results

### ? LEGACY

Other docs in docs/ may need review for consistency with unified architecture.

## RECOMMENDATIONS

### IMMEDIATE ACTIONS

1. **✓ DO NOT DELETE hamiltonian.py** - It's the core of the system!

2. **UPDATE benchmark_suite.py** to use unified architecture:
   ```python
   # Old (hybrid mode):
   h = HeliumHamiltonian(max_n=10, Z=2, kinetic_scale=-0.0625)

   # New (unified):
   mol = MoleculeHamiltonian(
       nuclei=[(0,0,0)],
       nuclear_charges=[2],
       max_n=10,
       kinetic_scale=-0.10298808
   )
   ```

3. **ADD chemistry tests to benchmark_suite.py**:
   - He: Target -2.903 Ha, unified architecture
   - H-: Target -0.527 Ha, kinetic_scale=+2.789
   - H3+: Target -1.343 Ha, with V_NN

4. **ARCHIVE unused research code** to `old_research_archive/`:
   - holographic_analysis.py
   - hyperfine_contact.py
   - proton_radius_improved.py
   - fundamental_constants.py
   - muonic_hydrogen.py (unless validated)
   - diagnostic scripts (diagnose_*.py, find_*.py)

5. **VERIFY AtomicSolver** - Does it duplicate MoleculeHamiltonian?
   - If yes: deprecate and use MoleculeHamiltonian
   - If no: keep and document use case

### ARCHIVE STRUCTURE

Proposed:
```
old_research_archive/
├── specialized_solvers/
│   ├── holographic_analysis.py
│   ├── hyperfine_contact.py
│   ├── muonic_hydrogen.py
│   ├── dirac_hamiltonian.py (if unused)
│   └── fundamental_constants.py
├── diagnostics/
│   ├── diagnose_hydride.py
│   ├── diagnose_kinetic_scale.py
│   ├── find_hydride_solution.py
│   └── final_hydride_scan.py
└── README.md (explaining what's archived and why)
```

## SUMMARY

**CLEAN & WORKING:**
- ✓ geovac/lattice.py (unified)
- ✓ geovac/hamiltonian.py (unified) **ESSENTIAL!**
- ✓ demo/chemistry_lab.py (all tests pass)

**NEEDS WORK:**
- ⚠️ tests/benchmark_suite.py (update to unified architecture)
- ⚠️ Add chemistry tests to benchmark suite
- ⚠️ Archive unused research code

**KEY POINT:** `hamiltonian.py` contains `MoleculeHamiltonian` - the CORE class we just refactored! It is absolutely essential and actively used. DO NOT DELETE!
