# GeoVac Project Guidelines

## 🌍 Project Context
**Name:** GeoVac (The Geometric Vacuum)
**Mission:** Spectral Graph Theory approach to Computational Quantum Chemistry. The discrete graph Laplacian is a dimensionless, scale-invariant topology (homologous to the unit S³) that is *mathematically equivalent* to the Schrödinger equation via Fock's 1935 conformal projection, as proven in Paper 7. This equivalence is exploited computationally to replace expensive continuous integration with highly efficient O(N) sparse matrix eigenvalue problems.
**Core Theory:** Papers 0-1, 7 in `papers/core/` (graph Laplacian mechanics, universal constant, dimensionless vacuum)
**Conjectures:** Papers 2-5 in `papers/conjectures/` (emergent spacetime, alpha derivation, holography — exploratory)

---

## 📂 Directory Structure Rules

**ROOT DIRECTORY MUST REMAIN CLEAN.** Only essential files allowed.

### **Core Directories**
- `geovac/` → Core package source code (`__init__.py`, `hamiltonian.py`, `lattice.py`, `dynamics.py`, etc.)
- `papers/` → Theory papers, split into two tiers:
  - `papers/core/` → **Defensible foundations.** Strictly testable O(N) sparse graph Laplacian implementations (Paper 0: geometric packing & universal constant, Paper 1: spectral graph theory & eigenvalue methods, Paper 7: dimensionless vacuum & Schrodinger recovery)
  - `papers/conjectures/` → **Theoretical physics explorations.** Speculative extensions beyond computational QC (Paper 2: alpha derivation, Paper 3: holography, Paper 4: universality, Paper 5: geometric vacuum synthesis, FAQ)
- `tests/` → Unit tests and validation suites (pytest format)

### **Development Directories**
- `debug/` → Scratchpad for analyzing one-off physics issues
  - `debug/plots/` → Generated plots and figures
  - `debug/data/` → Generated data files (.txt, .csv, .json)
  - `debug/*.py` → Debug analysis scripts
- `demo/` → Polished, customer-facing demonstration scripts
  - `demo/demo_h2.py` → H2 molecule example
  - `demo/demo_h2_dirac.py` → Relativistic example
- `benchmarks/` → Performance tracking and regression tests
  - `benchmarks/scripts/` → Benchmark suite scripts
  - `benchmarks/figures/` → Performance plots
  - `benchmarks/BENCHMARKS.md` → Performance documentation

### **Documentation**
- `docs/` → Project documentation, validation reports, analysis
  - `docs/releases/` → Release notes and version summaries
  - `docs/archive/` → Historical/deprecated documentation
  - `docs/*.md` → Current validation reports and analysis
- `old_research_archive/` → Legacy code/theory. **Reference only** - do not import from here.
  - **Archive Rule:** If any archived module is needed by tests or the application, it MUST be moved back into the appropriate source directory (`geovac/`, `ADSCFT/`, `tests/`, etc.). Never import directly from the archive.

### **Theoretical Research** (Retained & Tested)
- `ADSCFT/` → AdS/CFT correspondence and holographic duality research
  - `ADSCFT/boundary/` → Boundary theory (CFT/graph implementation)
  - `ADSCFT/bulk/` → Bulk theory (AdS/geometric embedding with 3D coordinates)
  - `ADSCFT/muonic_hydrogen.py` → Muonic hydrogen solver (mass-independence tests)
  - `ADSCFT/holographic_analysis.py` → Holographic entropy, central charge, spectral dimension
  - `ADSCFT/tests/` → Validation tests for AdS/CFT predictions
  - **Status:** Retained and tested via `tests/advanced_benchmarks.py`
  - **Note:** Isolated from main `geovac/` package. Import as `from ADSCFT import ...`
  - **Purpose:** Bridge graph topology (current) to geometric embedding (needed for fine structure, detailed contact geometry)

### **Files Allowed in Root**

✅ **Allowed:**
- `README.md` - Project overview and quick start
- `LICENSE` - MIT license
- `CHANGELOG.md` - Version history
- `CLAUDE.md` - This file (AI guidelines)
- `setup.py` - Package installation
- `.gitignore`, `.github/` - Git configuration

❌ **NOT Allowed:**
- `.py` scripts → Move to `debug/`, `demo/`, `tests/`, or `benchmarks/`
- `.png`/`.pdf` files → Move to `debug/plots/` or `papers/figures/`
- `.txt`/`.csv`/`.json` data → Move to `debug/data/`
- Documentation `.md` files → Move to `docs/`

---

## 📛 File Naming Conventions

### **Python Scripts**
- `demo_*.py` → Demo scripts (goes in `demo/`)
- `validate_*.py` → Validation scripts (goes in `debug/`)
- `test_*.py` → Unit tests (goes in `tests/`)
- `benchmark_*.py` → Benchmarks (goes in `benchmarks/scripts/`)
- `debug_*.py` → Debug analysis (goes in `debug/`)

### **Documentation**
- `RELEASE_NOTES_v*.md` → Version release notes (`docs/releases/`)
- `RELEASE_SUMMARY_v*.md` → Version summaries (`docs/releases/`)
- `*_VALIDATION_REPORT.md` → Validation reports (`docs/`)
- `*_ANALYSIS.md` → Technical analysis (`docs/`)
- `README.md` → Only in root

### **Data & Figures**
- `*.png`, `*.pdf`, `*.svg` → Plots (`debug/plots/` or `papers/figures/`)
- `*.txt`, `*.csv`, `*.json` → Data (`debug/data/`)
- `*.tex`, `*.bib`, `*.aux`, `*.log` → Papers (`papers/core/` or `papers/conjectures/`)

---

## ⚡ Coding Standards

### **1. Sparse vs Dense: Context-Dependent**

- **Hamiltonian and CI matrices (N > 100):** Always use `scipy.sparse` (csr_matrix,
  coo_matrix). Never densify.
- **Hot-loop lookup tables (ERI, h1 in direct CI):** Use dense NumPy arrays
  when the array fits comfortably in memory (n_spinorb ≤ ~300 at nmax=5 → 300²
  floats = 720 KB). `scipy.sparse._validate_indices` overhead (~24µs/call) is
  prohibitive at 100K+ lookups per assembly.

Rule of thumb: sparse for the physics matrix (N_SD × N_SD), dense for the
orbital-index lookup tables (n_spinorb × n_spinorb or n_spatial⁴).

```python
# ✅ GOOD — large CI matrix
from scipy.sparse import csr_matrix
H = csr_matrix((n_sd, n_sd))

# ✅ GOOD — small lookup table in tight loop
H1_dense = H1_sparse.toarray()  # n_spatial × n_spatial, used 100k+ times
eri_4d = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))  # < 100 MB

# ❌ BAD — densifying a large CI matrix
H = H.toarray()  # N_SD × N_SD can be 500k × 500k
```

### **2. Type Hints Required**
Use Python type hinting for all function signatures.

```python
# ✅ GOOD
def compute_ground_state(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    ...

# ❌ BAD
def compute_ground_state(self, n_states=1):
    ...
```

### **3. Physical Constants**
Import from `geovac.constants` or define at module top. Do not hardcode magic numbers.

**Exception:** `-1/16` is the universal topological constant (can be used directly).

```python
# ✅ GOOD
from geovac.constants import HBAR, C, ALPHA

UNIVERSAL_KINETIC_SCALE = -1/16  # Topological constant

# ❌ BAD
energy = 137.036 * ...  # What is 137.036?
```

### **4. Vectorization Over Loops**
Avoid Python loops for graph operations; use NumPy masking/vectorization.

```python
# ✅ GOOD
mask = (n_values >= 1) & (l_values < n_values)
states_filtered = states[mask]

# ❌ BAD
states_filtered = []
for state in states:
    if state.n >= 1 and state.l < state.n:
        states_filtered.append(state)
```

---

## 🔬 The Dimensionless Vacuum Principle (CRITICAL)

This principle governs all development on the GeoVac codebase. It was formally proven via 18/18 symbolic proofs (`tests/test_fock_projection.py`, `tests/test_fock_laplacian.py`) and documented in Paper 7.

### **The Principle**
The discrete graph Laplacian has **no intrinsic physical scale**. It is a pure, dimensionless combinatorial topology homologous to the unit three-sphere S³. The continuous Schrodinger equation — including the 1/r Coulomb potential and dimensionful energy levels E_n = -1/(2n²) — emerges from projecting this dimensionless topology into flat R³ coordinates via Fock's 1935 stereographic projection.

### **Mathematical Equivalence vs. Physical Priority**
The 18 symbolic proofs establish a *mathematical equivalence* between the discrete graph, the unit S³, and the Schrödinger equation — not a proof that the graph is physically more fundamental. Per Paper 7 Section VI.B, these are equivalent representations under conformal projection. Claims about ontological priority are interpretive, not proven. Code and documentation should reflect this distinction.

### **Prime Directive**
**Never attempt to modify the discrete graph Laplacian to artificially recover continuous differential terms (like 1/r or ∇²).** The graph is an exact, dimensionless S³ topology. The Schrodinger equation is merely its flat-space projection. If the graph eigenvalues do not match expected physics, the issue is in the projection or the energy-shell constraint — not in the graph itself.

### **Key Mathematical Facts**
1. **Conformal factor** Ω = 2p₀/(p² + p₀²) always produces the **unit** S³, regardless of p₀. The intrinsic geometry has no scale.
2. **Eigenvalues** of the Laplace-Beltrami operator on unit S³ are pure integers: λ_n = -(n² - 1). No physical dimensions.
3. **Energy** enters only through the energy-shell constraint p₀² = -2E, which acts as a "stereographic focal length." Energy is the coordinate penalty for flattening curved topology.
4. **The 1/r Coulomb potential** is not an input force law — it is the coordinate distortion created by the stereographic projection (chordal distance identity).

### **Topological Integrity Test**
The symbolic validation suite (`tests/test_fock_projection.py` + `tests/test_fock_laplacian.py`) is the **foundational topological integrity check** for the entire framework. These 18 tests must:
- **Never be broken or bypassed** by any code change
- **Always pass** before any release
- Be run alongside physics benchmarks: `pytest tests/test_fock_projection.py tests/test_fock_laplacian.py -v`

---

## 🧪 Workflow Protocols

### **1. The "Theory Check" Rule**

Before implementing new physics (e.g., relativistic corrections, new potentials):

1. **Check `papers/`** to see if the derivation exists
2. **If code contradicts paper** → Flag it and ask user
3. **If we change physics in code** → Prompt user to update papers/README

**Authoritative source:** The core papers (`papers/core/`) — Papers 0, 1, 6, 7, 10, 11, 12, and 13 — are the authoritative theoretical source. If the README or any documentation conflicts with these papers, the papers take precedence. Flag the conflict to the user rather than silently resolving it in favor of the README.

**Example:**
```
You: "I want to add nuclear recoil corrections"
AI: "Let me check papers/Paper_4_Universality.tex for the formula..."
AI: "Found! Section 3.2 has the derivation. Implementing now."
```

### **2. The "Benchmarking" Rule**

After **any** modification to `hamiltonian.py`, `lattice.py`, or `solver.py`:

1. Run topological integrity: `pytest tests/test_fock_projection.py tests/test_fock_laplacian.py -v`
2. Run validation: `pytest tests/advanced_benchmarks.py`
3. Verify **18/18 symbolic proofs pass** (The Topological Foundation)
4. Verify **H2+** < 0.1% error (The Topological Control)
5. Verify **H2 Full CI** < 1.0% error (The Accuracy Control)
6. Report any speed regression > 10%

**Rationale:** Symbolic proofs validate the mathematical foundation. H2+ validates topology. H2 Full CI validates physics.

### **3. The "Clean Room" Rule**

- Generated plots → `debug/plots/` or `papers/figures/`
- Generated data → `debug/data/`
- Documentation → `docs/`
- **Never leave scripts in root**

If you create `test_something.py` in root during development:
```bash
mv test_something.py debug/
```

---

## 🧠 Memory Bank (Do Not Forget)

### **Core Physics**
- **Universal Constant:** Kinetic scale is `-1/16` exactly (validated: H, He+, H2+, muonic H)
- **H2 Correlation:** ~17% mean-field error is **physical** (missing correlation), not a bug
- **Full CI:** Requires `method='full_ci'` with `V_ee` (repulsion) and `V_en` (cross-attraction)
- **Basis Limit:** `max_n=10` is our "High Res" standard (~20s runtime)
- **Mass-Independence:** All holographic properties (d_s, c, topology) are mass-invariant

### **Validation Benchmarks** (Critical Tests)
- **Topological Integrity:** 18/18 symbolic proofs must pass → validates S³ conformal geometry (FOUNDATIONAL)
- **H2+ (Single-electron):** Must be < 0.1% error → validates topology
- **He (Hyperspherical):** Must be < 0.1% error → validates multi-electron natural geometry
- **H2 Full CI:** Must be < 1.0% error → validates physics
- **H2 Neumann V_ee:** 92.4% D_e → validates algebraic integral approach
- **Hyperspherical tests:** 20/20 must pass → validates angular + adiabatic + radial solvers
- **Muonic H:** Energy ratio = 206.77, topology identical → validates mass-independence
- **Spectral Dimension:** d_s ≈ 1.8-2.0, mass-independent → validates holography
- **Central Charge:** c ≈ 0.057, mass-independent → validates CFT connection

### **Theory Papers**
- **Core (`papers/core/`):**
  - **Paper 0:** Geometric packing framework and universal constant K = -1/16
  - **Paper 1:** Spectral graph theory foundations and eigenvalue methods
  - **Paper 7:** The Dimensionless Vacuum — formal proof that the graph Laplacian is a scale-invariant unit S³ topology and the Schrodinger equation is its flat-space projection (18/18 symbolic proofs)
  - **Paper 10:** Nuclear lattice — rovibrational spectra from SU(2) algebraic chains
  - **Paper 11:** Molecular Fock projection — prolate spheroidal lattice for diatomics (H₂⁺ at 0.70%)
  - **Paper 12:** Algebraic two-electron integrals — Neumann expansion for V_ee (H₂ at 92.4% D_e); cusp diagnosis
  - **Paper 13:** Hyperspherical lattice — two-electron atoms as coupled channel graphs (He at 0.05%); first non-trivial fiber bundle; ab initio molecular spectroscopy
- **Conjectures (`papers/conjectures/`):**
  - **Paper 2:** Fine structure constant (α⁻¹) derivation (geometric ansatz)
  - **Paper 3:** Holographic entropy, spectral dimension, central charge
  - **Paper 4:** Mass-independence, universality, muonic hydrogen
  - **Paper 5:** Comprehensive geometric vacuum framework (synthesis)

### **Old Research Insights** (Reference Only)
- **Fine Structure:** Old method used symplectic plaquettes (0.15% error) vs. our graph impedance (exploratory)
- **Proton Radius:** Old method optimized contact factors (80% match) vs. our theoretical prediction (25% match)
- **Note:** Different physical approaches - use for insights, not direct implementation

### **AdS/CFT Research** (Retained & Tested)
- **Boundary theory:** Graph/CFT - quantum states as nodes, adjacency matrix
- **Bulk theory:** AdS/geometric - 3D coordinates (x,y,z), symplectic areas, gauge fields
- **Holographic tools:** Muonic hydrogen, holographic entropy, spectral dimension, central charge (in `ADSCFT/`)
- **Goal:** Bridge from graph topology to geometric embedding for:
  - Fine structure constant (symplectic plaquettes — 0.0045% error achieved)
  - Proton radius (3D contact geometry — 100% agreement achieved)
- **Status:** Retained in `ADSCFT/`, tested via `tests/advanced_benchmarks.py`
- **Reference:** Old research used ParaboloidLattice with 3D coordinates (0.15% error on α)

### **Known Limitations**
- Mean-field only (no electron correlation beyond Full CI for H2)
- Spectral dimension plateau detection needs improvement
- Fine structure extraction is exploratory (algorithm incomplete)
- Proton radius prediction needs full hyperfine calculation
- AdS/CFT correspondence incomplete (boundary implemented, bulk in development)

---

## 📚 Documentation Workflow

### **When Creating Validation Reports**
1. Run tests, collect data → `debug/data/results.txt`
2. Generate plots → `debug/plots/validation_plot.png`
3. Write report → `docs/VALIDATION_REPORT_v0.3.4.md`
4. Update `CHANGELOG.md` with key findings
5. If release-worthy → Create `docs/releases/RELEASE_NOTES_v0.3.4.md`

### **When Creating New Demos**
1. Write polished script in `demo/demo_feature.py`
2. Include docstring with usage example
3. Add to `README.md` under "Examples"
4. Ensure minimal dependencies (numpy, scipy only if possible)

### **When Debugging Physics**
1. Create `debug/debug_issue_name.py` with descriptive name
2. Save outputs to `debug/data/` or `debug/plots/`
3. Document findings in `debug/FINDINGS.md`
4. If resolved → Move summary to `docs/archive/`

---

## 🚀 Common Commands

### **Testing & Validation**
```bash
# Run all tests
pytest tests/

# Run topological integrity proofs (MUST pass before any release)
pytest tests/test_fock_projection.py tests/test_fock_laplacian.py -v

# Run advanced benchmarks (muonic H, holography, etc.)
pytest tests/advanced_benchmarks.py -v

# Validate universal constant
python debug/validate_universal_constant.py

# Check version validation
python debug/validate_v0.3.2.py
```

### **Demos**
```bash
# H2 molecule demo
python demo/demo_h2.py

# Dirac relativistic demo
python demo/demo_h2_dirac.py
```

### **Development**
```bash
# Check theory papers
ls papers/*.tex
evince papers/Paper_4_Universality.pdf  # (or your PDF viewer)

# Review current validation status
cat docs/COMPLETE_VALIDATION_REPORT_v0.3.3.md

# Check performance benchmarks
cat benchmarks/BENCHMARKS.md
```

### **Package Management**
```bash
# Install package (editable mode)
pip install -e .

# Build package
python setup.py sdist bdist_wheel

# Test installation
python debug/test_install.py
```

---

## 🤖 AI Assistant Guidelines

### **Before Making Code Changes**
1. ✅ Check `papers/` for theoretical basis
2. ✅ Review `docs/COMPLETE_VALIDATION_REPORT_v*.md` for current status
3. ✅ Check `CHANGELOG.md` for recent changes
4. ✅ **Ask user** if unsure about physics interpretation

**Authoritative source:** The core papers (`papers/core/`) — Papers 0, 1, 6, 7, 10, 11, 12, and 13 — take precedence over README and all other documentation. If a conflict is found, flag it to the user rather than silently resolving it in favor of the README.

### **When Writing New Code**
1. Read relevant sections from `papers/Paper_*.tex`
2. Follow existing patterns in `geovac/` modules
3. Add type hints to **all** function signatures
4. Use sparse matrices for N > 100 (always)
5. Validate against benchmarks before finishing

### **When Generating Files**
1. **Scripts** → Save to correct directory (`debug/`, `demo/`, `tests/`)
2. **Documentation** → Save to `docs/` (not root!)
3. **Data/Plots** → Save to `debug/data/` or `debug/plots/`
4. **Use descriptive names** (e.g., `VALIDATION_REPORT_v0.3.3.md`, not `report.md`)
5. **Include metadata** at top: date, status, version

### **Documentation Standards**
- Include date and status at top of document
- Use clear section headers with emoji (optional but nice)
- Link to source files when referencing code
- Use tables for comparative data
- Include "Next Steps" section in reports

### **When Cleaning Up**
1. **Never delete** files without user confirmation
2. **Move, don't delete** → Use `mv` not `rm`
3. **Preserve git history** → Use `git mv` for tracked files
4. **Update links** in README and docs after moving files
5. **Test after cleanup** → Run `pytest tests/` to ensure nothing broke

---

## 🎯 Quick Reference

### **"Where does this file go?"**
| File Type | Location | Example |
|:---|:---|:---|
| Core module | `geovac/` | `hamiltonian.py` |
| Unit test | `tests/` | `test_lattice.py` |
| Demo script | `demo/` | `demo_h2.py` |
| Debug script | `debug/` | `debug_bridge.py` |
| Validation script | `debug/` | `validate_*.py` |
| Generated plot | `debug/plots/` | `convergence.png` |
| Generated data | `debug/data/` | `results.txt` |
| Documentation | `docs/` | `VALIDATION_REPORT.md` |
| Release notes | `docs/releases/` | `RELEASE_NOTES_v0.3.3.md` |
| Theory paper | `papers/` | `Paper_4.tex` |

### **"What's the error tolerance?"**
| Test | Max Error | Purpose |
|:---|:---:|:---|
| Symbolic proofs (18 tests) | 0 failures | **Topological foundation** |
| H (hydrogen) | < 0.1% | Basic validation |
| He+ (helium ion) | < 0.1% | Z-scaling check |
| H2+ (ionized H2) | < 0.1% | **Topological control** |
| He (hyperspherical) | < 0.1% | **Multi-electron control** |
| H2 Full CI | < 1.0% | **Accuracy control** |
| H2 Neumann V_ee | 92.4% D_e | Algebraic integral accuracy |
| Muonic H energy ratio | < 0.01% | Mass-independence |
| Speed regression | < 10% | Performance control |
| V_ee S³ overlap (1s-1s, 1s-2s, 2s-2s) | < 0.01% | Topological integrity |
| Direct CI vs matrix CI (any system) | < 1e-8 Ha | Algorithmic consistency |
| Hyperspherical He (20 tests) | 0 failures | Angular + adiabatic + radial |
| LiH binding energy (CP-corrected) | report only | BSSE >> D_e at nmax=3 |

### **"Which paper has the physics I need?"**
| Topic | Paper | Location | Tier |
|:---|:---:|:---|:---:|
| Universal constant -1/16 | Paper 0 | Sec 2 | Core |
| Graph Laplacian method | Paper 1 | Sec 3 | Core |
| O(V) quantum dynamics | Paper 6 | All | Core |
| Rabi oscillations | Paper 6 | — | Core |
| Delta-kick spectroscopy | Paper 6 | — | Core |
| AIMD / Langevin thermostat | Paper 6 | — | Core |
| Fine structure α⁻¹ | Paper 2 | Sec 4-6 | Conjecture |
| Spectral dimension d_s | Paper 3 | Sec 4 | Conjecture |
| Holographic entropy S | Paper 3 | Sec 5 | Conjecture |
| Central charge c | Paper 3 | Sec 6 | Conjecture |
| Mass-independence | Paper 4 | Sec 3-4 | Conjecture |
| Muonic hydrogen | Paper 4 | Sec 5 | Conjecture |
| Contact geometry | Paper 4 | Sec 5 | Conjecture |
| Comprehensive framework | Paper 5 | All | Conjecture |
| **Dimensionless vacuum** | **Paper 7** | **All** | **Core** |
| Schrodinger recovery proof | Paper 7 | Sec 4 | Core |
| S³ conformal geometry | Paper 7 | Sec 3 | Core |
| V_ee on S³ (node overlap) | Paper 7 | Sec VI | Core |
| Slater F⁰ master formula | Paper 7 | Sec VI.B | Core |
| Node vs edge V_ee warning | Paper 7 | Sec VI.C | Core |
| Excitation-driven Direct CI | K&H 1984 | `direct_ci.py` | Core |
| Neumann V_ee expansion | Paper 12 | Sec III-V | Core |
| Prolate spheroidal CI (H₂) | Paper 12 | Sec VI | Core |
| Cusp diagnosis (7.6% gap) | Paper 12 | Sec VII | Core |
| **Hyperspherical coordinates** | **Paper 13** | **Sec II** | **Core** |
| Angular eigenvalue (Gaunt) | Paper 13 | Sec III | Core |
| Adiabatic potential curves | Paper 13 | Sec IV | Core |
| Fiber bundle structure | Paper 13 | Sec VII | Core |
| Natural geometry hierarchy | Paper 13 | Sec VIII | Core |
| Ab initio spectroscopy | Paper 13 | Sec IX | Core |

---

## 📖 Version History

- **v2.5** (Mar 15, 2026): Added Papers 12, 13 to reference table; hyperspherical solver, Neumann V_ee, fiber bundle, ab initio spectroscopy; qubit encoding benchmarks
- **v2.4** (Mar 7, 2026): Corrected sparse/dense rule to context-dependent; direct CI validated as O(N_SD × N_connected)
- **v2.3** (Mar 6, 2026): Added excitation-driven Direct CI (Knowles-Handy 1984) references and algorithmic consistency tolerance
- **v2.2** (Mar 6, 2026): Added V_ee S³ density-overlap (Paper 7 Sec VI) to theory reference table
- **v2.1** (Feb 27, 2026): Clarified mathematical equivalence vs. physical priority framing; added Paper 6 to core reference table; established core papers as authoritative source over README
- **v2.0** (Feb 22, 2026): Added Dimensionless Vacuum Principle, Paper 7, topological integrity tests
- **v1.0** (Feb 14, 2026): Complete rewrite with directory structure, file naming, workflows
- **v0.1** (Feb 11, 2026): Initial version (basic guidelines only)

---

**Last Updated:** February 27, 2026
**Status:** Active (official project guidelines)
**Compliance:** All current files organized per these standards
