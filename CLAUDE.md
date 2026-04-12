# GeoVac Project Guidelines

## 🌍 Project Context
**Name:** GeoVac (The Geometric Vacuum)
**Mission:** Topological Quantum Chemistry Engine. We solve the wave function using Spectral Graph Theory (Sparse Graph Laplacians) rather than continuous integration.
**Theory Foundation:** Papers 0-5 in `papers/` directory (LaTeX source of truth)

---

## 📂 Directory Structure Rules

**ROOT DIRECTORY MUST REMAIN CLEAN.** Only essential files allowed.

### **Core Directories**
- `geovac/` → Core package source code (`__init__.py`, `hamiltonian.py`, `lattice.py`, etc.)
- `papers/` → **Theory Source of Truth.** LaTeX papers (Paper_0 through Paper_5). Read here first when unsure about physics.
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

### **Theoretical Research** (Exploratory - May Not Be in Production)
- `ADSCFT/` → AdS/CFT correspondence and holographic duality research
  - `ADSCFT/boundary/` → Boundary theory (CFT/graph implementation)
  - `ADSCFT/bulk/` → Bulk theory (AdS/geometric embedding with 3D coordinates)
  - `ADSCFT/tests/` → Validation tests for AdS/CFT predictions
  - **Status:** Theoretical/exploratory work - less precise than core methods
  - **Note:** Isolated from main `geovac/` package. May not be included in releases.
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
- `*.tex`, `*.bib`, `*.aux`, `*.log` → Papers (`papers/`)

---

## ⚡ Coding Standards

### **1. Sparse Matrix First**
Always use `scipy.sparse` (csr_matrix, coo_matrix). **Never densify** matrices unless N < 100.

```python
# ✅ GOOD
from scipy.sparse import csr_matrix
H = csr_matrix((n_states, n_states))

# ❌ BAD
H = H.toarray()  # Only if n_states < 100
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

## 🧪 Workflow Protocols

### **1. The "Theory Check" Rule**

Before implementing new physics (e.g., relativistic corrections, new potentials):

1. **Check `papers/`** to see if the derivation exists
2. **If code contradicts paper** → Flag it and ask user
3. **If we change physics in code** → Prompt user to update papers/README

**Example:**
```
You: "I want to add nuclear recoil corrections"
AI: "Let me check papers/Paper_4_Universality.tex for the formula..."
AI: "Found! Section 3.2 has the derivation. Implementing now."
```

### **2. The "Benchmarking" Rule**

After **any** modification to `hamiltonian.py`, `lattice.py`, or `solver.py`:

1. Run validation: `pytest tests/advanced_benchmarks.py`
2. Verify **H2+** < 0.1% error (The Topological Control)
3. Verify **H2 Full CI** < 1.0% error (The Accuracy Control)
4. Report any speed regression > 10%

**Rationale:** H2+ (single-electron) validates topology. H2 Full CI validates physics.

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
- **H2+ (Single-electron):** Must be < 0.1% error → validates topology
- **H2 Full CI:** Must be < 1.0% error → validates physics
- **Muonic H:** Energy ratio = 206.77, topology identical → validates mass-independence
- **Spectral Dimension:** d_s ≈ 1.8-2.0, mass-independent → validates holography
- **Central Charge:** c ≈ 0.057, mass-independent → validates CFT connection

### **Theory Papers** (Source of Truth)
- **Paper 0:** Geometric packing framework and universal constant discovery
- **Paper 1:** Spectral graph theory foundations and eigenvalue methods
- **Paper 2:** Fine structure constant (α⁻¹) derivation (exploratory)
- **Paper 3:** Holographic entropy, spectral dimension, central charge
- **Paper 4:** Mass-independence, universality, muonic hydrogen
- **Paper 5:** Comprehensive geometric vacuum framework (synthesis)

### **Old Research Insights** (Reference Only)
- **Fine Structure:** Old method used symplectic plaquettes (0.15% error) vs. our graph impedance (exploratory)
- **Proton Radius:** Old method optimized contact factors (80% match) vs. our theoretical prediction (25% match)
- **Note:** Different physical approaches - use for insights, not direct implementation

### **AdS/CFT Research** (Theoretical/Exploratory)
- **Current:** Boundary theory (graph/CFT) - quantum states as nodes, adjacency matrix
- **Missing:** Bulk theory (AdS/geometric) - 3D coordinates (x,y,z), symplectic areas, gauge fields
- **Goal:** Bridge from graph topology to geometric embedding for:
  - Fine structure constant (needs symplectic plaquettes, not graph impedance)
  - Proton radius (needs 3D contact geometry, not simplified formula)
- **Status:** Isolated in `ADSCFT/` - may not be production-ready
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
| H (hydrogen) | < 0.1% | Basic validation |
| He+ (helium ion) | < 0.1% | Z-scaling check |
| H2+ (ionized H2) | < 0.1% | **Topological control** |
| H2 Full CI | < 1.0% | **Accuracy control** |
| Muonic H energy ratio | < 0.01% | Mass-independence |
| Speed regression | < 10% | Performance control |

### **"Which paper has the physics I need?"**
| Topic | Paper | Section |
|:---|:---:|:---|
| Universal constant -1/16 | Paper 0 | Sec 2 |
| Graph Laplacian method | Paper 1 | Sec 3 |
| Fine structure α⁻¹ | Paper 2 | Sec 4-6 |
| Spectral dimension d_s | Paper 3 | Sec 4 |
| Holographic entropy S | Paper 3 | Sec 5 |
| Central charge c | Paper 3 | Sec 6 |
| Mass-independence | Paper 4 | Sec 3-4 |
| Muonic hydrogen | Paper 4 | Sec 5 |
| Contact geometry | Paper 4 | Sec 5 |
| Comprehensive framework | Paper 5 | All |

---

## 📖 Version History

- **v1.0** (Feb 14, 2026): Complete rewrite with directory structure, file naming, workflows
- **v0.1** (Feb 11, 2026): Initial version (basic guidelines only)

---

**Last Updated:** February 14, 2026
**Status:** Active (official project guidelines)
**Compliance:** All current files organized per these standards
