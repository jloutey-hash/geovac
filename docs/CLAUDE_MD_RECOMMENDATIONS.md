# Recommendations for CLAUDE.md Updates

**Date:** February 14, 2026
**Status:** Proposed changes to improve project guidelines

---

## 📝 Recommended Updates

### **1. Add Missing Directory: `docs/`**

**Current CLAUDE.md** mentions:
- `geovac/`, `papers/`, `old_research_archive/`, `tests/`, `debug/`, `demo/`, `benchmarks/`

**Missing:**
- `docs/` - For project documentation (reports, analysis, release notes)

**Recommended addition:**
```markdown
- `docs/` → Project documentation, validation reports, release notes.
  - `docs/releases/` → Historical release notes
  - `docs/archive/` → Deprecated/historical documentation
```

---

### **2. Update Directory Structure Section**

**Current:**
```markdown
## 📂 Directory Structure Rules
**ROOT DIRECTORY MUST REMAIN CLEAN.** Do not create temp scripts here.

- `geovac/` → Core package source code (`__init__.py`, `hamiltonian.py`).
- `papers/` → **Theory Source of Truth.** If you are unsure about physics, read here first.
- `old_research_archive/` → Old code/theory. Use for reference, but do not import from here.
- `tests/` → Unit tests (pytest).
- `debug/` → Scratchpad for analyzing one-off physics issues (e.g., `debug_bridge_scaling.py`).
- `demo/` → Polished, customer-facing scripts (e.g., `demo_h2.py`).
- `benchmarks/` → Performance tracking scripts.
```

**Recommended:**
```markdown
## 📂 Directory Structure Rules
**ROOT DIRECTORY MUST REMAIN CLEAN.** Only essential files (README, LICENSE, setup.py, CLAUDE.md, CHANGELOG.md).

### **Core Directories**
- `geovac/` → Core package source code (`__init__.py`, `hamiltonian.py`, `lattice.py`, etc.).
- `papers/` → **Theory Source of Truth.** LaTeX source files (Paper_0-5). Read here first when unsure about physics.
- `tests/` → Unit tests and validation suites (pytest format).

### **Development Directories**
- `debug/` → Scratchpad for analyzing physics issues. Subdirectories:
  - `debug/plots/` → Generated plots and figures
  - `debug/data/` → Generated data files (.txt, .csv, .json)
- `demo/` → Polished, customer-facing demonstration scripts (e.g., `demo_h2.py`).
- `benchmarks/` → Performance tracking and regression tests.
  - `benchmarks/scripts/` → Benchmark suite scripts
  - `benchmarks/figures/` → Performance plots

### **Documentation**
- `docs/` → Project documentation, validation reports, analysis.
  - `docs/releases/` → Release notes and version summaries
  - `docs/archive/` → Historical/deprecated documentation
- `old_research_archive/` → Legacy code/theory. **Reference only** - do not import.

### **Files Allowed in Root**
✅ **Allowed:**
- README.md, LICENSE, CHANGELOG.md, CLAUDE.md
- setup.py, pyproject.toml, requirements.txt
- .gitignore, .github/

❌ **NOT Allowed:**
- .py scripts (move to debug/ or demo/)
- .png files (move to debug/plots/ or papers/figures/)
- .txt/.csv data files (move to debug/data/)
- Documentation files (move to docs/)
```

---

### **3. Add Section: File Naming Conventions**

**New section to add:**
```markdown
## 📛 File Naming Conventions

### **Python Scripts**
- `demo_*.py` → Demo scripts (goes in `demo/`)
- `validate_*.py` → Validation scripts (goes in `debug/`)
- `test_*.py` → Unit tests (goes in `tests/`)
- `benchmark_*.py` → Benchmarks (goes in `benchmarks/scripts/`)
- `debug_*.py` → Debug analysis (goes in `debug/`)

### **Documentation**
- `RELEASE_NOTES_v*.md` → Version release notes (goes in `docs/releases/`)
- `RELEASE_SUMMARY_v*.md` → Version summaries (goes in `docs/releases/`)
- `*_VALIDATION_REPORT.md` → Validation reports (goes in `docs/`)
- `*_ANALYSIS.md` → Technical analysis (goes in `docs/`)

### **Data Files**
- `*.png`, `*.pdf`, `*.svg` → Plots (goes in `debug/plots/` or `papers/figures/`)
- `*.txt`, `*.csv`, `*.json` → Data (goes in `debug/data/`)
- `*.tex`, `*.bib` → Papers (goes in `papers/`)
```

---

### **4. Update Memory Bank Section**

**Current:**
```markdown
## 🧠 Memory Bank (Do Not Forget)
- **Universal Constant:** The kinetic scale is `-1/16` exactly.
- **H2 Correlation:** The ~17% error in Mean Field is **Physical** (missing correlation), not a bug.
- **Full CI:** Requires `method='full_ci'` and includes `V_ee` (repulsion) and `V_en` (cross-attraction).
- **Basis Limit:** `max_n=10` is our current "High Res" standard (~20s runtime).
```

**Recommended additions:**
```markdown
## 🧠 Memory Bank (Do Not Forget)

### **Core Physics**
- **Universal Constant:** The kinetic scale is `-1/16` exactly (validated across H, He+, H2+).
- **H2 Correlation:** The ~17% error in Mean Field is **Physical** (missing correlation), not a bug.
- **Full CI:** Requires `method='full_ci'` and includes `V_ee` (repulsion) and `V_en` (cross-attraction).
- **Basis Limit:** `max_n=10` is our current "High Res" standard (~20s runtime).

### **Validation Benchmarks**
- **H2+ (Single-electron):** Must be < 0.1% error (topological control)
- **H2 Full CI:** Must be < 1.0% error (accuracy control)
- **Mass-independence:** Muonic/electronic ratio must be 1.0000 for all holographic properties

### **Theory Papers** (Source of Truth)
- **Paper 0:** Geometric packing and universal constant
- **Paper 1:** Spectral graph theory foundations
- **Paper 2:** Fine structure constant (exploratory)
- **Paper 3:** Holographic entropy and central charge
- **Paper 4:** Mass-independence and universality
- **Paper 5:** Geometric vacuum framework (comprehensive)

### **Old Research Archive** (Insights)
- **Fine structure:** Old method used symplectic plaquettes (0.15% error)
- **Proton radius:** Old method optimized contact factors (80% match)
- **Note:** Different physical approaches - use for reference, not direct implementation
```

---

### **5. Add Section: Documentation Workflow**

**New section to add:**
```markdown
## 📚 Documentation Workflow

### **When Creating Validation Reports**
1. Run tests and collect data → save to `debug/data/`
2. Generate plots → save to `debug/plots/`
3. Write report → save to `docs/` with descriptive name
4. Update CHANGELOG.md with key findings
5. If release-worthy, create `docs/releases/RELEASE_NOTES_v*.md`

### **When Creating New Demos**
1. Write polished script in `demo/`
2. Include docstring with usage example
3. Add to README.md under "Quick Start" or "Examples"
4. Ensure it runs successfully with minimal dependencies

### **When Debugging Physics**
1. Create script in `debug/` with descriptive name (e.g., `debug_bridge_scaling.py`)
2. Save outputs to `debug/data/` or `debug/plots/`
3. Document findings in `debug/FINDINGS.md` or similar
4. If resolved, consider moving findings to `docs/archive/`
```

---

### **6. Update Common Commands Section**

**Current:**
```markdown
## 🚀 Common Commands
- **Run Tests:** `pytest tests/`
- **Run Demo:** `python demo/demo_h2.py`
- **Validate Constant:** `python debug/validate_universal_constant.py`
```

**Recommended:**
```markdown
## 🚀 Common Commands

### **Testing & Validation**
- **Run Tests:** `pytest tests/`
- **Run Advanced Benchmarks:** `pytest tests/advanced_benchmarks.py -v`
- **Validate Universal Constant:** `python debug/validate_universal_constant.py`

### **Demos**
- **H2 Molecule:** `python demo/demo_h2.py`
- **Dirac Hamiltonian:** `python demo/demo_h2_dirac.py`

### **Development**
- **Check Theory:** Review `papers/Paper_*.tex` for physics equations
- **Review Validation:** See `docs/COMPLETE_VALIDATION_REPORT_v*.md`
- **Check Performance:** See `benchmarks/BENCHMARKS.md`

### **Package Management**
- **Install (editable):** `pip install -e .`
- **Build package:** `python setup.py sdist bdist_wheel`
- **Test install:** `python debug/test_install.py`
```

---

### **7. Add Section: AI Assistant Guidelines**

**New section to add:**
```markdown
## 🤖 AI Assistant Specific Guidelines

### **Before Making Changes**
1. **Check papers/** for theoretical basis
2. **Check docs/COMPLETE_VALIDATION_REPORT** for current status
3. **Check CHANGELOG.md** for recent changes
4. **Ask user** if unsure about physics interpretation

### **When Writing New Code**
1. **Read relevant theory** from papers/ first
2. **Follow existing patterns** in geovac/ modules
3. **Add type hints** to all function signatures
4. **Use sparse matrices** (scipy.sparse) for N > 100
5. **Validate against benchmarks** before committing

### **When Generating Documentation**
1. **Save to correct directory** (docs/, not root)
2. **Use clear, descriptive names** (e.g., `VALIDATION_REPORT_v0.3.3.md`)
3. **Include date and status** at top of document
4. **Link to source files** when referencing code

### **When Cleaning Up**
1. **Never delete** without user confirmation
2. **Move, don't delete** historical docs → `docs/archive/`
3. **Preserve git history** (git mv, not rm + add)
4. **Update links** in README and other docs after moving files
```

---

## 🎯 Priority Rankings

### **High Priority** (Implement immediately)
1. ✅ Add `docs/` to directory structure
2. ✅ Add file naming conventions
3. ✅ Expand allowed/not-allowed root files list

### **Medium Priority** (Nice to have)
4. ⚠ Add documentation workflow section
5. ⚠ Expand common commands
6. ⚠ Update memory bank with validation benchmarks

### **Low Priority** (Optional)
7. ℹ️ Add AI assistant guidelines (helpful but not critical)

---

## 📄 Proposed Updated CLAUDE.md

See attached: `CLAUDE_UPDATED.md` (complete rewrite with all recommendations)

---

**Status:** Ready for review
**Impact:** Improves clarity, completeness, and maintainability
**Breaking Changes:** None (purely additive)
