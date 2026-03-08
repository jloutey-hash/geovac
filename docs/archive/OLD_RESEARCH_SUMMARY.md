# Old Research Archive - Quick Summary

**Date:** February 13, 2026
**Purpose:** Quick reference for insights from previous implementations

---

## 🎯 Key Findings

### Fine Structure Constant (α⁻¹ = 137.036)

**Old Implementation: 0.15% error ✓**

```python
# Method: Symplectic plaquette areas (NOT graph impedance)
S_n = sum of geometric plaquette areas in 3D lattice
P_n = √[(2πn)² + δ²]  # Helical photon path!
κ_n = S_n / P_n

# For n=5:
S_5 = 4325.83
P_5 = 31.567  (with helical pitch δ = 3.081)
κ_5 = 137.04  ✓
```

**Current Implementation: 96% error ❌**

```python
# Method: Graph Laplacian effective resistance
Z_em = mean(L†_ii + L†_jj) for dipole transitions
Z_em = 5.03  ❌
```

**Conclusion:** Different physics! Old method uses geometric areas, new method uses graph resistance. Keep current method as "exploratory" or implement symplectic approach.

---

### Proton Radius Puzzle (Δr_p = 0.034 fm)

**Old Implementation: 80% match ✓**

```python
# Method: Optimize contact factors from experimental data
C_e = 0.6658  (fitted to E_HFS = 5.87 μeV)
C_μ = 0.5000  (fitted to E_HFS = 182.7 meV)
Δr_p = 0.043 fm  ✓
```

**Current Implementation: 25% match ❌**

```python
# Method: Assume contact factors from theory
C_e = 2/3  (theoretical)
C_μ = 1/2  (theoretical)
Δr_p = 0.011 fm  ❌ (4× too small)
```

**Conclusion:** Need to implement full hyperfine calculation with contact factor optimization (not just simple scaling).

---

## 📁 Key Files

### Fine Structure
- `old_research_archive/src/hydrogen_u1_impedance.py`
- `old_research_archive/src/compute_alpha.py`

### Proton Radius
- `old_research_archive/MUONIC_HYDROGEN_REPORT.md`
- `old_research_archive/src/muonic_hydrogen_analysis.py`

---

## 🔧 Recommended Fixes

### For Fine Structure (v0.4.0)
- **Option 1:** Implement symplectic plaquette method (major effort)
- **Option 2:** Keep as "exploratory" (current status is correct)

### For Proton Radius (v0.3.4)
- **Implement:** Full hyperfine calculator with contact factor optimization
- **Expected improvement:** 25% → 80% match

---

See [INSIGHTS_FROM_OLD_RESEARCH.md](INSIGHTS_FROM_OLD_RESEARCH.md) for full technical details.
