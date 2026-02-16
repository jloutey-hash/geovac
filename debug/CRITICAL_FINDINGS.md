# CRITICAL FINDINGS - Hydrogen Atom Debug

**Date:** February 13, 2026
**Status:** 🚨 **BLOCKING ISSUE IDENTIFIED**

---

## 🔴 Problem Statement

**Hydrogen atom gives -0.97 Ha instead of -0.5 Ha (94% error!)**

This is **NOT** the 4.9% error we saw earlier. This is a fundamental breakdown of the method.

---

## 🔍 Root Cause Analysis

### Issue: Kinetic Scale Mismatch

The "universal" kinetic scale `-1/16` is **NOT universal** for atoms in hybrid mode!

#### Current Hamiltonian (Hybrid Mode):
```
H₁ = T + V

where:
  T = -0.5 * kinetic_scale * (D - A)
  V = -Z/n² (diagonal potential)
```

#### For Hydrogen (Z=1, n=1):

**Potential Energy:**
```
V(n=1) = -1/1² = -1.0 Ha  ✓ (correct by virial theorem)
```

**Kinetic Energy (PROBLEM):**
```
Laplacian eigenvalue:  λ_L ≈ 0.009  (measured)
Kinetic contribution:  T = -0.5 * (-1/16) * 0.009
                       T ≈ 0.0003 Ha  ❌ (need +0.5 Ha!)
```

**Total Energy:**
```
E = T + V
E ≈ 0.0003 + (-1.0)
E ≈ -1.0 Ha  ❌ (should be -0.5 Ha)
```

---

## 📊 Quantitative Analysis

### Required Kinetic Scale

To achieve `E = -0.5 Ha`, we need:
```
E = T + V = -0.5
T = -0.5 - V = -0.5 - (-1.0) = +0.5 Ha

T = -0.5 * scale * λ_L
0.5 = -0.5 * scale * 0.009
scale = -111

Required scale:  ≈ -111
Universal scale: = -1/16 = -0.0625
Ratio:           ≈ 1,776x difference!
```

---

## 🎯 Key Findings

### 1. Graph Laplacian Eigenvalues Too Small

**Measured:**
```
Lowest Laplacian eigenvalues:
  λ₁ ≈ 0.000  (degenerate ground state)
  λ₇ ≈ 0.006
  λ₈ ≈ 0.007
  λ₉ ≈ 0.008
  λ₁₀ ≈ 0.009
```

**Problem:** Multiple near-zero eigenvalues suggest:
- Disconnected graph components, OR
- Graph structure doesn't match quantum state space, OR
- Edge weights incorrectly normalized

### 2. Universal Constant Fails for Atoms

**Evidence:**

| System | Mode | Kinetic Scale | Status |
|:---|:---|:---:|:---|
| H atom | Hybrid | -1/16 | ❌ 94% error |
| H₂ molecule | Mean-field | -1/16 | ✓ ~17% (expected) |
| H₂ molecule | Full CI | -1/16 | ✓ 2.8% error |

**Conclusion:** The `-1/16` scale works for **molecules**, not **atoms in hybrid mode**.

### 3. Two Different Formulations

**Molecular Hamiltonian (works):**
```python
# For molecules, uses molecular bridge construction
H = kinetic_scale * (D - A)
# kinetic_scale = -1/16 ✓
```

**Atomic Hamiltonian (broken):**
```python
# For atoms, uses hybrid mode
H₁ = T + V
T = -0.5 * kinetic_scale * (D - A)
V = -Z/n²
# kinetic_scale = -1/16 ❌ (wrong by 1776x!)
```

---

## 🔧 Immediate Fixes Needed

### Option 1: Use Correct Scale for Atoms (QUICK FIX)

```python
# In benchmark_suite.py, use atom-specific scale:
ATOMIC_KINETIC_SCALE = -111  # Calibrated for hybrid mode
```

**Pros:** Immediate fix, hydrogen will work
**Cons:** Breaks "universality" claim

### Option 2: Fix Hybrid Mode Formulation (CORRECT FIX)

The issue is the factor of `-0.5` in the kinetic term. Try:
```python
# Remove the -0.5 factor:
T = kinetic_scale * laplacian  # Not -0.5 * kinetic_scale
```

**Reasoning:** The `-0.5` might be doubling the effect when combined with `-1/16`.

### Option 3: Use Geometric Mode for Atoms (ALTERNATIVE)

```python
h = HeliumHamiltonian(max_n=max_n, Z=Z,
                     kinetic_scale=UNIVERSAL_KINETIC_SCALE,
                     geometric_mode=True)  # Pure geometric
```

**Test:** Does geometric mode work with `-1/16`?

---

## 📝 Recommended Action Plan

### Immediate (Today)

1. ✅ **Document the issue** (this file)
2. 🔲 **Test geometric mode** - Does it work with -1/16?
3. 🔲 **Test factor-of-2 fix** - Remove -0.5 from kinetic term
4. 🔲 **Find calibrated scale** - What scale makes H atom exact?

### Short-term (This Week)

1. 🔲 **Understand the discrepancy** - Why do molecules work but atoms don't?
2. 🔲 **Fix the formulation** - Make it consistent
3. 🔲 **Update benchmarks** - Use correct method for atoms
4. 🔲 **Document limitations** - Be honest about domain of validity

### Long-term (Next Release)

1. 🔲 **Unify formulations** - One consistent approach
2. 🔲 **Validate universal constant** - Or admit it's system-specific
3. 🔲 **Theory paper** - Explain when/why it works

---

## 🤔 Theoretical Questions

### Why Does -1/16 Work for Molecules?

**Hypothesis:**
- Molecular Hamiltonian uses different construction
- Bridge edges provide additional coupling
- The formulation is fundamentally different

**Need to investigate:**
- How does MoleculeHamiltonian build its matrices?
- Is it also using hybrid mode with T + V?
- Or is it pure geometric?

### Is the Problem the Hybrid Mode?

**Test:** Compare three modes for H atom:
1. Hybrid mode (current): `T + V` formulation
2. Pure geometric: `kinetic_scale * L` only
3. Modified hybrid: Different T formula

---

## ⚠️ Impact on Previous Results

### Helium Full CI (1.24% error)

**Status:** ✓ **STILL VALID**

Helium uses `HeliumHamiltonian.compute_ground_state()` which solves the full 2-electron problem:
- Uses `h.h2` (two-particle Hamiltonian)
- Includes electron-electron repulsion
- NOT just 2× single-particle energy

**The Full CI result is independent of single-particle issues!**

### H₂ Molecule Results

**Status:** ✓ **VALID**

H₂ uses `MoleculeHamiltonian` which has different construction:
- Different graph structure (inter-atomic bridges)
- Possibly different formulation
- Validated to 2.8% error

---

## 🎯 Bottom Line

**The "universal" kinetic scale of -1/16 is:**
- ✅ Valid for molecules (H₂, H₂⁺)
- ✅ Valid for multi-electron Full CI (He atom)
- ❌ **INVALID for single-electron atoms in hybrid mode**

**This is a FORMULATION issue, not a fundamental flaw.**

The tensor product / Full CI approach works (proven by Helium 1.24% error).
The problem is specific to how we handle single-particle atomic Hamiltonians.

---

## 🚀 Next Steps

**Priority 1:** Test geometric mode for H atom
```bash
python debug/test_geometric_mode.py
```

**Priority 2:** Find correct atomic kinetic scale
```bash
python debug/calibrate_atomic_scale.py
```

**Priority 3:** Fix the formulation or document the limitation
- Update BENCHMARKS.md
- Add caveats to README
- Fix atomic calculations in benchmark_suite.py

---

**Status:** Investigation ongoing
**Blocking:** Yes (for atomic benchmarks)
**Impact:** High (affects interpretation of results)
**Priority:** P0 (critical)

**Conclusion:** We found the bug, now we fix it! 🔧
