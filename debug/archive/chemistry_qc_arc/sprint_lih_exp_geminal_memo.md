# LiH π/M=16 marriage: Slater geminal e^{−γr} + the linexp c-offset (2026-09-26)

Reference throughout: M=16 π natural-orbital CI truncated at |c|>1e-2 (9 dets, nact 64,
E_trunc −8.025184). Analytic path: the O(nnz) sparse marriage (`lih_gvee_sparse.all_pieces_sparse`),
EXACT-Neumann path. VMC: the saved multi-geminal linear-method samples (`mg_lm_s1..8`, guide
|Ψ_CI(1+c(F−F̄))|² with the linexp0.5 analytic c) plus two new component-resolved VMC runs.

## How γ was set (geovac/ untouched)

`debug/lih_exp_geminal_run.py` reads `geovac/lih_r12ci/kernels.py`'s source, substitutes the
`GAM = 0.5` line (and optionally the `NXI, NETA, NPHI` grid line), and execs it as the module
`geovac.lih_r12ci.kernels` **before** any consumer is imported. A post-import monkeypatch would
not reach the consumers: gVee, hT, gT, triangle, hVee, phase2 and lih_pi_m16 all bind `GAM` by
`from .kernels import GAM` at import, and gVee builds `Kf`/`_KS` at import. The driver then
asserts that every namespace sees the new γ (KG, P2, M, GV, HT, GT, TR, HV, plus the `GV._KS`
keys). The same mechanism overrides triangle.py's `MMAX`/`LMAX`. The hard-coded φ-quadratures
(`build_kernel_gen` nphi 28, mode/smooth kernels 48) and the triangle rank (30) are patched
through `__defaults__` or wrappers.

**Sanity gate PASS:** linexp γ=0.5 through the patched path gives E_R12 = −8.042181268…,
c = +0.221588, and every piece identical to the banked `lih_pi_m16_energy_trunc1e-2.npz` values.
The patch changes nothing else.

## Job 1: exp geminal, γ scan (trunc 1e-2)

| γ | E_R12 (Ha) | corr vs E_trunc (mHa) | c | F̄ | σ² | cusp est. |
|---|---|---|---|---|---|---|
| 0.5 | −8.039002 | −13.82 | −0.22225 | 1.77882 | 0.16379 | 0.13 |
| 1.0 | −8.046711 | −21.53 | −0.24731 | 0.79490 | 0.07735 | 0.31 |
| **1.5** | **−8.048138** | **−22.95** | −0.23887 | 0.45722 | 0.04739 | **0.46** |
| 2.0 | −8.047573 | −22.39 | −0.22869 | 0.30183 | 0.03511 | 0.58 |
| linexp 0.5 (ref) | −8.042181 | −17.00 | +0.22159 | 3.49790 | 0.19337 | 0.25 |

Pieces:

| γ | h_T | h_Vne | h_Vee | g_T | g_Vne | g_Vee |
|---|---|---|---|---|---|---|
| 1.0 | +0.780522 | −0.988955 | +0.295473 | +0.989330 | −1.661178 | +0.324386 |
| 1.5 | +0.766392 | −0.930927 | +0.260626 | +0.827791 | −1.074195 | +0.220144 |
| 2.0 | +0.713469 | −0.848975 | +0.233404 | +0.771142 | −0.839108 | +0.178578 |

The γ trend is smooth, and **the minimum on the reference is near γ ≈ 1.5** (−8.04814), which
the VMC scan did not sample (it had 0.25/0.5/1/2).

**Cusp check.** For the multiplier 1 + c(F−F̄), the delivered coalescence cusp is
c·f′(0)/(1 + c(f(0) + ⟨F_rest⟩ − F̄)), estimated with ⟨F_rest⟩ ≈ (5/6)F̄ (the five spectator
pairs at their mean). Kato's antiparallel value is ½. The results:

- exp: 0.13 at γ=0.5, 0.31 at γ=1, **0.46 at γ=1.5**, 0.58 at γ=2.
- linexp: 0.20 at γ=0.25, 0.25 at γ=0.5, 0.21 at γ=1.

The energy optimum sits where the delivered cusp crosses Kato's ½. That fits the multi-geminal
finding that the gap was geminal *shape* (cusp), not geminal count. The linexp form r·e^{−γr}
delivers only about half the cusp at any γ, because f(0)=0 forces the whole cusp through a
single small c.

**Numerical convergence at γ=2 (the short-range worst case):**

| setting | E_R12 | ΔE vs base (mHa) |
|---|---|---|
| base 72×44, nphi 28/48, MMAX 4, rank 30, LMAX 34 | −8.047573 | 0 |
| grid 96×60, nphi 40/64 | −8.047586 | −0.013 |
| MMAX 8, rank 60, LMAX 46 | −8.047573 | +0.0004 |

**Gate PASS** (< 0.1 mHa). linexp γ=0.5 moves by −0.004 mHa on the fine grid, and γ=0.25 by
−0.002 mHa with every knob raised.

**Gate vs VMC (same reference), ingredient by ingredient, block jackknife over 200 blocks:**

| geminal | E analytic | E VMC | ΔE (mHa) | σ | c analytic | c VMC | c σ |
|---|---|---|---|---|---|---|---|
| exp 0.5 | −8.039002 | −8.039039 ± 0.39 | −0.04 | −0.10 | −0.22225 | −0.22321 | −1.5 |
| exp 1.0 | −8.046711 | −8.046798 ± 0.30 | −0.09 | −0.29 | −0.24731 | −0.24728 | +0.1 |
| exp 2.0 | −8.047573 | −8.047812 ± 0.29 (lin) / −8.047773 ± 0.20 (direct) | −0.24 / −0.20 | −0.81 / −1.0 | −0.22869 | −0.22834 | +1.7 |

At exp γ=2, F̄, σ², h, g and D all agree within 1.3σ. On the fine grid the energy gap is
−0.23 mHa.

**JOB 1 VERDICT: GO.** At γ=2 the analytic value sits 0.20–0.24 mHa (0.8–1.0σ) above the
same-reference VMC, the γ trend is consistent, and grid and truncation errors are ≤ 0.013 mHa.
The analytic exp geminal at γ=1.5 is the best trunc-1e-2 marriage energy: **−8.048138**,
6.0 mHa below linexp and 0.85 mHa below the same-reference Padé VMC (−8.04729). It still sits
above the 2-geminal linear-method VMC (−8.0490).

## Job 2: the linexp c offset

### Step 1: estimator artefact? No.

Decomposing the saved samples (`lih_c_offset_analysis.py`) for linexp0.5:

| quantity | VMC | analytic | σ |
|---|---|---|---|
| F̄ | 3.49816 ± 0.00049 | 3.49790 | +0.5 |
| σ² | 0.19317 ± 0.00048 | 0.19337 | −0.4 |
| h (col) | −0.07677 ± 0.00035 | −0.07671 | −0.2 |
| g_tot | −1.20967 ± 0.0039 | −1.20898 | −0.2 |
| D = g − Eσ² | 0.34386 ± 0.00116 | 0.34616 | **−2.0** |
| c (gen. eigen) | 0.22327 ± 0.00044 | 0.22159 | **+3.85** |

- The right-eigenvector c equals h_col/(Eσ² − g). h_col and D fluctuate together, so their ratio
  has a smaller relative error (0.2%) than either input. That correlation is why a 2σ shift in D
  shows up as a 3.85σ shift in c.
- Jackknife bias correction on the full ratio changes c by < 1e-6, and a bootstrap
  (2000 resamples) matches. Finite-sample ratio bias is therefore ruled out.
- The plug-in c from the symmetrized 2×2 built from VMC elements is 0.2244 ± 0.0019 (+1.45σ):
  consistent, just noisier.
- Seed scatter: χ²/dof is 12.6/7. Inflating the error by √(χ²/dof) gives 3.4σ.
- An **independent replicate** (10 new seeds, `cvmc_s101..110`) gives c = 0.22287 ± 0.00037
  (+3.6σ) and D offset −0.0033 ± 0.0009.

The offset reproduces across runs, so it is not a statistical fluke.

### Step 2: analytic numerics? No.

linexp0.5 under each refinement:

| setting | E_R12 | c | g_Vee |
|---|---|---|---|
| base 72×44/28/48, MMAX 4, rank 30, LMAX 34 | −8.042181 | 0.221588 | 0.6464884 |
| grid 96×60, nphi 40/64 | −8.042185 | 0.221601 | 0.6464878 |
| MMAX 8, rank 60 | −8.042181 | 0.221586 | 0.6464906 |
| LMAX 46 | −8.042181 | 0.221588 | 0.6464884 |

linexp0.25 with every knob raised: g_Vee 0.932859 vs 0.932862 at base. LMAX changes nothing
at all, because the reference pair densities are finite Legendre polynomials in η, so the
Neumann series terminates. **Discretization is excluded as the cause at the 1e-5 level.**

### Step 3: which piece?

I ran new VMC recording every component (`lih_c_offset_vmc{,2}.py`; guide as before;
nburn 4000, nsweep 9000, 1000 walkers). v2 has 12 seeds and a basis {Ψ, F_lx0.25Ψ, F_lx0.5Ψ},
with both kinetic estimators.

| linexp0.25 | VMC | analytic | diff | σ |
|---|---|---|---|---|
| g_T (non-sym) | 2.6391 ± 0.0047 | 2.6342 | +0.0049 | +1.0 |
| g_T (sym) | 2.6321 ± 0.0031 | 2.6342 | −0.0021 | −0.7 |
| g_Vne | −4.9320 ± 0.0072 | −4.9273 | −0.0047 | −0.65 |
| **g_Vee** | **0.93964 ± 0.00115** | **0.93286** | **+0.00678** | **+5.9** |
| h_T / h_Vne / h_Vee | — | — | — | +0.6 / +0.6 / +0.6 |
| D | 0.76049 ± 0.00123 | 0.75357 | +0.00692 | +5.6 |
| c | 0.15976 ± 0.00023 | 0.16146 | −0.00170 | −7.4 |
| E | −8.044806 | −8.044831 | +0.025 mHa | +0.1 |

| linexp0.5 | VMC (v2) | VMC (v1) | analytic | σ (combined) |
|---|---|---|---|---|
| g_Vee | 0.64387 ± 0.0012 | 0.64443 ± 0.0012 | 0.64649 | −2.8 |
| D offset | −0.0020 ± 0.0009 | −0.0033 ± 0.0009 | — | — |

- **The whole D offset is g_Vee** in both linexp cases: +0.0068 against +0.0069 at γ=0.25, and
  −0.0023 against −0.0026 at γ=0.5.
- g_Vee's VMC estimator is the simplest one in the problem: the weighted mean of (F−F̄)²V_ee,
  with no kinetic operator and finite variance.
- **Substituting the VMC g_Vee into the analytic 2×2 reproduces the VMC c exactly:**
  linexp0.5 goes 0.22159 → 0.22308, and linexp0.25 goes 0.16146 → 0.16003.
- The energy shift from that substitution is small: −0.11 mHa (γ=0.5) and +0.17 mHa (γ=0.25).
  That is why the energies agreed while c did not.

Where inside g_Vee? The analytic g_Vee is assembled as F2Vee − 2F̄·FVee + F̄²·Vee
(`lih_gvee_sparse.py:492`), a cancellation of about 200× at γ=0.25 (182.15 → 0.93). The offset
is a relative 4e-5 inconsistency among the raw moments, so the raw moments individually agree
with VMC (F2Vee 182.1453 ± 0.033 vs 182.1481; FVee and Vee within 0.6σ;
`lih_c_offset_vee_split.py`).

- Vee is not the culprit: it is consistent with E0_grid = E_trunc to 6e-8.
- FVee is disfavoured. An FVee error ε shifts h_Vee by ε and g_Vee by −2F̄ε. Producing +0.0068
  at γ=0.25 needs ε ≈ −4.6e-4, but h_Vee agrees with VMC (+3.0e-4 ± 4.7e-4).
- **Most likely locus: F2Vee (the two-f-edge + Coulomb classes: triangle, bridge, TREE4,
  pair).** The error is structural, not discretization. Its sign flips with γ: analytic low by
  0.0068 at γ=0.25, high by 0.0023 at 0.5, and within 2σ at 1.0 and 2.0. None of the exp
  geminals shows an offset above 1.2σ in D, apart from exp0.5 at −2.9σ with naive errors.
- The 2-MO precedent fits: that engine's g_Vee was +0.53834 against VMC +0.53605 (+0.0023),
  accepted at the time as "0.4%, within MC".

**JOB 2 VERDICT: a real analytic defect, localized to the piece g_Vee (most likely its F2Vee
raw moment). It is not estimator bias and not discretization.** It is not fixed: the raw-moment
error is 4e-5 relative, far below what a class-by-class VMC split can resolve (raw class errors
are about 0.01–0.03), so pinning the class needs an analytic-side audit of the F2Vee enumerator
rather than more sampling. The suggested next step is to evaluate each F2Vee topology class
(PAIR / share-vertex tree / TREE4 / BRIDGE / TRI) by direct 12-D importance MC on the
linexp0.25 reference, where the offset is largest (5.9σ). Impact on published numbers: ≤ 0.2 mHa
on linexp energies. For exp geminals there is no significant signal: the γ=2 D offset is
+1.2σ, and swapping it moves E by < 0.1 mHa.

## Recommendation on trunc 1e-3

Run it, **exp geminal at γ = 1.5** (about 2.5 h, about 9 GB, one run).

- At 1e-2, exp1.5 beats linexp0.5 by 6.0 mHa.
- The richer reference absorbed about 2 mHa of the linexp gain (corr −17.0 → −14.9 mHa). A
  similar partial absorption would put exp1.5 at 1e-3 near −8.048 to −8.049. That is below the
  full-CI Padé VMC (−8.047) and would make it the corpus-best rigorous LiH energy.
- If budget allows only one run, γ=1.5. γ=2 is the cross-check that has a VMC twin.
- Carry a ±0.2 mHa systematic for the open g_Vee defect.

## Files

- Drivers:
  - `debug/lih_exp_geminal_run.py` (γ, grid, nphi, MMAX, LMAX and rank overrides; propagation asserts)
  - `debug/lih_exp_geminal_compare.py`
  - `debug/lih_c_offset_{analysis,D,perseed_old,seedchi2,components,components2,vee_split}.py`
  - `debug/lih_c_offset_vmc{,2}.py`
- Data:
  - `debug/data/lih_exp_geminal_*_1e-02*.npz` (13 runs)
  - `debug/data/cvmc_s101..110.npz`, `debug/data/cvmc2_s201..212.npz`
  - logs in `debug/data/expgem/`
- **Patch for the PM:** none. The defect is localized but not fixed.
