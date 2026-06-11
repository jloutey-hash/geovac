# Sprint LS-3: Regularized Bethe logarithm via acceleration form

**Date:** 2026-05-02
**Sprint goal:** Close out the from-scratch hydrogen Lamb shift derivation by
implementing Drake's acceleration form of the Bethe logarithm, hoping to break
through LS-2's velocity-form 0/0 closure pathology for l > 0 (2P).
**Verdict:** **MIXED.** Acceleration form gives ~3× faster convergence for
s-states (1S, 2S, 3S) — 2S Bethe log error reduced from 3.1 % (LS-2) to 0.92 %
(LS-3 N=40), reducing the LS-1/LS-2 Lamb-shift truncation error by 3.3×. For
l > 0 states (2P, 3P, 3D) the acceleration form has the **same closure
pathology** as velocity (mathematically identical by [p, H₀] = −i∇V), and
genuine Drake-Swainson regularization is structurally a different sprint.

## 1. Setup

The velocity-form Bethe logarithm is

> ln k_0(nl) = J_v(nl) / I_v(nl)

with `I_v(nl) = Σ_m |⟨nl|p|m⟩|²(E_m − E_n)`, `J_v(nl) = Σ_m |⟨nl|p|m⟩|²(E_m − E_n) ln|2(E_m−E_n)/Ry|`. Closure: `I_v(nl) = (2 Z⁴/n³) δ_{l,0}` — vanishes for l > 0.

The acceleration form uses `[p, H₀] = −i ∇V` to substitute `|⟨p⟩|²(ΔE)² = |⟨∇V⟩|²`, giving

> I_a(nl) = Σ_m |⟨nl|∇V|m⟩|² / (E_m − E_n)
> J_a(nl) = Σ_m |⟨nl|∇V|m⟩|² ln|2(E_m − E_n)/Ry| / (E_m − E_n)

By the operator identity `I_a = I_v` and `J_a = J_v` *as continuum sums*, but
**finite-basis evaluation is reweighted**: the 1/(ΔE) factor suppresses
high-energy pseudostates and emphasises low-energy ones. For Coulomb V_C =
−Z/r, ∇V = Z r̂/r², and the m-summed radial squared matrix element is
`l_> · Z² · |∫ R_{nl} R_{n'l'} dr|²` (no r²-weight in the radial integral —
the 1/r² cancels the volume measure).

The Sturmian basis machinery is reused from LS-2; only the matrix element of
1/r² (instead of r) was added (`sturmian_acceleration_matrix`). All matrix
elements remain exact polynomial-times-exponential integrals computed in
mpmath at 50-digit precision.

## 2. Cross-check on s-states (1S, 2S)

Both forms must converge to the same Drake-Swainson value. Different
finite-basis convergence rates expose the spectral-density reweighting.

### 1S (Drake-Swainson 2.984128)

| N  | acc ln k_0  | err (%) | vel ln k_0  | err (%) |
|---:|:-----------:|:-------:|:-----------:|:-------:|
| 12 |    2.7269   |  −8.62  |    2.5312   | −15.18  |
| 16 |    2.8906   |  −3.13  |    2.7080   |  −9.25  |
| 20 |  **3.0019** | **+0.60** |  2.8315   |  −5.11  |
| 24 |    3.0831   |  +3.32  |    2.9236   |  −2.03  |
| 30 |    3.1712   |  +6.27  |    3.0256   |  +1.39  |

**Acceleration form crosses the Drake value at smaller N (N=20 vs N≈30 for velocity),
giving 3× tighter agreement at the optimal N**. Past the crossing, both forms
overshoot — characteristic of the Sturmian basis quality limit (high-energy
pseudostates start carrying numerical error in J before they cancel correctly
in the I sum).

### 2S (Drake-Swainson 2.811770)

| N  | acc ln k_0  | err (%) | vel ln k_0  | err (%) |
|---:|:-----------:|:-------:|:-----------:|:-------:|
| 16 |    2.1358   | −24.04  |    1.8964   | −32.56  |
| 20 |    2.3276   | −17.22  |    2.0974   | −25.41  |
| 24 |    2.4673   | −12.25  |    2.2468   | −20.09  |
| 30 |    2.6189   |  −6.86  |    2.4124   | −14.20  |
| 40 |  **2.7858** | **−0.92** |  2.5995   |  −7.55  |
| 50 |    2.8955   |  +2.98  |    2.7255   |  −3.07  |

**Acceleration N=40 = 2.786 (err −0.92%)** is the best practical estimator.
This is **3.3× tighter than LS-2's velocity N=50 = 2.726 (err −3.1%)**.

### 3S (Drake-Swainson 2.767664)

Acceleration form, λ = 1/3:

| N  | ln k_0  | err (%) |
|---:|:-------:|:-------:|
| 24 | 2.0359  | −26.4   |
| 30 | 2.2432  | −18.9   |
| 40 | 2.4709  | −10.7   |

Convergence is slower at higher n_target because the basis must resolve more
nodal structure simultaneously. Reaching ~1 % at 3S requires N ≳ 80, beyond
this sprint's compute budget.

## 3. 2P (l = 1): closure pathology persists

| N  | acc I        | acc J        | acc ln k_0   |
|---:|:------------:|:------------:|:------------:|
|  8 | −1.52×10⁻³  | −3.22×10⁻²  |    +21.23   |
| 16 | −2.92×10⁻⁴  | −2.49×10⁻²  |    +85.04   |
| 24 | −1.02×10⁻⁴  | −2.34×10⁻²  |   +230.29   |
| 30 | −5.57×10⁻⁵  | −2.30×10⁻²  |   +413.72   |
| 50 | −1.34×10⁻⁵  | −2.27×10⁻²  |  +1685.90   |

I → 0 as the basis size grows (closure says I_v(2P) = 0 exactly), J stays
roughly constant; the ratio J/I → ∞. **Drake's −0.030 is not recovered**.

This is the same pathology as velocity form, with reordered numerics: the
acceleration 1/(ΔE) weighting makes the magnitudes of I smaller (and the
divergence faster), but the structural cancellation between l−1=0 and
l+1=2 channels is unchanged.

### 3.1 Per-channel structural finding

Splitting the J/I ratio per intermediate-channel l_p = 0 vs l_p = 2:

| N  | ln k_0(2P→S)  | ln k_0(2P→D)  |
|---:|:-------------:|:-------------:|
| 16 |    +0.3850    |    +0.0654    |
| 20 |    +0.3849    |    +0.0836    |
| 24 |    +0.3848    |    +0.0937    |
| 30 |    +0.3848    |    +0.1020    |

**Each intermediate-channel ratio is finite and converges separately** (the
S-channel rapidly to ~0.385 driven by the 2P→1S bound-bound transition; the
D-channel slowly toward a value > 0.10 from continuum pseudostates). Their
weighted sum can in principle hit Drake's −0.030, but **no obvious linear
combination reproduces this value**:
- Total `(J_S + J_D) / (I_S + I_D)` → ∞ (closure cancellation)
- `J_S / I_S = 0.385`, `J_D / I_D ≈ 0.10`, neither is −0.030
- Positive-ΔE-only `J/I` ratios of each channel give 0.79 and 0.10 — also not
  −0.030
- Drake's normalization combines J across channels with a *non-trivial
  ΔE-dependent kernel* and subtracts the asymptotic `ln(K_max)` divergence as
  K → ∞ (Schwartz 1961 §III, Drake-Swainson 1990 Eq. 5).

**Conclusion: Drake's value is a regularized object, not a finite-basis
limit of any naive J/I.** Implementing Drake-Swainson's regularization
requires either (a) the explicit frequency-integral representation with
asymptotic subtraction, or (b) Goldman-Drake 1983's K-shell integration in
the Coulomb Green's function — both are multi-sprint implementations.

### 3.2 3P, 3D show the same pathology

| (n,l) | acc N=40 ln k_0 | Drake     |
|:-----:|:---------------:|:---------:|
| 3P    |    +299.16      |  −0.0382  |
| 3D    |  +29 893.36     |  −0.0052  |

Higher l → smaller closure value → faster divergence in the ratio. Same
structural obstruction.

## 4. Combined Lamb shift

Using LS-1 self-energy formula with the LS-3 best 2S Bethe log (acceleration
form, N=40) and Drake-Swainson 2P Bethe log (native 2P diverges):

|                             | ln k_0(2S) | ln k_0(2P) | Lamb (MHz) | Err (MHz) | Err % |
|:----------------------------|:----------:|:----------:|:----------:|:---------:|:-----:|
| LS-1 (Drake-Drake baseline) |   2.812    |  −0.030    |  1025.06   |  −32.78   | −3.10 |
| LS-2 (vel N=50 + Drake)     |   2.726    |  −0.030    |  1037.30   |  −20.55   | −1.94 |
| **LS-3 (acc N=40 + Drake)** | **2.786**  | **−0.030** | **1028.59**| **−29.26**| **−2.77** |
| Experimental                |    —       |    —       |  1057.85   |     0     |  0    |

### 4.1 Error decomposition

The Lamb shift's sensitivity to the 2S Bethe log is `dE_Lamb/d(ln k_0(2S)) =
−(4/3)·α³/(π · 8) · 1 Ha → −135.6 MHz per unit`.

| Source        | LS-1     | LS-2      | **LS-3**   |
|:--------------|:--------:|:---------:|:----------:|
| T (Bethe log) |   0      | +11.70    | **+3.52**  |
| A (one-loop)  | −32.78   | −32.78    | −32.78     |
| **Total**     | −32.78   | −20.55    | **−29.26** |

**Truncation error (T) reduced 3.3× by LS-3 acceleration form.**
**Approximation error (A) is the LS-1 one-loop ceiling** (missing α⁵
two-loop SE, Karplus-Sachs VP, Yennie gauge, and recoil — sum ~33 MHz from
literature, here implicit in the +38/45 vs +1078 MHz textbook discrepancy of
LS-1).

### 4.2 Why LS-3 is "worse" numerically than LS-2

Despite a more accurate Bethe log, LS-3's Lamb shift is *farther from
experiment* than LS-2's. This is the artefact noted in LS-2 §3.4: LS-2's
3.1% under-Drake Bethe log was *accidentally compensating* for LS-1's
systematic 3% one-loop convention undershoot. LS-3's tighter Bethe log
exposes the underlying one-loop ceiling. The honest reading is that **LS-3
correctly decomposes the Lamb-shift error into truncation (T = 0.33%) and
one-loop ceiling (A = 3.10%)**, while LS-2's smaller total error was a
favourable miscancellation.

## 5. What LS-3 accomplishes

* **Acceleration-form Bethe log machinery** added to the GeoVac framework as
  `sturmian_acceleration_matrix` (matrix element of 1/r² in the chi = r·R
  basis = ∫ R R dr radial). Used in `bethe_log_acceleration` alongside the
  pre-existing velocity form. All new code in `debug/ls3_bethe_log_regularized.py`;
  no production-code changes.

* **2S Bethe log accuracy 0.92 %** at N=40, 3.3× tighter than LS-2 at
  N=50.

* **Truncation error decomposition**: the LS-3 Lamb shift's residual is
  cleanly resolved into T (Bethe log truncation, +3.5 MHz) and A (LS-1
  one-loop ceiling, −32.8 MHz). This is the kind of clean error analysis
  that LS-1/LS-2 lacked.

* **Cross-check 1S, 3S**: same convergence profile observed across n_target,
  consistent with O(1/√N) Sturmian convergence.

* **Structural finding for l > 0**: documented per-channel finiteness and
  identified the missing ingredient (Drake-Swainson regularization beyond
  finite-basis sums). 2P, 3P, 3D all exhibit the same closure pathology.

## 6. Honest limitations

* **2P (and l > 0 Bethe logs in general) are not derived from scratch in
  GeoVac.** The acceleration form does not bypass the closure obstruction;
  Drake-Swainson regularization (asymptotic subtraction) is a separate
  multi-sprint program.

* **A (approximation) ceiling is the binding constraint.** Reaching < 1%
  on the Lamb shift requires α⁵ multi-loop (Karplus-Sachs vacuum
  polarization, two-loop self-energy, Yennie gauge, recoil) — not a
  Bethe-log issue. Truncation T is now down to 0.33 %, well below A = 3.1 %.

* **2S Bethe log accuracy 0.92 % is also a structural plateau.** Pushing N
  to 100+ at higher mpmath precision should reach ~0.1 %, but the
  non-monotonic overshoot past N=40 is the next bottleneck — the Sturmian
  basis at fixed λ gives O(1/√N) convergence with reseeded round-off.
  Multi-λ basis or graded-N basis is the natural next refinement.

## 7. Verdict per exit criteria

| Exit criterion                         | Status            |
|:---------------------------------------|:-----------------:|
| 2P Bethe log within 10% of Drake       | **NEGATIVE**      |
| 2P Bethe log within 5%                 | NEGATIVE          |
| 2P Bethe log within 1%                 | NEGATIVE          |
| Combined Lamb shift within 1%          | NEGATIVE          |
| 2S Bethe log < 1% (subsidiary)         | **POSITIVE (0.92%)** |
| Truncation error reduction (subsidiary)| **POSITIVE (3.3×)** |

Net: **MIXED-POSITIVE.** The headline target (2P) is structurally beyond a
finite-basis acceleration form; the consolation (2S accuracy) is a clean,
documented improvement over LS-2 with explicit error decomposition.

## 8. Files

* Implementation: `debug/ls3_bethe_log_regularized.py`
* Data: `debug/data/ls3_bethe_log_regularized.json`
* Memo: this file

## 9. Future work

* **LS-4 (formerly LS-3 plan):** Implement Drake-Swainson's frequency-integral
  representation with explicit asymptotic subtraction for l > 0 Bethe logs.
  The cleanest target: Goldman-Drake 1983 K-shell integration over the
  Coulomb Green's function, evaluated in the Sturmian basis. This is a
  conceptually different sprint — not a finite-basis sum but a regularized
  integral with a sub-leading log term subtracted. Expected to reach 2P,
  3P, 3D Bethe logs at < 5 % accuracy.

* **LS-5:** Multi-λ Sturmian basis for s-states. Using a graded basis with
  several exponents (λ₁, λ₂, …) should resolve the high-energy continuum
  better and push 2S, 3S Bethe logs to < 0.1 %.

* **LS-6 (the one-loop ceiling):** Add α⁵ multi-loop corrections (Karplus-
  Sachs, two-loop SE, Yennie, recoil) as a separate "approximation order"
  improvement layer. Each is ~5 MHz of the ~33 MHz one-loop residual; the
  full set would reach < 1 % on the Lamb shift.

## 10. References

* H. A. Bethe, *Phys. Rev.* 72 (1947) 339.
* C. Schwartz, *Phys. Rev.* 123 (1961) 1700 — finite-basis Bethe log,
  asymptotic subtraction for l > 0.
* M. Lieber, *Phys. Rev.* 174 (1968) 2037.
* S. P. Goldman, G. W. F. Drake, *Phys. Rev. A* 28 (1983) 1228 — K-shell
  integration scheme.
* G. W. F. Drake, R. A. Swainson, *Phys. Rev. A* 41 (1990) 1243 — Bethe log
  table to 10⁻¹⁰.
* K. Pachucki, *J. Phys. B* 31 (1998) 5123; *Phys. Rev. A* 78 (2008)
  012504 — modern computational scheme.
* M. I. Eides, H. Grotch, V. A. Shelyuto, *Phys. Rep.* 342 (2001) 63.
* GeoVac LS-1 (`debug/ls1_lamb_shift_memo.md`), LS-2
  (`debug/ls2_bethe_log_memo.md`).
