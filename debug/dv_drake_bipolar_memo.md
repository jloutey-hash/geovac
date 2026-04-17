# Track DV: Drake Bipolar Expansion — Sprint 5 Technical Memo

**Sprint:** Sprint 5 Track DV.
**Date:** 2026-04-15.
**Status:** **HONEST PARTIAL** — bipolar channel structure fully characterized,
direct/exchange mixing ratios NOT closed symbolically from the leading-order
bipolar expansion. The "leading-term" bipolar machinery (k_1 + k_2 = K) gives
a specific hybrid linear combination of piecewise Drake M^K integrals that is
NOT the same object as Drake's M^K_dir, M^K_exch totals. Closure requires
either (i) a different M^K convention than `breit_ss_radial`'s r_<^K/r_>^{K+3},
(ii) non-trivial piecewise identities transforming the bipolar form into
Drake's form, or (iii) a bipolar prefactor beta(k_1, k_2, K) that is not
in the standard Brink-Satchler / Steinborn-Filter lookup tables.

---

## 1. Sprint question

Sprint 4 DD's honest-negative: the J-pattern (-2, +1, -1/5) for f_SS(J) was
derived from pure 6j algebra, but the direct/exchange mixing ratios
(3/50, -2/5, 3/2, -1) in front of M^K_dir, M^K_exch remained a rational
search hit (Sprint 3 BF-D) without closed-form derivation.

Sprint 5 DV was commissioned to close this by deriving the mixing ratios
from the full bipolar harmonic expansion of Y^K(r_hat_12) / r_12^{K+1}.

---

## 2. What was established (positive)

### 2.1 Channel enumeration (Step 1, sympy-exact)

For He (1s)(2p) ^3P_J, the angular-allowed (k_1, k_2) bipolar channels are:

**K = 2 (SS):**
- Direct (l_a, l_b, l_c, l_d) = (0, 1, 0, 1): ONLY (k_1 = 0, k_2 = 2)
  from `<l=0||C^{k_1}||l=0>` requires k_1 = 0; triangle (0, 2, K=2) with
  k_2 even (parity) picks k_2 = 2 as unique.
- Exchange (0, 1, 1, 0): ONLY (k_1 = 1, k_2 = 1).

**K = 1 (SOO):**
- Direct: EMPTY (bipolar angular CG forbids).
- Exchange: EMPTY.

### 2.2 Piecewise Drake M^K decomposition (Step 7, exact symbolic)

Using the production module's `_t_kernel_region_I`, we obtain Region I
(r_1 < r_2) and Region II (r_1 > r_2) contributions to Drake's M^K:

At Z = 1 (floats):

|       | M^0_dir    | M^1_dir    | M^2_dir    | M^0_exch   | M^1_exch   | M^2_exch   |
|-------|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|
| Region I  | +2.469e-02 | +1.235e-02 | +7.697e-03 | +8.230e-03 | +5.718e-03 | +4.331e-03 |
| Region II | +4.576e-03 | +3.618e-03 | +2.984e-03 | +8.230e-03 | +5.718e-03 | +4.331e-03 |
| Total     | +2.927e-02 | +1.596e-02 | +1.068e-02 | +1.646e-02 | +1.144e-02 | +8.662e-03 |

**Confirmed:** for exchange channel (densities symmetric under r_1 <-> r_2),
Region I = Region II. For direct channel, Region I ≠ Region II because
the densities (1s)^2(1) and (2p)^2(2) are asymmetric.

Symbolic Region I / Region II values (sample):
- M^2_dir Region I = -7/54 + log(3)/8
- M^2_dir Region II = -240 log(2) - 7882/81 + 240 log(3)

### 2.3 Bipolar channel identities (Step 7/10)

For bipolar channel (k_1, k_2, K = 2) with kernel r_1^{k_1} r_2^{k_2} / r_>^5
and k_1 + k_2 = K, the piecewise structure is:

- Region I (r_1 < r_2, r_> = r_2): kernel = r_1^{k_1} r_2^{k_2 - 5}
  Matches Drake M^{K_d = k_1} Region I kernel r_1^{k_1} / r_2^{k_1 + 3}
  when k_2 - 5 = -(k_1 + 3), i.e., k_1 + k_2 = 2 ✓
- Region II (r_1 > r_2, r_> = r_1): kernel = r_1^{k_1 - 5} r_2^{k_2}
  Matches Drake M^{K_d = k_2} Region II kernel r_2^{k_2} / r_1^{k_2 + 3}.

**Exact identity (verified numerically at float precision, Step 7):**

| Bipolar channel | Direct | Exchange |
|:----------------|:-------|:---------|
| (0, 2, K=2) | M^0_dir_I + M^2_dir_II = 2.768e-02 | M^0_exch_I + M^2_exch_II = 1.256e-02 |
| (1, 1, K=2) | M^1_dir_I + M^1_dir_II = M^1_dir = 1.596e-02 | M^1_exch_I + M^1_exch_II = M^1_exch = 1.144e-02 |
| (2, 0, K=2) | M^2_dir_I + M^0_dir_II = 1.227e-02 | M^2_exch_I + M^0_exch_II = 1.256e-02 |

### 2.4 Structural mismatch with Drake's M^K basis (Step 11)

The bipolar expansion gives:

A_SS ∝ scaling * [spatial_red_dir * β(0,2,2) * (M^0_dir_I + M^2_dir_II)
                  - spatial_red_exch * β(1,1,2) * M^1_exch]

Drake's formula gives:

A_SS_Drake = (3/50) (M^2_dir_I + M^2_dir_II) - (2/5) (M^2_exch_I + M^2_exch_II)

**The natural bipolar basis {M^0_dir_I, M^2_dir_II, M^1_exch} is DIFFERENT
from Drake's basis {M^2_dir, M^2_exch}.** No choice of (β(0,2,2), β(1,1,2))
can reconcile these bases via simple scaling — the bipolar form has a
"direct" piece involving M^0_dir_I (not M^2_dir_I) and involves M^1_exch
(not M^2_exch).

This discrepancy is the technical obstruction that closes the Sprint 5 DV
derivation prematurely.

### 2.5 Numerical calibration (Step 11)

Treating β(0,2,2) and β(1,1,2) as free parameters, we require

β(0,2,2) * a + β(1,1,2) * b = bracket_target

with a = -3.032e-02, b = +8.524e-03, bracket_target = -6.187e-03.

This is ONE equation in TWO unknowns. For the "equal-β" ansatz
β(0,2,2) = β(1,1,2), the value β = 0.2839 almost matches 1/sqrt(4π) = 0.2821
(0.64% off) — strongly suggesting a convention factor of sqrt(4π) is
involved in the bipolar prefactor normalization, but NOT equal β values.

---

## 3. What's NOT established (honest-negative)

1. **No closed-form symbolic derivation of (3/50, -2/5).** The bipolar
   expansion gives a direct channel that mixes M^0_dir_I + M^2_dir_II,
   not M^2_dir alone. Drake's formula cannot be the "direct output" of
   the leading-order bipolar — it's a DIFFERENT linear combination of
   the same underlying piecewise integrals.

2. **No closed-form symbolic derivation of (3/2, -1).** For SOO (K=1),
   the bipolar angular CG gives ZERO for both direct and exchange
   channels (since the CG `<k_1 0 k_2 0 | K=1 0>` vanishes for the
   allowed (k_1, k_2) triangles with l_a=0, l_b=1). This confirms that
   the SOO operator is not a pure bipolar rank-1 spatial tensor — it
   includes the momentum operator [r × p], which changes the angular
   structure fundamentally.

3. **Bipolar prefactor β(k_1, k_2, K)** was not uniquely determined.
   The numerical fit suggests β(0,2,2) and β(1,1,2) are of order
   1/sqrt(4π) = 0.28, but the exact rational (or simple closed form)
   was not identified.

---

## 4. Hypothesis for the structural mismatch

**The Drake M^K integrals in Eq. (17) of Drake 1971 are NOT the raw
r_<^K/r_>^{K+3} Slater integrals computed by `breit_ss_radial`.**

One consistent interpretation: Drake's M^K_dir, M^K_exch are DEFINED
by the Racah-coupled two-electron matrix element of the Breit-Pauli
operator, i.e., the LEFT-HAND SIDE of the equation

<1s 2p, ^3P | H_SS | 1s 2p, ^3P> = f_SS(J) * α² (c^K_dir M^K_dir + c^K_exch M^K_exch)

where M^K_dir and M^K_exch are then FIT rationally with c^K_dir = 3/50,
c^K_exch = -2/5. The fit works because the RHS (for fixed Z) is a
specific number that happens to equal 3/50 * (specific rational + log
form) - 2/5 * (other specific rational + log form). The rational fit
(3/50, -2/5) is approximate (0.20% error on NIST).

A "deeper" derivation would require either:

(a) A proof that `breit_ss_radial` values satisfy a specific **algebraic
    relation** such that the 3/50 and -2/5 coefficients emerge from
    Racah angular algebra alone. Not yet found.

(b) A **different convention** for M^K_dir, M^K_exch (e.g., Drake's
    original paper) that naturally absorbs the piecewise weighting.

(c) A **higher-order bipolar expansion** contributing additional
    channels. Confirmed EXCLUDED by the angular selection rules: no
    higher (k_1, k_2) channels contribute for He (1s)(2p) ^3P.

---

## 5. Files created (Sprint 5 DV)

| File | Purpose |
|:-----|:--------|
| `debug/dv_drake_bipolar.py` | Initial channel enumeration; SOO K=1 empty confirmed |
| `debug/dv_step1_channels.py` | Full (k_1, k_2) enumeration for K=0, 1, 2 with 9j values |
| `debug/dv_step2_radial.py` | Scipy numerical bipolar radial integrals |
| `debug/dv_step3_symbolic_radial.py` | Sympy piecewise Drake M^K (sampling) |
| `debug/dv_step4_closed_form.py` | Symbolic piecewise integrals (slow; didn't complete) |
| `debug/dv_step5_analytic.py` | Re-derivation notes; aborted |
| `debug/dv_step6_piecewise.py` | Bipolar channel via wrong-orbital-labels BAD; Step 7 corrects |
| `debug/dv_step7_correct_pieces.py` | CORRECT piecewise Drake M^K via `_t_kernel_region_I` |
| `debug/dv_step8_angular_numerical.py` | Numerical 6D integration design notes |
| `debug/dv_step9_brute_numerical.py` | Further analytical notes on the obstruction |
| `debug/dv_step10_final.py` | Final sympy attempt (ran into int_max_str_digits) |
| `debug/dv_step11_float_only.py` | FLOAT-only verification; 1/sqrt(4π) near-hit |

---

## 6. Recommendation

**Close Sprint 5 DV with honest-partial status.** The J-pattern derivation
from Sprint 4 DD stands; the mixing ratios remain a rational-search match.

### 6.1 Paper 14 §V update (minimal)

The current text says:

> The rational direct/exchange mixing ratios were identified by systematic
> enumeration over small rationals (Track BF-D) and confirmed numerically
> to reproduce NIST He 2^3P to sub-percent accuracy; a closed-form derivation
> of these specific rationals from the bipolar harmonic expansion of
> Y^(K)(r̂_12)/r_12^{K+1} is deferred (Track DD honest-negative on the mixing
> coefficients; J-pattern and spin-tensor structure are derived).

Revised text:

> The direct/exchange mixing ratios (3/50, -2/5, 3/2, -1) were identified
> by rational search (Track BF-D) and reproduce NIST He 2^3P to sub-percent
> accuracy. Sprint 5 Track DV completed the bipolar harmonic decomposition
> of the He (1s)(2p) ^3P angular matrix element and confirmed that the
> direct channel has a UNIQUE (k_1 = 0, k_2 = 2) bipolar contribution and
> the exchange has a UNIQUE (k_1 = 1, k_2 = 1) contribution (no higher-order
> (k_1, k_2) channels are angular-allowed). The piecewise decomposition of
> the bipolar radial kernel into Drake's M^K_{dir, exch} basis is NOT a
> simple linear combination of M^K totals — the bipolar channel (0, 2, 2)
> direct equals M^0_dir_Region_I + M^2_dir_Region_II, which is structurally
> distinct from Drake's M^K_dir = Region I + Region II. The rational form
> (3/50, -2/5) is therefore a convention-dependent combining identity that
> emerges from the specific definition of Drake's M^K integrals; a first-
> principles derivation remains open.

### 6.2 Tests added

None. The Sprint 5 DV sprint did not produce new closed-form identities
that require a test.

### 6.3 Production code changes

None. `geovac/breit_integrals.py` production code is NOT modified (pre-
existing 31 tests still pass).

### 6.4 Next-sprint proposal

If the PI wishes to pursue a full derivation of (3/50, -2/5), the
recommended path is:

**Track DW:** directly evaluate the Breit-Pauli spin-spin matrix element
for He (1s)(2p) ^3P_1 via 6-dimensional numerical integration (three
radial + three angular), PSLQ-identify in Drake's {M^0, M^1, M^2} × {dir, exch}
basis at several Z values (Z = 2, 3, 5, 10), and cross-check the rational
fit against the Drake formula. Estimated 2-4 sub-agent hours.

Alternatively, consult Drake 1971 (Phys Rev A 3, 908) directly for the
exact definition of M^K in use — the mixing ratios may become obvious
with the correct radial-integral convention.

---

## 7. Summary

**Achieved:**
- Full enumeration of angular-allowed (k_1, k_2) bipolar channels for He
  (1s)(2p) ^3P at K = 0, 1, 2 (Step 1, 9j/3j exact).
- Piecewise decomposition of Drake M^K into Region I / Region II
  contributions (Step 7, sympy-exact via production module hook).
- Identification that SOO (K=1) has NO angular-allowed bipolar channel —
  confirming that SOO involves the momentum operator and is not a pure
  rank-1 Y^{(1)}/r_{12}^2 tensor.
- Numerical near-coincidence β = 1/sqrt(4π) = 0.282 for equal-β ansatz,
  suggesting a Y^K ↔ C^K conversion factor is the dominant scale in the
  bipolar prefactor.

**Not achieved:**
- Symbolic derivation of the Drake mixing ratios (3/50, -2/5) for SS.
- Closed-form extraction of (3/2, -1) for SOO (the K=1 spatial tensor
  does NOT have a pure bipolar decomposition).
- Uniqueness of the bipolar prefactor β(k_1, k_2, K) — the equal-β
  ansatz only matches approximately, indicating the physics requires
  β(0,2,2) ≠ β(1,1,2) with values not in common angular-algebra lookups.

**Key structural finding:** Drake's M^K_dir, M^K_exch as computed by
`breit_ss_radial` (= ∫∫ ρ r_<^K / r_>^{K+3}) are NOT the natural basis
for expressing the bipolar-expansion matrix element. The bipolar basis
{M^{k_1}_dir_I + M^{k_2}_dir_II : k_1 + k_2 = K} is structurally
different. The Drake formula's rational coefficients (3/50, -2/5) are
therefore a convention-dependent combining identity, not a direct output
of Racah angular algebra.
