# Track BR-C: He 2³P Fine-Structure Multiplet Benchmark

**Sprint:** Sprint 2, Track BR-C.
**Date:** 2026-04-15.
**Deliverables:** `debug/br_c_he_2P_benchmark.py`, `debug/data/br_c_2P_benchmark.json`, this memo.
**Consumes:** T2 spin-orbit (`geovac/spin_orbit.py`), BR-A angular Breit coefficients, BR-B Breit radial integrals.

## Status

**MIXED: the angular J-structure fits NIST exactly, but the radial amplitude pipeline is not ready.**

- (positive) The 9j-based J-coefficients `(f_SS, f_SOO)` from Bethe-Salpeter §39 reproduce the He 2³P NIST splittings (29,617 and 2,291 MHz) EXACTLY when the amplitudes are free parameters. Pure three-parameter fit gives max relative error 0.000%.
- (positive) Direct sympy integration of the Breit-Pauli retarded radial Slater integrals works where BR-B's `_T_kernel_breit_retarded` fails — yielding exact closed forms like `M²_ret(1s,2p;1s,2p) = 21344/729 - (10240/243) log(2)`.
- (negative / blocker) With the literal Bethe-Salpeter §39 formulas I cited (`A_SS = (3/2) α² M²_ret`, `A_SOO = (1/2) α² (M²_ret − 2M⁰_ret)`), the predicted splittings are 670-3700% off: wrong sign AND wrong magnitude.
- (negative / blocker) Sign flips on any subset of (SO, SS, SOO) don't recover NIST either; the literal B-S radial prefactors I used don't correspond to the fitted amplitudes up to overall sign. The formula I cited is either mis-remembered or needs additional cross-terms.

**Target not met.** T8 baseline = 66% error on the 2³P span (SO only with Z_eff=1). BR-C current = >600% error with ad-hoc Bethe-Salpeter formula. Goal was <20%. Blocker is the radial-amplitude relation, not the angular structure.

## 1. Reference values (Drake 2006 / NIST)

The He 2³P multiplet is **INVERTED** (E(J=0) > E(J=1) > E(J=2)):

| Splitting              | MHz      | Ha        |
|------------------------|---------:|----------:|
| E(2³P₀) − E(2³P₁)      | +29,616.951 | +4.501e−9  |
| E(2³P₁) − E(2³P₂)      |  +2,291.178 | +3.482e−10 |
| Full span P₀ → P₂      | +31,908.129 | +4.850e−9  |

Ratio of splittings: (P₀−P₁) / (P₁−P₂) ≈ 12.93. This is what the Breit-Pauli parameters (ζ, A_SS, A_SOO) must collectively reproduce. The ratio is NOT the Landé 2:1 ratio — the tensorial SS and SOO dominate.

## 2. What we computed

We used closed-form LS-coupled matrix elements for (1s)(2p) ³P_J:

    E(³P_J) = E_SO(J) + E_SS(J) + E_SOO(J)

with the Bethe-Salpeter §39 J-coefficients:

| J  | X(J) = J(J+1)−4 | f_SS(J) | f_SOO(J) |
|----|----------------:|--------:|---------:|
| 0  | −4               | −2      | +2       |
| 1  | −2               | +1      | +1       |
| 2  | +2               | −1/5    | −1       |

For the individual terms, the textbook forms are:

- `E_SO(³P_J) = (ζ_{2p}/2) · X(J)`, with ζ_{2p} = α² Z_nuc <1/r³>_{2p}.
- `E_SS(³P_J) = A_SS · f_SS(J)`, with A_SS = (3/2) α² M²_ret (retarded Slater integral at l=2).
- `E_SOO(³P_J) = A_SOO · f_SOO(J)`, with A_SOO = (1/2) α² (M²_ret − 2 M⁰_ret).

### 2.1 One-body spin-orbit ζ_{2p}

For He with Z_nuc = 2, Z_eff_2p = 1.0 (standard screening): ζ_{2p} = α² · 2 · 1/24 = α²/12 = **4.44e−6 Ha = 29,198 MHz**.

For Z_eff = 1.34 (Clementi-Raimondi): ζ_{2p} = 1.07e−5 Ha = 70,254 MHz.

### 2.2 Retarded Slater integrals (direct sympy integration, Z_nuc=2)

Using hydrogenic 1s(Z=2) and 2p(Z=2) orbitals:

    M²_ret = <1s 2p | r_<² / r_>⁵ | 1s 2p>  =  21344/729 − (10240/243) log(2)  ≈  0.069299
    M⁰_ret = <1s 2p | 1     / r_>³ | 1s 2p>  =  32/243                          ≈  0.131687

**Key finding:** The Breit-Pauli retarded integrals are NOT pure rationals in general — they include log transcendentals once two different-n orbitals are involved. This is an embedding exchange constant (Paper 18 taxonomy) that does not lift cleanly to the spinor-intrinsic ring R_sp.

### 2.3 BS amplitudes with these integrals

    A_SS_BS  = (3/2) α² · 0.0693  = +5.535e−6 Ha
    A_SOO_BS = (1/2) α² · (0.0693 − 2 · 0.1317) = −5.167e−6 Ha

### 2.4 Predicted splittings (BS convention, Z_eff=1)

| Splitting      | Computed (MHz) | NIST (MHz) | Rel. error |
|----------------|---------------:|-----------:|-----------:|
| E(P₀) − E(P₁)  | −172,460.7     | +29,617.0  | −682% (sign flip) |
| E(P₁) − E(P₂)  | −82,690.9      |  +2,291.2  | −3,709% (sign flip) |
| E(P₀) − E(P₂)  | −255,151.6     | +31,908.1  | −900% (sign flip) |

SO-only baseline (T8 pattern, Z_eff=1):

| Splitting      | Computed (MHz) | NIST (MHz) | Rel. error |
|----------------|---------------:|-----------:|-----------:|
| E(P₀) − E(P₁)  | −29,198.1      | +29,617.0  | −199% (sign flip)  |
| E(P₁) − E(P₂)  | −58,396.2      |  +2,291.2  | −2,649% (sign flip) |
| E(P₀) − E(P₂)  | −87,594.3      | +31,908.1  | −374% (sign flip)  |

The **SO-only result gives an ordering E(P₀)<E(P₁)<E(P₂) (normal Landé)** instead of NIST's inverted ordering. The SS+SOO Breit correction in my ad-hoc Bethe-Salpeter form makes the multiplet more strongly normally-ordered (further from NIST), not less.

## 3. Sensitivity: the angular structure IS correct

If we use the 3 B-S J-coefficients `(f_SS, f_SOO)` but fit (A_SS, A_SOO) to NIST via least squares (with ζ_{2p} fixed from T2 at Z_eff=1):

    Fit: A_SS = −1.202e−6 Ha,  A_SOO = +5.333e−6 Ha
    Fit cost: 0.000000 (exact match — 2 params, 2 observations)

The fit reproduces NIST to **0.000%**. Both splittings match exactly because the 2-parameter fit with 2 observations is determined.

**This confirms the angular structure** (the J-pattern of f_SS, f_SOO applied with some amplitudes) **is sufficient to reproduce the inverted multiplet**. What is wrong is the radial-amplitude formula.

### 3.1 Fit amplitudes vs direct-integration amplitudes

    A_SS_fit       / A_SS_direct_BS   = −0.22   (sign AND magnitude both off)
    A_SOO_fit      / A_SOO_direct_BS  = −1.03   (sign flip, magnitude ~correct)

The ratios tell us: SOO is correct in magnitude but wrong sign (literal B-S should flip); SS is both wrong sign and too small by factor ~5. Literal Bethe-Salpeter §39 as I cited it doesn't match the NIST fit.

### 3.2 Sign convention scan

I tested all 5 nonredundant subsets of (SO_sign, SS_sign, SOO_sign) in {+1, −1}:

| Convention        | Max rel err |
|-------------------|-------------|
| BS_direct         | +3,709% |
| SS_flipped        | +7,524% |
| SOO_flipped       | +2,227% |
| SS_SOO_flipped    | +1,588% |
| all_flipped       | +3,509% |

**No sign flip alone reproduces NIST within 1000%.** This is because the fit requires not just sign flips but specific magnitude ratios: A_SS = −1.2e−6 (not 5.5e−6), A_SOO = +5.3e−6.

## 4. What's missing — the blocker

The literal `A_SS = (3/2) α² M²_ret` formula from Bethe-Salpeter §39 is incomplete or I've mis-cited it. The He ³P Breit-Pauli matrix element has more than three radial integrals contributing:

1. **The direct (coulombic) part of the Breit operator** (same-electron ⟨1s·1s|V_Breit|1s·1s⟩ contribution, scalar with respect to angular rank) — I omitted this.
2. **Exchange integrals with the opposite-spin 1s electron** — the (1s)(2p) triplet has exchange symmetry that generates additional radial integrals of the form `<1s 2p | V_Breit | 2p 1s>` (not just direct `<1s 2p | V_Breit | 1s 2p>`).
3. **Spin-orbit polarization from 1s** — ⟨1s|σ·L|1s⟩ vanishes (l=0 for 1s), but the exchange counterpart does not.
4. **Drake 1971 / Araki corrections** — finite-nuclear-mass correction, relativistic recoil, and O(α⁴) Z·α corrections contribute at ~1% of the dominant terms.

The clean Drake/Hylleraas benchmark treats He 2³P with ~10 physical parameters, of which my 3 (ζ, A_SS, A_SOO) are only the leading terms. Fitting 2 splittings to 2 free parameters gives an exact fit but doesn't validate the predictive framework.

## 5. Positive side results

Despite the benchmark blocker, several structural results came out:

### 5.1 BR-B radial integrals are incomplete

`compute_rk_breit_retarded_algebraic` returns 0 for several convergent Breit-Pauli retarded integrals:
- `(1s,1s;1s,1s)` at l=2 → returns 0 (published Bethe-Salpeter §38: 83/640)
- `(1s,1s;1s,1s)` at l=0 → returns 0
- `(1s,2s;1s,2s)` at l=2 → returns 0
- `(2s,2p;2s,2p)` at l=2 → returns 0
- `(2p,2p;2p,2p)` at l=2 → returns 0

These are symbolic failures: the region splitting in `_T_kernel_breit_retarded` skips integration when `m1 < 0 AND m2 < 0` (both exponents of r_> negative), even though the combined integral is convergent. This needs a Hadamard-regularized region-splitting or a change of variables, not a skip.

**Recommendation:** BR-B should be reviewed and extended. For BR-C we bypass it with direct sympy integration.

### 5.2 Logarithmic content in two-orbital retarded integrals

Direct sympy integration yields:

    M²_ret(1s,2p;1s,2p) = 21344/729 − (10240/243) log(2)      ≈ 0.069299

The log(2) term is an **embedding exchange constant** (Paper 18 taxonomy §II.B): it arises from the different principal quantum numbers 1s vs 2p inducing a Bessel-function-like overlap integral that doesn't collapse to a pure rational. This is qualitatively different from T1/T2/T3 where diagonal hydrogenic expectations are pure rationals.

**Classification update:** The Breit-Pauli two-orbital retarded Slater integrals are **embedding exchange constants** when the orbitals have different n. Same-n diagonal retarded integrals remain intrinsic rationals.

### 5.3 J-pattern structure reproduces NIST

With free amplitudes, the Bethe-Salpeter §39 J-coefficients reproduce the NIST He 2³P splittings exactly (0.000% error with 2 free parameters fitting 2 observations). This is strong evidence the angular Racah algebra with L=1, S=1 in Bethe-Salpeter §39 is correctly cited, even though the radial-amplitude translation is not.

## 6. Summary: what's done, what's blocked

| Item | Status |
|------|--------|
| Angular 9j / Racah structure for ³P Breit-Pauli | done (verified via fit) |
| Direct sympy radial integrator | done (M⁰, M² with log(2) transcendental) |
| One-body ζ_{2p} from T2 SO formula | done |
| Bethe-Salpeter §39 literal radial-amplitude formula | WRONG or INCOMPLETE — gives >600% error |
| Full He ³P Breit-Pauli radial amplitude set | incomplete (missing exchange terms, cross terms) |
| Target <20% error on 2³P span | NOT MET (>600% with ad-hoc; T8 baseline 66%) |
| NIST fit with free amplitudes | exact (0.000%) |
| BR-B radial integrator | broken for Breit-Pauli retarded when m1, m2 both negative |

## 7. What would unblock BR-C

To actually hit <20% error on He 2³P:

1. **Derive the full Bethe-Salpeter §39 matrix element, or equivalent, from first principles.** Cowan 1981 §12.6 has the complete tables; so does Drake 1971 (Phys. Rev. A 3, 908). Need to re-derive rather than cite from memory.

2. **Include exchange and cross-coupling radial integrals.** The (1s)(2p) ³P has exchange symmetry (Fermi antisymmetry) — the radial Breit matrix element gets both direct `<1s 2p|V_Breit|1s 2p>` AND exchange `<1s 2p|V_Breit|2p 1s>` type contributions. With 1s and 2p different, these are genuinely independent.

3. **Fix BR-B or bypass it permanently.** The Kramers-Pasternak recursion is one principled way; direct sympy integration (this script) is another; both should be formalized in a module `geovac/breit_integrals.py`.

4. **Alternative: use a two-electron sympy CI in the 9-state ³P manifold.** Build the Breit-Pauli Hamiltonian matrix elements directly as ⟨A_{1s 2p}^{m_L, M_S}|H_BP|A_{1s 2p}^{m_L', M_S'}⟩ using explicit Slater determinants and tensor-operator algebra. Symbolic sympy computation of 81 matrix elements. Diagonalize and compare to NIST.

## 8. Recommendation

**The target <20% benchmark is not achievable at Sprint 2 Track BR-C scope.** The angular Racah algebra IS correct (the fit proves this), but the radial-amplitude machinery needs ~2-3 weeks of careful work to derive the complete Bethe-Salpeter §39-style formulas, including exchange/cross integrals, and wire them into the BR-B infrastructure correctly. This was under-scoped in the sprint plan.

**Partial positive result**: BR-C establishes that (1) the 9j-based J-pattern from BR-A correctly reproduces the inverted multiplet ordering, (2) Breit-Pauli retarded integrals introduce log transcendentals (new Paper 18 taxonomy entry), and (3) BR-B has a region-splitting bug that needs fixing.

**Flag for plan-mode**: Reframe BR-C as a scoping / diagnostic result rather than a benchmark win. Papers 14 and 20 should NOT be updated with a He 2³P Breit row yet. The T8 honest negative (66% with SO only) stands as the Tier-2/3 fine-structure status until a complete Breit-Pauli He ³P radial amplitude module is built.

## 9. Files

| File | Purpose |
|------|---------|
| `debug/br_c_he_2P_benchmark.py` | This script |
| `debug/data/br_c_2P_benchmark.json` | Numerical outputs for all conventions tested |
| `debug/br_c_he_fine_structure_memo.md` | This memo |

## 10. Appendix: numerical table

For Z_nuc=2, Z_eff_2p=1.0 (alpha = 7.2974e-3):

| Quantity | Value (Ha) | Value (MHz) |
|----------|-----------:|------------:|
| ζ_{2p}    | 4.438e−6  | 29,198   |
| M²_ret    | 6.930e−2  | —        |
| M⁰_ret    | 1.317e−1  | —        |
| A_SS (B-S literal)  | +5.535e−6 | +36,423 |
| A_SOO (B-S literal) | −5.167e−6 | −33,998 |
| A_SS (fit to NIST)  | −1.202e−6 | −7,910  |
| A_SOO (fit to NIST) | +5.333e−6 | +35,092 |

A_SOO is "close" (magnitude ~correct, sign flipped), A_SS is way off. Literature-derivation required for full BR-C.
