# Door 2 — BW Wedge Entropy Re-read (forcing-catalogue forward-run)

**Date:** 2026-06-01
**Driver:** `debug/door2_bw_entropy_reread.py`, `debug/door2_constant_pslq_tighten.py`
**Data:** `debug/data/door2_bw_entropy_reread.json`, `debug/data/door2_constant_pslq_tighten.json`
**Source crack:** `docs/forcing_catalogue.md` Door 2 ("BW wedge entropy mis-filed as failed area-law")
**Prior record:** `debug/bh_phase0_entanglement_entropy_memo.md` (filed NEGATIVE: "area-law REJECTED, R²=0.83")

---

## VERDICT: PARTIAL (clean exact structure beneath the negative; F-theorem link RULED OUT)

The forcing catalogue was **right that a positive result was shelved under the
negative**, and **wrong about what it is**. The shelved positive is *better* than
"slope ≈ 2 ≈ Cardy–Calabrese": the entropy has a fully **closed-form asymptotic**
with **slope provably exactly 2** and an **exact closed-form constant**. But the
suspected link to the Paper 50 F-theorem F-coefficient is **cleanly ruled out** —
the constant lives in a different transcendental ring, and the slope is a boundary
*dimension*, not a central charge.

Mapped to the decision gate:
- **slope provably → 2:** YES (proved analytically, not fitted).
- **constant has a clean closed form:** YES, `C_inf = coth(1) − log(2 sinh 1)`.
- **stateable link to the F-theorem:** NO — ruled out (different ring; this is the
  *informative* half of the result, consistent with what Paper 50 already says).

So: not a full "door" (no F-link), not a "wall" (the negative was incompletely
filed — there IS exact structure). **PARTIAL, leaning vindication-of-method.**

---

## What the entropy actually is (exact)

The wedge KMS state is `ρ_W = e^{−K_α}/Z` at canonical BW (β = 2π baked into the
integer spectrum of `K_α`). On the truncated full-Dirac triple:

- eigenvalue `k = 2m+1`, `m = 0 … n_max−1` (odd integers 1,3,…,2n_max−1);
- degeneracy `g_m = (n_max − m)(n_max − m + 1)`;
- per-state Boltzmann weight `e^{−(2m+1)}`.

Entropy in closed form (`S = log Z + ⟨K⟩`, both exact sums):

```
Z      = Σ_m g_m e^{−(2m+1)}
⟨K⟩    = Σ_m g_m (2m+1) e^{−(2m+1)} / Z
S(n)   = log Z + ⟨K⟩
```

This closed form reproduces the production `geovac/modular_hamiltonian.py` numerics
to machine precision (diff ≤ 4×10⁻¹⁶ at n = 2,5,7,10). The entropy is therefore an
**analytical object** — no fitting is needed at all, and the entire "curve-fit"
exercise is replaced by an exact asymptotic.

---

## (1) Slope: provably EXACTLY 2 (the 1.963 was a small-n artifact)

Resumming the geometric series for `n_max → ∞` (the upper limit `n−1 → ∞` adds only
`O(e^{−2n})`):

```
g_m = n² + (1−2m) n + m(m−1)
Z   = n²·A0 + n·A1 + A2,   A0 = Σ e^{−(2m+1)} = e^{−1}/(1−e^{−2}) = 1/(2 sinh 1)
log Z = 2 log n + log A0 + (A1/A0)/n + O(1/n²)
⟨K⟩  = B0/A0 + O(1/n),     B0 = Σ (2m+1) e^{−(2m+1)} = e^{−1}(1+e^{−2})/(1−e^{−2})²
```

The `n²` leading term **forces the slope = 2 exactly**. Numerically:

| window | least-sq slope | R² |
|:---|---:|---:|
| [2..7]   (original) | 1.9435 | 0.99991 |
| [8..20]            | 1.9962 | 0.9999996 |
| [30..200]          | 2.00026 | 1.0000000 |
| [350..2000]        | 2.00005 | 1.0000000 |

Local log-derivative `dS/d(log n)` climbs monotonically 1.896 → 2.0000 across
n ∈ [2, 2000]. The original `1.963` was the small-n value of a sequence converging
to exactly 2. **This is a degeneracy log: `S → log(n_eq) = log(n²) = 2 log n`,** where
`n_eq = n(n+1)` is the equator-shell multiplicity. The "2" is `dim(S²_equator)`,
exact.

## (2) Constant: exact closed form `C_inf = coth(1) − log(2 sinh 1)`

The `n → ∞` intercept of `S − 2 log n` is

```
C_inf = log A0 + B0/A0 = −log(2 sinh 1) + coth(1) = 0.4584487433681903606…
```

- residual of the closed form vs the resummed value: `2×10⁻⁸¹` (dps=80);
- PSLQ on `{C_inf, coth1, log2, log(sinh1)}` returns `[−1, 1, −1, −1]`, i.e.
  `C_inf = coth(1) − log2 − log(sinh1) = coth(1) − log(2 sinh 1)`, residual exactly 0;
- equivalent purely in `e`: `C_inf = (e²+1)/(e²−1) − log((e²−1)/e)`.

**The original fit's `0.540 / 0.56` was wrong** — a small-n contamination. With the
slope biased low (1.963), the fitted intercept was biased high (0.540) to compensate.
Forcing the true slope 2 and reading the `n → ∞` intercept gives `0.4584…`.

## (3) Cardy–Calabrese reading: FORCED ANALOGY, not an identity

A 1+1D Cardy–Calabrese entanglement log is `S = (c_eff/6) log L`. Here the slope is
`dS/d(log n_max) = 2`. Forcing the CC form gives `c_eff = 6·2 = 12`, a number with
no independent derivation and no relation to any framework central charge. The
construction is **not** a 1+1D CFT interval — it is a degeneracy log of the lowest
`K_α` shell on a (2+1)D wedge. The slope is a boundary **dimension**, not a central
charge. **The Cardy–Calabrese label is a forced analogy and should not be used.**

## (4) F-theorem link: RULED OUT (precise statement)

This is the load-bearing part of the audit, and it agrees with what Paper 50 §5
already proved (Prop 4.3 / `prop:F_vs_S_decomposition`).

- **Different ring.** The F-coefficients live in the master-Mellin M2/M3 ring
  `ℚ·log2 ⊕ ℚ·ζ(3)/π²` (Paper 50 §3, KPS bit-exact). PSLQ of `C_inf` against
  `{1, log2, ζ(3)/π²}` at maxcoeff 10⁸ returns **null** — `C_inf` is **not** in the
  F-theorem ring. Instead `C_inf` lives in the **Boltzmann-base ring** generated by
  `e` (via `coth 1`, `sinh 1`), tied to the integer spacing `Δk = 2` of the `K_α`
  spectrum at β = 2π. It is a thermal/spectral-spacing constant, not a geometric
  CFT transcendental.

- **Different object.** The F-coefficient is extracted from `ρ_W` by the
  **spectral-zeta** operation `−½ζ′_Δ(0)` (a property of the *operator* `Δ`/`D`).
  The wedge entropy is extracted from the **same `ρ_W`** by `−Tr[ρ log ρ]` (a
  property of the *state*). Paper 50's whole point (the three-observables-on-one-
  state structure) is that these two operations on the same modular state extract
  **categorically different** continuum data: the universal F-coefficient (operator
  side) vs. the wedge boundary dimension (state side). This re-read **confirms and
  sharpens** that: the state-side observable is `2 log n_max + C_inf` with the
  boundary dimension 2 (an integer) and a thermal constant `C_inf` (in the e-ring),
  while the operator-side observable is `F_s, F_D` (in the log2/ζ(3) ring). The two
  do not meet.

**Statement for synthesis:** *The BW wedge entropy is the state-side complement of
the F-theorem F-coefficient only in the trivial sense that both are computed from
the same modular state `ρ_W`. They are NOT two readings of one number: the entropy
extracts the boundary dimension `dim(∂W) = 2` (a degeneracy log, e-ring constant),
while the F-coefficient extracts the 3D free energy (log2/ζ(3)/π² ring). The
suspected entanglement↔F-coefficient identity does not hold.*

---

## Curve-fit audit (per `docs/curve_fit_audit_memo.md`)

| Check | Result |
|:---|:---|
| Free parameters | **ZERO** for the asymptotic claim. Slope 2 and constant `coth1−log(2 sinh1)` are *derived* from the closed-form degeneracy + Boltzmann weight, not fitted. The original 2-parameter `(a,b)` fit is replaced by a 0-parameter resummation. |
| Selection bias | None. No search over candidate constants; `C_inf` falls out of the geometric resummation (`A0 = 1/(2 sinh 1)`, `B0/A0 = coth 1`). PSLQ only *confirms*. |
| Alternative fits | log-quadratic had marginally higher R² on the small-n panel, but its quadratic-in-log coefficient → 0 as the window moves up — it was over-fitting small-n curvature. True form is single-log `2 log n + C_inf + O(1/n)`. |
| Robustness | Slope 1.896 → 2.0000 monotone across n ∈ [2,2000]; constant `S−2 log n → C_inf` with sign-correct `O(1/n)` tail (predicted coeff +0.687/n). |
| Independent cross-check | Closed form vs `modular_hamiltonian.py`: diff ≤ 4×10⁻¹⁶. |

---

## Tier classification (curve-fit memo scheme)

- Slope = 2: **Tier A** (derived; the `n²` leading degeneracy forces it).
- `C_inf = coth(1) − log(2 sinh 1)`: **Tier A** (closed-form resummation; PSLQ confirms, does not discover).
- "slope = 2 = dim(S²)": **Tier A** as a dimension statement (degeneracy log of `n_eq ~ n²`); **NOT** a Cardy–Calabrese central-charge claim.
- F-theorem ↔ entropy identity: **ruled out** (PSLQ-disjoint rings; structural reason in Paper 50 Prop 4.3).

---

## Recommended record update (NOT applied here — report for synthesis)

The §3 / F-11 "negative-with-structure" filing is **upgraded** but stays honest:
the area law was correctly rejected; what was missed is that the *replacement*
scaling is exact, not approximate. Suggested reframing of the entry:

> *BW wedge entropy `S(n_max) = 2 log n_max + C_inf + O(1/n)`, slope provably
> exactly 2 (= `dim(∂W) = dim S²`, degeneracy log of the equator shell
> `n_eq = n_max(n_max+1)`), constant `C_inf = coth(1) − log(2 sinh 1) = 0.45845`
> (e-ring, thermal-spacing). Area law `S ~ n²` is correctly rejected; the
> log scaling is exact. NOT a Cardy–Calabrese central-charge log, and NOT the
> entanglement-side complement of the Paper 50 F-coefficient (different ring,
> confirmed by null PSLQ vs `{1, log2, ζ(3)/π²}`; consistent with Paper 50
> Prop 4.3).*

Paper 50 §5 (`sec:wedge_kms`) could optionally be tightened: it currently reports
the small-n fit `S ≈ 1.94 log n + 0.56` and asserts slope "≈ 2". The exact
asymptotic `S = 2 log n + (coth 1 − log(2 sinh 1)) + O(1/n)` is a strictly better
statement (slope *exactly* 2, constant in closed form) and would make the
boundary-dimension claim rigorous rather than fitted. Flagged for PI; not applied.

---

## Meta (for the forcing catalogue)

Door 2's prediction — "a negative tested against one hypothesis hides structure" —
is **confirmed in form but corrected in content**. The shelved positive is real and
exact (closed-form slope + constant), but it is *not* the Cardy–Calabrese / F-theorem
link the catalogue guessed. The disciplined re-read (exact asymptotic instead of a
fit) was decisive: it both *upgraded* the positive (fit → exact) and *killed* the
over-reach (no F-link). This is exactly the value of running the audit before
asserting the identity.
