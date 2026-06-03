# Sprint M4 memo — Paper 50 §7 S⁵ multiplicity correction

**Date:** 2026-06-02
**Sprint:** M4 (Paper 50 §7 multiplicity fix)
**Driver:** `debug/math_sprint_m4_s5_multiplicity.py`
**Raw output:** `debug/data/math_sprint_m4_output.txt`
**Verdict:** CLEAN POSITIVE — bug fully characterized, §7 salvageable with a focused multi-line patch.

---

## 1. Confirmed: 4× factor error in the scalar multiplicity

The audit's M4 finding is verified.

| n | Paper 50 stated mult `(2n+4)(n+1)(n+2)(n+3)/6` | KPS-correct mult `(2n+4)(n+1)(n+2)(n+3)/24` | ratio |
|:-:|:--:|:--:|:--:|
| 0 |  4 |  1 | 4 |
| 1 | 24 |  6 | 4 |
| 2 | 80 | 20 | 4 |
| 3 | 200 | 50 | 4 |
| 4 | 420 | 105 | 4 |
| 5 | 784 | 196 | 4 |

The Paper 50 multiplicity is exactly **4×** the textbook S⁵ Laplacian multiplicity. The textbook formula (Bär 1996; Camporesi–Higuchi 1996) is

$$g_n^{(d)} \;=\; \frac{(2n + d - 1)(n + d - 2)!}{n!\,(d - 1)!}.$$

At $d = 5$, this gives $g_n^{(5)} = (2n+4)(n+1)(n+2)(n+3)/24$, with $g_0 = 1$ and $g_1 = 6$ (the $d+1 = 6$-dimensional vector representation of $SO(6)$).

The Paper 50 statement (L713) drops the $4!$ in the denominator in favour of $3!$ — a single character-level typo (`6` instead of `24`).

## 2. Verified: the downstream F-coefficient error is exactly 4×

Computing $F_s^{S^5} = -\tfrac12 \zeta'_{\Delta_{\text{conf}}}(0)$ via the same Hurwitz-zeta machinery as Paper 50, with each multiplicity:

| Multiplicity | $F_s^{S^5}$ | Reference |
|:--|:--:|:--|
| Paper 50 `/6`  | $-0.02297196544658850\ldots$ | matches Paper 50 §7.1 stated value |
| KPS-correct `/24` | $-0.005742991361647126\ldots$ | matches KPS Table 1 d=5 entry $-5.74 \times 10^{-3}$ |

Ratio: exactly 4.0 (machine precision at 80 dps).

## 3. PSLQ-verified corrected closed form

At 80 dps, PSLQ on the corrected $F_s^{S^5}$ against the basis $\{\log 2,\; \zeta(3)/\pi^2,\; \zeta(5)/\pi^4\}$ returns the integer relation $[256, 2, 2, -15]$ (residual $1.33 \times 10^{-72}$). This rationalizes to:

$$\boxed{F_s^{S^5} \;=\; -\frac{\log 2}{128} \;-\; \frac{\zeta(3)}{128\,\pi^2} \;+\; \frac{15\,\zeta(5)}{256\,\pi^4} \;\approx\; -5.74 \times 10^{-3}.}$$

The structure is identical to Paper 50's stated form, but every numerical coefficient is divided by 4 (Paper 50: $-\tfrac{\log 2}{32} - \tfrac{\zeta(3)}{32\pi^2} + \tfrac{15\,\zeta(5)}{64\pi^4}$).

Equivalently, Paper 50's PSLQ relation `[32, -2, -2, 15]` on the un-renormalized $\zeta'_{S^5}(0)$ should become `[64, -2, -2, 15]` (or equivalently the SAME `[32, -2, -2, 15]` relation, but applied to $\zeta'_{S^5}(0)/4$, i.e. with the correct prefactor $1/64$ rather than $1/16$).

## 4. The Dirac analog (Theorem `thm:dirac_S5`) is CORRECT as written

Paper 50's Dirac multiplicity on S⁵ is $(n+1)(n+2)(n+3)(n+4)/12$. Re-deriving from the Camporesi–Higuchi 1996 spinor multiplicity on round $S^d$ odd,

$$d_n^{\mathrm{Weyl,\,S^d}} \;=\; 2^{[d/2]-1} \cdot \binom{n+d-1}{n},$$

at $d = 5$ this gives $2 \cdot (n+1)(n+2)(n+3)(n+4)/24 = (n+1)(n+2)(n+3)(n+4)/12$ — bit-identical to Paper 50.

Numerical re-derivation:
$$D'_{\mathrm{Dirac}}(0)^{S^5} \;=\; -\frac{3\log 2}{128} \;-\; \frac{5\,\zeta(3)}{128\,\pi^2} \;-\; \frac{15\,\zeta(5)}{256\,\pi^4} \;\approx\; -0.02163$$
matches the Paper 50 closed form to residual $9.9 \times 10^{-83}$.

The Weyl-Dirac sector is unaffected by the scalar `/24` fix.

## 5. Theorem `thm:log3_cancellation` (log 3 cancellation) SURVIVES

Theorem 7.3 is about the cancellation of $(2/3)^s$ terms in the half-integer Hurwitz expansion on the **Dirac side**. The cancellation coefficient sum

$$\frac{1}{12}\left[\frac{81}{16} - \frac{5}{2}\cdot \frac{9}{4} + \frac{9}{16}\cdot 1\right] \;=\; \frac{1}{12}\cdot \frac{81-90+9}{16} \;=\; 0$$

depends only on the Dirac multiplicity polynomial $(u^2 - 9/4)(u^2 - 1/4)/12$. The scalar `/24` fix does not affect the Dirac sector, so the log 3 cancellation theorem stands verbatim.

## 6. Proposition `prop:dual_basis_S5` (dual-basis non-extension) SURVIVES

The structural conclusion — that the $(F_s, F_D)$ pair on $S^5$ cannot project to a single M-engine axis because the M-ring is 3-dimensional but the (scalar, Dirac) plane is 2-dimensional — is preserved under the fix. The DISPLAYED LINEAR-ALGEBRA EQUATIONS need updating, but the inconsistency conclusion is unchanged.

Old (using the WRONG $F_s = 4\times$ the correct value):
$$-4a - 5b = 0 \;\Longrightarrow\; a = -5b/4, \quad 60a - 15b = 0 \;\Longrightarrow\; a = b/4.$$

Corrected (using $F_s^{S^5} = -\log 2/128 - \zeta(3)/(128\pi^2) + 15\zeta(5)/(256\pi^4)$):

- $a F_s + b F_D$ coefficient of $\log 2$: $-(a + 3b)/128$
- $a F_s + b F_D$ coefficient of $\zeta(3)/\pi^2$: $-(a + 5b)/128$
- $a F_s + b F_D$ coefficient of $\zeta(5)/\pi^4$: $15(a - b)/256$

Isolate $\log 2$ (set the other two coefficients to zero):
$$a + 5b = 0 \;\Longrightarrow\; a = -5b, \quad a - b = 0 \;\Longrightarrow\; a = b.$$
Inconsistent unless $b = 0$. **Non-extension holds.**

The same kind of inconsistency arises for any other choice of "isolated" axis; the proposition statement is unchanged in substance.

## 7. The §7 header sentence (L709) does NOT need framing changes

> "The bit-exact match of Sec.~\ref{sec:partition} extends to round $S^{5}$ via the same Hurwitz-zeta machinery."

This sentence remains TRUE after the fix. The corrected closed form `-log(2)/128 - zeta(3)/(128 pi^2) + 15 zeta(5)/(256 pi^4)` IS bit-exact to KPS Table 1 d=5 (and to Beccaria–Tseytlin 2017 review eq 2.6). The PSLQ residual is $1.3 \times 10^{-72}$. The bit-exact framing extends faithfully to $S^5$ — Paper 50's claim was just stated with the wrong numerical scale.

## 8. Proposed paper patches (locked-in)

| Loc | Old | New |
|:--|:--|:--|
| L713 | scalar multiplicity `(2n+4)(n+1)(n+2)(n+3)/6` | `(2n+4)(n+1)(n+2)(n+3)/24` |
| Eq. `eq:Fs_S5` (L727) | `-log(2)/32 - zeta(3)/(32 pi^2) + 15 zeta(5)/(64 pi^4)`, `~ -0.02297` | `-log(2)/128 - zeta(3)/(128 pi^2) + 15 zeta(5)/(256 pi^4)`, `~ -5.74e-3` |
| L739 | "re-indexed multiplicity structure $v^2(v^2-1)/3$ at $v = n+2$" | "re-indexed multiplicity structure $v^2(v^2-1)/12$ at $v = n+2$" |
| L743 | prefactor `2/3` (= 2 * 1/3 from `/6 -> /3`) on Hurwitz derivatives, and `1/3` on the geometric tail sum | `2/12 = 1/6` on derivatives, `1/12` on the tail sum |
| L749 | `log(2)/16 + zeta(3)/(16 pi^2) - 15 zeta(5)/(32 pi^4)` (i.e. $\zeta'_{S^5}(0)$ value implied by `/6`) | `log(2)/64 + zeta(3)/(64 pi^2) - 15 zeta(5)/(128 pi^4)` (`/4` of the above; same PSLQ relation `[32, -2, -2, 15]` applies to `4·` the corrected value) |
| L833 | `-4a - 5b = 0 => a = -5b/4` | `a + 5b = 0 => a = -5b` (coefficient is `-(a + 5b)/128`, the `-4` and `-5` should both be `-1`, ratio 1:5 preserved as 1:5, but the displayed eqs cited `-4`/`-5` reflect the WRONG `/6` mult) |
| L834 | `60a - 15b = 0 => a = b/4` | `15a - 15b = 0 => a = b` (the `60` came from `15 * 4` after the wrong `/6` mult; the correct coefficient is `15(a - b)/256`, ratio 15:15 = 1:1) |

**No change needed** to: L709 ("bit-exact match extends to S⁵" header), Theorem `thm:dirac_S5` (correct as written), Theorem `thm:log3_cancellation` (Dirac side, unaffected), Proposition `prop:dual_basis_S5` structural statement (only the displayed equations).

## 9. Decision gate

**SALVAGEABLE with a focused patch.** The Hurwitz-zeta machinery is correct, the Dirac sector is correct, the structural theorems (`thm:log3_cancellation`, `prop:dual_basis_S5`) survive, the bit-exact-to-KPS framing of §7 is preserved. Only the scalar multiplicity prefactor (`/6` → `/24`) and its three immediate numerical downstream values need updating. Estimated patch size: ~7 LaTeX-line edits, no §7 restructure.

---

## Files

- `debug/math_sprint_m4_s5_multiplicity.py` — driver
- `debug/data/math_sprint_m4_output.txt` — full numerical output
- `debug/math_sprint_m4_memo.md` — this memo

## Cross-references

- Paper 50 §7.1, file `papers/group1_operator_algebras/paper_50_cft3_partition_function.tex`, lines 706–810
- KPS 2011, arXiv:1105.4598 (Klebanov–Pufu–Safdi, "F-Theorem without Supersymmetry"), Table 1 d=5
- Camporesi–Higuchi 1996, gr-qc/9505009, Eq. (3.14) (Dirac spinor multiplicity on S^d)
- Beccaria–Tseytlin 2017 review for S⁵ scalar F-coefficient closed form
- Audit: `debug/audit_critical_issues_review.md` §M4; `debug/review_paper50.md`
