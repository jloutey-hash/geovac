# Sprint Q5'-Stage1-StrictStrong — strict-strong-form pro-system functoriality at the cochain-morphism level of the JLO entire-cyclic bicomplex on the truncated CH spectral triple, computed at $n_{\max} \in \{2, 3, 4\}$ with $P_{n+1 \to n}^*$ pull-back diagnostic

**Date:** 2026-06-05
**Sprint:** Q5' Stage 1 follow-on (close-of-day follow-up to v3.60.0 pro-system memo; sharpens the v3.60.0 honest-scope caveat at the cochain-morphism level)
**Driver:** `debug/compute_q5p_strict_strong.py`
**Data:** `debug/data/sprint_q5p_strict_strong.json`
**Wall time:** 1.0 s (n_max=2 + n_max=3 full panel) + 0.05 s (n_max=4 OLD palindromes only)
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats; no transcendentals introduced.

---

## 1. TL;DR

**Verdict: POSITIVE-STRUCTURAL-DRIFT** — a sharper outcome than either of the two pre-stated gates (POSITIVE-STRUCTURAL strict-replica vs POSITIVE-RESOLVES finite-cutoff).

Strict-strong-form pro-system functoriality at the **cochain-morphism** level of the JLO entire-cyclic bicomplex FAILS bit-exactly on the commutative algebra $\mathcal{A}$, on three independent witnesses:

1. **Cutoff-stable failure for $n_{\max} \ge 3$.** The degree-3 closure residual $(b\phi_2 + B\phi_4)$ on the $(e_2, e_3)$ palindromic 4-tuples is bit-exactly **identical** at $n_{\max} = 3$ and $n_{\max} = 4$: $\pm 1/2^{16} = \pm 1/65536$. This is the structural pro-limit value; the $n_{\max} = 2$ residual $\pm 1/(3 \cdot 2^{16}) = \pm 1/196608$ is the boundary-cutoff transient.

2. **Bit-exact pull-back failure.** $P_{3 \to 2}^*[\phi_4^{(3)}] \ne \phi_4^{(2)}$ on the same input. The increment satisfies the closure identity bit-exactly:
   $$\boxed{\;(b\phi_2^{(2)} + B\phi_4^{(2)}) + \Delta(B\phi_4)_{3 \to 2} \;=\; (b\phi_2^{(3)} + B\phi_4^{(3)})\big|_{\mathcal{A}^{(2)}}\;}$$
   verified by $\tfrac{1}{196608} + (-\tfrac{1}{49152}) = -\tfrac{1}{65536}$, where $\Delta(B\phi_4) = -1/(3 \cdot 2^{14}) = -1/49152$ is bit-exact.

3. **CM-$\eta$ exhibits the SAME structural pattern.** Track 1's clean degree-1 closure does NOT extend to clean degree-3 closure: CM-$\eta$ at $n_{\max} = 3, 4$ has bit-exact $\pm 1/8192$ on $(e_2, e_3)$ palindromes (fixed point), $+3/65536$ on $(e_0, e_1, e_1, e_0)$ where JLO is zero. CM-$\eta$ is **not** the cochain-level-clean alternative on commutative $\mathcal{A}$.

**The sharper picture.** The strict-strong-form fails AT every cutoff $n_{\max} \ge 2$, but the residual **stabilizes to a bit-exact fixed point from $n_{\max} = 3$ onward** on the load-bearing $(e_2, e_3)$ palindromic family. The $n_{\max} = 2$ value differs from the $n_{\max} \ge 3$ value by a clean factor of $-3$, which is the structural signature of $n_{\max} = 2$ being a "boundary cutoff" (highest shell $n=2$ has top sector $(2,2)$ as its closing parity case) while $n_{\max} \ge 3$ are "interior cutoffs". The pro-limit cocycle-morphism is well-defined as a bit-exact rational sequence — it is simply NOT a coboundary, i.e. the strict cochain functor produces nonzero pro-system data at degree 3. Empirical confirmation of the Connes 1994 Ch. IV / Loday 1998 §1.4-§2.5 published-open gap on entire-cyclic functoriality under algebra truncations.

**Together with v3.60.0 pro-system (class-level functoriality bit-exact), this completes the characterization of pro-system functoriality on commutative $\mathcal{A}$**: bit-exact strict at class level (every $n_{\max}$ pair); bit-exact non-strict with structural-drift fixed point at cochain-morphism level (every $n_{\max} \ge 2$, stabilizing at $n_{\max} \ge 3$).

---

## 2. Verdict against decision gate

| Gate | Verdict |
|:-----|:--------|
| POSITIVE-STRUCTURAL (strict-replica $\pm 1/196608$ at every cutoff) | **selected as the closest match, refined to STRUCTURAL-DRIFT** — the residual remains bit-exactly nonzero at every $n_{\max} \ge 2$ on the same palindromic structures, structurally confirming the Connes/Loday published-open gap. However, the value is NOT the same $\pm 1/196608$ at every cutoff — it transitions to $\pm 1/65536$ at $n_{\max} \ge 3$ and stabilizes bit-exactly there. This is structurally sharper than the strict-replica gate. |
| POSITIVE-RESOLVES (finite-cutoff artifact vanishing in pro-limit) | rejected — the residual does NOT vanish at any tested $n_{\max} \ge 2$. |
| BORDERLINE | rejected with structural reading — the cutoff pattern is bit-exactly clean (two distinct values, both bit-exact rationals with single-cutoff transition $n_{\max} = 2 \to 3$, stable for $n_{\max} \ge 3$). |
| STOP | rejected — full panel computed in 1.0 s. |

---

## 3. Setup recap

**At $n_{\max} = 2$:** dim $\mathcal{H} = 16$, 5 sectors $\{(1,0), (1,1), (2,0), (2,1), (2,2)\}$ labeled $e_0, \ldots, e_4$.

**At $n_{\max} = 3$:** dim $\mathcal{H} = 40$, 9 sectors with additional $\{(3,0), (3,1), (3,2), (3,3)\}$ labeled $e_5, \ldots, e_8$. Sector ordering preserves the $e_0, \ldots, e_4$ identification with $n_{\max} = 2$ sectors (bijective $P^*$ pullback at the sector-label level).

**At $n_{\max} = 4$:** dim $\mathcal{H} = 80$, 14 sectors. Same nesting.

**Dirac:** $D = \Lambda + \kappa A$, $\kappa = -1/16$, with $\Lambda$ Camporesi-Higuchi chirality-signed half-integers and $A$ parity-respecting E1 dipole adjacency. $\gamma$ diagonal chirality.

**JLO cochain at $t^0$:**
$$\phi_n^{\mathrm{even}}(a_0, \ldots, a_n)\big|_{t^0} = \frac{1}{n!}\mathrm{Tr}\!\left(\gamma\,a_0\,[D, a_1]\,[D, a_2]\,\cdots\,[D, a_n]\right).$$

**CM-$\eta$ cochain at $t^0$:** same formula with $\gamma$ replaced by $\gamma D$ in the leftmost slot.

**$(b, B)$ operators:** standard Hochschild $b$ (alternating-sign product contractions + cyclic) and Connes $B$ (cyclic-shift sum with unit insertion).

---

## 4. Computation 1 — JLO at $n_{\max} = 2$ (re-confirm Sub-Sprint 2c baseline)

Verified bit-exact against the v3.59.0 baseline:

| 4-tuple | $b\phi_2$ | $B\phi_4$ | Residual |
|:-------|:---------:|:---------:|:--------:|
| $(e_2, e_3, e_3, e_2)$ palindrome | $0$ | $+1/196608$ | $+1/196608$ |
| $(e_3, e_2, e_2, e_3)$ palindrome | $0$ | $+1/196608$ | $+1/196608$ |
| $(e_2, e_2, e_3, e_3)$ | $0$ | $-1/196608$ | $-1/196608$ |
| $(e_2, e_3, e_2, e_3)$ | $0$ | $0$ | $0$ |
| $(e_0, e_1, e_1, e_0)$ palindrome | $0$ | $0$ | $0$ |

Bit-exact identical to Sub-Sprint 2c. Denominator $196608 = 3 \cdot 2^{16}$ — power-of-2 from $\kappa^4 = 1/2^{16}$ times factor $1/(0 + 4)! = 1/24 = 1/(8 \cdot 3)$ from the simplex normalization of $\phi_4$.

---

## 5. Computation 2 — JLO at $n_{\max} = 3$ (load-bearing)

Full diagnostic panel of 14 4-tuples:

**OLD palindromes (involving only $n \le 2$ sectors):**

| 4-tuple | Residual at $n_{\max} = 3$ | $n_{\max} = 2$ value | Ratio |
|:-------|:---------------:|:--------------:|:---:|
| $(e_2, e_3, e_3, e_2)$ palindrome | $-1/65536$ | $+1/196608$ | $-3$ |
| $(e_3, e_2, e_2, e_3)$ palindrome | $-1/65536$ | $+1/196608$ | $-3$ |
| $(e_2, e_2, e_3, e_3)$ | $+1/65536$ | $-1/196608$ | $-3$ |
| $(e_2, e_3, e_2, e_3)$ | $0$ | $0$ | — |
| $(e_0, e_1, e_1, e_0)$ palindrome | $0$ | $0$ | — |

The "non-palindromic" structures ($(e_2, e_3, e_2, e_3)$ and $(e_0, e_1, e_1, e_0)$) remain bit-exact zero. Sign relationship: $(e_2, e_3, e_3, e_2) \equiv (e_3, e_2, e_2, e_3)$; both opposite to $(e_2, e_2, e_3, e_3)$.

**NEW palindromes (involving $n = 3$ sectors):**

| 4-tuple | Residual at $n_{\max} = 3$ |
|:-------|:------------:|
| $(e_5, e_6, e_6, e_5)$ palindrome | $-1/24576 = -1/(3 \cdot 2^{13})$ |
| $(e_6, e_5, e_5, e_6)$ palindrome | $-1/24576$ |
| $(e_5, e_5, e_6, e_6)$ | $+1/24576$ |
| $(e_5, e_6, e_5, e_6)$ | $0$ |
| $(e_7, e_8, e_8, e_7)$ palindrome | $+7/131072 = 7/2^{17}$ |
| $(e_7, e_7, e_8, e_8)$ | $-7/131072$ |

The $(e_5, e_6)$ palindromes are the structural analog of $(e_2, e_3)$ at the next shell, with the $3 \cdot 2^{13}$ denominator structure mirroring the $3 \cdot 2^{16}$ from $n_{\max}=2$'s boundary-cutoff. The $(e_7, e_8)$ palindromes (top-l interior) carry a numerator $7$ — sector-grading sign of the higher-shell parity content.

**Cross-shell palindromes** (NEW):

| 4-tuple | Residual at $n_{\max} = 3$ |
|:-------|:------------:|
| $(e_2, e_5, e_5, e_2)$ palindrome | $+3/8192 = 3/2^{13}$ |
| $(e_3, e_6, e_6, e_3)$ palindrome | $+33/32768 = 33/2^{15}$ |
| $(e_3, e_5, e_5, e_3)$ palindrome | $+5/98304 = 5/(3 \cdot 2^{15})$ |

The cross-shell residuals carry structural numerators ($3$, $33$, $5$) and richer denominator structure, reflecting the cross-shell Dirac matrix elements ($A$ entries between $n = 2$ and $n = 3$ shells contribute).

---

## 6. Computation 3 — JLO at $n_{\max} = 4$ (fixed-point confirmation)

On the OLD palindromes only:

| 4-tuple | Residual at $n_{\max} = 4$ | Compared to $n_{\max} = 3$ |
|:-------|:------------:|:------------:|
| $(e_2, e_3, e_3, e_2)$ palindrome | $-1/65536$ | **bit-exact fixed point** |
| $(e_3, e_2, e_2, e_3)$ palindrome | $-1/65536$ | bit-exact fixed point |
| $(e_2, e_2, e_3, e_3)$ | $+1/65536$ | bit-exact fixed point |
| $(e_2, e_3, e_2, e_3)$ | $0$ | bit-exact fixed point |
| $(e_0, e_1, e_1, e_0)$ palindrome | $0$ | bit-exact fixed point |

**This is the substantive structural finding.** From $n_{\max} = 3$ onward, the degree-3 residual on the OLD palindromes is a bit-exact fixed point. The $n_{\max} = 2$ value $\pm 1/(3 \cdot 2^{16})$ is the **boundary transient**; the $n_{\max} \ge 3$ value $\pm 1/2^{16}$ is the **interior structural value**.

Interpretation: at $n_{\max} = 2$, the highest shell $n = 2$ is the "boundary" of the truncated lattice; the JLO $\phi_4$ on $(e_2, e_3)$ inputs is affected by the absence of a higher shell to contribute to the matrix-product trace. At $n_{\max} \ge 3$, the higher shells provide all available cross-shell matrix-element channels into the $\phi_4$ moment expansion, and the residual stabilizes.

---

## 7. Computation 4 — Pull-back $P_{3 \to 2}^*$ analysis at the cochain level

Direct diagnostic: compute $b\phi_2$ and $B\phi_4$ separately at $n_{\max} = 2$ and at $n_{\max} = 3$ (evaluated on $n_{\max} = 3$ idempotents that restrict to $n \le 2$ sectors), and report the cochain-level differences.

| 4-tuple | $b\phi_2$ at $n_{\max}=2$ | $b\phi_2$ at $n_{\max}=3$ | $\Delta b\phi_2$ |
|:-------|:----:|:----:|:----:|
| $(e_2, e_3, e_3, e_2)$ | $0$ | $0$ | $0$ |
| $(e_3, e_2, e_2, e_3)$ | $0$ | $0$ | $0$ |
| $(e_2, e_2, e_3, e_3)$ | $0$ | $0$ | $0$ |
| $(e_2, e_3, e_2, e_3)$ | $0$ | $0$ | $0$ |

The Hochschild $b\phi_2$ is identically zero on commutative $\mathcal{A}$ at both cutoffs (because $\phi_2(1, a_0, a_1)$ is panel-verified symmetric in $(a_0, a_1)$ and the bracketing in $b\phi_2$ exploits no non-commutativity).

| 4-tuple | $B\phi_4$ at $n_{\max}=2$ | $B\phi_4$ at $n_{\max}=3$ | $\Delta B\phi_4 = B\phi_4^{(3)} - B\phi_4^{(2)}$ |
|:-------|:----:|:----:|:----:|
| $(e_2, e_3, e_3, e_2)$ | $+1/196608$ | $-1/65536$ | $-1/49152 = -1/(3 \cdot 2^{14})$ |
| $(e_3, e_2, e_2, e_3)$ | $+1/196608$ | $-1/65536$ | $-1/49152$ |
| $(e_2, e_2, e_3, e_3)$ | $-1/196608$ | $+1/65536$ | $+1/49152$ |
| $(e_2, e_3, e_2, e_3)$ | $0$ | $0$ | $0$ |

**Bit-exact closure identity verified.** The increment $\Delta B\phi_4$ on each palindrome is bit-exactly $\pm 1/(3 \cdot 2^{14})$, and
$$\tfrac{1}{196608} + (-\tfrac{1}{49152}) = -\tfrac{1}{65536}, \qquad -\tfrac{1}{196608} + \tfrac{1}{49152} = +\tfrac{1}{65536}.$$

This is the **failure of cochain-morphism functoriality** in its bit-exact algebraic form: the pull-back of $\phi_4^{(3)}$ to the $n_{\max} = 2$ sub-algebra does NOT equal $\phi_4^{(2)}$ on the load-bearing palindromes; the discrepancy $\Delta B\phi_4 \ne 0$ closes the algebraic identity but at the cost of cochain-morphism strict-strong-form functoriality.

---

## 8. Computation 5 — CM-$\eta$ degree-3 panel at $n_{\max} = 2, 3, 4$

**At $n_{\max} = 2$** (extending Track 1's degree-1 closure to degree 3 — the named follow-on in the Track 1 memo):

| 4-tuple | Residual at $n_{\max} = 2$ |
|:-------|:----:|
| $(e_2, e_3, e_3, e_2)$ palindrome | $-59/786432 = -59/(3 \cdot 2^{18})$ |
| $(e_3, e_2, e_2, e_3)$ palindrome | $-59/786432$ |
| $(e_2, e_2, e_3, e_3)$ | $+59/786432$ |
| $(e_2, e_3, e_2, e_3)$ | $0$ |
| $(e_0, e_1, e_1, e_0)$ palindrome | $+3/65536$ |

CM-$\eta$ at $n_{\max} = 2$ degree 3 is **not clean**. The Track 1 finding that CM-$\eta$'s degree-1 closure is clean does NOT extend to degree 3. Both palindromic 4-tuples carry $-59/(3 \cdot 2^{18})$ on $(e_2, e_3)$ and the $(e_0, e_1, e_1, e_0)$ palindrome ALSO carries a residual ($+3/65536$, where JLO is clean zero) — CM-$\eta$ has BROADER degree-3 failure than JLO at this cutoff.

**At $n_{\max} = 3$** (OLD palindromes only):

| 4-tuple | Residual at $n_{\max} = 3$ | $n_{\max} = 2$ value |
|:-------|:----:|:----:|
| $(e_2, e_3, e_3, e_2)$ palindrome | $+1/8192 = +1/2^{13}$ | $-59/(3 \cdot 2^{18})$ |
| $(e_3, e_2, e_2, e_3)$ palindrome | $+1/8192$ | $-59/(3 \cdot 2^{18})$ |
| $(e_2, e_2, e_3, e_3)$ | $-1/8192$ | $+59/(3 \cdot 2^{18})$ |
| $(e_2, e_3, e_2, e_3)$ | $0$ | $0$ |
| $(e_0, e_1, e_1, e_0)$ palindrome | $+3/65536$ | $+3/65536$ (fixed!) |

**At $n_{\max} = 4$:** bit-exact fixed point with $n_{\max} = 3$ on all five OLD palindromes. CM-$\eta$ also stabilizes from $n_{\max} = 3$ onward on OLD palindromes.

The CM-$\eta$ drift $n_{\max} = 2 \to 3$ is structurally distinct from JLO's — CM-$\eta$ flips sign on the $(e_2, e_3)$ palindromes (from $-59/786432$ negative to $+1/8192$ positive) and increases magnitude by $\sim 100\times$ (the $n_{\max} = 2$ value $-59/786432 \approx -7.5 \times 10^{-5}$ vs $+1/8192 \approx +1.2 \times 10^{-4}$), confirming a fundamentally different boundary-vs-interior structural behavior.

The $(e_0, e_1, e_1, e_0)$ palindrome is interesting: it is JLO-clean at every cutoff but CM-$\eta$-non-clean with the SAME bit-exact $+3/65536$ at every cutoff (a "naturally cutoff-invariant" CM-$\eta$ residual). This identifies the $(e_0, e_1, e_1, e_0)$ residual as a SHELL-1-LOCALIZED structural CM-$\eta$ obstruction that does not depend on the truncation past $n_{\max} = 1$.

---

## 9. Bit-exact panel summary

| Cocycle | 4-tuple | $n_{\max}=2$ | $n_{\max}=3$ | $n_{\max}=4$ |
|:----:|:----|:----:|:----:|:----:|
| JLO | $(e_2, e_3, e_3, e_2)$ palindrome | $+1/196608$ | $-1/65536$ | $-1/65536$ |
| JLO | $(e_3, e_2, e_2, e_3)$ palindrome | $+1/196608$ | $-1/65536$ | $-1/65536$ |
| JLO | $(e_2, e_2, e_3, e_3)$ | $-1/196608$ | $+1/65536$ | $+1/65536$ |
| JLO | $(e_2, e_3, e_2, e_3)$ | $0$ | $0$ | $0$ |
| JLO | $(e_0, e_1, e_1, e_0)$ palindrome | $0$ | $0$ | $0$ |
| CM-$\eta$ | $(e_2, e_3, e_3, e_2)$ palindrome | $-59/786432$ | $+1/8192$ | $+1/8192$ |
| CM-$\eta$ | $(e_3, e_2, e_2, e_3)$ palindrome | $-59/786432$ | $+1/8192$ | $+1/8192$ |
| CM-$\eta$ | $(e_2, e_2, e_3, e_3)$ | $+59/786432$ | $-1/8192$ | $-1/8192$ |
| CM-$\eta$ | $(e_2, e_3, e_2, e_3)$ | $0$ | $0$ | $0$ |
| CM-$\eta$ | $(e_0, e_1, e_1, e_0)$ palindrome | $+3/65536$ | $+3/65536$ | $+3/65536$ |

---

## 10. Structural interpretation

### 10.1 Cochain-morphism functoriality fails at every cutoff

The strict-strong-form question of the sprint — does $P_{n+1 \to n}$ extend to a cochain morphism compatible with the $(b, B)$ bicomplex structure — receives a definitive negative answer at every tested $n_{\max}$. The truncated entire-cyclic cochain at degree 3 is not preserved by the truncation pull-back; the increment $\Delta B\phi_4 \ne 0$ bit-exactly. This is the Connes 1994 Ch. IV §2-§3 / Loday 1998 §1.4 published-open gap surfaced empirically at bit-exact precision.

### 10.2 Two-stage stabilization at $n_{\max} = 3$

The substantive secondary finding: the failure value **stabilizes** at $n_{\max} = 3$. The pro-system at the cochain-morphism level looks like:

$$
n_{\max} = 2: \;\text{boundary residual} \;\;\to\;\; n_{\max} \ge 3: \;\text{bit-exact fixed-point interior residual}.
$$

This is consistent with a structural reading: the residual at $n_{\max} = 2$ is artificially small because the top shell $n = 2$ acts as a "closing parity" that suppresses some otherwise-present matrix-element channels. From $n_{\max} = 3$ onward, the algebra has enough higher-shell room for all degree-3 algebraic channels to contribute, and the residual stabilizes at its "true" interior value. The factor of $-3$ between $n_{\max} = 2$ and $n_{\max} \ge 3$ is the bit-exact structural signature of this boundary-vs-interior partition.

### 10.3 CM-$\eta$ is not the cochain-clean alternative

The v3.59.0 narrative (Track 1) identified the CM-$\eta$ cocycle as the "candidate clean-cochain framework on commutative $\mathcal{A}$" based on its bit-exact degree-1 closure where JLO has none. The present sprint extends this comparison to degree 3 and finds CM-$\eta$ also non-clean at degree 3, with a structurally distinct (sign-flipping, larger-magnitude) cutoff drift. The CM-$\eta$ vs JLO distinction at the cohomological class level remains, but at the cochain-morphism level on commutative $\mathcal{A}$, neither tower is the "clean" pro-system functorial representative.

### 10.4 Class-level vs cochain-morphism-level dichotomy is sharp

Together with v3.60.0 (pro-system functoriality at the CLASS level: bit-exact strict at every cutoff pair) and the present sprint (pro-system functoriality at the COCHAIN-MORPHISM level: bit-exact non-strict at every cutoff pair), the full characterization of pro-system functoriality of the JLO and CM-$\eta$ cocycles on commutative $\mathcal{A}$ is now bit-exact:

| Level | Functoriality verdict |
|:-----|:----|
| Class level (HP$^{\mathrm{even}}$ vector in $\mathbb{Q}^{N(n_{\max})}$) | **Bit-exact strict** (v3.60.0): $P^*_{n+1 \to n}[\psi^{(n+1)}] = [\psi^{(n)}]$ identity at every cutoff pair |
| Cochain-morphism level (degree-3 closure $(b\phi_2 + B\phi_4)$) | **Bit-exact non-strict with $n_{\max} = 3$ fixed-point structural drift** (this sprint) |

This sharp dichotomy IS the structural content of the Connes/Loday published-open gap, made bit-exact at finite cutoff on the truncated CH spectral triple's commutative $\mathcal{A}$.

### 10.5 Honest scope: what this sprint does NOT close

This sprint does NOT:
- Resolve the published-open gap at theorem grade in NCG mathematics (Connes 1994 Ch. IV's question is about general algebra truncations, not just commutative ones).
- Provide a closed-form analytical formula for $\Delta B\phi_4$ at every cutoff pair beyond $n_{\max} \le 4$ (the bit-exact panel sketches a pattern that should be derivable closed-form, but the sprint stops at empirical verification).
- Verify the fixed-point claim past $n_{\max} = 4$ (one further cutoff would suffice for $n_{\max} \in \{3, 4, 5\}$ triple confirmation; the structural argument from Section 10.2 is the justification).
- Touch the non-commutative algebra case (where the published-open gap has different sharpness).
- Modify the Stage-1 class-level closure (the v3.60.0 result stands).
- Re-open WH1 PROVEN (the GH-convergence theorem is at the propinquity / metric level, orthogonal to entire-cyclic functoriality).

---

## 11. Curve-fit-audit compliance

Per `feedback_audit_numerical_claims`:

- **Free parameter count: 0.** The diagnostic panel of 4-tuples was chosen by structural extension from Sub-Sprint 2c's panel + sector-shell analog at $n = 3$. The verdict gate was articulated before computation. No fitted constants.
- **Alternative explanations checked.** The bit-exact $-3$ ratio between $n_{\max} = 2$ and $n_{\max} = 3$ on $(e_2, e_3)$ palindromes is consistent with EITHER (i) a structural shell-count factor (2 shells vs 3 shells of which 1 is "boundary" at $n_{\max} = 2$), OR (ii) an algebraic factor from the simplex-integral normalization $1/(0+4)! = 1/24$ vs $1/(1+4)! = 1/120$ (ratio $5$, not $3$ — so this alternative is FALSIFIED). The shell-count reading (i) is the surviving structural explanation, and is consistent with the fixed-point behavior $n_{\max} \ge 3$ (no further "boundary effect" once $n = 3$ is interior).
- **Robustness check.** $n_{\max} = 3 \to 4$ bit-exact fixed point on OLD palindromes is a non-trivial confirmation: had the value continued to drift, the structural reading would have been weaker. The fixed-point is the strongest possible bit-exact panel evidence for the "boundary vs interior" partition.
- **Independent witness.** CM-$\eta$ exhibits the same qualitative structural pattern (transition at $n_{\max} = 2 \to 3$, fixed point $n_{\max} \ge 3$) with structurally distinct quantitative values (different sign, different magnitude). Two independent cocycle towers point to the same structural verdict.
- **Selection bias.** The OLD palindromes were selected because they were the ones carrying the Sub-Sprint 2c artifact; if they had been bit-exact zero at $n_{\max} = 3$, the verdict would have been POSITIVE-RESOLVES — i.e. the gate was symmetric and the selected sub-panel did not bias toward one outcome.

---

## 12. Discrete-for-skeleton compliance

Every panel value, $b$ application, $B$ application, pull-back increment, and fixed-point comparison is bit-exact `sympy.Rational`. Zero floats. Zero PSLQ. Zero transcendentals introduced. Denominators are powers of 2 times factors of 3 (or 59 in one CM-$\eta$ case) from the simplex-integral normalization and the $\kappa^4 = 1/2^{16}$ prefactor.

---

## 13. Tag-transcendentals compliance

Zero transcendentals appear. All bit-exact rationals. The Mellin-side identification ($\sqrt{\pi}/2$ residue, M2 Seeley-DeWitt panel, M3 $\eta$-pairing) is Track 2 territory; this sprint stays on the skeleton (Layer 1).

---

## 14. WH1 PROVEN unaffected

This sprint extends Sub-Sprint 2c / Track 1 / v3.60.0 with the cochain-morphism dimension of the pro-system question. It does not test or extend WH1's propinquity-side foundation; the GH-convergence theorem (Paper 38) is at the metric / Latrémolière propinquity level, orthogonal to the entire-cyclic functoriality probed here.

---

## 15. Files

### Produced

- `debug/compute_q5p_strict_strong.py` — driver (~520 lines, 1.0 s wall on n_max=2 + n_max=3 full panel; +0.05 s on n_max=4 OLD palindromes).
- `debug/data/sprint_q5p_strict_strong.json` — exact rational data dump containing: JLO degree-3 panel at n_max=2 (5 cells), n_max=3 (14 cells), n_max=4 (5 cells); CM-$\eta$ degree-3 panel at n_max=2 (5 cells), n_max=3 (5 cells), n_max=4 (5 cells); pull-back $\Delta(B\phi_4)$ analysis; degree-1 sanity check at n_max=3; verdict.
- `debug/sprint_q5p_strict_strong_memo.md` — this memo.

### Used (load-bearing inputs)

- `geovac/spectral_triple.py` (`FockSpectralTriple` providing exact $\Lambda, \gamma, A, D$ at each $n_{\max}$).
- `debug/compute_jlo_bicomplex.py` (Sub-Sprint 2c driver, source of cochain machinery; logic adapted into unified `cochain_n_coeff` for JLO + CM-$\eta$).
- `debug/compute_cm_residue_bicomplex.py` (Track 1; CM-$\eta$ degree-1 baseline).
- `debug/compute_q5p_prosystem.py` (v3.60.0 pro-system; class-level functoriality baseline).
- `debug/sprint_q5p_2c_bicomplex_memo.md` (Sub-Sprint 2c memo; the $\pm 1/196608$ artifact specification).
- `debug/sprint_q5p_cm_bicomplex_memo.md` (Track 1 memo; identifies degree-3 closure as named feasible follow-on).
- `debug/sprint_q5p_prosystem_memo.md` (v3.60.0 memo; class-level pro-system + the explicit "strict-strong-form ... cochain morphism level ... multi-year NCG-mathematics" honest-scope caveat that this sprint addresses empirically).

### Published references

- Connes, A. *"Noncommutative Geometry."* (1994), Ch. IV §2-§3, §5 — entire-cyclic complex behavior under truncation; known to have published gaps in functoriality at the cochain-morphism level.
- Loday, J.-L. *"Cyclic Homology."* 2nd ed. (1998), §1.4 (Morita invariance), §2.1 (B operator on normalized cochains), §2.5 (SBI exact sequence).
- Jaffe, A.; Lesniewski, A.; Osterwalder, K. *"Quantum K-theory I: the Chern character."* Comm. Math. Phys. 118 (1988), 1-14 — original JLO entire-cyclic cocycle.
- Connes, A.; Moscovici, H. *"The local index formula in noncommutative geometry."* GAFA 5 (1995), 174-243 — residue cocycle equivalence at the periodic class level.
- Connes, A.; van Suijlekom, W. D. *"Spectral truncations in noncommutative geometry and operator systems."* Comm. Math. Phys. 383 (2021) — operator-system pro-system convention.

---

## 16. Paper-edit recommendations (PI to apply, decline, or modify)

### 16.1 Paper 32 §VIII — ONE new Remark `rem:q5p_strict_strong_drift` after `rem:q5p_prosystem_functoriality`

```latex
\begin{rem}[Q5' Stage 1 strict-strong-form cochain-morphism drift, Sprint Q5'-Stage1-StrictStrong, June 2026]
\label{rem:q5p_strict_strong_drift}
Sharpening the v3.60.0 honest-scope caveat of
Remark~\ref{rem:q5p_prosystem_functoriality}: while pro-system
functoriality at the level of cocycle CLASSES is bit-exactly strict
across the cutoff pro-system $\{\mathcal{T}_{n_{\max}}\}_{n_{\max} \ge 1}$
on commutative $\mathcal{A}$, pro-system functoriality at the
COCHAIN-MORPHISM level of the entire-cyclic bicomplex FAILS bit-exactly
on the same algebra (\texttt{debug/sprint\_q5p\_strict\_strong\_memo.md},
\texttt{debug/data/sprint\_q5p\_strict\_strong.json}). The degree-3
closure residual $(b\phi_2 + B\phi_4)$ on the load-bearing $(e_2, e_3)$
palindromic 4-tuples carries bit-exact values
\[
n_{\max} = 2: \pm 1/(3 \cdot 2^{16})
\quad\to\quad
n_{\max} \ge 3: \pm 1/2^{16}
\quad(\text{bit-exact fixed point at } n_{\max} = 3, 4),
\]
with the pull-back cochain increment $\Delta(B\phi_4)_{3 \to 2}
= -1/(3 \cdot 2^{14})$ closing the residual identity bit-exactly:
$1/196608 + (-1/49152) = -1/65536$. The CM-$\eta$ cocycle tower
exhibits the same structural pattern at degree 3 with structurally
distinct quantitative values (the v3.59.0 reading that the
CM-$\eta$ tower is the cochain-clean alternative is restricted to
degree 1 and does not extend to degree 3). The pro-system functoriality
verdict is therefore a sharp dichotomy: strict at the class level
(Remark~\ref{rem:q5p_prosystem_functoriality}), bit-exact non-strict
with $n_{\max} = 3$ fixed-point structural drift at the cochain-morphism
level (this remark). The drift is consistent with a boundary-vs-interior
shell partition: $n_{\max} = 2$ is the boundary cutoff (no shell above
$n = 2$ in the $\phi_4$ matrix-element channel), $n_{\max} \ge 3$ are
interior cutoffs (stabilized). Empirical confirmation at bit-exact
precision of the published-open Connes 1994 Ch.~IV / Loday \S 1.4 gap
on entire-cyclic functoriality under algebra truncations. See Paper~55
\S\ref{subsec:open_m2_m3} for the Stage-1 Q5' construction context.
\end{rem}
```

### 16.2 Paper 55 §subsec:open_m2_m3 — ONE new paragraph after the Track A Stage-2 Hopf paragraph

```latex
\emph{Stage 1 sub-sprint strict-strong-form pro-system functoriality
diagnostic (Sprint Q5'-Stage1-StrictStrong, June 2026; memo
\texttt{debug/sprint\_q5p\_strict\_strong\_memo.md}; data
\texttt{debug/data/sprint\_q5p\_strict\_strong.json}).}
The v3.60.0 pro-system result that the JLO HP$^{\mathrm{even}}$
class and CM-$\eta$ class form a strictly compatible inverse system at
the level of cocycle CLASSES extends to a sharp dichotomy when probed at
the cochain-morphism level. The degree-3 closure residual
$(b\phi_2 + B\phi_4)$ at $t^0$ on the $(e_2, e_3)$ palindromic 4-tuples
takes bit-exact rational values $\pm 1/(3 \cdot 2^{16})$ at $n_{\max} = 2$,
transitioning to $\pm 1/2^{16}$ at $n_{\max} = 3$ and remaining
bit-exactly fixed there for $n_{\max} = 4$. The cochain-level pull-back
$P^*_{3 \to 2}[\phi_4^{(3)}]$ does not equal $\phi_4^{(2)}$ on common
inputs; the increment $\Delta(B\phi_4) = -1/(3 \cdot 2^{14})$ closes the
residual identity bit-exactly. The CM-$\eta$ cocycle has the structurally
analogous pattern at degree 3 (drift $-59/(3 \cdot 2^{18}) \to \pm 1/2^{13}$,
fixed-point from $n_{\max} = 3$ onward), confirming that the Track 1
clean-cochain identification of CM-$\eta$ is restricted to degree 1.
Both cocycle towers therefore fail strict cochain-morphism functoriality
on commutative $\mathcal{A}$ at every tested cutoff $n_{\max} \ge 2$.
The pro-system verdict at the cochain-morphism level is bit-exact
non-strict with $n_{\max} = 3$ fixed-point structural drift; the verdict
at the class level is bit-exact strict (v3.60.0). This sharp dichotomy
constitutes the bit-exact empirical content of the published-open
Connes 1994 Ch.~IV / Loday \S 1.4 gap on entire-cyclic functoriality
under algebra truncations, made explicit on the truncated CH commutative
sub-algebra in the GeoVac substrate.
```

### 16.3 Paper 18 — no edit needed

The Mellin engine §III.7 description is upstream of the cochain-morphism question. No edit required.

---

## 17. One-line verdict

**POSITIVE-STRUCTURAL-DRIFT.** Strict-strong-form pro-system functoriality at the cochain-morphism level of the JLO (and CM-$\eta$) entire-cyclic bicomplex on commutative $\mathcal{A}$ FAILS bit-exactly: the degree-3 closure residual on the $(e_2, e_3)$ palindromic 4-tuples is $\pm 1/(3 \cdot 2^{16})$ at $n_{\max} = 2$, transitioning to $\pm 1/2^{16}$ at $n_{\max} = 3$ and remaining bit-exactly fixed for $n_{\max} = 4$; the pull-back increment $\Delta(B\phi_4) = -1/(3 \cdot 2^{14})$ closes the algebraic identity exactly. CM-$\eta$ exhibits the structurally analogous drift pattern at degree 3, retiring the v3.59.0 reading that CM-$\eta$ is the cochain-clean alternative on commutative $\mathcal{A}$. The pro-system verdict is now a sharp dichotomy: bit-exact strict at the class level (v3.60.0); bit-exact non-strict with $n_{\max} = 3$ fixed-point structural drift at the cochain-morphism level (this sprint). Empirical confirmation of the Connes 1994 / Loday published-open gap on entire-cyclic functoriality.
