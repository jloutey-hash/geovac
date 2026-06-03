# Sprint M2 — Paper 18 Eq.(38) D(4) substitution error: resolution memo

**Date:** 2026-06-02
**Driver:** `debug/math_sprint_m2_d4_substitution.py`
**Data:** `debug/data/math_sprint_m2_d4.json`
**Audit triggers:** `debug/audit_critical_issues_review.md` §M2; `debug/review_paper18.md` H2 / M1.

---

## 1. Error specification

Paper 18 §IV Theorem 1 part (2) currently presents (at L1196–1206):

```
Eq.(38) of audit / paper Eq.(\ref{eq:dirac_dirichlet_zeta3}):

  D_{g^Dirac}(4)
    = sum_{m >= 1} 2 m (m+1) / m^4
    = 2 zeta_R(2) + 2 zeta_R(3)
    = pi^2/3 + 2 zeta_R(3)   ≈ 5.694                          ❌ WRONG
```

This conflicts with:

- **Paper 18 Eq.(30)** itself (the *correct* Hurwitz decomposition derived two paragraphs earlier in Theorem 1 part (1)):
  $$D(s) = 2\,(2^{s-2} - 1)\,\zeta_R(s-2) \;-\; \tfrac{2^s - 1}{2}\,\zeta_R(s).$$
  At $s = 4$ this gives $D(4) = 6\zeta_R(2) - \tfrac{15}{2}\zeta_R(4) = \pi^2 - \pi^4/12$.

- **Paper 28 Eq.(43)** and Table 2: $D(4) = \pi^2 - \pi^4/12 \approx 1.752$.

- The Camporesi–Higuchi spectrum on unit $S^3$: $|\lambda_n| = n + 3/2$, $g_n = 2(n+1)(n+2)$. Paper 18's Eq.(38) accidentally **substituted integer $m$ for the half-integer eigenvalues** $n + 3/2$, dropping the $3/2$ shift.

---

## 2. Confirmed correct $D(4)$ value at high precision

All three methods (direct truncated CH sum, Eq.(30) Hurwitz form, closed form $\pi^2 - \pi^4/12$) **agree bit-exact to 60 decimal places**:

```
D(4) = 1.7521801482558221824644632758173918645030672678503...
```

residuals: $|$direct $-$ closed$| \sim 1.6 \times 10^{-61}$, $|$direct $-$ Eq.(30)$| \sim 1.6 \times 10^{-61}$.

Paper 18's wrongly-claimed value $2\zeta(2) + 2\zeta(3) \approx 5.694$ differs from $D(4)$ by $\approx 3.94$ — gross, not a rounding glitch.

**The number $5.694$ is itself correct** as the value of the wrongly-written sum $\sum_{m \ge 1} 2m(m+1)/m^4$ — what is wrong is its identification with $D(4)$. The corrected substitution is given in §4 below.

---

## 3. Cross-validation against Paper 28 Table 2 (even $s$)

| $s$ | Direct CH sum (60 dps) | Closed form | Match |
|---|---|---|---|
| 2 | DIVERGES (linear growth, $\sim 19999$ at $N=10^4$) | Analytic continuation: $-\tfrac{3}{2}\zeta_R(2) = -\pi^2/4$ via Eq.(30) | n/a (CH sum classically divergent; Eq.(30) gives reg. value) |
| 4 | $1.7521801482558221824644632758\ldots$ | $\pi^2 - \pi^4/12$ | ✓ to 60 dps |
| 6 | $0.4233905588239978444727961078\ldots$ | $\pi^4/3 - \pi^6/30$ | ✓ to 60 dps |
| 8 | $0.1653628948026883332468912438\ldots$ | $2\pi^6/15 - 17\pi^8/1260$ | ✓ to 60 dps |

(s=2 is a separate analytic-continuation subtlety; Paper 28 Theorem 1 deals with this via $\zeta(0) = -1/2$.)

This **confirms Paper 28's Theorem 5 / Table 2 line-by-line** at >60 dps and demonstrates **the *Hurwitz decomposition* (Paper 18 Eq.(30)) is internally consistent**; only the explicit substitution at Eq.(38) is wrong.

---

## 4. Where $\zeta_R(3)$ actually appears on the CH spinor bundle

This is the substantive structural finding.

Applying Eq.(30) at each integer $s$:

| $s$ | $D(s)$ closed form via Eq.(30) | $\zeta_R(3)$ present? |
|---|---|---|
| 2 | $-\tfrac{3}{2}\zeta_R(2) = -\pi^2/4$ | No |
| 3 | $2 \cdot (2^1 - 1) \cdot \zeta_R(1) - \tfrac{7}{2}\zeta_R(3)$ | **$\zeta_R(1)$ DIVERGES** ⇒ structurally inadmissible |
| 4 | $6\zeta_R(2) - \tfrac{15}{2}\zeta_R(4) = \pi^2 - \pi^4/12$ | No (purely $\pi^{\text{even}}$) |
| 5 | $14\,\zeta_R(3) - \tfrac{31}{2}\,\zeta_R(5)$ | **YES — first clean appearance** |
| 6 | $30\zeta_R(4) - \tfrac{63}{2}\zeta_R(6) = \pi^4/3 - \pi^6/30$ | No |
| 7 | $62\,\zeta_R(5) - \tfrac{127}{2}\,\zeta_R(7)$ | No ($\zeta_R(5), \zeta_R(7)$ instead) |
| 8 | $126\zeta_R(6) - \tfrac{255}{2}\zeta_R(8) = 2\pi^6/15 - 17\pi^8/1260$ | No |
| 9 | $254\,\zeta_R(7) - \tfrac{511}{2}\,\zeta_R(9)$ | No |

**Verified at 50+ dps:**
$$D(5) = 14\,\zeta_R(3) - \tfrac{31}{2}\,\zeta_R(5) = 0.75641643951208613746016922107627026582505384048734\ldots$$
direct CH sum vs closed-form residual: $3.9 \times 10^{-61}$.

**Structural conclusion.** The CH first-order spectral zeta $D(s)$ **does** have a clean, bit-exact $\zeta_R(3)$ witness — but it lives at $s = 5$, not at $s = 4$. At $s = 4$ the half-integer eigenvalues conspire to produce **purely $\pi^{\text{even}}$ content**, exactly because the squared spectrum at $s = 4$ is the squared-Dirac spectral zeta at $s = 2$, which by Paper 28 T9 must be in $\mathbb{Q}[\pi^2]$.

Paper 18 §IV.B Table (Eq.~\ref{eq:hurwitz_discriminant}) was already in the parity discriminant: "At even~$s$ both terms are $\pi^{\mathrm{even}}$; at odd $s \ge 3$ the terms are $\mathbb{Q}$-linearly independent of $\pi$" — this discriminant is correct. The only edit needed is **at $s = 4$ specifically** (Eq.~\ref{eq:dirac_dirichlet_zeta3}), which should not claim $\zeta_R(3)$ content.

---

## 5. Diagnosis: what Eq.(38) should look like with correct half-integer eigenvalues

The correct half-integer derivation:

$$D(4) \;=\; \sum_{n \ge 0} \frac{2(n+1)(n+2)}{(n + 3/2)^4} \;=\; 16 \sum_{n \ge 0} \frac{2(n+1)(n+2)}{(2n+3)^4}.$$

Substitute $m = 2n+3$ so $m$ runs over odd integers $\ge 3$, with $(n+1)(n+2) = (m^2-1)/4 \cdot (1/2)\cdot ... $ — cleanly: $2(n+1)(n+2) = (m-1)(m+1) \cdot \tfrac{1}{2} \cdot 2 = (m^2-1)/2$. Then

$$D(4) \;=\; 16 \cdot \tfrac{1}{2} \sum_{m \ge 3,\, m\text{ odd}} \frac{m^2 - 1}{m^4} \;=\; 8 \biggl[ \sum_{m \ge 1,\, m\text{ odd}} \frac{1}{m^2} - \sum_{m \ge 1,\, m\text{ odd}} \frac{1}{m^4} \biggr]$$

(the $m=1$ term contributes $0/2 = 0$ thanks to the half-integer shift — the *singular* term automatically cancels: this is the structural difference from Paper 18's wrong $m$-sum-starting-at-1 form). Then with $\lambda(s) := \sum_{m\,\text{odd}\ge 1} 1/m^s = (1 - 2^{-s})\zeta_R(s)$:

$$D(4) \;=\; 8\,\bigl[\lambda(2) - \lambda(4)\bigr] \;=\; 8\,\bigl[\tfrac{3}{4}\zeta_R(2) - \tfrac{15}{16}\zeta_R(4)\bigr] \;=\; 6\zeta_R(2) - \tfrac{15}{2}\zeta_R(4) \;=\; \pi^2 - \tfrac{\pi^4}{12}.$$

Numerical confirmation at 40 dps: closed-form, $\lambda$-rewriting, and direct CH sum all agree at relative residual $\sim 10^{-61}$ (`task5_corrected_eq38_diagnostic` in JSON).

---

## 6. Proposed paper-text patches

### Patch (a): replace Paper 18 Eq.(38) entirely

**Location:** `papers/group3_foundations/paper_18_exchange_constants.tex`, L1195–1221 (paragraph "Odd-zeta content from first-order spectral operators" through to "$\zeta_R(3)$ is structurally present on $S^3$ without participating in $K$.").

**Replacement.** Remove the entire `\begin{equation}\label{eq:dirac_dirichlet_zeta3}` block and the two paragraphs that interpret it as a $\zeta(3)$ source at $s = 4$. Replace with a corrected statement that:

1. Gives $D(4) = \pi^2 - \pi^4/12$ (purely $\pi^{\text{even}}$, matching Paper 28 Eq.(43)).
2. Re-identifies the first $\zeta_R(3)$ appearance on the CH spinor bundle at $s = 5$:
   $$D(5) \;=\; 14\,\zeta_R(3) - \tfrac{31}{2}\,\zeta_R(5).$$
3. Notes that the "odd-zeta content from first-order spectral operators" claim survives — it just first appears at odd $s$, not at the packing exponent $s = 4$.

Exact LaTeX replacement (drop-in for L1199–1221):

```latex
The Dirac operator on unit $S^3$ (Camporesi--Higuchi~\cite{camporesi_higuchi1996}
spectrum $|\lambda_n| = n + 3/2$, degeneracy
$g_n^{\text{Dirac}} = 2(n+1)(n+2)$) introduces a structurally new
transcendental: odd-zeta values $\zeta_R(2k{+}1)$ enter the Dirac
Dirichlet series $D(s) = \sum_n g_n^{\text{Dirac}}/|\lambda_n|^s$ at
odd integer $s \ge 5$.  The Hurwitz decomposition
(Eq.~\eqref{eq:hurwitz_discriminant} below) gives
\begin{equation}\label{eq:dirac_zeta3_witness}
  D(5) \;=\; 14\,\zeta_R(3) \;-\; \tfrac{31}{2}\,\zeta_R(5),
\end{equation}
the first integer~$s$ at which $\zeta_R(3)$ appears with a non-divergent
coefficient (at $s = 3$ the term proportional to $\zeta_R(s{-}2) = \zeta_R(1)$
diverges).  At the packing exponent $s = 4$, by contrast,
\begin{equation}\label{eq:dirac_d4_value}
  D(4) \;=\; 6\,\zeta_R(2) - \tfrac{15}{2}\,\zeta_R(4)
        \;=\; \pi^2 - \frac{\pi^4}{12},
\end{equation}
purely $\pi^{\text{even}}$---consistent with Theorem~\ref{thm:T9} below
(squaring the spectrum maps every odd-power Dirichlet sum to an even
power).  The half-integer eigenvalue $|\lambda_n| = n + 3/2$ is essential:
naive substitution of an integer~$m$ for the eigenvalue would predict an
even-$s$ odd-zeta witness, contradicting Theorem~1 part~(2) and
Paper~28~\cite{loutey_paper28} Theorem~5 / Eq.~(43).
The shift~$3/2$ is what cancels the singular $m=1$ term in
the odd-integer rewrite of the spectral sum.

This odd-zeta content (entering at odd $s$ on the spinor bundle) is
categorically different from calibration-tier content.  Calibration
constants carry even-zeta or $\pi^{\text{even}}$ content arising from
second-order Riemannian operators through Jacobi-$\theta$ inversion;
the odd-zeta $\zeta_R(3)$ enters through the first-order Dirac operator
via half-integer Hurwitz-zeta identities at \emph{odd} $s$, and is
$\mathbb{Q}$-linearly independent of $\zeta_R(2)$ by Ap\'ery's
theorem.  Within Paper~2's $K = \pi(B + F - \Delta)$, the odd-zeta tier
is a byproduct rather than a summand: $\zeta_R(3)$ is structurally
present on $S^3$ without participating in $K$.
```

This **(i)** removes the wrong Eq.(38) sum, **(ii)** replaces it with the correct $s = 5$ $\zeta_R(3)$ witness (Eq.~\ref{eq:dirac_zeta3_witness}), **(iii)** adds the corrected $s = 4$ value Eq.~\ref{eq:dirac_d4_value} explicitly, and **(iv)** flags the half-integer-vs-integer substitution as the load-bearing structural point (turning the previous error into a pedagogical aside).

### Patch (b): §IV motivic-weight loop-order paragraph (L1446–1473)

The current text says: "at $s = 4$, $D(4) = 2\zeta_R(2) + 2\zeta_R(3)$, reproducing Eq.~(\ref{eq:dirac_dirichlet_zeta3})." This sentence must be deleted and replaced with the $s = 5$ witness.

**Drop-in replacement at L1457–1464:**

```latex
The even/odd Hurwitz spectral zeta of the first-order Dirac operator
(Eq.~\ref{eq:hurwitz_discriminant} of Theorem~\ref{thm:operator_order})
provides a sharp computational check.
This exact formula separates the even and odd content at each~$s$:
at $s = 4$, both terms are $\pi^{\text{even}}$
(Eq.~\ref{eq:dirac_d4_value}: $D(4) = \pi^2 - \pi^4/12$); the first
odd-zeta witness on the CH spinor bundle is
Eq.~\eqref{eq:dirac_zeta3_witness}: $D(5) = 14\,\zeta_R(3) -
\tfrac{31}{2}\,\zeta_R(5)$.  At even~$s$ the spectral sum is purely
$\pi^{\text{even}}$; at odd $s \ge 5$ the $\zeta_R(s{-}2)$ and
$\zeta_R(s)$ terms are independently transcendental.
```

### Patch (c): §IV "motivic correspondence ⇒ QED loop order" argument

This argument **survives the correction**, but its sharpness changes. The current text reads:

> "The T9 theorem (Eq.~\ref{eq:zeta_d2}) guarantees that $\zeta_R(3)$ cannot enter through any single trace of a squared propagator $D^{-2s}$---i.e., at one loop.  It enters at two loops, where irreducible products of first-order propagators $D^{-1} \otimes D^{-1}$ produce odd motivic weight."

The structural claim "$\zeta_R(3)$ enters first-order Dirac sectors but not squared-Dirac sectors" is **correct and untouched** by M2 — it follows directly from T9 (the $D^2$ spectral zeta is always $\pi^{\text{even}}$). What the M2 correction sharpens is the GeoVac-internal **integer-$s$ location** of the witness: not the packing exponent $s = 4$, but the next odd integer $s = 5$.

The motivic-weight reading is preserved because $D(5)$ contains $\zeta_R(3)$ at weight 3 — Eq.~\eqref{eq:dirac_zeta3_witness} is the GeoVac realization of Zagier's $d_3 = 1$ statement on the CH spinor bundle. The two-loop Petermann–Sommerfield $\zeta(3)$ is matched at the **motivic-weight level** (both are weight-3 transcendentals from first-order propagators), but the *specific* sum identity that was claimed at the packing exponent has to migrate to $s = 5$.

**Recommended additional sentence after the loop-order argument** (insert after L1473):

```latex
The natural CH-internal realization of the weight-3 motive
$\mathfrak{M}(3) = \mathbb{Q}(-3)$ is Eq.~\eqref{eq:dirac_zeta3_witness}
at $s = 5$ (not at the packing exponent $s = 4$, which is forced into
$\mathbb{Q}[\pi^2]$ by Theorem~\ref{thm:T9}).  The Petermann--Sommerfield
two-loop $g{-}2$ coefficient~\cite{petermann1957,sommerfield1957} is
matched at the motivic-weight level: both are weight-3 periods entering
through irreducible products of first-order Dirac propagators on $S^3$.
```

This **preserves the §IV argument's content** (motivic correspondence ⇒ loop-order discriminant ⇒ Petermann–Sommerfield witness) while moving the *literal GeoVac equation* from $s = 4$ to $s = 5$.

### Patch (d): paragraph at L1496–1503 mentioning "$D_{\text{even}} + D_{\text{odd}} = D(4) = \pi^2 - \pi^4/12$"

This paragraph is **already correct** — it uses the right $D(4)$ value. No change needed.

### Patch (e): §IX summary

The §IX summary repeats the $\zeta_R(3)$ on the spinor bundle claim at L2707 and L2799. These references are at the *motivic-tier* level, not at a specific integer $s$, and should be updated to reference Eq.~\eqref{eq:dirac_zeta3_witness} instead of the (deleted) Eq.~\eqref{eq:dirac_dirichlet_zeta3}.

---

## 7. Cross-validation against Paper 28 verified identities

| Identity | Paper 28 location | M2 driver result | Verdict |
|---|---|---|---|
| $D(4) = \pi^2 - \pi^4/12$ | Eq.(43), Theorem 5 / Table 2 | bit-exact to 60 dps via three independent methods | ✓ Paper 28 is the correct source |
| $D(5) = 14\zeta(3) - 31\zeta(5)/2$ | Table 2 | bit-exact to 60 dps direct vs closed | ✓ |
| $D(6) = \pi^4/3 - \pi^6/30$ | Table 2 | bit-exact to 60 dps | ✓ |
| $D(8) = 2\pi^6/15 - 17\pi^8/1260$ | Table 2 | bit-exact to 60 dps | ✓ |
| $D_{\text{even}}(4) - D_{\text{odd}}(4) = 8(\beta(4) - G)$ | Eq.~\eqref{eq:D_diff_closed} | not in this sprint (Paper 28 §IV already verifies this) | ✓ (audit-verified, see review_paper18.md row 4) |
| T9: $\zeta_{D^2}(s) \in \mathbb{Q}[\pi^2]$ | Eq.~\eqref{eq:T9} | structural — not directly tested here | ✓ (audit row 4 already verified) |

**Bottom line.** Paper 28's verified-bit-exact identities ARE the canonical source. The M2 fix is to migrate Paper 18's Eq.(38) from a wrongly-substituted "integer-$m$" form (which produced the right number $\approx 5.694$ but identified it with the wrong observable) to the correct half-integer CH form (which gives $D(4) = \pi^2 - \pi^4/12$). The $\zeta_R(3)$ "witness" on the CH spinor bundle then sits at $s = 5$ (clean, Apéry-verified, $\mathbb{Q}$-linearly independent of $\pi$).

---

## 8. Decision gate verdict

> **DECISION GATE (per sprint brief):** M2 is mostly computational. If the corrected Eq.(38) and a clean ζ(3) source emerge, the §IV argument can be salvaged with patches. If the spinor-bundle ζ(3) source turns out structurally absent, the §IV argument may need a content rewrite.

**Verdict: SALVAGE WITH PATCHES.** A clean ζ(3) source exists on the CH spinor bundle — Eq.~\eqref{eq:dirac_zeta3_witness}: $D(5) = 14\zeta_R(3) - \tfrac{31}{2}\zeta_R(5)$, verified bit-exact at 50 dps. The §IV motivic-weight loop-order argument survives the correction: the discriminant (1st-order spinor ⇒ odd-zeta vs 2nd-order ⇒ $\pi^{\text{even}}$) is unchanged. Only the *specific integer $s$* where the GeoVac internal witness sits migrates from $s = 4$ (purely $\pi^{\text{even}}$) to $s = 5$.

Recommended patches (a)+(b)+(c)+(e) above are mechanical drop-in replacements. After application, Paper 18 §IV becomes internally consistent with Paper 28 Table 2 / Eq.(43), removes the Eq.(38) substitution error, and preserves the load-bearing motivic-weight–loop-order correspondence.

---

## 9. Files

- driver: `debug/math_sprint_m2_d4_substitution.py` (~360 lines)
- data: `debug/data/math_sprint_m2_d4.json`
- this memo: `debug/math_sprint_m2_memo.md`
- target paper file: `papers/group3_foundations/paper_18_exchange_constants.tex` (patches not applied yet — proposed text above; PI sign-off recommended before in-place edits per CLAUDE.md §13.8 since Paper 18 is an Always-load paper and the motivic-weight argument is load-bearing for the bound-state QED arc / Paper 36 cross-reference).
