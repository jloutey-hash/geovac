# Math sprint M3 — Paper 28 Table 1 resolution

**Date:** 2026-06-02
**Driver:** `debug/math_sprint_m3_table1.py`
**Data:** `debug/data/math_sprint_m3_table1.json`
**Verdict:** RESOLVED. Paper 28 Table 1 is wrong at all four rows; the
correct values follow directly from Theorem 1 / Eq. (5) and have been
verified three independent ways at 100 dps with all cross-checks closing at
≥ 60 dps.

---

## 1. The correct closed forms

Convention (Paper 28 §III, Eq. (5)):

$$
\zeta_{D^2}(s) \;=\; 2^{2s-1}\bigl[\lambda(2s-2) - \lambda(2s)\bigr],
\qquad
\lambda(2k) = (1 - 2^{-2k})\,\zeta_R(2k),
$$

with the mode-sum definition

$$
\zeta_{D^2}(s) \;=\; \sum_{n\ge 0} \frac{g_n}{|\lambda_n|^{2s}},
\qquad g_n = 2(n+1)(n+2),
\qquad |\lambda_n| = n + \tfrac{3}{2},
$$

equivalently $\zeta_{D^2}(s) = 2\,\zeta_H(2s-2, \tfrac{3}{2}) - \tfrac{1}{2}\zeta_H(2s, \tfrac{3}{2})$ (Hurwitz-zeta form, also already in Paper 28 Eq. just above (6)).

At $s = 1, 2, 3, 4$:

| $s$ | $\zeta_{D^2}(s)$ closed form | Compact form | Numerical value (50 dps) |
|:---:|:-----------------------------|:-------------|:-------------------------|
| 1 | $-\dfrac{\pi^{2}}{4}$ | $-\dfrac{\pi^{2}}{4}$ | $-2.4674011002723396547086227499690377838284248518102$ |
| 2 | $\pi^{2} - \dfrac{\pi^{4}}{12}$ | $\pi^{2} - \dfrac{\pi^{4}}{12}$ | $\phantom{-}1.7521801482558221824644632758173918645030672678503$ |
| 3 | $\dfrac{\pi^{4}(10 - \pi^{2})}{30}$ | $\dfrac{\pi^{4}}{3} - \dfrac{\pi^{6}}{30}$ | $\phantom{-}0.4233905588239978444727961078210404543352776215269$ |
| 4 | $\dfrac{\pi^{6}(168 - 17\pi^{2})}{1260}$ | $\dfrac{2\pi^{6}}{15} - \dfrac{17\pi^{8}}{1260}$ | $\phantom{-}0.1653628948026883332468912438027971951115029362837$ |

### Cross-checks

| Cross-check | Worst residual across $s\in\{1,...,4\}$ |
|:------------|:----------------------------------------|
| direct mode-sum (8000 terms + Hurwitz tail) vs Hurwitz-zeta form | $\le 4.6\times 10^{-97}$ |
| Hurwitz-zeta form vs Theorem 1 / Eq. (5) | $\le 1.9\times 10^{-99}$ |
| Theorem 1 / Eq. (5) vs sympy symbolic (substituting $\zeta_R(2k)$ Euler–Bernoulli) | $\le 4.5\times 10^{-60}$ (sympy `N(..., 60)`) |

All three independent computations agree at 100 dps. Theorem 1 itself is
verified as the closed form (the existing Paper 28 verification note
"sympy for $s=1,...,10$; agreement to $10^{-60}$ at $s=2,...,20$" stands —
only the table values were wrong).

### Sign convention

This convention is the **standard heat-kernel/spectral-zeta** convention
$\zeta_{D^2}(s) = \mathrm{Tr}\,|D^2|^{-s} = \sum_n g_n/\lambda_n^{2s}$
with $g_n > 0$ throughout, and **no** eta-invariant or Tate-twist sign
factor is needed. The sign of $\zeta_{D^2}(1)$ is genuinely negative
because $\lambda(0) = (1-1)\zeta_R(0) = 0$ in the closed form, so the
$s=1$ value is

$$
\zeta_{D^2}(1) = 2\bigl[\lambda(0) - \lambda(2)\bigr]
              = -2\lambda(2)
              = -2 \cdot \tfrac{3}{4} \cdot \tfrac{\pi^2}{6}
              = -\tfrac{\pi^2}{4}.
$$

This matches the audit's "should be $-\pi^2/4$". The negative sign is a
consequence of the $\lambda$-difference structure; no other convention
choice is needed.

---

## 2. Proposed exact replacement for Table 1 in Paper 28

Replace the four body rows of `tab:D2_values` in
`papers/group5_qed_gauge/paper_28_qed_s3.tex` lines 295–298 with:

```latex
1 & $-\pi^2/4$ \\
2 & $\pi^2 - \pi^4/12$ \\
3 & $\pi^4/3 - \pi^6/30$ \\
4 & $2\pi^6/15 - 17\pi^8/1260$ \\
```

(Caption and headers unchanged.)

A factored alternative for $s=3, 4$ if a more compact rendering is
preferred:

```latex
3 & $\pi^4(10 - \pi^2)/30$ \\
4 & $\pi^6(168 - 17\pi^2)/1260$ \\
```

Either form is correct; the un-factored form matches Table 2's house style
more closely.

---

## 3. Root-cause diagnosis

Comparing each original row to the (newly-derived) correct values:

| Original row $s$ | Original value | Matches $\zeta_{D^2}(s')$ for $s' = ?$ |
|:----------------:|:---------------|:--------------------------------------|
| 1 | $\pi^{2} - \pi^{4}/12$ | $s' = 2$ (bit-exact) |
| 2 | $2\pi^{4}/3 - 2\pi^{6}/15$ | none in $\{1,\ldots,6\}$ |
| 3 | $16\pi^{6}/15 - 16\pi^{8}/63$ | none in $\{1,\ldots,6\}$ |
| 4 | $32\pi^{8}/9 - 256\pi^{10}/495$ | none in $\{1,\ldots,6\}$ |

The $s=1$ row is a **clean off-by-one in $s$**: the original entry is
**exactly** the correct $\zeta_{D^2}(2)$ value, which is also Paper 28
Table 2's published $D(4)$ value (verified: Paper 28 line 386:
`4 & $\pi^2 - \pi^4/12$`). This confirms the audit's observation —
*"`$\pi^2 - \pi^4/12$` should be `$D(2 \cdot 2) = D(4) = \zeta_{D^2}(2)$`,
not `$\zeta_{D^2}(1)$`"*.

The $s=2, 3, 4$ rows do **not** match $\zeta_{D^2}(s')$ for any
$s' \in \{1,\ldots,6\}$ — they are not a clean off-by-one (or off-by-two,
or off-by-any-shift) of the genuine values. Examining the structure:

- $\zeta_{D^2}(3) = 2\pi^4/3 - 2\pi^6/15$ has the same leading rational
  prefactor $2/3$ as the original $s=2$ row, but the wrong $\pi$-exponents.
- The progression $\pi^4 \to \pi^6 \to \pi^8 \to \pi^{10}$ in the second
  column of the original rows matches the *correct* $\pi$-exponents for
  $\zeta_{D^2}(2), \zeta_{D^2}(3), \zeta_{D^2}(4), \zeta_{D^2}(5)$
  respectively (shifted by one), but the *first* column $\pi$-exponents
  are also shifted, and the rational prefactors are independently wrong.

**Most plausible root cause**: a hand-copy of an earlier draft of Table 1
that used a different normalization convention or a different prefactor
power, combined with the $s$-column being mis-aligned. The cleanest
hypothesis is that an old draft tabulated $2\zeta_{D^2}(s+1)$ (or
equivalently $4^s[\lambda(2s) - \lambda(2s+2)]$) at $s = 1, 2, 3, 4$ and
the labels did not get updated when the convention was re-fixed to the
current Theorem 1 form. Whatever the exact transcription history, the
practical takeaway is that the $s=1$ row provably equals the correct
$s=2$ value, and the remaining three rows are inconsistent with any
single transformation of the genuine values — they should simply be
replaced.

**No theorem, proof, equation, or other table in Paper 28 is affected.**
The Hurwitz-zeta closed form (line 191–192), Theorem 1 / Eq. (5), the
Theorem 2 proof (which uses Eq. (5) in the same line of argument), and
all of Tables 2–7 are independent of the wrong Table 1 entries. The
verification note on line 283–285 ("symbolic identity confirmed in sympy
for $s = 1, \ldots, 10$. Numerical agreement to $10^{-60}$ at
$s = 2, \ldots, 20$") refers to Eq. (5), not the table, and remains
correct.

---

## 4. Apply-the-patch checklist

For the PI applying the patch:

1. Open `papers/group5_qed_gauge/paper_28_qed_s3.tex`.
2. Replace lines 295–298 with the four rows in §2 above.
3. No other edits to Paper 28 are required (Theorem 1, the proof, Table 2,
   Remark `rem:zetaD2_no_fe`, and all downstream sections are unaffected).
4. Verify by recompiling and spot-checking row $s=2$ against Table 2 row
   $s=4$: they must match (both equal $\pi^2 - \pi^4/12 \approx 1.7522$).

---

## 5. Files produced

- `debug/math_sprint_m3_table1.py` (driver)
- `debug/data/math_sprint_m3_table1.json` (corrected values + diagnostic)
- `debug/math_sprint_m3_memo.md` (this memo)
