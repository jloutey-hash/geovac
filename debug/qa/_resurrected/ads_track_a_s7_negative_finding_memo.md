# Track AdS-A S^7 Sprint — Negative Finding on Ring Conjecture

**Date:** 2026-05-25 (continuation post-v3.2.1)
**Sprint:** Testing Paper 50 Remark 7.5 conjecture on round S^7.
**Verdict:** Conjecture FALSIFIED in its simple form. S^7 scalar zeta'(0) is NOT in {log 2, ζ(3)/π², ζ(5)/π⁴, ζ(7)/π⁶} ring at PSLQ tolerance 10^-30 to 10^-80 with maxcoeff up to 10^16.

---

## §0. The framework's value

Numerical computation via Hurwitz expansion at $\mathrm{dps} = 200$, $k_{\max} \in \{40, 80, 100\}$ all converge to the same value:

$$\zeta'_{S^7,\text{conf}}(0) = -0.00159466155346957386175492433977\ldots$$

Convergence is fast: identical to displayed 30 digits across $k_{\max} \ge 40$, indicating the Hurwitz series is rapidly convergent (geometric ratio $1/4$) and the value is correct to many more digits than displayed.

So $F^{S^7}_{\text{scalar}} = -\tfrac{1}{2}\zeta'(0) \approx +0.000797$.

---

## §1. PSLQ search across multiple bases

Tried bases and PSLQ tolerance/maxcoeff combinations:

| Basis | PSLQ result |
|:------|:-----------:|
| $\{\log 2, \zeta(3)/\pi^2, \zeta(5)/\pi^4, \zeta(7)/\pi^6\}$ | no relation found |
| above + $\{1\}$ | no relation found |
| $\{\log 2, \zeta(3), \zeta(5), \zeta(7), \pi^2, \pi^4, \pi^6\}$ | no relation found |
| above + $\{\zeta(9)/\pi^8\}$ | no relation found |
| above + $\{\log(\pi)\}$ | no relation found |
| above + $\{\log 3\}$ | no relation found |

All trials at tolerance $10^{-30}, 10^{-50}, 10^{-80}$ and maxcoeff $10^8, 10^{12}, 10^{16}$ returned no integer relation.

---

## §2. Possible interpretations

**(a) Coefficients are larger than 10^16.** Unlikely for a closed-form CFT partition function value; standard KPS-type formulas have small-integer coefficients (S^3 was $[8, 2, -3]$, S^5 was $[32, -2, -2, 15]$).

**(b) Additional transcendentals are needed.** The full ring might include $\log 3$, $\log 5$, Catalan $G$, Dirichlet $\beta(2k)$ values, or other less-obvious constants.

**(c) The conjectured Remark 7.5 ring is wrong.** The (d-1)/2+1 dimensional ring with only $\log 2$ and odd-zeta-over-even-pi might not extend to higher odd $d$. The structural cancellations that worked on S^3 (no log p for any odd prime p) and S^5 (log 3 cancelled via the (1, -5/2, 9/16) combination) might fail at S^7 where the multiplicity polynomial (1, -5, 4 with /360) is structurally different.

**(d) Computational error.** Possible but unlikely given the Hurwitz expansion's robustness on S^3 and S^5. The verification of $\zeta_R'(-2), \zeta_R'(-4), \zeta_R'(-6)$ against the standard closed forms passed to 200+ digits.

---

## §3. What this means for Paper 50

The Paper 50 Remark 7.5 conjecture should be **honestly weakened** to:

**Refined statement:** On round $S^3$ and $S^5$, the framework's spectral-zeta partition function values are bit-exactly in the (d-1)/2+1-dim ring $\{\log 2, \zeta(3)/\pi^2, \ldots, \zeta(d-2)/\pi^{d-3}\}$. The extension to $S^d$ for $d \ge 7$ is **open**; preliminary computation at $S^7$ does not match this simple ring at PSLQ tolerance 10^-30 with maxcoeff $10^{16}$.

This is itself a paper-worthy structural observation — the simple pattern breaks (or at least, its closed form becomes complex enough to require an enlarged basis or different structural treatment) at $S^7$.

---

## §4. Open follow-on directions

1. **Compute the S^7 value via an INDEPENDENT method** (e.g., direct numerical Σ + Euler-Maclaurin asymptotic correction, or KPS-style Barnes-zeta closed form). If the value matches, then the issue is genuinely PSLQ basis incompleteness.

2. **Locate the published closed form for S^7 free CFT partition function.** Beccaria-Tseytlin and related literature should have this. Compare.

3. **Look for product terms** like $\log 2 \cdot \zeta(3)/\pi^2$ that linear PSLQ wouldn't find.

4. **Test a Dirichlet L-value basis** (Catalan $G$, $\beta(2)$, $\beta(4)$, $L(s, \chi_{-4})$ at various $s$). These appear naturally from vertex-parity sums in QED (Paper 28); they might appear in S^7 spectral data via Hurwitz at higher shifts.

5. **Honest write-up:** add a Remark 8.1 (or footnote) to Paper 50 §7 documenting the negative finding at $S^7$ and stating the open question.

---

## §5. Files

- `debug/ads_track_a_s7_scalar_partition_function.py` — main script
- `debug/ads_s7_pslq_search.py` — wider PSLQ basis search
- `debug/data/ads_track_a_s7_scalar_partition_function.json` — result (no PSLQ relation)
- THIS memo
