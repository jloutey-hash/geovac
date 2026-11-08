# Sprint: Forcing-catalogue forward-run (Directions 1/2/3) + verify-the-verifier

**Date:** 2026-06-03
**Type:** Leader-brief forward-run — 3 sub-agent dispatches + main-session re-verify + 3 paper corrections
**Verdict:** Door 1 ladder GRADUATES (+ corrected a released-paper false negative); seam → relative theorem of necessity (placed); Paper 54 two-body → STOP (paper vindicated, one equation fixed).

Drivers/data: `debug/door1_s11_graduation.py`+`.json`; `debug/seam_packing_scoping_memo.md`;
`debug/paper54_recompute_both_constructions.py`+`debug/data/paper54_recompute_both.json`;
`debug/paper54_two_body_forward_scoping_memo.md` (with correction banner).

---

## 1. The arc in one paragraph

A Leader Strategic Brief picked three forward-run Directions off the Forcing
Catalogue; three sub-agents executed them in parallel. A main-session
verify-the-verifier pass then **overturned a sub-agent claim in each of the two
paper-touching Directions** — once finding a real false-negative the agent had
accidentally surfaced (Paper 50), once catching a phantom "erratum" that would
have corrupted a correct paper (Paper 54). Net: three released/active papers
corrected or extended, all compiling clean; one new structural result placed.

## 2. Direction 1 — Door 1 F-coefficient ladder graduation

**Verdict: GO on the ladder, STOP on the recursion constant.**

- **(A) GO.** The conformal-scalar and Weyl-Dirac F-coefficients
  $F = -\tfrac12\zeta'_\Delta(0)$ on round $S^d$ ($d$ odd) lie bit-exactly in the
  ring $\mathcal{R}_d = \{\log 2\}\cup\{\zeta(2k{+}1)/\pi^{2k}: k=1..(d{-}1)/2\}$
  (top atom $\zeta(d)/\pi^{d-1}$) — verified at $S^7, S^9, S^{11}$, 260 dps,
  frozen basis, decoy-controlled (Catalan $G$, $\log 3$, unscaled $\zeta(3)$ all
  $\to 0$). No out-of-ring transcendental → WH2 grid intact, no 4th Mellin
  sub-mechanism. S¹¹ scalar closed form (sample):
  $-\tfrac{7}{131072}\log2 - \tfrac{3897}{45875200}\tfrac{\zeta3}{\pi^2}
   - \tfrac{485}{8257536}\tfrac{\zeta5}{\pi^4} + \tfrac{609}{1310720}\tfrac{\zeta7}{\pi^6}
   + \tfrac{425}{262144}\tfrac{\zeta9}{\pi^8} + \tfrac{1023}{524288}\tfrac{\zeta11}{\pi^{10}}$.

- **(B) STOP.** The cross-species recursion
  $\mathrm{scalar} - m\,\mathrm{Dirac} = c_d\,\mathrm{Dirac}_{S^{d-2}}$ closes at
  $S^5$ only (top-atom ratio $2^{-(d-5)/2}$ hits the needed $m{=}1$ at $d{=}5$
  alone). The borderline $\{c_5{=}1, c_7{=}c_9{=}1/2\}$ was an artifact of a
  non-frozen Dirac normalization (a hidden binary free parameter). Curve-fit
  audit caught it.

**Verify-the-verifier — the false negative.** The S¹¹ result CONTRADICTED Paper
50's existing remark, which documented the $S^7$ scalar value
$\zeta'(0)\approx-0.00159466$ as having **no PSLQ relation** in
$\{\log2,\zeta3/\pi^2,\zeta5/\pi^4,\zeta7/\pi^6\}$. Independent main-session check
(driver's value fn + my own `mpmath.pslq`, 200 dps): the value is correct
(three computations agree) and the relation IS genuine —
$-\tfrac1{512}\log2 - \tfrac{41}{15360}\tfrac{\zeta3}{\pi^2} + \tfrac5{1024}\tfrac{\zeta5}{\pi^4}
+ \tfrac{63}{2048}\tfrac{\zeta7}{\pi^6}$, reconstruction error $2.4\times10^{-202}$.
The old "negative" was a **false negative** from an under-resolved search
(30-digit value vs maxcoeff $10^{16}$; the real relation has coeffs $\sim10^4$ and
needs $\gtrsim200$ dps to lock). My initial JSON inspection truncated the
reconstruction-error exponents with `[:12]` (showing "4.37" for "4.37e-262"),
which briefly looked like the relation was spurious — corrected by the
from-scratch PSLQ.

**Paper edit (Paper 50 §7).** Rewrote `rem:general_d_conjecture`: ladder continues
through $S^{11}$ (closed forms + erratum correcting the false $S^7$ negative);
fixed a ring-spec typo ($\zeta(d-2)/\pi^{d-3}\to\zeta(d)/\pi^{d-1}$); preserved the
still-correct dual-basis non-extension (Prop `dual_basis_S5`) — the actual
structural-specificity content, untouched. 17 pp, two-pass clean, 0 warnings.

## 3. Direction 2 — forced/free seam scoping

**Verdict: NO-GO (and the NO-GO is the prize).** Relative to the Paper 0 packing
construction, the inner-factor residue ($N_{\rm gen}$, inner KO-dim) is
packing-inaccessible. Two established facts compose: Door 4 seam (inner ring ⊥
forced outer ring) + Paper 0 §VII (packing is kinematic — labels + topology, no
operators/real-structures/couplings/multiplicities). $N_{\rm gen}$ is a
Hilbert-space multiplicity, inner KO-dim a $(J_F,\gamma_F)$ datum — both of the
type packing doesn't emit. The single escape (a packing-native factor-3
multiplicity) is empty: Paper 0's inventory is $2\ell+1$ (grows with $\ell$, not
generation-universal) + one spin-$\mathbb{Z}_2$. Explains why H1/W3/Koide all
failed (wrong side of the seam).

**Paper edit (Paper 32 §VIII.C).** Placed as a `\paragraph{}` after Door 4f
(which had named this exact "deep wall"), framed as a **relative** theorem of
necessity with the named sibling-axiom escape — house style, not a formal
theorem (consistent with Door 4 honest-scope). 62 pp, two-pass clean. The
initial `\begin{theorem}` draft was downgraded to paragraph form and removed.

## 4. Direction 3 — Paper 54 two-body forward scoping

**Verdict: STOP — and the sub-agent's "erratum" was the artifact.**

The agent flagged Paper 54 Prop 2's $n{=}3$ connected fraction (32.4%) as a stale
wrong number, recommending replacement with 71.4%. Main-session recompute
(`paper54_recompute_both_constructions.py`) ran BOTH gauge-field constructions:

| Quantity | Published | DOUBLE sum $\Sigma_i\Sigma_j$ | DIAGONAL $\Sigma_i$ (eq:A_full as written) |
|---|---|---|---|
| Prop 2 conn n=2 / n=3 | 76.7% / 32.4% | **76.7% / 32.4%** ✓ | 84.3% / 71.4% |
| Thm 2 Gaunt n=3 | 100% | **100%** ✓ | 97.3% |
| Thm 2 pure k=0 | both | **k0=100% both** ✓ | k2=2.5% at n=3 |
| Prop 3 Pearson n=2 / n=3 | 0.58 / 0.41 | **0.5818 / 0.4080** ✓ | 0.42 / 0.097 |
| Prop 3 sign n=2 / n=3 | 68% / 54% | **68% / 54%** ✓ | 52% / 55% |
| Table I n=3 (ε=.1/.5/1) | .8384/.4102/.1854 | **.8384/.4102/.1854** ✓ | .9907/.9544/.9112 |

The **double sum reproduces every published number bit-for-bit.** The paper's
numbers are correct and self-consistent. The agent had computed the **diagonal**
sum — `eq:A_full` as literally written — a *different* construction, producing
its phantom erratum. The real defect is the opposite and far smaller:
`eq:A_full` was mis-transcribed as the diagonal single-sum; it should be the
double sum $\Sigma_i\Sigma_j$, which is also the mathematically correct inner
fluctuation for the product algebra $\mathcal{A}\otimes\mathcal{A}$.

**Paper edit (Paper 54).** Corrected `eq:A_full` to the double sum + one sentence
pinning all results to it. All numbers stand. 4 pp, two-pass clean, 0 warnings.
Added a correction banner to `paper54_two_body_forward_scoping_memo.md` so its
withdrawn erratum can't mislead. The STOP verdict (don't pursue an
implementation sprint — radial weights are the Fock conformal-factor wall,
Pearson 0.58/0.41) **stands**; only the numerical sub-claim was wrong.
Literature: Bochniak–Sitarz 2022 (arXiv:2201.03839, spectral-action two-geometry
interaction = bimetric, never Coulomb) independently corroborates the negative;
Vanhecke 1999 (math-ph/9902029) is product-triple prior art. Citations not yet
added to the paper (optional follow-on).

## 5. Methodology lesson (durable)

A sub-agent reproducing a paper's number can silently use a **different
construction** than the paper (diagonal vs double-sum gauge field; under- vs
adequately-resolved PSLQ), producing a phantom erratum or a phantom negative.
Before a reproduction drives an edit to a released paper, confirm the
construction/precision matches the paper's — this is the verify-the-verifier
discipline of the audit cycle, and it fired productively in **both** directions
this session. (Captured here rather than a new memory file: MEMORY.md is over its
size limit, and this is an instance of an already-baked-in lesson.)

## 6. Honest scope

**Theorem grade / bit-exact:**
- Paper 50 S⁷/S⁹/S¹¹ in-ring closure: PSLQ bit-exact, 260 dps (driver) +
  independent 200 dps (main session). The $S^7$ correction is verified three
  ways (value) and two ways (relation).
- Paper 54 number reproduction: bit-exact match of all published values by the
  double-sum construction (float64 precision).

**Structural / paper-level (not theorem-grade):**
- The seam "relative theorem of necessity" (Paper 32 §VIII.C) is a structural
  argument composing two established facts; explicitly **relative** to the
  current Paper 0 construction, with a named (un-handled) sibling-axiom escape.
  It is NOT a bit-exact computation and is framed as a `\paragraph{}`, matching
  the Door 4 honest-scope ("nothing at theorem grade — structural analysis").

**Numerical observation:**
- The Door 1 (B) recursion STOP: $c_d$ is a single-rung $S^5$ coincidence, not a
  dimensional constant. The $2^{-(d-5)/2}$ top-atom ratio is the keep-worthy
  by-product.

**Named open follow-ons:**
- Paper 54: optionally add the Bochniak–Sitarz / Vanhecke citations (verified;
  support the central negative). Not applied.
- Paper 50: optionally formalize the S⁷/S¹¹ ladder as a `tests/` test
  (currently verified via debug driver, matching Paper 50's existing S³/S⁵
  convention — no `test_paper50` exists).
- Direction 1 next rung (S¹³) would further harden the ladder; not run.
- The seam sibling-axiom direction (force $N_{\rm gen}$ from a new primitive)
  remains multi-year with no present handle — recorded, deferred.

**Verification / prohibitions:**
- No `geovac/` code modified → 18/18 topological proofs and accuracy benchmarks
  untouched, no regression risk.
- Hard prohibitions (§13.5) clean: no fitted params, no geometry-hierarchy
  change, no negative-result deletion, Paper 2 conjectural label untouched. The
  Paper 50 change is a §13.8-authorized correction of a claim contradicted by
  new evidence.
- Paper 32's 8 compile warnings are pre-existing (revtex `Note1` / cross-ref
  debt); this sprint's paragraph added no `\label`/`\ref` and introduced none.

## 7. Files

**Created:** `debug/door1_s11_graduation.py`+`.json`+`_memo.md`;
`debug/seam_packing_scoping_memo.md`;
`debug/paper54_full_sum_gauge_probe.py`(+data), `debug/paper54_radial_under_fullsum_probe.py`(+data);
`debug/paper54_recompute_both_constructions.py`+`debug/data/paper54_recompute_both.json`;
`debug/paper54_two_body_forward_scoping_memo.md`; this memo.

**Modified:** `papers/group1_operator_algebras/paper_50_cft3_partition_function.tex`
(§7 remark rewrite); `papers/group3_foundations/paper_54_tensor_product_two_body.tex`
(eq:A_full); `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
(§VIII.C seam paragraph).

**Removed:** `debug/seam_theorem_of_necessity_draft.tex` (superseded by the placed paragraph).
