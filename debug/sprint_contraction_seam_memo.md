# Sprint memo: the contraction seam behind Paper 60's two-centre metric

**Date:** 2026-09-12 · **Version:** v5.11.4 · **Trigger:** PI conversational
question — where does the two-centre overcompleteness enter *mathematically*,
and can the decompactification/transcendental pattern be brought to bear on it?

Probes `debug/p60_{contraction_seam,window_constant}_probe.py`; data
`debug/data/p60_{contraction_seam,window_constant}.json`; scan
`debug/lit_scan/contraction_seam_e3_memo.md`; backing
`tests/test_paper60_{contraction_window,mcentre_orders}.py`.

---

## 1. The question, and what was already settled

Overcompleteness does **not** enter through the basis. v5.11.2 measured the
Bessel deficit of a displaced Sturmian against the one-centre span and found it
plateaus at 0.380 / 0.696 / 0.907 for `kR = 1/2/4`, flat over `N = 16..256`.
It enters at **one point of the continuous dual variable**: the Shibuya–Wulfman
operator is multiplication by the translation phase `e^{ip·R}` on the Fock
sphere, and at `p = 0` that phase is 1 for every `R` and every `k`. That is
Ron–Shen's `||sigma||_inf = 1`, attained.

The PI's question was whether the corpus's decompactification vocabulary
(compact skeleton, released axis, transcendental toll) applies here. It does,
and the payoff is a **classification of the two π's**, not a removal of the wall.

## 2. Measured (five identities plus two follow-ups)

| test | statement | result |
|:--|:--|:--|
| A | sine basis = S³ zonal Gegenbauer in half-density form, `sin(a chi) = sin(chi) C^(1)_(a-1)(cos chi)` | exact, 4.9e-14 |
| B | Mehler–Heine `n^-1 C^(1)_(n-1)(cos(z/n)) -> j0(z)` | exponent −2.014 |
| C | `(1/4pi) int e^{ip.R} dOmega = j0(pR)` | 7.2e-15 |
| D | near-null direction's spread in the Fock polar angle | `n·rms(theta)` → ~3.13, rising |
| E | weight-independence of the conditioning law | see §4 |
| F1 | Richardson of `n·rms(theta)` | **3.14158** vs π = 3.14159 |
| F2 | `n^2 min<theta^2>` over the band, **no Bessel present** | **9.86949** vs π² = 9.86960 |
| F3 | `1 - sigma_max = (kR)^2 <theta^2>/24` on that direction | 0.3% at n=320, kR=0.5..4 |

So the conditioning law is an identity in the contraction window:

    1 - sigma_max  =  <1 - j0>  =  (kR)^2/24 · <theta^2>  =  (kR)^2/24 · pi^2/n^2

## 3. The result that was kept: two π's, opposite sides of the seam

**The `pi^2` of `eq:sigma_law` carries no Bessel content.** It is the KMS
constant `c_1`, fixed by the symbol's zero order together with the finite
section. Paper 18 tier: calibration, M2 (`pi^2·Q` half). **Corrected by C23
run #3 (same day):** the Dirichlet-eigenvalue reading, the extremal
Wirtinger–Sobolev problem, and the independence of `c_alpha` from `b` are all
Böttcher–Widom's, in the source already cited; and F2 is the `b == 1` case of
that theorem, i.e. a change of representation, **not** an independent route.
The independent-route claim below is withdrawn.

**The `(2pi)^-1/2` and `pi/4` of `eq:chirp_decay` are continuum-side** —
normalisation and branch phase of a Bessel asymptotic at the `p -> infinity`
pole, the decompactified direction. Paper 18 tier: calibration, M2 (`sqrt(pi)·Q`
half).

**WITHDRAWN, same day (C23 run #3).** This section originally drew an
operational corollary — that a truncation-side price is matrix-level and a
preconditioner reaches it, while a continuum-side price is symbol-level and no
congruence can touch it. Both halves are wrong. The v5.11.0 preconditioner is
built FROM the symbol (matching polynomial chosen to share the symbol's zero,
Serra), and preconditioning IS a congruence, replacing the symbol by `f/g`.
The surviving statement is about the symbol on both sides and turns on the
KIND of feature: a banded congruence cancels a finite-order zero but cannot
alter a decay class — a second-order zero at one pole, a `j^-5/4` envelope at
the other. **Open and load-bearing:** whether the two poles are independent
features or two faces of one is unsettled; `debug/lit_scan/``toeplitz_finite_section_memo.md` records both conclusions at different
points, and the split adopted here is the weaker, measured one.

## 4. The one genuinely new thing: the law is carried by the translation

Deform the metric by any smooth positive radial weight `W` on the Fock sphere —
intra block = finite section of `W`, cross block = that of `W j0`. The
generalised symbol is the quotient `W j0 / W = j0`, which does not see `W`.
Measured at `kR = 2`, `n = 160`:

| W | collapse `(1-sigma_max)(n/kR)^2` | exponent |
|:--|--:|--:|
| 1 (SW reference) | 0.4072 | −1.979 |
| 1 + 0.8 cos χ | 0.4123 | −2.006 |
| 2 + sin χ | 0.4107 | −1.989 |
| e^−χ | 0.4164 | −2.011 |
| **control: 1 + cos χ** (vanishes at χ=π) | **0.828** | −1.967 |

against π²/24 = 0.4112. **The control is the load-bearing half** — without it
an insensitive pipeline would pass. Independently re-measured by the scan agent
across eight weights, same verdict.

**Scope, stated in-paper because the overstatement is close by.** This is
**not** `V_0`-independence: a position-space-local `V_0` acts on momentum space
by *convolution*, not multiplication, and leaves this class entirely. What is
shown is narrower — the law is carried by the *translation phase*, not by any
momentum-space weight multiplying it. Consistent with the one measured
out-of-class case (v5.10.15: the L² overlap degrades at the same exponent with a
constant ≈1.4× larger, which is what a non-multiplication metric should look
like). **Whether any position-space `V_0` preserves the constant is OPEN.**

## 5. C23 run #2 — three verdicts, and one defect it found in live text

Scan memo: `debug/lit_scan/contraction_seam_e3_memo.md`.

- **C1, the SO(4)→E(3) contraction reading: PRIOR ART.** Díaz Martín &
  Pacharoni (arXiv:1807.03904) state it in exactly the Gelfand-pair form;
  lineage Inönü–Wigner (PNAS 39, 1953, 510), **Clerc, Studia Math. 57 (1976) 27**
  (the named Mehler–Heine theorem on symmetric spaces), Dooley–Rice.
  **Therefore NOT written into the paper.** It is expository only, available to
  the PI as a framing sentence *with* attribution. Two technical cautions if it
  is ever used: the degeneracy sits at `chi = pi` (the *antipodal* end) while
  Mehler–Heine is classically stated at `chi -> 0`, and the antipodal limit
  carries a parity factor, `n^-1 U_{n-1}(cos(pi - z/n)) -> (-1)^{n+1} j0(z)` —
  it converges along parities, not outright. Stated bare it is false.
- **C2, the "trivial character" reading: PRIOR ART, and already cited** (flat
  limit; Driscoll–Fornberg + Barthelmé–Usevich are in the paper). A relabel.
- **C3, the weight-independence of §4: ABSENT.** The one thing here that is ours.

**The scan's best catch was a live over-claim, not an attribution.** Paper 60
read "spanned by a fixed vector set that does not move with geometry" and "a
fixed rank-(M−1) rotation removes it", unqualified. The null *space* is
geometry-independent; the *rates* are not. Re-derived and re-measured locally
before editing (C23 hard rule):

`j0(pd) = 1 - (pd)^2/6 + O(p^4)`, so the order-`p^2` form on `1⊥` is `P D2 P`
with `(D2)_ij = d_ij^2`. For collinear centres
`d_ij^2 = h^2(i^2 1^T + 1(j^2)^T - 2 x x^T)` and `P` annihilates the two outer
terms from both sides, leaving `-2h^2 P x x^T P` — **rank one**.

| geometry | orders in p | rank(P D2 P) |
|:--|:--|--:|
| collinear M=3 | (2, 4) | 1 |
| collinear M=4 | (2, 4, 6) | 1 |
| equilateral M=3 | (2, 2) | 2 |
| bent water-like M=3 | (2, 2) | 2 |
| tetrahedral M=4 | (2, 2, 2) | 3 |

(mpmath dps=60, two independent p-ratios.) **Water's `A_1` is bent, hence
full-rank — which is why its measured table holds.** A *linear* polyatomic is
not, and the single `tri(1,2,1)` reaches only its one order-2 direction. The
lever is now scoped in-paper to M=2 and non-collinear M=3; the collinear case is
open and not claimed. External: Batenkov–Demanet–Goldman–Yomdin
(arXiv:1809.00658), exponent controlled by maximal cluster size, our M=2 their
ℓ=2 — abstract verified at source before the bibitem was added.

## 6. Owed / PI items

1. **New Paper 34 projection candidate — PI-gated, NOT added.** §III covers
   "Wigner D-matrix rotation between molecular centers" but nothing covers
   *translation* between centres, which is where `j0` enters. Under the
   tag-transcendentals STOP rule a candidate projection is flagged for §VIII
   review rather than written into §III. Flagging it here.
2. **The §3 M2 tagging is itself a new `[SYMBOLIC]` classification claim** and
   carries its own C23 trigger-2 (at-authorship) obligation. Unscanned.
3. **Shibuya–Wulfman (1965) remains UNVERIFIABLE** (Royal Society 403, no 1965
   preprint), as the 2026-09-11 scan also left it. So "do SW themselves give a
   group-theoretic reading of their integrals?" is still unanswered. Avery's
   canon likewise paywalled.
4. **Do not add a Serra citation** for the Toeplitz-pencil ratio-symbol result
   underlying §4 — verified numerically and by a monotonicity argument here, but
   no primary was opened.

## 7. Process note

The probe's own framing had to be corrected mid-sprint: the agent brief said the
exponent would survive a change of `V_0` "while the CONSTANT is not", which the
measurement contradicted in the *favourable* direction — the constant survives
too, for every weight non-vanishing at `p = 0`. The paper carries the corrected
(stronger) form. Worth recording because the wrong version reached a dispatched
agent before it reached a test.

## 7b. C23 run #3 -- the at-authorship trigger paid for itself immediately

Memo `debug/lit_scan/c23_run_003_m2_tagging_memo.md`. Run on the tagging
paragraph the day it was written, per the v5.11.2 scope change.

| claim | verdict |
|:--|:--|
| T1 `c_1 = pi^2` as a Dirichlet eigenvalue | **PRIOR ART** -- Boettcher-Widom, the source already cited |
| T2 band-limited second-moment minimum | **PRIOR ART** -- Wirtinger branch, NOT Slepian/Landau/Pollak |
| T3 the two-sided taxonomy | **ABSENT** as a taxonomy, but its corollary was wrong |

**The scan's value was a defect in prose written hours earlier.** The
removability corollary was false in both halves (see Sec. 3, withdrawn).
Verified locally before editing. Corrected at the owner and swept to all three
dependents.

One reported defect did NOT hold: the inline `c_1` attribution was flagged as
C20 bibitem-less, but the bibitem exists and the sentence cites it. Not
"fixed". One did: `bottcher_widom2005` was arXiv-only; published coordinates
added after Crossref verification, **without** the series volume number, which
was search-level only.

*Open direction the scan could not close:* Serra-Capizzano's *Practical Band
Toeplitz Preconditioning and Boundary Layer Effects* was unreachable (Springer
IDP redirect) and is the source most likely to state a truncation-vs-symbol
separation in the literature's own words. T3's ABSENT verdict should be re-run
against it before being leaned on.

## 7c. Primaries pass -- both targets moved a claim

Memo `debug/lit_scan/primaries_sw1965_toeplitz_pencil_memo.md`.

| target | verdict | consequence |
|:--|:--|:--|
| Shibuya-Wulfman 1965 | **SECONDARY-QUOTED** (abstract + full reference list; body still unread) | the **translation** identification is prior art on three counts; only the **symbol** survives as ours |
| Toeplitz pencil / ratio symbol | **REACHED** | Ahmad et al., *Numer. Algorithms* **78**(3) 867-893 (2018), Eq. (22): an identity for the tau/DST-I algebra. Our weight-independence result is prior art, re-tiered |

Plus: `eq:sigma_law`'s ~1% residue is **dominated by the `n -> n+1` grid
convention** (-0.99% -> +0.26% at `n=160`, `kR=2`), re-measured locally before
editing; an `O(1/n)` term survives in both conventions.

**Two reported defects did not hold** (the `avery2004` volume is already right;
the inline `c_1` attribution does cite). Recording this because the pattern
across three scans is consistent: **the attribution findings are reliable, the
"live citation defect" findings are roughly half right**, and every one must be
checked against the current file before acting.

### OWED, and it threatens a standing claim

**Monkhorst & Jeziorski, "No linear dependence or many-center integral problems
in momentum space quantum chemistry", J. Chem. Phys. 71, 5268 (1979).**
Unopened. A 1979 title of that exact shape bears directly on
`memory/avery_method_and_prior_art_gaps.md` items 2 and 3 -- "there is NO prior
art for secular-matrix norm growth" and "conditioning is a blind spot in the
whole Avery canon" -- which are load-bearing for Paper 60's novelty framing and
have already survived three scans. **Read it before repeating either claim.**

Also owed: Serra-Capizzano's *Practical Band Toeplitz Preconditioning and
Boundary Layer Effects* (Springer IDP redirect, unreachable), the source most
likely to state a truncation-vs-symbol separation in the literature's own
words. C23 run #3's T3 ABSENT verdict should be re-run against it.

## 7d. The antipodal parity, and a gap this session opened

Withdrawing the independent-route claim (Sec. 7b) left "the minimiser is the
Dirichlet ground state in the band index" in the paper with **no backing test**
-- a coverage gap created by this session. Closing it caught an omission in the
claim itself.

The identification holds only with the antipodal parity:
`c_a ~ (-1)^(a+1) sin(pi a/(n+1))`, correlation 0.9999998 at n=320, while the
**unalternated mode is exactly orthogonal** (1e-7). Not an approximation of the
right answer -- its complement. Mechanism: the basis index is `chi`, the
Dirichlet mode is natural in `theta`, `sin(a chi) = (-1)^(a+1) sin(a theta)`.

**Same factor, second appearance.** C23 run #2 flagged exactly this parity as
the caution on the contraction reading of `j0`. It arrived here from a
different direction entirely, which is mild evidence the two readings are
describing one object.

Paper sentence made precise; guards added asserting both halves and fire-tested
three ways (drop the parity in either guard; feed the second band mode).

## 7e. Monkhorst-Jeziorski impact set, enumerated BEFORE the verdict

Built while the read was running, so that whichever way it goes the sweep is
already scoped. Four loci carry the exposure, and only two are genuinely at
risk.

| locus | what it says | exposure |
|:--|:--|:--|
| `memory/avery_method_and_prior_art_gaps.md` items 2-3 | "no prior art for secular-matrix norm growth"; "conditioning is a blind spot in the whole Avery canon" | **HIGH** — the direct target; already SUSPENDED pending the read; auto-loads every session |
| Paper 60 `sec:obstruction` (~L186) | shared-scale Coulomb Sturmians "grow linearly dependent" at common `k` | **MEDIUM** — a measured `L^2` fact; needs a scoping clause only if M-J's claim is about a DIFFERENT inner product |
| Paper 60 `sec:quantum` (~L673) | "no existing quantum algorithm combining (i) Sturmian basis, (ii) quantum eigenvalue routine, (iii) isoenergetic inversion" | **LOW** — scoped to *quantum algorithms*; a 1979 classical paper cannot reach it |
| group2 synthesis L737-740 | "measures rather than assumes its cost"; "the sharpest ... basis-independent conditioning" | **LOW** — measurements and a comparative, not novelty claims |

**The resolution the corpus should test first**, because it is the one its own
recent work predicts: GeoVac's linear-dependence statement is about the **L^2**
overlap, while the momentum-space method's natural inner product is the
**V_0-weighted** one — and v5.10.13/v5.10.15 established that the SW matrix IS
`V_0`, with the L^2 metric cancelling identically in the metric-free posing. If
M-J's "no linear dependence" is a statement in the momentum-space metric, both
claims can be true at once and the corpus already owns the reason.

*But that resolution is incomplete as it stands*, and the gap is the
interesting part: GeoVac ALSO measures the SW metric itself — i.e. `V_0` — to be
ill-conditioned at two centres (`cond ~ n^2`). So if M-J are making a
`V_0`-metric claim, either they are single-centre, or their construction avoids
the `p = 0` degeneracy in a way this corpus has not identified. **That second
branch would be a lever rather than a correction**, which is why the read is
worth its cost regardless of the attribution verdict.

## 7f. Monkhorst-Jeziorski 1979 -- verdict, and the resolution

**SECONDARY-QUOTED.** Abstract + record verified (Crossref, re-verified by the
PM); two-page body closed, zero repository copies, UNREAD. Correct DOI
`10.1063/1.438337`, JCP **71**(12), 5268-5269. Memo
`debug/lit_scan/monkhorst_jeziorski_1979_memo.md`.

| suspended claim | verdict |
|:--|:--|
| "no prior art for secular-matrix norm growth" | **SURVIVES**, narrowed; rests on abstract + ref list + page count |
| "conditioning is a blind spot in the whole Avery canon" | **RETRACTED as phrased** |

**The predicted resolution (Sec. 7e) was WRONG.** The PM expected an
`L^2`-vs-`V_0` metric distinction. It is the *same pencil*; the degeneracy is in
their matrix too. The difference is that they never **invert**: the overlap
enters only multiplicatively, so a near-null direction gives a spurious branch
rather than amplified error, and the basis is exactly orthonormal in the metric
used. **They keep the degeneracy out of the denominator.** GeoVac inverts
because a block-encoding wants a standard Hermitian eigenproblem -- so the
conditioning exposure belongs to the ENCODING REQUIREMENT, not the basis or the
metric. This is a better answer than the one that was predicted, and it is the
single most useful thing the whole prior-art arc produced.

**Open lever:** determinantal re-posing removes the conditioning multiplier and
reinstates the outer nonlinear scale search `eq:secular` exists to eliminate.
Unpriced.

**Impact set (Sec. 7e) verified against the verdict:** the memory items were the
real exposure and both moved; Paper 60 `sec:quantum`'s novelty claim was
insulated exactly as predicted (it is scoped to *quantum* algorithms, which a
1979 classical note cannot reach); `sec:obstruction` gained the resolution
rather than needing a scoping clause. Pre-enumeration held.

---

## 8. Follow-on items (2026-09-12, PI-directed)

### 8.1 Paper 34 candidate (e) logged in §VIII -- DONE, still PI-gated

Verified first that the gap is real: the string "translat" appears **nowhere**
in Paper 34, and the only Shibuya-Wulfman mentions (§III.22, multipole
expansion) name the SW *basis* expansion as a **contrast** -- something that
does NOT terminate -- never as a projection. §III.11 covers Wigner-`D`
*rotation* between centres and explicitly preserves rationality up to
`Q[sqrt2,sqrt3,sqrt6]`.

Candidate (e) is now logged in §VIII's open-question list alongside (a)-(d),
**not** written into §III, per the tag-transcendentals STOP rule. It carries:
the two-sided M2 signature (truncation prices `pi^2`, symbol asymptotics price
`sqrt(pi)`); the caution that the contraction reading is prior art and would be
cited rather than claimed, with the antipodal parity caveat; and a clause
reconciling it with the Sprint-3 "structurally complete" verdict -- that verdict
is scoped to the master Mellin engine's *accounting*, and (e) is M2 on both
sides, so what it adds is a **slot, not a mechanism**.

Promotion to a §III entry remains a PI call. group6 gates all PASS.

### 8.2 A gate for the splice class -- BUILT, PROPOSED, not adopted

`debug/qa/check_prose_continuity.py` (proposed **C24**). The 2026-09-12
self-catch was that an applier's anchor ended mid-sentence, and **every gate
passed** -- C10 compiles references, not prose. The class is not new: the
criteria document records that C20's own registration block was spliced into
the middle of the C19 sentence. Script-driven paragraph insertion is the
corpus's standard editing mechanism, so the class is systematic.

Detects two independent signals, each on its own: a prose paragraph that ends
without terminal punctuation, and one that opens with a lowercase ordinary
word. Neighbour-aware, so prose running INTO a display equation and prose
resuming AFTER one are both exempt, as is prose introducing a displayed
theorem.

**Measured: 0 findings across all 70 papers in 9 scopes**, while firing on both
halves of the real defect when it is re-planted into Paper 60.

*Two corrections it needed, both found by `debug/qa/_prose_gate_probe.py`, which
re-plants the real defect and requires the gate to fire.* The first draft
reported PASS on three scopes with zero findings and **would have shipped as a
gate that examines nothing**:

1. Environment tracking counted `\begin{document}`, which never closes until
   the last line, so every paragraph in the body was skipped as "inside an
   environment".
2. The conjunction was the wrong shape. The draft required one seam to BOTH end
   unterminated and be followed by a lowercase opener, on the reasoning that
   this keeps false positives near zero. But when a paragraph is spliced into
   the middle of a sentence, the two halves land at **opposite ends of the
   insertion**. The real defect triggered neither half.

Refinement then took the corpus from 107 findings to 0 without losing the
detection: inline `\begin{smallmatrix}` was disqualifying whole prose
paragraphs; preamble macro blocks were being judged as prose; trailing
`\checkmark` left whitespace the closer-strip did not remove.

**Adding a QA criterion is a gate change**, hence a minor version and a PI
decision (the precedent is C23 at v5.11.0). Proposed, not adopted; the script
and its probe stand on their own until then.
