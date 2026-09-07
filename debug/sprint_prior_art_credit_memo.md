# Sprint memo — prior-art credit + the DLMF spheroidal identity

**Date:** 2026-09-06 · **Version:** v5.10.7 · **Trigger:** PI conceptual question
**Canonical memo for this sprint** (§13.11 rule 1). Agent deliverables: `debug/lit_scan/*.md` (4).

---

## 1. What was asked

The PI, reasoning from the v5.10.2 decompactification arc, asked whether
de-compactifying a "polar coordinate" axis appears elsewhere in quantum
physics, and then: *"one wonders how this could be an unexplored, unifying
principle."* Four adversarial external-literature scans were dispatched at PI
direction, on the standing `debug/lit_scan/` convention (assume prior art
exists; verdicts COLLISION / ADJACENT / BACKGROUND).

**Answer to the framing question: it is not unexplored, and it is not ours.
The principle is owned; the accounting is not.**

---

## 2. The four scans

| Scan | Memo | Headline verdict |
|:--|:--|:--|
| Abelian-residue principle | `abelian_residue_principle_memo.md` | NOT novel — Curie (1894) + Bethe (1929) + Pontryagin, composed. The abelian/non-abelian *asymmetry* alone is unstated as a general principle. |
| Compactness ⇒ discreteness | `compactness_discreteness_memo.md` | Named in mathematics (Pontryagin, compact resolvent, Peter–Weyl); unnamed folk theorem in physics; **stated causally by LQG** and, for time, by Dolce. |
| Deformation metrics / separability | `decompactification_metrics_memo.md` | Deforming labels = textbook special-function theory *with a monotonicity theorem*; decay-length front = Herring; quantum defect ↔ SO(4) already said. |
| Information-theoretic QM | `info_theoretic_compactness_memo.md` | WH7 anticipated in substance (Dolce, Rovelli); WH8 anticipated as a position (QBism). No reconstruction *derives* discrete outcome menus. |

Cross-scan convergence: each independently landed on the same small residue —
**the transcendental accounting, the re-targeting to sparsity, and two
unwritten bridges.**

---

## 3. PM verification record

Load-bearing citations were verified by the PM against primary sources, not
taken from agent memos (memory rule; the corpus has twice inherited citation
defects from explorer memos). **Two agent claims and three of the PM's own
in-conversation statements needed correction.**

### Verified as claimed
| Item | Route | Result |
|:--|:--|:--|
| Curie's principle, K̃ = K ∩ G | IUCr pamphlet 18 §8, fetched | Verbatim: "greatest common subgroup of the symmetry group of the crystal without the influence and of the symmetry group of the external influence" |
| Bethe 1929, Mulliken 1930/1932, Herring 1962, Herring–Flicker 1964, Erikson–Hill 1949 | Crossref `works/{doi}` | All six resolved; titles/volumes/pages/years match |
| Liu & Noui, CQG 34, 135008 | arXiv HTML full text | Sentence confirmed, **Introduction** |
| Dolce, arXiv:1707.00677 | arXiv abstract | Verbatim confirmed |
| Rovelli RQM Postulate 1 | SEP entry (rev. 2025-02-04) | Verbatim confirmed, both halves |
| Maz'ya–Shubin, Ann. of Math. 162, 919 | Crossref | Title is literally "Discreteness of spectrum and positivity criteria…" |
| Fall & Kondo 2604.27125; Lax 2604.27339 | arXiv abstracts | **Both real** — checked because 2026 IDs are the usual fabrication shape |

### Corrections made during verification
1. **Liu & Noui quote nearly reported as unverified.** The sentence is *not* in
   the abstract, and a PDF text-extraction pass failed to find it. Only the
   arXiv HTML full text settled it. It is real, and the paper states it as an
   *established reading in its field*, not as its own result — which makes the
   collision stronger, not weaker.
2. **Agent's Paper 18 §III flag was mis-aimed.** It recommended "compact
   resolvent, not compact manifold" *for §III*. §III's actual mathematics is
   Peter–Weyl (compact group) and Selberg (compact quotient), both correct as
   written. The real exposure is the CLAUDE.md §1.7 slogan "discreteness IS
   compactness" — a biconditional — and §III's section title. Fixed at the
   right level.
3. **PM error — "the Erikson–Hill constant."** Not an attested name. The paper
   is real (PR 75, 29 (1949)); write "the second-order integral of Erikson and
   Hill."
4. **PM error — the three-centre non-integrability bridge, as pitched, is
   refutable.** Knauf & Taimanov (arXiv:math/0401202 = Math. Ann. 331, 631
   (2005)) construct independent integrals of **Gevrey class g > 1** for n ≥ 3
   at positive energies; only real-analytic ones fail. The wall survives —
   Gevrey integrals are not Killing tensors and give no separation of variables
   — but the qualifier *analytically* is load-bearing. Also: the integrable
   planar case is **Euler's two** centres; the Darboux family runs 2 → 4,
   skipping 3.
5. **PM error — the quantum-defect proposal is mostly prior art.** Krug &
   Buchleitner (EPL 49, 176 (2000); arXiv:physics/9911064, titled *"Residual
   Symmetries…"*) state that a non-hydrogenic core destroys Runge–Lenz
   conservation. Only the specifically Fock-projection framing is unclaimed.

---

## 4. New result — the η-equation is the spheroidal equation

**[SYMBOLIC + MEASURED]** Captured in Paper 58, Sec. "The continuous side",
new paragraph *"The continuity is analytic, with a bounded rate."*

Matching the m = 0, homonuclear η-equation of Paper 11,
`d/dη[(1−η²)G'] + (−A + c²η²)G = 0`, against DLMF 30.2.1,
`d/dz[(1−z²)w'] + (λ + γ²(1−z²) − μ²/(1−z²))w = 0`, term by term:

> **γ² = −c², λ^m_n(γ²) = c² − A**

The framework's angular separation constant *is* a spheroidal eigenvalue.
Three consequences, each verified against the production solver:

| DLMF | Becomes | Measured |
|:--|:--|:--|
| (30.3.3) λ(0) = n(n+1) | A(0) = −n(n+1) | exact to 1e-12, all m ≤ 2, n ≤ 4 |
| §30.3(ii) analytic in γ² | coupling **real-analytic** in c² | — (structural) |
| (30.3.4) −1 < dλ/dγ² < 0 | **0 < dA/d(c²) < 1** | [0.14, 0.87], 9 modes × 40 pts, **0 violations** |
| (30.3.1) ordering | doublets at large c² = **separated-atom limit** | gap ~ c²e^{−2c}; c²=140 within 1e-7 vs between O(10) |

**Both ends of the correlation diagram are now pinned.** (30.3.3) fixes the united-atom end to the Legendre
label; the doublet collapse fixes the separated-atom end. The gap law is measured as ~c²·e^{−2c} — local
slope −1.69 → −1.89 approaching −2 from above, prefactor c^2.07 at R²=0.999. **A naive single-exponential
fit returns −1.67, which is the prefactor contaminating the slope, not a decay rate** — the backing test
asserts the two constraints separately for that reason. Since c ∝ R this is Herring's exchange-splitting
class, i.e. today's front citation arriving at the opposite end of the diagram. Reported as measured on
this branch; the splitting quantity is λ, not an energy, so no literature asymptotic is claimed.

**Why it matters.** v5.10.2 established "Boolean label, continuous coupling" by
measurement. This upgrades the second half from *measured smooth* to
*real-analytic with a two-sided rate bound*: l-decompactification can neither
stall nor outrun c². It also pins the R → 0 end of the correlation diagram to a
special-function identity rather than an assertion.

**Honesty note (carried in-paper).** DLMF states (30.3.4) for the real variable
γ² with no explicit domain restriction; our branch is γ² = −c² < 0, the oblate
branch. The bound is therefore reported as **verified numerically there**, not
cited as proved there.

**Independent-route by-product.** `molecular_sturmian._angular_sep_const`'s
docstring has asserted `A = -obl_cv(m, n, c)` since the module was written and
was never tested. Now verified against scipy at **2.3e-13** across the grid.

Backing `tests/test_paper58_dlmf_spheroidal.py` (28 tests, 0.8 s). Driver
`debug/dlmf_spheroidal_crosscheck.py`. Five guards, all fire-tested against
**subject** plants (not test-side edits): c²-sign flip → rate bound fires;
wrong Legendre diagonal → united-atom limit fires; basis detune → scipy route
fires.

*Incidental:* the H-matrix build is duplicated between `_angular_sep_const` and
`_angular_eigenvector` (fire_test reported "anchor occurs 2x"). Harmless; not
refactored.

---

## 5. Instrument fix — C21 read bibliographies as claims

Adding the Maz'ya–Shubin bibitem made C21 FAIL: the page range **"919--942"**
matched the retired `cah_rel_n2_pauli` value 942. A bibliography entry's
volume, pages and year are metadata about *someone else's* paper and can never
be an assertion of this corpus. `check_retired` now skips `thebibliography`
spans.

Fire-tested both directions, which is the point — the risk of the fix is
over-broadening: a retired value planted in Paper 18 **body** text still FIRES;
the same value in a bibitem stays silent; the restored tree PASSes.

---

## 6. Verification summary

- 12 deterministic gates; C16 (trunk/group2/group3) PASS, C17 (same) PASS,
  C21 PASS post-fix, C13/C14/C15/C19/C20/C22 PASS.
- Papers 58 (14 pp) and 18 (30 pp): 3-pass compile, **0 errors, 0 undefined**.
- Regression slice 48 passed (18 symbolic S³ proofs + new DLMF file + v5.10.2
  front/η tests). Registry self-test 17/17.
- Working tree confirmed free of fire-test mutations.

---

## 6b. Second instrument pass — and a guard disarmed by a dict literal

Two further C21 collisions, both from this session's own new prose:

1. `exponent` matched **inside the word "exponential"**, colliding the new
   log-gap slope −1.69 with retired `exp_lambda_4pt`. Fixed with a
   `spheroidal|DLMF` forbid rather than by dropping a detection anchor;
   verified the retired value still fires in a genuine 1-norm context.
2. **Latent bug found while fixing (1):** `exp_lambda_4pt` appeared in
   `RETIRED` **twice**, as `1.69` and `1.690` — the same float. The dict
   literal silently collapsed them and the later entry's overwrite
   **discarded the earlier one's `alpha|lambda` anchors with no error**.
   The "guard that cannot fail" class, arriving through a dict literal
   rather than an assertion. Entries merged; a new self-test scans the
   *source* (the duplicate is gone by import time) and names the offending
   literals. Verified against a genuine planted duplicate — clean assertion
   failure. 57 literal keys, 56 distinct: one had been silently lost.

*Process note:* two of this session's "discrimination checks" initially
reported the guard asleep because the **plant anchor did not match** and
nothing was actually planted. Both were caught and redone. A plant that
fails to apply reads exactly like a passing gate.

## 7. Owed / PI calls

1. **WH7 claim-level restatement (PI call).** The interval-vs-circle and
   real-vs-imaginary-time conflation is a genuine hole in the hypothesis *as
   written*. Status line updated; the claim and falsifier were not touched
   (§1.7 governance requires PI direction).
2. ~~CLAUDE.md §1.7 organizing observation~~ — **DONE** (PI-directed, same
   session): corrected to "compactness forces discreteness; the converse is
   false", with the compact-resolvent form, the Paper 24 counterexample, and
   credit to Liu–Noui / Rovelli / Dolce.
3. ~~Paper 58's headline vs its name~~ — **DONE** (PI-directed, same session):
   two scope statements added to the abstract rather than a retitle, since the
   paper carries a DOI and is cited as `loutey_paper58`.
4. **Two unwritten bridges**, both unclaimed and unstarted: Lüscher ↔ quantum
   defect theory; physics compactness ↔ periods-as-compact-volumes
   (Kontsevich–Zagier: every period is the volume of a compact semi-algebraic
   set).
5. **Do not use the Efimov "log-radial compactification" reading** until the two
   compactness notions in it are separated (periodicity mod ln λ₀ sets the
   *ratio*; a finite log interval sets the *count*). Pazy (PRE 102, 022136)
   already maps Efimov to Bloch states.
6. Unverified items in the agent memos are labelled UNVERIFIED there and must
   not enter a bibliography without a primary-source check. Nothing UNVERIFIED
   was cited in this sprint.
