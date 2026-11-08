"""STEP 8: chronicle the run. CHANGELOG entry + CLAUDE.md Sec.2 one-liner + version."""
import io

ENTRY = """## [v5.10.16] - 2026-09-11

**`/qa paper_60` DELTA verification = DEFECTS, remediated. NOT certified.** Five dimensions, unseeded, against criteria frozen with PI approval. Deterministic layer 12/12 green. **2 LARGE + 25 MATERIAL-SMALL + 5 upgrade candidates — and zero mathematical defects.** Five adversarial passes, one of which independently re-derived every keystone and fire-tested the guards it was handed on an 8-plant sample of its own construction, found nothing wrong with the arithmetic. Every finding was prose, attribution, staleness, or coverage. That is the same distribution as the v5.4.4 arc and the same lesson: *the arithmetic was never wrong; the misses are all in what got said about what was measured.*

### The LARGE was in CLAUDE.md, and it was mine

`CLAUDE.md:119` asserted four retired values — `K^0.84`, "local slope 0.906 by K=340", "the nuclear diagonal T0 is the sublinear block (K^0.70)", "T' is SUPERlinear (K^1.05)" — with no supersession marker, **in the file loaded into every session and every subagent dispatch**, two lines below the v5.10.12 bullet that supersedes it. §13.11 rule 9 ("status updates replace, never append") was not applied: four new §2 bullets were added for v5.10.12–15 and the v5.10.10 bullet was left standing. Every subagent dispatched this week read the inverted mechanism as current; that they still returned correct answers is luck, not process. Replaced; superseded text relocated verbatim to the frontier archive.

### Three instrument findings, fixed before any content edit

**(1) No C16 entry could reach the high-traffic surfaces.** `files` lists were papers-and-docs only, so `CLAUDE.md`, `claims_register.md`, `code_architecture.md` and the tracked `geovac/` docstrings were unreachable. **(2) `p60-sublinear-as-regime` stated its own CORRECTION in retired values** — it said the fix was "T^0 at K^0.70 (stable), T' at K^1.05 (SUPERlinear)", both retired the next day. A gate whose correction text is wrong teaches the wrong answer to everyone who trusts it. **(3) Its `exempt_if_nearby` contained `nuclear diagonal|T^0`** — the vocabulary of the retired mechanism — so a locus asserting "the nuclear diagonal T^0 is the sublinear block" **exempted itself by its own wrong mechanism**. That is the authored-exemption failure declared fixed-as-a-class on 2026-09-08, recurring.

Two new entries (`p60-tprime-superlinear`, `p60-l2-metric-diverges`) use `exempt_if_nearby = (?!)` and the standardized per-entry marker only. That choice is forced rather than stylistic: **CLAUDE.md:119 sits two lines from a bullet containing both "WITHDRAWN" and "2026-09-08"**, so any plausible retirement vocabulary in a ±5-line window would have exempted the very LARGE the entry exists to catch.

Both proven to discriminate **two ways — 6/6 fire on retired wording, 8/8 silent on corrected wording.** The silent direction earned its keep: an early draft of my own pattern fired on the *denial* ("Neither $T'$ grouping is superlinear") at two loci, i.e. on the corrected text. A guard that fires on the correction is worse than no guard — it makes the right answer unwritable, and the cheapest escape is to reword the correction until the gate goes quiet.

### The one wrong number

The s-only span deficit read `4.43` mHa at K=136. The driver's own output gives **4.40376 at K=136 and 4.43413 at K=105** — the free value quoted at one basis and the locked value at another, a mismatched pair in the sentence carrying the floor's re-attribution to the scale lock. It survived because it was an **unregistered literal**, invisible to C21 — and because the declared-debt table written into the DoD this session asserted "Every one is currently CORRECT" of a list that did not include it. *An enumeration offered as complete is a stronger claim than the literals it lists.* Both halves now registered as a matched pair.

### Two upgrades: the backing proved more than the prose

**`eq:W_diagonal` is exact, not measured.** Three cases: angular orthogonality; hermiticity of `(T−E)` on the two Sturmian equations giving `(Q_mu − Q_nu)·W = 0`; and disjoint n-multisets in the degenerate branch, which first occurs at `n_max = 35` against a largest computed basis of 17. The identity is unconditional. And it is **the same generalized-eigenvalue orthogonality the paper already invokes molecularly** ("automatic for eigenvectors of a common Sturmian problem with distinct beta") — the atomic case is that statement's own special case, measured numerically for want of noticing it was already argued eighty lines later.

**The variational bound is proved by inertia, and proves more.** `H(lambda) + ½lambda²S = lambda(lambda·1 − M)` **exactly** — verified to 1.4e-17, machine precision, an algebraic identity rather than a fit. With S ≻ 0, Sylvester's law of inertia gives the number of pencil roots below `−½lambda²` as `#{k : lambda_k(M) > lambda}`: zero at `lambda_max` (the bound), exactly k at `lambda_k` — the **root-by-root correspondence** every excited-state number in §4 depends on and that the test only asserted numerically. Confirmed against a direct pencil solve at k=0,1,2,3.

### The sweep

Four live zombies fixed — including two docstrings in tracked `geovac/`, one of them (`_whiten`: "ill-conditioned by construction ... this truncation is required, not cosmetic") **introduced by commit c435653, the same commit that corrected the claim next door**, and independently false since the truncation provably never fires (0 of 78/105/290 directions dropped). Sixteen chronicle loci given the standardized marker. The six-locus "ill-conditioned" cluster swept claim-wide — its sharpest instance called cond 5.8 ill-conditioned while §2 called cond 56.3 ordinary, same paper, factor of ten, opposite verdicts. `reverses to superlinear` corrected at owner + citer + register (no measured full-shell point is superlinear: the total exponent rises only 0.867→0.911). The floor/bracket contradiction resolved — 6.44 is a full-range fit that falls outside its own windowed bracket [6.47, 6.62], and only the ground state was inconsistent because 2¹S quotes a measured endpoint while the ground state quoted a fit, which is exactly what the paper's own stated policy forbids. The Acknowledgments' four surrendered theorems reclaimed. Two over-tight identity pins loosened from 1e-8 — *tighter than the quantity's spread across legitimate domains*, so raising the box rule would have failed them, the same anti-pattern that retired this paper's `cond(S)` guard — and fire-tested in both directions (fires at 5e-7, accepts 1e-9).

### Owed

**The second LARGE is a coverage gap, and it is owed as its own pass:** six abstract-level `[MEASURED]` families (the span-deficit pair, the free-scale matched set, the K=452 state pair, posing-cost roots 2–3, the state-prep overlap, the floor brackets) have driver-only backing in the prunable `debug/` tree. The reviewer independently reproduced every one, so the exposure is regression protection rather than correctness — and §9 requires guard-writing to be separate, separately-reviewed work, so it is not bundled here.

### Also

`cor:dual_p0` said "no shared p0 exists **that simultaneously represents both atomic length scales**" — a representation failure — which the auto-loaded `memory/polyatomic_state_of_play.md` rendered as "prove no shared-p0 Sturmian basis **exists**". Paper 60 *builds* one. Same unscoped-guardrail class that commit `3a3903a` fixed in Paper 58; the sweep had stopped at Paper 58, so `paper_fci_molecules`' flat "mathematically falsified as a viable discrete framework" was still live and is now scoped.

"""

c = io.open("CHANGELOG.md", encoding="utf-8").read()
A = "## [v5.10.15] - 2026-09-09"
assert A in c and "v5.10.16" not in c
io.open("CHANGELOG.md", "w", encoding="utf-8").write(c.replace(A, ENTRY + A, 1))
print("CHANGELOG: v5.10.16 added")

m = io.open("CLAUDE.md", encoding="utf-8").read()
m = m.replace("**Version:** v5.10.15 (September 9, 2026)",
              "**Version:** v5.10.16 (September 11, 2026)", 1)
B = "- **Molecular Sturmian build: exact, and it ends at a known wall (2026-09-09, v5.10.15):**"
NEW = ("- **/qa paper_60 DELTA = DEFECTS, remediated (2026-09-11, v5.10.16):** 2 LARGE + 25 SMALL, "
       "ZERO mathematical. The LARGE was this file's own Sec.2, and C16 could not reach it. "
       "See CHANGELOG v5.10.16.\n")
assert B in m and "v5.10.16):**" not in m
io.open("CLAUDE.md", "w", encoding="utf-8").write(m.replace(B, NEW + B, 1))
print("CLAUDE.md: version -> v5.10.16, Sec.2 one-liner added (29 words)")
