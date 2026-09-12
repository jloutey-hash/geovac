"""STEP 7: the definition-of-done, including the part of it that was wrong.

The claims reviewer found that the criteria FROZEN FOR THIS RUN carry the
mechanism the paper withdrew on 2026-09-07:

    C8.1  "... driven by overlap ill-conditioning."
    C10   dispatch line: "... L2 ill-conditioning growth; F0=5/8; ..."

So a reviewer following the DoD was instructed to audit a headline that the
owner document had already retracted -- the same class the /qa paper_61 run
named on 2026-09-08 ("a `.done.md` RATIFYING a retired claim"), recurring in a
record I rewrote and froze this session with PI approval.

The delegation convention adopted this session does not reach it, and that is
worth recording precisely:  the convention delegates VALUES to their owning gate
(C21 key / C17 family), and this is a MECHANISM stated in prose.  Prose
mechanisms have no owning gate, so nothing was delegated and nothing was
checked.  That is a measured gap in a convention one day old, not a
counterargument to it.

Also records the run and the two upgrades, and adds C8.15/C8.16 so the two
newly-proved results are criteria rather than folklore.
"""
import io

D = "docs/qa/paper_60.done.md"
d = io.open(D, encoding="utf-8").read()

# ---- 1. C8.1's retired mechanism.
OLD_C81 = ("1. **Naive L2 inflation [MEASURED].** Shared-scale Coulomb-Sturmian is L2-non-orthogonal;\n"
           "   JW/Löwdin LCU 1-norm inflates λ∼Q^3.33 vs Q^1.19 hydrogenic; driven by overlap\n"
           "   ill-conditioning.")
NEW_C81 = ("1. **Naive L2 inflation [MEASURED].** Shared-scale Coulomb-Sturmian is L2-non-orthogonal;\n"
           "   JW/Löwdin LCU 1-norm inflates λ∼Q^3.33 vs Q^1.19 hydrogenic. **Mechanism corrected\n"
           "   2026-09-11:** this line read \"driven by overlap ill-conditioning\" — the\n"
           "   characterization §2 of the owner WITHDREW on 2026-09-07 (converged cond(S) is\n"
           "   ordinary, ~0.12·K). The inflation is driven by the *density* of Löwdin's\n"
           "   S^{-1/2} and the resulting spread of the coefficient distribution, not by\n"
           "   numerical instability. Asserting ill-conditioning here = MATERIAL; C16\n"
           "   `p60-l2-metric-diverges` guards it.")
assert OLD_C81 in d, "C8.1 locus not found"
d = d.replace(OLD_C81, NEW_C81, 1)

# ---- 2. the C10 dispatch line that sends a reviewer at the same retired headline.
OLD_C10 = "  L2 ill-conditioning growth; F0=5/8; single-config He variational −2.847/−2.84766;"
NEW_C10 = ("  L2 overlap NON-ORTHOGONALITY growth (cond is modest — do NOT audit this as\n"
           "  ill-conditioning, that reading is withdrawn); F0=5/8; single-config He\n"
           "  variational −2.847/−2.84766;")
assert OLD_C10 in d, "C10 dispatch locus not found"
d = d.replace(OLD_C10, NEW_C10, 1)

# ---- 3. the two upgrades become criteria.
ANCHOR = "## Un-delegated literals — DECLARED DEBT (2026-09-11)"
NEW_CRIT = """## C8.15 / C8.16 — the two results upgraded from MEASURED to proved (2026-09-11)

**C8.15 — `eq:W_diagonal` is [SYMBOLIC], not measured.** The one-body Coulomb
metric is exactly diagonal, by three cases: angular orthogonality when
`l_μ ≠ l_ν`; hermiticity of `(T − E)` on the two Sturmian equations giving
`(Q_μ − Q_ν)·W_μν = 0` when the `l` agree and the roots differ; and disjoint
n-multisets in the degenerate branch, which first occurs at `n_max = 35` (the
largest basis computed is 17). The `5×10⁻¹¹` entrywise agreement is a check on
the implementation, **not the evidence for the claim** — stating it as the
evidence is now an UNDERCLAIM and MATERIAL. The equation is stated at unit scale;
at general λ the diagonal is `λ R_ν`, and dropping that qualifier is MATERIAL
(it makes the labelled equation false by a factor λ).

**C8.16 — the variational bound is proved by inertia, and proves more.** Because
W is diagonal, `H(λ) + ½λ²S = λ(λ𝟙 − M)` **exactly** (verified to 1.4e-17 —
machine precision, i.e. an algebraic identity, not a fit). S ≻ 0, so by
Sylvester's law of inertia the number of pencil roots below `−½λ²` equals
`#{k : λ_k(M) > λ}`: zero at `λ_max` (the bound), exactly k at `λ_k` (the
**root-by-root correspondence**, on which every excited-state number in §4
depends). Verified against a direct pencil solve at k=0,1,2,3. The paper must
NOT revert to asserting "E_iso is the lowest root" as an unproved intermediate —
it is a consequence, not a premise.

"""
assert ANCHOR in d and "C8.15" not in d
d = d.replace(ANCHOR, NEW_CRIT + ANCHOR, 1)

# ---- 4. record the run.
BANNER = "# Paper 60 (Generalized-Sturmian Secular Equation) — `/qa` profile"
RUN = """
> **RUN 2026-09-11 — DELTA verification, unseeded: DEFECTS, remediated. NOT certified.**
> Five dimensions. Deterministic layer 12/12. Citations CLEAN. **2 LARGE + 25
> MATERIAL-SMALL + 5 upgrade candidates, and zero mathematical defects** — every
> finding was prose, attribution, staleness or coverage. Five adversarial passes,
> one of which re-derived every keystone independently and fire-tested the guards
> it was handed, found nothing wrong with the arithmetic.
>
> **LARGE 1 (corpus surface).** `CLAUDE.md:119` asserted four retired values
> (`K^0.84`, "local slope 0.906", "T⁰ is the sublinear block at K^0.70", "T′ is
> SUPERlinear at K^1.05") with no supersession marker, in the file loaded by every
> session and every subagent dispatch — two lines below the bullet that supersedes
> it. §13.11 rule 9 was not applied: four newer bullets were appended and this one
> left standing. Replaced; superseded text relocated to the frontier archive.
>
> **LARGE 2 (coverage).** Six abstract-level `[MEASURED]` families have
> driver-only backing in prunable `debug/`. Independently reproduced by the
> reviewer, so the exposure is regression protection rather than correctness.
> **OWED as its own pass** (§9: guard-writing is separate, separately-reviewed work).
>
> **The one wrong number in the paper:** the s-only span deficit read `4.43` at
> K=136; the driver's own output gives `4.40376` there and `4.43413` at K=105 — a
> mismatched pair, in the sentence carrying the floor's re-attribution. It was an
> *unregistered* literal, so C21 was blind to it, and it was absent from the
> declared-debt table that asserted its list was all correct. Now registered as a
> matched pair (`p60_span_deficit_sonly_locked` / `_free`).
>
> **Instrument findings (three), all fixed before any content edit:** no C16 entry
> could reach `CLAUDE.md`, `claims_register.md`, `code_architecture.md` or
> `geovac/sturmian_variational.py`; `p60-sublinear-as-regime` stated its own
> CORRECTION in values retired the next day; and its `exempt_if_nearby` contained
> `nuclear diagonal|T^0` — **the vocabulary of the retired mechanism**, so a locus
> asserting "the nuclear diagonal T⁰ is the sublinear block" exempted itself. Two
> new entries (`p60-tprime-superlinear`, `p60-l2-metric-diverges`) use the
> standardized per-entry marker and `exempt_if_nearby = (?!)`, because CLAUDE.md:119
> sits two lines from a bullet containing both "WITHDRAWN" and "2026-09-08" — any
> plausible nearby vocabulary would have exempted the LARGE. Both proven to
> discriminate two ways: 6/6 fire on retired wording, 8/8 silent on corrected
> wording, including two loci where an early draft of my own pattern fired on the
> *denial* ("Neither T′ grouping is superlinear") and would have made the right
> answer unwritable.
>
> **Two upgrades** (C8.15/C8.16 below) — the backing proved more than the prose on
> both results this week rests on.
>
> Full sweep: 4 live zombies fixed, 16 chronicle loci marked, 6-locus
> "ill-conditioned" cluster swept claim-wide, `reverses to superlinear` corrected at
> owner + citer + register, the floor/bracket contradiction resolved, the
> Acknowledgments' four surrendered theorems reclaimed, and the `cor:dual_p0`
> overstatement fixed in the auto-loaded memory and in `paper_fci_molecules`.

"""
assert BANNER in d and "RUN 2026-09-11" not in d
d = d.replace(BANNER, BANNER + "\n" + RUN, 1)

io.open(D, "w", encoding="utf-8").write(d)
print("paper_60.done.md updated:")
print("  C8.1  retired 'ill-conditioning' mechanism corrected")
print("  C10   dispatch line no longer sends reviewers at a withdrawn headline")
print("  C8.15 eq:W_diagonal [SYMBOLIC] + scale convention")
print("  C8.16 variational bound by inertia + root-by-root correspondence")
print("  RUN   2026-09-11 DELTA recorded")
