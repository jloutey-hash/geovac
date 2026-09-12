"""Close the coverage LARGE in the DoD, the CHANGELOG, and CLAUDE.md Sec.2.

Six of the seven families are now BACKED-SOUND against tracked `geovac/` code.
The seventh -- the floor bracket -- is PARTIAL by design, and the DoD says which
half is backed and which is not, because a record that rounds PARTIAL up to done
is the `.done.md`-ratifying-a-retired-claim class in a different costume.
"""
import io

# ---------------------------------------------------------------- DoD
D = "docs/qa/paper_60.done.md"
d = io.open(D, encoding="utf-8").read()

OLD = """> **LARGE 2 (coverage).** Six abstract-level `[MEASURED]` families have
> driver-only backing in prunable `debug/`. Independently reproduced by the
> reviewer, so the exposure is regression protection rather than correctness.
> **OWED as its own pass** (§9: guard-writing is separate, separately-reviewed work)."""
NEW = """> **LARGE 2 (coverage) — CLOSED 2026-09-11**, as its own pass per §9. Six
> abstract-level `[MEASURED]` families had driver-only backing in prunable
> `debug/`; all six are now recomputed from tracked `geovac/` code in
> `tests/test_paper60_resource_ladder.py` (7 tests, all `@pytest.mark.slow`,
> 14 min). Four were fire-tested against the specific wrong answer each excludes
> — including the K=105/K=136 substitution that was this run's one wrong number,
> which FIRES. The seventh family, the floor **bracket**, is **PARTIAL and
> declared**: the claim form (fit-from-below vs Shanks-from-above) is backed on
> the spdf ladder K=74–244; the endpoint values [6.47, 6.62] / [1.647, 1.676]
> need the full ladder to K=452 (~20 min) and remain driver-backed, recorded as
> PARTIAL in `docs/claim_test_matrix.md`.
>
> Two things the writing found that the review had not. **(i)** The
> bracket-direction claim is SECTOR-SPECIFIC: written first on the cheap s-only
> ladder, it failed, because there the windowed fit *falls* (4.3098 → 4.3059 →
> 4.3035) where on the spdf ladder it *rises*. A cheap proxy in the wrong sector
> would have reported the claim backed while measuring something that behaves
> oppositely. **(ii)** A seventh locus of the "ill-conditioned" cluster, in
> `docs/claim_test_matrix.md` and in a backing test **named**
> `test_paper60_l2_overlap_illconditioned_grows` — the test-asserts-the-zombie
> sub-flavour catalogued in v4.43.5. Its assertions were sound and are unchanged;
> the name and framing were not. C16's pattern was widened twice to reach the
> bare adjective (post- and pre-nominal), and re-proved silent on the paper's own
> denial."""
assert OLD in d, "DoD LARGE-2 locus not found"
d = d.replace(OLD, NEW, 1)
io.open(D, "w", encoding="utf-8").write(d)
print("paper_60.done.md: LARGE 2 marked CLOSED, with the PARTIAL half named")

# ---------------------------------------------------------------- CHANGELOG
C = "CHANGELOG.md"
c = io.open(C, encoding="utf-8").read()
OLD_OWED = """### Owed

**The second LARGE is a coverage gap, and it is owed as its own pass:** six abstract-level `[MEASURED]` families (the span-deficit pair, the free-scale matched set, the K=452 state pair, posing-cost roots 2–3, the state-prep overlap, the floor brackets) have driver-only backing in the prunable `debug/` tree. The reviewer independently reproduced every one, so the exposure is regression protection rather than correctness — and §9 requires guard-writing to be separate, separately-reviewed work, so it is not bundled here."""
NEW_CLOSED = """### The coverage LARGE, closed as its own pass

Six abstract-level `[MEASURED]` families — the span-deficit pairs, the free-scale matched set, the K=452 state pair, the posing-cost ladder, the state-preparation overlaps — had their only backing in `debug/p60_*.py`, a tree §9 prunes by design. All six now recompute from tracked `geovac/` code in **`tests/test_paper60_resource_ladder.py`** (7 tests, all `@pytest.mark.slow`, 14 min total). Every number was reproduced independently *before* the test was written, including the K=452 pair at 431 s.

Four guards were fire-tested against the specific wrong answer each excludes: planting the K=105 value where K=136 belongs — **this run's one wrong number** — FIRES; so do "whitening is free" (`Hh = H`), a state-independent posing cost (all roots → root 0), and the naive L²-amplitude substituted for the S-metric overlap.

**The floor bracket is PARTIAL, and declared as such.** Its claim form (windowed fit from below, Shanks from above) is backed on the spdf ladder K=74–244; the endpoint values [6.47, 6.62] and [1.647, 1.676] need the full ladder to K=452 (~20 min) and stay driver-backed. `docs/claim_test_matrix.md` carries that as a PARTIAL row rather than a silent omission.

**Two findings the writing produced that the review had not.** *The bracket direction is sector-specific.* Written first on the cheap s-only ladder, the test failed — there the windowed fit FALLS (4.3098 → 4.3059 → 4.3035) where on the spdf ladder it RISES (6.3887 → 6.4192 → 6.4388, reproducing the driver exactly). The approach direction belongs to the sector, not the extrapolator, so a cheap proxy would have certified the claim while measuring something that behaves oppositely. *And a seventh "ill-conditioned" locus* — in the claim matrix, and in a backing test **named** `test_paper60_l2_overlap_illconditioned_grows`, the test-asserts-the-zombie sub-flavour catalogued in v4.43.5. Its assertions were always sound (cond(S) 3.0→32.2, shared-scale λ inflating faster than hydrogenic) and are untouched; the name and framing asserted the withdrawn reading. Renamed, reframed, and C16 widened twice — post- and pre-nominal — then re-proved silent on the paper's own denial."""
assert OLD_OWED in c, "CHANGELOG owed-section locus not found"
c = c.replace(OLD_OWED, NEW_CLOSED, 1)
c = c.replace("## [v5.10.16] - 2026-09-11", "## [v5.10.17] - 2026-09-11", 1)
io.open(C, "w", encoding="utf-8").write(c)
print("CHANGELOG: Owed section replaced by the closure; entry -> v5.10.17")

# ---------------------------------------------------------------- CLAUDE.md
M = "CLAUDE.md"
m = io.open(M, encoding="utf-8").read()
m = m.replace("**Version:** v5.10.16 (September 11, 2026)",
              "**Version:** v5.10.17 (September 11, 2026)", 1)
OLD_B = ("- **/qa paper_60 DELTA = DEFECTS, remediated (2026-09-11, v5.10.16):** 2 LARGE + 25 SMALL, "
         "ZERO mathematical. The LARGE was this file's own Sec.2, and C16 could not reach it. "
         "See CHANGELOG v5.10.16.")
NEW_B = ("- **/qa paper_60 DELTA = DEFECTS, remediated (2026-09-11, v5.10.17):** 2 LARGE + 25 SMALL, "
         "ZERO mathematical. The LARGE was this file's own Sec.2, and C16 could not reach it. "
         "Coverage LARGE closed: 6 abstract families now tracked-backed. See CHANGELOG v5.10.17.")
assert OLD_B in m
m = m.replace(OLD_B, NEW_B, 1)
io.open(M, "w", encoding="utf-8").write(m)
print("CLAUDE.md: version -> v5.10.17, Sec.2 bullet updated")
