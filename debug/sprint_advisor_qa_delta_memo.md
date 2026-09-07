# Sprint memo — the /advisor pass and its group2/group3 surgery, closed by a /qa DELTA (2026-09-06, v5.10.6)

## 1. Origin and shape

Started from a PI "broad strokes" question — where does GeoVac go from here, and is
there anything worth taking to James Avery or the wider community? It became a full
arc: a strategic synthesis, four external-literature scans, a new **thesis-advisor
skill** (`/advisor`), an advisor pass on Papers 58/59/60, the paper surgery that pass
drove, and a `/qa` DELTA that certified the diff.

## 2. The strategic read (the "point of it all")

The load-bearing synthesis, offered to the PI: **"discreteness is compactness"
(Paper 18) turned into a measured, continuous phenomenon** by the chemistry arc — you
can watch compactness release axis-by-axis as a bond forms (m stays compact, l's SO(4)
degeneracy dies), and the crossing is a continuous front (R\* ≈ 2.14 n/Z), not a wall.
The through-line across the whole program: **GeoVac reads standard chemistry objects
(Sturmian ERIs, corresponding orbitals, compound matrices, O(4) correlation) at the
discrete/continuum boundary — genus, transcendence weight, LCU 1-norm, a measured
decompactification front — a lens the field has the objects for but does not apply.**

Four adversarial lit scans (memos in `debug/lit_scan/`) calibrated novelty:
- **3-center ERI = elliptic Bessel moment / sunrise bridge → OPEN** (genuine frontier;
  the chemistry↔elliptic-Feynman literatures are disjoint; Avery 2013 is the closest
  adjacent, stops at numerics). Reaches a *new* community (Brown/Kleinschmidt periods).
- **Atomic metric-free Sturmian encoding, ‖M‖₁~K⁰·⁸⁴ → OPEN** (Koridon 2021's N¹·³⁴ is
  the best prior lever, still superlinear; Babbush 2019 is only a *wording* collision).
- **Projector-angle bonding → ADJACENT/COLLISION** — principal angles, the commutator,
  and compound matrices are all mature named machinery (Amos–Hall, King, Löwdin,
  Burton, Halmos, Böttcher–Spitkovsky). What is ours is narrow: the R-continuous front
  reading. This drove the citation-hygiene pass.
- **SO(4)-breaking / decay-length front → ADJACENT (framing) / OPEN (the constants)**;
  Solov'ev's hidden-crossings program is the adjacent prior art.

## 3. The `/advisor` skill (new, PI-requested)

Created `.claude/commands/advisor.md`: a thesis-committee-chair review of one specific
paper — is the contribution clear and defensible, is it honestly positioned, what is
the single hostile-examiner question it can't yet answer, where to push. It **advises,
never edits** (§13.10), carries no gates/PASS verdict (distinct from `/qa`), and catches
framing in **both** directions (overclaim *and* undersell). Ran on 58/59/60 (one Opus
advisor each, grounded in the relevant lit memo). Verdict on all three:
**DEFENSIBLE AFTER NAMED FIXES.** Briefs saved in `debug/advisor/paper{58,59,60}_advisory.md`.

## 4. The surgery the advisor drove

- **Citation hygiene (58/60/32):** ~14 mature prior-art references added where the
  paper presented standard objects without attribution (principal angles, compound
  matrices, two-projections theory, SO(4)-breaking, the geometric-mean combining rule).
- **D_e gate closed (P58):** 1.961 eV is *correct* — Huang et al. 2010 (15815 cm⁻¹,
  Crossref-verified), re-cited off the anachronistic Huber–Herzberg attribution. (A
  self-caught error: I first entered fabricated author initials, caught it as
  fabrication, pulled the authoritative list from Crossref.)
- **P59 → Paper 61 split** (PI-directed): the number-theory tower (modular/cosmic-Galois
  + Bessel-moment period algebra) moved verbatim to a new companion, Paper 61 (group3),
  aimed at the periods/amplitudes community; P59 keeps the chemistry spine + a pointer.
  Verified a *move, not a copy*; both compile clean, zero undefined refs.
- **P60 coherence → P58 relocate** (PI-directed): the M2 signed-coherence-collapse
  result + its test + registry keys moved to P58's decompactification arc (its natural
  home); P60 points to it. No loss.
- **P60 restructure:** abstract now leads with the atomic metric-free win; the Koridon
  clause's axis-mismatch (my own earlier edit) fixed; Gaussian ratio softened;
  "molecular payoff" → "equivalent-center payoff".
- **NaH R_eq migration (P58):** the "owed" well-minimum is now test-backed. The finding
  that made it worth doing: 3.736 is the driver's *reproducible* coarse-grid (0.25 a₀)
  parabolic value — NOT a corruption (I ran the canonical driver; it emits 3.736); a
  finer 0.02 grid gives the true minimum 3.72, a negligible grid artifact. PI call:
  pin 3.736, note 3.72. Migration test fire-tested (planted 3.900 → FAIL).
- **LiH 30-digit:** found *already* test-backed — never genuinely owed.

## 5. The `/qa` DELTA (unseeded) on {32, 58, 59, 60, 61}

**VERDICT: DEFECTS → REMEDIATED.** 5 LLM reviewers + 12 deterministic gates.
- **Deterministic layer CLEAN** (all 12 gates; C13/C14/C15/C20 re-run on group2+group3
  because the trunk-scope default would have missed the changed papers — the
  gate-scope-audit rule earned its keep).
- **Findings: 1 LARGE + 6 SMALL, every one introduced by this session's changes:**
  - LARGE **F1** — `group3_foundations_synthesis.tex` credited Paper 61's entire tower
    (Legendre/Γ(2), 66-digit certification, cosmic-Galois) to "Paper 59" and never cited
    Paper 61. **The claim-impact reviewer caught it** — invisible to every gate and
    byte-diff; the split's bookkeeping swept P59/P61 + central registries but not the
    *citer* documents. Textbook proof-of-value for the semantic-diff protocol.
  - SMALL — P58 cross-ref pointing out to P60 for now-in-P58 content; **kim_gordon1972
    → actually F. T. Smith** (Crossref); **prosser_hagstrom1968 wrong title** in P58+P60
    (Crossref); P61 cites owed in group2 synthesis + Paper 56; 6 stale rows in
    `topic_to_paper_lookup.md`.
- **Remediated gate-first**, each fix re-verified (Crossref for the two citation
  defects; cite/ref consistency on all edited files; C19 clean; the new bibitems all
  resolve). Two items left with documented disposition: the `[PANEL-VERIFIED]` tag
  (conservative, not an overclaim) and two `test_paper60_coherence_front.py` refs in a
  dated CHANGELOG entry + a transient track log (historically accurate for their date).

## 6. Honest scope

- **Theorem-grade closed this sprint:** none new. The split *relocated* existing proofs
  (irreducibility stays in P59; the modular/period-algebra results move to P61
  unchanged). No new theorem was proved.
- **Structural / process:** the `/advisor` skill (a durable process artifact); the
  P59→P61 split (organizational — a new paper, its own bibliography, cross-cites); the
  P60→P58 coherence relocation.
- **Numerical observation:** NaH R_eq = 3.736 a₀ is a **coarse-grid** value (fine grid
  3.72, negligible vs the +4–5% minimal-basis error); D_e = 1.071 eV reproduces exactly;
  the migration test is a *reproduction pin*, not a new result. All disclosed in Table II.
- **Corrections caught (not new science, but load-bearing):** my own fabricated Huang
  author initials (→ Crossref); kim_gordon→Smith and the Prosser title (inherited from
  explorer memos, caught by the citation reviewer, Crossref-verified).
- **Named open follow-ons:**
  1. **Paper 61 cert OWED** — its content post-dates the last /qa; it needs its own full
     cert pass (it inherited the tower's owed status; not a regression).
  2. **A clean re-delta** is the formal precondition for the eventual full certifying run
     (this sprint closed with a DEFECTS-then-remediated delta; the remediation is
     verified at the deterministic + citation + consistency level).
  3. **`test_paper59_* → test_paper61_*` file renames** deferred (docstring/matrix
     repointing was done; the filename convention is a separate cosmetic pass).
  4. **The elliptic frontier** (closed form via elliptic polylogarithms; the scoped
     PSLQ-against-Γ₁(N) probe) remains the research direction that reaches a new
     community — now housed in Paper 61.
  5. **Avery / Brown externalization** stays PI-gated (unchanged; do not draft/prep).

## 7. Files

- **New:** `.claude/commands/advisor.md`; `papers/group3_foundations/paper_61_bessel_moment_periods.tex`;
  `debug/advisor/paper{58,59,60}_advisory.md`, `p59_split_bookkeeping.md`;
  `debug/lit_scan/*.md` (5 scans); `debug/nah_pes_reeq_verify.py`;
  `tests/test_paper58_coherence_front.py` (renamed from `test_paper60_*`).
- **Edited:** papers 58, 59, 60, 32, 56; group2 + group3 syntheses; `papers/INDEX.md`;
  `CLAUDE.md §6`; `docs/claim_test_matrix.md`; `docs/topic_to_paper_lookup.md`;
  `tests/test_paper58_nah_ladder.py` (+migration test), `test_paper58_headline_numbers.py`,
  10 `test_paper59_*.py` (docstring repointing).
- **Page/size diffs:** Paper 59 shrank (tower removed, 1675→~1601 lines); Paper 61 new
  (8 pp); P58/P60 grew slightly (relocated paragraph / lead sentence). P58 + P61
  compile-verified two-pass clean (exit 0, 0 errors, 0 undefined refs); all edited files
  pass C19 + cite/ref consistency.
