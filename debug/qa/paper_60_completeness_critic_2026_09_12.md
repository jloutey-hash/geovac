# Paper 60 — FULL `/qa` certification run, COMPLETENESS-CRITIC pass (2026-09-12)

**Read-only.** Nothing was modified. This pass names what NO reviewer examined.
It does not re-adjudicate anything the six dimensions did examine.

**Method.** Structural inventory of the target (16 labelled equations, 67 inline
tiers, 6 tabulars, 46 bibitems, 1 caption, 1 appendix, 7 `geovac/sturmian_*.py`
modules, 20 backing test files), then a criteria walk C1–C23 measuring *surface*
rather than verdict, then a programmatic test of the C16 registry against the
claims it is named as guarding. Every load-bearing finding below was confirmed by
executing the gate's own code, not by eye — `grep -F` misreports in this shell and
two candidate findings were **retracted** when re-checked that way (see §5).

---

## 1. Coverage gaps — most consequential first

### G1. The withdrawn `eq:blowup` mechanism is LIVE in the module that owns it, and its named C16 guard cannot reach it

`geovac/sturmian_l2_encoding.py:20-25`:

```
Paper 60 eq:blowup claims this naive encoding inflates lambda as:
    hydrogenic : lambda ~ Q^1.19
    sturmian   : lambda ~ Q^3.33  (the shared-scale overlap ill-conditions,
                 and Loewdin whitening spreads that ill-conditioning into the
                 transformed integrals -- the motivation for the isoenergetic
                 reformulation, which removes the metric for atoms).
```

C8 headline #1 declares exactly this MATERIAL: *"this line read 'driven by overlap
ill-conditioning' — the characterization §2 of the owner WITHDREW on 2026-09-07…
Asserting ill-conditioning here = MATERIAL; C16 `p60-l2-metric-diverges` guards it."*

It does not guard it. Executed against the file, **0 of the 63 C16 entries fire**;
all four `p60-*` entries return `False`. The entry's `files` list *does* include
`geovac/sturmian_*.py`, so the file is scanned — the **pattern** misses. The
2026-09-11 widening added `overlap[^.\n]{0,60}is ill-conditioned`,
`ill-conditioned[^.\n]{0,30}overlap` and the literal `ill-conditions with basis`;
the docstring uses the bare verb (`overlap ill-conditions`) and the gerund
(`ill-conditioning`), neither of which is an adjective and neither of which is
followed by `with basis`.

The provenance makes this sharper. `debug/qa_p60_fix_illcond_test_name.py` shows the
2026-09-11 remediation found this claim in a *test* docstring worded
`"the overlap ill-conditions with basis size"`, fixed that one locus, and added a
pattern matching **that exact wording**. The pattern was written to the locus, not to
the claim — the locus-by-locus class this record has now caught four times — and the
module docstring one file over was never reached. It is worse than a stale sentence:
it states the withdrawn mechanism as *"the motivation for the isoenergetic
reformulation"*, i.e. as the reason for the paper's central result.

### G2. The group2 synthesis re-claims `eq:sigma_law` as derived — C8#21 MATERIAL

`papers/synthesis/group2_quantum_chemistry_synthesis.tex:737-738`:

> "the conditioning grows polynomially, **with a derived band-limited law**"

C8#21: *"`eq:sigma_law` is Kac-Murdock-Szego [PRIOR ART]. Not derived here… **Re-claiming
the asymptotic = MATERIAL.**"* The strings `Kac`, `Murdock`, `Szeg` appear **0 times**
in the synthesis. The owner paper was corrected on 2026-09-11 (P60 L852-861 carries the
[PRIOR ART] attribution); the citer was not. Classic owner-corrected / citer-stale,
§9's retraction→dependents shape — and no C16 entry exists for it (see G4), so no
pattern could catch it. Reviewer 6 ran "10 withdrawn-reading hunts"; those are keyed to
retired *phrases*, and the citer restated the claim in its own word — "derived".

### G3. The synthesis carries NOTHING of v5.11.0–v5.11.3, and its closing clause is the pre-breach reading

Counts in the group2 synthesis: `precondition` **0**, `Toeplitz` **0**, `breach` **0**,
`amplitude floor` **0**, `direct block-encoding` **0**, `overcomplete` **0**,
`circulant` **0**, `Serra` **0**, `DST` **0**, `chirp` **0**, `null space` **0**.

Its Paper-60 block ends (L741-743):

> "…a symmetry-unique heavy atom reinstates the *growth*, **confining that lever to
> homonuclear-diatomic-like systems**."

v5.11.0 breached exactly that case (band-Toeplitz preconditioner bounds cond for the
diatomic **and** water `A_1` — C8#16), and v5.11.3 took the metric penalty n³→n
(C8#19). C9 is a **GATING** dimension for this target. Reviewer 6 audited 8 loci that
exist; it could not, by construction, report on a paper-level result that has no locus.
**Absence is not compliance.**

### G4. Two more live instances of the C8#21 re-claim, in tracked code — and C16 has ZERO coverage of headlines 15–27

`geovac/sturmian_sigma_law.py:14-27` states the asymptotic as internal work —
*"the band-limited concentration rate … i.e. the conditioning exponent is exactly 2
(asymptotically; finite windows fit lower slopes such as the N^1.85 / N^1.97 reported
in Paper 60…)"* — with `Sprint chronicle: CHANGELOG v4.103.0` and **no KMS attribution**.
`tests/test_paper60_sigma_law.py:3` — *"Backs the **derived** conditioning law"*. Both
files predate the 2026-09-11 re-attribution and were not swept.

Executed against the registry, **no C16 entry fires on any of the six withdrawn or
re-attributed readings** pre-registered today:

| C8 headline | withdrawn/re-tiered reading | C16 guard |
|---|---|---|
| #20 | "overcompleteness is the price of one-centre completeness" (frames reading) | **NO GUARD** |
| #24 | the removability corollary (truncation-side reachable / continuum-side untouchable) | **NO GUARD** |
| #26 | the translation law claimed as novel (Ahmad et al. prior art) | **NO GUARD** |
| #21 | `eq:sigma_law` re-claimed as derived | **NO GUARD** |
| #23 | stationary-phase origin for the `π/4` | **NO GUARD** |
| #27 | the translation identification claimed as ours | **NO GUARD** |

The DoD's own standing caution (L658-660) says C16 entries are *owed* for two of these.
It is six. The named backing for #21/#22/#23 is `tests/test_paper60_kms_attribution.py`
— a test, which pins the mathematics in one file and cannot sweep the corpus for
re-claims. G2 and G4 are what that gap costs.

### G5. `sturmian_secular.py`'s grid note warrants its domain with the two quantities Appendix A proves are blind to the defect

`geovac/sturmian_secular.py:51-54`:

> "Grid note: the radial grid (`R_MAX`, `N_GRID`) is calibrated so the
> cumulative-trapezoid (O(dr^2)) Slater potentials are converged — `(5/8) Q` and
> `E(1s^2) = -2.84766` to 5-6 digits. **Do not change it without re-validating the gate
> numbers.**"

`R_MAX: float = 60.0` (L76). Appendix A of the owner paper, L1587-1589:

> "Equal-$n$ configurations and $1s^{2}$ are unaffected, so single-configuration checks
> — including $E(1s^{2})=-729/256$ and the $(5/8)Q$ unit-charge integral — **are blind to
> it**."

So the module offers as its convergence warrant precisely the pair the paper proves
cannot see the truncation bias that caused two successive failed remediations — and
instructs future editors to re-validate against them. The DoD's top watch-note predicted
this class in words (*"the module's own calibration cannot see it"*); nobody read the
module's own calibration statement against it. Related: **no code anywhere enforces
`eq:boxrule`**; `R_MAX = 60.0` is below `3 n_max²` for every `n_max ≥ 5`, and the tests
that need a converged domain monkey-patch it by hand.

### G6. Three labelled equations are outside every enumeration, one of them the paper's own gatekeeper

| label | line | test | claim-matrix | C21 key | C17 family | named in DoD |
|---|---|---|---|---|---|---|
| `eq:boxrule` | 1571 | **none** | **none** | none | none | **never** |
| `eq:T0_closed` | 500 | **none** | none (claims_register only) | none | none | **never** |
| `eq:pw` | 280 | yes | yes | — | — | never |

`eq:boxrule` is the equation the **branch-defining criterion** rests on — *"a run must
verify the domain satisfies `R_MAX ≥ 3n_max²` (App. A) before accepting any exponent."*
The criterion by which every exponent in the paper is accepted is itself unbacked,
untracked and unnamed. `eq:T0_closed` is flagged in-paper (L508) as *"the one leg of this
section that is derived rather than [measured]"* — a SYMBOLIC claim with no test.

### G7. C21 examines ZERO of the molecular half; three registry keys are cited from no `.tex` at all

29 `\gvq` annotations in the paper, **highest at L670**. By region:

| region | lines | decimal literals | `\gvq` |
|---|---|---|---|
| abstract | 23–131 | 22 | 5 |
| obstruction | 180–271 | 12 | 1 |
| atomic | 323–701 | 116 | 23 |
| **molecular** | 731–1095 | **46** | **0** |
| **resource** | 1096–1415 | **79** | **0** |
| **manyelectron** | 1416–1500 | **8** | **0** |
| **conclusion** | 1501–1558 | **5** | **0** |
| **Appendix A** | 1561–1590 | **3** | **0** |

**141 decimal literals past L670, none annotated** — the whole v5.11.x result set plus
the entire resource-table family. And `p60_weighted_collapse_control`,
`p60_window_richardson_pi2`, `p60_window_rms_richardson_pi` are registered but appear in
**no `.tex` file in the repository** (only `CHANGELOG.md`), so C21's annotation check has
nothing to check for them either. The DELTA's own lesson — an unregistered literal is one
C21 is blind to — holds over half this paper.

### G8. `check_cert_staleness.py` ignores the `extra` scope member for every single-paper target

`debug/qa/check_cert_staleness.py:74-77` (the `paper_(\d+)` branch) globs the paper file
and returns. The `GROUPS` branch (L57-60) adds its synthesis; the `trunk` branch (L62-72)
adds the group3 synthesis with a comment saying omitting it *"was the same scope gap as
C19's, one level up."* **That fix was never propagated to the single-paper branch**, even
though `qa_scopes.py` gives `paper_58/59/60/61` explicit `extra` synthesis members.

Measured: the group2 synthesis has **5 commits since 2026-08-18**, including the one that
added the Paper-60 block (2026-09-06) and the DELTA remediation (2026-09-11). The banner
in `paper_60.done.md` — which the record itself tells you to trust over its own prose by
re-running this script — reports **"1 `.tex` changed"**. The staleness instrument reports
zero change on the GATING dimension, for all four single-paper certs.

### G9. The Acknowledgments is a completeness claim frozen one day before thirteen headlines landed

L1544-1557 credits the Averys externally and then enumerates *"**four** results we claim
as our own"*: `eq:scale_lock` + `eq:W_diagonal`, `eq:no_selection`, the `eq:general_v0`
identification, and the variational bound. C8 headlines **15–20** add six more results the
paper claims (the preconditioner lever, the water `A_1` transfer, the M-centre scoping,
the amplitude floor, the direct block-encoding, the one-direction result, plus C8#27's
*"what survives as ours is the SYMBOL"*), and 2026-09-12 added roughly twenty external
lineages (Kac–Murdock–Szegő, Slater–Koster, Löwdin, DLMF, Böttcher–Widom, Ahmad et al.,
Serra, Batenkov et al., Jaffard, Gröchenig–Leinert, Barthelmé–Usevich, Driscoll–Fornberg,
Jordan, Björck–Golub, Hartman–Wintner, Eijkhout–Vassilevski, Rokob et al.,
Monkhorst–Jeziorski, Halmos, Loring, Klappenecker, Böttcher–Spitkovsky) whose provenance
the paragraph does not reflect. The record notes the Acknowledgments was actively
maintained on 2026-09-11 ("the four surrendered theorems reclaimed") — and then not
updated. Reviewer 4 covered "Acknowledgments"; checking four named theorems is not
checking whether four is still the number. Exactly the lesson this DoD wrote on
2026-09-11 about its own debt table: *"an enumeration offered as complete is a stronger
claim than the literals it lists."*

### G10. The transcendental-tagging chain is half-present, and its reference is invisible to every gate

The M2 tagging paragraph (L938-976) is careful and correct as far as it goes, but:

- it names **Paper 18** and **no Paper 34 projection**. CLAUDE.md §6: *"tag every
  transcendental against **both**"*; §4: *"anonymous transcendentals are not allowed in
  production code or papers."* `Paper 34` / `paper34` appear **0 times** in the paper.
- `Paper~18` at **L939** and **L974** are **bare prose references with no `\cite` and no
  bibitem**. C11 checks bibitem titles and years; a bare `Paper~NN` is invisible to it.
  Paper 18 is not in the `paper_60` deterministic scope. So the paper's one dependency on
  a trunk-tier document is reachable by no gate, and no reviewer beat covered
  transcendental tagging — the DoD lists it as W10 and assigns it to no dimension. It is
  a standing memory rule (`feedback_tag_transcendentals`) and a §4 prime-directive item.

*(The `rests on:` discipline was applied in `docs/claim_test_matrix.md:552`, which records
the Paper-18 §"Compactness as the source of discreteness" dependency, and that section
does still exist. Credit where due — the gap is in the paper and the gates, not the matrix.)*

### G11. `tab:resource` was re-derived and enumerated, but never asked whether it is still the paper's resource position

`tab:resource` (L1131-1151, the paper's only float and only caption) prices four metrics
by `κ` / `d_inv` / qubits, where `d_inv` **is** the metric penalty. `eq:amplitude_floor`
(L1304) and `eq:ratio_symbol` (L1342) sit in the *same section*, later, and C8#19 states
the metric penalty now scales as `n` rather than `n^3`. The table carries no
preconditioned row and no caveat in its caption. Reviewer 1 re-derived its values;
reviewer 3 enumerated every row; reviewer 4 covered it. All three are value-level. Nobody
asked the claim-level question.

### G12. "All 46 bibitems verified at source" is unreconciled with the corpus's own standing UNVERIFIABLE caps

Three separate records declare cited sources unread or unverifiable:
`monkhorst_jeziorski1979` (C8#27: *"that paper's two-page body is UNREAD (closed, no
repository copy)"*), `bottcher_spitkovsky2010` (`docs/qa/c23_run_001_paper_60.md`:
*"remains cited-but-unread and flagged"*), and the three Avery book/thesis sources
(2026-08-18 cert: *"UNVERIFIABLE … inaccessible primary sources"*). Either those caps were
lifted during this run or the dimension's summary overstates. Nobody reconciled it, and
the caps are load-bearing prose in the paper (C8#27: *"Dropping that cap = MATERIAL-SMALL"*).

### G13. C23 runs #2 and #3 exist only in the prunable tree

`docs/qa/` holds `c23_run_001_paper_60.md`. Run #2 is recorded only in
`debug/sprint_contraction_seam_memo.md` + CHANGELOG; run #3 only in
`debug/lit_scan/c23_run_003_m2_tagging_memo.md` + CHANGELOG. The DoD cites "C23 run #3"
as the provenance for C8#24's prior-art credit. Also open and tracked nowhere:
`debug/sprint_contraction_seam_memo.md:233` — *"C23 run #3's T3 ABSENT verdict should be
re-run against it."*

### Regions checked and found genuinely covered (recorded so the next run need not re-walk)

Every `\section` is inside one of the two claims beats; the 67 tier labels partition
cleanly across them; the conclusion carries none; there are **no figures and no footnotes**;
the five bare tabulars beyond `tab:resource` are inside `sec:resource` and `sec:molecular`.
`\date{August 18, 2026}` is stale against a 2026-09-12 last-commit, but so is every
sibling (P58 Aug 10 / P59 Aug 16 / P61 Aug 16), and §6's correct-in-place directive makes
git the version record — **convention, not defect.**

---

## 2. Unexercised criteria — UNMEASURED, not passed

| criterion | surface in this target | verdict |
|---|---|---|
| **C5** (K-label half of the hard prohibitions) | **0** occurrences of K = π(B+F−Δ) | **UNMEASURED.** The fitted-parameter half does have surface. |
| **C6** (discrete-vs-continuum precision) | the word "graph" appears **once** in 1815 lines; the paper makes no discrete-graph claim | **UNMEASURED — no surface.** |
| **C7** (trunk-dependent status) | no trunk paper is cited by bibitem; the sole trunk dependency (Paper 18's M2 taxonomy) is bare prose outside the scope | **UNMEASURED**, and the one real dependency is outside the gate (G10). |
| **C12** (K-label cleanliness, deterministic) | ran, examined nothing | **UNMEASURED.** |
| **C16** (retracted-claims) | exercised on the pre-2026-09-12 surface only | **Headlines 15–27: UNMEASURED.** 0/63 entries fire on any of six withdrawn readings (G4). |
| **C17** (headline-number registry) | 2 families, both atomic/molecular exponents | **Partly unmeasured.** The DoD's own note (L292-295) declares SEVEN owed on freeze; five were never created, and thirteen new headlines have none. |
| **C19**, **C20** | ran clean, but named **0 times** in the DoD | **Unassigned** — no dimension owns them. |
| **C21** (numeric consistency) | atomic half only; **0** annotations past L670 | **Molecular half UNMEASURED** (141 literals, G7). |
| **C22** (test-claim backing) | named **0 times** in the DoD; `check_test_claim_backing.py --gate paper_60` **exits 2 — the script takes no `--gate`** | **NOT EXERCISED.** Structurally cannot see G4: its check C keys off the C16 registry (empty for 15–27) and it inspects claim-matrix rows, not test docstrings. |
| **C23** (inverse citation) | run #1 documented; named in the DoD only as provenance inside two C8 headlines | **Unassigned as a dimension**; runs #2/#3 unrecorded in `docs/qa` (G13). |
| **Transcendental tagging** (§4 + `feedback_tag_transcendentals`, DoD W10) | present for Paper 18, **absent for Paper 34** | **UNMEASURED** — no dimension was assigned it. |

---

## 3. Self-referentially stale records

1. **`docs/qa/paper_60.done.md:665`** — the seeding plan still reads *"≥1 catchable by each
   EXERCISED dimension (code / prose / citation; **C9 not exercised**)"*. The same file
   corrects that premise **twice** (scope block L275-285; dimensions block L327-334:
   *"C9 is therefore a GATING dimension"*). The retired premise survived in a third place.
   This is the class the prompt asks for, in this record, again.
2. **`docs/qa/paper_60.done.md:292`** — *"the headline-number registry currently has **NO
   Paper-60 families**"*. False since 2026-08-18; two exist. The seven it says *"must be
   ADDED on freeze"* were never completed and the note carries no supersession marker.
3. **`docs/qa/paper_60.done.md:300`** — the C1 dimension still directs the code reviewer to
   *"RUN `pytest tests/test_paper60_sturmian.py --slow` (17 tests)"*. There are now **13**
   `test_paper60_*.py` files and 7 `test_sturmian*.py`. A reviewer obeying the DoD
   literally runs 1 of 20. (Both code reviewers went well beyond it — the DoD, not the run,
   is what is broken.)
4. **`docs/qa/paper_60.done.md:320-326`** — the C4 dimension enumerates 15 cites; the paper
   has **46 bibitems**. Thirty-one external sources sit outside the declared citation scope.
5. **`docs/qa/paper_60.done.md:423-428`** — W2, the #2-ranked watch-note, navigates to
   *"Abstract (line 62)"* and *"body (line 409)"*. L62 is now "only the one-dimensional
   radial content need be shipped in"; L409 is the scale-lock posing. Dead navigation on a
   HIGH watch-note.
6. **`docs/qa/paper_60.done.md:389-413`** — the declared-debt table is correct for what it
   enumerates and frozen at 2026-09-11, hence silent on the entire v5.11.x literal surface.
   Its own lesson now applies to it a second time.
7. **`debug/qa/qa_scopes.py:134-135`** — *"Papers 59 and 60 have no synthesis footprint
   (their DoDs put C9 out of scope), so their scope is the paper alone"* stands immediately
   **above** the comment block that corrects it. §13.11 rule 9 (replace, never append) —
   the same shape as the `CLAUDE.md:119` LARGE this record caught on 2026-09-11, in the
   file that defines this target's scope.
8. **`geovac/sturmian_secular.py:51-54`** — the grid note's convergence warrant, invalidated
   by the owner paper's own Appendix A (G5).

---

## 4. What a reader should NOT conclude from this run's clean dimensions

The six dimensions established, credibly, that **the numbers reproduce** and that **the
prose that exists matches its declared tier**. Two independent code passes re-derived
essentially every measured quantity and fire-tested 13 planted defects; two claims passes
enumerated every section, table row and tier label; citations and synthesis loci were
walked. That is real and it is the hard part.

It is not a certification that:

- **the paper's newest results are guarded.** C16 has zero coverage of headlines 15–27 and
  C17 has no family for any of them; three of those headlines record a *withdrawn* reading,
  which this corpus names its most frequent defect class. A clean C16 run here means the
  blocklist found no instance of the claims it knows about — it knows about four, from
  before today.
- **the synthesis reflects the paper.** C9 is GATING and it passed on the loci that exist,
  while the entire v5.11.0–v5.11.3 arc has no locus at all and the block's closing clause is
  the pre-breach reading. A dimension that audits present text cannot report absent text.
- **the code says what the paper says.** Tests passing is orthogonal to docstrings being
  true. Three live re-claims (G1, G4) sit in module and test docstrings that every test run
  imports and no assertion touches — including one that presents a withdrawn mechanism as
  the motivation for the paper's central result.
- **the staleness banner is measuring the scope.** It measures one file of two (G8).
- **the DoD is a usable map.** Five of its navigational statements are stale (§3), one of
  them re-asserting a premise the same file corrects twice.

The honest reading of a clean run against this DoD: **the mathematics is sound and the
instruments are behind the mathematics.** That is the same asymmetry §9 records —
*"every substantive mathematical claim of that cycle survived independent re-derivation,
while nearly every guard written to protect one did not."*

---

## 5. Retracted during this pass (recorded so they are not re-raised)

- **"Paper 60 L1450 claims content Paper 58 does not contain."** `grep -niF` returned empty
  for `HeH` and `many-electron` in Paper 58; a ripgrep/Python re-check found both (P58
  L716, 808, 810, 968, 972, 979, 990, 991, 1005, 1006). The cross-reference, including the
  named section *"The continuous side: the decompactification front"* (P58 L792), is
  **supported**. `grep -F` is unreliable in this shell; every finding above was re-verified
  in Python.
- **"The water exponent is a mismatched 1.96/1.97 pair"** (the shape of the DELTA's
  `4.43`/`4.40` defect). The paper states both explicitly and correctly at L1286-1287:
  *"The raw column reproduces the $N^{1.97}$ of the probe above independently
  ($N^{1.96}$ here)."* **Not a defect.**
- **"`\date{August 18, 2026}` is stale."** Corpus convention; all four sibling papers match.
  **Not a defect.**
