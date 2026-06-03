# Follow-on register

**Last updated:** 2026-06-03 (post-v3.45.3, post-Paper-41-marcolli-vs-closure)

This file is the durable list of carved-out follow-ons from recent sprints that
the PI has not yet authorized for execution. New items are added at the top of
each section; items get crossed out (or removed) when completed. Each entry
should be self-contained enough that a fresh PM session can pick it up cold,
including the originating sprint reference for context.

---

## A. Substantive research sprints (PI decision required)

### A1. Inhomogeneous Robertson–Walker extension of Mixed-Tate test
**Source:** Sprint Mixed-Tate Test (v3.45.3, 2026-06-03; memo
`debug/sprint_mixed_tate_test_memo.md` §"Open questions" #1).
**Estimated effort:** 2–4 weeks.
**Question:** does the Fathizadeh–Marcolli mixed-Tate classification extend
to the GeoVac $S^3$ sector with a non-trivial scaling factor $a(t)$ (the
actual R–W substrate, not the static-$S^3$ sub-case proved in v3.45.3)?
**Why interesting:** F–M proved the result on continuum R–W. The GeoVac
discrete substrate with time-dependent $a(t)$ would need its own
re-derivation, and the discrete-vs-continuum gap is itself informative.

### A2. $S^5$ Bargmann–Segal mixed-Tate analog
**Source:** Sprint Mixed-Tate Test (v3.45.3) §"Open questions" #2.
**Estimated effort:** mini-sprint (1–2 weeks).
**Question:** does the same pure-Tate refinement hold on $S^5$ (where the
analogous CC expansion has three power-law terms including $R^2$, per
Paper 51 `rem:two_term_uniqueness`)?
**Status:** $S^5$ retains the pure-Tate sub-ring at the volume-normalized
SD level; explicit closed forms not yet written out.

### A3. M3 cyclotomic-mixed-Tate verdict at level 4 over $\mathbb{Z}[i]$
**Source:** Sprint Mixed-Tate Test (v3.45.3) §"Open questions" #3.
**Estimated effort:** parallel mini-sprint (1–2 weeks).
**Question:** verify that the M3 sub-mechanism (vertex-parity Hurwitz,
Catalan $G$, $\beta(s)$) lives in cyclotomic mixed-Tate at level 4 over
$\mathbb{Z}[\zeta_4] = \mathbb{Z}[i]$ against the Brown–Goncharov
cyclotomic mixed-Tate classification.
**Why interesting:** would close the arithmetic classification of the
master Mellin engine at all three sub-mechanisms (M1 trivial, M2 pure-Tate
over $\mathbb{Q}$ via v3.45.3, M3 cyclotomic-MT-level-4 via this sprint).

### A4. Inner-factor mixed-Tate extension (Sprint H1 / Yukawa Dirichlet ring)
**Source:** Sprint Mixed-Tate Test (v3.45.3) §"Open questions" #4.
**Estimated effort:** multi-week.
**Question:** the almost-commutative extension $\mathcal{A}_{GV} \otimes
(\mathbb{C} \oplus \mathbb{H})$ (Paper 32 §VIII.C) introduces Yukawa
coupling content classified in Paper 18 §IV's "inner-factor input data"
tier. Whether these inner-factor SD coefficients remain mixed-Tate over
$\mathbb{Q}$ is open.
**Dependency:** would benefit from completing A3 first (M3 mechanism
classified) to scope cleanly.

### A5. Fathizadeh–Marcolli full paper read for precise relative-motive structure
**Source:** Sprint Mixed-Tate Test (v3.45.3) §"Open questions" #5.
**Estimated effort:** 1–2 days reading + writeup.
**Why:** v3.45.3's verdict relied on the F–M abstract + search-surfaced
quotes. A full F–M read would let us write down the specific affine
complement of quadrics + hyperplanes that the GeoVac $S^3$ SD coefficients
reduce to (essentially trivial — a point — in our case, since there are
no continuous $a^{(n)}(t)$ derivatives).
**Why interesting:** would let us state the F–M inheritance corollary as
an explicit motive-equivalence rather than the current text-level
inheritance statement.

---

## B. Collaboration outreach (PI decision required; user has flagged
   as a separate later comprehensive push)

### B1. Garofalo / Urbach group (Bonn / DESY)
**Source:** Phase 4 long-tail lit survey (v3.45.2, 2026-06-03; memo
`debug/lit_survey_phase4_longtail_memo.md`); now cited in Paper 30
(Garofalo et al. arXiv:2311.15926 + Jakobs et al. arXiv:2503.03397).
**Why:** they have built SU(2) Hamiltonian lattice gauge theory using
$S_3$ symmetric-group partitionings — same Wilson Hamiltonian on the
same group manifold as Paper 30, from a quantum-computing-digitization
motivation. Two communities (math-phys QFT-on-graphs vs lattice-QCD-for-QC)
not aware of each other's work.

### B2. Brian Hall (Notre Dame) on coherent states for $S^3$ / Bargmann–Segal
**Source:** Phase 4 long-tail lit survey (v3.45.2); now cited in Paper 24
(Higgs–Pickrell arXiv:2503.23549; Hall's coherent-states-on-spheres
papers also cited).
**Why:** Hall's group's open question for $d \ge 3$ on SO(d)-invariant HO
on d-spheres is exactly Paper 24's $d=5$ case. Natural conversation partner.

### B3. Stetcu / Sarma-Stevenson nuclear VQE
**Source:** Phase 4 long-tail lit survey (v3.45.2); now cited in Paper 23
(Sarma–Stevenson arXiv:2510.02124, Pillai et al. arXiv:2509.08948).
**Why:** direct head-to-head competitors for Paper 23's deuteron 16-qubit /
592-Pauli headline.

### B4. Hekkelman / van Suijlekom / Latrémolière / Farsi cluster
**Source:** Phase 2 connection-finding lit survey (v3.45.2; memo
`debug/lit_survey_phase2_connections_memo.md`).
**Why:** 4 directly-adjacent papers in 2024 alone (Hekkelman-McDonald
arXiv:2412.00628, van Suijlekom arXiv:2409.02773, Farsi-Latrémolière
arXiv:2404.00240, Farsi-Latrémolière-Packer Adv. Math. 437). GeoVac math.OA
arc (Papers 38–49) is in citation distance.

---

## C. Mechanical / technical follow-ons (PM-actionable, low-substantive)

### C1. Paper 14 figure regeneration
**Source:** v3.45.3 corpus cleanup. Replaced 2 missing-PNG figure includes
with `\fbox` TODO placeholders.
**Action:** run `benchmarks/qubit_encoding/pauli_comparison.py`, capture the
two output PNGs (`pauli_scaling.png`, `eri_density.png`), and place them at
`papers/group4_quantum_computing/paper_14_figures/`. The generator currently
outputs to `benchmarks/qubit_encoding/` so a path fix or manual copy is
required. After placement, replace the `\fbox` placeholders in Paper 14
with the original `\includegraphics` lines.
**Estimated effort:** 30 min (assuming the generator script's dependencies
work in the current environment).
**Blocker risk:** uses `openfermion` and other packages — may need a venv
setup if those aren't currently installed.

### C2. Paper 14 §sec:ft_gaussian BLISS-THC comparison table
**Source:** Phase 3 physics audit (v3.45.2; memo
`debug/lit_survey_phase3_physics_memo.md`). Caesura et al. 2025 cite added
but full BLISS-THC fault-tolerant comparison table not built.
**Action:** populate a comparison row in §sec:ft_gaussian using the
Caesura et al. resource estimates vs GeoVac composed numbers from Paper 20.
**Estimated effort:** 1–2 days (requires reading Caesura et al. carefully
for the right comparison metrics).

---

## D. Closed in v3.45.3 (no longer follow-ons)

- ~~Mixed-Tate period test sprint~~ → POSITIVE-with-pure-Tate-refinement;
  Paper 32 Corollary applied; Paper 18 forward pointer closed; Paper 28 T9
  footnote added; Paper 51 rem:cc_uncanny cross-ref added; CLAUDE.md WH2
  status updated.
- ~~Paper 14 dcolumn shorthand bug~~ → fixed (preamble macro added).
- ~~Paper 14 missing figure compile error~~ → fixed (placeholders;
  see C1 above for actual figure regeneration).
- ~~Paper 30 remark environment undefined~~ → fixed (preamble macro added).
- ~~Paper 41 marcolli_vs2014 orphan bibitem~~ → closed (added "Relation to
  Marcolli–vS gauge networks" paragraph at end of related-work section
  citing all three: Marcolli-vS 2014, Perez-Sanchez 2024, Perez-Sanchez 2025;
  fixes the orphan AND closes the Phase 2 recommendation to add
  Perez-Sanchez follow-ups across Papers 25/30/41 — 25 and 30 already had
  them, 41 didn't).

## E. Closed in v3.45.2 (no longer follow-ons)

- ~~Comprehensive 4-phase literature audit~~ → done; memos in
  `debug/lit_survey_phase{1,2,3,4}_*_memo.md`.
- ~~Citation verification gate~~ → done; 18% hallucination rate caught,
  6 mismatches corrected; table in `debug/citation_verification_table.md`.
- ~~21 verified citations added across 14 papers~~ → done.
- ~~5 framing edits (Papers 35, 36, 48, 50, 51)~~ → done.
- ~~Pre-existing Paper 51 `\nmax` macro bug~~ → fixed.
- ~~Pre-existing Paper 36 `\to` in `\text{}` bug~~ → fixed.

---

## Maintenance notes

- Add a row to the appropriate section when a new follow-on is carved out
  by a sprint that does not authorize execution.
- Move completed rows to section D / E / etc. with a strikethrough and
  brief note rather than deleting (preserves the audit trail).
- Reference the originating sprint memo for each entry so the next PM
  session can find context cold.
- This file is referenced from CLAUDE.md §2 most-recent-sprint entry; keep
  it in sync.
