# Paper 59 → Paper 61 structural split: central-bookkeeping owed to the PM

**Date:** 2026-09-06
**What happened:** Paper 59's number-theory tower (the two periods/amplitudes-facing
sections) was split into a new companion, **Paper 61**, per PI direction. The
mechanical edits to the two `.tex` files are done and both compile clean (see
"Verification" below). This file lists the central-registry / index / doc updates
that were deliberately **NOT** done in the split (they are PM / gated-artifact
territory) and that you now owe.

---

## 1. Files changed by the split (for your CHANGELOG / §2 one-liner)

- `papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex` — EDITED
  (two sections removed → pointer paragraph; 8 cross-refs converted; 11 bibitems
  removed; `\bibitem{loutey_paper61}` added and cited). 1675 → 1601 lines; 49 → 39
  bibitems.
- `papers/group3_foundations/paper_61_bessel_moment_periods.tex` — CREATED (new
  companion, Paper 61; 8 pp; 19 bibitems).

**Moved content (verbatim except ref conversions):** the two sections
`\section{The modular structure and its cosmic-Galois periods}` (`sec:modular`)
and `\section{The Bessel-moment period algebra: master family and quadratic
relations}` (`sec:bessel_algebra`). Everything else — Secs 1–6, the F12 section
`sec:f12` (chemistry-facing, STAYS), Scope, Conclusion — stayed in P59.

---

## 2. `papers/INDEX.md` (NOT touched — you own it)

**(a) Add a Group-3 row for Paper 61.** Group 3 table currently ends at Paper 57
(lines ~46–49). Suggested row (match the existing style):

```
| 61 `paper_61_bessel_moment_periods.tex` | ACTIVE | Number-theory tower of the three-center ERI (companion to Paper 59): Legendre/Γ(2) modular family, in-domain CM Γ-value fibres, X₀(2) home + conductor-4 (ℤ[i], G=β(2)) arithmetic, resurgent-Lambert reading; Bessel-moment period algebra — Wronskian W₀D⁻², quadratic relations B=πΩ forcing Sp₄(ℤ); Eisenstein/CM period, not a cusp-form L-value; cosmic-Galois bridge ↔ Paper 56; closed form open |
```

**(b) Revise the Paper 59 row (line 86).** Its current description advertises
content that MOVED: *"period on the Legendre/Γ(2) modular family with CM Γ-value
fibres (cosmic-Galois bridge, ↔ Paper 56)."* Recommend trimming P59's row to the
chemistry/genus-jump/irreducibility core (which stays) and appending *"; the
modular / cosmic-Galois structure and the Bessel-moment period algebra are now in
the companion Paper 61."*

---

## 3. `CLAUDE.md` §6 (NOT touched — PM-editable, you own it)

- **Folder table (§6):** `group3_foundations/` count **11 → 12**. `group2_quantum_chemistry/`
  is unchanged at **11** (Paper 59 stays in group2; the split created a NEW group3 paper).
- **§6 "Load on topic" list:** the periods/foundations line reads
  `foundations/periods → 18, 54, 55, 56, 57`. Consider adding **61** (it is a
  periods/cosmic-Galois paper). Paper 59 remains under `chemistry solvers`.
- No §6 status-flag change is required (neither P59 nor P61 is a keystone/guardrail).

---

## 4. Numeric registry (`debug/qa/numeric_registry.py`) — NOT touched

**Finding: no registry key's owning locus moved.** I searched the registry (and
`check_retracted_terms.py`) for the moved constants — the T2 value
`0.395355765901713964325…` (66 digits), the period pairing `{-π, 0, 2π}` /
`B = πΩ`, the diagonal `A = J(¼,¼,1)/4 = 0.0790540321…`, `K(½)=Γ(¼)²/(4√π)`, the
disc-8 period, the co-area / θ₃² values — and **none of them are registered.**
The registry's keys are chemistry/composed-scaling and Papers 2/7/22/32/38/40
constants; there is no Paper-59 Bessel-moment entry to re-point.

**So there is nothing to re-point, but flag for your judgment:** these load-bearing
moved numbers (esp. the 66-digit **T2** value, backed by
`benchmarks/certified_reference/` + `test_paper59_t2_value.py`, and **B = πΩ**) are a
pre-existing C21 registry gap. They now live in **Paper 61**. If/when they are
registered, the owning locus is Paper 61 (`sec:modular` / `sec:bessel_algebra`),
not Paper 59. Recommend registering at least the T2 value and the `{-π,0,2π}`
pairing when the numeric surface of the periods papers is next swept.

---

## 5. `tests/test_paper59_*` backing the MOVED sections — LIST ONLY (not moved, not modified)

Per instruction these were **left in place and untouched** (filenames keep the
`paper59` stem even though the material is now in Paper 61). Tests whose claims are
now in **Paper 61**:

| test file | backs (moved-section claim) |
|:--|:--|
| `tests/test_paper59_t2_value.py` | T2 66-digit value, K⁻⁷ tail law, digit-19 anchor correction (`sec:modular`) |
| `tests/test_paper59_corner_sigma2.py` | ρ=σ² spectral corner, collinear value, guarded PSLQ-negative (`sec:modular`) |
| `tests/test_paper59_theta_chi4.py` | θ₃² χ₋₄ theta series, Γ(4) level / conductor-4 (`sec:modular`) |
| `tests/test_paper59_bd_jacobian.py` | modular-lambda (Borwein) Jacobian dλ/dτ (`sec:modular`) |
| `tests/test_paper59_coarea_reduction.py` | co-area reduction to X(2)-locus, Φ(ρ)=Φ(1/ρ)/ρ² symmetry (`sec:modular`) |
| `tests/test_paper59_gamma0_2.py` | X₀(2) fold / involution (`sec:modular`) |
| `tests/test_paper59_resurgence_corner.py` | Gevrey/Borel-radius resurgence, corner power counts (`sec:modular`) |
| `tests/test_paper59_resurgent_skeleton.py` | resurgent-skeleton "5-for-5" pattern, algebraic Stokes (`sec:modular`) |
| `tests/test_paper59_bessel_moment_algebra.py` | master-period Wronskian, period pairing −π,0,2π (`sec:bessel_algebra`) |
| `tests/test_routeC_momentum.py` | (partly) intersection-form/planes symbolic + diagonal-A tests cited by `sec:bessel_algebra`/`sec:modular`; also backs KEPT momentum reduction |
| `tests/test_paper59_four_subspace_lambda.py` | four-subspace λ / cross-ratio (`sec:modular`) — on disk; not cited by exact filename in prose but topically backs moved content |

**Shared (backs BOTH kept and moved):**
- `tests/test_paper59_euclidean_dictionary.py` — cited in KEPT P59 (Euclidean
  dictionary, Secs 2/4) **and** in moved `sec:modular` (Euclidean-interval Stokes
  location). Leave in place; backs both papers.

**Decision for you:** whether to rename any of these to `test_paper61_*` and add a
`docs/claim_test_matrix.md` row mapping them to Paper 61. I did **not** rename or
touch any test (out of scope per the split instruction). The `test_routeC` bibitem
STAYS in P59; its annotation was edited so the two clauses that back moved material
now read "backing the companion Paper 61".

---

## 6. Bibliography split (already applied — for your audit)

- **Moved to P61 (cited only by moved sections), 11:** `beilinson_levin1994`,
  `broadhurst_dorigoni_characters2025`, `broedel_eisenstein2018`,
  `brown_cosmic2017`, `brown_levin2011`, `brown_mmv2014`, `fantini_rella2024`,
  `loutey_paper35`, `loutey_paper56`, `manin2006`, `mcspirit_rolen2025`.
- **Duplicated (cited by both; kept in P59 AND copied to P61), 7:** `bbbg2008`,
  `broadhurst2008`, `broadhurst2016`, `broadhurst_dorigoni2026`,
  `fresansabbahyu2023`, `zhou_wick2018`, `zhou_wronskian2018`.
- **New:** `loutey_paper59` added to P61 (self-title of P59); `loutey_paper61`
  added to P59 (title of the companion).

---

## 7. Cross-references converted (enumerated)

**In P59 (kept) — refs to now-moved labels → prose "the companion (Paper 61)":**
1. `sec:obstruction` closing sentence: `Sec.~\ref{sec:modular}` → "developed in the companion (Paper 61 \cite{loutey_paper61})".
2. `sec:literature` (arithmetic obs): `theta series θ₃² of Sec.~\ref{sec:modular}` → "of the companion (Paper 61)".
3. `sec:literature` (Reconciliation): stale `Sec.~4` (value-level searches) → "the companion (Paper 61)".
4. `sec:scope`: `Section~\ref{sec:bessel_algebra} localizes it` → "The companion (Paper 61) localizes it".
5. `sec:scope`: `(Sec.~\ref{sec:bessel_algebra})` (Ω identified in closed form) → "(developed in the companion, Paper 61)".
6. `sec:scope`/Conclusion tail: `(Sec.~\ref{sec:bessel_algebra})` → "(companion, Paper 61)".
7. `test_routeC` bibitem: `backing Sec.~\ref{sec:bessel_algebra})` → "backing the companion Paper 61)".
8. `test_routeC` bibitem: `backing Sec.~\ref{sec:modular})` → "backing the companion Paper 61)".

**In P61 (moved) — refs to labels that stayed in P59 → prose "Paper 59":**
- `Eq.~\eqref{eq:curve}` → "the genus-one curve of Paper 59" (×1)
- `Sec.~\ref{sec:obstruction}` → "Paper 59" (×6; and one parenthetical `(Sec.~\ref{sec:obstruction})` after "formally self-adjoint" was dropped as redundant)
- `connection~\eqref{eq:pf}` → "connection of Paper 59" (×2)
- `operator~\eqref{eq:pf}` → "operator of Paper 59" (×1)
- `Eq.~\eqref{eq:pf}` → "the Picard–Fuchs equation of Paper 59" (×2)
- `fibre~\eqref{eq:besselmoment}` → "fibre of Paper 59" (×1)
- `Eq.~\eqref{eq:K0}` → "the K₀ identity of Paper 59" (×1)
- stale `Sec.~4` (shadow reading, in `sec:bessel_algebra`, points to its sibling
  `sec:modular`) → `Sec.~\ref{sec:modular}` (resolves inside P61).

**Untouched (resolve within their own file):** P61's `\ref{sec:modular}`,
`\ref{sec:bessel_algebra}`, `\eqref{eq:lambda_rho}`, `\eqref{eq:kw}`,
`\eqref{eq:planepi}`; and P61's `\cite{loutey_paper35}`, `\cite{loutey_paper56}`,
`\cite{brown_cosmic2017}` (other GeoVac/external papers, bibitems moved with them).

---

## 8. Verification (already run)

- **pdflatex (MiKTeX), 2 passes each:** both PDFs build. **Zero** undefined
  references, **zero** undefined citations, no fatal errors (only the standard
  benign `nameref` "definition of \label has changed" notice). P61 = 8 pp, P59 = 9 pp.
- **Cite ↔ bibitem:** P59 — 39 bibitems, all 39 cited, no undefined `\cite`, no dup
  keys. P61 — 19 bibitems, all 19 cited, no undefined `\cite`, no dup keys.
- **Ref ↔ label:** no undefined `\ref`/`\eqref` in either file. (P59 has 7
  pre-existing unused labels — `eq:laplace`, `eq:modpf`, `eq:exactform`,
  `eq:inmodule`, `sec:genus`, `sec:literature`, `sec:f12` — these were never
  cross-referenced before the split; not introduced here.)

---

## 9. Items flagged for your review (judgment calls I made)

- **Stale `Sec.~4` references.** P59 contained two hard-coded `Sec.~4` refs that do
  NOT match P59's actual §4 (`sec:reduction`, which has no "shadow"/"value-level"
  content). The referenced content (value-level integer-relation searches, the
  D→0-shadow reading) lives in the moved `sec:modular`. I therefore treated the
  kept-side one (in `sec:literature`) as a cross-boundary ref → "the companion
  (Paper 61)", and the moved-side one (in `sec:bessel_algebra`) as an intra-P61 ref
  → `\ref{sec:modular}`. If you intended `Sec.~4` to mean literal §4, revisit.
- **Companion `\date{August 16, 2026}`** — I matched P59's date (the material is
  P59's), not the split date (today, 2026-09-06). Change if you want the companion
  dated at creation.
- **Companion title** — used the PI's suggested form: "The Bessel-Moment Period
  Algebra of the Three-Center Coulomb Integral: Modular Structure and Cosmic-Galois
  Placement".
- **Pointer placement** — the single pointer paragraph sits at the end of
  `sec:literature` (where the two sections were), before the F12 section, so it
  reads as the closing note of the placement discussion. No new section heading was
  added for it.
