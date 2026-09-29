# Retired Papers Register

Every `.tex` in `papers/archive/` has a row here. The register is the thing the
**C24** gate reads; a file in the archive with no row, or a row missing its
required fields, fails the gate.

**Why this file exists.** Archiving a paper removes it from the working set and
from the Zenodo manifest. Its *content* survives (the Zenodo deposit holds the
PDF bytes, and the `.tex` stays in `papers/archive/`), but its *approaches*
become invisible to anyone starting new work. That is the same failure mode
CLAUDE.md §3 exists to prevent one level down, and the same answer applies:
record the attempt with trigger terms, so the existing "check before you start"
discipline picks it up without anyone having to remember a new place to look.

---

## Retirement classes

| Class | Meaning | Re-derivation risk | Triggers required? |
|:------|:--------|:-------------------|:-------------------|
| **SUCCESSOR-COVERED** | A live document occupies the same topic and supersedes this one. | None — the successor is what a search finds. | No |
| **CLOSED** | No successor. The work stopped, was retracted, or was overtaken. | **Real** — nothing live occupies the topic. | **Yes** |
| **CLOSED-VALID** | No successor, and the content is *not* refuted — it was orphaned, not wrong. | Real, and the right response is **reuse**, not re-derivation. | **Yes** |

The distinction matters. A superseded draft needs no trigger probe because its
successor already answers the search. A *closed* paper leaves a hole in the
topic map, and that hole is exactly where a future sprint re-derives something.

---

## Register

| File | Class | Retired | Reason | Trigger terms (`;`-separated) |
|:-----|:------|:--------|:-------|:------------------------------|
| `paper_3_holography.tex` | CLOSED | pre-2026-05-22 | Holographic machinery partially retracted. | holographic duality; AdS/CFT; holographic hydrogen; bulk-boundary correspondence; paraboloid lattice |
| `Paper_4_Universality.tex` | CLOSED | pre-2026-05-22 | Proton-radius claims overtaken by measurement; the surviving kernel was absorbed into Papers 34/35. | holographic central charge; central charge c = 0.045; lepton mass independence; universal holographic behavior |
| `Paper_5_Geometric_Vacuum.tex` | CLOSED | pre-2026-05-22 | Early synthesis, superseded in substance but not by a single named document. | information impedance; emergent spacetime; metric tensor from graph topology |
| `Paper_6_Quantum_Dynamics.tex` | CLOSED-VALID | pre-2026-05-22 | Valid tool paper; citation-orphan in the corpus. **Not refuted — reuse it rather than rebuild it.** | real-time quantum dynamics; Crank-Nicolson propagator; time propagation; broadband molecular spectroscopy; O(V) dynamics |
| `paper_10_nuclear_lattice.tex` | CLOSED | pre-2026-05-22 | Early nuclear draft. Note: the *nuclear* program was later rebuilt independently as Paper 23; this paper's molecular-vibration content was not. | molecular vibration lattice; rotational spectra from graph; nuclear degrees of freedom on S3 |
| `paper_18_exchange_constants_v1.tex` | SUCCESSOR-COVERED | pre-2026-05-22 | Superseded by the current Paper 18. | — |
| `paper_21_geometric_vacuum_synthesis.tex` | SUCCESSOR-COVERED | pre-2026-05-22 | Superseded by the two group syntheses. | — |
| `paper_43_lorentzian_extension.tex` | SUCCESSOR-COVERED | 2026-09-28 | Both headlines lost in the 2026-09-28 /qa trunk cycle (PI decision): the finite four-witness period closure is the automatic integer-spectrum closure of the posited K = 2J_z (identifies no Bisognano--Wichmann flow), and the Pythagorean 1/pi^2 is the posited beta = 2pi/kappa_g squared, not an M1 output; its order-zero/one remark wrongly called growing failures vanishing artifacts. **Successor: Paper 42 Sec. 11** (`sec:lorentzian_closure`), which carries the surviving content with honest tiers: Krein carrier + Lorentzian Dirac (11.1), (4,6) axiom audit + BBB universal-axiom failure (11.2), Lorentzian wedge (11.3), Krein-carrier period closure [SYMBOLIC] + six-witness by construction (11.4), HS orthogonality as a chirality-parity identity (11.7), M3 convention finding (11.8). | Lorentzian four-witness; Krein period closure; Krein-carrier period closure; Pythagorean orthogonality; Pythagorean HS-orthogonality; H_local signature-independence; BBB universal axiom failure |
| `paper_46_strong_form_lorentzian_propinquity.tex` | CLOSED-VALID | 2026-09-14 | Strong-form construction descoped: the full-Krein operator-norm Lipschitz seminorm is degenerate, so the closed-form bound evaluates a rate formula rather than a distance. **Survives: Lemma 3.2, the degeneracy diagnosis** — reuse it, do not re-derive it. | strong-form Lorentzian propinquity; full-Krein operator-norm Lipschitz seminorm; strong-form Lipschitz seminorm |
| `paper_47_two_rate_hybrid_convergence.tex` | CLOSED-VALID | 2026-09-14 | Propinquity arrow descoped with the Paper 45 substrate. **Survives: the norm-resolvent arrow and the three-carrier identification.** | two-rate hybrid convergence; de-compactification gap; compact-to-non-compact resolvent comparison |
| `paper_48_krein_ms_bridge.tex` | CLOSED-VALID | 2026-09-14 | Metric-level theorems open pending a repair whose own verdict is that it weakens the claim to convention. **Survives: the bridge's categorical design** (the contravariant functor and the F2-mismatch resolution via thermal time). | Krein-MS bridge; pinned proper quantum metric space; covered Lorentzian pre-length space; F2 mismatch |
| `paper_49_oslpls_strong_form_bridge.tex` | CLOSED-VALID | 2026-09-14 | Λ inheritance descoped. **Survives: the cocycle-deficit / TICI algebra.** | OSLPLS; operator-system Lorentzian pre-length space; cocycle deficit; TICI algebra |

---

## Note: the Lorentzian tail (Papers 46--49, archived 2026-09-14)

These four were archived together as one decision. The model they rest on --
a Lorentzian quantum metric via K⁺ compression of the Krein Dirac -- was
withdrawn by the Paper 45 annihilation theorem, and the identified repair path
(Toeplitz temporal compressions) carries its own recorded verdict that it
weakens the claim to convention. The register lists them as CLOSED-VALID
rather than CLOSED because **each has surviving content that is not refuted**;
the trigger terms exist so that a future sprint reuses those parts instead of
rebuilding the whole arc to reach them.

**Papers 43, 45 and 50 were deliberately NOT archived** at that time. Paper 45 is descoped
but load-bearing, with five healthy dependents including Paper 38, the WH1
keystone. Papers 43 and 50 were ACTIVE with seven and three external citers.
(Paper 43 was archived later, 2026-09-28, as SUCCESSOR-COVERED by Paper 42 Sec. 11,
after both of its headlines were withdrawn; see its row above.)
**Papers 52 and 53 were also not archived**: they are DRAFT with no dependents,
which is unfinished rather than supplanted, and §9 says explicitly that zero
dependents is not grounds.

---

## Adding a row

A retirement is a **PI decision**. C24 nominates candidates; it never decides,
and nothing is moved without direction. When a paper is retired:

1. Move the `.tex` to `papers/archive/`. **Archive, never delete** — the repo
   copy is what future work greps. Recovering a pruned artifact has already
   cost this project real time (`memory/feedback_resurrect_pruned_artifacts.md`).
2. Add a row here with its class, reason, and — for `CLOSED` and
   `CLOSED-VALID` — trigger terms specific enough not to match the whole corpus.
3. Keep its `papers/INDEX.md` row as a tombstone so the paper stays findable by
   name.
4. Stamp every document that cited it, the same way a retracted *claim* stamps
   its dependents (§9 retraction → dependents rule).

**Trigger terms are the load-bearing field.** Write the phrase a future sprint
would actually use when proposing the same thing, not the paper's title. If a
term matches broadly across the live corpus, it is too generic to be useful —
C24 reports its hit count so over-broad terms are visible.
