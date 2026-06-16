# group3 QA follow-up: angular-ERI-density naming/value + Coulomb/HO layer count

**Date:** 2026-06-15
**Scope:** ANALYSIS ONLY — no paper/doc/code edited. Recommendation spec for the
two open group3 integrity items: (A) the angular-ERI-density naming/value mismatch
across the corpus, and (B) the Coulomb/HO "layer count" inconsistency in Paper 31.

---

## PART A — Angular ERI density

### A.0 The two densities (recap, both correct under their own rule)

| convention | rule | l_max=1 | 2 | 3 | 4 | 5 |
|:--|:--|:--|:--|:--|:--|:--|
| **Global-M_L** (Coulomb 1/r₁₂) | m_a+m_b = m_c+m_d | 14.84% | 8.52% | 6.06% | 4.83% | 3.99% |
| **Pair-diagonal** D_pd (axial / q=0 per pair) | m_a=m_c AND m_b=m_d | 7.81% | 2.76% | 1.44% | 0.90% | 0.62% |

Both are arithmetically correct (reproduced exactly below). They count DIFFERENT
selection rules. The integrity problem is **naming**: the bare symbol `D(l_max)`
is used for the global numbers in Paper 22's Theorem-3 table, and for the
pair-diagonal numbers everywhere downstream (Paper 31, synthesis Table I,
CLAUDE.md §1, claims register).

### A.1 KEY NEW FINDING — the PRODUCTION qubit pipeline realizes PAIR-DIAGONAL

This is the load-bearing fact the earlier review (`debug/review_paper22.md`) did
not establish; it only traced the *diagnostic* `angular_zero_count`, not the
*value* path that builds the actual Hamiltonians.

**Trace:** every production molecule's ERI tensor is built by
`geovac/composed_qubit.py::_build_eri_block` (called for LiH/BeH₂/H₂O at lines
740, 1076, 1088, 1100, 2065, 2075, 2089, 2657, 2667, 2682). `_build_eri_block`
populates its angular table exclusively from `_ck_coefficient` (line 360).

`_ck_coefficient` (line 184) uses `q = mc - ma` and evaluates
`(la k lc; -ma, q, mc)`. The 3j bottom row then sums to
`-ma + (mc-ma) + mc = 2(mc-ma)`, which is nonzero unless **ma = mc**; `_wigner3j`
returns 0 for any nonzero m-sum (line 147). **Result: `_ck_coefficient` zeroes
every m_a ≠ m_c coupling → the realized ERI tensor is pair-diagonal.**

Empirically confirmed:
- `_ck_coefficient(1,1, 1,-1, 2)` → **0.0** (production drops it)
- correct global Gaunt 3j(1 2 1; -1 2 -1) = **0.447** (genuinely nonzero)
- pure-angular density realized by the production `_ck_coefficient` value path:
  **7.81 / 2.76 / 1.44 / 0.90 / 0.62%** at l_max=1..5 — i.e. PAIR-DIAGONAL.

So the **51×–1712× Pauli advantage and the O(Q^2.5) scaling are measured on
pair-diagonal Hamiltonians.** 1.44% is therefore the *honest realized GeoVac
number*; 6.06% is the *generic full-Coulomb-basis bound*.

There is an internal inconsistency WITHIN composed_qubit.py worth noting (not
load-bearing for the headline): the *counting* helper
`estimate_cross_center_eri_count` (line 1503+) uses correct global logic
(explicit `|q|≤k` + global m-conservation, NOT the 3j m-sum check), so it would
count the GLOBAL density, while the *value* path drops to pair-diagonal. The
values win — they define the Hamiltonian that gets JW-transformed and counted.

`casimir_ci.py::_gaunt_ck` uses `q = m1 - m2` (correct global: bottom row sums to
0 identically) — it is the global-convention sibling, used by the atomic FCI path,
not the composed qubit path.

### A.2 Is pair-diagonal physically appropriate for the GeoVac composed basis?

**Basis fact:** `_enumerate_states` (line 119) builds COMPLEX spherical harmonics,
m = −l..+l (line 134). It is NOT a real/axial or m-decoupled basis. So
pair-diagonal is NOT auto-justified by the basis construction.

**Physics:** under the true Coulomb 1/r₁₂ rule a global-allowed, pair-off-diagonal
element such as ⟨p₊₁ p₋₁ | p₋₁ p₊₁⟩ is genuinely nonzero (angular factor 0.24,
computed). It is the "angular-correlation-transfer" channel where the two
electrons exchange units of m about the axis. The production code **omits these
elements**. So the production ERI tensor is a *strict physical subset* of the full
Coulomb tensor (q=0 per-pair sub-block).

**Verdict on A.2:** 1.44% is NOT the generic-Coulomb-basis number; it is the
density of the **q=0 pair-diagonal sub-block that GeoVac actually builds**. It is
honest as "the density GeoVac realizes," but it is *over-optimistic* if quoted as
"the angular ERI density of the Coulomb interaction" (the generic figure is
6.06%, ~4.2× denser). The QA seed-key control **C-A** already encodes the
accepted truth: *"6.06% (Coulomb rule) / 1.44% (pair-diagonal) … CORRECT
(MEASURED)"* — i.e. both are correct under their own label, and neither may be
flagged material. The only defect is the **missing convention tag** downstream.

**Open physics question (flag to PI, not resolvable here):** whether the q=0
pair-diagonal restriction is a *deliberate modeling choice* for the composed
fiber-bundle (each block treated as axially decoupled about its own center) or a
*latent convention bug* in `_ck_coefficient`. If deliberate, 1.44% is the right
realized headline and the corpus should say so explicitly. If a bug, fixing
`_ck_coefficient` to `q = ma - mc` would (a) densify every production Hamiltonian
toward 6.06% and (b) raise all Pauli counts / re-price the 51×–1712× advantage.
**This is a PI decision with code+benchmark consequences; do NOT fix silently.**

### A.3 Every place the density appears

| # | Location | Current text | Convention quoted | Issue |
|:--|:--|:--|:--|:--|
| 1 | Paper 22 abstract (~l.23–32) | lists 14.84/8.52/6.06/4.83/3.99 as `D(l_max)` "Coulomb selection rule"; then pair-diag as "stricter" 2.76/1.44 | global = D; pd = D_pd | **OK / best-in-corpus** — both named, conventions explicit |
| 2 | Paper 22 Thm 3 table (~l.254–266) | columns `D(l_max)` = global, `D_pd(l_max)` = pair-diag | both, tagged | **OK** — the authoritative definition |
| 3 | Paper 22 Thm-3 caption (~l.268–279) | "D = physical Coulomb (M_L); D_pd = stricter pair-diagonal … only for axially-symmetric / m-decoupled bases" | both, tagged | **OK** |
| 4 | Paper 22 §IV Table II + abstract "97.24% vanish" | pair-diagonal magnitude (2.76% density / 97.24% zeros) presented as the Coulomb verification | pd, **untagged** in abstract | abstract "97.24%" is the pd figure; global is 91.48% — tag needed |
| 5 | Paper 22 §V spinor tables | both conventions, `d_sc^FG` reproduces global D bit-exactly | both, tagged | **OK** |
| 6 | **Paper 31 Eq. eq:eri_density (l.314–324)** | `D(l_max)` = **7.81/2.76/1.44/0.90/0.62** (pair-diagonal) | **pd under the name "D(l_max)"** | **MISMATCH** — same symbol as Paper 22's global column |
| 7 | **synthesis Table I `tab:eri_density` (l.500–520)** | `D(l_max)` = 7.81/2.76/**1.44**/0.90/0.62; caption "Universal across Coulomb, HO, WS, SQW, Yukawa" | **pd, untagged + caption conflates** | **MISMATCH** — pd numbers, named D, caption asserts physical-Coulomb universality |
| 8 | **CLAUDE.md §1.6 Phase 4 (l.43)** | "ERI density … verified at l_max=3 as **1.44%**" | **pd, untagged** | **MISMATCH** — bare 1.44% |
| 9 | CLAUDE.md §2 key-results (l.205) | "ERI density depends only on l_max, not V(r). Universal across potentials." | no number | OK (no value) |
| 10 | **docs/claims_register.md row 3** | "1.44% at l_max = 3" | **pd, untagged** | **MISMATCH** — bare 1.44% |
| 11 | Paper 31 §III prose (l.35, abstract) | "ERI density 1.44% at l_max=3" as canonical universal-sector content | **pd, untagged** | **MISMATCH** |

### A.4 RECOMMENDED single naming+value scheme

Adopt ONE convention-tagged vocabulary, corpus-wide:

- `D(l_max)` ≡ **global-M_L Coulomb density** = 14.84/8.52/**6.06**/4.83/3.99%.
  This matches Paper 22's authoritative Theorem-3 column and the physically
  correct 1/r₁₂ selection rule.
- `D_pd(l_max)` ≡ **pair-diagonal (q=0) density** = 7.81/2.76/**1.44**/0.90/0.62%.
  Tag it explicitly as "the q=0 per-pair sub-block, = the density the GeoVac
  composed pipeline realizes."

**Headline rule:** when a document needs ONE number for the GeoVac realized
Hamiltonian (sparsity/Pauli context), quote **D_pd = 1.44%** AND name it
"pair-diagonal (realized)". When it needs the generic-Coulomb-basis bound, quote
**D = 6.06%** and name it "global-M_L (Coulomb)". Never quote either under the
bare unqualified phrase "the angular ERI density."

This keeps both numbers (both are real and both are useful — 1.44% is the actual
GeoVac advantage, 6.06% is the honest generic bound) while making `D` vs `D_pd`
mean exactly one thing each everywhere.

### A.5 Per-document spec (NO edits applied)

- **Paper 22** — already the cleanest. Two mechanical fixes only:
  (i) abstract "97.24% vanish" → tag as the pair-diagonal (q=0) figure and add
  "91.48% under global M_L"; (ii) §III Definition prose currently describes the
  GLOBAL rule but the table column D_pd is computed under pair-diagonal — add one
  clause "(and, for D_pd, m_a=m_c, m_b=m_d)". **PM-mechanical.**
- **Paper 31 Eq. eq:eri_density (l.314–324)** — the numbers are pair-diagonal but
  the symbol is `D(l_max)`. Either (a) relabel the symbol to `D_pd(l_max)` and add
  "(pair-diagonal; the density the composed pipeline realizes; global-M_L Coulomb
  is ~4× denser — Paper 22 Table)", or (b) swap to the global 6.06 numbers under
  `D`. **Recommend (a)** — it is the realized GeoVac figure and the paper's whole
  point is the realized sparsity advantage. **PM-mechanical relabel; (a) vs (b) is
  a small framing call — default (a).**
- **synthesis Table I (l.500–520)** — relabel column `D(l_max)` → `D_pd(l_max)`,
  fix the caption: drop "fraction … for which the Gaunt coefficient does not
  vanish" (that's the global rule) and replace with "pair-diagonal (q=0) density,
  = the composed-pipeline-realized figure; global-M_L Coulomb density is given in
  Paper 22 Table" — and KEEP the "universal across 5 potentials" clause (that part
  is true for both conventions). **PM-mechanical.**
- **CLAUDE.md §1.6 (l.43)** — "verified at l_max=3 as 1.44%" →
  "verified at l_max=3 as 1.44% (pair-diagonal / realized; 6.06% under the global
  Coulomb M_L rule)". §1.6 is a **HARD-NO PM-edit zone** (Project Phase) →
  **PI-decision / PI must apply.**
- **claims_register.md row 3** — "1.44% at l_max = 3" →
  "1.44% (pair-diagonal, realized) / 6.06% (global Coulomb) at l_max=3".
  **PM-mechanical** (register is PM-editable).
- **Paper 31 §III/abstract prose "1.44%"** — same tag as Eq. **PM-mechanical.**

**The one PI-decision item in Part A:** whether `_ck_coefficient`'s pair-diagonal
restriction is intended physics (keep 1.44% as realized headline, tag it) or a bug
to fix (densify to 6.06%, re-benchmark Pauli counts). The naming fixes above are
correct EITHER way; only the code question needs the PI.

---

## PART B — Coulomb/HO layer count

### B.1 Owning paper (Paper 24) = the current authority

Paper 24 §sec:asymmetry_layer4 (l.616–649) enumerates **SIX layers** explicitly:

1. spectrum-computing role of L₀ (Coulomb-specific)
2. calibration π (Coulomb-specific)
3. Wilson lattice gauge with natural matter coupling (Coulomb-specific)
4. modular-Hamiltonian Pythagorean orthogonality (Coulomb-specific)
5. spectral-action gravity termination — Paper 51 (Coulomb-specific)
6. master Mellin engine chirality-parity selection rule (new, June 2026)

So **the current correct count per the owning paper is SIX.** (Note: Paper 24 has
its own internal debt — a stale footnote at l.437–438 still says "refined to five
layers", and the subsection heading at l.603 reads "Fourth layer" though the list
runs to six. That is Paper 24's own cleanup, separate from Paper 31, but worth a
row in the same fix-pass.)

### B.2 Paper 31 is internally self-contradictory on the count

| location | says | count |
|:--|:--|:--|
| abstract (l.56–58) | "from two layers to **three**" | 3 |
| §five-layers title (l.659) + text (l.664) | "from two layers … to **five** layers" | 5 |
| §five-layers body list (l.670–700) | enumerates Layer 1, 2, 3 only | 3 |
| Table tab:five_layers (l.704–724) | **five** rows above the divide | 5 |
| conclusion (l.1388–1391) | "**five** layers, not two" then lists 3 | 5/3 |

Paper 31 never reaches six; it predates Paper 24's Layer 6 (and arguably Layer 4/5
were added to Paper 24 after Paper 31's "five" table was written — the table rows 4
and 5 are "Modular-Hamiltonian Pythagorean (Paper 24 §V)" and "Gravity termination
(Paper 51 G4-5a)", matching Paper 24 layers 4 and 5).

### B.3 The dangling \cite{paper51}

There is **no `\cite{paper51}`** in Paper 31 (grep: no matches). The reference is
the **plain text** "Gravity termination (Paper~51 G4-5a)" in row 5 of
tab:five_layers (l.717). Paper 31's bibliography has **no `\bibitem{paper51}`**.
So this is not a dangling `\cite` producing a `[?]`; it is an *uncited
plain-text reference* to Paper 51. (Task statement's "dangling \cite{paper51}" is
slightly off — it is an un-bibitem'd textual mention, which is the lighter
problem.) Paper 24's footnote correctly cites Paper 51 by name too.

### B.4 RECOMMENDED fix for Part B (NO edits applied)

**Reconcile Paper 31 to SIX, matching the owning paper (Paper 24).** Concretely:

1. **abstract (l.56–58)** — "from two layers to three" → "from two layers to six"
   and update the parenthetical to list all six (or, lighter: "to a multi-layer
   stratification — six layers as of Paper 24"). **PM-mechanical.**
2. **§five-layers** — retitle "The Five-Layer Sharpening" → "The Six-Layer
   Sharpening"; expand the body list (l.670–700) from 3 to the full 6, importing
   verbatim from Paper 24 §asymmetry_layer4 (single source of truth). Update
   §"Why five and not two" (l.730) accordingly. **PM-mechanical** (pure
   synchronization to the owning paper; no new physics).
3. **tab:five_layers** — rename label/caption to six-layer; add the two missing
   rows (Layer 6 master-Mellin chirality; Layer 1 is already row "L₀ computes
   spectrum"). Re-count: table currently has 5 Coulomb-specific rows = layers
   1,2,3,4,5; add layer 6. **PM-mechanical.**
4. **conclusion (l.1388–1391)** — "five layers, not two" → "six layers, not two"
   and list all six. **PM-mechanical.**
5. **\cite{paper51}** — add a `\bibitem{paper51}` (Paper 51, the spectral-action
   gravity paper) and convert the plain-text "Paper~51 G4-5a" in the table to
   `\cite{paper51}`. **PM-mechanical.** (Alternative: if the PI prefers Paper 31
   stay scoped to the asymmetry as known at its writing, drop rows 4/5/6 and
   revert to a clean THREE-layer statement consistent with the abstract — but that
   contradicts the owning paper's current six and loses information. **Recommend
   sync-up to six, not down to three.**)
6. **Paper 24 self-debt** (same pass) — fix the stale "five layers" footnote
   (l.437–438) → six, and the "Fourth layer" subsection heading (l.603) →
   "Higher layers" or "Layers 4–6". **PM-mechanical** (owning paper internal
   consistency).

**Part B has no PI-decision content** — it is pure synchronization of Paper 31 (and
Paper 24's stale footnote/heading) to Paper 24's authoritative six-layer list, plus
one missing bibitem. All PM-mechanical.

---

## Summary of dispositions

| item | fix | owner |
|:--|:--|:--|
| A — Paper 22 abstract "97.24%" + §III def tag | tag conventions | PM-mechanical |
| A — Paper 31 Eq./prose `D`→`D_pd` (1.44% realized) | relabel + tag | PM-mechanical |
| A — synthesis Table I relabel + caption | relabel + caption | PM-mechanical |
| A — claims_register row 3 dual-number | add tag | PM-mechanical |
| A — CLAUDE.md §1.6 "1.44%" tag | add tag | **PI (Phase section, PM-locked)** |
| A — is `_ck_coefficient` pair-diagonal intended or a bug? | keep+tag vs fix+re-benchmark | **PI-decision (code)** |
| B — Paper 31 three/five → SIX everywhere | sync to Paper 24 | PM-mechanical |
| B — Paper 31 add \bibitem{paper51} + \cite | add bibitem | PM-mechanical |
| B — Paper 24 stale "five" footnote + "Fourth layer" heading | sync internal | PM-mechanical |

**Naming scheme to adopt corpus-wide:** `D(l_max)` = global-M_L Coulomb
(6.06% @ l_max=3); `D_pd(l_max)` = pair-diagonal / realized (1.44% @ l_max=3).
Never quote either unqualified.
