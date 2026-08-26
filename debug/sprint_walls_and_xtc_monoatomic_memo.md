# Sprint: `/walls` pass + atomic-xTC monoatomic sweep (2026-08-24)

Umbrella memo for one session that landed two related deliverables: a new negative-results
refinement tool (`/walls`), and — as its first "ride the recommendation" follow-through — the
completion of the atomic xTC angular-sparsity result across the whole monoatomic library.
Versions v5.0.10 (tool) and v5.0.11 (sweep).

---

## A. `/walls` — the negative-results refinement gate

### A.1 Motivation
`/qa` systematically refines the *positive* corpus (papers/claims) but nothing refines the
*negative* corpus — the walls (CLAUDE.md §3 + `debug/` dead-ends + `memory/` meta-findings), of
which the project has dozens. Yet two of the project's best *positive* structural insights —
`multi_focal_wall_pattern` and `two_kinds_of_sparsity` — were born from clustering negatives until
a shared mechanism fell out. In the skeleton picture the walls are the survey stakes on the
forced/free boundary, so refining them sharpens the skeleton map. `/walls` is the machine that does
that on purpose.

### A.2 What was built
- **`.claude/commands/walls.md`** — the skill. Four moves: (1) RE-VERIFY each wall vs current state
  (STANDING / SOFTENED / BREACHED / MIS-SCOPED / SUPERSEDED; no bare verdicts — every status cites
  the current-state evidence that decides it); (2) CLASSIFY HARD / SOFT / OPEN-LEANING; (3)
  CLUSTER + MINE the shared mechanism (CRYSTALLIZED / FORMING / SINGLETON); (4) PROMOTE earned
  crystallizations, PI-gated. PI-invoked genuine-trigger command, no memory rule (§13.9a).
- **`docs/walls/register.md`** — the living register (refinement layer atop the append-only §3
  ledger). Bootstrapped by working the **CHEM-ACCURACY** cluster fully; seven further clusters named
  and queued.
- **CLAUDE.md §13.9b** — `/walls` row added (7→8 commands; genuine-trigger group 4→5). *PI-directed
  edit of a normally-PM-locked §13* (Josh: "add walls to Claude.md").
- Memory: `walls_pass.md` + MEMORY.md index line. CHANGELOG v5.0.10.

### A.3 First-run finding (the CHEM-ACCURACY crystallization)
Working the cusp/sparsity/composition family produced a single sharper statement — the **two
compounding walls**:
- **Wall A (accuracy — the e-e cusp):** a Layer-2 basis limit; breachable in principle by explicit-r₁₂
  machinery.
- **Wall B (sparsity — non-orthogonality):** GeoVac has only symmetry sparsity, which dies at two
  centers; every attempt to keep it while adding molecular coupling relocates a non-orthogonality cost
  (dual-basis theorem). HARD.
- **Compounding:** they co-locate favorably *only at one center* — for atoms Wall A is cheaply
  breachable (xTC) while Wall B is absent; for molecules breaching Wall A re-triggers Wall B (the
  geminal is ~93% collinear → κ(S)≈200). So chemical accuracy and native sparsity co-locate at one
  center and separate at two. Operational: atoms GO; molecular *accuracy* STOP; molecular
  *sim-structure* GO.

This SOFTENED the "cheap cusp loses sparsity" wall (breached for atoms by xTC) and is flagged as a
**PI-promotion candidate** (proposed home: a sharpening of `two_kinds_of_sparsity` and/or a Paper 14
remark). Not written into a paper — promotion is PI-gated.

---

## B. Atomic xTC — monoatomic library fleshed out (v5.0.11)

### B.1 The move
`/walls` flagged the atomic lane as GO. The tracked angular-support engine
(`geovac/xtc_angular_sparsity.py`) is radial-free and EXACT (integer block counts), and the fill-in
depends only on the reference density's angular multipoles — so "flesh out all the monoatomics" is a
cheap, exact sweep over open-shell angular character. Driver: `debug/xtc_monoatomic_sweep.py`
(pre-registered prediction, class table, validations, per-atom Z=1–56 map); data
`debug/data/xtc_monoatomic_sweep.json`.

### B.2 Result — the Unsöld law
xTC fill-in is **0 at every basis iff the reference density is spherical**, which by Unsöld's theorem
means **closed subshells, half-filled high-spin subshells, and s-open shells**. This strictly
generalizes the v5.0.9 "s-reference → 0" to include the half-filled shells (N 2p³, Cr/Mn 3d⁵, …).

Census across Z=1–56: **29 spherical** (exact 0 fill-in), **16 open-p**, **11 open-d**.

| Basis | spherical | open-p | open-d | Coulomb density |
|:--|--:|--:|--:|--:|
| s+p | 0 | 0 | 0 | 14.8% |
| s+p+d | 0 | 4 | 4 | 8.5% |
| s+p+d+f | 0 | 96 | 96 | 6.1% |
| s+p+d+f+g | 0 | 264 | **268** | 4.8% |

The transition metals (open-d) carry an additional **L′=4** density multipole (open-p has only L′=2);
it makes no difference through s+p+d+f (identical 0/4/96) and first bites at s+p+d+f+g (+4 blocks of
390,625). Fill-in never exceeds ~0.15% of blocks; support density stays ~5–6% throughout.

### B.3 Validations (all exact, passed)
- m-filling independence of support (two C 2p² fillings → identical 4).
- core support-invariance (p² bare vs +[Ne] core → identical).
- p² ≡ p⁴ (both 96); half-filled p³ ≡ d⁵ ≡ 0.
- Unsöld law tested *generally* over 13 (l,k) cases, both directions.

### B.4 Correction surfaced
The prior Paper 14 paragraph listed **N** as an open-p (fill-in) reference; N's 2p³ is half-filled →
spherical → 0 fill-in. Corrected in the paper.

### B.5 What was banked
- `geovac/xtc_angular_sparsity.py`: `open_subshell(l,k)`, `is_spherical(l,k)` promoted (tracked).
- `tests/test_paper14_xtc_angular_sparsity.py`: 13 → **30 fast** (Unsöld law general test, half-filled,
  open-d, g-divergence). Green.
- Paper 14 `sec:tc_atomic_sparsity`: `tab:xtc_fillin` extended (3 classes × 4 bases); spherical
  statement generalized to the Unsöld law + 29/16/11 census; L′=4 divergence + Paper 16 periodicity.
  Compiles clean (0 undefined refs).
- `docs/claim_test_matrix.md` +2 rows; `docs/walls/register.md` CHEM-ACCURACY (A)-lane updated.
- Visual artifact: periodic table Z=1–56 coloured by xTC class (the ½-filled islands visible).

---

## §6. Honest scope

**Closed at exact (theorem-adjacent) grade — the strongest tier here:**
- The Unsöld law (fill-in = 0 ⟺ spherical reference density) and all fill-in counts
  (0/4/96/264 open-p, 0/4/96/268 open-d) are **EXACT** — angular support is radial-independent, so
  these are integer facts, not tolerances. Tested generally over (l,k), both directions.
- The 29/16/11 census is exact given standard ground-state configurations (reference chemistry data).

**Structural (mechanism):** the monopole collapse — a spherical reference density carries only the
L′=0 multipole, forcing the correlator leg to the monopole and collapsing the non-abelian vertex.

**MEASURED (unchanged from v5.0.9, NOT swept this sprint):** the radial 1-norm ratios (0.76–0.91×)
and the 279-Pauli count are grid-engine measurements on a single-common-k Coulomb–Sturmian basis with
one spin-independent geminal (classical non-Hermitian PoC). These are resource numbers, not accuracy.

**Named open follow-ons:**
1. **Per-atom accuracy** — the radial/1-norm/energy side of xTC across the library is NOT swept; this
   is the honest open frontier in the atomic lane (the `/walls` register flags it).
2. **Molecular (2-center)** — the l-sparsity does not survive to two centers; xTC preserves only
   m-conservation sparsity there ([[two_kinds_of_sparsity]]). Atomic-only scope stands.
3. **CHEM-ACCURACY two-walls crystallization** — PI-promotion candidate (register), not yet a paper
   claim.

**Hard-prohibition check (§13.5):** no fitted parameters; no geometry-hierarchy change; no negative
result deleted (§3 untouched; register annotates); K=π(B+F−Δ) not touched. The one normally-locked
edit (§13.9b, part of §13) was made under **explicit PI direction**.
