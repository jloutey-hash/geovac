# Scoping memo: `/walls` architecture-swap audit — composed → prolate-native diatomic

**Status:** EXECUTED + PI-adjudicated 2026-09-19 (v5.14.6). Results committed to
`docs/walls/register.md` § "Architecture-swap audit: composed → prolate-native diatomic"
— that section is the durable record; this file is the (transient) spec. Original spec
below.

**Status (original):** SPEC for a PI-authorized scoped `/walls` pass (PI: "go for it", 2026-09-19).
Reads `docs/failed_approaches_ledger.md`; writes PROPOSED refinements to
`docs/walls/register.md` only (never the ledger — append-only, §13.5). PI adjudicates
every MIS-SCOPED and every guardrail-flagged verdict before it is committed.

## Why this pass exists

Positive claims carry five-plus staleness gates (C16/C17/C21/C24 + the `cited_by`/
`rests-on` dependency graph). Negative results carry deletion-protection (append-only,
"don't re-derive") but **no staleness mechanism and no dependency edge** naming which
architectural components a wall's proof used. So a wall proven on architecture X can
silently stop applying once X is swapped, and nothing surfaces it. The corpus's move
away from composed toward prolate-native diatomics (v5.13+ Paper 12 work) is exactly
such a swap. This pass retroactively fills the missing `rests-on:` edge for the
affected rows and asks, per wall: does its proof's load-bearing component survive the
swap?

## PRIMARY GATE — where the audit has purchase (PI caution, 2026-09-19)

The audit only bites where an architecture *choice* exists.

- **Atom / diatomic:** prolate-native, NOCI, and exact-integral routes are genuine
  alternatives to composed → a composed-proven wall MAY no longer apply → **IN SCOPE.**
- **Triatomic+ (3+ centers):** composed is the *sole* architecture (prolate spheroidal
  holds exactly two foci; there is no alternative). A composed-proven wall there is a
  **permanent constraint, not a stale artifact.** Re-checking it is wasted motion.
  → **OUT OF SCOPE; STANDS by default.**

Consequence: verdicts are per **(wall × system-class)**. The same wall — e.g. "PK is
the irreducible cost of composed-geometry factorization" — is MIS-SCOPED for diatomics
(prolate-native needs no PK) and STANDING for triatomics (composition is forced). The
pass fills only the atom+diatomic column; the triatomic column is pre-set STANDS.

## INTEGRITY RULE — non-negotiable, the whole risk of the exercise

**MIS-SCOPED means the OLD negative was proven on machinery the target does not use, so
it does not INFORM the new program. It does NOT mean the approach works.** Reactivation
of any re-scoped approach requires a fresh positive test — never the mere removal of the
old wall (`memory/feedback_validate_before_reducing.md`). Every MIS-SCOPED verdict
carries an explicit anti-laundering note. A pass that quietly resurrects dead ends has
failed, regardless of how many rows it "clears."

## Target architecture: prolate-native diatomic, componentwise

| Component | Composed (old) | Prolate-native diatomic (target) |
|---|---|---|
| Geometry | fiber-bundle per electron group | prolate spheroidal, two foci |
| Core | PK pseudopotential | explicit, or absent (H2/H2+) |
| Electron space | factorized, incompatible coords | one shared coordinate system |
| Basis | hydrogenic per-n / Z_eff | Laguerre x Legendre (Paper 12) |
| V_ee | composed / cross-block | Neumann, closed-form {E1, gamma, ln} |
| Orthogonality | Lowdin across centers | native |
| Ceiling | polyatomic (walled) | two foci = DIATOMIC, hard stop |

## Classification rubric — per candidate row

1. `rests-on:` — architectural component(s) the ORIGINAL proof used. Vocabulary:
   COMPOSED / PK / NESTED / CONCAT / LOWDIN / HYDROGENIC-PER-N /
   ADIABATIC-HYPERSPHERICAL / GEOMETRY-2FOCI / TWO-BODY-CUSP / SYMMETRY-SPARSITY /
   PROLATE-NATIVE. (This is the field the corpus has never had for negatives.)
2. `system-class:` atom / diatomic / triatomic+  (classify only atom + diatomic)
3. `retained-by-target?:` yes / no / partial
4. `verdict:`
   - **STANDING** — rests on a component the target keeps; still binds the new program.
   - **MIS-SCOPED** — proven on a dropped component; does not inform the new program
     (+ mandatory anti-laundering note).
   - **OPEN-LEANING** — premise genuinely changed, answer untested; needs a fresh test.
   - **STANDS-BUT-ORTHOGONAL** — still true, but constrains a different goal (usually
     qubit-resource sparsity, not bond length).
   - **SUPERSEDED** — already overtaken.
5. `guardrail?:` — if the row belongs to Papers 8-9 / FCI-M / Track DF, FLAG for PI;
   never silently re-scope a guardrail negative.
6. `anti-laundering-note:` — mandatory on every MIS-SCOPED.

## Candidate set (atom + diatomic slice ONLY)

Provisional triage; the pass verifies each row against its actual ledger text.

**Tier 1 — STANDS and binds the prolate program (do NOT relax; the real constraints).**
Recent prolate-native + geometry-level:
- Two-block per-block radial exponent, Paper 12 H2 prolate CI (v5.14.2) — PROLATE-NATIVE.
- Explicit r12 (James-Coolidge) on the prolate 2e basis / FD (FD leg re-scoped v5.13.10;
  the "analytical rises, FD falls, meet ~86%" lesson is live) — PROLATE-NATIVE.
- Pairing explicit r12 with the Neumann algebraic V_ee (v5.13.9) — PROLATE-NATIVE.
- Cusp treatments (3) — TWO-BODY-CUSP (basis/geometry-independent). STANDING in
  substance; re-scope only the hyperspherical (alpha, theta12) *parametrization*.
- Three-center genus jump / two-foci ceiling — GEOMETRY-2FOCI. STANDING; caps the
  program at diatomics.

**Tier 2 — candidate MIS-SCOPED (diatomic applications only).**
- PK modifications (6) — PK. Anti-laundering note applies hard: a prolate 4e LiH is
  still a hard 4-electron two-center CI; PK negatives do not make it easy, they just
  do not bear on it.
- Inter-group antisymmetry (3) — COMPOSED factorization; prolate 2e has one coordinate
  system so antisymmetry is native. (A point FOR the target.)
- Coupled composition / l_max-via-2D — COMPOSED.
- Nested LiH x3 (single-center, charge-center, heterogeneous) + LCAO concatenation —
  NESTED / CONCAT. GUARDRAIL-FLAGGED (Papers 8-9, FCI-M, Track DF).
- W1e / second-row-binding cluster (~10: kernel-shape, Pauli-orthogonality, mean-field
  J-K, Schmidt, [Ne] correlation, basis-enlargement, DMRG-FCIDUMP, explicit-core-HF,
  off-diagonal cross-block h1, NaH Z_orb) — COMPOSED heteronuclear binding.
- Lowdin retrofit / non-orthogonal encoding / Sturmian-CI 1-norm — likely
  STANDS-BUT-ORTHOGONAL (bind the qubit-RESOURCE claim, not bond length).

**Tier 3 — OPEN-LEANING (premise changed, answer untested).**
- Full N-electron radial solvers (adiabatic / coupled-channel / 2D) — ADIABATIC-
  HYPERSPHERICAL; prolate diatomic uses Neumann-CI, a different solver. The "angular
  basis is the bottleneck" lesson may or may not travel.
- Hydrogenic-per-n / completeness (k_n = Z/n) — the completeness lesson could travel to
  a prolate multi-exponent set; untested there.

## Out of scope
- Triatomic+ rows (polyatomic coupling x3, gerade lever, PK-on-triatomic, H2O/BeH2
  composition) — composed is the sole architecture; walls STAND (PI caution).
- Non-chemistry clusters (QED / gravity / periods / alpha / nuclear).
- Positive-claim staleness (owned by C16/C17/C21).

## Output
Proposed entries under a new `docs/walls/register.md` section
"ARCHITECTURE-SWAP: composed -> prolate-native diatomic", one per row with the six
rubric fields. MIS-SCOPED and guardrail-flagged rows are presented to the PI for
adjudication BEFORE commit. The filled `rests-on:` column is the seed of the general
negative-side dependency mechanism (Part 2 of the 2026-09-19 discussion); if it earns
its keep here, that is the evidence for making the field standing.

## How it runs
One read+classify pass over the ~30 atom/diatomic chemistry rows. PI-invoked `/walls`
(scoped); the PM drafts, the PI adjudicates. No `geovac/` code touched; no ledger row
edited.
