# Tier 2 Track T5 — Verdict Memo

**Sprint:** Dirac-on-S³ Tier 2, Track T5.
**Date:** 2026-04-15.
**Status:** Complete. Certifier landed; T3 symbolic Hamiltonian certifies
cleanly; Paper 18 drop-in drafted; 25/25 new tests + 51/51 D1 regression
tests pass.

## 1. Summary

T5 extends the Tier 1 π-free rational certificate
(`geovac.dirac_s3.verify_pi_free`) from the scalar Camporesi–Higuchi
spectrum to the full Breit–Pauli spin-orbit Hamiltonian built by T2.
The new certifier `verify_spinor_pi_free` walks the sympy expression
tree of each coefficient and rejects any transcendental content outside
the Paper 18 spinor-intrinsic ring

    R_sp := Q(α²)[γ] / (γ² + (Zα)² − 1)       (Eq. 1 of the drop-in)

with `{Z, α, γ}` as the registered free-symbol set.

## 2. Does T3's symbolic Hamiltonian pass?

**Yes, cleanly.** When `build_so_hamiltonian_block(n_max, Z=Z_sym,
alpha=alpha_sym)` is run with symbols kept free (no numerical binding),
every matrix element lies in `R_sp`. Tested for n_max ∈ {1, 2, 3, 4},
every l ≥ 1 shell, and every κ branch (κ<0 and κ>0 both accepted).

Specific diagonal eigenvalues audited:

| State | κ | Symbolic expression | Ring membership |
|:------|:-:|:--------------------|:----------------|
| 2p_{3/2} | −2 | +Z⁴α²/96 | Q(α²)·Q → ✓ |
| 2p_{1/2} | +1 | −Z⁴α²/48 | Q(α²)·Q → ✓ |
| 3d_{5/2} | −3 | +Z⁴α²/540 | Q(α²)·Q → ✓ |
| 3d_{3/2} | +2 | −Z⁴α²/360 | Q(α²)·Q → ✓ |
| 1s_{1/2} | −1 | 0 (Kramers) | Integer(0) → ✓ trivially |

All agree with the T2 design memo §Hydrogen fine-structure benchmark.
No coefficient surfaced any π, ζ, log, E₁, or unregistered symbol.

## 3. Obstruction report

**None at the Breit–Pauli leading order.** Every coefficient closed as
`rational × α² × hydrogenic radial seed`, which is exactly the ring the
drop-in proposes. γ is present in T1 only as a reserved symbol (not
bound to any T2 coefficient), but the certifier accepts it and its
explicit `sqrt(1-(Z·α)²)` form — it is ready for Tier 3+ when
Martínez-y-Romero radial recursions are wired in.

**Potential future obstructions (flagged, not hit):**

- **α⁴ Darwin + mass-velocity.** When/if these enter, they will add
  `α⁴` powers to `R_sp`. The ring is closed under arbitrary powers of
  α², so no taxonomic change is needed — the certifier accepts them.
- **α³ Lamb-shift / QED corrections.** If these enter, `R_sp` expands
  to `Q(α)[γ]/(γ² + (Zα)² − 1)` (i.e., α itself, not just α²). The
  certifier's registered symbol set already includes bare α; the
  taxonomic classification remains spinor-intrinsic.
- **Breit-interaction two-body corrections.** Introduce new radial
  integrals but same (α², γ) transcendental content; stays in `R_sp`.

None of these break the subtier; all populate it further.

## 4. The deliberate contamination test

Test `test_rejects_pi_in_so_block` injects `sp.pi * alpha_sym**2 / 96`
into a real T2-shaped dict (replacing a legitimate κ≠−1 coefficient)
and confirms the certifier raises `SpinorTaxonomyError` with a π-aware
diagnostic. Six other negative controls cover: bare π, π², ζ(3), log(Z),
expint (E₁), and an unregistered free symbol. All six raise cleanly.

## 5. Drop-in as drafted

`docs/paper18_spinor_subtier_proposal.tex` contains:

1. A new `\paragraph{Spinor-intrinsic content from first-order
   operators on a spinor bundle (new subtier).}` that sits immediately
   after the Tier-1 odd-zeta paragraph in §IV.
2. A closed-form derivation trail: Eq. (BP_SO) → Eq. (BP_SO_matrix) →
   Eq. (gamma_def), with the ring definition stated explicitly.
3. Four itemized criteria (source / mechanism / dependency / algebraic
   status) distinguishing the subtier from calibration, embedding,
   flow, and odd-zeta.
4. An updated operator-order × bundle structural table
   (Table `tab:bundle_operator_tiers`):

   | Tier | Examples | Operator order | Bundle |
   |:---|:---|:---:|:---|
   | Intrinsic | π^(d/2) Weyl | — | scalar |
   | Calibration | κ = −1/16, F = π²/6 | 2nd | scalar on curved space |
   | Embedding | e^a E₁(a), 1/r₁₂ | — | graph → R³ projection |
   | Flow / composition | µ(ρ,R) piecewise, PK | — | multi-geometry |
   | Odd-zeta (Tier 1 new) | ζ(3), ζ(5), … | 1st | Dirac spectrum, scalar |
   | **Spinor-intrinsic (Tier 2 NEW)** | **α², γ = √(1−(Zα)²)** | **1st** | **spinor bundle** |

5. A 2×2 operator-order × bundle summary table that makes the structural
   parallel to Paper 24 explicit. The lower-left cell ("2nd-order on
   spinor bundle") is flagged as not yet encountered — Tier 3+ may or
   may not populate it.
6. A "Certification and taxonomy wiring" closing paragraph
   cross-referencing the T5 certifier module and test file.

## 6. Recommendation to T6

**Apply the Paper 18 drop-in as part of the sprint commit.**

Rationale:

- **Structurally clean.** Every claim in the drop-in is backed by an
  exact sympy identity + unit test in `tests/test_spinor_certificate.py`.
- **No conflicts with existing Paper 18 content.** The drop-in is
  strictly additive: a new `\paragraph` + a new table, adjacent to (not
  replacing) the Tier-1 odd-zeta paragraph. The existing "Parallel to
  Paper 24" paragraph remains valid as a discussion of the odd-zeta
  tier; the new spinor-intrinsic subtier gets its own "Parallel to
  Paper 24" sub-paragraph at the end of the drop-in.
- **Completes the Tier 1/Tier 2 taxonomic parallel.** Without T5's
  drop-in, the Dirac-on-S³ sprint results are only half-captured in
  Paper 18. Tier 1 added the spectrum-side transcendental (ζ(3)); T5
  captures the coupling-side transcendental (α², γ).
- **Mechanically certified.** `verify_spinor_pi_free` is not a
  documentation artifact — it is an enforceable post-commit gate that
  can be invoked on any future Hamiltonian to confirm it stays in the
  ring. This matches Paper 24's pattern of coupling a structural claim
  to an executable certifier.

**Minor edit flagged for T6's editorial pass:**

- The drop-in's Table `tab:bundle_operator_tiers` duplicates structure
  with the existing `tab:taxonomy` earlier in §IV. Consider whether to
  keep both (they emphasize different axes — `tab:taxonomy` is per-Level
  in the GeoVac hierarchy, while the new table is operator-order × bundle)
  or merge. Current recommendation: **keep both**; they address different
  reader questions ("where does this Level's transcendental come from?"
  vs. "how is the transcendental content structurally classified?").

## 7. No new transcendentals surfaced beyond T1/T2/T3 memos

The T1 memo flagged α and γ as reserved symbols; T2 memo confirmed α²
enters the Breit–Pauli prefactor. T3 memo noted the composed pipeline
inherits the same transcendental content. T5 confirms this: no
coefficient outside the `{α, γ, rational, hydrogenic radial seed}` set
appears anywhere in the T2/T3 symbolic path. The one flagged item from
T3 memo §8 (the PK quadrature is still grid-based) is orthogonal to the
taxonomic certification: PK is a composition/flow exchange constant,
already classified in the existing taxonomy, and its numerical
quadrature is a production-code implementation detail, not a new
transcendental tier.

## 8. Files

- `geovac/spinor_certificate.py` — new module (the T5 deliverable).
- `tests/test_spinor_certificate.py` — 25 tests, all passing.
- `docs/paper18_spinor_subtier_proposal.tex` — drop-in for §IV, ready
  for T6.
- `docs/tier2_t5_verdict.md` — this memo.

## 9. Guardrails observed

- Did NOT modify `papers/core/paper_18_exchange_constants.tex`.
- Did NOT modify any other paper.
- Did NOT modify `geovac/dirac_s3.py`, `geovac/dirac_matrix_elements.py`,
  `geovac/spin_orbit.py`, `geovac/composed_qubit.py`,
  `geovac/composed_qubit_relativistic.py`, `geovac/molecular_spec.py`.
  (The certifier is a new module; everything upstream remains locked.)
- Did NOT classify α² as "calibration" (it is a new subtier, per the
  PI-recorded Tier 2 Decision 6).
- Did NOT introduce new transcendentals — only classified existing
  T1/T2/T3 content.
- No TC, no Dirac Fock projection, no graph Dirac operator, no S⁵/S⁷.
