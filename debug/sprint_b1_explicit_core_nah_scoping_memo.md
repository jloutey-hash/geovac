# Sprint B.1 — Explicit-Core NaH HF (scoping memo)

**Date:** 2026-06-07.
**Sprint goal:** Test the decisive Thread-B hypothesis: does the NaH chemistry wall close when the [Ne] core is un-frozen and treated explicitly at HF level?
**Decision gate (from Plan agent A/B scoping):** HF on explicit-core NaH PES exhibits an interior minimum on R ∈ [3.0, 4.5] bohr.  D_e accuracy is secondary; topology of binding is the question.
**Estimated effort:** 4 days (Plan agent estimate); revised after scoping → ~3-5 days.

---

## §1 What "explicit-core NaH" means structurally

Current NaH spec (`molecular_spec.hydride_spec(Z=11)`):
- Core treatment: `'frozen'` → uses `FrozenCore(11)` to compute an effective potential, freezes the 10 core electrons (1s² 2s² 2p⁶) as a screening density
- Valence: 2 electrons (Na 3s + H 1s) in a bond block at `Z_eff ≈ 2.2` (Clementi–Raimondi screened exponent)
- Total encoded electrons: 2 (the 10 core electrons live in the ecore constant + V_cross_nuc(R))

Explicit-core NaH:
- Core treatment: un-frozen, 10 core electrons in explicit core block at `Z_center = 11` (bare Na nucleus)
- Valence: 2 electrons (Na 3s + H 1s) in a bond block at the same `Z = 11` (no screening — the core electrons screen via the encoded Hamiltonian itself)
- Total encoded electrons: 12
- Total qubits at max_n=2: 5 core spatial orbitals + 5 valence spatial orbitals (3s + bond + ...) ≈ 20-30 qubits

The qubit count and Hilbert-space dimension grow substantially. At HF level (mean-field), we never construct the full FCI Hilbert space — we only need to diagonalize a Roothaan SCF equation iteratively, which is cheap even for 12 electrons.

## §2 What needs to change in production code

The Plan agent estimated 350 lines.  Mapping that to concrete edits:

### §2.1 `geovac/molecular_spec.py` (~80 lines)

Add `core_method='explicit'` override for normally-frozen-core elements:

```python
def hydride_spec(
    Z: int,
    R: Optional[float] = None,
    max_n: int = 2,
    include_pk: bool = True,
    core_method: str = 'pk',           # existing
    force_explicit_core: bool = False, # NEW
) -> MolecularSpec:
    ...
    if force_explicit_core and core_type == 'frozen':
        # Build the spec as if the element were explicit-core.
        core_type = 'explicit'
        n_core_e = _NEON_CORE_SIZE_FOR_Z[Z]   # 10 for Z=11, 12 for Z=12 (no), etc.
        # Core block: hydrogenic 1s, 2s, 2p_{-1}, 2p_0, 2p_{+1} at Z_center = Z (bare)
        # Valence block: 3s at Z_center = Z (bare), bond partner H
```

The structural change: the existing `'explicit'` path assumes a 1s² core (n_electrons=2, max_n=2 generating one 1s orbital).  For [Ne] core (n_electrons=10, max_n=2 generating five spatial orbitals: 1s, 2s, 2p_{-1,0,+1}), the core block needs to be parameterized differently.

A clean way: add `core_n_electrons` and `core_max_n` to `OrbitalBlock.block_type='core'`.  At max_n=2 the orbital count is already 5 spatial = 10 spin-orbitals = exact [Ne] occupancy.  The block_type='core' machinery already handles arbitrary `n_electrons` at the dataclass level; the bug is just that `_FIRST_ROW_CORE_ENERGY` hardcodes the He-like 2-electron core energy, and the spec factory hardcodes 1s² populations.

### §2.2 `geovac/composed_qubit.py` (~50 lines)

The composed builder iterates over blocks and assigns orbital indices, builds h1 from kinetic + V_ne to `Z_nuc_center`, and builds eri tensor block-diagonally.  For explicit-core NaH:

- Core block h1: kinetic + V_ne(Z=11) for each of the 5 spatial orbitals (1s, 2s, 2p×3).  The diagonal values are the bare hydrogenic eigenvalues `-Z²/(2n²) = -60.5, -15.125, -15.125, -15.125, -15.125` Ha.  Sum of the 10 core spin-orbitals occupancy gives total bare-core energy `≈ -181.5 Ha`.  This is more negative than the [Ne] core's true energy `-128.5 Ha` because the electron-electron repulsion is not yet accounted for — that's the eri tensor's job, which adds back `+53 Ha` of pair repulsion.
- Bond block h1: kinetic + V_ne(Z=11) for valence orbitals (3s + 1s_H).  Now the bond block's `Z_center=11` (bare), not `Z_eff=2.2`.
- Cross-block eri: must include `(core | core | core | core)`, `(core | core | bond | bond)`, `(core | bond | core | bond)`, etc. — the existing balanced builder only has within-block eri + cross-V_ne, so cross-block eri is a new piece of infrastructure for HF.

The complication: existing `balanced_coupled.compute_cross_block_eri` is designed for valence-only cross-block, with 4-index integrals on a single l_max=2 multipole expansion.  For core ↔ valence cross-block eri, we need higher angular content (1s ↔ 3s coupling involves much larger overlap region) and the multipole expansion may need higher L_max.

### §2.3 HF SCF driver in `debug/` (~150 lines)

The actual SCF iteration is straightforward once h1 and eri are in hand:

```python
# Roothaan SCF for 12-electron NaH
F = h1 + sum_{rs} P[r,s] * (eri[p,q,r,s] - 0.5 * eri[p,r,q,s])
# Diagonalize F → new orbitals
# Build new P from occupied orbitals
# Iterate until P converges
E_HF = trace(P @ (h1 + F)) / 2 + ecore
```

This is ~50 lines of NumPy.  No exotic dependencies.  Use direct eigh + density-matrix mixing for convergence.

Driver runs on R panel {2.5, 3.0, 3.566, 4.0, 5.0, 6.0} bohr.

### §2.4 Reference comparison

Compare to:
- Continuous Level 4 reference: doesn't exist for NaH (`ComposedDiatomicSolver` has no `NaH_ab_initio` classmethod).  Need to add one.  ~30 lines in `composed_diatomic.py`.
- CCCBDB experimental: R_eq = 3.566 bohr, D_e = 0.075 Ha.

## §3 What the decision gate actually tests

The Plan agent's gate: "interior minimum on R ∈ [3.0, 4.5] bohr."  This is qualitative — does HF predict NaH is bound?

Three possible outcomes:

**Outcome 1 — HF binds NaH at the right R.**  Frozen-core projection was the culprit.  Path forward: implement FCI on the explicit-core Hamiltonian for higher accuracy (multi-week sprint), and the entire chemistry-engineering thread A becomes "make the explicit-core path practical and accurate" rather than "fix the frozen-core path."

**Outcome 2 — HF binds NaH but at the wrong R (e.g., still over-attractive at small R).**  Frozen-core is partial; another wall lurks at HF level (likely the bond-block fixed Z_eff or the cross-block ERI shape).  This is the case where A.3 (R-adaptive Z_eff) is the natural next step, now applied to the explicit-core spec.

**Outcome 3 — HF still does not bind NaH on explicit core.**  The wall is genuinely structural at HF level.  This rules out frozen-core AND multi-determinant correlation as the source.  The remaining hypothesis is the *basis structure* itself (B.4 bonding-orbital Pauli repulsion against the same-center core, or a deeper issue with the multipole expansion's representation of bond formation).  In this case, the chemistry-engineering thread is essentially exhausted and we pivot to NCG reformulation (Bratteli networks).

Outcomes 1 and 3 are decisive for the strategic direction.  Outcome 2 is the "ambiguous" case that requires another sub-sprint to disambiguate, but informs the next step uniquely (A.3-on-explicit-core).

## §4 First concrete piece (today's scope)

Before committing the 3-5 day implementation, two diagnostics — DONE 2026-06-07, see `debug/b1_basis_sanity_diagnostic.py`.

**4a) Hydrogenic bare-core energy spectrum.**  Bare Coulomb 11-electron energy at Z=11 hydrogenic basis is ~-1200 Ha; real Na atom is -162 Ha; HF needs to recover ~+1000 Ha of e-e repulsion.  This is consistent with the standard "bare hydrogenic at high Z is far from physical" pattern; HF SCF self-consistency closes the gap if the basis has enough variational freedom.

**4b) Na 3s orbital extent.**  Hydrogenic 3s mean radius:
- At bare Z=11: 1.23 bohr (inside [Ne] 2p extent of 0.45 bohr) — wrong region
- At Z_eff=2.2 (Clementi–Raimondi screened): 6.14 bohr — physically right
- HF self-consistency must push the 3s outward via mixing with higher-n.

**Key constraint surfaced by (4b):  max_n MUST be ≥ 3.**  The Plan agent's outline assumed max_n=2, but at max_n=2 the basis has only `n=1` and `n=2` shells — no 3s orbital exists at all.  The 3s is the valence orbital for Na; B.1 is meaningless without it.  At max_n=3:
- 14 spatial orbitals (1 + 4 + 9 across n=1,2,3) → 28 spin-orbitals
- 12 electrons → ~43% filling, comparable to STO-3G NaH
- HF iteration cost trivial (O(14^4) ≈ 4×10^4 ops per iteration)
- Cross-block ERI build cost grows ~16× vs max_n=2 (one-time per R; still feasible)

**Revised time estimate:** 4-5 days (was 4 days).  Single architectural change required vs. Plan agent's outline: extend `force_explicit_core` to also enforce `max_n=3` automatically for frozen-core elements being un-frozen, so the user doesn't have to remember the constraint.

**Outcome of §4 diagnostic: GO with revised parameters.**  The hydrogenic basis at bare Z=11 is a valid HF starting basis provided max_n≥3.  HF SCF should converge.  No structural reason to abandon B.1 before implementation.

## §5 Strategic note

The Bratteli-network deep-dive is running in parallel (background agent dispatched 2026-06-07).  When that returns:
- If it says Bratteli accommodates infinite-dim GeoVac CH triples cleanly, then Track CD becomes a theorem AND we have a structural framework for explicit-core spec construction.  The B.1 implementation becomes a *test* of the Bratteli construction's predictions.
- If it says Bratteli doesn't fit, the engineering thread is still our only path and B.1 must succeed on its own.

Run B.1 implementation in parallel with the Bratteli reading — they don't compete for code-time, only PI-attention.

## §6 Files

### To be created (during B.1 implementation)
- `debug/sprint_B1_explicit_core_nah_driver.py` — HF SCF + PES sweep
- `debug/sprint_B1_explicit_core_nah_memo.md` — sprint results

### To be modified (production code)
- `geovac/molecular_spec.py` — `force_explicit_core` kwarg, extended `_FIRST_ROW_CORE_ENERGY` or new `_NEON_CORE_ENERGY` table
- `geovac/composed_qubit.py` — cross-block eri for core↔valence at HF integrals
- `geovac/composed_diatomic.py` — add `NaH_ab_initio` classmethod for continuous reference (small, ~30 lines)

### Diagnostic-first targets (this turn)
- `debug/b1_basis_sanity_diagnostic.py` — quick check of (4a) and (4b) above before any production-code work
