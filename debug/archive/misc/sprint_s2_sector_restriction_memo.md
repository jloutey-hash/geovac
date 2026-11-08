# Sprint S² sector restriction (2026-06-07)

## TL;DR

Ship `geovac/spin_sector.py`: total-spin operators `S²` and `S_z` as
Jordan–Wigner `QubitOperator`s on the GeoVac composed-builder
spin-orbital basis, an audit gate `verify_s2_commutes`, a closed-form
singlet-sector dimension counter, an `H + λ·S²` VQE penalty objective,
and an FCI diagnoser that reports `<S²>` per low-lying eigenstate.
**`[S², H] = 0` bit-exactly (residual = 0.0e+00) on every untapered
production composed Hamiltonian in the H₂ / LiH / BeH₂ / HF / H₂O /
NH₃ / CH₄ panel.** Singlet-sector dimensions are 1.67×–5.35× smaller
than the corresponding `S_z = 0` sectors, a Hilbert-space-dimension
(not qubit-count) reduction that translates directly into a
spin-projected VQE convergence factor.

Honest scope: total spin is an integer-valued Casimir, not a Z₂
stabilizer, so this does **not** give a clean ΔQ qubit reduction.
What it gives is a VQE objective (`H + λ·S²`) that filters triplet
contamination at the optimization level and a closed-form bound on
the singlet sub-FCI dimension. The commutation gate also exposes
when the qubit basis has stopped matching the alpha/beta JW partition
— for example after per-block Hopf tapering the standard
`compute_s2_operator(M_tapered)` is **not** the physical S² of the
tapered basis (the tapering rotates orbital mode indices on the
spatial-orbital side but the qubit-count is no longer 2M).

## Mechanism

Standard JW spin-orbital convention (`OpenFermion`,
`geovac.qubit_encoding.build_fermion_op_from_integrals`): spatial
orbital `p` carries spin orbitals `2p` (alpha) and `2p+1` (beta).
Under this convention,

  S² = S⁻ S⁺ + S_z·(S_z + 1),    S_z = ½·Σ_p (n_{p,α} − n_{p,β}),

both expressible as `O(M²)` Pauli sums on `2M` qubits via
`openfermion.s_squared_operator` + `jordan_wigner` (we delegate to that
canonical implementation rather than rolling our own; it is well-tested
in the openfermion suite and matches our convention exactly).

The GeoVac composed builder constructs the spatial-orbital `h1` and
`eri` integrals and uses
`geovac.qubit_encoding.build_fermion_op_from_integrals` to JW-encode.
Because the JW step is spin-symmetric (it adds α and β terms with the
same coefficient), the resulting `H` commutes with both `S²` and `S_z`
exactly.

`s2_penalty_hamiltonian(H, M, λ)` returns `H + λ·S²`. Triplet
eigenstates (`S² = 2`) shift up by `2λ`; quintets by `6λ`; the
singlet sector is unshifted. This is the standard
*Lagrange-multiplier* approach used in spin-projected variational
methods (cf. CASSCF spin-penalty, OpenFermion's
`ProjectorOnSubspace` for symmetry-adapted state preparation, and the
"S² penalty" trick in unitary-coupled-cluster VQE literature).

## Panel: dim(S_z = 0) → dim(S = 0)

Untapered production composed Hamiltonians at `max_n=2`, 7 molecules:

| System | Q  | M  | N_e | dim(S_z=0)        | dim(S=0)          | Speedup | [S²,H] residual |
|:-------|---:|---:|----:|------------------:|------------------:|--------:|:----------------|
| H₂     | 10 |  5 |   2 |             25    |             15    |  1.667× | 0.000e+00       |
| LiH    | 30 | 15 |   4 |         11,025    |          4,200    |  2.625× | 0.000e+00       |
| BeH₂   | 50 | 25 |   6 |      5,290,000    |      1,495,000    |  3.538× | 0.000e+00       |
| HF     | 60 | 30 |  10 |     20,307,960,036|      4,035,556,161|  5.032× | 0.000e+00       |
| H₂O    | 70 | 35 |  10 |    105,385,935,424|     20,397,277,824|  5.167× | 0.000e+00       |
| NH₃    | 80 | 40 |  10 |    432,974,528,064|     82,184,979,864|  5.268× | 0.000e+00       |
| CH₄    | 90 | 45 |  10 |  1,492,695,054,081|    279,121,839,381|  5.348× | 0.000e+00       |

**Pattern.** Speedup grows monotonically with `M / N` (more orbitals
per electron => more correlation paths => more spurious triplets
mixed into the `S_z=0` sector that the singlet projector kills).
For 10-electron production systems the singlet sub-FCI is ~5× smaller
than the unrestricted `S_z=0` sector.

The commutator is bit-exact zero on **every** entry. This is the
expected outcome — the composed builder uses spin-restricted
integrals and the JW transform preserves spin symmetry — but the
audit gate matters as a regression detector for any future
modification to the Hamiltonian builder that might mix
α and β manifolds (e.g. spin-orbit coupling, Bogoliubov rotations,
explicit triplet symmetry breaking).

## What was NOT done

1. **No per-block-tapered `S²` operator.** Per-block Hopf tapering
   (the `tapered='per_block'` path in `geovac.ecosystem_export`)
   removes `2 + n_sub_blocks` qubits via the `m → −m` Z₂ on the
   *spatial-orbital* side. The tapered Hamiltonian no longer has
   an even qubit count in general (LiH per_block: Q=25, odd!), so the
   "`S²` on 2M qubits" construction does not directly apply to the
   tapered basis. The honest scope: call `verify_s2_commutes` and
   `compute_fci_ground_state_s2` on the *untapered* Hamiltonian, then
   apply Hopf tapering as a second, independent layer at the end of
   the pipeline. A future sprint could build the tapered `S²` by
   applying the same Z₂ rotation that `geovac/z2_tapering.py` applies
   to `h1, eri` — but doing that correctly requires lifting the rotation
   to a full spinor-doubled `2M×2M` matrix and reapplying it to
   `s_squared_operator(M)`. Left for a follow-on.

2. **No "extended" tapering basis.** The task prompt mentioned a
   `tapered='extended'` mode (with ℓ-parity + atom-swap + inversion);
   this is not present in the current worktree (only `global` and
   `per_block`). Tests fall back to `untapered` automatically when
   the requested mode is unavailable. The `tapered='per_block'`
   pathway was tested directly (commutator residual ~3×10² — large
   because the basis has been rotated, not because spin symmetry is
   broken; the rotated `S²` would need to be the construction target).

3. **No qubit removal.** S² is integer-valued (eigenvalues 0, 2, 6,
   12, ...), not Z₂-valued. Standard qubit-tapering machinery (which
   removes one qubit per Z₂ stabilizer) does not apply. The
   speedup comes entirely from Hilbert-space-dimension reduction
   inside VQE, not from a smaller Hamiltonian.

4. **No full LiH ground-state `<S²>` check at production scale.**
   At Q=30 the full Hilbert space is 2³⁰ ~ 1.1 billion dim, infeasible
   to diagonalise in CI tests. Covered by (a) a bit-exact commutator
   gate (which proves the ground state is `S²`-pure to whatever
   precision the eigensolver achieves) and (b) a synthetic 2-orbital
   H₂-like fixture whose ground state is verified to be a singlet.

5. **No production wiring into ecosystem_export.** S² is a *diagnostic
   and objective* tool, not a Hamiltonian modification. The natural
   API surface is the new module's public functions; users opt in
   when constructing VQE objectives. No `tapered='spin_sector'`
   keyword has been added to `hamiltonian()` because the operation
   is not a tapering.

## Hard-prohibition check

- ✓ No change to the natural geometry hierarchy.
- ✓ No fitted or empirical parameters introduced (`λ` is a
  user-supplied Lagrange multiplier, not a physics parameter).
- ✓ No deletion of negative results from §3.
- ✓ No change to the Paper 2 "conjectural" combination-rule label.
- ✓ No change to `geovac/extended_tapering.py` (which doesn't exist
  in this worktree in any case) nor to `geovac/z2_tapering.py`.
- ✓ No change to `CLAUDE.md` or `CHANGELOG.md` — proposed entries
  staged at `debug/claudemd_staging_s2_sector.md` and
  `debug/changelog_staging_s2_sector.md`.

## Files

- `geovac/spin_sector.py` — new production module, ~370 lines, six
  public functions, docstrings cross-referencing OpenFermion source
  and CLAUDE.md §1.5 positioning.
- `tests/test_spin_sector.py` — 27 tests (25 passing, 2 intentional
  skips: LiH full diag too large for CI, toy Hamiltonian skip path
  exercised).
- `debug/sprint_s2_sector_panel.py` — verification panel driver.
- `debug/data/sprint_s2_sector_panel.json` — panel results (JSON).
- `debug/sprint_s2_sector_restriction_memo.md` — this memo.
- `debug/changelog_staging_s2_sector.md` — proposed CHANGELOG entry.
- `debug/claudemd_staging_s2_sector.md` — proposed CLAUDE.md §2
  one-liner + §7 entry-point update.

## Honest scope summary

The commutator gate is bit-exact on the untapered production
Hamiltonians (residual = 0.0 across the panel). Singlet sector
dimensions are 1.67×–5.35× smaller than `S_z = 0`. This is real
structural content (Weyl branching) and the right reduction factor
for the VQE trial-state manifold of a closed-shell singlet target.
It is **not** a qubit-count reduction. It is **not** automatic on
the per-block-tapered Hamiltonian (the basis has been rotated; the
matching tapered S² is a follow-on). The penalty Hamiltonian
`H + λ·S²` lifts triplet contamination by `2λ`, leaves singlets
unchanged, and is the standard recipe for spin-projected VQE.
