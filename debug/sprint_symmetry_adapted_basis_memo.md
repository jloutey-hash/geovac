# Sprint Symmetry-Adapted Basis — Hidden Z₂ via Z₂³ irrep decomposition

**Date**: 2026-06-07
**Worktree**: `agent-a3ef4376fd00df85a`
**Verdict**: **POSITIVE** — per-sub-block particle-conservation Z₂'s found in the
hidden-symmetry scan, surviving GF(2) linear-independence audit and producing
additional sprint-scale qubit savings (+24 qubits across the 6-molecule panel)
above the v3.89.0 extended (Hopf × ℓ-parity × atom-swap × inversion) tapering.

## TL;DR

Combined rotation $U_{\rm total} = U_{\rm swap} U_{\rm Hopf}$ acts on a
chemistry-Hamiltonian orbital index space as a single orthogonal map; assigning
each rotated orbital a triple-Z₂ irrep label
$(p_{\rm Hopf},\,p_\ell,\,p_{\rm swap})\in\{\pm1\}^3$ produces a strict
block-diagonalisation of $h_1$ and the four-index ERI tensor (max
off-sector entry $< 10^{-10}$ on all 6 panel molecules). Scanning each sector
block for accidental connected-component sub-structure surfaces the
**per-sub-block particle-conservation Z₂** of the standard composed builder
(`composed_qubit.py` ERIs vanish across sub-blocks by construction) as a
genuinely new — but linearly-dependent on existing stabilizers in some
combinations — class of stabilizers. After GF(2) audit + commutator gate, the
surviving hidden Z₂'s save 1 (BeH₂, H₂O — atom-swap merges sub-blocks) to 8
(CH₄) additional qubits per molecule.

The combined rotation is orthogonal to $< 2\times 10^{-16}$ and
$U_{\rm Hopf}$ commutes with $U_{\rm swap}$ bit-exactly. Spectrum preservation
is verified at machine precision on LiH max\_n=1 (full diagonalisation
$|E_{\rm ext}-E_{\rm eph}|<2\times10^{-15}$).

## Module + tests

- `geovac/symmetry_adapted_basis.py` — new module, additive on top of
  `extended_tapering.py`. Public API:
    * `build_symmetry_adapted_rotation(spec, nuclei=None)` —
      returns `(U_total, sector_labels, parity_hopf, parity_swap, orbital_table, diagnostics)`.
    * `decompose_hamiltonian_by_sector(spec, nuclei=None, builder='composed')` —
      returns dict keyed by `(p_Hopf, p_ell, p_swap)` with intra-sector
      `h1_block`, `eri_block`, `orbital_indices`. Includes a `('__meta__',)`
      entry carrying `U_total`, `h1_rot`, `eri_rot`, off-sector residuals.
    * `find_hidden_z2_in_sector(spec, nuclei=None, builder='composed')` — scans
      each sector's intra-sector $h_1+\mathrm{ERI}$ for connected components;
      emits a candidate Z-string per component-as-anti and audits against the
      full rotated $H$.
    * `extended_plus_hidden_tapered_from_spec(spec, ...)` — runs v3.89.0
      extended tapering, then appends validated hidden stabilizers and re-runs
      the GF(2) drop + commutator audit + final `taper_off_qubits`.

- `tests/test_symmetry_adapted_basis.py` — 13 regression tests covering
  orthogonality, commutation, block-diagonality, sector counts, candidate
  audit, qubit savings, and spectrum preservation.

- `debug/sprint_symmetry_adapted_basis_panel.py` — verification panel driver
  that writes `debug/data/sprint_symmetry_adapted_basis_panel.json`.

## Verification panel

| Molecule | M  | $n_{\rm sect}$ | max off-$h_1$ | max off-ERI | $n_{\rm cand}$ | $n_{\rm valid}$ | $Q_{\rm ext}$ | $Q_{\rm eph}$ | hidden save | Pauli reduction |
|:--------:|:--:|:--------------:|:-------------:|:-----------:|:--------------:|:---------------:|:-------------:|:-------------:|:-----------:|:---------------:|
| LiH      | 15 | 3              | 5.1e-17       | 2.4e-17     | 9              | 9               | 22            | 20            | **2**       | 1.9%            |
| HF       | 30 | 3              | 3.5e-16       | 8.6e-17     | 18             | 18              | 46            | 41            | **5**       | 7.0%            |
| BeH₂     | 25 | 6              | 4.9e-16       | 1.9e-17     | 15             | 3               | 41            | 40            | 1           | 0.4%            |
| H₂O      | 35 | 6              | 2.4e-14       | 3.9e-17     | 21             | 3               | 59            | 58            | 1           | 1.0%            |
| NH₃      | 40 | 3              | 2.2e-16       | 4.4e-17     | 24             | 24              | 62            | 55            | **7**       | 9.6%            |
| CH₄      | 45 | 3              | 2.0e-16       | 4.8e-17     | 27             | 27              | 70            | 62            | **8**       | 7.8%            |
| **Total**|    |                |               |             |                |                 | **300**       | **276**       | **24**      | 0.4 – 9.6%      |

Block-diagonality residuals are at the float64 floor. Atom-swap molecules
(BeH₂, H₂O) save less because the joint-sub-block merging built into
`extended_tapering._build_joint_sb_map` already absorbs the per-sub-block
factorisation that the hidden-Z₂ scan exploits — only sectors with parity
across the swap pair contribute new (independent) candidates.

## Sector decomposition example: BeH₂

BeH₂ has 6 populated sectors of $Z_2^3$ (atom-swap doubles the Hopf×ℓ
3-sector LiH count):

```
sector (+1, +1, +1): dim=6, orbs=[0, 1, 5, 6, 10, 11]   (1s, 2s; symm under swap)
sector (-1, -1, +1): dim=3, orbs=[2, 7, 12]            (2p_{-1}; antisym in m)
sector (+1, -1, +1): dim=6, orbs=[3, 4, 8, 9, 13, 14]  (2p_0, 2p_{+1}; ell-odd)
sector (+1, +1, -1): dim=4, orbs=[15, 16, 20, 21]      (1s, 2s; ANTISYM swap)
sector (-1, -1, -1): dim=2, orbs=[17, 22]              (2p_{-1}; antisym in both m and swap)
sector (+1, -1, -1): dim=4, orbs=[18, 19, 23, 24]
max off-sector h1:  4.86e-16  (bit-exact zero)
max off-sector eri: 1.93e-17  (bit-exact zero)
```

The unpopulated cells $(\pm1,\pm1,\pm1)\setminus$ {6 above} correspond to
orbitals that don't exist in the basis (e.g. there's no $m=0$ orbital with
$\ell$-odd and $m$-antisymmetric simultaneously).

## Structural finding (substantive content)

**Per-sub-block particle conservation is the connected-component signature in
the Z₂³ decomposition.** Each sector block of LiH (no atom-swap) decomposes
into **3 connected components** — one per sub-block (Li\_core, LiH\_bond\_center,
LiH\_bond\_partner). This is the conservation law that

$$
N_{\rm sub\text{-}block} = \sum_{i \in {\rm sub\text{-}block}}\,(n_{i,\alpha}+n_{i,\beta})
$$

is separately preserved when the ERIs are intrinsically block-diagonal in the
sub-block label. The Hopf-per-block stabilizer in v3.52.0 captures only the
$m_\ell\to -m_\ell$ reflection within each sub-block, NOT the per-sub-block
electron number — these are independent Z₂'s in spin-orbital number space,
and the GF(2) audit confirms 2 of the 9 LiH candidates are linearly
independent of the existing alpha/beta/Hopf-per-block/ℓ-parity-per-block
stabilizers.

This is the chemistry-side analog of the Z₂ Hopf-U(1) sprint finding (Z₂ is
the only sub-action of a continuous U(1) that commutes with a real-integer
adjacency): per-sub-block particle conservation is a $U(1)^{n_{\rm sb}}$
group, but the only sub-action that produces a Z-string stabilizer compatible
with the standard JW encoding is its Z₂ subgroup that flips the parity of one
component vs the rest.

## Honest scope

- Hidden Z₂'s found here are **structurally equivalent to per-sub-block
  particle-number parity**, not a "new" symmetry of the underlying physics.
  However, they are independent linearly of the v3.52.0 / v3.89.0 stabilizers
  in the GF(2) sense (2-of-9 for LiH; up to 8-of-27 for CH₄) and produce
  measurable qubit-count savings.
- Spectrum preservation is verified at machine precision on LiH max\_n=1
  (full diag, $\Delta E < 2\times10^{-15}$ Ha). At production max\_n=2 the
  tapered operators are too large for full diag in this sprint; the audit
  gate ($\|[H,P]\|\le 10^{-10}$) is the load-bearing correctness witness.
- BeH₂ / H₂O show low hidden savings (1 each) because the swap rotation
  mixes sub-blocks and the existing `_build_joint_sb_map` mechanism in
  `extended_tapering.py` already absorbs the per-pair-of-sub-blocks
  factorisation. The headline savings concentrate on the non-swap panel
  members (LiH, HF, NH₃, CH₄).
- The detection mechanism (connected components of intra-sector $h_1+$ERI)
  is general — it would surface any further accidental block-structure in
  the rotated Hamiltonian, not just the per-sub-block one. On this panel,
  per-sub-block is the only signature found.

## Gates passed

- Combined rotation orthogonal: $\max|UU^T-I|<2\times10^{-16}$ on all panel
  members.
- $U_{\rm Hopf}$ and $U_{\rm swap}$ commute: bit-exact zero.
- Rotated $h_1$ block-diagonal across $Z_2^3$ sectors: $<2\times10^{-14}$
  (worst case H₂O; structural artefact of multi-center accumulated error,
  still 10⁴× below the 1e-10 audit gate).
- Rotated ERI block-structured: max off-sector $< 9\times10^{-17}$.
- Every emitted hidden-Z₂ candidate passes the commutator audit against
  the full rotated $H$ to $<10^{-10}$.
- Existing v3.89.0 extended\_tapering regression: **12/12 PASS** (unchanged).
- New test suite: **13/13 PASS** in 9.7 s.

## Files

- `geovac/symmetry_adapted_basis.py` (new, ~550 lines)
- `tests/test_symmetry_adapted_basis.py` (new, 13 tests, 9.7 s wall)
- `debug/sprint_symmetry_adapted_basis_panel.py` (driver)
- `debug/data/sprint_symmetry_adapted_basis_panel.json` (panel data, 6 rows)
- `debug/sprint_symmetry_adapted_basis_memo.md` (this file)
- `debug/changelog_staging_symadapt.md`
- `debug/claudemd_staging_symadapt.md`

## Where this lives in the framework

This module sits in the Paper 14 §sec:hopf\_tapering arc. Compared to the
v3.52.0 per-sub-block Hopf tapering and the v3.89.0 extended (Hopf × ℓ-parity
× atom-swap × inversion) tapering, this sprint provides:

1. A **unifying basis-rotation framing** that makes the three known Z₂'s
   manifest as orbital labels (an Z₂³ irrep tuple).
2. A **systematic detection mechanism** for further Z₂'s via connected-
   component analysis within each irrep block.
3. A **production module** that adds the validated hidden stabilizers on
   top of the extended pipeline without modifying v3.89.0 code.

Paper 14 §sec:hopf\_tapering would benefit from a single short paragraph
naming the Z₂³ irrep decomposition and the (1.9 – 9.6%) additional Pauli
reduction it produces above the v3.89.0 baseline on the non-swap panel
members. No new theorems are needed for this sprint; the structural content
(per-sub-block particle conservation as the symmetry that produces the
hidden Z₂'s) is already implicit in the composed builder's block-diagonal
ERI structure.
