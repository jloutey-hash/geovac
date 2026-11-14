# Pre-registration LB: rebuild of the LiH STO-3G qubit baseline

Registered 2026-11-13, before any computation. PI-directed ("We can build LiH. Let's do it.").

## Why

Papers 14, 19, 20 and the group-4 synthesis quote a raw Jordan-Wigner LiH STO-3G baseline
(registry `lih_sto3g_raw` 631 terms, `lih_sto3g_lambda` 16.48 Ha, `lih_sto3g_qwc` 179 groups,
Q = 12, R = 3.015 bohr). The script that produced it is not on file, so the group count could not
be re-measured under the v6.7.24 deterministic QWC tie-break. This record rebuilds it from the
project's own Gaussian integral engine and pins it with a test.

## Construction (fixed before computing)

- Basis: STO-3G for Li (1s, 2s, 2p) and H (1s), standard exponents and contraction coefficients,
  Cartesian p, 6 spatial functions. Li at the origin, H at +R on z, R = 3.015 bohr.
- Integrals: `geovac/noci_engine.py` McMurchie-Davidson engine (`integral_set_md`).
- Orbitals: closed-shell RHF canonical MOs, all 6 kept (no frozen core, no active space).
  The degenerate pi pair is rotated so one MO is pure p_x and the other pure p_y (the operator,
  and so its term count, depends on the orientation inside a degenerate pair; the axis-aligned
  choice is the standard one).
- Qubit operator: `geovac.qubit_encoding.build_fermion_op_from_integrals` +
  `exact_jordan_wigner` (the v6.7.23 drop-free path). Terms counted with the identity.
- Groups: `geovac.measurement_grouping.count_qwc_groups` (v6.7.24 tie-break).

## Independent-route control

The same pipeline at R = 1.45 Angstrom is compared with OpenFermion's shipped
`H1-Li1_sto-3g_singlet_1.45.hdf5` (integrals from a different code). Pass: RHF energy within
1e-6 Ha of the file's `hf_energy`, and the same term count. If the control fails, LB1-LB5 are
VOID, not scored.

## Predictions

- **LB1:** RHF energy at R = 3.015 bohr is -7.8620 +/- 0.0005 Ha.
- **LB2:** 631 Pauli terms, identity included (the registered value).
- **LB3:** 1-norm, identity included, 16.48 +/- 0.05 Ha.
- **LB4:** QWC groups between 172 and 182 inclusive. (Not a prediction that 179 reproduces.)
- **LB5:** the "2.8x fewer QWC groups" sentence (composed LiH 64 groups) survives: groups / 64
  rounds to 2.8, i.e. 176 <= groups <= 182.

## Scoring (2026-11-13, `debug/lih_sto3g_rebuild.py`, `debug/data/lih_sto3g_rebuild.json`)

**Control PASSED.** At 1.45 A our RHF energy is -7.8625677855 Ha against the file's
-7.8625677857 (difference 1.7e-10); both operators give 631 terms, 1-norm 16.45629 Ha
(difference 8e-8) and 178 QWC groups.

| | registered | measured | verdict |
|---|---|---|---|
| LB1 RHF energy | -7.8620 +/- 0.0005 | -7.862009 | HELD |
| LB2 terms | 631 | 631 | HELD |
| LB3 1-norm | 16.48 +/- 0.05 | 16.4768 (12.3423 without the identity, identity -4.1344) | HELD |
| LB4 QWC groups | 172-182 | 177 | HELD |
| LB5 "2.8x" survives | 176-182 | 177 (177/64 = 2.77) | HELD |

The published 179 is the round-off insertion order: replacing the tie-break with the size-only
sort on this same build gives exactly 179. The groups are unchanged under three 1e-13 relative
perturbations (177, 177, 177). The identity coefficient and the non-identity 1-norm match the
2026-11-04 registry provenance to the printed digits, so that recompute was this construction.

**A construction note was wrong.** The record said the operator depends on the orientation inside
the degenerate pi pair. It does not: a rotation inside the pair is a rotation about the bond axis,
a symmetry of LiH, and the MO integrals agree to 1e-16 after a 30 degree rotation. The axis-aligned
choice is harmless, not required.

**Verdict: GO.** Registry `lih_sto3g_qwc` 179 -> 177 (179 registered as retired); Papers 14, 19,
20, the group-4 synthesis, README, claims register and one-pager updated. Guard
`tests/test_lih_sto3g_baseline.py` (2 tests), fire-tested: the size-only order fails both
(179 vs 178 at the control, 179 vs 177 at 3.015), a 1% error in one Li 2p contraction coefficient
fails both (control energy off by 2.1e-5 Ha).
