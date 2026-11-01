# Pre-registered prediction: Dong's gauge group of the true Camporesi-Higuchi window (WG)

**Registered 2026-10-31** (PI: "Love your plan. Let's do it", the research-thread ranking that put this first),
before any dimension of the window's gauge algebra is computed.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 32 `prop:odd` and the scope paragraph after `thm:h1_yukawa_non_selection` (the YC wall); Paper 38
(the true window, dim O = 5, 30, 91); record `docs/predictions/2026-10-29_ccm_yukawa_coefficient_window.md` (YC);
CHANGELOG v5.21.27 (the lit scan: Dong, Choi, Verhoeven) and v5.21.30 (the edge-crossing guard).
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope; no molecular encoding.

## Why

The YC wall says the compressed fluctuation on the true window is not gauge invariant, because a gauge function u
acts on the window only if it commutes with the spectral projection (Choi). That says the LOCAL gauge group is
gone. It does not say what survives. Dong (SIGMA 18 (2022) 060, Lemma 4.1, Prop. 4.2, Cor. 4.3) answers this for
the truncated circle: the gauge group of an operator system E, G(E) = {U unitary in C*(E) : U E U* = E}, is
U(1) × (U(1) ⋊ Z₂) for the Toeplitz system at every n, i.e. a phase, the rotations and the reflection: only the
isometries survive. Nobody has computed it for the true window of S³. This record does.

## Method (fixed before scoring)

- O^CH: the compressions to the true window (complete Dirac eigenspaces, |λ| ≤ n_max + 1/2) of the band-limited
  multipliers, built by the tracked Paper 38 guard's Route A (`tests/test_p38_true_ch_window.py::_window`);
  dim H = 4, 16, 40, 80 and dim O = 5, 30, 91, 204 at n_max = 1, 2, 3, 4. C*(O) = M_N (Verhoeven 2026, Thm 4.10).
- Scored object: the Lie algebra of the identity component, g = {X anti-Hermitian : [X, e] ∈ O for all e ∈ O}.
  Linear in X. Because O is *-closed, the complex solution space has complex dimension = real dim g. Computed as the
  null space of the positive semidefinite Gram form G = Σ_k L_k†(1 − P_O)L_k over an orthonormal basis e_k of O,
  L_k X = [X, e_k]. An eigenvalue of G counts as zero if below 1e-9 × λ_max(G).
- Validity rule: a point is VOID if the spectrum of G has an eigenvalue between 1e-9 and 1e-4 of λ_max (no clean
  gap).
- Driver: `debug/gauge_group_window.py` (`control`, `control_complex`, `inclusion` run; `score` not run).

## What I know before comparing (disclosed)

- **Controls (measured, non-scored):**
  - Dong's own case, the Toeplitz system Toep_n, n = 2..7: dim g = 2 at every n (Dong: U(1) × U(1) identity
    component), null eigenvalues ≤ 7e-16, smallest nonzero eigenvalue 0.80 to 2.0.
  - The same systems conjugated by a random unitary (complex matrices): dim g = 2 at every n, null ≤ 3.3e-15, gap
    ≥ 0.80. Added because a conjugation bug in the basis routine passed the real control and was caught by the
    inclusion check below; fixed before any scored run.
- **Inclusion (measured; a lower bound, not the scored number):** seven known anti-Hermitian operators lie in g at
  n_max = 1, 2, 3, with relative residual ≤ 1.6e-15, and they are linearly independent (rank 7): the identity (the
  phase), the three left translations (spin-j matrices on m) and the three right translations with spin
  (J_R on m′ plus σ/2), all of which commute with D to 2e-15 and so preserve the window. So dim g ≥ 7.
- No dimension of g on the window has been computed.
- What I expect and why: the commutant of O is the scalars (Paper 38, rem. on the Dirac kernel; Choi), so nothing
  beyond the phase acts trivially; the isometries of round S³ that preserve the spin structure and D form
  Spin(4) = SU(2)_L × SU(2)_R, of dimension 6. Dong's circle has exactly phase + isometries. I expect the same.

## Predictions

**WG0 (controls; measured, PASS).** Toeplitz and conjugated Toeplitz: dim g = 2, n = 2..7.

**WG1 (only the isometries and a phase survive).** At n_max = 2, 3 and 4: dim g = 7, with a clean gap.
Since the seven known generators lie in g (rank 7), dim g = 7 means g is exactly u(1) ⊕ su(2)_L ⊕ su(2)_R.
- Falsifier: dim g > 7 at any of the three (a symmetry of the truncation that is not an isometry: interesting),
  or dim g < 7 (would contradict the inclusion check: a bug).

**WG2 (reported, not predicted).** dim g at n_max = 1, where the window is four-dimensional and an accidental
larger symmetry is plausible.

**WG3 (reported, not predicted).** Whether an orientation-reversing element (the analogue of Dong's reflection)
lies in G, tested by solving the linear intertwining problem U e = ι(e) U for the inversion g ↦ g⁻¹ on functions.
Reported only if computed; not gated.

## Decision

- **GO:** WG0 passes and WG1 holds at all three n_max.
- **BORDERLINE:** WG1 holds at two of three.
- **STOP:** WG1 fails at two or more, or any scored point is VOID.

## Reading if GO, stated in advance

On the true window the gauge group's identity component is U(1) · Spin(4): an overall phase and the isometries,
the S³ analogue of Dong's circle. The local gauge group C(S³, U(1)) has no remnant beyond its constants. This
sharpens the YC wall from "local gauge invariance fails" to "exactly the global isometries survive". It does not
repair the instrument and changes no paper claim on its own.

## Disposition (2026-10-31)

Runs after the record was committed (2e3ba9c2): `python debug/gauge_group_window.py score 1 2 3 4` and
`reflection 1 2 3`. Rows in `debug/data/gauge_group_window.jsonl`, log `debug/data/gauge_group_window_score.log`.

| n_max | dim H | dim O | dim g | largest null eigenvalue | smallest nonzero eigenvalue | λ_max(G) | known 7 in g |
|:--:|:--:|:--:|:--:|:--|:--|:--|:--:|
| 1 (reported) | 4 | 5 | 7 | 4.2e-16 | 1.000 | 4.000 | rank 7 |
| 2 | 16 | 30 | 7 | 2.4e-15 | 0.899 | 6.137 | rank 7 |
| 3 | 40 | 91 | 7 | 1.4e-14 | 0.832 | 7.929 | rank 7 |
| 4 | 80 | 204 | 7 | 2.7e-14 | 0.769 | not printed | rank 7 |

- **WG0: PASS** (as registered).
- **WG1: HELD** at n_max = 2, 3, 4. dim g = 7, so g = u(1) ⊕ su(2)_L ⊕ su(2)_R exactly. Validity: no eigenvalue in
  the void band at n_max ≤ 3 (smallest nonzero / λ_max = 0.25, 0.15, 0.11). At n_max = 4 λ_max was not printed; the
  band would need λ_max > 7,700 for 0.769 to fall in it, against 4.0, 6.1, 7.9 at n_max = 1, 2, 3.
- **WG2 (reported):** dim g = 7 at n_max = 1 as well; no accidental symmetry on the smallest window.
- **WG3 (reported):** the inversion g ↦ g⁻¹ is implemented. For each of three sign conventions of the inversion on
  matrix coefficients (Condon-Shortley sign, no sign, plain index swap), the intertwining equation U e = ι(e) U has a
  one-dimensional solution space, unitary to 2e-14, and that U ANTICOMMUTES with D (to 7e-14; ‖[U, D]‖ = 6.0, 18.3,
  38.9), as an orientation reversal must. Control: the identity map gives the scalars only (null dimension 1, U
  commuting with D). Not shown: that these are the only components.
- **Gate: GO.**

Reading (as stated in advance, now held): on the true window the identity component of Dong's gauge group is
U(1) · Spin(4), the phase and the isometries of round S³, and the group also contains an orientation-reversing
element that anticommutes with D, the analogue of Dong's reflection of the circle. The local gauge group C(S³, U(1))
leaves nothing beyond its constants. This sharpens the YC wall from "local gauge invariance fails" to "exactly the
global isometries (and a reflection) survive". It does not repair the instrument.

## Erratum (2026-10-31, after the disposition)

The registration's control disclosure says a conjugation bug in the basis routine "passed the real control and was
caught by the inclusion check". That cause is wrong. Writing the guard showed that conjugating the basis vectors
changes nothing here: Toep_n and the true window's O are both closed under entrywise conjugation (and transposition),
measured to 5e-15 at n_max = 1, 2, 3, so the routine's conjugation was harmless for every system this record touches.
The inclusion check had failed for a different reason, a 0/0 relative residual for the identity, which commutes with
all of O; that was fixed (zero commutator returns zero) before any scored run. The conjugated-Toeplitz control is
therefore not a test of the conjugation; it stays as a second control of the Gram form on complex matrices. No
scored number, verdict or disposition changes.

