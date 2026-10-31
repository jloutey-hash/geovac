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
