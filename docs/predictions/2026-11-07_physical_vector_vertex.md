# Pre-registered prediction: the physical vector-photon vertex on S^3 (VV, stage A)

**Registered 2026-11-07** (PI: "Ok, do stage A"), before any scored case is computed.
**Commit status:** committed on its own, before the driver exists.
**Rests on:** Paper 28 sec. "Vertex selection rule" (eq:vertex_rule is a vertex model, the scalar-harmonic triangle;
the physical vertex is stated to allow q = n1+n2+1 and to couple the ground state at q = 1, by a Killing-spinor
witness) and Theorem 4 (Sigma(0) = 0, re-tiered to the model in v6.5.0). Owed computation "P28 vector-photon
vertex" in `docs/qa/group5.debt.md`.
**Guardrail / ledger check:** the sec. 3 rows "sigma-vertex and direction-resolved vector QED on Fock graph
(VQ-1..VQ-5)" and "co-exact mode q-labeling" are about photons on the *graph*. This record is about the *continuum*
vertex on round S^3 (Camporesi-Higuchi spinors, transverse vector harmonics, Clifford multiplication), so it is not
inside their scope. It does not try to put the vector photon on the graph.

## Why

Paper 28's spectral sums couple spinor level n1 to n2 through photon level q under
|n1 - n2| <= q <= n1 + n2 with n1 + n2 + q odd. That is the triangle for scalar harmonics. Theorem 4's ground-state
zero, and the parity discriminant, follow from it. Paper 33's 1 + 6 + 1 partition scores eight selection rules in
the same model, and its "Dirac spinor phase" rule rests on a diagonal vertex V(a, a) = 0 that the code
(`geovac/vector_qed.py`) imposes by hard-code. Stage A asks what the physical vertex
int_{S^3} <psi_1, A . psi_2> dvol actually allows.

## The rule (derived before any computation)

Spin(4) = SU(2)_L x SU(2)_R acts by isometries, and Clifford multiplication by a 1-form is equivariant, so a vertex
can be nonzero only if rho_1* (x) V (x) rho_2 contains a singlet. Content:
- Spinors at level n (eigenvalue of D = +-(n + 3/2), (n+1)(n+2) states each sign): sign + is ((n+1)/2, n/2), sign -
  is (n/2, (n+1)/2) (labelling of which sign is which is a convention; the rules below are stated for it).
- Transverse (coclosed) vector harmonics at level q >= 1, Hodge eigenvalue (q+1)^2: component A = ((q+1)/2, (q-1)/2)
  and component B = ((q-1)/2, (q+1)/2), q(q+2) states each.
SU(2) is multiplicity-free, so each (sign pair, component) has at most one invariant: one reduced amplitude.
Triangle and integrality in L and R give:

**Same sign (++; -- is the mirror, A <-> B):** n1 + n2 + q odd, and
- component A open for |n1 - n2| + 1 <= q <= n1 + n2 + 1;
- component B open for |n1 - n2| + 1 <= q <= n1 + n2 - 1.
Inside the model's window the count W in {0, 1, 2} agrees with the model; the one new channel is q = n1 + n2 + 1
(one component).

**Opposite sign (+-):** n1 + n2 + q **even**, with the L/R triangles
- A: |n1 + 1 - n2| <= q + 1 <= n1 + n2 + 1 and |n1 - n2 - 1| <= q - 1 <= n1 + n2 + 1;
- B: |n1 + 1 - n2| <= q - 1 <= n1 + n2 + 1 and |n1 - n2 - 1| <= q + 1 <= n1 + n2 + 1.
The model admits no triple with n1 + n2 + q even, so every opposite-sign channel is outside it.

**Ground state (n_ext = 0):** same sign, q = n_int + 1, one component, every n_int >= 0; opposite sign, q = n_int,
component B (for ext sign +), every n_int >= 1. The model has none of these.

The rule is a necessary condition. That every open (sign pair, component, n1, n2, q) has a nonzero reduced amplitude
is the generic expectation, not derived; SA5 tests it.

## Method (fixed before scoring)

- **Driver:** `debug/vector_vertex_s3.py` (written after this record; sha256 recorded in the outcome). Writes
  `debug/data/vector_vertex_s3.json`. Memo `debug/sprint_vector_vertex_s3_memo.md`.
- **Exact arithmetic only.** Zero and nonzero are exact (sympy), so there is no numerical floor to set.
- **Route 1 (group theory):** amplitudes as products of SU(2) Clebsch-Gordan coefficients over the content above,
  with the spinor and vector harmonics built from the left-trivialized frame (functions (j, j) tensored with the
  spin-1/2 or spin-1 frame index). Range n1, n2 <= 3, q <= 7, both signs, both components, all m.
- **Route 2 (coordinates, independent of route 1's CG tables):** S^3 in R^4; spinor eigenfunctions and coclosed
  1-form eigenforms built as polynomials in (x1..x4), checked against D and the Hodge Laplacian directly; vertex
  integrals by exact monomial integration over S^3. Range n1, n2 <= 2, q <= 5.
- **Scored quantity:** for each (sign pair, component, n1, n2, q), whether the vertex is zero for every choice of m
  labels (closed) or nonzero for some choice (open). Plus the diagonal amplitude of SA4.

**Controls (exact; failing any one makes the record VOID, not a miss):**
- C1. D spectrum +-(n + 3/2), multiplicity (n+1)(n+2) per sign, n <= 3 (route 1); explicit D-check n <= 2 (route 2).
- C2. Transverse vector harmonics: coclosed, Hodge eigenvalue (q+1)^2, 2q(q+2) states, q <= 7 (route 1 by
  dimension; route 2 explicitly, q <= 5).
- C3. Route 1 restricted to eq:vertex_rule reproduces `geovac.qed_anomalous_moment.get_photon_channels` channel
  counts on n1, n2 <= 3 (the model is what the code says it is).
- C4. Routes 1 and 2 agree on the open/closed pattern, and on squared amplitudes summed over m, wherever both are
  computed. A disagreement means a construction error and voids the record.

## What I know before comparing (disclosed)

Nothing computed. Paper 28 already asserts, by one witness (a Killing spinor against a left-invariant coframe
1-form), that the ground state couples to itself at q = 1; SA3's n_int = 0 clause is that witness. The other clauses,
the opposite-sign rule and the diagonal amplitude, have not been computed anywhere in the corpus that I know of.

## Predictions

- **SA1 (same-sign rule).** The open same-sign set on route 1's range is exactly the set above. In particular
  q = n1 + n2 + 1 is open (one component) for every n1, n2 <= 3, and nothing outside the stated windows is open.
- **SA2 (opposite-sign rule).** The open opposite-sign set is exactly the set above, every member has n1 + n2 + q
  even, and the set is nonempty (e.g. n1 = 0, n2 = 1, q = 1, component B). So the model's parity clause is not a
  rule of the physical vertex once both signs of D are summed.
- **SA3 (ground state).** For n_ext = 0 the open channels are exactly: same sign at q = n_int + 1 (n_int = 0..3) and
  opposite sign at q = n_int (n_int = 1..3), one component each. Hence Sigma(0) in eq:self_energy with the physical
  channel count is a sum of positive terms with at least one present: nonzero. Theorem 4 does not survive.
- **SA4 (diagonal vertex).** For the ground state a (n = 0, either sign) there is a q = 1 transverse harmonic and an
  m with <a | A . a> != 0. So V(a, a) = 0 is not a property of the physical vertex.
- **SA5 (no accidental zeros).** Every rep-allowed (sign pair, component, n1, n2, q) on route 1's range has a nonzero
  amplitude for some m labels.

Scoring: each prediction HELD or MISSED on its stated set; any mismatch is reported as an extra (a channel open that
the rule closes: the derivation is wrong) or an accidental zero (a channel the rule allows, computed closed).

## Reported, not scored

- Sigma_phys(n_ext), n_ext = 0..3, in eq:self_energy form with the physical channel count and both internal signs,
  at cutoffs n_max = 3, 5, 8, 12 (convergence as observed); and the same with squared reduced amplitudes in place of
  the count.
- For each of Paper 33's eight rules: holds, fails, or not decidable at stage A, with the triple that decides it.
- Normalization, alpha, and the flat-space limit are stage B and are not attempted.
