# VV stage A: the physical vector-photon vertex on S^3 (outcome memo)

**Record:** `docs/predictions/2026-11-07_physical_vector_vertex.md` (commit 974e3d45, before any computation).
**Date:** 2026-11-07. **Verdict:** controls C1-C4 all PASS; **SA1-SA5 all HELD**, no extras, no accidental zeros.

**Artifacts:**
- `debug/vector_vertex_s3.py`, the driver. sha256 `3e1284b16f9332f47ffd21995921c50ee777e8f74681c8ccb881a26182a72ef4`.
- `debug/vector_vertex_s3_route1.py` (sha256 `eee1123d...cad8bf`).
- `debug/vector_vertex_s3_route2.py` (sha256 `daa2296f...dcc1`).
- `debug/data/vector_vertex_s3.json`.

The full run takes about 30 s and uses exact arithmetic throughout.

## Method as executed

**Route 1 (group theory).** S^3 = SU(2) with the left-invariant frame X_i = -2i R_i. Functions are the Wigner D^j_{m m'}; right translation acts on m'.

- **Spinors.** f (x) chi, with D = sum (i sigma_i) X_i + 3/2 = 2 sigma.R + 3/2. The spin connection gives (1/4) eps_ijk gamma_i gamma_j gamma_k = 3/2. The eigenstates are the CG couplings (j, 1/2) -> J.
- **1-forms.** f (x) e_mu, with *d = -2 S.R - 2. The constant is forced by *d(grad) = 0.
- **Vertex.** V = int <psi_1, c(A) psi_2> dg over Haar measure of volume 1, with c(A) = sum a_i sigma_i. It factorizes exactly as
  V = sqrt((2j1+1)(2jA+1)(2j2+1))/(2j1+1) * <jA mA j2 m2|j1 m1> * Rt(M1,MA,M2).
  Here Rt is the frame-side contraction of four CGs and a sigma element, computed by brute force over all m. Numbers live in the multi-quadratic field, so zero means exactly zero.
- **Range.** n1, n2 <= 3 and q <= 7, both signs and both components: 896 (sign pair, component, n1, n2, q) cells.
- **Closed form.** The brute-force sum over m of |V|^2 equals
  (2jA+1)(2j2+1) * 6(2j1+1)(2JA+1)(2J2+1)(2J1+1) {jA j2 j1; 1 1/2 1/2; JA J2 J1}^2
  on all 896 cells (0 mismatches). It is used only to extend the reported Sigma tables beyond route 1's range.

**Route 2 (coordinates in R^4; imports neither route 1 nor geovac).**

- **Spinors.** S^+ = C^2 of R^4 restricted to S^3. The tau_a are (1, i sigma_k), and gamma_a gamma_b on S^+ is tau_a tau_b^dag. The candidate operator is D = sum_{a<b} tau_a tau_b^dag L_ab + c0, with c0 = +-3/2 from the Gauss formula. Clifford multiplication is c(v) = xslash vslash.
- **Eigenspinors.** Exact nullspaces of D - lambda over Q(i) on S^+ (x) H_k for k <= 3, with H_k the harmonic polynomials.
- **1-forms.** Polynomial 4-vectors of degree q, with (*dA)_a = eps_abcd x_b d_c A_d. Eigenforms are the exact nullspaces of *d -+ (q+1) over Q, q <= 5. Tangency and coclosedness are verified afterwards, not imposed.
- **SU(2) x SU(2) content.** Measured by the Casimirs of the two su(2) factors of the total so(4) generators G_ab = L_ab + (1/2) tau_a tau_b^dag (spinors) or L_ab + E_ab (vectors). They are calibrated on the scalar x1, which has content (1/2, 1/2).
- **Vertex.** Each eigenspace gets an exact Gram-Schmidt orthogonal basis. The full tensor V_ijk is computed for every basis vector, with exact monomial moments on S^3 (volume 1). No Schur shortcut is used.
- **Range.** As fixed in the record: n1, n2 <= 2 and q <= 5, giving 360 cells.

**Sign convention (stated as the task asked).** "First" SU(2) factor means the factor in which the D = +(n+3/2) spinors carry the larger spin (n+1)/2. This is the record's convention.

- In route 1, the first factor is the right-translation SU(2), frame index included (label J). The second factor is the left SU(2) (label j).
- In route 2 it is the factor with Casimir index 0. The measured content of every space matched the record's pattern: + = ((n+1)/2, n/2) and - = (n/2, (n+1)/2) for n <= 2. Components A and B are identified by measured content.
- Incidentally, A carries *d = -(q+1) in both routes.

**Mirror.** I checked that swapping the convention (sign -> -sign, A <-> B) maps the record's rule onto itself, for all sign pairs. There were 0 mismatches over 896 cells, including the (-,+) case, which is the transpose of the stated (+,-) rule. Swapping the external sign in the Sigma sums gives identical exact values.

**Rule transcription.** The record's rule as written agrees with the generic triangle-plus-integrality test on the content, cell by cell (0 mismatches).

## Controls

- **C1 (D spectrum): PASS.**
  - Route 1: 4 S.R + 3/2 acting on the coupled states gives exactly +-(n+3/2), with multiplicity (n+1)(n+2), for n <= 3.
  - Route 2: c0 = -3/2 gives the symmetric spectrum on S^+ (x) H_k, k = 0..3. The eigenspace dimensions (2 | 6+2 | 12+6 | 20+12) exhaust 2(k+1)^2. The choice c0 = +3/2 fails, so exactly one constant passes.
  - Every route-2 eigenspinor was checked directly against D.
  - D commutes with all six so(4) generators on every space.
- **C2 (transverse vector harmonics): PASS.**
  - Route 1: -2 S.R - 2 gives -+(q+1) on J = j +- 1, with q(q+2) states each, q <= 7. It gives 0 on J = j, the exact forms.
  - Route 2: q(q+2) eigenforms per sign for q <= 5. All are tangent and divergence-free, and (*d)^2 = (q+1)^2 holds. *d also annihilates tangential gradients, which checks the formula independently of orientation.
- **C3 (model = code): PASS.** On the 28 vertex_rule triples with n1, n2 <= 3, the route-1 same-sign open-component count equals `len(get_photon_channels(...))`, called with the content convention `qed_anomalous_moment` uses (jL = (n+1)/2, jR = n/2 for +; mirrored for --). Also, `vertex_allowed` is false on every triple outside the model window.
- **C4 (routes agree): PASS.** All 360 common cells agree on open/closed, and on the exact sum over m of |V|^2 (0 mismatches). There are 70 open cells, including the new channels. Examples:
  - ++|A|0,0,1: S = 6.
  - +-|B|0,1,1: S = 6.
  - ++|A|1,1,1: S = 12.
  - ++|B|1,1,1: S = 2.
  - +-|A|1,1,2: S = 16/3.

## Scoring (route 1's range, as written)

| Prediction | Verdict | Mismatches |
|:--|:--:|:--|
| SA1 same-sign rule | **HELD** | 0 extras, 0 accidental zeros (++ and --). q = n1+n2+1 has exactly one open component for every n1, n2 <= 3, both signs |
| SA2 opposite-sign rule | **HELD** | (+,-): 0 extras, 0 accidental zeros. Every open opposite-sign cell has n1+n2+q even. The set is nonempty, and the witness +-\|B\|0,1,1 is open. (-,+) scored against the transposed rule (reported only): 0 / 0 |
| SA3 ground state | **HELD** | 0 mismatches, with n_ext in either slot and both external signs. Open channels: same sign at q = n_int+1, opposite sign at q = n_int (n_int >= 1), one component each. So Sigma(0) > 0 with the physical count: Theorem 4 does not survive for the physical vertex |
| SA4 diagonal vertex | **HELD** | Route 1: for each n = 0 state a, the sum over photon m of \|<a\|c(A)\|a>\|^2 at q = 1 equals 1. The open component is A for D=+ and B for D=- (the other component gives 0). Route 2 gives 1 for every basis state of both signs |
| SA5 no accidental zeros | **HELD** | Every rule-allowed cell in all 896 has a nonzero amplitude |

Census: 168 of 896 cells are open: 44 ++, 44 --, 40 +-, 40 -+. Of the 88 same-sign open cells, 32 lie outside the model's window, all at q = n1+n2+1. Inside the window the counts equal the model's W (that is C3).

Also verified in route 1: every nonzero Rt satisfies M1 = MA + M2, and the left CG conserves m_L. So m is conserved in both factors separately.

## Reported, not scored

### Sigma_phys(n_ext)

Conventions:
- External sign +, psi_1 = external. Both internal signs are summed, with the per-sign degeneracy g_s = (n+1)(n+2).
- q runs over everything that is open, which means q <= n_ext + n_int + 1.
- "count" is the eq:self_energy form with W -> number of open components, and d_q^T = q(q+2).
- "amp" replaces W g d by the sum over m_int and m_gamma of |V|^2 at a fixed external state. With Haar volume 1 and unnormalized coupling, this is the true unnormalized one-loop sum.
- mu_q = q(q+2) follows Paper 28 eq:hodge1_spectrum. The variant mu = (q+1)^2 is the Hodge eigenvalue of the coclosed family.
- Exact rationals are in the JSON.

| n_ext | n_max | count, mu=q(q+2) | count, mu=(q+1)^2 | amp, mu=q(q+2) | amp, mu=(q+1)^2 |
|:-:|:-:|--:|--:|--:|--:|
| 0 | 3 | 0.959741 | 0.786627 | 0.266936 | 0.206902 |
| 0 | 5 | 1.072368 | 0.895897 | 0.270242 | 0.210107 |
| 0 | 8 | 1.157450 | 0.979710 | 0.271503 | 0.211348 |
| 0 | 12 | 1.214382 | 1.036219 | 0.271924 | 0.211766 |
| 1 | 3 | 2.330562 | 1.999638 | 0.188676 | 0.156053 |
| 1 | 5 | 2.668441 | 2.326779 | 0.192115 | 0.159364 |
| 1 | 8 | 2.923689 | 2.578093 | 0.193401 | 0.160628 |
| 1 | 12 | 3.094486 | 2.747601 | 0.193827 | 0.161050 |
| 2 | 3 | 2.919154 | 2.639586 | 0.111840 | 0.099594 |
| 2 | 5 | 3.482286 | 3.182037 | 0.115517 | 0.103077 |
| 2 | 8 | 3.907700 | 3.600438 | 0.116846 | 0.104379 |
| 2 | 12 | 4.192361 | 3.882883 | 0.117278 | 0.104807 |
| 3 | 3 | 3.145440 | 2.941209 | 0.070350 | 0.065355 |
| 3 | 5 | 3.933825 | 3.691546 | 0.074454 | 0.069067 |
| 3 | 8 | 4.529404 | 4.276156 | 0.075848 | 0.070424 |
| 3 | 12 | 4.927929 | 4.671420 | 0.076290 | 0.070862 |

Observed convergence:
- The count version converges slowly, roughly like 1/n_max: the terms go as g/lambda^4 ~ 1/n^2 with a bounded number of open q.
- The amplitude version converges fast.
- Sigma_phys(0) is nonzero at every cutoff.
- The amplitude-weighted Sigma decreases with n_ext. That is the opposite of the "increases monotonically" statement Paper 28 makes about the model's count. The count version does increase.

### Paper 33's eight rules against the physical vertex

| # | Rule | Stage-A verdict | Deciding cell |
|:-:|:--|:--|:--|
| 1 | Gaunt/CG triangle sparsity | **holds** in its Spin(4) form: every open cell satisfies both SU(2) triangles, with 0 extras in 896 | e.g. ++\|*\|0,3,1 closed (both components) |
| 2 | Vertex parity / GS structural zero | **fails** | ++\|A\|0,0,1 open (S = 6). Also +-\|B\|0,1,1 is open with n1+n2+q even |
| 3 | SO(4) channel count W vanishes on disallowed cells | **holds** as the representation-theory statement (0 extras, 0 accidental zeros). **Fails** as the model's W table: W_phys = 1 where W_model = 0 | ++\|A\|0,0,1; every q = n1+n2+1 |
| 4 | Delta m conservation | **holds**, in both factors separately (checked on every nonzero Rt) | none fails |
| 5 | Spatial parity (E1: same-parity states decouple, photon parity -1) | **fails** in the literal form: a state couples to itself (SA4). The antipodal element of Spin(4) gives photon parity (-1)^(q+1), not -1. Its selection content is the triangle integrality already in rules 1/8 | ++\|A\|0,0,1 diagonal, sum over m_A of \|V(a,A,a)\|^2 = 1 |
| 6 | Ward identity | **holds**, in Ward-Takahashi form: <u, c(df) v> = (lambda_u - lambda_v)<u, f v> exactly. Checked on 384 element pairs (f harmonic of degree 1, 2; n <= 1), 56 with nonzero LHS. Pure-gauge photons decouple between degenerate states | operator identity [D, f] = c(df) |
| 7 | Furry (as scored in Paper 33: diagonal alpha_z vertex V(a,a) = 0) | the **scoring mechanism fails**: V(a,a) is nonzero (SA4). Furry's theorem proper (C-odd loops) is **undecidable at stage A**: it needs C on S^3 spinors and three-vertex traces | ++\|A\|0,0,1 diagonal |
| 8 | SO(4) triangle \|n1-n2\| <= n_gamma <= n1+n2-2 | **fails** | ++\|A\|0,0,1 (q = 1 > n1+n2-2 = -2); every q = n1+n2+1 |

## Findings outside the record (flag to PM/PI)

1. **The W that feeds eq:self_energy in code is not the model's channel count.** `geovac/qed_self_energy._so4_channel_count` couples content ((n1+1)/2, n1/2) to (n2/2, (n2+1)/2) (mixed signs) and checks the triangle inequalities without the integrality clause.
   - Inside the model window with n <= 3, it disagrees with `get_photon_channels` and with the physical same-sign count at 9 triples. For example, (1,1,1): code W = 0 against physical 2; (1,2,2): code W = 1 against physical 2.
   - The units it does count are representation-theoretically empty, because the triangle sums are half-integer. For example, in (1,2,2) both triangles of component B have half-integer sums.
   - With the integrality check added, it would give W = 0 on every vertex_rule triple, because opposite-sign cells need n1+n2+q even.
   - So "the model is what the code says it is" holds for `get_photon_channels` (C3 passed), but `qed_self_energy.py` implements a third object.
   - Theorem 4's Sigma(0) = 0 is unaffected: it rests on the window alone. But any printed model Sigma(n_ext >= 1) values inherit this W.
2. **Paper 28's mu_q = q(q+2) is not the Hodge eigenvalue (q+1)^2 of the family the photon lives in.** Paper 28 already flags this labelling issue. Both variants are reported above.

## What in the record I think is mis-posed

- **The opposite-sign rule is stated only for (+,-).** (-,+) is its transpose, which also equals the mirror. SA2 was scored on (+,-) as written; (-,+) is reported, also 0/0.
- **SA3 does not say which vertex slot n_ext occupies.** It was scored in both slots; there are no mismatches either way.
- **C3 does not fix the content passed to `get_photon_channels`.** I used the convention `qed_anomalous_moment` itself uses. The other module (finding 1) would not pass, so C3's "the model is what the code says" is true only for that function.
- **"Sigma_phys ... with the physical channel count and both internal signs" leaves the per-sign degeneracy and mu_q open.** I used (n+1)(n+2) per sign and both mu conventions.
- **SA4's "an m" holds only with the components summed.** For a D=- ground state, component A gives zero and only B opens.

None of these changes a verdict.

## What stage B needs

- **Coupling normalization.** Fix e and the photon field normalization: the transverse photon propagator on S^3 in a stated gauge, with the coexact-mode measure, and the physical Hodge eigenvalue (q+1)^2 rather than q(q+2). Connect the Haar-volume-1 amplitudes to the radius R.
- **Reduced amplitudes.** The 9j closed form above gives them for all levels. Stage B needs their large-n asymptotics, to define Sigma with regularization (the count version converges only like 1/n_max) and to take the flat-space limit R -> infinity at fixed physical momentum.
- **Gauge sector.** Include the exact and longitudinal modes and the ghosts (Feynman gauge). The Ward identity checked here is the input to that.
- **Spin structure for Furry.** A charge-conjugation operator on S^3 spinors, to decide rule 7 properly.
- **Paper updates (PM/PI).** The SA verdicts here would re-price Paper 28 Theorem 4's physical reading and Paper 33's rules 2, 3, 5, 7, 8 for the physical vertex. No paper was edited.
