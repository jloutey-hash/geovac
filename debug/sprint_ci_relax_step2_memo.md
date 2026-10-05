# LiH residual, Step 2: relaxing the CI coefficients alongside the geminal (diagnostic, 2026-10-05)

**Verdict: BORDERLINE.** Relaxing the CI coefficients alongside the e^{-1.5 r12} geminal (ansatz A)
gains **2.7 +/- 0.1 mHa** (LM, correlated) / **1.9 +/- 0.4 mHa** (independent re-sampled VMC) on the
full 90-determinant pi reference. A geminal-dressed copy of every determinant (ansatz B) gains
**3.8 +/- 0.15 / 2.9 +/- 0.4 mHa**. Both are real (>> 3 sigma) and both sit in the 2-5 mHa band. They
close about 10-18% of the 21.5 mHa residual. Frozen CI coefficients are a contributor, not the main
lever. Best variational energy here is (B) = **-8.0524 +/- 0.0003 Ha**, i.e. 18.1 mHa above exact (-8.0705).

Scripts: `debug/ci_relax_step2_lm.py` (LM sampler + re-sampling VMC) and `debug/ci_relax_step2_analyze.py`
(correlated generalized eigenproblem + jackknife). Data: `debug/data/cirelax_*.{npz,log}`.
Nothing in geovac/, tests/, papers, CHANGELOG, CLAUDE.md or docs/ was touched.

## 1. Method

The basis is the 2N functions {Phi_I, F Phi_I}, where Phi_I = dA[a_I] dB[b_I] runs over the kept
determinants of the truncated M=16 pi natural-orbital reference, and F = sum_{i<j} exp(-1.5 r_ij).
All three ansatzes are subspaces of this basis, so ONE set of sampled matrices serves all of them
(exact correlated sampling).

- (C) span{Psi_CI, F Psi_CI}. This is the frozen 2x2. (F - Fbar)Psi_CI spans the same space.
- (A) span{Phi_I, F Psi_CI}, N+1 parameters.
- (B) span{Phi_I, F Phi_I}, 2N parameters.

Sampling density: P = |g Psi_CI|^2 + eps sum_I |Phi_I|^2. Here g = 1 + c(F - Fbar) is the frozen
marriage multiplier (c = -0.2269, Fbar = 0.4560, from `lih_pi_m16_energy_exp1.5.npz`). The eps term
carries 5% of the mass. It bounds every ratio X_k/sqrt(P), so S and H have finite variance at the
nodes of Psi_CI and at the zeros of g. Estimators are the standard non-symmetric LM ones:
S_kl = <X_k* X_l / P> and H_kl = <X_k* (H X_l) / P>. The local action uses
H(F Phi) = F H Phi - 1/2 (sum_i lap_i F) Phi - grad F . grad Phi, with analytic orbital Laplacians.

The eigenproblem is solved by canonical orthogonalisation of S, then eig of the projected H, taking
the lowest real root. Errors come from a jackknife over 50 block groups. The optimum vectors were then
re-sampled with plain VMC on |Psi_opt|^2 using the standard local energy, fresh seeds, and burn-in
>= 3000 sweeps (9-det) / 4000 sweeps (90-det). The corpus lesson is that 90 determinants need >= 4000
burn-in sweeps. Threads were pinned to 1.

**Instrument checks** (before any production run):
- The local energy of Psi_CI reproduces `lih_vmc.local_energy_analytic` to 2e-14.
- The local energy of (1 + cF)Psi_CI reproduces `lih_vmc.local_energy_linear` to 2e-14.
- H(F Phi_I) matches a finite-difference Laplacian to a relative error of 2e-8.
- **LM half-sample check:** the two halves of the 90-det seeds give A-C = -2.68 / -2.65 mHa and
  B-C = -3.88 / -3.80 mHa. No sample-size bias is visible at 180 parameters.

## 2. Results

Statistics:
- 9-det LM: 12 seeds x 1000 walkers x 3000 production sweeps (3000 burn-in), thin 5.
- 90-det LM: 10 seeds x 1000 walkers x 5000 production sweeps (4000 burn-in), thin 5.
- Re-sampled VMC: 2 seeds per vector (9-det, 8000 sweeps total, 3000 burn-in) or 4 seeds per vector
  (90-det, 10000 sweeps total, 4000 burn-in), 1000 walkers each.

| reference | quantity | LM (correlated) | re-sampled VMC (independent) |
|---|---|---|---|
| 9-det (trunc 1e-2) | E(bare Psi_CI) | -8.02525 +/- 0.00067 (E_trunc -8.025184) | -- |
| | **E(C)** frozen 2x2 | -8.04817 +/- 0.00030 | -8.04820 +/- 0.00056 |
| | **E(A)** CI free | -8.05047 +/- 0.00024 | -8.05066 +/- 0.00039 |
| | **E(B)** dressed dets | -8.05067 +/- 0.00025 | -8.05042 +/- 0.00044 |
| | E(A)-E(C) | **-2.31 +/- 0.12 mHa** | -2.46 +/- 0.68 mHa |
| | E(B)-E(C) | -2.51 +/- 0.14 mHa | -2.23 +/- 0.71 mHa |
| | E(B)-E(A) | -0.20 +/- 0.05 mHa | -- |
| 90-det (trunc 1e-3) | E(bare Psi_CI) | -8.02874 +/- 0.00045 (E_trunc -8.028633) | -- |
| | **E(C)** frozen 2x2 | -8.04897 +/- 0.00022 | -8.04950 +/- 0.00030 |
| | **E(A)** CI free | -8.05164 +/- 0.00017 | -8.05144 +/- 0.00025 |
| | **E(B)** dressed dets | -8.05281 +/- 0.00016 | -8.05239 +/- 0.00027 |
| | E(A)-E(C) | **-2.67 +/- 0.12 mHa** | -1.94 +/- 0.39 mHa (-2.40 vs analytic C) |
| | E(B)-E(C) | **-3.84 +/- 0.15 mHa** | -2.88 +/- 0.40 mHa (-3.34 vs analytic C) |
| | E(B)-E(A) | -1.17 +/- 0.07 mHa | -0.94 +/- 0.37 mHa |

The LM column is a finite-sample eigenvalue of the projected problem. The re-sampled column is an
unbiased variational energy of that fixed vector. The re-sampled values are slightly higher, by
0.2-0.4 mHa for A and B, which is the usual LM optimism plus noise. The re-sampled (C) sits 0.46 mHa
below the analytic value, a 1.5 sigma deviation, so the re-sampled gains are probably slightly understated.
Taking all estimates together: **(A) is worth 2.0-2.7 mHa, (B) is worth 2.9-3.8 mHa.**

**What the relaxation does.** The relaxed (A) coefficients are scaled-down copies of the frozen ones
(90-det: double 0.1307 -> 0.1199, 0.0649 -> 0.0588; 9-det: 0.0313 -> 0.0194). The overlap with the
reference is |<d_A|c_ref>| = 0.9998, and the geminal coefficient c grows from -0.205 to -0.229. This is
the textbook double-counting correction: the geminal already carries part of the pair correlation
that the CI doubles were mimicking, so the CI de-emphasises them. The 9-det gain (2.3) is close to the
90-det gain (2.7), which means the effect is concentrated in the few leading dets. (B) adds a further
1.2 mHa on 90 dets (0.2 on 9). That increment is the determinant-dependent geminal amplitude.

## 3. Sanity (step 3): PASS

- 9-det (C): LM -8.048166 +/- 0.30 mHa and re-sampled -8.04820 +/- 0.56 mHa, against analytic
  -8.048141. The deviations are -0.03 and -0.06 mHa.
- 90-det (C): LM -8.048972 +/- 0.22 mHa against analytic -8.049045 (+0.07 mHa, 0.3 sigma) and the
  VMC twin -8.04921 +/- 0.14 (+0.9 sigma). The re-sampled value -8.04950 +/- 0.30 is 1.5 sigma below
  the analytic value.
- The bare references reproduce E_trunc within 0.1-0.2 sigma.

The instrument is calibrated.

## 4. What an analytic (A) needs from the RI-free engine (sketch, not built)

The (N+1)-dimensional problem in {Phi_I} ∪ {(F - Fbar)Psi_CI} needs four blocks:

1. **<Phi_I|H|Phi_J>.** Available: this is the M=16 CI Hamiltonian (h1/eri in the NO basis), with S = 1.
2. **<FPsi|H|FPsi>, <FPsi|FPsi>.** Available: these are the existing frozen 2x2 g and sigma^2 pieces
   (with Fbar). This is the expensive O(nnz) g_Vee part, and it is computed once.
3. **S_{I,F} = <Phi_I|(F - Fbar)|Psi_CI>.** This is a Fbar-type one-f-edge moment with a
   *transition* coefficient tensor T_I = (1/4) sum e_I[a,b] CS[a',b'] E E in place of the
   reference's own CS (x) CS. `pi_m16`'s machinery separates cleanly here:
   - `contract_diagram_m` builds per-electron factor tensors over the active ordered pairs, and
     these do not depend on T.
   - `RefM.contract` / `traces_coo` contract them against T at O(nnz) cost.
   - Every pair (a_I, a') generated by a kept determinant against the reference is already in the
     reference's active ordered-pair set, because T(CS(x)CS) already contains all cross-det pairs.

   So S_{I,F} for all I = the existing diagram factors + N cheap re-contractions with transition COO
   tensors. **Missing:** a `build_T` variant taking two different coefficient matrices (bra e_I,
   ket CS). Its output is not bra/ket symmetric, so the ordered-pair COO must not be symmetrised.
   This is a few lines of code.
4. **H_{I,F} = <Phi_I|H(F - Fbar)|Psi_CI>** (the h-type transition elements).
   - **V_ne and V_ee parts:** the same diagrams as the existing h_Vne / h_Vee, with one f edge and
     at most one Coulomb edge. They involve only chains and disjoint pairs, never the triangle or
     the 4-body bridge (those enter only g). Each is a re-contraction with T_I. Available modulo the
     transition-T builder.
   - **Kinetic part: genuinely new.** The existing h_T uses an IBP collapse (Part B -> scalar Yukawa
     Ybar), which relies on bra = ket. For a transition, -<Phi_I|grad F . grad Psi> keeps an
     antisymmetric current term (phi_p grad phi_q - phi_q grad phi_p). One clean route stays scalar.
     Write H(F Psi) = F (H Psi) - 1/2 (lap F) Psi - grad F . grad Psi and integrate the last term by
     parts onto the bra side. That gives <(grad Phi_I . grad Psi + Phi_I lap Psi) F>.
     - This needs the existing gradient base (G_pq, 'grad') plus a **new ordered-pair base
       L_pq = phi_p lap phi_q**. It is not symmetric in p,q. The orbital Laplacians are already on
       the grid (`pi_m16` `nlap`, used for the T gate) but are not exposed as a RefM base.
     - The Yukawa term -1/2 <Phi_I (lap F) Psi> is the existing kernel with T_I.

   Equivalently, the kinetic (F H Psi) piece can be assembled from the 'rhov', 'lap' and Coulomb bases.

**Cost.** All factor tensors are T-independent, so the whole transition column is roughly ONE extra
h-type evaluation (no g) plus N O(nnz) contractions. At 90 dets that is minutes, against the ~2.5 h
g build. Analytic (B) is different: it needs <Phi_I|F H F|Phi_J> for N^2 pairs. These are g-type
elements with transition tensors on both sides, i.e. the g_Vee triangle/bridge machinery re-contracted
N^2 times with non-symmetric T_IJ. That is feasible in principle, since the factors are again
T-independent, but it is about N^2 = 8100 O(nnz) contractions of the dense TREE4 fields. Expensive,
not blocked.

## 5. Ledger check

`docs/failed_approaches_ledger.md` has no prior CI-coefficient-relaxation attempt for the LiH
marriage. The only mention is step 1's own (uncommitted) note naming frozen coefficients as a
candidate (L337). Related but distinct prior rows:
- "Multi-geminal expansion / quadratic F^2 term for the LiH r12 marriage gap" (2026-09-26): extra
  geminals saturate at +1.4 mHa, F^2 is worth 0.013 mHa. That linear-method scan also kept Psi_CI frozen.
- Step 1 (r12 t^2, STOP): 0.3-3 mHa.

## 6. Reading: what the residual is

After (B), the best variational energy is -8.0524, still **18.1 mHa** above exact. The tested levers
on the pi reference are:

| lever | worth |
|---|---|
| more geminals | 1.4 mHa |
| F^2 | 0.01 mHa |
| r12 t^2 / e-e-n | 0.3-3 mHa |
| CI relaxation | 2-3.8 mHa |

Each is small. None of them changes the one-particle space. The most plausible dominant remainder
is **the orbital basis together with the single short-range linear-geminal form**. The evidence:
1. The M=16 FCI is 41.4 mHa above exact. Enlarging the orbital set saturates at an orbital ceiling of
   -8.032, which is the radial/angular limit of *this* primitive family. In the geminal-free limit,
   the full orbital lever is only ~3 mHa beyond M=16.
2. The geminal dresses this fixed orbital space. Relaxation shows the geminal and the CI doubles
   overlap: the doubles shrink by 10-40%. What neither can supply is correlation outside both spans.
   This is mainly the Li 1s^2 core-core radial/angular correlation (the additive-F12 core estimate is
   30 mHa for an isolated Li+ cusp, of which the marriage plus relaxation captures only part) and
   core-valence correlation needing tight polarization functions that the M=16 core (two extra 1s
   exponents) does not carry.
3. A single gamma = 1.5 serves the core pair (natural range ~1/Z) and the valence pair alike. The
   9-det -> 90-det comparison shows relaxation saturating, not growing, with reference size.

**Next cheapest discriminator (suggested, not run).** Repeat (B) in the same LM harness with a
two-range basis {Phi_I, F_{1.5} Phi_I, F_{g_core} Phi_I}, plus core-only/valence-only split
geminals. It is the same sampler with a few extra basis columns. If that also stalls at ~1-2 mHa,
the residual is the orbital basis (core-correlating functions), and the move is a core-polarization
enrichment of the primitive set, not more r12 machinery.

## 7. Summary numbers

- Full 90-det pi reference, gamma = 1.5:
  - E(C) = -8.04897 +/- 0.00022 (LM) / -8.04950 +/- 0.00030 (re-sampled); analytic -8.049045.
  - E(A) = -8.05164 +/- 0.00017 / -8.05144 +/- 0.00025.
  - E(B) = -8.05281 +/- 0.00016 / -8.05239 +/- 0.00027.
- Gains: A-C = 2.67 +/- 0.12 mHa (correlated) and B-C = 3.84 +/- 0.15 mHa (correlated); the
  independent re-sampled values are 1.9 +/- 0.4 and 2.9 +/- 0.4.
- Verdict: **BORDERLINE (2-5 mHa).** An analytic (A) is cheap: one transition-T builder plus one new
  ordered-pair Laplacian base. It is worth ~2.5 mHa on the corpus-best energy, taking -8.049045 to
  about -8.0515, and is variational. It does not close the 21.5 mHa.
