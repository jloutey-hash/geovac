# Paper 60 [OPEN] item: is there a preconditioner for the sigma-basis SW metric? (2026-10-09)

**Verdict: STOP for a preconditioner of the s-sector kind. BORDERLINE for the one lever that does work, which removes directions instead of preconditioning them.**
The mechanism is exact. On the sigma basis the cross-centre block C is a finite section of the
translation operator, and that operator is *unitary* on the axial (m = 0) sector in the SW metric.
So the limiting two-centre Gram operator is [[I, U], [U*, I]], whose spectrum is {0, 2}, each with
infinite multiplicity. The metric is the Gram matrix of *two complete orthonormal bases of the same
space*, so its near-null directions are genuine redundancy and grow with the basis. A symbol zero
is not the cause.
Consequences:
- A centre-block-diagonal preconditioner provably gains nothing.
- Structured per-l, per-shell and band-Toeplitz preconditioners gain constant factors at best.
- Fixed-rank deflation does not bound the conditioning.
- Canonical orthogonalisation (Avery's "phase out the smallest") holds cond <= 1e6 with no
  measurable H2+ accuracy cost through N = 8 (36 functions per centre). At any *fixed* threshold,
  though, convergence stalls once N exceeds a threshold-dependent value.

Guardrail note: this work evaluates and conditions a two-centre Sturmian molecular metric. It does
not propose a single-centre encoding, so the Papers 8–9 guardrail does not apply. One point is
adjacent and flagged rather than pursued. The redundancy reading says the one-centre axial set is
complete in this metric, so in the limit the B functions add nothing that A cannot represent. At
finite N, however, the B functions carry exactly the content A converges to only slowly. This memo
makes no claim about single-centre expansions.

Drivers: `debug/p60_sigma_preconditioner_diag.py`, `debug/p60_sigma_preconditioner_n8.py`.
Data: `debug/data/p60_sigma_preconditioner{,_n8}.json` (+ `.log`). Builder: tracked
`geovac/sturmian_sigma_basis.py` (k = 1, R = kR; 160x160 prolate quadrature). Nothing outside
`debug/` was modified.

---

## 1. Mechanism (exact, with measured confirmation)

**The metric is a global inner product.** The Sturmian equation turns the SW metric entry
`S_{t't} = (n_t/k) <chi_t'|1/r_{a(t)}|chi_t>` into `<chi_t'|(p^2 + k^2)/(2k^2)|chi_t>`. This is one
inner product, independent of which centre a function sits on. In it:

1. **Each centre's set is orthonormal and complete.** Orthonormality is the intra-centre identity,
   measured to 2e-12. Completeness of the fixed-k Coulomb-Sturmian set is standard (Fock / Rotenberg /
   Avery), and restricting to m = 0 gives the axial sector.
2. **The B set is the image of the A set under translation T_R.** T_R commutes with p^2, so it is
   unitary in this inner product. With R along z it also preserves m = 0. Hence `C_ab = <chi_a, T_R chi_b>`
   is the N-truncation of a **unitary** operator U on the axial sector.
3. **The gerade/ungerade split.** Let `Pi = diag((-1)^l)` (parity about A). Then `C Pi` is the
   matrix of `T_R P_A`, which is inversion about the bond midpoint: an involution J with J^2 = 1 and
   spectrum {+1, -1}. The gerade and ungerade blocks are exactly `I + C Pi` and `I - C Pi`. In the
   limit they become `I ± J`, which are twice the projectors onto J's eigenspaces. Each has an
   **infinite-dimensional kernel**.

**Contrast with the s-sector.** The l = 0 subspace is *not* invariant under T_R. The compression of
U to it is the multiplication operator with symbol j0(s cot(chi/2)), and |j0| = 1 only at p = 0.
That gives one soft direction with a symbol zero of known order, which tri(1,2,1) matches. On the
sigma basis the limit operator itself is singular on half the space. The finite S_N therefore
approaches an operator with an infinite-dimensional kernel, and its soft eigenvalues are the
truncation's approach to that kernel. Their rate is exponential in the number of shells, the way a
concentration (Slepian-type) compression approaches 0 and 1. KMS behaviour (polynomial) does not
appear.

**Measured confirmation** (`spectra` block of the data):

| kR | N (fns/centre) | frac sigma_k(C) > 0.9 | #eig(S) < 1e-2 (= #eig > 2 - 1e-2) | #mu(C Pi) < -0.99 / > +0.99 |
|---|---|---|---|---|
| 1 | 2 (3) → 7 (28) | 0.33, 0.50, 0.60, 0.67, 0.71, 0.75 | 1, 3, 6, 9, 14, 19 | 0/1 → 8/11 |
| 3 | 2 (3) → 7 (28) | 0.00, 0.33, 0.40, 0.47, 0.52, 0.57 | 0, 1, 3, 5, 8, 13 | 0/0 → 5/8 |
| 8 | 2 (3) → 7 (28) | 0.00, 0.00, 0.10, 0.13, 0.19, 0.25 | 0, 0, 0, 1, 3, 5 | 0/0 → 2/3 |

- The fraction of singular values near 1 rises monotonically toward 1, as the unitary limit predicts.
- The spectrum of S is exactly symmetric about 1 (spec = 1 ± sigma_k), so half of S's eigenvalues go soft.
- `C Pi` is symmetric to 2e-16, and its eigenvalues pile up at **both** −1 (gerade-soft) and +1
  (ungerade-soft), as the involution predicts.
- Involution test on the {1s, 2s, 2p} block at kR = 3: `||((C Pi)^2)_low − I||_max` =
  0.73, 0.63, 0.50, 0.45, 0.39, 0.35, 0.31, 0.27 for N = 2..9. It decreases monotonically and slowly.
  The slow approach is expected: strong convergence needs J chi (an off-centre function) to be
  represented in the A basis, which is the same slow single-centre convergence that drives the
  Bessel deficit to zero.

## 2. Raw conditioning and the precision floor

| kR | N=2 (3) | 3 (6) | 4 (10) | 5 (15) | 6 (21) | 7 (28) |
|---|---|---|---|---|---|---|
| 1: cond S | 5.1e3 | 1.7e8 | 3.6e13 | floor | floor | floor |
| 1: cond gerade | 10.0 | 1.1e4 | 4.3e8 | 1.0e14 | floor | floor |
| 3: cond S | 18.3 | 5.1e3 | 1.2e7 | 1.1e11 | >2e11 (floor) | floor |
| 3: cond gerade | 2.3 | 21 | 7.4e3 | 2.1e7 | 2.1e11 | floor |
| 8: cond S | 1.23 | 4.1 | 76 | 9.2e3 | 3.6e6 | 3.3e9 |
| 8: cond gerade | 1.11 | 1.70 | 3.8 | 53 | 6.9e3 | 3.1e6 |

**Precision floor.** Matrix entries are quadrature-limited at about 2.5e-12. Float64 eigensolver
error, about eps·||S|| ≈ 4e-16, is far smaller. To locate the floor, the lowest three eigenvalues
were recomputed with a different quadrature (176 x 320) and compared:

| eigenvalue size | agreement between the two quadratures |
|---|---|
| ~1e-8 to 1e-11 | 5–6 digits (e.g. kR=3 N=5: 1.808788e-11 vs 1.808778e-11) |
| ~1e-14 | 1% (kR=1 N=4: 5.504e-14 vs 5.472e-14) |
| ~1e-15 | noise (kR=3 N=6: 8.1e-16 vs 1.5e-15) |

So cond values up to ~1e13 are reliable, and anything above ~1e14 is marked "floor". Negative
"cond" entries in the JSON are floor artefacts. The kR=3, N=6 second eigenvalue, 9.50e-12, is
stable, so the true cond there exceeds 2e11. Raising the Gauss–Laguerre order further is not
possible with `laggauss`, because the weights overflow exp(t) beyond ~178 nodes. No conclusion
below rests on a value under the floor: every lever comparison is also made at N ≤ 5 (kR = 1, 3)
and N ≤ 7 (kR = 8), where the numbers are clean.

Growth is exponential in basis size. log10 cond(S) rises about 0.8 per added function per centre
at kR = 3 and about 0.42 at kR = 8. The per-shell factor itself increases with N (kR = 8: ×3.3, ×19,
×120, ×390, ×910).

## 3. Soft-subspace characterisation (gerade block, kR = 3, eigenvalues < 1e-2)

| N | #soft | l-weight of soft vectors (l=0,1,2,3,…) | uniform l-share | max single-l purity per vector |
|---|---|---|---|---|
| 4 | 1 | .24 .55 .19 .01 | .40 .30 .20 .10 | .55 |
| 5 | 2 | .28 .50 .18 .04 .01 | .33 .27 .20 .13 .07 | .64, .36 |
| 6 | 3 | .16 .39 .18 .18 .08 .01 | .29 .24 .19 .14 .10 .05 | .70, .44, .49 |
| 7 | 5 | .16 .33 .19 .21 .09 .02 .00 | .25 .21 .18 .14 .11 .07 .04 | .75, .52, .57, .38, .35 |

- **There is not one soft direction per l.** Every soft vector mixes several l (purity 0.35–0.75).
  The weight is tilted toward l = 1 and toward the outer shells (n-weight at N = 7: n = 3..6 carry
  0.19/0.21/0.24/0.18).
- **Each fixed-l diagonal block of the gerade metric is well-conditioned.** The worst single-l
  block has cond 2.7 → 7.7 over N = 4..7 at kR = 3, and ≤ 53 even at kR = 1, N = 7. All the
  softness is **inter-l**.
- This is what the mechanism predicts. Midpoint inversion re-expanded about A couples every l, so
  the kernel of `I ± J` is not aligned with any (n, l) or multipole channel.

## 4. Levers tested (cond after preconditioning; whitening-invariant)

**(a) Centre-block-diagonal preconditioners: no gain, provably.**
Lemma: for any SPD `P = diag(P_A, P_B)`, `cond(P^{-1}S) ≥ cond(S)`.
Proof: take `z± = (±u_k, v_k)`, where (u_k, v_k) is the top singular pair of C. Both vectors have the
same P-norm, `u^T P_A u + v^T P_B v`. Their S-norms are 2(1 ± sigma_max). The ratio of the two
Rayleigh quotients is therefore (1 + sigma_max)/(1 − sigma_max) = cond(S). ∎
This holds for the s-sector too. There, tri(1,2,1) works only because it acts in the g/u frame,
which is not centre-diagonal. Measured: per-l tri(1,2,1) on each centre gives ratios 1.1–2.6 (≥ 1
everywhere), and random block-diagonal SPD gives 1.3–3.6.

**(b) On the gerade and ungerade blocks (kR = 3):**

| N | raw g | l-block Jacobi | n-shell Jacobi | per-l tri(1,2,1) in n | s-sector tri on l=0 only | ideal deflation r=1 | r=3 | r=5 |
|---|---|---|---|---|---|---|---|---|
| 3 | 21 | 16.5 | 15.8 | 52 | 42 | 2.2 | 1.7 | 1.0 |
| 4 | 7.4e3 | 6.4e3 | 6.4e3 | 1.8e4 | 2.0e4 | 7.4 | 2.1 | 1.8 |
| 5 | 2.1e7 | 1.8e7 | 2.1e7 | 5.6e7 | 7.8e7 | 906 | 4.2 | 2.0 |
| 6 | 2.1e11 | 1.8e11 | 2.2e11 | 6.1e11 | 1.0e12 | 1.5e6 | 187 | 3.1 |
| 7 | (floor) | – | – | – | – | 1.0e10 | 1.8e5 | 56 |

The ungerade block behaves the same way. At kR = 3, N = 5: raw 1.1e11; l-Jacobi 4.8e10; n-Jacobi
9.5e10; per-l tri 1.4e11; l = 0 tri 6.7e10. At kR = 1, l-Jacobi gains at most about 10x (u, N = 4:
3.6e13 → 3.5e12). At kR = 8 every structured lever is within 1.5x of raw, or worse.

- **l-block and n-shell Jacobi** are the strongest preconditioners local in l or in n, since each
  uses the exact diagonal blocks. They gain constant factors of ≤ 10x and leave the exponential
  rate untouched. They cannot do more, because by §3 those blocks are already near-identity.
- **Per-l band-Toeplitz tri(1,2,1) in n**, the direct generalisation, and **the s-sector
  preconditioner on l = 0** both make things *worse* by 2–5x. They impose a symbol zero at an
  index location where the sigma-basis metric has none.
- **Ideal fixed-rank deflation** removes the exact lowest r eigenvectors. By Courant–Fischer it is
  the best any rank-r known-subspace deflation can do. At fixed r it still grows exponentially
  (r = 1: 2.2 → 1.0e10). The rank needed to reach cond ≤ 100 grows with the basis: gerade at kR = 3
  needs 0, 0, 1, 2, 4, 5 for N = 2..7, and kR = 1 needs 0, 1, 2, 4, 6, 9. This is the same picture
  as the water Woodbury result, where flattening needed r ∝ basis size.

**(c) What *any* bounded preconditioner would need.**
Suppose `cond(P^{-1/2} S P^{-1/2}) ≤ K`. Then `cP ≤ S ≤ cKP` in Loewner order. Courant–Fischer
then gives `c·λ_j(P) ≤ λ_j(S) ≤ cK·λ_j(P)` for every j, so P must itself contain about half its
spectrum at exponentially small eigenvalues, aligned with S's near-kernel to within K. That
near-kernel is the image of the translation operator, which has dimension ∝ basis size and mixes
every l (§3). A preconditioner whose structure does not already encode T_R — per centre, per l,
per shell, band-Toeplitz — cannot meet this requirement, and the measurements agree. Encoding the
near-kernel means computing the SVD of C, which is canonical orthogonalisation in another form.
This is the STOP.

## 5. Canonical orthogonalisation: the BORDERLINE lever (H2+, R = 2 bohr, fixed-R Galerkin)

Directions with lambda < tau·lambda_max are dropped, and k is solved self-consistently. By
construction cond(kept) ≤ 1/tau. Error is measured against −1.102634214 Ha.

| N (fns total) | tau = 1e-2 | 1e-3 | 1e-4 | 1e-6 | 1e-8 | 1e-12 (≈ undropped) |
|---|---|---|---|---|---|---|
| 3 (12) | 9.12e-5 [1] | 9.12e-5 [1] | 9.12e-5 [0] | 9.12e-5 | 9.12e-5 | 9.12e-5 |
| 4 (20) | 4.64e-5 [3] | 4.64e-5 [2] | 2.06e-5 [1] | 2.06e-5 [1] | 2.06e-5 | 2.06e-5 |
| 5 (30) | 1.87e-4 [6] | 5.44e-6 [4] | 5.44e-6 [4] | 5.44e-6 [2] | 5.44e-6 [1] | 5.44e-6 |
| 6 (42) | 5.31e-4 [9] | 1.62e-6 [8] | 1.62e-6 [7] | 1.61e-6 [5] | 1.51e-6 [3] | 1.51e-6 [1] |
| 7 (56) | 1.32e-6 [13] | 1.32e-6 [12] | 1.10e-6 [11] | 4.49e-7 [8] | 4.49e-7 [6] | 4.49e-7 [3]* |
| 8 (72) | 4.60e-6 [19] | 4.60e-6 [18] | 2.71e-6 [15] | 1.49e-7 [12] | 1.49e-7 [9] | (1e-10: 1.49e-7 [8]) |

[number dropped]. *At or near the quadrature floor.

- **At tau = 1e-6 (cond ≤ 1e6, flat in N by construction) the loss is unmeasurable through N = 8.**
  The error keeps converging at the undropped rate, ≈ 3–3.6x per shell, to 1.49e-7 Ha, while 12 of
  72 directions are dropped. The undropped metric at these N is beyond 1e15.
- **Any fixed tau eventually stalls convergence, and the error can even go non-monotone.**
  tau = 1e-2 stalls at N ≈ 5 (error ~2–5e-4, *worse* than N = 3). tau = 1e-3 and 1e-4 stall at
  N ≈ 7–8 (N = 8 is worse than N = 7). To keep converging, tau must shrink as N grows. On these
  points, cond 1e3 suffices for ~5e-6 Ha and 1e6 for ~1.5e-7 Ha. Two points do not establish a
  law, so no exponent is quoted.
- Applying the lever needs the eigendecomposition (or the SVD of C) classically. For the
  one-electron metric that is cheap. It is the same sparsity-destroying operation already excluded
  as a *structural* lever (locality row, Paper 60 Sec. resource). Canonical orthogonalisation is
  therefore a conditioning remedy, not a preconditioner of the s-sector kind.

## 6. Does the large-R lever survive on sigma?

**Only as a delay.** At fixed N, cond falls steeply with kR. At N = 5: kR = 1 gives floor
(> 1e14), kR = 3 gives 1.1e11, kR = 8 gives 9.2e3. At every kR, though, it still grows
exponentially in N (kR = 8: 1.23 → 3.3e9 over N = 2..7, gerade 1.1 → 3.1e6). Larger separation
postpones the onset, the number of shells needed before A's span reaches into B's. It does not
bound the growth. Under the unitary-limit mechanism this is forced: at any R the limit operator is
[[I, U], [U*, I]]. On the s-sector, by contrast, larger R only lowers the constant of a polynomial
law (1 − sigma_max ≈ (kR)^2 π^2 / (24 n^2)).

## 7. Replacement sentence for Paper 60 (~line 1191)

Replace "whether a preconditioner exists for the $\sigma$ basis is open" with:

> and on the $\sigma$ basis no preconditioner of that kind exists: in the SW metric each center's
> axial Sturmian set is complete and orthonormal and the other center's is its image under a
> translation, which is unitary in that metric, so the two-center metric is a finite section of
> $\begin{psmallmatrix}I&U\\U^\ast&I\end{psmallmatrix}$ with $U$ unitary --- an operator with
> spectrum $\{0,2\}$ of infinite multiplicity, whose near-null directions mix every $\ell$ and
> multiply with the basis.  A center-block-diagonal preconditioner then provably leaves
> $\mathrm{cond}(S)$ unimproved, and per-$\ell$ or per-shell (including band-Toeplitz)
> preconditioners gain only constant factors;\ what does work is Avery and Avery's remedy of
> phasing out the smallest eigenvectors, which at threshold $10^{-6}$ holds the conditioning at
> $10^{6}$ with no measured loss of $\mathrm{H}_2^+$ accuracy through $36$ functions per center,
> although any fixed threshold eventually stalls convergence.

Tier: [MEASURED] for the numbers. The unitary-limit identification and the block-diagonal lemma
are short derivations, and each needs its own backing test before entering the paper (equation
gate, §13.4a). The guard should be written as a separate pass. Suggested guard: for N = 2..6 at
kR = 3, the fraction of sigma_k > 0.9 increases monotonically, and every centre-block-diagonal SPD
P gives cond(P^{-1}S)/cond(S) ≥ 1 − 1e-9. The wrong answer it rejects is "a per-centre / per-l
preconditioner bounds the sigma metric". A related walls-register update (conditioning row:
"OPEN on a sigma basis" → STANDING for preconditioning, with canonical orthogonalisation as the
BORDERLINE remedy) is a PI/PM call.
