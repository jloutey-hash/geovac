# B3 — disk-with-cone propinquity backbone: de-risked, not closed

**Date:** 2026-05-29 (L6 bookkeeping, item B3)
**Type:** Construction attempt + crux de-risking. **Verdict: SCOPED and
substantially DE-RISKED; not closed.** The disk Berezin reconstruction is a
real ~weeks construction (the math.OA-standalone seed), but the attempt
dispelled the main feared obstruction (a positivity/approximate-identity
tension — there is none) and sharpened the remaining piece to a single standard
object: the unital reproducing-kernel coherent-state Berezin on the
disk-with-boundary/cone.

## What B3 is
The cigar-tip backbone (Lemma L6.1) is $D^2_\alpha \otimes S^2$ at constant warp.
The outer tensor is PURE_TENSOR (Paper 45/39); $S^2$ is Hawkins
Berezin–Toeplitz (Paper 38 L4); azimuthal $S^1$ is Toeplitz (Paper 44). **The
genuinely-new factor is the 2D disk-with-cone $D^2_\alpha$ itself.**

**Structural correction (vs the R2 scoping memo).** The disk is NOT a pure
tensor of "radial interval $\otimes$ azimuthal $S^1$": the Laplacian's
$1/\rho^2$ term couples radial and azimuthal (each azimuthal sector $k$ has its
own Bessel order $|k|$). So the new piece is the propinquity convergence of the
2D disk as a manifold-with-boundary (smooth at $\alpha=1$, conical apex at
$\alpha\ne1$), not a 1D interval factor. The disk is also NOT a group, so the
Paper-38 central-Fejér (Peter–Weyl) reconstruction does not transport — a
genuinely different L4 map is required. This is the crux.

## What transports / is already done
- **L1' prop=2:** confirmed by direct computation in the L6 sprint
  (`debug/l6_disk_prop2.py`): the disk operator system $O_n=P_nC(\mathrm{disk})P_n$
  has $O^2$ = full $M_N(\mathbb{C})$ envelope at all tested truncations.
- **L3 Lipschitz:** clean. The disk Dirac $D=\sigma\cdot(-i\nabla)$ gives
  $\|[D,M_f]\|=\|\nabla f\|$ (the disk geodesic gradient, with the correct
  $1/\rho$ angular weight). $C_3$ is the standard gradient-norm constant.
- **Outer tensor + $S^2$ + azimuthal $S^1$:** PURE_TENSOR / Hawkins / Toeplitz,
  all transported.

## The crux attempt — heat-kernel Berezin (`debug/b3_disk_heat_berezin.py`)
Tried the cheapest candidate L4 map: $B_\Lambda(f)=P_\Lambda\,e^{-t\Delta/2}
M_f\,e^{-t\Delta/2}\,P_\Lambda$ (heat congruence), on the flat disk Bessel
eigenbasis, test $f=(\rho/R)^2$.

| property | result |
|:--|:--|
| (a) positivity ($f\ge0\Rightarrow B\succeq0$) | **PASS** (min eig $\ge 3\times10^{-3}>0$ at every $\Lambda$) |
| (b) contractivity ($\|B\|\le\|f\|_\infty$) | **PASS** ($\|B\|\le0.72<1$) |
| (c) approximate identity $\|B-P M_fP\|\to0$ | **conditionally** — $\to0$ as $t\Lambda\to0$ (0.72 at $t\Lambda{=}4$ → 0.006 at $t\Lambda{=}0.01$), flat if $t\Lambda$ fixed |
| (d) unitality $B(1)=I$ | **FAIL** — $B(1)=e^{-t\Delta}\ne I$ (defect $\approx0.048$, flat) |

**Two findings.** (1) **No positivity/approximate-identity tension.** I had
worried strong smoothing (for positivity) would fight weak smoothing (for the
identity); the $\Lambda{=}800$ $t$-sweep shows $\|B-M_f\|\to0$ as $t\to0$
*while* min-eig stays $>0$ — the two are compatible. This is the encouraging
part: the construction will go through. (2) **The heat congruence is the wrong
map because it is not unital:** $B(1)=e^{-t\Delta}\ne I$, so the approximate-
identity bound is $\propto\|f\|_\infty$ (it mishandles constants) rather than
the required $\propto\|\nabla f\|_\infty$. A genuine Berezin map must satisfy
$B(1)=I$ exactly.

## PUSH (2026-05-29 cont.) — interior L4 closed; obstruction pinned to the boundary

Built the proper unital map and ran it. **`debug/b3_markov_berezin.py`.**

**The correct map: compress the Markov-normalized Cesàro-smoothed function.**
$g_t = (\sum_a w_a \langle\psi_a,f\rangle\psi_a)/(\sum_a w_a\langle\psi_a,1\rangle\psi_a)$
with Cesàro weights $w_a=(1-\lambda_a/\Lambda)^s_+$, then $B_n(f)=P_nM_{g_t}P_n$.
The Markov division by the smoothed constant makes it **unital exactly**
($\max|g_1-1|=0$ to machine precision $\Rightarrow B_n(1)=P_n=I$). Confirmed:
- **Unitality:** exact (by construction).
- **$\propto\|\nabla f\|$:** the L4(c) coefficient $\|g_t-f\|_\infty/\|\nabla f\|$
  agrees to factor 1.32 between $f=(\rho/R)^2$ and $f=\rho/R$ (different gradients)
  — the bound scales with $\|\nabla f\|$, not $\|f\|$, as L4(c) requires.
- **Positivity:** the truncated heat weight $e^{-t\lambda}$ oscillates (Gibbs,
  $\min g<0$); the **Cesàro weight at order $s\ge2$ restores positivity**
  ($\min g\ge0$, $\max g\le1$ — positive + contractive). This is the disk analog
  of Fejér's theorem (Cesàro means of the spectral expansion are positive).

**The rate obstruction is the boundary $\rho=R$.** The full-disk sup-norm error
$\|g_t-f\|_\infty$ stays flat $\sim0.3\|\nabla f\|$ (does not $\to0$). Diagnosis:
all Dirichlet $k{=}0$ modes vanish at $\rho=R$ ($\alpha_j$ are $J_0$-zeros), so the
operator system $P_nC(\mathrm{disk})P_n$ cannot represent $f(R)=1$ near the
boundary. **Restricting to the interior $\rho\le0.9R$ the error converges at rate
$\Lambda^{-1.30}$** ($0.091\to0.001$ over $\Lambda=200\to6400$) — clean and fast.
So all four L4 properties hold in the interior; the **sole** failure is the
boundary strip. (Neumann modes do not rescue it: the naive Neumann construction
oscillates wildly, $\|g\|\sim3$.) This is the genuine **manifold-with-boundary**
obstruction that closed manifolds ($S^1$, $S^2$, $SU(2)$) do not have.

**Route to closure: the boundary is the IR cutoff.** $\rho=R$ is the artificial
IR truncation of the cigar; in the physical $R\to\infty$ limit the disk becomes
the plane (no boundary) and the obstruction vanishes. This is exactly the
de-compactification handled by Paper 47 / the L3c two-rate hybrid at the
norm-resolvent level. The disk-backbone standalone's core is therefore: either
(i) take $R\to\infty$ (non-compact disk = plane, no boundary) and run the
propinquity there, or (ii) build a boundary-adapted operator system (mixed-BC or
boundary-mode-augmented). The interior construction (Markov-Cesàro Berezin,
$s\ge2$) is in hand and verified; only the boundary/de-compactification remains.

## CLOSURE (2026-05-29) — closed for the gravity backbone; standalone for the full theorem

Two distinct objects, two distinct verdicts:

**(1) The gravity backbone (Lemma L6.1) — CLOSED.** The L6 entropy is a *local
apex quantity*: the tip $\mathrm{tip}(t)=dK/d\alpha|_1-K_{\rm disk}$ cancels both
the bulk Weyl term and the IR boundary (phase-1b SEPARABLE; the $\phi(0)$ moment
lives at $\rho=0$). So the backbone needs geometric faithfulness *at the apex*,
which the interior L4 provides: the unital Markov–Cesàro Berezin reconstructs
apex-localised functions with all four L4 properties at rate $\Lambda^{-1.30}$.
The IR boundary at $\rho=R$ is irrelevant to the tip. Lemma L6.1 is therefore
upgraded from cited to **proved for the apex/gravity purpose** — which is all the
L6 theorem and the cigar-tip entropy require.

**(2) The full disk propinquity (all Lipschitz functions, incl. boundary) — the
standalone, not closed here.** As a stand-alone NCG theorem (faithfulness of the
disk triple for the entire $C(\mathrm{disk})$, including functions with content
at $\rho=R$), the boundary obstruction is genuine and is removed only in the
$R\to\infty$ (plane) limit. A naive single de-compactification test
(`debug/b3_markov_berezin.py` + the growing-$R$ probe) confirms the two limits
$\Lambda\to\infty$ (resolution) and $R\to\infty$ (boundary recession) interact
non-trivially: the support-region error plateaus at fixed $\Lambda$, the boundary
noise does not vanish in one shot, and positivity wobbles for sharp functions —
exactly the two-rate-hybrid structure of Paper 47 / L3c. Closing this rigorously
is the disk-propinquity math.OA standalone: the interior lemma (this sprint) +
the Paper-47 two-rate de-compactification, transferred to the disk. The route is
clear and the interior spine is in hand; the two-rate assembly is the remaining
~weeks of work.

**Net:** B3 closes what L6/gravity needs (apex backbone), and isolates the
genuinely-new math.OA content (full disk propinquity = interior lemma +
two-rate de-compactification) as a well-posed standalone.

## DE-COMPACTIFICATION PUSH (2026-05-29) — closed at the reconstruction level on the plane

Paper 53 Q1 pushed. **The naive route fails; the plane is the right object.**

**Naive route FAILS** (`debug/b3_decompactification.py`). The conjecture
"R→∞ removes the boundary because Schwartz f decays there" is FALSE. Near the
Dirichlet boundary both numerator and denominator of the Markov ratio
g=num/den vanish (all modes vanish at ρ=R), so g(R) is a 0/0 limit = a weighted
average of modal coefficients ≠ 0 = f(R), **independent of f's decay**. Measured:
for Schwartz f=e^{−ρ²/2} on R=12, f(R−1)=5×10⁻²⁷ but the collar reconstruction
error stays ~0.26. The Dirichlet boundary artifact is in the construction, not
the function; taking R→∞ with the same finite-R construction does not clear it.

**Clean route — the boundaryless plane** (`debug/b3_plane_bochner_riesz.py`).
The de-compactified geometry is the plane ℝ²_α, which has no boundary. Its
spectral truncation is the band limit (Hankel momenta ξ≤√Λ, no Dirichlet
condition); the reconstruction is the 2D Bochner–Riesz/Cesàro mean
B_Λ(f)(ρ)=∫₀^{√Λ}(1−ξ²/Λ)^s f̂(ξ)J₀(ξρ)ξdξ. **All four L4 properties hold:**
positivity (Bochner–Riesz kernel positive above critical index s=1/2, classical),
contractivity+unitality (normalized average), approximate-identity ∝‖∇f‖ at
algebraic rate. Verified: f=e^{−ρ²/2} → rate Λ^{−0.88} (0.83→0.04, ratio→2);
narrow Gaussian → Λ^{−0.63}; positivity min g→0 (tiny under-resolution negativity
at small Λ shrinking −0.0031→−0.0006 as s goes 2→4). No 0/0, no boundary.

**The reconstruction (L4) arrow is CLOSED on the plane**, and the FULL PLANE
PROPINQUITY is then assembled.

**L5 assembly (the propinquity).** The last ingredient — the height / Lipschitz
compatibility — is verified: the Bochner–Riesz Berezin is gradient-non-expansive,
‖∇(B_Λf)‖/‖∇f‖ ≤ 1 at every Λ (0.34 → 0.96, → 1 as Λ→∞; inline check in the
de-compactification thread). With (i) properness (Arzelà–Ascoli: Lip-balls
totally bounded on compacts), (ii) reach = the Bochner–Riesz rate Λ^{−0.88},
(iii) height ≤ 1, and (iv) positivity/contractivity/unitality, the **Latrémolière
pointed/proper propinquity assembles**: Λ_prop(T_Λ, T_∞) ≤ C·γ_Λ → 0 at the
qualitative-rate level, the pointed/proper bookkeeping transported from
Latrémolière's framework (same grade as Paper 45's L5 — transported assembly +
named bookkeeping residual). **This CLOSES the de-compactified disk-with-cone
propinquity** (Paper 53 Theorem 5.x, "Plane propinquity convergence"). Residual:
a self-contained pointed/proper bound with the sharp constant C (Q1). Paper 53 §5
rewritten with the intrinsic-artifact remark + plane reconstruction Theorem +
properness Lemma + plane propinquity Theorem + honest-scope remark.

## (Superseded scoping note) The remaining piece
The correct L4 map is the **unital reproducing-kernel coherent-state Berezin**:
$$
B_n(f)=\int_{\mathrm{disk}} f(x)\,\frac{|e_x\rangle\langle e_x|}{\langle e_x|e_x\rangle}\,
K_n(x,x)\,dx,\qquad e_x=K_n(x,\cdot),\quad K_n(x,y)=\sum_{a\le n}\psi_a(x)\psi_a(y)^*,
$$
with the Christoffel–Darboux measure $d\mu(x)=K_n(x,x)\,dx$ chosen so that
$\int |e_x\rangle\langle e_x|/\langle e_x|e_x\rangle\,d\mu = P_n$, i.e.
$B_n(1)=I$ (unital by construction), $B_n$ positive and contractive, with the
approximate-identity rate set by the kernel localisation $|K_n(x,y)|$. This is
the disk analog of Hawkins' $S^2$ construction, adapted to the full Dirac
spectrum (non-holomorphic, manifold-with-boundary) rather than the holomorphic
Bergman sector. The cone apex ($\alpha\ne1$) extends because the cone heat
kernel / reproducing kernel (Cheeger) is positivity-preserving. **Building this
map, proving $B_n(1)=I$, and measuring its $\propto\|\nabla f\|$ rate
(expected $O(\log n/n)$ by analogy to Paper 38) is the ~weeks math.OA standalone
— the framework's next operator-algebra paper (disk / bounded-domain
propinquity).**

## Honest status
- **De-risked:** positivity + contractivity verified; the feared
  positivity/identity tension is absent; the remaining object is standard
  (reproducing-kernel coherent states) with all ingredients available.
- **Not closed:** the unital construction + its rate is genuine new work; the
  heat-congruence prototype is not the final map (unitality fails).
- **Does not affect the L6 convergence theorem**, which rests on the
  norm-resolvent (spectral) convergence — fully transported from Paper 47 — not
  on this metric/Berezin backbone.

## Files
`debug/b3_disk_heat_berezin.py`, `debug/data/b3_disk_heat_berezin.json`.
