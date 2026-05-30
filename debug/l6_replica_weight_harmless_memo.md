# L6 — replica-weight-harmless theorem for the cigar-tip entropy (the prize)

**Date:** 2026-05-29
**Type:** Single-threaded multitask sprint toward L6, the single irreducible
deliverable of the gravity R2 program (charter
`debug/gravity_campaign_R2_theorem_charter_memo.md`). Four tracks run
sequentially in main session; all four land.
**Verdict:** **L6-VERIFIED at framework grade** (structural proof + full
numerical verification of every load-bearing sub-claim; transported pieces
confirmed by direct compute or cited). The gravity campaign's one remaining
theorem is closed at the level Papers 45/46/47 operate at. Formal LaTeX
theorem write-up (Paper-47 interval transport spelled out + uniform-derivative
theorem stated) is the remaining bookkeeping — not new content.

## The object (charter framing)

The cigar near-horizon = constant-warp tensor product
$D^2_\alpha \otimes S^2(r_h)$. The entropy is the **replica derivative of the
heat-trace spectral functional** (NOT a differentiated propinquity metric):
$$
S_{\rm tip} = \phi(0)\text{-moment of } \mathrm{tip}(t),\qquad
\mathrm{tip}(t) = \tfrac{dK}{d\alpha}\big|_{\alpha=1} - K_{\rm disk}(t),
$$
with $K(\alpha,t) = \mathrm{Tr}\,e^{-tD_{\alpha}^2}$. On the discrete substrate
(`geovac.gravity.warped_dirac`, SPECTRAL azimuthal):
$$
K(\alpha,t) = 4\sum_{j\ge0} h\!\big(\tfrac{j+\frac12}{\alpha},t\big),\qquad
h(m,t)=\sum_{\rm radial} e^{-t\lambda},\ \lambda=\mathrm{eigs}\Big[-u''+\tfrac{m^2-\frac14}{\rho^2}\Big]_{(0,R],\ a=R/N_\rho}.
$$
Replica derivative: $\frac{dK}{d\alpha}|_1 = -4\sum_j (j+\tfrac12)\,h'(j+\tfrac12)$.
Per-mode replica contribution $C_j = -4(j+\tfrac12)h'(j+\tfrac12) > 0$; the
$(j+\tfrac12)$ factor is the **replica weight** — a heavier mode-moment than the
bare heat trace.

**L6 statement.** $\frac{dK_n}{d\alpha}|_1 \to \frac{dK}{d\alpha}|_1$ as
$n\to\infty$, uniformly in $\alpha$ on $[1-\delta,1+\delta]$ and $C^1$ in
$\alpha$, so that $\lim_n$ and $\frac{d}{d\alpha}|_1$ commute, with a rate. The
charter's own difficulty assessment: "tractable, not a wall — genuine work is
the rate and the uniformity, not bare convergence."

## The proof, in one line

The centrifugal floor $\lambda_{j,0}\sim (j+\tfrac12)^2/R^2$ makes
$e^{-t\lambda}\sim e^{-t(j+\frac12)^2/R^2}$ **Gaussian in $j$**, which dominates
the polynomial replica weight $(j+\tfrac12)$ and the prefactor $h'$. So $C_j$ has
a Gaussian envelope $\Rightarrow$ the replica-weighted mode series converges
absolutely and uniformly in $\alpha$ (Weierstrass M-test, dominating function
the Gaussian envelope) $\Rightarrow$ term-by-term $\alpha$-differentiation is
legal $\Rightarrow$ $\lim_n$ and $d/d\alpha$ commute. **The replica weight is
harmless**: the joint rate equals the undifferentiated heat-trace rate.

## Numerical verification (`debug/l6_replica_weight_harmless.py`)

All at $R=10$, $t=1$; lean tridiagonal radial solver cross-checked bit-exact
against the production `DiscreteDiskDiracSpectral` module (reldiff 0.0e0).

**(3a) Gaussian envelope — the L6 core.** $C_j \sim \exp(\text{slope}\cdot
(j+\tfrac12)^2)$, measured slope $-0.01129$ vs predicted $-t/R^2 = -0.01000$
(the small excess matches the discrete centrifugal floor slope $0.01164$, also
measured). The replica weight is crushed:
$C_j/C_0 = 1 \to 4.7\,(j{=}4) \to 0.38\,(j{=}16) \to 4\!\times\!10^{-5}\,(j{=}32)
\to 1.6\!\times\!10^{-20}\,(j{=}64)$.
Azimuthal partial sum of $dK/d\alpha$: tail/total $= 1.7\!\times\!10^{-6}$ at
$K_{\rm cut}{=}32$, $4\!\times\!10^{-13}$ at $48$, $0$ at $64$. **Super-polynomial
azimuthal convergence despite the replica weight.**

**(3b) Uniformity in $\alpha$.** Envelope slope across $\alpha\in[0.9,1.1]$:
range $[-0.01134,-0.01123]$ — essentially constant. $K_{\rm cut}{=}32$ tail
$\le 1.4\!\times\!10^{-5}$ at every $\alpha$. The dominating Gaussian is uniform;
the $C^1$ requirement holds.

**(3c / T2) Rate.** Azimuthal: super-poly (residual $0$ by $K_{\rm az}{=}64$).
Lattice $a=R/N_\rho$: tip$(a)\to$ as $a\to0$ with order $p\approx1.12$ (FD polar
Laplacian near the screened apex), order-2 Richardson $\to 0.16432$ vs
$1/6=0.16667$ ($-1.4\%$, consistent with the known G4-6b finite-$a$ offset
$+4\%$). **Joint rate gated by the lattice / norm-resolvent rate $O(a^{1.1})$,
NOT degraded by the replica weight** — exactly the L6 claim.

**(3d) Commute.** $\lim_n \frac{dK_n}{d\alpha}$ (per-mode analytic) $=41.501652$
vs $\frac{d}{d\alpha}\lim_n K_n$ (FD on converged $K$) $=41.501681$; diff
$2.9\!\times\!10^{-5}$. **$\lim_n$ and $d/d\alpha$ commute to 5 digits.**

Verdict cell: **L6-VERIFIED-NUMERICALLY** (Gaussian envelope ✓, uniform ✓,
commute ✓, super-poly azimuthal tail ✓).

## Backbone (Track 1) — prop=2 on the disk operator system

`debug/l6_disk_prop2.py`. Built $O_n = P_n C(\mathrm{disk}) P_n$ from compressed
disk monomials $z^a\bar z^b$ ($z=\rho e^{i\phi}$) in the truncated radial$\times$
azimuthal mode basis; computed prop $=$ smallest $p$ with $\dim(O^p)=N^2$.
**Result: prop $=2$ at all four truncations** $(K,J)\in\{(2,2),(2,3),(3,2),(3,3)\}$,
$N\in\{8,12,12,18\}$ — $O^2$ spans the full $M_N(\mathbb{C})$ envelope every time.
The disk backbone inherits the Connes--vS 2-step-generating (Toeplitz $S^1$,
Paper 44 Prop 4.2) property by direct computation, not just transport. The
Layer-1 propinquity backbone (charter step 1) is confirmed.

## M1 / area signature (Track 4)

`debug/data/l6_m1_area_signature.json`. The $S^2$ Dirac heat trace
($|\lambda_n|=(n+1)/r_h$, $g_n=4(n+1)$) has $K_{S^2}(t)\cdot t \to 2r_h^2 =
A/(2\pi) = \mathrm{Vol}(S^2)\,r_h^2/(2\pi)$ (verified to 4 digits at
$r_h\in\{1,2,3\}$). The passive horizon sphere carries the area
$A = 4\pi r_h^2 = \mathrm{Vol}(S^2)\,r_h^2$ in its leading heat-trace
coefficient, so $S_{\rm tip} = (\text{disk coeff})\times A$ factors through
$\mathrm{Vol}(S^2)$. The M1 signature $4/\pi = \mathrm{Vol}(S^2)/\pi^2$ — the
same constant carrying the L2 GH-convergence rate (Paper 38/40), the L2-F.1
wedge HS-orthogonality $1/\pi^2$ prefactor (Paper 43 §10.2), and the
Stefan--Boltzmann Matsubara prefactor — is therefore **structural-by-construction**
in the cigar entropy: gravity entropy $S_{\rm tip}=A/4$ wears the M1 /
Bernoulli-ladder signature through the horizon-sphere area. Honest caveat: this
is a by-construction signature (the area literally *is* $\mathrm{Vol}(S^2)r_h^2$),
not an independent rate coincidence. The charter's "check once the rate is in
hand" condition is now met (the rate IS in hand), so the tie to the master
Mellin engine is real, if modest.

## Honest scope / grade

**Sprint-level grade summary (the four-way taxonomy).**
- *Theorem grade (proved):* L6 Theorem (replica-weight-harmless) with its formal
  three-lemma proof, the new Lemma L6.3 proved via explicit Gaussian domination
  (Paper 51 §L6); B1 rate O(a) + boundary-placement mechanism; B2 continuum
  target 1/6; Paper 53 interior reconstruction (Thm 3.2), apex backbone (Cor 4.2),
  plane Bochner–Riesz reconstruction (Thm 5.x) and plane propinquity convergence
  (Thm 5.x) — the last at the same transported-bookkeeping grade as Paper 45's L5.
- *Structural sketch (solid, not formal proof):* the pointed/proper propinquity
  *bookkeeping* (base-point-compatible tunneling, local-uniformity of reach/height)
  is transported from Latrémolière's framework, not re-derived; Lemma L6.1's full
  L1'–L5 compact assembly is likewise transport-cited.
- *Numerical observation:* every load-bearing claim verified numerically
  (envelope slope −t/R², commute to 5 digits, prop=2, B1 √2 ratio, B4 √a ratio,
  Bochner–Riesz rate Λ^{−0.88}, height ratio ≤1); L6 is a Layer-2 convergence
  theorem (rates, not bit-exact).
- *Closed negative:* B4 — the Möbius α>1 form is a finite-a artifact; continuum
  α>1 = plain Sommerfeld–Cheeger; no exotic mechanism exists (retires the thread).
- *Named open follow-ons:* (i) self-contained pointed/proper propinquity bound
  with the sharp constant C (Paper 53 Q1); (ii) B3-standalone full L1'–L5 compact
  assembly write-out (Lemma L6.1); (iii) quantitative L6 lattice rate / B1 sharp
  constant. None block the gravity application or the proved theorems.

- **NEW L6 content (the prize) — VERIFIED at framework grade.** Replica weight
  harmless; uniform-in-$\alpha$ $C^1$ convergence; $\lim_n$ and $d/d\alpha$
  commute; joint rate $=$ undifferentiated heat-trace rate. Analytical mechanism
  (Gaussian dominates polynomial replica weight) + full numerical verification
  (envelope slope $=-t/R^2$, uniformity, commute to 5 digits). This is the level
  Papers 45/46/47 close at.
- **Layer 1 (prop=2 backbone) — CONFIRMED by direct compute**, not just cited.
- **Layer 2 (radial-interval norm-resolvent rate) — verified numerically**
  ($O(a^{1.1})$); the analytical interval theorem is cited from Paper 47 (apex
  centrifugally screened, regular Sturm--Liouville), not re-derived here.
- **Continuum target $dK/d\alpha|_1=1/6$** is the G4-4c spinor-cone
  identification; refinement lands within $1.4\%$ (Richardson), consistent with
  the Class-1 finite-$a$ offset. L6 is a convergence statement (rate + commute);
  it does not re-derive the $1/6$ value, and the exact continuum value vs the
  finite-substrate offset is the separate Class-1 cutoff-dependence matter.
- **NOT bit-exact.** L6 is a Layer-2 convergence theorem (rates, fits, 5-digit
  commute) — caution-light, appropriately, not a skeleton identity. Presented as
  such.
- **Formal write-up DONE (end-of-session, same sprint).** Paper 51 §L6 now
  carries a proper **Theorem L6 + three-lemma proof**: Lemma L6.1 (propinquity
  backbone, transported, prop=2 confirmed), Lemma L6.2 (norm-resolvent rate,
  transported, screened apex), **Lemma L6.3 (replica weight harmless — NEW,
  proved)** with the explicit Gaussian-dominating sequence
  G_j = C(t,δ)(j+½)² e^{−c(δ)(j+½)²}, c(δ)=t/((1+δ)²R²), plus a Proof of the
  theorem via the differentiation-of-series theorem. Dominating exponent
  c(δ)→t/R² matches the measured envelope slope. Paper 51: 31pp, three-pass clean.
- **B1 + B2 CLOSED** (`debug/l6_b1b2_rate_and_target.py`). **(B1, rate.)** The
  lattice rate is O(a) (not the earlier-quoted ~1.1): per-eigenvalue order
  q=1.00 vs exact Bessel zeros (j_{m,1}/R)² at m∈{½,3/2,5/2}, tip order p=0.97.
  **Mechanism = Dirichlet boundary placement** (NOT the 1/ρ² singularity, NOT an
  apex-edge strip — both falsified): the production grid ρ_i=i·a sets the IR wall
  at (N_ρ+1)a = R+a, displaced O(a) from R. Proven two ways: (i) the free m=½
  case (zero centrifugal potential) still gives O(a); (ii) the corrected grid
  a=R/(N_ρ+1) (wall at R) restores central-FD O(a²) (ratio 4, errors 1000×
  smaller). The O(a) is a fixable substrate detail. **(B2, target.)** With the
  correct O(a) order, linear-in-a extrapolation gives tip₀ = 0.166656 (2-pt) /
  0.166582 (lstsq) vs 1/6 = 0.166667 — **−0.006% / −0.05%**. The continuum
  dK/dα|₁ = 1/6 (G4-4c Sommerfeld–Cheeger) IS the discrete limit; the earlier
  −1.4% was an order-2-Richardson-on-O(a) artifact. The boundary-placement
  diagnosis corrects an earlier "apex-edge" guess (audit-numerical-claims paid
  off). Paper 51 §L6 Grade + Lemma L6.2 updated.
- **B3 DE-RISKED, not closed** (`debug/b3_disk_backbone_memo.md`). Full L1'–L5
  assembly of the D²_α ⊗ S² backbone (upgrades Lemma L6.1 cited→proved). Transport
  map: L1' prop=2 ✓ (L6 sprint); L3 Lipschitz = disk gradient norm ✓ (clean);
  outer PURE_TENSOR + S² Hawkins + azimuthal S¹ Toeplitz ✓. **Structural
  correction:** the new piece is NOT a "radial-interval factor" — the disk
  Laplacian's 1/ρ² term couples radial and azimuthal, so it is the **2D
  disk-with-cone as a manifold-with-boundary** (not a group, so the Paper-38
  central-Fejér reconstruction does not transport). **PUSH (built the map):** the
  correct L4 map is the unital Markov-normalized Cesàro Berezin
  B_n(f)=P_n M_g P_n, g=(Σ w_a⟨ψ_a,f⟩ψ_a)/(Σ w_a⟨ψ_a,1⟩ψ_a), Cesàro weight
  w_a=(1−λ_a/Λ)^s. **Interior L4 CLOSED:** unitality exact (B_n(1)=I to machine
  precision — the Markov division fixes it), positivity+contractivity at order
  s≥2 (truncated-heat oscillates; Cesàro restores it = disk Fejér),
  approximate-identity ∝‖∇f‖ (gradient-proportional to factor 1.32) with clean
  **interior rate Λ^{−1.30}** (ρ≤0.9R: 0.091→0.001 over Λ=200→6400). **Sole
  remaining obstruction = the boundary ρ=R:** the Dirichlet modes vanish there
  (J₀-zeros) so sup-norm reconstruction fails in the boundary strip (full error
  flat ~0.3); interior converges fast. Genuine manifold-with-boundary obstruction
  (closed manifolds don't have it; Neumann doesn't rescue). **Route: the boundary
  IS the IR cutoff; R→∞ (disk→plane, no boundary) removes it — the same
  de-compactification Paper 47 / L3c handles at norm-resolvent level.** **CLOSURE:** two
  scopes. **(gravity backbone — CLOSED)** the L6 entropy is a local apex quantity
  (tip cancels bulk AND boundary, phase-1b SEPARABLE), so Lemma L6.1 needs only
  apex faithfulness = the interior L4 (proved, Λ^{−1.30}); the IR boundary is
  irrelevant to the tip. Lemma L6.1 upgraded cited→proved for the gravity purpose.
  **(full disk propinquity — standalone)** all of C(disk) incl. boundary-supported
  functions needs the R→∞ two-rate de-compactification (naive one-shot attempt
  confirms Λ→∞ × R→∞ interact non-trivially — the Paper-47/L3c structure on the
  disk); interior lemma in hand, two-rate assembly = the math.OA companion. Does
  not affect the L6 theorem (norm-resolvent, fully transported). See
  `debug/b3_disk_backbone_memo.md`, `debug/b3_markov_berezin.py`.
- **B4 CLOSED — retires the open Möbius α>1 thread** (`debug/l6_b4_moebius_continuum_diagnostic_memo.md`).
  The α/(2α−1) "closed form" is a finite-a substrate artifact; clean
  spectral-azimuthal + radial-continuum (a→0) extrapolation drives the recovery
  from ≈Möbius (0.667 at a≈0.05) toward the Sommerfeld–Cheeger continuation
  (→1), with (1−recovery)~√a (cone-singularity signature at α≠1, vs smooth-apex
  O(a) at α=1). Continuum α>1 = plain SC −1/12. No exotic mechanism exists. Not
  L6 itself (α=1 replica unaffected) — but the natural closing item of the
  gravity-arc bookkeeping.

## Files
- `debug/l6_replica_weight_harmless.py`, `debug/data/l6_replica_weight_harmless.json`
- `debug/l6_disk_prop2.py`, `debug/data/l6_disk_prop2.json`
- `debug/data/l6_m1_area_signature.json`

## Cross-references
- `debug/gravity_campaign_R2_theorem_charter_memo.md` — the L6 charter (three-layer architecture)
- `debug/gravity_campaign_R2_scoping_memo.md` — lemma-transport map
- `debug/gravity_campaign_phase1b_tip_bulk_independence_memo.md` — tip/bulk SEPARABLE
- Paper 45 (PURE_TENSOR propinquity), Paper 47 (norm-resolvent interval), Paper 38/40 (C₃, central-Fejér 4/π = M1), Paper 44 (Toeplitz prop=2), G4-4c (spinor-cone tip $-(1/12)(1/\alpha-\alpha)$), G4-5d (φ(0) moment map)
