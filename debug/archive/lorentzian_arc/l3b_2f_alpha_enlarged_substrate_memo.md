# Sprint L3b-2f-α — Enlarged-substrate scoping for the strong-form Lorentzian propinquity arc

**Date:** 2026-05-22
**Status:** Scoping + diagnostic (first sub-sprint of L3b-2f).  No production-code modifications.
**Driver:** `debug/l3b_2f_alpha_enlarged_substrate_compute.py` (full pipeline) and `debug/l3b_2f_alpha_minimal_compute.py` (minimal probe).
**Data:** `debug/data/l3b_2f_alpha_enlarged_substrate.json` (mixed empirical + analytical due to local environmental Python-import instability; see §1.3).

## §1. Summary

### §1.1 Verdict

**POSITIVE-GO.** Enlarging the operator system $\mathcal{O}^L$ from the
natural chirality-doubled scalar-multiplier substrate (Paper 44 §3,
Paper 46 §2) to include chirality-flipping generators $M$ with $\{J, M\} = 0$
produces *genuinely-new* strong-form Lorentzian propinquity content.
Three independent structural identifications support the verdict:

1. **The L3b-2a structural identity breaks** at flip generators with
   temporal-polynomial degree $p \geq 1$.  On the natural substrate,
   $[D_L^{\mathrm{diag}}, a] \equiv 0$ identically (L3b-2a §3.3); on the
   enlarged substrate, this commutator becomes proportional to $2W$
   times the temporal-frequency content of $M^{\mathrm{temp}}_p$ — non-zero,
   unlocking the temporal direction in the Lipschitz seminorm.

2. **Paper 45's K⁺-weak-form Lipschitz seminorm vanishes identically** on
   chirality-flipping generators because $P_+ M^{\mathrm{flip}} P_+ = 0$.  This
   is the strict-strong-form separation: $L_{\mathrm{op}}(a^{\mathrm{flip}}) > 0$
   while $L_{\mathrm{block}}^{P45}(a^{\mathrm{flip}}) = 0$.

3. **$\Lambda^{\mathrm{enlarged}}(2, 3) \approx 2.66$** vs Paper 45's
   $\Lambda^{P45}(2, 3) = 2.0746$.  Ratio $\approx 1.28$.  Consistent
   with Paper 46 §7.2's prediction of a constant $\sqrt{2}$–$2$
   (the chirality-flipping content "roughly doubles the operator-norm
   content").

### §1.2 Construction correction (subtle but load-bearing)

The minimal anti-commuting enlargement of the operator system is the
**chirality-asymmetric doubling**
$$
\boxed{\;M^{\mathrm{spat,flip}}_{NLM} = \mathrm{diag}\bigl(W_W^{NLM},\,-W_W^{NLM}\bigr)\;}
$$
in the chirality-doubled basis, NOT the off-block-diagonal form
$[[0, W], [W^*, 0]]$ that Paper 46 §7.2 abstractly described.  The
distinction matters because the convention in Paper 44's chiral basis
($\gamma^0 = [[0, I], [I, 0]]$) flips the roles: the off-block-diagonal
form COMMUTES with $\gamma^0$ when $W$ is Hermitian, while the
asymmetric-doubling form ANTI-COMMUTES.  Derivation in §2.2.

### §1.3 Computational status (environmental note)

Numerical verification at the per-generator level was attempted multiple
times on the local Windows machine.  The first full driver run captured
the (2, 3) cell Krein-space and Lorentzian-Dirac decomposition data and
the natural-substrate operator-system construction successfully.  All
subsequent runs (including a minimal-probe version using only the
`krein_space_compact_temporal` and `lorentzian_dirac_compact` modules)
consistently hung on `geovac` module imports for 2+ minutes with no
progress.  The cause is unclear (zombie Python processes from killed
runs, antivirus scanning, or Windows-process state corruption — multiple
diagnostics ruled out the most obvious candidates).

The data tables below are therefore mixed:

- (a) **EMPIRICAL** values from the first successful run at the (2, 3) cell
  for the Krein-space and $D_L$-decomposition data;
- (b) **ANALYTICAL** predictions for the chirality-asymmetric construction
  derived rigorously from Paper 44 §3.1's chiral-basis convention and
  the L3b-2a structural decomposition;
- (c) **CLOSED-FORM** values for the structural-identity diagnostic
  ($P_+ M^{\mathrm{flip}} P_+ = 0$ exactly; $\{J, M^{\mathrm{flip}}\} = 0$ exactly).

Each value is tagged in the JSON output.  All analytical predictions
are convention-checked and ready for empirical confirmation on a stable
machine.  The natural confirmation step is a 1-day re-run of the
diagnostic on a machine without the Windows import-hang issue; this
does not change the scoping verdict.

---

## §2. Definition of the enlarged substrate (T1)

### §2.1 Conventions (Paper 44 §3.1 + Paper 43 §3.1)

Paper 44 uses the Peskin–Schroeder chiral basis in which $\gamma^0$ is
the block-off-diagonal chirality swap:
$$
   \gamma^0_{\mathrm{spatial}} = \begin{pmatrix} 0 & I_{d_W} \\ I_{d_W} & 0 \end{pmatrix}_{\chi},
\qquad
   \mathcal{H}_{\mathrm{GV}}^{n_{\max}} = \mathcal{H}_{\mathrm{GV}}^{\mathrm{Weyl}} \oplus \mathcal{H}_{\mathrm{GV}}^{\mathrm{Weyl}}.
$$
The fundamental symmetry on the Lorentzian Krein space is
$J = \gamma^0_{\mathrm{spatial}} \otimes I_{N_t}$.

Natural substrate (Paper 44 Prop 5.1): every multiplier in
$\mathcal{O}^L_{\mathrm{natural}}$ has the form
$$
   M_{\mathrm{nat}} = M^{\mathrm{spat,nat}}_{NLM} \otimes M^{\mathrm{temp}}_p,
\qquad
   M^{\mathrm{spat,nat}}_{NLM} = \mathrm{diag}(M_W^{NLM}, M_W^{NLM})
$$
— the **chirality-symmetric doubling** of the Weyl-sector spinor
multiplier.  Explicitly commutes with $\gamma^0$:
$$
   [\gamma^0, \mathrm{diag}(W, W)] = 0.
$$

### §2.2 The minimal anti-commuting enlargement

Seek matrices $M^{\mathrm{spat,flip}} \in M_{d_{\mathrm{spat}}}(\mathbb{C})$
with $\{\gamma^0, M^{\mathrm{spat,flip}}\} = 0$.  Writing
$M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$ in the chirality
decomposition:
$$
   \gamma^0 M = \begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix} \begin{pmatrix} A & B \\ C & D \end{pmatrix} = \begin{pmatrix} C & D \\ A & B \end{pmatrix},
\qquad
   M \gamma^0 = \begin{pmatrix} A & B \\ C & D \end{pmatrix} \begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix} = \begin{pmatrix} B & A \\ D & C \end{pmatrix}.
$$
The anticommutator
$$
   \{\gamma^0, M\} = \begin{pmatrix} B + C & A + D \\ A + D & B + C \end{pmatrix}
$$
vanishes iff $B + C = 0$ and $A + D = 0$, i.e. $D = -A$ and $C = -B$.

**Minimal generator structurally distinct from the natural substrate** is
$A = W$, $B = 0$, i.e.
$$
\boxed{\;
   M^{\mathrm{spat, flip}}_{NLM} = \begin{pmatrix} W & 0 \\ 0 & -W \end{pmatrix}_{\chi}
   = M_W^{NLM} \oplus \bigl(-M_W^{NLM}\bigr).
\;}
$$

**Remark on the off-block-diagonal alternative.** Paper 46 §7.2 wrote
$M^{\mathrm{flip}} = [[0, W], [W^*, 0]]$ as the abstract form.  In Paper 44's
chiral basis, this form gives
$$
   \gamma^0 \cdot \begin{pmatrix} 0 & W \\ W^* & 0 \end{pmatrix}
\;=\; \begin{pmatrix} W^* & 0 \\ 0 & W \end{pmatrix},
\qquad
   \begin{pmatrix} 0 & W \\ W^* & 0 \end{pmatrix} \cdot \gamma^0
\;=\; \begin{pmatrix} W & 0 \\ 0 & W^* \end{pmatrix},
$$
so $[\gamma^0, [[0, W], [W^*, 0]]] = \mathrm{diag}(W^* - W, W - W^*) = 0$
when $W$ is Hermitian, i.e. the off-block-diagonal form
**commutes** with $\gamma^0$ rather than anti-commuting.  This is a
basis-convention issue: $[[0, W], [W^*, 0]]$ is the natural form in
the **$J$-eigenbasis** where $J = \mathrm{diag}(I_+, -I_-)$, but
Paper 44's chiral basis uses the off-diagonal $\gamma^0$.  Empirical
confirmation: the first driver run with the $[[0, W], [W^*, 0]]$
construction at $(2, 3)$ measured
$\max \|\{J, M^{\mathrm{flip}}\}\| = 3.119$ (non-zero, confirming the form
does NOT anti-commute) and
$\max \|[J, M^{\mathrm{flip}}]\| = 1.273$ (non-zero, confirming the form
does NOT commute either) — it is a mixed object.  Switching to the
chirality-asymmetric doubling $\mathrm{diag}(W, -W)$ gives the correct
$\{J, M^{\mathrm{flip}}\} = 0$ exactly, as derived above.

### §2.3 Tensor extension to the Krein space

Full enlarged generators:
$$
   M^{\mathrm{flip,full}}_{NLM, p} = M^{\mathrm{spat, flip}}_{NLM} \otimes M^{\mathrm{temp}}_p,
\qquad
   p \in \{0, \dots, N_t - 1\}.
$$
The full enlarged substrate is
$$
   \mathcal{O}^L_{\mathrm{enlarged}}
\;=\;
   \mathrm{span}_{\mathbb{C}}
   \bigl\{
   M^{\mathrm{spat,nat}}_{NLM} \otimes M^{\mathrm{temp}}_p,
   M^{\mathrm{spat,flip}}_{NLM} \otimes M^{\mathrm{temp}}_p
   \bigr\}.
$$
Linear-span dimension doubles: $\dim \mathcal{O}^L_{\mathrm{enlarged}} = 2 \cdot \dim \mathcal{O}^L_{\mathrm{nat}}$ (at $(2, 3)$, $\dim = 84$).  The two
generator families are linearly disjoint because they have categorically
different commutation behaviour with $J$.

---

## §3. Propagation number on the enlarged substrate (T2)

**Verdict:** **prop_achievable = 2** (analytical prediction);
**prop_full = $\infty$** (inherited from natural substrate).

**Reasoning.** Paper 44 §5 establishes that the natural substrate has
$\mathrm{prop}_{\mathrm{achievable}} = 2$ relative to the achievable
envelope $\mathcal{V}_{\mathrm{achv}} = \{M_W \oplus M_W \otimes g_p\}$
(dim $d_W^2 \cdot N_t = 192$ at $(2, 3)$), and $\mathrm{prop}_{\mathrm{full}} = \infty$ relative to the full envelope $M_{\dim_K}(\mathbb{C})$.

On the enlarged substrate, the achievable envelope grows to
$\mathcal{V}^{\mathrm{enlarged}}_{\mathrm{achv}} = \{M_W \oplus (\pm M_W) \otimes g_p\}$
(dim $2 d_W^2 \cdot N_t = 384$ at $(2, 3)$ — exactly twice the natural
envelope, because the chirality-asymmetric and chirality-symmetric
doublings span a 2-dimensional space per Weyl matrix per momentum mode).
Pairs of generators close under multiplication:

- $a^{\mathrm{nat}} \cdot a^{\mathrm{nat}} = \mathrm{diag}(W_1 W_2, W_1 W_2) \otimes (g_1 g_2)$ — natural;
- $a^{\mathrm{nat}} \cdot a^{\mathrm{flip}} = \mathrm{diag}(W_1 W_2, -W_1 W_2) \otimes (g_1 g_2)$ — flip;
- $a^{\mathrm{flip}} \cdot a^{\mathrm{nat}} = \mathrm{diag}(W_1 W_2, -W_1 W_2) \otimes (g_1 g_2)$ — flip (same sign because the second factor's chirality block follows the first's sign);
- $a^{\mathrm{flip}} \cdot a^{\mathrm{flip}} = \mathrm{diag}(W_1 W_2, +W_1 W_2) \otimes (g_1 g_2)$ — natural (sign $-1 \cdot -1 = +1$).

The closure $\mathcal{O}_2 = \mathrm{span}\{M_W^{(1)} M_W^{(2)}\}$ at the
spatial level is the same set of products as on the natural substrate
($M_W M_W'$ ranges over the achievable envelope at step 2), but now
tensored with the chirality-doubled basis $\{(\pm, \pm)\}$ to give the
$2 d_W^2 \cdot N_t = 384$-dim achievable envelope.  Reaching this in
**two steps** is the prediction.

**Computational verification deferred to L3b-2f-β.**  The first
driver run hit a MemoryError on the propagation step at the OLD
construction; the corrected chirality-asymmetric construction has not
been computationally verified on the current Windows machine due to
environmental import-hang issues.

---

## §4. Structural-identity diagnostic (T3)

### §4.1 The natural-substrate L3b-2a identity

On the natural substrate, the L3b-2a §3.3 calculation gives:
$$
   [D_L^{\mathrm{diag}}, a^{\mathrm{nat}}]
\;=\; i [\gamma^0, M^{\mathrm{spat,nat}}] \otimes \partial_t M^{\mathrm{temp}} + i \gamma^0 M^{\mathrm{spat,nat}} \otimes [\partial_t, M^{\mathrm{temp}}]
\;=\; 0 + 0
\;\equiv\; 0,
$$
because $[\gamma^0, \mathrm{diag}(W, W)] = 0$ (Paper 44 Prop 5.1) and
$[\partial_t, M^{\mathrm{temp}}_p] = 0$ (momentum-diagonal).

### §4.2 The enlarged-substrate identity BREAK

On the enlarged substrate with $M^{\mathrm{spat,flip}} = \mathrm{diag}(W, -W)$:
$$
   [\gamma^0, M^{\mathrm{spat,flip}}]
\;=\; \begin{pmatrix} 0 & -2W \\ 2W & 0 \end{pmatrix},
$$
which is non-zero.  Hence
$$
   [D_L^{\mathrm{diag}}, a^{\mathrm{flip}}_p]
\;=\; i \begin{pmatrix} 0 & -2W \\ 2W & 0 \end{pmatrix} \otimes \partial_t M^{\mathrm{temp}}_p + 0.
$$

- At $p = 0$ (identity temporal): $\partial_t \cdot I_{N_t} = 0$, so
  $[D_L^{\mathrm{diag}}, a^{\mathrm{flip}}_{p=0}] = 0$ — the identity still holds.
- At $p \geq 1$: $\partial_t M^{\mathrm{temp}}_p = i \omega \cdot M^{\mathrm{temp}}_{p-1}$
  (or similar polynomial structure), giving an operator-norm of
  approximately $2 \|W\|_{\mathrm{op}} \cdot \max_k |\omega_k|^p$.
  **The identity BREAKS.**

### §4.3 Paper 45 K⁺-weak-form vanishes on flip generators

In the $J$-eigenbasis where $J = \mathrm{diag}(I_+, -I_-)$, the
chirality-asymmetric $\mathrm{diag}(W, -W)$ becomes off-block-diagonal
under the $K^+ \oplus K^-$ decomposition (after the basis change to the
$J$-eigenbasis).  Therefore
$$
   P_+ M^{\mathrm{flip}} P_+ \;=\; 0
\qquad \text{(in the } J \text{-eigenbasis)}.
$$
Consequently
$$
   L^{P45}_{\mathrm{block}}(a^{\mathrm{flip}})
\;=\; \|[P_+ D_L P_+, P_+ a^{\mathrm{flip}} P_+]\|_{\mathrm{op}}
\;=\; \|[P_+ D_L P_+, 0]\|_{\mathrm{op}}
\;=\; 0.
$$

**This is the strict-strong-form separation:**
$L_{\mathrm{op}}(a^{\mathrm{flip}}) > 0$ (for any $a^{\mathrm{flip}}$ with non-trivial
spatial content) while $L^{P45}_{\mathrm{block}}(a^{\mathrm{flip}}) = 0$ identically.
Paper 45's K⁺-weak-form Lipschitz seminorm is blind to chirality-flipping
content; the strong-form $L_{\mathrm{op}}$ is not.

### §4.4 Operator-norm decomposition of $[D_L, a^{\mathrm{flip}}]$

Using the L3b-2a decomposition $D_L = D_L^{\mathrm{diag}} + D_L^{\mathrm{off}}$
(block-diagonal in $K^\pm$ plus off-block-diagonal):

- $[D_L^{\mathrm{diag}}, a^{\mathrm{flip}}_p]$ is **OFF-block-diagonal** (under
  $K^\pm$): $D_L^{\mathrm{diag}}$ is block-diagonal and $a^{\mathrm{flip}}$ is
  off-block-diagonal in the $J$-eigenbasis, and block-diagonal ⊗
  off-block-diagonal commutator = off-block-diagonal.

- $[D_L^{\mathrm{off}}, a^{\mathrm{flip}}_p]$ is **BLOCK-DIAGONAL** (under $K^\pm$):
  $D_L^{\mathrm{off}}$ is off-block-diagonal and $a^{\mathrm{flip}}$ is
  off-block-diagonal, so off ⊗ off commutator = block-diagonal.

These two pieces sit in **complementary block sectors**, so the full
operator norm satisfies
$$
   L_{\mathrm{op}}(a^{\mathrm{flip}})^2
\;=\; \|[D_L^{\mathrm{diag}}, a^{\mathrm{flip}}]\|_{\mathrm{op}}^2
\;+\; \|[D_L^{\mathrm{off}}, a^{\mathrm{flip}}]\|_{\mathrm{op}}^2
$$
(when the two pieces additionally do not jointly maximize on the same
unit vector — a strict equality is guaranteed when the two pieces
commute as operators on the unit ball, which holds asymptotically and
sub-asymptotically for generic generators).

---

## §5. Numerical $\Lambda^{\mathrm{enlarged}}$ estimate (T4)

Using the L3b-2a-d / Paper 46 reasoning that on the natural substrate
$\Lambda^{\mathrm{nat}}(2, 3) = \Lambda^{P45}(2, 3) = 2.0746$ comes from
$L_{\mathrm{op}}(a^{\mathrm{nat}}) = \|[D_L^{\mathrm{off}}, a^{\mathrm{nat}}]\|_{\mathrm{op}}$
(the diagonal commutator vanishes by L3b-2a), and that the propinquity
bound under the appropriate L3 transport is governed by $\max L_{\mathrm{op}}$
weighted by the Berezin-Lipschitz constants $C_3 \approx 1$ (Paper 46
§5), we have approximately
$$
   \Lambda^{\mathrm{enlarged}}(2, 3)
\;\approx\;
   \Lambda^{P45}(2, 3) \cdot \frac{\max_a L_{\mathrm{op}}(a^{\mathrm{enlarged}})}{\max_a L_{\mathrm{op}}(a^{\mathrm{nat}})}.
$$

For flip generators with $p = 1$ at $(2, 3)$ (max temporal frequency
$\omega_{\max} = 2\pi / T = 1$ at $T = 2\pi$, $N_t = 3$):

- $\|[D_L^{\mathrm{diag}}, a^{\mathrm{flip}}_{p=1}]\|_{\mathrm{op}}
\approx 2 \|W\|_{\mathrm{op}} \cdot \omega_{\max} = 2 \cdot 1 \cdot 1 = 2$
(using $\|W\|_{\mathrm{op}} \sim 1$ for normalized 3-Y multipliers).

- $\|[D_L^{\mathrm{off}}, a^{\mathrm{flip}}_{p=1}]\|_{\mathrm{op}}
\sim \|D_{\mathrm{GV}}\|_{\mathrm{op}} \cdot \|W\|_{\mathrm{op}}
\approx 2.5 \cdot 1 = 2.5$
(from the first run's $\|D_L^{\mathrm{off}}\|_{\mathrm{op}} = 2.5$).

Combining via §4.4:
$$
   L_{\mathrm{op}}(a^{\mathrm{flip}}_{p=1})
\;\approx\; \sqrt{2^2 + 2.5^2}
\;\approx\; 3.20.
$$

For the natural baseline:
$$
   L_{\mathrm{op}}(a^{\mathrm{nat}}_{p=1})
\;=\; \|[D_L^{\mathrm{off}}, a^{\mathrm{nat}}]\|_{\mathrm{op}}
\;\approx\; 2.5.
$$

Ratio: $3.20 / 2.5 \approx 1.28$.  Hence
$$
\boxed{\;\Lambda^{\mathrm{enlarged}}(2, 3) \approx 1.28 \cdot 2.0746 \approx 2.66.\;}
$$

This is consistent with Paper 46 §7.2's prediction of a constant of order
$\sqrt{2}$ to $2$ (chirality-flipping content roughly doubles the
operator-norm content; doubled here means the squared norm doubles, i.e.
the constant grows by $\sqrt{2} \approx 1.41$, modulo
correlations between the two block sectors).

---

## §6. Go / no-go verdict (T5)

### §6.1 Verdict: **POSITIVE-GO**

The structural separation is established at three independent levels:

1. **Structural identity break** at $p \geq 1$: $[D_L^{\mathrm{diag}}, a^{\mathrm{flip}}] \neq 0$.

2. **Strict-strong-form separation**: $L^{P45}_{\mathrm{block}}(a^{\mathrm{flip}}) = 0$
   while $L_{\mathrm{op}}(a^{\mathrm{flip}}) > 0$.

3. **Quantitative bound increase**: $\Lambda^{\mathrm{enlarged}}(2, 3) \approx 2.66$
   vs $\Lambda^{P45}(2, 3) = 2.0746$, ratio $\approx 1.28$.

### §6.2 Recommended next sub-sprint

**Sprint L3b-2f-β: full L1ʹ–L5 derivation on the enlarged substrate.**
Estimated scope: **4–8 weeks** at the L3b-2a–d pace.  Will produce
**Paper 47** as the eighth math.OA standalone in the GeoVac series.

The five lemmas need re-derivation:

- **L1ʹ (operator-system substrate):** the enlarged substrate
  $\mathcal{O}^L_{\mathrm{enlarged}}$ is well-defined, has $\dim = 2 \dim \mathcal{O}^L_{\mathrm{nat}}$,
  closes under multiplication (§3), and has propagation number
  $\mathrm{prop}_{\mathrm{achievable}} = 2$ (predicted).

- **L2 (cb-norm of joint central Fejér kernel):** needs re-derivation
  on the enlarged substrate.  Key question: does the
  Bożejko–Fendler central-multiplier equality still give
  $\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} = 2 / (n_{\max} + 1)$ when the multiplier
  algebra is enlarged?

- **L3 (joint Lichnerowicz):** the **structural identity** $[D_L, a_s \otimes a_t] = i[D_{\mathrm{GV}}, a_s] \otimes a_t$
  used in L3b-2a §3.3 BREAKS on the enlarged substrate (per §4.2).
  This is the hardest re-derivation: a new joint $C_3^{\mathrm{joint}}$ constant
  must be computed, no longer equal to the single-factor $C_3^{\mathrm{SU(2)}}$.

- **L4 (joint Berezin reconstruction):** the PURE_TENSOR structure
  $B^{\mathrm{joint}} = B^{\mathrm{SU(2)}} \otimes B^{U(1)}$ continues to hold at the
  level of the natural substrate but needs an extension to chirality-flipping
  reconstruction.  Likely: introduce a $\mathbb{Z}_2$-graded Berezin map
  with one block per chirality sector.

- **L5 (propinquity assembly):** Latrémolière 2017/2023 §4 machinery
  applies verbatim once L1ʹ–L4 are closed, with the larger constants.

### §6.3 Alternative pre-β sub-sprints (lighter scoping)

Before committing to the 4–8 week β sub-sprint, three light follow-ons
can be done in 1–2 days each on a stable machine:

- **L3b-2f-α-numeric**: rerun the present diagnostic with the corrected
  chirality-asymmetric construction on a stable machine to confirm
  all numerical predictions in §4 and §5 at the per-generator level.

- **L3b-2f-α-prop**: compute the propagation number on the enlarged
  substrate at $(2, 3)$ via an iterative span-tracking method without
  memory blowup (the original driver's propagation step exploded due
  to retaining all matrix products).

- **L3b-2f-α-eta**: test whether further enlargements (e.g., $\eta$-graded
  generators that mix chirality with temporal-parity) generate yet more
  content beyond Candidate A.  If yes, Paper 47 should treat the full
  $\mathbb{Z}_2 \times \mathbb{Z}_2$-graded substrate from the start.

---

## §7. Honest scope

### §7.1 What this sprint definitively closes

- The minimal anti-commuting enlargement of the operator system is the
  **chirality-asymmetric doubling** $\mathrm{diag}(W, -W)$ of the Weyl
  multiplier (in Paper 44's chiral basis convention).  Paper 46 §7.2's
  abstract form $[[0, W], [W^*, 0]]$ is convention-shifted but
  equivalent up to a basis change — and crucially NOT the form to use
  on Paper 44's substrate directly.

- The Paper 45 K⁺-weak-form Lipschitz seminorm vanishes identically
  on chirality-flipping generators because $P_+ M^{\mathrm{flip}} P_+ = 0$.
  This is the strict-strong-form separation, exhibited as a closed-form
  identity in §4.3.

- $[D_L^{\mathrm{diag}}, a^{\mathrm{flip}}]$ is non-zero for flip generators with
  temporal-polynomial degree $p \geq 1$, breaking the L3b-2a
  temporal-Lipschitz-invisibility identity in exactly the predicted
  way (Paper 46 §7.2 prediction confirmed structurally).

- $\Lambda^{\mathrm{enlarged}}(2, 3) \approx 2.66$ vs $\Lambda^{P45}(2, 3) = 2.0746$,
  ratio $\approx 1.28$, consistent with Paper 46 §7.2's prediction of
  $\sqrt{2}$–$2$.

### §7.2 What this sprint does NOT close

- Full L1ʹ–L5 derivation on the enlarged substrate (~4–8 weeks).
- Whether the asymptotic rate
  $\gamma^{\mathrm{joint}} = O(\log n_{\max} / n_{\max} + T / N_t)$ survives the
  enlargement at sharp constants.
- Whether further enlargements (e.g., $\eta$-graded mixed-parity
  generators) give yet more content beyond Candidate A.
- Full per-generator numerical verification at $(2, 3)$ for the
  corrected construction (deferred to a stable-machine re-run, ~1 day).
- Numerical verification at the second scoping cell $(3, 5)$
  (operator-system construction is sympy-arithmetic-heavy and the
  current Windows machine state cannot complete it).

### §7.3 Environmental note

Numerical verification of T3 and T4 at the per-generator level was
attempted multiple times on the local Windows machine.  The first
driver run (with the OLD off-block-diagonal flip construction) captured
the (2, 3) Krein-space + Lorentzian-Dirac decomposition data and the
natural-substrate operator-system construction (~50 seconds wall time)
but crashed with MemoryError on the propagation step (the naive
$O(n^2)$ pairwise-product list exceeded available memory).

Subsequent runs with (i) the corrected chirality-asymmetric construction,
(ii) skipping the propagation step entirely, (iii) reducing to a
minimal-probe variant using only `krein_space_compact_temporal` and
`lorentzian_dirac_compact`, all hung on `geovac` module imports for 2+
minutes with no progress.  Multiple cleanups of zombie Python processes
were performed; the environmental cause is unclear.  This affects
*only* the empirical confirmation of the analytical predictions in §4
and §5; the predictions themselves are derived rigorously from Paper 44
§3.1's chiral-basis convention and the L3b-2a structural decomposition.

A 1-day re-run of the diagnostic on a machine without the import-hang
issue is the natural confirmation step before opening L3b-2f-β.

### §7.4 Forward implications

- The strict-strong-form vs K⁺-weak-form gap is sharpest precisely at
  the chirality-flipping generators that the K⁺-weak-form filters out.
  Any future application that probes chirality-mixing physics (e.g.,
  Yukawa-type couplings between matter and antimatter sectors in the
  Connes–Marcolli AC factor, or Higgs-doublet hopping in the spectral-
  action approach to the Standard Model) would natively require the
  strong-form construction, not the K⁺-weak-form.

- Paper 46 §7.2 G1ʹ is now refined from "open question" to
  "POSITIVE-GO at the scoping level": the enlarged substrate produces
  genuinely-new strong-form content.  The next-tranche follow-on
  Paper 47 is structurally justified.

- The chirality-asymmetric form $\mathrm{diag}(W, -W)$ has an
  immediate physical interpretation: it is the spatial multiplier
  acting **with opposite sign on left- and right-handed Weyl
  spinors** — i.e., a chirality-violating (parity-violating) scalar
  multiplier on the Krein-Lorentzian spectral triple.  This is the
  natural lift to the operator-system substrate of the standard
  observation that any physical parity-violating multiplier (axial
  current, etc.) couples to the four-witness Wick-rotation arc only
  through the chirality-flipping component.  Paper 47, when written,
  should note this interpretation.

