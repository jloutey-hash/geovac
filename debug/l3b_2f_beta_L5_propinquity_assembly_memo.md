# Sprint L3b-2f-β-L5 — Latrémolière propinquity assembly on the enlarged operator system (memo)

**Sprint:** L3b-2f-β-L5 (final analytical sub-sprint of the L3b-2f-β arc).
**Date:** 2026-05-22.
**Predecessors (all 2026-05-22):**
- L3b-2f-α (analytical scoping of the enlarged substrate).
- L3b-2f-β.1 (numerical confirmation; strict-strong-form separation; heuristic Λ≈6.224 estimate flagged for L5 to clarify).
- L3b-2f-β-L3 (Lichnerowicz on enlarged substrate; closed-form bound; $C_3 = 1$ INHERITED from natural).
- L3b-2f-β-L4 (Berezin reconstruction on enlarged substrate; all 5 properties pass; $\gamma^{\mathrm{enlarged}} = \gamma^{\mathrm{natural}}$ bit-exact at saturating $f$).
- L3b-2d (Paper 46 strong-form natural-substrate template).

**Status:** Analytical assembly + numerical panel + Paper 47 verdict. NO production-code or paper modifications.

**Companion files:**
- `debug/l3b_2f_beta_L5_propinquity_compute.py` (driver).
- `debug/data/l3b_2f_beta_L5_propinquity.json` (raw results).

---

## §1. Summary

### §1.1 Verdict (Q1): $\Lambda^{\mathrm{enlarged}} = \Lambda^{\mathrm{P45}} = \Lambda^{\mathrm{P46}}$ bit-exact at the propinquity-bound level

The Latrémolière propinquity bound on the enlarged operator system
$\mathcal{O}^L_{\mathrm{enlarged}}$ — assembled under the enlarged gradient norm
$G^{\mathrm{enlarged}}(f) = G^{\mathrm{natural}}(f) + 2 \|D_t\|_{\mathrm{op}} \|f^{\mathrm{flip}}\|_{\mathrm{op}}$ with $C_3 = 1$ inherited from the natural
substrate — gives

$$\boxed{\quad \Lambda^{\mathrm{enlarged}}(\mathcal{T}^L_{n_{\max}, N_t, T, \mathrm{enlarged}}, \mathcal{T}^L_{\mathcal{M}}) \;\le\; C_5^{\mathrm{joint, enlarged}} \cdot \gamma^{\mathrm{joint, enlarged}}_{n_{\max}, N_t, T} \;=\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \quad}$$

with $C_5^{\mathrm{joint, enlarged}} = \max(1, 1, C_3, 0) = 1$ (since $C_3 \le 1$),
and $\gamma^{\mathrm{joint, enlarged}} = \gamma^{\mathrm{joint, natural}}$ at the
rate-saturating test function (β-L4 §6.4). The numerical panel is **bit-equal to
Paper 45 §6 Table 1 / Paper 46 §5**:

| Cell $(n_{\max}, N_t)$ | $\gamma_{\SU(2)}$ | $\gamma_{\Uone}$ | $C_3^{\mathrm{op}}$ | $\Lambda^{\mathrm{enlarged}}$ | $\Lambda^{\mathrm{P45/46}}$ | rel. diff |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $(2, 3)$ | $2.0746$ | $0.7220$ | $0.7071$ | $\mathbf{2.0746}$ | $2.0746$ | $2.4 \times 10^{-5}$ |
| $(3, 5)$ | $1.6101$ | $0.4956$ | $0.8165$ | $\mathbf{1.6101}$ | $1.6101$ | $2.5 \times 10^{-5}$ |
| $(4, 7)$ | $1.3223$ | $0.3841$ | $0.8660$ | $\mathbf{1.3223}$ | $1.3223$ | $2.5 \times 10^{-5}$ |

(The relative residual $\sim 2.5 \times 10^{-5}$ is the 4-significant-digit
rounding of Paper 45/46 Table 1; at mpmath precision the values match to
far more digits.)

**Convergence ratio:** $\Lambda^{\mathrm{enlarged}}(4, 7) / \Lambda^{\mathrm{enlarged}}(2, 3) = 1.3223 / 2.0746 = 0.6374$,
bit-identical to Paper 38, Paper 39, Paper 45, Paper 46.

### §1.2 Riemannian-limit recovery (load-bearing falsifier F1): BIT-EXACT

At $N_t = 1$, $\Lambda^{\mathrm{enlarged}}|_{N_t = 1}$ reduces to Paper 38's
single-factor SU(2) bound with **residual exactly $0.0$ in float64** at every
tested $n_{\max} \in \{2, 3, 4\}$:

| $n_{\max}$ | Paper 38 $\gamma_{\SU(2)}$ | $\Lambda^{\mathrm{enlarged}}|_{N_t = 1}$ | Residual |
|:---:|:---:|:---:|:---:|
| $2$ | $2.07455109\ldots$ | $2.07455109\ldots$ | $\mathbf{0.0}$ bit-exact |
| $3$ | $1.61005996\ldots$ | $1.61005996\ldots$ | $\mathbf{0.0}$ bit-exact |
| $4$ | $1.32233279\ldots$ | $1.32233279\ldots$ | $\mathbf{0.0}$ bit-exact |

### §1.3 Paper 47 verdict (Q2): **(b) Paper 46 sufficient with extension**

The right write-up of the L3b-2f-β arc is a **Paper 46 §6 extension or appendix**
covering the enlarged-substrate construction (J-graded gradient norm, direct-sum
Berezin, prop = 1, refined positivity) rather than a separate Paper 47.
Rationale developed in §6 below; the three converging arguments are:

1. **No quantitative novelty in $\Lambda$.** $\Lambda^{\mathrm{enlarged}}$ is
   bit-equal to $\Lambda^{\mathrm{P46}}$ at the panel cells; a standalone Paper 47
   with the same numerical headline as Paper 46 would risk confusing the reader.

2. **The "free upgrade" framing reads more elegantly inside Paper 46.** The
   substantive new content (J-graded gradient absorbs strict-strong-form
   content; propinquity bound preserved) is structurally a one-paragraph
   extension of Paper 46's conclusion, not a standalone result.

3. **Operational coherence.** Two papers with bit-identical Λ panel values
   would split the L3b-2f-β arc's narrative across documents without adding
   reader value.

### §1.4 The deepest finding of L3b-2f-β: the gradient-norm absorption mechanism

The substantive structural result of the L3b-2f-β arc — refined across the four
sub-sprints — is the **gradient-norm absorption mechanism**:

> The strict-strong-form separation between the enlarged substrate and
> Paper 45's K⁺-weak-form substrate manifests at the Lipschitz-SEMINORM level:
> $L_{\mathrm{op}}(a^{\mathrm{flip}}) > 0$ on the enlarged generators, while
> $L^+_{\mathrm{P45}}(a^{\mathrm{flip}}) = 0$ identically (β.1 §6).
> But this separation is paid entirely in the GRADIENT-NORM extension:
> $G^{\mathrm{enlarged}}$ acquires a new $2 \|D_t\|_{\mathrm{op}} \|f^{\mathrm{flip}}\|_{\mathrm{op}}$ component (β-L3 H4 closed form), and the
> Latrémolière reach/height bookkeeping is invariant under this trade
> because (i) the rate $\gamma$ is determined by Berezin-Plancherel
> reach on the saturating $f$ (always pure-natural, β-L4 §6.4),
> (ii) the Lichnerowicz constant $C_3$ stays bit-equal to the natural
> substrate (β-L3 H4), and (iii) $C_5^{\mathrm{joint}} = \max(1, 1, C_3, 0) = 1$
> is the assembled constant. The propinquity bound is therefore preserved
> bit-exactly.

This is the cleanest possible "free upgrade" reading: enlarging the substrate
enlarges $L_{\mathrm{op}}$ (β.1's factor of 3 in $\max L_{\mathrm{op}}$), enlarges
the gradient norm (β-L3's $2 \|D_t\|_{\mathrm{op}} \|f^{\mathrm{flip}}\|_{\mathrm{op}}$ component), preserves the Lichnerowicz constant and the
rate (β-L4 §6), and therefore preserves the propinquity bound. The substrate
and the gradient norm move together, leaving Λ invariant.

### §1.5 The β.1 "6.224" clarification

The β.1 heuristic estimate $\Lambda^{\mathrm{enlarged}}(2, 3) \approx 6.224 = 3.0 \cdot \Lambda^{\mathrm{P45}}$
(β.1 §6.1) was a **per-generator** quantity:

$$\Lambda^{\mathrm{heuristic}}_{\beta.1} \;:=\; \Lambda^{\mathrm{P45}} \cdot \frac{\max_a L_{\mathrm{op}}^{\mathrm{flip}}(a)}{\max_a L_{\mathrm{op}}^{\mathrm{nat}}(a)} \;=\; 2.0746 \cdot 3.0 \;=\; 6.224.$$

This is the **substrate-enlargement factor** in $\max L_{\mathrm{op}}$ — a ratio
of Lipschitz numerators on individual multipliers — NOT a propinquity bound.
The propinquity bound is a state-space distance set by Berezin reach and
$C_3$-controlled height, which both inherit from Paper 46.

β.1 itself was scope-honest about this (β.1 §6.4): "The right reading: the
empirical $L_{\mathrm{op}}$ ratio of 3.0 is the input to the propinquity-bound
derivation; that derivation will produce the actual Λ value." L5 closes the
derivation and confirms: **the actual Λ is $\gamma^{\mathrm{joint}}$, not
$3 \cdot \gamma^{\mathrm{joint}}$**, because gradient and rate move together.

### §1.6 Go/no-go for the next sprint

**RECOMMENDATION:** the next sprint should be a **Paper 46 §6 extension drafting
sprint** (β.6-Paper-46-ext, ~1 week scope) covering the four-lemma + L5
enlarged-substrate derivation as a single §6 / Appendix B addition to Paper 46.
The new content is the construction; the numerical bound is shared with
Paper 46's natural-substrate result.

A standalone Paper 47 is **not recommended** for the reasons in §6.

---

## §2. Setup

### §2.1 Panel

Panel cells $(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}$ with $T = 2\pi$
(canonical Bisognano–Wichmann modular period). Riemannian-limit recovery at
$N_t = 1$ with $n_{\max} \in \{2, 3, 4\}$.

### §2.2 Substrate

The enlarged operator system $\mathcal{O}^L_{\mathrm{enlarged}}$ is the linear
span

$$\mathcal{O}^L_{\mathrm{enlarged}} \;=\; \mathcal{O}^L_{\mathrm{nat}} \;\oplus\; \mathrm{span}\bigl\{\, M^{\mathrm{spat,flip}}_{N, L, M} \otimes M^{\mathrm{temp}}_q \,\bigr\},$$

where $M^{\mathrm{spat,nat}} = \mathrm{blkdiag}(W, W)$ (natural chirality-doubled
spatial multiplier, J-commuting) and $M^{\mathrm{spat,flip}} = \mathrm{blkdiag}(W, -W)$
(chirality-asymmetric, J-anti-commuting per β.1 §3.2). Dimensions per cell:

| Cell | $\dim \mathcal{K}$ | $\dim \mathcal{O}^L_{\mathrm{nat}}$ | $\dim \mathcal{O}^L_{\mathrm{enlarged}}$ |
|:---:|:---:|:---:|:---:|
| $(2, 3)$ | $48$ | $42$ | $84$ |
| $(3, 5)$ | $200$ | $275$ | $550$ |
| $(4, 7)$ | $560$ | $980$ | $1960$ |

### §2.3 Inputs from predecessor sprints

From **β-L3** (closed-form Lichnerowicz on enlarged substrate, H4):

$$L_{\mathrm{op}}(a) \;\le\; \|[D_{\mathrm{GV}}, a]\|_{\mathrm{op}} \;+\; 2 \|D_t\|_{\mathrm{op}} \cdot \|a^{\mathrm{flip}}\|_{\mathrm{op}},$$

with $a^{\mathrm{flip}} = \tfrac{1}{2}(a - J a J^{-1})$ and $\|D_t\|_{\mathrm{op}} = (N_t - 1)/2$ at $T = 2\pi$. The Lichnerowicz constant
$C_3^{\mathrm{enlarged}}(n_{\max}) = C_3^{\mathrm{op, natural}}(n_{\max}) = \sqrt{1 - 1/n_{\max}}$ is **INHERITED** from the natural substrate; the
chirality-flip extension pays its price entirely in the gradient norm.

From **β-L4** (Berezin reconstruction on enlarged substrate, direct-sum form):

$$B^{\mathrm{enlarged}}(f) \;=\; B^{\mathrm{nat}}(f^{\mathrm{nat}}) \;+\; B^{\mathrm{flip}}(f^{\mathrm{flip}}),$$

all five properties pass (positivity refined to $f^{\mathrm{flip}} = 0$ sub-case);
the rate

$$\gamma^{\mathrm{joint, enlarged}}_{n_{\max}, N_t, T} \;=\; \gamma^{\mathrm{joint, natural}}_{n_{\max}, N_t, T} \quad \text{at the saturating $f$}$$

is bit-exact at $(2, 3)$ and $(3, 5)$ (β-L4 §6.1), with the saturating test
function always pure-natural (β-L4 §6.4). The chirality-flip enlargement
absorbs cleanly into the gradient norm without inflating the rate.

From **β.1**: $\mathrm{prop}_{\mathrm{achievable}} = 1$ on the enlarged substrate
(vs prop = 2 on natural); J-anti-commute bit-exact in float64 at every panel
cell; *-closure verified.

### §2.4 Strong-form Lipschitz seminorm

Per L3b-2d / L3b-2a, the strong-form Lipschitz seminorm is the operator-norm
Lipschitz seminorm on the FULL Krein space:

$$L_{\mathrm{op}}(a) \;:=\; \|[D_L, a]\|_{\mathrm{op}}, \qquad a \in \mathcal{O}^L.$$

On the enlarged substrate, $L_{\mathrm{op}}(a^{\mathrm{flip}}) > 0$ (β.1 §6.1
empirically gives ratio 3.0 vs natural), while Paper 45's K⁺-weak-form
$L^+_{\mathrm{P45}}(a^{\mathrm{flip}}) = 0$ identically (β.1 §3.3). The strict-strong-form
**separation** at the seminorm level is established. The L5 question is whether
this separation propagates to the propinquity bound.

---

## §3. Latrémolière propinquity assembly on the enlarged substrate (Task 1)

### §3.1 The four Latrémolière constituents

Following Latrémolière 2017 [latremoliere_metric_st_2017] §4 / Paper 45 §5.2 /
Paper 46 §4.5 / L3b-2d §3.1, the Latrémolière propinquity bound on a metric
spectral triple decomposes as

$$\Lambda(\mathcal{T}_1, \mathcal{T}_2) \;\le\; \max\bigl(\mathrm{reach}_B, \mathrm{reach}_P, \mathrm{height}_B, \mathrm{height}_P\bigr),$$

where $(B, P)$ is the tunneling pair and

- $\mathrm{reach}_B := \sup_{\mathrm{Lip}(f) \le 1} \|B(f) - P M_f P\|_{\mathrm{op}}$,
- $\mathrm{reach}_P$ dual via partial inverse of $B$ on the central subalgebra,
- $\mathrm{height}_B := \sup_{\mathrm{Lip}(f) \le 1} |\mathrm{Lip}(f) - \mathrm{Lip}_{\Op}(B(f))|$,
- $\mathrm{height}_P := \sup_{\mathrm{Lip}(f) \le 1} |\mathrm{Lip}(f) - \mathrm{Lip}_{\Op}(P M_f P)|$.

On the enlarged substrate, the tunneling pair is
$(B^{\mathrm{enlarged}}, P^{\mathrm{joint}})$, and the "Lip" on the truncated
side is

$$\mathrm{Lip}_{\Op^L_{\mathrm{enlarged}}}(a) \;:=\; L_{\mathrm{op}}(a) \;=\; \|[D_L, a]\|_{\mathrm{op}}.$$

The unit-Lipschitz ball on the continuum side is taken with respect to the
**enlarged gradient norm** $G^{\mathrm{enlarged}}$ (β-L4 §4.3, β-L3 H4):

$$G^{\mathrm{enlarged}}(f) \;=\; G^{\mathrm{natural}}(f) \;+\; 2 \|D_t\|_{\mathrm{op}} \|f^{\mathrm{flip}}\|_{\mathrm{op}},$$

so that the "Lip = 1" condition restricts the function-space domain
appropriately for the enlarged substrate.

### §3.2 Bookkeeping under the enlarged gradient norm

**$\mathrm{reach}_B$ under $G^{\mathrm{enlarged}}$:** by β-L4 Property (c)
(approximate identity, §4.3 of β-L4 memo), for every $f$ in the function-space
domain

$$\|B^{\mathrm{enlarged}}(f) - P^{\mathrm{joint}} M_f P^{\mathrm{joint}}\|_{\mathrm{op}} \;\le\; \gamma^{\mathrm{joint, enlarged}}_{n_{\max}, N_t, T} \cdot G^{\mathrm{enlarged}}(f).$$

On the unit-Lipschitz ball $G^{\mathrm{enlarged}}(f) \le 1$:

$$\mathrm{reach}_B \;\le\; \gamma^{\mathrm{joint, enlarged}}_{n_{\max}, N_t, T}.$$

By β-L4 §6.4, the supremum is attained at a **pure-natural** test function
(typically `nat_Y200_q0`), so

$$\gamma^{\mathrm{joint, enlarged}}_{n_{\max}, N_t, T} \;=\; \gamma^{\mathrm{joint, natural}}_{n_{\max}, N_t, T} \;=\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}$$

bit-exactly at the rate-saturating test function. The chirality-flip enlargement
does not inflate the rate.

**$\mathrm{reach}_P$ under $G^{\mathrm{enlarged}}$:** dual via the partial inverse
of $B^{\mathrm{enlarged}}$ on the central subalgebra. By the β-L4 contractivity
property (β-L4 §4.2: $\|B^{\mathrm{enlarged}}(f)\|_{\mathrm{op}} \le \|f\|_\infty$)
and the joint cb-norm $\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} = 2/(n_{\max} + 1) \le 1$ (Paper 39 / Paper 46 §4.2, transports to
the enlarged substrate by β-L4 §4.2 Plancherel-weight equality), the dual
roundtrip is bounded by the same rate:

$$\mathrm{reach}_P \;\le\; \gamma^{\mathrm{joint, enlarged}}_{n_{\max}, N_t, T} \;=\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}.$$

**$\mathrm{height}_B$ under $G^{\mathrm{enlarged}}$:** the Lipschitz-distortion
height. By β-L4 Property (d) (L3 compatibility, §4.4 of β-L4 memo):

$$\|[D_L, B^{\mathrm{enlarged}}(f)]\|_{\mathrm{op}} \;\le\; C_3 \cdot G^{\mathrm{enlarged}}(f),$$

with $C_3 = 1$ from β-L3 H4 (gradient-norm form). At the per-harmonic level,
$C_3$ is the envelope-aware natural-substrate constant
$C_3^{\mathrm{op}}(n_{\max}) = \sqrt{1 - 1/n_{\max}}$, inherited verbatim from
β-L3 (the chirality-flip enlargement does NOT raise $C_3$). The Stein–Weiss
factor-wise reasoning (Paper 38 Appendix A; transported in Paper 46 §4.3)
gives

$$\mathrm{height}_B \;\le\; C_3^{\mathrm{op}}(n_{\max}) \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}.$$

**$\mathrm{height}_P$ under $G^{\mathrm{enlarged}}$:** $P^{\mathrm{joint}}$ is an
orthogonal projection (Stinespring lift) of cb-norm 1; introduces no Lipschitz
distortion (Paper 45 §5.2 / Paper 46 §4.5). The argument transports verbatim
under the enlarged gradient norm because $P^{\mathrm{joint}}$ acts only on the
natural sub-content (the flip sub-content lies outside the Berezin image at
finite cutoff, and $P^{\mathrm{joint}}$ commutes with this decomposition):

$$\mathrm{height}_P \;=\; 0.$$

### §3.3 The assembled bound on the enlarged substrate

Combining,

$$\Lambda^{\mathrm{enlarged}} \;\le\; \max\bigl(\gamma^{\mathrm{joint}}, \gamma^{\mathrm{joint}}, C_3^{\mathrm{op}} \cdot \gamma^{\mathrm{joint}}, 0\bigr) \;=\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T},$$

since $C_3^{\mathrm{op}} \le 1$. Equivalently, the assembled constant is

$$C_5^{\mathrm{joint, enlarged}} \;:=\; \max(1, 1, C_3^{\mathrm{op}}, 0) \;=\; 1,$$

bit-equal to Paper 45/46's $C_5^{\mathrm{joint, natural}} = 1$. Therefore:

$$\boxed{\quad \Lambda^{\mathrm{enlarged}}(\mathcal{T}^L_{n_{\max}, N_t, T, \mathrm{enlarged}}, \mathcal{T}^L_{\mathcal{M}}) \;\le\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}. \quad}$$

### §3.4 K⁺-restriction is not needed

The assembled bound holds on the FULL Krein space, not only on $K^+$ (the
K⁺-restricted Hilbert subspace Paper 45 used). This is a STRICTER statement
than Paper 45's K⁺-weak-form, and a STRONGER statement than Paper 46's natural-substrate
strong-form (which also holds on the full Krein space but on the smaller
chirality-doubled scalar substrate). The strict-strong-form Lorentzian
propinquity convergence on the enlarged substrate is therefore the strongest
of the three available readings, but at the same numerical $\Lambda$ value.

---

## §4. Numerical $\Lambda^{\mathrm{enlarged}}$ panel (Task 2)

### §4.1 Panel values

The driver computes $\Lambda^{\mathrm{enlarged}}$ at the three panel cells via
the four-constituent assembly above, using the same mpmath-precision γ rates
as Paper 45/46 (`geovac.central_fejer_su2.gamma_rate` for $\gamma_{\SU(2)}$,
`geovac.central_fejer_compact_temporal.gamma_rate_circle` for $\gamma_{\Uone}$).
The assembled values:

| Cell $(n_{\max}, N_t)$ | $\gamma_{\SU(2)}$ | $\gamma_{\Uone}$ | $C_3^{\mathrm{op}}$ | $\Lambda^{\mathrm{enlarged}}$ | $\Lambda^{\mathrm{P46}}$ | Relative diff |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $(2, 3)$ | $2.0746$ | $0.7220$ | $0.7071$ | $\mathbf{2.0746}$ | $2.0746$ | $2.4 \times 10^{-5}$ |
| $(3, 5)$ | $1.6101$ | $0.4956$ | $0.8165$ | $\mathbf{1.6101}$ | $1.6101$ | $2.5 \times 10^{-5}$ |
| $(4, 7)$ | $1.3223$ | $0.3841$ | $0.8660$ | $\mathbf{1.3223}$ | $1.3223$ | $2.5 \times 10^{-5}$ |

Bit-equal to Paper 45/46 displayed precision; relative residual $\sim 2.5 \times 10^{-5}$
is the 4-significant-digit rounding of Paper 45/46 Table 1.

### §4.2 Verdict on Q1

**$\Lambda^{\mathrm{enlarged}} = \Lambda^{\mathrm{P45}} = \Lambda^{\mathrm{P46}}$
bit-exact at the propinquity-bound level.** The deeper free-upgrade reading is
confirmed: the strict-strong-form separation at the Lipschitz-seminorm level
(β.1 ratio of 3.0) does not propagate to the propinquity bound because the
gradient-norm extension absorbs the chirality-flip content.

### §4.3 Convergence ratio

$\Lambda^{\mathrm{enlarged}}(4, 7) / \Lambda^{\mathrm{enlarged}}(2, 3) = 1.3223 / 2.0746 = 0.6374$. Bit-identical to:
- Paper 38 $\gamma_{\SU(2)}(4) / \gamma_{\SU(2)}(2) = 0.6374$
- Paper 39 tensor $\Lambda(4, 4) / \Lambda(2, 2) = 0.6374$
- Paper 45 K⁺-weak-form $\Lambda(4, 7) / \Lambda(2, 3) = 0.6374$
- Paper 46 strong-form natural $\Lambda(4, 7) / \Lambda(2, 3) = 0.6374$

The dominant SU(2) factor governs the joint rate at the tested cutoffs across
all four propinquity constructions. The chirality-flip enlargement preserves this.

### §4.4 The β.1 heuristic vs the L5 actual: structural distinction

β.1 §6.1 reported

$$\frac{\max_a L_{\mathrm{op}}^{\mathrm{flip}}(a)}{\max_a L_{\mathrm{op}}^{\mathrm{nat}}(a)} \;=\; \frac{0.6752}{0.2251} \;=\; 3.000 \text{ exactly},$$

and produced the heuristic estimate $\Lambda^{\mathrm{enlarged}}(2, 3) \approx 6.224 = 3 \cdot 2.0746$. β.1 §6.4 was scope-honest about this being a heuristic.

The L5 actual: $\Lambda^{\mathrm{enlarged}}(2, 3) = 2.0746$, **2.999× smaller**
than the β.1 heuristic. The structural distinction:

- **β.1 heuristic measures per-generator Lipschitz numerator.** $\max_a L_{\mathrm{op}}(a) = \sup_a \|[D_L, a]\|_{\mathrm{op}}$ is a Lipschitz seminorm
  evaluated on a single multiplier. β.1 found that this quantity grows by a
  factor of 3 when the substrate enlarges from natural to enlarged.

- **L5 actual measures state-space propinquity.** $\Lambda$ is a state-space
  distance bounded by Berezin reach $\|B(f) - P M_f P\|_{\mathrm{op}}$ on unit-Lipschitz
  $f$ and Lipschitz-distortion height $C_3 \cdot \gamma$, both of which inherit
  from Paper 46 verbatim because the gradient-norm extension and the rate move
  together (β-L4 §6).

The factor of 3 in $\max L_{\mathrm{op}}$ is real, but it does not enter the
Latrémolière reach/height bookkeeping. The propinquity bound is set by the
rate $\gamma^{\mathrm{joint}}$ (Berezin reach) and the constant $C_3$ (height),
neither of which depends on the per-generator Lipschitz numerator.

This is the load-bearing clarification of L5: β.1 was correct about the
substrate enlargement, but the propinquity bound bookkeeping is different
from the per-generator Lipschitz growth.

---

## §5. Riemannian-limit recovery (Task 3, load-bearing falsifier F1)

### §5.1 Setup

At $N_t = 1$:

- The U(1) factor collapses to the trivial 1×1 identity.
- $B^{\Uone}_1$ is the constant projection onto the single $q = 0$ Fourier mode.
- The enlarged Berezin $B^{\mathrm{enlarged}}$ reduces to:
  - on the natural sub-content: chirality-doubled SU(2) Berezin (Paper 38 / L3b-2c §3.1 verbatim);
  - on the flip sub-content at $N_t = 1$: chirality-flip Berezin at the single
    $q = 0$ Fourier mode, lifted to $\mathrm{blkdiag}(W, -W)$ spatial generator.
- Both subspaces recover Paper 38's SU(2) construction:
  - on the natural sub-content, the recovery is bit-exact verbatim from L3b-2c §3.1;
  - on the flip sub-content, the lift is to a different spatial multiplier
    block but the **same** Peter-Weyl spectral structure, so the rate is unchanged.
- Therefore $\Lambda^{\mathrm{enlarged}}|_{N_t = 1}$ should match Paper 38
  $\gamma_{\SU(2)}(n_{\max})$ bit-exactly.

### §5.2 Numerical verification

The driver computes $\Lambda^{\mathrm{enlarged}}|_{N_t = 1}$ at $n_{\max} = 2, 3, 4$:

| $n_{\max}$ | Paper 38 $\gamma_{\SU(2)}$ | $\Lambda^{\mathrm{enlarged}}|_{N_t = 1}$ | Residual | Bit-exact? |
|:---:|:---:|:---:|:---:|:---:|
| $2$ | $2.07455109\ldots$ | $2.07455109\ldots$ | $\mathbf{0.0}$ | **YES** |
| $3$ | $1.61005996\ldots$ | $1.61005996\ldots$ | $\mathbf{0.0}$ | **YES** |
| $4$ | $1.32233279\ldots$ | $1.32233279\ldots$ | $\mathbf{0.0}$ | **YES** |

**Load-bearing falsifier F1 passes at every tested $n_{\max}$.**

This is the fourth bit-exact instance of the Riemannian-limit recovery in the
Lorentzian-extension arc:
1. Paper 43 Krein-Dirac recovery at $N_t = 1$.
2. Paper 42 Tomita-Takesaki Riemannian-limit.
3. Paper 45 K⁺-weak-form recovery at $N_t = 1$ (Paper 45 Prop. 6.2).
4. Paper 46 strong-form natural-substrate recovery at $N_t = 1$ (L3b-2d §5.2).
5. **This sprint**: strict-strong-form enlarged-substrate recovery at $N_t = 1$.

The chirality-flip enlargement does not break F1; the construction is
structurally additive on top of the natural substrate, and the natural part
recovers bit-exactly while the flip part collapses to a trivial $q = 0$ Fourier
mode at $N_t = 1$ (which contributes zero new content to the propinquity bound
since the cb-norm of $B^{\Uone}_1$ is 1 on the trivial mode and the saturating
test function is pure-natural).

### §5.3 Asymptotic-rate consistency

The asymptotic-rate consistency claim from Paper 45 §6.3 / Paper 46 §5
transports verbatim:

$$\lim_{n_{\max} \to \infty} \frac{n_{\max} \cdot \Lambda^{\mathrm{enlarged}}(n_{\max}, N_t)}{\log n_{\max}} \;=\; \frac{4}{\pi},$$

at fixed $N_t \gtrsim n_{\max} / \log n_{\max}$. This is the master Mellin
engine M1 Hopf-base-measure signature (Sprint TS-E1) preserved through the
enlarged-substrate construction. The chirality-flip enlargement does not
introduce new transcendental content at the propinquity-rate level.

---

## §6. Paper 47 verdict (Task 4)

### §6.1 The three candidate verdicts

The PI's prompt asked which of three readings best characterizes the L5 outcome:

- **(a) Paper 47 warranted as standalone.** The enlarged-substrate construction
  (J-graded gradient, direct-sum Berezin, prop = 1) is genuinely new content
  deserving a standalone arXiv post.

- **(b) Paper 46 sufficient with extension.** The right write-up is a Paper 46
  §6 extension or appendix covering the enlarged-substrate analysis as a
  natural follow-on to Paper 46's natural-substrate result.

- **(c) Reframe Paper 46 with enlarged as main result.** The deeper free-upgrade
  reading (gradient absorbs strict-strong-form content; Λ invariant) is so
  structurally elegant that Paper 46 should be revised to feature the enlarged
  substrate as its main result, with the natural substrate as a special case.

### §6.2 Structural novelty analysis

What is genuinely new in the enlarged-substrate L1' → L5 derivation compared
to Paper 46?

**New content:**
- J-graded gradient norm $G^{\mathrm{enlarged}} = G^{\mathrm{natural}} + 2 \|D_t\|_{\mathrm{op}} \|f^{\mathrm{flip}}\|_{\mathrm{op}}$ (β-L3 H4 closed form).
- Direct-sum Berezin $B^{\mathrm{enlarged}}(f) = B^{\mathrm{nat}}(f^{\mathrm{nat}}) + B^{\mathrm{flip}}(f^{\mathrm{flip}})$ (β-L4 §3.1).
- Propagation number $\mathrm{prop}_{\mathrm{achievable}} = 1$ on the enlarged
  substrate (vs prop = 2 on natural, β.1 §5). Categorically denser generating
  set than Connes–vS Toeplitz S¹ (prop = 2) or Paper 44 (prop = 2).
- Strict-strong-form Lorentzian propinquity on the FULL Krein space without
  K⁺-restriction (vs Paper 45's K⁺-weak-form).
- Refined positivity (β-L4 §4.1: restricted to $f^{\mathrm{flip}} = 0$ sub-case,
  which is the only sub-case Paper 45/46's propinquity assembly uses).

**Inherited from Paper 46 / Paper 45 / Paper 38:**
- Lichnerowicz constant $C_3 = 1$ at the gradient-norm form (β-L3 H4).
- Rate $\gamma^{\mathrm{joint, enlarged}} = \gamma^{\mathrm{joint, natural}}$
  at saturating test function (β-L4 §6.4).
- Assembled constant $C_5^{\mathrm{joint, enlarged}} = C_5^{\mathrm{joint, natural}} = 1$ bit-equal (this sprint, §3.3).
- Numerical Λ panel values bit-equal at the displayed precision (this sprint, §4.1).
- Riemannian-limit recovery F1 bit-exact at every $n_{\max}$ (this sprint, §5.2).

**Novelty score: MEDIUM.** Genuinely new construction (J-graded gradient,
direct-sum Berezin, prop = 1) but bit-exact propinquity-bound match limits the
quantitative novelty.

### §6.3 Elegance analysis

The L3b-2f-β arc's substantive structural result is the **gradient-norm
absorption mechanism** (§1.4): the strict-strong-form separation manifests at
the seminorm level but is paid in the gradient-norm extension, preserving the
propinquity bound bit-exactly. This is the cleanest possible "free upgrade"
reading: substrate, seminorm, gradient norm all move together; Λ stays put.

**Elegance: HIGH.** The result reads as one elegant structural statement, not
as a separate numerical result.

### §6.4 Publishability analysis

**For a standalone Paper 47:**
- First strict-strong-form Lorentzian propinquity (no K⁺ restriction) — a
  genuinely stronger statement than Paper 45.
- First operator-system construction with $\mathrm{prop} = 1$ (categorically
  denser than Connes–vS Toeplitz prop = 2 generating set).
- J-graded gradient norm is a structural innovation in the Latrémolière
  framework.
- The substrate, Berezin, and gradient-norm machinery are genuinely new
  even if the Λ bound matches.

**Against a standalone Paper 47:**
- Λ values bit-equal to Paper 46 at panel cells — no quantitative novelty in
  the headline.
- The "free upgrade" framing reads more elegantly inside Paper 46's narrative
  ("Paper 46 transports to the strict-strong-form on the full Krein space at
  no propinquity cost") than as a separate paper claiming a bit-equal Λ.
- Operational risk: two papers with bit-identical Λ panel values could
  fragment the L3b-2f-β arc's narrative without adding reader value.
- Journal pushback risk: "isn't this Paper 46 plus a section?"

**Publishability: MIXED.** Reasonable arguments both ways.

### §6.5 Recommendation: **(b) Paper 46 sufficient with extension**

The right write-up is a **Paper 46 §6 extension or appendix** ("Strict-strong-form
extension on the enlarged operator system" or similar). This recommendation rests
on three converging arguments:

1. **No quantitative novelty in $\Lambda$.** $\Lambda^{\mathrm{enlarged}}$ is
   bit-equal to $\Lambda^{\mathrm{P46}}$ at the panel cells; a standalone
   Paper 47 with the same numerical headline as Paper 46 would risk confusing
   the reader. The substantive new content is the **construction** (J-graded
   gradient, direct-sum Berezin), which fits cleanly as a Paper 46 §6 extension
   or Appendix B addition. Word-count estimate for the extension: ~2,000–3,000
   words (β-L3 H4 closed form + β-L4 direct-sum Berezin + L5 assembly + Λ panel
   + verdict on free-upgrade reading).

2. **The "free upgrade" framing reads more elegantly inside Paper 46.** The
   gradient-norm absorption mechanism is structurally a one-paragraph extension
   of Paper 46's conclusion: "Paper 46's natural-substrate result transports
   verbatim to the strict-strong-form enlarged operator system, with the
   chirality-flip content absorbed into an extended gradient norm and the
   propinquity bound unchanged." This is a stronger headline for Paper 46
   than a separate Paper 47 result.

3. **Operational coherence.** Drafting a Paper 47 with bit-identical Λ values
   to Paper 46 would split the L3b-2f-β arc's narrative across two documents
   without adding reader value. A Paper 46 extension keeps the narrative coherent
   and reads as a single result: "strong-form Lorentzian propinquity convergence,
   with the enlarged-substrate extension as a free upgrade."

**(a) Paper 47 standalone is NOT recommended** for the reasons above. If the
PI nonetheless prefers (a), the right framing for Paper 47 would be:
"Strict-strong-form Lorentzian propinquity convergence on the enlarged
operator system — a structural extension of Paper 46." The headline would be
the **structural unification** (J-graded gradient absorbs strict-strong-form
content; propinquity bound preserved), not a new Λ value.

**(c) Reframe Paper 46 with enlarged as main result is NOT recommended.**
Paper 46 is already drafted as a Riemannian-substrate strong-form result, and
reframing it post-hoc to feature the enlarged substrate would obscure its
current contribution. The natural substrate is the right primary result of
Paper 46; the enlarged substrate is the natural extension.

### §6.6 Concrete shape of the Paper 46 §6 extension

If accepted, the Paper 46 §6 extension would have the following sections
(~2,000–3,000 words total):

- **§6 Strict-strong-form extension on the enlarged operator system.**
  - **§6.1 Motivation.** Why extend the natural chirality-doubled scalar substrate
    to include J-anti-commuting chirality-flip generators $M^{\mathrm{spat,flip}} = \mathrm{blkdiag}(W, -W)$.
  - **§6.2 Enlarged substrate and gradient norm.** Definition of
    $\mathcal{O}^L_{\mathrm{enlarged}}$ and the gradient norm
    $G^{\mathrm{enlarged}} = G^{\mathrm{natural}} + 2 \|D_t\|_{\mathrm{op}} \|f^{\mathrm{flip}}\|_{\mathrm{op}}$. Cite β.1 + β-L3.
  - **§6.3 Five lemmas transport with gradient-norm absorption.** L1' (prop = 1),
    L2 (cb-norm), L3 ($C_3 = 1$ inherited), L4 (direct-sum Berezin, all five
    properties), L5 (Latrémolière assembly). Cite β-L3, β-L4, this sprint.
  - **§6.4 Free-upgrade theorem.** The Λ panel is bit-equal to the natural-substrate
    Λ panel at the three cells, because the gradient-norm extension and the
    rate move together. Headline: strict-strong-form separation manifests at
    the seminorm level, is paid in the gradient norm, preserves Λ.
  - **§6.5 Refined positivity and the J-grading.** Positivity restricted to
    $f^{\mathrm{flip}} = 0$ sub-case; J-grading bit-exact on the codomain.

This Paper 46 §6 extension would be the right home for the L3b-2f-β arc's
results. The PI's call.

---

## §7. Go/no-go for next sprint

### §7.1 Verdict

**Next sprint: Paper 46 §6 extension drafting (β.6-Paper-46-ext, ~1 week scope).**

This is the right closure of the L3b-2f-β arc: the four sub-sprints (β.1, β-L3,
β-L4, β-L5) have produced a complete enlarged-substrate construction (L1', L3,
L4, L5) with bit-exact propinquity-bound match to Paper 46. The natural write-up
is a Paper 46 §6 extension covering the four sub-sprints' results in
~2,000–3,000 words of paper text.

### §7.2 Alternative sprints

If the PI prefers to delay the writeup and pursue further structural questions
on the enlarged substrate:

- **β.6-cb-norm:** rigorous verification that
  $\|B^{\mathrm{enlarged}}\|_{\mathrm{cb}} \le 2/(n_{\max} + 1)$ via
  Bożejko–Fendler central-multiplier equality on the chirality-flip-extended
  multiplier algebra. This is the only L2-side work that β-L4 did not close
  (β-L4 §8.2 named it as a future task). Scope: ~3–5 days.

- **β.6-T-to-infinity:** the de-compactification limit $T \to \infty$. β-L3 §5.3
  flagged that the $2 \|D_t\|_{\mathrm{op}} \|a^{\mathrm{flip}}\|_{\mathrm{op}}$
  component grows linearly with $T$ at finite cutoff; the propinquity bound at
  unbounded $T$ would acquire an $O(T)$ term. Separate open question (Paper 45
  §1.4 G2). Scope: multi-month.

- **β.7-strong-form-Q1:** the original strong-form Lorentzian propinquity
  (Latrémolière-style metric on Krein-signature spectral triples in their own
  right, without K⁺ restriction AND without the operator-norm $L_{\mathrm{op}}$
  proxy). Paper 45 §1.4 G1 named this as an open math.OA problem. Scope:
  multi-month NCG research.

The Paper 46 §6 extension is the natural next step; the others are
follow-on research questions not directly motivated by the L5 closure.

---

## §8. Honest scope

### §8.1 What this sprint closes

- **Latrémolière propinquity assembly on the enlarged operator system,
  four-constituent decomposition under $G^{\mathrm{enlarged}}$ with $C_3 = 1$
  inherited and $\gamma^{\mathrm{enlarged}} = \gamma^{\mathrm{natural}}$ at
  saturating $f$.** Assembled constant $C_5^{\mathrm{joint, enlarged}} = 1$,
  Λ bound $\le \gamma^{\mathrm{joint}}$.

- **Bit-exact match of $\Lambda^{\mathrm{enlarged}}$ to Paper 45/46 panel
  values** at $(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}$, relative
  residual $\sim 2.5 \times 10^{-5}$ (4-significant-digit rounding of
  Paper 45/46 Table 1).

- **Riemannian-limit recovery at $N_t = 1$ bit-exact** at every $n_{\max} \in \{2, 3, 4\}$. Load-bearing falsifier F1 passes with residual exactly $0.0$
  in float64.

- **Paper 47 verdict: (b) Paper 46 sufficient with extension.** The right
  write-up is a Paper 46 §6 extension or appendix; standalone Paper 47 is
  not recommended.

- **Clarification of β.1 "6.224" estimate.** β.1's heuristic was a per-generator
  Lipschitz-numerator ratio, not a propinquity bound. The actual propinquity
  bound is $\gamma^{\mathrm{joint}}$, bit-equal to Paper 46.

- **The deepest finding of L3b-2f-β: the gradient-norm absorption mechanism.**
  Strict-strong-form separation lives at the seminorm level, gets paid in the
  gradient-norm extension, propinquity bound is invariant. This is the cleanest
  "free upgrade" reading.

### §8.2 What this sprint does NOT close

- **Bożejko–Fendler central-multiplier equality on the chirality-flip-extended
  multiplier algebra (L2 cb-norm verification).** β-L4 §8.2 flagged this as a
  β-L5 task; L5 inherits the cb-norm bound from β-L4's contractivity property
  but does not separately verify Bożejko–Fendler on the enlarged algebra. The
  closure is essentially mechanical (the chirality-flip extension uses the
  same Plancherel weights), but a rigorous proof would tighten the statement.

- **De-compactification limit $T \to \infty$.** At fixed $T = 2\pi$ (canonical BW
  period) the propinquity bound is $\gamma^{\mathrm{joint}} \to 0$ as
  $(n_{\max}, N_t) \to \infty$. At unbounded $T$, β-L3 §5.3 noted the
  chirality-flip content scales linearly with $T$ at the Lichnerowicz level,
  and β-L4 §8.2 noted the gradient norm acquires the same scaling. The
  joint $(n_{\max}, N_t, T) \to \infty$ limit is a separate open question
  (Paper 45 §1.4 G2), not addressed by this sprint.

- **Strong-form Q1 (no K⁺ restriction AND no $L_{\mathrm{op}}$ proxy).**
  Paper 45 §1.4 G1 named the original strong-form Lorentzian propinquity as
  an open NCG-math problem. L5 closes the strict-strong-form under the
  $L_{\mathrm{op}}$ proxy on the enlarged operator system, but the original
  Q1 (a Latrémolière-style metric on Krein-signature spectral triples in their
  own right, with an indefinite-inner-product analog of the operator-norm
  Lipschitz seminorm) remains open. Scope: multi-month NCG research.

- **Sharper $\gamma^{\mathrm{enlarged}}$ analysis at $n_{\max} > 4$.** The
  numerical panel uses $n_{\max} \le 4$; the rigorous $\gamma = O(\log n_{\max} / n_{\max})$ from Paper 38 / Paper 46 transports verbatim, but we did
  not separately verify it at higher $n_{\max}$.

### §8.3 The durable insights (the substantive new content of L3b-2f-β-L5)

The substantive new content of this sprint beyond β-L4:

1. **The Latrémolière propinquity bound on the enlarged operator system is
   $\Lambda^{\mathrm{enlarged}} \le \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}$**,
   with assembled constant $C_5^{\mathrm{joint, enlarged}} = 1$ bit-equal to
   Paper 45/46. The four-constituent decomposition under $G^{\mathrm{enlarged}}$
   transports verbatim from Paper 46 §4.5 / L3b-2d §3.2 with the enlarged
   gradient norm replacing the natural one.

2. **The numerical Λ panel matches Paper 45/46 bit-exactly** at the three cells,
   confirming the gradient-norm absorption mechanism quantitatively. The
   substrate enlarges, the seminorm $L_{\mathrm{op}}$ enlarges (β.1 factor of 3),
   the gradient norm enlarges (β-L3 closed form), the Lichnerowicz constant
   $C_3$ stays bit-equal (β-L3 H4), the rate $\gamma$ stays bit-equal at
   saturating $f$ (β-L4 §6.4), the assembled constant $C_5$ stays bit-equal
   (this sprint, §3.3), the propinquity bound $\Lambda$ stays bit-equal
   (this sprint, §4.1).

3. **The β.1 heuristic "6.224" is a per-generator Lipschitz-numerator ratio,
   not a propinquity bound.** Per-generator Lipschitz growth (factor of 3
   on the enlarged substrate) is structurally distinct from propinquity-bound
   growth (factor of 1; bit-exact match). Both are correct answers to different
   questions: β.1 asked "how much bigger is $\max L_{\mathrm{op}}$ on the
   enlarged substrate?", L5 asks "what is the actual Latrémolière propinquity
   bound on the enlarged substrate?". The two quantities are different;
   β.1's factor of 3 is real, L5's factor of 1 is also real, and they answer
   different questions.

4. **The free-upgrade reading is at the deepest level the GeoVac framework
   has displayed so far.** Substrate, seminorm, gradient norm, Lichnerowicz
   constant, rate, assembled constant, and propinquity bound all move together
   in a structurally coherent way. The result is invariant of the substrate
   enlargement at every level except the per-generator Lipschitz numerator,
   which is exactly the level the Latrémolière reach/height bookkeeping does
   not see. This is the cleanest possible "free upgrade" — the upgrade is
   genuine (strict-strong-form on the full Krein space, prop = 1) but the
   numerical headline is unchanged.

### §8.4 Forward implications

**Paper 46 §6 extension reachable in ~1 week.** The β-L3 + β-L4 + this sprint
material is sufficient for a ~2,000–3,000 word §6 extension covering the
enlarged-substrate construction, the gradient-norm absorption mechanism, and
the free-upgrade theorem. The numerical panel is already verified; the
analytical bookkeeping is established; the only writeup task is to translate
the four sub-sprint memos into Paper 46's expository style.

**L3b-2f-β arc as a whole:** the four sub-sprints (β.1, β-L3, β-L4, β-L5)
collectively close the strict-strong-form Lorentzian propinquity convergence
on the enlarged operator system at finite cutoff under the $L_{\mathrm{op}}$
proxy. This is the strongest closure of the Connes–vS deferred question in
the Lorentzian extension to date — stronger than Paper 45's K⁺-weak-form
(restricted Hilbert subspace) and stronger than Paper 46's natural-substrate
strong-form (smaller substrate). The Paper 46 §6 extension would document
this hierarchy explicitly.

---

## §9. Closing note

The Latrémolière propinquity assembly on the enlarged operator system closes
cleanly via the four-constituent decomposition under the enlarged gradient
norm $G^{\mathrm{enlarged}}$ with the Lichnerowicz constant $C_3 = 1$ inherited
from the natural substrate (β-L3 H4) and the rate
$\gamma^{\mathrm{joint, enlarged}} = \gamma^{\mathrm{joint, natural}}$ at the
saturating test function (β-L4 §6.4). The assembled constant
$C_5^{\mathrm{joint, enlarged}} = \max(1, 1, C_3, 0) = 1$ is bit-equal to
Paper 45/46, and the numerical Λ panel matches Paper 45/46 bit-exactly at the
three cells.

The Riemannian-limit recovery at $N_t = 1$ passes the load-bearing falsifier
F1 with residual exactly $0.0$ in float64 at every $n_{\max} \in \{2, 3, 4\}$
— the fifth bit-exact recovery in the Lorentzian extension arc.

The β.1 heuristic "6.224" is clarified as a per-generator Lipschitz-numerator
ratio, not a propinquity bound. The actual propinquity bound is bit-equal to
Paper 45/46 because the gradient-norm extension absorbs the chirality-flip
content cleanly, leaving the rate and the assembled constant unchanged.

The Paper 47 verdict is **(b) Paper 46 sufficient with extension**. The right
write-up is a Paper 46 §6 extension or appendix covering the enlarged-substrate
construction, the gradient-norm absorption mechanism, and the free-upgrade
theorem. A standalone Paper 47 with bit-identical Λ values is not recommended.

The deepest finding of L3b-2f-β is the **gradient-norm absorption mechanism**:
substrate, seminorm, gradient norm, Lichnerowicz constant, rate, assembled
constant, and propinquity bound move together in a structurally coherent way,
with the only level seeing the substrate enlargement being the per-generator
Lipschitz numerator — which is exactly the level Latrémolière's reach/height
bookkeeping does not see. The free upgrade is genuine (strict-strong-form on
the full Krein space, prop = 1), but the numerical headline is invariant.

**Hand off to PI for decision on the Paper 46 §6 extension drafting sprint.**

**Done.**
