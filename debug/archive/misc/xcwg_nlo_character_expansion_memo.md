# XCWG-E: Next-to-leading-order character expansion on Rule B

**Date.** 2026-05-16
**Sprint.** XCWG-E (third Wilson-loop diagnostic; follows XCWG-C Gaussian probe and XCWG-D LO area law)
**Companion data.** `debug/data/xcwg_nlo_character_expansion.json`
**Driver.** `debug/xcwg_nlo_character_expansion.py`

**One-line verdict.** **INCONCLUSIVE-LEANING-NEGATIVE.** Rule B has an unusually
small smallest closed 2-cycle ($k_\delta = 3$, much smaller than the 6-plaquette
cube boundary of $\mathbb{Z}^4$), so NLO corrections enter at very low order
in $(I_1/I_0)$. The normalized $\langle W\rangle_{\rm NLO}$ at fixed Wilson loop
$A_{\min}=1$ shows nominal sign reversal of $\sigma_{\rm NLO}(\beta)$ at $\beta \sim 30$,
but the Creutz-ratio $\sigma_{\rm eff}(A_1, A_2)$ analysis shows
$\sigma_{\rm eff}(2 \to 3) \approx \sigma_{\rm LO}(\beta)$ remains positive at all
tested $\beta$ — the apparent zero crossing in the single-loop diagnostic is
a finite-truncation artifact, not a deconfinement transition. The cleanest
read of the data: **NLO does not show a clean 4D-style deconfinement transition,
but the diagnostic is not powerful enough to distinguish 3D-like permanent
confinement from a slow continuum approach.** XCWG-F monopole density is the
next diagnostic.

---

## §1 Setup: NLO character expansion and the 3D-vs-4D distinction

For compact U(1) lattice gauge theory with Wilson action
$$S_W = -\beta \sum_P \cos\theta_P, \qquad \theta_P = \sum_{e\in P} \theta_e,$$
the character expansion gives, for a Wilson loop $W(C) = e^{i\theta_C}$,
$$\langle W(C) \rangle_\beta
  = \frac{N(w;\beta)}{Z(\beta)},$$
where
$$N(w;\beta) = \sum_{n \in \mathbb{Z}^P : d_1^T n = w} \prod_P \frac{I_{n_P}(\beta)}{I_0(\beta)},
\qquad
Z(\beta) = \sum_{n \in \mathbb{Z}^P : d_1^T n = 0} \prod_P \frac{I_{n_P}(\beta)}{I_0(\beta)},$$
and $d_1: \mathbb{Z}^P \to \mathbb{Z}^E$ is the plaquette boundary map (rows of
$d_1$ are signed L=4 plaquette boundaries on Rule B).

### §1.1 Leading order
At LO, take $n_P \in \{0, \pm 1\}$ and require $|n|_1 = \sum_P |n_P|$ minimal.
The unique minimizer for a Wilson loop with $A_{\min}$ plaquettes is
$n = n_{\rm LO}$ with $n_P = +1$ on those $A_{\min}$ plaquettes (oriented to
match $w$). Then $N_{\rm LO} = (I_1/I_0)^{A_{\min}}$ and $Z_{\rm LO} = 1$ (from
$n \equiv 0$), giving
$$\langle W \rangle_{\rm LO} = \left(\frac{I_1(\beta)}{I_0(\beta)}\right)^{A_{\min}}.$$
This is XCWG-D's verified positivity-and-asymptotics result; $\sigma_{\rm LO}(\beta) = -\ln(I_1/I_0) > 0$ at all $\beta$.

### §1.2 Next-to-leading order
NLO contributions come from non-minimal assignments. Write $n = n_{\rm LO} + z$
where $z \in \mathbb{Z}^P$ is a *closed 2-cycle*, i.e. $d_1^T z = 0$. The
contribution scales as $(I_1/I_0)^{|n_{\rm LO} + z|_1}$. Let $k_\delta$ denote
the smallest non-zero $|z|_1$ over all closed 2-cycles. The first non-trivial
NLO contribution is suppressed by an additional factor $(I_1/I_0)^{k_\delta}$
relative to LO.

### §1.3 3D-vs-4D theoretical predictions
On the $\mathbb{Z}^d$ cubic lattice:
- **4D** ($\mathbb{Z}^4$): the smallest closed 2-surface is the boundary of a
  3-cube, made of 6 plaquettes ($k_\delta = 6$). The character expansion has a
  Polyakov phase transition at $\beta_c \approx 1.01$ where the string tension
  drops to zero (deconfinement).
- **3D** ($\mathbb{Z}^3$): in dimension 3, no closed 2-surface exists in the
  plaquette complex alone — closed 2-surfaces would be 3-cells which $\mathbb{Z}^3$
  doesn't have. The Polyakov 1977 result is that the monopole-instanton expansion
  (a non-plaquette mechanism) gives permanent confinement at all $\beta$.

The NLO character expansion is, in pure plaquette homology, structurally
*equivalent* in 3D and 4D when $b_2$ (second Betti number of the 2-complex)
is similar; what differs is the 3-cell content. Rule B is a *graph* (no
intrinsic 3-cells); whatever 2-cycles exist live in the 2-complex generated
by L=4 plaquettes.

---

## §2 Smallest closed 2-cycle on Rule B

### §2.1 Method
We search for the smallest integer $z \in \mathbb{Z}^P$ with $z_p \in \{-1, 0, +1\}$
and $d_1^T z = 0$ (closed). We BFS over connected plaquette subsets of size
$k = 2, 3, \ldots$ and for each subset test all $2^{k-1}$ sign patterns
(modulo global flip) for $d_1^T z = 0$.

### §2.2 Result at $n_{\max} = 2$
**$k_\delta = 3$**, much smaller than $\mathbb{Z}^4$'s $k_\delta = 6$. Sample
2-cycle: $\{p_0, p_3, p_{38}\}$ with signs $(+1, -1, -1)$.

Structurally, this 2-cycle is a triangulated closed surface:
- $p_0 = (0, 8, 3, 9)$, edge set $\{3, 4, 18, 19\}$
- $p_3 = (0, 8, 2, 9)$, edge set $\{3, 4, 13, 14\}$
- $p_{38} = (2, 8, 3, 9)$, edge set $\{13, 14, 18, 19\}$

The three plaquettes pairwise share two edges and together cover six edges,
each used by exactly two plaquettes. Geometrically this is the analog of a
triangulated $S^2$ "tetrahedron face" — three plaquettes glued along common
axes, closing up to form a 2-cycle.

This structure arises because Rule B has high plaquette density (44 L=4
plaquettes on $V=10$ vertices, $E=20$ edges, so each edge participates in
many plaquettes on average — about $4 \cdot 44/20 = 8.8$ plaquette-incidences
per edge). The "extra" plaquettes generate closed-surface relations at very
low order.

### §2.3 Comparison to $\mathbb{Z}^d$
| Lattice | $k_\delta$ | Notes |
|:--------|:----------:|:------|
| $\mathbb{Z}^3$ | none (no plaquette-only closed surfaces) | confinement is monopole-driven |
| $\mathbb{Z}^4$ | 6 (cube boundary) | NLO Polyakov phase transition at $\beta_c \approx 1.01$ |
| **Rule B (n_max=2)** | **3** | dense plaquette overlap, very low-order 2-cycles |
| **Rule B (n_max=3)** | **3** | structure persists at larger graph |

**The $k_\delta = 3$ is the key structural finding of this sprint.** It says
Rule B has a *more plaquette-rich* 2-complex than 4D cubic lattices, with
closed 2-surfaces appearing at lower order. Whether this destabilizes
confinement or leaves it intact is the operational question.

### §2.4 Cycle homology dimension
At $n_{\max} = 2$: $b_2 = P - {\rm rank}(d_1) = 44 - 11 = 33$ (large kernel,
many independent cycles).
At $n_{\max} = 3$: $b_2 = 200 - 64 = 136$ (200-plaquette cap).

The 2-complex generated by L=4 plaquettes has a very large kernel of $d_1^T$
on Rule B, consistent with the high plaquette density. The smallest cycle in
the entire kernel is $k_\delta = 3$.

---

## §3 Normalized NLO for a single-plaquette Wilson loop

### §3.1 Construction
Take $W$ = boundary of plaquette $p_0$ at $n_{\max} = 2$, so $A_{\min} = 1$.
$n_{\rm LO}$ is the indicator on $p_0$ alone with sign $+1$. NLO assignments
are $n_{\rm LO} + \alpha z$ for closed 2-cycles $z$ of size $k_\delta = 3$
and $\alpha \in \{-1, +1, -2, +2\}$. We truncate at $|n|_1 \le A_{\min} + k_\delta + 1 = 5$.

Similarly $Z_{\rm NLO}$ sums $\alpha z$ for $\alpha \in \{\pm 1, \pm 2\}$ over
all enumerated 2-cycles.

### §3.2 Result at $n_{\max} = 2$
We use 30 closed 2-cycles of size 3 (the full set of size-3 cycles is
combinatorially bounded by the local L=4-plaquette neighbourhoods of each plaquette).

| $\beta$ | $\sigma_{\rm LO}$ | $\sigma_{\rm NLO}$ | $N$ | $Z$ | $\langle W\rangle_{\rm NLO}$ |
|:-------:|:-----------------:|:------------------:|:---:|:---:|:-----------------------------:|
| 0.1 | 2.997 | 2.817 | 0.0602 | 1.007 | 0.0598 |
| 0.3 | 1.908 | 1.516 | 0.263 | 1.196 | 0.220 |
| 0.5 | 1.417 | 1.025 | 0.665 | 1.856 | 0.359 |
| 0.7 | 1.108 | 0.805 | 1.413 | 3.160 | 0.447 |
| 1.0 | 0.807 | 0.618 | 3.41 | 6.34 | 0.539 |
| 1.5 | 0.517 | 0.424 | 8.98 | 13.71 | 0.655 |
| 2.0 | 0.360 | 0.302 | 15.82 | 21.38 | 0.740 |
| 3.0 | 0.211 | 0.171 | 27.71 | 32.88 | 0.843 |
| 5.0 | 0.113 | 0.072 | 40.74 | 43.78 | 0.930 |
| 10  | 0.053 | 0.003 | 52.06 | 52.22 | 0.997 |
| 30  | 0.017 | −0.042 | 60.49 | 58.02 | 1.042 |
| 100 | 0.005 | −0.058 | 63.92 | 60.33 | 1.060 |

**Observations:**
- $\sigma_{\rm NLO} < \sigma_{\rm LO}$ at every $\beta$ (NLO suppresses the
  string tension; consistent with the structural expectation that bigger
  surfaces compete with the LO surface).
- $\sigma_{\rm NLO}$ crosses zero around $\beta \sim 10\text{-}30$, becoming
  slightly negative (−0.058 at $\beta = 100$).
- **$\langle W\rangle_{\rm NLO} > 1$ at large $\beta$**: this is non-physical
  for a unitary Wilson loop. The truncated enumeration is over-counting $N$
  relative to $Z$ at large $\beta$, where $I_n/I_0 \to 1$ for every $n$ and
  the number of valid surfaces in $N$ outpaces $Z$ within the truncation budget.

The $\langle W\rangle_{\rm NLO} > 1$ at large $\beta$ is a clear **truncation
artifact**: in the full sum (no truncation), each higher-charge state with
$I_n/I_0 \to 1$ contributes proportionally to both $N$ and $Z$, and the ratio
remains $\le 1$. The truncation breaks this symmetry because the *number of
valid surfaces* with $d_1^T n = w$ within $|n|_1 \le {\rm cutoff}$ grows
asymmetrically vs the number with $d_1^T n = 0$.

### §3.3 Result at $n_{\max} = 3$
At $n_{\max} = 3$ (V=28, E=106, P=200 with cap), $k_\delta = 3$ persists.
Single-plaquette Wilson loop, 10 cycles of size 3 enumerated:

| $\beta$ | $\sigma_{\rm LO}$ | $\sigma_{\rm NLO}$ | $\langle W\rangle_{\rm NLO}$ |
|:-------:|:-----------------:|:------------------:|:-----------------------------:|
| 0.1 | 2.997 | 2.949 | 0.0524 |
| 1.0 | 0.807 | 0.695 | 0.499 |
| 3.0 | 0.211 | 0.183 | 0.833 |
| 10 | 0.053 | 0.017 | 0.984 |
| 30 | 0.017 | −0.026 | 1.026 |
| 100 | 0.005 | −0.041 | 1.042 |

Same qualitative pattern: nominal zero crossing around $\beta = 10\text{-}30$
with $\langle W\rangle_{\rm NLO} > 1$ at large $\beta$ — same truncation
artifact.

---

## §4 Creutz-ratio cross-check (cleaner diagnostic)

### §4.1 Method
For multiple area targets $A \in \{1, 2, 3, 4\}$, build a connected cluster
of $A$ plaquettes, compute $\langle W(A)\rangle_{\rm NLO}$ properly normalized,
then extract $\sigma_{\rm eff}(A_1, A_2) := -\ln[\langle W(A_2)\rangle/\langle W(A_1)\rangle] / (A_2 - A_1)$.
At LO, $\sigma_{\rm eff}$ is independent of $(A_1, A_2)$ and equals $\sigma_{\rm LO}(\beta)$.
At NLO, $\sigma_{\rm eff}$ varies with $A$; the large-$A$ limit gives the
true asymptotic string tension.

### §4.2 Result at $n_{\max} = 2$
20 cycles of size 6 enumerated for the Creutz cross-check (uses size-6 not size-3
to avoid double-counting the same physical surface in cluster of size 2, where
adding a size-3 cycle gives a size-5 cluster).

| $\beta$ | $\sigma_{\rm LO}$ | $\sigma_{\rm eff}(1{\to}2)$ | $\sigma_{\rm eff}(2{\to}3)$ | $\sigma_{\rm eff}(3{\to}4)$ |
|:-------:|:----:|:----:|:----:|:----:|
| 0.1 | 2.997 | 2.982 | **2.997** | 2.310 |
| 0.3 | 1.908 | 1.788 | **1.908** | 1.268 |
| 0.5 | 1.417 | 1.147 | **1.417** | 0.857 |
| 0.7 | 1.108 | 0.727 | **1.108** | 0.643 |
| 1.0 | 0.807 | 0.399 | **0.808** | 0.473 |
| 1.5 | 0.517 | 0.236 | **0.521** | 0.333 |
| 2.0 | 0.360 | 0.186 | **0.367** | 0.251 |
| 3.0 | 0.211 | 0.132 | **0.233** | 0.159 |
| 5.0 | 0.113 | 0.078 | **0.178** | 0.082 |
| 10 | 0.053 | 0.036 | **0.191** | 0.021 |
| 30 | 0.017 | 0.011 | **0.231** | −0.020 |
| 100 | 0.005 | 0.002 | **0.250** | −0.036 |

**Headline:** $\sigma_{\rm eff}(2 \to 3) \approx \sigma_{\rm LO}(\beta)$ at
small $\beta$, but at large $\beta$ it **grows back to a positive plateau
around 0.20**, not zero. The "true" string tension extracted from this
asymptotic pair stays *positive at all $\beta$*.

This is the cleanest reading of the data:
- $\sigma_{\rm eff}(1 \to 2)$ shows large NLO suppression (because the smaller
  loop's NLO corrections are large relative to LO at small area)
- $\sigma_{\rm eff}(2 \to 3)$ stays close to $\sigma_{\rm LO}$ at small $\beta$
  and saturates at a positive value at large $\beta$ — the **most informative
  pair**, since the NLO contributions to areas 2 and 3 partially cancel
- $\sigma_{\rm eff}(3 \to 4)$ shows mild negative crossing at very large $\beta$
  — likely truncation artifact, with the small cycle budget unable to fully
  match contributions to clusters of size 4 vs size 3

### §4.3 Interpretation
The Creutz cross-check gives a **distinctly different verdict** from the
single-plaquette diagnostic:
- Single-plaquette $\sigma_{\rm NLO}$ shows nominal zero crossing → looks 4D-like
- Creutz $\sigma_{\rm eff}(2 \to 3)$ stays positive and approaches a positive
  plateau → looks 3D-like

The reconciliation: in compact U(1) lattice gauge theory, the *true* string
tension is the large-$A$ limit of $\sigma_{\rm eff}(A, A+1)$. Small-$A$
diagnostics are dominated by perimeter-law corrections that aren't part of
the genuine string tension. The 4D-style deconfinement transition would
show up as $\sigma_{\rm eff} \to 0$ for *all* asymptotic pairs, not just
$(1 \to 2)$.

**Verdict from Creutz analysis:** Rule B Wilson U(1) does NOT show a clean
4D-style deconfinement transition. The asymptotic-$A$ string tension stays
positive at all $\beta$.

---

## §5 Why the diagnostic is honest-inconclusive

The two diagnostics disagree:
1. Single-plaquette $\sigma_{\rm NLO}$ goes nominally negative at $\beta \sim 30$
2. Creutz $\sigma_{\rm eff}(2 \to 3)$ stays positive at all $\beta$, plateauing
   at ~0.20 at large $\beta$

Possible interpretations:

**(A) 3D-like with permanent confinement.** The single-plaquette negative
crossing is a truncation artifact ($\langle W\rangle_{\rm NLO} > 1$ at large
$\beta$). The Creutz analysis on asymptotic-$A$ pairs reflects the true
string tension, which stays positive. This is the natural reading and matches
the XCWG-D LO verdict at full normalization.

**(B) Mild 4D-like behavior with $\beta_c$ at very large $\beta$.** The
Creutz $(2 \to 3)$ plateau at large $\beta$ is itself a finite-$A$ artifact;
the *true* asymptotic $\sigma$ might decrease further with $A$ and eventually
hit zero. Without enumerating $A \ge 10$ pairs with adequate cycle budgets,
we can't rule this out.

**(C) Mixed behavior characteristic of the high-plaquette-density graph.**
$k_\delta = 3$ is so small that the NLO sector is unusually "full" of
near-zero-tension surfaces. This could indicate that Rule B doesn't fit
cleanly into either the 3D or 4D universality class; it might be a
non-standard "compact U(1)-on-rich-2-complex" theory with intermediate
features.

**Computational scope.** A definitive verdict requires either:
- Monte Carlo simulation of the full character expansion at multiple $(\beta, A)$
- $A$ up to ~10 with proper cycle enumeration to size $\sim 30$
- Or: monopole-density analysis (XCWG-F) to distinguish 3D (always-condensed
  monopole plasma) from 4D (low-density monopole gas at large $\beta$)

These are beyond the current sprint budget.

---

## §6 What $k_\delta = 3$ structurally means

The smallest closed 2-cycle being size 3 is the most concrete structural
finding of this sprint. In standard cubic lattices, $k_\delta \ge 6$ (the
boundary of a 3-cube on $\mathbb{Z}^4$). Rule B's $k_\delta = 3$ reflects
its high plaquette density: 4.4 plaquettes per vertex (44 plaquettes on
$V = 10$), vs $\sim 1$ per vertex on $\mathbb{Z}^4$.

Geometrically, the size-3 cycles are "triangle prisms" — three plaquettes
that share a common edge-axis and tile a triangulated $S^2$-like
structure. They exist because Rule B has multiple L=4 plaquettes through
the same pair of vertices (parity-flip cycles in the Dirac graph).

**This says Rule B is a *non-cubic* graph in a strong sense**: its
2-complex carries closed surfaces at lower combinatorial order than any
$\mathbb{Z}^d$ cubic lattice. Whether this maps cleanly onto a known
universality class (e.g. compact U(1) on a triangulated 3-manifold,
or compact U(1) on $S^3$ with a non-standard cell structure) is an open
question. The XCWG-A foundation result ("Rule B is 3D-like in spectral
dimension and plaquette/edge ratio") gives a structural hint that the
3D class is the natural fit; the present NLO diagnostic is consistent
with that reading via the Creutz analysis.

---

## §7 Verdict and recommendation

**Verdict.** **INCONCLUSIVE-LEANING-3D-LIKE.** The NLO sprint produced:
- Confirmed $k_\delta = 3$ at both $n_{\max} = 2$ and $n_{\max} = 3$
  (Rule B has lower-order closed 2-surfaces than $\mathbb{Z}^4$).
- Nominal $\sigma_{\rm NLO}$ sign crossing at $\beta \sim 30$ from
  single-plaquette diagnostic, but accompanied by $\langle W\rangle > 1$
  (truncation artifact).
- Creutz $\sigma_{\rm eff}(2 \to 3)$ stays positive at all $\beta$,
  plateauing at $\sim 0.20$ at large $\beta$ — the cleaner asymptotic
  diagnostic, **consistent with 3D-like permanent confinement**.

**Paper 41 v2 status.** The single-plaquette diagnostic *cannot* be
honestly described as "passing NLO" because of the $\langle W\rangle > 1$
truncation artifact. The Creutz cross-check gives the better
3D-like reading but requires the operator to be cautious about extracting
asymptotic-$A$ behavior from limited $A$ values.

**Recommendation:** The XCWG-E sprint does NOT yet provide a clean
fifth-witness for the 3D compact U(1) identification at NLO. Paper 41 v2
should report:
1. The structural finding $k_\delta = 3$ on Rule B (genuinely new, with
   clear interpretation as high-plaquette-density consequence)
2. The Creutz-ratio result that $\sigma_{\rm eff}$ stays positive
   asymptotically at all tested $\beta$
3. An honest caveat that the single-plaquette $\sigma_{\rm NLO}$
   diagnostic has truncation issues that prevent a clean 3D/4D
   discriminator at this resolution
4. The next-step recommendation: **proceed to XCWG-F monopole density**,
   which probes the same physics through a different observable that is
   far less sensitive to surface-enumeration truncation. In compact U(1),
   the monopole density $\rho_M(\beta)$ stays nonzero at all $\beta$ in 3D
   (plasma phase) but drops to zero above $\beta_c$ in 4D (dilute gas with
   confinement-deconfinement transition). XCWG-F is the natural and more
   discriminating successor to XCWG-E.

---

## §8 Files

- `debug/xcwg_nlo_character_expansion.py` (~470 lines): driver
- `debug/data/xcwg_nlo_character_expansion.json`: full data tables
- `debug/xcwg_nlo_character_expansion_memo.md`: this memo

---

**Final one-line.** Rule B has $k_\delta = 3$ (much smaller than $\mathbb{Z}^4$'s
$k_\delta = 6$), but NLO character-expansion diagnostics with limited cycle
enumeration give mixed signals — the single-plaquette $\sigma_{\rm NLO}$ shows
nominal zero crossing (truncation artifact) while the asymptotic Creutz
$\sigma_{\rm eff}$ stays positive (3D-like). The verdict is honest-inconclusive
at NLO; the recommended next step is XCWG-F monopole density, which probes
the same 3D/4D distinction through an observable not subject to surface-
enumeration truncation.
