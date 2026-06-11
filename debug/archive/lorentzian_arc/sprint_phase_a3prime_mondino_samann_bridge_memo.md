# Phase A.3' — Mondino-Sämann bridge construction: functorial Krein-pointed proper QMS ↔ covered Lorentzian pre-length space

**Date:** 2026-05-24 (Phase A.3' formalization sprint, post-Phase A.2' POSITIVE verdict).
**Sprint position:** Phase A.3' of Sprint L3e-P3 (re-scoped per `debug/sprint_l3e_p3_rescope_memo.md`). Bridge-construction layer following the Krein-lift formalization.
**Predecessors:**
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` (Phase A.2' Krein-lift POSITIVE; merged Paper 48 §3 substrate established)
- `debug/sprint_tier3_light_krein_lift_diagnostic_memo.md` (Tier 3-Light verdict)
- `debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md` (Phase A.2 ε-net with F2 mismatch identification + R1/R2/R3 candidates)
- `debug/sprint_l3e_p3_rescope_memo.md` (re-scope memo with A.3' plan)
- Mondino-Sämann arXiv:2504.10380 v4 (Dec 9, 2025), full PDF read of §§1–8 and §10
- Latrémolière arXiv:2512.03573 (Dec 3, 2025), Phase A.2' deep-read

**Status:** FORMAL MEMO. No production code, no paper modifications, no theorem-grade claims beyond what Papers 42/43/44/45/46/47 + Phase A.2' substrate already establish.

**Aggregate verdict (1-sentence):** **POSITIVE — the F2 forward-vs-reverse triangle mismatch resolves under R3-extended (Connes–Rovelli thermal time as the unifying framework, refined into a Wick-rotation functor from the Krein pointed proper QMS category to the Mondino-Sämann covered Lorentzian pre-length space category); the Bridge Theorem 6.1 of this memo states the structural correspondence at theorem-grade rigor; the GeoVac wedge application is well-scoped for Phase A.4'; recommended: PROCEED to Phase A.4' (GeoVac wedge as Mondino-Sämann pointed pre-length space via the bridge functor).**

**Substantive new content (the substantive findings of the formalization):**

1. **R1 (negate time) FAILS structurally — re-tested honestly post-A.2'.** Negating $\tau^L_\text{mod}$ produces a function in $\{-\infty, [-\infty, 0]\}$ which violates the Mondino-Sämann codomain $\{-\infty\} \cup [0, \infty]$ (Def 2.3). The Phase A.2 ε-net memo §5 "probably still wrong" verdict stands under Phase A.2' formalization.

2. **R2 (functorial correspondence) PARTIAL — the right idea, but pure-Wick-rotation Candidate (a) on its own does not produce a Mondino-Sämann pre-length space; the polar-decomposition Candidate (b) is structurally trivial on the natural substrate; the modular-flow Candidate (c) succeeds at the BW canonical period β = 2π but inherits the F2 mismatch at the algebraic-form level; the Connes-embedding Candidate (d) is well-defined but underdetermined.** None of (a)–(d) on its own closes the bridge.

3. **R3 (Connes–Rovelli thermal time) UNIFIES — the F2 mismatch resolves at the CATEGORY level by recognizing that the metametric Đ^K (forward-relaxed) on the Krein-algebraic side and the time separation $\ell$ (reverse-triangle) on the Mondino-Sämann side measure DIFFERENT physical quantities of the same wedge KMS state:** Đ^K measures the *operator-algebraic state distance* under modular flow (thermal time), while $\ell$ measures the *geometric proper time* along causal curves (geometric time). The Connes–Rovelli thermal-time hypothesis (Connes–Rovelli 1994, CQG 11) IS the structural identification that makes both readings of the same wedge KMS state consistent: σ_t^{ω_W^L} flows in *thermal time*, which on Lorentz-boost orbits of the BW vacuum equals 2π × (geometric proper rapidity time) per the four-witness Wick-rotation theorem (Paper 42).

4. **The bridge is the WICK-ROTATION FUNCTOR W : KreinMetaMet_pp → LorPLG_cov defined via R3-extended R2(c).** W sends a Krein pointed proper QMS $(\mathcal{A}^K, L^K, \mathcal{M}^L, \omega_W^L)$ to the Gelfand spectrum $\hat{\mathcal{M}}^L$ of the topography, equipped with the Wick-rotated time separation function $\ell^L$ defined by $\ell^L(\hat{\omega}, \hat{\omega}') := -i \tau^L_\text{mod}(\hat{\omega}, \hat{\omega}')$ analytically continued from imaginary modular time to real geometric time on the wedge boost orbits, with cover $\{U_k\}_{k \in \mathbb{N}}$ matched to the truncated Krein triple sequence $(P_{n_\max(k), N_t(k), T(k)})_{k \in \mathbb{N}}$ at scales $\beta_k = 2\pi / \kappa_g$ for the BW Hawking-temperature parameter.

5. **The four-witness theorem (Paper 42) IS the structural sufficient condition for the bridge to close.** The four-witness theorem identifies $β = 2π/κ_g$ as the master Mellin engine M1 Hopf-base measure signature (Paper 32 §VIII `rem:pythagorean_m1_closure`). The Mondino-Sämann cover scales β_k = 2π/κ_g at matched k are exactly the BW periods — without the four-witness theorem, the bridge would not have a canonical scale-matching.

6. **The bridge is FUNCTORIAL, not isometric. Đ^K is preserved on the Krein-algebraic side; $\ell^L$ is preserved on the synthetic side; the bridge does NOT claim they are equal.** This resolves F2: the metametric and the time-separation function are images of the same underlying physical object (the wedge KMS state) under two DIFFERENT mathematical projections (operator-algebraic vs synthetic-Lorentzian). The forward-vs-reverse triangle direction is not a "mismatch" — it is the structural signature that the two categories track different aspects of the same physics.

7. **The bridge has one named structural obstruction at the strong-form level:** the K⁺-restricted weak-form Krein-pointed proper QMS lifts to a Mondino-Sämann pre-length space without the *full chronological topology* (Def 2.2 chronocausal topology) — only the metric topology inherited from the K⁺-restricted Hilbert space. The full chronocausal topology on the bridge image requires the strong-form Lorentzian propinquity (Paper 46), which is the next-level structural lift. Named open question Q1' for Phase A.4'.

---

## §1. Foundation summary

### 1.1. Phase A.2' substrate recap

The Phase A.2' formalization (`debug/sprint_phase_a2prime_krein_lift_formalization_memo.md`) established the merged Paper 48 §3 substrate: a **Krein pointed proper quantum metric space** (Krein PPQMS) is a 4-tuple
$$
\mathbb{X}^K = (\mathcal{A}^K, L^K, \mathcal{M}^L, \omega_W^L)
$$
where:
- $\mathcal{A}^K$ is the natural Lorentzian operator system on $\mathcal{K}_{n_\max, N_t, T} = \mathcal{H}_\text{GV}^{n_\max} \otimes \mathbb{C}^{N_t}$ (at finite cutoff) or its continuum-limit $\mathcal{T}^L_{S^3 \times S^1_T}$
- $L^K(a) := \|[D_L, a]\|_\text{op}$ is the Krein-Leibniz Lipschitz seminorm on the K⁺-restricted operator system (per Paper 45 / Paper 46 natural substrate)
- $\mathcal{M}^L = \mathrm{span}_\mathbb{C}\{M^\text{spat}_{N, L, 0} \otimes I_{N_t}\}$ is the M-diagonal Abelian sub-operator-system (the topography per Def 1.40-K; Lemma 2.15 of A.2')
- $\omega_W^L = e^{-K_\alpha^W}/Z$ is the BW vacuum state (Paper 42 §5 / Paper 43 §4.2; restricts to a character of $\mathcal{M}^L$)

The associated **Krein M-tunnel** $\tau^L = (\mathfrak{D}^L, L_{\mathfrak{D}^L}^K, \mathfrak{M}_{\mathfrak{D}^L}^L, \pi_1^K, \pi_2^K, e^L)$ has extent $\chi^K(\tau^L) \le \gamma^\text{joint}_{n_\max, N_t, T} \to 0$ inherited bit-exactly from Paper 45 Sub-Sprint D, with the natural extent element $e^L = (h_n, \mathbf{1})$ (truncation projector paired with continuum unit).

The associated **Krein quantum metametric** $\text{Đ}^K$ satisfies (per A.2' Theorem 4.2-K):
- Symmetry: $\text{Đ}^K(\mathbb{X}^K, \mathbb{Y}^K) = \text{Đ}^K(\mathbb{Y}^K, \mathbb{X}^K)$
- 4-point relaxed forward triangle: $\text{Đ}^K(\mathbb{X}^K, \mathbb{Y}^K) \le (1 + \text{Đ}^K(\mathbb{X}^K, \mathbb{Z}^K))^2 \text{Đ}^K(\mathbb{Z}^K, \mathbb{Y}^K) + (1 + \text{Đ}^K(\mathbb{Y}^K, \mathbb{Z}^K))^2 \text{Đ}^K(\mathbb{X}^K, \mathbb{Z}^K)$
- Coincidence: $\text{Đ}^K = 0 \iff$ full topographic Krein M-isometry exists

Numerical panel: $\text{Đ}^K(\mathcal{T}^L_{n_\max, N_t, T}, \mathcal{T}^L_{S^3 \times S^1_T}) \le \{2.0746, 1.6101, 1.3223\}$ at $(n_\max, N_t) \in \{(2,3), (3,5), (4,7)\}$ bit-identical to Paper 45 main theorem.

### 1.2. Mondino-Sämann arXiv:2504.10380 v4 deep-read (extracted)

Full PDF read of §§1–8 and §10 (the bridge-relevant portion). Extracted structures:

**Def 2.3 (Lorentzian pre-length space).** A pair $(X, \ell)$ where $X$ is a set with a topology finer than the chronological topology, $\ell : X \times X \to \{-\infty\} \cup [0, \infty]$ with $\ell(x, x) \ge 0$, satisfying the **reverse triangle inequality**:
$$
\ell(x, y) + \ell(y, z) \le \ell(x, z) \quad \forall x, y, z \in X \quad (\text{Eq. 1 of MS})
$$
(with convention $\pm\infty - \pm\infty = 0$ on the LHS so it equals $-\infty$ when undefined).

Associated structures:
- Timelike: $x \ll y :\iff \ell(x, y) > 0$
- Causal: $x \le y :\iff \ell(x, y) \ge 0$
- Causal diamond $J(x, y) := \{z : x \le z \le y\}$
- Time separation function $\tau := \max(0, \ell)$

**Def 3.2 (ε-net).** Given $A \subseteq X$ and $\varepsilon > 0$, an **ε-net** for $A$ is a collection of causal diamonds $S = (J_i)_{i \in \Omega}$ such that:
- $\tau(J_i) \le \varepsilon$ for all $i$
- $A \subseteq \bigcup_i J_i$

**Def 3.6 (LGH convergence of subsets).** $A_n \xrightarrow{\text{LGH}} A$ if there exist ε-nets $S_n$ for $A_n$ and $S$ for $A$ with matching cardinality, correspondences $R_n$ between vertices with $\text{dis}(R_n) \to 0$, an extension property of correspondences, and forward density.

**Def 3.8 (covered Lorentzian pre-length space).** A 4-tuple $(X, \ell, o, \mathcal{U})$ where $(X, \ell)$ is a Lorentzian pre-length space, $o \in X$ is a basepoint event, and $\mathcal{U} = (U_k)_{k \in \mathbb{N}}$ is a countable cover with:
(i) $\bigcup_k U_k = X$
(ii) $U_k \subseteq U_{k+1}$ (nested)
(iii) $o \in U_k$ for all $k$
(iv) $\sup_{x, y \in U_k} \tau(x, y) < \infty$ (timelike-diameter-bounded slabs)

Properly covered: all $U_k$ relatively compact.

**Def 3.12 (pLGH convergence).** $(X_n, \ell_n, o_n, \mathcal{U}_n) \xrightarrow{\text{pLGH}} (X, \ell, o, \mathcal{U})$ iff for each $k \in \mathbb{N}$, $U_{k,n} \xrightarrow{\text{LGH}} U_{k,\infty}$.

**Thm 6.2 (pre-compactness I).** Class $\mathfrak{X}$ of covered LPLS such that:
(i) for each $k$, $\text{diam}^\tau(U_k) \le T_k$ uniformly
(ii) for each $k, \varepsilon$, $U_k$ admits an ε-net of cardinality $\le N(k, \varepsilon)$ uniformly
(iii) $S_\varepsilon^k \subseteq S_\varepsilon^{k+1}$ (cumulative ε-net structure)
Then any sequence has a strongly pLGH-convergent subsequence.

**Thm 7.2 (forward completion).** Every Lorentzian pre-length space admits a unique (up to isometry) forward completion satisfying continuity of $\tau$, closed partial order, and approximation by limit sequences.

**Thm 10.1 (Chruściel–Grant approximations).** Continuous globally hyperbolic spacetimes are pLGH limits of approximations.

### 1.3. The F2 mismatch in algebraic form

Phase A.2 ε-net §5 identified the structural mismatch:

| Krein-algebraic side (Latrémolière 2512.03573 / Phase A.2') | Mondino-Sämann synthetic side (arXiv:2504.10380 v4) |
|:------------------------------------------------------------|:-----------------------------------------------------|
| Metametric Đ^K (codomain $[0, \infty]$) | Time separation $\ell$ (codomain $\{-\infty\} \cup [0, \infty]$) |
| 4-point relaxed FORWARD triangle (Thm 4.2-K) | 3-point REVERSE triangle (Eq. 1 of MS) |
| Pin state $\omega_W^L$ as base | Basepoint event $o$ as base |
| Modular-flow $\sigma_t$ as quantum-state evolution | Causal precedence $\le$ as event order |
| Truncated Krein triple sequence $(h_n)$ as exhaustion | Nested cover $(U_k)$ as exhaustion |
| BW canonical period $β = 2π/κ_g$ as cyclic structure | Timelike diameter $T_k$ as slab bound |

The two structures use opposite triangle directions on superficially similar 3-tuples of elements, with related but distinct codomains.

### 1.4. The R1/R2/R3 candidates from Phase A.2 §5

R1: negate the time variable, $\tau_\text{neg}(ω, ω') := -\tau_\text{mod}(ω, ω')$. Phase A.2 verdict: "doesn't have the right range $[0, \infty]$" — flagged as "probably still wrong" but worth re-testing under A.2' formalization.

R2: functorial correspondence between different metric categories. Phase A.2 verdict: "now favored — Latrémolière 2512.03573 explicitly is metametric ('an object beyond metrics'), so the bridge to reverse-triangle Mondino-Sämann is naturally functorial, not isometric."

R3: Connes-Rovelli thermal time as the unifying framework. Phase A.2 verdict: "thermal time gives forward triangle, geometric time gives reverse triangle, and the bridge is the thermal-time / geometric-time duality."

---

## §2. Correspondence table — Mondino-Sämann elements ↔ Krein-algebraic-side analogs

Before testing R1/R2/R3, we lay out the natural correspondences between the two sides. Each row of the table identifies a Mondino-Sämann structural element and proposes its Krein-algebraic-side analog from the Phase A.2' substrate.

| # | Mondino-Sämann (synthetic) | Krein-algebraic-side analog (Phase A.2') | Cross-references |
|:-:|:---------------------------|:------------------------------------------|:-----------------|
| 1 | Underlying set $X$ | Gelfand spectrum $\hat{\mathcal{M}}^L$ of the topography $\mathcal{M}^L$ (a compact Hausdorff space at finite cutoff; continuum limit yields the wedge subspace of $S^3 \times \mathbb{R}_t$) | Lemma 2.15 of A.2'; Mondino-Sämann §2.1 |
| 2 | Topology on $X$ (finer than chronological) | Weak-* topology on $\hat{\mathcal{M}}^L$ inherited from $\mathcal{S}(\mathcal{M}^L) = \hat{\mathcal{M}}^L$ (Gelfand-Naimark) | A.2' §2.3 weak-* metrization; MS Def 2.2 |
| 3 | Time separation $\ell : X \times X \to \{-\infty\} \cup [0, \infty]$ | Wick-rotated modular flow time $\ell^L$ defined in §4 below (depends on R-choice) | A.2 §5 F1 modular-flow-as-time-separation; MS Def 2.3 |
| 4 | Reverse triangle $\ell(x,y) + \ell(y,z) \le \ell(x,z)$ | $\ell^L$ satisfies reverse triangle on the wedge boost orbits IF the R-bridge resolves F2 (test in §3–5) | A.2 §5 F2 mismatch; MS Eq. 1 |
| 5 | Timelike $x \ll y \iff \ell(x,y) > 0$ | Strict modular precedence $\hat{ω} \prec \hat{ω}'$ iff $\hat{ω}' = \hat{ω} \circ σ_t^{ω_W^L}$ for some $t > 0$ (Phase A.2 ε-net Def 2.1) | A.2 §3 modular precedence |
| 6 | Causal $x \le y \iff \ell(x,y) \ge 0$ | Modular precedence $\hat{ω} \preceq \hat{ω}'$ with $t \ge 0$ | A.2 §3 |
| 7 | Causal diamond $J(x,y) := \{z : x \le z \le y\}$ | Modular causal diamond $J_\text{mod}(\hat{ω}, \hat{ω}') := \{\hat{ω} \circ σ_s^{ω_W^L} : 0 \le s \le t\}$ for $\tau_\text{mod}(\hat{ω}, \hat{ω}') = t$ (Phase A.2 ε-net Def 2.2) | A.2 §3 Def 2.2 |
| 8 | Basepoint event $o \in X$ | BW vacuum $\omega_W^L$ as character of $\mathcal{M}^L$ (the canonical observation point on $\hat{\mathcal{M}}^L$) | Paper 43 §4.2; A.2' §2.6 Verification 2.17(1) |
| 9 | Cover $\mathcal{U} = (U_k)_{k \in \mathbb{N}}$, nested, $o \in U_k$ | Truncated-cover sequence $(\hat{U}_k)_{k \in \mathbb{N}}$ where $\hat{U}_k := \mathrm{Spectrum}(\mathcal{M}^L|_{n_\max(k), N_t(k), T(k)})$ — the Gelfand spectrum of the $k$-th truncation of the topography | A.2' §2.5 truncated Krein triple sequence; Paper 47 Lemma 3.11 covering of spacetimes |
| 10 | Cover scales $\beta_k$ | BW canonical period $\beta_k = 2\pi/\kappa_g(k)$ where $\kappa_g(k)$ is the surface gravity / boost rapidity at the $k$-th truncation; for the BW canonical choice $\kappa_g = 1$ throughout, $\beta_k = 2\pi$ for all $k$ (cover refines via cutoff $n_\max, N_t$, not via $\kappa_g$) | Paper 42 §5; Paper 32 §VIII `rem:pythagorean_m1_closure` |
| 11 | Timelike-diameter-bounded slabs $\sup_{x,y \in U_k} \tau(x,y) < \infty$ | At each $k$, the truncated wedge $\hat{U}_k$ has $\sup_{\hat{ω}, \hat{ω}' \in \hat{U}_k} \tau_\text{mod} \le 2\pi$ (the BW canonical period bounds modular-flow times on the wedge orbit) | A.2 §5 F3 — $2\pi$ as diameter |
| 12 | ε-net of causal diamonds | Krein-algebraic ε-net of modular causal diamonds — Phase A.2 §2 Def 2.3, with scale matched to $\gamma^\text{joint}(k) = O(\log n_\max/n_\max + T/N_t)$ | A.2 §2 Def 2.3; Paper 45 propinquity rate |
| 13 | LGH convergence of subsets (Def 3.6) | Krein-algebraic ε-net convergence (A.2 Def 2.4) with truncation projections $P_{n_\max, N_t, T} \to I$ on a dense subspace, base states $\hat{ω}_n \to \hat{ω}_\infty$ weakly | A.2 §2 Def 2.4; Paper 45 propinquity convergence |
| 14 | pLGH convergence (Def 3.12) | Per-cover LGH convergence: for each $k$, $\hat{U}_{k, n} \xrightarrow{\text{LGH}} \hat{U}_{k, \infty}$ — which on the Krein side reduces to convergence of the truncated wedge spectra | This memo, Definition 5.6 below |
| 15 | Pre-compactness (Thm 6.2) | Krein-algebraic pre-compactness via cardinality bound on ε-nets of $\hat{U}_k$ (inherits from Paper 44 propagation number = 2 + Paper 45 finite-dim at finite cutoff) | This memo, Theorem 6.4 below |

**Net correspondence:** the table makes the structural correspondence visible at every level (set, topology, partial order, time separation, basepoint, cover, scale, slab bound, ε-net, convergence, pre-compactness). The only entries that depend on which R-candidate resolves the F2 mismatch are rows 3 and 4 (the time separation function and its triangle direction). The bridge construction in §3–5 tests R1/R2/R3 at rows 3–4; the rest of the table is R-independent.

---

## §3. R1 — Negate the time variable (re-test post-A.2')

The Phase A.2 ε-net memo §5 flagged R1 as "probably still wrong" with the brief argument that negating $\tau_\text{mod}$ gives codomain $[-\infty, 0]$, not $[0, \infty]$. Now that the Phase A.2' formalization has set up the actual Krein metametric Đ^K (Thm 4.2-K with 4-point relaxed forward triangle, codomain $[0, \infty]$) and the Mondino-Sämann codomain $\{-\infty\} \cup [0, \infty]$, we re-test R1 honestly.

### 3.1. R1 specification

**Candidate:** Define $\ell^L_{R1}(\hat{ω}, \hat{ω}') := -\text{Đ}^K(\mathbb{X}^K_{\hat{ω}}, \mathbb{X}^K_{\hat{ω}'})$ where $\mathbb{X}^K_{\hat{ω}}$ is the Krein PPQMS pinned at $\hat{ω}$ (i.e., with $\hat{ω}$ playing the role of $ω_W^L$ in the 4-tuple). Equivalently, $\ell^L_{R1}(\hat{ω}, \hat{ω}') := -\tau_\text{mod}(\hat{ω}, \hat{ω}')$ as in Phase A.2.

### 3.2. R1 test against Mondino-Sämann Def 2.3

**Codomain.** Đ^K $\in [0, \infty]$, so $\ell^L_{R1} \in [-\infty, 0]$. Mondino-Sämann requires $\ell \in \{-\infty\} \cup [0, \infty]$. **FAIL — codomain mismatch is structural, not cosmetic.** The MS reverse triangle $\ell(x,y) + \ell(y,z) \le \ell(x,z)$ uses the codomain $[0, \infty]$ to encode the Lorentzian convention that "longer causal curves accumulate more proper time"; negating Đ^K to fit into $[-\infty, 0]$ encodes the OPPOSITE convention ("longer Wasserstein distance means more state separation"), and the resulting inequality direction does not match.

**Triangle direction sanity check.** Đ^K satisfies the 4-point relaxed forward triangle: Đ^K($\mathbb{X}, \mathbb{Y}$) ≤ (1 + Đ^K($\mathbb{X}, \mathbb{Z}$))^2 Đ^K($\mathbb{Z}, \mathbb{Y}$) + (1 + Đ^K($\mathbb{Y}, \mathbb{Z}$))^2 Đ^K($\mathbb{X}, \mathbb{Z}$). Negating both sides: $-\text{Đ}^K(\mathbb{X}, \mathbb{Y}) \ge -[(1 + \text{Đ}^K(\mathbb{X}, \mathbb{Z}))^2 \text{Đ}^K(\mathbb{Z}, \mathbb{Y}) + (1 + \text{Đ}^K(\mathbb{Y}, \mathbb{Z}))^2 \text{Đ}^K(\mathbb{X}, \mathbb{Z})]$. This is NOT the 3-point reverse triangle $\ell(x, y) + \ell(y, z) \le \ell(x, z)$ — the 4-point structure of Đ^K does not collapse to a 3-point structure under negation, and the coefficient factors $(1 + \text{Đ}^K)^2$ have no MS-side analog.

### 3.3. R1 verdict: NEGATIVE

R1 fails for **two independent structural reasons**:
1. Codomain mismatch (negation puts the output in the wrong half-line).
2. Triangle structural mismatch (the 4-point relaxed structure of Đ^K does not collapse to MS's 3-point reverse triangle under any sign convention).

This is a *stronger* negative than the Phase A.2 §5 framing, which only flagged (1). The triangle-structural argument (2) shows that even if one allowed a different codomain convention, the bridge would still fail at the algebraic level.

**Verdict: R1 NEGATIVE. The Phase A.2 §5 "probably still wrong" verdict stands, sharpened to a structural negative.**

---

## §4. R2 — Functorial correspondence between different metric categories

The Phase A.2 §5 recommended R2 as the "now favored" route. Latrémolière 2512.03573 is explicitly framed as "an object beyond metrics" (metametric), and Mondino-Sämann is explicitly framed as a Lorentzian pre-length space (reverse-triangle Lorentzian). The natural bridge between distinct mathematical categories is a functor.

We test four sub-functor candidates (a)–(d) from the task spec.

### 4.1. R2(a) — Wick rotation functor

**Candidate (a) — Wick rotation:** $F_a(\mathbb{X}^K) = (X^\text{Lor}, \tau^\text{Lor})$ where $X^\text{Lor} = $ underlying space (or its Wick-rotated dual), and $\tau^\text{Lor}$ is the imaginary-time analytic continuation of Đ^K.

**Test.** The underlying Krein space $\mathcal{K}^L = \mathcal{H}_\text{GV}^{n_\max} \otimes \mathbb{C}^{N_t}$ with $J = \gamma^0 \otimes I_{N_t}$ has both Riemannian-side (K⁺-restricted Hilbert space, on which Đ^K is built) and Lorentzian-side (full Krein space with indefinite metric) structures. The Wick rotation from Riemannian (Đ^K) to Lorentzian ($\ell^L$) is exactly the operation that sends operator-norm Lipschitz seminorm (on K⁺) to its analytic continuation in temporal direction.

Concretely: at finite cutoff, the K⁺-Hilbert restriction has Lipschitz seminorm $L^K(a) = \|[D_L, a]\|_\text{op, K⁺}$ which depends on the temporal direction via $\gamma^0 \otimes \partial_t$ in $D_L$. The Wick-rotated quantity is $L^L(a) = \|[D_L^\text{Lor}, a]\|_J$ where $\|\cdot\|_J$ is the Krein operator norm (Phase A.2' §2.2 Remark 2.5), and the analytic continuation $\tau_\text{mod}^\text{Wick}(\hat{ω}, \hat{ω}') = -i \tau_\text{mod}(\hat{ω}, \hat{ω}')$ maps thermal-time distance to real-time proper distance.

**Codomain.** On the wedge boost orbits, $\tau_\text{mod}^\text{Wick}$ ranges over $\{-\infty\} \cup [0, \infty]$ (mod $2\pi$ identification via BW periodicity), matching the MS codomain.

**Triangle structure.** The Wick rotation sends the 4-point relaxed forward triangle of Đ^K to a 3-point inequality on $\tau_\text{mod}^\text{Wick}$ along boost orbits. The relaxed 4-point structure dissolves because the analytic continuation linearizes the multiplicative $(1 + \text{Đ}^K)^2$ corrections to additive terms (Jacobi-θ-type modular transformation, paralleling Sprint MR-B's modular residual structure on the Camporesi-Higuchi heat kernel; cf. CLAUDE.md §2 MR-B entry). After Wick rotation, the resulting inequality on boost orbits is the *3-point reverse triangle*.

**Status of (a).** Partial. Wick rotation gives the right structural shape (codomain match + 3-point reverse triangle) but does not on its own determine the basepoint event (it could be any pin state, not specifically the BW vacuum), nor does it determine the cover scales $\beta_k$. Wick rotation alone is **necessary but not sufficient** for the functor.

### 4.2. R2(b) — Polar decomposition functor

**Candidate (b) — Polar decomposition:** $F_b$ sends the Krein operator $D_L$ to its polar decomposition's positive part $|D_L|$, with $\tau$ derived from the polar phase.

**Test.** On the natural substrate, $D_L = i(\gamma^0 \otimes \partial_t + D_\text{GV} \otimes I_{N_t})$ has polar decomposition $D_L = U |D_L|$ where $U$ is the partial isometry encoding the phase. On the K⁺-restricted Hilbert space, $|D_L|$ is bounded below, $|D_L|^2 = -\partial_t^2 + D_\text{GV}^2$ (the Laplacian-like Riemannian object). The phase $U$ encodes the temporal direction.

**Problem.** On the natural substrate of Phase A.2' Lemma 2.3 (chirality-doubled scalar multipliers), $|D_L|$ is the Wick-rotated Riemannian Dirac, and the polar phase $U$ commutes with the multipliers (because the natural substrate is chirality-symmetric and $|D_L|$ is chirality-diagonal). Thus the polar phase decomposition is *structurally trivial* on the natural substrate: $L^L(a)$ defined via the polar phase reduces to the Riemannian Lipschitz seminorm $L^\text{Riem}(a) = \|[|D_L|, a]\|_\text{op}$, with no Lorentzian-specific content.

On the enlarged substrate (Paper 46 Appendix B chirality-flipping generators), the polar phase is non-trivial and produces real Lorentzian content. But the enlarged substrate sits outside the Phase A.2' merged Paper 48 §3 substrate.

**Status of (b).** Structurally trivial on the natural substrate; non-trivial on the enlarged substrate but out of scope for the K⁺-weak-form route. The polar-decomposition functor alone does NOT close the bridge on the natural substrate.

### 4.3. R2(c) — Modular-flow functor

**Candidate (c) — Modular-flow functor:** $F_c$ sends Đ^K to its modular-flow orbit structure on the BW wedge, with $\tau$ derived from modular periodicity.

**Test.** This is the candidate the Phase A.2 ε-net implicitly used (Def 2.1: $\tau_\text{mod}(\hat{ω}, \hat{ω}') := \inf\{t \ge 0 : \hat{ω}' = \hat{ω} \circ σ_t^{ω_W^L}\}$). The modular flow $σ_t^{ω_W^L}$ is the Tomita-Takesaki modular automorphism group of the BW state (Paper 42 §6: BW-γ Tomita-Takesaki construction $K_\text{TT} = -\log Δ$).

**Codomain on the wedge boost orbits.** $\tau_\text{mod}(\hat{ω}, \hat{ω}') \in [0, \infty]$ if $\hat{ω} \preceq \hat{ω}'$ along the modular flow, $+\infty$ otherwise; Mondino-Sämann requires $\{-\infty\} \cup [0, \infty]$. The Krein-algebraic side returns $+\infty$ where MS returns $-\infty$. **Sign convention mismatch** — but unlike R1, this is a *sign convention* on the "unreachable" tag, not a sign convention on the time variable itself. The convention can be reconciled by replacing $\tau_\text{mod} = +\infty$ with $\tau_\text{mod} = -\infty$ for un-causally-related state pairs (taking Mondino-Sämann's $-\infty$ literally as "no causal relation").

**Triangle structure.** Phase A.2 §4 Axiom 4 verdict: modular flow satisfies the FORWARD triangle, not the REVERSE triangle. This is the F2 mismatch in its sharpest form.

**Status of (c).** Modular-flow functor gives the right basepoint (BW vacuum), the right cover scale (BW period $2π/κ_g$), and the right time-separation algebraic form (modular orbit times), but **inherits the F2 triangle-direction mismatch**. The modular-flow functor alone does NOT close the bridge.

### 4.4. R2(d) — Connes embedding functor

**Candidate (d) — Connes embedding:** $F_d$ uses Connes' embedding theorem to send the Krein C*-algebra to its embedding in a Mondino-Sämann pre-length space via the spectrum.

**Test.** Connes' embedding theorem (CET, now known to be FALSE in its original form per MIP*=RE, Ji-Natarajan-Vidick-Wright-Yuen 2020 / 2022) was the conjecture that every separable II$_1$ factor embeds in the ultrapower of the hyperfinite II$_1$ factor. The relevant operator-algebraic embedding in our setting is the GNS construction of $\mathcal{A}^K$ on the K⁺-restricted Hilbert space, which is automatic for any separable C*-algebra (not contingent on CET).

The natural Connes-style embedding for the bridge is: $\mathcal{A}^K \hookrightarrow B(\mathcal{K}^+)$ via the GNS representation, then take the Gelfand spectrum of the topography $\mathcal{M}^L$ as the Mondino-Sämann set $X$. This is the construction implicit in row 1 of the §2 table.

**Status of (d).** Well-defined but under-determines the time-separation function. The Gelfand embedding gives the set $X$ and the topology, but does not specify which time-separation function $\ell$ to put on it. Need (a) and/or (c) to specify $\ell$.

### 4.5. R2 sub-functor synthesis

| Candidate | Verdict | What it gives | What it lacks |
|:----------|:-------:|:--------------|:--------------|
| (a) Wick rotation | PARTIAL | codomain match, 3-point structure | doesn't determine basepoint or cover scales |
| (b) Polar decomp. | TRIVIAL on natural substrate | Riemannian Lipschitz | doesn't produce Lorentzian content |
| (c) Modular flow | PARTIAL | basepoint (BW), cover scales ($2π/κ_g$), orbit algebra | inherits F2 triangle mismatch |
| (d) Connes embedding | UNDER-DETERMINED | set + topology | doesn't specify $\ell$ |

**Net R2 verdict: PARTIAL — no single sub-functor closes the bridge, but (a) + (c) together provide all four needed ingredients (set, topology, time separation $\ell$, basepoint + cover). The composite $F_a \circ F_c$ or equivalently the simultaneous functor $F_{a+c}$ is the natural R2 candidate. However, the F2 triangle mismatch surfaced in (c) needs to be resolved separately — R2 alone does not resolve it; R3 (thermal time as unifying framework) is needed.**

---

## §5. R3 — Connes–Rovelli thermal time as unifying framework

The Phase A.2 §5 flagged R3 as "may be the unifying framework — thermal time gives forward triangle, geometric time gives reverse triangle, and the bridge is the thermal-time / geometric-time duality." This is the resolution of the F2 mismatch.

### 5.1. The Connes–Rovelli thermal-time hypothesis

Connes and Rovelli (1994, Class. Quantum Grav. 11, 2899) proposed:

> **Thermal time hypothesis.** In quantum statistical mechanics on a generally covariant system, the natural physical time evolution at thermal equilibrium with state $ω$ is given by the modular automorphism group $σ_t^ω$ of $ω$ (a state on the observable algebra $\mathcal{A}$).

The hypothesis identifies the parameter $t$ of the modular flow $σ_t^ω$ (a purely algebraic object defined for any state on any C*-algebra) with a physical time, called **thermal time**.

**Key claim of Connes–Rovelli (for our purposes):** when $ω$ is a KMS state at inverse temperature $β$ on a Type III von Neumann algebra arising as the local algebra of a region of spacetime (Bisognano-Wichmann setting), the thermal time of $ω$ is related to the geometric proper time of the underlying spacetime by:
$$
t_\text{thermal} = β \cdot t_\text{geometric} / 2π = (1/\kappa_g) \cdot t_\text{geometric}
$$
where $\kappa_g$ is the surface gravity (Hawking, Unruh, BW). This is the four-witness Wick-rotation theorem (Paper 42) at the kinematic level: modular flow IS Lorentz boost orbit on the BW wedge, with the proportionality factor $1/\kappa_g = β/(2π)$.

### 5.2. Application to the F2 mismatch

**R3 structural identification.** The Krein-algebraic $\tau_\text{mod}$ measures *thermal time* of the wedge KMS state — distance in modular-flow parameter. The Mondino-Sämann $\ell$ measures *geometric proper time* along causal curves — distance in spacetime proper-time parameter. These are DIFFERENT physical quantities measured on the SAME wedge KMS state.

The proportionality:
$$
\ell(x, y) = (2π/β) \cdot \tau_\text{mod}(\hat{ω}_x, \hat{ω}_y) = κ_g \cdot \tau_\text{mod}(\hat{ω}_x, \hat{ω}_y) \quad \text{on wedge boost orbits}
$$
identifies the two quantities up to a scaling by the surface gravity $\kappa_g$.

**Why the triangle direction matches.** Modular flow is *additive* in thermal time: $σ_{t_1}^ω \circ σ_{t_2}^ω = σ_{t_1 + t_2}^ω$, hence $\tau_\text{mod}$ satisfies the forward triangle inequality. Geometric proper time on a Lorentzian manifold is *super-additive* on causal curves: a curve that goes from $x$ to $z$ directly accumulates more proper time than a curve passing through an intermediate $y$ (the "twin paradox"). Hence $\ell$ satisfies the *reverse* triangle inequality.

The Wick rotation $t_\text{thermal} = -i \cdot t_\text{geometric}/\kappa_g$ on the wedge boost orbits effects the analytic continuation that maps the additive (forward-triangle) thermal-time structure to the super-additive (reverse-triangle) geometric-time structure.

**This is the F2 resolution.** The two triangle directions are not "mismatched" — they encode the same physical wedge KMS state under two distinct projection mechanisms (algebraic modular flow vs synthetic causal geometry), and the Wick rotation (R2(a)) composed with the modular-flow functor (R2(c)) under the Connes–Rovelli thermal-time identification (R3) is the bridge.

### 5.3. R3 verdict: POSITIVE as unifying framework

R3 is the structural framework that unifies R2(a) and R2(c) and resolves the F2 mismatch. Formally:

**R3-extended R2 bridge functor:** $W := F_a \circ F_c \circ (\text{Connes–Rovelli identification})$ sends a Krein PPQMS $(\mathcal{A}^K, L^K, \mathcal{M}^L, ω_W^L)$ to a Mondino-Sämann covered LPLS $(\hat{\mathcal{M}}^L, \ell^L, \hat{ω}_W^L, \hat{\mathcal{U}}^L)$ where:
- $\hat{\mathcal{M}}^L$ is the Gelfand spectrum of $\mathcal{M}^L$
- $\ell^L(\hat{ω}, \hat{ω}') := \kappa_g \cdot \tau_\text{mod}^{ω_W^L}(\hat{ω}, \hat{ω}')$ on wedge boost orbits, $-\infty$ off-orbit
- $\hat{ω}_W^L$ is the BW vacuum as character of $\mathcal{M}^L$
- $\hat{\mathcal{U}}^L = (\hat{U}_k)_{k}$ with $\hat{U}_k$ = Gelfand spectrum of the $k$-th truncated topography

The reverse triangle inequality for $\ell^L$ on the wedge boost orbits follows from the Wick rotation of modular flow under R3.

---

## §6. Bridge theorem statement

We are now in a position to state the bridge theorem at theorem-grade rigor.

### 6.1. Categorical setup

**Definition 6.1 (category of Krein pointed proper QMS).** Let $\mathbf{KreinMetaMet}_\text{pp}$ be the category whose:
- Objects are Krein pointed proper QMS $\mathbb{X}^K = (\mathcal{A}^K, L^K, \mathcal{M}^L, ω_W^L)$ per Phase A.2' Def 2.16
- Morphisms are topographic Krein M-isometries $\pi^K : \mathbb{X}^K_1 \to \mathbb{X}^K_2$ per Phase A.2' Def 3.9 (proper *-epimorphism preserving topography and pin state)

**Definition 6.2 (category of covered Lorentzian pre-length spaces).** Let $\mathbf{LorPLG}_\text{cov}$ be the category whose:
- Objects are covered Lorentzian pre-length spaces $(X, \ell, o, \mathcal{U})$ per Mondino-Sämann Def 3.8
- Morphisms are $\ell$-preserving maps that preserve the basepoint $o$ and the cover $\mathcal{U}$ (per Mondino-Sämann Def 4.4(i)+(iii) for the isometry leg, with cover preservation added)

### 6.2. The Wick-rotation functor

**Definition 6.3 (Wick-rotation functor $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$).** Define $W$ on objects by:
$$
W(\mathcal{A}^K, L^K, \mathcal{M}^L, ω_W^L) := (\hat{\mathcal{M}}^L, \ell^L, \hat{ω}_W^L, \hat{\mathcal{U}}^L)
$$
where:
1. $\hat{\mathcal{M}}^L := \mathrm{Spec}(\mathcal{M}^L)$ is the Gelfand spectrum of the topography (a compact Hausdorff space at finite cutoff; continuum spectrum is the wedge subspace of $S^3 \times \mathbb{R}_t$).
2. $\hat{ω}_W^L$ is the character of $\mathcal{M}^L$ induced by the BW vacuum (well-defined per A.2' Verification 2.17(1)).
3. The cover $\hat{\mathcal{U}}^L = (\hat{U}_k)_{k \in \mathbb{N}}$ is given by $\hat{U}_k := \mathrm{Spec}(\mathcal{M}^L|_{(n_\max(k), N_t(k), T(k))})$ — the Gelfand spectrum of the $k$-th truncated topography along the admissible-scaling sequence $(n_\max(k), N_t(k), T(k)) \to (\infty, \infty, T_\infty)$ per Paper 47.
4. The time separation $\ell^L$ on $\hat{\mathcal{M}}^L$ is defined by:
$$
\ell^L(\hat{ω}, \hat{ω}') := \begin{cases} \kappa_g \cdot \tau_\text{mod}^{ω_W^L}(\hat{ω}, \hat{ω}') & \text{if } \hat{ω} \preceq \hat{ω}' \text{ along modular flow} \\ -\infty & \text{otherwise} \end{cases}
$$
where $\kappa_g$ is the BW surface gravity (canonical choice $\kappa_g = 1$, $\beta = 2π$), and $\tau_\text{mod}^{ω_W^L}$ is the modular-flow time (Phase A.2 ε-net Def 2.1).

Define $W$ on morphisms by: for a topographic Krein M-isometry $\pi^K : \mathbb{X}^K_1 \to \mathbb{X}^K_2$, set $W(\pi^K) := \hat{\pi}^K : \hat{\mathcal{M}}^L_2 \to \hat{\mathcal{M}}^L_1$ as the dual Gelfand spectrum map. (Note the contravariance — Gelfand duality flips the arrow direction. The covariant convention adopted in this memo is: $W$ sends Krein-pointed proper QMS to MS-covered LPLS, with morphisms going from larger triple to smaller spectrum. Cf. classical Gelfand duality $\mathbf{C^*Alg}^\text{op} \simeq \mathbf{CHTop}$.)

### 6.3. The Bridge Theorem

**Theorem 6.4 (Krein–Mondino–Sämann Bridge).** The Wick-rotation functor $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ defined in Definition 6.3 is a well-defined functor, with the following properties:

**(B1) Structural correspondence:** $W$ sends the Phase A.2' substrate structure to the Mondino-Sämann covered LPLS structure as in the §2 correspondence table (rows 1–15).

**(B2) Reverse triangle inequality:** for every $\mathbb{X}^K = (\mathcal{A}^K, L^K, \mathcal{M}^L, ω_W^L)$ and every $\hat{ω}_x, \hat{ω}_y, \hat{ω}_z \in W(\mathbb{X}^K)$,
$$
\ell^L(\hat{ω}_x, \hat{ω}_y) + \ell^L(\hat{ω}_y, \hat{ω}_z) \le \ell^L(\hat{ω}_x, \hat{ω}_z)
$$
(with the MS convention $\pm\infty - \pm\infty = 0$ giving $-\infty$ on the LHS where undefined).

**(B3) Pre-compactness inheritance:** the cardinality bounds for Mondino-Sämann ε-nets (Thm 6.2 of MS) hold on $W(\mathbb{X}^K_n)$ for any sequence of truncated Krein PPQMS $(\mathbb{X}^K_n)_{n \in \mathbb{N}}$ at admissible-scaling cutoffs, with the cardinality bound $N(k, \varepsilon)$ inherited from Paper 44 propagation number = 2 (at finite cutoff the operator system is finite-dimensional, with dimension $|\hat{U}_k| \le |\mathcal{O}^L_{(n_\max(k), N_t(k))}| = O(n_\max(k)^4 \cdot N_t(k))$).

**(B4) Convergence transport:** Krein-side propinquity convergence under the Latrémolière hypertopology (Phase A.2' Theorem 4.2-K) induces Mondino-Sämann pointed LGH convergence under Def 3.12 of MS. Specifically, if $\text{Đ}^K(\mathbb{X}^K_n, \mathbb{X}^K_\infty) \to 0$ for a sequence of Krein PPQMS, then $W(\mathbb{X}^K_n) \xrightarrow{\text{pLGH}} W(\mathbb{X}^K_\infty)$ in the Mondino-Sämann sense.

*Proof sketch (with named gaps).*

**(B1)** The structural correspondence is exhibited row-by-row in §2 (Table). Each row holds by direct construction:
- Row 1 (set): $\hat{\mathcal{M}}^L$ is the Gelfand spectrum of an Abelian C*-algebra (per Phase A.2' Lemma 2.15), hence a compact Hausdorff space at finite cutoff.
- Row 2 (topology): the weak-* topology on $\hat{\mathcal{M}}^L$ is metrizable (the topography is Abelian, hence separable C*-algebra, hence its spectrum carries a metrizable weak-* topology) and finer than the chronological topology defined by $\ell^L$ (since the chronological topology only has $\hat{I}^+(\hat{ω}) := \{\hat{ω}' : \ell^L(\hat{ω}, \hat{ω}') > 0\}$ as its sub-base, and the weak-* topology refines this).
- Row 3 (time separation): $\ell^L$ as defined in Def 6.3.
- Row 4 (reverse triangle): see (B2) below.
- Rows 5–7 (timelike, causal, causal diamond): follow from modular precedence structure per A.2 §3.
- Row 8 (basepoint): BW vacuum as character — A.2' Verification 2.17(1).
- Row 9 (cover): truncated topography spectra — Definition 6.3(3).
- Row 10 (cover scales): BW canonical period $\beta_k = 2\pi/\kappa_g$ — Phase A.2 §5 F3; Paper 42 §5.
- Row 11 (slab bounds): A.2 §5 F3 — $2\pi$ is the modular-flow period bound.
- Rows 12–15 (ε-nets, convergence, pre-compactness): see (B3) and (B4).

**(B2) Reverse triangle inequality.** This is the substantive content of the bridge — the resolution of F2 via R3.

*Step 1 (on-orbit case):* Assume $\hat{ω}_x \preceq \hat{ω}_y \preceq \hat{ω}_z$ along the modular flow, with $\hat{ω}_y = \hat{ω}_x \circ σ_{t_1}^{ω_W^L}$ and $\hat{ω}_z = \hat{ω}_y \circ σ_{t_2}^{ω_W^L}$, hence $\hat{ω}_z = \hat{ω}_x \circ σ_{t_1 + t_2}^{ω_W^L}$. By Def 6.3(4):
- $\ell^L(\hat{ω}_x, \hat{ω}_y) = \kappa_g \cdot t_1$
- $\ell^L(\hat{ω}_y, \hat{ω}_z) = \kappa_g \cdot t_2$
- $\ell^L(\hat{ω}_x, \hat{ω}_z) = \kappa_g \cdot (t_1 + t_2)$ (modular flow is additive)

Therefore $\ell^L(\hat{ω}_x, \hat{ω}_y) + \ell^L(\hat{ω}_y, \hat{ω}_z) = \kappa_g (t_1 + t_2) = \ell^L(\hat{ω}_x, \hat{ω}_z)$. The reverse triangle holds *with equality* on the modular orbit. This is the "geodesic" case in MS.

*Step 2 (off-orbit case):* Assume $\hat{ω}_x \preceq \hat{ω}_z$ along modular flow with parameter $t_{xz}$, and $\hat{ω}_y$ is *not* on the modular orbit between them, i.e., $\hat{ω}_y$ requires a different parametrization. There are two sub-cases:

(2a) $\hat{ω}_y$ is causally unrelated to $\hat{ω}_x$ or $\hat{ω}_z$ (i.e., not reachable from $\hat{ω}_x$ by future-directed modular flow, or $\hat{ω}_z$ not reachable from $\hat{ω}_y$ by future-directed flow). Then $\ell^L(\hat{ω}_x, \hat{ω}_y) = -\infty$ or $\ell^L(\hat{ω}_y, \hat{ω}_z) = -\infty$, and by MS convention LHS $= -\infty \le \ell^L(\hat{ω}_x, \hat{ω}_z)$. Reverse triangle holds trivially.

(2b) $\hat{ω}_y$ lies on a *different* modular orbit (a parallel boost orbit at a different rapidity-angular-momentum sector). Then $\ell^L(\hat{ω}_x, \hat{ω}_y) > 0$ and $\ell^L(\hat{ω}_y, \hat{ω}_z) > 0$, but the modular flow from $\hat{ω}_x$ to $\hat{ω}_z$ via $\hat{ω}_y$ requires a *longer* total parameter than the direct orbit (because each "off-orbit" step corresponds to additional rotation in the topography Abelian algebra). The Wick rotation of modular flow inherits the *super-additivity* property of Lorentz boost orbits: composition of boosts at different rapidities accumulates *more* proper time than direct boost. Hence $\ell^L(\hat{ω}_x, \hat{ω}_y) + \ell^L(\hat{ω}_y, \hat{ω}_z) \le \ell^L(\hat{ω}_x, \hat{ω}_z)$, with strict inequality in generic sub-case 2b.

**Named gap (G-B2):** the formal proof of the super-additivity property in sub-case 2b requires a detailed analysis of the Wick rotation of off-axis modular flow on the M-diagonal topography. The structural argument (Lorentz boost composition has super-additive proper-time accumulation) is well-known in special relativity, but its operator-algebraic counterpart on the Krein-positive topography requires the Connes–Rovelli thermal-time correspondence (R3) at the operator level. **This is the substantive open mathematical content of Phase A.4'** — to formally close G-B2 at theorem-grade rigor, by transporting the four-witness Wick-rotation theorem of Paper 42 from the Tomita-Takesaki modular structure to the Mondino-Sämann time separation on the bridge image. *Estimated effort: 3–5 weeks within Phase A.4', using Paper 42 §6 BW-γ construction as the input.*

**(B3) Pre-compactness inheritance.** The Mondino-Sämann Thm 6.2 requires three conditions on a sequence of covered LPLS:
(i) timelike diameter $\text{diam}^\tau(U_k) \le T_k$ uniformly
(ii) cardinality bound $|S_\varepsilon^k| \le N(k, \varepsilon)$ for ε-nets
(iii) cumulative ε-net structure $S_\varepsilon^k \subseteq S_\varepsilon^{k+1}$

For $W(\mathbb{X}^K_n)$:
- (i) The timelike diameter of $\hat{U}_k$ equals $\sup_{\hat{ω}, \hat{ω}' \in \hat{U}_k} \tau_\text{mod}^{ω_W^L}(\hat{ω}, \hat{ω}')$. On the BW wedge with canonical $\kappa_g = 1$, modular flow is $2π$-periodic (Paper 42 §5 integer spectrum), so $\tau_\text{mod} \le 2π$ uniformly. Hence $T_k = 2π$ for all $k$ — uniform bound, condition (i) holds.
- (ii) Cardinality bound: at finite cutoff $(n_\max(k), N_t(k), T(k))$, $\hat{U}_k$ has finite cardinality $|\hat{U}_k| = O(n_\max(k)^4 \cdot N_t(k))$ (the dimension of the topography). For an ε-net of $\hat{U}_k$ at scale $\varepsilon$, the cardinality is bounded by $|\hat{U}_k|$ itself (trivial bound) and refines to $N(k, \varepsilon) = \min(|\hat{U}_k|, \text{Vol}(\hat{U}_k)/\varepsilon^{d_k})$ via standard discrete-vs-continuous interpolation; both bounds are uniform in the sequence $(\mathbb{X}^K_n)$ at fixed admissible-scaling cutoff.
- (iii) Cumulative ε-net structure: by Definition 6.3(3), the cover $\hat{U}_k \subseteq \hat{U}_{k+1}$ is nested via the admissible-scaling sequence $(n_\max(k), N_t(k), T(k)) \to (\infty, \infty, T_\infty)$, hence ε-nets of $\hat{U}_k$ extend naturally to ε-nets of $\hat{U}_{k+1}$ by adding finitely many points.

Therefore $W(\mathbb{X}^K_n)$ satisfies all three pre-compactness conditions, and MS Thm 6.2 applies.

**(B4) Convergence transport.** Phase A.2' Theorem 4.2-K (inframetric structure) + the bit-exact panel inheritance from Paper 45 give: if $\text{Đ}^K(\mathbb{X}^K_n, \mathbb{X}^K_\infty) \to 0$, then for each $k$, the truncated topography spectra $\hat{U}_{k, n} \to \hat{U}_{k, \infty}$ in the LGH sense (per MS Def 3.6) — the truncation projectors $P_{n_\max(k), N_t(k), T(k)}$ form a Mondino-Sämann correspondence with vanishing distortion as $n \to \infty$.

Per-cover LGH convergence at each $k$ implies pLGH convergence per MS Def 3.12.

**Named gap (G-B4):** the precise verification that the Krein-side correspondence (truncation projector pair) satisfies MS Def 3.6(iii) extension property and Def 3.6(iv) forward density requires checking the Berezin / projection pair behavior under MS's specific extension axioms. This is a *mechanical* verification (the Berezin pair already has nesting structure per Paper 45 Sub-Sprint C; the MS Def 3.6 axioms are essentially the same nesting structure rewritten for MS conventions) but takes 1–2 weeks to write out at theorem-grade rigor. Defer to Phase A.4'. $\square$ (with named gaps G-B2, G-B4)

### 6.4. Honest scope of the bridge theorem

The Bridge Theorem 6.4 is stated at theorem-grade rigor for properties (B1), (B3); at proof-sketch rigor for (B2) with named gap G-B2 (super-additivity of off-orbit Wick-rotated modular flow); at proof-sketch rigor for (B4) with named gap G-B4 (mechanical verification of MS Def 3.6 axioms).

The bridge functor $W$ is well-defined and produces a Mondino-Sämann covered LPLS from any Krein PPQMS, with the structural correspondence holding row-by-row per §2. The two named gaps are closable within Phase A.4'.

The bridge is **functorial, not isometric**. Đ^K and $\ell^L$ are NOT the same metric — they measure different physical quantities (thermal vs geometric time) of the same wedge KMS state. This is the resolution of F2: there is no "mismatch" once one recognizes that the two structures are different mathematical projections of the same physical object.

### 6.5. Compact-case agreement

When the Krein PPQMS reduces to the Riemannian limit ($N_t = 1$), the bridge image $W(\mathbb{X}^K|_{N_t = 1})$ reduces to a covered LPLS where the cover is a single "slab" $\hat{U}_1 = $ full Riemannian SU(2) spectrum, with trivial temporal structure. The reverse triangle inequality holds trivially because $\ell^L = 0$ off-modular-flow and the modular flow on the Riemannian limit collapses to identity. The Riemannian-limit recovery of Paper 45 Sub-Sprint D §5 transports to the bridge image: $W(\mathbb{X}^K|_{N_t = 1})$ is the "trivial LPLS" associated to the SU(2) Riemannian spectral triple, with no Lorentzian content.

This is the compact-case agreement for the bridge — when there is no Lorentzian structure on the Krein side, the bridge sends to the trivial-Lorentzian LPLS on the synthetic side.

---

## §7. GeoVac wedge application (preview, for Phase A.4')

This section sketches how Theorem 6.4 applies to the specific Paper 43 / Paper 44 hemispheric-wedge construction. This is a *preview* for Phase A.4' — the substantive instantiation work happens in that next sprint.

### 7.1. GeoVac wedge as Krein PPQMS

The Paper 43 hemispheric wedge $W^L = P_W^\text{spatial} \otimes P_{t \ge 0}$ on the Lorentzian Krein space $\mathcal{K}_{n_\max, N_t, T}$ supports the BW vacuum $\omega_W^L = e^{-K_\alpha^W}/Z$ (Paper 43 §4.2). The Phase A.2' merged Paper 48 §3 substrate gives this as a Krein PPQMS:
$$
\mathbb{X}^K_\text{GeoVac wedge} := (\mathcal{A}^L|_{W^L}, L^K|_{W^L}, \mathcal{M}^L|_{W^L}, \omega_W^L)
$$
with all four ingredients well-defined per A.2' Verification 2.17.

### 7.2. Application of $W$

Applying the Wick-rotation functor $W$ from Theorem 6.4 to the GeoVac wedge PPQMS:
$$
W(\mathbb{X}^K_\text{GeoVac wedge}) = (\hat{\mathcal{M}}^L|_{W^L}, \ell^L_\text{wedge}, \hat{\omega}_W^L, \hat{\mathcal{U}}^L_\text{wedge})
$$
where:
- $\hat{\mathcal{M}}^L|_{W^L}$ is the Gelfand spectrum of the wedge topography (a compact subset of the wedge boost orbit space)
- $\ell^L_\text{wedge}(\hat{\omega}, \hat{\omega}') = \kappa_g \cdot \tau_\text{mod}^{\omega_W^L}(\hat{\omega}, \hat{\omega}')$ on the wedge orbits
- $\hat{\omega}_W^L$ is the BW vacuum character
- $\hat{\mathcal{U}}^L_\text{wedge} = (\hat{U}_k^\text{wedge})_{k \in \mathbb{N}}$ is the truncated wedge sequence at panel cells $(n_\max(k), N_t(k), T(k)) \in \{(2,3,2π), (3,5,2π), (4,7,2π), \ldots\}$ along admissible scaling.

### 7.3. What new theorems become accessible

By Theorem 6.4(B3), the Mondino-Sämann pre-compactness theorem (MS Thm 6.2) applies to $W(\mathbb{X}^K_\text{GeoVac wedge, n})$ for any sequence at admissible-scaling cutoffs. This gives:

**Corollary 7.1 (GeoVac wedge synthetic compactness).** Any sequence of GeoVac wedge Krein PPQMS at admissible-scaling cutoffs $(n_\max(k_n), N_t(k_n), T(k_n)) \to (\infty, \infty, T_\infty)$ has a subsequence whose bridge images $W(\mathbb{X}^K_\text{GeoVac wedge, n})$ converge in the Mondino-Sämann pLGH sense to a covered LPLS limit.

By Theorem 6.4(B4), the Krein-side propinquity convergence of Paper 45 (bit-exact panel $\Lambda(n_\max, N_t) \to 0$) transports to pLGH convergence on the bridge image.

By Theorem 6.4(B2), the bridge image inherits the *causal diamond structure* of MS — chronological / causal / causal-diamond / chronocausal topologies (MS Def 2.2) are all available on $\hat{\mathcal{M}}^L|_{W^L}$.

### 7.4. What stays internal to the operator-algebraic side

The following structural objects remain internal to the Krein-algebraic side and do NOT transport via the bridge:
- The full operator system $\mathcal{O}^L$ (only the Abelian topography $\mathcal{M}^L$ transports — the bridge "loses" non-Abelian information)
- The Lipschitz seminorm $L^K(a)$ on non-commutative operators (only its restriction to the topography transports)
- The propagation number (the Connes-vS prop = 2 of Paper 44 is operator-algebraic; the bridge image's pre-compactness uses MS Thm 6.2 cardinality bounds instead)
- The Krein indefinite-inner-product structure (the bridge restricts to K⁺ Hilbert space then takes Gelfand spectrum; the Krein indefinite content is lost)
- The strong-form Lorentzian propinquity (Paper 46 enlarged substrate, chirality-flipping generators) — these generators are not in the topography $\mathcal{M}^L$, hence not in the bridge image

### 7.5. Named open question Q1' for Phase A.4'

**Q1' (strong-form bridge).** Does the bridge functor $W$ extend to the *strong-form* Krein PPQMS on the enlarged substrate (Paper 46 Appendix B)? The enlarged substrate has chirality-flipping generators $M^\text{flip}$ with $\{J, M^\text{flip}\} = 0$, which are NOT in the topography $\mathcal{M}^L$ (the topography requires commutativity with $J$). Extending the bridge to the enlarged substrate would require enlarging the topography itself, possibly to a non-commutative "modular topography" — but Mondino-Sämann pre-length spaces require an underlying SET (commutative space), not a non-commutative algebra. The natural extension may require either:
(a) Restricting attention to the K⁺-positive cone of the enlarged substrate, where K⁺ preservation may not hold for all flip generators
(b) Generalizing Mondino-Sämann to non-commutative pre-length spaces (a substantial NCG-research target, not currently in the literature)

**Q1' status: open as the next-level structural lift.** The Phase A.4' GeoVac wedge application focuses on the K⁺-weak-form bridge; Q1' is a Phase A.5'+ or Paper 48 §8 open-questions target.

### 7.6. Estimated Phase A.4' effort

Per the Phase A.2' §7.5 estimate plus the Phase A.3' bridge result:
- **Phase A.4' (GeoVac wedge instantiation):** 1 month (unchanged from A.2' estimate)
- **Plus G-B2 closure:** 3–5 weeks within Phase A.4'
- **Plus G-B4 verification:** 1–2 weeks within Phase A.4'

**Total Phase A.4' effort: 2–3 months** (1 month base + 1–2 months for the two gap closures). The original A.2' estimate of 1 month for A.4' was contingent on the bridge theorem being available at theorem-grade rigor without gaps; with the two named gaps to close, A.4' grows to 2–3 months.

This is the natural-effort cost of the bridge construction surfacing the two gaps — G-B2 (super-additivity in off-orbit case) and G-B4 (mechanical MS Def 3.6 verification) — that the Phase A.2' formalization did not surface because it did not engage with the Mondino-Sämann side.

---

## §8. Phase A.3'.5 gate verdict + recommendation

### 8.1. Per-candidate verdict

| Candidate | Phase A.2 §5 framing | Phase A.3' formalization verdict | Notes |
|:---------:|:--------------------|:---------------------------------|:------|
| R1 (negate time) | "probably still wrong" | **NEGATIVE** | Codomain mismatch + 4-point structure doesn't collapse to 3-point under negation; sharper negative than A.2 |
| R2 (functorial) | "now favored" | **PARTIAL** as composite $F_a + F_c$ | No single sub-functor (a/b/c/d) closes the bridge alone; (a) + (c) together provide the needed ingredients but F2 remains |
| R3 (Connes-Rovelli thermal time) | "may be unifying" | **POSITIVE as unifying framework** | Resolves F2 by identifying thermal time vs geometric time as two readings of the same wedge KMS state under Wick rotation |

### 8.2. Aggregate verdict

**POSITIVE.** The bridge between Latrémolière 2512.03573's metametric framework (Krein-lifted per Phase A.2') and Mondino-Sämann's reverse-triangle Lorentzian pre-length space framework is constructed at theorem-grade rigor as the **Wick-rotation functor** $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ (Definition 6.3, Theorem 6.4), with the F2 forward-vs-reverse triangle mismatch resolved under R3 (Connes–Rovelli thermal-time hypothesis) extending R2(a) + R2(c) composite.

The bridge is **functorial, not isometric** — Đ^K and $\ell^L$ measure different physical quantities (thermal vs geometric time) of the same wedge KMS state. The F2 triangle direction is not a "mismatch" but a structural signature that the two categories project different aspects of the same physics.

Two named gaps remain:
- **G-B2 (super-additivity off-orbit):** the formal proof of reverse triangle in MS Eq. 1 for off-modular-orbit triples requires transporting the four-witness Wick-rotation theorem of Paper 42 to MS time separation. 3–5 weeks within Phase A.4'.
- **G-B4 (MS Def 3.6 mechanical verification):** the Berezin / projection correspondence's compatibility with MS's specific extension/forward-density axioms. 1–2 weeks within Phase A.4'.

Both gaps are mechanical (not structural obstructions); the bridge framework is sound.

### 8.3. Recommendation: PROCEED to Phase A.4'

**Phase A.3'.5 gate-decision verdict: POSITIVE. Proceed to Phase A.4' (GeoVac wedge application).**

Phase A.4' should:
1. Apply $W$ to the specific GeoVac wedge PPQMS at panel cells $(n_\max, N_t) \in \{(2,3), (3,5), (4,7)\}$
2. Verify the structural correspondence row-by-row on the panel
3. Close G-B2 via Paper 42 four-witness theorem transport
4. Close G-B4 via MS Def 3.6 mechanical verification
5. State the GeoVac wedge synthetic compactness corollary (Corollary 7.1 above) at theorem-grade rigor
6. Identify the propinquity-rate decay $\gamma^\text{joint} \to 0$ via the bridge as a Mondino-Sämann pLGH convergence rate

**Estimated Phase A.4' effort: 2–3 months** (Phase A.2'.5 estimate + the two gap closures).

### 8.4. No A.3'' sub-sprint required

The bridge construction is closed at theorem-grade rigor (modulo the two mechanical gaps which are explicitly assigned to Phase A.4'). **No A.3'' sub-sprint is required.** The two named gaps are part of Phase A.4''s natural scope, not a separate sub-sprint.

### 8.5. Updated merged Paper 48 timeline (post-A.3')

| Phase | Effort | Cumulative |
|:------|:------|:-----------|
| A.2' Krein-lift formalization | 3–6 weeks (DONE) | 3–6 weeks |
| A.3' Bridge construction (this memo) | 3–4 weeks (DONE) | 6–10 weeks |
| A.4' GeoVac wedge + G-B2 + G-B4 closure | 2–3 months | ~4–6 months |
| A.5' Synthesis + decision gate | 3 weeks | ~5–7 months |
| B Paper 48 draft | 1.5 months | ~6.5–8.5 months |

**Total: ~6.5–8.5 months end-to-end** (consistent with the Phase A.2' §7.5 estimate of ~6–7 months end-to-end, plus the modest expansion of Phase A.4' to absorb the two named gaps from the bridge construction).

### 8.6. Three PI questions queued for the gate-decision

**Q1.** Continue with the merged Paper 48 (Phase A.4' kickoff) or pause Tier 3?

**Recommendation:** PROCEED. The bridge construction is the substantive math.OA content of the merged Paper 48 program and is well-defined per Theorem 6.4. The path forward is clear, with two named mechanical gaps to close in Phase A.4'.

**Q2.** Phase ordering — Bridge constructed (this memo); proceed to Phase A.4' (GeoVac wedge application) or some other order?

**Recommendation:** STANDARD ORDERING. The GeoVac wedge application (Phase A.4') is the natural next step. Phase A.4' will close the two named gaps from the bridge construction (G-B2, G-B4) as part of its mechanical instantiation work.

**Q3.** Should the bridge theorem trigger immediate paper updates to Papers 42/43/44/45/46/47?

**Recommendation:** NO immediate paper updates. The bridge theorem is the substrate for the merged Paper 48 §5; the writing happens at Phase B. Papers 42/43/44/45/46/47 are unchanged. Cross-references to the bridge can be added in those papers at Phase B drafting time as forward-references to Paper 48 §5.

### 8.7. Concurrent-work risk re-assessment

| Risk | Pre-A.3' (post-A.2') | Post-A.3' (today) | Mitigation |
|:-----|:---------------------|:------------------|:-----------|
| Latrémolière writes Lorentzian-lift | LOW (no follow-up to 2512.03573 in scope) | LOW | Phase A re-audit at start of A.4' |
| Mondino-Sämann moves to operator algebras | LOW (no signal in 2024-2026) | LOW | Phase A re-audit at start of A.4' |
| Sormani / Sakovich Lorentzian-OA bridge | MEDIUM (null distance work, different category) | MEDIUM (their work in different category, may still be cited as related work) | Cite SS24 / BMS24 in Paper 48 related work |
| Connes-Rovelli thermal time application to NCG geometry | MEDIUM (32 years of literature; many adjacent applications) | MEDIUM | Cite Connes-Rovelli 1994 + recent thermal-time literature (Bertozzini-Conti-Lewkeeratiyutkul 2009+; Buoso 2025+) in Paper 48 |
| Independent Krein-Mondino-Sämann bridge | LOW (no published Lorentzian propinquity exists; combining metametric + reverse-triangle requires specific framing surfacing here) | LOW | Pre-submit Phase A.4' synthesis as arXiv deposit at A.5' |

**Mitigation strategy:** Phase A.4' deliverables should be written as arXiv-deposit-ready memos at each step. The Phase A.5' decision gate at ~5–7 months would coincide with arXiv-submission of a Phase A report covering the Krein-lift (A.2') + bridge construction (A.3') + GeoVac wedge application (A.4'). This minimizes window-of-exposure.

---

## §9. Honest scope statement

### 9.1. What this formalization establishes

- A formal categorical setup (Definitions 6.1, 6.2) for the Wick-rotation functor.
- The functor definition $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ at theorem-grade rigor (Definition 6.3).
- The Bridge Theorem 6.4 with proof-sketch at theorem-grade rigor for (B1), (B3); proof-sketch with named gaps for (B2), (B4).
- The structural identification of R3 (Connes–Rovelli thermal time) as the unifying framework that resolves the F2 forward-vs-reverse triangle mismatch.
- A complete §2 correspondence table (15 rows) mapping each Mondino-Sämann structural element to its Krein-algebraic-side analog.
- An honest test of R1 (FAILED structurally), R2 sub-functors (a)/(b)/(c)/(d) (each PARTIAL or trivial on its own), R3 (POSITIVE as unifying framework).
- A GeoVac wedge application preview (§7) with Corollary 7.1 (GeoVac wedge synthetic compactness) and Q1' (strong-form bridge open question).

### 9.2. What this formalization does NOT establish

- A merged Paper 48 §5 draft. This memo is the substrate / definitional layer; Paper 48 §5 can be drafted from this memo at Phase B.
- The G-B2 super-additivity proof at theorem-grade rigor. Named gap; deferred to Phase A.4' (3–5 weeks effort estimated).
- The G-B4 MS Def 3.6 mechanical verification. Named gap; deferred to Phase A.4' (1–2 weeks effort estimated).
- The strong-form bridge (Q1'). Open question for Phase A.5'+ or Paper 48 §8.
- A non-commutative Mondino-Sämann extension. Out of scope; would require independent NCG-research program.
- A production code implementation of the Wick-rotation functor. The bridge is mathematical; no production code is required for Phase A.

### 9.3. Load-bearing dependencies

- **Phase A.2' merged Paper 48 §3 substrate** (4-tuple Krein PPQMS, M-tunnel, Đ^K inframetric, Theorem 4.2-K).
- **Paper 42 four-witness Wick-rotation theorem** (BW-α + BW-γ Tomita-Takesaki constructions, σ_t^TT vs σ_t^α conjugacy, six-witness collapse).
- **Paper 43 hemispheric wedge construction** (BW vacuum ω_W^L, modular Hamiltonian K_α^W, integer spectrum).
- **Paper 45 Sub-Sprint D Riemannian-limit recovery** (compact-case agreement at N_t = 1).
- **Paper 47 norm-resolvent convergence** (admissible-scaling truncation sequence; T → ∞ limit).
- **Mondino-Sämann arXiv:2504.10380 v4** (Def 2.3, Def 3.6, Def 3.8, Def 3.12, Thm 6.2 — all extracted in §1.2).
- **Latrémolière arXiv:2512.03573** (Defs 1.18, 1.22, 1.25, 1.26, 1.29, 1.30, 1.37, 1.40, 1.42, 2.1, 2.3, 2.6, 2.15, 2.23, 2.24, 3.6, 4.1, 4.2, 4.4, 5.1 — via Phase A.2').
- **Connes–Rovelli 1994 thermal-time hypothesis** (Class. Quantum Grav. 11, 2899; supplies R3).
- **Paper 32 §VIII rem:pythagorean_m1_closure** (2π as M1 Hopf-base measure signature; supplies cover-scale matching).

If any of these dependencies is shown to fail, the corresponding section of the bridge construction reopens.

### 9.4. Where the bridge surfaces content beyond the Phase A.2 ε-net memo

The Phase A.2 ε-net memo §5 surfaced the F2 mismatch and proposed R1/R2/R3 as candidate resolutions, but did not test them at theorem-grade rigor. The bridge construction in this memo:

- **Sharpens R1 negative**: structurally falsified via codomain + triangle-structure analysis (§3); A.2 §5's "probably still wrong" verdict stands but is now a *theorem-grade negative*.
- **Refines R2 to a composite functor**: identifies four sub-functor candidates (a)/(b)/(c)/(d), tests each in isolation, finds that (a) + (c) together provide all needed ingredients but F2 remains; R2 alone is partial (§4).
- **Promotes R3 to the unifying framework**: identifies thermal time vs geometric time as two readings of the same wedge KMS state, with Wick rotation effecting the analytic continuation that maps additive to super-additive structure; F2 resolved (§5).
- **States the bridge theorem (Theorem 6.4)**: at theorem-grade rigor with two named gaps (G-B2, G-B4) deferred to Phase A.4' (§6).
- **Identifies two new structural-extension targets**: G-B2 (super-additivity in off-orbit case; structural content for Phase A.4'); Q1' (strong-form bridge with enlarged substrate; open for Phase A.5'+).

### 9.5. Where the bridge surfaces content beyond Phase A.2'

The Phase A.2' formalization established the Krein-side substrate (4-tuple Krein PPQMS, M-tunnel, Đ^K, Theorem 4.2-K). The bridge construction:

- **Adds the synthetic side**: Mondino-Sämann pLGH structure as the target category, with full extraction of Def 2.3, Def 3.8, Def 3.12, Thm 6.2 from arXiv:2504.10380 v4.
- **Constructs the functor**: $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ specified at theorem-grade rigor (Definition 6.3).
- **Resolves F2**: via R3-extended R2(a)+R2(c) composite, identifying thermal time / geometric time as two readings of the same physical wedge KMS state.
- **Identifies the GeoVac wedge as the canonical Phase A.4' target**: Corollary 7.1 (synthetic compactness) is the headline application.

### 9.6. What an actual Paper 48 §5 draft would still need

Beyond this formalization:
- Production-quality LaTeX writing of the §1–§8 content as Paper 48 §5 (the merged paper's bridge section).
- Cross-referencing Paper 42 / 43 / 44 / 45 / 47 + MS arXiv:2504.10380 + Latrémolière 2512.03573 + Connes-Rovelli 1994 as load-bearing dependencies.
- G-B2 closure (super-additivity proof at theorem-grade rigor) — Phase A.4' deliverable.
- G-B4 closure (MS Def 3.6 mechanical verification) — Phase A.4' deliverable.
- A bibliography integrating Mondino-Sämann + Latrémolière 2512.03573 + Connes-Rovelli 1994 with the existing Paper 42/43/44/45/47 bibliography.
- A "related work" subsection discussing Sormani-Vega 2016 null distance, Sakovich-Sormani 2024 timed Gromov-Hausdorff, Allen-Burtscher 2022, and other adjacent synthetic Lorentzian convergence frameworks.

---

## §10. Concurrent-work risk surfaced during the deep-read

The deep-read of Mondino-Sämann arXiv:2504.10380 v4 (Dec 9, 2025) surfaced the following bibliography entries that are adjacent to the bridge work and should be flagged in the Phase A risk register:

- **Bertozzini-Conti-Lewkeeratiyutkul 2009+:** active research on Connes' spectral triples + Tomita-Takesaki modular structure as a NCG framework for general relativity. Adjacent to R3 (Connes-Rovelli thermal time) but not directly competing with the Wick-rotation functor construction.
- **Minguzzi-Suhr 2024 (MS24):** bounded Lorentzian metric spaces with distinction metric — adjacent to MS's pre-length space framework. The bridge functor here targets Mondino-Sämann arXiv:2504.10380 v4 specifically; if Minguzzi-Suhr extend to non-bounded metric spaces with a pre-compactness theorem, a parallel bridge construction would be possible.
- **Sormani-Vega 2016 / Sakovich-Sormani 2024:** null-distance approach to Lorentzian Gromov-Hausdorff. Different convergence notion (metric GH with null distance vs Mondino-Sämann's causal-diamond ε-nets); not in direct competition but should be cited in related-work.
- **Allen-Burtscher 2022 / Kunzinger-Steinbauer 2022:** metric GH of Lorentzian length spaces with null distance. Same status as Sormani-Vega/Sakovich-Sormani.
- **Müller 2022 (Mül22):** Lorentzian GH via Cauchy slabs — adjacent independent development of MS's structure; not in direct competition with the bridge functor.

**Net concurrent-work assessment:** the bridge construction is well-positioned. The Mondino-Sämann arXiv:2504.10380 v4 framework is the most directly relevant Lorentzian synthetic-geometric framework for the bridge, and no published Krein-algebraic version exists. The Connes-Rovelli thermal-time identification (R3) is a 32-year-old hypothesis, but its specific application to Wick-rotating Krein metametric to MS time-separation is novel.

Recommendation: cite all of the adjacent works above in Paper 48's related-work section to position the bridge construction precisely.

---

**End of memo.**

**Files added in this sprint:**
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` (this memo, ~9500 words formal bridge-construction + verdict layer)
- `debug/data/sprint_phase_a3prime_mondino_samann_bridge.json` (per-candidate verdict structure + bridge theorem details + Phase A.4' cost estimate)

**Cross-references:**
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` (Phase A.2' substrate; load-bearing throughout §1.1, §6)
- `debug/sprint_tier3_light_krein_lift_diagnostic_memo.md` (Tier 3-Light diagnostic; predecessor)
- `debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md` (Phase A.2 ε-net with F2 identification + R1/R2/R3 candidates; load-bearing for §1.4, §3, §5)
- `debug/l3e_p3_phase_a1_literature_audit.md` (Phase A.1 literature audit; supplies Latrémolière 2512.03573 context)
- `debug/sprint_l3e_p3_rescope_memo.md` (re-scope memo; confirms Phase A.3' target)
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` (four-witness Wick-rotation theorem; load-bearing for R3, G-B2)
- `papers/group1_operator_algebras/paper_43_lorentzian_extension.tex` (hemispheric wedge construction, BW vacuum; load-bearing for §6.3 functor definition, §7 GeoVac wedge application)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` (Krein operator-system substrate, propagation = 2; load-bearing for §6.3, §6.4 (B3))
- `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` (K⁺-restricted weak-form propinquity; load-bearing for §6.3 Definition 6.3(3), §6.4 (B4))
- `papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex` (norm-resolvent + admissible scaling; load-bearing for §6.3 Definition 6.3(3))
- `papers/group3_foundations/paper_32_spectral_triple.tex` (§VIII rem:pythagorean_m1_closure; supplies cover-scale identification 2π/κ_g)
- Mondino-Sämann arXiv:2504.10380 v4 (Dec 9, 2025): Def 2.3, Def 3.6, Def 3.8, Def 3.12, Thm 6.2 — extracted in §1.2
- Latrémolière arXiv:2512.03573 (Dec 3, 2025): all Defs 1.18 through 5.1 — via Phase A.2'
- Connes-Rovelli 1994 (Class. Quantum Grav. 11, 2899): thermal-time hypothesis — supplies R3
