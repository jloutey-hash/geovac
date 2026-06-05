# Sprint memo: Propinquity-aware Trotter bounds for GeoVac composed Hamiltonians

**Date:** 2026-06-04
**Track:** Bridge math.OA propinquity machinery (Paper 38) to quantum-simulation cost (Paper 14/20)
**Goal:** Derive an operator-algebraic Trotter error bound from the propinquity bound $\Lambda(T_{n_{\max}}, T_{S^3}) \le C_3 \cdot \gamma_{n_{\max}}$. Numerically compare against naive Suzuki-Trotter for LiH, BeH$_2$, H$_2$O at production $n_{\max} = 2$.

## 1. Set-up

Let $H_\infty$ denote the (hypothetical) GeoVac composed Hamiltonian in the $n_{\max} \to \infty$ continuum limit on the round $S^3$ Camporesi-Higuchi spectral triple $T_{S^3}$. Let $H_{n_{\max}}^{\text{composed}}$ denote the production-grade truncated qubit Hamiltonian from `geovac.ecosystem_export.hamiltonian(...)` at finite Fock cutoff $n_{\max}$ on the truncated triple $T_{n_{\max}}$.

A quantum simulator runs a first-order Suzuki-Trotter product formula approximating
$$
e^{-i H_{n_{\max}} t} \approx U_{\text{Trotter}}^{r}(H_{n_{\max}}, t) := \prod_{k=1}^{r} \prod_{j} e^{-i c_j P_j t/r},
$$
where $H_{n_{\max}} = \sum_j c_j P_j$ is the Jordan-Wigner Pauli decomposition. We want to bound
$$
\bigl\| e^{-i H_\infty t} - U_{\text{Trotter}}^{r}(H_{n_{\max}}, t) \bigr\|_{\text{op}} \le \epsilon.
$$

By the triangle inequality,
$$
\epsilon \le \underbrace{\bigl\| e^{-i H_\infty t} - e^{-i H_{n_{\max}} t} \bigr\|}_{=: \epsilon_{\text{trunc}}}
       + \underbrace{\bigl\| e^{-i H_{n_{\max}} t} - U_{\text{Trotter}}^{r} \bigr\|}_{=: \epsilon_{\text{Trotter}}}.
$$
Naive Suzuki-Trotter analysis budgets the full $\epsilon$ to $\epsilon_{\text{Trotter}}$ and silently assumes $\epsilon_{\text{trunc}} = 0$. The propinquity from Paper 38 supplies the missing ingredient for budgeting $\epsilon_{\text{trunc}}$ honestly.

## 2. Truncation bound via propinquity

Paper 38, Lemma L3 ($C_3 = 1$) plus Lemma L4(d) gives, for any test function $f \in C^\infty(S^3)$,
$$
\bigl\| [D_{\mathrm{CH}}, M_f - B_{n_{\max}}(f)] \bigr\|_{\text{op}} \le \gamma_{n_{\max}} \cdot \|\nabla f\|_{L^\infty},
$$
and via Paper 38, Lemma L4(c),
$$
\bigl\| M_f - P_{n_{\max}} M_f P_{n_{\max}} \bigr\|_{\text{op}} \le \gamma_{n_{\max}} \cdot \|\nabla f\|_{L^\infty}.
$$
These are single-particle bounds on the round $S^3$.

The GeoVac composed Hamiltonian is an antisymmetric sum of single-particle and two-particle operators (`geovac/composed_qubit.py`), each factoring as a single-particle multiplier/Dirac truncation acting on one fermion line. By the standard $N_a$-body norm bound on fermionic 1RDM + 2RDM operators (Childs-Su-Tran-Wiebe-Zhu 2021 §IV), per-electron truncation propagates linearly:
$$
\bigl\| H_\infty - H_{n_{\max}}^{\text{composed}} \bigr\|_{\text{op}} \le N_{\text{active}} \cdot \gamma_{n_{\max}} \cdot L_H,
$$
where $L_H$ is the effective Lipschitz constant of the dominant single-particle multiplier. We use $L_H = \text{mean}\{|c_j|\}$ (mean Pauli coefficient) as a falsifier-safe representative of the typical multiplier strength. The max coefficient is dominated by the identity term and gives a vastly looser bound; the median is closer to the true Lipschitz contribution but is overly aggressive. The mean lands between the two.

Combining with Duhamel ($\| e^{-iAt} - e^{-iBt}\| \le t \|A-B\|$):
$$
\boxed{\;\;\epsilon_{\text{trunc}} \le t \cdot N_{\text{active}} \cdot \gamma_{n_{\max}} \cdot L_H.\;\;}
$$
This is a **propinquity-derived truncation bound**: the central-Fejer kernel rate $\gamma_{n_{\max}}$ (Paper 38 Lemma L2) appears as the *Lipschitz-distortion* control on the multiplier-Dirac difference. It vanishes as $n_{\max} \to \infty$ at rate $\gamma_{n_{\max}} \sim (4/\pi)\log n / n$, with **leading constant $4/\pi = \text{Vol}(S^2)/\pi^2$ — the M1 (Hopf-base measure) signature of Paper 18 §III.7's master Mellin engine**. The Trotter convergence rate carries the same transcendental signature that the framework's heat-kernel expansions imprint on its bulk observables. This is, on its own, a clean structural connection between propinquity geometry and Trotter cost — even though it does not produce a tighter $r$.

## 3. Trotter step bound

Commutator-tightened first-order Suzuki-Trotter (Childs-Su-Tran-Wiebe-Zhu 2021 Eq. 48):
$$
\epsilon_{\text{Trotter}} \le \frac{t^2}{r} \cdot \alpha_{\text{comm}},
\qquad \alpha_{\text{comm}} := \sum_{j<k:\,P_j\, \text{anticomm}\, P_k} |c_j| |c_k|.
$$
The classical 1-norm bound replaces $\alpha_{\text{comm}}$ with $\lambda^2/2$ where $\lambda = \sum_j |c_j|$.

## 4. Joint optimization

Split the error budget evenly: $\epsilon_{\text{trunc}}, \epsilon_{\text{Trotter}} = \epsilon/2$.

**Truncation feasibility:** $n_{\max}$ must satisfy $t N_{\text{active}} \gamma_{n_{\max}} L_H \le \epsilon/2$. If not, **no number of Trotter steps achieves the target error at this $n_{\max}$**.

**Trotter step count (propinquity-aware):**
$$
\boxed{\;\;r_{\text{prop}} = \Bigl\lceil \frac{2 t^2 \alpha_{\text{comm}}}{\epsilon} \Bigr\rceil
        \quad \text{(given truncation feasible)}.\;\;}
$$

Compare to:
- $r_{\text{naive}} = \lceil t \lambda / \sqrt{2\epsilon} \rceil$ (1-norm bound, no truncation budget)
- $r_{\text{comm}} = \lceil \sqrt{t^2 \alpha_{\text{comm}} / \epsilon} \rceil$ (commutator bound, no truncation budget)

## 5. Numerical results

All computations: `python debug/sprint_trotter_propinquity.py`, JSON at `debug/data/sprint_trotter_propinquity.json`, console log at `debug/data/sprint_trotter_propinquity_console.log`. Simulation time $t = 1.0$ atomic units. Propinquity values from `geovac.gh_convergence.compute_propinquity_bound`: $\gamma_{2} = 2.0746$, $\gamma_3 = 1.6101$, $\gamma_4 = 1.3223$.

### 5.1 Production panel at $n_{\max} = 2$, $\epsilon = 10^{-3}$

| Mol     | Q  | n_terms | $\lambda$ | $L_H$ (mean) | $\alpha_{\text{comm}}$ | $\epsilon_{\text{trunc}}$ bound | $r_{\text{naive}}$ | $r_{\text{comm}}$ | $r_{\text{prop}}$ | Feasible? |
|:--------|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| LiH     | 30 | 334 | 32.59 | 0.098 | 6.97  | 0.40 | 729 | 84  | N/A | NO |
| BeH$_2$ | 50 | 556 | 66.04 | 0.119 | 17.34 | 0.99 | 1,477 | 132 | N/A | NO |
| H$_2$O  | 70 | 778 | 358.91 | 0.461 | 188.41 | 7.66 | 8,026 | 435 | N/A | NO |

The truncation budget $\epsilon/2 = 5 \times 10^{-4}$ is **overshot by 800x for LiH, 2000x for BeH$_2$, 15,000x for H$_2$O** at the production cutoff. The propinquity-aware $r_{\text{prop}}$ is therefore N/A at production parameters for all three molecules.

### 5.2 n_max feasibility sweep for LiH at $\epsilon = 10^{-3}$, $t = 1$ a.u.

| $n_{\max}$ | $\gamma_{n_{\max}}$ | $\epsilon_{\text{trunc}}$ bound | Feasible ($\le 5{\times}10^{-4}$)? | $r_{\text{prop}}$ |
|---:|---:|---:|:---:|---:|
| 2     | 2.0746 | 0.4049 | NO | — |
| 3     | 1.6101 | 0.3142 | NO | — |
| 4     | 1.3223 | 0.2581 | NO | — |
| 5     | 1.1302 | 0.2206 | NO | — |
| 10    | 0.6724 | 0.1312 | NO | — |
| 20    | 0.3868 | 0.0755 | NO | — |
| 50    | 0.1800 | 0.0351 | NO | — |
| 100   | 0.0586 | 0.0114 | NO | — |
| 200   | 0.0337 | 0.0066 | NO | — |
| 500   | 0.0158 | 0.0031 | NO | — |
| 1000  | 0.0088 | 0.0017 | NO | — |
| 5000  | 0.0022 | 0.0004 | **YES** | 13,934 |

**The feasibility crossover for LiH is at $n_{\max} \approx 5000$**, where $\gamma_{n_{\max}}$ has shrunk by $\sim 1000\times$ relative to production. At this cutoff $r_{\text{prop}} = 13,934$, compared to naive $r_{\text{naive}} = 729$ and commutator $r_{\text{comm}} = 84$. The propinquity-aware bound at feasibility is $\sim 19\times$ MORE expensive than naive Suzuki-Trotter and $\sim 166\times$ MORE expensive than commutator-tightened Trotter.

### 5.3 BeH$_2$ and H$_2$O scale even worse

BeH$_2$ ($N_{\text{active}} = 4$) needs $n_{\max} \sim 10{,}000$ for $\epsilon = 10^{-3}$ feasibility; H$_2$O ($N_{\text{active}} = 8$, $L_H = 0.46$) needs $n_{\max} \sim 10^5$. Building these cutoffs is utterly out of computational reach (would require billions of basis vectors).

## 6. Diagnostic of the loose bound

The propinquity bound massively over-estimates the actual truncation error. Empirically, LiH at $n_{\max} = 2$ gives FCI energy at $\sim 2\%$ error vs exact (Paper 20 Table I), which corresponds to a Hamiltonian-norm truncation error on the order of $|E_{\text{FCI}} - E_{\text{exact}}|/N_{\text{shells}} \sim 0.01$ Ha — three orders of magnitude smaller than the propinquity-derived 0.40 upper bound.

Three sources of slack:

1. **$L_H$ is a single-particle Lipschitz constant**, but the composed Hamiltonian's truncation error is dominated by the leading-shell ERI contributions, which are typically much smaller (the Gaunt-restricted ERI density is $\sim 1.4\%$ at $\ell_{\max} = 2$, Paper 22).
2. **$\gamma_{n_{\max}}$ has a leading constant $4/\pi \approx 1.27$**, but the relevant test functions for the composed Hamiltonian are concentrated on the lowest few harmonics, where the actual Stein-Weiss constant is much smaller than the all-functions sup.
3. **The Duhamel + linearity step assumes the composed Hamiltonian decomposes into commuting single-particle blocks**, but in reality the cross-electron ERIs couple Fock shells in a way that creates partial cancellation (the basis-intrinsic sparsity from Gaunt selection rules, Paper 22, Paper 14).

A sharper $L_H$ and a sharper $\gamma_{n_{\max}}$ would each shrink the bound by $\sim 10\times$ — together $\sim 100\times$ — bringing LiH feasibility crossover from $n_{\max} \sim 5000$ to $n_{\max} \sim 50$. Even then, the propinquity-aware bound would not beat the commutator bound on $r$.

## 7. Honest scope statement

1. **The propinquity bound is a single-particle Lipschitz-distortion bound.** It is genuinely smaller as $n_{\max} \to \infty$ (the $\gamma_{n_{\max}} \to 0$ statement from Paper 38), but at all production-relevant cutoffs the structural upper bound $N_{\text{active}} \cdot \gamma_{n_{\max}} \cdot L_H$ is loose by ~3 orders of magnitude relative to actual FCI truncation error.

2. **The bound DOES NOT beat naive Suzuki-Trotter on Trotter steps.** Even where feasible ($n_{\max} = 5000$ for LiH), $r_{\text{prop}}$ is 19x larger than $r_{\text{naive}}$ and 166x larger than $r_{\text{comm}}$.

3. **What the bound DOES supply is a feasibility certificate.** It names the minimum $n_{\max}$ at which a target $\epsilon$ is achievable, separate from Trotter cost. This is genuinely new operator-algebraic content: propinquity is the natural object for this question, even though the bound is too loose to be useful in practice at production parameters.

4. **Named obstruction (loose $L_H$).** Using $L_H = \text{mean}\{|c_j|\}$ is conservative. A sharper $L_H$ from the actual multipole structure of the composed Hamiltonian (Paper 19) would shrink the bound by $\sim 10\times$. Still does not close the 1000x gap to make production $n_{\max} = 2$ feasible at $\epsilon = 10^{-3}$.

5. **Named obstruction (single-particle-to-many-body).** The Duhamel-plus-linearity step is loose by $\sim 2\times$ for exchange-correlation terms. Tightening would require extending Paper 39's tensor-product propinquity to antisymmetric $N$-electron operators — multi-month NCG work, not sprint-scale. Even with a factor 2 improvement, the bound is structurally insufficient.

6. **Named obstruction (Lipschitz constants in L2 not sharp).** Paper 38 Lemma L2 quantitative rate is $\gamma_{n_{\max}} = (4/\pi)\log n_{\max}/n_{\max} + O(1/n_{\max})$ with a UNIFORM upper bound $6\log n_{\max}/n_{\max}$. For the specific multipliers entering the composed Hamiltonian (which are NOT all of $C^\infty(S^3)$ but only the low-harmonic component), a much sharper rate is plausible. This is the L2 Stein-Weiss "Track C" sharpening flagged in `geovac/gh_convergence.py`. A factor 10 improvement is plausible but not yet executed.

## 8. Decision: STOP

**STOP because** the propinquity bound on $\|H_\infty - H_{n_{\max}}\|_{\text{op}}$ is structurally too loose — at the production cutoff $n_{\max} = 2$, the bound overshoots empirical FCI error by 3 orders of magnitude and overshoots the eps/2 truncation budget by $10^3$ to $10^4$ across LiH/BeH$_2$/H$_2$O. The first molecule reaches feasibility only at $n_{\max} \approx 5000$ (computationally inaccessible), and at that point requires $r_{\text{prop}} = 13{,}934$ vs $r_{\text{naive}} = 729$ — the propinquity-aware bound is uniformly LOOSER than naive Suzuki-Trotter when both are evaluated honestly.

The structural finding stands as a NEGATIVE: the operator-algebraic propinquity, which gives a beautiful convergence theorem for spectral triples (Paper 38), does not transfer to a useful Trotter cost certificate at production parameters via the canonical Duhamel-plus-linearity lift. The gap is concentrated in three places: (i) sharper $L_H$, (ii) sharper $\gamma_{n_{\max}}$ on the low-harmonic component, (iii) per-Pauli-term local propinquity in place of the global single-particle bound. Each is a clean future research target, but none is sprint-scale and the combined improvement needed (~$10^4$) is multi-year.

The headline content suitable for a Paper 20 footnote is the *negative structural finding* itself: "propinquity is the natural object for this question, but the bound is currently too loose to compete with Suzuki-Trotter at production parameters; tightening is named here as a multi-step research program rather than a near-term result."

## 9. Proposed Paper 20 update (small footnote, NEGATIVE FINDING)

**Target paper:** Paper 20 (resource benchmarks). Better fit than Paper 14 (which is structural-sparsity).

**Target section:** Insert as a single-paragraph footnote attached to the existing §III.B "1-norm and Trotter step counts" subsection (or its equivalent). DO NOT add a full subsection — the finding is a NEGATIVE, and structural negatives are better captured as honest footnotes than as standalone subsections.

**Proposed footnote text** (PM should apply after audit; do NOT edit .tex directly per sprint instructions; avoids merge conflict with Z2-tapering track):

```latex
\footnote{
The GH-convergence theorem of \cite{GeoVac_Paper38}
(propinquity $\Lambda(T_{n_{\max}}, T_{S^3}) \le C_3 \gamma_{n_{\max}}$ with
$\gamma_{n_{\max}} \sim (4/\pi)\log n_{\max}/n_{\max}$) gives a
single-particle Lipschitz-distortion bound on the truncation error
$\|H_\infty - H_{n_{\max}}\|_{\text{op}} \le N_{\text{active}}
\gamma_{n_{\max}} L_H$. Combining with the standard commutator-tightened
Trotter bound \cite{ChildsSu2021} via a split error budget produces a
\emph{feasibility certificate}: the minimum $n_{\max}$ at which a
target $\epsilon$ is achievable, separate from the Trotter step count.
At production $n_{\max} = 2$ and target $\epsilon = 10^{-3}$ the
truncation budget is overshot by three to four orders of magnitude
across LiH, BeH$_2$, H$_2$O; LiH first becomes feasible at $n_{\max}
\approx 5000$ with $r_{\text{prop}} = 13{,}934$ steps vs
$r_{\text{naive}} = 729$ for the same molecule under standard 1-norm
analysis. The structural conclusion is that the propinquity bound is
currently too loose to compete with Suzuki-Trotter at production
parameters; tightening would require a sharper Lipschitz constant
$L_H$ (a multipole-decomposition refinement of Paper~19) and a sharper
Stein-Weiss rate on the low-harmonic component (the L2 "Track C"
sharpening of Paper~38) and is named here as a multi-step research
program. The $4/\pi$ constant identifies the M1 (Hopf-base measure)
signature of the master Mellin engine of \cite{GeoVac_Paper18}, so
sharpening the propinquity-Trotter bridge connects to the broader
master-Mellin-engine program; this is a positive structural reading
of an otherwise negative finding.
}
```

## Verdict

**STOP because** the propinquity-derived truncation bound is uniformly looser than naive Suzuki-Trotter at production parameters (3–4 orders of magnitude over the $\epsilon/2$ budget for LiH/BeH$_2$/H$_2$O at $n_{\max} = 2$, $\epsilon = 10^{-3}$; LiH's feasibility crossover is at $n_{\max} \approx 5000$ where $r_{\text{prop}}$ is 19× the naive count), and the named sharpening paths (sharper $L_H$ via multipole decomposition; sharper $\gamma_{n_{\max}}$ via L2 Stein-Weiss "Track C") are each multi-step research programs rather than sprint-scale tightening — the structural finding is genuinely new (single-particle propinquity → many-body Trotter via Duhamel + $N_{\text{active}}$-linearity carries the $4/\pi$ M1 Hopf-base measure signature) but the numerical bound is not useful for production simulation.
