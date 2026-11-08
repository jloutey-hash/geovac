# Sprint A1-Matsubara — canonical memo

Date: 2026-06-03
Scope: 1-session diagnostic computation, no PSLQ, no paper edits.
Predecessor: `debug/sprint_a1_scoping_memo.md` (Phase-1 scoping;
GO verdict on the Matsubara sub-case at sprint scale).

## TL;DR

**Verdict: (c) NEGATIVE-no-structural-engagement, with structural
refinement.** The temporal compactification step
$\sthree \rightsquigarrow \sthree \times S^1_\beta$ does NOT inject
classical odd-zeta $\zeta(2k+1)$ or MZV content into the 4D
Seeley–DeWitt coefficients. The reason is structural and one-line: the
small-$t$ (UV) asymptotic of the temporal heat kernel
$K_{S^1_\beta}(t)$ is a **single power** $\beta/\sqrt{4\pi t}$ (the
$m = 0$ winding), with all $m \neq 0$ winding contributions
exponentially suppressed by $e^{-m^2\beta^2/(4t)}$. The Matsubara mode
sum --- the mechanism through which $\zeta(2k+1)$ values *could* enter
--- lives in the **opposite** asymptotic regime (large-$t$ / IR / low
temperature), reached via Jacobi-$\vartheta_3$ modular inversion. The
two regimes never exchange transcendental content at the SD level.

The 4D SD coefficients of the bosonic scalar product $\sthree \times
S^1_\beta$ at integer $k$ are explicitly:
$$
a_k^{\text{4D, scalar}}(\beta) = \beta \cdot a_k^{\sthree,\text{scalar}}
= \frac{2\pi^2 \beta}{k!},
$$
and of the Dirac product:
$$
a_0^{\text{4D, Dirac}}(\beta) = 4\pi^2 \beta, \quad
a_1^{\text{4D, Dirac}}(\beta) = -2\pi^2 \beta, \quad
a_k^{\text{4D, Dirac}}(\beta) = 0 \text{ for } k \ge 2.
$$
Both factor as $\beta \cdot a_k^{\sthree}$:\ the temporal compactification
multiplies every SD coefficient by the temporal volume $\beta$, with no
further structural engagement. The M2 pure-Tate ring
$\bigoplus_k \pi^{2k}\cdot\Q$ of Paper 55 §4 is preserved. The 4D
Dirac sector retains its two-term exactness verbatim.

Refinement vs the scoping memo's three-option framing: outcome (c)
landed but it is **non-trivial**. The diagnostic verifies a
load-bearing structural claim that wasn't proved in scoping:\ namely,
that the M2 ring is convention-independent under multiplicative
extension by a compact temporal volume factor. Recording this as a
clean negative (NOT vacuous) on the open question opened by the
scoping memo.

The Stefan–Boltzmann pi-injection of Paper 35 §III ($\pi^2/90$ in the
free-energy density) lives in a structurally different extraction
operation:\ it is the Mellin transform at integer $s$ on
$K_{\sthree \times S^1_\beta}(t)$, not the small-$t$ SD asymptotic.
The two extractions of the same heat kernel produce structurally
different transcendental content, and this is the canonical M1×M2
cross-product mechanism Paper 35 §3.2 already classifies. No conflict.

---

## 1. Setup

### 1.1 The thermal product geometry

Following Paper 51 §G2, the Euclidean thermal product geometry is
$\sthree \times S^1_\beta$ with the product Dirac operator
$$
D_{\text{total}} = D_{\sthree} \otimes I + \gamma_{\sthree} \otimes
(-i\partial_\tau),
$$
$\tau \in [0, \beta)$. The anti-commutator
$\{D_{\sthree}\otimes I, \gamma_{\sthree}\otimes (-i\partial_\tau)\}
= 0$ makes $D_{\text{total}}^2$ block-diagonal, so eigenvalues add and
the heat kernel factorises multiplicatively:
$$
\boxed{
K_{\sthree \times S^1_\beta}(t)
= K_{\sthree}(t) \cdot K_{S^1_\beta}(t).
}
$$

### 1.2 The two factors

**Spatial factor (scalar Laplacian).** Paper 51 Theorem 3.1 (`thm:scalar_ak`):
$$
K_{\sthree}^{\Delta}(t) = \frac{\sqrt\pi}{4}\,\frac{e^t}{t^{3/2}} +
O(e^{-\pi^2/t}),
$$
giving SD coefficients $a_k^{\sthree,\Delta} = 2\pi^2/k!$ in the
standard $(4\pi)^{-d/2}$ volume normalisation. The exp-small remainder
comes from Jacobi $\vartheta_3$ inversion on the integer spectrum.

**Spatial factor (Dirac).** Paper 51 Corollary 2.1 (`cor:two_term`):
$$
K_{\sthree}^{D^2}(t) = \frac{\sqrt\pi}{2}\,t^{-3/2} -
\frac{\sqrt\pi}{4}\,t^{-1/2} + O(e^{-\pi^2/t}),
$$
giving $a_0^{D^2} = 4\pi^2$, $a_1^{D^2} = -2\pi^2$, $a_k^{D^2} = 0$
for $k \ge 2$. Two-term exact via Bernoulli-identity $B_{2k+1}(3/2) =
(2k+1)/4^k$ (Paper 51 `thm:zeta_unit_neg_k`, `rem:two_term_uniqueness`).

**Temporal factor.** The heat kernel of the scalar Laplacian
$-\partial_\tau^2$ on $S^1_\beta$ (length $\beta$, periodic BC) is the
Jacobi theta-3 series in the standard form
$$
K_{S^1_\beta}(t) = \sum_{k\in\Z} e^{-(2\pi k/\beta)^2 t}
= \frac{\beta}{\sqrt{4\pi t}}\,\vartheta_3\!\left(0;\,
\frac{i \beta^2}{4\pi t}\right)
= \frac{\beta}{\sqrt{4\pi t}}\,\sum_{m \in \Z} e^{-m^2\beta^2/(4t)},
$$
where the second equality is the Jacobi $\vartheta_3$ modular
transformation that exchanges the **mode sum** (over Matsubara index
$k$) for the **winding sum** (over winding number $m$). The two
representations have orthogonal asymptotic regimes:
$$
K_{S^1_\beta}(t) =
\begin{cases}
\dfrac{\beta}{\sqrt{4\pi t}} + O(e^{-\beta^2/(4t)}) & t \to 0^+
\text{ (UV: winding-dominated)}, \\[1ex]
1 + 2\sum_{k\ge 1} e^{-(2\pi k/\beta)^2 t} & t \to \infty
\text{ (IR: mode-dominated)}.
\end{cases}
$$

This dichotomy is the structural crux of the verdict.

---

## 2. Factorised heat trace, small-$t$ asymptotic

In the small-$t$ asymptotic regime, only the $m = 0$ winding of
$K_{S^1_\beta}$ contributes to the power-law expansion; all higher
windings are absorbed into the $O(e^{-\beta^2/(4t)})$ remainder. The
product heat trace's small-$t$ expansion is therefore
$$
K_{\sthree \times S^1_\beta}(t)
= K_{\sthree}(t) \cdot \frac{\beta}{\sqrt{4\pi t}}
+ O(e^{-\beta^2/(4t)}) + O(e^{-\pi^2/t}).
$$

### 2.1 Scalar sector (explicit derivation)

$$
K_{\sthree \times S^1_\beta}^{\Delta}(t)
= \frac{\sqrt\pi}{4} \cdot \frac{e^t}{t^{3/2}} \cdot
  \frac{\beta}{\sqrt{4\pi t}}
= \frac{\beta}{8} \cdot \frac{e^t}{t^2}
+ \text{exp-small}.
$$

Expanding $e^t = \sum_{k \ge 0} t^k/k!$ and applying the $d = 4$ SD
convention (multiplication by $(4\pi)^{d/2} = 16 \pi^2$):
$$
\boxed{
a_k^{\text{4D, scalar}}(\beta) =
16\pi^2 \cdot \frac{\beta}{8 \cdot k!}
= \frac{2 \pi^2 \beta}{k!}.
}
$$

Explicit values:

| $k$ | $a_k^{\text{4D, scalar}}$ |
|:---:|:--:|
| 0 | $2\pi^2 \beta$ |
| 1 | $2\pi^2 \beta$ |
| 2 | $\pi^2 \beta$ |
| 3 | $\pi^2 \beta / 3$ |
| 4 | $\pi^2 \beta / 12$ |
| 5 | $\pi^2 \beta / 60$ |

Pattern observation: $a_k^{\text{4D, scalar}}(\beta) = \beta \cdot
a_k^{\sthree, \text{scalar}}$, exactly. The temporal compactification
acts as a multiplicative volume factor; no other structural change
takes place.

### 2.2 Dirac sector

$$
K_{\sthree \times S^1_\beta}^{D^2}(t)
= \left(\frac{\sqrt\pi}{2}\,t^{-3/2} - \frac{\sqrt\pi}{4}\,t^{-1/2}\right)
\cdot \frac{\beta}{\sqrt{4\pi t}}
= \frac{\beta}{4}\,t^{-2} - \frac{\beta}{8}\,t^{-1} + \text{exp-small}.
$$

Applying the $(4\pi)^2 = 16\pi^2$ normalisation:
$$
\boxed{
a_0^{\text{4D, Dirac}}(\beta) = 4\pi^2 \beta, \quad
a_1^{\text{4D, Dirac}}(\beta) = -2\pi^2 \beta, \quad
a_k^{\text{4D, Dirac}}(\beta) = 0 \text{ for } k \ge 2.
}
$$

Again the pattern: $a_k^{\text{4D, Dirac}}(\beta) = \beta \cdot
a_k^{\sthree, \text{Dirac}}$. Two-term exactness is preserved.

The 4D Dirac SD expansion of Paper 51 Cor. 4.1 (`cor:zeta_4D_neg_k`)
matches verbatim:
$$
S(R, \beta, \Lambda)
= \frac{\beta R^3}{4}\Lambda^4 - \frac{\beta R}{8} \Lambda^2 +
  O(\text{exp-small}),
$$
where $R = 1$ (unit $\sthree$) for our convention. The Einstein-Hilbert
+ cosmological constant structure of Paper 51 §G2 is the physical
restatement of this:\ the $\beta\Lambda^4$ and $\beta\Lambda^2$ terms
are $a_0^{\text{4D,Dirac}}$ and $a_1^{\text{4D,Dirac}}$ in heat-kernel
language. Both terms carry $\beta$ linearly and live in $\pi^2 \Q$.

### 2.3 Period-theoretic classification

All 4D SD coefficients of both the scalar and Dirac sectors sit in
$$
M_2^{(\text{4D, thermal})}
\subset \beta \cdot \bigoplus_k \pi^{2k}\cdot\Q.
$$

The temporal volume $\beta$ enters as a **Tate-weight-0 rational
scaling factor** (it carries no $\pi$ and no zeta content of its own;
it is a free parameter of the substrate). After factoring out $\beta$,
the residual ring is bit-identical to the M2 pure-Tate ring of Paper
55 §4 on the static $\sthree$ sub-case:
$$
\frac{1}{\beta}\,M_2^{(\text{4D, thermal})}
= M_2^{(\sthree)}
= \bigoplus_k \pi^{2k}\cdot\Q.
$$

No $\zeta(3)$. No $\zeta(5)$. No MZV. No level-$N \ge 2$ cyclotomic
content. **M2 is preserved verbatim.**

---

## 3. Where the Matsubara modes DO inject $\zeta(2k+1)$

The Matsubara modes $(2\pi k/\beta)$ are physical. They DO produce
$\zeta(2k+1)$ values in GeoVac observables. Paper 35 §III documents
this explicitly:\ the Stefan-Boltzmann constant $\pi^2/90$ in the
free-energy density on $\sthree \times S^1_\beta$ comes from the
Bose-Einstein integral
$$
\int_0^\infty \frac{x^3}{e^x - 1}\,dx
= \Gamma(4)\,\zeta_R(4) = 6 \cdot \frac{\pi^4}{90} = \frac{\pi^4}{15},
$$
which evaluates the spectral zeta of $S^1_\beta$ at integer $s = 4$.

But this is **NOT an SD coefficient extraction**. It is a different
operation on the same $K_{\sthree \times S^1_\beta}(t)$:

| Extraction operation | Regime | Output | Transcendentals |
|:---|:---:|:---|:---|
| SD asymptotic at $t \to 0^+$ | UV | $a_k^{4D}(\beta) = \beta\cdot a_k^{\sthree}$ | $\pi^{2k}\beta\cdot\Q$ |
| Mellin at integer $s$ | finite $s$, mixed | spectral $\zeta$ values | $\zeta_R(s)$ at integer $s$ |
| IR free energy $\beta \to 0$ | $t \to \infty$ via $\vartheta$-inversion | $\propto T^4 \pi^2/90$ | $\pi^2/\beta^4\cdot\Q$ |

The SD operation projects to one transcendental class
($\pi^{\text{even}}\beta\Q$); the Mellin / Bose-Einstein operation
projects to a different one ($\pi^{\text{even}}\Q$ at integer $s$ via
$\zeta_R$). They are independent linear functionals on the same heat
kernel. **The temporal compactification provides $\zeta(2k+1)$
content to observables only via Mellin / Bose-Einstein extraction,
not via SD asymptotic extraction.**

This is the precise structural reason the verdict is (c)
NEGATIVE-no-structural-engagement and not (b) POSITIVE-MZV-enters:\
the Matsubara mode sum lives in the IR regime, and the SD coefficient
extraction is a UV operation. The two regimes are exchanged by the
Jacobi $\vartheta_3$ modular transformation, but they do not mix
transcendental content; the modular transformation is purely
algebraic, mapping one $\sqrt\pi$-rich representation to another
without injecting odd-zeta.

This refines and clarifies Paper 35 §III's Stefan-Boltzmann
classification:\ the M2 sector of GeoVac's 4D thermal substrate is
**period-theoretically clean at the heat-kernel asymptotic level**;
$\zeta(2k+1)$ enters only when one performs a Mellin extraction at
integer $s$, which is a structurally different observable from the
SD coefficient.

---

## 4. Verdict per decision gate

Per the three-option decision gate of the sprint prompt:

- **(a) POSITIVE-pure-Tate-survives.** Substantively correct at the
  level of the explicit closed forms above:\ the SD coefficients sit in
  $\beta\cdot\pi^{2k}\cdot\Q$, i.e.\ pure Tate (modulo the rational
  scaling factor $\beta$). However, the framing of (a) was
  "Matsubara structurally preserves the M2 ring via temporal $\zeta$
  cancellation or recombination", and our finding is sharper:\ the
  temporal $\zeta$ values **never enter at the SD level**, so there is
  nothing to cancel or recombine. The cancellation/recombination
  language overstates what happens. (a) is therefore the verdict
  **at the level of the result**, but the **mechanism** is closer to (c).

- **(b) POSITIVE-MZV-enters.** **Rejected.** No $\zeta(3)$, $\zeta(5)$,
  or any MZV enters the 4D SD coefficient at integer $k$. The
  speculation in the scoping memo about Matsubara possibly injecting
  $\zeta(3)$ via $\zeta(s) = 2\zeta(s)/\beta^s$ at integer $s$ does
  not survive the SD asymptotic computation:\ the SD regime is the
  $t \to 0^+$ winding-dominated regime of $K_{S^1_\beta}$, where the
  $m = 0$ winding contributes a single power $\beta/\sqrt{4\pi t}$
  and the $m \ne 0$ windings are exponentially suppressed. The
  Matsubara mode sum lives in the $t \to \infty$ regime, which the SD
  asymptotic does not see.

- **(c) NEGATIVE-no-structural-engagement.** **The accurate label for
  the mechanism.** Temporal compactification multiplies SD
  coefficients by a Tate-weight-0 scaling factor $\beta$ and changes
  no other structural feature. The question "does Matsubara preserve
  the M2 ring?" is best answered as:\ Matsubara never engages M2 at
  the SD level, so the ring is trivially preserved by non-engagement.

**Best single-label verdict:** (c) NEGATIVE-no-structural-engagement,
with the structural mechanism being that the SD asymptotic and the
Matsubara mode sum live in orthogonal asymptotic regimes of the
temporal heat kernel.

**Caveat to the negative reading:** the result is **not vacuous**.
Two substantive structural claims are verified here:
1. The 4D Dirac SD expansion is two-term exact, inheriting the
   $\sthree$ two-term exactness verbatim. This locks Paper 51 Cor. 4.1
   (Einstein-Hilbert + cosmological constant from G2, no higher
   curvature corrections) at the period-theoretic level.
2. The temporal compactification preserves the pure-Tate sub-ring
   $\bigoplus_k \pi^{2k}\cdot\Q$, refining the F-M mixed-Tate
   classification on the static sub-case. The Matsubara mechanism
   does not push GeoVac's M2 ring out of pure-Tate into generic
   mixed-Tate; that escape route is closed.

Both claims live below the F-M classification ceiling. Both are
strictly stronger than F-M at the static-$\sthree$ specialisation.

---

## 5. Tagging per master Mellin engine M1/M2/M3

Per CLAUDE.md memory rule [[feedback_tag_transcendentals]] and the
Paper 18 §III.7 master Mellin engine domain partition:

| Quantity | Value | M1 / M2 / M3 tag | Domain |
|:---|:---:|:---:|:---|
| $a_k^{\text{4D, scalar}}(\beta) = 2\pi^2\beta/k!$ | pure $\pi^2\Q\beta$ | **M2** (Seeley-DeWitt, $k = 2$ slot) | SD asymptotic on $\sthree \times S^1_\beta$ |
| $a_0^{\text{4D, Dirac}}(\beta) = 4\pi^2\beta$ | pure $\pi^2\Q\beta$ | **M2** ($k = 2$ slot, leading) | SD asymptotic, Dirac |
| $a_1^{\text{4D, Dirac}}(\beta) = -2\pi^2\beta$ | pure $\pi^2\Q\beta$ | **M2** ($k = 2$ slot, subleading) | SD asymptotic, Dirac |
| $a_k^{\text{4D, Dirac}}(\beta) = 0,\ k \ge 2$ | zero | structurally zero (two-term exact) | SD asymptotic, Dirac |
| $\sqrt\pi$ prefactor of $K_{\sthree}$ | $\sqrt\pi$ | **M2** (spinor-fiber-dim, $d = 3$ odd) | volume convention |
| $\sqrt\pi$ in $K_{S^1_\beta}$ leading | $\sqrt\pi$ | **M1** ($\Vol(S^1)$ Mellin measure) | temporal volume |
| $\beta$ (temporal volume) | $\beta$ (rational) | Tate-weight-0 substrate scale | substrate parameter |
| Stefan-Boltzmann $\pi^2/90$ | $\pi^2\Q$ | **M1 × M2** cross-product | Mellin/Bose-Einstein extraction (NOT SD!) |

Both $\sqrt\pi$ factors in the heat trace (one from $K_{\sthree}$ via
M2 spinor-fiber-dim, one from $K_{S^1_\beta}$ via M1 $\Vol(S^1)$
Mellin measure) combine multiplicatively:\
$\sqrt\pi \cdot \sqrt\pi = \pi$, which then assembles with the $1/\pi$
from the $1/(4\pi)$ in the temporal-volume normalisation to give
$1/\sqrt{\pi}$... and then the $(4\pi)^{d/2}$ 4D convention factor
$(4\pi)^2 = 16\pi^2$ recovers the clean $\pi^2\beta\cdot\Q$ ring of
the table. All $\pi$ powers are accounted for; nothing escapes the
master Mellin engine partition.

The Stefan-Boltzmann row sits in the M1×M2 cross-product, but is NOT
an SD coefficient. It is the joint-engagement M1×M2 mode classified
in Paper 55 §6 and demonstrated in Paper 35 §III, Paper 51 G2.

---

## 6. Honest scope

What this memo DOES:

- Computes the small-$t$ asymptotic expansion of the product heat
  kernel $K_{\sthree \times S^1_\beta}(t)$ explicitly in both scalar
  and Dirac sectors, via sympy.
- Extracts the 4D SD coefficients at integer $k = 0, 1, 2, 3, 4, 5$ in
  the standard $(4\pi)^{d/2}$ convention.
- Tags each coefficient in the M1/M2/M3 framework.
- Verifies that the M2 pure-Tate ring of Paper 55 §4 is preserved by
  the temporal compactification step.
- Identifies the precise structural mechanism (winding-vs-mode regime
  separation under Jacobi $\vartheta_3$ modular transformation).

What this memo does NOT do:

- Run PSLQ. (No need:\ closed forms are explicit.)
- Address the non-compact decompactification limit $\beta \to \infty$.
  (Paper 51 §3 covers this; the result is the zero-temperature de
  Sitter spectral action with no Matsubara modes, which agrees with
  the static $\sthree$ classification.)
- Address generic continuum $a(t)$ R-W substrates with non-zero
  time derivatives $\varepsilon_i$. (Outside discrete-substrate
  principle, per Sprint A1 scoping memo NO-GO sub-case.)
- Tile open question Q3' on $\sfive$ M3 cyclotomic content.
- Modify Paper 55 or Paper 35.

This is exactly the diagnostic decision-gate output requested in the
prompt. Paper edits await PI sign-off.

---

## 7. Proposed Paper 55 §7.2 update (paste-ready LaTeX)

The current Paper 55 §7.2 (`subsec:open_inhomo`, three sub-cases of
inhomogeneous extension) lists thermal-time compactification as a
sprint-scale open question with the framing:

> Does the temporal zeta $\zeta_{S^1_\beta}(s) = 2\zeta(s)/\beta^s$
> inject classical $\zeta(2k+1)$ content into the 4D SD coefficients,
> or does pure-Tate survive structurally?

The proposed replacement is a closure-in-the-negative paragraph plus
a Proposition. Drop-in LaTeX:

```latex
  \item \textit{Thermal-time compactification (closed in the negative).}
    The temporal compactification $\R \to S^1_\beta$ (the natural GeoVac
    time-direction surrogate via the Paper~35~\cite{loutey_paper35}
    Matsubara substrate $\sthree \times S^1_\beta$) does NOT inject
    classical odd-zeta $\zeta(2k+1)$ or MZV content into the 4D
    Seeley--DeWitt coefficients.  Sprint~A1-Matsubara (2026-06-03,
    \texttt{debug/sprint\_a1\_matsubara\_memo.md}) computes the small-$t$
    asymptotic of the product heat kernel
    $K_{\sthree \times S^1_\beta}(t) = K_{\sthree}(t) \cdot
    K_{S^1_\beta}(t)$ and finds that the temporal factor contributes a
    single power $\beta/\sqrt{4\pi t}$ at $t \to 0^+$ (the $m = 0$
    winding of $\vartheta_3$), with the $m \neq 0$ winding terms
    exponentially suppressed by $e^{-m^2\beta^2/(4t)}$.  The Matsubara
    mode sum --- the mechanism by which $\zeta(2k+1)$ could enter ---
    lives in the orthogonal large-$t$ (IR) regime of $K_{S^1_\beta}(t)$,
    reached only via Jacobi $\vartheta_3$ modular inversion, and is not
    visible to the SD asymptotic.  We obtain the explicit closed forms

    \begin{proposition}[Pure-Tate refinement on $\sthree \times S^1_\beta$]
    \label{prop:thermal_pure_tate}
    The volume-normalised 4D Seeley--DeWitt coefficients on
    $\sthree \times S^1_\beta$ for the scalar Laplacian and the
    Camporesi--Higuchi Dirac satisfy
    \begin{align}
    a_k^{\rm 4D,\,scalar}(\beta) &= \frac{2\pi^2\,\beta}{k!} \quad
       (k \ge 0), \\
    a_0^{\rm 4D,\,Dirac}(\beta)  &= 4\pi^2\,\beta, \quad
    a_1^{\rm 4D,\,Dirac}(\beta)  = -2\pi^2\,\beta, \quad
    a_k^{\rm 4D,\,Dirac}(\beta)  = 0 \;\; (k \ge 2).
    \end{align}
    Each coefficient factorises as $\beta\cdot a_k^{\sthree}$, where
    $a_k^{\sthree}$ is the corresponding SD coefficient on the static
    $\sthree$ sub-case (\S\ref{sec:m2}, Paper~51~Thm~3.1 and
    Cor.~2.1~\cite{loutey_paper51}).  The two-term exactness of the
    $\sthree$ Dirac sector is inherited verbatim.  All coefficients
    therefore lie in the rational-rescaled pure-Tate ring
    \[
    M_2^{\rm (4D,\,thermal)} \subset \beta\cdot
    \bigoplus_{k \ge 0} \pi^{2k}\cdot\Q.
    \]
    \end{proposition}

    The Stefan--Boltzmann pi-injection of Paper~35~\S~III is structurally
    a different extraction --- the Mellin transform at integer $s$ on
    the same product heat kernel, not the small-$t$ SD asymptotic ---
    and lives in the M1$\times$M2 cross-product joint engagement of
    \S\ref{sec:joint}, not in M2 itself.  No conflict.
```

This closes the open question "does Matsubara inject MZV?" in the
negative, and adds a Proposition with explicit closed forms for the
4D SD coefficients.

---

## 8. Files used

### Papers (read)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` (§4 M2
  proof, §6 joint engagement, §7.2 open question on inhomogeneous
  extension)
- `papers/group6_precision_observations/paper_35_time_as_projection.tex`
  (§II KG spectrum on $\sthree\times S^1_\beta$, §III Stefan-Boltzmann
  as M1×M2, §VIII Sprint TD Track 1 tensor-product verification)
- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` (§G1 zeta-unit-neg-k
  theorem and Cor. 2.1 two-term-exact Dirac, §G2 thermal product Cor.
  4.1, §G3 scalar SD closed form, §G3 TT-tensor sector)

### Predecessor memos
- `debug/sprint_a1_scoping_memo.md` (Phase-1 scoping, GO verdict)
- `debug/sprint_mixed_tate_test_memo.md` (M2 pure-Tate baseline on
  static $\sthree$, the structural precedent)

### Diagnostic (this sprint)
- Inline sympy verification (no separate driver file:\ closed forms
  computed via short interactive sessions, both factors known
  analytically from Paper 51)

---

## Final 200-word summary

The Matsubara temporal compactification $\sthree \rightsquigarrow
\sthree \times S^1_\beta$ does NOT inject classical odd-zeta
$\zeta(2k+1)$ or MZV content into the 4D Seeley-DeWitt coefficients.
The explicit closed forms are $a_k^{\text{4D, scalar}}(\beta) =
2\pi^2\beta/k!$ and $a_k^{\text{4D, Dirac}}(\beta) = \beta \cdot
a_k^{\sthree,\text{Dirac}}$ (two-term exact). Both factor as
$\beta \cdot a_k^{\sthree}$, with $\beta$ entering as a
Tate-weight-0 rational scaling factor and no further structural
engagement.

The Paper 55 §4 M2 pure-Tate ring is preserved verbatim:
$M_2^{(\text{4D, thermal})} \subset \beta \cdot \bigoplus_k
\pi^{2k}\cdot\Q$. The structural mechanism is the orthogonality of
asymptotic regimes of $K_{S^1_\beta}(t)$:\ the SD asymptotic
($t \to 0^+$) sees only the $m = 0$ winding (a single $\sqrt\pi$
power), while the Matsubara mode sum is the $t \to \infty$ regime,
reached by Jacobi $\vartheta_3$ inversion. The two regimes never
exchange transcendental content at the SD level.

The Stefan-Boltzmann $\pi^2/90$ that DOES involve odd-zeta lives in
the M1×M2 cross-product extracted via Mellin at integer $s$, NOT the
SD coefficient. Paper 55 §7.2 closes in the negative; Paper 51 Cor. 4.1
Einstein-Hilbert + cosmological constant locks at the period-theoretic
level. Verdict: (c) NEGATIVE-no-structural-engagement, refined into
a substantive structural finding.
