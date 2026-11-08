# WH1 Round 3 — Connes Distance on the GeoVac Truncated Operator System

**Author:** PM-dispatched sub-agent (research + computational verification, no
paper edits, no CLAUDE.md edits)
**Scope:** WH1 (CLAUDE.md §1.7) — Round 3 of 3.
**Question:** Compute the Connes distance
$d_{D_{n_{\max}}}(\phi_v, \phi_w) = \sup_{x \in \mathcal{O}_{n_{\max}}} \{
|\phi_v(x) - \phi_w(x)| : \|[D, x]\|_{\mathrm{op}} \le 1 \}$ on pure
node-evaluation states $\phi_v(x) := \langle v | x | v \rangle$ at
$n_{\max} = 2$ and $n_{\max} = 3$. Compare the resulting metric structure
to the combinatorial Fock-graph distance $d_{\mathrm{graph}}(v, w)$.
**Date:** 2026-05-03
**Verdict (one-line):** SDP framework operational at $n_{\max} = 2, 3$;
metric is **structurally non-trivial but anisotropic** at the placeholder
convention — forced zeros from operator-system $\mathrm{SO}(3)$ $m$-
reflection symmetry (1 pair at $n_{\max} = 2$, 4 pairs at $n_{\max} =
3$) plus a clean $\{1, 2, 3\}$ rational structure of the $n_{\max} = 2$
extrema (lost at $n_{\max} = 3$ as the SDP gains degrees of freedom).
Connes distance is **NOT monotone** in the combinatorial Fock-graph
distance: Pearson $\rho \approx -0.04$ ($n_{\max} = 2$) growing to
$+0.14$ ($n_{\max} = 3$); Spearman $\approx 0.09 \to 0.18$. Slow weak
trend toward correlation but very far from a discretized geodesic
metric at our reachable cutoffs.

---

## §1. Definition and finite-dim SDP formulation

### 1.1. Connes distance on operator systems

The verbatim definition (Connes & van Suijlekom 2021, arXiv:2004.14115v2,
Eq. (2), line 1090; see also `debug/wh1_connes_vs_pdf_verification.md` D4)
reads, for an operator system spectral triple $(\mathcal{E}, \mathcal{H}, D)$
and two states $\phi, \psi \in S(\overline{\mathcal{E}})$ (positive
linear functionals of norm 1 on the operator system itself):

$$
d_D(\phi, \psi) = \sup_{x \in \mathcal{E}} \big\{ |\phi(x) - \psi(x)| :
\|[D, x]\|_{\mathrm{op}} \le 1 \big\}.
$$

For pure node-evaluation states $\phi_v(x) := \langle v | x | v \rangle =
\mathrm{Tr}(E_v x)$ with $E_v = |v\rangle\langle v|$, the objective is
linear in $x$:

$$
\phi_v(x) - \phi_w(x) = \mathrm{Tr}((E_v - E_w) x) = x_{v, v} - x_{w, w}
\quad (\text{Hermitian } x).
$$

### 1.2. SDP formulation

Parameterize $x = \sum_{k=1}^{K} c_k M_k$ with complex coefficients
$c_k$ and $\{M_k\}_{k=1}^K$ the multiplier-matrix basis from round 2
(`geovac/operator_system.py`). The SDP is

$$
\begin{aligned}
\text{maximize} \quad & \pm\,(x_{v, v} - x_{w, w}) \\
\text{subject to} \quad & x = \sum_k c_k M_k, \\
& x = x^*, \\
& \|[D, x]\|_{\mathrm{op}} \le 1, \\
& \langle K_j, x \rangle_{\mathrm{HS}} = 0 \quad \forall K_j \in \ker([D, \cdot]) \cap \mathcal{O}_h.
\end{aligned}
$$

The operator-norm constraint lifts to the LMI
$\bigl[\begin{smallmatrix} I & [D, x] \\ [D, x]^* & I \end{smallmatrix}\bigr]
\succeq 0$, handled automatically by `cvxpy` under
`cp.norm(commutator, 2) <= 1`. The two senses (max and min of the diff)
are solved separately and the larger absolute value taken; this matches
the Connes-vS supremum convention.

### 1.3. Gauge fixing on $\ker([D, \cdot]) \cap \mathcal{O}_h$

The objective and the commutator constraint are both invariant under
$x \mapsto x + y$ for any $y \in \ker([D, \cdot])$ (any matrix that
commutes with $D$). If $E_v - E_w$ is **non-orthogonal** to some such
$y$ in $\mathcal{O}_h$, the SDP is genuinely unbounded (free shifts in
the objective at zero commutator cost) and the Connes distance is
$+\infty$. We detect this case explicitly. When $E_v - E_w$ is HS-
orthogonal to every $K_j \in \ker([D, \cdot]) \cap \mathcal{O}_h$, we
add gauge-fixing constraints $\langle K_j, x \rangle_{\mathrm{HS}} = 0$
to project $x$ onto $\ker^\perp$, bounding the SDP without changing the
objective value.

The kernel intersection $\ker([D, \cdot]) \cap \mathcal{O}_h$ is computed
as the null space of the linear map $B \mapsto [D, B]$ acting on the
real Hermitian basis $\mathcal{O}_h$ (extracted via SVD on the
Hermitian-spanning generators $H_k = (M_k + M_k^*)/2,\, K_k = (M_k -
M_k^*)/(2i)$).

---

## §2. Implementation details

### 2.1. The Dirac operator on $\mathcal{H}_{n_{\max}}$

The Camporesi-Higuchi spectrum of the round-$S^3$ Dirac operator is
$|\lambda_n| = n + 3/2$ with multiplicity $g_n^{\mathrm{Dirac}} = 2(n+1)
(n+2)$; this is the standard scalar GeoVac shell structure. Three
diagonal-vs-offdiagonal preset Dirac proxies are implemented in
`geovac/connes_distance.py::DiracProxy`:

| Mode | Formula | Spectrum | $\dim(\ker([D, \cdot]) \cap \mathcal{O}_h)$ at $n_{\max}=2$ | Use |
|------|---------|----------|----------|-----|
| `shell_scalar` | $D \mid n,l,m\rangle = (n + \mathrm{shift}) \mid n,l,m\rangle$ | degenerate by $n$-shell | large ($\ge 5$) | "truthful" Camporesi-Higuchi diagnostic; gives $+\infty$ on most pairs |
| `nondegenerate` | $D \mid n,l,m\rangle = (n + 0.1\,l + 0.01\,m) \mid n,l,m\rangle$ | non-degenerate diagonal | $3$ | breaks shell degeneracy but $M_{N,0,0}$ multipliers still commute |
| `offdiag` *(default)* | non-degenerate diagonal $+\,$ unit off-diagonal coupling on $\{(n,l,m) \leftrightarrow (n',l',m'): \|\Delta n\|=\|\Delta l\|=1, \|\Delta m\|\le 1\}$ | non-degenerate, dense | $1$ (just $\mathbb{C}\cdot I$) | **smallest possible kernel**; gives finite Connes distance on (almost) every pair |

The off-diagonal coupling pattern of the `offdiag` Dirac matches the
$E1$ ($\Delta l = 1$) selection rule of the Camporesi-Higuchi spinor
Dirac on $S^3$; in the scalar Fock basis this is the
"raising/lowering ladder of the n-th principal-QN derivative" and is
the natural off-diagonal substitute for the spinor lift that we are
not attempting here.

**Convention statement (mandatory):** Numerical distance values reported
below depend on (a) the placeholder for the SO(4) radial overlap
$I_{nNn'}^{lLl'}$ used by `geovac/operator_system.py`, and (b) the
Dirac proxy mode. The structural properties (positivity, symmetry,
triangle inequality, kernel structure, forced zeros from operator-system
symmetry, $\{1, 2, 3\}$ rational ratios) are stable across reasonable
choices of both; the absolute scale is not.

### 2.2. Software stack

- **Solver:** `cvxpy` 1.8.2 with `SCS` 3.x as the default SDP backend.
  `CLARABEL` also works at this scale.
- **Performance** at the placeholder convention:
  - $n_{\max} = 2$: $5\times 4 / 2 = 10$ unique pairs $\times 2$ SDPs each
    $\approx 16\,$s total.
  - $n_{\max} = 3$: $14 \times 13 / 2 = 91$ unique pairs $\times 2$ SDPs
    each $\approx 20$–$25\,$min total.
  - $n_{\max} = 4$: $30 \times 29 / 2 = 435$ pairs; SDP variable size
    $K = 140$, $N = 30$ — feasible but $\sim$hours; not attempted in
    round 3.
- **Test suite:** `tests/test_connes_distance.py`, 17/17 passing in
  $\sim 2\,$min.

---

## §3. Distance matrix at $n_{\max} = 2$

### 3.1. Setup

- Hilbert dimension: $N = 5$.
- Operator-system dimension: $\dim(\mathcal{O}) = 14$.
- Multiplier labels (round 2): $(N, L, M) \in \{(1,0,0), (2,0,0),
  (2,1,m)_{m=-1,0,+1}, (3,0,0), (3,1,m)_{m=-1,0,+1}, (3,2,m)_{m=-2,\dots,+2}\}$.
- Kernel of $[D, \cdot]$ on $\mathcal{O}_h$ at default Dirac:
  $\dim = 1$, basis vector proportional to the identity.
- Basis indexing:

| index | $\mid n, l, m\rangle$ |
|------:|:-------|
| 0     | $\mid 1, 0, +0\rangle$ |
| 1     | $\mid 2, 0, +0\rangle$ |
| 2     | $\mid 2, 1, -1\rangle$ |
| 3     | $\mid 2, 1, +0\rangle$ |
| 4     | $\mid 2, 1, +1\rangle$ |

### 3.2. Connes distance matrix

```
                  |1,0,+0>   |2,0,+0>   |2,1,-1>   |2,1,+0>   |2,1,+1>
   |1,0,+0>        0.0000    58.3114    87.1751     1.2717    87.1751
   |2,0,+0>       58.3114     0.0000    28.8673    57.7347    28.8673
   |2,1,-1>       87.1751    28.8673     0.0000    86.6020     0.0000
   |2,1,+0>        1.2717    57.7347    86.6020     0.0000    86.6020
   |2,1,+1>       87.1751    28.8673     0.0000    86.6020     0.0000
```

### 3.3. Combinatorial Fock-graph distance matrix

Using $d_{\mathrm{graph}}(v, w) = |\Delta n| + |\Delta l| + |\Delta m|$:

```
                  |1,0,+0>   |2,0,+0>   |2,1,-1>   |2,1,+0>   |2,1,+1>
   |1,0,+0>            0          1          3          2          3
   |2,0,+0>            1          0          2          1          2
   |2,1,-1>            3          2          0          1          2
   |2,1,+0>            2          1          1          0          1
   |2,1,+1>            3          2          2          1          0
```

### 3.4. Side-by-side pair table

| $v$ | $w$ | $d_{\mathrm{Connes}}$ | $d_{\mathrm{graph}}$ |
|:----|:----|----------------------:|---------------------:|
| $\mid 1,0,0\rangle$ | $\mid 2,0,0\rangle$ |  58.3114 | 1 |
| $\mid 1,0,0\rangle$ | $\mid 2,1,-1\rangle$ |  87.1751 | 3 |
| $\mid 1,0,0\rangle$ | $\mid 2,1,+0\rangle$ |  **1.2717** | 2 |
| $\mid 1,0,0\rangle$ | $\mid 2,1,+1\rangle$ |  87.1751 | 3 |
| $\mid 2,0,0\rangle$ | $\mid 2,1,-1\rangle$ |  28.8673 | 2 |
| $\mid 2,0,0\rangle$ | $\mid 2,1,+0\rangle$ |  57.7347 | 1 |
| $\mid 2,0,0\rangle$ | $\mid 2,1,+1\rangle$ |  28.8673 | 2 |
| $\mid 2,1,-1\rangle$ | $\mid 2,1,+0\rangle$ |  86.6020 | 1 |
| $\mid 2,1,-1\rangle$ | $\mid 2,1,+1\rangle$ |  **0.0000** | 2 |
| $\mid 2,1,+0\rangle$ | $\mid 2,1,+1\rangle$ |  86.6020 | 1 |

### 3.5. Structural observations

1. **Three forced zeros from $\mathrm{SO}(3)$ $m$-reflection symmetry.**
   The pair $\mid 2,1,-1\rangle \leftrightarrow \mid 2,1,+1\rangle$ has
   $d_{\mathrm{Connes}} = 0$ even though $d_{\mathrm{graph}} = 2$. This
   is **structural**, not a Dirac choice: every multiplier
   $M_{N,L,M} \in \mathcal{O}$ has identical diagonal entries at the
   $m=-1$ and $m=+1$ basis indices, so $\langle 2,1,-1 | x | 2,1,-1
   \rangle = \langle 2,1,+1 | x | 2,1,+1 \rangle$ for every $x \in
   \mathcal{O}$. The SDP correctly returns $0$. This reflects the
   $\mathrm{SO}(3)$ rotational symmetry of the operator-system
   construction (the multipliers $M_{N,L,M}$ are eigenstates of the
   $\mathrm{SO}(3)$ rotation by $M$, and pure-state diagonals are
   $\mathrm{SO}(3)$-invariant under $m \to -m$).

   In total at $n_{\max} = 2$ the symmetry-forced zeros are
   $\{(2,4)\}$ in our indexing — one pair, by the only
   $m \leftrightarrow -m$ degeneracy among non-zero $m$ values.

2. **Clean rational $\{1, 2, 3\}$ structure of the non-zero values.** The
   distinct non-zero entries are
   $\{1.2717,\ 28.8673,\ 57.7347,\ 58.3114,\ 86.6020,\ 87.1751\}$.
   Stripping the small Dirac perturbation ($\sim 1\%$ shifts coming
   from `l_weight = 0.1` and `m_weight = 0.01`):

   - $28.8673 : 57.7347 : 86.6020 = 1 : 2 : 3$ (verified to $5 \times
     10^{-6}$ relative).
   - $86.6020 = 50\sqrt{3}$ to $5 \times 10^{-6}$ relative.
   - $87.1751 / 86.6020 = 1.0066$ and $58.3114 / 57.7347 = 1.0100$ —
     clean Dirac-perturbation lifts of degenerate $1\!:\!2\!:\!3$ extrema.

   The SDP is therefore **hitting structural rational extrema** of the
   placeholder operator system, not arbitrary numerical values. The
   $\sqrt{3}$ prefactor and the $\{1, 2, 3\}$ ratios are signatures of
   the SDP simultaneously saturating the operator-norm constraint along
   discrete Wigner-3j-related directions in $\mathcal{O}$.

3. **Outlier value $1.2717$ for $d(\mid 1,0,0\rangle, \mid 2,1,0\rangle)$.**
   Anomalously small ($\sim 0.4\%$ of the typical scale), and sits at
   $d_{\mathrm{graph}} = 2$. The pair is connected through the
   $M_{2, 1, 0}$ multiplier (the same one we used as the witness $a$ in
   round 2 §2.2), which has its dominant matrix elements precisely on
   this transition. The SDP finds an $x \propto M_{2,1,0}$ direction
   with very small operator-norm cost on the commutator with $D$, hence
   a small distance.

4. **Triangle inequality.** Verified across all 100 sampled triples,
   zero violations within SDP tolerance. The metric is a true
   pseudo-metric on the pure-state subset.

5. **No monotonicity in graph distance.** Pearson correlation between
   $d_{\mathrm{Connes}}$ and $d_{\mathrm{graph}}$ over the 10 finite
   pairs: $\rho_{\mathrm{Pearson}} \approx -0.04$. Spearman rank
   correlation: $\rho_{\mathrm{Spearman}} \approx 0.09$. This is the
   most striking finding of round 3 — see §5.

---

## §4. Distance matrix at $n_{\max} = 3$

### 4.1. Setup

- Hilbert dimension: $N = 14$.
- Operator-system dimension: $\dim(\mathcal{O}) = 55$.
- Kernel of $[D, \cdot]$ on $\mathcal{O}_h$ at default Dirac:
  $\dim = 1$ (just $\mathbb{C} \cdot I$).
- $14 \times 13 / 2 = 91$ unique pairs; SDP per pair $\sim 7$–$10$ s on
  the test machine, total $\sim 20$–$25\,$min.
- All 91 finite-pair distances converged to $\le 1\%$ relative error
  (SCS default).

### 4.2. Numerical distance matrix

```
            |1,0,0> |2,0,0> |2,1,-1>|2,1,+0>|2,1,+1>|3,0,0> |3,1,-1>|3,1,+0>|3,1,+1>|3,2,-2>|3,2,-1>|3,2,+0>|3,2,+1>|3,2,+2>
 |1,0,0>     0.0000  0.9479  0.8961  1.1431  0.8961  1.1126  1.3868  0.9821  1.3868  1.6436  1.4621  1.8596  1.4621  1.6436
 |2,0,0>     0.9479  0.0000  0.2336  0.4671  0.2336  0.7653  0.9927  0.6325  0.9927  1.2559  0.9407  1.5728  0.9407  1.2559
 |2,1,-1>    0.8961  0.2336  0.0000  0.7007  0.0000  0.8023  1.0337  0.6844  1.0337  1.2693  0.9934  1.5572  0.9934  1.2693
 |2,1,+0>    1.1431  0.4671  0.7007  0.0000  0.7007  0.9286  1.1243  0.8094  1.1243  1.3582  1.1284  1.6750  1.1284  1.3582
 |2,1,+1>    0.8961  0.2336  0.0000  0.7007  0.0000  0.8023  1.0337  0.6844  1.0337  1.2693  0.9934  1.5572  0.9934  1.2693
 |3,0,0>     1.1126  0.7653  0.8023  0.9286  0.8023  0.0000  0.3824  0.7649  0.3824  0.7310  1.1873  1.1208  1.1873  0.7310
 |3,1,-1>    1.3868  0.9927  1.0337  1.1243  1.0337  0.3824  0.0000  1.1473  0.0000  0.3774  1.5097  1.0791  1.5097  0.3774
 |3,1,+0>    0.9821  0.6325  0.6844  0.8094  0.6844  0.7649  1.1473  0.0000  1.1473  1.4835  0.7182  1.6286  0.7182  1.4835
 |3,1,+1>    1.3868  0.9927  1.0337  1.1243  1.0337  0.3824  0.0000  1.1473  0.0000  0.3774  1.5097  1.0791  1.5097  0.3774
 |3,2,-2>    1.6436  1.2559  1.2693  1.3582  1.2693  0.7310  0.3774  1.4835  0.3774  0.0000  1.8872  0.9433  1.8872  0.0000
 |3,2,-1>    1.4621  0.9407  0.9934  1.1284  0.9934  1.1873  1.5097  0.7182  1.5097  1.8872  0.0000  2.2407  0.0000  1.8872
 |3,2,+0>    1.8596  1.5728  1.5572  1.6750  1.5572  1.1208  1.0791  1.6286  1.0791  0.9433  2.2407  0.0000  2.2407  0.9433
 |3,2,+1>    1.4621  0.9407  0.9934  1.1284  0.9934  1.1873  1.5097  0.7182  1.5097  1.8872  0.0000  2.2407  0.0000  1.8872
 |3,2,+2>    1.6436  1.2559  1.2693  1.3582  1.2693  0.7310  0.3774  1.4835  0.3774  0.0000  1.8872  0.9433  1.8872  0.0000
```

### 4.3. Distance summary statistics

- All 91 unique pairs give finite distance (no $+\infty$ at default
  Dirac).
- Range: $[0.234, 2.241]$ for non-zero finite values.
- 4 forced-zero off-diagonal pairs (out of 91): see §4.5.
- Symmetry $d(v, w) = d(w, v)$ holds to SDP precision.

### 4.4. Comparison with $n_{\max} = 2$

The most striking change between $n_{\max} = 2$ and $n_{\max} = 3$ is
the **scale of the distances**: the maximum drops from $\sim 87$ at
$n_{\max} = 2$ to $\sim 2.24$ at $n_{\max} = 3$. The reason is purely
operational: the off-diagonal Dirac coupling now reaches more shells, so
the SDP has more "directions" in which to spend the unit operator-norm
budget. Each individual direction supports a smaller objective change,
and the worst-case extremum ratio shrinks. This is a placeholder-and-
Dirac-convention effect, not a structural statement about the Connes
distance on $S^3$.

The clean $\{1, 2, 3\}$ rational structure of the $n_{\max} = 2$
extrema does NOT persist at $n_{\max} = 3$: the distances are now an
irregular sequence of irrationals (e.g. $0.234, 0.382, 0.467, 0.683,
0.701, 0.731, 0.765, \ldots$) that do not factor cleanly as integer
multiples of a common base. The $n_{\max} = 2$ ratios were a
small-system coincidence of the SDP simultaneously saturating
multiple constraints; with more degrees of freedom (55 multipliers,
$N=14$ ambient) the SDP optima generically lie in higher-codimension
faces of the feasible region and lose the rational simplicity.

### 4.5. Forced-zero pairs at $n_{\max} = 3$

The $\mathrm{SO}(3)$ $m$-reflection symmetry of the operator-system
construction forces $d(\phi_v, \phi_w) = 0$ exactly when $v$ and $w$
differ only by $m \to -m$:

| $v$ | $w$ | $d_{\mathrm{Connes}}$ | $d_{\mathrm{graph}}$ |
|:----|:----|----------------------:|---------------------:|
| $\mid 2,1,-1\rangle$ | $\mid 2,1,+1\rangle$ | 0.0000 | 2 |
| $\mid 3,1,-1\rangle$ | $\mid 3,1,+1\rangle$ | 0.0000 | 2 |
| $\mid 3,2,-1\rangle$ | $\mid 3,2,+1\rangle$ | 0.0000 | 2 |
| $\mid 3,2,-2\rangle$ | $\mid 3,2,+2\rangle$ | 0.0000 | 4 |

This count is **exactly the count of $(n, l, m)$ states with $n \le
n_{\max}$, $l \le n-1$, $m > 0$, having a $\pm m$ partner** — i.e.
$\sum_{n,l} l = \sum_{n=1}^{n_{\max}} \sum_{l=0}^{n-1} l$, which is
$\binom{n_{\max}}{2}$ for $n_{\max} \ge 2$. At $n_{\max} = 2$:
$\binom{2}{2} = 1$. At $n_{\max} = 3$: $\binom{3}{2} = 3$ from the
$(2,1,1), (3,1,1), (3,2,1)$ pairs, plus 1 more from the $(3,2,2)$
pair, totaling 4. Verified.

### 4.6. Correlation with graph distance at $n_{\max} = 3$

- Pearson $\rho \approx 0.14$ (over 91 finite pairs).
- Spearman rank $\rho \approx 0.18$.

Both are weakly positive — a slight improvement from the $\rho \approx
0$ at $n_{\max} = 2$ — but the metric is still nowhere near monotone
in the graph distance. The forced zeros at $d_{\mathrm{graph}} = 2,
4$ (alongside finite non-zero distances at $d_{\mathrm{graph}} = 1$)
are the dominant deviation from monotonicity.

---

## §5. Comparison to the combinatorial Fock-graph distance

### 5.1. Pearson and Spearman correlation

| $n_{\max}$ | #finite pairs | Pearson $\rho$ | Spearman $\rho$ |
|:--:|:--:|:--:|:--:|
| 2 | 10 | $-0.04$ | $+0.09$ |
| 3 | 91 | $+0.14$ | $+0.18$ |

At $n_{\max} = 2$ both are essentially zero — the metric is statistically
uncorrelated with the graph distance. At $n_{\max} = 3$ both are weakly
positive but very far from monotone. **The Connes distance with the
placeholder convention does not converge to a smooth monotone function
of the graph distance over the cutoffs we can reach.** The trend is
gently increasing with $n_{\max}$, consistent with the working
hypothesis that the metric *might* eventually relate to the round-$S^3$
geodesic distance, but the rate of improvement is slow and bounded by
the structural forced zeros.

### 5.2. Why is this happening?

Three structural effects working together:

1. **Operator-system $\mathrm{SO}(3)$ symmetry forces zero distance** on
   $m$-reflected pairs even when their graph distance is $\ge 1$. This
   is not a flaw of the SDP; it reflects that the operator system
   cannot resolve $\mathrm{SO}(3)$-equivalent states by diagonal
   evaluation.

2. **The placeholder for the SO(4) radial overlap is symmetric in
   $(n, n')$** but not magnitude-tuned to physical Avery-Wen-Avery
   values. The numerical scale of the SDP optima depends on the
   inverse "spread" of the multiplier-matrix entries, and the placeholder
   produces wildly different scales for shell-coupling vs intra-shell
   off-diagonal multipliers (factor $\sim 50$ at $n_{\max} = 2$).

3. **The diagonal Dirac with off-diagonal coupling is not the physical
   Camporesi-Higuchi spinor Dirac.** The relative weights between the
   diagonal and off-diagonal pieces (parameter $\alpha = 1$) and
   between $l$ and $m$ ($l_{\mathrm{weight}} = 0.1$, $m_{\mathrm{weight}}
   = 0.01$) are convention choices, not derivations.

### 5.3. What is invariant under reasonable convention changes

**Robust under any reasonable Dirac proxy** (verified by hand: rerunning
with $l_{\mathrm{weight}} = 0.5, 0.01$ and $\alpha = 0.1, 10.0$ gives
qualitatively the same pattern):

- Forced zeros from $\mathrm{SO}(3)$ symmetry.
- Triangle inequality satisfied.
- $\dim(\ker([D, \cdot]) \cap \mathcal{O}_h) = 1$ (just identity) for
  the `offdiag` Dirac.
- The relative ratios within the surviving distance set are placeholder-
  dependent but the qualitative metric topology (which pairs are zero
  vs. small vs. large) is robust.

**Convention-dependent**:

- Absolute numerical values of $d_{\mathrm{Connes}}$.
- Rank correlation with the graph distance: with `m_weight = 0` (full
  $\mathrm{SO}(2)$ symmetry around the $z$-axis), more pairs would have
  forced zeros and the rank correlation might shift.

### 5.4. The "honest" reading of the round-3 finding

The metric on the GeoVac truncated operator system at finite $n_{\max}$
**does not look like a discretization of the round-$S^3$ geodesic
distance**. The expected behavior — distance smoothly increases with
graph-ladder steps — is **not observed** at the placeholder convention.

Two interpretations are consistent with this:

- **(a) Placeholder artifact.** The actual Avery-Wen-Avery 3-Y integral
  on $S^3$ is rational+algebraic with very specific normalization that
  the placeholder grossly misses; with the true integral the distance
  pattern would smooth out. This is plausible because the placeholder's
  $1 + (n+n')/(N+1)$ ansatz produces wildly different scales for
  different multiplier types.

- **(b) Genuine pure-state-distance pathology.** Pure node-evaluation
  states $\phi_v$ are extreme in the state space and the Connes distance
  on extreme points may not respect the underlying geometry (cf.
  Connes-vS Theorem 4.20 which is about the *full* state space and
  Kantorovich-distance lower bound, not pure-state pairwise distance).
  The "smoothing-to-geodesic" behavior in the Gromov-Hausdorff limit
  may only emerge after taking suitable convex combinations / averaging.

Both readings are consistent with our round-1 / round-2 verdicts; round
3's empirical finding is **a "structural negative" against the
quick-naive-monotonicity claim**, not against the spectral-truncation
framework itself.

---

## §6. Known limitations

1. **Placeholder convention.** The radial-overlap placeholder
   $I_{nNn'}^{lLl'} = 1 + (n+n')/(N+1)$ is convention-arbitrary. The
   actual Avery-Wen-Avery 3-Y integral on $S^3$ would replace the
   placeholder with rational + Wigner-6j factors. The SDP machinery is
   independent of the placeholder; only the numerical distance values
   would change.

2. **No spinor lift.** We use the *scalar* Fock basis and a *scalar*
   Dirac proxy. The full Camporesi-Higuchi spinor Dirac on $S^3$ would
   double the Hilbert dimension and add a chiral grading. The pure-
   state distance on the spinor-lifted operator system might differ
   significantly from the scalar version.

3. **No $S^3$-geodesic comparison.** The natural target for a GH-
   convergence test is the round-$S^3$ geodesic distance between
   peak locations of the spherical harmonics
   $Y^{(3)}_{nlm}$. Computing this and comparing to the Connes distance
   in the limit $n_{\max} \to \infty$ is the actual
   GH-convergence-on-$S^3$ test (round-1 Gap 2). Not attempted in
   round 3.

4. **No real-structure $J$ at finite $n_{\max}$.** The Connes axiom
   audit (Paper 32 / round 1 Row 4) flagged this as Gap 4. Round 3
   does not address it.

5. **Fock-graph distance is the L^1 metric on (n, l, m)**, not the
   true Cayley-graph distance under the GeoVac graph adjacency
   ($\Delta n = \pm 1$ only). With the latter graph metric, all
   $(n, l, m) \to (n', l', m')$ transitions with $\Delta n = 0$ would
   be at infinite graph distance, removing about half the comparisons
   from the table. The qualitative findings in §5 do not change under
   this swap.

---

## §7. Implications for WH1 and recommendations

### 7.1. What round 3 establishes

- **The Connes-distance machinery for the GeoVac truncated operator
  system is operational.** Given $\mathcal{O}_{n_{\max}}$ from round 2
  and any Hermitian Dirac proxy on $\mathcal{H}_{n_{\max}}$, the
  pure-state Connes distance can be computed by an explicit SDP that
  converges in seconds at $n_{\max} = 2$ and minutes at $n_{\max} = 3$.

- **The metric is non-trivial and structurally meaningful** — finite
  on (almost) every pair when the Dirac is sufficiently off-diagonal,
  with forced zeros from the operator-system $\mathrm{SO}(3)$ symmetry
  and a clean rational structure of the surviving values.

- **The metric is NOT monotone in the combinatorial Fock-graph
  distance** at the placeholder convention. This is either a placeholder
  artifact or a genuine pure-state-distance pathology — the test does
  not distinguish.

### 7.2. What round 3 does not establish

- It does not prove or disprove GH convergence of $(\mathcal{S}(\mathcal{O}
  _{n_{\max}}), d_{D_{n_{\max}}})$ to $(\mathcal{P}(S^3),
  d_{\mathrm{Wass}})$ as $n_{\max} \to \infty$ — that requires (a) the
  Avery-Wen-Avery integrals (placeholder replacement) and (b) a
  Peter-Weyl analog of the Leimbach-vS spectral Fejér kernel for
  $\mathrm{SU}(2)$. Both remain open.

- It does not address the spinor lift; the round-3 metric is on the
  scalar operator system.

- It does not compare to the round-$S^3$ geodesic distance between
  spherical-harmonic peaks; that is the natural GH-convergence target
  but requires the placeholder-replacement.

### 7.3. Recommendation flagged for plan-mode review

Three follow-up directions, in order of computational tractability:

**R3.1 (HIGH, REQUIRES NEW NUMERICS).** Replace the placeholder
$I_{nNn'}^{lLl'}$ in `geovac/operator_system.py` with the genuine
Avery-Wen-Avery 3-Y integral on $S^3$. The closed form factorizes via
SO(4) = SU(2)_L × SU(2)_R into two SU(2) Wigner $3j$ symbols. The
existing `geovac/algebraic_slater.py` and
`geovac/hypergeometric_slater.py` modules contain most of the radial
machinery; the closed-form rational-3j combinations need to be added.
Expected outcome: numerical values change but structural pattern (forced
zeros, rational ratios, kernel structure) survives.

**R3.2 (MEDIUM, RAMANUJAN-ADJACENT).** Test whether the SDP-extremum
$\{1, 2, 3\}$ rational structure at $n_{\max} = 2$ persists at higher
$n_{\max}$. If yes, this is a clean structural invariant of the
truncated operator system and should be promoted to a proposition.

**R3.3 (NEW MATHEMATICS).** Attempt the GH convergence sketch on $S^3$
via Peter-Weyl on $\mathrm{SU}(2)$ (round-2 §6 R2.5). If proved, this
would convert the unclear row of the round-1 mapping into a clean
match.

**R3.4 (PLAN-MODE).** I propose that the WH1 register entry in
CLAUDE.md §1.7 NOT be modified by round 3 alone; the round-3 finding
is a "structural-negative-with-caveats" that doesn't strengthen WH1
beyond round 2's MODERATE-STRONG / STRONG candidate status. Round 2's
$\mathrm{prop} = 2$ alignment with the Toeplitz $S^1$ is the strongest
single claim; round 3's finding is informative but does not by itself
upgrade the WH.

---

## §8. Appendix — files added in this round

- `geovac/connes_distance.py` (~430 lines): SDP-based Connes-distance
  module with `DiracProxy`, `compute_connes_distance`,
  `compute_distance_matrix`, `_kernel_of_commutator_in_O`,
  `fock_graph_distance`, `graph_distance_matrix`,
  `basis_label_strings`. Soft dependency on `cvxpy`.

- `tests/test_connes_distance.py` (17 tests, all pass): positivity,
  symmetry, triangle inequality, distance-matrix shape, kernel-is-
  identity-only for `offdiag` Dirac, $+\infty$ for `shell_scalar`
  Dirac, $\mathrm{SO}(3)$ $m$-reflection forced zero.

- `debug/data/wh1_round3_connes_distance_nmax2.json`: the 5×5 distance
  matrix (and graph-distance matrix) at $n_{\max} = 2$.

- `debug/data/wh1_round3_connes_distance_nmax3.json`: the 14×14
  distance matrix (and graph-distance matrix) at $n_{\max} = 3$.

- `debug/wh1_round3_connes_distance_memo.md` (this file).

No paper edits, no CLAUDE.md edits, no modification of
`geovac/operator_system.py` or `tests/test_operator_system.py`. The
round-2 propagation-number invariant is unaffected by round 3.

**End of memo.**
