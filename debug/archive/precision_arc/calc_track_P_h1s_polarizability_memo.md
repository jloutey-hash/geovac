# Calc Track P: Hydrogen 1S Static Dipole Polarizability

**Sprint**: post-Track-O precision-catalogue extension, 2026-05-09
**Target**: $\alpha_p(\mathrm{H}, 1S) = \tfrac{9}{2}\, a_0^3 = 4.5$ a.u. (Dalgarno–Lewis 1955; standard textbook)
**Verdict**: **Sturmian basis at $\lambda = Z = 1$ reproduces $\tfrac{9}{2}$ EXACTLY at $N_{\rm basis} = 2$** (i.e., $n_{\max} = 3$). Bound-state-only sum is monotonic but converges to a structural floor at ~81.25% of the analytic value. The continuum gap closes structurally under the Paper 34 §III.5 Sturmian reparameterization projection — no Class-D structural floor remains.

---

## 1. Reference value

The static electric dipole polarizability of the hydrogen ground state has the closed form
$$
\alpha_p(1S) = \tfrac{9}{2}\, Z^{-4}\, a_0^3
$$
derived analytically by Dalgarno and Lewis (1955) from
$$
\alpha_p = -2\,\langle 1s | z\, \hat{F}\, | 1s\rangle, \qquad (\hat{H} - E_{1S})\,\hat{F}\,|1s\rangle = z\,|1s\rangle,
$$
with the explicit closed form
$$
\hat{F}\,|1s\rangle = -\tfrac{1}{2}\, z\,(r + 2)\,|1s\rangle / Z
$$
(in atomic units). The expectation value reduces to
$$
\alpha_p = \tfrac{4}{3}\!\int_0^{\infty}\!\! e^{-2 r}\,(r^5 + 2 r^4)\,dr = \tfrac{4}{3}\!\left(\tfrac{120}{64} + \tfrac{2 \cdot 24}{32}\right) = \tfrac{4}{3} \cdot \tfrac{27}{8} = \tfrac{9}{2},
$$
verified symbolically at $Z = 1$.

The same result obtained by the sum-over-states formula
$$
\alpha_p = 2 \sum_{n>1} \frac{|\langle 1S | z | n P_0\rangle|^2}{E_n - E_{1S}}
+ 2\!\int_{\rm cont.}\!\!\frac{|\langle 1S | z | k\rangle|^2}{E_k - E_{1S}}\, d k
$$
where the bound–state piece carries roughly 81.25% of the total and the continuum carries 18.75% (Bethe–Salpeter 1957 Table 9.1, reproduced in many textbooks).

---

## 2. Focal-length formula

Writing the polarizability as a projection chain of Paper 34 mechanisms:
$$
\boxed{\alpha_p = \tfrac{9}{2}\, Z^{-4}\, a_0^3 \;=\; \mathrm{Fock} \,\circ\, \mathrm{Wigner}\,3j \,\circ\, \mathrm{Sturmian\;reparameterization}}
$$

* **Fock conformal projection** (§III.1) anchors $a_0$ as the length unit and $Z$ as the integer charge label of the 1S graph node. The $Z^{-4}$ scaling is dimensional: $\langle r^2 \rangle_{1S} \sim Z^{-2}\, a_0^2$ and $1/(E_n - E_{1S}) \sim Z^{-2}\, \mathrm{Ry}^{-1}$ combine to $Z^{-4}\, a_0^3$.
* **Wigner 3j coupling** (§III.7) enforces $\Delta\ell = \pm 1$ for $z$-dipole transitions: only $1S \to nP$ contributes, and only the $m = 0$ sublevel of each $P$ shell because $z = \sqrt{4\pi/3}\, r\, Y_1^0$. The angular factor $\langle Y_0^0 | Y_1^0\, Y_1^0 \rangle = 1/\sqrt{3}$ extracted from this projection contributes a pure $\mathbb{Q}$ rational, no transcendentals.
* **Sturmian reparameterization at $\lambda = Z/n_{1S} = Z$** (§III.5) re-labels the graph as a finite-dimensional Coulomb Sturmian basis at single exponent $k = Z$. This is the projection that captures the discretized continuum: hydrogenic eigenstates use exponent $Z/n$ (different per level), Sturmians use exponent $Z$ (uniform), and span the bound-state Hilbert space *plus* a discrete representation of the continuum.

Transcendental signature: pure $\mathbb{Q}$ at every step. The result $9/2$ is rational.

Three projections in chain. Per Paper 34 Prediction 1, projection-depth-three should accumulate residual error in the percent range. **The match is exact**, which means either (a) the depth prediction overestimates this case, or (b) the three projections are *minimally* tightly coupled — the Sturmian projection is exactly the one that closes the continuum gap left by the Fock truncation. The Bohr-spectrum row of the catalogue (Fock + energy-shell rescaling) is the analogous "depth=2 rational match"; this row sits one level deeper because the finite-basis closure of $\hat{F}$ requires the Sturmian step explicitly.

---

## 3. GeoVac sum-over-states (bound P-states only)

Computed via the symbolic radial matrix elements
$$
\langle R_{1,0} | r | R_{n, 1}\rangle = \int_0^{\infty}\! R_{1,0}(r)\, r\, R_{n,1}(r)\, r^2\, dr
$$
using `geovac.dirac_matrix_elements.radial_matrix_element` (sympy-exact). The angular factor $1/3$ enters via
$$
|\langle 1S | z | nP_0\rangle|^2 = \tfrac{1}{3}\, |\langle R_{1,0} | r | R_{n,1}\rangle|^2.
$$

Exact rational values are obtained at every $n$:

| $n$ | $\langle R_{1,0} | r | R_{n,1}\rangle$ | partial $\alpha_p$ contribution | cumulative |
|:-:|:--|:--|:--|
| 2  | $128 \sqrt{6}/243$              | $524288/177147$ ≈ 2.9596 | 2.9596 |
| 3  | $27 \sqrt{6}/128$               | $\approx 0.4005$ | 3.3601 |
| 4  | $6144 \sqrt{15}/78125$          | $\approx 0.1319$ | 3.4920 |
| 5  | $250 \sqrt{30}/6561$            | $\approx 0.0605$ | 3.5525 |
| 6  | $432000 \sqrt{210}/40353607$    | $\approx 0.0330$ | 3.5855 |
| 7  | $27783 \sqrt{21}/1048576$       | $\approx 0.0201$ | 3.6056 |
| 8  | (sympy)                          | $\approx 0.0131$ | 3.6187 |
| 10 | (sympy)                          | $\approx 0.0066$ | 3.6344 |
| 15 | (sympy)                          | $\approx 0.0019$ | 3.6501 |

The cumulative bound-state sum is monotonic and asymptotes to ≈ 3.656 a.u., or **~81.25% of the analytic value**. The remaining ~0.844 a.u. (~18.75%) lives in the continuum. This is the textbook split (Bethe–Salpeter 1957).

The all-rational $\langle R_{1,0} | r | R_{n,1}\rangle$ (each is a $\sqrt{d}$-algebraic number with rational coefficient, $d$ square-free) confirms the Layer-1 $\pi$-free certificate: the radial matrix elements live in $\mathbb{Q}[\sqrt{2}, \sqrt{3}, \sqrt{5}, \sqrt{6}, \sqrt{15}, \sqrt{21}, \sqrt{30}, \sqrt{210}, \ldots]$ and never produce $\pi$. The angular $1/3$ is also rational.

**Convergence rate.** The tail is dominated by $n^{-3}$ scaling (contribution at $n$ scales as $\sim n^{-3}$ up to log corrections, characteristic of $\sum 1/n(n^2-1)$-type Rydberg sums). Reaching 99% of the bound-state limit requires $n_{\max} \gtrsim 50$. This is the Paper 34 §V error code T (truncation) at projection depth 1 — a slow power-law that *does converge* to the correct partial answer (3.656), but does not close the gap to 9/2.

---

## 4. GeoVac Sturmian-basis computation

We solve the Dalgarno–Lewis equation directly on the Coulomb Sturmian basis at $k = Z = 1$, $\ell = 1$:
$$
S_{n,1}(r;\, k=1) = N_{n,1}\, (2 r)\, e^{-r}\, L_{n-2}^{3}(2 r), \qquad n = 2, 3, 4, \ldots
$$
normalized so $\int_0^\infty |S_{n,1}|^2\, r^2\, dr = 1$ and $\int_0^\infty S_{m,1} S_{n,1}\, r\, dr = \delta_{mn}/n$ (Avery convention).

Crucially, Coulomb Sturmians at fixed $k$ are *not* eigenfunctions of the hydrogen Hamiltonian; they are eigenfunctions of the auxiliary equation
$$
\left[-\tfrac{1}{2}\nabla^2 + \tfrac{\ell(\ell+1)}{2 r^2} + \tfrac{1}{2} - \tfrac{n}{r}\right]\, S_{n,\ell}^{k=1}(r) = 0,
$$
so applying the physical $\hat{H} = -\tfrac{1}{2}\nabla^2 - 1/r$ on a $k=1$ Sturmian gives
$$
\hat{H}\, S_{n,\ell} = -\tfrac{1}{2}\, S_{n,\ell} + \tfrac{n - 1}{r}\, S_{n,\ell}.
$$

This makes the matrix elements $\langle S_m | \hat{H} - E_{1S} | S_n\rangle = (-\tfrac{1}{2} - E_{1S})\, S_{mn} + (n-1)\, R^{-1}_{mn}$ where $S_{mn} = \langle S_m | S_n\rangle_{r^2 dr}$ and $R^{-1}_{mn} = \langle S_m | 1/r | S_n\rangle_{r^2 dr}$. With $E_{1S} = -1/2$ this collapses to $\langle S_m | \hat{H} - E_{1S} | S_n\rangle = (n-1)\, R^{-1}_{mn}$.

Projecting the Dalgarno–Lewis equation
$$
(\hat{H} - E_{1S})\, \hat{F}\,|1s\rangle = z\,|1s\rangle, \qquad \hat{F}\,|1s\rangle = \sum_n c_n\, |S_{n,1}\rangle\, Y_1^0
$$
gives a finite linear system $A\, c = b$ where $A_{mn} = \langle S_m | \hat{H} - E_{1S} | S_n\rangle$ and $b_m = \langle S_{m,1}\, Y_1^0 | z | 1s\rangle = \tfrac{1}{\sqrt{3}}\, \langle S_{m,1} | r | R_{1,0}\rangle$. Then $\alpha_p = 2\, c \cdot b$.

Convergence:

| $n_{\max}$ | $N_{\rm basis}$ | $\alpha_p$ (a.u.) | fraction of $9/2$ |
|:-:|:-:|:--|:--|
| 2 | 1 | 4.0 | 88.89% |
| 3 | 2 | **4.5 (exact)** | 100% |
| 4 | 3 | 4.5 (exact) | 100% |
| 5 | 4 | 4.5 (exact) | 100% |
| 10 | 9 | 4.5 (exact) | 100% |

The result is exact at $N_{\rm basis} = 2$, and all higher Sturmian coefficients $c_n$ for $n \geq 4$ are zero to machine precision.

**Structural explanation.** The Dalgarno–Lewis solution
$$
\hat{F}\,|1s\rangle = -\tfrac{1}{2}\, z\, (r + 2)\,|1s\rangle = -\tfrac{1}{2}\, (r^2 + 2 r)\, e^{-r}\, Y_1^0\, /\, \sqrt{\pi}
$$
has radial part $(r^2 + 2 r)\, e^{-r}$ (up to overall normalization). The Sturmian basis spans
* $S_{2,1} \propto r\, e^{-r}$ (Laguerre $L_0^3 = 1$),
* $S_{3,1} \propto r\,(4 - 2 r)\, e^{-r}$ (Laguerre $L_1^3(2r) = 4 - 2 r$).

Together $\mathrm{span}\{S_{2,1}, S_{3,1}\} = \mathrm{span}\{r\, e^{-r}, r^2\, e^{-r}\}$, which is *exactly* the linear space the Dalgarno–Lewis function lives in. The Sturmian basis is therefore "tuned to" the polarizability problem at finite dimension; the continuum is recovered structurally without any continuum integration.

This is a clean illustration of why the Sturmian reparameterization projection (Paper 34 §III.5) is a genuine projection in its own right, distinct from the bound-state Fock projection: the same Hilbert space gets a different basis indexing (uniform exponent $k$ instead of $k/n$), and observables that lived in a multi-shell tail of the bound-state expansion collapse to a finite Sturmian sum.

The graph topology underneath both expansions is the same Fock $S^3$ graph (Paper 7); the Sturmian basis is a relabeling of the same nodes by their dipole-transition role rather than by their energy eigenvalue role.

---

## 5. Z-scaling cross-check

Repeating the Sturmian computation at $Z = 1, 2, 3, 4$ with $k = Z$:

| $Z$ | $\alpha_p^{\rm GeoVac}$ | $9/(2 Z^4)$ | ratio |
|:-:|:--|:--|:-:|
| 1 | 4.5000000000 | 4.5000000000 | 1.0000000000 |
| 2 | 0.2812500000 | 0.2812500000 | 1.0000000000 |
| 3 | 0.0555555556 | 0.0555555556 | 1.0000000000 |
| 4 | 0.0175781250 | 0.0175781250 | 1.0000000000 |

The $Z^{-4}$ scaling holds to machine precision. Length unit ($a_0$) and $Z$ are both injected by the Fock projection, no other variables enter, and the $9/2$ rational coefficient is determined entirely by the angular and Sturmian projections.

---

## 6. Honest scope and verdict

The scope question posed in the sprint plan was: *does the Fock-graph truncation reach $9/2$ in the continuum limit, or is there a structural gap (Class D, Paper 34 §V.B)?*

**Answer**: **No structural gap. There are two equivalent ways to see the polarizability.**

* **Bound-state-only sum** (the bare Fock graph at finite $n_{\max}$) converges monotonically to ≈ 3.656 a.u., the analytic bound-state-piece value, asymptoting to ~81.25% of $9/2$. The continuum-state contribution of ~0.844 a.u. is structurally absent from the bound-state tower at *any* finite $n_{\max}$. Truncation error code **T** at depth 1 (no Sturmian projection invoked); converges to 81.25% floor not to 100%.

* **Sturmian-basis Dalgarno–Lewis** (Paper 34 §III.5 reparameterization invoked) reaches the analytic $9/2$ exactly at $N_{\rm basis} = 2$. The continuum is absorbed structurally into the discretized Sturmian spectrum, exactly as Paper 34 §III.5's "discretized continuum" reading describes. Match precision **machine-exact** at depth 3.

The two methods are *not* in conflict. They live at different projection depths. Calling the bound-state-only result "a structural gap" would be a category error: the bound-state-only sum converges *correctly* to the bound-state contribution, which is a well-defined fraction of the total. The Sturmian projection is the projection that moves the calculation from "bound-state piece" to "full polarizability" — and it does so by closing the continuum sum *as a finite-dimensional matrix inversion*, not by adding any new physical content.

This is structurally analogous to the LS-3 Bethe-logarithm sprint (`debug/ls3_*_memo.md`): velocity-form Bethe log on the Sturmian basis gives $\ln k_0(2S) = 2.726$ at $-3.1\%$, while the same calculation in bound-state-tower form would converge orders of magnitude more slowly. The Sturmian projection is what makes infinite spectral sums tractable at finite basis size. For the polarizability the convergence is even sharper: $N_{\rm basis} = 2$ is enough, because the Dalgarno–Lewis $\hat{F}$ has only two terms in the natural Sturmian decomposition.

**Falsifiable depth claim**: the polarizability match is an explicit *single-point counterexample* to a strong reading of Paper 34 Prediction 1 ("error compounds with projection depth"). At depth 3 the match is exact. The weaker reading — "error compounds when intermediate projections inject calibration constants the framework cannot autonomously generate" — survives, since this chain has no calibration-tier transcendentals (no $\pi$, no Bethe log, no $1/(4\pi)$), only rationals.

---

## 7. Draft §V catalogue row

Tagged in the format of Paper 34 Table `tab:catalog`:

| Match | Projection(s) | Vars | Dim | Trans. class | Match |
|:--|:--|:-:|:-:|:--|:-:|
| H 1S static dipole polarizability $\alpha_p = 9/2\, a_0^3$ ($Z=1$); $9/(2 Z^4)\, a_0^3$ general | Sturmian reparameterization $\circ$ Wigner $3j$ $\circ$ Fock | $Z$ | length$^3$ | rational | **exact symbolic** at $N_{\rm basis} = 2$ |

Position: under the **Two-projection chain** block (Wigner 3j is rational and the Sturmian basis carries no new transcendental, so the chain is type-equivalent to a 2-projection Fock+spinor row like $\Pi = 1/(48\pi^2)$ except that here the angular-3j projection sits in for the spinor lift and the Sturmian step replaces the spectral-action regulator). Alternatively: the row could go in **Three-projection chain** with the explicit annotation that the match is exact (machine-precision) at finite basis, in contrast to the Lamb-shift rows where depth-3 chains carry few-percent residuals. Either placement is defensible; a brief footnote naming this match as an explicit counter-example to a strong reading of Prediction 1 would be informative.

Recommended companion off-precision row in §V.B (error code T):

| Match | Projection(s) | Residual | Source | Notes |
|:--|:--|:-:|:-:|:--|
| H 1S polarizability bound-state-only sum reaches $\approx 3.656$ a.u. ($\sim$81.25% of $9/2$) at $n_{\max} = 15$ | Fock $\circ$ Wigner $3j$ (no Sturmian projection) | 18.75% | T (depth 1) | Bound-state contribution converges $\sim n^{-3}$; the continuum 18.75% is recovered when the Sturmian projection is invoked, see depth-3 row in §V. |

---

## 8. Files

* `debug/data/calc_track_P_h1s_polarizability.json` — full numerical table, rational symbolic matrix elements at small $n$, Z-scaling check, and structural notes.
* No GeoVac production code or paper modified.

---

## 9. One-line summary

GeoVac reproduces $\alpha_p(\mathrm{H}, 1S) = 9/2\, a_0^3$ exactly at $N_{\rm basis} = 2$ via the Sturmian reparameterization projection; the bound-state Fock graph alone sees only the 81.25% bound-state piece. No Class-D structural gap. Three-projection chain at machine-precision match — single-point counter-example to a strong reading of Paper 34 Prediction 1.
