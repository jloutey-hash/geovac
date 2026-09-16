# The Tao–McCurdy–Rescigno method, and what it says about Paper 12

**Canonical memo.** PI request 2026-09-14: *"Let's look into the Tao–McCurdy–Rescigno
method"* — follow-up to the 2026-09-13 accuracy-path scan
(`debug/lit_scan/sturmian_accuracy_paths_memo.md`), which surfaced TMR as the one
published chemical-accuracy result in GeoVac's own Level-2 coordinate system.

**Headline: Paper 12's 92.4% plateau is the σ-only ($m_1=m_2=0$) restriction, not the
electron–electron cusp.** Two independent routes agree, one external and verified at
source, one internal and already in the corpus. Paper 12's stated diagnosis
("the cusp requires non-analytic $r_{12}^{1/2}$, $r_{12}\ln r_{12}$ terms that no
polynomial prolate spheroidal basis can represent") is contradicted by both.

---

## 1. The method

L. Tao, C. W. McCurdy, T. N. Rescigno, "Grid-based methods for diatomic quantum
scattering problems. III. Double photoionization of molecular hydrogen in prolate
spheroidal coordinates," *Phys. Rev. A* **82**, 023423 (2010),
DOI 10.1103/PhysRevA.82.023423. **VERIFIED at source** — OSTI accepted manuscript
(`https://www.osti.gov/servlets/purl/1051651`), text-extracted with `pdftotext -layout`;
all quotations below are verbatim from that text. (APS, AIP, Wiley, ScienceDirect all
403 the fetcher; OSTI PURLs and arXiv work.)

**Coordinates.** Prolate spheroidal $(\xi,\eta,\varphi)$ — the same system as Papers 11
and 12.

**One-particle basis.** Product of a FEM/DVR function in $\xi$ with an *ordinary
spherical harmonic* $Y_{lm}$ in $(\cos^{-1}\eta, \varphi)$:
$\phi_{ilm}(\mathbf r) = \chi_i(\xi)\,Y_{lm}(\eta,\varphi)$. The two-electron expansion
is $\Psi = \sum C_{ijl_1l_2m}\,\phi_{il_1m}(\mathbf r_1)\,\phi_{jl_2,M-m}(\mathbf r_2)$
— **$m$ is summed**, with only the total $M$ fixed. For the ground state $M=0$, so
$m_2=-m_1$ and the $\pi^2$, $\delta^2$, … configurations are all present.

**Electron repulsion.** Not a quadrature and not a closed-form Green's function:
they use the *Poisson* route. The Green's function for Poisson's equation in prolate
spheroidal coordinates, expanded in spherical harmonics in $\eta$, gives a
one-dimensional ODE in $\xi$ per $(l,m)$ channel, solved by inverting the DVR
representation $[T_{lm}]^{-1}$ of that operator. Consequence: the two-electron integral
is **diagonal in the radial DVR indices** ($\delta_{ii'}\delta_{jj'}$ in their final
expression) — $N_\xi^2$ nonzeros instead of $N_\xi^4$.

**Angular couplings.** Pure Wigner $3j$ (Gaunt), done analytically, with the
$(\xi^2-\eta^2)$ Jacobian handled by the $Y_l^m Y_2^0$ recoupling identity. The sum
**terminates exactly**: $l_{\max} = \min(l_1+l_1'+2,\,l_2+l_2'+2)$, "the latter identity
following from the properties of the $3-j$ symbols." This is the same exact-termination
structure as Paper 15's split-region Legendre expansion.

**Parameters and result.** Grid in $\xi$: 8 real intervals of length 8.0 from 1.0 to
65.0, 14th-order DVR each, plus one complex ECS element (scattering only). Angular:
$l_{\max}=6$. Verbatim:

> "we retained terms up to lmax = 6. The ground-state H2 energy, obtained using only
> the first two real elements, was **-1.17442 hartree**, in excellent agreement with the
> accurate value of -1.17447 hartree results of Wolniewicz [26]. We note that
> calculations in spherical coordinates reported in ref. [9] using a single-center
> expansion with lmax = 7 give a target energy of **-1.16908 hartree**, which is 0.15 eV
> higher than the accurate value."

So: **0.05 mHa** (99.97% of $D_e$) at $R=1.4$, from ~27 radial DVR functions and 7
angular values, with **no explicit $r_{12}$ and no non-analytic basis functions**. Their
own single-centre control at $l_{\max}=7$ is 5.4 mHa — a ~110× coordinate-choice
penalty at comparable truncation order.

---

## 2. What Paper 12 actually does

Verified in `papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex` and
`geovac/neumann_vee.py`.

**Basis** (Eq. `basis_function`):
$\phi_i = \mathcal N_i e^{-\alpha(\xi_1+\xi_2)}\xi_1^{j}\xi_2^{k}\eta_1^{l}\eta_2^{m} + (1\leftrightarrow2)$.
**There is no $\varphi$ dependence anywhere.** The label `m` in that expression is the
power of $\eta_2$, not an azimuthal quantum number. One common $\alpha$, variationally
optimized; powers $j,l \le 3$; $N=72$ maximum.

**Kernel** (§ before Eq. `neumann_sigma`): "For $^1\Sigma_g^+$ states of homonuclear
diatomics, the wavefunction has $m=0$ symmetry. The azimuthal integration … onto the
$m=0$ component." The module docstring says it outright: *"The Coulomb kernel expands as
(for σ states, m=0 only)."*

**The error.** A $^1\Sigma_g^+$ state has **total** $M=m_1+m_2=0$. That does not make
each electron $m_i=0$. The $\pi^2$ ($m_1=+1,m_2=-1$), $\delta^2$ ($\pm2$), … configurations
are all $^1\Sigma_g^+$ and all carry $M=0$. Projecting the kernel onto its $m=0$ Neumann
component **and** using a $\varphi$-independent basis removes every one of them. What
remains is left–right and in–out correlation; angular correlation is structurally absent.

This also explains, correctly, the paper's own observation that the plateau is not
fixable by more basis functions — "additional σ-type basis functions with higher $\xi$
and $\eta$ powers cannot access the missing physics." That sentence is right. The
attribution that follows it is not.

---

## 3. The two routes that settle it

**Route 1 (external, verified).** TMR: same coordinates, polynomial-in-$\eta$ basis, no
$r_{12}$, no non-analytic terms, $m\ne0$ included → 0.05 mHa. If the cusp needed
$r_{12}^{1/2}$ and $r_{12}\ln r_{12}$ basis functions to get past 13 mHa, this
calculation could not exist.

**Route 2 (internal, already measured — Paper 15).** Paper 15 ran the controlled
experiment in its own geometry and reported it in
Table `tab:extended_convergence` and the surrounding text:

> "a constant ${\sim}6.6$ percentage-point additive offset from $\pi$ channels,
> independent of $l_{\max}$"

and in Table (comparison), at $l_{\max}=4$: σ-only 87.0%, σ+π **94.1%** (+7.1 points).

**Arithmetic.** $D_e^{\rm exact}=0.174475$ Ha. Paper 12: $E=-1.161304$, $D_e=0.161304$,
= 92.45%; gap = **13.17 mHa**. Paper 15's π offset, 6.6 pp of $D_e$, = **11.5 mHa** —
**87% of Paper 12's entire gap**, and the residual 1.0 pp (1.7 mHa) is the size of the
δ channels plus basis slack. Adding π to Paper 12 predicts ≈ 99.0%.

Two geometries, one external and one internal, agree that the missing quantity is
angular correlation, not cusp representability.

---

## 4. Consequences (all PI calls — §13.4 gate 4, and one touches §13.5)

1. **Paper 12 abstract + `sec:gap`.** The cusp diagnosis is the paper's closing
   argument and it is wrong as stated. What survives intact: the Neumann $V_{ee}$ *is*
   exact, it *does* beat numerical quadrature by 12–20 points at every basis size, and
   the plateau *is* a basis-space limit not an integration error. Only the identification
   of *which* limit changes: σ-only, not cusp non-analyticity. Also affected by the
   §9 summary-surface rule: abstract, `sec:gap`, and the "missing geometry" subsection,
   which uses the cusp reading to motivate the move to hyperspherical coordinates.
2. **Paper 15's comparative claim is not like-for-like, and the matched comparison
   REVERSES it (measured 2026-09-14).** "Exceeding the 92.4% achieved by prolate
   spheroidal CI" and "This confirms that the hyperspherical coordinates' cusp
   resolution provides a genuine advantage" compare Level 4 **σ+π** (94.1%) against
   Paper 12 **σ-only** (92.4%) and read the 1.7-point difference as a property of the
   coordinates. Run matched, prolate spheroidal wins: σ+π in Paper 12's own basis is
   **99.15%** (Sec. 5c) against Level 4's 94.1% at l_max=4 and 96.0% at l_max=6 with a
   Schwartz cusp correction — and TMR reach **99.97%** in the same coordinates. The
   sentence is not merely unsupported; its conclusion is the wrong way round.
3. **§5 hierarchy — flag only, do not touch (§13.5 hard prohibition).** Level 4 exists
   for H₂ partly because the cusp is a coordinate singularity in prolate spheroidal.
   TMR reaching 0.05 mHa there without cusp-adapted functions weakens that rationale.
   The hierarchy is PI-owned; this memo flags, nothing more.
4. **The Neumann/Poisson equivalence is worth stating.** TMR's Poisson solve and Paper
   12's $P_l(\xi_<)Q_l(\xi_>)$ Neumann kernel are the same mathematical object; Paper 12
   uses the *closed-form* Green's function and is therefore the **more** algebraic of
   the two. The 13 mHa versus 0.05 mHa difference is entirely the one-particle basis.
   This is "exact ≠ accurate" (memory `polyatomic_state_of_play` §1) demonstrated inside
   the corpus's own Level 2.

---

## 5. THE EXPERIMENT — RUN 2026-09-14. Diagnosis confirmed.

Two independent implementations, both validated against things computed by
other people or other machinery.  Prediction stated in the previous revision of
this memo was **98-99% of D_e**; measured **98.95%** in the native basis and
**99.10-99.42%** in a Gaussian basis.

### 5a. Independent route first: Gaussian-basis FCI

`debug/p12_m_channel_probe.py`, `debug/p12_m_channel_validate.py`.  Cartesian
Gaussians through the corpus's own McMurchie-Davidson engine
(`geovac/noci_engine.py`), where the sigma / pi / delta split is exact and the
restriction is imposed by orbital selection.  Touches none of the
prolate-spheroidal machinery, so it is a genuinely independent route
(memory rule `feedback_independent_route_crosscheck`).

Controls, all with answers known independently of this corpus:

| control | result | expected |
|:--|:--|:--|
| one H atom, 8s | -0.499888 | -0.5, basis-limited |
| H2 at R = 20 bohr | -0.999265 | -1.0, basis-limited |
| H2 RHF, 8s3p1d | -1.133270 | -1.13363 (HF limit), 0.36 mHa |

Measurement (H2, R = 1.4011):

| basis | sigma only | + \|m\| = 1 | + \|m\| = 2 | d(\|m\| >= 1) |
|:--|--:|--:|--:|--:|
| 6s2p  | -1.159791 (91.58%) | -1.170557 (97.75%) | -- | 10.77 mHa |
| 8s3p  | -1.160900 (92.22%) | -1.171839 (98.49%) | -- | 10.94 mHa |
| 8s3p2d| -1.161103 (92.34%) | -1.172907 (99.10%) | -1.173463 (99.42%) | 12.36 mHa |

**The sigma-only ceiling converges to -1.16110 (92.34%); Paper 12's sigma-only
value is -1.161304 (92.45%).  They agree to 0.2 mHa in two completely
unrelated bases.**  Paper 12 is therefore *at* its sigma-only basis-set limit --
its plateau is a real ceiling, and the ceiling is the restriction, not the
cusp.  Adding the azimuthal channels recovers 12.36 of the 13.17 mHa gap (94%);
the rest is this probe's own basis incompleteness (it too stops at 99.42%).

### 5b. The native experiment: general-m Neumann in Paper 12's own basis

`debug/prolate_ci_general_m.py` (module + validation), `debug/p12_mu_table.py`
(table).  The one-electron factor gains the azimuthal quantum number,

    u_{j,l,mu}(xi,eta) = xi^j eta^l (xi^2-1)^{mu/2} (1-eta^2)^{mu/2} e^{-a xi}

paired as m1 = +mu, m2 = -mu through cos(mu (phi1-phi2)) -- the M = 0,
^1Sigma_g^+ combination.  mu = 0 is Paper 12's basis exactly.  S, T and V_ne stay
**exact** (every integrand is a polynomial in xi and eta times e^{-2 a xi},
evaluated against the A_n recurrence and elementary eta moments) and are
diagonal in mu.  V_ee uses the full Neumann expansion, whose e^{i m dphi} factor
is what couples different mu -- the coupling Paper 12 discards.

**Validation (the load-bearing part).**
- The mu = 0 V_ee matrix reproduces the corpus's exact
  `geovac.neumann_vee.compute_vee_matrix_neumann` to **1.06e-9 relative**,
  elementwise, on all 27 x 27 entries.  So the ordered-xi quadrature, the eta
  moments, the Jacobian expansion and every prefactor are right.
- The mu = 0 energy at (j,l) = (2,2), N = 27 is **-1.160960** against Paper 12's
  published **-1.160961** -- 1e-6, against a number produced by completely
  different (exact recurrence) machinery.

**A numerical landmine found and fixed.**  The Neumann sum terminates exactly,
but only if the termination is *imposed*.  Integrating by parts m times (the
boundary terms vanish because (1-eta^2)^s has a zero of order s >= m at
eta = +-1) shows the eta moment is exactly zero for l > Q + 2s - m, and parity
kills (Q + l - m) odd.  Left to floating point, the ~1e10 Legendre-derivative
coefficients leave a ~1e-6 residue which the ~1e20 xi integral amplifies: at
l_neumann = 18 the energy came out **-3.2e8 Ha**.  With the rule imposed
explicitly the result is bit-identical from l_neumann = 8 to 24.  This is the
same guard the existing m = 0 code carries, generalised.

### 5c. Result

`debug/p12_mu_table.py` / `debug/prolate_ci_general_m.py`.  H2 at R = 1.4011,
Neumann V_ee (exact; the sum is capped at the proven cutoff l <= Q + 2s - m),
alpha scanned, and the generalised eigenproblem solved by **canonical
orthogonalisation** (see below).

| (j,l) | mu<=0 N (kept) | E | D_e% | mu<=1 N (kept) | E | D_e% |
|:--|--:|--:|--:|--:|--:|--:|
| (1,1) |   6 | -1.129445 | 74.19 |  12 | -1.142433 | 81.63 |
| (2,1) |  12 | -1.132371 | 75.87 |  24 | -1.145635 | 83.47 |
| (2,2) |  27 | -1.160960 | 92.25 |  54 | -1.172667 | 98.96 |
| (3,2) |  46 | -1.161155 | 92.37 |  92 (76) | -1.172722 | 99.00 |
| (3,3) |  72 (65) | -1.161246 | 92.42 | 144 (115) | -1.172881 | **99.09** |

Paper 12's own sigma-only column at the same truncations: -1.129606, -1.132373,
-1.160961, -1.161162, -1.161304.  Reproduced to **161, 1.6, 1.0, 6.7 and
58 uHa**.

**At the largest basis the azimuthal channels are worth +11.64 mHa,
92.42% -> 99.09% of D_e.**  The predicted band (98-99%) is met, and the two
independent routes land on top of each other: **99.09%** here, **99.10%** in the
Gaussian basis (Sec. 5a).

**Why this is the channel and not just a bigger basis.**  Paper 12's own
convergence table is the control: its sigma space is saturated -- N = 27 -> 46
-> 72 buys 92.2% -> 92.4% -> 92.4%, i.e. **0.34 mHa for 2.7x the functions**.
The Gaussian probe says the same independently: a full d shell moved the
sigma-only space 0.22 mHa.  Growth *along the sigma axis* does nothing; opening
the azimuthal axis is worth 11.6 mHa.  That is a channel effect, not a count
effect.

**A third finding, small but worth recording: Paper 12's largest basis is
linearly dependent.**  Measured cond(S) = **2.56e14** at (j,l) = (3,3), mu = 0 --
the N = 72 basis whose energy the paper headlines -- and 2.0e16 once mu = 1
doubles it, past double precision.  A raw `eigh(H, S)` there is not reliable:
at N = 144 it returned -79 Ha, wildly non-variational.  Under canonical
orthogonalisation (7 of 72 directions dropped at mu = 0, 29 of 144 at mu <= 1)
everything is stable and variational.  Paper 12's published -1.161304 sits
**58 uHa below** the conditioned value, i.e. its last two digits are
linear-dependence contamination rather than physics.  This changes no
conclusion -- it is 0.4% of the gap under discussion -- but the paper quotes six
decimals and only four are meaningful.

### 5d. What this does and does not establish

**Does:** Paper 12's 92.4% plateau is the sigma-only restriction.  Its
"the cusp requires non-analytic r12^{1/2}, r12 ln r12 terms" diagnosis is wrong
as the explanation of that plateau -- no r12 factors and no non-analytic
functions appear anywhere in either calculation above, and both clear 98%.

**Does not:** this is not yet a production Hamiltonian.  The |m| = 2 sector was
not obtained here at all -- the d^4 Q_l cancellations near xi = 1 defeat the
quadrature -- so the delta contribution (0.46-0.56 mHa) rests on the Gaussian
route alone.  S, T, V_ne are exact,
but the ordered xi integral in V_ee is computed by a graded-panel spectral
quadrature rather than by the closed-form A_l/B_l/X_l recurrences the m = 0 path
uses.  The corpus's "quadrature-free" property therefore is **not** yet extended
to mu > 0.  Doing so means generalising those three auxiliary tables to
associated Legendre functions (P_l^m, Q_l^m) -- a well-defined derivation sprint,
now with a validated numerical reference to check it against, which is the right
order to do it in.

## 6. The follow-on (was Sec. 5: the experiment this licenses)

**Extend `geovac/neumann_vee.py` to general $m$ and rerun Paper 12's H₂ CI.**

The general-$m$ Neumann expansion in prolate spheroidal coordinates is standard:
$$\frac{1}{r_{12}} = \frac{2}{R}\sum_{l}\sum_{m}(-1)^m(2l+1)\!\left[\frac{(l-|m|)!}{(l+|m|)!}\right]^{2}\!
P_l^{|m|}(\xi_<)\,Q_l^{|m|}(\xi_>)\,P_l^{|m|}(\eta_1)P_l^{|m|}(\eta_2)\,e^{im(\varphi_1-\varphi_2)}$$
The current module implements only the $m=0$ term.

**What already exists.** Paper 11 §"π, δ states ($m\ne0$)" has the associated-Legendre
$\eta$ machinery (tridiagonal/pentadiagonal recurrences for $P_l^{|m|}$) and the
associated-Laguerre radial basis, with the $m\ne0$ content reduced to the single
transcendental seed $e^aE_1(a)$ (Track J). `prolate_spheroidal_lattice.py`,
`prolate_scf.py` and `two_center_eri.py` all carry associated-Legendre code. The
one-electron half of this is built.

**What is new.** $Q_l^{|m|}(\xi)$ moments (the $B_l$ analogue at general $m$), basis
functions carrying $e^{im\varphi}$, and the $m$-sum in the assembly.

**Prediction, stated before the run (falsifier).** With $|m|\le1$ at the same
$(j,l)\le3$, $N\approx$ 200: **98–99% of $D_e$**, i.e. $E \approx -1.172$ to $-1.173$ Ha.
Below 97% means something other than the π channels is also binding; above 99.5% at this
basis size would be suspicious.

**Why it is worth doing.** It would give the corpus a **quadrature-free, discrete-angular,
Gaunt-coupled, non-overcomplete two-electron molecular Hamiltonian at (or near) chemical
accuracy** — the qubit-encoding anchor the accuracy scan said the framework lacks, built
from two components it already owns. Honest cap, unchanged: diatomic, two electrons.
Three centres remain Paper 59's genus-1 wall.

**Second-order opportunity, not required.** TMR's DVR diagonality is *distance* sparsity
in the radial index — the kind memory `two_kinds_of_sparsity` says GeoVac does not have.
A DVR radial basis is permitted by §4 (numerical radial amplitude within a fixed channel
structure) but forfeits closed forms. Worth a separate assessment, not this sprint.

---

## 7. Verification ledger

**Verified at source by the PM:** TMR 2010 (OSTI manuscript, full text, all quotations
above); Paper 12 `.tex` Eq. `basis_function`, `neumann_sigma` paragraph, `tab:convergence`,
`sec:gap`; `geovac/neumann_vee.py` docstring and `m`-usage; Paper 15
`tab:extended_convergence`, the σ-only/σ+π comparison table, and the "6.6 percentage-point"
sentence; Paper 11 §"π, δ states".

**Not verified:** Wolniewicz's −1.17447 is taken from TMR's citation of it, not read at
source (the corpus separately uses Kolos–Wolniewicz −1.174475 at $R=1.4011$; consistent).
The σ/π/δ decomposition of H₂ correlation energy was searched for in the literature and
**not** cleanly found — the 6.6 pp figure used above is the corpus's own Paper 15
measurement, which is why Route 2 is stated as internal.

**Not done:** no calculation was run this session. §5's prediction is a prediction.
