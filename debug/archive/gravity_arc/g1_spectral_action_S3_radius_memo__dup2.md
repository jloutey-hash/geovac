# Sprint G1 — Spectral action on $S^3_R$ parameterized by radius

**Date:** 2026-05-28
**Path:** Gravity Path 1 — minimum-cost test of whether the framework's spectral action has Connes–Chamseddine gravity-like structure on a one-parameter family of GeoVac geometries.
**Verdict:** POSITIVE-WITH-NUANCE. The asymptotic spectral action on $S^3_R$ has EXACTLY two terms with NO power-law corrections — substantive new structural finding. The two coefficients have opposite signs and the asymptotic admits a formal extremum at $u_{\rm crit} \cdot \Lambda = O(1)$, the GeoVac analog of CC's de-Sitter radius selection. However, the formal extremum sits in the IR regime where the asymptotic does NOT faithfully represent the QM exact spectral action. The QM exact sum is monotonically increasing in $u = \Lambda R$ and has no literal extremum. Under the CC spectral-action principle (which treats the asymptotic AS the classical action), GeoVac reproduces CC's Einstein–Hilbert + cosmological constant structure with a bonus exactness theorem; under a literal-QM reading, there is no minimum.

## 1. Setup

Camporesi–Higuchi Dirac spectrum on $S^3$ of radius $R$:
$$|\lambda_n(R)| = \frac{n + 3/2}{R}, \qquad g_n = 2(n+1)(n+2)$$

(Full 4-component Dirac, matching the convention in `geovac/dirac_s3.py` and Paper 28.) The spectral zeta on $S^3_R$ scales cleanly:
$$\zeta_R(s) := \sum_n g_n\, |\lambda_n(R)|^{-2s} = R^{2s}\, \zeta_{\rm unit}(s)$$

where $\zeta_{\rm unit}(s) = \sum_n 2(n+1)(n+2)\,(n+3/2)^{-2s}$.

Paper 28's two-term exactness theorem (Jacobi $\vartheta$ inversion on the half-integer-shifted Dirac spectrum) states
$$\mathrm{Tr}\,e^{-tD^2}\Bigm|_{\rm unit} = \frac{\sqrt{\pi}}{2}\,t^{-3/2} - \frac{\sqrt{\pi}}{4}\,t^{-1/2} + O\!\left(e^{-\pi^2/t}\right)$$

with all higher Seeley–DeWitt coefficients $a_k = 0$ for $k \ge 2$.

## 2. The substantive new finding: $\zeta_{\rm unit}(-k) = 0$ for all $k \ge 0$

**Spectral-zeta translation of two-term exactness.** Via Mellin inversion,
$$\mathrm{Tr}\,e^{-tD^2}\bigm|_{\rm unit} = \frac{1}{2\pi i} \int \Gamma(s)\,\zeta_{\rm unit}(s)\,t^{-s}\,ds$$

The expansion of $\mathrm{Tr}\,e^{-tD^2}$ at small $t$ picks up residues from $\Gamma(s)\zeta_{\rm unit}(s)$ in the right half-plane. The two power-law terms in two-term exactness correspond to simple poles of $\zeta_{\rm unit}(s)$ at $s = 3/2$ (residue 1) and $s = 1/2$ (residue $-1/4$). Two-term exactness states there are *no other power-law contributions* — in particular, no Taylor terms $t^k$ for $k \ge 0$.

But $\Gamma(s)$ has simple poles at $s = 0, -1, -2, \dots$ with residues $(-1)^k/k!$. These would contribute Taylor terms $\zeta_{\rm unit}(-k) \cdot \frac{(-1)^k}{k!}\,t^k$ to $\mathrm{Tr}\,e^{-tD^2}$ UNLESS $\zeta_{\rm unit}(-k) = 0$ for all non-negative integers $k$.

**Symbolic verification (this sprint, $k = 0,\dots,5$):**
$$\zeta_{\rm unit}(-k) = 2\zeta_H(-2k-2, 3/2) - \tfrac{1}{2}\zeta_H(-2k, 3/2)$$

Using $\zeta_H(-n,a) = -B_{n+1}(a)/(n+1)$:
$$\zeta_{\rm unit}(-k) = -\frac{2 B_{2k+3}(3/2)}{2k+3} + \frac{B_{2k+1}(3/2)}{2(2k+1)}$$

This vanishes for $k = 0, 1, 2, 3, 4, 5$ — exact rational zero in sympy at every $k$. The equivalent identity
$$4(2k+1)\,B_{2k+3}(3/2) = (2k+3)\,B_{2k+1}(3/2)$$

is a Bernoulli-polynomial closed form. So $\zeta_{\rm unit}(-k) = 0$ is *equivalent* to two-term exactness on the heat-kernel side — they are two sides of the same theorem.

## 3. Exact two-term spectral action

For cutoff $f$ with Mellin transform $\phi(s) = \int_0^\infty f(x)\,x^{s-1}\,dx$:
$$S(R,\Lambda) := \mathrm{Tr}\,f(D^2/\Lambda^2) = \frac{1}{2\pi i}\int \phi(s)\,(\Lambda R)^{2s}\,\zeta_{\rm unit}(s)\,ds$$

The contour shift through both poles of $\zeta_{\rm unit}$ and through the poles of $\phi(s)$ at non-positive integers gives:
- $s = 3/2$ residue: $\phi(3/2)\cdot(\Lambda R)^3 \cdot 1$
- $s = 1/2$ residue: $\phi(1/2)\cdot(\Lambda R) \cdot (-1/4)$
- $s = -k$ residue from $\phi$: $\mathrm{Res}(\phi,-k)\cdot(\Lambda R)^{-2k}\cdot \zeta_{\rm unit}(-k) = 0$ (since $\zeta_{\rm unit}(-k) = 0$)

**Result:**
$$\boxed{\;S(R,\Lambda) = \phi(3/2)\,(\Lambda R)^3 - \tfrac{1}{4}\,\phi(1/2)\,(\Lambda R) + O(\text{exp small in }(\Lambda R)^2)\;}$$

This is EXACT modulo exp-small corrections, for ANY cutoff function $f$ whose Mellin transform is meromorphic with poles only at non-positive integers. No power-law corrections to any order. The two-term form is structurally forced by Paper 28's exactness theorem.

**Dimensional reading.** $\phi(3/2)\cdot(\Lambda R)^3$ is the cosmological-constant slot (Volume$(S^3_R) = 2\pi^2 R^3$). $-\tfrac{1}{4}\phi(1/2)\cdot(\Lambda R)$ is the Einstein–Hilbert slot ($R_{\rm scalar}\cdot V \propto R$ on $S^3_R$). The two-term spectral action thus has the exact form of Einstein–Hilbert + cosmological constant — the structural target of CC spectral action gravity.

## 4. Numerical verification (three cutoffs)

`debug/g1_spectral_action_S3_radius.py`, `debug/data/g1_spectral_action_S3_radius.json`.

| Cutoff | $\phi(3/2)$ | $\phi(1/2)$ | $u_{\rm crit}$ | $S_{\rm asymp}(u_{\rm crit})$ |
|--------|-------------|-------------|----------------|---------|
| Gaussian $e^{-x}$ | $\sqrt{\pi}/2 \approx 0.886$ | $\sqrt{\pi} \approx 1.772$ | $1/\sqrt{6} \approx 0.408$ | $-0.121$ |
| Polynomial $e^{-x^2}$ | $\tfrac{1}{2}\Gamma(3/4) \approx 0.613$ | $\tfrac{1}{2}\Gamma(1/4) \approx 1.813$ | $\approx 0.497$ | $-0.150$ |
| Sharp $\Theta(1-x)$ | $2/3$ | $2$ | $1/2$ | $-1/6$ |

Formal extremum derived from $dS_{\rm asymp}/du = 0$: $u_{\rm crit} = \sqrt{B/(3A)}$ where $A = \phi(3/2)$, $B = \tfrac{1}{4}\phi(1/2)$. All three cutoffs give $u_{\rm crit}$ in $[0.41, 0.50]$ — all $O(1)$ as predicted.

**QM exact vs asymptotic agreement** (Gaussian):

| $u = \Lambda R$ | $S_{\rm QM}$ | $S_{\rm asymp}$ | rel diff |
|---|---|---|---|
| 0.3 | $5.6 \times 10^{-11}$ | $-0.109$ | $10^9$ (huge) |
| 0.5 | $4.9 \times 10^{-4}$ | $-0.111$ | 225 |
| 1.0 | $0.445$ | $0.443$ | $4 \times 10^{-3}$ |
| 1.5 | $2.326$ | $2.326$ | $3 \times 10^{-8}$ |
| 2.0 | $6.204$ | $6.204$ | $1 \times 10^{-15}$ |
| 3.0 | $22.60$ | $22.60$ | $1 \times 10^{-37}$ |
| 10.0 | $881.80$ | $881.80$ | $3 \times 10^{-51}$ |

The asymptotic matches the QM exact sum to MACHINE PRECISION for $u \ge 2$ (Gaussian). For $u \le 0.7$ the asymptotic is in *disagreement* with the QM exact: the QM sum is tiny (modes Boltzmann-suppressed by $\exp(-(3/2)^2/u^2)$) while the asymptotic continues linearly into negative values.

**Regime boundary** (smallest $u$ where $|S_{\rm QM} - S_{\rm asymp}|/|S_{\rm asymp}| < 0.01$):

| Cutoff | Boundary $u^*$ | $u_{\rm crit}$ | Extremum inside regime? |
|--------|---------|---|---|
| Gaussian | 1.0 | 0.41 | NO |
| Polynomial | 2.0 | 0.50 | NO |
| Sharp | 7.0 | 0.50 | NO |

**The formal extremum lies in the IR (small-$u$) regime where the asymptotic does NOT represent the QM action.** This is the key honest finding of the sprint.

## 5. Structural reading

**Two readings of "does the spectral action select a preferred geometry?":**

### Reading A — Literal QM action

$S_{\rm QM}(u) = \sum_n g_n\,f(\lambda_n^2/\Lambda^2)$ is the *literal* quantum action. For all three cutoffs, $S_{\rm QM}(u)$ is monotonically increasing from 0 (no modes pass the cutoff at $u \to 0$) to $A u^3$ (UV asymptotic, $u \to \infty$). **No extremum.**

Under this reading, GeoVac's spectral action does NOT select a preferred radius.

### Reading B — CC spectral-action principle

The CC spectral-action principle (Chamseddine–Connes 1997, Connes–Marcolli 2008) treats the *asymptotic expansion* as the effective classical action for the metric, regardless of whether the extremum sits at finite $u$ or in the deep IR. Under this reading:

$$S_{\rm classical}(R,\Lambda) := \phi(3/2)\,(\Lambda R)^3 - \tfrac{1}{4}\phi(1/2)\,(\Lambda R)$$

has a unique extremum at $u_{\rm crit} = \sqrt{B/(3A)} = O(1)$, which is the GeoVac analog of CC's de-Sitter radius selection.

This is the standard reading in NCG-style gravity: the spectral action's role is to *encode* Einstein–Hilbert + cosmological constant at the level of the asymptotic local invariants. The "extremum picks out a metric" claim is about the asymptotic, not the literal QM sum.

**Under Reading B, GeoVac reproduces CC-style spectral-action gravity, with a STRONGER STRUCTURAL RESULT: the asymptotic terminates exactly at two terms with no higher-curvature corrections.** This is a consequence of Paper 28's two-term exactness theorem, which is specific to the half-integer-shifted Dirac spectrum on $S^3$. CC's spectral action on a generic 4-manifold has higher-curvature corrections at every order; GeoVac on $S^3_R$ has *none* beyond Einstein–Hilbert.

## 6. What this DOES and DOES NOT establish

**Reached:**
- The asymptotic spectral action on $S^3_R$ has EXACTLY two terms in the Mellin-residue expansion, with NO power-law corrections at any order.
- The two coefficients have opposite signs (for all three natural cutoffs tested), giving a formal extremum at $u_{\rm crit} = O(1)$.
- The structural identity $\zeta_{\rm unit}(-k) = 0$ for $k \ge 0$ (the spectral-zeta side of two-term exactness) — substantive new statement, verified symbolically.
- The framework's spectral action has the structural form of CC Einstein–Hilbert + cosmological constant.
- Under the CC spectral-action principle, this gives a formal de-Sitter-radius selection $R_{\rm crit} \sim 1/\Lambda$.

**Not reached:**
- A literal QM-level extremum of the spectral action (the QM sum is monotonic).
- Derivation of Einstein equations from variation over a one-parameter family — the radius $R$ varies the *size* of $S^3$, not the *shape* of the metric. Real Einstein-equation tests require variation over $g_{\mu\nu}$ at fixed topology.
- A prediction of the cosmological-constant scale. As in CC's spectral action, the cutoff function moments $\phi(3/2)$ and $\phi(1/2)$ are external inputs, not fixed by GeoVac. This is the same gap that gives the standard cosmological-constant problem in CC.
- Coupling to matter dynamics. Matter fields would couple to the geometry through $D \to D + \omega$ inner fluctuations (Marcolli–vS lineage), which would source the spectral action — not done here.

**Scope:** This sprint tested the structural form of the spectral action on the family $\{S^3_R : R > 0\}$. It is the cheapest possible test that the spectral-action machinery does what it's supposed to on the GeoVac substrate. Two-term exactness gives a cleaner result than CC's continuum spectral action (no higher-curvature corrections), but inherits the same separation between "spectral action principle" (Reading B) and "literal QM minimum" (Reading A).

## 7. Forward paths

Three natural extensions, in order of cost:

**G2 — Higher-curvature corrections on $S^3$ × $S^1_\beta$.** Sprint TD Track 4 already reproduced $T_H = 1/(8\pi M)$ from the cigar geometry. The spectral action on $S^3 \times S^1_\beta$ is the natural two-parameter family to test next: does the two-term exactness extend? Does the Stefan–Boltzmann free-energy structure emerge from heat-kernel expansion? Reachable, 2–4 weeks.

**G3 — Linearized graviton spectrum on $S^3$.** Lichnerowicz Laplacian on symmetric 2-tensors. Continuum spectrum known. Test if GeoVac discrete substrate carries a graviton mode, see if it converges to the continuum result with Paper 38–style propinquity rate. 2–4 weeks.

**G4 — Cigar parameterized by $M$.** $T_H = 1/(8\pi M)$ already in hand. Computing the spectral action across the $M$-family and chasing Bekenstein–Hawking $S = A/4$ is the natural multi-month extension. Paper 48 §11 (Bousso–Casini–Fisher–Maldacena bulk-dual) is the gravity-side anchor. Multi-month.

Path 1's structural lesson — "two-term exactness gives EXACT Einstein–Hilbert + cosmological constant, with no higher-curvature corrections at any order" — likely persists in some form across these extensions, since it's driven by the half-integer Dirac spectrum which is preserved across the family.

## 8. Cross-references in the GeoVac corpus

- **Paper 28 §4** (two-term exactness theorem) — spectral-zeta translation of two-term exactness is now $\zeta_{\rm unit}(-k) = 0 \forall k \ge 0$, a new identity supplementing the heat-kernel statement.
- **Paper 18 §III.7** (master Mellin engine) — the M2 mechanism (Seeley–DeWitt, $k=2$). This sprint is M2 on $S^3$.
- **Paper 32 §VIII** (case-exhaustion theorem) — the two-term form is the framework's M2 signature; the structural exactness on $S^3$ is specific to the half-integer Dirac spectrum.
- **Paper 50 §8** (M2 ↔ holographic Weyl anomaly) — companion finding on the AdS/CFT side; the same M2 mechanism that runs in CC-style gravity runs in holographic anomaly coefficients.
- **CLAUDE.md §1.7 WH5** — $\alpha$ as projection constant; cosmological constant inherits the same calibration-data status (external input via cutoff function moments).
- **memory/external_input_three_class_partition.md** — gravitational coupling is plausibly Class 1 calibration data, like Yukawas and $\alpha$.

## 9. Files produced

- `debug/g1_spectral_action_S3_radius.py` (~340 lines) — driver
- `debug/data/g1_spectral_action_S3_radius.json` — structured numerical results
- `debug/g1_spectral_action_S3_radius_memo.md` — this memo
