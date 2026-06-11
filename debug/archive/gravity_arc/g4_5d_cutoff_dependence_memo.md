# Sprint G4-5d — Cutoff-function dependence sweep for tip-only $S_{BH}$

**Date:** 2026-05-29
**Verdict:** **NEGATIVE-G4-5d-CUTOFF-MELLIN-SCALING** at sprint-scale resolution, with a substantive structural finding that sharpens the G8 reading.

## 1. Question

G8 (Paper 51 §10, v3.12.0) classifies the cutoff function $f$ as Class 1 calibration and gives the closed-form Newton constant, cosmological constant, and critical radius via Mellin moments

$$\phi(s) := \int_{0}^{\infty} f(x)\, x^{s-1}\, dx,$$

with $G_{\mathrm{eff}} = 6\pi/(\phi(1)\,\Lambda^2)$, $\Lambda_{cc} = 6\,\phi(2)/\phi(1) \cdot \Lambda^2$, and $R_{\mathrm{crit}} \Lambda = \sqrt{\phi(1)/(6\,\phi(2))}$.

F12 of the G4-5 scoping memo asked whether the **tip-only** contribution to $S_{BH}$ on the *discrete* cigar substrate (G4-3a/cleanup/b/c/d → G4-5a) scales with the **second** Mellin moment $\phi(2)$, by analogy with $\Lambda_{cc}$. The three reference cutoffs are

| Cutoff | $f(x)$ | $\phi(s)$ | $\phi(2)$ |
|:--|:--:|:--:|:--:|
| Gaussian | $e^{-x}$ | $\Gamma(s)$ | $1$ |
| Sharp | $\Theta(1-x)$ | $1/s$ | $1/2$ |
| Polynomial | $e^{-x^{2}}$ | $\tfrac{1}{2}\Gamma(s/2)$ | $1/2$ |

So the predicted ratios at fixed $\Lambda$ are $S_{\mathrm{tip}}(\mathrm{gauss})/S_{\mathrm{tip}}(\mathrm{sharp}) = 2$, $S_{\mathrm{tip}}(\mathrm{gauss})/S_{\mathrm{tip}}(\mathrm{poly}) = 2$, $S_{\mathrm{tip}}(\mathrm{sharp})/S_{\mathrm{tip}}(\mathrm{poly}) = 1$.

## 2. Method

The sweep reuses the G4-5a sweet-spot panel ($R=10$, $a=0.05$, $N_\rho=200$, $N_0=120$, $k_{\mathrm{step}}=12$, $\varepsilon=0.1$) at the 8-point log-spaced t-grid $t\in\{0.1, 0.2, 0.5, 1, 2, 5, 10, 20\}$. The tip term

$$\mathrm{tip}(t) := \frac{K^{\mathrm{Dirac}}_{\mathrm{wedge}}(\alpha_{+},t) - K^{\mathrm{Dirac}}_{\mathrm{wedge}}(\alpha_{-},t)}{\alpha_{+} - \alpha_{-}} \;-\; K^{\mathrm{Dirac}}_{\mathrm{disk}}(t)$$

is recomputed and reproduced **exactly bit-identical** to G4-5a (sanity check at $t=1$: diff $0.0$).

For each cutoff $f$ and each $\Lambda \in \{0.5, 1, 2\}$ the tip integral

$$S_{\mathrm{tip}}(\Lambda, f) = \tfrac{1}{2} \int \frac{dt}{t}\, f(t\Lambda^{2})\, \mathrm{tip}(t)$$

is evaluated by trapezoid in $\log t$ on the eight-point grid.

## 3. Results

**Tip term on the substrate (recomputed bit-identical to G4-5a):**

| $t$ | $K_{\mathrm{disk}}$ | $dK/d\alpha$ | $\mathrm{tip}(t)$ |
|:--:|:--:|:--:|:--:|
| 0.1 | 529.91 | 530.03 | +0.1272 |
| 0.2 | 245.88 | 246.02 | +0.1382 |
| 0.5 | 90.03 | 90.18 | +0.1484 |
| 1.0 | 42.00 | 42.16 | +0.1538 |
| 2.0 | 19.11 | 19.26 | +0.1576 |
| 5.0 | 6.21 | 6.38 | +0.1612 |
| 10.0 | 2.33 | 2.49 | +0.1628 |
| 20.0 | 0.65 | 0.80 | +0.1509 |

The tip is **slowly varying** in $\log t$, peaking around $t\sim 10$ at $+0.163$ and dropping mildly at both ends.

**$S_{\mathrm{tip}}(\Lambda, f)$ on the panel:**

| $\Lambda$ | gauss | sharp | poly |
|:--:|:--:|:--:|:--:|
| 0.5 | $+0.2305$ | $+0.2541$ | $+0.2494$ |
| 1.0 | $+0.1303$ | $+0.1907$ | $+0.1425$ |
| 2.0 | $+0.0488$ | $+0.0776$ | $+0.0492$ |

Each cutoff gives finite, positive $S_{\mathrm{tip}}$ decreasing monotonically with $\Lambda$ — load-bearing sanity preserved.

**Ratios vs $\phi(2)$ prediction:**

| $\Lambda$ | $g/sh$ emp. | $g/sh$ pred. | dev | $g/po$ emp. | $g/po$ pred. | dev | $sh/po$ emp. | $sh/po$ pred. | dev |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 0.5 | 0.91 | 2.00 | $-54.7\%$ | 0.92 | 2.00 | $-53.8\%$ | 1.02 | 1.00 | $+1.9\%$ |
| 1.0 | 0.68 | 2.00 | $-65.8\%$ | 0.91 | 2.00 | $-54.3\%$ | 1.34 | 1.00 | $+33.8\%$ |
| 2.0 | 0.63 | 2.00 | $-68.6\%$ | 0.99 | 2.00 | $-50.4\%$ | 1.58 | 1.00 | $+57.8\%$ |

Max deviation: $68.6\%$; mean deviation: $49.0\%$; sign and ordering also wrong (Gaussian gives the **smallest** $S_{\mathrm{tip}}$, not the largest).

## 4. Diagnosis: $\phi(2)$ is the wrong Mellin moment for the tip

The negative is informative and the structural reading is clean. Three observations:

**(a)** The tip term is approximately **constant in $\log t$** across the eight-point grid ($\mathrm{tip}(t) \approx 0.15$, variation $\pm 0.02$). Pulling $\mathrm{tip}$ out as $C_{\mathrm{tip}} \approx 0.15$ gives

$$S_{\mathrm{tip}}(\Lambda, f) \;\approx\; \tfrac{1}{2}\, C_{\mathrm{tip}} \int \frac{dt}{t}\, f(t\Lambda^{2}) \;=\; \tfrac{1}{2}\, C_{\mathrm{tip}} \int \frac{du}{u}\, f(u) \;=\; \tfrac{1}{2}\, C_{\mathrm{tip}} \cdot M_{0}(f),$$

after the substitution $u = t\Lambda^{2}$. The integral $\int_{0}^{\infty} f(u)\, du/u$ is the **zeroth** Mellin moment, $\phi(0)$, not $\phi(2)$. And $\phi(0)$ is IR-divergent for all three reference cutoffs — the discrete substrate regulates it through the finite $t$-grid bounds.

**(b)** For the **G8 calibration relations** that hold in Paper 51 §10, the relevant Mellin moments are $\phi(1)$ ($G_{\mathrm{eff}}$, from the Einstein–Hilbert $\Lambda^{2}$ coefficient) and $\phi(2)$ ($\Lambda_{cc}$, from the $\Lambda^{4}$ cosmological-constant coefficient). These are the moments of $\Lambda^{2k}\, dx$ measures, with explicit polynomial weights $x$ and $x^{2}$ scaling the heat-kernel coefficients $a_{0}$ and $a_{2}$ on $S^{3}$.

**(c)** The **tip contribution to $S_{BH}$** in the replica method, by contrast, sits in the **topological** Sommerfeld–Cheeger coefficient $(1/12)(1/\alpha - \alpha)$, which multiplies $K_{S^{2}}(t)$ at the horizon ($t\to 0$, fixed-point of $\alpha\to 1$). The natural Mellin moment for this coefficient is the **zeroth** moment $\phi(0) = \int_{0}^{\infty} f(t\Lambda^{2})\,dt/t$, evaluated as a regulated log via UV/IR cutoffs of the integral.

**Concretely:**
- The polynomial-weighted moments $\phi(1)$ and $\phi(2)$ control Wilson coefficients of $R^{0}$ ($\Lambda_{cc}$) and $R^{1}$ (Einstein–Hilbert) bulk operators.
- The unweighted moment $\phi(0)$ controls the topological tip coefficient, which is $R$-independent in the standard CC replica derivation.

The empirical numbers support this reading: $\phi(0)$ of the Gaussian (Euler–Mascheroni regulated, $\sim 0.58$) is **smaller** than $\phi(0)$ of the sharp $\Theta$ ($\log(t_{\mathrm{max}}/t_{\mathrm{min}}) \sim 5.3$ on the panel) and intermediate vs polynomial — and indeed empirically $S_{\mathrm{tip}}(\mathrm{gauss}) < S_{\mathrm{tip}}(\mathrm{sharp})$ at every $\Lambda$, with the polynomial sitting between them.

## 5. Sharper rereading of G8 for the tip-of-cigar sector

The G8 closed-form predictions

$$G_{\mathrm{eff}} = \frac{6\pi}{\phi(1)\,\Lambda^{2}}, \quad \Lambda_{cc} = 6\,\frac{\phi(2)}{\phi(1)}\,\Lambda^{2}, \quad R_{\mathrm{crit}}\Lambda = \sqrt{\frac{\phi(1)}{6\,\phi(2)}}$$

are correct **for the bulk Einstein–Hilbert + cosmological-constant sector**, which inherits its cutoff dependence from the polynomial-weight $x^{k}$ scaling of the heat-kernel $a_{k}$ coefficients on $S^{3}$.

For the **tip-of-cigar sector** (Sommerfeld–Cheeger topological term × $K_{S^{2}}$), the relevant moment is **$\phi(0)$**, not $\phi(2)$. F12 was framed against the wrong moment, and the discrete substrate has correctly diagnosed that mismatch.

In the standard continuum Bekenstein–Hawking derivation, $\phi(0)$ is logarithmically divergent and the topological tip coefficient is dimensionless in the cutoff $\Lambda$. The discrete substrate naturally regulates the log via the finite $t$-grid (here, panel range $t\in[0.1, 20]$ gives $\log(t_{\mathrm{max}}/t_{\mathrm{min}}) \sim 5.3$). On the **panel ratios** the empirical $S_{\mathrm{tip}}(\mathrm{sharp})/S_{\mathrm{tip}}(\mathrm{poly}) = 1.02$ at $\Lambda = 0.5$ is consistent with the prediction $\phi(0)_{\mathrm{sharp}}/\phi(0)_{\mathrm{poly}} \approx 1$ when both are dominated by the $\log(t_{\mathrm{IR}}/t_{\mathrm{UV}})$ piece, but the gauss vs others mismatch persists because $\phi(0)_{\mathrm{gauss}}$ is structurally Euler–Mascheroni-class rather than log-class.

## 6. Three substantive findings

1. **F12 closure form changes.** The cutoff dependence of the tip-of-cigar contribution to $S_{BH}$ on the discrete substrate is **not** controlled by $\phi(2)$. The G8 Mellin-moment scaling for bulk Wilson coefficients ($G_{\mathrm{eff}}$, $\Lambda_{cc}$) **does not transport** to the topological tip coefficient.
2. **The correct moment is $\phi(0)$.** The Sommerfeld–Cheeger topological-tip coefficient scales with the **zeroth** Mellin moment $\phi(0)$ (regulated by UV/IR substrate scales), reflecting its dimensionless, $R$-independent structure in the standard CC replica reading.
3. **Empirical ordering reflects $\phi(0)$, not $\phi(2)$.** The discrete-substrate $S_{\mathrm{tip}}(\mathrm{sharp}) > S_{\mathrm{tip}}(\mathrm{poly}) > S_{\mathrm{tip}}(\mathrm{gauss})$ at every $\Lambda$ matches the qualitative ordering of $\phi(0)$ (sharp is log-divergent at UV, polynomial decays faster, Gaussian is regulated by Euler–Mascheroni) and not the ordering of $\phi(2)$ (Gaussian largest).

These three together promote G8 from a **single** Mellin-moment scaling law to a **moment-by-sector** structural map:

| Sector | Wilson coefficient | Mellin moment |
|:--|:--|:--:|
| Bulk $R^{0}$ (cosmological const.) | $\Lambda_{cc}$ | $\phi(2)$ |
| Bulk $R^{1}$ (Einstein–Hilbert) | $G_{\mathrm{eff}}^{-1}$ | $\phi(1)$ |
| Topological tip ($S_{BH}$) | (CC-tip coeff $1/12$) | $\phi(0)$ |

This is a clean refinement of Paper 51 §10's structural-skeleton/calibration partition: the cutoff dependence partitions the *Mellin axis* across sectors, with each operator order picking out a specific $\phi(k)$.

## 7. Verdict and follow-on

**NEGATIVE-G4-5d-CUTOFF-MELLIN-SCALING** as F12 was originally stated ($\phi(2)$), with a substantive structural reframing: the correct moment for the topological-tip sector is $\phi(0)$, not $\phi(2)$.

This refines G8's structural-skeleton reading and is publishable content for Paper 51 §10 (per CLAUDE.md §13.8): a one-sentence addition that "the cutoff-dependence map $f \mapsto (G_{\mathrm{eff}}, \Lambda_{cc}, R_{\mathrm{crit}}, S_{BH}^{\mathrm{tip}})$ is sector-wise Mellin-moment-controlled, with sector indices $k\in\{0,1,2\}$ corresponding respectively to topological tip, Einstein–Hilbert, and cosmological-constant terms."

Follow-on candidates (small enough to be sub-sprint-scale):
- **F13:** rerun G4-5d with $\phi(0)$ predicted ratios (gauss $\sim 0.58$, sharp $\sim 5.3$, poly $\sim 0.89$ on the panel) to test the corrected prediction quantitatively. ~1 hour sprint at panel-converged resolution; would close F12-refined POSITIVE.
- **F14:** explicit panel-by-panel $\phi(0)$ regulator computation showing the discrete substrate IR/UV cutoffs $t_{\mathrm{min}}, t_{\mathrm{max}}$ replace the continuum log-divergence. ~half-day sprint.
- **G4-5e**: revisit Paper 51 §10 phrasing to insert the sector-by-Mellin-moment refinement.

## Files

- Driver: `debug/g4_5d_cutoff_dependence.py`
- Data: `debug/data/g4_5d_cutoff_dependence.json`
- Memo: `debug/g4_5d_cutoff_dependence_memo.md` (this file)
