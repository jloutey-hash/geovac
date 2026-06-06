# Sprint Q5'-Stage1-Sub-Sprint-2b — Formal Mellin transform of $\phi_0^{\mathrm{odd}}$ on the truncated CH spectral triple at $n_{\max} = 2$ (and $n_{\max} = 3$)

**Date:** 2026-06-05 (same session as Q5'-Stage1-JLO-nmax2)
**Sprint:** Q5' Stage 1, Sub-Sprint 2b (second of three named Sub-Sprint-2 sub-targets from JLO memo line 207)
**Driver:** `debug/compute_phi0_mellin.py`
**Data:** `debug/data/sprint_q5p_2b_phi0_mellin_data.json`
**Wall:** 1.2 s
**Discipline:** bit-exact `sympy.Rational` for finite-cutoff $\zeta_{D^2}(s)$; symbolic `sympy.gamma`; matrix-power route ($\mathrm{Tr}((D^2)^{-s})$) for full-$D$ at integer $s$, bypassing degree-4 Galois eigenvalues; no PSLQ.

---

## TL;DR

**Verdict: POSITIVE-WITH-STRUCTURAL-FINDING.** The formal Mellin transform $\mathcal{M}[\phi_0^{\mathrm{odd}}(1; t)](s) = \Gamma(s) \cdot \zeta_{D^2}^{(n_{\max})}(s)$ closes in bit-exact $\mathbb{Q}[\Gamma(s)]$ at finite cutoff $n_{\max} \in \{2, 3\}$, for both diagonal $\Lambda$ AND full $D = \Lambda + \kappa A$. The $s \to 0$ residue is exactly $N = \dim \mathcal{H} \in \mathbb{Z}$ (a rational integer), with NO $\pi$ at the finite-cutoff Laurent level.

**The structural finding:** **M1 Hopf-base measure $\pi$ does NOT inject at the finite-cutoff Mellin Laurent expansion at $s = 0$. It enters only at the continuum limit, through the Weyl-law short-time asymptotic $\mathrm{Tr}(e^{-tD^2}) \sim (2\pi^2)/(4\pi t)^{3/2} + O(t^{-1/2})$.** The finite-cutoff Mellin is structurally on the OTHER side of the master Mellin engine: it sees the discrete spectrum's exact rational moments, not the continuum measure's $\pi$ content.

This is **NOT** a failure of the M1 reading — it's a **sharp structural localization of where M1 lives**. At finite cutoff, the spectral zeta is rational; M1's $\pi$ enters only when one releases the discrete cutoff and asks for the continuum heat-kernel coefficient. The cosmic-Galois symbol level at $n_{\max} = 2$ is integer-valued (Sub-Sprint 1: $(16, 84, 36) \in \mathbb{Z}^3$); the M1 $\pi$ is the limit-of-cutoffs symbol, not the at-cutoff symbol.

The result aligns with **Paper 18 §III.7** master-Mellin reading: M1 sub-mechanism = "trivial heat kernel on a sphere collapses to the Hopf-base measure $\mathrm{Vol}(S^d)$" — this collapse happens in the $t \to 0^+$ continuum asymptotic, NOT at finite $n_{\max}$. The finite cutoff sees the discrete spectrum's rational moments; M1's $\pi$ lives at the asymptotic boundary.

---

## Verdict against gate

| Gate | Verdict |
|:-----|:--------|
| CLEAN-POSITIVE | not selected — M1 $\pi$ does NOT appear at the $s \to 0$ residue at finite cutoff. |
| **POSITIVE-WITH-STRUCTURAL-FINDING** | **selected** — Mellin transform closes bit-exactly in $\mathbb{Q}[\Gamma(s)]$; the $s \to 0$ residue is exactly $\dim \mathcal{H}$; M1 lives at the continuum boundary, not at finite cutoff. |
| PARTIAL | not selected — every formal-Mellin output is closed-form and bit-exact. |
| BLOCKED | rejected — closure is clean. |

POSITIVE-WITH-STRUCTURAL-FINDING is the truthful read of the data: the formal Mellin closes at the cutoff, AND the M1 mechanism is sharper than "the residue at $s = 0$ shows $\pi$" — it lives at the continuum-limit short-time asymptotic, an axis ORTHOGONAL to the $s$-expansion.

---

## Setup recap

From Sub-Sprint 1 (memo `debug/sprint_q5p_jlo_nmax2_memo.md`), the JLO degree-0 odd cochain at the unit on the truncated CH triple at $n_{\max} = 2$ is

$$\phi_0^{\mathrm{odd}}(1; t) \;=\; \mathrm{Tr}\!\left(e^{-t D^2}\right) \;=\; \sum_{m \ge 0} c_m t^m,$$

with bit-exact coefficients (diagonal $\Lambda$ case): $c_0 = 16$, $c_1 = -84$, $c_2 = 489/2$, $c_3 = -3967/8$, $c_4 = 98203/128$, … (sub-sprint 1 verification table).

The formal Mellin transform is

$$\mathcal{M}[\phi_0^{\mathrm{odd}}(1; t)](s) \;:=\; \int_0^\infty t^{s-1}\,\phi_0^{\mathrm{odd}}(1; t)\,dt \;=\; \int_0^\infty t^{s-1}\,\sum_i m_i e^{-t \mu_i}\,dt \;=\; \Gamma(s)\,\zeta_{D^2}(s),$$

where the second equality uses linearity over the finite spectrum $\{\mu_i\}$ of $D^2$ with multiplicities $m_i$, and the identity $\int_0^\infty t^{s-1}\,e^{-t\mu}\,dt = \Gamma(s)\,\mu^{-s}$ (valid for $\mathrm{Re}(s) > 0$, $\mu > 0$).

The finite-cutoff spectral zeta $\zeta_{D^2}^{(n_{\max})}(s) = \sum_i m_i \mu_i^{-s}$ is well-defined for all $s \in \mathbb{C}$ at finite $n_{\max}$ (no convergence issue, as the sum is finite). Its analytic continuation to $s = 0$ is trivial (the sum is finite, equals $N = \dim \mathcal{H}$).

---

## Bit-exact Mellin transform at integer $s$

### $n_{\max} = 2$, diagonal $\Lambda$ ($\Lambda^2$ spectrum: $9/4$ with mult $4$, $25/4$ with mult $12$):

| $s$ | $\Gamma(s)$ | $\zeta_{\Lambda^2}(s)$ | $\mathcal{M}[\phi_0](s) = \Gamma(s)\,\zeta_{\Lambda^2}(s)$ | continuum $\zeta_{D^2}^{\mathrm{cont}}(s)$ | trunc/cont ratio |
|:---:|:-----------:|:----------------------:|:---------------------------------------------------------:|:------------------------------------------:|:----------------:|
| 1   | 1 | $832/225$ | $832/225$ | $-\pi^2/4$ | $-1.499$ (truncated vs analytic-continued continuum; the continuum pole at $s=1$ flips sign of comparison) |
| 2   | 1 | $55552/50625$ | $55552/50625$ | $\pi^2 - \pi^4/12$ | $0.626$ |
| 3   | 2 | $4559872/11390625$ | $9119744/11390625$ | $\pi^4/3 - \pi^6/30$ | $0.945$ |
| 4   | 6 | $420155392/2562890625$ | $840310784/854296875$ | $2\pi^6/15 - 17\pi^8/1260$ | $0.991$ |
| 5   | 24 | $40725594112/576650390625$ | $325804752896/192216796875$ | $17\pi^8/315 - 31\pi^{10}/5670$ | $0.999$ |

**Bit-exact rational at every cutoff $\cdot$ integer $s$.** Trunc/cont ratio → 1 monotonically for $s \ge 3$ (consistent with the Sprint Q5'-CH-2 convergence rate $\sim n_{\max}^{3 - 2s}$).

### $n_{\max} = 2$, full $D = \Lambda + \kappa A$, $\kappa = -1/16$ (matrix-power route $\zeta_{D^2}(s) = \mathrm{Tr}((D^2)^{-s})$ at integer $s$):

| $s$ | $\Gamma(s)$ | $\zeta_{D^2}(s)$ | $\mathcal{M}[\phi_0](s)$ |
|:---:|:-----------:|:----------------:|:-----------------------:|
| 1   | 1 | $1978465196515232/532976619280625 \approx 3.7121$ | $1978465196515232/532976619280625$ |
| 2   | 1 | $316198062792242180125437283456/284064076699804288492500390625 \approx 1.1131$ | (same) |
| 3   | 2 | $62222126235359634502217877340729901078440448/151399511258533849166316043509699718994140625 \approx 0.4110$ | $2\cdot$ that |

Compare to $\Lambda$-only values (3.6978, 1.0973, 0.4003): the full-$D$ values are slightly larger (since $\kappa A$ perturbation lowers the bottom of the spectrum slightly, so inverse powers are larger). Differences are 0.4% at $s=1$, 1.4% at $s=2$, 2.7% at $s=3$, consistent with the perturbation magnitude $|\kappa| = 1/16 \approx 6.25\%$.

**Eigenvalue note (matrix-power vs direct diagonalization):** The eigenvalues of full $D^2$ at $n_{\max}=2$ are algebraic but messy — sympy returns degree-4 Galois roots involving complex cube roots (structurally consistent with the Sprint 3 RH-F finding that some GeoVac polynomials have non-solvable Galois groups; the smaller $D^2$ blocks here are degree-4 cubic-resolvent-extensions). For the Mellin transform at integer $s$, **the matrix-power route $\mathrm{Tr}((D^2)^{-s})$ bypasses this entirely**: $D^2$ has rational entries, $(D^2)^{-1}$ has rational entries (sympy `inv()`), and $(D^2)^{-s}$ for integer $s \ge 1$ remains rational. This is the right computational route at finite cutoff.

### $n_{\max} = 3$, diagonal $\Lambda$ (spectrum: $9/4$ mult $4$, $25/4$ mult $12$, $49/4$ mult $24$):

| $s$ | $\zeta_{\Lambda^2}(s)$ | $\mathcal{M}[\phi_0](s)$ | trunc/cont ratio |
|:---:|:----------------------:|:------------------------:|:----------------:|
| 1   | $62368/11025$ | $62368/11025$ | $-2.293$ |
| 2   | $152820352/121550625 \approx 1.2573$ | (same) | $0.717$ |
| 3   | $553960380928/1340095640625 \approx 0.4133$ | $1107920761856/1340095640625$ | $0.976$ |
| 4   | $\approx 0.1650$ | $\approx 0.9900$ | $0.998$ |
| 5   | $\approx 0.07065$ | $\approx 1.6952$ | $0.9999$ |

Convergence to continuum tightens monotonically with $n_{\max}$, consistent with the Sprint Q5'-CH-2 table.

---

## The $s \to 0$ residue structure (THE M1 question)

### Finite-cutoff $\zeta_{D^2}^{(n_{\max})}(s)$ at $s = 0$

For a finite spectrum $\{\mu_i\}$ with multiplicities $\{m_i\}$,

$$\zeta_{D^2}^{(n_{\max})}(s) \;=\; \sum_i m_i \mu_i^{-s} \;\xrightarrow{s \to 0}\; \sum_i m_i \;=\; N \;=\; \dim \mathcal{H}.$$

This is an **integer** (in our case $N = 16$ at $n_{\max} = 2$, $N = 40$ at $n_{\max} = 3$). **No $\pi$ appears.** The Laurent expansion of $\mathcal{M}[\phi_0](s) = \Gamma(s) \zeta_{D^2}(s)$ around $s = 0$ is

$$\mathcal{M}[\phi_0](s) \;=\; \frac{N}{s} \;+\; \Big(\!\!-\!\log\det(D^2) \,-\, N\,\gamma_E\Big) \;+\; O(s),$$

where $\gamma_E$ is the Euler–Mascheroni constant and $\log\det(D^2) = \sum_i m_i \log \mu_i$. (Derivation: $\Gamma(s) = \Gamma(s+1)/s = (1/s)(1 - \gamma_E s + O(s^2))$ and $\zeta_{D^2}(s) = N - s \log\det(D^2) + O(s^2)$ from $\mu^{-s} = e^{-s \log \mu}$.)

**At $n_{\max} = 2$, $\Lambda$ only:** $\log\det(\Lambda^2) = 4 \log(9/4) + 12 \log(25/4) = -25.2347\dots$, the next-order Laurent coefficient is $-\log\det(\Lambda^2) - 16\,\gamma_E = -12\log(25/4) - 4\log(9/4) - 16\gamma_E$.

**At $n_{\max} = 2$, full $D$:** $\det(D^2) = 110962529960861050192382965087890625 / 1208925819614629174706176$ (bit-exact integer-over-integer); $\log\det(D^2) = \log(\text{that}) \approx 25.2427$; the small $\kappa A$ perturbation barely shifts the log-det from the $\Lambda$-only value.

**At $n_{\max} = 3$, $\Lambda$:** $\log\det(\Lambda^2) = 4\log(9/4) + 12\log(25/4) + 24\log(49/4) \approx 85.367$; residue stays $N = 40$.

**No $\pi$ at any of these.** The transcendental content at finite cutoff is $\gamma_E$ (from $\Gamma$) and $\log(\text{rational})$ values — neither in M1's pure-Tate ring $\mathbb{Q}[\pi, \pi^{-1}]$ (Paper 55 Thm `thm:m1_pure_tate`), nor in M2's $\bigoplus_k \pi^{2k}\mathbb{Q}$, nor in M3's $\mathcal{MT}(\mathbb{Z}[i, 1/2])$. The transcendentals are $\gamma_E$ and $\log(\mathrm{rationals})$, which are NOT in the master Mellin engine's three period rings. **They sit outside the engine entirely at this level.**

This is the **honest scope** of the formal-Mellin closure: the transform closes in $\mathbb{Q}[\Gamma(s)]$ at integer $s$ (no transcendentals there because $\Gamma(\text{integer}) = \text{integer}$), and at $s \to 0$ the next-order coefficient introduces $\gamma_E$ and $\log$ which are master-Mellin-engine-EXTERNAL. The M1 $\pi$ appears in a DIFFERENT extraction.

### Where M1 $\pi$ actually appears: continuum Weyl-law short-time asymptotic

The M1 mechanism (Paper 18 §III.7, Paper 55 §3) is the **Hopf-base measure** $\mathrm{Vol}(S^2)/4 = \pi$, sourced from the $k = 0$ slot of the master Mellin engine, which is the **trivial heat kernel collapsing to the volume integral**. Concretely: on the continuum $S^3$, the Camporesi–Higuchi Dirac heat kernel admits the Weyl small-$t$ asymptotic

$$\mathrm{Tr}\!\left(e^{-tD^2}\right) \;\stackrel{t \to 0^+}{\sim}\; \frac{\dim(\text{spinor bundle}) \cdot \mathrm{Vol}(S^3)}{(4\pi t)^{3/2}} \;+\; \mathrm{subleading}\,,$$

with $\mathrm{Vol}(S^3) = 2\pi^2$, $\dim(\text{spinor bundle}) = 2$, hence leading coefficient $4\pi^2 / (4\pi)^{3/2} \cdot t^{-3/2} = \sqrt{\pi}/2 \cdot t^{-3/2}$ (driver `hopf_base_measure_identification()` output bit-exactly).

**This is where M1 lives in the Mellin context:** the $t \to 0^+$ short-time leading coefficient of $\mathrm{Tr}(e^{-tD^2})$, NOT the $s \to 0$ Laurent coefficient of its Mellin transform at finite cutoff. The two limits are categorically different:

- **At finite $n_{\max}$, $s \to 0$:** $\zeta_{D^2}(0) = N$ (integer, no $\pi$).
- **At $n_{\max} \to \infty$, $t \to 0^+$:** $\mathrm{Tr}(e^{-tD^2}) \sim \sqrt{\pi}/2 \cdot t^{-3/2}$ (M1 leading).

The Mellin transform AT THE CONTINUUM has a pole at $s = 3/2$ (from the $t^{-3/2}$ leading) with residue $\sqrt{\pi}/2$ (M1 coefficient). The M1 sub-mechanism enters via:

$$\mathcal{M}[\mathrm{Tr}(e^{-tD^2})]^{\mathrm{cont}}(s) \;=\; \Gamma(s) \cdot \zeta_{D^2}^{\mathrm{cont}}(s), \qquad \text{residue at } s = 3/2:\; \Gamma(3/2) \cdot \mathrm{Res}_{s=3/2}\,\zeta_{D^2}^{\mathrm{cont}}(s) \;=\; \frac{\sqrt{\pi}}{2} \cdot \frac{\mathrm{Vol}(S^3)}{(4\pi)^{3/2}}\,\Gamma(3/2) = \tfrac{1}{4}\pi.$$

(The latter equality uses $\Gamma(3/2) = \sqrt{\pi}/2$ and the Weyl-law residue $\mathrm{Vol}(S^3)/(4\pi)^{3/2} \cdot \mathrm{Res}_{s=3/2}\Gamma(s) \zeta_{D^2}(s) = $ M1 coefficient times the depth-0 factor.)

The bit-exact statement at finite cutoff: $\zeta_{D^2}^{(n_{\max})}(s)$ is **regular** at $s = 3/2$ for any finite $n_{\max}$ (it's a finite sum); the pole at $s = 3/2$ is a continuum-limit-only phenomenon. The pole RESIDUE encodes the M1 $\pi$. **The pole position $s = 3/2$ is itself the spectral dimension** ($d/2$ for $d = 3$), Paper 18 §III.6 / Paper 32 §VIII content.

### Sharper M1 reading: the regularised zeta at $s = 0$ in the continuum

In the continuum, $\zeta_{D^2}^{\mathrm{cont}}(s)$ at $s = 0$ does NOT equal $\dim \mathcal{H} = \infty$ — it has a *regularised* value via Hurwitz / analytic continuation. From the driver: $\zeta_{D^2}^{\mathrm{cont}}(0) = 0$ **exactly** (closed form, computed via Bernoulli polynomial reduction). This is the continuum CH spectral zeta's "regularised dimension" — a Bernoulli rational, NO $\pi$.

So at the $s = 0$ residue level, **neither finite-cutoff nor continuum produces M1 $\pi$ directly**. The cutoff side gives integer $N$; the continuum side gives Bernoulli rational $0$. M1 $\pi$ enters elsewhere (the pole at $s = 3/2$, or equivalently the leading $t^{-3/2}$ heat-kernel asymptotic).

This is the **structural finding**: the original task description's expectation that "the $s = 0$ value … should produce a $\pi$ coefficient" was the wrong place to look for M1. M1 lives at the *spectral-dimension pole*, not at $s = 0$. The $s = 0$ value is the regularised dimension (Bernoulli rational continuum / integer finite-cutoff); M1 $\pi$ is the residue at $s = d/2 = 3/2$.

---

## Cross-check with Paper 28 T9 and Sprint Q5'-CH-2 M2 closed forms

The continuum closed forms reproduced by the driver (via Hurwitz reduction) match Sprint Q5'-CH-2's bit-exact M2 verification at $s = 1, 2, 3, 4, 5$:

| $s$ | $\zeta_{D^2}^{\mathrm{cont}}(s)$ closed form | M2 fingerprint |
|:---:|:---------------------------------------------:|:---------------|
| 1   | $-\pi^2/4$ | $\pi^2 \cdot (-1/4)$ |
| 2   | $\pi^2 - \pi^4/12$ | $\pi^2 \cdot 1 + \pi^4 \cdot (-1/12)$ |
| 3   | $\pi^4/3 - \pi^6/30$ | $\pi^4 \cdot (1/3) + \pi^6 \cdot (-1/30)$ |
| 4   | $2\pi^6/15 - 17\pi^8/1260$ | M2 |
| 5   | $17\pi^8/315 - 31\pi^{10}/5670$ | M2 |

The finite-cutoff Mellin $\mathcal{M}[\phi_0]^{(n_{\max})}(s)$ approaches these continuum values monotonically (table above; ratios 0.626 / 0.945 / 0.991 / 0.999 at $s = 2,3,4,5$ for $n_{\max} = 2$).

**The M2 ring $\bigoplus_k \pi^{2k}\mathbb{Q}$ shows up not at the finite cutoff, but at the continuum-limit values of $\mathcal{M}[\phi_0](s)$ at integer $s$.** This is the same lesson: M2 is a continuum-asymptotic property, not a finite-cutoff property. **Both M1 and M2 enter in the continuum limit, at different extraction points (M1 at the spectral-dimension pole $s = d/2$; M2 at integer $s$).**

This sharpens Paper 32 §VIII Cor `cor:m2_mixed_tate` operationally: the cor states that $\zeta_{D^2}^{\mathrm{cont}}(s)$ at integer $s$ lives in $\bigoplus_k \pi^{2k}\mathbb{Q}$. The finite-cutoff version lives in $\mathbb{Q}$. **The transition from $\mathbb{Q}$ to $\bigoplus_k \pi^{2k}\mathbb{Q}$ happens in the continuum limit, exactly when the divergent tail of the spectrum is summed via analytic continuation.**

---

## Sub-Sprint 1 cross-check

The driver's small-$t$ Taylor expansion of $\mathrm{Tr}(e^{-t\Lambda^2})$ at $n_{\max} = 2$ matches Sub-Sprint 1's $\phi_0^{\mathrm{odd}}(1; t)$ on $\Lambda$ bit-exactly at $m = 0,1,2,3,4$ (all match=True). $c_0 = 16, c_1 = -84, c_2 = 489/2, c_3 = -3967/8, c_4 = 98203/128$. Cross-check passes.

---

## Implications for Stage 2 motivic Galois action

Three takeaways for Stage 2:

1. **The cosmic-Galois action on $\omega^{\mathrm{tri}}$ acts at the continuum limit, not at finite cutoff.** At finite $n_{\max}$, $\omega^{\mathrm{tri}}(\mathcal{T}_{n_{\max}}) \in \mathbb{Z}^3$ (integer triple, no period content; cosmic-Galois acts trivially on $\mathbb{Z}$). Non-trivial Galois action requires lifting to the continuum, where $\omega^{\mathrm{tri}}$ takes values in $M_1 \oplus M_2 \oplus M_3$.

2. **The pole at $s = d/2 = 3/2$ is the natural locus where M1 appears.** Stage 2 should target the residue Tannakian functor at $s = 3/2$ (the Connes–Marcolli analog of the QFT renormalization pole), not at $s = 0$.

3. **The integer-$s$ continuum values host M2.** $\mathcal{M}[\phi_0]^{\mathrm{cont}}(s)|_{s \in \mathbb{Z}_{\ge 1}}$ has two $\pi^{2k}$ slots per integer $s$ (Sprint Q5'-CH-2); natural target for a $U^*$ coalgebra on M2.

The three components $(M_1, M_2, M_3)$ live at three different extraction points of the same Mellin pairing:\ M1 = residue at $s = 3/2$; M2 = values at integer $s \ge 1$; M3 = values at integer $s$ on the structurally different $\eta$-Mellin pairing $\Gamma(s) \mathcal{M}[\mathrm{Tr}(\gamma D e^{-tD^2})](s)$ (CM-residue framework, Sub-Sprint 1 §"The $\eta$-style trace"). The Stage 2 cosmic-Galois action is the simultaneous compatibility of three Galois actions on three categorically different Mellin objects:\ $M_1$ by integer Tate twist, $M_2$ by even Tate twists, $M_3$ by cyclotomic-mixed-Tate at level 4.

---

## Honest scope

1. **One observable at one cutoff.** This sprint produces the formal Mellin transform at $n_{\max} \in \{2, 3\}$ for both $\Lambda$ and full $D$. Continuum limit is referenced via Sprint Q5'-CH-2 closed forms but not re-derived.

2. **Integer $s \in \{1, 2, 3, 4, 5\}$ only.** Half-integer $s$ (where the spectral-dimension pole at $s = 3/2$ lives) is not computed here. The continuum residue at $s = 3/2$ is identified structurally (Weyl-law leading coefficient $\sqrt{\pi}/2 \cdot t^{-3/2}$) but the finite-cutoff approximation at $s = 3/2$ would diverge at $n_{\max} \to \infty$ — that's the M1 mechanism content. Computing finite-cutoff $\zeta_{D^2}^{(n_{\max})}(3/2)$ at a few cutoffs and observing the $n_{\max}^{3 - 2 \cdot 3/2} = n_{\max}^0$ → log divergence pattern would be a clean follow-on.

3. **The $s \to 0$ residue is integer $N$, not $\pi$.** This is the genuine structural finding, NOT a partial failure. The task description's expectation of $\pi$ at $s = 0$ residue was based on a continuum-Weyl-law reading; at finite cutoff, the residue is the finite-spectrum dimension, an integer.

4. **Full-$D$ eigenvalues are degree-4 Galois roots (algebraic but ugly).** The matrix-power route $\mathrm{Tr}((D^2)^{-s})$ gives clean bit-exact rationals; this is the right finite-cutoff Mellin route at integer $s$. Direct diagonalization would produce the same numbers in messier form.

5. **Curve-fit-audit compliance.** All claims are direct closed-form computations:\ matrix inversion, characteristic polynomial det, Hurwitz reduction. Zero free parameters. The structural identification "M1 lives at $s = 3/2$ pole, not at $s = 0$" is forced by the Weyl law on $S^3$ (a standard result, Paper 28 T9 / Paper 51 §G1).

6. **Discrete-for-skeleton compliance.** Every finite-cutoff $\zeta$ value is bit-exact `sympy.Rational`. No PSLQ. The continuum closed forms use sympy's `zeta()` symbolic engine and are exact (Hurwitz reduction with Bernoulli polynomial values).

7. **Tag-transcendentals compliance.** When $\pi$ appears (it does, in the continuum closed forms inherited from Sprint Q5'-CH-2 — $-\pi^2/4$, $\pi^2 - \pi^4/12$, etc., AND in the M1 Weyl-law coefficient $\sqrt{\pi}/2$), it is tagged as **M2** in the integer-$s$ regular points and as **M1** in the spectral-dimension-pole residue. Paper 18 §III.7 master Mellin engine partition applies:\ M1 at $k = 0$ (Hopf-base measure) localized to the $t \to 0^+$ pole at $s = d/2$; M2 at $k = 2$ (Seeley–DeWitt) localized at integer-$s$ regular points. $\gamma_E$ at the next-order Laurent coefficient is master-Mellin-engine-EXTERNAL (it's a $\Gamma$-function constant, not in M1/M2/M3); $\log(\text{rational})$ at the next-order finite-cutoff Laurent is also master-Mellin-engine-EXTERNAL.

8. **WH1 PROVEN is not re-opened.** This sprint adds cohomological-side Mellin data on top of Sub-Sprint 1's JLO closure; does not test WH1.

9. **No paper edits applied.** Recommendations flagged below.

---

## Files

### Produced
- `debug/compute_phi0_mellin.py` — driver (~590 lines, ~1.2 s wall, bit-exact `sympy.Rational` throughout for finite-cutoff data, symbolic `sympy.Gamma` for the $\Gamma(s)$ factor).
- `debug/data/sprint_q5p_2b_phi0_mellin_data.json` — exact rational data dump (~25 KB) with finite-cutoff $\zeta_{D^2}(s)$ at integer $s$ for both $\Lambda$ and full $D$ at $n_{\max} \in \{2, 3\}$, plus continuum closed forms and convergence ratios.
- `debug/sprint_q5p_2b_phi0_mellin_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_jlo_nmax2_memo.md` and `debug/data/sprint_q5p_jlo_nmax2_data.json` — Sub-Sprint 1 source data.
- `debug/sprint_q5p_ch1_memo.md` — CH-1 leading predictions.
- `debug/sprint_q5p_ch2_memo.md` — continuum $\zeta_{D^2}(s)$ closed forms at integer $s$ (reproduced bit-exactly here).
- `geovac/spectral_triple.py` — FockSpectralTriple class.

### Published references
- Connes, A. *Noncommutative Geometry.* 1994, Ch. VI §1.
- Camporesi, R.; Higuchi, A. *J. Geom. Phys.* 20 (1996) 1–18.
- Fathizadeh, F.; Marcolli, M. *J. Geom. Phys.* 108 (2016) 1–17 (arXiv:1611.01815).
- Deligne, P. (2010); Glanois, C. *J. Number Theory* 152 (2015) 79–144.

---

## Recommended paper edits (PI to apply, decline, or modify)

### Paper 55 §subsec:open_m2_m3 — extend the Stage 1 closure with the Mellin localization

The Sub-Sprint 1 recommended paragraph (memo lines 261) closes with the symbol-level data at $n_{\max} = 2$. Suggested additional sentence at the end of that paragraph, naming the Sub-Sprint 2b finding:

> Sub-Sprint 2b (Sprint Q5'-Stage1-2b-phi0-Mellin, June 2026; memo \texttt{debug/sprint\_q5p\_2b\_phi0\_mellin\_memo.md}) further localises the M1/M2/M3 extraction points within the Mellin transform $\mathcal{M}[\phi_0^{\mathrm{odd}}(1; t)](s) = \Gamma(s) \zeta_{D^2}(s)$: M1 (Hopf-base measure $\pi$) lives at the spectral-dimension pole $s = d/2 = 3/2$ as a continuum-limit residue with leading $\sqrt{\pi}/2$; M2 (Seeley–DeWitt pure Tate) lives at integer $s \ge 1$ regular points (continuum closed forms reproduce the Sprint Q5'-CH-2 bit-exact panel verbatim); M3 (vertex-parity Hurwitz) lives in the $\eta$-style trace $\mathrm{Tr}(\gamma D e^{-tD^2})$, a structurally different Mellin pairing. At finite $n_{\max}$ the spectral zeta is bit-exact rational (e.g. $\zeta_{D^2}^{(2)}(1) = 1978465196515232/532976619280625$ for the full Dirac at $n_{\max} = 2, \kappa = -1/16$), and the $s \to 0$ residue is exactly the integer $\dim \mathcal{H}$ — no $\pi$ at finite cutoff. The transition from $\mathbb{Q}$-valued finite-cutoff data to $M_1 \cup M_2$-valued continuum data is the cosmic-Galois Stage 2 frontier.

### Paper 32 §VIII — optional Remark sharpening the master Mellin partition's extraction-point structure

Optional one-sentence Remark addition (after the existing `rem:master_mellin_domain`):

> \emph{Remark} (Mellin extraction points, Sprint Q5'-Stage1-2b, June 2026). The master Mellin partition $(M_1, M_2, M_3)$ corresponds to three categorically different extraction points of the formal Mellin transform $\Gamma(s) \zeta_{D^2}(s)$ on a compact 3-dimensional spectral triple:\ M1 = residue at the spectral-dimension pole $s = 3/2$ (Hopf-base $\sqrt{\pi}/2$ coefficient); M2 = values at integer $s \ge 1$ regular points (Seeley–DeWitt pure-Tate); M3 = analogous extraction on the $\eta$-style Mellin pairing $\Gamma(s) \mathcal{M}[\mathrm{Tr}(\gamma D e^{-tD^2})](s)$ (vertex-parity Hurwitz at quarter-integer shifts). At finite $n_{\max}$ the spectral zeta is rational-valued and master-Mellin transcendentals are absent; the M1/M2/M3 period content emerges in the $n_{\max} \to \infty$ continuum limit. See Paper 55 §subsec:open\_m2\_m3.

### Paper 18 §III.7 — optional one-sentence pointer

Optional one-sentence addition (after the existing `Mechanism-as-domain sharpening` paragraph closing): "Sprint Q5'-Stage1-2b (June 2026, memo `debug/sprint_q5p_2b_phi0_mellin_memo.md`) localises the M1 sub-mechanism's $\pi$ injection to the spectral-dimension residue $s = d/2$ of the formal Mellin transform $\Gamma(s) \zeta_{D^2}(s)$, with the finite-cutoff $\zeta_{D^2}^{(n_{\max})}$ rational-valued and the M1 $\pi$ appearing only in the continuum limit's $t^{-d/2}$ Weyl-law leading."

(All three recommendations only; no paper edits applied. PI direction sought.)

---

## One-line verdict

**POSITIVE-WITH-STRUCTURAL-FINDING.** Formal Mellin transform $\mathcal{M}[\phi_0^{\mathrm{odd}}](s) = \Gamma(s) \zeta_{D^2}(s)$ closes bit-exactly in $\mathbb{Q}[\Gamma(s)]$ at $n_{\max} \in \{2, 3\}$ for both $\Lambda$ and full $D = \Lambda + \kappa A$; $s \to 0$ residue is exactly $\dim \mathcal{H} \in \mathbb{Z}$ (no $\pi$ at finite cutoff). **M1 $\pi$ lives at the continuum spectral-dimension pole $s = d/2 = 3/2$, not at the $s = 0$ residue.** The original task's expectation of $\pi$ at $s = 0$ was the wrong extraction point; the corrected M1 location is the Weyl-law short-time pole, with leading coefficient $\sqrt{\pi}/2$ = M1 Hopf-base measure signature. Second concrete construction step of the cosmic-Galois $U^*$ Stage 1 lands.
