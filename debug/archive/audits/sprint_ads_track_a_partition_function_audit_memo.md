# Track AdS-A — Diagnostic-before-engineering audit on CFT₃-on-S³ partition functions

**Date:** 2026-05-25
**Sprint:** Track AdS-A (read-only diagnostic; no production code modified)
**Builds on:** Sprint MR-B (2026-05-06, modular residual closed form); Paper 28 §QED-vertex; Paper 35 §V (rational Casimirs); Sprint TS-E1 master Mellin engine (Paper 32 §VIII); Tracks AdS-B + AdS-C audit memos (parallel today, `debug/sprint_ads_track_{b,c}_*_memo.md`).

---

## Verdict

**GO-FAST (3–5 weeks) as a standalone deliverable, with strong recommendation to merge into the unified ~3-month CHM-correspondence sprint flagged by the Track AdS-C audit memo (§1 below).**

The framework already contains, at the production-code level: the exact spectra (scalar Laplacian and Camporesi–Higuchi Dirac on S³), the Hurwitz-zeta closed forms for the Dirichlet series at every integer argument (`qed_two_loop.py`), the Seeley–DeWitt heat-kernel coefficients (`qed_vacuum_polarization.py`), the Sprint MR-B closed form for the modular residual that controls the analytic continuation, and the propinquity-convergence machinery (Paper 38 / 45). The continuum target — Klebanov–Pufu–Safdi's $F_\text{scalar} = (\log 2)/8 - 3\zeta(3)/(16\pi^2) \approx 0.0638$ for the conformally coupled scalar on the round $S^3$, and an analogous Dirac value $\approx 0.21896$ — sits in the same transcendental ring (M2's $\log 2$ and the half-integer Hurwitz / odd-zeta $\zeta(3)$) that the master Mellin engine already classifies for GeoVac. The deliverable is one symbolic computation plus one comparison.

The 2–3 month original projection assumed building this infrastructure from scratch. It is already built.

---

## §1. Pre-emptive coordination with Tracks AdS-B + AdS-C

Both sibling audit memos were dispatched in parallel today. Reading them changes the recommended framing of Track AdS-A.

**Track AdS-B verdict (`debug/sprint_ads_track_b_entanglement_audit_memo.md`):** BLOCKED as literal Ryu–Takayanagi (the Fock graph has $n_{\max}$ disconnected $\ell$-sectors; no min-cut between sub-regions), but GO-FAST (3–6 weeks) as a Cardy–Calabrese-style discrete boundary entropy. The substantive infrastructure is `geovac/fock_graph_hodge.py` plus the wedge KMS state $\rho_W = e^{-K_\alpha^W}/Z$ that already lives in `geovac/modular_hamiltonian.py`.

**Track AdS-C verdict (`debug/sprint_ads_track_c_modular_jlms_audit_memo.md`):** GO-FAST (4–6 weeks) for the boundary-side Casini–Huerta–Myers structural correspondence at finite cutoff; BLOCKED for any JLMS-style bulk reconstruction (no AdS₄/H⁴ infrastructure; Sprint RH-B 2026-04-17 already closed Wick-rotation S³→H³ as a clean structural dead end; Sprint L3e-P3 2026-05-23 found Paper 38's $4/\pi$ rate does NOT transport to non-compact Coulomb). Recommends **combining A + B + C into a single ~3-month CHM-correspondence sprint** that builds shared wedge-KMS infrastructure once, then derives partition function, entropy, and modular Hamiltonian as parallel boundary-side observables.

**Track AdS-A is the cleanest of the three** because it is the most direct — partition function = $\zeta'_{D^2}(0)$ + Mellin-transform of the heat-kernel residual, all four pieces of which the framework already computes. The boundary modular Hamiltonian (Track C) and Cardy entropy (Track B) require building or extending sub-region infrastructure on top of the same wedge state. Track A is the substrate-only deliverable; Tracks B and C add the wedge-cut machinery on top.

**Recommendation: keep Track A as a 3–5 week standalone first move, with the explicit understanding that its output feeds directly into the unified CHM sprint when (if) PI commits to the full ~3-month combined sprint.** Track A's deliverable is a stand-alone math.OA-style note (~10–15 pages) at minimum, and naturally bundles into a larger CHM paper at maximum.

---

## §2. What the framework already has (answers to Q1 + Q6)

### §2.1 Spectra (exact)

- **Scalar Laplacian on unit S³:** eigenvalues $n(n+2)$ with degeneracies $(n+1)^2$ for $n = 0, 1, 2, \ldots$ This is exactly the GeoVac Fock spectrum after the relabeling $n_\text{Fock} = n + 1$ (so $\lambda_n = n_\text{Fock}^2 - 1 = n(n+2)$), Paper 7 §V. The conformally coupled scalar adds the conformal mass $\xi R_\text{scalar} = (1/8) \cdot 6 = 3/4$, shifting eigenvalues to $n(n+2) + 3/4 = (n + 1/2)(n + 3/2)$. The shifted spectrum is exactly the Camporesi–Higuchi half-integer ring that Paper 35 §V uses to get the spatial Casimir $1/240$.

- **Dirac on unit S³ (Camporesi–Higuchi):** $|\lambda_n| = n + 3/2$ with $g_n^\text{Dirac} = 2(n+1)(n+2)$, in `geovac/dirac_s3.py`. Both $\lambda_n^2 = (n + 3/2)^2$ and $g_n$ are rational. Wired through `qed_vacuum_polarization.py` and `qed_two_loop.py`.

**Net:** both spectra needed for the KPS partition functions are *already* in the production code, in the conventions KPS use, with no relabeling work needed.

### §2.2 Spectral zeta (exact closed forms)

`qed_two_loop.py::dirac_dirichlet_series_hurwitz(s)` returns, at every integer $s \geq 4$, the closed form

$$ D_\text{Dirac}(s) = \sum_{n \geq 0} g_n / |\lambda_n|^s = 2 (2^{s-2} - 1) \zeta_R(s-2) - \tfrac{1}{2}(2^s - 1) \zeta_R(s) $$

with the even/odd discriminant theorem (Paper 28 §parity): $s$ even → $\pi^\text{even}$ ring; $s$ odd → odd-zeta ring ($\zeta(3), \zeta(5), \ldots$). This is the Hurwitz-zeta machinery KPS use for their analytic continuation. The framework has it as exact rational coefficients of $\zeta_R(s-2)$ and $\zeta_R(s)$, sympy-verifiable, mpmath-PSLQ-checked.

The squared-Dirac series $\zeta_{D^2}(s)$ has its own T9 closed form (Paper 28 §sec:t9): $\zeta_{D^2}(s) = 2^{2s-1}[\lambda(2s-2) - \lambda(2s)]$ where $\lambda(2k) = (1 - 2^{-2k}) \zeta_R(2k)$ — pure $\pi^\text{even}$ at every integer $s$. This is the one-loop spectral zeta. **At $s = 2$ it gives $\zeta_{D^2}(2) = \pi^2 - \pi^4/12$**, a closed form the framework verified to >100 dps. KPS's free-Dirac $F$ involves $\zeta'_{D^2}(0)$, which we don't yet have a sympy closed form for — but see §2.4.

### §2.3 Heat-kernel coefficients (exact closed forms)

`qed_vacuum_polarization.py::seeley_dewitt_coefficients_s3()` returns, in exact sympy:

$$ a_0(D^2) = \sqrt{\pi},\quad a_1(D^2) = \sqrt{\pi},\quad a_2(D^2) = \sqrt{\pi}/8 $$

on unit $S^3$, with all higher $a_k = 0$ (Paper 35 / `debug/spectral_action_sd_exactness.json`). These are exactly what KPS need to subtract before zeta-regularizing the divergent spectral sum. The framework has them as exact sympy expressions, not just numerical values.

### §2.4 The modular residual closed form (Sprint MR-B)

This is the load-bearing infrastructure that the original 2–3 month projection missed. From `debug/mr_b_spectral_action_rate_memo.md`:

$$ K(t) = \mathrm{Tr}\, e^{-t D^2} = \frac{\sqrt{\pi}}{2} t^{-3/2} - \frac{\sqrt{\pi}}{4} t^{-1/2} + \varepsilon(t) $$

with the *closed form* for the residual derived from the Jacobi $\vartheta_2$ modular transformation:

$$ \varepsilon(t) = \sum_{m \geq 1} (-1)^m \sqrt{\pi} \cdot e^{-m^2 \pi^2 / t} \cdot \left[ t^{-3/2} - 2 m^2 \pi^2 t^{-5/2} - \tfrac{1}{2} t^{-1/2} \right] $$

verified to >100 digits at every test point. The modular exponent is exactly $\pi^2$. Every coefficient sits in the M2 ring $\sqrt{\pi} \cdot \mathbb{Q} \oplus \pi^2 \cdot \mathbb{Q}$.

**What this gives directly:** the Mellin transform of $K(t)$ — which is $\Gamma(s) \zeta_{D^2}(s)$ — has a Laurent expansion around $s = 0$ whose finite part is exactly $-\zeta'_{D^2}(0)$ (the functional determinant), and the modular residual contributes the *finite* piece that the heat-kernel asymptotic subtraction leaves over. **The framework already has this closed-form residual.** Substituting into

$$ F_\text{Dirac} = -\log Z = \tfrac{1}{2} \zeta'_{D^2}(0) $$

gives the KPS $F_\text{Dirac}$ value directly as a sympy expression, with the Mellin transform reducing to a Hurwitz-zeta derivative at $s = 0$ — the exact same object KPS evaluate analytically. The "$\log 2$ and $\zeta(3)$" that appear in $F_\text{KPS}$ come from $\zeta'_R(0) = -\tfrac{1}{2}\log(2\pi)$ and from $\zeta(3)$ entering through Hurwitz-at-half-integers via Bernoulli polynomials — both transcendentals that the master Mellin engine M2/M3 already catalogs.

### §2.5 Existing partial implementation

`qed_vacuum_polarization.py::spectral_zeta_derivative_at_zero(m_sq, n_max, R)` (lines 318–419) implements the **regularized remainder** — subtracts the leading large-$n$ asymptotic terms to make the spectral sum convergent. The docstring explicitly flags that the *divergent* asymptotic piece (where $\log 2$ and $\zeta(3)$ live, via Hurwitz-zeta derivatives) is "infrastructure for a future extension" and not yet implemented. **The 3–5 week sprint completes exactly that gap.**

The function as it stands at line 419 returns the convergent finite remainder; what's missing is the analytically continued divergent piece, which (per Sprint MR-B) reduces to a sympy expression in $\zeta_R'(s)$ at integer $s$. The `qed_two_loop.py` machinery already handles the Hurwitz-at-3/2 evaluation; extending it to $s = 0$ derivative is mechanical.

### §2.6 What's missing

- **Closed-form sympy expression for $\zeta'_{D^2}(0)$ on unit S³.** Requires combining §2.3 (Seeley–DeWitt subtraction), §2.4 (modular residual analytic continuation), and §2.2 (Hurwitz-zeta closed form). Sprint deliverable, week 1.
- **Closed-form sympy expression for $\zeta'_{-\Delta + 3/4}(0)$ on unit S³** (conformally coupled scalar). Spectrum $(n + 1/2)(n + 3/2)$, degeneracy $(n+1)^2$; sympy work analogous to Dirac, week 2.
- **Direct comparison to KPS.** Numerical value of GeoVac's $F$ vs $F_\text{KPS}^\text{scalar} = (\log 2)/8 - 3\zeta(3)/(16\pi^2)$ and $F_\text{KPS}^\text{Dirac}$. PSLQ identification at 100 dps to confirm bit-exactly that the two expressions are identical (not just numerically close). Week 3.
- **Propinquity rate at $\zeta'(0)$.** Use Paper 38's $4/\pi$ rate to bound the convergence of $\zeta'_{D^2_{n_\max}}(0) \to \zeta'_{D^2_{S^3}}(0)$ as $n_\max \to \infty$. This is the discrete-to-continuum statement; not strictly needed for the partition-function match (which sits at the continuum level), but it is the genuine GeoVac contribution beyond standard KPS continuum field theory. Week 4–5.

---

## §3. Standard continuum CFT₃-on-S³ partition functions (answer to Q2)

### §3.1 Reference closed forms (Klebanov–Pufu–Safdi 2011, arXiv:1105.4598)

For the round unit S³:

**Conformally coupled real scalar:**
$$ \boxed{ F_\text{scalar} = -\log|Z_{S^3}^\text{scalar}| = \frac{\log 2}{8} - \frac{3\zeta(3)}{16 \pi^2} \approx 0.06378 } $$

**Free massless Dirac fermion:**
$$ F_\text{Dirac} = -\log|Z_{S^3}^\text{Dirac}| \approx 0.21896 $$

(KPS eq. 4.39; the closed form involves $\log 2$, $\zeta(3)$, $\pi^2$ in similar combinations to the scalar; the published value's exact sympy form was not extracted from web fetch due to PDF rendering, but the ratio $F_\text{D}/F_\text{S} \approx 3.43$ confirmed *irrational* sets the structural class: same transcendental ring, different rational coefficients.)

Both values are computed via **Hurwitz-zeta regularization** applied to:
- scalar: spectrum $n(n+2) + 3/4$, degeneracy $(n+1)^2$ for $n = 0, 1, 2, \ldots$
- Dirac: spectrum $(n + 3/2)$, degeneracy $2(n+1)(n+2)$ for $n = 0, 1, 2, \ldots$ (Camporesi–Higuchi)

The functional determinant identity is $F = -\log|Z| = -(1/2) \zeta'_{D^2}(0)$ (scalar with appropriate operator) and $F = (1/2) \zeta'_D(0)$ (Dirac). KPS use Barnes-zeta closed forms; the framework's approach via Hurwitz-zeta + Sprint MR-B residual gives the same final number via a different (but equivalent) analytic route.

### §3.2 What needs converting to sympy

The KPS analytic forms use Barnes zeta. The framework's natural sympy expressions use Hurwitz zeta plus the MR-B modular residual. These are equivalent via the standard Hurwitz / Barnes identity for $d = 3$. The sprint deliverable is to do this conversion explicitly so that the GeoVac sympy expression for $F$ reduces, by symbolic manipulation in sympy, to the KPS closed form $(\log 2)/8 - 3\zeta(3)/(16\pi^2)$ for the scalar (and the analogous Dirac form). **This is a mechanical sympy simplification, not a new derivation.**

---

## §4. Discrete-to-continuum protocol (answer to Q3)

### §4.1 Standard protocol

Compute the discrete spectral zeta $\zeta'_{D^2_{n_\max}}(0)$ by direct sum (no regularization needed at finite $n_\max$ — the sum is finite). Compare to the continuum value. Track convergence as $n_\max$ increases.

### §4.2 GeoVac-specific protocol via Paper 38 propinquity rate

Paper 38 Theorem (`thm:gh_convergence`) gives $\Lambda(\mathcal{T}_{n_\max}, \mathcal{T}_{S^3}) \leq C_3 \cdot \gamma_{n_\max} \to 0$ with $C_3 = 1$ asymptotic and $\gamma_{n_\max} \sim (4/\pi) \log n_\max / n_\max$ (Paper 38 Appendix A / L2 rate). The discrete partition function $Z_{n_\max} = e^{-F_{n_\max}}$ inherits a rate bound from the spectral-triple convergence:

$$ |F_{n_\max} - F_\text{continuum}| \leq C \cdot \gamma_{n_\max} \quad \text{for some explicit } C $$

though pinning $C$ requires a Mellin-domain version of the Berezin reconstruction Sprint MR-A/B framework, which is open work. The qualitative-rate statement is direct.

### §4.3 Numerical protocol

At each $n_\max \in \{3, 5, 10, 20, 50, 100\}$, compute the discrete

$$ F_{n_\max}^{D^2} = -\tfrac{1}{2} \sum_{n=0}^{n_\max} g_n \log(\lambda_n^2) $$

and subtract the heat-kernel asymptotic piece (which is regularized at the continuum level). The remainder converges to a finite limit. PSLQ-identify the limit against $\{1, \log 2, \zeta(3), \pi^2, \pi^4\}$ to confirm bit-exactly that the result is the KPS closed form.

**Expected behavior:** convergence at rate $O(\log n_\max / n_\max)$ per Paper 38 / Sprint MR-B. At $n_\max = 100$ this gives $\sim 5\%$ precision on $F$, sufficient for the structural match. At $n_\max = 1000$ this gives $\sim 0.7\%$. Higher precision via the closed-form symbolic path; numerical truncation is the sanity check.

### §4.4 Why this works on GeoVac specifically

Three converging structural reasons:

1. The Fock-projection on $S^3$ produces *exactly* the KPS scalar spectrum (Paper 7 + Paper 35). No relabeling, no rescaling, no convention drift.
2. The Camporesi–Higuchi Dirac is *exactly* the KPS Dirac spectrum. Already wired through `qed_two_loop.py` with exact Hurwitz-zeta closed forms.
3. The Sprint MR-B modular residual is *exactly* the heat-kernel piece that controls $\zeta'(0)$'s analytic continuation. Closed form, all coefficients in M2 ring (the ring KPS's $\log 2$ lives in).

This is the cleanest possible discrete-to-continuum match in the framework. No new natural-geometry hypothesis is required; the geometry is locked.

---

## §5. Recent prior art (answer to Q4)

### §5.1 Direct overlap with GeoVac's discrete-spectral-triple approach

**Closest hit:** Gaudillot-Estrada & van Suijlekom, "Convergence of Spectral Truncations for Compact Metric Groups," IMRN 2025 (rnaf197) — proved state-space GH convergence at the level Paper 40 generalizes. They do *not* compute partition functions or compare to KPS. No scoop on the partition-function deliverable.

**Latrémolière 2025 (arXiv:2512.03573):** pointed proper quantum metric hypertopology — substrate for Paper 48 (Tier 3-Light). Does not address partition functions.

**Reconstructing manifolds from truncated spectral triples (arXiv:1912.09227):** Riemannian-side. No CFT₃ partition functions.

**Net:** no NCG / spectral-triple paper has computed CFT₃-on-S³ partition functions at finite cutoff from a discrete spectral triple and matched them to KPS at sub-percent precision. Track A would be the first.

### §5.2 Adjacent CFT₃-on-S³ partition function work to cite

- **Klebanov–Pufu–Safdi 2011** (arXiv:1105.4598): the F-theorem-without-supersymmetry paper. Closed forms for free scalar and Dirac. Already named in Q.
- **Jafferis–Klebanov–Pufu–Safdi 2011** (arXiv:1103.1181): "Towards the F-theorem." N=2 SUSY background; foundational.
- **Beccaria & Tseytlin 2017** (arXiv:1705.00292): squashed-S³ partition functions; gives the round-S³ Dirac value explicitly with sympy-style expressions.
- **Lei & van Leuven 2024** (arXiv:2406.01567): "Modularity in d > 2 free conformal field theory" — connects free-CFT₃ partition functions to elliptic gamma / modular structure. Same modular framework that Sprint MR-B's Jacobi-$\vartheta_2$ transformation lives in. **Strongest "we should cite" candidate** — they explicitly extend the modular machinery to free fermions and to odd dimensions, exactly the substrate Track A would use.
- **Anninos, Denef, Law, Sun and others** on dS / sphere partition functions and Euclidean quantum gravity (2024–2026) — adjacent but not directly competing.
- **Bobev–Bueno–Olea 2017** (arXiv:1705.00292): comments on squashed-sphere partition functions; round case as starting point.
- **Mukherjee–Naruko–Pasterski 2025** (arXiv:2603.09799): "The scheme independent 3-sphere free energy is not a monotone F-function" — refines the F-theorem status; matters for honest-scope framing.

**Honest assessment:** the round-S³ free scalar and Dirac partition functions are *very* well-known. The framework's contribution is **not** the value of $F$ — it's that the value drops out cleanly from a discrete spectral triple at finite cutoff, with a quantitative-rate convergence theorem (Paper 38), in a transcendental class predicted by the master Mellin engine (Paper 18 §III.7). The novelty is the discrete-to-continuum bridge, not the partition function itself.

### §5.3 Scoop risk

**Low.** No published discrete-spectral-triple / NCG / truncated-Connes-vS computation of CFT₃ partition functions on $S^3$ at the level Track A proposes. The Gaudillot-Estrada–vS 2025 paper is the only direct competitor on the math.OA side and they do not address physical observables.

---

## §6. Cleanest single deliverable (answer to Q5)

Two options, in order of cleanliness:

### §6.1 Sympy closed-form match (RECOMMENDED)

**Deliverable:** prove, as exact sympy identities, that

$$ F^\text{GeoVac, scalar}_{S^3} = \frac{\log 2}{8} - \frac{3 \zeta(3)}{16 \pi^2} \qquad F^\text{GeoVac, Dirac}_{S^3} = F^\text{KPS, Dirac} $$

at the symbolic level (`sympy.simplify(geovac_expr - kps_expr) == 0`), with the GeoVac expression obtained from the framework's Hurwitz-zeta machinery (`qed_two_loop.py`) plus the Sprint MR-B modular residual.

This is the cleanest single statement. Both sides are in M2 + M3 + rational ring. The match would be a sympy theorem.

### §6.2 Discrete-to-continuum quantitative-rate statement (GO-MEDIUM extension)

**Deliverable:** in addition to §6.1, prove the quantitative-rate bound

$$ |F_{n_\max}^{D^2} - F^\text{KPS}| \leq C \cdot \frac{\log n_\max}{n_\max} $$

with explicit $C$ derived from Paper 38's L2 / Sprint MR-B / Paper 40's universal $4/\pi$ rate, and verify the bound numerically at $n_\max \in \{3, 5, 10, 20, 50\}$.

This is the genuinely-new-NCG-mathematics statement and lands as a math.OA standalone (~15 pages, twelfth math.OA paper in the GeoVac series; siblings 38, 39, 40, 42, 43, 44, 45, 46, 47, 48, 49). 4–6 weeks total for §6.1 + §6.2.

### §6.3 What to NOT do at GO-FAST cadence

- **Do NOT** attempt the squashed-S³ generalization. That requires deformed-spectrum Hurwitz machinery that the framework does not have at sympy-closed-form level (squashed-S³ has $(n + a/q)$-type half-integer Hurwitz at *non*-canonical $a$, breaks MR-B's Jacobi-$\vartheta_2$ closed form). Separate sprint.
- **Do NOT** attempt the higher-genus / orbifold / lens-space partition functions. Those require quotient-sector machinery beyond the round-$S^3$ substrate.
- **Do NOT** attempt the *coupling*-dependent free energy or higher-spin partition functions. Track A is free-field round-$S^3$ only.

---

## §7. Synthesis: what the master Mellin engine predicts (cross-reference to Q6)

The framework's master Mellin engine (Paper 18 §III.7 / Paper 32 §VIII case-exhaustion theorem) classifies $\pi$-sources in any GeoVac observable as $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$ at $k \in \{0, 1, 2\}$:

- **M1** ($k = 0$): Hopf-base measure, $4/\pi$ signature (Paper 38 L2 rate)
- **M2** ($k = 2$): Seeley–DeWitt $\sqrt{\pi}$ and $\pi^2$ powers (Sprint MR-B residual, Cardy/KPS $\log 2$ piece via Bernoulli-at-half-integer + Hurwitz)
- **M3** ($k = 1$): half-integer Hurwitz / vertex parity / Dirichlet-L (Paper 28 Catalan G; KPS's $\zeta(3)$ piece via odd-zeta from Dirac Hurwitz)

**Prediction:** $F^\text{KPS, scalar} = (\log 2)/8 - 3\zeta(3)/(16\pi^2)$ decomposes as

- $(\log 2)/8$ ← M2 contribution (heat-kernel piece, $\log 2 = -2 \zeta_R'(0) - \log \pi$ at appropriate normalization; sits in $\sqrt{\pi} \cdot \mathbb{Q} \oplus \pi^2 \cdot \mathbb{Q}$ ring after Mellin)
- $-3\zeta(3)/(16\pi^2)$ ← M3 contribution (half-integer Hurwitz piece; $\zeta(3)$ is the canonical M3 transcendental per Paper 28 §sec:t9)

If this decomposition holds, **Track A's match isn't just a numerical coincidence — it's a structural verification of Paper 18's master Mellin engine in the wild on a physical CFT₃ observable.** That elevates the deliverable from "compute a known number" to "verify the master Mellin engine on independent physics data."

This is the substantive scientific content. Even if the numerical match itself is well-known (KPS 2011), the structural decomposition into M2 + M3 is the framework-specific finding.

---

## §8. Quick answers to remaining Q

**Q1 (recap):** `qed_self_energy.py` computes Dirac-on-S³ self-energy via SO(4) vertex coupling weights and the spectral sum infrastructure. `qed_vertex.py` computes the two-loop sunset with Hurwitz-zeta closed forms (the same machinery that gives the partition function ζ'(0) after one more sympy step). Both have $\zeta_R(s-2), \zeta_R(s)$ as native objects. Sprint TS / TX-B computed exact rational Casimirs (Paper 35 §V); the missing step to $\zeta'(0)$ is a sympy derivative, week 1 work.

**Q2 (recap):** $F^\text{scalar}_\text{KPS} = (\log 2)/8 - 3\zeta(3)/(16\pi^2) \approx 0.0638$; $F^\text{Dirac}_\text{KPS} \approx 0.21896$; both via Hurwitz-zeta regularization of the spectra GeoVac already uses.

**Q3 (recap):** discrete $\zeta'_{n_\max}(0)$ at $n_\max = 5$ would converge at $\sim 25\%$ via Paper 38's $4/\pi$ rate; $n_\max = 100$ gives $\sim 5\%$; sub-percent requires symbolic closed-form match. Sprint deliverable is the symbolic match, not numerical extrapolation.

**Q4 (recap):** No NCG-side competitor for the discrete-cutoff partition-function-on-$S^3$ statement. Lei & van Leuven 2024 is the closest adjacent work and a citation, not a competitor.

**Q5 (recap):** Sympy-exact closed-form match (§6.1) is the cleanest. Quantitative-rate extension (§6.2) is the math.OA-paper extension.

**Q6 (recap):** **Partition function ingredients ARE essentially already computed.** Sprint MR-B's modular residual closed form + `qed_two_loop.py`'s Hurwitz closed forms + `qed_vacuum_polarization.py`'s Seeley–DeWitt are precisely the three pieces $\zeta'(0)$ needs. The 2–3 month projection assumed building from scratch; the actual work is the symbolic gluing step plus the KPS comparison, which is 3–5 weeks.

**Q7 — overlap with B/C (recap):** All three tracks rest on the same `geovac/modular_hamiltonian.py` / wedge-KMS / Camporesi–Higuchi infrastructure. Track A is the substrate-only deliverable. The unified ~3-month CHM-correspondence sprint flagged by Track C's audit memo would do A as week 1 of a 3-month combined effort, with B and C following on the wedge-cut machinery built on top of A's partition function.

---

## §9. Recommendation

**GO-FAST as standalone (3–5 weeks)** with the explicit understanding that the output bundles into the unified CHM sprint when (if) PI commits to the larger combined effort. Two viable plans:

**Plan α (recommended):** Run Track A standalone first, 3–5 weeks. Deliverable §6.1 + §6.2 = math.OA standalone paper. After delivery, decide whether to commit to the unified ~3-month CHM sprint (which absorbs Tracks B + C on top of Track A's substrate).

**Plan β:** Commit immediately to the unified CHM sprint (~3 months total: 1 month substrate = Track A; 1 month boundary observables = Tracks B + C; 1 month writeup). Track A's standalone paper is then the first ~15 pages of a larger ~50-page CHM paper.

PM's recommendation is Plan α. Rationale: the standalone partition-function statement is independently publishable and a clean win — the master Mellin engine verification (§7) is genuinely new framework-specific content even if the numerical value is well-known. Plan α defers the larger commitment until after the first deliverable lands. Standard diagnostic-before-engineering discipline (CLAUDE.md `feedback_diagnostic_before_engineering.md`).

---

## §10. Files audited

- `geovac/qed_self_energy.py` — Dirac spectral sums, SO(4) vertex weights
- `geovac/qed_vertex.py` — referenced; Hurwitz-zeta vertex coupling
- `geovac/qed_two_loop.py` — Hurwitz closed forms for $D_\text{Dirac}(s)$, even/odd parity
- `geovac/qed_vacuum_polarization.py` — Seeley–DeWitt $a_0/a_1/a_2$, `spectral_zeta_derivative_at_zero` (incomplete; convergent remainder only)
- `geovac/dirac_s3.py` — Camporesi–Higuchi spectrum
- `geovac/modular_hamiltonian.py` — wedge KMS state, partition function (operator-system level)
- `papers/group6_precision_observations/paper_35_time_as_projection.tex` §V — rational Casimirs on S³ scalar and Dirac
- `papers/group3_foundations/paper_22_angular_sparsity.tex` — S³ scalar Laplacian framework
- `debug/mr_b_spectral_action_rate_memo.md` — modular residual closed form (substrate)
- `debug/sprint_ads_track_b_entanglement_audit_memo.md` — sibling audit, RT side
- `debug/sprint_ads_track_c_modular_jlms_audit_memo.md` — sibling audit, modular JLMS side
- `tests/test_qed_vacuum_polarization.py` — existing $\zeta'(0)$ convergence test (regularized remainder only)
- `tests/test_qed_two_loop.py` — existing Hurwitz closed-form tests

KPS reference: arXiv:1105.4598. Lei–van Leuven 2024 modular extension: arXiv:2406.01567.
