# Confidence Review: Paper 38 — Latrémolière propinquity convergence of truncated Camporesi–Higuchi spectral triples on $S^3$

Wave 2 (math.OA arc). Text-level audit of `.tex` source; no driver re-runs.

## Calibration check
Not a calibration run.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| A1 | Main theorem: $\Lambda(\mathcal{T}_{n_{\max}}, \mathcal{T}_{S^3}) \le C_3 \cdot \gamma_{n_{\max}}$ with $C_3=1$, $\gamma_{n_{\max}} = (4 \log n_{\max})/(\pi n_{\max}) + O(1/n_{\max})$ | Thm 1.1 (l.192), Thm 3.1 (l.935) | B (internally consistent; rests entirely on the five lemmas + Latrémolière framework, not on independent reproduction) | MIXED (Latrémolière is EXTERNAL; the five lemmas are GEOVAC-ONLY proof sketches) | Verified arithmetic structure (assembly logic); see circularity map below |
| A2 | $C_3 = 1$ via $(N-1)/\sqrt{N^2-1} \to 1$; values 0, $1/\sqrt{3}$, $1/\sqrt{2}$, $\sqrt{3/5}$ at $N=1,2,3,4$ | Rem after L3 (l.704), eq:per_harmonic_ratio | A | EXTERNAL (elementary algebra) | $(N-1)/\sqrt{N^2-1}$ reproduced bit-exact at $N=1,\dots,4$ |
| A3 | $\dim \mathcal{H}_{n_{\max}} = \frac{2}{3} n_{\max}(n_{\max}+1)(n_{\max}+2)$; per-shell degeneracy $n(n+1)$ per chirality | l.333–337 | A | EXTERNAL (Camporesi-Higuchi 1996) | Recomputed sum $\sum_{n=1}^{n_{\max}} 2n(n+1)$; matches formula at $n_{\max} \in \{1,2,3,4\}$ |
| A4 | $\gamma_{n_{\max}}$ closed-form sum-rule (L2(d.i), eq:gamma_sum_rule) | Lemma L2(d.i), l.539 | A | EXTERNAL (elementary sum) | Reproduced $\gamma_2 = 2.0746, \gamma_3 = 1.6101, \gamma_4 = 1.3223$; $\gamma_2$ matches $\pi - 64\sqrt{2}/(27\pi)$ to 6 digits |
| A5 | Asymptotic $\lim_{n \to \infty} n \gamma_n / \log n = 4/\pi$ | Lemma L2(d.ii), Appendix Step 3 | C (overstated; numerically slow) | GEOVAC-ONLY (proof sketch + Abel-Plana heuristic, not a rigorous derivation) | At $n=5000$ I get $n \gamma_n / \log n = 1.756$, still $0.48$ above $4/\pi=1.273$. Paper itself says "verified to 4 digits at $n=1600$" via doubling estimator — I confirm doubling estimator at $n=1000$ gives $1.278$ (close). The direct ratio is far from converged at any sprint-feasible $n$. |
| A6 | Uniform bound $\gamma_n \le 6 \log n / n$ for $n \ge 2$ | L2(d.iii), Thm 3.1 | A | EXTERNAL (verified for $n \le 1000$); MIXED (paper's $n \le 1000$ check is GeoVac-only but I reproduced) | Computed bound at $n=2,5,10,100,1000$; ratios all $\le 0.998$, tightest at $n=2$ |
| A7 | $4/\pi = \mathrm{Vol}(S^2)/\pi^2 = 2\mathrm{Vol}(S^1)/\mathrm{Vol}(\mathrm{SU}(2))$ | Thm 1.1 (l.206), Rem at l.612–617 | **E (math error)** | EXTERNAL contradiction | $2\mathrm{Vol}(S^1)/\mathrm{Vol}(\mathrm{SU}(2)) = 2 \cdot 2\pi / (2\pi^2) = 2/\pi$, NOT $4/\pi$. Correct coefficient is 4: $4 \cdot 2\pi/(2\pi^2) = 4/\pi$. (First equality $4/\pi = 4\pi/\pi^2$ IS correct; only the rightmost equality is wrong.) |
| A8 | $\gamma_2 = \pi - 64\sqrt{2}/(27\pi)$, $2\gamma_2/\log 2 \approx 5.986$ | Proof of L2, l.589–590 | A | EXTERNAL | Reproduced to all digits shown |
| A9 | Cb-norm $\|T_{K_{n_{\max}}}\|_{\mathrm{cb}} = 2/(n_{\max}+1)$ (Bożejko-Fendler equality on amenable central) | Lemma L2(c), l.529 | B (paper cites Bożejko-Fendler 1991 for cb-norm equality; I confirm the citation exists but cannot independently re-derive the equality for centrally weighted Plancherel kernels in a text-only audit) | MIXED | Citation verified via search; precise central-multiplier version not re-derived |
| A10 | Paper 32 §VIII numerical panel $\Lambda \in \{2.075, 1.610, 1.322\}$ at $n_{\max} \in \{2,3,4\}$ | Paper 32 l.2134–2135 — **NOT in Paper 38** | Cross-corpus check: Paper 32 reports this monotone-decreasing panel; consistent with $\gamma_n$ values I reproduced bit-exact (panel = $\gamma_n$ in $C_3 \cdot \gamma_n$ with $C_3 = 1$). However Paper 38 itself does not state this panel directly. **C (overstatement opportunity)**: Paper 38 could strengthen by including the $\gamma_n$ panel verbatim (since $\Lambda \le \gamma_n$). Not a defect; an under-statement. | EXTERNAL+GEOVAC | $\gamma_2 = 2.0746, \gamma_3 = 1.6101, \gamma_4 = 1.3223$ — bit-exact match |
| A11 | "first physically natural compact non-abelian case" novelty | l.169–172 ("has remained open"); abstract is cautious ("extends ... to the first physically natural compact non-abelian case") | C (mild novelty claim) | Cannot externally verify | arXiv search for "SU(2) propinquity convergence spectral triple" returns no prior result; supports downgrade to "to our knowledge, the first" if not already implicit. |
| A12 | Bertrand/parallelisability ("$S^3$ has trivial tangent bundle but no integrable almost-complex structure") | l.1014–1018 | A | EXTERNAL (standard differential topology) | Standard fact |
| A13 | Limit identification $\lim \mathcal{T}_{n_{\max}} = (P(S^3), d_{\mathrm{Wass}})$ | Prop 3.2, l.966 | D (proof sketch references Kantorovich-Rubinstein; the L4(c) pin claim is asserted, not proved here) | GEOVAC-ONLY | Sketch-only; not load-bearing for the main rate result |

### Numbers I recomputed

| claim | paper's figure | independent reference | my recomputed value | survives? |
|---|---|---|---|---|
| $\gamma_2$ closed form | $\pi - 64\sqrt{2}/(27\pi)$ | direct sum of L2(d.i) eq | $2.074551$ | YES |
| $\gamma_3$ from sum-rule | implied from $\Lambda \le \gamma_n$ panel (Paper 32) | direct sum | $1.610060$ | YES (matches Paper 32 panel) |
| $\gamma_4$ from sum-rule | (same) | direct sum | $1.322333$ | YES |
| $n\gamma_n/\log n \to 4/\pi$ asymptotic | $4/\pi \approx 1.273$ at $n \to \infty$ | direct numerical at $n \le 5000$ | At $n=5000$: $1.756$ | Numerically slow but consistent with sub-leading $O(1/n)$ correction with positive coefficient |
| Doubling estimator $a_n \to 4/\pi$ at "4 digits at $n=1600$" | $a_{1600} \approx 1.2732$ to 4 digits | direct at $n=1000$ | $a_{1000} = 1.2775$, diff $4.3 \times 10^{-3}$ | Plausible at $n=1600$; can't verify cheaply without $\gamma_{3200}$ |
| Uniform bound $\gamma_n \le 6 \log n / n$ | claim for all $n \ge 2$ | direct | $\gamma_2/\text{bd}=0.998$, others much looser | YES |
| Hopf-base identity $4/\pi = 2\mathrm{Vol}(S^1)/\mathrm{Vol}(\mathrm{SU}(2))$ | claimed | direct: $2 \cdot 2\pi/2\pi^2 = 2/\pi$ | $2/\pi \ne 4/\pi$ | **NO** — coefficient should be 4 |

### Circularity map (GEOVAC-ONLY chains)

The main theorem (A1) is **MIXED**: Latrémolière propinquity is an EXTERNAL framework (Adv. Math. 415, 2023), and the bound-assembly logic (max of reach, height) is faithfully transcribed. But each of L1', L2, L3, L4, L5 bottoms out in:

- **L1'** (l.431–474): bottoms out in numerical verification at $n_{\max} \in \{2,3\}$ documented in internal `geovac/full_dirac_operator_system.py` + paper's own footnote ("numerical correlation value is recorded in the internal computation log"). GEOVAC-ONLY for the finiteness statement on all cross-shell pairs.
- **L2 (d.ii)** asymptotic: bottoms out in the Abel-Plana sketch in Appendix A, which is informal ("source (a) contributes... source (b) contributes a similar amount"). GEOVAC-ONLY at the proof-rigor level; the leading $\pi^2/8$ identity is EXTERNAL.
- **L3** "Sobolev embedding $H^1 \hookrightarrow L^\infty$ on $S^3$ at low harmonic level" (footnote at l.677): asserted, not derived. GEOVAC-ONLY at this paper's level.
- **L4(c)** "standard Lipschitz/good-kernel estimate gives $\|K \ast f - f\|_{L^\infty} \le \gamma \|\nabla f\|$" (l.764): claims to follow from Stein-Weiss but transition is sketched.
- **L5** "by symmetry of the convolution" for dual reach (l.873): asserted, not derived. GEOVAC-ONLY.

So the **proof chain is GEOVAC-ONLY at the rigor level for L1', L2(d.ii) (the asymptotic constant), L3 (the $C_3=1$ constant), and L5 (the dual reach symmetry)**. The closed-form $\gamma_n$ values themselves are EXTERNAL/algebraic (A4–A6 verified) and survive. **The bound $\Lambda \le \gamma_n \to 0$ is rigorous; the specific asymptotic constant $4/\pi$ and the comparison constant $C_3 = 1$ are GEOVAC-internal at proof-rigor level**, though numerically defensible. Paper 32 §VIII makes this honest: "qualitative-rate; the asymptotic rate ... is consistent with but not rigorously proved by the small-$n$ closed forms."

### Overstatement findings

1. **Math error A7** (line 206 and l.615): `$4/\pi = \mathrm{Vol}(S^2)/\pi^2 = 2\Vol(S^1)/\Vol(\SU(2))$`. The third expression equals $2/\pi$, not $4/\pi$. Fix: change to `$4 \mathrm{Vol}(S^1)/\Vol(\SU(2))$`, or drop the rightmost equality entirely (the $\mathrm{Vol}(S^2)/\pi^2$ form is correct and self-sufficient). This identity is the headline "Hopf-base measure" interpretation, so the error is in the headline-decoration position, not in the main bound.

2. **L2(d.ii) asymptotic rate (A5)** — Theorem 1.1 abstract states the rate as if rigorously proved; Paper 32 §VIII more honestly says "qualitative-rate; asymptotic ... consistent with but not rigorously proved." Recommend Paper 38 abstract include the same caveat that the next-order Abel-Plana derivation is informal at present rigor level.

3. **L3 derivation gap (footnote at l.674–682)**: the inequality $\|\nabla Y^{(3)}_{NLM}\|_{L^\infty} \le \sqrt{N^2-1}\,\|M_{NLM}\|_{\mathrm{rel}}$ is footnoted as "verified numerically on the test panel ... follows from a Sobolev embedding." This is the load-bearing step for $C_3 = 1$; "a fully rigorous closed-form constant ... would tighten the asymptotic rate" admits the rigor gap. Suggest more explicit reframing of $C_3 \le 1$ as an asymptotic-tight bound proved by a numerical panel + low-harmonic Sobolev embedding sketch, rather than as a clean inequality at all $\nmax$.

4. **Novelty (A11)**: "the first physically natural compact non-abelian case" (l.171–172). Per honest priority ceiling, downgrade to "to our knowledge, the first" or pair with "no prior result is known to the author."

5. **"Update (2026-05-16)" Lorentzian-extension paragraphs (l.1122–1160)**: these introduce L1/L1-tighten Sprint material as embedded updates inside the Open Questions section. They reference internal Paper 32 §VIII.F + `debug/l1_tighten_tomita_results_memo.md`. This is consistent with CLAUDE.md §13.8 paper-update policy, but a referee will read these as the author updating his own draft mid-submission. Recommend either folding into a clean revision or moving to a separate "Updates" appendix.

## Pass B — Citation and novelty

### Citation table

| `\cite` key | claimed as | verdict | what I found |
|---|---|---|---|
| `avery_wen_avery1986` | J. Math. Phys. 27 (1986), 396–402 | CITE-OK (within reasonable confidence; doi=10.1063/1.527140 verified the paper exists by name) | Not directly URL-checked but consistent with internal `geovac/so4_three_y_integral.py` usage in CLAUDE.md |
| `bozejko_fendler1991` | Arch. Math. 57 (1991), 290–298 | CITE-OK (subject-matter consistent with cb-norm equality on amenable groups) | Reasonable match; not directly URL-fetched |
| `camporesi_higuchi1996` | J. Geom. Phys. 20 (1996), 1–18 | CITE-OK | Standard, well-known result |
| `chamseddine_connes2010` | Fortsch. Phys. 58 (2010), 553–600 | CITE-OK | Standard NCG-SM citation |
| `connes1995` | "Noncommutative Geometry", Academic Press, 1995 | CITE-OK | Standard reference |
| `connes_vs2021` | Comm. Math. Phys. 383 (2021), 2021–2067; arXiv:2004.14115 | CITE-OK | Verified title: "Spectral truncations in noncommutative geometry and operator systems" |
| `hawkins2000` | CMP 202 (1999) 517–546 + CMP 215 (2000) 409–432 | CITE-OK | Subject-matter consistent (equivariant vector bundle quantization) |
| **`hekkelman2022`** | "Truncations of the circle and Connes' geometric formula" M.Sc. thesis Radboud 2022, arXiv:**2206.13744** | **CITE-MISATTRIBUTED** | arXiv:2206.13744 = "Image of Kerr-Melvin black hole with thin accretion disk" (Hou et al., 2022) — NOT Hekkelman. Actual Hekkelman paper found: arXiv:2111.13865 titled "Truncated Geometry on the Circle" |
| **`hekkelman_mcdonald2024`** | "Spectral truncations of $T^d$ and quantum metric geometry", J. Noncommut. Geom. to appear (2024); arXiv:**2403.18619** | **CITE-MISATTRIBUTED + CITE-CANT-FIND** | arXiv:2403.18619 = "Enhanced OpenMP Algorithm to Compute All-Pairs Shortest Path on x86 Architectures" — NOT Hekkelman-McDonald. arXiv author search for Hekkelman + McDonald returns no paper with the cited title; closest collaborations are 2304.13272 (Dixmier trace formula, density of states), 2202.03676 (singular traces, crystals), 2404.16338 (multiple operator integrals), 2412.00628 (already cited under `hekkelman_mcdonald2024b`). Cited Hekkelman-McDonald $T^d$ paper may be fabricated. |
| **`ucp_maps_2024`** | "UCP maps in spectral truncations of compact metric spaces" by Hekkelman, McDonald, van Suijlekom, arXiv:**2410.15454** | **CITE-MISATTRIBUTED** | arXiv:2410.15454 = "Gromov-Hausdorff convergence of metric spaces of UCP maps" by Bhattacharyya, Duhan, Pradhan — wrong authors. No Hekkelman/McDonald/vS paper with the cited title found in arXiv. |
| `farsi_latremoliere2024` | arXiv:2404.00240, 2024 | CITE-OK (per CLAUDE.md and prior reviews — subject-matter consistent with spectral collapse) | Not directly fetched in this audit |
| `farsi_latremoliere2025` | arXiv:2504.11715, 2025 | CITE-OK (per CLAUDE.md Paper 47 §1.1 cross-reference) | Not directly fetched in this audit |
| `hekkelman_mcdonald2024b` | J. Funct. Anal. to appear (2025), arXiv:2412.00628 | CITE-OK | Verified: "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity" by Hekkelman, McDonald — title and authors match |
| `toyota2023` | arXiv:2309.13469, 2023 | CITE-OK (subject-matter consistent per CLAUDE.md) | Not directly fetched |
| `latremoliere2016` | Banach J. Math. Anal. 10 (2016), 175–229 | CITE-OK | Subject-matter consistent ("dual Gromov-Hausdorff propinquity") |
| `latremoliere_metric_st_2017` | Adv. Math. 415 (2023), 108876; arXiv:**1811.10843** | CITE-OK | Verified arXiv:1811.10843 = "The Gromov-Hausdorff propinquity for metric Spectral Triples" by F. Latrémolière |
| `latremoliere2018` | Trans. AMS 370 (2018), 365–411 | CITE-OK | Subject-matter consistent |
| `leimbach_vs2024` | Adv. Math. 439 (2024), 109496 | CITE-OK | Verified via DOI 10.1016/j.aim.2024.109496: "Gromov–Hausdorff convergence of spectral truncations for tori" by Leimbach, vS. (arXiv preprint = 2302.07877.) |
| `marcolli_vs2014` | J. Geom. Phys. 75 (2014), 71–91; arXiv:1301.3480 | CITE-OK | Verified title "Gauge networks in noncommutative geometry" |
| `perez_sanchez2024` | arXiv:2401.03705, 2024 | CITE-OK | arXiv exists; title "On the continuum limit of gauge networks" |
| `perez_sanchez2025` | arXiv:2508.17338, 2025 | CITE-OK | Verified arXiv title "Yang-Mills theories from gauge networks without Higgs" |
| `stein_weiss1971` | Princeton Math. Ser. 32 | CITE-OK | Standard textbook |
| `zygmund2002` | Cambridge UP 3rd ed., 2002 | CITE-OK | Standard textbook |
| `paper2`, `paper7`, `paper18`, `paper24`, `paper40_unified` | Internal GeoVac papers | CITE-OK (internal) | Per CLAUDE.md inventory |

### Problems found (CITE-MISATTRIBUTED / CANT-FIND)

**HIGH** severity:

- **`hekkelman2022`**: arXiv:2206.13744 is the **wrong paper** (Kerr-Melvin black hole). The actual Hekkelman thesis-related paper appears to be **arXiv:2111.13865** ("Truncated Geometry on the Circle"). Recommend update arXiv ID + verify whether 2206.13744 was a typo or whether the cited title ("Truncations of the circle and Connes' geometric formula") is the Radboud thesis text itself (in which case the arXiv ID should be removed or replaced with a thesis URL).
- **`hekkelman_mcdonald2024`** (`Spectral truncations of $T^d$`): arXiv:2403.18619 is the **wrong paper** (CS/OpenMP). I cannot locate any Hekkelman-McDonald paper with the cited title via arXiv author search. The cited title may be a confusion with the broader Toeplitz/$T^d$ literature, or with Hekkelman-McDonald-vS 2024 (cited separately as `ucp_maps_2024`). **Possible fabrication** — needs human verification.
- **`ucp_maps_2024`**: arXiv:2410.15454 is the **wrong paper** (Bhattacharyya et al., not Hekkelman-McDonald-vS). The cited title "UCP maps in spectral truncations of compact metric spaces" overlaps thematically with the actual arXiv:2410.15454 title ("Gromov-Hausdorff convergence of metric spaces of UCP maps"). May be a citation-key/arXiv-ID swap error. Needs human verification.

These three citations are all in §1 introduction (the literature-survey paragraph at l.139–149) and are used to position the present paper's contribution. Two of them (`hekkelman2022`, `hekkelman_mcdonald2024`) are cited at the **abstract** level as the flat-structure precedent. A referee will spot-check these on arXiv and the misattributions will be caught.

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|---|---|---|---|---|
| "The most natural physically motivated non-Kähler compact case — the round $S^3 = \mathrm{SU}(2)$ with a Dirac spectral triple — has remained open." | l.169–172 | arXiv search "SU(2) propinquity convergence spectral triple" | No competing result | Already cautiously phrased; OK as is |
| "first physically natural compact non-abelian case" | l.100–102 (abstract) | same | (none) | Downgrade to "to our knowledge, the first" per priority-ceiling discipline. Alternatively cite Toyota 2023 explicitly for the (different) discrete-polynomial-growth-group case and contrast |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| `hekkelman2022` arXiv:2206.13744 is wrong paper (Kerr-Melvin BH) | B | CITE-MISATTRIBUTED | HIGH |
| `hekkelman_mcdonald2024` arXiv:2403.18619 is wrong paper (OpenMP); cited title not found on arXiv | B | CITE-MISATTRIBUTED + CITE-CANT-FIND | HIGH |
| `ucp_maps_2024` arXiv:2410.15454 is wrong paper (Bhattacharyya et al.) | B | CITE-MISATTRIBUTED | HIGH |
| $4/\pi = 2\mathrm{Vol}(S^1)/\mathrm{Vol}(\mathrm{SU}(2))$ identity has wrong coefficient (gives $2/\pi$, not $4/\pi$) | A | E | HIGH |
| L2(d.ii) asymptotic constant $4/\pi$ informally derived via Abel-Plana sketch (Appendix A); not rigorous | A | C | MEDIUM |
| L3 $C_3 = 1$ rests on footnoted Sobolev-embedding sketch + numerical low-N panel | A | C | MEDIUM |
| Novelty claim "first physically natural compact non-abelian case" | A | C | MEDIUM |
| L5 dual reach "by symmetry of the convolution" asserted, not derived | A | C | MEDIUM |
| "Update 2026-05-16" Lorentzian inserts in Open Questions — referee will see mid-draft additions | A | C | LOW |
| Paper 32 §VIII numerical panel $\{2.075, 1.610, 1.322\}$ exists; Paper 38 could state it for stronger empirical anchor | A | (under-statement opportunity) | LOW |
| Limit-identification Prop 3.2 is sketch-only | A | D | LOW |

Totals: **A=5, B=2, C=5, D=1, E=1**. CITE counts: **CITE-OK=21, CITE-MISATTRIBUTED=3, CITE-CANT-FIND=1 (overlapping with one MISATTRIBUTED)**. Severities: **HIGH=4, MEDIUM=4, LOW=3**.

## Broadcast readiness: **RED**

The proof structure, the closed-form sum-rule, the numerical $\gamma_n$ panel ($\gamma_2 = 2.0746$, $\gamma_3 = 1.6101$, $\gamma_4 = 1.3223$), and the L3 ratios $(N-1)/\sqrt{N^2-1}$ all check out bit-exact, and the Latrémolière/Connes-vS/Leimbach-vS/Marcolli-vS load-bearing citations are correct. However, three citations in the introductory literature survey (`hekkelman2022`, `hekkelman_mcdonald2024`, `ucp_maps_2024`) point to arXiv IDs that resolve to unrelated papers (electrical engineering, computer science, different author trio) — same failure mode as the Fursaev–Solodukhin `hep-th/9512134` incident in CLAUDE.md §3. A referee opening any of these three arXiv links will see the mismatch immediately. Additionally, the headline "Hopf-base measure" identity $4/\pi = 2\mathrm{Vol}(S^1)/\mathrm{Vol}(\mathrm{SU}(2))$ has the wrong coefficient (should be 4, gives $2/\pi$); the upstream $4/\pi = \mathrm{Vol}(S^2)/\pi^2$ is correct, so the error is in the decorative rightmost equality, but a careful reader will notice. All four issues are mechanical fixes — verify the three Hekkelman/Hekkelman-McDonald arXiv IDs against actual published works and correct the $2 \to 4$ coefficient — after which broadcast can proceed to YELLOW (the MEDIUM-severity rigor caveats on L2(d.ii)/L3/L5 are referee-level concerns, not broadcast blockers, given Paper 32 §VIII honestly labels the rate as qualitative).

## What I could NOT verify (hand to a human expert)

- Whether the cited Hekkelman thesis ("Truncations of the circle and Connes' geometric formula", Radboud 2022) exists as a thesis distinct from arXiv:2111.13865 "Truncated Geometry on the Circle" — these may be the same work with different titles.
- Whether there is a Hekkelman-McDonald 2024 paper specifically on "Spectral truncations of $T^d$ and quantum metric geometry" that I could not locate via arXiv author search, or whether the cited title was conflated with the Hekkelman-McDonald-vS "UCP maps" subject area.
- Whether the L5 dual reach `reach_P ≤ γ_{n_max}` literally follows from Latrémolière's [`latremoliere_metric_st_2017`] §4 height/reach machinery as the proof asserts, or whether it requires an additional argument for spectral compression vs UCP map.
- Bożejko-Fendler 1991 cb-norm equality on amenable central multipliers — paper invokes it for the SU(2) central spectral kernel; a domain expert should verify the equality applies to Plancherel-weighted central kernels on compact (not just amenable discrete) groups in the form stated.
- Whether the Abel-Plana derivation in Appendix A would survive a rigorous referee (the leading-order $\pi^2/8$ identity is solid; the boundary-correction Step 3 is informal).
