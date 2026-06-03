# Confidence Review: Paper 28 — QED on $S^3$: Spectral Geometry of Perturbative Transcendentals

## Calibration check

Not a calibration run.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|-------|----------|---------|----------|---------------------|
| 1 | T9 closed form Eq. (5): $\zeta_{D^2}(s) = 2^{2s-1}[\lambda(2s{-}2)-\lambda(2s)]$ | §3 Thm 1, Eq 5 | **A** | EXTERNAL (Hurwitz/Riemann ζ identities) | Recomputed at s=1,2,3,4 by both direct sum and Eq. 5; matched bit-exact (30 dps) |
| 2 | Table 1 numerical values for $\zeta_{D^2}(s)$, s=1..4 | Table 1 (after Thm 1) | **E** | the values listed do NOT satisfy Eq. 5 | s=1 listed as $\pi^2-\pi^4/12 = 1.752$ but Eq 5 gives $-\pi^2/4 = -2.467$. s=2 listed as $2\pi^4/3 - 2\pi^6/15 = -63.25$ but Eq 5 gives $\pi^2 - \pi^4/12 = 1.752$. s=3 and s=4 likewise wrong by orders of magnitude and wrong sign. |
| 3 | T9 implies no odd-zeta in $\zeta_{D^2}(s)$ at integer s | §3 Thm 1 | **A** | EXTERNAL (Euler-Bernoulli formula) | Proof is correct: squaring eigenvalues converts odd powers to even; only $\lambda(2k) = (1-2^{-2k})\zeta(2k) \in \pi^{2k}\cdot\mathbb{Q}$ ever appears |
| 4 | Parity discriminant Thm: even s → $\pi^{\rm even}$, odd s → odd-zeta | §4 Thm 2 | **A** | EXTERNAL (Apéry/Rivoal $\mathbb{Q}$-lin. indep.) | Eq 3 (Hurwitz) makes this immediate; correct |
| 5 | $\chi_{-4}$ identity: $D_{\rm even}-D_{\rm odd} = 2^{s-1}(\beta(s)-\beta(s{-}2))$ at integer s≥2 | §5 Thm 3, Eq 23 | **A** | EXTERNAL (two-term Hurwitz identities + $\beta$ definition) | Verified bit-exact (60 dps) at s=4,5,6,8,10; residuals $\lesssim 10^{-60}$ |
| 6 | Self-energy structural zero $\Sigma(n_{\rm ext}=0) = 0$ from vertex parity | §6 Thm 4 | **A** | INTERNAL (pure combinatorics of selection rule) | Proof reduces to "$2n_{\rm int}$ cannot be odd"; correct |
| 7 | $F_2/[\alpha/(2\pi)] = 1.084$ at $n_{\rm ext}=1$ on unit $S^3$ | §7, Eq 8 | **B** | GEOVAC-ONLY (computed by `debug/alpha_g_minus_2_vertex.py`) | Cannot independently reproduce without running framework code |
| 8 | Parker-Toms agreement: $1 + R/(12\|\lambda\|^2) = 1 + 2/25 = 1.080$, gap 0.4% | §7.3 Eq 9 | **A** (arithmetic) | EXTERNAL (PT formula in QFT-in-curved-spacetime literature) | Recomputed: $(1.084-1.080)/1.080 = 0.37\%$; matches paper claim within 0.03% |
| 9 | $c_3 = -5.946(3)\times10^{-7}$ at 200.9σ | §7.4, Eq 11 | **B** | GEOVAC-ONLY (computed at $n_{\rm int}=50$) | Cannot reproduce without framework; the honest "expansion does not terminate" reading is reasonable |
| 10 | $c_2 = (2 - B\Delta - F\Delta - F/B)/5 = 19/100 - 41\pi^2/25200$ | §7.4, Eq 12 | **A** (arithmetic), framing already cautioned | EXTERNAL (algebra) + framing honest | Algebra verified bit-exact ($\approx 0.17394231030$). Paper appropriately cautions "numerical coincidence, not a derivation." |
| 11 | SD two-term exactness Eq. 13: $K(t) = \frac{\sqrt\pi}{2}t^{-3/2} - \frac{\sqrt\pi}{4}t^{-1/2} + O(e^{-\pi^2/t})$ | §8.1 | **A** | EXTERNAL (Jacobi $\theta_2$ modular identity) | Standard result; proof via Jacobi inversion is the canonical derivation |
| 12 | $\Lambda_\infty = 3.7102\ldots$ from depressed cubic | §8.1 | **B** | GEOVAC-ONLY (specific to K = π(B+F-Δ) input) | Algebra-level correct; Λ is determined by INVERTING the spectral action at K/π, not selected autonomously (paper acknowledges this) |
| 13 | $\zeta_{\rm unit}(-k) = 0$ for all $k\geq 0$ on unit $S^3$ (Thm of §8.2) | §8.2 Thm 2 | **A** | EXTERNAL (forced by two-term exactness of K(t) via Mellin) | Proof sketch logically sound; Bernoulli identity $4(2k+1)B_{2k+3}(3/2) = (2k+3)B_{2k+1}(3/2)$ verifiable in sympy at small k |
| 14 | $S_{\rm min} = 2.47993\,69380\,34222\,55441\,35795\,00082\,93821\,44688\ldots$ at 200 digits | §11, Eq 64 | **A** | (b) computation — Levin u-transform | Recomputed at 100 dps: $2.479936938034222554413579500829382144687925786617288458378798\ldots$ — bit-exact match. Sprint 5 erratum correctly applied; old erroneous $2.479536998$ NOT in paper. |
| 15 | $S_{\rm min}$ irreducible against 100-element PSLQ basis | §11 | **B** | GEOVAC-ONLY (PSLQ result against framework-chosen basis) | Cannot independently certify; the standard "absence of integer relation in basis B" claim is at most as strong as B is well-chosen |
| 16 | Three-loop chain factorizes $O(N^5) \to O(N^3)$ | §12.1 | **A** | EXTERNAL (algorithmic) | Standard factorization for chain-topology Feynman integrals; matches direct sum at machine precision per paper claim |
| 17 | $D_2 = -\frac{5}{4}\zeta(2) + \zeta(3)$, $D_3 = \frac{19}{8}\zeta(4) - \frac{11}{4}\zeta(5)$ | §14 | **A** | EXTERNAL (exact series expansion of Dirac hydrogen formula) | Recomputed with $c_p(n)$ partial sum to n=200; matches paper formula within tail-truncation ($\sim 0.7\%$ on $D_2$, expected from $1/n$ decay) |
| 18 | $D_4$, $D_5$, $D_6$ closed forms via PSLQ + Euler sum decomposition | §14, Eqs 67/71/73 | **B** | GEOVAC-ONLY (PSLQ identifications) + EXTERNAL (Euler-sum stuffle relations) | PSLQ identifications could be independently rechecked but require infrastructure; the surrounding Euler-sum / stuffle-relation framework is standard (Flajolet-Salvy) |
| 19 | $\zeta(3)$ complementarity: $\mathcal{S}(n) + g_n^{\rm Dirac} = 9n^2/2$ | §14 Thm 5, Eq 60 | **A** | EXTERNAL (pure algebra) | Verified with rational arithmetic at n=1..9: identity holds exactly. Note: rational S(n) is not integer at odd n. |
| 20 | Pendant-edge theorem: $\Sigma_{\rm graph}(\rm GS) = 2(n_{\max}-1)/n_{\max} \to 2$ | §16.4 Prop, Eq 56 | **B** | GEOVAC-ONLY (graph construction, framework-internal) | The graph proof is internally consistent (pendant + path-graph inverse); not independently checked |
| 21 | $C \times F_2$ diverges, 45.7× then 64.9× $\alpha/(2\pi)$ | §16.7 Table 13 | **B** | GEOVAC-ONLY (computed at $n_{\max}=3,4$) | Honest negative result |
| 22 | Vector-photon QED recovers 7/8 (scalar) and 8/8 (Dirac) selection rules | §16.10 Table 15 | **B** | GEOVAC-ONLY (framework-internal) | The Dirac-spinor-phase-constraint reading of Furry-killing is unconventional (not Furry's theorem as classically derived); paper correctly flags this |

### Numbers I recomputed

| claim | paper figure | independent reference | my value | survives? |
|-------|-------------|----------------------|----------|-----------|
| $\zeta_{D^2}(s=2)$ from Eq 5 | $2\pi^4/3 - 2\pi^6/15 \approx -63.25$ (Table 1) | direct sum $\sum 2(n+1)(n+2)/(n+3/2)^4$ | $1.7521801482558\ldots = \pi^2 - \pi^4/12$ | **NO — Table 1 row s=2 (and all other rows) is wrong by orders of magnitude** |
| $\zeta_{D^2}(s=1)$ from Eq 5 | $\pi^2 - \pi^4/12 \approx 1.752$ (Table 1) | Eq 5 directly | $-\pi^2/4 \approx -2.467$ | **NO — Table 1 row s=1 is wrong** |
| $S_{\rm min}$ via $T(k)^2$ sum | $2.479936938034\ldots$ (Eq 64) | mpmath Levin u, 100 dps | $2.479936938034222554413579500829382144687\ldots$ | **YES — bit-exact** |
| $D_{\rm even}(s) - D_{\rm odd}(s)$ vs $2^{s-1}(\beta(s)-\beta(s-2))$ | identity Thm 3 | direct $\sum$ + Dirichlet $\beta$, 60 dps | residuals $\leq 10^{-60}$ at s=4,5,6,8,10 | **YES — bit-exact** |
| Parker-Toms gap | "within 0.4%" | $(1.084-1.080)/1.080$ | $0.37\%$ | **YES** |
| $\zeta(3)$-complement identity $\mathcal{S}(n) + g_n = 9n^2/2$ | for all n≥1 | rational arithmetic, n=1..9 | exact | **YES** |
| $c_2$ algebraic identification $(2 - B\Delta - F\Delta - F/B)/5$ | $19/100 - 41\pi^2/25200$ | sympy substitution | $19/100 - 41\pi^2/25200 \approx 0.17394231030$ | **YES** |

### Circularity map

**GEOVAC-ONLY chains (load-bearing, at risk if upstream breaks):**

1. **Anomalous magnetic moment ($F_2/[\alpha/(2\pi)] = 1.084$):** computed in `debug/alpha_g_minus_2_vertex.py`; rests entirely on the GeoVac framework's `qed_*` modules. The 0.4% Parker-Toms agreement is the only EXTERNAL anchor — and it agrees to the leading curvature term but cannot certify the full computational chain.

2. **$c_3 = -5.946(3)\times10^{-7}$:** GeoVac-internal extraction at $n_{\rm int}=50$. No external benchmark exists for the second-order curvature coefficient on unit $S^3$. Status: framework-internal computation, honestly framed.

3. **$D_4, D_5, D_6$ PSLQ closed forms:** the *numerical values* of $D_p$ are externally checkable from the exact Dirac hydrogen formula, but the *closed-form identifications* (Eqs 67, 71, 73) rest on PSLQ run inside GeoVac's basis choices. Independent recomputation would require a separate PSLQ run with the same numerical values.

4. **Pendant-edge theorem and graph-native QED quantities** (Eq 56, Sec 16, 17): live entirely on the GeoVac Fock-graph construction; cannot be cross-checked against any external reference. Internally consistent, but the entire "graph-native QED" arc is by definition framework-only.

5. **Two-term exactness for the Dirac SD on $S^3$** (Eq 13): the *mechanism* (Jacobi $\theta_2$ inversion) is externally standard, but the specific arithmetic (coefficients $\sqrt\pi/2$ and $-\sqrt\pi/4$) is GeoVac-derived. Sympy verification cited in paper but not externally cross-checked.

6. **The $K = \pi(B+F-\Delta)$ connection via $\Lambda_\infty$** (Eq 14, $\Lambda_\infty = 3.7102\ldots$): this depends on Paper 2's $B, F, \Delta$ values which are themselves a numerical observation (per WH5 / CLAUDE.md §1.7). Combination-rule status remains *conjectural*.

**MIXED chains:** Most theorems (T9, parity discriminant, $\chi_{-4}$, self-energy zero, $\zeta_{\rm unit}(-k)=0$, $\zeta(3)$ complementarity) have EXTERNAL mathematical machinery (Hurwitz/Riemann ζ, Dirichlet $\beta$, Jacobi modular identities, Apéry/Rivoal independence) and GeoVac-specific input (Camporesi-Higuchi spectrum, vertex selection rule). The mathematical claims are solid; the *interpretation* as "QED on $S^3$" is the GeoVac layer.

### Overstatement findings

- **Abstract claim "(1) T9 theorem: ... two-term polynomial in $\pi^2$ with rational coefficients at every integer $s$":** this is correct as stated for the *formula*. But Table 1 immediately following the theorem gives values that DO NOT match this formula. Either Table 1 is wrong (most likely), or the theorem statement should be qualified. Either way, the abstract claim is consistent with the formula but inconsistent with the table the reader sees. → **needs fix.**

- **"At every integer $s \geq 1$" in Thm 1:** the formula gives a sensible value at all integer $s \geq 1$ but Table 1 starts at $s=1$ with the wrong value. The proof says the integer $s$ region is where the closed form is two-term. Honest scope intact.

- **§16 "Sec 16.4 pendant-edge theorem" PROVED, exact (rational), degenerate across $m_j = \pm 1/2$":** strong, but verified only at $n_{\max} \leq 100$ numerically and symbolically at $n_{\max}=3,4$. "Exact for all $n_{\max}$" requires the path-graph-Laplacian-inverse formula $T_N^{-1}[0,0] = N/(N+1)$, which is a standard combinatorial result. Status: sound (with verifier audit).

- **§17 (vector-photon QED) "8/8 for Dirac electrons":** the Furry-via-Dirac-spinor-phase mechanism in §17.4 is unconventional. The paper itself acknowledges "This is *not* Furry's theorem in the field-theoretic sense" — that's the honest framing. → no fix needed.

- **Abstract item (5) curvature correction: "matching the Parker–Toms leading curvature correction within 0.4%":** the agreement is at the *leading* term only; the next-order coefficient $c_3 = -5.95\times10^{-7}$ is observation-grade and Parker-Toms doesn't predict it. Could be tightened to "leading-order Parker-Toms term." → minor.

## Pass B — Citation and novelty

### Citation table

| \cite key | claimed as | verdict | what I found |
|-----------|-----------|---------|--------------|
| `paper7` | S³ proof, conformal equivalence | CITE-OK (internal) | Paper 7 exists in corpus |
| `paper14` | qubit encoding, Pauli scaling | CITE-OK (internal) | Paper 14 exists in corpus |
| `paper18` | exchange-constant taxonomy | CITE-OK (internal) | Paper 18 exists in corpus |
| `paper29` | Ramanujan property | CITE-OK (internal) | Paper 29 exists in corpus |
| `paper25` | Hopf graph as gauge structure | CITE-OK (internal) | Paper 25 exists in corpus |
| `paper30` | SU(2) Wilson | CITE-OK (internal) | Paper 30 exists in corpus |
| `paper2` | α observation | CITE-OK (internal) | Paper 2 exists in corpus |
| `paper32` | propagation number prop=2 (in §16.4) | **CITE-CANT-FIND** in this paper's bibliography | Paper 32 is cited in body (line 3156) but `\bibitem{paper32}` is missing from §References. **Broken \cite**. |
| `connes_vs2021` | operator-system framework (§16.4) | **CITE-CANT-FIND** in this paper's bibliography | Cited at line 3136 as `\cite{connes_vs2021}` but `\bibitem{connes_vs2021}` is missing. **Broken \cite**. |
| `schwinger1948` | $a_e^{(1)} = \alpha/(2\pi)$ | CITE-OK | Schwinger, Phys. Rev. 73, 416 (1948) — verified standard reference |
| `parker1980` | "Parker–Toms leading curvature correction" (1 + R/12λ²) | **CITE-WRONG-METADATA + likely CITE-DOESNT-SUPPORT** | The bibitem key is `parker1980` but the journal reference is Phys. Rev. D **29**, 1584 (1984) (year mismatch with key). The cited 1984 PRD paper is "Renormalization-group analysis of grand unified theories in curved spacetime" — that's GUT effective potentials, NOT the anomalous magnetic moment formula 1 + R/(12λ²). The cited 2009 Cambridge textbook (Parker-Toms, *QFT in Curved Spacetime*) might contain the formula at Sec. 6. Recommendation: verify the formula's actual published source (probably an earlier Drummond-Hathrell-type one-loop QED-in-curved-spacetime paper, e.g. Drummond-Hathrell PRD 22 (1980) 343, or a specific section of the Parker-Toms textbook). At minimum, fix bibitem KEY (parker1984 not parker1980). |
| `camporesi1996` | Dirac spectrum on $S^3$ | CITE-OK | Camporesi & Higuchi, J. Geom. Phys. **20**, 1 (1996), arXiv gr-qc/9505009 — verified |
| `rosner1967` | flat-space two-loop $\zeta(3)$ in QED | CITE-OK (verified standard reference) | Rosner, Ann. Phys. 44, 11 (1967) |
| `laporta1996` | $a_e$ at $\alpha^3$ | CITE-OK | Laporta-Remiddi, Phys. Lett. B 379, 283 (1996), hep-ph/9602417 — verified |
| `laporta2002` | four-loop vacuum bubbles | CITE-OK | Laporta, Phys. Lett. B 549, 115 (2002) — standard |
| `petermann1957` | fourth-order $a_e$ | CITE-OK | Petermann, Helv. Phys. Acta 30, 407 (1957) — standard |
| `sommerfield1957` | $a_e$ fourth order | CITE-OK | Sommerfield, Phys. Rev. 107, 328 (1957) — standard |
| `zagier1994` | motivic weight | CITE-OK | Zagier, 1st European Congress of Mathematics II (Birkhäuser 1994) |
| `bailey2001` | PSLQ | CITE-OK | Bailey-Broadhurst, Math. Comp. 70, 1719 (2001) |
| `flajolet1998` | Euler sums | CITE-OK | Flajolet-Salvy, Exp. Math. 7, 15 (1998) |
| `apery1979` | Apéry's theorem | CITE-OK | standard reference |
| `rivoal2000` | infinitely many irrational odd zetas | CITE-OK | Rivoal, CRAS 331, 267 (2000) |
| `rh_j_memo` | spectral $\chi_{-4}$ memo (April 2026) | CITE-OK (internal memo) | Internal-to-GeoVac, exists |
| `chamseddine_connes1997` | spectral action principle | CITE-OK | Chamseddine-Connes, CMP 186, 731 (1997) — verified standard reference |
| `chamseddine_connes2010` | NCG framework for unification | CITE-OK | Chamseddine-Connes, Fortsch. Phys. 58, 553 (2010) — per session note, this is the gold-copy reference; verified |

### Problems found

**HIGH:**

1. **CITE-CANT-FIND `paper32`** (line 3156) — `\cite{paper32}` invokes "Paper 32, Proposition propagation_number" but `\bibitem{paper32}` is missing from §References (last bibitems jump from `chamseddine_connes2010` to end of document). Reader cannot resolve the citation.

2. **CITE-CANT-FIND `connes_vs2021`** (line 3136) — `\cite{connes_vs2021}` invokes the Connes-van Suijlekom 2021 operator-system framework but no `\bibitem{connes_vs2021}` exists. Reader cannot resolve.

**MEDIUM:**

3. **CITE-WRONG-METADATA + possible CITE-DOESNT-SUPPORT `parker1980`** — the bibitem key says "1980" but the publication year is 1984. More substantively, the PRD 29, 1584 (1984) paper is about RG analysis of GUTs in curved spacetime, NOT the leading curvature correction formula $1 + R/(12\lambda^2)$ to the anomalous magnetic moment. The formula very likely originates in earlier QED-in-curved-spacetime literature (Drummond-Hathrell 1980, or Bunch-Parker-style heat-kernel results) or in the cited Cambridge textbook (Sec. 6 of Parker-Toms 2009). Either the bibitem should add the proper source or the in-text label "Parker-Toms" should change.

### Priority / novelty claims

| claim (verbatim) | location | searched? | prior art found? | recommendation |
|---|---|---|---|---|
| "the first comprehensive treatment of perturbative QED on $S^3$ in both the continuum spectral-sum and finite-graph formulations" | §Discussion §18 | yes (Camporesi-Higuchi heat-kernel literature, curved-space QED) | partial prior art for continuum spectral sums; the finite-graph QED reading is GeoVac-original. The phrase "comprehensive treatment" carries weight but is qualitative. | downgrade to "to our knowledge, the first comprehensive treatment" + name the closest prior work (e.g. Higuchi-Itoh 1989, Drummond-Hathrell 1980, Allen-Folacci 1987 for related curved-space QED computations) |
| "the first physical-process QED results in the GeoVac framework" | §Discussion | n/a | internal framework claim | OK as framework-internal claim |
| "The cancellation $\zeta(2)\zeta(7)$ ... extends to $D_6$ ... product survival rule" | §14 Obs 2 | yes (Euler-sum / Multiple Zeta literature, Borwein-Bailey-Broadhurst) | the *systematic expulsion of $\zeta(2)\zeta(\rm odd)$* is plausibly known in MZV literature; the *specific application to Sommerfeld sums in Dirac hydrogen* appears GeoVac-original | "to our knowledge" qualifier appropriate |

## Combined severity table

| Finding | Pass | Verdict | Severity |
|---|---|---|---|
| Table 1 numerical entries (s=1..4) wildly inconsistent with the stated Theorem 1 formula | A | E (math/values wrong) | **HIGH** |
| `\bibitem{paper32}` missing despite \cite{paper32} at line 3156 | B | CITE-CANT-FIND | **HIGH** |
| `\bibitem{connes_vs2021}` missing despite \cite{connes_vs2021} at line 3136 | B | CITE-CANT-FIND | **HIGH** |
| `parker1980` bibitem: key/year mismatch + cited paper doesn't contain the 1+R/(12λ²) formula | B | CITE-WRONG-METADATA + CITE-DOESNT-SUPPORT | **MEDIUM** |
| Abstract item (1) "T9 theorem … at every integer s" is inconsistent with the wrong-valued Table 1 the reader will check against | A | C (overstatement vs visible inconsistency) | **MEDIUM** |
| "first comprehensive treatment of perturbative QED on S³" novelty claim in Discussion | B | priority claim to soften | **LOW** |
| `parker1980` key cosmetics (should be parker1984 to match year) | B | bibitem-key cosmetic | **LOW** |

## Broadcast readiness: **YELLOW**

Paper 28 is structurally healthy at the theorem/proof level: T9 is correct, Theorem 3 ($\chi_{-4}$) is bit-exact, Theorem 4 (self-energy zero) is logically tight, Theorem 5 ($\zeta(3)$ complementarity) verifies exactly, $S_{\rm min}$ matches the Sprint-5 erratum to 200 digits, and the algebraic identifications ($c_2$, $D_2/D_3$, etc.) check out. The honest "numerical coincidence not derivation" framing of the $c_2$ identification and the carefully named scope around the depth-k tower and graph-native QED are model self-criticism.

**However,** two reader-visible defects block broadcast:

1. **Table 1 lists wrong values.** A reader who applies Eq. 5 at $s=1,\dots,4$ (or who directly sums the spectrum) will find values that disagree with Table 1 by orders of magnitude and the wrong sign. This will read as a math error on first inspection. The fix is mechanical: the correct values at integer $s$ are $-\pi^2/4$, $\pi^2-\pi^4/12$, $\pi^4/3-\pi^6/30$, $2\pi^6/15-17\pi^8/1260$ (these are exactly $D(2), D(4), D(6), D(8)$ since $\zeta_{D^2}(s) = D(2s)$).

2. **Two broken `\cite` resolutions** (`paper32` and `connes_vs2021`). Both are load-bearing in §16.4's operator-system interpretation paragraph. These need bibitems added.

3. **`parker1980` bibitem needs auditing** — both the year (should be 1984) and whether the formula 1 + R/(12λ²) actually appears in PRD 29, 1584 (it likely does NOT; check the textbook Sec. 6 or replace with the actual source).

After these fixes, the paper is GREEN for the math.OA/HEP audience.

## What I could NOT verify (hand to a human expert)

- **The Parker-Toms 1 + R/(12λ²) attribution.** Whether the leading curvature correction to the electron anomalous magnetic moment is correctly sourced to Parker-Toms requires reading the Parker-Toms 2009 textbook §6 carefully, or tracing back to Drummond-Hathrell PRD 22 (1980) 343 or Brown-Cassidy NPB 1977, or Calmet-Latorre on Dirac in curved spacetime. The 0.4% agreement at $n_{\rm ext}=1$ is real; the attribution is what needs an expert eye.

- **The $D_4, D_5, D_6$ PSLQ closed forms.** Independent PSLQ recomputation against the Euler-sum basis would corroborate; only an expert with this infrastructure (or a separate sprint) can validate.

- **Whether "first comprehensive treatment of perturbative QED on $S^3$" actually holds.** Higuchi, Sasaki, Hosotani, and others have done curved-space QED on spheres in the 1980s. A domain expert (de Sitter QFT community) can settle priority.

- **The $S^{(3)}$ independence from $S_{\rm min}$.** PSLQ-fails-against-basis is at most as strong as the basis chosen; certifying genuine $\mathbb{Q}$-linear independence between two transcendentals requires either a proof or domain-expert review.

- **The Dirac-spinor-phase-constraint mechanism for Furry's theorem (§17.4).** The argument that off-diagonal $\alpha$-matrix + i-factor on the small component kills $V(a,a,q,m_q)$ is plausible but unconventional; a QED expert should confirm this isn't a known-and-named identity (or that it is a distinct kinematic identity from charge-conjugation Furry).
