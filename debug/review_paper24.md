# Confidence Review: Paper 24 — The Bargmann–Segal Lattice: π-Free Discretization of the Harmonic Oscillator on S^5

Reviewer: Confidence Reviewer (combined CONFIDENCE_AUDITOR + CITATION_CHECKER)
Date: 2026-06-01
Pass: Wave 1 re-fire, text-level audit of `.tex` source
Source file: `papers/group3_foundations/paper_24_bargmann_segal.tex` (1368 lines)

## Calibration check

Not a calibration run; no known-honest answer was provided.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | Node count $|V(N_{\max}=5)|=56$ | Eq.(11), §III.A | **A** | EXTERNAL | Recomputed: $\sum_{N=0}^5 (N{+}1)(N{+}2)/2 = 1{+}3{+}6{+}10{+}15{+}21 = 56$. ✓ |
| 2 | Edge count $|E(N_{\max}=5)|=165$ | abstract + §III.B + §III.B Thm 1 | **A** | EXTERNAL (combinatorial enumeration) | Enumerated by hand from the stated SU(3) dipole rules ($\Delta N{=}\pm 1$, $\Delta l{=}\pm 1$, $\Delta m_l{=}0,\pm 1$): 165 directed edges. ✓ |
| 3 | Bit-exact π-free at $N_{\max}=5$ (zero irrationals encountered) | Thm 1 / §III.B | **B** | GEOVAC-ONLY for the numerical certificate; MIXED for the algebraic statement | Ran `verify_pi_free(N_max=5)` from `geovac/nuclear/bargmann_graph.py`. Output bit-exactly matches the certificate quoted in the paper (56 nodes, 165 edges, `irrationals_encountered: []`, `pi_free: True`). The algebraic statement (Thm 1 proof) rests on standard rational closed forms for the radial HO matrix elements and rank-1 Wigner $3j$ coefficients, which I verified independently (see #6 below). The numerical certificate by itself is internal consistency — GeoVac code reporting on GeoVac code. The Theorem (proof) itself is external-grade. |
| 4 | 17/17 unit tests pass | §III.B Numerical Verification | **B** | GEOVAC-ONLY | `pytest tests/test_bargmann_graph.py` → `17 passed in 0.46s`. Confirms paper's claim. Still internal-consistency only. |
| 5 | SU(3) shell degeneracy $(N{+}1)(N{+}2)/2$ matches $(N,0)$ Weyl dim | Eq.(7), §II.B | **A** | EXTERNAL | $\dim(N,0) = (N{+}1)(0{+}1)(N{+}0{+}2)/2 = (N{+}1)(N{+}2)/2$ from the SU(3) Weyl dimension formula. ✓ at $N=0,\ldots,5$: $1,3,6,10,15,21$. |
| 6 | Radial reduced matrix elements: $|\langle n_r,l{+}1|r|n_r,l\rangle|^2 = n_r+l+3/2$ and $|\langle n_r{+}1,l{-}1|r|n_r,l\rangle|^2 = n_r+1$ | Eqs.(15)–(16) | **A** | EXTERNAL | Numerical integration with `scipy` HO radial wavefunctions (16 $(n_r,l)$ pairs for each formula) reproduces both to 6 digits. Matches Moshinsky/Wybourne. |
| 7 | $B$-Bargmann measure normalization $d\mu = \pi^{-3} e^{-|z|^2} d^6 z$ integrates to 1 on $\mathbb{C}^3$ | Eq.(2) | **A** | EXTERNAL | $\pi^{-3}\cdot(\sqrt{\pi})^6 = 1$ identically. ✓ |
| 8 | HO Hamiltonian → $\hbar\omega(\hat N + 3/2)$ under Bargmann | Eq.(5) | **A** | EXTERNAL | Standard Bargmann 1961 / Hall 1994; reproduced in every textbook (e.g. Folland). |
| 9 | Borel–Weil identification of $(N,0)$ irreps with holomorphic sections of $\mathcal O(N)\to\mathbb{CP}^2$ | §II.C | **A** | EXTERNAL | Standard Borel–Weil–Bott theorem. ✓ |
| 10 | Theorem 2 (Coulomb/HO asymmetry, structural origin of calibration π) | §IV.A, Thm 2 | **B** then **C** | MIXED | The CAUSAL CHAIN (six numbered steps) is supported by classical inputs at each link (Weyl asymptotics, Borel–Weil, Folland–Stein, Paper 18). The theorem's status is openly hedged in the paper's own next paragraph: "This is a structural description rather than a stand-alone theorem in the strict mathematical sense." That self-hedge mitigates the **C** flag, but the heading "Theorem" remains stronger than the content. The HONESTY of the immediate follow-up paragraph keeps me from rating this **C**, but a reader who only skims the box gets a stronger claim than the body delivers. |
| 11 | Theorem 3 (Bargmann–Segal HO rigidity) | §V, Thm 3 | **A** for the forward direction and Step 1 (Jauch–Hill); **A** for Steps 2,3 too | MIXED (classical inputs are external, the framing as "structural dual of Paper 23 Fock rigidity" rests on Paper 23) | Each of the three classical inputs is well-established: Jauch–Hill 1940 for the biconditional uniqueness of SU(3) symmetry; Borel–Weil; Bargmann–Segal–Hall. The combination as written is sound. |
| 12 | "Five layers" Coulomb/HO asymmetry footnote | §IV footnote, line 432–449 | **E** (internal inconsistency) | GEOVAC-ONLY (internal taxonomy) | The footnote says "refined to **five layers**" then enumerates **four** items (i, ii, iii, iv), and item (iv) is labelled "(\textbf{universal})" — not part of the asymmetry. The §V.subsec "Fourth layer" then enumerates five layers (1 through 5), with Layer 4 = Pythagorean orthogonality and Layer 5 = gravity termination. The conclusion of §V at line 699 refers to "four-layer Coulomb/HO asymmetry." Three different enumerations co-exist in the paper. (See Overstatement section for the suggested fix.) |
| 13 | "Layer 4: Pythagorean orthogonality" bit-exact at every tested $(n_{\max}, N_t)$ across 18-cell six-witness panel | §V Layer 4 paragraph | **B** | GEOVAC-ONLY | This is a downstream claim about Paper 32/42/43 sprint output (referenced memos exist in `debug/`). Re-verification would require running the cited sprint compute; the paper itself flags "honest scope" — "A formal proof of the subspace decomposition at general $(n_{\max},N_t)$ is a named follow-on; the bit-exact numerical evidence at $n_{\max}\in\{1,\ldots,5\}$ … is strong evidence the orthogonality is structural, not coincidental." Self-flagged appropriately. |
| 14 | Layer 5 (NEW): spectral-action gravity termination on $S^{2m+1}$ → only $S^3$ ($m{=}1$) gives pure Einstein | §V Layer 5 paragraph | **B** | GEOVAC-ONLY (Paper 51) | This claim rests entirely on Paper 51 and the Bernoulli identity for odd polynomials at half-integer shifts. I did not independently re-derive it for this review; flagged as B (rests on a not-yet-externally-verified GeoVac paper). |
| 15 | Coulomb/HO Dirichlet series: $\sum_{N\ge 1} g_N^{HO}/N^6 = \tfrac12[\zeta(4)+3\zeta(5)+2\zeta(6)]$ "with no clean $\zeta(2)=\pi^2/6$ projection" | §V.D, line 1056–1062 | **A** | EXTERNAL | Algebra is correct: $g_N=(N^2+3N+2)/2$, so the sum splits exactly as claimed. Confirmed numerically via mpmath. ✓ The "no clean $\zeta(2)$" framing is also correct — there is no $\zeta(2)$ term in the sum. |
| 16 | $\mathbb{CP}^2$ Laplacian spectrum $\lambda_k = 4k(k+2)$ with degeneracy $(k+1)^3$ | §V.D, line 1003 | **A** | EXTERNAL | Standard result (Berger–Gauduchon–Mazet 1971); reproduced by my calculation: $\dim(k,k) = (k{+}1)^2(2k{+}2)/2 = (k{+}1)^3$. ✓ |
| 17 | Paper 2's $K=137.036064$ matches $1/\alpha$ to $8.8\times 10^{-8}$ (re-cited in §VI.A) | §VI.A | **C** at the precision claim (Paper 2 issue, re-cited here) | EXTERNAL recomputation against CODATA | $K = \pi(42+\pi^2/6-1/40) = 137.0360644…$, CODATA 2018 $\alpha^{-1} = 137.035999084$, so the relative error is **$4.77\times 10^{-7}$**, not $8.8\times 10^{-8}$. The $8.8\times 10^{-8}$ figure is the Paper 2 claim being re-stated here; the discrepancy is upstream. Paper 24 is OK in re-citing what Paper 2 claims, but if Paper 2 is wrong by ~5× the figure quoted in Paper 24's §VI.A then Paper 24 inherits the overstatement. **Recommendation: cross-check with Paper 2 audit; if confirmed there, soften here from "matches … to $8.8\times 10^{-8}$" → "matches … to within $\sim 5\times 10^{-7}$" or whatever the correct figure is.** |
| 18 | "Pi-free graph at every finite $N_{\max}$" (abstract / §III.B / conclusion) | repeated throughout | **A** | EXTERNAL (Theorem 1 proof is sound) | The proof in §III.B is solid: products of two rationals are rational, no transcendental enters at any finite $N_{\max}$. ✓ |
| 19 | $\beta_1 = 110$ (cycle classes / Wilson-loop holonomies) at $N_{\max}=5$ | §V.E (gauge-structure corollary) | **B** | GEOVAC-ONLY | $\beta_1 = E - V + c = 165 - 56 + 1 = 110$. The arithmetic checks out (computed by hand). The interpretation as "carrying Wilson-loop holonomies" is internal-framework reading. |
| 20 | The HO graph has no Wilson SU(3) structure (Corollary in §V.E) | §V.E paragraph "Reduced content (no Wilson SU(3) structure)" | **A** | EXTERNAL (representation theory) | Argument is sound: link variables would need to live in a fixed rep on every edge for a Wilson lattice gauge theory; the $(N,0)$ tower has different reps on different shells, so the dipole operator is a Clebsch–Gordan intertwiner not a group element. Confirmed by lattice gauge theory standards. ✓ |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value | Survives? |
|---|---|---|---|---|
| Node count at $N_{\max}=5$ | 56 | — (pure combinatorics) | 56 (computed in `debug/`) | YES |
| Edge count at $N_{\max}=5$ | 165 | — (pure combinatorics) | 165 (enumerated) | YES |
| Bargmann measure normalization | $\int \pi^{-3}e^{-|z|^2}d^6z = 1$ | Multivariate Gaussian | $1.000\ldots$ (machine) | YES |
| HO radial reduced ME #1 | $n_r+l+3/2$ | scipy direct integration | matches to 6 digits across 16 cases | YES |
| HO radial reduced ME #2 | $n_r+1$ | scipy direct integration | matches to 6 digits across 12 cases | YES |
| SU(3) $(N,0)$ Weyl dim | $(N+1)(N+2)/2$ | SU(3) dim formula | matches at $N=0..5$ | YES |
| $\mathbb{CP}^2$ Laplacian first 6 eigenvalues | $\{0,12,32,60,96,140,\ldots\}$ | Berger-Gauduchon-Mazet | $4k(k+2)$ matches | YES |
| $(k,k)$ SU(3) degeneracy | $(k+1)^3$ | SU(3) dim formula | matches | YES |
| Paper 2 $K=137.036064$ matches $1/\alpha$ to $8.8\times 10^{-8}$ | $8.8\times 10^{-8}$ | CODATA 2018 $\alpha^{-1}=137.035999084$ | **$4.77\times 10^{-7}$** | **NO** (off by 5×; but this is a Paper 2 issue re-cited here) |
| $\pi$-free certificate at $N_{\max}=5$ | True | `verify_pi_free` | matches exactly | YES (internal consistency) |
| 17 unit tests pass | 17 | `pytest tests/test_bargmann_graph.py` | 17 passed | YES (internal consistency) |
| HO Dirichlet series $\sum g_N/N^6$ | $\tfrac12[\zeta(4)+3\zeta(5)+2\zeta(6)]$ | mpmath direct | $3.11389631…$ matches | YES |

### Circularity map

**EXTERNAL-bottomed chains (solid):**
- Node count and edge count at $N_{\max}=5$ → pure combinatorics from the SU(3) dipole rules + the stated $(N,0)$ branching to SO(3).
- Theorem 1 (π-free certificate): proof rests on **rationality of radial HO reduced matrix elements** (external — Moshinsky/Wybourne formulas verified) + **rationality of squared rank-1 Wigner $3j$** (external — Edmonds/Varshalovich closed forms). Solid.
- Theorem 3 (HO rigidity): three external classical inputs (Jauch–Hill 1940 biconditional; Borel–Weil; Bargmann–Segal–Hall). The combination is the contribution of this paper; the inputs are not. Solid.
- Algebra of the Dirichlet series $\sum g_N/N^6$ and the $\mathbb{CP}^2$ Laplacian spectrum: external.

**GEOVAC-ONLY-bottomed chains (house-of-cards risk if upstream fails):**
- The "Layer 4: Pythagorean orthogonality" structural claim → rests on Papers 32, 42, 43 sprint outputs. Paper itself flags this honestly ("named follow-on" for formal proof).
- The "Layer 5: spectral-action gravity termination" claim → rests entirely on Paper 51 (and the Bernoulli identity for odd polynomials at half-integer shifts; that piece is external but the spectral-action machinery is GeoVac-internal).
- The "five layers" framing of the Coulomb/HO asymmetry → enumeration is internal-taxonomy, not externally established. (Internal inconsistency, see Finding 12 above.)
- 17/17 unit tests pass and `verify_pi_free` returns the certificate quoted: GeoVac testing GeoVac. Confirms the code does what the paper says it does; does NOT confirm the math is right (that's covered by the external chain above).
- Re-cited Paper 2 $K = 137.036064$ figure with claimed precision $8.8\times 10^{-8}$ → rests on Paper 2's framing. My recomputation against CODATA 2018 gives $4.77\times 10^{-7}$ instead.

**MIXED chains:**
- Theorem 2 (structural origin of calibration π) is a combination of external classical results (Weyl asymptotics, Borel–Weil, Folland–Stein) and the GeoVac-internal Paper 18 taxonomy. Paper 24 itself hedges this honestly: "This is a structural description rather than a stand-alone theorem in the strict mathematical sense."

### Overstatement findings

1. **The "five-layer" footnote (§IV footnote, line 432–449) versus the §V "Fourth layer" subsection.** Three counts of the layers co-exist in the paper:
   - Footnote: "refined to **five layers**" then enumerates four items (i, ii, iii, iv), of which (iv) is labelled "(\textbf{universal})" — not part of the asymmetry.
   - §V.subsec "Fourth layer": Layers 1–5 with Layer 4 = Pythagorean orthogonality and Layer 5 = gravity termination.
   - Line 699: "four-layer Coulomb/HO asymmetry is consistent across spectrum computation, calibration, gauge-matter coupling, and modular structure."

   **Suggested fix:** Make the footnote and §V agree. Easiest fix: change the footnote to "refined to **five layers** across Sprint 5 Track S5, Papers 30/ST-SU3, Paper 51 (gravity termination), and Sprint L2-F.1 (modular-Hamiltonian Pythagorean orthogonality)" and add Layer 5 in the footnote enumeration. Then change the §V conclusion "four-layer" → "five-layer" for consistency. CLAUDE.md §6 Paper 24 entry already lists four layers (does NOT include gravity termination as Layer 5); a corpus-level decision is needed on whether to canonize 4 or 5.

2. **Paper 2 precision figure $8.8\times 10^{-8}$.** This is re-cited from Paper 2 §VI.A. My CODATA recomputation gives $4.77\times 10^{-7}$, ~5× larger. This may be a Paper 2 issue rather than a Paper 24 issue, but Paper 24 inherits any overstatement. Recommend the Paper 2 reviewer cross-check (Paper 2 was in Wave 1 already and the team verified Paper 2 corrections; please check if this specific figure is in the corrected file).

3. **The phrase "first key structural difference" (§III.C line 419).** Mild — the paper writes "This is the first key structural difference between the Coulomb and the harmonic oscillator discretizations, and it motivates the asymmetry analysis of Section IV." Fine as written; not an overstatement, just flagging that the phrase exists.

4. **Theorem 2 labelling.** As discussed under Finding 10, the paper itself acknowledges this is "a structural description rather than a stand-alone theorem in the strict mathematical sense." The honesty mitigates the overstatement. Alternative wording: "Structural observation 2" or "Proposition 2 (structural)" would calibrate the heading to the body. LOW severity.

## Pass B — Citation and novelty

### Citation table

| `\cite` key | claimed as | verdict | what I found |
|---|---|---|---|
| `bargmann1961` | Bargmann transform, integral kernel | CITE-OK | V. Bargmann, *Comm. Pure Appl. Math.* **14**, 187 (1961). Confirmed via Wiley/scirp/Wikidata. DOI 10.1002/cpa.3160140303. |
| `hall1994` | Bargmann–Segal transform for compact Lie groups | CITE-OK | B. C. Hall, *J. Funct. Anal.* **122**, 103 (1994). Confirmed via Springer/ResearchGate; paper title and content match. |
| `segal1963` | Mathematical problems of relativistic physics; Boulder Summer Seminar 1960 | CITE-OK | I. E. Segal, AMS, 1963; based on a 1960 Boulder seminar. Confirmed via Wiley reviews + AMS catalog. |
| `fock1935` | Fock's 1935 stereographic projection of hydrogen | CITE-OK | V. Fock, *Z. Phys.* **98**, 145 (1935). Well-known landmark; spot-verified. |
| `jauchhill1940` | Biconditional uniqueness of $\mathrm{SU}(3)$ dynamical symmetry of the 3D HO | CITE-OK | J. M. Jauch and E. L. Hill, *Phys. Rev.* **57**, 641 (1940). Confirmed via APS. |
| `follandstein1974` | Estimates for the $\bar\partial_b$ complex and analysis on the Heisenberg group | CITE-OK | G. B. Folland and E. M. Stein, *Comm. Pure Appl. Math.* **27**, 429 (1974). Confirmed; DOI 10.1002/cpa.3160270403. |
| `borelweil` | Borel–Weil theorem (cited via Serre exposition + Fulton-Harris ch. 23) | CITE-OK | Standard textbook content; verifiable in any algebraic geometry / representation theory reference. The Serre 1954 Séminaire Bourbaki exposé 100 is well-known and the Fulton-Harris textbook chapter is the canonical introduction. |
| `wybourne1974` | Wybourne, *Classical Groups for Physicists*, Wiley-Interscience, 1974, Ch. 16 | CITE-OK | Standard reference book. The 1974 publication and Wiley-Interscience are correct. Did not verify Ch. 16 specifically holds the cited material, but the book itself exists and is the standard reference for SU(3) HO. |
| `iachello2006` | F. Iachello, *Lie Algebras and Applications*, Springer, 2006, Ch. 4 | CITE-OK | Real book by Iachello, Springer 2006 series Lecture Notes in Physics; standard reference. |
| `edmonds` | Edmonds, *Angular Momentum in Quantum Mechanics*, Princeton, 1957 | CITE-OK | Standard reference for closed-form $3j$ symbols. |
| `varshalovich` | Varshalovich, Moskalev, Khersonskii, *Quantum Theory of Angular Momentum*, World Scientific, 1988 | CITE-OK | Standard reference. |
| `lichnerowicz1963` | Lichnerowicz, *Spineurs harmoniques*, C. R. Acad. Sci. Paris **257**, 7–9 (1963) | CITE-OK | Confirmed via Wikipedia + multiple secondary sources. Foundational paper for the Lichnerowicz identity. |
| `moshinsky_smirnov_1996` | Moshinsky-Smirnov, *The Harmonic Oscillator in Modern Physics*, Harwood Academic, 1996 | CITE-OK | Standard textbook for HO brackets and the lab-to-relative transformation; book exists with the cited authorship and publisher. |
| `berger_gauduchon_mazet` | Berger, Gauduchon, Mazet, *Le Spectre d'une Variété Riemannienne*, LNM **194**, Springer, 1971 | CITE-OK | Confirmed via Springer (ISBN 9783540054375, 266 pages, LNM 194, 1971). The cited Ch. 5 for the $\mathbb{CP}^n$ Fubini–Study spectrum is correct. |
| `paper0`, `paper7`, `paper18`, `paper22`, `paper23`, `paper25`, `paper27_entropy`, `paper43` | Internal GeoVac papers | n/a (internal) | All exist in the corpus and are correctly listed in CLAUDE.md §6. |
| `paper2_alpha` | Paper 2 — fine-structure constant conjecture | n/a (internal) | Tagged "[Conjectural.]" in the bibitem, which is appropriate per CLAUDE.md §13.5. |
| `paper14`, `paper32`, `rowe1985` | DEFINED IN BIB BUT NEVER CITED in body | CITE-OK (entries are valid) but stale | Three unused bibitems. LOW severity cleanup. |

### Problems found

None of CITE-MISATTRIBUTED / CITE-DOESNT-SUPPORT / CITE-CANT-FIND severity.

The only citation issue is **three unused bibitems** (`paper14`, `paper32`, `rowe1985`) — cosmetic, LOW.

### Priority / novelty claims

| Claim (verbatim or close paraphrase) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "The Bargmann–Segal lattice is bit-exactly $\pi$-free at every finite $N_{\max}$" | abstract, Thm 1, conclusion | The π-free framing as a notable result is GeoVac-internal (it's specific to the GeoVac exchange-constant program). The mathematical fact that HO matrix elements are rational is itself standard. | N/A — this is a framework-internal observation, not a priority claim against the math literature. The cautious reading is: "the rational-arithmetic content of the HO matrix elements is well known; what is new in Paper 24 is the framing of this as an 'exchange constant taxonomy' result." | NO change needed. Paper is not claiming mathematical priority on the rationality. |
| "Calibration $\pi$ is Coulomb-specific … no analog for the harmonic oscillator's linear complex-analytic projection" (abstract) | abstract + Thm 2 | This framing is GeoVac-internal (Paper 18 taxonomy). The mathematical content (HO has linear spectrum, $S^3$ Laplacian has quadratic; one needs spectral-zeta regularization, the other doesn't) is well-established in spectral geometry. | N/A — framework-internal taxonomy. | NO change needed. |
| Theorem 3 (Bargmann–Segal rigidity for the 3D HO) | §V | The biconditional uniqueness of SU(3) symmetry (Step 1) is Jauch–Hill 1940. The Borel–Weil identification (Step 2) is standard. The Bargmann–Segal–Hall unitary equivalence (Step 3) is standard. The COMBINATION as the "HO analog of Fock projection rigidity" is the GeoVac-specific framing. | Not a priority claim against the math literature — it's a synthesis with classical inputs. | NO change needed; the proof correctly attributes each ingredient. |
| "First manifold-WITH-BOUNDARY carrier in the GeoVac propinquity series" — NOT in this paper, only in Paper 53. Skipped. | n/a | n/a | n/a | n/a |
| "Coulomb/HO asymmetry … sharpens Paper 2's framing" (§VI.A) | §VI.A | This is a framing-sharpening claim, not a priority claim. | n/a | NO change needed. |

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---|---|---|---|
| 1 | Internal inconsistency in layer count (footnote = "five layers" listing four, §V = five-layer enumeration, line 699 = "four-layer"). Mutually contradictory statements within Paper 24. | A | E (internal inconsistency in text) | MEDIUM |
| 2 | Paper 2 precision figure $8.8\times 10^{-8}$ re-cited from Paper 2 §VI.A. My CODATA recomputation gives $4.77\times 10^{-7}$, ~5× larger. Issue is upstream in Paper 2 but Paper 24 inherits it. | A | C (overstatement re-cited) | MEDIUM |
| 3 | Theorem 2 labelling: paper itself acknowledges in the next paragraph that this "is a structural description rather than a stand-alone theorem in the strict mathematical sense." Calibration mismatch between Theorem heading and content. | A | C (mild) | LOW |
| 4 | Three unused bibitems: `paper14`, `paper32`, `rowe1985`. | B | CITE-OK (entries valid) but stale | LOW |
| 5 | All other load-bearing citations (Bargmann 1961, Hall 1994, Jauch-Hill 1940, Folland-Stein 1974, Lichnerowicz 1963, Borel-Weil, Berger-Gauduchon-Mazet, Moshinsky-Smirnov, Edmonds, Varshalovich, Wybourne, Iachello, Segal 1963) | B | CITE-OK | none |
| 6 | π-free certificate (numerical) verifies bit-exactly via `verify_pi_free`. 17 unit tests pass. | A | B (internal consistency) — but the THEOREM 1 PROOF is external-grade A. | none |
| 7 | Node count 56, edge count 165 at $N_{\max}=5$: independently verified. | A | A | none |
| 8 | All radial HO matrix element formulas, SU(3) Weyl dimensions, $\mathbb{CP}^2$ Laplacian spectrum, Bargmann measure normalization, Dirichlet series identity: independently verified. | A | A | none |
| 9 | Theorem 3 (HO rigidity) proof is sound; the three classical inputs are correctly cited. | A | A | none |

## Broadcast readiness: **YELLOW**

Paper 24's mathematical core is solid. The π-free certificate (Theorem 1) is correctly stated and correctly proved; the rationality argument rests on standard rational closed forms (radial HO matrix elements + rank-1 $3j$) that I verified independently. The numerical certificate (56 nodes, 165 edges, zero irrationals) reproduces bit-exactly via the reference implementation. The HO rigidity theorem (Theorem 3) is correctly assembled from three established classical inputs (Jauch–Hill 1940, Borel–Weil, Bargmann–Segal–Hall). All load-bearing external citations check out.

Two MEDIUM issues prevent a GREEN verdict:

1. **The "layer count" of the Coulomb/HO asymmetry is internally inconsistent**: the §IV footnote says "five layers" but enumerates four (one of which is universal, not Coulomb-specific); the §V "Fourth layer" subsection lists five layers including a "(new)" Layer 5; line 699 of the conclusion says "four-layer." A reader cannot tell whether the asymmetry has four layers or five from reading Paper 24 alone. CLAUDE.md §6 Paper 24 entry says four layers and does not mention gravity termination as Layer 5; corpus-level decision needed.

2. **The Paper 2 precision figure $8.8\times 10^{-8}$ is re-cited but is ~5× tighter than my CODATA recomputation** ($4.77\times 10^{-7}$). This is fundamentally a Paper 2 issue, but Paper 24 §VI.A re-states the figure and inherits any overstatement. The Paper 2 audit team should reconcile, and Paper 24 should match whatever the corrected figure is.

Plus several LOW-severity items: three unused bibitems and a heading-vs-content calibration on Theorem 2.

The framework-internal layers (Layer 4 Pythagorean orthogonality, Layer 5 gravity termination) rest on GeoVac papers (32, 42, 43, 51) rather than external sources — the paper flags this honestly ("named follow-on" for formal proof of Layer 4) and the reading remains B (internally consistent, externally unverified). I am not flagging these as defects, but a reader should know they are GeoVac-internal claims, not externally established theorems.

## What I could NOT verify (hand to a human expert)

- The numerical eigenvalues of the 12-sector $m_l$-quotient Laplacian ($\{0, 2.22, 4.87, 4.94, 11.65, 12.11, 18.74, 28.09, 33.58, 50.61, 63.77, 99.43\}$). These are the output of a specific computation in `debug/data/s5_graph_spectrum.json`; I did not re-derive them from scratch.
- The closed-form structure of the Layer 4 Pythagorean residual $\|H_{\mathrm{local}} - D_W^L\|_F^2 = \kappa_g^2 \cdot S(n_{\max})/(4\pi^2) + D(n_{\max})$. Verified at the qualitative level (matches Paper 43 §10.2 closed form per CLAUDE.md), not re-derived here.
- The Layer 5 claim about spectral-action gravity termination ($S^{2m+1}$ uniquely gives $(m+1)$ power-law terms) — rests on Paper 51 and the Bernoulli identity for odd polynomials at half-integer shifts. The Bernoulli identity is external; the spectral-action embedding is GeoVac-internal.
- The Paper 2 precision claim ($8.8 \times 10^{-8}$) — should be re-checked by the Paper 2 reviewer against CODATA. My quick recomputation gave $4.77\times 10^{-7}$, but I did not exhaustively check alternative $\alpha$ source values (e.g., older CODATA or a non-CODATA recommendation).
- Whether "to our knowledge" priority claims elsewhere in the corpus (Paper 38, 45, etc.) are correct — Paper 24 itself does not make priority claims, so this is out of scope here.

---

## Summary for dispatcher

- **Paper:** 24 — The Bargmann–Segal Lattice: π-Free Discretization of the Harmonic Oscillator on $S^5$
- **Verdict:** YELLOW
- **Pass A counts:** A=10, B=4, C=2 (one re-cited from Paper 2), D=0, E=1 (internal layer-count inconsistency)
- **Pass B counts:** CITE-OK = 14 load-bearing citations (Bargmann 1961, Hall 1994, Segal 1963, Fock 1935, Jauch-Hill 1940, Folland-Stein 1974, Borel-Weil, Wybourne 1974, Iachello 2006, Edmonds, Varshalovich, Lichnerowicz 1963, Moshinsky-Smirnov 1996, Berger-Gauduchon-Mazet) + 8 GeoVac-internal `cite` keys; 0 CITE-WRONG-METADATA; 0 CITE-MISATTRIBUTED; 0 CITE-DOESNT-SUPPORT; 0 CITE-CANT-FIND; 3 unused bibitems (`paper14`, `paper32`, `rowe1985`)
- **Severity totals:** HIGH = 0; MEDIUM = 2 (layer-count inconsistency + Paper 2 precision re-citation); LOW = 2 (Theorem 2 heading-vs-content, unused bibitems)
- **Top finding:** The Coulomb/HO asymmetry "layer count" is internally inconsistent — footnote says "five layers" but enumerates four (one of which is universal, not part of the asymmetry); §V "Fourth layer" subsection lists five with a "(new)" Layer 5 (gravity termination from Paper 51); conclusion says "four-layer." Paper 24 needs a single canonical layer enumeration that agrees with CLAUDE.md §6.
