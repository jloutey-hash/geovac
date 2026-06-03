# Confidence Review: Paper 1 — The Geometric Atom: Quantum Mechanics as a Packing Problem

**Reviewer:** GeoVac Confidence Reviewer (combined CONFIDENCE_AUDITOR + CITATION_CHECKER)
**Mode:** Wave 1 re-fire, text-level only (.tex source); no computational drivers re-run beyond formula-level verification
**Paper file:** `papers/group3_foundations/paper_1_spectrum.tex` (dated Feb 22, 2026)
**Date of review:** 2026-06-01

---

## Calibration check

Not a calibration run. The "Berry phase retracted v1.2.0" note in CLAUDE.md §6 alerted me to look for residual unretracted Berry-phase claims; cross-referencing `debug/qa_sprint/berry_phase_reconciliation.md` (the 2026-03-15 reconciliation memo) and `debug/berry_phase_convergence/results.md` confirms the retraction wording in §III.A of the current Paper 1 text is consistent with the reconciliation memo's required actions.

---

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|---|---|---|---|---|
| 1 | "Transition operators $T_\pm, L_\pm$ reproduce the exact Rydberg spectrum $E_n = -1/(2n^2)$ from operator eigenvalues" | Abstract; §II.A | A | EXTERNAL (Barut-Kleinert 1967; well-known SO(4,2) result) | Standard result, externally verified |
| 2 | "Berry phase $\theta = 0$ identically for plaquettes built from $T_\pm, L_\pm$ (real positive Biedenharn-Louck CGs)" | §III.A; Eq. (8); Eqs. (5)–(7) | A | EXTERNAL (SU(2)+SU(1,1) CG coefficients are real positive in the BL convention) | Recomputed: $\sqrt{(n-l)(n+l+1)/4} \cdot \sqrt{(l-m)(l+m+1)} \cdot \sqrt{(n+1-l)(n+1+l)/4} \cdot \sqrt{(l+m+1)(l-m)} > 0$ for any valid plaquette; `arg(positive real) = 0` |
| 3 | "Erratum: $k = 2.113 \pm 0.015$ exponent retracted; replaced by $\theta = 0$ identically" | §III.A erratum paragraph | A | GEOVAC-ONLY for prior text; EXTERNAL for the retraction | `debug/qa_sprint/berry_phase_reconciliation.md` documents the issue (paper placeholder, no reproducing code); the current paper text properly retracts |
| 4 | Log-holonomy $\Theta(n) = -2\ln((n+1)/n)$ with power-law exponent $k=1.0$, $R^2=1.0$ | §III.B; Eq. (10); Appendix A.B | A | EXTERNAL (elementary log identity); MIXED if "k=1.0 exactly" is read as a numerical fit | Symbolic: $\ln(1/(n+1)^2) - \ln(1/n^2) = -2\ln((n+1)/n)$ identically (sympy verified residual = 0); asymptote is $-2/n + 1/n^2 + O(n^{-3})$ (NB: sign — see overstatement finding below) |
| 5 | "At $n_{\max} = 10$ (385 nodes): $\lambda_{2s} = 0.0202$, $\lambda_{2p} = 0.0238$, $\Delta E = 0.0035$, 16% relative" | §III.C; Eqs. (12)-(14); abstract; Table I; Fig. 2 caption | B | GEOVAC-ONLY (computed from internal Lanczos diagonalization) | Per the instruction "text-level only, do NOT re-run computational drivers" I did not re-execute the eigenvalue solve. Per debug/qa_sprint/berry_phase_reconciliation.md and Section IV.A flow, the 385-node count is consistent with $\sum_{n=1}^{10} n^2 = 385$. Values not externally verified. |
| 6 | "Convergence sweep $n_{\max}=5\to30$: 13% → 16% (peak) → 0.3% → 0.005%; oscillatory then $\to <0.01\%$" | §III.D | B | GEOVAC-ONLY | Same caveat as #5. The qualitative claim "spectral aliasing decays in continuum limit" is structurally correct for a graph Laplacian truncation; the specific numbers were not externally verified. |
| 7 | "$D_{(2,0,0)} = 0.854$, $D_{(2,1,0)} = 3.416$, ratio $\approx 4\times$" | §III.B; Table I | B | GEOVAC-ONLY | Not re-verified (text-only review). |
| 8 | "Operators commute exactly $[T_+, L_+] = 0$" | §IV.C item 2 | A | EXTERNAL ($SU(2)$ and $SU(1,1)$ act on disjoint quantum numbers in this basis; their generators commute) | Algebraically: $L_\pm$ raises $m$ at fixed $(n,l)$; $T_\pm$ raises $n$ at fixed $(l,m)$. The two raises act on independent labels in the basis ordering used here, hence $[T_+, L_+] = 0$ as operators on the tensor product. |
| 9 | "Gearing ratio $\|L_+\|/\|T_+\| \to \approx 1.77$" | §IV.C item 1 | B | GEOVAC-ONLY | Not re-verified. |
| 10 | "$\alpha$ cannot emerge from a single-particle electron lattice; requires photon coupling" | Abstract; §IV.C; §VI | A | EXTERNAL (standard physics: $\alpha = e^2/(4\pi\epsilon_0\hbar c)$ measures e-γ coupling) | Standard, externally true. The paper's framing is honest and explicit. |
| 11 | "We conjecture that $\alpha$ emerges as the geometric impedance — the ratio of action densities between the electron paraboloid and the photon fiber" | §VI "Outlook" paragraph | D | GEOVAC-ONLY (conjecture, properly labeled) | The paper labels this a conjecture explored in the companion paper. Properly hedged. |
| 12 | "$\Delta E_{\text{Lamb}} \sim \alpha^5 \cdot mc^2$" | §IV.C, scaling list | C | EXTERNAL (Bethe scaling is $\alpha^5 mc^2 / n^3 \cdot \ln(1/\alpha) \cdot Z^4$ for hydrogen; the $\alpha^5$ leading exponent is correct, but it is conventionally accompanied by a Bethe log factor) | This is a conventional informal scaling list; the omission of the Bethe log is standard textbook practice for an order-of-magnitude scaling list. Borderline but acceptable. No fix recommended. |

### Numbers I recomputed

| claim | paper's figure | independent reference | my recomputed value/error | survives? |
|---|---|---|---|---|
| $\theta = 0$ for all plaquettes | $\theta = 0$ identically | Direct evaluation at $(n,l,m)=(2,1,-1)$: $T_+ = 1.0$, $L_+ = \sqrt{2}$, $T_- = \sqrt{2}$, $L_- = \sqrt{2}$; product $= 2\sqrt{2} > 0$; $\arg = 0$ | bit-exact match | YES |
| Log-holonomy formula | $\Theta(n) = -2\ln((n+1)/n)$ | sympy symbolic simplification of $\ln(1/(n+1)^2) - \ln(1/n^2)$ | residual = 0 (exact identity) | YES |
| Leading asymptote | $\Theta \sim 2/n$ (Eq. 10 RHS, unsigned) | sympy series at $n \to \infty$: $-2/n + 1/n^2 - 2/(3n^3) + \ldots$ | sign mismatch: should be $-2/n$, not $+2/n$ | Numerical magnitude survives; sign in Eq. (10) is wrong (LOW severity, see Overstatement) |
| 2s/2p split at $n_{\max}=10$ | $\Delta E = 0.0035$, 16% relative | NOT recomputed (text-only review) | n/a | UNVERIFIED — flagged as GEOVAC-ONLY (verdict B) |
| Per-plaquette $k$ convergence to 1.0 | "$k=1.0$ exactly, $R^2=1.0$" | `debug/berry_phase_convergence/results.md` shows finite-$n_{\max}$ fits give $k \in \{0.899, 0.928, 0.955, 0.970, 0.976, 0.980\}$ for $n_{\max} \in \{10,20,50,100,150,200\}$; convergence is monotonic to 1 from below | The "exactly" qualifier is correct in the analytical / continuum limit but slightly overstated for the finite-$n_{\max}$ numerical fit; the analytical formula's exponent is 1 exactly | YES (interpreted as analytical statement); MEDIUM-LOW finding (see Overstatement) |

### Circularity map

**GEOVAC-ONLY chains identified:**

- **Quantitative 2s/2p splitting numbers (claims 5, 6, 7):** computed from `geovac/lattice.py` + `geovac/hamiltonian.py` Lanczos diagonalization. The structural claim "splitting is a discretization artifact that vanishes in the continuum limit" is externally supportable (standard graph-Laplacian truncation behaviour); the specific numerical values are internal. **No corroborating external reference is cited for these numerical values, and none is plausible — they are specific to this lattice construction.** Per CLAUDE.md §6, Paper 1 is the originating spectrum paper; no other paper independently regenerates these values. RISK: low for the qualitative claim (well-understood phenomenon); medium for the specific numerical values being reproduced exactly by a third party without the code.

- **Conjecture: $\alpha$ as geometric impedance (claim 11):** rests on the existence of a companion paper (`\cite{companion_alpha}`). The companion paper exists (Paper 2, `papers/group5_qed_gauge/paper_2_alpha.tex`) but uses an entirely different mechanism (cubic equation in $\alpha$, Hopf-fibration spectral invariants $B+F-\Delta$). The "geometric impedance" framing in Paper 1 §VI does NOT match what Paper 2 actually does. See CITE-WRONG-METADATA finding in Pass B.

**EXTERNAL anchors (solid):**

- Fock 1935 SO(4) dynamical symmetry (claim 1's foundation)
- Barut-Kleinert 1967 SO(4,2) noncompact algebra and Biedenharn-Louck CGs (claims 1, 2, 8)
- Berry 1984 phase definition $\theta = \arg(\prod \cdots)$ (claim 2's framework)
- Chung 1997 graph Laplacian formalism $L = D - A$ (claim 3's framework)
- Standard textbook QED scaling laws for $\alpha$ corrections (claim 10)

### Overstatement findings

1. **(MEDIUM, §III.B Eq. 10):** The exact log-holonomy `$\Theta(n) = -2\ln((n+1)/n) \sim 2/n$` mixes a NEGATIVE quantity $-2\ln((n+1)/n) < 0$ with a POSITIVE asymptote `~ 2/n`. The correct asymptote is `~ -2/n`. The downstream `$\Theta \propto n^{-1}$` framing is correct (a scaling exponent). Suggested fix: change Eq. (10)'s rightmost expression to `$\sim -\frac{2}{n}$` for consistency with the leading minus sign. (Note: Appendix A.B repeats `$\Theta(n) = -2\ln((n+1)/n)$` correctly but does not re-state the asymptote.)

2. **(LOW, abstract):** "...the edge-weight log-holonomy $\Theta(n) = -2\ln((n+1)/n) \sim n^{-1}$ measures the curvature of the weight function across the lattice." This is correct as a scaling statement. No fix.

3. **(LOW, §III.B):** "The per-plaquette log-holonomy decays as $n^{-1}$ exactly ($k = 1.0$, $R^2 = 1.0$)." The analytical asymptote is exactly $n^{-1}$, but `debug/berry_phase_convergence/results.md` Table 1 shows the *fitted* exponent for finite $n_{\max}$ is sub-unity ($0.899, 0.928, \ldots, 0.980$ for $n_{\max}=10,\ldots,200$), converging to 1 only in the limit. The current text reads as if a fit gave $k=1.0$ with $R^2=1.0$. Suggested fix: rewrite as "The leading-order asymptote is $n^{-1}$ exactly (the next-order correction is $+n^{-2}$); a power-law fit to the analytical formula across $n \in [1, 200]$ recovers $k \to 1$ monotonically from below, consistent with the $+n^{-2}$ subleading term." This is also LOW because the analytical statement is unambiguous.

4. **(MEDIUM, §III.C and abstract):** The abstract says "spectral aliasing from lattice truncation that oscillates and decays to $<0.01\%$." This is a strong empirical claim (the 0.005% at $n_{\max}=30$ value) presented as established. The 2s/2p split numbers in §III.C and §III.D are NOT externally cross-checked (no independent code re-runs the sweep). They rest only on internal GeoVac code. Suggested softening: add a one-line caveat in §III.D that the convergence sweep is computed within the GeoVac framework and consistent with general graph-Laplacian truncation behaviour, but the specific numerical values $\{0.13, 0.16, 0.003, 0.00005\}$ are framework-internal. Not blocking; the qualitative structural claim is well-founded.

5. **(MEDIUM, §V.A "wave-like" vs "particle-like" language):** The paper invokes "wave-particle duality" as a heuristic ("This duality is suggestive: the algebraic operators provide exact spectral content, while the graph topology introduces finite-size corrections..."). The conclusion paragraph already states "The algebra-geometry duality is suggestive of a connection to wave-particle complementarity, though we do not claim it establishes one." This is properly hedged — accept as-is.

6. **(LOW, §VI "Outlook"):** The "geometric impedance" $\alpha$ conjecture mismatches what Paper 2 actually does (Paper 2 is Hopf-fibration spectral cubic, NOT a helical photon gauge fiber with symplectic ratio). The body honestly labels it a conjecture in a companion paper, so the framing is hedged — but a reader who follows the citation will find a paper with a different mechanism. See CITE-WRONG-METADATA in Pass B.

---

## Pass B — Citation and novelty

### Citation table

| \cite key | claimed as | verdict | what I found (URL / arXiv / DOI) |
|---|---|---|---|
| barut1967 | "$SO(4,2)$ \cite{barut1967,fock1935}" — dynamical symmetry group of hydrogen | CITE-OK | Barut & Kleinert, Phys. Rev. **156**, 1541 (1967), DOI [10.1103/PhysRev.156.1541](https://link.aps.org/doi/10.1103/PhysRev.156.1541). Title "Transition Probabilities of the Hydrogen Atom from Noncompact Dynamical Groups" matches. |
| fock1935 | "$SO(4,2)$ \cite{barut1967,fock1935}" — dynamical symmetry group of hydrogen | CITE-OK | V. Fock, Z. Phys. **98**, 145 (1935), "Zur Theorie des Wasserstoffatoms". [NASA ADS](https://ui.adsabs.harvard.edu/abs/1935ZPhy...98..145F), [Springer Link](https://link.springer.com/article/10.1007/BF01336904). Note: Fock's original paper actually establishes $SO(4)$ symmetry (the rotation group of a 4-sphere in momentum space), not $SO(4,2)$. $SO(4,2)$ is the larger dynamical (spectrum-generating) algebra added later (Barut and others). The grouped citation "$SO(4,2)$ \cite{barut1967,fock1935}" reads naturally as "we use $SO(4,2)$ following these two foundational references," and is a standard citation pattern. CITE-OK with the minor sharpening note. |
| biedenharn1981 | "Biedenharn-Louck calculus for $SU(2)$ and $SU(1,1)$ Clebsch-Gordan coefficients" | CITE-OK | Biedenharn & Louck, *Angular Momentum in Quantum Physics*, Encyclopedia of Math. and Its Applications **8** (Addison-Wesley, 1981). Cambridge later republication 1984; original Addison-Wesley publication 1981 [Amazon DE](https://www.amazon.de/Angular-Momentum-Quantum-Physics-Encyclopedia/dp/0521302285), [Google Books](https://books.google.com/books/about/Angular_Momentum_in_Quantum_Physics.html?id=8HWQwdXqcXAC). Carruthers is the series editor, not a co-author of this volume (some catalog listings list all three). Cited as covering $SU(2)$ CGs (standard) and $SU(1,1)$ CGs (also covered by this reference series). CITE-OK. |
| berry1984 | Berry phase definition $\theta = \arg(\prod \cdot)$ in Appendix A | CITE-OK | Berry, Proc. R. Soc. Lond. A **392**, 45 (1984), DOI [10.1098/rspa.1984.0023](https://royalsocietypublishing.org/doi/10.1098/rspa.1984.0023). The cited use (defining geometric phase from holonomy of adiabatic eigenstates around a closed loop) is exactly what the Berry 1984 paper establishes. |
| chung1997 | Spectral graph theory framework $L = D - A$ | CITE-OK | F. Chung, *Spectral Graph Theory*, CBMS Regional Conference Series in Mathematics **92** (AMS, 1997). Verified via [AMS](https://www.ams.org/books/cbms/092/). |
| companion_alpha | "...we conjecture that $\alpha$ emerges as the geometric impedance — the ratio of action densities between the electron paraboloid and the photon fiber. This hypothesis is explored in a companion paper, which introduces a helical photon gauge fiber and computes the symplectic coupling ratio." | CITE-WRONG-METADATA / CITE-DOESNT-SUPPORT | The bibitem text reads "The Geometric Atom: Deriving the Fine Structure Constant from Lattice Helicity (companion paper, 2026)." A paper with that title does NOT exist in the corpus. The companion $\alpha$-paper that does exist is `papers/group5_qed_gauge/paper_2_alpha.tex`, titled "The Fine Structure Constant from Spectral Geometry of the Hopf Fibration." Paper 2's mechanism is a cubic in $\alpha$ with coefficient $K = \pi(B + F - \Delta)$ built from Hopf-fibration spectral invariants — NOT a "helical photon gauge fiber" / "symplectic coupling ratio" / "geometric impedance" framework. Paper 1's §VI Outlook framing and bibitem title both anticipate a Paper 2 manuscript that was reorganized in subsequent sprints. The cited work exists under a different title and pursues a different mechanism — citation should be updated. |
| paper13 | "...adiabatic Berry connection from eigenstates of a parameter-dependent Hamiltonian (as in the hyperspherical $P_{\mu\nu}(R)$ connection of Ref. \cite{paper13})" | CITE-OK | Paper 13, "Hyperspherical Lattice: Two-Electron Atoms as Coupled Channel Graphs," `papers/group2_quantum_chemistry/paper_13_hyperspherical.tex`. The Hellmann-Feynman / nonadiabatic coupling matrix $P(R)$ is in fact the central technical object of Paper 13. CITE-OK. |

### Problems found

**CITE-WRONG-METADATA / CITE-DOESNT-SUPPORT — `companion_alpha`:**

Paper 1 line 397 lists the companion as:
> J. Loutey, "The Geometric Atom: Deriving the Fine Structure Constant from Lattice Helicity," (companion paper, 2026).

The actual companion exists as Paper 2: "The Fine Structure Constant from Spectral Geometry of the Hopf Fibration" (`papers/group5_qed_gauge/paper_2_alpha.tex`). The mechanism described in Paper 1 §VI (helical photon gauge fiber, symplectic impedance, action-density ratio) is NOT what Paper 2 implements (Paper 2's mechanism is the cubic $\alpha^3 - K\alpha + 1 = 0$ with $K = \pi(B + F - \Delta)$ from Hopf-fibration spectral invariants). Furthermore, CLAUDE.md §13.5 protects the "conjectural" framing of Paper 2's combination rule, and CLAUDE.md §1.7 WH5 records Paper 2 has been moved to Observations status (2026-05-02). Recommended fix at MEDIUM severity:

1. Update the bibitem title to match the actual Paper 2 title.
2. Either (a) update §VI "Outlook" paragraph to describe Paper 2's actual mechanism (Hopf-fibration spectral cubic), or (b) drop the specific "geometric impedance" / "helical gauge fiber" / "symplectic coupling" predictive framing and keep only a generic "explored in a companion paper" sentence.

Severity: MEDIUM. Does not block broadcast on its own — the Paper 1 abstract's central claims (algebraic exactness, Berry-phase erratum, log-holonomy) are independent of this citation. But a domain expert following the citation will see a mismatch and may discount the surrounding framing. This is the single most important pre-broadcast item to fix.

### Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|---|---|---|---|---|
| (none — the paper makes no "first in the literature" or "novel" claim that I could find) | — | — | — | — |

The paper is structurally well-hedged on novelty: it does not assert priority on the SO(4,2) dynamical group framing, on graph Laplacian techniques, or on Berry phase definitions, all of which are cited to standard references. The "algebra-geometry duality" framing is presented as suggestive, not as a new theorem. No novelty-claim flags arise.

---

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---|---|---|---|
| F1 | `companion_alpha` bibitem title and §VI "geometric impedance" framing do not match actual Paper 2 (which is now Hopf-fibration spectral cubic, in Observations status) | B | CITE-WRONG-METADATA / CITE-DOESNT-SUPPORT | **MEDIUM** |
| F2 | Eq. (10) `$\sim 2/n$` asymptote has wrong sign (should be `$\sim -2/n$`) — cosmetic; magnitude/scaling claims downstream are correct | A | C (overstatement / sign typo) | **LOW** |
| F3 | §III.B "$k = 1.0$, $R^2 = 1.0$" wording reads as a numerical fit but is actually an analytical asymptote; finite-$n_{\max}$ fits give sub-unity values approaching 1 | A | C (mild) | **LOW** |
| F4 | §III.C and §III.D 2s/2p splitting numbers are GEOVAC-ONLY (not externally reproducible without internal code); structural claim is OK, specific numbers are not externally checked | A | B | **LOW-MEDIUM** (not blocking; recommend a one-line caveat) |
| F5 | §IV.C $\Delta E_{\text{Lamb}} \sim \alpha^5 \cdot mc^2$ omits Bethe log factor (conventional textbook informality) | A | C borderline | **LOW** (no fix needed) |
| F6 | Figure 1 (lattice visualization), Fig. 2 (eigenvalue bar chart), Fig. 3 (log-holonomy plot) are `[Placeholder]` boxes, not actual figures | A | (style) | **LOW-MEDIUM** for broadcast — three of three figures are placeholders. The captions still contain numerical statements ($\lambda_{2s}, \lambda_{2p}$, etc.) that read as if a figure exists. A pre-broadcast reader who skims the figures will see only placeholders, which is awkward. Recommend replacing with real figures OR replacing the fbox+text with a "Figure deferred to revised version" notice. |
| F7 | Berry-phase erratum (k = 2.113 retracted) is properly handled in §III.A; cross-checked against `debug/qa_sprint/berry_phase_reconciliation.md` and `debug/berry_phase_convergence/results.md` | A | A (positive) | n/a |

Severity totals: HIGH = 0, MEDIUM = 1, LOW = 5 (plus 1 LOW-MEDIUM cosmetic).

A–E verdict counts (Pass A claim inventory): A = 5, B = 4, C = 2, D = 1, E = 0.

Pass B citation-verdict counts: CITE-OK = 6 (barut1967, fock1935, biedenharn1981, berry1984, chung1997, paper13), CITE-WRONG-METADATA = 1 (companion_alpha; also borderline CITE-DOESNT-SUPPORT), CITE-MISATTRIBUTED = 0, CITE-CANT-FIND = 0.

---

## Broadcast readiness: YELLOW

Paper 1's load-bearing claims (algebraic exactness of the Rydberg spectrum, Berry-phase erratum, log-holonomy formula) are externally sound. All five external citations verify cleanly with correct authors/venues/years. The Berry-phase retraction is properly handled in the current text and matches the reconciliation memo.

The blocking-to-broadcast issue is F1: the `companion_alpha` bibitem points to a paper with a different title and a different mechanism than the §VI "Outlook" paragraph describes. A domain expert following the citation will read the Paper 2 abstract and find a Hopf-fibration spectral cubic, not the helical photon gauge fiber / symplectic impedance / geometric impedance framework Paper 1 promises. This is fixable in a single edit (update bibitem title + tighten §VI to describe Paper 2's actual mechanism, or to drop the specific helical/symplectic predictive language). Not a math error and not GREEN-blocking on its own, but a clean broadcast warrants the fix.

The three placeholder figures (F6) are cosmetic-LOW for an arXiv broadcast but read as unfinished if a reviewer flips to the figures expecting data. Recommend either generating real figures or rewording the placeholder boxes to a clear "figure deferred" notice.

The sign typo in Eq. (10) (F2) and the "$k=1.0, R^2=1.0$" wording (F3) are LOW cosmetic.

Recommended action: apply F1 + F6 fixes before broadcast (MEDIUM + cosmetic-medium); batch F2/F3 into a later cleanup pass.

---

## What I could NOT verify (hand to a human expert)

- The 2s/2p eigenvalue splitting numerical values ($\lambda_{2s}=0.0202$, $\lambda_{2p}=0.0238$, $\Delta E = 0.0035$ at $n_{\max}=10$) and the convergence sweep values (13% → 16% → 0.3% → 0.005% at $n_{\max} = 5/10/20/30$): per the instruction "text-level only, do NOT re-run computational drivers," I did not re-execute `geovac.lattice.GeometricLattice` + Lanczos. The qualitative claim (oscillatory decay to zero in the continuum limit) is structurally sound; the specific numerical values are framework-internal and externally unchecked. A human expert with hands-on access to the codebase could verify these in minutes.
- Whether Paper 1's intended companion paper has further unwritten content the bibitem still anticipates (Paper 2's actual evolution is documented in CLAUDE.md but a human PI is the right arbiter of whether the §VI Outlook reframe is desirable or whether the bibitem should simply be updated to match the existing Paper 2 title).

---

## Sources (web-verified citations)

- Barut, A. O. & Kleinert, H. (1967). [Phys. Rev. **156**, 1541 — Transition Probabilities of the Hydrogen Atom from Noncompact Dynamical Groups](https://link.aps.org/doi/10.1103/PhysRev.156.1541)
- Fock, V. (1935). [Z. Phys. **98**, 145 — Zur Theorie des Wasserstoffatoms (NASA ADS)](https://ui.adsabs.harvard.edu/abs/1935ZPhy...98..145F)
- Biedenharn, L. C. & Louck, J. D. (1981). [*Angular Momentum in Quantum Physics*, Encyclopedia of Math. and Its Applications 8 (Cambridge later editions)](https://books.google.com/books/about/Angular_Momentum_in_Quantum_Physics.html?id=8HWQwdXqcXAC)
- Berry, M. V. (1984). [Proc. R. Soc. Lond. A **392**, 45 — Quantal phase factors accompanying adiabatic changes](https://royalsocietypublishing.org/doi/10.1098/rspa.1984.0023)
- Chung, F. R. K. (1997). [*Spectral Graph Theory*, CBMS Regional Conference Series in Mathematics 92 (AMS)](https://www.ams.org/books/cbms/092/)
