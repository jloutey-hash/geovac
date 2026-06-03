# Confidence Review: Paper 0 — Angular Momentum as Information Geometry: Isotropic Packing and the Universal Angular Cross-Section

(Wave 1 re-fire — combined Pass A + Pass B; text-level only.)

## Calibration check

Not a calibration run. A prior AUDITOR-only Pass A exists at `debug/audit_paper0.md`. Some fixes from that earlier audit have already been applied in the current `.tex`: §V line 508 now reads "96.0%" (the 94.1% stale figure is gone). This Wave 1 re-fire runs both Pass A and Pass B against the current `.tex` and inherits one open item (Ω(0) = 2 normalization) from the earlier audit for explicit resolution.

## Pass A — Content audit

### Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|-------|----------|---------|----------|---------------------|
| 1 | "Shell $k$ has angular capacity $2k-1$, isomorphic to $Y_{\ell m}$ with $\ell=k-1$" | Abstract; §III; Eqs. (3),(5),(6) | A | EXTERNAL (SO(3) rep theory; $\dim V_\ell = 2\ell+1$ standard) | Algebraic identity from Eq. (1): $A_k = \pi k^2 d_0^2 - \pi(k-1)^2 d_0^2 = \pi d_0^2(2k-1)$; division by $\sigma_0 = \pi d_0^2/2$ gives $N_k = 2(2k-1)$; the factor $2k-1 = 2(k-1)+1 = 2\ell+1$ is the standard SO(3) irrep dim. |
| 2 | $\sum_{k=1}^n (2k-1) = n^2$ matches hydrogen $n^2$ degeneracy | §IV.A Eq. (7) | A | EXTERNAL (textbook arithmetic) | Sum of consecutive odd integers, induction. Hydrogen $n^2$ degeneracy is Fock/Bargmann 1935-36. |
| 3 | Factor of 2 ↔ orientability of compactified packing manifold ↔ spin multiplicity | Abstract; §IV.B | C | MIXED (Z₂ doubling under $S^2$ orientability is math; identification with spin-½ is interpretive) | Body explicitly hedges as "structural correspondence, not a derivation"; abstract reads more strongly. |
| 4 | $K_{\mathrm{vac}} = -1/16$ derivable from Fock projection: $c^2(n,l) = (1/16)[1-l(l+1)/(n(n+1))]$, equals $1/16$ for $l=0$; equivalently $1/\Omega^4(0)$ with $\Omega(0)=2$ | §VI.C lines 631-638 | **B with internal-inconsistency caveat** | GEOVAC-ONLY plus DLMF Gegenbauer recurrence | Cross-corpus check against Paper 18 §VII (lines 2632-2650) and Paper 7 §III (lines 161-209): the formula $c^2(n,l) = (1/16)[1 - l(l+1)/(n(n+1))]$ is internally consistent. The "Ω(0) = 2" claim is conditional on $p_0 = 1$ (Paper 7 line 200: "$p_0 = 1/n$"; ground state $n=1 \Rightarrow p_0 = 1$). Paper 7 line 206 says explicitly "At $\mathbf{p}=0$: $\Omega = 2/p_0$." So Ω(0) = 2 only at the ground-state energy shell. Paper 18 itself uses TWO inconsistent definitions: line 173 has $\Omega = (1+p^2/p_0^2)^{-1}$ (which would give $\Omega(0) = 1$), but line 2364 has $\Omega = 2p_0/(p^2+p_0^2)$ (gives $\Omega(0) = 2/p_0$). The internal-corpus inconsistency is Paper 18's; Paper 0 inherits it. |
| 5 | "ground-state eigenvalue of $H$ converges to $-0.5$ Ha" as $n_{\max} \to \infty$ | §VI.C line ~643 | B | GEOVAC-ONLY (Paper 18 §VII; Paper 7 conformal equivalence) | Uses $H = K_{\mathrm{vac}} L$ alone (no node weights $W$). Per Paper 18 §III line 165-168, the *full* per-shell spectrum $E_n = -Z^2/(2n^2)$ requires $H = \kappa\mathcal{L} + W$. Paper 0's claim is correctly scoped to the *ground-state* only — but the abstract's "spectrally converge to the hydrogen eigenvalues" (plural) is broader than what the body supports without $W$. |
| 6 | "Paper 15 ... recovers 96.0\% of the H$_2$ dissociation energy" | §V line 508 | A | EXTERNAL within-corpus | Cross-corpus check: Paper 15 abstract, intro, table caption, and conclusion all give 96.0%. CLAUDE.md §2 best-results table confirms. **Fix from prior audit already applied.** |
| 7 | "errors ranging from <0.1\% for hydrogen to ~6\% for LiH equilibrium geometry" | §IX Conclusions line 802 | C (minor) | MIXED | CLAUDE.md §2: LiH best is 5.3% via composed at l_max=2 (Paper 17). "~6%" rounds up; conservative but slightly imprecise. |
| 8 | Five-level hierarchy table (S³, prolate spheroidal, hyperspherical, mol-frame hypersp., composed) | §V Table II | A/B | MIXED | Geometries themselves are external/standard (Fock 1935 for S³; prolate spheroidal is textbook; hyperspherical is standard). GeoVac-internal interpretive layer is "universal $(\ell,m)$ angular cross-section." |
| 9 | "the lattice is the invariant; the continuum geometries are projections" | §V.A (combinatorial invariant subsection) | C | GEOVAC-ONLY interpretive | Reads as assertion; the very next paragraph correctly hedges ("we state this as a mathematical observation, not a metaphysical claim"). Internal tension with CLAUDE.md §1.5 dual-description rhetoric rule. |
| 10 | "Paper 7 provides 18 symbolic proofs for the $S^3$ embedding" | §V.A; §VI.C; conclusions | B | GEOVAC-ONLY (Paper 7 + GeoVac test suite) | Internal-consistency framing only — these tests check GeoVac code against GeoVac claims. Honest broadcast framing would say "18 internal symbolic test cases." |
| 11 | Step 1 placement: "two points ... separated by the fundamental distance $d_0$" placed "on an isotropic shell of radius $r_1 = d_0$" | §II.B Step 1 | B (geometric wording) | minor inconsistency | Two points on a circle of radius $d_0$ have a separation that depends on placement; "separation $d_0$" and "radius $d_0$" coincide only for a 60° arc. The fundamental area $\sigma_0 = \pi d_0^2/2$ is the disk area divided by 2 (the two points). Wording could clarify that $d_0$ is the shell radius and $\sigma_0$ is area-per-state from the initialization disk. |

### Numbers I recomputed

| Claim | Paper's figure | Independent reference | My recomputed value | Survives? |
|---|---|---|---|---|
| H$_2$ dissociation recovery (Paper 15) | 96.0% | Paper 15 abstract/intro/conclusion all 96.0%; CLAUDE.md §2 table | 96.0% | YES |
| $\sum_{k=1}^n (2k-1) = n^2$ | $n^2$ | Textbook | $n^2$ exactly (induction) | YES |
| Annular area $A_k = \pi d_0^2(2k-1)$ (Eq. 1) | Eq. 1 | Elementary geometry | $\pi k^2 d_0^2 - \pi(k-1)^2 d_0^2 = \pi d_0^2(2k-1)$ | YES |
| $N_k = 2(2k-1)$ (Eq. 2) | Eq. 2 | Elementary | $A_k/\sigma_0 = (2k-1)/(1/2) = 2(2k-1)$ | YES |
| Cumulative count = $2k^2$ (Table I) | $2k^2$ | Arithmetic | $\sum_{k=1}^n 2(2k-1) = 2 n^2$ | YES |
| $c^2(n,l) = (1/16)[1 - l(l+1)/(n(n+1))]$ | §VI.C | Paper 18 §VII (line 2640) verbatim match | identical | YES (internal consistency) |
| $\kappa = -1/16 = -1/\Omega^4(0)$ with $\Omega(0) = 2$ | §VI.C | Paper 7 line 163, 206: $\Omega(0) = 2/p_0$, equals 2 only when $p_0 = 1$ | $\Omega(0) = 2/p_0$; equals 2 at ground-state $p_0 = 1$ only | survives WITH unstated $p_0 = 1$ caveat; abstract claim is unconditional |
| LiH error ~6% | "~6%" | CLAUDE.md §2: 5.3% best (Paper 17) | 5.3% | survives as conservative round-up |

### Circularity map

GEOVAC-ONLY chains explicitly listed:

1. **$K_{\mathrm{vac}} = -1/16$ kinetic scale derivation.** Paper 0 §VI.C → Paper 18 §VII (v2.26.1) → Paper 7 §III (Fock 1935 + chordal-distance identity). The Fock paper itself is EXTERNAL; the *identification* of $1/\Omega^4(0)$ as the kinetic scale, and the $\Omega(0) = 2$ assignment, are internal corpus content. The Gegenbauer/Chebyshev three-term recurrence yielding $c^2 = 1/16$ for $l=0$ is mathematically external (DLMF 18.9), but application here is GeoVac-internal.

2. **"Ground state → −0.5 Ha as $n_{\max} \to \infty$."** Paper 0 §VI.C → Paper 18 §VII Eq. (kappa_derived). Numerical verification ($\lambda_{\max}(5)=6.62$, $\lambda_{\max}(10)=7.61$, $\lambda_{\max}(20)=7.90 \to 8$) is GeoVac code testing GeoVac predictions.

3. **"Paper 7 ... 18 symbolic proofs"** — pure internal consistency (tests of GeoVac code against GeoVac code per CLAUDE.md §9).

4. **Universal angular cross-section across Levels 1–5.** Paper 0 §V → Papers 7, 11, 13, 15, 17. The *separation of variables* in each geometry is standard external math; the "this same $(\ell,m)$ is the universal fiber" framing is GeoVac's interpretive contribution.

5. **5.3% LiH accuracy.** Paper 0 §IX → Paper 17 → GeoVac composed-geometry solver. Per CLAUDE.md §1.5 benchmarking rule, Paper 17 does include external comparison; Paper 0's restatement is downstream.

### Overstatement findings

| Exact phrase | Suggested honest replacement |
|---|---|
| (Abstract) "construction yields shell capacities $2k-1$ for shell~$k$" | "construction yields per-shell *angular* capacities $2k-1$" — Table I shows shell 1 has 2 states, not 1; the $2k-1$ is the angular factor after dividing out the global 2. |
| (Abstract) "spectrally converge to the hydrogen eigenvalues" | "ground-state eigenvalue spectrally converges to the hydrogen ground-state energy" OR cross-reference to Paper 18 §III that the full per-shell spectrum needs node weights $W$. |
| (Abstract last sentence) "making it the foundational motif of the framework" | "we treat it as a foundational motif of the framework" — minor framing softening; current phrasing self-assigns status. |
| (§V.A) "The lattice is the invariant; the continuum geometries are projections." | "Equivalently, the lattice may be read as the invariant and the continuum geometries as projections; the converse reading is equally consistent with the mathematics (cf. CLAUDE.md §1.5 dual-description rule)." |
| (§VI.C) "$\Omega(0) = 2$ is the stereographic conformal factor at the projection center" | "$\Omega(0) = 2/p_0 = 2$ at the ground-state energy shell $p_0 = 1$" — make the implicit $p_0$ choice explicit; cross-reference Paper 7 line 163. |
| (§VI.C) "18 symbolic proofs validate the conformal geometry" | "18 internal symbolic test cases validate the conformal geometry" — flag internal-consistency framing. |
| (§IX Conclusions) "~6\% for LiH equilibrium geometry" | "~5\% (specifically 5.3\%) for LiH equilibrium geometry" — match CLAUDE.md §2. |

## Pass B — Citation and novelty

### Citation table

| `\cite` key | Claimed as | Verdict | What I found |
|---|---|---|---|
| `fock1935` | Fock, "Zur Theorie des Wasserstoffatoms," Z. Phys. 98, 145 (1935) | **CITE-OK** | Z. Phys. 98, 145-154 (March 1935). DOI 10.1007/BF01336904. NASA/ADS 1935ZPhy...98..145F. |
| `bargmann1936` | Bargmann, "Zur Theorie des Wasserstoffatoms: Bemerkungen ...," Z. Phys. 99, 576 (1936) | **CITE-OK** | Z. Phys. 99, 576-582 (March 1936). DOI 10.1007/BF01338811. |
| `barut1967` | Barut & Kleinert, "Transition Probabilities of the Hydrogen Atom from Noncompact Dynamical Groups," Phys. Rev. 156, 1541 (1967) | **CITE-OK** | Phys. Rev. 156, 1541 (Apr 25, 1967). DOI 10.1103/PhysRev.156.1541. |
| `biedenharn1981` | Biedenharn & Louck, *Angular Momentum in Quantum Physics: Theory and Application*, Addison-Wesley, Reading, MA, 1981 | **CITE-OK** | Confirmed: Addison-Wesley, Reading, 1981 (Encyclopedia of Mathematics and its Applications, vol. 8). |
| `bohr1913` | Bohr, "On the Constitution of Atoms and Molecules," Phil. Mag. 26, 1 (1913) | **CITE-OK** | Phil. Mag. (Series 6) vol. 26, no. 151, pp. 1-25 (July 1913). DOI 10.1080/14786441308634955. |
| `paper7`, `paper11`, `paper13`, `paper15`, `paper17` | Internal GeoVac references | **CITE-OK (internal)** | All files present at the cited paths in `papers/group2_quantum_chemistry/` (Papers 11, 13, 15, 17) and `papers/group3_foundations/` (Paper 7). |

### Problems found

None for Pass B. All five external bibitems verified against authoritative sources (Springer / APS / Taylor & Francis); years, volumes, page numbers, and titles match. No CITE-MISATTRIBUTED / CITE-CANT-FIND / CITE-DOESNT-SUPPORT findings.

### Priority / novelty claims

| Claim (verbatim) | Location | Searched | Prior art found? | Recommendation |
|---|---|---|---|---|
| "What is new here is the specific packing algorithm that produces the $(\ell, m)$ labeling from geometric–combinatorial reasoning independent of the operator formalism" | §I, intro, line 119-122 | Web search not feasible at confidence the literature is fully indexed; the claim is modest ("what is new here is X" framing, not "first in the literature") | None found in cursory search; but the search ceiling rule applies (cannot confirm novelty). | Already honestly framed — "what is new here is" reads as the author's contribution, not as an absolute priority claim. **No change needed.** |
| (No "first in the literature" or "first published" claims in Paper 0.) | — | — | — | Paper 0 does not make explicit priority claims; this is good practice. |

## Combined severity table

| # | Finding | Pass | Verdict | Severity |
|---|---|---|---|---|
| 1 | Ω(0) = 2 unconditional in §VI.C / abstract; actually conditional on $p_0 = 1$ (Paper 7 line 163, 206) | A | C (overstatement) **with internal-corpus inconsistency in Paper 18** | **MEDIUM** |
| 2 | Abstract "shell capacities $2k-1$" ambiguous re factor-of-2 | A | C | MEDIUM |
| 3 | Abstract "spectrally converge to the hydrogen eigenvalues" (plural) broader than $H = K_{\mathrm{vac}} L$ delivers | A | C | MEDIUM |
| 4 | "The lattice is the invariant; the continuum geometries are projections" assertion vs CLAUDE.md §1.5 dual-description rule | A | C | MEDIUM |
| 5 | "18 symbolic proofs" framing without internal-consistency caveat | A | B (label needed) | LOW |
| 6 | "~6\%" LiH (actual best is 5.3%) | A | C | LOW |
| 7 | "foundational motif of the framework" mildly self-promoting | A | C | LOW |
| 8 | Step 1 wording: "separated by $d_0$" on circle of "radius $d_0$" | A | B | LOW |
| 9 | (All external citations checked clean) | B | CITE-OK ×5 | — |

**Severity totals:** HIGH = 0, MEDIUM = 4, LOW = 4.

**Pass A verdict-type totals:** A = 3 (claims 1, 2, 6), B = 4 (claims 4, 5, 10, 11), C = 4 (claims 3, 7, 9, + abstract/overstatement findings), D = 0, E = 0.

**Pass B citation totals:** CITE-OK = 5 (external) + 5 (internal Paper-*); CITE-WRONG-METADATA = 0; CITE-MISATTRIBUTED = 0; CITE-DOESNT-SUPPORT = 0; CITE-CANT-FIND = 0.

## Broadcast readiness: **YELLOW**

The core mathematical claim — uniform-density isotropic packing produces $(\ell, m)$ labeling via the elementary annular-area identity $A_k = \pi d_0^2(2k-1)$ — is solid (verdict A) and the structural correspondences ($n$ as cumulative grouping; factor of 2 from orientability) are honestly framed in the body. All five external citations check clean. The earlier audit's HIGH finding (94.1% → 96.0%) is already fixed.

What blocks GREEN: four MEDIUM-severity content findings, all in the same category — abstract/headline language is slightly broader than what the body actually delivers. Three are abstract-precision tightenings (shell capacity / "hydrogen eigenvalues" plural / "foundational motif"); one is the unstated $p_0 = 1$ caveat behind the "$\Omega(0) = 2$" claim, which is a Paper 18 corpus-internal inconsistency Paper 0 inherits (Paper 18 itself uses $\Omega = (1+p^2/p_0^2)^{-1}$ in §III and $\Omega = 2p_0/(p^2+p_0^2)$ in §VIII — two different conformal-factor normalizations under the same symbol). The "lattice is the invariant" framing tension (§V.A vs the very next paragraph) is also MEDIUM. After these are softened, Paper 0 goes to GREEN. No HIGH findings. No CITE-* problems.

## What I could NOT verify (hand to a human expert)

1. **The Ω(0) = 2 normalization convention.** Paper 7 line 163 defines $\Omega(\mathbf{p}) = 2p_0/(p^2 + p_0^2)$, and line 206 confirms $\Omega(0) = 2/p_0$. Paper 18 §III line 173 defines $\Omega = (1 + p^2/p_0^2)^{-1}$ which would give $\Omega(0) = 1$. Paper 18 §VIII line 2364 reverts to $\Omega = 2p_0/(p^2 + p_0^2)$. The "$1/\Omega^4(0) = 1/16$" requires $\Omega(0) = 2$, which selects the Paper 7 / §VIII convention AND requires $p_0 = 1$ (ground-state energy shell). Recommended fix: in *Paper 18*, harmonize the Ω-definition between §III and §VIII to a single form (the $\Omega = 2p_0/(p^2+p_0^2)$ convention is the one used in the kinetic-scale derivation), and in *Paper 0* §VI.C add a half-sentence noting the ground-state energy-shell convention ($p_0 = 1$) under which $\Omega(0) = 2$. This is a paper-18-and-7 issue Paper 0 inherits; a domain expert in conformal-projection geometry could rule on which convention is most standard.

2. **Independent reproduction of $\lambda_{\max}(n_{\max}) \to 8$.** Paper 18 §VII reports the convergence, but this is GeoVac code testing GeoVac predictions. A spectral graph theorist could verify in ~10 minutes by computing the largest eigenvalue of $L = D - A$ on the $\{(n,\ell,m): 1 \le n \le n_{\max}, 0 \le \ell < n, -\ell \le m \le \ell\}$ vertex set with the edge set defined in Paper 0 §VI.B.

3. **The 18 symbolic proofs in Paper 7.** External verification of at least one of the 18 (e.g., the chordal-distance identity, Eq. 6 of Paper 7) before broadcast.

4. **"Factor of 2 from orientability of $S^2$" → spin-½.** Body labels this "structural correspondence, not a derivation"; honest. A mathematical physicist may want to see the explicit (non-)connection to spin-structure theory on $S^2$ (the obstruction is the $w_2$ Stiefel-Whitney class) stated.

5. **The novelty of the packing-algorithm framing itself.** Paper 0 phrases this as "what is new here is" (an authorial contribution claim) rather than "first in the literature" (a priority claim). The phrasing is honest and does not require search-based confirmation. But a domain expert (e.g., in geometric / information-theoretic foundations of QM, à la Wheeler / Hardy / Caticha) would be the right judge of priority.

## Recommended edits (prioritized)

**MEDIUM** (errata batch):
- Abstract: "construction yields shell capacities $2k-1$" → "construction yields per-shell *angular* capacities $2k-1$".
- Abstract: "spectrally converge to the hydrogen eigenvalues" → "the ground-state eigenvalue spectrally converges to the hydrogen ground-state energy; the full per-shell spectrum requires node weights, addressed in Paper 18 §III."
- §V.A: soften "The lattice is the invariant; the continuum geometries are projections" to the dual-reading form.
- §VI.C: add "$\Omega(0) = 2/p_0$, equal to 2 at the ground-state energy shell $p_0 = 1$" to make the implicit convention explicit. Recommend a one-sentence note that the Ω-definition harmonization issue in Paper 18 §III vs §VIII is a known cross-paper inconsistency to be cleaned up; flag to PI rather than silently resolving.

**LOW** (cleanup batch):
- §IX Conclusions: "~6%" → "~5% (specifically 5.3%)" for LiH.
- §II.B Step 1: minor wording cleanup ($d_0$ as shell radius vs separation).
- §VI.C: "18 symbolic proofs" → "18 internal symbolic test cases".
- Abstract last sentence: "the foundational motif of the framework" → "a foundational motif of the framework".

## Summary line for dispatcher

Paper 0 — **YELLOW** — Pass A: A=3, B=4, C=4, D=0, E=0; Pass B: CITE-OK=10, all other CITE-* = 0; HIGH=0, MEDIUM=4, LOW=4. **Top finding:** the abstract/§VI.C claim "$\Omega(0) = 2$" is unconditional but requires $p_0 = 1$ (ground-state energy shell) per Paper 7 lines 163, 206, and inherits an unresolved Ω-definition inconsistency between Paper 18 §III (Ω(0)=1) and §VIII (Ω(0)=2/p_0) — a corpus-level harmonization is needed before broadcast.
