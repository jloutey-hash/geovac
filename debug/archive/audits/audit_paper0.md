# Confidence Audit: Paper 0 — Angular Momentum as Information Geometry: Isotropic Packing and the Universal Angular Cross-Section

## Calibration check

This is not a calibration run (no stated known-honest answer was provided). The audit proceeds against the corpus as live work, with cross-corpus check (Paper 18, Paper 15, Paper 17, CLAUDE.md §2) and exact-text verification.

## Claim inventory + verdicts

| # | Claim | Location | Verdict | Rests on | Evidence I produced |
|---|-------|----------|---------|----------|---------------------|
| 1 | "Shell $k$ has angular capacity $2k-1$, isomorphic to $Y_{\ell m}$ with $\ell=k-1$" | Abstract; §III; Eqs. (3),(5),(6) | A | EXTERNAL (SO(3) rep theory; $\dim V_\ell = 2\ell+1$ is standard) | Direct algebraic identity from Eq. (1) annular-area formula. The annular area is $\pi d_0^2(2k-1)$ by elementary geometry; divided by $\sigma_0 = \pi d_0^2/2$ gives $N_k = 2(2k-1)$. The angular factor $2k-1 = 2(k-1)+1 = 2\ell+1$ is the standard SO(3) multiplicity. Identity, not coincidence. |
| 2 | Cumulative count $\sum_{k=1}^{n}(2k-1) = n^2$ matches hydrogen $n^2$ degeneracy | §IV.A; Eq. (7) | A | EXTERNAL (elementary arithmetic identity) | $\sum_{k=1}^{n}(2k-1) = n^2$ is the textbook sum of consecutive odd integers. Hydrogen $n^2$ degeneracy is established (Schrödinger 1926, Fock 1935). Paper is honest that the grouping is a structural correspondence, not derived. |
| 3 | "factor of 2 ... orientability of compactified packing manifold ($S^2$)" → spin multiplicity | Abstract; §IV.B | C/D | MIXED (Z₂ doubling under $S^2 \cong \mathbb{C}^*$ orientability is mathematical; identification with spin is interpretive) | The paper explicitly states this is a *structural correspondence*, not a derivation. The text is appropriately hedged. However, the abstract's bare phrase "a factor of 2 arises from the orientability of the compactified packing manifold ($S^2$), corresponding to spin multiplicity" reads as derivation to a hostile outsider; needs slight softening (verdict C, minor). |
| 4 | $K_{\mathrm{vac}} = -1/16$ derivable from Fock projection: $c^2(n,l) = (1/16)[1-l(l+1)/(n(n+1))]$, equals $1/16$ for $l=0$; equivalently $1/\Omega^4(0)$ with $\Omega(0)=2$ | §VI.C | B (cross-paper internal consistency) | GEOVAC-ONLY (rests on Paper 7 + Paper 18 §VII v2.26.1 derivation) | Verified against Paper 18 §VII lines 2632–2650 verbatim: same formula, same equivalence, same conformal-factor reading. Paper 0's text matches Paper 18 v2.26.1's "Fock weight derivation" exactly. The Gegenbauer/Chebyshev recurrence used to obtain $c^2 = 1/16$ is established (DLMF 18.9). Listed B not A because the *application* to hydrogen energy calibration is GeoVac-internal; the mathematics behind the recurrence is external. |
| 5 | "ground-state eigenvalue of $H$ converges to $-0.5$ Ha" as $n_{\max} \to \infty$ | §VI.C (line ~644) | B | GEOVAC-ONLY (Paper 18 §VII Eq. eq:kappa_derived; Paper 7 conformal equivalence) | This claim uses $H = K_{\mathrm{vac}}\mathcal{L}$ (no node weights, only the kinetic scale). Paper 18 §VII derives the same: $\kappa \lambda_{\max} = -1/2$ Ha with $\lambda_{\max} \to 8$. Internally consistent. **Subtlety:** Paper 18 §III line 167 explicitly states "Per-shell energies $E_n = -Z^2/(2n^2)$ emerge from the full graph Hamiltonian $H = \kappa\mathcal{L} + W$ ... not from $\kappa$ alone." Paper 0's claim is about the *ground-state* only, which $H = \kappa\mathcal{L}$ does deliver (as Paper 18 confirms). So the claim is correctly scoped. **Recommend** adding a one-line note that this is the ground state only; the full $E_n = -Z^2/(2n^2)$ requires the node-weight term $W$, per Paper 18 §III. |
| 6 | "Paper 15 ... recovers 94.1\% of the H₂ dissociation energy using this geometry" | §V, last paragraph of Level 4 entry (line ~509) | E | EXTERNAL contradiction | **Stale figure.** Paper 15 explicitly states 96.0\% in its abstract (line 44), introduction (line 106), Table caption (line 690), and concluding sentence (line 884). CLAUDE.md §2 best-results table also gives 96.0\%. 94.1\% is from an earlier draft of Paper 15. Fix: replace 94.1\% with 96.0\%. |
| 7 | "errors ranging from ${<}\,0.1\%$ for hydrogen to ${\sim}\,6\%$ for LiH equilibrium geometry" | §IX Conclusions | C (minor) | MIXED | Hydrogen ${<}0.1\%$ is correct per CLAUDE.md §5. LiH best result is 5.3\% (CLAUDE.md §2, Paper 17 §1477). "$\sim 6\%$" is slightly imprecise but conservative (rounded up). Recommend: "${\sim}\,5\%$" or "5.3\%" for accuracy. |
| 8 | Five-level natural-geometry hierarchy table (S³, prolate spheroidal, hyperspherical, mol-frame hypersp., composed) | §V Table II | A/B | MIXED | The geometries themselves are external/established (Fock 1935 for S³; prolate spheroidal is textbook; hyperspherical is standard). Their application to GeoVac and the "universal $(\ell,m)$ angular cross-section" framing rests on Papers 7, 11, 13, 15, 17 (GEOVAC-ONLY interpretive layer). Honest framing in body. |
| 9 | "lattice is the invariant; the continuum geometries are projections" | §V.A (combinatorial invariant subsection) | C/B | GEOVAC-ONLY (interpretive) | The text appropriately follows with "We state this as a mathematical observation, not a metaphysical claim ... Whether one regards the discrete or continuum description as more fundamental is a choice of interpretation." This matches CLAUDE.md §1.5 rhetoric rule. The hedging is correct. Slight tension: "the lattice is the invariant" reads like an assertion. **Recommend** rewording to "Equivalently, the lattice can be read as the invariant and the continuum geometries as projections; the converse reading is equally consistent with the mathematics." |
| 10 | "Paper 7 ... 18 symbolic proofs validate the conformal geometry" | §VI.C; conclusions | B | GEOVAC-ONLY (rests on Paper 7 + GeoVac test suite) | Internal consistency only — these are tests of GeoVac code against GeoVac claims. Honest broadcast framing would say "verified by 18 internal symbolic test cases." |
| 11 | "Step 1 ... two points ... separated by the fundamental distance $d_0$" placed "on an isotropic shell of radius $r_1 = d_0$" | §II.B Step 1 | B (internal geometry) | EXTERNAL geometric inconsistency (minor) | Two points on a circle of radius $d_0$ have maximum separation $2 d_0$ (diametrically opposed), not $d_0$. The text equates "separated by $d_0$" with the circle radius being $d_0$. To get separation $d_0$ with two points on a circle of radius $d_0$, one needs a 60° arc — but Step 1's text says "minimum needed to define a distance," not specifying placement on the circle. The fundamental area $\sigma_0 = \pi d_0^2 / 2$ that follows uses the area of the disk (radius $d_0$) divided by 2 (the two points). This is consistent IF $\sigma_0$ is defined as "area-per-state from the initialization disk," not "area-per-state with $d_0$ being a separation." Recommend: minor wording cleanup to clarify that $d_0$ is the shell radius and $\sigma_0 = (\text{disk area})/N_{\mathrm{init}}$, rather than describing $d_0$ as a "separation." |

## Numbers I recomputed

| claim | paper's figure | independent reference | my recomputed value/error | survives? |
|-------|----------------|----------------------|---------------------------|-----------|
| H₂ dissociation energy recovery (Paper 15) | 94.1% (Paper 0 §V) | Paper 15 abstract, intro, Table 1, conclusion all state 96.0% | 96.0% per cross-corpus | **NO — figure is stale; fix to 96.0%** |
| Cumulative $\sum_{k=1}^n(2k-1)$ matches $n^2$ | $n^2$ | Standard arithmetic | $\sum_{k=1}^n (2k-1) = n^2$ exactly by induction | YES |
| Annular area $A_k = \pi d_0^2 (2k-1)$ | (Eq. 1) | Elementary geometry | $\pi k^2 d_0^2 - \pi (k-1)^2 d_0^2 = \pi d_0^2 (2k-1)$ | YES |
| $N_k = 2(2k-1)$ | (Eq. 2) | Elementary | $A_k / \sigma_0 = \pi d_0^2(2k-1) / (\pi d_0^2 / 2) = 2(2k-1)$ | YES |
| $\kappa = -1/16$ from $\kappa \lambda_{\max} = -1/2$ Ha and $\lambda_{\max} = 8$ | (implicit via Paper 18 §VII) | Paper 18 §VII Eq. eq:kappa_derived | $\kappa = (-1/2)/8 = -1/16$ | YES (internal) |
| $\Omega(0) = 2$ at $p = 0$ for Fock projection | (§VI.C) | Fock 1935; Paper 7 | $\Omega(p) = (1 + p^2/p_0^2)^{-1}$, so $\Omega(0) = 1$ unless normalization is $\Omega = 2/(1+p^2/p_0^2)$ — wait | Need to verify (see "What I could NOT verify") |
| Best results: H < 0.1%, LiH ~6% | "${\sim} 6\%$" | CLAUDE.md §2; Paper 17 = 5.3% | 5.3%, not 6% | survives at one-sig-fig precision (5.3 rounds to ~5, not 6); minor overstatement |

## Circularity map

The GEOVAC-ONLY chains, stated explicitly:

1. **$K_{\mathrm{vac}} = -1/16$ as kinetic scale calibration.** Chain: Paper 0 §VI.C → Paper 18 §VII (v2.26.1 derivation) → Paper 7 (S³ conformal equivalence + Fock 1935). The Fock paper is external; the *interpretive identification* of $1/\Omega^4(0)$ as the GeoVac kinetic scale is internal. The Gegenbauer recurrence giving $c^2 = 1/16$ for $l=0$ is mathematically external (DLMF 18.9), but its application here is GeoVac-internal.

2. **"Ground state $\to -0.5$ Ha as $n_{\max} \to \infty$."** Chain: Paper 0 §VI.C → Paper 18 §VII Eq. eq:kappa_derived. The graph $\lambda_{\max} \to 8$ claim has numerical verification in Paper 18 (5: 6.62; 10: 7.61; 20: 7.90) but that is GeoVac code against its own predictions. No independent reproduction.

3. **"18 symbolic proofs validate the conformal geometry."** Chain: Paper 0 → Paper 7 → GeoVac test suite (`tests/test_fock_projection.py`, `tests/test_fock_laplacian.py`). Pure internal consistency; the tests check GeoVac claims against GeoVac code.

4. **The universal angular cross-section across Levels 1–5.** Chain: Paper 0 §V → Papers 7, 11, 13, 15, 17. The *separation of variables* in each natural geometry (S³, prolate spheroidal, hyperspherical, mol-frame hypersp., composed fiber) is standard external math. The "this same $(\ell, m)$ cross-section is the universal fiber" framing is GeoVac's interpretive claim — defensible but interpretive.

5. **The "5.3\% LiH" accuracy.** Chain: Paper 0 conclusion → Paper 17 → GeoVac composed-geometry solver. Per CLAUDE.md §1.5 benchmarking rule, the LiH value should be (and is, in Paper 17) compared against an external reference. But Paper 0's restatement is downstream of GeoVac internal code.

These chains are not *wrong* — they are appropriately scoped — but they should be flagged in any broadcast as "internal verification, not externally reproduced."

## Overstatement findings

| Exact phrase | Suggested honest replacement |
|--------------|------------------------------|
| (Abstract) "construction yields shell capacities $2k-1$ for shell $k$" | "construction yields per-shell angular capacities $2k-1$ for shell $k$ (with a global factor of 2 giving total occupancy $2(2k-1)$ per shell, addressed below)" — disambiguate "shell capacity" since Table I shows shell 1 = 2 states, not 1. |
| (§V, line ~509) "Paper 15 recovers 94.1\% of the H$_2$ dissociation energy using this geometry" | "Paper 15 recovers 96.0\% of the H$_2$ dissociation energy using this geometry" — **HIGH severity: stale figure contradicted by Paper 15's own abstract.** |
| (§IX Conclusions) "errors ranging from ${<}\,0.1\%$ for hydrogen to ${\sim}\,6\%$ for LiH equilibrium geometry" | "errors ranging from ${<}\,0.1\%$ for hydrogen to ${\sim}\,5\%$ (specifically 5.3\%) for LiH equilibrium geometry" — minor. |
| (§V.A) "The lattice is the invariant; the continuum geometries are projections." | "Equivalently, one may read the lattice as the invariant and the continuum geometries as projections; the converse reading is equally consistent with the mathematics, as discussed below." Better aligns with CLAUDE.md §1.5 dual-description rhetoric rule. |
| (Abstract) "spectrally converge to the hydrogen eigenvalues" | "ground-state eigenvalue spectrally converges to the hydrogen ground-state energy" — the *full* hydrogen spectrum $E_n = -Z^2/(2n^2)$ requires node weights $W$ (Paper 18 §III line 167), not the graph Laplacian alone. Calling the convergence "to the hydrogen eigenvalues" (plural) is slightly more than the body delivers without $W$. |
| (Abstract, last sentence) "making it the foundational motif of the framework" | "we treat it as a foundational motif of the framework" — minor framing softening; the current phrasing reads as self-assigned status. |

## Severity table

| Finding | Type | Severity |
|---------|------|----------|
| H₂ figure 94.1% → 96.0% (item 6) | E (wrong number) | **HIGH** (would embarrass in front of an expert who knows Paper 15) |
| Ground-state-only claim needs node-weight footnote (item 5) | C (overstatement / abstract→body gap) | MEDIUM |
| "shell capacities $2k-1$" abstract ambiguity (item 1 / overstatement #1) | C | MEDIUM |
| "lattice is the invariant" framing (item 9) | C | MEDIUM |
| "Paper 15 has 18 symbolic proofs" internal-consistency framing (item 10) | B (label needed) | LOW |
| "spin from orientability" abstract bareness (item 3) | C | LOW (body hedges adequately) |
| LiH "~6%" vs 5.3% (item 7) | C | LOW |
| Step 1 wording about "separated by $d_0$" on a circle of radius $d_0$ (item 11) | B (minor geometric wording) | LOW |
| "foundational motif" rhetoric (overstatement #6) | C | LOW |

**Severity totals:** HIGH = 1, MEDIUM = 3, LOW = 5

**Verdict-type totals:** A = 2 (claims 1, 2), B = 4 (claims 4, 5, 10, 11), C = 5 (claims 3, 6 overstatements, 7, 9, plus abstract softening items), D = 0 (one calibration item flagged below), E = 1 (claim 6: stale H₂ figure)

## Broadcast readiness: YELLOW

The paper's core mathematical claim — that uniform-density isotropic packing produces the $(\ell, m)$ angular labeling via an elementary annular-area identity — is solid and externally verified (A). The structural correspondences ($n$ as cumulative grouping, factor of 2 from orientability) are honestly framed in the body as correspondences-not-derivations.

Three issues block GREEN:

1. **One genuine wrong number** (verdict E): Paper 0 §V cites Paper 15's H₂ result as 94.1\%, but Paper 15 itself reports 96.0\% in four separate places. An expert reader cross-referencing the series will catch this immediately. Fix is a one-character edit (94 → 96, 1 → 0).

2. **Abstract→body precision gap** (verdict C, multiple): the abstract uses "shell capacities $2k-1$" (which is the *angular* capacity post-factorization, not the total $N_k = 2(2k-1)$ in Table I) and asserts "spectrally converge to the hydrogen eigenvalues" (which requires $H = \kappa\mathcal{L} + W$ per Paper 18 §III, not the graph Laplacian alone delivered in Paper 0 §VI.C).

3. **One framing assertion** ("the lattice is the invariant; the continuum geometries are projections," §V.A) reads more strongly than CLAUDE.md §1.5's dual-description rhetoric rule allows. The very next paragraph correctly hedges with "Whether one regards the discrete or continuum description as more fundamental is a choice of interpretation, not a consequence of the mathematics." Mild internal tension.

After (1) is fixed and (2)–(3) are softened, the paper goes to GREEN. The core packing→$(\ell, m)$ result is well-defended.

## What I could NOT verify (hand to a human expert)

1. **$\Omega(0) = 2$ at $p = 0$.** Paper 0 §VI.C and Paper 18 §VII both state "$1/16 = 1/\Omega^4(0)$ where $\Omega(0) = 2$ is the stereographic conformal factor at $p = 0$." But Paper 18 also writes $\Omega = (1 + p^2/p_0^2)^{-1}$, which gives $\Omega(0) = 1$, not 2. There may be a normalization convention difference (e.g. $\Omega = 2/(1 + p^2/p_0^2)$ in the Fock-stereographic literature), or the "$\Omega(0) = 2$" may refer to a different conformal factor (e.g. the $p_0$-dependent prefactor that maps the $\mathbb{R}^3$ measure to the $S^3$ measure). Hand to expert: verify the precise normalization that makes $1/\Omega^4(0) = 1/16$ consistent. This is a paper-18-and-7 issue, not a paper-0-only issue, but Paper 0 inherits it.

2. **Independent reproduction of $\lambda_{\max} \to 8$.** Paper 18 §VII reports $\lambda_{\max}(5) = 6.62$, $\lambda_{\max}(10) = 7.61$, $\lambda_{\max}(20) = 7.90$. I did not re-run this computation. A spectral graph theorist could verify in 10 minutes; recommended for broadcast.

3. **Conformal-projection equivalence (Paper 7's 18 symbolic proofs).** These are tests of GeoVac code; an external referee should re-derive at least one of the 18 symbolic identities independently before broadcast.

4. **The "factor of 2 from orientability of $S^2$" identification with spin-$\tfrac{1}{2}$ multiplicity.** Paper 0 honestly labels this a "structural correspondence, not a derivation," but a mathematical physicist may want to see the connection (or lack thereof) to spin-structure theory on $S^2$ stated explicitly. Mild D verdict — needs a domain expert to weigh.

## Recommended edits (prioritized)

**HIGH** (fix before any broadcast):
- §V, line ~509: `94.1\%` → `96.0\%` (Paper 15 figure).

**MEDIUM** (fix in errata batch):
- Abstract: disambiguate "shell capacities $2k-1$" → "per-shell angular capacities $2k-1$".
- Abstract: soften "spectrally converge to the hydrogen eigenvalues" → "ground-state eigenvalue spectrally converges to the hydrogen ground-state energy" OR add a footnote that the full spectrum $E_n$ requires the node-weight term $W$ (per Paper 18 §III).
- §V.A: soften "The lattice is the invariant; the continuum geometries are projections" to a dual-reading statement.
- §VI.C: add a one-line note that "ground-state eigenvalue ... converges to $-0.5$ Ha" uses $H = K_{\mathrm{vac}}\mathcal{L}$ alone; the full per-shell spectrum $E_n$ requires node weights $W$ (cross-reference Paper 18 §III).

**LOW** (cleanup batch):
- §IX Conclusions: `~6%` → `~5%` or `5.3%` for LiH.
- §II.B Step 1: minor wording on $d_0$ as shell radius vs separation.
- §VI.C: replace "18 symbolic proofs validate" with "18 internal symbolic test cases validate" to signal internal-consistency.
- Abstract last sentence: "foundational motif of the framework" → "a foundational motif of the framework."
- Hand to expert: verify the $\Omega(0) = 2$ normalization convention (this is a Paper-7/18 issue inherited by Paper 0).
