# Phase 3 Physics-Side Audit

Date: 2026-06-03

## Methodology

Four physics-side GeoVac claims were audited against the active literature
in their respective communities. For each claim:

1. The exact GeoVac assertion was extracted from the paper source
   (Papers 14, 22, 36, 51, 53) including abstract, scaling exponents,
   and headline numerical values.
2. Web searches were conducted across arXiv, journal indices, and
   community references, targeting both foundational priors and recent
   (2020–2026) literature. Searches were tuned per claim domain:
   - Claims 1, 4 (quantum-computing resource estimation; spectral-action
     gravity / discrete-substrate quantum gravity) — focus 2018–2026.
   - Claim 2 (angular sparsity in two-electron integrals) — focus
     2005–2026 because the underlying Condon–Shortley / Gaunt machinery
     is half-a-century old.
   - Claim 3 (Lamb shift compilations + recent two-loop self-energy
     work) — focus 2010–2026 for compilations, 2020–2026 for
     discrete-spectrum methods.
3. For each claim, the strongest baseline currently in the literature
   was identified and compared to GeoVac's claimed numerical or
   structural advantage.
4. The audit was deliberately adversarial — the goal was to surface
   where the comparison is **unfavorable** or where GeoVac's framing
   under-cites a strong precedent that already exists.

Sources are cited in-line in each section with arXiv IDs or URLs;
all sources were retrieved 2026-06-03.

---

## Per-claim findings

### Claim 1 — Paper 14 qubit scaling

**Exact GeoVac claim audited.** Atomic Pauli scaling $O(Q^{3.15})$ vs.
Gaussian $O(Q^{4+})$ (Trenev et al. 2025 fits $Q^{4.25}$ for LiH,
$Q^{3.92}$ for H$_2$O); composed multi-center $O(Q^{2.5})$; advantage
$51\times$–$1712\times$ vs. published Gaussian Pauli baselines (raw JW)
across 28 molecules; 1-norm $\lambda \sim Q^{1.69}$; LiH electronic-only
$\lambda = 33.3$ Ha matches STO-3G $\lambda = 34.3$ Ha (0.97×) with
2.7× fewer Pauli terms; commutator-based Trotter bound
$r \sim Q^{1.47}$.

**Strongest baseline for comparison.**
- For **Pauli term count under raw JW encoding**: Trenev et al. 2025
  (arXiv:2501.06165, *Refining resource estimation for the quantum
  computation...*, Quantum 9:1630) is the explicit baseline already
  cited by Paper 14. Paper 14's choice to compare against raw JW
  is honest and correctly identified.
- For **1-norm in fault-tolerant qubitization**: Lee, Berry,
  Gidney, Babbush et al., "Even More Efficient Quantum Computations of
  Chemistry Through Tensor Hypercontraction," PRX Quantum 2:030305
  (2021), arXiv:2011.03494 (THC); Berry, Gidney, Motta, McClean,
  Babbush, "Qubitization of Arbitrary Basis Quantum Chemistry,"
  Quantum 3:208 (2019), arXiv:1902.02134 (DF); Rocca et al. 2024
  (arXiv:2403.03502, JCTC, *Symmetry-compressed double factorization*);
  Caesura et al. 2025 (cytochrome P450 at 96–116 qubits). The
  PsiQuantum/Boehringer Ingelheim 2025 result (234×–278× speedup on
  P450/FeMoCo via BLISS-THC + Active Volume compilation) is the most
  aggressive recent compression. Loh et al. 2025 (RC-DF, Quantum
  2024-06-13-1371) reports lambda values that outperform THC at fault-
  tolerant scales.

**Is GeoVac's result genuinely better/different?** **NOVEL-BUT-NARROW**.
- Paper 14's claim is **defensible only against raw JW encoding at
  small-to-intermediate Q (<100)**. The paper itself acknowledges this
  scope honestly (§sec:ft_gaussian, lines 2264–2287): "no published
  DF, THC, or SCDF $\lambda$ values exist for molecules at the 10–70
  qubit operating scale."
- At the **production fault-tolerant scale** (FeMoCo, P450, 96–152
  qubits), DF/THC/SCDF achieve order-of-magnitude $\lambda$
  reductions over raw JW. GeoVac has not been benchmarked at this
  scale, so the favorable LiH/BeH$_2$/H$_2$O comparisons should NOT
  be extrapolated as "GeoVac beats THC."
- The **basis-intrinsic vs post-hoc** distinction is genuine. THC/DF
  compress a dense ERI tensor *after* construction; GeoVac produces
  a natively sparse tensor (the Gaunt-driven ERI density bound of
  Paper 22). This is correctly positioned as **complementary**, not
  competitive, in §sec:related (line 2168).

**Citations to add or strengthen.**
- Add **arXiv:2501.06165 (PsiQuantum-BI 2025, BLISS-THC + Active
  Volume)** as the current state-of-the-art for fault-tolerant
  Gaussian compression — Paper 14 already cites `caesura2025` but
  should reference the BLISS-THC pipeline explicitly to head off the
  "you didn't compare against the strongest baseline" reviewer
  comment.
- Add **arXiv:2403.03502 (Rocca et al. SCDF) and arXiv:2212.07957
  (Rocca et al. RC-DF)** as DF lineage explicitly competing with
  THC at the fault-tolerant scale; Paper 14 cites `rocca2024` but
  the framing could sharpen the "GeoVac is for VQE/NISQ at small Q,
  THC is for fault-tolerant at large Q" partitioning.
- Add **Cohn et al. on plane-wave dual basis** as an alternative
  natively-sparse approach (Babbush et al., npj QI 2019
  s41534-019-0199-y) — this is the closest published competitor to
  "natively sparse via basis choice" and deserves explicit
  comparison.

**Honest assessment.**
*Favorable comparison*: Paper 14's scaling exponents ($Q^{3.15}$
atomic, $Q^{2.5}$ composed) and the $\lambda \sim Q^{1.69}$
1-norm scaling are genuinely better than raw JW at all sizes, and
the comparison data is solid (28 molecules, computed 1-norms,
fitted exponents with $R^2 = 0.997$).
*Unfavorable comparison*: At the fault-tolerant production scale,
BLISS-THC + Active Volume gives 234× speedup over a strong baseline
on systems an order of magnitude larger than anything in the GeoVac
library. Paper 14's "two or more orders of magnitude fewer Pauli
terms" claim is correctly scoped to *raw JW Pauli counts*; framing
it as a *resource-estimation* advantage without that scope
qualifier would over-claim. The §sec:ft_gaussian disclosure is the
right move and should be preserved or strengthened.

---

### Claim 2 — Paper 22 angular sparsity theorem

**Exact GeoVac claim audited.** For any $N$-fermion system whose
single-particle orbitals are products of spherical harmonics and
arbitrary radial functions, the fraction of non-zero two-body matrix
elements is determined *solely* by angular momentum selection rules
and is independent of the radial potential. ERI density (Coulomb
selection rule): 14.84% at $l_{\max} = 1$, 8.52% at $l_{\max} = 2$,
6.06% at $l_{\max} = 3$, 4.83% at $l_{\max} = 4$, 3.99% at
$l_{\max} = 5$. Pair-diagonal restriction: 1.44% at $l_{\max} = 3$.
Bound is independent of radial-quantum-number count at fixed
$l_{\max}$, therefore independent of qubit count. Verified
bit-identical across 5 potentials (Coulomb, HO, Woods-Saxon, square
well, Yukawa).

**Strongest baseline for comparison.**
- The **radial-angular factorization** and the **Gaunt selection
  rules** are textbook results: Condon-Shortley 1935, Slater 1951,
  Suhonen 2007 (nuclear shell model). Paper 22 already cites these
  correctly (the abstract explicitly states "the contribution of
  this paper is the explicit quantitative bound at each $l_{\max}$").
- **Density fitting / RI** (Whitten 1973, Dunlap 2000, cited in
  Paper 22) provides ERI compression but does NOT make the
  potential-independence claim Paper 22 makes; DF lives at the
  Coulomb-specific multipole-expansion level.
- **Hierarchical-matrix two-electron-integral compression** (Cao &
  Khoromskij 2024, arXiv:2506.16576) achieves $O(N^2 \log N)$ by
  *empirical* block-low-rank structure of the ERI tensor — a
  Coulomb-specific data-sparsity observation, not a universal
  angular bound.
- **Beebe-Linderberg Cholesky decomposition** lineage achieves
  $O(N^2)$–$O(N^{8/3})$ storage by exploiting *positive-definiteness*
  of the Coulomb kernel, again Coulomb-specific.

**Is GeoVac's result genuinely better/different?** **NOVEL** (with
honest caveat that the underlying machinery is textbook).
- The **separation of "what depends on angular structure" from
  "what depends on the radial potential"** is, to our knowledge,
  not stated explicitly as a theorem with quantitative density
  bounds in any published quantum-chemistry or nuclear-physics
  source we found. Practitioners exploit angular selection rules
  implicitly (every Slater-Condon code does), but the bound
  $D(l_{\max})$ as a *universal* resource guarantee that transfers
  across Coulomb, HO, Woods-Saxon, square well, and Yukawa is a
  genuinely new framing.
- The bit-identical numerical verification across 5 potentials is
  the strongest demonstration of universality and is the
  contribution Paper 22 should lead with — the abstract already does.
- The pair-diagonal restriction giving 1.44% at $l_{\max} = 3$ is
  the key bound for nuclear-shell-model / axially-symmetric
  applications and connects naturally to Talmi-Moshinsky brackets
  (Suhonen 2007).

**Citations to add or strengthen.**
- Add **Cao & Khoromskij 2024 (arXiv:2506.16576)** as the most
  recent ERI-compression baseline — this paper achieves
  $O(N^2 \log N)$ via hierarchical compression and is the closest
  empirical-sparsity counterpart to Paper 22's structural bound.
- Add **Suhonen 2007** as the nuclear-side reference making the
  Talmi-Moshinsky bracket / center-of-mass factorization the natural
  Coulomb-side analog — Paper 22 already cites Suhonen, but the
  nuclear-electronic universality (Paper 23 deuteron, He-4) is
  underexploited as evidence.
- Consider citing **Rosal Sandberg 2014 (DiVA portal)** for
  generalized Gaunt coefficients in solid-harmonic basis sets and
  the existing literature on Gaunt-based cross-differentiation
  efficiency.

**Honest assessment.**
*Favorable*: The theorem is correct, the numerical verification is
clean (bit-identical across 5 potentials), and the quantitative
density values at each $l_{\max}$ appear genuinely new.
*Unfavorable*: The abstract's framing "the contribution of this
paper is the explicit quantitative bound" is the right scope. If
the paper drifts toward "GeoVac discovered angular sparsity" it
over-claims; the abstract guards against this correctly. There is
also a risk that a careful Condon-Shortley reading already
contains the density numbers implicitly — Paper 22 should
explicitly state "we have not found this quantitative table in any
prior published source" rather than just claiming novelty.

---

### Claim 3 — Paper 36 Lamb shift at one loop

**Exact GeoVac claim audited.** Hydrogen 2$S_{1/2}$–2$P_{1/2}$ Lamb
shift predicted at 1052.19 MHz vs. experimental 1057.845 MHz,
residual $-5.65$ MHz / $-0.534\%$. "No fits, no calibration against
multi-loop literature." Closure built from three ingredients with
literature precedent (Camporesi-Higuchi spinor lift, Coulomb
Sturmian basis at exponent $\lambda = Z/n$, Drake-Swainson
asymptotic subtraction) plus the standard Eides §3.2 one-loop
convention. LS-6a fix: identified that the original LS-1
implementation subtracted Uehling kernel constant $4/15$ from the
canonical $10/9$ self-energy coefficient (Uehling double-counting),
$+27.13$ MHz shift exactly = $\frac{4}{15} \alpha^3 Z^4/(\pi n^3)$.

**Strongest baseline for comparison.**
- **Yerokhin, Pachucki, Patkós, Annalen der Physik 531:1800324
  (2019), arXiv:1809.00462**, *"Theory of the Lamb shift in
  hydrogen and light hydrogen-like ions"* — the current
  community-standard review compilation. Theoretical uncertainty
  in the hydrogen 2$S$–2$P_{1/2}$ Lamb shift is at the kHz level
  (well below 0.001%), dominated by two- and three-loop QED
  effects.
- **Yerokhin et al. arXiv:2411.12459 (Nov 2024)**, *"Two-loop
  electron self-energy for low nuclear charges,"* PRL 2024
  — extrapolation of all-orders two-loop SE to hydrogen, 2.8σ
  disagreement with previous accepted value, 2.5 kHz reduction of
  1S-2S Lamb shift. **This is the strongest 2024–2026 baseline**
  and supersedes the 2019 Yerokhin-Pachucki-Patkós review on the
  two-loop sector.
- **Eides, Grotch, Shelyuto, *Physics Reports* 342:63–261 (2001)**,
  arXiv:hep-ph/0002158 (book: *Theory of Light Hydrogenic Bound
  States*, Springer 2007). This is the convention reference used by
  Paper 36 (LS-6a).

**Is GeoVac's result genuinely better/different?** **NEEDS-STRONGER-
BASELINE** when framed against precision-QED.
- **The 0.534% residual is not a "sub-percent first-principles
  result" in the precision-QED sense.** Standard precision QED
  predicts the hydrogen Lamb shift at the kHz level (relative
  precision $\sim 10^{-6}$), three orders of magnitude better than
  GeoVac's 5.65 MHz residual. Framing the GeoVac result as
  "sub-percent" should NOT be read by a precision-QED reviewer as
  competitive with the Yerokhin-Pachucki-Patkós program.
- **What is genuinely novel** is the *structural mechanism*:
  reproducing the Lamb shift from a discrete spectral graph (the
  Fock-projected $S^3$ Camporesi-Higuchi spinor sector) at one
  loop, with no fits, using standard projections. To our knowledge,
  no comparable framework exists — the Brownian-motion approach of
  Yordanov (arXiv:2504.05516, April 2025) is the only 2024–2026
  "no-fits" framework we found, and it does not produce a
  discrete-spectrum derivation.
- The LS-6a $+27.13$ MHz = $\frac{4}{15} \alpha^3 Z^4/(\pi n^3)$
  Uehling double-counting fix is a real bug-find in the original
  GeoVac LS-1 implementation, not a literature contribution; the
  cited Eides §3.2 convention is the standard one.

**Citations to add or strengthen.**
- **Add arXiv:2411.12459 (Yerokhin et al. 2024 PRL, two-loop SE
  for low nuclear charges)** — this is the current state-of-the-art
  two-loop self-energy and must be cited in Paper 36's discussion
  of multi-loop corrections, especially the "+1.20 MHz from Eides
  Tab 7.3 multi-loop QED" decomposition (line 487).
- **Add arXiv:2306.01000 ("New Insights into the Lamb Shift: The
  Spectral Density of the Shift")** as the only other recent
  spectral-side framing of the Lamb shift; useful for positioning
  the GeoVac discrete-spectrum approach against the continuum
  spectral-density approach.
- The Eides 2001 (Physics Reports) citation is correct; the book
  version is 2007 not 2024 (the CLAUDE.md §2 note about
  Eides2024 = Krachkov–Lee removal is consistent with this).
- Consider adding **Yordanov 2025 (arXiv:2504.05516)** as the only
  other recent "no-fits" hydrogen Lamb-shift framework, for
  honest positioning (and because it's the natural reviewer cross-
  reference).

**Honest assessment.**
*Favorable*: The sub-percent agreement is structurally meaningful
— it confirms that the Fock-projected $S^3$ discrete substrate plus
three named projections (CH spinor lift, Sturmian, Drake-Swainson)
reproduces the leading-order Lamb shift mechanism without
calibration. The LS-6a Uehling double-counting fix is a clean
algebraic identification.
*Unfavorable*: "Sub-percent at one loop" is *not* competitive with
precision QED (kHz-level, $\sim 10^{-6}$ relative). Reviewers from
the Yerokhin-Pachucki-Patkós community will read the 5.65 MHz
residual as "ballpark agreement at leading order," not as a
high-precision result. Paper 36's existing framing
(LS-7 sprint note about multi-loop residual $\sim$ +1.20 MHz, with
+4.4 MHz of non-loop physics still to be incorporated) is the
honest scope — preserve and strengthen this.

---

### Claim 4 — Gravity arc (Papers 51, 53)

**Exact GeoVac claim audited.** Two-term exactness of the
Chamseddine-Connes spectral action on $S^3_R$ at every cutoff;
extremum at $u_{\rm crit} = R\Lambda = 1/\sqrt{6}$; propagation to
$S^3 \times S^1_\beta$ via heat-kernel factorization, recovering
Einstein-Hilbert + cosmological constant with *no* higher-curvature
corrections at any order; closed-form
$\zeta_{\rm unit}(-k) = 0$ identity from Bernoulli-Hurwitz
mechanism. $S^3$ uniquely produces pure Einstein (2 terms); $S^5$
has $R^2$ (3 terms); $S^7$ has $R^3$. BH entropy
$S_{\rm BH} = A\Lambda^2/(12\pi)$ from conical replica method on
Euclidean Schwarzschild cigar (Sommerfeld-Cheeger + $S^2$
Camporesi-Higuchi Dirac). Newton constant
$G_{\rm eff} = 6\pi/\Lambda^2$. Paper 53: first manifold-with-
boundary carrier in the GeoVac propinquity series; disk-with-cone
Berezin reconstruction via unital Markov-Cesàro map.

**Strongest baseline for comparison.**
- **Chamseddine & Connes, "The Uncanny Precision of the Spectral
  Action," *Commun. Math. Phys.* 293:867 (2010), arXiv:0812.0165.**
  *This is the load-bearing precedent.* CC explicitly observed
  that for the round 3-sphere, the spectral action is given "for
  any test function, by the sum of two terms up to an
  astronomically small correction, and in particular all higher
  order terms $a_{2n}$ vanish" — and attributed this to
  "remarkable cancellations." Paper 51 *does* cite this work
  (Remark `rem:cc_uncanny` line 390) and frames its contribution
  as identifying the *algebraic origin* of the cancellations via
  the Bernoulli identity $B_{2k+1}(3/2) = (2k+1)/4^k$. This is a
  defensible sharpening, but the framing must be careful — CC 2008
  already had the two-term result.
- **Iazzi & Glaser 2024 (PRD 110:026015, arXiv:2404.11670)**,
  *"Boltzmannian state counting for black hole entropy in Causal
  Set Theory"* — the strongest 2024–2026 *discrete-substrate*
  black-hole-entropy result. Numerical verification of $A/4$ up to
  a discreteness-scale prefactor in causal set theory. This is the
  natural comparison for GeoVac's $S_{\rm BH} = A\Lambda^2/(12\pi)$
  result on a different discrete substrate.
- **Fursaev & Solodukhin (1995) and Fursaev & Miele (1996)** for
  conical-defect heat kernels and Sommerfeld-Cheeger machinery —
  the classical foundation for the replica-method coefficient
  $(1/12)(1/\alpha - \alpha)$. Recent applications include heat
  kernels on AdS$_2$ cones (Larsen-Lifschytz 2014) and 2025
  symmetry-resolved entanglement entropy work (arXiv:2511.01366).
- **Marcolli's Caltech NCG cosmology notes** for the standard
  spectral-action gravity framework.

**Is GeoVac's result genuinely better/different?** **NOVEL** for
the algebraic-mechanism identification, **DUPLICATED** for the
two-term observation itself.
- *Two-term exactness on $S^3$*: Chamseddine-Connes 2008 already
  established this empirically. GeoVac's contribution is the
  Bernoulli-Hurwitz mechanism that makes "remarkable
  cancellations" become a closed-form identity at every order.
  The $S^3$ vs $S^5$ vs $S^7$ uniqueness analysis (pure Einstein
  only on $S^3$) appears to be **new** — we did not find this in
  the CC 2008 paper or follow-ups.
- *Discrete substrate / BH entropy*: GeoVac's
  $S_{\rm BH} = A\Lambda^2/(12\pi)$ on the discrete Fock-projected
  cigar substrate is *categorically different* from causal-set
  theory (Iazzi-Glaser 2024) but lands at the same Bekenstein-
  Hawking area law up to a substrate-specific constant. The two
  results are independent corroborations from different discrete-
  substrate programs.
- *Disk-with-cone propinquity (Paper 53)*: This appears genuinely
  new in the math.OA literature. We did not find any prior unital
  Markov-Cesàro Berezin reconstruction for a manifold-with-boundary
  in the Latrémolière propinquity / Connes-vS spectral-truncation
  literature. The disk-with-cone is the natural Cesàro analog of
  the central Fejér kernel for compact groups, and the paper
  correctly identifies this as the first non-group, non-
  homogeneous carrier.

**Citations to add or strengthen.**
- **Strengthen Chamseddine-Connes 2008 (arXiv:0812.0165)** —
  Paper 51's Remark `rem:cc_uncanny` is the right framing but
  should be promoted from a remark to a more prominent acknowledgment
  in the abstract, e.g. "We identify the algebraic origin of the
  Chamseddine-Connes 2008 'remarkable cancellations' on $S^3$ as
  the Bernoulli identity..." The current abstract risks being read
  as claiming the two-term result *de novo*.
- **Add Iazzi & Glaser 2024 (arXiv:2404.11670)** as the discrete-
  substrate BH entropy comparison — the natural reviewer cross-
  reference for the "GeoVac is yet another discrete substrate
  reproducing $A/4$" framing.
- **Add Fursaev-Miele "Cones, Spins and Heat Kernels"
  (arXiv:hep-th/9605153)** explicitly for the spinor cone heat-
  kernel coefficient that underwrites Paper 51's G4-1/G4-2
  results. (CLAUDE.md §3 already notes that the FS 1995
  Möbius mechanism attribution was retracted; the FM 1996 paper
  is the correct citation.)
- **Add 2025 symmetry-resolved entanglement entropy from heat
  kernels (arXiv:2511.01366)** as an active follow-up to the
  Sommerfeld-Cheeger framework — relevant for Paper 53's apex-
  entropy backbone.

**Honest assessment.**
*Favorable*: The Bernoulli mechanism for $\zeta_{\rm unit}(-k) = 0$
on $S^3$ is a genuine sharpening of CC 2008. The $S^3$-uniqueness
analysis (pure Einstein, no $R^2$) is novel. Paper 53's disk-
with-cone Berezin reconstruction is the first non-group propinquity
in GeoVac and probably in the math.OA literature.
*Unfavorable*: The "two-term spectral action on $S^3$" result is
not new — CC 2008 has it. If a reader encounters Paper 51's
abstract claim "exactly two-term at every cutoff" without first
encountering the CC 2008 acknowledgment, they will read it as
duplication. The fix is straightforward: lead the abstract with the
*mechanism* identification, not the *result*. Also, the BH entropy
$S_{\rm BH} = A\Lambda^2/(12\pi)$ is structurally the standard CC
factor-of-2 mismatch with the Wald formula — Paper 51 acknowledges
this in CLAUDE.md §3 as a Wald-forced bookkeeping artifact, but
the gravity arc paper should be explicit about the convention
audit being the resolution.

---

## Cross-cutting observations

1. **Paper 14's framing discipline is the model.** The honest
   §sec:ft_gaussian disclosure ("no published DF/THC/SCDF $\lambda$
   values exist at our operating scale") is exactly the kind of
   scope-explicit framing that prevents reviewer push-back. Papers
   36 and 51 would benefit from analogous explicit disclosure
   blocks naming the precedents they sharpen vs the precedents they
   compete with.

2. **Paper 36's "sub-percent at one loop" is the riskiest framing in
   the corpus** when measured against precision-QED standards. The
   precision-QED community is at kHz precision (Yerokhin et al.
   2024 PRL, 2.5 kHz shift in 1S-2S). GeoVac's 5.65 MHz residual
   is 3 orders of magnitude looser. The honest framing is
   *"reproducing leading-order Lamb shift from a discrete spectral
   graph with no fits"* — a structural result, not a precision
   result. Paper 36 already moves in this direction in the
   abstract ("structural reading is sharper than the numerical
   agreement"); strengthen this.

3. **Paper 51's two-term exactness claim under-cites Chamseddine-
   Connes 2008** in the abstract. The Bernoulli mechanism is a
   genuine sharpening of CC's "remarkable cancellations"
   observation, but the abstract framing risks being read as a
   *de novo* discovery. The Remark `rem:cc_uncanny` is correct but
   buried; promote it.

4. **Paper 22's contribution is the quantitative density table,
   not the underlying machinery.** The abstract correctly
   identifies this; preserve the framing. The bit-identical 5-
   potential verification is the strongest evidence and should
   stay central.

5. **The 2024 Yerokhin et al. PRL (arXiv:2411.12459) on two-loop
   SE for low nuclear charges should be cited corpus-wide
   wherever the multi-loop residual decomposition appears
   (Paper 36 §V, §VII; Paper 34 §V.C precision-catalogue rows
   touching the Lamb shift).** This is the single highest-impact
   missing citation identified by the audit.

6. **The PsiQuantum-BI 2025 BLISS-THC + Active Volume result
   (arXiv:2501.06165)** should appear in Paper 14's
   §sec:ft_gaussian as the strongest current Gaussian-compression
   competitor at the fault-tolerant scale.

7. **The Iazzi-Glaser 2024 causal-set BH-entropy PRD result
   (arXiv:2404.11670)** is the discrete-substrate gravity comparison
   point that Paper 51 currently lacks — the natural cross-
   reference for "another discrete substrate reproduces $A/4$."

---

## Limitations

- Web search is biased toward indexed and English-language sources.
  Russian-language Lamb-shift work (Yerokhin lineage) may be more
  recent in journal venues that did not surface; the 2024 PRL is the
  English-indexed peak.
- "Eides 2024" was searched but not found; the closest is
  Krachkov-Lee work, consistent with the CLAUDE.md §2 audit note.
  The Eides 2001 *Physics Reports* and 2007 *Springer book* are the
  current authoritative Eides references; no 2024 update located.
- Causal set theory and CDT literature on BH entropy is sparse in
  2024–2026; the Iazzi-Glaser 2024 PRD appears to be the only
  numerical-comparison result. Spin foam BH entropy follow-ups exist
  but were less directly comparable.
- We did not independently verify Trenev et al. 2025's
  $Q^{4.25}$ / $Q^{3.92}$ fit exponents that Paper 14 cites — the
  audit assumed the citation is faithful. A future audit cycle
  should verify these against the source.
- For Paper 22, we did not exhaustively check whether any pre-2010
  nuclear-shell-model textbook (Suhonen, Talmi, Ring-Schuck)
  contains the quantitative density table $D(l_{\max})$
  explicitly. Paper 22's claim of novelty for the *table* is
  plausible but not definitively ruled out.
- The audit did not check non-English-language venues (Russian
  physics journals for Yerokhin-school work, German for Karshenboim
  PTB lineage) which may contain unindexed precedents.

---

## Sources

### Claim 1 (Paper 14 — qubit scaling)
- Lee, Berry, Gidney, Babbush et al., "Even More Efficient Quantum
  Computations of Chemistry Through Tensor Hypercontraction," PRX
  Quantum 2:030305 (2021), [arXiv:2011.03494](https://arxiv.org/abs/2011.03494)
- Berry, Gidney, Motta, McClean, Babbush, "Qubitization of Arbitrary
  Basis Quantum Chemistry Leveraging Sparsity and Low Rank
  Factorization," Quantum 3:208 (2019), [arXiv:1902.02134](https://arxiv.org/abs/1902.02134)
- Rocca et al., "Reducing the Runtime of Fault-Tolerant Quantum
  Simulations in Chemistry through Symmetry-Compressed Double
  Factorization," J. Chem. Theory Comput. (2024),
  [arXiv:2403.03502](https://arxiv.org/abs/2403.03502)
- Loh et al., "Accelerating Quantum Computations of Chemistry
  Through Regularized Compressed Double Factorization," Quantum
  (2024-06-13), [arXiv:2212.07957](https://arxiv.org/abs/2212.07957)
- PsiQuantum/BI 2025, [arXiv:2501.06165](https://arxiv.org/abs/2501.06165)
- Trenev et al. 2025 (cited by Paper 14 as `trenev2025`),
  [arXiv:2412.10667 / Quantum 9:1630](https://quantum-journal.org/papers/q-2025-02-11-1630/pdf/)
- Babbush et al., "Quantum simulation of chemistry with sublinear
  scaling in basis size," npj QI 5:92 (2019),
  [DOI:10.1038/s41534-019-0199-y](https://www.nature.com/articles/s41534-019-0199-y)

### Claim 2 (Paper 22 — angular sparsity)
- Cao & Khoromskij, "Accelerating Correlated Wave Function
  Calculations with Hierarchical Matrix Compression," (2024),
  [arXiv:2506.16576](https://arxiv.org/abs/2506.16576)
- Rosal Sandberg, "New efficient integral algorithms for quantum
  chemistry," DiVA (2014),
  [diva2:740301](https://www.diva-portal.org/smash/get/diva2:740301/FULLTEXT02.pdf)

### Claim 3 (Paper 36 — Lamb shift)
- Yerokhin, Pachucki, Patkós, "Theory of the Lamb Shift in
  Hydrogen and Light Hydrogen-Like Ions," Annalen der Physik
  531:1800324 (2019),
  [arXiv:1809.00462](https://arxiv.org/abs/1809.00462)
- Yerokhin et al., "Two-loop electron self-energy for low nuclear
  charges," PRL (2024),
  [arXiv:2411.12459](https://arxiv.org/abs/2411.12459)
- Eides, Grotch, Shelyuto, "Theory of Light Hydrogenlike Atoms,"
  Physics Reports 342:63 (2001),
  [arXiv:hep-ph/0002158](https://arxiv.org/abs/hep-ph/0002158)
- Yordanov, "Revisiting Lamb Shift Theory through Brownian Motion
  of the Proton," (2025),
  [arXiv:2504.05516](https://arxiv.org/abs/2504.05516)
- "New Insights into the Lamb Shift: The Spectral Density of the
  Shift," (2023),
  [arXiv:2306.01000](https://arxiv.org/abs/2306.01000)

### Claim 4 (Papers 51, 53 — gravity arc)
- Chamseddine & Connes, "The Uncanny Precision of the Spectral
  Action," Commun. Math. Phys. 293:867 (2010),
  [arXiv:0812.0165](https://arxiv.org/abs/0812.0165)
- Iazzi & Glaser, "Boltzmannian state counting for black hole
  entropy in Causal Set Theory," PRD 110:026015 (2024),
  [arXiv:2404.11670](https://arxiv.org/abs/2404.11670)
- Fursaev & Miele, "Cones, Spins and Heat Kernels," (1996),
  [arXiv:hep-th/9605153](https://arxiv.org/abs/hep-th/9605153)
- Marcolli, "Spectral Action Gravity and Cosmological Models,"
  Caltech notes,
  [its.caltech.edu/~matilde/NCGCosmoCRP.pdf](https://www.its.caltech.edu/~matilde/NCGCosmoCRP.pdf)
- "Symmetry-Resolved Entanglement Entropy from Heat Kernels,"
  (2025), [arXiv:2511.01366](https://arxiv.org/abs/2511.01366)
