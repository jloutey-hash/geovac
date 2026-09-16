# Sturmian / exponential-basis paths to chemical accuracy — literature scan (2026-09-13)

**Canonical memo for this scan.** Working trails: `sturmian_accuracy_paths_{A,B,C,D}.md`
(one per axis; search queries, hits, dead ends). PI question (2026-09-13): *"Let's send
out a few agents to see if there's any paths among the Sturmian folks that may help
lead us to chemical accuracy."*

Method: four parallel Opus scans, web-enabled, each capped at 12 verified sources and
25 web calls, each carrying the corpus walls (scale lock, Löwdin / l-selection,
three-centre genus 1, PK, the failed ledger) and a GO / BORDERLINE / STOP gate per
candidate. Axes: **A** molecular Sturmian calculation record; **B** explicit
correlation + hyperspherical; **C** STO / B-function integral technology + production
ETO codes; **D** two-focus (spheroidal) and multi-scale Sturmians.

The PM re-verified at source the two claims that bear on a recommendation or contradict
the corpus (Tao–McCurdy–Rescigno 2010 via the OSTI accepted manuscript, text-extracted;
Hoggan arXiv:1010.5425 abstract) and checked Paper 8's Remark, Paper 12's basis
definition and CLAUDE.md §3.5 against the agents' readings.

---

## 1. Verdict

**No published exponential-type-basis calculation on a molecule (two or more centres)
reaches chemical accuracy through closed-form integrals.** The field's wall, read off
all four axes at once:

| centres | best correlated reach | who |
|:--|:--|:--|
| 1 | 4 electrons (Hy-CI; 4-electron integral is the published bottleneck) | Sims–Hagstrom |
| 2 | 2 electrons, ~1e-15 Ha on H₂; Pachucki's own reach claim is "an arbitrary diatomic" | Pachucki 2010/2012 |
| 3+ | nothing analytic; accuracy bought by a locality device (ADF PARI-MP2, MAD < 1 kcal/mol) or stochastically (Caffarel zero-variance MC, CH₄ near-FCI, 3 µHa) | Förster et al. 2020; Caffarel 2019 |

Two corollaries the PI should hold:

- **GeoVac is not behind the field at three centres.** The genus-1 wall of Paper 59
  is where everyone stops analytically; the ETO community never even posed the
  analytic-class question (§4 below).
- **"Symmetry sparsity only" is the binding difference.** Every ETO route to
  production molecular accuracy buys it with distance sparsity (pair-atomic fitting,
  Schwarz screening) — a sparsity kind the framework does not have and π-free closed
  forms cannot supply.

---

## 2. Candidate paths, ranked

All four keep the **diatomic cap**; none reaches polyatomics.

### Path 1 — Rebuild the Paper 12 prolate-spheroidal two-electron CI at literature resolution (axis D; highest confidence)

**External fact, VERIFIED at source (PM, OSTI PDF text):** Tao, McCurdy, Rescigno,
*Phys. Rev. A* **82**, 023423 (2010), DOI 10.1103/PhysRevA.82.023423 — "we retained
terms up to lmax = 6. The ground-state … was −1.17442 hartree, in excellent agreement
with the accurate value of −1.17447 hartree results of Wolniewicz [26]. We note that
calculations in spherical coordinates reported in ref. [9] using a single-center
expansion with lmax = 7 give a target energy of −1.16908 hartree." Same coordinate
system as Paper 12, FEM/DVR radial basis, electron–electron interaction by a
Neumann-type expansion with closed-form angular factors.

**Corpus state (verified):** Paper 12's basis is Hylleraas-type in prolate spheroidal
coordinates, ξ^j ξ^k η^l η^m e^{−α(ξ₁+ξ₂)}, **one common α optimized variationally**,
maximum powers j = 3, l = 3 at N = 72; the Neumann result "plateaus at 92.4%" of D_e
(≈ 13 mHa short) and the abstract diagnoses "a one-electron basis completeness limit:
the electron-electron cusp requires non-analytic terms (r₁₂^{1/2}, r₁₂ ln r₁₂) that no
polynomial prolate spheroidal basis can represent."

**Contradiction:** a polynomial-in-η, DVR-in-ξ basis at l_max = 6 reaches 0.05 mHa on
the same problem. The 13 mHa plateau is therefore **radial incompleteness (one α,
low powers) plus angular truncation at 3**, not a cusp wall; the cusp is the slow
partial-wave tail everyone lives with, ~110× smaller than Paper 12's residual at
comparable order.

**What it keeps:** discrete (l, m) labels, Gaunt-algebraic angular couplings, the
Neumann V_ee (Paper 12), Paper 11's algebraic Laguerre radial machinery, no
overcompleteness by construction (one two-focus set). No guardrail fires.

**Cost / deltas:** (i) angular order 3 → 6; (ii) a radial set that genuinely spans —
several exponents or a DVR-like Laguerre set — instead of one α. **Gate:** reproduce
−1.1744 within 0.5 mHa; anything above 2 mHa at l_max = 6 says the corpus's radial
machinery, not the geometry, is the limit.

**Honest cap:** diatomic, two electrons. The payoff is a chemically accurate,
quadrature-free, discrete-angular, non-overcomplete two-electron molecular Hamiltonian
as the qubit-encoding anchor in place of a 92.4% D_e result — not polyatomics.

### Path 2 — Gill's Coulomb resolution on the GeoVac basis (axis C; cheapest decisive test)

**External fact, VERIFIED at source (PM, arXiv abstract):** Hoggan, arXiv:1010.5425
(2010): "Coulomb resolutions provide an excellent approximation that reduces these
integrals to a sum of one-electron overlap-like integral products that each involve
orbitals on at most two centers. Such two-center integrals are separable in prolate
spheroidal co-ordinates."

**Why it is interesting:** it attacks the three-centre two-body block (Poly-2 wall)
without surrendering what the corpus owns — the factors land in exactly the two-centre
class where weight-one, π-free closed forms in {exp, E₁, log, γ} exist, and Gaunt
selection survives in the one-electron factors (unlike a Löwdin retrofit).

**Gate the corpus can run from prior measurement:** the resolution needs an auxiliary
set complete in the molecular metric, and v5.11.2 measured the one-centre set as NOT
complete there (Bessel deficit plateau 0.38 / 0.70 / 0.91 at kR = 1 / 2 / 4). Convert
that deficit into a resolution error on a single (XY|XZ) against the Paper 59
momentum-space reference (0.20494172). Above ~1e-6 the path closes on the corpus's own
data; below, the block opens through closed forms already in hand. Expected outcome:
closes. Worth the day because the scalar is decisive either way.

**Honest cap:** it is an approximation (10–20 resolution terms to chemical accuracy in
Hoggan's SCF tests); exactness and Lindemann decidability are forfeited for the
three-centre block.

### Path 3 — Hy-CI restriction (≤ one r₁₂ per configuration) at two centres (axis B)

The only device in the explicitly correlated literature that buys correlated accuracy
while **capping the integral order**: it never generates the three-body operator the
corpus measured as non-collapsing under Gaunt/6j, and it stays Hermitian (the ledgered
TC attempts were non-Hermitian). The hard integral is published — Pachucki, *PRA* **86**,
052514 (2012), two-centre two-electron exponential integrals with integer r₁₂ powers —
and the corpus already has the one-centre instance working (R12-CI He, 0.80 mHa, ~6
functions). Untraveled step = the second centre, not the ansatz.

**Honest cap:** the r₁₂ block has no Gaunt closure — a dense correction bolted onto a
sparse Hamiltonian. Deliverable is chemical accuracy on diatomics as *validation*, not a
sparsity claim; inherits the field's 2e / 2-centre ceiling.

### Path 4 — Split-shell / two-scale Sturmians at two centres (axis A)

Both canons locate the He floor in the same place (one exponent per Goscinskian shell;
Avery's split-shell 1s1s′ mechanism; Herbst's "basis sets with multiple ks?" listed as
unanswered) and neither tested the fix in a molecule. GeoVac's atomic free-per-shell
λ already works. At two centres the pure-number table in s = kR becomes 2-D in
(k₁R, k₂R) — bounded cost — but the locked posing's cancellation
(β_μ − β_ν)⟨Φ_μ|V₀|Φ_ν⟩ = 0 assumes one β per configuration, which a split shell
breaks.

**Honest cap:** most probable outcome is a *quantified* restatement of the known trade
(accuracy ⇔ metric returns), not an escape. Measured nowhere by anyone.

### Recorded non-paths

- **F12/R12** is not an exponential-basis precedent: the Slater-type geminal is
  Gaussian-fitted and orbitals + CABS are Gaussian. Do not cite it as such.
- **Transcorrelated / Jastrow** — ledgered (three-body operator; non-Hermitian).
- **Mixed-scale per-centre k (k_A ≠ k_B)** — **no literature record at all**; every
  implementation reached uses one shared exponent. Paper 60's "Avery's documented
  route" should read "documented as a direction".
- **Avery's product resolution at doubled exponent 2k** ("like exact density fitting",
  *Adv. Quantum Chem.* **70**, 265 (2015)) — nameable, not walled, but mixes l in the
  auxiliary index and buys decidability, not accuracy.
- **ADF PARI-MP2** (Förster, Franchini, van Lenthe, Visscher, *JCTC* **16**, 875 (2020),
  VERIFIED body) — GO at MP2 level, but the enabling move is locality sparsity, which
  the framework structurally lacks.
- **Caffarel zero-variance MC** (arXiv:1906.04515 / *JCP* 2019, VERIFIED full text) —
  the correct way to de-bias a Gaussian evaluator; keeps the symmetry-zero pattern,
  forfeits determinism and π-freeness. Runner-up for bias removal, not structure.

---

## 3. Corpus corrections surfaced (route each; none is silent)

| # | Locus | Defect | Status |
|:--|:--|:--|:--|
| 1 | CLAUDE.md §3.5 guardrail row | "it binds (Avery SW closed forms; Herbst-Avery-Dreuw)" credits HAD 2019 with molecular binding; the paper is atoms-only HF. Paper 8's Remark is correct ("there for atoms"). | **FIXED 2026-09-13** (mechanical, one clause) |
| 2 | Paper 12 abstract + §Convergence | Cusp diagnosis of the 92.4% plateau contradicted by TMR 2010 (0.05 mHa, same coordinates, l_max = 6). | **OWED — PI call** (Remark citing TMR; soften abstract/conclusion per the summary-surface rule) |
| 3 | Paper 15 abstract | "exceeding the 92.4% … (Paper 12)" is true but the published two-focus ceiling is 0.05 mHa; the "angular basis is the bottleneck" reading needs the coordinate-choice caveat (110× between spherical l_max=7 and prolate l_max=6 in TMR). | OWED |
| 4 | Paper 60 §molecular | "Avery's documented route" for mixed scale → "documented as a direction" (no worked mixed-scale molecular calculation exists). | OWED |
| 5 | Paper 58/59 + memory | "Prolate spheroidal has two foci; Slater functions have no product theorem" is asserted unbacked; Hoggan states it verbatim — cite. | OWED (bibitem) |
| 6 | Polyatomic wall framing | Poly-2 three-centre block is a **closed-form** wall, not an **accuracy** wall (Caffarel reaches 2e-9 integrals stochastically). Say "no closed form and no cheap deterministic route". | memory FIXED; papers OWED |
| 7 | `noci_engine` evaluator | Molpro's SMILES manual warns STO-9G integrals "may have not sufficient accuracy for post-HF calculations" — the corpus's evaluator is the same construction. Check n_gauss adequacy on one post-HF number. | OWED diagnostic |
| 8 | Paper 14/17 benchmarking rule | Strongest LiH baseline is ECG, 0.3 cm⁻¹ (Tung, Pavanello, Adamowicz, *JCP* **134**, 064117 (2011), UNVERIFIED at source), not STO-3G. | OWED |
| 9 | Paper 60 bibliography | Klahn & Bingel, *IJQC* **11**, 943 (1977) overcompleteness taxonomy — citation debt (UNVERIFIED abstract). | OWED |
| 10 | Paper 59 | Prior-art test of "connection unmade": **survives**. Weniger, Safouhi, Hoggan full-text sweeps: 0 hits for transcendence / Appell / Lauricella / Meijer; "elliptic" in that literature = elliptical coordinates. Keep the hedge; do not strengthen. | no change |

---

## 4. Verification ledger

**PM-verified at source:** TMR 2010 (OSTI accepted manuscript, pdftotext; sentences
quoted above); Hoggan arXiv:1010.5425 (abstract); Paper 8 Remark `rem:overlap_imposed`;
Paper 12 Eq. (basis_function) + §Convergence (j = 3, l = 3, N = 72, common α);
CLAUDE.md §3.5 row.

**Agent-verified at source (full text or abstract):** Herbst–Avery–Dreuw PRA 99, 012512
(full text, A); Avery & Avery JPCA 113, 14565 (abstract, A/D); Pachucki PRA 82, 032509
and PRA 86, 052514 (B); Padhy arXiv:1609.00269 (B); Tolstikhin–Watanabe–Matsuzawa PRL
74, 3573 (metadata, B); Weniger arXiv:0811.3406 (full text, C); Slevinsky–Safouhi
arXiv:1905.13537 (abstract, C); Molpro SMILES manual (C); Förster et al. JCTC 16, 875
(body, C); Caffarel arXiv:1906.04515 (full text, C); Lehtola IJQC 119, e25944 (D);
Mitnik–López–Ancarani Mol. Phys. 119, e1881179 (D); Kereselidze et al. Mol. Phys. 113,
3471 (abstract, D).

**UNVERIFIED (must be reached at the publisher before entering any paper):**
Avery & Avery Mol. Phys. 110, 1593 (2012); Aquilanti et al. THEOCHEM 2004; Klahn–Bingel
IJQC 1977; Sims–Hagstrom 2004/2015; Ten-no 2004/2012; Cohen et al. JCP 151, 061101;
Tung–Pavanello–Adamowicz JCP 134, 064117; Duchon–Dumont-Lepage–Gazeau JCP 76, 445 (DOI);
Vanne–Saenz J. Phys. B 2004; Safouhi J. Comput. Phys. 176, 1 (content).

Access note (from the trails): APS, AIP, Wiley, T&F, ScienceDirect and PubMed 403 the
fetcher; Crossref REST, EuropePMC REST, OSTI PURLs and arXiv work; binary PDFs saved by
the fetcher are recoverable with `pdftotext -layout`.

---

## 5. Answer to the PI's question, in one paragraph

There is no hidden Sturmian trick that yields chemical accuracy on molecules while
keeping the framework's sparsity; the literature confirms the corpus's own wall map,
and the ETO community reached production accuracy only by importing locality or
stochastics. What the scan did find is that **the framework's own Level-2 geometry has
a published chemical-accuracy result the corpus never cited**, and that Paper 12's
diagnosis of its 92.4% plateau is wrong: the two-focus partial-wave basis reaches
0.05 mHa on H₂ at l_max = 6 when the radial part genuinely spans. That is a
diatomic-only, two-electron result, but it is the one path here that keeps discrete
angular labels, algebraic couplings and no overcompleteness — and it is half built.
The Coulomb-resolution scalar test is the second thing worth a day, mainly because it
is decisive on the corpus's own data.
