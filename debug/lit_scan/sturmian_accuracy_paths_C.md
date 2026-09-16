# Literature scan, axis C: STO / B-function / momentum-space multi-centre integral
# technology, and production STO codes

Date: 2026-09-13. Read-only scan. No papers, code or tests edited.
Caps observed: 25 web calls (at cap), 7 sources verified at source (cap 12).

Assigned question: does any exponential-type-orbital (ETO) code reach chemical
accuracy on a **correlated (post-HF) molecule**, how does it handle three- and
four-centre integrals, and what does that cost? Plus: has anyone in the STO /
B-function literature ever stated an **elliptic or higher-transcendence class**
for three-centre integrals (the Paper 59 prior-art test)?

Scope discipline: the molecular-Sturmian record, explicit correlation (R12/F12,
Hy-CI), and two-centre spheroidal Sturmians are other agents' axes and are not
treated here. (Note: the shared WebFetch cache in this session contains PDFs
those agents pulled -- Herbst's Coulomb-Sturmian talk, Mitnik-Lopez-Ancarani
prolate GSF, Ruiz-Sims-Padhy Hy-CI. I did not read them for this memo.)

---

## 1. Result table

| Path / code | Best demonstrated accuracy + system | 3-/4-centre integrals | Enabling move | What GeoVac forfeits | Verdict | Source |
|:--|:--|:--|:--|:--|:--|:--|
| **B functions + Fourier method** (Weniger, Steinborn, Grotendorst, Homeier) | Integral technology only. No correlated molecule reported. | Reduce to a **highly oscillatory semi-infinite spherical-Bessel integral**; Weniger states explicitly that quadrature "can become prohibitively difficult because of the highly oscillatory nature of the integrands" | Fourier transform of a B function is exceptionally simple; STOs are finite linear combinations of B functions | Closed forms (replaced by accelerated quadrature); pi-free | **BORDERLINE** | Weniger, *The strange history of B functions or how theoretical chemists and mathematicians do (not) interact*, arXiv:0811.3406, Int. J. Quantum Chem. (2009), DOI 10.1002/qua.22014 — **VERIFIED** (full text) |
| **Safouhi / Berlu / Slevinsky nonlinear transformations** (SD-bar, D-bar, S-transform + double-exponential + Sinc; rational minimax) | High *pre-determined* accuracy on 3-centre nuclear attraction and 4-centre Coulomb integrals. No correlated molecule. | Same semi-infinite Bessel integral; S-transform to a sine integral, DE transformation to a bi-infinite fast-decaying integrand, Sinc quadrature; or rational minimax approximants in double precision | Convergence acceleration / nonlinear sequence transformations. Weniger's own verdict: "This may well be the currently most promising approach for the evaluation of complicated molecular multicenter integrals of exponentially decaying functions." | Closed forms; pi-free (quadrature nodes); exactness | **BORDERLINE** | Slevinsky & Safouhi, *Compact Formulae for Three-Center Nuclear Attraction Integrals Over Exponential Type Functions*, arXiv:1905.13537; J. Math. Chem. (2022) DOI 10.1007/s10910-022-01362-7 — **VERIFIED** (abstract). Safouhi, *Efficient and rapid numerical evaluation of the two-electron, four-center Coulomb integrals...*, J. Comput. Phys. **176**, 1–19 (2002) — bib **VERIFIED** via Weniger ref [53]; content **UNVERIFIED** |
| **SMILES** (Fernandez Rico, Lopez, Ramirez, Ema, Zorrilla, Ishida), shipped inside **Molpro** | 1- and 2-centre: 12 decimals, few microseconds/integral. 3- and 4-centre: **8 decimals, hundreds of microseconds**. Callable from HF, CCSD(T), MULTI, MRCI. | 3-centre (AB|AC) and 4-centre (AB|CD) **by STO-nG Gaussian expansion** (discrete Gauss transform), NGSSTO default 9, max 30 | Gaussian expansion used as an *evaluator* — exactly GeoVac's noci_engine pattern | Exactness. Molpro's own manual warns: integrals from STO-9G expansions "may have not sufficient accuracy for post-HF calculations, specially with high quality basis sets" | **STOP as stated** (the corpus already runs this evaluator and knows its bias); BORDERLINE if the bias is removed | Molpro manual, section SMILES — **VERIFIED** (documentation). Fernandez Rico et al., *Efficiency of the algorithms for the calculation of Slater molecular integrals in polyatomic molecules*, J. Comput. Chem. **25**, 1987–1994 (2004) — bib **VERIFIED** via Molpro docs; accuracy figures **UNVERIFIED at source** (Wiley 403, PubMed cookie wall; figures from two independent search summaries) |
| **ADF: PARI-MP2 and SOS double hybrids with STOs** (Förster, Franchini, van Lenthe, Visscher) | **MAD < 1 kcal/mol** vs DF-MP2/CBS on S66 dimerisation and on 152 GMTKN30 conformer points; **0.30 kcal/mol** vs CCSD(T)/CBS on HEAVY28. Production code. | **Four-centre integrals are never formed.** Pair-atomic resolution of the identity: AO pair products expanded in **Slater-type auxiliary fit functions centred only on the two atoms of the pair**; plus distance screening (DCAB, DCAC multipole, NHF) | **Density fitting = locality/distance sparsity** | Closed forms; pi-free; *and it requires a sparsity kind GeoVac does not have* | **GO at MP2 level** (correlated, molecular, sub-kcal/mol) | Förster, Franchini, van Lenthe, Visscher, *A Quadratic Pair Atomic Resolution of the Identity Based SOS-AO-MP2 Algorithm Using Slater Type Orbitals*, J. Chem. Theory Comput. **16**(2), 875–891 (2020), DOI 10.1021/acs.jctc.9b00854 — **VERIFIED** (body) |
| **Zero-variance Monte Carlo ERIs + near-FCI with STOs** (Caffarel) | **Be**: exFCI −14.6180020(1) vs Ema et al. FCI −14.618002 (agreement to 1e−7). **CH4** (VB1 STO): E_HF −40.21485042(7), exFCI −40.413651(3), i.e. **3 microhartree statistical error on a 5-centre molecule**; cf. cc-pVDZ GTO exFCI −40.392975. **Cyanine [H2N(CH)NH2]+**, 90 STOs, 8.3e6 two-electron integrals: exFCI statistical error 2e−4 Ha ≈ 0.1 kcal/mol, stated as "sub-chemical accuracy" | **No analytic integration at all.** Multi-centre ERIs over *arbitrary* orbitals by zero-variance MC: a Gaussian expansion supplies the deterministic control variate, correlated sampling computes the exact-minus-approximate difference, so the Gaussian bias is *removed*, not inherited. Absolute integral error ~2e−9 achieved | Zero-variance MC estimator + universal Gaussian sampling | Closed forms; determinism; pi-free. **Keeps symmetry sparsity** — integrals that vanish by selection rule are simply never sampled | **GO** (the only correlated ETO molecule found at true chemical accuracy with genuinely multi-centre integrals) | Caffarel, *Evaluating two-electron-repulsion integrals over arbitrary orbitals using Zero Variance Monte Carlo: Application to Full Configuration Interaction calculations with Slater-type orbitals*, arXiv:1906.04515; J. Chem. Phys. (2019), DOI 10.1063/1.5114703 — **VERIFIED** (full text). Volume/page **UNVERIFIED** |
| **Coulomb resolution (Gill) applied to ETOs; STOP package** (Hoggan, Berlu) | H2 SCF total energy **−1.1284436 Ha** vs an HF-limit estimate −1.1336296 Ha (l = 0 only, 5.2 mHa short). (H2)2: total −2.256998 Ha, well depth 0.069 kcal/mol at 6 au; 4-centre integrals in **milliseconds** (STOP 12 ms, POISON 10 ms, OVERLAP 2 ms). NMR shielding tensors at CPHF level are the actual application. | **3- and 4-centre two-electron integrals are reduced to sums of products of one-electron overlap-like integrals involving orbitals on at most TWO centres**, each separable in prolate spheroidal coordinates. "may require 10 or even 20 terms to converge to chemical accuracy"; Schwarz screening for distant pairs | **Resolution of the identity on the Coulomb operator itself** (1/r12 = sum |phi_i><phi_i| with <f_i|1/r12|f_j> = delta_ij). Claimed SCF cost n^2 rather than n^4 | Exactness of the two-electron integral (the resolution is truncated); **but it PRESERVES the two-centre closed forms and prolate separability that GeoVac already owns** | **BORDERLINE** on demonstrated accuracy; **most interesting structurally** | Hoggan, *How specific exponential type orbitals recently became a viable basis set choice in NMR shielding tensor calculation*, arXiv:1010.5425 (2010) — **VERIFIED** (full text). Companion: Hoggan, *Four-center Slater-type orbital molecular integrals without orbital translations*, Int. J. Quantum Chem. (2010), DOI 10.1002/qua.22213 — **UNVERIFIED** (vol/pages not reached) |

### Caution on the Hoggan H2-dimer numbers
The H2 calculation is SCF; the quoted 0.057–0.069 kcal/mol "Van der Waals well"
at 6–6.4 au cannot be a dispersion well at Hartree-Fock level, where dispersion
is absent by construction. Either a correlated step is implied but not stated in
the text read, or the feature is a basis-set artefact. **Do not cite this as a
correlated ETO result.** The only defensible accuracy datum in that paper is the
SCF total energy, 5.2 mHa above the HF limit with an l = 0 basis.

---

## 2. Paper 59 prior-art test

**Claim under test** (papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex,
L60): "This places a molecular integral inside the modern theory of Bessel
moments and two-mass elliptic Feynman integrals; a targeted literature search
finds that connection unmade in either field."

**Result: none found.** No statement of an elliptic, hypergeometric, or any other
*transcendence class* for three-centre exponential-type integrals appears in the
axis-C literature scanned. What the literature says instead is weaker and
categorically different — it is a statement of *absence of closed form* and of
*numerical* difficulty, never of analytic class:

- Hoggan (arXiv:1010.5425), verbatim: "the infinite series arising when
  Hartree-Fock two-electron integrals **that do not possess closed forms (three
  and four center terms)** are evaluated converge much more slowly when the
  negative alpha functions are used."
- Hoggan, on why Gaussians won, verbatim: "The essential advantage they had over
  exponential basis sets was **the simple product theorem for gaussians** on two
  different atomic centers. This allows all the two-electron integrals, including
  three- and four-center terms to be expressed as single-center two-electron
  integrals." And immediately after: "The corresponding relationship for
  exponential type orbitals generally led to **infinite sums**..."
- Weniger (arXiv:0811.3406), verbatim: "The key problem of this approach is that
  evaluation by numerical quadrature can become prohibitively difficult because
  of the **highly oscillatory nature of the integrands**." His forward-looking
  sentence names convergence acceleration, not analysis: Safouhi's approach "may
  well be the currently most promising approach for the evaluation of complicated
  molecular multicenter integrals of exponentially decaying functions."
- Safouhi/Slevinsky's own framing is numerical throughout: the enabling device in
  the 2019/2022 three-centre paper is "rational minimax approximants that
  minimize the maximum error on the interval of evaluation" — an approximation-
  theory object, not a period.
- The two-centre side is where "elliptic" *does* appear in this literature, and
  only as **elliptical (prolate spheroidal) coordinates**, not elliptic
  integrals: e.g. the standard remark that with unequal exponents one has "two
  elliptical exponents alpha and beta", so the Neumann expansion becomes an
  infinite series. This is a coordinate-system word collision and must not be
  mistaken for prior art. (Source: search-level summary of the three-centre
  Coulomb-exchange literature; the ScienceDirect original returned 403 —
  **UNVERIFIED**.)
- On the mathematics side the Bessel-moment/elliptic literature is well
  developed and entirely disjoint from chemistry: Bailey, Borwein, Broadhurst &
  Glasser, *Elliptic integral evaluations of Bessel moments*, J. Phys. A **41**,
  205203 (2008) (arXiv:0801.0891) — **UNVERIFIED at source here**, and already
  covered on the Feynman side by debug/lit_scan/elliptic_eri_feynman_memo.md.
  Its applications are named as QFT, lattice Green functions (hexagonal, diamond,
  cubic) and condensed matter. No molecular-integral application appears.

**Verdict:** the Paper 59 sentence survives this axis. Strengthened, in fact,
with a mechanism for *why* the bridge is unmade: the ETO community never posed
the analytic-class question at all. From Steinborn's 1970s programme onward the
remaining semi-infinite Bessel integral was treated as a **quadrature problem**
to be beaten with sequence transformations, so the object whose period Paper 59
identifies was never asked what it was. Paper 59's existing hedge ("an absence
from a targeted search, not a proof") remains the right tier and should not be
strengthened on this evidence.

Searched (axis C): Weniger + Steinborn B functions Fourier multicentre; Safouhi
SD-bar three-centre nuclear attraction; four-centre B-function Coulomb
extrapolation; "elliptic integral" + three-centre + Slater/exponential;
molecular integrals + complete elliptic integrals + three-centre two-electron;
Weniger asymptotic analysis + semi-infinite Bessel + "no closed form"; Bessel
moment + elliptic curve + quantum chemistry + three-centre. Plus full-text
keyword sweeps for ellipt / hypergeom / closed form / transcend / Appell /
Lauricella / Meijer over the Weniger and Hoggan full texts (1 hit for "ellipt",
a reference to Agmon on elliptic *equations*; 2 for "hypergeom", both the
terminating 1F1 form of the reduced Bessel function; 0 for transcend / Appell /
Lauricella / Meijer).

---

## 3. Contradictions with corpus statements (5 bullets)

1. **No contradiction — direct corroboration of the product-theorem reason.**
   The corpus's stated reason Gaussians won ("Slater functions have no product
   theorem") is Hoggan's verbatim reason too. This should be *cited*, not
   re-derived: the corpus currently asserts it without external backing.
2. **The corpus's Gaussian evaluator has a named production precedent AND a
   named production caution.** SMILES/Molpro does exactly what noci_engine does
   (STO-nG expansion of 3-/4-centre ERIs), and Molpro's own manual warns it "may
   have not sufficient accuracy for post-HF calculations." This is external
   support for the corpus's "exact != accurate" orientation, but it is also a
   warning aimed squarely at the corpus's own evaluator: at STO-9G the bias is a
   *post-HF-relevant* bias, not a rounding detail.
3. **"GeoVac has symmetry sparsity only, no locality sparsity" is confirmed as
   the binding difference, not merely a difference.** Every ETO path that
   actually reaches production accuracy on molecules buys it with a locality
   device — PARI + DCAB/DCAC/NHF screening in ADF, Schwarz screening in the
   Coulomb-resolution work. None buys it with closed forms. The corpus is
   correct that density fitting would be a *new sparsity kind*; the literature
   adds that it is, empirically, the *only* kind that has worked at scale for
   exponential bases.
4. **One corpus framing is too strong in one direction.** The Poly-2 record
   treats the three-centre two-body block as the wall. Caffarel's result shows
   the wall is specifically a *closed-form* wall, not an *accuracy* wall: with a
   stochastic evaluator the same integrals reach 2e−9 absolute and support a
   near-FCI CH4 at 3 microhartree statistical error. The corpus should say
   "no closed form, and no deterministic cheap route" rather than anything that
   reads as "not computable to chemical accuracy".
5. **Nothing found contradicts the three-bases map, the two-centre closed-form
   inventory, or the l-sparsity-lost-at-two-centres theorem.** The ETO
   literature independently confirms two-centre separability in prolate
   spheroidal coordinates as the thing that works and the thing that stops at
   two centres.

---

## 4. Single most promising untraveled path

**Gill's Coulomb resolution, applied to the GeoVac basis.** It is the only device
found that attacks the three-centre two-body block *without* surrendering what
GeoVac owns. The move is a resolution of the identity on the operator rather than
on the density: write 1/r12 = sum_i |phi_i><phi_i| in a set chosen so that the
Coulomb operator is the identity matrix in it, and every three- and four-centre
two-electron integral becomes a finite sum of products of one-electron
overlap-like integrals, each involving orbitals **on at most two centres** —
which is precisely the class for which the corpus already has weight-one,
pi-free closed forms in {exp, E1, log, gamma}, and precisely the class that is
separable in prolate spheroidal coordinates. It converts the missing product
theorem into a truncation, and Hoggan reports 10–20 terms to chemical accuracy
with Schwarz screening for distant pairs. It also preserves symmetry sparsity
(the one-electron factors carry the same Gaunt selection rules) rather than
destroying it the way a Löwdin retrofit does.

The honest obstruction, which the corpus is already positioned to evaluate and
which should be the first gate: a Coulomb resolution needs an **auxiliary set
complete in the relevant metric**, and v5.11.2 measured exactly that for the
one-centre set and found it wanting — overcompleteness is one direction only,
with the Bessel deficit plateauing at 0.38–0.91, so the one-centre set is *not*
complete in the molecular metric. The resolution's truncation error is therefore
not obviously controllable with GeoVac's existing one-centre functions. The
cheap decisive experiment is to take the measured deficit and ask what
resolution error it implies for a single (XY|XZ) integral against the corpus's
own momentum-space reference value — a scalar test, no new machinery. If the
error is above ~1e−6 the path closes on the corpus's own prior measurement; if
it is below, the three-centre two-body block opens through two-centre closed
forms the corpus already has.

Runner-up, for a different purpose: Caffarel's **zero-variance correlated
sampling** is the correct way to remove the bias from the Gaussian evaluator the
corpus already runs, at the price of determinism and 4800 cores for hours. It is
not a structural advance and it forfeits pi-freeness, but it is the demonstrated
route to a chemically-accurate correlated molecule in a genuine exponential
basis, and it leaves the qubit-Hamiltonian sparsity *pattern* untouched because
symmetry-zero integrals are never sampled.

---

## 5. Search trail

Web calls: 25 (at cap). Searches (10): Weniger/Steinborn B functions Fourier
multicentre STO; Safouhi SD-bar three-centre NAI; four-centre B-function Coulomb
Safouhi extrapolation; STO correlated MP2/CCSD "chemical accuracy" multicentre;
ADF STO MP2 double hybrid PARI benchmark kcal/mol; "elliptic integral"
three-centre molecular integrals STO closed form; Rico/Lopez/Ramirez/Ema SMILES
four-centre ERI; Weniger asymptotic semi-infinite Bessel "no closed form";
molecular integrals "complete elliptic integral" three-centre two-electron;
Bouferguene STOP / Talman numerical Fourier STO; Fernandez Rico + FCI/CI/MP2
microhartree; Bessel moment + elliptic curve + quantum chemistry.

Fetches (15), of which succeeded: arXiv abs 0811.3406 (abstract + bib);
arXiv PDF 0811.3406 (binary -> extracted locally with pypdf, 13 pp);
arXiv PDF 1010.5425 (binary -> extracted locally, 23 pp); arXiv abs 1905.13537;
arXiv abs 1906.04515; arXiv PDF 1906.04515 (binary -> extracted locally, 11 pp);
PMC7027358 (ADF PARI-MP2, body + full bib, two fetches); Molpro SMILES manual.

Dead ends: ScienceDirect S0009261405010109 (403); Wiley jcc.20131 (403);
PubMed 15473010 (cookie wall); Semantic Scholar search page (empty render).
Note for future scans: WebFetch on an arXiv PDF returns an unreadable binary to
the summarizing model but **saves the PDF locally**; extracting it with pypdf
and grepping the text is far more reliable than the fetch summary, and costs no
web call. This recovered three of the four load-bearing full texts in this memo.

Working extracts (transient, in the session tool-results dir, not in the repo):
axisC_extract.txt (Weniger + Hoggan), caffarel.txt.
