# Literature scan — Axis D: two-centre / spheroidal Sturmians and multi-scale exponential bases, ACCURACY side

Date: 2026-09-13. Read-only scan. No papers, code or tests edited.
Caps observed: 23 web calls (limit 25), 12 sources (limit 12).

Assigned question: does a two-focus or multi-scale exponential basis (a) keep a
discrete angular label with algebraic couplings, (b) avoid the overcompleteness of
two one-centre sets, and (c) reach chemical accuracy on a diatomic or larger — and
what does a third centre cost?

---

## 1. Headline table

| # | Path | Best demonstrated accuracy + system + reference | Angular label discrete? couplings algebraic? | Overcompleteness avoided how? | Third-centre cost | Verdict | Source |
|:--|:-----|:-----|:-----|:-----|:-----|:-----|:-----|
| D1 | Prolate-spheroidal FEM/DVR, partial waves Y_lm in the eta variable, Neumann expansion for 1/r12 | **H2 two-electron ground state -1.17442 Ha vs exact -1.17447 (Wolniewicz) = 0.05 mHa, at l_max = 6.** Same paper: a single-centre spherical expansion at l_max = 7 gives -1.16908 (5.4 mHa, "0.15 eV higher"). | YES (l, m; M = m1+m2 conserved). Angular overlaps are CLOSED FORM in (l, m) — their Eq. (15) is a three-term algebraic expression for the (xi^2 - eta^2) volume-element matrix element; the Neumann expansion of 1/r12 gives Gaunt-type angular products exactly. RADIAL part is a numerical Poisson solve on the DVR grid (sparse, diagonal in radial DVR indices). | Structurally: ONE two-focus basis, and the DVR is orthogonal by construction. Not a remedy — the problem does not arise. | Dies. Prolate spheroidal has exactly two foci; the Neumann expansion is a two-centre object. Paper says it "opens the way to calculations on more complicated **diatomic** targets". | **BORDERLINE-high** (accuracy reached; angular closure algebraic, radial numerical; 1e and 2e only) | Tao, McCurdy, Rescigno, "Grid-based methods for diatomic quantum scattering problems. III. Double photoionization of molecular hydrogen in prolate spheroidal coordinates", Phys. Rev. A **82**, 023423 (2010), DOI 10.1103/PhysRevA.82.023423 — **VERIFIED** (accepted-manuscript full text read at OSTI 1051651; journal/vol/page from the PRA landing page, article body read) |
| D2 | HelFEM: partial-wave prolate-spheroidal finite elements, chi_nlm = B_n(mu) Y_l^m(nu, phi) | **Sub-microhartree ROHF-limit energies for 70 diatomics, periods 1–4**, reproducing / correcting published numerical-HF limits with orders of magnitude fewer parameters than x2dhf. l_max needed is LARGE: Cl2- needs l_sigma = 38, l_pi = 29. | YES. "in the partial wave approach angular integrals can be evaluated analytically in closed form"; **Gaunt coefficients are precomputed and stored**; the Neumann-expansion angular integrals "can be performed analytically within the partial wave expansion". Radial = FEM (15th-order Lobatto), numerical. | Explicit: "**because the finite element basis set is never ill-conditioned**, symmetric orthonormalization is used…". One two-focus set; no linear dependence. | Dies. HF/DFT on **diatomics only**; the partial-wave method is limited to two centres. | **BORDERLINE** (accuracy is the HF/DFT limit, not a correlated energy; no post-HF has been published on this basis) | Lehtola, "Fully numerical Hartree–Fock and density functional calculations. II. Diatomic molecules", Int. J. Quantum Chem. **119**, e25944 (2019), DOI 10.1002/qua.25944, arXiv:1810.11653 — **VERIFIED** (full arXiv text read) |
| D3 | Generalized Sturmian Functions (GSF) in prolate spheroidal coordinates | **H2+ 1sigma_g: E = -1.1026346 vs benchmark -1.1026342 (7 digits) with 4 angular + 6 radial functions**; excited states 1sigma_u, 2sigma_g, 3sigma_u, 3d to 6–7 digits; heteronuclear HHe2+ -2.2506056, HLi3+ -4.7501126; R_eq found at 1.99704 a0, E_tot -0.602635. | Discrete: the eta equation is a separation-constant (A) eigenproblem — a discrete quasi-angular label. "**The matrices of the generalised eigenvalue problem … are easily constructed as they are all analytical**" for the ANGULAR equation. The RADIAL side is "completely numerical" (predictor-corrector ODE generation of the basis + 10^4-point quadrature for matrix elements). | ONE two-focus set (basis functions are functions of xi and eta, not of two centres). | Dies at two foci. Paper's own roadmap is H2 and N2 as **two-electron / quasi-two-electron diatomic** targets — no polyatomic claim. | **BORDERLINE** (one-electron only; smallest basis in the field) | Mitnik, López, Ancarani, "Generalized Sturmian Functions in prolate spheroidal coordinates", Mol. Phys. **119**, e1881179 (2021), DOI 10.1080/00268976.2021.1881179, arXiv:2006.06616 — **VERIFIED** (full text read) |
| D4 | Coulomb Sturmians **derived in** spheroidal coordinates (separation of variables + direct solution of the 1D equations) | **H2+ -1.102614 vs -1.1026342 = 0.02 mHa; HLi3+ -4.750111 vs -4.7501126**, with 10 basis functions per nucleus. | Discrete (spheroidal separation constant). Couplings: NOT established from the sources reached — the abstract claims "completeness and good convergence"; no algebraic-closure statement reached. UNVERIFIED. | NOT avoided in the form used: "10 basis functions **per nucleus**" (per Mitnik's comparison), i.e. still two sets. The abstract's selling point is the *similarity* of one- and two-centre orbitals in spheroidal coordinates, not the elimination of a second set. | Dies at two foci. | **BORDERLINE / STOP** (one-electron; keeps the per-nucleus duplication the corpus already priced) | Kereselidze, Chkadua, Defrance, "Coulomb Sturmians in spheroidal coordinates and their application for diatomic molecular calculations", Mol. Phys. **113**(22), 3471–3479 (2015), DOI 10.1080/00268976.2015.1036146 — **VERIFIED** (abstract at KCL Pure; the two energies read in Mitnik's Table 5 at source). Companion chapter Kereselidze & Ogilvie, Adv. Quantum Chem. (2018), DOI 10.1016/S0065-3276(18)30003-0 — **UNVERIFIED** (Semantic Scholar page returned empty). |
| D5 | Avery many-centre Coulomb Sturmians with molecular V0 = sum_A Z_A/r_A; Goscinskian isoenergetic configurations; Shibuya–Wulfman integrals | **H2+ -1.10220 (0.43 mHa error) and HLi3+ -4.75011 with 10 Coulomb Sturmians per nucleus** (read in Mitnik Table 5). For **N-electron molecules the 2009 paper proposes a method** — "A method is proposed for using isoenergetic configurations formed from many-center Coulomb Sturmians…"; no N-electron molecular energies reached in the sources scanned. | Discrete (n, l, m) per centre; SW integrals are algebraic via Fock hyperspherical harmonics (this is the corpus's own Paper 60 machinery). | **NOT avoided** — this is exactly the two-one-centre-sets configuration the corpus measured (cond ~ (1/error)^1.8, one near-dependent direction). | Generalizes to many centres in principle (SW integrals are many-centre by construction) — the failure is not the third centre, it is that the shared-k requirement is baked in. | **STOP for accuracy** (no demonstrated N-electron molecular accuracy; the one-electron numbers are *worse* than every two-focus route above) | Avery & Avery, "Can Coulomb Sturmians Be Used as a Basis for N-Electron Molecular Calculations?", J. Phys. Chem. A **113**(52), 14565–14572 (2009), DOI 10.1021/jp9040502 — **VERIFIED at abstract level** (abstract via PubMed 19807119 listing + publisher listing; the numbers via Mitnik Table 5 read at source). Avery, "Many-center Coulomb Sturmians and Shibuya–Wulfman integrals", Int. J. Quantum Chem. (2004), DOI 10.1002/qua.10820 — **VERIFIED at abstract level**. |
| D6 | **Mixed-scale Sturmians: per-centre k_A != k_B** | **NONE FOUND.** Every implementation reached states a single shared exponent: "the CS exponent k is identical for each basis function"; the molecular Sturmian "automatic scaling" adjusts one k with R, it does not give per-centre k. | n/a | n/a — this is the proposed *remedy*, untested | n/a | **UNTRAVELED** (no evidence for or against; not a ledger entry, not walled) | Herbst, Dreuw, Avery et al. lineage: "Quantum chemistry with Coulomb Sturmians: Construction and convergence … at Hartree–Fock level", Phys. Rev. A **99**, 012512 (2019), DOI 10.1103/PhysRevA.99.012512, arXiv:1811.05777 — **VERIFIED at abstract/summary level only** (shared-k statement from search summary of the paper, body not read). Treat the shared-k claim as **UNVERIFIED-at-body**. |
| D7 | Novosadov's Sturmian alternative to LCAO (as transmitted by Duchon–Dumont-Lepage–Gazeau) | **One-electron many-centre only.** The paper gives two Sturmian techniques for the energy eigenvalues of a many-centre **one-electron** system; the first is Novosadov's (also revisited by Monkhorst–Jeziorski), the second is closer to LCAO and shows "excellent convergence" on **H2+**. No numerical values reached; no molecular N-electron result. | Not reached. | Not reached. | Formally many-centre (that is the technique's point) — but never carried past one electron in the record reached. | **STOP** (one-electron; and the first technique is the Monkhorst–Jeziorski thread the corpus has already scanned) | Duchon, Dumont-Lepage, Gazeau, "On two Sturmian alternatives to the LCAO method for a many-center one-electron system", J. Chem. Phys. **76**(1), 445–447 (1982) — **VERIFIED at listing level** (AIP article-abstract listing; full abstract page returned HTTP 403). **DOI UNVERIFIED — not reached, do not invent one.** Novosadov's own Russian-language originals: **UNVERIFIED, not reached.** |
| D8 | Ponomarev–Somov / Solov'ev two-centre quasiradial–quasiangular functions | Asymptotic (large-R) and WKB/phase-integral expressions for the two-Coulomb-centre problem; agreement with numerics "when R is greater than the shell size". Not a variational basis with a competitive energy. | Discrete (separation constants), but the content is asymptotic expansion, not a matrix method. | n/a | n/a | **STOP** (no competitive accuracy; one-electron; asymptotic regime) | Thread verified only at review/summary level — Ponomarev & Somov phase-shift work and the WKB two-Coulomb-centre paper, Theor. Math. Phys. (2017), DOI 10.1134/S0040577917030047 — **UNVERIFIED** (search-level only) |
| D9 | B-spline CI in prolate spheroidal coordinates (Vanne–Saenz) | Ground and excited states of **H2**, including doubly-excited autoionizing states, "favourable comparison with literature values in all cases"; no digit count reached. | Prolate spheroidal, molecular symmetry fully accounted for; angular treatment not reached. | One two-focus set. | Dies at two foci. | **BORDERLINE, unpriced** | Vanne & Saenz, J. Phys. B (2004) — **UNVERIFIED** (volume/pages/DOI not reached; description from a search summary only) |

Benchmarks used by the above, for the record: Scott, Aubert-Frécon, Grotendorst, Chem. Phys. **324**, 323 (2006), DOI 10.1016/j.chemphys.2005.10.031 (H2+ reference A and E); Madsen & Peek, At. Data Nucl. Data Tables **2**, 171 (1971); Bian, Phys. Rev. A **90**, 033403 (2014), DOI 10.1103/PhysRevA.90.033403; Wolniewicz (H2 -1.17447 Ha) — all **VERIFIED as cited in the bodies read**, not independently opened.

---

## 2. Answers to the three assigned questions

**(a) Discrete angular label with algebraic couplings — YES, and it is standard practice.**
The two-focus literature does *not* abandon the (l, m) label. HelFEM and the
Tao–McCurdy–Rescigno FEM/DVR both expand the eta (quasi-angular) coordinate in
ordinary spherical harmonics Y_l^m, keep m as an exactly conserved axial label, and
evaluate every angular integral in closed form from **Gaunt coefficients** — Lehtola
precomputes and stores the Gaunt tables, and both note that the **Neumann expansion
of 1/r12 in prolate spheroidal coordinates integrates analytically inside a partial-wave
expansion**. This is, structurally, GeoVac Paper 11 x Paper 12. The novelty in the
literature is entirely on the radial side (finite elements / DVR / numerically generated
Sturmians), which is where all three keep numerical quadrature.

**(b) Overcompleteness — the two-focus route does not have the problem at all.**
It is not remedied, it is absent. Lehtola states it flatly: the finite-element basis in
prolate spheroidal coordinates "is never ill-conditioned". The reason is structural and
matches the corpus's own v5.11.2 finding from the other side: overcompleteness is a
property of *two one-centre sets* spanning overlapping regions, not of exponential bases
per se. Move to one set on two foci and the near-dependent direction has no way to form.
Note the contrast with D4/D5, which stay on per-nucleus sets and therefore keep it.

**(c) Chemical accuracy — reached, twice, but only at 1 and 2 electrons or only at HF.**
- 2-electron correlated: **0.05 mHa on H2 at l_max = 6** (D1). That is 30x inside chemical accuracy.
- N-electron: **sub-microhartree, but only to the ROHF limit**, for 70 diatomics through period 4 (D2). No correlated post-HF calculation on a two-focus partial-wave basis was found.

**Third centre:** every path dies, and for the reason the corpus already owns — prolate
spheroidal has exactly two foci, and the Neumann expansion is a two-centre kernel. Lehtola
is explicit that the partial-wave approach is limited to diatomics. Nothing in this scan
offers a third-centre route; nothing contradicts Paper 59's genus-1 result.

---

## 3. Contradictions / re-pricings against corpus statements

1. **Paper 15's "the per-channel angular basis is the convergence bottleneck" is
   geometry-specific, not an l_max law.** Tao–McCurdy–Rescigno reach 0.05 mHa on H2 at
   **l_max = 6**, and in the *same paper* report that a single-centre spherical expansion
   at **l_max = 7** gives 5.4 mHa. Same molecule, same order of angular truncation, ~110x
   difference. So l_max = 6 is not intrinsically a coarse angular basis; the mol-frame
   hyperspherical choice is what makes it one. This re-prices the bottleneck statement,
   it does not falsify a claim.

2. **The corpus's own two-electron prolate-spheroidal CI is 260x off the published
   ceiling in the same coordinate system.** Paper 12 reports 92.4% of exact D_e
   (E = -1.161304 Ha at N = 72, i.e. **13.2 mHa** from -1.17447), plateauing "regardless of
   basis size". The published FEM/DVR number in the identical coordinate system with the
   identical Neumann kernel is 0.05 mHa. The two visible differences: Paper 12 goes only to
   **l = 3** where McCurdy uses l_max = 6, and Paper 12 uses a **single variationally
   optimized orbital exponent alpha** on a Laguerre set where McCurdy uses a complete radial
   DVR. Paper 12 attributes the plateau to "grid truncation error" — that attribution is
   worth re-examining, since a two-element 14th-order DVR suffices on the published side.
   (Flag, not a correction: I did not re-run anything.)

3. **"Overcompleteness is ONE DIRECTION" (v5.11.2) is consistent with, and gets a
   mechanism from, this scan.** The direction exists because there are two sets. HelFEM's
   never-ill-conditioned statement is the control experiment: one set on two foci, no
   direction. The corpus's band-Toeplitz preconditioner (v5.11.0) treats the symptom; the
   two-focus basis removes the cause.

4. **Avery's many-centre Coulomb Sturmians are, on the one-electron benchmark, the
   *least* accurate of the two-focus family** (0.43 mHa on H2+ with 10 functions/nucleus,
   vs 0.02 mHa for spheroidal Sturmians and 7 digits for prolate GSF with 10 functions
   total). Useful calibration for Paper 60's framing: the shared-k many-centre CS route is
   already near chemical accuracy at one electron, but it is not the accuracy leader even
   there, and its N-electron molecular application remains a 2009 *proposal*.

5. **Per-centre k_A != k_B has no literature record at all.** Every implementation reached
   uses one shared exponent; "generalized Shibuya–Wulfman integrals for mixed scales" did
   not surface as a published, applied object. The corpus's statement that this is "Avery's
   documented route" should be read as *documented as a direction*, not as a demonstrated
   method — I found no worked mixed-scale molecular calculation. (Caveat: the Herbst–Dreuw
   body was not read; shared-k is VERIFIED only at summary level.)

---

## 4. The single most promising untraveled path

**Rebuild GeoVac's two-electron prolate-spheroidal CI (Papers 11 + 12) at the angular and
radial resolution the published two-focus literature actually uses, and take the 13.2 mHa
plateau to the 0.05 mHa ceiling.**

This is the one path in axis D that is (i) demonstrated to reach well inside chemical
accuracy on a correlated diatomic, (ii) built on a discrete (l, m) label with Gaunt-algebraic
angular couplings — so it does not touch the prime directive, (iii) free of overcompleteness
by construction rather than by preconditioning, and (iv) **already half-built in the corpus**:
Paper 11 supplies the algebraic radial solver (sigma states quadrature-free; m != 0 reduced to
the single transcendental seed e^a E1(a)) and Paper 12 supplies the Neumann V_ee. No guardrail
fires: this is not a single-centre encoding (Papers 8–9 do not apply), not graph concatenation
(FCI-M does not apply), and not nested hyperspherical (Track DF does not apply). What the
literature says is missing from the corpus's version is specific and cheap to test: l = 3 ->
l_max = 6, and a single variational exponent -> a radial set that actually spans (McCurdy needs
only two elements of 14th-order DVR; GeoVac's own spectral Laguerre with several exponents, or
a genuine cross-n set, is the algebraic analogue). If the plateau moves, Paper 12's
"grid truncation" attribution is wrong and the corpus gains a correlated diatomic two orders of
magnitude better than its current best. If it does not move, the corpus has isolated a real
structural ceiling of the algebraic-radial choice and can say so with a published control to
compare against.

The honest cap, stated up front: **this is a diatomic-only program.** Two foci is two foci,
and the third centre is the genus-1 wall of Paper 59. Nothing in this scan changes that. The
value is not a route to polyatomics; it is that the corpus would then hold a *chemically
accurate, quadrature-free, discrete-angular, non-overcomplete* two-electron molecular
Hamiltonian — which is a far better anchor for the qubit-encoding story than a 92.4% D_e result,
and a far better platform from which to price what the third centre costs.

Secondary, much weaker: **per-centre k_A != k_B (D6)** is genuinely untraveled and is the only
remedy in scope aimed at the *accuracy* side of the two-one-centre-sets configuration. But it
has no literature to inherit, it does not remove overcompleteness (two sets remain), and the
generalized Shibuya–Wulfman integrals it needs are not a published applied object. Recommend
it only if the PI wants to stay on the many-centre CS architecture for encoding reasons.

---

## 5. Search trail

**Web calls: 23. Sources verified: 12.**

Queries run (in order), and what each yielded:

1. `Sturmian functions prolate spheroidal coordinates two-center H2+ expansion basis` — HIT: the Mitnik/López/Ancarani GSF-prolate family (arXiv 2006.06616, Mol. Phys. 2021, IAFE PDFs).
2. `Avery generalized Sturmians molecular V0 weighting potential Goscinskian configurations different scaling parameters` — HIT: the Avery & Avery 2009 JPCA "Can Coulomb Sturmians…" paper and the Generalized Sturmians book; established that N-electron molecular use is a *proposal*.
3. FETCH arXiv:2006.06616 abstract — abstract only; no accuracy numbers. Prompted the PDF route.
4. FETCH `users.df.uba.ar/dmitnik/publications/sturprolates.pdf` — WebFetch could not parse the PDF, but saved it locally; `pdftotext -layout` then gave the complete article. This was the highest-yield move of the scan (Tables 1–7, the "all analytical" angular statement, the "completely numerical" radial statement, and the comparison Table 5 against Avery and Kereselidze).
5. `Kereselidze Chkadua Defrance "prolate spheroidal" Coulomb Sturmian ... 2015` — HIT: KCL Pure records for both the 2015 Mol. Phys. paper and a 2016/2018 companion.
6. `Avery "Can Coulomb Sturmians Be Used as a Basis..."` — HIT: PubMed 19807119, ACS listing, DOI 10.1021/jp9040502; abstract text via search summary.
7. FETCH KCL Pure record for Kereselidze 2015 — full abstract, volume, pages, DOI. No accuracy numbers in the record (they came from Mitnik's Table 5 instead).
8. FETCH Semantic Scholar page for the Kereselidze–Ogilvie Adv. Quantum Chem. chapter — **DEAD END** (page returned empty content).
9. `Duchon Dumont-Lepage Gazeau JCP 1982 76 445 ... Novosadov` — HIT: exact title and the statement that technique 1 is Novosadov's and was revisited by Monkhorst–Jeziorski; one-electron only.
10. `Lehtola HelFEM ... Kobus x2dhf ...` — HIT: arXiv 1810.11653, the IJQC 2019 paper, the HelFEM repo, and the Lehtola 2019 review.
11. FETCH arXiv:1810.11653 abstract — confirmed basis form chi_nlm = B_n(mu) Y_l^m and the 68-diatomic validation; not enough detail on Gaunt/conditioning.
12. FETCH AIP article-abstract page for JCP 76, 445 — **DEAD END, HTTP 403.** DOI therefore left UNVERIFIED rather than guessed.
13. FETCH arXiv PDF 1810.11653 — WebFetch could not parse; saved locally; `pdftotext` gave the full text. Yielded the Gaunt-precomputation statement, the "never ill-conditioned" statement, the l_sigma = 38 / l_pi = 29 partial-wave counts for Cl2-, and the diatomic-only scope.
14. `Avery Shibuya-Wulfman ... different exponents per nucleus` — partial: confirmed the SW parameter is S = k|X_a0 - X_a| with a **single** k; no mixed-scale applied work.
15. `Vanne Saenz H2 prolate spheroidal B-spline CI` — weak hit; confirmed the method exists for two-electron H2 incl. doubly-excited states, but no volume/pages/DOI reached. Left UNVERIFIED.
16. `Tao McCurdy Rescigno prolate spheroidal FEM-DVR ... PRA 79 012719` — HIT: PRA 79 012719 (2009, paper I), PRA 80 013402 (paper II), PRA 82 023423 (paper III, two-electron), plus a free OSTI PDF of paper III.
17. FETCH core.ac.uk PDF of paper I — **DEAD END, HTTP 403.**
18. FETCH PubMed 19807119 — **DEAD END** (cookie wall; only the notice was returned).
19. FETCH OSTI 1051651 PDF — WebFetch could not parse; saved locally; `pdftotext` gave the full text. **Highest-value hit of the scan:** the l_max = 6 / -1.17442 vs -1.17447 H2 number, the single-centre l_max = 7 / -1.16908 control in the same paper, the Y_lm angular basis, the closed-form Eq. (15) angular overlap, the Neumann expansion, and the Poisson-solve radial treatment.
20. `Ponomarev Somov two-centre ... Solov'ev` — thread is asymptotic/WKB, not a competitive variational basis. Classified STOP; left UNVERIFIED at source.
21. `fully numerical correlated diatomic ... MP2 CCSD full CI ... finite element` — **important negative:** no published post-HF correlated calculation on a two-focus partial-wave basis was found. The literature summary itself notes the partial-wave method "is not limited to the Hartree-Fock level" but that no MP2/CCSD/FCI implementation is evident.
22. `Coulomb Sturmian ... "different values of k" ... per atom` — **important negative:** confirmed the shared-exponent convention ("the CS exponent k is identical for each basis function") and the molecular "automatic scaling" of a single k with R. No mixed-scale molecular calculation found.
23. (local, not a web call) grep of the corpus for `Neumann expansion` and the Paper 11/12 two-electron statements — established the 92.4% D_e / -1.161304 Ha / l = 3 / single-alpha baseline that the recommendation is measured against.

**Dead ends worth recording for the next scanner:** ACS, PubMed and core.ac.uk all refuse WebFetch (403 / cookie wall). The reliable pattern in this scan was: fetch the PDF URL even when WebFetch reports it cannot parse it, then run `pdftotext -layout` on the file WebFetch saves to the tool-results directory. Three of the four load-bearing sources were recovered this way.

**Not reached, deliberately (budget):** Rotenberg's 1970 Adv. At. Mol. Phys. Sturmian review; Judd; Ovchinnikov–Rakitin; Novosadov's Russian-language originals; the Kereselidze 2015 body (only the abstract and Mitnik's quotation of its numbers); Vanne–Saenz bibliographic details; the Herbst–Dreuw PRA 99 body. All tagged UNVERIFIED above.
