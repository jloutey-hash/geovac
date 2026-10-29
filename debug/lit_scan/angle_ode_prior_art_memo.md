# Prior art: an ODE in the bond angle for (XY|XZ), and how many terms guessing needs (2026-10-28)

Scope: literature search before the AO (angle-operator) sprint. Already credited in Paper 59 and not re-reported here: Pachucki PRA 2009/2010, Lesiuk-Moszynski arXiv:1209.0839, Puchalski-Kedziera-Pachucki PRA 80 032521, Laporta 2000, Chetyrkin-Tkachov 1981.
How each source was read: **FT** = full text read (PDF, pdftotext); **ABS** = abstract or publisher page only; **SEC** = known only from a search snippet or a citing source. Downloads are in `debug/lit_scan/dl_angle_*/`, excluded from git.

## Verdicts

| # | Question | Verdict |
|---|---|---|
| 1 | DE in a geometric parameter (angle, or distance to a third centre) for 3- or 4-centre Slater/Yukawa integrals | **NOT FOUND** |
| 2 | Position-space IBP with exponential/Yukawa lines and fixed centres | **PARTIAL** (two-centre only, and done through momentum space) |
| 3a | BKKS: thousands of Mellin moments | **FOUND** (FT) |
| 3b | Rule for how many terms guessing needs; trading order for degree | **FOUND** (FT) |
| 3c | Approximants hide the true equation at finite data | **PARTIAL** (FT; the exact-arithmetic version only) |
| 4 | Closed form or ODE in the bond angle for 3-centre nuclear attraction | **NOT FOUND** |

## 1. Geometric-parameter DEs for multicentre integrals: NOT FOUND

- Pachucki, arXiv:0905.1635 (FT). The DE is in the internuclear distance and is obtained by mapping to momentum space. Abstract: "the master integral with the single negative power of all distances can be obtained from the second order differential equation in r, the distance between nuclei." Sec. II: "We use the method of integration by parts identities [2] ... in the momentum space representation of the integral G". Here G is the one-centre tetrahedron integral with six massive lines (Eq. 8). After a Laplace transform, the internuclear distance appears as one more mass. The text also says: "The solution of this differential equation is presented in the work of Fromm and Hill ... However, they have not used the differential equation approach". None of arXiv:0905.1635, 1209.1258 or 1307.6065 (all FT) uses the strings "three-cent", "polyatom" or "triatom" (grep).
- Our own reading, not a statement in any source: the trick works because a single distance Laplace-transforms into a mass. A bond angle is not a distance, so the same trick does not supply an angle ODE.
- Özdoğan and Ruiz, arXiv:1209.3755 (SEC). Three-centre two-electron STO integrals, analytic, "first for linear conformation"; the triangular case was deferred. Not checked for any recursion in the angle.
- Searches for Yukawa three-centre integrals, H3+ explicitly correlated exponentials, and "three-center" + "integration by parts" found nothing relevant (Gaussian recursions run over Cartesian indices, not the angle).

## 2. Position-space IBP with fixed centres: PARTIAL

- The nearest prior art is the Pachucki line above. Its IBP is done in momentum space, not position space, and it covers two centres only.
- Chetyrkin, Kataev, Tkachov, *Nucl. Phys. B* **174**, 345 (1980), the Gegenbauer x-space technique (ABS/SEC). It works in coordinate space for **massless** diagrams. Per a citing source, "only the lines, attached to the so-called root vertex (x = 0) in coordinate space, may have arbitrary powers."
- Schröder and Meyer, arXiv:2511.17349 (ABS). Coordinate-space perturbation theory with one mass scale; Bessel-function propagators; three-point correlators via a Gegenbauer expansion. The abstract mentions neither IBP nor DEs.
- Loebbert, Müller, Münkler, *PRD* **101**, 066006 (2020), arXiv:1912.05561, and Loebbert, Miczajka, Müller, Münkler, *SciPost Phys.* **11**, 010 (2021) (ABS). Yangian differential constraints for conformal and massive one-loop integrals in dual-coordinate variables, i.e. DEs in cross-ratios with fixed external points. The propagators are 1/(x²+m²), not e^{-mr}/r, and only one loop is treated. This is an analogue only.
- Structural observation (ours, check before using): (XY|XZ) is a 3-D position-space two-loop graph. The integrated vertices are r1 and r2, the fixed external points are X, Y, Z, the propagators are e^{-mr}/r (with m-derivatives for 1s), and 1/r12 is massless. In this language a geometric ODE is a DE in the external-point invariants. **No paper found** reduces such a graph with fixed external points by IBP.

## 3. Guessing from many coefficients

**3a. FOUND (FT).** Blümlein, Kauers, Klein, Schneider, *Comput. Phys. Commun.* **180**, 2143 (2009), arXiv:0902.4091.
- Abstract: "For the most complicated subproblem 5114 moments were needed in order to produce a recurrence of order 35 whose coefficients have degrees up to 938."
- Why so many moments are needed (Sec. 2): "the number of unknowns in the linear system determines the number of initial terms qi that have to be known a priori in order to detect the recurrence. According to our experience, the size of the linear system is minimized when the order l and the degree d are approximately balanced."
- The minimal-order recurrence has the smallest coefficients but "the maximal number of unknowns in the corresponding linear system, and vice versa". The system for C^(3)_{2,q,CF^3} "would require approximately 11Tb of memory", so the authors work modulo primes. The moments are **exact rationals**: moment N = 5114 has "13888 numerator and 13881 denominator digits".
- Sec. 5 notes that direct computation could reach only N ≈ 16, so the method needs far more moments than exist ab initio.
- Related (SEC only): Behring et al., arXiv:1908.03779, polarized γ_qg^(2), "4000 moments were generated ... 2640 moments turned out to be sufficient".

**3b. FOUND (FT).** Kauers and Koutschan, "Guessing with Little Data", ISSAC 2022, pp. 83-90, arXiv:2202.07966.
- With (r+1)(d+1) unknowns, "If r and d are chosen such that (r+1)(d+2) ≤ N+2, this homogeneous linear system is not expected to have a (nontrivial) solution. If it does have a solution, this is interpreted as evidence ... The more the number of equations exceeds the number of variables, the stronger is the evidence that the solution is not just noise".
- On trading order for degree: "an equation of slightly higher order may have substantially lower degree, so that we may find an equation using r+1 and d/2, for example. This phenomenon is well understood [14]". Ref. [14] is Jaroschek, Kauers, Chen, Singer, "Desingularization explains order-degree curves for Ore operators", ISSAC 2013. We saw it only as a reference entry.
- Tests that separate a true guess from a wrong one (from Bostan-Kauers, FPSAC 2009): short integers; integrality when unrolled; "nice singularities"; nilpotent p-curvature.

Bostan, Boukraa, Guttmann, Hassani, Jensen, Maillard, Zenine, *J. Phys. A* **42**, 275209 (2009), arXiv:0904.1601 (FT). This is the clearest statement of the cost and the fix:
- "The minimal order linear ODE may contain a very large number of apparent singularities and can thus only be determined from a very large number of series coefficients (generally speaking (q+1)(D0+1) terms are needed). Other (non-minimal order) linear ODEs, because they carry polynomials of smaller degrees, may require fewer series coefficients".
- "In the case of χ̃(5) this optimal linear ODE can be discovered from some 7400 terms while the minimal order linear ODE requires almost 49100 terms."
- The minimal order (29) is then *inferred* from several non-minimal ODEs by their "ODE formula", which gives the number of required terms as linear in Q and D (App. B). Spurious exponents of a non-minimal ODE are filtered by the rationality of the exponents.

**3c. PARTIAL (FT).** Guttmann, arXiv:1604.08254 (differential approximants): "such ODEs implicitly define all coefficients in the generating function, but if N terms are used in the construction of the ODE, all terms of order z^{N+1} and beyond will be approximate, unless the exact ODE is discovered". All of the sources above work in exact or modular arithmetic. **We found no source** that analyses guessing from *floating-point* coefficients, where at finite precision near-null vectors of an under-determined system look like true equations. The transferable rules are: require the system to be well over-determined; prefer non-minimal operators to cut the number of terms needed; check every guess with the a-posteriori tests above.

## 4. Bond-angle ODE for 3-centre nuclear attraction: NOT FOUND

- Bouferguene and Jones, *JCP* **109**, 5718 (1998) (ABS): a convergence analysis of the STO addition theorem that picks a series accelerator. No angle ODE appears in the abstract.
- Barker and Eyring, *JCP* **21**, 912 (1953) (ABS): s-orbital three-centre attraction expressed through auxiliary functions. Its recurrences run over integer indices, not over the angle.
- Slevinsky and Safouhi, arXiv:1905.13537 (ABS): compact analytic formulae with minimax approximants. The abstract mentions no recursion in the geometry.
- Bağcı and Hoggan, arXiv:1503.03690 (SEC): three-centre nuclear attraction integrals "generally lack analytically closed forms" (snippet paraphrase).
- Jain, Jones, Weatherford, Hoggan, *IJQC* (2004), doi 10.1002/qua.20089 (SEC): truncated series, not a closed form.
- None of these sources gives an ODE in cos γ.

## Bottom line

The angle-ODE target looks unclaimed in the sources we could check. The guessing budget is well documented: unknowns = (order+1)(degree+1); the minimal operator is the most expensive to guess; non-minimal operators need several times fewer terms (7400 vs 49100 for χ̃(5)). An ~140-coefficient float fit failing is consistent with this literature. The 1000+ estimate is plausible for a minimal operator, but a non-minimal ansatz (higher order, lower degree) can cut it.
