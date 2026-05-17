# Bekenstein-Hawking Entropy Sprint — Scoping Memo

**Date:** 2026-05-16
**Status:** SCOPING ONLY — no production code written, no closure attempted.
**Verdict:** GO-WITH-PREREQUISITES, estimated 8–14 weeks for first-pass closure; longer if Newton-constant identification is treated rigorously.
**Builds on:** Sprint TD Track 4 (`debug/sprint_td_track4_memo.md`, T_H reproduction via M1 on the τ-circle); Track D Bisognano–Wichmann reading; Sprint Unruh-pendant (four-witness theorem).

---

## §1. Infrastructure already in place

The temperature side is closed. Sprint TD Track 4 used `geovac/thermal_tensor_triple.py::matsubara_spectrum()` verbatim with $\beta_{\rm cigar}=8\pi M$, reproducing $T_H=1/(8\pi M)$ from the M1 Hopf-base measure mechanism. Paper 32 §VIII `rem:bisognano_wichmann_reading` and Paper 34 §V.B Hawking row both anchor this as machinery-witness with `M` external input. The four-witness Wick-rotation theorem (Hawking+Sewell+BW+Unruh) supplies the bridge to Lorentzian-side modular flow at the level of correlation-function periodicity.

What is reusable for the entropy side: (a) `matsubara_spectrum()` and `modular_residual_thermal_tensor()` factorization $\mathrm{Tr}\,e^{-sD^2}=K_{\rm spatial}(s)\cdot K_{\rm temporal}(s)$; (b) the master Mellin engine partition (M1/M2/M3) and the case-exhaustion theorem (Paper 32 §VIII Thm `thm:pi_source_case_exhaustion`); (c) the Seeley–DeWitt machinery in `geovac/qed_vacuum_polarization.py::seeley_dewitt_coefficients_s3()` — but **only on S³**, with the closed-form $a_0, a_1, a_2$ for the squared Dirac on the round 3-sphere via Branson–Gilkey / Vassilevich coefficients.

**The entropy side needs a different Seeley–DeWitt construction (on the horizon S² as a 2-manifold) and a Newton-constant identification, neither of which exists in the framework.**

## §2. Horizon S² spectral data — the load-bearing missing piece

Search results: there is **no S² Dirac spectrum, no S² Seeley–DeWitt machinery, no 2-manifold spectral triple construction anywhere in `geovac/`**. The S² content that exists is (a) Hopf-base content in Paper 25 §IV ($G^{S^2}$ quotient of the Fock graph with Laplacian spectrum $\{0,0,0,1,3,6\}$ at $n_{\max}=3$); (b) S² Gaunt integrals in `geovac/so4_three_y_integral.py` and `geovac/operator_system.py` (Wigner-3j angular factors). Both are *graph-quotient* / *Gaunt-coupling* constructions, not the continuum spectral-triple data needed for Connes–Chamseddine on the horizon.

On a 2-manifold, the canonical Seeley–DeWitt coefficients for $D^2$ on a spin-2 fermion are $a_0 = \mathrm{Area}/(4\pi)$ (dimensional, from $(4\pi t)^{-d/2}$ at $d=2$) and $a_2(D^2) = (1/(48\pi))\int R \,d\mathrm{vol} - (1/(8\pi))\int E\,d\mathrm{vol}$, with Lichnerowicz $E = R/4$. For round S² of radius $r_h$, $R = 2/r_h^2$ and $\mathrm{Area}=4\pi r_h^2$. These need to be re-derived from `qed_vacuum_polarization.py`'s framework with $d=3 \to d=2$ — a 1–2 week module-level refactoring task, but conceptually unproblematic (closed-form Branson–Gilkey on a constant-curvature 2-sphere).

## §3. Area-law prefactor accounting

The Connes–Chamseddine spectral action on a 4-manifold expands as $S_{CC}[\Lambda] = \mathrm{Tr}\,f(D^2/\Lambda^2) \sim \Lambda^4 a_0 + \Lambda^2 a_2 + a_4 + O(\Lambda^{-2})$ with $f$-moments $f_4,f_2,f_0$ depending on the cutoff function. The horizon-area term is the leading bulk Einstein–Hilbert piece via the Chamseddine–Connes identification $f_2\Lambda^2/(96\pi^2)\int R \,d\mathrm{vol}_4 \leftrightarrow (1/(16\pi G))\int R\,d\mathrm{vol}_4$, giving $\Lambda^2 \sim 1/G$ (Chamseddine-Connes 1997, J. Geom. Phys. 22). Wald-formula evaluation at the bifurcation surface then reads $S = A/(4G)$. The numerical match that would constitute closure: framework computes the spectral-action coefficient on the Euclidean Schwarzschild saddle, identifies the area term via Wald, and recovers $S_{BH} = A/(4G)$ with a single calibration choice for $\Lambda^2$ in terms of $G$. **One free parameter ($\Lambda^2\leftrightarrow G$) is unavoidable** — Newton's constant is calibration input on the GeoVac side, exactly parallel to how the cigar's `M` is external in Track 4.

## §4. Named obstacles

(i) **No horizon-Dirac construction (the dominant gap, 3–6 weeks).** Need to construct $D_{S^2_{r_h}}$ from scratch on the horizon, then assemble the full spectral triple $\mathcal{T}_{\rm horizon}=(\mathcal{A},\mathcal{H},D)$ for Connes–Chamseddine evaluation. This is mathematically standard but framework-internally absent.

(ii) **Newton-constant identification (1–2 weeks, but calibration not derivation).** $\Lambda^2 \leftrightarrow 1/G$ is one input parameter — accept it cleanly as a §V.D convention exposure or as the gravity-side W3 inner-factor analog. Honest scope per CLAUDE.md §1.5.

(iii) **Wick-rotation factor accounting (1 week).** Track D's structural-correspondence-not-literal-identification verdict applies: Euclidean spectral action on the cigar gives a saddle-point free energy $F$; Wick rotation reads $S = \beta F - (\beta\partial_\beta - 1)F$ to extract entropy. The framework computes Euclidean side; Lorentzian-side identification is the Wick chain.

(iv) **Non-compact $r$-direction (1 week, regularization choice).** Track 4 explicitly flagged the cigar's $r\in[2M,\infty)$ as outside the master Mellin engine's compact-spectrum scope. For the entropy computation the horizon S² IS compact, so the obstacle reduces to IR-regularizing the $r$-integration in the saddle evaluation (standard QFT-in-curved-spacetime: finite spacelike box, take limit).

(v) **Master Mellin engine M2 transferability (structural, partial).** Sprint TD Track 4 §6 explicitly stated "M2 needs case-by-case work for non-compact / non-Coulomb manifolds." Branson–Gilkey on S² is closed-form; the integration over horizon vs. non-compact bulk is the genuine new content.

## §5. Sprint length estimate

| Phase | Task | Weeks |
|:------|:-----|:-----:|
| 1 | S² Dirac + Seeley–DeWitt on the horizon (module refactor of `qed_vacuum_polarization.py` to $d=2$, Branson–Gilkey closed forms) | 1–2 |
| 2 | Spectral action evaluation at Schwarzschild saddle (Euclidean cigar geometry, area term extraction via Wald) | 2–3 |
| 3 | Newton-constant calibration $\Lambda^2\leftrightarrow 1/G$ + bookkeeping for $f_2$ cutoff moment | 1–2 |
| 4 | IR regularization of non-compact $r$-direction + horizon saddle | 1–2 |
| 5 | Verification (sympy residual zero on $S_{BH}=A/(4G)$ after one-parameter calibration) + Paper updates (32 §VIII, 35 §VIII, 34 §V.B Bekenstein row) | 1–2 |
| 6 | Buffer (literature alignment with Chamseddine–Connes 2010 black-hole spectral-action literature, Marcolli–vS lineage check) | 2–3 |
| **Total** | | **8–14** |

Estimate dominated by (1)+(2); (4)+(5) are routine once (1)+(2) land. If the agent insists on operator-level (not calibration-level) derivation of $G$ from internal spectral data, that becomes an unbounded W3-class problem — refuse.

## §6. Feasibility verdict

**GO-WITH-PREREQUISITES.** The horizon-Dirac construction (§4(i)) is the named precursor sprint; everything else is routine given that. The "GO" reading rests on three load-bearing observations: (a) infrastructure on the temperature side already factorized cleanly per Track 4; (b) Seeley–DeWitt machinery already exists for $d=3$ Dirac on S³, generalizing to $d=2$ on S² is a refactor not a new theory; (c) the Newton-constant calibration is one acknowledged-external parameter, parallel to `M` in Track 4 — accepting it preserves the structural-skeleton-scope discipline (CLAUDE.md §2 pattern crystallization).

The "WITH-PREREQUISITES" reading: this is **not** an entry-level 2-week sprint. Approximately 8–14 weeks to a clean first-pass closure. The result, if it lands, is a Paper-32 §VIII addendum or a standalone Paper 42 entry; it does NOT promote $S_{BH}$ to a GeoVac prediction — it is at best a verification that the master Mellin engine's M2 sub-mechanism on a 2-manifold reproduces the area law under standard Connes–Chamseddine bookkeeping, with one calibration parameter (Newton's constant). This is consistent with how the cigar's `M` was treated in Track 4: machinery-witness, error class C, no new physics.

**Not NO-GO** because no structural blocker requires full Lorentzian extension first — the Euclidean spectral action on the cigar is well-defined Riemannian geometry, and the Wick chain to Lorentzian entropy was already validated at the temperature level in Track D. The Lorentzian-NCG scoping memo (`debug/lorentz_boost_scoping_memo.md`) confirmed no published Lorentzian propinquity sequel as of May 2026 — but the BH entropy computation is *Euclidean*, so this is not blocking.

## §7. Recommendation

Defer until at least one of the following triggers: (i) a precision-physics observable becomes load-bearing for an entropy-side BH check (none currently identified); (ii) the master Mellin engine's M2 sub-mechanism on non-Coulomb manifolds becomes a structural target (it is currently flagged as future work in TD Track 4 §6); (iii) the PI wants a paper-scale deliverable in the gravity direction. Otherwise the chemistry-solver W1c-residual sprint, the BBB93/KTT screening upgrade, or the V.C autopsy backlog are higher-value uses of an 8–14 week window.

If pursued, scope as a single track owned by the PM with explicit pre-implementation diagnostic at week 0 (per `feedback_diagnostic_before_engineering.md`): can the S² Seeley–DeWitt module land in ≤2 weeks as the agent claims? If no, escalate before committing.
