# Paper 59 — carryforward (remediation + synthesis sprint owed before cert)

> **Origin: the 2026-08-17 first-cert FULL `/qa paper 59` run = FAIL (trustworthy).**
> Panel fully calibrated (6/6 seeds, 0/5 controls); deterministic C10–C18 all green.
> PI direction 2026-08-17: **scope the whole remediation as its own sprint + synthesis
> work** (do NOT fix inside the /qa run). This file is that scope. Cert record +
> compact findings live in `docs/qa/paper_59.done.md` (frozen DoD + change log);
> seed key in `debug/qa/paper_59_seed_key.json`.

**Cert path:** do Parts A–C → **DELTA-verification `/qa` run** (diff-scoped, seeded) →
once clean-delta, the **FULL certifying `/qa` run** (with the completeness-critic that
was deferred on the FAIL run) → PASS. Part D (synthesis) is parallel and PI-gated.

---

## Part A — Backing cluster (the substantive work: `tests/`)

The paper's *qualitative structural* claims are soundly backed (122/122 tests pass:
elliptic-vs-rational curve, exact-form g(1)=0, indicial exponents, irreducibility
architecture + the ρ=1/5 integer monodromy, λ=1−ρ family, CM Γ-values). The gap is that
several *quantitative headline numbers* are witnessed only in `debug/` research-stage
drivers or at float64 precision, and two headlines have **no** pytest backing. Pin them:

- **A1 [MATERIAL — C1/C8] Collinear value at 19 digits.** Paper §modular/§bessel_algebra
  claim `0.3953557659017139641` (~19 digits) via the Duffy `ρ=σ²` spectral corner, but no
  test witnesses the *assembled* T2 at 19 digits — `test_paper59_corner_sigma2.py` pins only
  the corner sub-triangle T1(δ=0.05) to ~13 digits, and `test_routeC_momentum.py::test_cosmic_galois_integrated_value`
  pins the full value only to ~6 (float64). **Add** a high-precision (mpmath) test that
  assembles corner + bulk and pins T2 to ~19 digits (port from `debug/routeC_fast_evaluator.py`
  / the corner_sigma2 driver). **Then reconcile `docs/claim_test_matrix.md` row 408** — it is
  STALE (says "~16 digits, independent certification ~6", cites only the ~6-digit test; must
  become 19 digits + the new corner backing). *Alternatively,* if the 19-digit assembly is not
  cheaply pinnable, down-state the paper's digit claim to what a test witnesses — but the
  matrix and paper must agree.
- **A2 [MATERIAL — C1] PSLQ-negative has no test.** The "guarded integer-relation search with
  a same-magnitude decoy control = NEGATIVE" headline (§modular + §bessel_algebra) is backed
  only by `debug/routeC_cosmic_galois_rung3b.py` (matrix row 408). **Add** a regression test
  that runs the guarded PSLQ + decoy and asserts the search-negative (height-grows /
  decoy-matched) at the reachable precision. Keep the tier honest: it pins a *search-negative
  at ~19 digits*, NOT a proof.
- **A3 [NIT] eq:K0 (6e-18) and eq:period (31-digit)** are witnessed only by float64
  `scipy.integrate.quad` at 1e-10 / 1e-9 (`test_single_dispersion_factor_is_besselK0`,
  `test_three_center_kernel_is_elliptic`). **Upgrade** to mpmath high-precision so the tested
  tolerance matches the quoted digits.
- **A4 [NIT] eq:pf residual ~1e-24.** No test plugs the actual master period N(D) into the
  L₄ operator and checks annihilation (existing L4 tests cover self-adjointness / exponents /
  factorization / monodromy, and `test_single_bessels_and_periods_are_not_D_solutions` only
  checks the *wrong* candidates are NOT annihilated). **Add** the direct L₄[N]≈0 residual test.
- **A5 [NIT] ρ=1/3 monodromy leg.** `test_L4_monodromy_numerically` hardcodes `w2=4 # ρ=1/5`;
  the paper's "same integer matrix for ρ=1/5 AND ρ=1/3 (topological invariant)" is asserted in
  a comment only. **Add** the ρ=1/3 numerical re-derivation.
- **A6 [NIT] Slater-FT 4.3e-7 / reduced-form 1.8e-6** (§momentum) have no assertion. **Pin** both.
- **A7 [NIT] tolerance/prose drift** on several existing tests (e.g. `test_modulus_pf_annihilates_periods`
  actual ~1e-31 vs quoted 1e-18; `test_cosmic_galois_cm_periods` actual ~4.6e-41 vs quoted ~1e-51;
  `test_modulus_source_is_in_module` 1e-25 tested vs 5e-44 quoted). Tighten tolerances to match
  the quoted digits, or state the tested tolerance in the matrix notes.

## Part B — Small paper + matrix fixes

- **B1 [MATERIAL — C3/C8/branch-crit-5] §scope line 779:** "the **proven** quadratic relations
  {−π,0,2π}" → "the **measured** (25-digit) quadratic relations {−π,0,2π}" (they are
  `[MEASURED]`; the intersection-form identification is `[OBSERVATION]`, "one step short of a
  theorem" — §bessel_algebra + matrix row 410 already tier it correctly).
- **B2 [NIT] provenance visibility:** the paper's `\bibitem{test_routeC}` cites only
  `test_routeC_momentum.py` + `test_two_center_eri_aabb.py`; add `test_paper59_bessel_moment_algebra.py`
  and `test_paper59_corner_sigma2.py` (they back §bessel_algebra + §modular). Add matrix rows
  if any A-item creates new tests.

## Part C — Citation NITs (`citation-reviewer` verified; cosmetic)

- `ozdogan2012` title "Slater-**type** orbitals" → "Slater orbitals" (actual).
- `fromm_hill1987` title "**the** three-electron integral" → "three-electron integral**s**".
- `broadhurst2016` — add the published journal home (*Commun. Number Theory Phys.* **10**,
  527–569, 2016) alongside the arXiv id (optional; arXiv-only is legitimate).
- Orphan bibitems: `loutey_paper18`, `loutey_paper58` are referenced as bare prose "Paper 18/58"
  (never `\cite`'d); `loutey_paper34` is fully orphaned (never mentioned). Either `\cite` them at
  their supporting claims or drop the unused bibitems.
- *(Not a defect — recorded for the record: the frontier amplitudes / number-theory / cosmic-Galois
  citation surface was verified clean against primary sources; the Katz FL-preserves-simplicity
  attribution — the load-bearing leg of the irreducibility proof — is topically sound.)*

## Part D — C9 synthesis (PI-gated; parallel to A–C)

> **STATUS 2026-08-17: DONE (PI-directed).** The fold was found **already present** in the
> working tree (uncommitted; added by a v4.85–88 sprint — last commit predates Paper 59):
> two blocks in `papers/synthesis/group3_foundations_synthesis.tex` — the §Tannakian
> "elliptic layer (Paper 59)" paragraph (~L1073) and the forward-look "elliptic layer of the
> cosmic Galois" subsection (~L1213). Both are current and tier-honest (they already carry the
> v4.86 L₄-irreducibility proof + the "corroborates rather than lifts the Hain–Brown negative"
> framing). One clause brought up to the paper's current placement: "candidate Γ(2) multiple
> modular value" → "candidate Eisenstein/CM **Bessel moment** (not a cusp-form L-value), of
> which a Γ(2) MMV is only the regular D→0 shadow." Consistency with Paper 56 `rem:paper59_cm`
> verified; synthesis compiles (13 pp). **D2 C9 surface now live** — a future `/qa group3` (or
> paper-59-recert) must exercise C9 on it. *Pre-existing, out-of-scope:* one latent natbib
> "Citation `Note1' undefined" at L661 (a footnote artifact, no literal `\cite{Note1}`; not
> introduced by the fold) — flag for a group3 pass.

- **D1** Paper 59 has **no synthesis footprint** today (why C9 was N/A this cert). PI direction
  2026-08-17: **fold Paper 59 into the group3-foundations synthesis**
  (`papers/synthesis/group3_foundations_synthesis.tex`) — its natural home via the Paper 56
  cosmic-Galois / Tannakian-substrate tie (Paper 56's `rem:paper59_cm` already lands the
  Q(i) CM-point corroboration there). Write current-state-only (per authoring-conventions
  rule 11): the genus-1 relocation, the proven L₄ irreducibility, the OPEN transcendental
  closed form, the Eisenstein/CM Bessel-moment period placement.
- **D2** Once folded, the group3-foundations synthesis gains a Paper-59 C9 surface → a future
  `/qa` (group3 or paper-59-recert) must exercise C9 on it.

---

## Change log
- 2026-08-17 — Created from the first-cert FAIL run (PI: scope as its own sprint + synthesis).
- 2026-08-17 — **Parts A–C DONE** (remediation sprint). New tests: A1
  `test_paper59_assembled_collinear_value` (assembles corner+bulk on an independent
  tiling; agrees with the headline 0.3953557659017139641 to ~12 digits at deg4/M384,
  cross-validated to ~16 at deg5); A2 `test_paper59_pslq_negative_multiprecision_guarded`
  (guarded weight-3 disc-4 + disc-8 multi-precision stability search-negative + decoy + positive control);
  A4 `test_paper59_L4_annihilates_the_physical_master_N` (direct eq:pf annihilation on
  N(D)=s_K, ~1e-18); A5 monodromy parametrized → ρ=1/5 AND ρ=1/3 both re-derived;
  A6 `test_momentum_true_slater_density_ft` (4.3e-7) + `test_momentum_angular_reduced_j0_form`
  (1.8e-6); A3 eq:K0 + eq:period upgraded float64→mpmath (dps=40, exact); A7 modulus_pf
  1e-10→1e-25, cm-periods dps=40→55 tol 1e-30→1e-50, in-module tolerance clarified in the
  matrix note. Paper: B1 §scope "proven"→"measured (25-digit)" quadratic relations; B2
  test_routeC bibitem expanded to cite both `test_paper59_*` files; C ozdogan/fromm_hill
  titles fixed (primary-source verified), broadhurst2016 journal home added (CNTP 10,
  527–569, 2016; verified), orphan bibitems P18/34/58 `\cite`d at their support points.
  `claim_test_matrix.md` rows 405/406/407/408 reconciled + new §momentum row. Paper compiles
  clean (9 pp, no undefined cites). Ready for the DELTA-verification `/qa` run.
- 2026-08-17 — **Part D DONE** (PI-directed "b"). Group3-foundations synthesis fold found
  already present (uncommitted v4.85–88 work) + brought current (Bessel-moment placement clause);
  Paper 56 tie verified; compiles (13 pp). D2 C9 surface now live for a future `/qa group3`. One
  pre-existing natbib `Note1` warning flagged (out of scope).
- 2026-08-17 — **FULL certifying `/qa` run = FAIL (trustworthy) → MATERIAL remediated.** Panel
  fully calibrated (7/7 seeds, 0/5 controls) incl. the C9 synthesis dimension. One MATERIAL: the A2
  PSLQ test under-witnessed its `[MEASURED]` headline (searched only weight-≤2 disc-4; paper claims
  weight-3 + disc-8 + height-grows). Remediated: extended to the multi-precision stability guard
  (`test_paper59_pslq_negative_multiprecision_guarded`) covering weight-3 + disc-8 + decoy + positive
  control; broedel co-author NIT fixed. Full record in `paper_59.done.md`. **Next: PI-invoked DELTA
  `/qa` to verify the fix → clean-delta → FULL certifying run → PASS.**
- 2026-08-17 — **DELTA run = DEFECTS → FALLBACK (PI-approved).** The extended multi-precision PSLQ
  test had a precision-management bug (the delta caught it); fixing it exposed that a decisive
  weight-3/disc-8 test is provably impossible at the value's ~19-digit precision (pslq needs ≥16 digits;
  the ~20-element basis needs ~40; extending the value = the ~40-digit frontier). Fallback: permanent
  test = the decidable **disc-4 weight-≤2** negative (`test_paper59_pslq_negative_disc4_guarded`); paper
  §modular/§bessel_algebra + C8 #7 + matrix re-scoped so weight-3/disc-8 read as bounded/driver-observed
  (decisive test needs ~40 digits). Compiles clean. **Next: a fresh DELTA `/qa` verifies the fallback.**
