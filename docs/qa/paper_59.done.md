<!-- CERT-STALENESS-BANNER -->
> ### ⚠ RE-CERTIFICATION OWED
> This record certifies the state as of **2026-08-21**. Since then **1 `.tex` changed**: paper_59_elliptic_bessel_moment.tex.
>
> **The CERTIFIED verdict below is therefore historical, not current.** Do not cite it as present-tense status.
>
> Re-measure rather than trusting this banner — it is itself a snapshot and will go stale the same way:
> `python debug/qa/check_cert_staleness.py --detail`
<!-- /CERT-STALENESS-BANNER -->

# Paper 59 (Elliptic Bessel Moment) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only the Paper-59 scope + deltas + watch-notes.

> **STATUS: RE-RUN 2026-09-07 — clean on this paper; NOT re-certified.**
> `/qa` (PI-invoked). Deterministic layer green on scope `paper_59`; no
> content defect found in the paper itself. One **instrument** gap closed:
> **C21 examined ZERO annotations here** — a PASS carrying no information.
> The salience report was no help, since it lists only multi-document
> numerals and its four "measurement-shaped" candidates on this paper were
> bibliography volume/page numbers (145, 151, 286, 376). The real quantity
> had to be read out of the text: the three-center ERI ground truth
> $(XY|XZ) = 0.204941722$, now registered and annotated (C21: 1 annotation
> checked). `\gvq` was undefined in this file and has been added — caught by
> C10 before it shipped. Done-record gates were running under `--gate group2`,
> which does not contain this paper; now `--gate paper_59`.
>
> *(superseded — historical)* **CERTIFIED ✅ 2026-08-17 — FULL certifying run #4 = PASS** (calibrated 8/8 sens,
> 0 FP; four dimensions incl. C9; completeness-critic clean; zero genuine non-seed MATERIAL).
> Cert arc: first-cert FULL = FAIL → DELTA #1 fallback = CLEAN-DELTA → FULL #2 = FAIL (PSLQ
> under-witness) → remediated → FULL #3 = FAIL (thin: three-master eq:pf coverage) → remediated →
> FULL #4 = **PASS**. Honest ceiling in the change log. (Original freeze note retained below for history.)
>
> **[FROZEN 2026-08-17 — first cert of a fresh target (FULL run).]** Paper 59
> was created at v4.82.0 (2026-08-16), after both the group2 cert (2026-06-28) and
> the v4.76.0 close-out freeze; it has never been through `/qa`. This run discharges
> the v4.84–v4.88 "Phase-4 re-review OWED" debt logged in CLAUDE.md §2. Framing
> decisions PI-approved 2026-08-17 (three-question freeze): (1) **C9 = N/A** — Paper 59
> has no synthesis footprint (see Scope); (2) branch-defining criterion **swapped** from
> group2 chemistry-benchmarking to **transcendence-tier + open-frontier honesty**;
> (3) proceed straight to the FULL first-cert run.

**Scope (single-paper):**
- **Paper 59** — `papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex`
  ("The Three-Center Electron-Repulsion Integral Is an Elliptic Bessel Moment:
  Genus One at the Third Center in Momentum-Space Slater Theory").
- **C9 (synthesis) = GATING (corrected 2026-09-07).** This record said until then:
  *"N/A this run. Paper 59 has no footprint in the group2 synthesis (verified: the only
  '59' hit there is an unrelated rotational constant $B_e=59.5\,\mathrm{cm}^{-1}$)."*
  **That premise was true when frozen (2026-08-17) and is now false.** The group2
  synthesis gained `\subsection{The transcendence frontier at the third center
  (Paper~59)}` (L666–681) plus abstract, open-questions and bibitem loci on 2026-09-06
  (v5.10.6). A scope exclusion inherited from a premise dies with the premise: C9 is a
  **gating dimension** for Paper 59, it was **not exercised** in the FULL run of
  2026-09-07, and that alone forced INCONCLUSIVE. The class is general — nothing in the
  apparatus audits a frozen DoD against the paper it certifies, so a stale premise reads
  from outside exactly like a passed criterion.
  Paper 56's `rem:paper59_cm` remark stays a lightweight **C7** check as before.
  - **STANDING FOLLOW-UP (PI direction 2026-08-17, NOT part of this run):** fold Paper 59
    into the **group3-foundations synthesis** (its natural home via the Paper 56
    cosmic-Galois / Tannakian-substrate tie). To be done as a separate reviewed task
    *after* this cert — adding synthesis content mid-`/qa` would re-open scope and would
    itself need a C9 pass.
- **Out of scope:** all other group2 papers (unchanged since the 2026-06-28 cert) and
  Paper 56/58 (certified). Trunk papers (0/1/7/14/18/34) canonical; in scope only where
  Paper 59 restates them (C7).

**Deterministic `--gate`:** `paper_59`.
  Note the scope is the SINGLE-PAPER one, not `group2`: `qa_scopes.py` deliberately excludes 58/59/60 from `group2` (they are their own cert targets), so the `--gate group2` this record carried until 2026-09-07 examined NOTHING for the paper it certifies.

## Dimensions exercised (ALL, one invocation — FULL run, first cert of a fresh target)

- **Code / test-backing (C1–C2)** — `code-reviewer` ×1 on Paper 59 (Sonnet-tier;
  2 seeds). Map + **RUN**: `tests/test_routeC_momentum.py` (momentum reproduction of
  $(XY|XZ)$; the single-dispersion $K_0$ identity; the elliptic-curve period),
  `tests/test_paper59_bessel_moment_algebra.py` (the rank-4 master periods, the
  Wronskian $W_0 D^{-2}$, the quadratic relations $\{-\pi,0,2\pi\}$, the L₄
  irreducibility battery / monodromy certificate), `tests/test_paper59_corner_sigma2.py`
  (the $\rho=\sigma^2$ spectral corner → the collinear value, since certified to 66
  digits via eq:kw + the guarded
  PSLQ-negative), and `tests/test_two_center_eri_aabb.py` (the two-center weight-one
  census the genus-zero baseline rests on). For each: does the test genuinely *prove*
  the claim (not tautological / weaker than prose / right-answer-wrong-reason)? Are the
  MEASURED residual tolerances (1.9e-14, 6e-18, 31-digit period, ~1e-24 PF, 5e-44
  in-module, 25-digit relations, integer monodromy) actually asserted?
- **Paper claims / prose (C3, C5, C6, C7, C8 + branch criterion)** — `claims-reviewer`
  ×1 on Paper 59 (Opus-tier; 2 seeds), enumeration-forced: every `[TIER]` tag, every
  headline number, every `[OPEN]`/`[OBSERVATION]` marker, every transcendental.
- **External citations (C4)** — `citation-reviewer` ×1 on Paper 59 (Sonnet-tier;
  2 seeds). **HIGH-FABRICATION-RISK SURFACE** — the largest, most specialized citation
  surface in the corpus (frontier amplitudes + number theory the PM has not
  independently grounded). Prioritize the **load-bearing** cites (below).
- **Synthesis faithfulness (C9)** — **GATING (corrected 2026-09-07; this line
  said N/A until then).** The v5.10.6 relocation gave the group2 synthesis a live
  Paper-59 surface, so the premise the N/A rested on is gone. See the C9 entry in
  the criteria block above for the full correction; the scope now carries
  `synthesis/group2_quantum_chemistry_synthesis.tex` and `--gate paper_59`
  resolves to 2 files.
- **Deterministic (C10–C18)** — the step-1 scripts, `--gate paper_59`.
- **Completeness-critic** ×1 (FULL run).

## Branch-defining criterion (SWAPPED, PI direction 2026-08-17): transcendence-tier + open-frontier honesty

This replaces group2's chemistry-benchmarking/guardrail-negative criterion (which is
vacuous for a differential-Galois / cosmic-Galois paper). It is a **sharpening of the
shared C3 (prose ≤ tier) + C8 (headline honesty) + C5 (no §3 negative suppressed)**, not
a new number. The reviewers must verify ALL of:

1. **SYMBOLIC PROOF means a proof, not a match.** The two strongest claims —
   **L₄ irreducibility** and the **fourth-order Picard–Fuchs equation** — carry
   `[SYMBOLIC PROOF]`. Irreducibility rests on TWO independent legs: (i) the
   Fourier–Laplace argument (FL of an irreducible rank-1 connection is simple, Katz 1990)
   and (ii) the explicit integer-monodromy certificate $M_0$. The prose may assert
   irreducibility as proven ONLY on the strength of these; a `[SYMBOLIC PROOF]` tag
   backed only by a numerical residual (a MEASURED fit) is an overclaim = MATERIAL. The
   self-adjointness ($L^*=L$), the exponents $\{0,1,1,2\}$, and the eigenring
   $E(L)=\mathbb{C}$ are the supporting SYMBOLIC/computed facts.
2. **Reducibility-negative ≠ closed-form-impossibility.** The paper proves
   **no closed form via factorization** (the operator is irreducible). It does NOT prove
   the object has no closed form: a **transcendental** closed form at the irreducible
   level (elliptic-polylog / $\Gamma(2)$ MMV) is explicitly `[OPEN]`. No sentence may
   collapse "irreducible ⇒ no factorization closed form" into "no closed form exists" or
   "the closed form is impossible" = MATERIAL.
3. **The finite closed form is NOT delivered.** The paper gives a *relocation* + a
   *diagnosis of hardness* + an *accurate/fast evaluator*, not a solution. No prose may
   imply the three-center integral is solved in closed form, or that a polyatomic energy
   is delivered.
4. **The PSLQ is NEGATIVE — never dressed up.** The guarded integer-relation search
   (with a same-magnitude decoy control) is **negative**: the value is *not* a low-height
   classical CM-Γ period at reachable precision; the relation height *grows* with
   precision and is *matched by the decoy*. This is "no low-height relation found,"
   NOT "proven no relation exists" and NOT a positive identification / transcendence
   proof. *(Updated 2026-08-21, v4.106.x delta remediation: the (k,w) refactorization
   certified 66 digits and the guarded searches are now DECISIVE through weight three
   including the Catalan-completed ring at height <=10 — height-bounded negatives, still
   never a positive ID or transcendence proof.)* Any framing of
   the PSLQ-negative as a positive ID, or as a proof of transcendental independence, or
   any surviving stale "digit-9 anchor-vs-composite discrepancy" (a **false alarm**,
   RectB GL under-resolution; anchor CONFIRMED per CLAUDE.md §2) = MATERIAL.
5. **OBSERVATION tiers stay observations.** "The bridge appears unmade" (literature
   novelty) is an *absence from a targeted search, not a proof of absence* (the paper
   says so — verify the hedge survives). The cosmic-Galois / mixed-elliptic-motive
   *placement* is an OBSERVATION, not a theorem. *(Goalpost updated 2026-08-21,
   v4.106.x delta remediation, superseding the frozen framing:)* the quadratic-relation
   intersection form was CLOSED SYMBOLICALLY in v4.97.0 (B = pi*Omega forced; ~~Galois in
   Sp4(Z)~~; backed by `test_intersection_form_*` legs) — [SYMBOLIC] prose is now correct,
   not an inflation; only the *placement* (Eisenstein/CM vs cusp-form) stays OBSERVATION.

   > **CORRECTED 2026-09-07 (`/qa paper_61`, two reviewers converging independently).**
   > The struck clause was FALSE, and this goalpost *ratified* it — which is why
   > it is graded LARGE: a `.done.md` is what the next certifying run measures
   > against, so a reviewer reading this line would have restored the defect on
   > the file's own authority. Correct form: **monodromy in Sp4(Z), differential
   > Galois in Sp4(C)**. Sp4(Z) is discrete, so a Zariski-closed subgroup of
   > GL_4(C) inside it is finite, forcing every solution algebraic — contradicting
   > the irregular singularity at infinity (Poincare rank 1) and the exponential
   > torus (C*)^2 that Paper 59:607 establishes. The backing variable is literally
   > `_L4_MONODROMY` and the test asserts `M0^T Omega M0 = Omega`, a monodromy
   > statement. **Nothing numeric changes** — B = pi*Omega and the integral
   > structure stand; only the group the containment is a statement ABOUT was
   > wrong. Registry entry: `p61-galois-in-sp4-z`.
6. **exact ≠ accurate (inherited from Paper 58 W1).** Closed-form / π-free / weight-one /
   momentum-native-evaluator content is a *decidability / diagnosis* result, NOT an
   accuracy improvement. §sec:scope must carry this; any "exact ⇒ better energy / more
   accurate molecule" reading = MATERIAL.

## Paper-59-specific watch-notes (the risk surface — ranked)

- **W1 — irreducibility overclaim [HIGHEST RISK; criterion 1/2].** The abstract +
  §obstruction + §scope + conclusion all restate the irreducibility result. Verify each
  restatement (a) attributes it to the FL/monodromy proof, not the numerical fits, and
  (b) keeps the transcendental closed form OPEN. The conclusion's "we prove irreducible …
  so it does not close by any factorization … whether it closes at the irreducible level
  … is the one frontier that remains" is the correct shape — hold every restatement to it.
- **W2 — the Katz FL attribution [load-bearing citation].** `katz1990` (Exponential Sums
  and Differential Equations) is cited for "FL is an exact auto-equivalence of holonomic
  D-modules on the affine line, preserving simplicity." This single citation carries the
  "in principle" leg of the headline irreducibility proof. It MUST be verified against the
  primary source. A wrong attribution here is MATERIAL (it is the difference between a
  proof and an assertion).
- **W3 — PSLQ-negative honesty [criterion 4].** §modular + §bessel_algebra both state the
  negative. Verify neither drifts to a positive claim. *(Canonical accounting updated
  2026-08-21, v4.106.x delta remediation:)* the certified value is
  **T2 = 0.395355765901713964325229296804847564…** (**66 digits**, decomposed certification,
  eq:kw; executable witness `test_kw_mpmath_witness` via `geovac/t2_kw.py` at ~21 digits);
  the pre-(k,w) anchor **0.3953557659017139641 is correct to 18 digits only** (digit 19
  superseded) and is guarded by the C17 family `t2-collinear-anchor`. The searches are
  DECISIVE through weight three incl. the Catalan-completed ring at h<=10 (64 working
  digits); height budgets scale 10^(D/n) in ring dimension (the old "~32-40 digits
  decides" was never true for the dim-20 ring).
- **W4 — genus grading [transcendental tagging, Paper 18/34].** The embedding tier is
  *genus-graded*: genus-0 $\{E_1,\ln,\gamma\}$ at two centers, genus-1 (elliptic) at
  three. Verify: the two-center weight-one negatives (no π, no dilog) are stated as
  genus-zero properties; the genus-1 claim is exactly one genus up (not "genus two");
  π enters via the period / quadratic relations; the Γ-values are CM-fibre periods; no
  anonymous transcendental. Paper 18 §Level-2 + Paper 34 carry the genus grading
  (CLAUDE.md §2) — C7 consistency.
- **W5 — modulus-map / CM-fibre claims [MEASURED / classical].** $\lambda(\tau(\rho))=1-\rho$
  (10⁻⁴¹), the Legendre/Γ(2) universal-family identification, and the Chowla–Selberg
  Γ-values (ρ=1/2 → Γ(1/4)²/(4√π); disc−8) are MEASURED/classical. Verify the honest
  deflation is intact: Γ(2) is the *universal* level (below the sunrise's Γ₁(6)); the
  content is the rational modulus map + in-domain CM fibres, "not an exotic level."
- **W6 — quadratic period relations [SYMBOLIC Wronskian + MEASURED relations].**
  $W(D)=W_0 D^{-2}$ is [SYMBOLIC] (Abel). $B[s_K,s_I]=-\pi$, $B[s_K,s_J]=0$,
  $B[s_I,s_J]=2\pi$ are [MEASURED, 25 digits]. The Broadhurst–Roberts/FSY identification
  is OBSERVATION. *(Renamed 2026-09-07: Broadhurst–**Mellit** names the
  DETERMINANT formulae; the QUADRATIC relations are Broadhurst–**Roberts**,
  proved by Fresán–Sabbah–Yu and independently by Zhou. Registry:
  `p61-broadhurst-mellit-quadratic`.)* *(Updated 2026-08-21: the intersection-form structure itself was
  closed symbolically in v4.97.0 — B = pi*Omega forced, [SYMBOLIC] is the correct tier.)*
  Verify tiers not inflated beyond these.
- **W7 — not-a-cusp-form-L-value argument [OBSERVATION].** The $\dim S_k(\Gamma(2))$
  counting (vanishes k=2,4; first 1 at k=6) + weight ≤ 3 + Eisenstein-flavoured pairing ⇒
  not a cusp-form L-value. This is a structural OBSERVATION, correctly hedged; verify it
  is not stated as a theorem and that the weight bookkeeping ("length ≤ 2 over X(2)",
  "motivic weight five") is internally consistent.
- **W8 — MEASURED validation residuals [C1/C2].** The paper's credibility rests on a long
  list of numerical residuals (1.9e-14, 4.3e-7, 1.8e-6, 6e-18, 31-digit period, 1e-24 PF,
  0 exact form, 5e-44 in-module, 2e-14→3e-30 VoP-approx, 10⁻¹⁸ modulus annihilation,
  19-digit L(0,ρ)=K(1−ρ), 10⁻¹⁷ masters, 25-digit relations, integer M₀, 10⁻⁴¹ λ,
  ~10⁻⁵¹ CM). Each MEASURED tag must be backed by a test that actually computes it; a
  `[MEASURED]` tag with no backing artifact is a coverage gap (C1) → raise if load-bearing.

## C8 headlines (enumerated, with tiers — the frozen goalposts)

1. **Momentum-space reduction [MEASURED].** $(XY|XZ)$ in momentum space (three phases,
   not three foci; eq:momentum exact) reproduces ground truth $0.204941722$ to
   **1.9×10⁻¹⁴** (closed-form Gaussian FT) / **4.3×10⁻⁷** (true Slater transform); angular
   integral → single $j_0$ (eq:angular); reduced eq:reduced to **1.8×10⁻⁶**.
2. **Genus jump [MEASURED / classical].** Single-scale $K_0$ identity eq:K0 (**6×10⁻¹⁸**);
   curve eq:curve elliptic (genus 1) for $c_1\neq c_2$, genus-0 on the diagonal; period
   eq:period = complete elliptic $K(m)$ (**31 digits**), diagonal → $\pi/(2\sqrt c)$.
3. **Picard–Fuchs [SYMBOLIC PROOF; MEASURED residual ~10⁻²⁴].** $N(D)$ = Laplace transform
   of the holomorphic elliptic differential (eq:laplace); explicit 4th-order PF eq:pf,
   char. roots = the 4 branch points; $D\ln D$ non-analyticity at $D=0$ [MEASURED].
4. **Elliptic-dilog route obstructed [SYMBOLIC PROOF; exact residual 0 / MEASURED 45
   digits].** $\mathcal{M}_\rho[Q^{-1/2}]$ = total derivative (eq:exactform, $g(1)=0$ ⇒ no
   tadpole); source closes only in-module (eq:inmodule, residual **5×10⁻⁴⁴**); VoP against
   periods only approximates (never saturates) ⇒ one level above the sunrise.
5. **L₄ irreducibility [SYMBOLIC PROOF, positive-control validated; MEASURED integer
   monodromy].** eq:pf irreducible over $\mathbb{C}(D)$ hence $\mathbb{Q}(\rho)(D)$: (i) FL
   of the irreducible rank-1 $Q^{-1/2}$ connection is simple, generic rank 4 = total drop
   (Katz 1990); (ii) monodromy around $D=0$ in the Lefschetz-thimble basis = exact integer
   unipotent $M_0$, same for $\rho=1/5,1/3$, single 2×2 Jordan block, no invariant
   coordinate subspace. Supporting: $L^*=L$ ⇒ Galois ⊆ $\mathrm{Sp}_4$; exponents
   $\{0,1,1,2\}$; eigenring $E(L)=\mathbb{C}$. ⇒ no closed form via factorization;
   transcendental closed form OPEN.
6. **Literature placement [OBSERVATION].** Three-center Slater integral (via Fock
   projection) = two-scale Bessel moment on the sunrise elliptic curve; bridge appears
   **unmade** — absence from a targeted search, not proof of absence.
7. **Modular structure + cosmic-Galois [MEASURED 10⁻⁴¹ / classical; OBSERVATION].**
   $\lambda(\tau(\rho))=1-\rho$ exact (**10⁻⁴¹**); Legendre/Γ(2) universal family, 3 cusps;
   CM fibres → Γ-values ($\rho=\tfrac12$ → $\Gamma(1/4)^2/4\sqrt\pi$; disc−8, **~10⁻⁵¹**).
   Collinear value **T2 = 0.395355765901713964325229296804847564…** (**66 digits**,
   decomposed certification via the (k,w) refactorization eq:kw, v4.104.0; the pre-(k,w)
   anchor 0.3953557659017139641 is correct to 18 digits only). Guarded PSLQ
   with decoy = **NEGATIVE**, now **DECISIVE through weight three** incl. the
   Catalan-completed corrected ring at height <=10 (64 working digits); disc-4/disc-8
   smaller rings at heights up to 1e12. *(Superseded framing retained below for history:
   re-scoped 2026-08-17 the wider legs were only "consistent" and a definitive test was
   projected to need
   ~40 digits — driver-observed, `debug/routeC_pslq_v2.py`.)*
   Cosmic-Galois placement [OBSERVATION].
8. **Bessel-moment period algebra [SYMBOLIC + MEASURED 25 digits + OBSERVATION].** 4 master
   periods (3 regular satisfy eq:pf; residual precision as stated in the paper,
   currently ~10⁻¹⁵ per-master); Wronskian $W(D)=W_0 D^{-2}$ [SYMBOLIC,
   Abel]; quadratic relations $B[s_K,s_I]=-\pi$, $B[s_K,s_J]=0$, $B[s_I,s_J]=2\pi$
   [MEASURED 25 digits]; Broadhurst–Roberts/FSY type + intersection-form = OBSERVATION /
   "one step short of a theorem"; NOT a cusp-form $L$-value (Eisenstein/CM period)
   [OBSERVATION].

## Seeding plan (worktree only; never touches the real corpus)

K ≈ 6 planted defects, ≥1 catchable by each gating dimension (code / claims / citation),
spanning the watch-notes; **tiered agents (code + citation = Sonnet) get 2 seeds each**,
claims (Opus) gets 2:
- **code ×2** (in a P59 backing test): (a) a loosened/tautological tolerance that no longer
  proves its MEASURED residual (e.g. monodromy-integer or 25-digit relation weakened to a
  trivially-passing check); (b) a right-answer/wrong-evaluation-space corruption.
- **claims ×2**: (a) an **irreducible ⇒ no-closed-form-exists** collapse (W1/criterion 2);
  (b) the **PSLQ-negative dressed as a positive ID / transcendence proof** (W3/criterion 4).
- **citation ×2**: (a) a **mis-attributed Katz FL claim** (W2 — the load-bearing leg);
  (b) a wrong-year/venue on one amplitudes cite (e.g. `fresansabbahyu2023` /
  `adamsbognerweinzierl2014`).

M ≈ 5 known-good controls that must NOT be flagged: the 1.9×10⁻¹⁴ momentum agreement, the
integer monodromy matrix $M_0$, $\lambda(\tau(\rho))=1-\rho$, the quadratic relations
$\{-\pi,0,2\pi\}$, the exponents $\{0,1,1,2\}$. Answer key →
`debug/qa/paper_59_seed_key.json`. Reviewers path-pinned to the worktree; forbidden to read
the real corpus.

## Change log
- 2026-08-17 — **DRAFTED + FROZEN** by PM, PI framing approved (three-question freeze):
  C9 = N/A (no synthesis footprint; fold-into-group3-foundations logged as standing
  follow-up), branch criterion swapped to transcendence-tier + open-frontier honesty,
  proceed to FULL first-cert run. First cert of a fresh target created post-close-out;
  discharges the v4.84–v4.88 "Phase-4 re-review OWED" debt. Inherits criteria.md C1–C18.
- 2026-08-17 — **FIRST-CERT FULL run = FAIL (trustworthy).** Panel FULLY CALIBRATED across
  all three gating dimensions: sensitivity **6/6** planted seeds caught (code S1 tautological
  Wronskian + S2 loosened B[K,I] tol; claims S3 irreducible⇒no-closed-form collapse + S4
  PSLQ-negative-as-transcendence-proof; citation S5 Fresán-Sabbah-Yu vol/year + S6 Katz
  vol/year), specificity clean (**0/5** controls false-flagged). Deterministic C10–C18 all
  GREEN. Genuine verified MATERIAL (real corpus, seeds excluded):
  (a) **[claims/branch-crit-5]** §scope L779 "the *proven* quadratic relations {−π,0,2π}"
  inflates [MEASURED, 25-digit] / "one step short of a theorem" → proven; contradicts
  §bessel_algebra + matrix row 410. (b) **[code/C1/C8]** collinear-value headline
  inconsistent + under-witnessed: paper §modular/§bessel_algebra say ~19 digits
  (0.3953557659017139641), `claim_test_matrix` row 408 STALE at ~16 digits / ~6-digit test,
  and NO test witnesses the assembled T2 at 19 digits (`test_paper59_corner_sigma2` pins only
  the corner sub-triangle ~13 digits). (c) **[code/C1]** PSLQ-negative headline has NO pytest
  backing (matrix row 408 = `debug/routeC_cosmic_galois_rung3b.py` only). NIT cluster
  (fix-on-sight, Paper-58 debug-backed precedent): eq:K0 6e-18 / eq:period 31-digit witnessed
  only at float64 quad; eq:pf ~1e-24 no L4-annihilation test; ρ=1/3 monodromy leg not
  numerically re-derived; Slater-FT 4.3e-7 / reduced 1.8e-6 unpinned; the two `test_paper59_*`
  files uncited in the paper; citation cosmetics (ozdogan/fromm_hill titles, broadhurst2016
  journal ref, orphan bibitems P18/34/58); P58-scope claims NIT. Completeness-critic DEFERRED
  to the certifying run (a FAIL does not need completeness closure; panel + matrix cross-check
  already enumerated exhaustively). Seed key `debug/qa/paper_59_seed_key.json`; worktree removed,
  no seed leaked (all 6 verified absent from real corpus; correct katz 124/1990 + FSY 17/2023
  intact). **Path to cert:** PI dispositions the backing cluster (pin the debug/ headline values
  + PSLQ + L4-residual + ρ=1/3 into tests/, reconcile matrix row 408) → remediate small tier
  words → DELTA-verification run → FULL certifying run (with completeness-critic).
- 2026-08-17 — **REMEDIATION Parts A–C DONE** (`docs/qa/paper_59.carryforward.md`). Backing
  cluster pinned into `tests/`: A1 assembled collinear value (independent tiling, ~12-digit
  witness of `0.3953557659017139641`, cross-validated ~16 at deg5), A2 guarded PSLQ-negative +
  decoy, A4 direct L₄[N]=s_K annihilation (~1e-18), A5 monodromy re-derived for ρ=1/5 AND ρ=1/3,
  A6 true-Slater-FT (4.3e-7) + reduced-form (1.8e-6), A3 eq:K0/eq:period float64→mpmath, A7
  tolerance drift fixed (modulus_pf 1e-25, cm-periods dps=55/1e-50, in-module clarified in the
  matrix). Paper: B1 §scope "proven"→"measured (25-digit)"; B2 bibitem cites both `test_paper59_*`;
  C citation titles fixed (ozdogan/fromm_hill, primary-source verified) + broadhurst2016 journal
  home + P18/34/58 orphan bibitems `\cite`d. `claim_test_matrix` rows 405/406/407/408 reconciled
  + new §momentum row. **Full slow P59 suite = 127 passed** (`test_routeC_momentum` +
  `test_paper59_{corner_sigma2,bessel_moment_algebra}` + `test_two_center_eri_aabb`); paper
  compiles clean (9 pp, no undefined cites). Ready for the PI-invoked **DELTA-verification `/qa`** run.
- 2026-08-17 — **Part D DONE** (PI-directed after A–C). The group3-foundations synthesis fold was
  found **already present** in the working tree (uncommitted v4.85–88 work): the §Tannakian
  "elliptic layer (Paper 59)" paragraph + the "elliptic layer of the cosmic Galois" forward-look
  subsection, both current (they carry the v4.86 L₄-irreducibility proof + the "corroborates
  rather than lifts" framing). Brought one clause current ("candidate Γ(2) MMV" → "candidate
  Eisenstein/CM Bessel moment, not a cusp-form L-value, of which a Γ(2) MMV is the regular D→0
  shadow"); Paper 56 `rem:paper59_cm` tie verified; compiles (13 pp). **D2: the group3-foundations
  synthesis now carries a live Paper-59 C9 surface** — a future `/qa group3` (or paper-59-recert)
  must exercise C9 on it. Pre-existing natbib `Note1` warning (L661 footnote artifact, not fold-
  introduced) flagged for that pass.
- 2026-08-17 — **SECOND FULL certifying run = FAIL (trustworthy); MATERIAL remediated.** Panel
  FULLY CALIBRATED across all four gating dimensions (C9 included this run per PI direction):
  **sensitivity 7/7** (code S1 dead-tol + S2 wrong-master `_sI`; claims S3 no-closed-form collapse
  + S4 PSLQ-as-transcendence-proof; citation S5 Katz vol/yr + S6 FSY vol/yr; synthesis S7
  lifts/fills overstatement), **specificity 0/5** controls false-flagged. Deterministic C5,C10–C18
  GREEN. **One verified MATERIAL (code/C1/C8):** the guarded PSLQ-negative regression test
  under-witnessed its `[MEASURED]` headline — it searched only weight-≤2 disc-4 (single precision),
  while the paper claims weight-3 (poly/Laurent) + disc-8-inclusive + height-grows-with-precision
  (that fuller search lived only in transient `debug/routeC_pslq_v2.py`). **REMEDIATED same day:**
  `test_paper59_pslq_negative_guarded_with_decoy` → `test_paper59_pslq_negative_multiprecision_guarded`
  — fits the raw period W=V·π/8 against weight-graded disc-4 + disc-8 monomials through weight 3, with
  a multi-precision **stability guard** (a genuine closure must be target-coeff-nonzero, low-height,
  AND survive a higher-precision residual re-check), a magnitude-matched decoy, and a **positive
  control** (planted 2π+3ϖ−P8 IS caught, so the negative is not vacuous). Diligence established the
  weight-3/disc-8 hits are precision-UNSTABLE (different vector each dps) and decoy-matched — a bounded
  ~19-digit exclusion (value cross-confirmed ~15–16 digits, best-est ~19; definitive test needs ~40).
  Paper prose stays accurate (W-polynomial fit subsumes the Laurent-in-V leg). NIT fixed on-sight:
  `broedel_eisenstein2018` dropped co-author Penante (verified 5 authors, arXiv:1803.10256).
  Completeness-critic DEFERRED (FAIL). Seed key `debug/qa/paper_59_certrun_seed_key.json`; worktree
  removed, no seed leaked. **Path to cert:** DELTA-verification `/qa` (diff-scoped, seeded — verifies
  the extended PSLQ test + broedel fix) → clean-delta → FULL certifying `/qa` (with completeness-critic)
  → PASS.
- 2026-08-17 — **DELTA-verification run = DEFECTS (trustworthy) → FALLBACK adopted (PI-approved).**
  Panel calibrated (code 2/2 seeds DA vacuous-decoy + DB trivial-control; citation 2/2 seeds DC year +
  DD arXiv; controls clean; broedel fix verified correct). Genuine non-seed defect (code): the extended
  multi-precision PSLQ test had a **precision-management bug** — `_cm_period_constants` forced
  `mp.mp.dps=dps+25`, so PSLQ ran at ~40 effective digits (not the intended grid) and the stability
  re-check was dead code. Fixing it exposed a **hard wall**: mpmath.pslq needs ≥16 digits, the value is
  trusted to only ~15–19 (extending = the ~40-digit collaboration frontier, a ~20 h computation, not a
  test), and the weight-3/disc-8 basis (~20 elements) is over-determined at ~19 digits — so a decisive
  weight-3/disc-8 regression test is **provably not achievable** at reachable precision (this IS the
  paper's own "~40 digits" wall). **Resolution (PI-approved fallback):** the permanent test asserts only
  the **decidable** disc-4 weight-≤2 negative (`test_paper59_pslq_negative_disc4_guarded` — no low-height
  ≤40 closure of W=V·π/8, height ~140–390, decoy-matched, + positive control; fast); Paper 59
  §modular/§bessel_algebra prose softened so the weight-3/disc-8 legs read as **bounded/driver-observed**
  ("consistent"; over-determined at ~19 digits; decisive test needs ~40), with disc-4 weight-≤2 the
  decisive leg. C8 headline #7 + matrix row re-scoped to match. Paper compiles clean (9 pp). Seed key
  `debug/qa/paper_59_delta_seed_key.json`; worktree removed, no leak. **Path to cert:** a FRESH DELTA
  `/qa` (verifies the fallback test + softened prose) → clean-delta → FULL certifying `/qa` → PASS.
- 2026-08-17 — **DELTA-verification run #2 = CLEAN-DELTA (trustworthy).** Verifies the PI-approved fallback
  (disc-4 permanent test + softened weight-3/disc-8 prose). Affected dimensions = **code + claims** only
  (citation + synthesis untouched by the fallback diff). Panel FULLY CALIBRATED: **sensitivity 4/4**
  (code S-code-1 W=V·π/4 wrong-object + S-code-2 positive-control removed→vacuous negative; claims
  S-claims-1 "transcendentally independent"=crit-4 dress-up + S-claims-2 "settles the open problem / no
  finite closed form"=crit-2 collapse), **specificity 0/5** controls false-flagged. Every MATERIAL finding
  mapped exactly to a planted seed; **zero genuine non-seed material defects** in the real fallback corpus.
  Bonus: the code reviewer independently re-ran PSLQ and CONFIRMED the real disc-4 negative is genuine and
  non-vacuous (W height 389, decoy 254, positive control height 3 = caught). Deterministic C5/C10–C18 all
  GREEN (C10: 9 pp, clean 2nd pass, no undefined cites). Seed key `debug/qa/paper_59_delta2_seed_key.json`;
  worktree removed, no seed leaked (π/8 + positive-control assert intact; neither claims seed present).
  **Precondition met — the next PI-invoked `/qa paper 59` fires the FULL certifying run (all dimensions,
  whole-paper enumeration, + the deferred completeness-critic) → PASS.**
- 2026-08-17 — **FULL certifying run (post-clean-delta) = FAIL (trustworthy).** All FOUR gating dimensions
  exercised (code / claims / citation / C9-synthesis on the Part-D group3-foundations fold — the DoD scope
  §C9=N/A line is superseded by Part D + the second-cert PI direction). Panel FULLY CALIBRATED:
  **sensitivity 8/8** (code S-code-A Wronskian tautology + S-code-B L4-uses-`_sI`-not-`_sK`; claims
  S-claims-B abstract "no finite closed form exists" crit-2 [caught MATERIAL] + S-claims-A cusp-form
  [OBSERVATION]→[SYMBOLIC PROOF]/"We prove" [caught, graded NIT]; citation S-cite-A katz vol 124→127 +
  S-cite-B FSY vol/yr 17/2023→15/2021; synthesis S-synth-A "corroborates→lifts/filled" + S-synth-B
  "open→now delivered"), **specificity 0** false-positives on control *values*. All four reviewers read the
  worktree (each caught worktree-only seeds). Deterministic C5/C10–C18 GREEN (group2 + synthesis).
  **VERIFIED GENUINE MATERIAL (code/C1 — several C8 `[MEASURED, N-digit]` headlines under-witnessed by the
  permanent tests; same class the paper already FAILED on twice — collinear 19-vs-12, PSLQ):**
  (a) **C8 #1 momentum 1.9×10⁻¹⁴** — `test_momentum_reproduces_three_center_eri` runs grid (90,32,32),
  achieves **1.45e-12**, asserts **<1e-9**; the 1.9e-14 appears only at grid (120,48,48)→2.7e-14 / (160³)→2.1e-15,
  run by NO test (PM-verified directly). (b) **C8 #8 quadratic relations "25 digits"** —
  `test_paper59_period_pairing_minus_pi_0_2pi` asserts **<1e-10** at dps=30 (achieves ~1e-17); 25 digits needs
  dps≥50, unreached (PM-verified). (c) supporting/same-class: **C8 #7 λ=1−ρ "10⁻⁴¹"** test asserts 1e-20 /
  achieves ~1e-31; in-module "45 digits" unreached at dps=30 (this one A7-dispositioned as matrix-note NIT).
  **CLEAN dimensions:** claims (real abstract correctly says "no closed form via factorization … we leave
  open"), citation (real katz=**124**/1990, FSY=**17**/541–602/**2023** both correct; 35/37 clean), synthesis
  (real fold says "open" + "corroborates rather than lifts … not filled") — every MATERIAL there was a seed.
  Two-way UPGRADE surfaced: `test_modulus_pf_annihilates_periods` achieves ~5e-32 vs the quoted 10⁻¹⁸ tag —
  could carry a tighter digit count. NITs: citation "Broadhurst–Mellit" vs the FSY-native "Broadhurst–Roberts"
  quadratic-relations naming **[APPLIED 2026-09-07 by `/qa paper_61` — logged here as a
  NIT on 2026-08-21, still unapplied when a citation reviewer had to
  re-find it and it was swept across 9 loci in 4 documents. A NIT that names a real defect
  is debt, not noise]**; the two `test_paper59_*` files well-cited. Completeness-critic DEFERRED (FAIL).
  Seed key `debug/qa/paper_59_certrun2_seed_key.json`; worktree removed, all 8 seeds verified absent from the
  real corpus. **Path to cert:** remediate the backing under-witnessing — for each C8 `[MEASURED, N-digit]`
  headline, either witness N in a permanent test (run the momentum test at the finer grid asserting ~1e-14;
  raise the period-pairing test to dps≥50 asserting ~1e-25; raise λ dps) OR down-state the headline digit to
  the witnessed value (and reconcile the abstract/§/C8/matrix) — then DELTA-verification `/qa` → clean-delta →
  FULL certifying `/qa` → PASS. Qualitative architecture (irreducibility/monodromy/PF/Legendre-Γ(2)/PSLQ
  methodology) is genuinely, non-tautologically proven — only the numeric-headline backing precision is short.
- 2026-08-17 — **REMEDIATION of the FULL-cert FAIL DONE (PI-directed "remediate").** The C8
  `[MEASURED, N-digit]` headline under-witnessing closed — every gap either WITNESSED in a permanent test
  or the tag RECONCILED to the regression witness (+ driver figure disclosed):
  (a) **C8 #1 momentum 1.9e-14** — added a converged-grid witness to `test_momentum_reproduces_three_center_eri`
  (grid (128,48,48), residual **7.8e-15**, asserts **<1.9e-14**; the coarse <1e-9 method check kept), ~1.7 s.
  (b) **C8 #8 quadratic relations "25 digits"** — `test_paper59_period_pairing_minus_pi_0_2pi` raised dps 30→**55**,
  assertions 1e-10→**1e-25** (worst residual ~2.3e-29). (c) **C8 #7 λ=1−ρ "10⁻⁴¹"** — `test_cosmic_galois_family_is_gamma2`
  raised dps 30→**50**, assertion 1e-20→**1e-41** (worst ~2.7e-51). (d) **C8 #4 in-module "45 digits"** — witnessing
  45 digits is genuinely expensive and the load-bearing claim is the saturation *discriminator*, so the paper
  (eq:inmodule) tag was RECONCILED `[MEASURED, 45 digits]`→`[MEASURED]` with prose "residual saturating flat …
  (below 10⁻²⁵ in the regression witness, 5×10⁻⁴⁴ at driver precision)". (e) **two-way UPGRADE** — eq:modpf tag
  tightened `10⁻¹⁸`→`10⁻²⁵` (underclaim; test witnesses <1e-25). **Collinear ~19 digits LEFT unchanged** — the
  frozen DoD (W3/C8 #7) canonicalizes the value to ~19 digits with the ~16-digit anchor CONFIRMED (cert-1
  dispositioned; deg5 ~16-digit re-derivation is >120 s, too slow to pin; test witnesses ~12 at deg4). `claim_test_matrix`
  rows 405/406/408/411 reconciled to the new witnesses. **Full P59 backing suite = 127 passed** (`--slow`:
  test_routeC_momentum + test_paper59_{corner_sigma2,bessel_moment_algebra} + test_two_center_eri_aabb); no dps
  contamination; paper compiles clean (9 pp, 0 undefined). No headline VALUE changed (so no C17/C16 registry edit) —
  only backing precision raised / tags reconciled. **Ready for the PI-invoked DELTA-verification `/qa`** (diff-scoped:
  the momentum fine-grid witness, the dps-raises, the two prose reconciliations) → clean-delta → FULL certifying `/qa`
  → PASS.
- 2026-08-17 — **FULL certifying run #3 (PI skipped the delta, direct FULL cert) = FAIL (trustworthy, thin) → MATERIAL remediated.**
  All four gating dimensions exercised (code/claims/citation/C9-synthesis), whole-paper enumeration. Panel FULLY
  CALIBRATED: **sensitivity 8/8** (code: λ-assert loosened 1e-41→1e-15 + Wronskian tautology; claims: §scope
  "no closed form" crit-2 collapse + "completing the theorem" crit-5; citation: adamsbognerweinzierl wrong-ID
  vol/yr/arXiv + brown_cosmic vol; synthesis: "resolves/filled" + "candidate→proven"), **specificity 0** FP on
  control values (every remediated headline — momentum <1.9e-14, period-pairing dps=55 <1e-25, in-module reconciled —
  marked SOUND). Deterministic C5/C10–C18 GREEN. **Genuine non-seed findings:**
  (D) **[MATERIAL, thin — code/C8 #8]** the paper's "the *three* regular masters satisfy eq:pf to ~10⁻¹⁷" was
  directly annihilation-tested for **s_K only**; claim VERIFIED TRUE (PM-checked: s_K 1.1e-18, s_I 3.4e-17,
  s_J 9.8e-16) but under-covered, and "~10⁻¹⁷" is best-case (s_J is ~10⁻¹⁵). **REMEDIATED same day:** new
  `test_paper59_all_three_regular_masters_satisfy_pf` (loops s_K/s_I/s_J through eq:pf, guarded <1e-14); paper
  bound softened `~10⁻¹⁷`→`~10⁻¹⁵` with the per-master residuals stated; matrix row updated.
  (A) **[NIT — DoD-dispositioned]** collinear ~19-digit value witnessed to ~12 in tests: the frozen DoD **W3/C8 #7
  canonicalizes ~19 with the 16-digit anchor CONFIRMED** and matrix row 409 discloses the ~12–16 witness — a frozen
  disposition; the paper correctly matches the DoD, so re-flagging = goalpost-move → left unchanged.
  NITs (fix-on-sight, non-blocking): eq:besselmoment "30 digits" no direct test (inherited via K0+Laplace);
  two untagged interpretive sentences (add [OBSERVATION]); two-way UPGRADE — eq:K0 6e-18 tag under-states its
  <1e-35 backing. **Verification:** new test passes (~8 s); C17 GREEN (no headline value changed — the 1e-17→1e-15
  is a residual-precision softening, not a registered family); paper compiles clean (9 pp, 0 undefined). Seed key
  `debug/qa/paper_59_certrun3_seed_key.json`; worktree removed, all 8 seeds verified absent from the real corpus
  (correct katz 124/1990 + FSY 17/2023 + adamsbognerweinzierl 1405.5640/55/2014 + brown_cosmic vol 11 intact).
  **Path to cert:** DELTA-verification `/qa` (diff-scoped: the new three-master test + the ~10⁻¹⁵ softening) →
  clean-delta → FULL certifying `/qa` (completeness-critic) → PASS. Qualitative architecture proven throughout;
  the residual gap was one plural-claim coverage detail on a verified-true C8 element.
- 2026-08-17 — **FULL certifying run #4 (PI: "paper 59 full") = PASS ✅ — CERTIFIED.** All four gating
  dimensions exercised (code/claims/citation/C9-synthesis), whole-paper enumeration + completeness-critic.
  Panel FULLY CALIBRATED: **sensitivity 8/8** (code: momentum-witness loosened 1.9e-14→1.9e-9 + three-master
  loop `_sK` thrice; claims: §modular "proves … transcendentally independent" crit-4 + §scope "now delivers a
  polyatomic energy" crit-3; citation: bbbg2008 vol 41→43 + blochkerrvanhove2015 vol 151→153; synthesis:
  block-1 "…locus is reached" + block-2 "a lift … established"), **specificity 0** FP on control values. All
  reviewers read the worktree. Deterministic C5/C10–C18 GREEN (group2 + synthesis). Completeness-critic surfaced
  only seeds + logged NITs (no new MATERIAL). **Zero genuine non-seed MATERIAL:** the code reviewer's largest
  non-seed flag (L4 factor-battery/eigenring prose "over ℚ(ρ)(D)" / "robust across ρ" vs ρ=1/5 tests) is a
  **corroborating-evidence NIT** — the load-bearing irreducibility is PROVEN by the ρ-generic Fourier–Laplace
  argument (Katz) + the two-ρ integer monodromy (both adequately backed; DoD C8 #5 frames the battery as
  "Supporting"), and the completeness-critic independently confirmed the ρ=1/5 battery is by-design. The
  re-flagged collinear ~19-digit value is **DoD-dispositioned** (W3/C8 #7 canonicalizes ~19, anchor CONFIRMED,
  matrix 409 discloses ~12–16 witness). **Fix-on-sight NITs applied:** eigenring tag "robust across pole windows
  and ρ" → "at ρ=1/5; robust across pole windows"; period-pairing "constant across the tested D and ρ" →
  "at D=1 across the tested ρ, and D-independent by the self-adjoint Lagrange identity". **Logged NITs (non-blocking):**
  eq:besselmoment "30 digits" + L(0,ρ)=K(1−ρ) "19 digits" have no direct test (inherited via the K0/Laplace/period
  chain); §5 L537 [OBSERVATION] tag under-states its FL+monodromy proof (two-way upgrade candidate); module
  docstring stale. Paper compiles clean (9 pp, 0 undefined); no headline VALUE changed (no C16/C17 edit). Seed key
  `debug/qa/paper_59_certrun4_seed_key.json`; worktree removed, all 8 seeds verified absent (real: "open"/"unestablished",
  vol 41/151, 1.9e-14, s_K/s_I/s_J).
  **HONEST CEILING:** PASS = survived the calibrated detectors for the seeded defect classes + the pre-registered
  criteria (transcendence-tier + open-frontier honesty, W1–W8, C1–C18). NOT provably perfect. Un-exercised this arc:
  a >~19-digit re-derivation of the collinear value (the ~40-digit collaboration frontier) and the weight-3/disc-8
  PSLQ decisive test (provably out of reach at ~19 digits) — both DoD-acknowledged as beyond scope. Standing debt:
  Paper 59 edits (v4.84–88 + this cert arc) still sit on `work/sparsity-boundary`; merge-to-`main` + any Zenodo
  Release remain PI-only. The C9 surface in the group3-foundations synthesis is now certified for Paper-59 content
  (a future `/qa group3` inherits it).

- 2026-08-19 — **FULL re-cert of the modular arc (v4.82–4.99, PI-directed to clear the Phase-4 OWED)
  = FAIL, calibration-backed.** Deterministic C11/C13/C14/C16/C17/C18 all PASS whole-paper (incl.
  C17: the v4.99.0 sec:modular edit's earlier value-conflict was self-caught + reverted pre-review, so
  the canonical 0.3953557659017139641 is intact). Panel = code-reviewer(Sonnet) + claims-reviewer(Opus)
  + citation-reviewer(Sonnet), worktree-seeded (K=6: code×2, claims×2, citation×2; M=5 controls).
  **CALIBRATION:** claims 2/2 seeds + 0/5 FP (CLEAN); citations 2/2 seeds + 0 FP (CLEAN); code C1 caught
  (→NIT: reviewer correctly argued the loosened err-tol is backstopped by the exact integer-matrix
  assertion), **C2 un-reachable (seeding artifact — planted in a test masked by F1's own import
  failure)**; code reviewer competence otherwise demonstrated (found F1/F2, re-derived Q'(x)Q'(-x)
  identity, built its own PSLQ positive control). **VERIFIED MATERIAL (code/C1 — the FAIL): F1** —
  `test_intersection_form_pi_and_planes_are_symbolic` (cited for eq:planepi [SYMBOLIC] B²=-π²) and
  `test_f12_kernel_swap_and_native_geminal` (cited for sec:f12 [MEASURED] 1e-14 / 5a/8), plus
  `test_intersection_form_thimble_block_structure`, `importlib`-load **UNTRACKED** transient
  `debug/routeC_intersection_form.py` / `debug/fock_f12_momentum_probe.py` → do not run from a clean
  checkout, while the file docstring falsely claims "Self-contained (no debug/ import)". Underlying math
  re-derived correct; permanent backing absent. Pre-existing (v4.82–4.98 arc). **NITs:** F2 (the -π,0,2π
  inline cite points to `test_intersection_form_period_cut_values` ~1e-10, not the 25-digit witness
  `test_paper59_period_pairing_minus_pi_0_2pi`); one proven-irreducibility clause tagged [OBSERVATION]
  (upgrade→[SYMBOLIC]); duplicate `sabbah2013`/`sabbah_stokes2013` bibitem key; `fresansabbahyu2023`
  printed vol/year (correct is ANT 17, 2023; arXiv ID fine); trivial `harris_michels1967` pages; two-way
  tag upgrades eq:K0 (6e-18→<1e-35 achieved), eq:modpf (1e-25→~1e-32). **PATH TO CERT:** remediate F1 —
  port the 3 tests self-contained (inline branch_point_proof/thimble_form + the f12 J_momentum/kernels,
  matching the policy the v4.99.0 tests already follow) OR move the driver logic to a tracked module —
  and fix the docstring; sweep the NITs; then a DELTA-verification /qa → clean-delta → certifying FULL run.
  No headline VALUE changed (no C16/C17 registry edit). Seed key `debug/qa/paper_59_seed_key.json`.
- 2026-08-19 — **REMEDIATION of the FULL-cert FAIL DONE (PI-directed "do it").** F1 closed: the 3
  debug-dependent tests (`test_intersection_form_thimble_block_structure`,
  `test_intersection_form_pi_and_planes_are_symbolic` [cited for eq:planepi SYMBOLIC B²=-π²],
  `test_f12_kernel_swap_and_native_geminal` [cited for sec:f12 MEASURED 1e-14/5a/8]) ported
  SELF-CONTAINED — the `debug/routeC_intersection_form.py` (`thimble`/`concomitant`/`thimble_form`/
  `branch_point_proof`) and `debug/fock_f12_momentum_probe.py` (`D`/kernels/`J_momentum`/
  `J_direct_samecenter`) logic inlined as `_if_*`/`_f12_*` helpers; the file docstring's
  "Self-contained (no debug/ import)" is now TRUE. **Full `test_routeC_momentum.py` = 36/36 GREEN**
  (was 33/36 with 3 FileNotFoundError). NITs swept: F2 (the -π,0,2π backing sentence now also cites
  the 25-digit witness `test_paper59_period_pairing_minus_pi_0_2pi`); duplicate `sabbah2013` bibitem
  merged into `sabbah_stokes2013` (arXiv:0912.2762 preserved, lone cite repointed). DEFERRED (harmless
  two-way, standing disposition): the [OBSERVATION]→[SYMBOLIC] upgrade on the proven-factorization
  clause L622 (already correctly hedged), and the eq:K0 6e-18 / eq:modpf 1e-25 conservative tags
  (frozen headlines; correct-but-conservative). Deterministic C11/C13/C14/C16/C17 all PASS; paper
  compiles clean (11 pp). No headline VALUE changed (no C16/C17 registry edit). **PATH TO CERT:**
  DELTA-verification `/qa` on this diff (code = the 3 ports; citation = the bibitem merge) → clean-delta
  → certifying FULL run. Remediation functionally verified (36/36 green) but not yet adversarially
  delta-reviewed.
- 2026-08-20 — **DELTA-VERIFICATION of the remediation → CLEAN-DELTA.** Diff-scoped /qa on the
  remediation (worktree `../geovac-qa-delta-p59`, now removed). Two Sonnet reviewers, each seeded
  with 2 blind diff-planted controls (`debug/qa/paper_59_delta_seed_key.json`). **Both dimensions
  CALIBRATED:** code-reviewer caught D1 (loosened tolerance `cross < 10**1` in the ported thimble
  test) + D2 (tautology `Bsq - Bsq == 0` in the ported plane-π test); citation-reviewer caught D3
  (wrong LNM volume 2011→2060) + D4 (fabricated test name `test_paper59_quadratic_relations_25digits`).
  4/4 seeds caught, 0 false positives. **Genuine findings (non-seed), verified against primary
  code/corpus:** (1) NIT — my helper-block insertion displaced `@pytest.mark.slow` onto a helper
  (`_if_thimble`) instead of `test_intersection_form_thimble_block_structure`; **FIXED** (decorator
  re-attached to the test; `pytest --collect-only -m slow` now lists the test, not the helper; test
  still passes under `--slow`). (2) The citation-reviewer flagged the F2 inline cite
  `test_paper59_bessel_moment_algebra.py::test_paper59_period_pairing_minus_pi_0_2pi` as "file does
  not exist in the worktree." **Verified against the real corpus (mandatory §6 reconcile): the finding
  is a worktree-sync artifact, NOT a paper-claim defect** — the file is present in `tests/` (a
  *permanent* home per §9, distinct from the transient `debug/` the F1 defect was about), the function
  exists, the 25-digit witness PASSES on the real corpus (11.7s, worst residual ~2.3e-29 at dps=55),
  C13 `--gate group2` PASSES, and the DoD (this file, prior B2 record) sanctions the citation. The
  reviewer saw "missing" only because the delta worktree synced the modified *tracked* files but not
  the *untracked* `tests/` files. **Genuine repo-state residual (PI hand-off, not a paper defect):**
  two cited backing tests — `tests/test_paper59_bessel_moment_algebra.py` (25-digit relations,
  DoD C8 #8) and `tests/test_paper59_corner_sigma2.py` (σ² corner) — are **untracked** (`??`), on the
  same footing as the entire uncommitted v4.85–v4.99 arc (the paper's own sec:bessel_algebra content
  is likewise absent from HEAD). They MUST be `git add`ed when the arc is committed, else the
  DoD-frozen 25-digit + σ² backing dangles on a clean checkout. **CONSEQUENCE FOR THE CERTIFYING FULL
  RUN:** its seed worktree must include these untracked `tests/` files — either fire it *after* the
  arc is committed (worktree from HEAD then has them), or have the worktree-sync copy untracked
  `tests/` files too. Deterministic C11/C13/C14/C16/C17/C18 all PASS on the real corpus; no headline
  VALUE changed. **VERDICT: CLEAN-DELTA** (both dimensions calibrated, 0 verified paper-claim defects,
  1 NIT fixed). Clean-delta is the precondition met; the certifying FULL run remains a PI-timed step.
- 2026-08-20 — **sec:modular RECONCILIATION (v4.100.0), post-DELTA.** After the CLEAN-DELTA above,
  a PI-directed two-track push (a: outer-corner numerics; b: Γ(2) pullback) reconciled into
  sec:modular. NO headline changed (T2=0.3953557659017139641 C8 #7 frozen; 25-digit relations
  C8 #8 unchanged). Three surgical sharpenings, all backing-UPGRADES not reversals: (1) resurgence
  [OBSERVATION] now carries a measured Gevrey-1/Borel-2 discriminator; (2) ring-correction scoped —
  the wt-3-requiring-G element is the D→0 *shadow*, the period-only negative = ring-incompleteness
  (non-classicality comes from the fibre's irregular twist, not the PSLQ); (3) outer wall
  characterized (only (0,0) fractional ρ^{3/2}; oscillatory corners integer-leading; complex
  off-axis singularity; independent momentum evaluator confirms ~13–14 digits then plateaus).
  New tracked backing `tests/test_paper59_resurgence_corner.py` (2 tests, self-contained; corner
  orders + resurgence signature) — cited inline; `test_L4_irregular_at_infinity` (tracked) also
  cited. Paper compiles 12 pp; deterministic C10–C18 all PASS on the real corpus. **NEW untracked
  cited backing to commit with the arc:** `tests/test_paper59_resurgence_corner.py` (joins the
  already-flagged `test_paper59_bessel_moment_algebra.py` + `test_paper59_corner_sigma2.py`).
  Compounds the standing Phase-4 re-review OWED; the certifying FULL run remains PI-timed and must
  see these untracked tests/ files (commit-first or worktree-sync).
- 2026-08-20 — **BD-pullback front integrated (v4.100.0, third T2 track), post-reconciliation.**
  The τ-plane pullback's weight-2 Jacobian is now the explicit closed form
  dλ/dτ=iπθ₂⁴θ₄⁴/θ₃⁴=iπλ(1−λ)θ₃⁴ (classical modular-lambda derivative; verified ~1e-32
  finite-difference + symbolic), upgrading sec:modular's prior bare "weight-2 λ′ Jacobian"
  assertion; and the BD closing step is localized to ONE named step (the Eichler integral of that
  Jacobian against the D=1 Bessel twist N(D); a(n) = the twist's Fourier–Whittaker coefficients,
  not Eisenstein divisor sums). NO headline changed (T2 value + 25-digit relations frozen). New
  tracked backing `tests/test_paper59_bd_jacobian.py` (self-contained; Jacobian identity + leading
  Lambert coeff). Paper compiles 12 pp; deterministic C10/C13/C17 re-run PASS. **NEW untracked cited
  backing to commit with the arc:** `tests/test_paper59_bd_jacobian.py` (joins
  `test_paper59_resurgence_corner.py` + `test_paper59_bessel_moment_algebra.py` +
  `test_paper59_corner_sigma2.py`). Compounds Phase-4 OWED.
- 2026-08-20 — **Eichler/co-area front integrated (v4.100.0, fourth T2 track).** Exact length-2→
  length-1 reduction T2=(16/π)∫₀¹Φ(ρ)dρ (co-area; folded by Φ(ρ)=Φ(1/ρ)/ρ²), reproducing the
  anchor to 2.2e-9 (re-run confirmed; Jacobian verified by hand). Corrects the third-track twist
  naming: the genuine twist is the scale-integrated two-mass Φ(ρ), N(D) its ρ→0 slice; leading
  log coeff A=J(¼,¼,1)/4 (genus-0 diagonal); cusp series asymptotic (reliable terms overshoot
  68%→51%) ⇒ convergent closure excluded, wall = Borel–Lambert resummation. NO headline changed
  (T2 frozen). New tracked backing `tests/test_paper59_coarea_reduction.py` (symmetry + A + coarse
  reduction, self-contained). Paper compiles 12 pp; C10/C13/C17 re-run PASS. **NEW untracked cited
  backing to commit with the arc:** `tests/test_paper59_coarea_reduction.py`.
