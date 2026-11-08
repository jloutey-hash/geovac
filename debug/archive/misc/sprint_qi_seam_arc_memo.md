# Sprint memo — the ℚ(i) seam arc: KMS/Hodge circle, seam audit, literature, Euclidean kinematics
Date: 2026-08-21 | Branch: work/sparsity-boundary | Version: v4.102.0
Umbrella memo (per /sprint-close: one umbrella + per-track sub-memos). Sub-memos, NOT superseded:
- Track 2: `debug/sprint_qi_seam_audit_memo.md` (adversarial seam audit)
- Track 3: `debug/sprint_kms_cm_literature_memo.md` (literature verification)
- Track 4: `debug/sprint_t2_euclidean_kinematics_memo.md` (kinematic dictionary)
This memo carries: Track 1 (main-session), the mirror-spin probe, the capture log, honest scope.

## 0. Arc origin
PI question after the v4.101.0 T2 close-out: "we had some stuff about time in the spectral
triple — is there maybe a connection?" — i.e. does the T2/Γ(2)/CM structure (Paper 59) connect
to the temporal/KMS structure (WH7, Paper 35, P56 Kramers)? Run as 1 main-session track + 3
parallel opus agents (PI directed parallelization).

## 1. Track 1 (main session) — the KMS circle IS the Hodge circle (12/12 exact)
Driver `debug/_kms_mt_torus.py`, promoted to `tests/test_paper56_kms_hodge_circle.py` (5/5).
On the j=1/2 doublet, with K = diag(2m_j) (the compact BW modular generator, WH7) and J the
Kramers/Hodge complex structure (P56 prop:hodge_cm_point):
- U e^{itK} U† = cos t·I + sin t·J (the flow, in the Kramers frame, is the Hodge circle);
- J = e^{i(π/2)K} transported — the flow's β/4 quarter-period point;
- quarter-periods = {I, J, −I, −J} ≅ μ₄ ⇒ conductor 4 = order of the flow's torsion point;
- e^{iπK} = −I on the spinor doublet vs +I on integer-l scalars (spin double cover visible
  inside the thermal circle at β/2; the scalar/spinor discriminant lives in the FLOW).
Deflations (stated in the paper remark): tori-conjugacy in SU(2) is automatic; both halves are
separately classical; the 2π-closure re-reads the already-registered compactness fact.
Terminology fix en route (Track-3 catch, verified): the computed group is the HODGE group
(special MT, MT∩SL₂), not the full Mumford–Tate group (which adds weight scalars,
Res_{ℚ(i)/ℚ}𝔾_m) — P56 prop renamed accordingly.

## 2. Track 2 verdict — CONVERGENT, not shared mechanism (sub-memo §1–§6)
Same arithmetic invariant (a ℚ-rational J², i.e. ℚ(i)/μ₄ data) on two unrelated ℚ²-carriers in
two different Tannakian categories; no comparison map. Route (a) Dirac: SPIN-forced (half-integer
shift × parity → quarter-integer Hurwitz → χ₋₄), mixed-Tate over ℚ(i). Route (b) T2: LEVEL-forced
— the audit's net-positive find: **θ₃² is the theta series of ℤ[i]** (Jacobi r₂(n)=4Σχ₋₄(d);
weight-1 level-4 Eisenstein series, L=4ζβ), so conductor-4 sits in the fibre period at EVERY
fibre; the Bessel twist makes √ρ=θ₄²/θ₃² physical ⇒ level Γ(2)→Γ(4)=ℚ(μ₄)=ℚ(i). The τ=i-fibre
corroboration leg is near-vacuous (any full-modulus family hits every CM fibre; physical fibre is
the ρ=1 cusp). Candidate mechanisms: metaplectic/spin = shared 4=2×2 PATTERN not cause (but
anchored: theta characteristics ARE spin structures, Atiyah 1971); KMS/μ₄ = real WH7↔P56 seam,
route-(a)-only (T2 has no thermal circle); real-forms/"Wick = base change to ℚ(i)" = PUN,
falsified by the compact-but-pure-Tate scalar sector. Upgrade targets: T-1 (natural J-equivariant
comparison map carrying QJ=I to the elliptic Riemann form) or T-2 (measured β(2) in T2 at ≥32
digits — blocked by the same second-cusp resummation as the closed form ⇒ **the seam question
and the T2 closed-form question coincide**).

## 3. Track 3 — literature (sub-memo, per-item verdict table)
NEW (4-framing search-negative): the flow-level BW-circle=Hodge-circle identification (with the
noted caveats); and the conductor-4 GAP — Catalan absent from the core Bessel-moment corpus
(Broadhurst 1604.03057, FSY 2006.02702, Zhou 1706.08308, BBBG 0801.0891 sit at 3/6/8/15), nearest
χ₋₄ = BD 2607.14020 §5.1 topological-string sector. NOT new (cite): CMR math/0501424 Thm 5.1
(KMS_∞ states LABELLED by CM points — not a flow identification); Chowla–Selberg/Gross–Deligne;
Hodge group of CM curve = norm-1 torus (classical); Angius–Volpato 2605.30418 (closest shape:
U(1)×U(1) R-symmetry flow realising a Hodge decomposition). Defects: connes_marcolli2004 usage in
P56 VERIFIED CORRECT (renormalization cosmic-Galois, not KMS); zhou_wick2018 must not be cited
for arithmetic Wick (never was); the MT/Hodge-group naming (fixed, §1).

## 4. Track 4 — Euclidean kinematics (sub-memo; dictionary validated 1e-30/1e-31)
The T2 fibre IS a 2D-Euclidean two-propagator correlator: masses aᵢ=ζ/√cᵢ (Feynman parameters =
Källén–Lehmann spectral masses, ρ = squared mass ratio), Euclidean times pᵢ=Dᵢ√cᵢ, twist
aᵢpᵢ=ζDᵢ family-invariant (D=1 = one Compton wavelength), j₀ = box smearing of the shared
spatial coordinate (third centre enters only through space). Stokes location = squared
complexified Euclidean interval |z*|=p₁²+|W|² (endpoint pinch; D≠1 rows discriminate, 35% at
D=0.5); no median ambiguity = Euclidean positivity; (0,0) corner = UV coincidence point where
the interval closes (explains the family-integration obstruction); NEW closed form for the three
oscillatory corners J=σ⁴(π/2)(7/e)²α(1−α)/|W|. Two-axis partition: mass split → elliptic,
time displacement → irregular (clean at on/off level, triangular at Stokes-data level);
corroborates WH7/Paper 35 (non-compact temporal displacement injects exponentials, not π);
signature honesty: the temporal reading is structural, not metric.

## 5. Mirror-spin probe (PI picture; main session, exact)
PI picture: pulling the spheroid's foci apart breaks the Paper-0 nodal pattern into two
mirror-reflection pieces "like spin up and spin down." Adjudicated: CORRECT with a precise home —
Paper 8's own J^(±)=½(L±A) (SO(4)=SU(2)×SU(2), both j=(n−1)/2). Verified exact on n=2 ((½,½)
rep): (i) parity/SWAP conjugates J^(+)↔J^(−) (the two factors are mirror twins); (ii) the
two-centre axis acts through A_z = J_z^(+)−J_z^(−) ("up minus down"), which MIXES the spherical
labels (A_z|2s⟩=|2p₀⟩ exactly) and is DIAGONAL in the two-spin weight basis (eigs m₊−m₋ =
{1,0,0,−1} — the Stark/parabolic ladder); (iii) spatial parity = −SWAP on n=2 (matches Paper 8's
Runge–Lenz-parity phase). Caveat: orbital pseudo-spins, NOT electron spin; the honest continuum
reading is the two Euclidean chiralities of SO(4) (Wick-rotated Lorentz), parity-swapped. Even-n
shells carry genuinely half-integer twins. NOT captured in a paper yet — Paper 8 remark drafted
in concept, PI-gated (Paper 8 is certified group2). Driver: inline session check (reproduced in
the proposed test if the remark is approved).

## 6. Honest scope
**Theorem grade:** none new. (All classical-algebra pieces are individually textbook; the exact
checks are identifications, not new theorems.)
**Exact / bit-exact:** Track 1's 12 checks (sympy, promoted to 5 pytest functions); the n=2
mirror-spin checks (6 sympy checks); the theta-series identities (4 pytest functions, machine-
exact / ≥30 digits); the eq:period c_max prefactor (verified both orderings, 1e-18).
**Measured:** the Euclidean dictionary (1e-30/1e-31/1e-18/1e-10 legs); z* interval form (0.6–5%
at 8 points); oscillatory-corner closed form (1e-4, σ² approach); the PSLQ/convention table.
**Adjudications (evidence-graded, not theorems):** seam = CONVERGENT; τ=i leg near-vacuous;
level-forcing on route (b); novelty = search-negative (4 framings), not proof of absence.
**Corrections applied to certified/in-flight papers:** P59 G-in-ring justification was in the
wrong (modulus) convention — in the paper's own parameter convention both cited integrals are
rational (∫K dm=2, ∫E dm=4/3); replaced with the convention-free theta-series reason. P59
eq:period prefactor fixed to 1/√c_max (printed form wrong by ×2.3 when c₁=c_max). Two
"coincident scales" mislabels fixed (ρ→0 and c₂→0 are the ONE-MASS limit; coincident is ρ=1).
P56 prop renamed Mumford–Tate → Hodge group (conclusion unaffected).
**Named open follow-ons:**
1. Seam upgrade tests T-1 (comparison map) / T-2 (β(2) in T2 at ≥32 digits, = the second-cusp
   resummation frontier). The seam is not independently testable until T2 is.
2. Paper 8 mirror-spin remark (PI-gated; drafted in concept, test sketched).
3. WH7 falsifier program unchanged; the μ₄/Hodge-circle re-reading strengthens the compactness
   leg's interpretation, adds no new falsifier.
4. Paper 34 candidate row (from Track 4): non-compact Euclidean-time displacement → exponential/
   irregular period (distinct from the compactification 2π entry) — not yet added to Paper 34.
**Dead end recorded (§3 row):** ℚ(i) seam as shared mechanism — all four candidate common causes
fail (pattern/one-sided/route-(a)-only/pun); real-forms "Wick = base change" falsified by the
compact-but-pure-Tate scalar sector.

## 7. Capture log (what went into papers, all compiled clean)
- **Paper 59** (12→14 pp over the arc): theta/χ₋₄ ring-justification fix + Γ(4) level promotion
  (sec:modular); kinematics paragraph (sec:reduction); light-cone/two-axis paragraphs
  (sec:modular); conductor-4-gap placement datum (sec:placement); eq:period prefactor fix; two
  one-mass-mislabel fixes; loutey_paper35 bibitem. Tests: test_paper59_theta_chi4.py (4),
  test_paper59_euclidean_dictionary.py (5).
- **Paper 56** (26→27 pp): Hodge-group rename (prop + proof line); rem:paper59_cm rewritten
  (convergent-not-common-cause; theta-series leg load-bearing; τ=i demoted; upgrade targets);
  new rem:kms_hodge_circle (flow=Hodge circle, μ₄, spinor sign; three qualifications; CMR +
  Angius–Volpato defensive cites); bibitems atiyah1971, connes_marcolli_ramachandran2005,
  angius_volpato2026. Tests: test_paper56_kms_hodge_circle.py (5).
- Both papers compound the standing Phase-4 re-review OWED.

## 8. Files
- Drivers: debug/_kms_mt_torus.py; agent drivers debug/t2_euclidean_dictionary.py,
  debug/_t2_dict_secC{,3}.py (+ data logs debug/data/_t2_dict_*.out).
- Tests (tracked, all green): tests/test_paper56_kms_hodge_circle.py,
  tests/test_paper59_theta_chi4.py, tests/test_paper59_euclidean_dictionary.py.
- Memos: this file + the three track sub-memos listed at top.
- Nothing committed; PI controls commit/release.
