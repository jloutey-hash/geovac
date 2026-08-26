# Period/transcendental-object census (base-rate denominator for the QI-SEAM)

Sprint: /aha Track 2b, 2026-08-21. Diagnostic only — no paper edits, no production-code changes.
Purpose: count every independently-constructed period/transcendental-identification object or
attempt in the corpus, classify by invariant class, and use the resulting population as the
base-rate denominator against which the P56-Dirac / P59-T2 ℚ(i)/conductor-4 recurrence
(`debug/sprint_qi_seam_audit_memo.md`) is judged ABOVE / AT / BELOW chance.

**Two populations are tracked separately** (collapsing them would misrepresent the finding):

- **Population A — blind numerical trials.** PSLQ / integer-relation / exhaustive-scan attempts
  that test a *given* numerical constant against a period basis, with no construction forcing
  the answer in advance. This is the natural "denominator" in the ordinary sense (Yukawa-PSLQ,
  W1e, S_min dossier, etc., exactly the examples named in the task brief).
- **Population B — forced/constructed mechanisms.** Theorems or proofs where a specific
  algebraic/geometric construction (vertex-parity projection, a Hodge complex structure, an
  elliptic modular family) *necessarily* produces period content of a determined class as a
  consequence of the construction, not a search. These are not "trials" in the base-rate sense —
  each is a narrow, purpose-built route — but they are the corpus's only source of non-pure-Tate
  content, so they matter for the "constrained landscape" reading (see verdict memo).

Class key: **pure-Tate** = ℚ[π,1/π] (M1) or ⊕π^{2k}·ℚ (M2); **disc-4** = ℚ(i)-CM-adjacent
(Γ(¼)-periods, Catalan G=β(2), β(4), χ₋₄ content); **disc-8** = ℚ(√-2)-CM-adjacent;
**MZV-L1/L2** = mixed-Tate, level ≤2, no χ₋₄ (log2, ζ(odd), classical Hoffman t-values);
**algebraic** = not a period at all (√(2k+1)-type); **rational** = ℚ, no transcendental content;
**unidentified** = attempted, no period-ring membership found at any tested precision/basis.

---

## Population A — blind numerical trials

| # | Object | Where defined | Class | Hit / attempt | N contrib. | Projection-route notes |
|---|--------|---------------|-------|----------------|:---:|-------------------------|
| A1 | SM Yukawa couplings, 9 fermions (e,μ,τ,u,d,s,c,b,t) × 3 transforms (yf, yf², ln yf) × 2 scales (M_Z, GUT) | `debug/sprint_yukawa_pslq_memo.md`; basis M1∪M2 pure-Tate, 9 generators, ceilings {10,100,1000} | — (external Class-1 data) | **CLEAN NEGATIVE**, 162/162 cells, 0 hits | 9 | Blind PSLQ. Basis was pure-Tate only (M3 excluded a priori by η-trivialization theorem), so this campaign could not have hit disc-4 even in principle — but it independently confirms zero low-coefficient pure-Tate structure either. |
| A2 | W1e NaH chemistry-correction terms, 11 distinct quantities (F4 PK barrier, F4/F6 FCI well depths, F5 Hartree J/K, experimental D_e, …) | `debug/sprint_w1e_period_class_memo.md` §1.1/§2.1; bases M1, M2, M3 (incl. disc-4 generators G,β(4),ζ(3),ζ(5)), INNER | — (external/algebraic-implicit) | **0/11** genuine outer-factor hits at audit ceiling 100 across ALL FOUR bases incl. M3; 3 M2 "hits" were rational-artifacts (filtered); 3 INNER "hits" were curve-fit-audit-failure-regime (filtered) | 11 | Blind PSLQ. This is the one campaign that explicitly tested against the disc-4 M3 basis (Catalan G, β(4), ζ(3), ζ(5)) and got **zero** hits — direct evidence against disc-4 being a promiscuous landing zone. |
| A3 | Δ = 1/40 mechanism-unification search (Paper 2 α ingredient) | `papers/group3_foundations/paper_18_exchange_constants.tex` L2721–2745, "Negative results on candidate mechanisms" | unidentified | 7 named mechanisms (Hopf-twist comparison, higher Casimir traces, Hopf quotient Laplacian, discrete/continuous S¹-fiber zeta, S⁵ spectral geometry, **3,228-candidate ζ-combination scan** over Riemann/Hurwitz zetas + Bernoulli numbers + Stieltjes constants) — **ALL NEGATIVE**, zero hits within 1e-25 of 1/40 | 1 (campaign; 3,228+ internal candidates) | Blind scan. Explicitly notes round-sphere spectral invariants "never" produce π² additively next to a rational — i.e. this campaign never even approaches disc-4 territory; it's testing pure-Tate/rational combinations only. |
| A4 | S_min (2-loop min-weighted CG sunset, Paper 28 eq:S_min) | `debug/smin_dossier_round1_memo.md`; old bases: 47-element (`smin_identification.py`) + 100-element (`smin_extended_pslq.py`), 15+1 "irreducibility" attempts, 150–200 dps | resolved: **mixed-Tate level 1/2** (ℚ[π²,ln2,ζ(3),ζ(5)]) | Old PSLQ runs: 16 clean FAILURES against those bases (later diagnosed as a **basis-coverage artifact** — the bases held only 2/8 of the weight-5 level-2 monomials). Resolved this round to an exact closed form via Hoffman-t-value reduction, residual 1.7e-129 | 1 | Landed at level-1/2 MZV — **not** disc-4. The old bases *did* include level-4 depth-1 Dirichlet-β monomials (i.e. disc-4-adjacent content) and still missed — S_min genuinely lives outside disc-4. |
| A5 | S_min^diff, depth-2 χ₋₄ analog (Sprint RH-P) | `papers/group1_operator_algebras/paper_29_ramanujan_hopf.tex` L899–918 | targeted disc-4 test | **35 PSLQ attempts across 7 basis strategies**, 100 dps, weight ≤8, specifically against χ₋₄ / β(s) / ζ(s) bases — **CLEAN NEGATIVE**, "numerically irreducible" | 1 (35 internal attempts) | This is the corpus's most **directly targeted** disc-4 attempt outside Paper 28/56/59 — a purpose-built χ₋₄ basis, aimed at exactly this class — and it missed. Strong evidence that landing disc-4 is not "easy" even when explicitly sought. |
| A6 | S^(3) (3-loop chain sunset, factorized form) | CHANGELOG 2026-06-11 (~L3098–3190); Paper 28/55 | mixed-Tate level 1/2 | 4 pending PSLQs **ACCEPT** against a complete level-2 weight-homogeneous basis (220–340 dps, residuals 1e-216..1e-339); no disc-4 basis element involved or needed | 1 | Landed classical level-2 MZV — disc-4 not reached, not tested for (basis had no χ₋₄ generator; the object's weight-graded structure ruled it out a priori). |
| A7 | Hain–Brown Sym² outer-factor period comparison (GeoVac Sym² periods vs the generic non-CM modular ring {E₄,E₆,ζ(3),ζ(5),Δ}) | `debug/sprint_hb_pslq_test_memo.md`, `debug/sprint_hb_eichler_kernel_memo.md`; Paper 56 §sec:open_g4_hodge, rem:cm_explains_hb | generic (non-CM) modular ring | **240-cell PSLQ panel** (20 periods × 4 bases × 3 precisions) — CLEAN NEGATIVE; repeated with a richer **22-generator Eichler-kernel-transformed basis**, second 240-cell panel — also CLEAN NEGATIVE | 1 (≈480 internal cells) | Blind search — but against the *wrong* ring. Paper 56 later explains the miss structurally: GeoVac's periods are ℚ(i)-CM-**pure** (disc-4), and disc-4-CM periods are degenerate/non-generic relative to the full Hain–Brown modular ring (which needs E₄,E₆,Δ). The miss is itself indirect evidence FOR the narrower disc-4 class, not evidence against period structure altogether. |
| A8 | Hydrogen 1S Bethe logarithm | CHANGELOG ~L10116 (Sprint TD-PSLQ-1) | unidentified | 27 PSLQ tests × 3 ceilings vs a frozen 64-form mechanical M1/M2/M3/ALG/depth-2 basis (**includes disc-4 M3 generators**) — clean NULL | 1 | Blind; basis explicitly included disc-4 content and still missed. |
| A9 | A_60(1S) spacetime-channel constant | CHANGELOG ~L10126 (Sprint TD-PSLQ-2) | unidentified | 30 PSLQ tests × 3 ceilings — clean NULL | 1 | Blind; same 64-form basis family, same miss. |
| A10 | Paper 59 integrated T2 observable (the collinear 3-center ERI value itself, as opposed to its fiber periods) | `papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex` L766–838; `debug/sprint_T2_*.md` | **UNDECIDED / open** | Guarded, decoy-controlled search: weight-≤2 disc-4 ring {1,π,K(½)} → **DECISIVE NEGATIVE**; weight-≤3 ring {π,K(½),1/K(½),G} (G=Catalan) → negative at reachable precision (~19 digits; would need ~32–40 to be conclusive) | 1 | The single live, still-open attempt at whether the **integrated observable** (not just its fiber) reaches disc-4. Currently a miss at reachable precision; independently characterized as a Γ(2) *resurgent Lambert series* (Gevrey-1, Borel radius 2), i.e. structurally NOT a finite classical period at all — the natural endpoint may be "never lands," not "landed and undetected." |
| A11 | Paper 59 diagonal cusp coefficient A = J(¼,¼,1)/4 | Paper 59 L913–921; `tests/test_paper59_diagonal_A.py` | unidentified (irreducible single-scale Bessel moment) | Guarded, decoy-controlled PSLQ vs weight-1 two-center ring {π, e^{-a}, K₀, K₁, E₁} at a∈{2,2√2} — **NO low-height closure** | 1 | Blind; no hit, disc-4 not reached at this sub-object either. |

**Population A total: N_blind = 30** independent objects/campaigns tested. **Confirmed disc-4 hits
in Population A: 0.** (A10 is the one still-open candidate and it currently leans negative at
reachable precision.) A2 and A5 are the two campaigns that explicitly built and tested a disc-4
basis and both returned clean negatives.

---

## Population B — forced/constructed mechanisms (not blind trials)

| # | Object | Where defined | Class | Route |
|---|--------|---------------|-------|-------|
| B1 | M1 (Hopf-base measure) | Paper 55 §M1; Paper 18 §taxonomy | pure-Tate ℚ[π,1/π] | Forced by the S¹→S³→S² Hopf-fibration measure Vol(S²)/4; witnessed across dozens of observables (state-space GH rate, gauge choice, temporal compactification, …). |
| B2 | M2 (Seeley–DeWitt heat kernel) | Paper 55 §M2; Sprint Mixed-Tate Test 2026-06-03 | pure-Tate ⊕π^{2k}·ℚ | Forced by the spectral-action heat trace on unit S³; √π cancels against (4π)^{3/2}, no ζ(3)/disc-4 content at any order. |
| B3 | M3 untwisted sub-sector (F-theorem F_s, F_D) | Paper 50 Thm 3.5; Paper 55 §"Canonical example 2" L1520–1565 | mixed-Tate **level 1** (log2, ζ(3)/π²) | M1×M3(untwisted) cross-product; uses the un-restricted MT(ℤ) sub-sector of M3, explicitly **not** the vertex-restricted disc-4 sub-sector. |
| B4 | **M3 vertex-parity Hurwitz identity** (Paper 28 Theorem 3 / RH-J.1) | Paper 28 "recent subsection"; cited widely (Papers 18, 29, 32, 55, 56) | **disc-4** (Catalan G=β(2), β(4)) | D_even(s) − D_odd(s) = 2^{s-1}(β(s) − β(s-2)). **Constructed, not searched** — the route is the vertex-parity projection γ_P=(−1)^n **explicitly applied** to the Dirac Dirichlet series. This is the corpus's origin point for disc-4 content. |
| B5 | M3 Galois-descent / period-map rank-2 | Paper 32 cor:m3_cyclotomic_mixed_tate; Paper 56 def:period_map, thm:injection_g4 | disc-4 (𝓜𝒯(ℤ[i,1/2]) level ≤4) | **Downstream of B4** — same vertex-parity mechanism, deepened via Deligne 2010 + Glanois 2015 import. NOT independently counted (same construction). |
| B6 | **Paper 56 Dirac Kramers-doublet CM point** | Paper 56 Prop. hodge_cm_point, L1673–1695; `tests/test_paper56_hodge_sl2.py` | **disc-4** (ℚ(i) exactly) | J²=−1 forces the Hodge/Mumford–Tate group of the fundamental rep to the norm-1 torus of ℚ(i). **Independent mechanism from B4** — algebraic complex-structure forcing, not vertex-parity Hurwitz. Bit-exact theorem. |
| B7 | **Paper 59 elliptic family CM-fiber sweep** (Rung 1+2, cosmic-Galois probe) | Paper 59 §sec:modular L707–736; `debug/routeC_cosmic_galois_rung{1,2}.py` | **disc-4 AND disc-8** | Rational modulus map λ(τ(ρ))=1−ρ (exact to 1e-41) makes the family the universal Legendre/Γ(2) family; physical domain ρ∈[0.36,2.78] happens to include τ=i (disc-4, ρ=½) and τ=i√2 (disc-8, ρ=2√2−2). Chowla–Selberg match ~1e-51. **Independent mechanism from B4 and B6** — level-forced via classical elliptic modular theory, not spin/vertex-parity, not Kramers-J. |
| B8 | Paper 59 native θ₃² χ₋₄ appearance | Paper 59 L811–828; `tests/test_paper59_theta_chi4.py` | disc-4 (level Γ(2)→Γ(4), field ℚ(μ₄)=ℚ(i)) | K(λ)=(π/2)θ₃² and θ₃² **is** the ℤ[i] theta series (Jacobi two-squares, r₂(n)=4Σχ₋₄(d)) at *every* fibre, not just CM points. **Downstream of / deepens B7** — same family, NOT independently counted. |
| B9 | Paper 38 GH-rate constant 4/π | Paper 55 §"Canonical example 3" L1566–1587 | pure-Tate M1 | Forced; = Vol(S²)/π²; no disc-4. |
| B10 | B = 42 (Casimir trace, Paper 2 α ingredient) | Paper 18 Thm thm:k_decomposition (B) | rational integer | Finite combinatorial SO(3) Casimir sum; not a period at all. |
| B11 | F = ζ_R(2) = π²/6 (Dirichlet identity, Paper 2 α ingredient) | Paper 18 eq:F_dirichlet | pure-Tate (M2-flavor) | D_{n²}(s=4) = ζ_R(2); forced Dirichlet-series identity at the packing exponent. |

**Population B, counting only independently-constructed mechanisms (collapsing B5→B4 and
B8→B7): 9 mechanisms** (B1, B2, B3, B4, B6, B7, B9, B10, B11). **Of these, 3 land disc-4-or-8**
(B4, B6, B7) = 3/9 ≈ 33% headline, but see verdict memo — the correct conditioning is on
*which* mechanisms are even capable of leaving pure-Tate territory at all (spin-doubling routes
and CM-elliptic-family routes), among which the disc-4 landing rate is 3/3 = 100%.

---

## Cross-check: Paper 34's 28-row projection dictionary

`papers/group6_precision_observations/paper_34_projection_taxonomy.tex`, Table
`tab:projection_axes`, L2886–2966. This is the corpus's own catalogue of every projection
(Layer-1 graph → Layer-2 observable) and its transcendental signature — a pre-existing,
independently-curated denominator.

Of the **28 rows**, exactly **1** (row 7, "Camporesi–Higuchi spinor," L2904–2906) carries
disc-4-adjacent content ("half-integer Hurwitz → odd ζ; Catalan G, β(4) at vertex") — and its
own column explicitly names the route as "at vertex," i.e. the same vertex-parity projection as
B4. Every other row is pure-Tate (π, π^{2k}, 2π — M1/M2 flavor, ≈8 rows), algebraic
(√(2k+1)-type Wigner-3j/6j content — not periods at all, ≈4 rows), ring-preserving over ℚ(α)
(≈6 rows), rational/trivial/"none" (≈7 rows), or explicitly flagged **PSLQ-disjoint from
M1∪M2∪M3** (the von Neumann-entropy apparatus-identity row, L2962–2965). **1/28 ≈ 3.6%** base
rate for disc-4 content across the framework's own master projection catalogue, and that one row
requires the identical vertex-parity construction as B4/B5.

---

## Summary counts

| Population | N | Confirmed disc-4/8 hits | Rate |
|---|---:|---:|---:|
| A (blind numerical trials) | 30 | 0 (1 open, leaning negative) | 0% (≤3.3% if A10 flips) |
| B (forced/constructed mechanisms) | 9 | 3 | 33% |
| B, conditioned on "mechanism engages spin-doubling or CM-elliptic structure" | 3 | 3 | 100% |
| Paper 34 projection dictionary (independent cross-check) | 28 | 1 | 3.6% |
