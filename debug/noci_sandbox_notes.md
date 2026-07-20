# NOCI sandbox notes (branch sandbox/noci — frozen-repo exploration, NO release path)

**Status:** sandbox only. Repo frozen at v4.76.0; push disabled (pushurl + pre-push hook).
No CHANGELOG, no version bump, no paper edits. This file is the preservation record.

**Thread:** the named open follow-on (a) of `debug/sprint_commutator_and_explorer_memo.md`
— "NOCI/VB-on-GeoVac feasibility (the only literature opening — relocates metric cost to a
device-measured many-body overlap S_IJ; scope before any build)." Prior art already in the
memo: Huggins et al. NJP 22, 073009 (2020); Baek et al. PRX Quantum 4, 030307 (2023).

## Step 1 — H2 machinery validation (`noci_h2_probe.py`, 2026-07-19)

STO-6G-fitted 1s orbitals (exactly renormalized), s-Gaussian closed-form integrals,
singlet spatial-pair NOCI. Anchor hit: covalent-only (Heitler-London) D_e = 3.158 eV @
R_eq = 1.645 a0 vs literature HL-STO zeta=1 ~3.14 eV @ ~1.64 a0. Span identity
(non-orth 3-config == Loewdin-orbital 3-config, same span) ~1e-15.
**Headline: 1 non-orthogonal covalent config binds; the SAME config on
Loewdin-orthogonalized orbitals is completely unbound.** MO det binds but fails
dissociation (classic); NOCI dissociates exactly. cond(S_config) <= 449.
Data: `debug/data/noci_h2_probe_results.json`.

## Step 2 — LiH, 4 electrons, genuine integrals (`noci_lih_probe.py`, 2026-07-19)

Basis {Li1s (6G fit, z=2.69), Li2s (STO-3G 2s shape, z=0.65), H1s (6G, z=1.0)};
exponents variational ON ATOMS ONLY (molecule is a prediction). General-N non-orthogonal
determinant machinery by explicit permutation expansion; S_IJ cross-checked vs det(M)
(worst 2e-16); complete-space identity FCI == Loewdin-FCI (2e-14); cond <= 15.6.

| dets | D_e | % of in-basis FCI binding |
|---|---|---|
| cov (2) | 0.857 eV | 64.3% |
| cov + Li+H- ionic (3) | 1.302 eV | **97.6%** |
| all 4 | 1.322 eV | 99.1% |
| FCI (15) | 1.334 eV @ R_eq 3.257 a0 | 100% |
| Loewdin cov (2) | unbound | -49% |
| Loewdin cov+ionH (3) | 0.200 eV | 15.0% |

NOCI-3 R_eq = 3.287 vs in-basis FCI 3.257 (0.9% truncation error on geometry).
Absolute: underbinds experiment (1.33 vs 2.51 eV) — expected, no Li 2p, no diffuse H.
Data: `debug/data/noci_lih_probe_results.json`.

## Reading

Non-orthogonality at the STATE level carries the chemistry: 3 fragment-native determinants
~= FCI; the same 3 after orbital orthogonalization keep 15% (2-det case: destroyed).
Compactness is a property of the representation, and qubit products only ever use truncated
spaces — this is the state-side face of the composition wall (Loewdin = forced commutation).

**Measurement ledger (NOQE-style):** n_det=3 → 6 S_IJ + 6 H_IJ per PES point; H stays at
native sparsity (corpus: composed LiH 334 Pauli @ Q30; Loewdin retrofit = 17.9x inflation,
v4.73.1). Trade in one line: ~12 overlap-type measurements per point buys escape from an
18x Pauli blowup AND restores binding the orthogonal truncation cannot express.

## Step 3a — surrogate-consistency audit (`noci_consistency_audit.py`, 2026-07-19)

Degrade the genuine LiH (S,h,g) GeoVac-builder-style, layer by layer; rerun ladder.
Collapse gate = E_min below genuine-FCI min (the -54 Ha Löwdin-retrofit failure mode).

| run | S | h | g | outcome |
|---|---|---|---|---|
| A | genuine | genuine | genuine | binds correctly (1.30/1.32 eV @ 3.26-3.29 a0) |
| B | genuine | cross-center off-diag zeroed (W1d mimic) | genuine | **catastrophic collapse** (E_min -11.3/-12.8 Ha, hole at small R) |
| C | genuine | genuine | cross-center overlap-density ERIs zeroed | **collapse** (-8.47/-8.56, spurious 15-18 eV attraction) |
| D | genuine | surr | surr | **collapse** (-8.75/-8.94) |
| E | **identity** | surr | surr | mild monotone over-attraction into small R, no interior min — **reproduces the corpus's balanced/composed failure signature** (NaH-style monotone descent) |

**Verdict: the hybrid shortcut is dead, in the worst possible way — genuine S bolted onto
surrogate h/g is WORSE than the orthonormal pretense (B/C/D collapse; E only drifts).**
Metric-operator mismatch is the failure mechanism (v4.73.1 confirmed in a controlled lab
where the exact answer is known). Consistency must be end-to-end: (S,h,g) from the same
genuine orbitals. Combined with the arc-closure line of the commutator memo ("the
l-sparsity lives in the BARE integral tensor; the only way to keep it is to NOT transform
the integrals"), the coherent roadmap is: genuine bare fragment integrals (SW / Avery
closed forms — sparse by Gaunt selection, untransformed) + NOCI at the state level +
device-measured S_IJ/H_IJ. No cheap bridge exists; the bridge IS the Avery machinery.
Caveat: builder-style surrogates are mimics on a textbook basis (structure-level, not
builder-reproducing). E's "R_eq=2.0" rows are grid-edge, i.e. no interior minimum.

## Step 3b — genuine-integral sparsity census (`noci_n3b_census.py`, 2026-07-19)

**Inventory (roadmap question answered):** `geovac/shibuya_wulfman.py` exposes
cross-center V_ne ONLY (same-center bra-ket pairs, off-center nucleus; plus
mismatched-lambda and multi-zeta variants) — no two-center overlap, kinetic, or ERIs.
Topos-3's exact Mulliken/Ruedenberg A_n/B_n machinery SURVIVES
(`debug/compute_topos3_exact_meet.py` + `tests/test_topos3_exact_meet.py`): genuine
two-center overlaps at arbitrary (n,l,m,Z) as exact rationals, zero-decidable
(Lindemann). `geovac/neumann_vee.py` is m=0-only Neumann V_ee (sigma states, the
spheroidal solvers) — not a general engine. **No two-center ERI engine exists anywhere
in the corpus.**

**New machinery (exact, validated):** extended the Topos-3 route with 1/r_A and 1/r_B
kernels (each cancels against the volume element (xi^2-eta^2); cancellation asserted)
→ genuine cross-center h via the hydrogenic eigen-trick, with the bra/ket assembly
identity (E_a S − Z_B I_B = E_b S − Z_A I_A) holding EXACTLY (rational equality) on
every entry — a machinery certificate. Four quadrature cross-checks at ~2e-16.

**Census (A: Z=3 at origin, B: Z=1 at R=3 zhat, hydrogenic n_max=2 per center, M=10):**

| object | genuine nonzero | builder nonzero | dense | genuine zeros are |
|---|---|---|---|---|
| S | 32 | 10 (identity) | 100 | m-rule only |
| h | 44 | 22 (no cross block = W1d) | 100 | m-rule only |
| g (ordered) | 2,944 (29.4%) | 214 (2.1%) | 10,000 | m-rule + Gaunt |

g at n_max=3 (M=28): genuine 114,280 (18.6%) vs builder 7,600 (1.2%) — **inflation
13.8x → 15.0x, rising with basis.** The genuine tensor is essentially m-rule-only
(per-side Gaunt feasibility cuts < 6% beyond the m-rule); the cross classes
(AA|AB + AB|AB + AB|BB) are 80% of it. Cross-S magnitudes 0.08–0.52; cross-h
magnitudes 0.19–0.80 Ha — **bond-scale, not perturbative** (an independent
quantification of W1d/W1e: the builder's architecturally-absent cross-block h entries
are load-bearing). Data: `debug/data/noci_n3b_census_results.json`.

**GATE VERDICT: STOP branch.** Gaunt l-selection does NOT survive on cross-center
terms — only the axial m-rule does (exact, structural; Topos-3's meet-is-m-grading,
now quantified on the full (S,h,g)). The genuine bare molecular tensor is l-dense in
every cross-center class. The roadmap's differentiator (i) — "untransformed bare-tensor
use keeps l-sparsity" — is FALSE at the molecular level: untransformedness protects
only the same-center (atomic) sector.

**The re-reading this buys (the census's real yield):** the corpus's 17.9x Löwdin
inflation (v4.73.1) and Track DF's 14x were never the price of orthogonalization —
they are the price of the genuine cross-center physics itself (~15x is the genuine
tensor's own inflation over the builder). "Löwdin destroys sparsity" and "the
molecular tensor is intrinsically l-dense cross-center" are the same fact seen twice;
Löwdin was paying the honest price. l-sparsity is a property of the atomic sector
only; no molecular representation keeps it while binding, because binding lives in
the cross-center classes.

What survives: N2's compactness (3 dets = 97.6% FCI binding) untouched; the m-rule is
exact (factor 3–5 at these bases); the (AA|BB) Coulomb class is low-rank via bipolar
multipole (SW-style); S/h are M^2 objects. A constant-factor story — not the corpus's
structural-scaling story.

**Thread disposition per the pre-registered gate: STOPPED honestly at N3b.** N4/N5/N6
not run as QC-product steps. Residual option (PI call only): N4-as-diagnosis — NaH
3-det NOCI with genuine integrals as a W1e closure-mechanism confirmation (would show
the wall is bindable cross-center physics, decoupled from any sparsity claim).

## Step 4 — N4-as-diagnosis: NaH BINDS (`noci_nah_probe.py`, 2026-07-19, PI-authorized)

PI authorized the diagnosis variant post-N3b ("Avery style, soup-to-nuts"); NO sparsity
claim. New machinery (`noci_md_engine.py`), all validated at machine precision before
use: McMurchie–Davidson s/p Gaussian engine (vs N1 closed forms 4e-16; p functions vs
FD center-derivatives ~1e-11); general-N Löwdin-cofactor non-orthogonal Slater–Condon
(vs N2 permutation machinery 4e-16; vs the stored N2 LiH ladder 2e-15 — N2's explicit
permutation expansion dies past ~6 electrons, this replaces it); bitstring FCI (vs
stored N2 FCI 2e-15); 6-Gaussian STO shape fits for 2p/3s (⟨fit|STO⟩ = 1.000000).

Setup: all-electron 12e NaH; Na {1s, 2s, 2p×3, 3s} + H {1s} (M=7 spatial); zetas
variational on the ISOLATED Na atom only (best = CR-like {10.63, 3.3, 3.44, 0.836};
E(Na) = −161.0845 Ha vs Clementi minimal-STO RHF anchor ≈ −161.12); H ζ=1 — the
molecular curve is a fragment prediction. Checks: rotational invariance (z vs 111 axis)
1.1e-13; non-orthogonal 91-det complete space = Löwdin bitstring FCI span identity
0.0e+00 bit-exact.

| method | R_eq (a0) | D_e | % of in-basis FCI binding |
|---|---|---|---|
| cov (2) | 3.954 | 0.686 eV | 58.4% |
| cov + Na⁺H⁻ ionic (3) | **3.736** | **1.071 eV** | **91.1%** |
| all 4 | 3.713 | 1.080 eV | 91.9% |
| FCI (91) | 3.595 | 1.175 eV | 100% |
| Löwdin cov (2) | — | unbound | −38.0% |
| Löwdin cov+ionH (3) | 3.658 | 0.136 eV | 11.6% |

**ALL FOUR PRE-REGISTERED GATES PASS.** G1 interior minimum ✓ (the corpus NEVER
produced one for NaH); G2 R_eq = 3.736 a0 = +4.8% vs experiment 3.566 ✓; G3 D_e =
1.071 eV sane vs experiment 1.961 ✓ (55% — minimal-basis underbinding, the same
pattern as N2's LiH at 53% of experiment); G4 compactness 91.1% ≥ 85% ✓. The in-basis
FCI R_eq = 3.595 a0 lands +0.8% from experiment.

**Diagnosis confirmed by construction.** The wall the corpus hit on NaH (W1e: monotone
descent through F1–F6, Schmidt, kwarg sweep, Sprint B.1 explicit-core HF, R3-B DMRG)
is the HAMILTONIAN SPECIFICATION — surrogate/missing cross-center physics — not the
correlation treatment: with genuine end-to-end integrals, the same 3-determinant
compactness that carried LiH binds NaH at ~5% geometry error. Both N2 readings
replicate at 12 electrons / second row / with p orbitals: non-orthogonal compactness
(91.1%) and orthogonalization-damages-only-truncated-spaces (Löwdin 3-det keeps 11.6%,
Löwdin 2-det destroyed — LiH was 15.0% / destroyed).

Avery-call relevance: natural centerpiece demo — the framework's hardest documented
failure, bound by the Averys' species of machinery (genuine two-center integrals over
Slater-type fragments) plus 3 fragment-native determinants.
Data: `debug/data/noci_nah_probe_results.json`.

## Honest gaps

Steps 1–2 and 4 use Gaussian-fitted STO integrals, NOT GeoVac-native integrals. Step 3a
closed the hybrid question (no metric patching, ever); Step 3b closed the sparsity
question (genuine cross-center tensor is l-dense; gate STOP). The g census is
structural (symmetry-rule counting in the complex-m basis; 'allowed' entries verified
generically nonzero only for S/h where the exact machinery decides — no numerical
two-center ERIs were computed). The real-harmonic builder count differs only by the
±m mixing. The builder comparator is the Gaunt-selected same-center superset
(convention-B-like); the production pair-diagonal A-convention is sparser still, so
the reported inflation factors are LOWER bounds.

## THE LONG PATH — roadmap to "does NOCI-on-GeoVac work?" (written 2026-07-19, pre-reset)

**The elimination argument (why this is the only road):** operator-side routes are all
closed by corpus theorems (Löwdin 17.9× Pauli, biorthogonal S⁻¹ just as dense, cross-block
h1 = 16× over-binding); hybrid metric-patching is closed by N3a (collapse); the bare
untransformed genuine tensor RETAINS l-sparsity (commutator memo arc-closure) and NOCI is
the one architecture that never transforms it. Every alternative is a documented dead end.

- **N3b — genuine-integral sparsity census.** [DONE 2026-07-19 — GATE = STOP; see
  Step 3b above. The genuine tensor densifies (l-dense cross-center, ~15x); the
  roadmap's sparsity payoff does not exist.] Inventory what
  `geovac/shibuya_wulfman.py` actually exposes (V_ne only? overlap? kinetic?) and whether
  Topos-3's Mulliken A_n/B_n closed-form code survives anywhere (grep debug/, tests/).
  Then: compute the GENUINE bare (S,h,g) tensor for a real GeoVac fragment pair
  (hydrogenic orbitals w/ l>0, two centers) and CENSUS its sparsity — how many nonzeros
  vs the surrogate tensor's? Gate: if the genuine bare tensor keeps O(same) sparsity
  (angular Gaunt selection surviving on cross-center terms), the roadmap holds; if it
  densifies, the NOCI route loses its sparsity payoff and the thread STOPS honestly.
- **N4 — NaH, the prize target.** [RUN 2026-07-19 as N4-as-diagnosis (PI-authorized,
  decoupled from sparsity claims) — **ALL GATES PASS**, see Step 4. The p-orbital
  engine blocker was closed by `noci_md_engine.py`.] Needs p
  orbitals on Na → either an s+p Gaussian engine (Obara–Saika recursion, a real build)
  or reuse of N3b genuine Sturmian integrals. Run the same 3-det ladder (covalent Na–H,
  ionic Na⁺H⁻). Gate: interior minimum, R_eq near exp ≈ 3.57 a0, D_e sane (exp ≈
  1.9–2.0 eV — VERIFY against a standard source before gating). Success = compact NOCI
  binds the framework's hardest documented failure with genuine integrals.
- **N5 — device-side NOQE resource table.** [NOT RUN — N3b gate STOP.] State-prep for fragment determinants over the
  native qubit register (Thouless rotation from an orthonormal computational reference —
  needs only the M² orbital-coefficient object, NOT the M⁴ transformed ERI tensor: this
  distinction is the whole ballgame). Hadamard-test counts: n_det(n_det+1)/2 × (1 + native
  Pauli terms) per PES point; shot-noise amplification via κ(S_config) (measured ≈ 15 at
  LiH). Deliverable: LiH/NaH resource table vs the Löwdin-inflated alternative.
- **N6 — verdict + the Avery intersection.** [RESOLVED by N3b on the negative branch:
  the honest negative is written into Step 3b; the Avery intersection remains real but
  as THEIR framework's home ground (standard molecular Sturmian QC), not as a
  GeoVac-sparsity joint program.] If N3b–N5 hold: "GeoVac fragment states +
  genuine Avery integrals + NOQE = sparse quantum chemistry that binds" — the natural
  joint program (their integral machinery × our encoding/labels × device-measured glue).
  If any gate fails, write the honest negative into these notes and stop.

**Standing constraints on the path (do not relearn these):** never orthogonalize anywhere
(N1/N2 kill); never mix metrics (N3a collapse — genuine S over surrogate h/g is the WORST
option); the "replacement-not-modification" cliff (v4.73.1) — the differentiators vs
standard Sturmian QC must remain (i) untransformed bare-tensor use, (ii) device-side
overlap payment, (iii) label-generated angular structure; if a step needs transformed
integrals, we have fallen off the cliff. Papers 8–9 guardrail not triggered (multi-center,
per-center exponents), but reload it before any basis unification idea.

## Resume recipe (post context-reset)

1. `git switch sandbox/noci` (work is COMMITTED there, local-only; main stays at v4.76.0).
2. Read this file top to bottom; probes: `debug/noci_{h2,lih}_probe.py`,
   `debug/noci_consistency_audit.py` (run from debug/ cwd; they import each other);
   data: `debug/data/noci_*_results.json`.
3. Freeze lockdown: pushurl = PUSH-DISABLED... + `.git/hooks/pre-push`; undo only on PI
   direction (`git config --unset remote.origin.pushurl` + delete hook). Sandbox rules:
   no pushes, no version bumps, no CHANGELOG, no paper edits.
4. Thread state (2026-07-19 end of day): N3b gate = STOP (sparsity payoff dead);
   N4-as-diagnosis RUN, ALL GATES PASS (NaH binds, +4.8% R_eq, 91.1% compactness —
   W1e localization confirmed by construction). N5/N6 remain not-run (QC-product
   steps, moot under the N3b STOP). Nothing queued. Open: Avery call prep owed when
   a date lands — N4 is the centerpiece demo. Drivers: `debug/noci_n3b_census.py`,
   `debug/noci_md_engine.py` (+ its validation suite), `debug/noci_nah_probe.py`.
5. Avery thread (separate but converging): availability email sent 2026-07-19; when a
   call date lands, build the one-page call prep. Primer for the PI:
   `debug/avery_framework_primer.md` (copied here from the temp scratchpad).
