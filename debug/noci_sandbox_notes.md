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

## Honest gaps

These probes use Gaussian-fitted STO integrals, NOT GeoVac-native integrals. Step 3a
closed the hybrid question (no metric patching, ever); the remaining gap is genuine
fragment integrals at GeoVac scale — which is the Avery machinery's home ground.

## THE LONG PATH — roadmap to "does NOCI-on-GeoVac work?" (written 2026-07-19, pre-reset)

**The elimination argument (why this is the only road):** operator-side routes are all
closed by corpus theorems (Löwdin 17.9× Pauli, biorthogonal S⁻¹ just as dense, cross-block
h1 = 16× over-binding); hybrid metric-patching is closed by N3a (collapse); the bare
untransformed genuine tensor RETAINS l-sparsity (commutator memo arc-closure) and NOCI is
the one architecture that never transforms it. Every alternative is a documented dead end.

- **N3b — genuine-integral sparsity census (next, no new math needed).** Inventory what
  `geovac/shibuya_wulfman.py` actually exposes (V_ne only? overlap? kinetic?) and whether
  Topos-3's Mulliken A_n/B_n closed-form code survives anywhere (grep debug/, tests/).
  Then: compute the GENUINE bare (S,h,g) tensor for a real GeoVac fragment pair
  (hydrogenic orbitals w/ l>0, two centers) and CENSUS its sparsity — how many nonzeros
  vs the surrogate tensor's? Gate: if the genuine bare tensor keeps O(same) sparsity
  (angular Gaunt selection surviving on cross-center terms), the roadmap holds; if it
  densifies, the NOCI route loses its sparsity payoff and the thread STOPS honestly.
- **N4 — NaH, the prize target.** The system the framework could NOT bind. Needs p
  orbitals on Na → either an s+p Gaussian engine (Obara–Saika recursion, a real build)
  or reuse of N3b genuine Sturmian integrals. Run the same 3-det ladder (covalent Na–H,
  ionic Na⁺H⁻). Gate: interior minimum, R_eq near exp ≈ 3.57 a0, D_e sane (exp ≈
  1.9–2.0 eV — VERIFY against a standard source before gating). Success = compact NOCI
  binds the framework's hardest documented failure with genuine integrals.
- **N5 — device-side NOQE resource table.** State-prep for fragment determinants over the
  native qubit register (Thouless rotation from an orthonormal computational reference —
  needs only the M² orbital-coefficient object, NOT the M⁴ transformed ERI tensor: this
  distinction is the whole ballgame). Hadamard-test counts: n_det(n_det+1)/2 × (1 + native
  Pauli terms) per PES point; shot-noise amplification via κ(S_config) (measured ≈ 15 at
  LiH). Deliverable: LiH/NaH resource table vs the Löwdin-inflated alternative.
- **N6 — verdict + the Avery intersection.** If N3b–N5 hold: "GeoVac fragment states +
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
4. Continue at N3b (cheap, decisive) unless the PI says otherwise. NaH (N4) blocked on
   the p-orbital engine or N3b integrals.
5. Avery thread (separate but converging): availability email sent 2026-07-19; when a
   call date lands, build the one-page call prep. Primer for the PI:
   `debug/avery_framework_primer.md` (copied here from the temp scratchpad).
