# Pre-registered prediction: the hydrogen ring, two bases and two label sets (HR)

**Registered 2026-11-01** (PI: "I like your refinement of the question with the hydrogen ring", then "Sounds like
we have some tasks that we can continue with go ahead and do so"), before any energy of the scored system is
computed in the Slater basis, and before any label measure of the scored system.
**Commit status:** committed on its own, before any computation it scores.
**Rests on:** Paper 62 (water: Slater DZP below cc-pVDZ at equal size; matched-accuracy record MA, GO); record HS
(hydrogen fluoride, GO); Paper 25 (two groupings of one shell; Dawson-Dullin credit); the ring prior-art memos
`debug/lit_scan/prior_art_ring_molecules_qc_memo.md` and `prior_art_tori_quantum_numbers_memo.md`.
**Guardrail check (Papers 8-9, FCI-M, Track DF):** not in scope as registered. Every basis here is atom-centred, one
centre per hydrogen. The ring-centred one-centre Sturmian arm suggested in the QC memo IS in scope of Papers 8-9
(single-centre molecular encoding) and is left out; adding it needs the PI's explicit acknowledgement.

## Why

The PI proposed a toroidal carbon molecule and predicted: "toroidal maybe has multiple sets of quantum numbers that
can be used differently. Kinda like the multiple gauge groups earlier." A carbon torus is out of reach, so the
stand-in is a ring of six hydrogens. Two questions:
- **Part A.** Does Paper 62's equal-size result (exponential functions beat Gaussians at the same size) hold on a
  ring of pure hydrogen? Water and HF drew most of their margin from the heavy atom. Here there is none.
- **Part B.** A ring has two label sets for one set of orbitals: the site j and the ring momentum k. They are a
  Fourier transform of each other, so the exact energy is the same. The prediction is that they are useful for
  different things. Ring momentum gives exact zeros and extra qubit symmetries. Site labels give locality: terms
  that fall off fast and can be dropped.

## Method (fixed before scoring)

- **System:** H6, regular hexagon (D6h), nearest-neighbour distance 1.0 A (near the FCI minimum), xy plane, 6
  electrons, singlet. Driver `debug/hring_study.py`; rows `debug/data/hring_study.jsonl`; exponents
  `debug/data/hring_exponents.json`.
- **Solver:** the validated RHF + spin-orbital CCSD(T) of `tests/test_water_ccsd.py`, all electron, RHF reference,
  run through a memory-lean copy in the driver (`lean_run`). It uses the same equations but builds only the
  spin-orbital integral blocks it needs, and it never forms the virtual-four block (that term is contracted from
  the spatial integrals). Checked before registration: on H10 STO-6G, H6 cc-pVDZ and H6 Slater [2s] it gives the
  same RHF, CCSD and CCSD(T) as the original to every printed digit. Each block agrees exactly on random
  integrals. The rewritten term agrees to 1e-13 on a scale of 400, and a planted factor error in it shows as 53.
- **Slater DZP [2s1p] per H, Paper 62's H rule unchanged:** 1s, 1s' = 1.212, 0.627 (H2 RHF at R = 1.4); 2p exponent
  optimized on the ring's CCSD with an 8/6-Gaussian evaluator (bounded Brent on log p in [0.3, 4.0], xatol 1e-3, at
  most 60 iterations; mode `pol sto`), scored at 12 Gaussians per function. 30 functions.
- **cc-pVDZ, two forms:** as published, and with its H p exponent re-optimized the same way (mode `pol gto`; the
  fair comparison). Matched-accuracy ladder as in HS: cc-pVDZ+TZpol (cc-pVDZ s shells with cc-pVTZ's two H p
  shells, 48 functions, our construction), aug-cc-pVDZ (54), cc-pVTZ (84). Basis data: the files already used by
  Paper 62 (`debug/data/basis/`), hydrogen entries only.
- **Qubit cost (Part A):** Jordan-Wigner over RHF canonical orbitals, no frozen core, `jw_pauli_stats` (identity
  counted in N_Pauli, excluded from the 1-norm), tolerance 1e-10. Mode `final <basis>`.
- **Labellings (Part B), in the Slater DZP space and, as a second space, in cc-pVDZ (published):** *site* = AOs
  with local-frame p functions (radial, tangential, z), then symmetric Lowdin orthogonalization; *ring* = the real
  discrete Fourier transform over sites of the site orbitals, per local function type (k = 0, 3 and cos/sin pairs
  for k = 1, 2). Same space, an orthogonal rotation. Measures (mode `labels <basis>`): JW Pauli count at tolerances
  1e-10, 1e-6, 1e-4, 1e-3, 1e-2; the 1-norm; and the number of Z-type Z2 symmetries, computed as the dimension of
  the F2 null space of the parity rows of every nonzero integral (no group theory is fed in).
- **VOID:** a polarization fit that ends without success or at its iteration cap voids the predictions that use
  that basis. The cc-pVTZ pipeline check below failing by more than 1e-5 Eh voids Part A.

## What I know before comparing (disclosed)

- **Pipeline checks (measured):** H10 ring STO-6G at 1.0 A, our FCI -5.4229584 against Stair-Evangelista's
  -5.422958 (a published value; the new string FCI `fci_strings` agrees with the dense `GE.fci_ground` to 3e-15 on
  an H4 ring). H2 cc-pVTZ FCI -1.172335 against the literature -1.172336.
- **cc-pVTZ on the scored ring, computed before registration as a pipeline check:** its CCSD is compared against
  Magoulas-Shen-Piecuch (arXiv:2111.13787, Table I): FCI -3.410069 and CCSD error +4.424 mEh, so CCSD -3.405645.
  Ours: RHF -3.255794, CCSD -3.405645 (-3.4056446), off by 3.7e-7 Eh, inside the 1e-5 bar. Part A is not void.
  Only CCSD was run here, not (T). This value was published anyway.
- **cc-pVDZ (published exponents) on the scored ring was computed during a capability probe, before this record
  was drafted:** RHF -3.249608, CCSD -3.384403, CCSD(T) -3.387147. The Slater DZP and the re-optimized cc-pVDZ have
  not been computed on the ring.
- **Solver error on a ring (measured, non-scored):** H6 in a Slater [2s] basis (12 functions), CCSD(T) -3.322047
  against FCI -3.322348: 0.30 mEh. H10 STO-6G: 0.71 mEh. Magoulas: CCSD error 4.4 mEh at cc-pVTZ.
- **H2 control (measured, non-scored; same rule, R = 1.4, FCI = CCSD for two electrons):** Slater [2s1p] (p
  re-optimized: 1.791) FCI -1.167511, RHF -1.133419; cc-pVDZ (p re-optimized: 0.742) FCI -1.163403, RHF -1.128756.
  Gap 4.11 mEh; RHF gap 4.66 mEh (share 1.13: correlation gives some of it back). cc-pVDZ+TZpol -1.165854 (+1.66),
  aug-cc-pVDZ -1.164608 (+2.90), cc-pVTZ -1.172335 (-4.82): the crossing lies between cc-pVDZ+TZpol and cc-pVTZ.
- **Label controls (measured, non-scored, all at 1.0 A):**

| Ring, basis | Functions | Z2 site / ring | Pauli 1e-10 site / ring | 1-norm site / ring | Pauli 1e-3 site / ring | Pauli 1e-2 site / ring |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| H6 STO-6G | 6 | 2 / 4 | 1,795 / 407 | 12.0 / 14.7 | 379 / 407 | 139 / 355 |
| H10 STO-6G | 10 | 2 / 4 | 14,171 / 2,735 | 25.1 / 44.6 | 751 / 2,695 | 291 / 1,299 |
| H6 Slater [2s] | 12 | 2 / 4 | 29,449 / 6,681 | 60.7 / 62.1 | 3,229 / 3,481 | 541 / 825 |
| H4 cc-pVDZ | 20 | 3 / 5 | 122,197 / 30,629 | 383.8 / 350.5 | 24,909 / 17,645 | 2,373 / 3,521 |

  FCI energies in the two labellings agree to 1e-14 wherever FCI was run (H6 and H10 STO-6G, H6 Slater [2s]).
  The 1-norm ordering flips once p functions are present, so the 1-norm is reported, not predicted. The truncation
  crossover moves from between 1e-4 and 1e-3 up to between 1e-3 and 1e-2 once p functions are present.
- **The Z2 counts are derived, not guessed.** Both labellings conserve N_alpha and N_beta (two Z2) and, with p
  functions, the reflection z -> -z (one more: p_z functions are odd). Ring orbitals also carry a definite parity
  under the half-turn (k even or odd) and under the mirror j -> -j (cos even, sin odd). Site orbitals do not: those
  operations move them to other sites. Derived: 3 site, 5 ring with p; 2 and 4 without.

## Predictions

**Part A (energies, CCSD(T); the H2 gap times three, the number of H2 units, is 12.3 mEh)**

**HR-A1 (equal size, the fair form).** Slater DZP lies below cc-pVDZ with re-optimized p by more than 6.2 mEh (half
of three H2 gaps).
- Falsifier: 6.2 mEh or less, or the wrong sign.

**HR-A2 (scale, not gating).** That gap lies in [6.2, 24.7] mEh (half to twice three H2 gaps).

**HR-A3 (mean-field share).** The RHF part of the HR-A1 gap is more than half of it (H2: 1.13).
- Falsifier: 0.5 or less.

**HR-A4 (matched accuracy).** cc-pVDZ+TZpol (48) and aug-cc-pVDZ (54) lie above the Slater DZP, and cc-pVTZ (84)
lies below it.
- Falsifier: either of the first two at or below it, or cc-pVTZ above it.

**HR-A5 (shape-blind Pauli count).** The Slater DZP's JW Pauli count over canonical orbitals equals cc-pVDZ's
exactly (30 functions, 60 qubits).
- Falsifier: any difference.

**Part B (labels; scored in the Slater DZP space and again in the cc-pVDZ space; a prediction holds only if it
holds in both)**

**HR-B1 (symmetries, derived).** The Z2 finder returns exactly 3 for site labels and exactly 5 for ring labels.
- Falsifier: any other count.

**HR-B2 (exact sparsity).** At tolerance 1e-10 the ring-label Pauli count is below the site-label count by more than
a factor of 3 (controls 4.0 to 5.2).
- Falsifier: a ratio of 3 or less.

**HR-B3 (locality).** At tolerance 1e-2 the site-label Pauli count is below the ring-label count (all four
controls).
- Falsifier: site at or above ring.

**HR-B4 (reported, not predicted).** The 1-norms in both labellings; the counts at 1e-6, 1e-4 and 1e-3, and where
the crossover falls; the fraction of exactly zero two-electron integrals.

**HR-A6 (reported, not predicted).** The 1-norms over canonical orbitals; the cost ratios of the matched-accuracy
bases; where between cc-pVDZ+TZpol and cc-pVTZ the crossing falls; the cc-pVDZ published-to-re-optimized change.

## Novelty status (set before computing; sources: the two ring scans of 2026-11-01)

- **HR-A1, HR-A3, HR-A4, HR-A5: NOT FOUND.** The QC scan found one Slater-against-Gaussian comparison on hydrogen
  rings: Reinhardt, Toulouse, Assaraf, Umrigar and Hoggan (arXiv:1107.3537), one function per H, at the
  Hartree-Fock and QMC level, about where symmetry breaking sets in. No correlated comparison at equal size above
  minimal basis was found, and no matched-accuracy ladder. The search was one scan of arXiv and journal indexes; it
  is repeated at scoring time.
- **HR-B1: PROBABLY KNOWN.** Using point-group Z2 symmetries to remove qubits is standard practice. The tapering
  method is credited to Bravyi, Gambetta, Mezzacapo and Temme (arXiv:1701.08213), and point-group tapering of
  molecules to Setia and coworkers (2020). Neither was re-read for this record. The count for this ring is derived
  above.
- **HR-B2: PROBABLY KNOWN.** Momentum conservation zeroing the integrals of a ring or periodic system is
  textbook: Hueckel annulenes, and the momentum-space form of lattice models. No count for a Coulomb ring in a
  correlated basis was searched.
- **HR-B3: PROBABLY KNOWN in idea, NOT FOUND as stated.** Stair and Evangelista (arXiv:2005.11349) report that the
  best basis for a 1D ring depends on the method: DMRG does best with localized orbitals, the other methods with
  delocalized ones. A Pauli-count crossover under truncation between site and ring labels was not found.

## Controls

- **Positive (method):** the H10 STO-6G FCI and the H2 cc-pVTZ FCI above; the cc-pVTZ CCSD check on the ring.
- **Positive (labels):** wherever FCI was run, the two labellings give the same energy, to 1e-14 (they span the
  same space).
- **Non-scored calibration:** the H2 control (Part A bounds) and the four label controls (Part B bounds).

## Decision

- **GO:** HR-A1, HR-A3, HR-A4, HR-A5 and HR-B1, HR-B2, HR-B3 hold. HR-A2 is scored but does not gate.
- **BORDERLINE:** HR-A1, HR-A4 and HR-B1 to HR-B3 hold, and one of HR-A3, HR-A5 fails.
- **STOP:** HR-A1, HR-A4 or any of HR-B1 to HR-B3 fails, or a VOID.

## Reading if GO, stated in advance

Part A: Paper 62's equal-size result does not need a heavy atom; on a pure-hydrogen ring with no core the
exponential basis still wins at equal size and equal Pauli count, and the Gaussian basis needs more than
aug-cc-pVDZ to match it. Part B: on a ring one set of orbitals carries two label sets that do different jobs. Ring
momentum is the label that gives exact selection rules and taperable qubits; the site label is the one that gives
locality and truncation. This is the planar-ring form of the PI's prediction. It says nothing yet about a torus,
where the tube winding is only approximately conserved (Encinosa-Mott); that is a separate test.

## Disposition

Scored 2026-11-01 after the scored modes ran (rows in `debug/data/hring_study.jsonl`, modes `final` and `labels`;
exponents `debug/data/hring_exponents.json`). Driver unchanged since this record's commit (313996c0; sha256 in its
message). **Decision: STOP.** HR-A1 and HR-A4 failed. Part B held in full. Not VOID: both polarization fits ended
with success, well under the iteration cap (Slater p = 1.5518 in 7 evaluations, cc-pVDZ p = 0.6085 in 8), and the
cc-pVTZ pipeline check passed (3.7e-7 Eh).

| Basis (H6, 1.0 A) | Functions | RHF | CCSD | CCSD(T) | JW Pauli | 1-norm |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| Slater DZP | 30 | -3.253852 | -3.389089 | -3.391958 | 288,959 | 910.8 |
| cc-pVDZ, p re-optimized | 30 | -3.248909 | -3.385643 | -3.388492 | 284,895 | 944.2 |
| cc-pVDZ, published | 30 | -3.249608 | -3.384403 | -3.387147 | 286,927 | 987.1 |
| cc-pVDZ+TZpol | 48 | -3.250136 | -3.392375 | -3.395593 | 1,812,729 | 3290.7 |
| aug-cc-pVDZ | 54 | -3.249807 | -3.389846 | -3.393029 | 3,038,807 | 2846.7 |
| cc-pVTZ | 84 | -3.255794 | -3.405645 | -3.409413 | 16,902,041 | 14679.4 |

**Part A**
- **HR-A1 MISSED.** The Slater DZP lies 3.47 mEh below the re-optimized cc-pVDZ, short of the 6.2 mEh bar. The
  sign is right.
- **HR-A2 MISSED (not gating).** 3.47 mEh is below the band [6.2, 24.7].
- **HR-A3 HELD.** The RHF gap is 4.94 mEh, 1.43 times the CCSD(T) gap. Correlation gives back 1.48 mEh, more than
  on H2 (share 1.13).
- **HR-A4 MISSED.** cc-pVDZ+TZpol (48 functions) lies 3.64 mEh below the Slater DZP, and aug-cc-pVDZ (54) lies
  1.07 mEh below it. cc-pVTZ lies below it, as predicted. Paper 62's water result, aug-cc-pVDZ still 8 mHa short,
  does not carry over to this ring.
- **HR-A5 MISSED as registered.** Over canonical orbitals the counts differ: Slater 288,959 against cc-pVDZ
  286,927 (published) and 284,895 (re-optimized). *Diagnosis, after scoring:* the difference is not the radial
  shape. In both fixed labellings of Part B, the Slater DZP and cc-pVDZ give identical counts at 1e-10 (site
  653,971, ring 149,591 in both). The two cc-pVDZ forms, with the same shape family, also differ from each other
  over canonical orbitals. The D6h ring has degenerate orbital pairs, and the RHF solver returns an arbitrary
  rotation within each pair, which moves integrals across the 1e-10 threshold. Water has no degenerate orbitals,
  which is why Paper 62's canonical counts matched. The registered prediction did not anticipate this. It stands
  as missed.
- **HR-A6 (reported).**
  - The Slater DZP has the lowest 1-norm of the four 30-function rows, 3.5% below the re-optimized cc-pVDZ.
  - The cc-pVDZ+TZpol, the cheapest basis that beats the Slater DZP, costs 6.3 times the Pauli terms and 3.6
    times the 1-norm.
  - Re-optimizing the cc-pVDZ p exponent gains 1.34 mEh.
  - The crossing lies between 30 and 48 functions, at or below the cc-pVDZ+TZpol rung, lower than on H2, where it
    sat between cc-pVDZ+TZpol and cc-pVTZ.

**Part B (both spaces)**
- **HR-B1 HELD.** Z2 counts: 3 site and 5 ring, in the Slater DZP space and in the cc-pVDZ space.
- **HR-B2 HELD.** Ring labels cut the 1e-10 count by a factor of 4.37 in both spaces (653,971 to 149,591).
- **HR-B3 HELD.** At 1e-2 the site count is lower: 4,267 against 7,987 (Slater), and 4,171 against 7,979 (cc-pVDZ).
- **HR-B4 (reported).**
  - 1-norms, site against ring: 764.1 against 769.4 (Slater) and 746.4 against 779.2 (cc-pVDZ).
  - Exactly zero two-electron integrals: 44.8% (site) and 87.3% (ring).
  - The crossover lies between 1e-4 and 1e-3 in both spaces. At 1e-3 site is already lower: 51,115 against
    70,851 (Slater), 50,395 against 63,031 (cc-pVDZ). This is lower than in the H4 cc-pVDZ control, where it lay
    between 1e-3 and 1e-2.

**Novelty at scoring time.** One rescan, a web search for a correlated Slater-against-Gaussian comparison on hydrogen
rings at equal size. It found Reinhardt et al. again (Hartree-Fock and QMC, minimal basis) and nothing correlated.
Status unchanged.

**Reading (STOP).** Part A: on a pure-hydrogen ring the exponential basis still wins at equal size, but by about
3.5 mEh, not the 12 mEh that three H2 units would give. One extra shell of Gaussian polarization closes the gap.
Paper 62's matched-accuracy margin is a water result. It does not carry over to the ring, so it cannot be read as
general to hydrogen-rich molecules. Part B holds as stated in advance: one orbital space, two label sets, two jobs.
Ring momentum gives the exact zeros and two extra taperable qubits. Site labels give the lower count once small
terms are dropped.
