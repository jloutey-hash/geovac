# Water: the many-electron barrier (2026-10-10, PI-directed, after v5.20.10)

Canonical memo for this arc. The three-centre integrals are evaluable (v5.20.10:
shared-centre channels, ~1e-9 on water). This arc asks what stands between an
exponential-orbital water calculation and an accurate one.

## Step 0: where water stands (measured)

Minimal Slater basis (7 orbitals: O 1s 7.66, 2s 2.25, 2p 2.23; H 1s 1.0),
R_OH 1.809, 104.5 deg; Poly-0 pipeline (`debug/poly0_three_center_burden.py`),
integrals via McMurchie-Davidson over Gaussian fits of the Slater orbitals.

| n_gauss (evaluator) | RHF | FCI | correlation |
|---|---|---|---|
| 8 | -75.65806 | -75.70547 | -0.04740 |
| 12 | -75.66142 | -75.70883 | -0.04741 |

- **The evaluator is not converged at the mHa level:** 3.4 mHa between 8 and 12
  Gaussians, the same O 1s cusp lag seen in the pairwise referee. Exact
  exponential-orbital integrals remove this.
- **Basis is the barrier, not correlation method.** Literature anchors (to verify
  before quoting in a paper): HF limit about -76.07, exact non-relativistic about
  -76.44. Minimal-basis RHF is about 0.4 Ha above the HF limit; FCI recovers only
  0.047 of the roughly 0.37 Ha correlation energy.

## Plan, with gates

1. **Exact integral tensor for Slater-basis water.** Every class has an origin
   choice where each density holds an orbital at the origin, or one density is
   one-centre (finite L): channel engine for all two-electron classes; one-electron
   classes exact. Gate: reproduces the n_gauss -> large Gaussian limit; removes the
   3.4 mHa drift.
2. **Grow the basis** (double-zeta + polarization Slater, ~25 functions). Correlated
   solver needed (FCI out of reach). Gate: a strong baseline in the SAME pipeline:
   cc-pVDZ at equal size (benchmarking rule).
3. **The GeoVac question:** at matched accuracy, do exact exponential orbitals give
   a smaller basis and/or a cheaper qubit Hamiltonian (Pauli count, 1-norm) than
   Gaussians? Prior art: Slater bases are known to be more compact per function;
   the qubit-resource comparison is the part with a chance of being new.

## Step 2: Slater DZP vs cc-pVDZ at equal size (24 functions)

**Solver (own; PySCF does not build on this Windows host).** RHF + spin-orbital
CCSD + (T), `tests/test_water_ccsd.py`. Validated against the Crawford
programming-project STO-3G water reference: (T) to 2e-12, CCSD correlation to
4e-10, RHF to 3e-8. CCSD == FCI for H2 to 1e-9. Eight planted wrong answers
caught. (Found on the way: the basis parser read STO-3G's SP shells as two s
functions; cc-pVDZ, which has no SP shells, was unaffected.)

**Fairness rule: cc-pVDZ's own design, with Slater primitives.**
- O: primitives 1s,1s',2s,2s' / 2p,2p',2p'' optimized on the O average-of-
  configuration HF energy, contracted exactly as cc-pVDZ is:
  s = {atomic 1s, atomic 2s, most diffuse s primitive}, p = {atomic 2p, most diffuse p}.
  Exponents 9.317, 6.765 | 2.710, 1.694 | 5.719, 2.846, 1.450; E_av = -74.768205.
  The atom energy comes from a new exact one-centre radial HF (closed-form
  Slater F^k/G^k). Cross-checks: He -729/256 exact; Ne radial = the 3D Gaussian-fit
  limit, approached from above (-4.9, -0.20, -0.09 mHa at 8, 12, 16 Gaussians).
- Negative result recorded: an *uncontracted* 3s2p optimized on the atom spends
  its second s function on the core (10.03 / 8.33 / 2.32), which leaves the valence
  s single-zeta. That is why cc-pVDZ uses general contraction, and why the Slater
  basis here does too.
- The first 3D open-shell atom energy (fractional-SCF orbitals) sat 10.7 mHa above
  the variational average energy. Its orbitals are not variational for E_av.
  Superseded by the radial code.
- H: 1s,1s' optimized on H2 RHF at R = 1.4: 1.212, 0.627.
- Polarization (O 3d, H 2p): optimized on **frozen-core** CCSD for both bases.
  First tried all-electron: the energy kept falling with a tighter O 3d, and the
  frozen-core scan showed this was valence, not core (core part ~3.4 mHa,
  nearly constant). Frozen core is still the fair yardstick, since cc-pVDZ is a
  valence-correlation basis. cc-pVDZ d/p re-optimized the same way:
  1.185 -> 1.312, 0.727 -> 0.605. Slater: O 3d 2.756, H 2p 1.579 (1D parabolas on
  frozen-core CCSD with an 8/6-Gaussian evaluator, `debug/data/water_sto_polgrid.jsonl`).
  Slater 3d/2p shapes are diffuse at a given zeta (<r> = 7/(2 zeta), 5/(2 zeta)), so
  the optimum exponents look large next to Gaussian ones.

R_OH 1.809 bohr, 104.5 deg; 24 functions each; CCSD(T) with O 1s frozen ("fc") and
all-electron ("ae").

| basis (24 fns) | RHF | CCSD(T) fc | CCSD(T) ae |
|---|---|---|---|
| cc-pVDZ (standard) | -76.026795 | -76.241025 | -76.243143 |
| cc-pVDZ, pol re-opt on fc CCSD | -76.027257 | -76.242201 | -76.244383 |
| **Slater DZP, cc-pVDZ design rule** (12 Gaussians/fn, 8 on pol) | **-76.055573** | **-76.281739** | -76.285720 |
| same, 16 Gaussians/fn (evaluator check) | -76.055644 | -76.281811 | -76.285792 |
| cc-pVTZ (58 fns, 116 qubits; anchor) | -76.057162 | -76.332191 | -- |

- **Slater DZP is 39.5 mHa below cc-pVDZ at frozen-core CCSD(T), at equal size**
  (48 qubits). 28.8 mHa is SCF: RHF is ~12 mHa above the HF limit (literature
  ~-76.067, to verify) vs ~41 mHa for cc-pVDZ. 12.1 mHa is extra correlation.
- Evaluator converged: 12 -> 16 Gaussians moves every energy by -0.07 mHa. The
  fitted basis is a legitimate Gaussian basis with exact integrals, so its RHF
  is a strict upper bound for that basis.
- **Place on the Gaussian ladder (cc-pVTZ anchor).** The Slater DZP RHF is 1.5 mHa
  above cc-pVTZ, so it covers ~94% of the DZ -> TZ SCF gap. Frozen-core CCSD(T) is
  40.8 mHa below cc-pVDZ and 50.4 mHa above cc-pVTZ, ~45% of the DZ -> TZ gap.
  The exponential shape buys mainly the mean-field part (cusp + tail). Correlation
  is bought by the number of correlating functions, and there the gain is smaller
  (12 mHa).
- Prior art (checked 2026-10-11; moderate web search, primary sources where
  reachable). The SCF half ("STO beats GTO per function") is textbook and not new.
  The design recipe has precedent: Ema, ..., Fernandez Rico, Paldus, J. Comput.
  Chem. 24, 859 (2003), doi:10.1002/jcc.10227 (STO valence sets for H-Ne with HF+CI-
  optimized exponents, atoms only); Forster & Visscher, JCTC 2021,
  doi:10.1021/acs.jctc.1c00308 (correlation-consistent-style STO sets for GW). The
  closest molecular STO-vs-cc-pVXZ paper (Fernandez Rico et al., IJQC 100, 221 (2004),
  doi:10.1002/qua.20088) compares forces, not correlated energies. No matched-size
  correlated (CCSD(T)) STO-vs-GTO head-to-head was found, nor the SCF-vs-correlated
  split of the DZ->TZ gap. The qualitative reading "orbital shape matters little for
  correlation" is attributed to Persson & Taylor, JCP 1996 (Gaussian-geminal paper),
  UNVERIFIED (full text blocked). Still to read before any claim: Persson-Taylor 1996,
  the Madrid group's later molecular STO-CI papers, and the water CBS anchor in
  Helgaker, Klopper, Koch, Noga, JCP 106, 9639 (1997), doi:10.1063/1.473863.

## Step 3: the qubit question (done 2026-10-11)

Jordan-Wigner, RHF canonical orbitals, interleaved spin, terms |c| > 1e-10,
identity counted in N_Pauli and excluded from the 1-norm. Driver
debug/water_qubit_resources.py; counter tests/test_water_ccsd.py::jw_pauli_stats
(closed-form per string class, matches openfermion's full JW on STO-3G water,
five planted errors all caught). External check: frozen-core cc-pVDZ = 107,382
terms, exactly Trenev et al. 2025 Table 5. Data debug/data/water_qubit_resources.json.

| basis | qubits (fc) | N_Pauli fc | 1-norm fc (Ha) | N_Pauli all-e (Q) | 1-norm all-e | CCSD(T) fc |
|---|---|---|---|---|---|---|
| Slater DZP | 46 | 107,382 | 706.3 | 128,793 (48) | 754.6 | -76.2817 |
| cc-pVDZ | 46 | 107,382 | 674.3 | 128,793 (48) | 717.0 | -76.2410 |
| cc-pVTZ | 114 | 3,983,652 | 7,653.5 | 4,281,643 (116) | 7,858.6 | -76.3322 |

- The Pauli count is identical in both 24-function bases. It is fixed by the
  orbital count and the C2v symmetry pattern (both bases share the irrep
  composition), not by the radial shape.
- The 1-norm is 4.7% higher for Slater DZP (fc). So the 41 mHa gain over
  cc-pVDZ costs about 5% in 1-norm and nothing in qubits or terms.
- This is not yet a matched-accuracy comparison. Slater DZP sits between cc-pVDZ
  and cc-pVTZ, and cc-pVTZ costs 37x the terms and 11x the 1-norm for the other
  ~55% of the gap. A fair matched-accuracy test needs a Gaussian basis of
  intermediate size reaching -76.28 (e.g. cc-pVDZ plus added functions), which
  has not been done.
- Not GeoVac's basis. These are tuned STOs, not fixed-exponent hydrogenic or
  Sturmian functions. The dense multicentre integrals give no Gaunt sparsity on
  either side, so Paper 14's sparsity mechanism does not enter here.

## Step 4: GeoVac's own bricks in the same test (2026-10-11, PI-directed)

The shells n <= 3 on O (14 functions) and n <= 2 on H (5 functions) give 24 functions,
the same count as cc-pVDZ. There are two GeoVac-native radial rules, and no
polarization or exponent is tuned to water in either:

- **Hydrogenic** (the production rule of composed_qubit / two_center_eri): decay Z/n
  with the bare nuclear charge, so the kit has zero parameters. Builder
  `hydrogenic_kit`; its test checks that every function is a one-electron eigenfunction
  (-Z^2/2n^2) and that same-l functions are orthogonal.
- **Coulomb Sturmian** (the theory's basis, i.e. Avery's): one shared exponent per atom.
  k_H = 1, which makes the H 1s exact. k_O = 3.781 minimizes the radial O-atom energy;
  the best that one shared k can do for the atom is -71.406 Ha, 3.4 Ha above the
  Slater primitives' -74.768.

Frozen-core CCSD(T), same geometry, same solver:

| kit | functions | RHF | CCSD(T) fc | N_Pauli fc (46 q) | 1-norm fc |
|---|---|---|---|---|---|
| STO-3G | 7 | -74.9629 | -75.0123 | -- | -- |
| Coulomb Sturmian (k per atom) | 24 | -74.2879 | -74.5603 | 107,382 | 848.8 |
| hydrogenic Z/n, bare Z | 24 | -75.8397 | -76.0452 | 107,382 | 510.8 |
| cc-pVDZ | 24 | -76.0268 | -76.2410 | 107,382 | 674.3 |
| Slater DZP (tuned) | 24 | -76.0556 | -76.2817 | 107,382 | 706.3 |
| cc-pVTZ | 58 | -76.0572 | -76.3322 | 3,983,652 (114 q) | 7,653.5 |

Reading:
- Neither GeoVac-native kit survives the test. The hydrogenic kit is 196 mHa above
  cc-pVDZ at equal size, and 236 mHa above the tuned Slater kit. The shared-k Sturmian
  kit, at 24 functions, is worse than STO-3G at 7. So the 41 mHa exponential-shape win
  of step 2 came from exponent tuning, not from the hydrogen-like shape.
- Mechanism, hydrogenic: bare Z puts the O n=2 shell at zeta = 4, about twice as compact
  as the tuned valence (1.45-2.85). The valence then has to live in the n=3 functions
  (zeta 2.67), which leaves almost nothing to correlate with. The core is fine
  (zeta 8 vs 7.7-9.3 tuned). This is the job Z_eff screening does in composed_qubit;
  it is not tested here.
- Mechanism, Sturmian: one k cannot serve a core (~8) and a valence (~2) at once. The
  remedy in the Sturmian literature is generalized or Goscinskian Sturmians with
  several scales, which is Avery's route and is not built here.
- The hydrogenic kit has the lowest 1-norm (511 vs 674, -24%), but it is at a much worse
  energy, so this is not a resource win. All 24-function kits share one Pauli count
  (symmetry-fixed), so radial shape does not move it.

Data: debug/data/water_ccsd.jsonl, debug/data/water_qubit_resources.json.

## Step 5: the two remedies (2026-10-11, PI-directed)

- **Hydrogenic with composed's screening.** O 1s at the bare Z = 8; O shells n = 2-3
  at Z_eff_valence = 6 (atomic_classifier, i.e. Z minus the 1s^2 core). Zero parameters.
  This maps composed's core/valence blocks onto one 24-function set; composed itself
  starts its valence blocks at n = 1, so the mapping is ours.
- **Two-scale Coulomb Sturmian.** O: a core 1s at kc plus Sturmians sharing kv
  (s n <= 2, p n <= 3, 3d) = 14 functions. (kc, kv) = (7.793, 2.281) minimize the radial
  O-atom energy (-74.5027). H: k = 1. This is a simplified relative of Avery's
  generalized Sturmians, not his isoenergetic construction.

| kit (24 functions) | RHF | CCSD(T) fc | 1-norm fc |
|---|---|---|---|
| Sturmian, one k per atom | -74.2879 | -74.5603 | 848.8 |
| hydrogenic, Z_eff = 6 for n >= 2 | -75.7466 | -75.9471 | 405.1 |
| hydrogenic, bare Z | -75.8397 | -76.0452 | 510.8 |
| Sturmian, two-scale | -75.8903 | -76.1113 | 620.0 |
| cc-pVDZ | -76.0268 | -76.2410 | 674.3 |
| Slater DZP (tuned) | -76.0556 | -76.2817 | 706.3 |

(All 24-function kits: 107,382 Pauli terms frozen-core, 46 qubits.)

- **Screening made the hydrogenic kit worse, by 98 mHa.** The diagnosis comes from the
  O atom (average of configuration, 3D, ng 8) in the kits as built:
  - bare -74.557, Z_eff -74.420;
  - the same kits with the polynomial contractions released give -74.757
    (radial, upper bound on the contracted energy).
  So **200-340 mHa is lost to the fixed hydrogenic node structure**. Each hydrogenic
  function carries a fixed radial polynomial (its nodes), set for orthogonality to the
  1s of the *same* charge. Bare Z keeps every shell mutually consistent; Z_eff = 6 puts
  the valence nodes in the wrong place relative to the Z = 8 core. The contraction, not
  the decay rate, is the binding constraint.
- **Sturmians do not pay this.** A shared-k Sturmian set spans the same space as its
  primitives, so there is no contraction penalty; the two-scale kit is the best
  GeoVac-native kit, 130 mHa above cc-pVDZ.
- **No GeoVac-native kit reaches cc-pVDZ at 24 functions.** The best is 130 mHa behind
  with 2 dials on O; the zero-parameter kits are 196 and 294 mHa behind.
- **Observation, not a claim:** across the 24-function kits other than the
  one-k Sturmian, the 1-norm rises with accuracy (405, 511, 620, 674, 706 Ha). A low
  1-norm at equal size here tracks a worse description, not a cheaper one.
