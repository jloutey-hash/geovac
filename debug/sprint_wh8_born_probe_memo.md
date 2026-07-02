# Sprint WH8 Step-1 — Born-measure probe (solo, 2026-07-01/02)

**PI-directed** (group4-recert pause; "run that sprint... solo it without subagents").
**Question:** does the graph's skeleton-side measure generate Born statistics?
**Verdict:** **NEGATIVE (as expected)** — the Born rule's Class-1 placement (external-input
three-class partition) is upgraded from classification to **tested negative**; **WH8
registered** (CLAUDE.md §1.7, PI direction).

## Prior art (current-state check — found mid-sprint)

Paper 34 §VIII already carried a **CLOSED structural position** (2026-05-26, sprint
`ahha_born_rule_attempt`): three derivation routes (§III.28+GNS; Pythagorean
orthogonality + preferred basis; M1 modular trace) all reduce to standard GNS/Gleason
inheritance — no GeoVac-specific reduction of Born to substrate primitives; framework
placed in the Wigner/calibration-tier camp, explicitly not improving on Gleason. That
closure was **derivation-direction and argumentative** (no test artifact). This sprint
adds the complementary **generation-direction, quantitative** leg + the frozen falsifier
the 05-26 closure lacked. (PM note: I initially told the PI "nobody ran it" — half
wrong; corrected in-session.)

## The probe (debug/wh8_born_measure_probe.py → debug/data/wh8_born_measure_probe.json)

The skeleton's only native measures are counting-type (rational, built from discrete
invariants: nodes, degeneracies, sectors). Three legs:

- **A — coincidence locus (exact).** Counting = Born **iff** the state is flat
  (uniform |amplitude| over its support). On the hydrogen n=2 multiplet (d=4):
  TV(flat)=0 exactly; amplitudes ∝(1,2,3,4) give TV = **1/3 exactly** (rational
  identity, pinned in exact Fractions). The coincidence is real — it is why
  "graph-native probability" survives casual inspection — and it is precisely the
  flat locus, nothing more.
- **B — physical divergence.** He graph-native CI (max_n=3, `slater_full`,
  zero parameters; N_SD=378; E0=−2.840439 Ha, consistent with the known series →
  0.19% at n_max=7): ground-state Born weight is **90.8% on 1s²**; TV vs
  uniform-on-support = **0.858**, vs uniform-on-basis 0.986, vs uniform-on-S_z=0
  sector 0.974. The most charitable skeleton candidate misses by O(1).
- **C — dynamical instability.** A flat state on the GS support evolves off the
  coincidence locus under the graph's own H: TV ~ t², measured small-t exponent
  **2.00**; TV(t=0.5)=6.1e−3. The locus is measure-zero AND not dynamically
  preserved — no state-independent skeleton rule can track Born through dynamics.

**Fence (no compute):** given the projection lattice (dim≥3), any countably-additive
non-contextual measure IS Born (Gleason 1957) — the import is unique; the residual
mystery is single outcomes, which no interpretation resolves without added structure.

## WH8 (registered, PI direction)

Born measure = the **exchange constant of the observation projection** (§III.28) —
the projection's irreducible imported content, exactly as π is the Coulomb
projection's import and 2π/β the temporal one's (WH7 linkage: the same compact
observer window that discretizes time injects |ψ|² at registration; menu-discreteness
= compactness of records; KMS diagonality = four-witness). Falsifier + status: see
CLAUDE.md §1.7 WH8.

## Artifacts

- Driver `debug/wh8_born_measure_probe.py`; data `debug/data/wh8_born_measure_probe.json`.
- Frozen falsifier `tests/test_wh8_born_probe.py` (3 passed; self-contained recompute,
  no debug/ import).
- Paper 34 §VIII Born entry: quantitative-converse addendum (tier-honest, §1.5-safe).
- CLAUDE.md: §1.7 WH8 registered; §3 dead-end row (counting measure as Born generator);
  §2 one-liner. Memory `external_input_three_class_partition.md` updated.

## Honest scope

This resolves the *discreteness* and *effective-diagonality* components of measurement
into skeleton-side structure and fences the *selection* component as a unique,
non-derivable, now-tested import. It does not touch the definite-outcome problem.
Papers stay under §1.5 throughout; the bold reading lives in the register only.
