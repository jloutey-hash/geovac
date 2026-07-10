# Sprint F4 — bonding-orbital Phillips-Kleinman closure attempt for W1e

**Date:** 2026-05-23 (same-day continuation of F3).
**Sprint position:** Named F3 follow-on attempting closure of the W1e
inner-region overattraction sub-layer. Per the algebraic-first sprint
mandate, gated on Step 1 diagnostic before any production wiring.
**Verdict line: STOP — algebraic Step 1 diagnostic confirms PK barrier
is structurally insufficient to close W1e (predicted closure ~0.1% of
the 4.37 Ha wall). Pivot to alternative W1e closure mechanisms.**
**Cross-references:** `debug/sprint_f3_cross_block_h1_memo.md` (W1d closed,
W1e named), `debug/sprint_w1c_full_arc_synthesis_memo.md` (full arc),
`geovac/phillips_kleinman_cross_center.py` (existing PK infrastructure
that F4 would have extended), `geovac/cross_block_h1.py` (F3 module),
CLAUDE.md §3 W1c-residual entries.

---

## §0. Executive summary + verdict

**Verdict line: STOP at Step 1 — bonding-orbital PK is structurally
insufficient to close W1e.**

The Step 1 algebraic diagnostic — computed in three complementary forms
(2×2 minimal subspace; full F3 h1 with explicit core overlaps; rank-1
PK injection FCI sensitivity test) — converges on a single quantitative
finding: the F3 bonding orbital at NaH max_n=2 has substantial overlap
with the Na [Ne] core (max |S| = 17.5% on Na 2s) but the resulting
Phillips-Kleinman barrier is only **+0.194 Ha** at the strict
absolute-PK reference $E_v = 0$. This is **22× smaller than the 4.37 Ha
wall depth** F4 was meant to close.

The FCI sensitivity test directly answers the gate question by
injecting a rank-1 PK shift $\Delta H = c \cdot |\text{bond}\rangle\langle
\text{bond}|$ on the F3 h1 at varying $c$. At the Step-1-predicted
$c = +0.19$ Ha, the FCI well depth changes by **0.05%** (4.3735 →
4.3712 Ha). Even at $c = 1.0$ Ha (5× the prediction) it closes only
14% of the wall. At $c = 2.0$ Ha (10× the prediction) it closes 43%
but saturates — still 33× experimental NaH $D_e \approx 0.075$ Ha.

The bonding orbital and the PK projection direction are well-separated
in the FCI's 2-electron singlet space at small R: the FCI ground state
puts significant weight on configurations OTHER than the bonding-pair
doubly-occupied state, and a rank-1 perturbation on the bonding orbital
does not lift those configurations. The W1e mechanism is structurally
distinct from a "bonding orbital orthogonal to core" reading — it
involves the 2-electron correlation channels in the FCI eigensolve,
not just the single-particle h1 eigenvalue.

Per the explicit Step 1 gate logic in the F4 sprint prompt:
- Strong gate (max |S| > 5% AND barrier > 1 Ha → PROCEED): FAILS
  (barrier 0.19 Ha < 1 Ha threshold).
- Stop gate (max |S| < 1% → STOP, not a Pauli problem): FAILS in the
  opposite direction (|S| = 17.5% > 1%).
- Intermediate gate (proceed but flag weak expected impact): formally
  triggered, but the FCI sensitivity test sharpens "weak" to
  "structurally insufficient" at the predicted magnitude.

The diagnostic-before-engineering rule (CLAUDE.md §1.7,
`feedback_diagnostic_before_engineering.md`) applies decisively here:
the FCI sensitivity test at Step 1 cost ~5 minutes and ruled out the
1-week F4 implementation sprint cleanly. The expected outcome of F4
Step 3 with full production wiring (extending
`phillips_kleinman_cross_center.py` to operate on the cross-block-h1-
constructed bonding orbital) would have been: D_e = 4.37 → 4.37 Ha
(unchanged to 0.05%), no internal equilibrium minimum, no closure of
W1e.

### Predicted (and now empirically rejected) D_e^F4: **+4.37 Ha**
(essentially unchanged from F3 baseline). No internal minimum.

---

## §1. Step 1 — bonding-vs-core overlap diagnostic

### §1.1 Setup

At NaH $R = R_\text{eq} = 3.566$ bohr with the full F3 stack (W1c
screening + multi-zeta Na 3s + cross-block h1), compute the bonding
orbital from h1 diagonalization and compute its overlap with each Na
[Ne] core orbital (1s, 2s, 2p).

Three diagnostics performed:

1. **2×2 minimal subspace** (`debug/sprint_f4_step1_diagnostic.py`):
   diagonalize h1 in {|Na 3s⟩, |H 1s⟩} subspace to obtain bonding
   coefficients, then compute overlaps with Na core.
2. **Full F3 h1** (`debug/sprint_f4_step1_full_h1_diag.py`):
   diagonalize the full 10×10 h1 from `build_balanced_hamiltonian` with
   the F3 stack enabled, identify the lowest eigenvector (bonding
   orbital), compute Na core s-orbital overlaps.
3. **Na 2p_z extension** (`debug/sprint_f4_step1_with_2p.py`): add the
   $m=0$ component of Na 2p core via axial quadrature.

### §1.2 Results

**2×2 minimal subspace (multi-zeta Na 3s):**
- $h_{12} = -1.371$ Ha, $S_{AB} = +0.397$
- Bonding eigvec: $c_\text{Na3s} = -0.542, c_\text{H1s} = -0.653$ at
  $E_\text{bond} = -1.232$ Ha
- Na core overlaps:
  - (1s): $|S| = 1.08\%$, $\Delta H^\text{PK}_\text{Ev=0} = +0.005$ Ha
  - (2s): $|S| = 6.02\%$, $\Delta H^\text{PK}_\text{Ev=0} = +0.010$ Ha
  - (2p): $|S| = 2.57\%$, $\Delta H^\text{PK}_\text{Ev=0} = +0.001$ Ha
- **Total PK barrier (2×2): +0.016 Ha**

**Full F3 h1 (10 orbitals):**
- Lowest h1 eigenvector at $E = -3.14$ Ha. Composition:
  - NaH_bond_center(1,0,0) = −0.620 (Na 3s)
  - NaH_bond_center(2,0,0) = +0.352 (Na 4s)
  - NaH_bond_partner(1,0,0) = −0.467 (H 1s)
  - NaH_bond_partner(2,0,0) = +0.522 (H 2s)
- Na core overlaps (sum of weighted single-orbital overlaps):
  - (1s): $|S| = 5.16\%$, $\Delta H^\text{PK}_\text{Ev=0} = +0.108$ Ha
  - (2s): $|S| = 17.52\%$, $\Delta H^\text{PK}_\text{Ev=0} = +0.086$ Ha
- **Total PK barrier (full h1, s-cores only): +0.194 Ha**

**Na 2p_z extension:**
- Sum of l=1, m=0 components on bonding: c_Na_2p_z = −0.0193,
  c_H_2p_z = +0.0286 (both small)
- $S_\text{bond,2pz} = -0.0032$ ($|S| = 0.32\%$)
- $\Delta H^\text{PK}_{2p_z} = +0.0000$ Ha (negligible)
- **Total PK barrier (s + 2p_z): +0.194 Ha** (essentially unchanged)

### §1.3 Mechanism

The full-F3 bonding orbital has SUBSTANTIAL overlap with the Na 2s core
(17.5%) — well above the 5% strong-gate threshold. The mechanism is the
cross-block h1's mixing of Na block_n=1 (physical n=3) AND block_n=2
(physical n=4) states with H block_n=1 (physical n=1) AND block_n=2
states. The Na 2s core (Bethe-Salpeter conv: Z_eff=6.4 with
Clementi-Raimondi zeta_2s=3.218) extends to ~0.5 bohr; the Na 3s
multi-zeta (physical mean radius 4.5 bohr) and the Na 4s component have
non-trivial amplitude at small r where the core lives. The cross-block
h1's energy-minimizing mixing pulls the bonding orbital INWARD slightly
toward the Na nucleus, giving it 17.5% overlap with the Na 2s core.

By contrast the predicted PK barrier weight ($E_v=0 - E_c = +2.797$ Ha
for Na 2s) times $|S|^2 = 0.031$ gives only 0.086 Ha. Even the Na 1s
core ($E_c = -40.5$ Ha) contributes only 0.108 Ha despite its large
energy magnitude because $|S|^2 = 0.0027$ is small.

### §1.4 FCI sensitivity test (rank-1 PK injection)

**Most important Step 1 result for the gate decision.** Rather than
wiring the full bonding-orbital PK extension into production code,
inject the predicted barrier directly as a rank-1 perturbation on h1:

$$\Delta H = c \cdot |\text{bond}\rangle\langle\text{bond}|$$

where $|\text{bond}\rangle$ is the lowest h1 eigenvector at each R, and
$c$ varies from 0 to 5.0 Ha. Run 2-electron FCI at $R \in \{2.0, 3.5,
10.0\}$ bohr for each $c$.

| $c$ (Ha) | $D_e$ (Ha) | reduction | fraction |
|---:|---:|---:|---:|
| +0.00 | +4.3735 | +0.0000 | 0.00% |
| +0.05 | +4.3730 | +0.0005 | 0.01% |
| +0.10 | +4.3724 | +0.0010 | 0.02% |
| **+0.20** | **+4.3712** | **+0.0023** | **0.05%** |
| +0.50 | +4.3648 | +0.0087 | 0.20% |
| +1.00 | +3.7714 | +0.6020 | 13.77% |
| +2.00 | +2.4998 | +1.8737 | 42.84% |
| +5.00 | +2.4823 | +1.8912 | 43.24% |

**At Step 1 predicted barrier $c = +0.19$ Ha, the FCI well depth
reduces by 0.05%.** The rank-1 PK shift on the bonding orbital is
structurally near-orthogonal to the 2-electron FCI's dominant
correlation channels at small R; the perturbation cannot lift the
configurations that produce the inner-region overattraction.

Even at saturation ($c \to \infty$), the closure fraction asymptotes
to ~43% — leaving the residual at ~2.5 Ha, still 33× experimental NaH
$D_e \approx 0.075$ Ha. This is the structural ceiling on what a
bonding-orbital PK can accomplish: it can shift the bonding orbital
out of the inner region, but it cannot remove the 2-electron
correlation channels that the FCI uses to find lower-energy
configurations involving other orbital combinations.

### §1.5 Step 1 verdict

Per the explicit F4 sprint gate logic:

| Criterion | Threshold | Actual | Verdict |
|:---|:---|---:|:---:|
| max |S| with at least one core | > 5% | 17.5% | PASS |
| AND predicted PK barrier | > 1 Ha | 0.194 Ha | FAIL |
| Strong gate | both above | — | **FAIL** |
| Intermediate gate | 1% < max |S| < 5% OR barrier intermediate | barrier intermediate | TRIGGERED |
| Stop gate | max |S| < 1% | 17.5% | NOT TRIGGERED |

Formally the gate triggers "proceed but flag weak expected impact."
But the FCI sensitivity test sharpens this to "structurally
insufficient" — at the predicted barrier magnitude, the closure is
0.05% of the wall, far below the 2× experimental window
$[0.0375, 0.150]$ Ha that would constitute closure.

**Decision: STOP the F4 sprint at Step 1.** Pivot to alternative W1e
closure mechanisms (see §5).

The diagnostic-before-engineering rule pays off: a ~5-minute Step 1
diagnostic ruled out a 1-week production extension cleanly. The
alternative would have been to extend `phillips_kleinman_cross_center`
to operate on the cross-block-h1 bonding orbital (~2 hours wiring +
1-week test cycle), produce a near-identical 4.37 Ha well depth
result, and only THEN discover the PK was insufficient.

---

## §2. (Not executed) Step 2 — Extend PK module + wire kwarg

Per the Step 1 STOP gate, Step 2 was not pursued. The intended
extension of `geovac/phillips_kleinman_cross_center.py` would have
been:

1. Accept arbitrary orbital combinations (not just single basis
   orbitals) — extend `compute_pk_cross_center_barrier` to take a
   bonding-orbital coefficient vector instead of single-orbital labels.
2. Wire `bonding_pk=False` kwarg in `balanced_coupled.build_balanced_hamiltonian`.
3. Add new module `geovac/bonding_orbital_pk.py` (~250 lines, s-s pairs only).
4. Add tests in `tests/test_bonding_pk.py` (~5 tests).

These are queued as named follow-ons if the W1e closure path is
reopened with a different operator architecture.

---

## §3. (Not executed) Step 3 — NaH FCI test with bonding-orbital PK

Per the Step 1 STOP gate, Step 3 was not pursued. The Step 1
FCI sensitivity test (§1.4) is a near-complete substitute: it
directly answers the question "what does the F3 stack + a rank-1
bonding-orbital PK barrier give as $D_e$ at NaH max_n=2?" by
injecting the predicted PK shift directly into h1.

Predicted Step 3 outcome at $c = +0.19$ Ha (matching the Step 1
prediction):
- $D_e^\text{F4} = 4.371$ Ha (vs F3 baseline 4.374 Ha)
- $R_\text{min} = 2.0$ bohr (unchanged)
- No internal equilibrium minimum
- Naturals signature unchanged ($\sim [1.999, 0.001]$ bonding)

This matches the F3 outcome to within 0.05% and would not have
satisfied the F4 binding criteria.

---

## §4. (Not executed) Step 4 — Mini-PES

Not pursued; conditional on Step 3 producing physically reasonable
D_e (it would not have).

---

## §5. Wall taxonomy update — W1e status after F4 Step 1

### §5.1 W1e sub-layers refined

The Step 1 diagnostic refines the W1e sub-layer characterization. W1e
is NOT a "bonding orbital not orthogonal to core" mechanism alone — the
bonding orbital DOES have substantial core overlap (17.5%), but the
PK-class barrier from that overlap is small (0.19 Ha) relative to the
wall depth (4.37 Ha).

The actual W1e mechanism is a 2-electron FCI correlation effect: the
FCI ground state at small R is NOT exclusively the doubly-occupied
bonding-orbital configuration. It puts significant weight on
configurations involving other orbitals — including, plausibly,
configurations with one electron in the bonding orbital and another
in a less-bonding orbital that has its own (possibly large) energy
lowering at small R. The rank-1 PK on the bonding orbital cannot
suppress those other configurations.

**Refined W1e characterization:** W1e is a 2-electron FCI correlation
wall, not a single-particle Pauli-orthogonality wall. The bonding
orbital is necessary for ANY bond formation (closed in F3) but the
2-electron correlation channels in the FCI eigensolve are what produce
the over-attraction at small R.

### §5.2 Candidate W1e closures (new diagnosis)

Three alternative mechanisms for W1e closure, ranked by structural
plausibility:

1. **Explicit frozen-core electrons** (highest priority): lift to a
   "4 + N_core" electron FCI rather than 2-electron, including the [Ne]
   core orbitals as additional FCI orbitals with appropriate occupancy
   constraints. The 2-body Coulomb repulsion between the bonding pair
   and core electrons acts as Pauli exclusion at the correlated level
   (not the single-particle level). This would address the actual
   mechanism (FCI correlation channels) rather than the PK-orthogonality
   approximation. **Implementation cost:** moderate, requires extending
   the FCI driver to handle occupancy constraints (~2-week sprint).
   **Expected closure:** physically motivated, could close the wall to
   within ~2× experimental.

2. **Strict Schmidt orthogonalization** in the basis-construction step:
   orthogonalize the heavy-atom-side valence orbitals against the
   [Ne] core orbitals BEFORE forming the cross-block h1. This is the
   mathematically definitive solution; the bonding orbital construction
   would proceed in an orthogonal basis from the start. **Implementation
   cost:** high, requires reworking the basis construction in
   `composed_qubit` and propagating through the cross-block h1
   integration (~3-4 week sprint). **Expected closure:** depends on
   whether the FCI correlation channels are also addressed by
   orthogonality (open).

3. **Bonding-orbital occupation projector** (low priority, possibly
   negative): replace the rank-1 PK barrier with a HIGHER-rank
   projector that lifts not just the bonding orbital but also
   correlation-relevant configurations. Structurally ad hoc; would
   amount to a fitted constraint.

The chemistry arc's diagnostic-before-engineering discipline (CLAUDE.md
§1.7) applies here too: before committing to candidate #1 or #2 above,
a Step 1-class diagnostic should test whether the predicted closure
mechanism actually addresses the FCI correlation channels. The rank-1
PK test in §1.4 IS that diagnostic for the PK family; an analogous
diagnostic for the explicit-core-electron candidate would inject the
core orbitals into the FCI sub-space and verify the 2-electron
correlation channels are addressed.

### §5.3 Sub-layer hierarchy update

| Sub-layer | Mechanism | Status |
|:---|:---|:---|
| W1c-cross-screening | Frozen-core $Z_\text{eff}(r)$ reduces cross-V_ne | CLOSED (Phase C-W1c) |
| W1c-multi-zeta-basis | Physical Na 3s replaces hydrogenic | CLOSED (F1-P1+P2) |
| W1d-cross-block-h1 | Off-diagonal h1 between centers | CLOSED (F3) |
| **W1e-FCI-correlation** | **2-electron correlation channels create inner-region overattraction; PK on bonding orbital insufficient** | **NEWLY DIAGNOSED (F4 Step 1)** |

W1e was originally named in F3 as "inner-region overattraction" with
the candidate mechanism "missing Pauli repulsion against frozen core."
F4 Step 1 sharpens this to "FCI correlation channels," shifting the
candidate closure path from PK-orthogonality to explicit-core-electrons
or Schmidt-basis-orthogonalization.

---

## §6. Recommended next sprint after F4

### Priority 1 — Sprint F5: explicit frozen-core electrons in FCI (~2 weeks)

The most structurally plausible W1e closure candidate. Extend the FCI
driver in `balanced_coupled` (or the dispatched driver in
`debug/sprint_f3_step3_fci.py`) to handle FCI with frozen-core
orbitals included as additional spatial orbitals with constrained
occupancy = 2. The cross-block 2-body Coulomb repulsion between the
bonding pair and core electrons acts as Pauli exclusion at the
correlated level.

**Decision gate:** an internal minimum at NaH max_n=2 with $R_\text{eq}$
within 1 bohr of 3.566 and well depth in $[0.0375, 0.150]$ Ha
closes W1e and the W1c-residual wall to within 2× experimental.

**Risk:** the explicit core electrons increase FCI dim by $\binom{10+N_\text{core}}{N_e+N_\text{core}}$
which at NaH ([Ne] core, 10 valence electrons in 5 valence orbitals at
max_n=2) becomes prohibitive without occupancy constraints. Need to
implement the occupancy-constraint FCI carefully.

### Priority 2 — Sprint F4-extend: Schmidt orthogonalization (~3-4 weeks)

Orthogonalize the heavy-atom-side valence orbitals against the [Ne]
core BEFORE the cross-block h1 step. Modifies `composed_qubit` and
propagates through `cross_block_h1`. Higher architectural lift but
mathematically definitive at the basis-construction level.

**Risk:** unclear whether Schmidt orthogonalization addresses the FCI
correlation channels OR just the PK orthogonality. Diagnostic needed
before commit.

### Priority 3 — Pivot scoping sprint: pause chemistry arc, consolidate

The chemistry arc has produced F1 → F2 → F3 → F4 in rapid succession
on 2026-05-23, each sprint refining the W1c-residual sub-layer
characterization. Natural pause point is now, before the F5 (or
F4-extend) implementation, to:

1. Update CLAUDE.md §3 with the F4 Step 1 negative finding
2. Update `papers/group2_quantum_chemistry/paper_17_composed_geometries.tex` §6.10 with the F4 STOP verdict
3. Update `papers/group2_quantum_chemistry/paper_19_coupled_composition.tex` §sec:w1c_residual with the W1e refinement
4. Pivot back to math.OA paper drafts (Papers 38/39/40/45/46/47 polish)
   if PI wants a documentation pass before further chemistry compute.

### Recommendation

**Default:** apply the documentation batch from F3 + F4 (paper edits +
CLAUDE.md updates + memory file), then PI decides between F5
(explicit-core FCI) or pivot to math.OA.

**If pursuing F5:** the explicit-core-electron approach is the most
physically motivated and addresses the actual mechanism W1e reveals
(FCI correlation channels). The 2-week sprint is feasible with the
existing FCI machinery; the main work is the occupancy constraint.

**If pivoting:** F4 Step 1 provides a clean "PK is structurally
insufficient; W1e is a correlation effect" narrative that completes
the W1c-residual story at F3-maturity. Math.OA arc has clear
deliverables ready.

---

## §7. Files

### Created (debug drivers)
- `debug/sprint_f4_step1_diagnostic.py` — 2×2 minimal-subspace
  bonding-vs-core overlap test.
- `debug/sprint_f4_step1_full_h1_diag.py` — full F3 10-orbital h1 with
  proper bonding-orbital identification and per-core overlap.
- `debug/sprint_f4_step1_with_2p.py` — extends to Na 2p_z core overlap
  via axial quadrature.
- `debug/sprint_f4_step1_fci_sensitivity.py` — rank-1 PK injection FCI
  sensitivity test (most decisive Step 1 result).

### Created (data)
- `debug/data/sprint_f4_step1_diagnostic.json` (2×2 result)
- `debug/data/sprint_f4_step1_full_h1_diag.json` (full-h1 result)
- `debug/data/sprint_f4_step1_with_2p.json` (2p_z addition)
- `debug/data/sprint_f4_step1_fci_sensitivity.json` (rank-1 PK FCI)
- `debug/data/sprint_f4_results.json` — consolidated results +
  verdict.

### Created (memo)
- `debug/sprint_f4_bonding_pk_memo.md` — this memo.

### NOT modified
- Production `geovac/` modules — no changes per the STOP gate.
- Tests — no changes; baseline regression preserved.
- Papers — per the sprint mandate, paper edits deferred to a
  comprehensive β sprint after the chemistry arc consolidation.
- CLAUDE.md — same deferral.

### Regression
- No production code modified; regression skipped per the explicit
  mandate. The Sprint F3 regression (146 + 1 skipped) preserved
  bit-exactly as-is.

---

**End of Sprint F4 memo. Verdict: STOP at Step 1 — PK structurally
insufficient.** F3 W1d closure stands; W1e re-diagnosed as a 2-electron
FCI correlation wall, not a single-particle Pauli wall. Candidate
closure paths refined: explicit frozen-core electrons (F5 priority 1) or
Schmidt orthogonalization (F4-extend priority 2). The diagnostic-before-
engineering discipline saved a ~1-week production sprint that would
have produced a near-null result. The bonding-orbital PK extension
infrastructure can be revisited if the actual W1e mechanism turns out
to be addressable by a higher-rank projector (currently judged unlikely
based on §5.2 analysis).
