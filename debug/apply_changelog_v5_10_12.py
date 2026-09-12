"""Prepend the v5.10.12 CHANGELOG entry and update the CLAUDE.md Sec.2 index."""
import io

ENTRY = """## [v5.10.12] - 2026-09-08

**Paper 60's accuracy floor is the scale lock, and it is a ground-state pathology.** Both halves of the previous attribution -- angular truncation, and "cannot reach chemical accuracy at any `K`" -- are retired. PI-directed diagnostic; not a `/qa` run. Drivers `debug/p60_{variational_probe,scale_scan,freescale_resource,posing_cost_by_state,excited_ladder,floor_windows,floor_models,stateprep_overlap}.py`; backing `tests/test_paper60_scale_lock.py`.

### The identity

The metric-free posing is not a separate construction. Written at a free global scale `λ`, the Goscinskian one-body Coulomb metric is exactly diagonal --- `W_μν = R_ν δ_μν`, new `eq:W_diagonal`, verified entrywise to 5e-11 --- so the kinetic matrix at unit scale is `T = 1 − S/2` and the variational problem over the same span is `H(λ)C = E·S·C`. Substituting `E = −λ²/2` cancels every `S` term **identically** and returns `eq:secular`. Hence, new `eq:scale_lock`:

    metric-free  ⟺  E = −λ²/2  ⟺  λ = p_κ = √(−2E)

The posing *is* the variational problem of its own span, evaluated at the one scale where the L² metric drops out. Confirmed to 2.3e-12 (`l_max=0`) and 7.1e-13 (`l_max=3`) by two disjoint assemblies sharing one primitive.

### What the lock costs, and what it does not

The lock is not the variational optimum. Freeing `λ` over the **identical** span reaches 1.28 mHa at K=130 (spdf) against 7.46 locked, and 0.15 mHa of the independently known He s-limit at K=136 against 4.43 locked. **The span was never the limitation.** The price of freeing it is the whole encoding advantage: `‖M‖₁ ~ K^0.72` becomes `‖H(λ*)‖₁ ~ K^1.95` and `‖S^{-1/2}HS^{-1/2}‖₁ ~ K^2.75`, with `cond(S) ~ K^0.94` --- a factor 5e3 at K=136. Metric-free posing, sublinear 1-norm and accuracy floor are **one fact, not three**, and the two cost risks the abstract advertises as structurally absent return together. Two corollaries: the variational bound is automatic (the isoenergetic root is the lowest root of `H(p_κ)C = E·S·C`), and the He `l≥4` partial-wave tail is 0.187 mHa (ground) / 0.008 mHa (2¹S) --- 3% and 0.5% of the floors --- so angular truncation was never a candidate.

### It is a ground-state pathology

`‖M‖₁` is a property of `M`, not of which root is extracted, so every ¹S state costs the same to encode and only the delivered accuracy differs. At K=452 the ground state sits **4.28×** above chemical accuracy and 2¹S sits **1.08×**. Reference-free, the posing cost at K=105 runs 4.212 / 0.983 / 0.323 / 0.125 mHa across the first four roots (reductions 4.3×, 3.0×, 2.6×). Rydberg gaps do not bind at chemical accuracy --- margins 254 / 27 / 8.6 / 3.9 --- so the query count ~2e5 is state-independent through 4¹S. State preparation is modest: L²-metric overlaps 0.992 / 0.798 / 0.864 / 0.889, so 2¹S is the *hardest* of the four and interior roots get easier again. One practical trap: the dominant configuration does not track the spectroscopic label (the root at He 3¹S is dominated by `(0,1,4)`), so a prep heuristic keyed to the physical label picks the wrong configuration.

### Prior art: the mechanism is Avery's, the price is ours

Relayed from Avery's canon (PI-relayed source with book access) and confirmed by a reference-free test here: helium's ground state needs **in-out radial correlation**, supplied variationally by split-shell `1s1s'` functions with two independent exponents. A Goscinskian `1s²` configuration puts both electrons at `n=1` and hence at one exponent --- **for any weighting potential**, since a separable `V₀` still assigns one `β_ν` per configuration. Excited configurations get two scales free from `n_a ≠ n_b`. The metric-free form is likewise **general in `V₀`** (subtracting the two Sturmian equations gives `(β_μ − β_ν)⟨Φ_μ|V₀|Φ_ν⟩ = 0`), derived independently here and by the relay. But **parameter-freeness does not survive a general `V₀`**, and a *constant* rescaling is a no-op since `β_ν Z_w = p_κ/R_ν` eliminates `Z_w`. The split to carry forward: **mechanism = Avery's; its price in qubits = ours.** `memory/avery_method_and_prior_art_gaps.md` updated.

**Also a re-instance of a documented wall:** §3's "Fock energy-shell self-consistency for He --- `k² = −2E` over-constrains 2-electron problem; single parameter insufficient" (Track DI Sprint 2) is this same finding reached from the opposite direction.

### Withdrawn

The comparison to "−2.90250 Ha with 102 optimized Coulomb--Sturmian configurations" is **withdrawn as a characterization of this method**. That figure (1.2 mHa short) matches our *scale-optimized* ladder (1.64 mHa at K=100) and not our *locked-scale* one (7.70 mHa at the same K), and which posing produced it could not be confirmed against the primary source. It was also never in the numeric registry, so no gate could have caught it. **If the locked-scale method really does reach 1.2 mHa at K≈102, the floor is wrong and the error is ours** --- recorded in memory so it does not evaporate.

### Honest limits

The floors are **bracketed, not pinned**, and the dominant uncertainty is not the one we looked for first. Across 21 fit windows the free-floor value drifts 1.3% (ground) / 0.4% (2¹S) and drifts *upward*; Shanks descends from above, bracketing [6.47, 6.62] and [1.647, 1.676] mHa. But *model family* is the real exposure: a `c + b/lnK` form, rejected at 340× worse RMS, puts the 2¹S floor **below** chemical accuracy at 0.82×. The paper therefore cites measured ladder endpoints throughout and states that "2¹S saturates above chemical accuracy" is a model-selection conclusion.

### Instruments, and what they caught

- Registry: 7 new entries plus 2 **derived** ratios, so the ×-chemical-accuracy figures recompute from the energies and cannot drift apart. `p60_energy_floor` provenance annotated with the mechanism.
- **C21 caught the author.** `\\gvq{p60_exc_gap_k202}{1.12}` keyed a *ratio* to an entry holding an *energy* --- the exact mis-keying class the gate exists for. Fixed by registering the ratios as derived.
- **C16 caught an unrevisited dependent.** `docs/qa/paper_60.done.md` was *ratifying* the retired claim (the ".done.md as re-infection vector" class named on 2026-09-08). Revisited and given a supersession note, then stamped --- not stamped blind.
- **The new backing test caught two prose errors written this same sprint**: "the ratio widens with `K`" (it narrows, 5.17 → 4.29; what widens is the absolute separation, which then flattens) and "threefold to fourfold per rung" (the third rung is 2.6×). Both were live at five loci including the synthesis and the registry provenance. This is the §9 guard-separation rule paying out inside one day.
- The group2 synthesis had restated the floor claim in its own words --- the owner-corrected / citer-stale shape §9 records as eight of ten recurring defects. Fixed in the same pass, with a C16 entry naming both loci.

"""

ch = io.open("CHANGELOG.md", encoding="utf-8").read()
anchor = "## [v5.10.11] - 2026-09-08"
assert anchor in ch and "v5.10.12" not in ch
ch = ch.replace(anchor, ENTRY + anchor, 1)
io.open("CHANGELOG.md", "w", encoding="utf-8").write(ch)
print("CHANGELOG: v5.10.12 entry prepended")

cm = io.open("CLAUDE.md", encoding="utf-8").read()
BULLET = ("- **P60 floor is the scale lock, not l_max (2026-09-08, v5.10.12):** "
          "metric-free iff lambda = p_kappa, so the span was never the limit; "
          "ground-state pathology (4.28x chem vs 1.08x for 2^1S). See CHANGELOG v5.10.12.\n")
a2 = "- **/qa 58/59/60 FULL + group1 DELTA = FAIL, remediated (2026-09-07, v5.10.10):**"
assert a2 in cm and "scale lock, not l_max" not in cm
cm = cm.replace(a2, BULLET + a2, 1)
cm = cm.replace("**Version:** v5.10.11 (September 8, 2026)",
                "**Version:** v5.10.12 (September 8, 2026)", 1)
io.open("CLAUDE.md", "w", encoding="utf-8").write(cm)
print("CLAUDE.md: Sec.2 one-liner added; version cursor -> v5.10.12")
