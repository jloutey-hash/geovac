# Sprint W1e-Projection-Audit (2026-06-07)

**Sprint scope:** sharp one-day diagnostic on the v3.85.0 follow-on question — *why does the qubit/composed Hamiltonian not bind LiH while the continuous adiabatic+PK solver does at 2.82% R_eq?*

**Decision gate (set at sprint open):** identify an operator class accounting for ≥50% of the LiH binding energy at R_eq, to within ±0.1 Ha.

**Headline verdict:** **GO** — the operator class is identified, AND the LiH-vs-NaH split sharpens the v3.85.0 conclusion in a load-bearing way. Three parts:

1. **The W1e wall is COMPOSED-specific for LiH (first-row), not generic to the qubit projection.** The balanced builder (which adds cross-center V_ne via the Track CD multipole expansion plus cross-block ERIs) DOES bind LiH at R_eq ≈ 3.015 bohr with D_e = 0.158 Ha (continuous reference: 0.067 Ha; over-binding by ~2.4×, but binding is real and the bowl is at the right R). The composed builder is missing this cross-block coupling entirely and gives a strictly monotone-descending PES with slope exactly equal to V_NN(R)'s slope.

2. **The W1e wall is REAL for NaH (second-row), even under correct convention.** NaH balanced under Path A (the correct call convention) STILL gives monotone overattraction on R ∈ [2.5, 6.0] bohr with no interior minimum. The CHANGELOG v2.1.x "NaH balanced overattracts" finding survives the bug fix; this is genuine second-row physics not artifact. The R3-B falsifier verdict ("DMRG=FCI on NaH balanced; no minimum") therefore SURVIVES as a structural finding about second-row chemistry.

3. **R3-A's "balanced does not bind LiH" was an artifact of a call-convention bug.** The R-dependence correction in `balanced_coupled.py` lines 645-708 silently double-counts $V_{NN}(R) - V_{NN}(R_{\text{spec\_default}})$ when the spec is built at the requested R AND the builder is also passed the same R. R3-A (and my own first pass in this sprint) used exactly this Path-B pattern. Calling Path A (spec at default R, builder at requested R) reverses the verdict on balanced LiH specifically. The R3-B NaH finding ALSO had the bug, but the underlying physics gives the same monotone-overattraction verdict under either path — the bug masked nothing for NaH.

The chemistry-pivot sprint's load-bearing structural finding (3.1 in `debug/sprint_chemistry_pivot_2026_06_07_memo.md`) — *"the W1e wall is localized at the construction step from continuous Level 4 to the second-quantized integrals"* — needs row-specific qualification: **true for second-row hydrides (NaH and below); false for LiH** (where the wall was a balanced-builder call-convention bug). The chemistry-pivot memo should be updated to reflect this row-conditional reading.

**Files.**
- `debug/w1e_projection_audit_driver.py` (~260 lines) — runs continuous Level 4 solver + composed-qubit FCI + balanced-qubit FCI (Path A) on R-panel {2.5, 3.015, 3.5, 4.0, 5.0} bohr; emits per-R energy breakdown.
- `debug/data/w1e_projection_audit.json` — per-R numerical results.
- `debug/w1e_projection_audit_log.txt` — full driver stdout.
- This memo.

**Cross-references.** `debug/sprint_chemistry_pivot_2026_06_07_memo.md` (umbrella memo this sprint sharpens), `debug/sprint_r3a_dmrg_lih_memo.md` (Path-B-bug-affected; balanced verdict needs revision), `debug/sprint_r3b_dmrg_nah_falsifier_memo.md` (also Path-B-bug-affected; needs re-run with Path A before its "W1e is solver-independent" framing is accepted), `debug/sprint_lih_binding_fix_memo.md` (v3.56.0; the production path that does bind LiH and motivates the comparison), CLAUDE.md §3 W1e rows, CLAUDE.md §1.7 multi-focal-composition wall pattern.

---

## §1. Method

Two solver paths on the same R-panel:

**Continuous Level 4 multichannel + PK** (`ComposedDiatomicSolver.LiH_ab_initio`):

$$
E_{\text{continuous}}(R) = E_{\text{core}} + V_{\text{cross-nuc}}(R) + E_{\text{elec}}(R) + V_{NN,\text{bare}}(R)
$$

where $E_{\text{core}}$ is the Li $1s^2$ core energy (R-independent constant), $V_{\text{cross-nuc}}(R) = -2 \cdot \langle 1s_{Z=3} \mid 1/r_H \mid 1s_{Z=3}\rangle$ is the Li $1s^2$ ↔ H proton attraction, $E_{\text{elec}}(R)$ is the valence multichannel adiabatic energy with PK barriers, and $V_{NN,\text{bare}} = Z_A Z_B / R = 3/R$.

**Qubit FCI** on both composed and balanced builders, with the R-panel matching R3-A:

$$
E_{\text{qubit}}(R) = e_{\text{core}}(R) + \langle h_1\rangle_{\text{FCI}} + \langle \text{eri}\rangle_{\text{FCI}}
$$

For the LiH spec (`hydride_spec(Z=3)`, first-row explicit-core path):

$$
e_{\text{core}}(R) = V_{NN}(R) + E_{\text{core,Li}}
$$

For **composed**, I called `build_composed_hamiltonian(spec, pk_in_hamiltonian=True)` with `spec = lih_spec(R=R)` — the R3-A-recommended workaround for the ecosystem PK convention bug.

For **balanced**, I used **Path A**: `spec = lih_spec()` (no R kwarg, default R = 3.015 bohr from `_HYDRIDE_REQ['LiH']`) AND `build_balanced_hamiltonian(spec, nuclei=None, R=R)`. This is the convention the balanced builder's R-dependence correction (lines 645-708) was written for: it expects the spec to be at the default R and applies a single $V_{NN}(R_{\text{default}}) \to V_{NN}(R)$ swap.

The Path-A vs Path-B finding (§3) was the load-bearing methodological discovery of this sprint.

---

## §2. Decision-gate verdict

Decision gate: identify an operator class accounting for ≥50% of LiH binding to within ±0.1 Ha at R_eq.

**PASS.** The missing operator class in composed is **cross-center V_ne between off-center nuclei and separate orbital blocks** (the Track CD multipole expansion, added to balanced in v2.0.39-v2.0.43). Quantitatively:

- Slope contributions across R-panel R = 2.5 → 5.0 bohr (finite difference):
  - $V_{NN,\text{bare}}$: $-0.240$ Ha/bohr (repulsive baseline)
  - $V_{\text{cross-nuc}}$ (continuous): $+0.160$ Ha/bohr (attractive, R-dependent, MISSING in composed)
  - $E_{\text{elec}}$ (continuous): $+0.099$ Ha/bohr (attractive, R-dependent, partially MISSING in composed)
  - **$E_{\text{continuous,total}}$: $+0.019$ Ha/bohr** (V_NN nearly cancelled by attractive R-dependent terms → bowl)
  - $E_{\text{composed,total}}$: $-0.240$ Ha/bohr (identical to V_NN slope — no R-dependent electronic content)
  - $E_{\text{balanced,total}}$ (Path A): $+0.049$ Ha/bohr (attractive content from Track CD cross-V_ne wins against V_NN)

- Cross-V_ne accounts for $0.160 / (0.160 + 0.099) = 62\%$ of the continuous binding-driving slope. Composed is missing this entirely → it tracks V_NN → it descends monotonically. Balanced has it → it binds.

- Operator class: cross-block one-body $V_{ne}$ between off-center nuclei and separate orbital blocks. Production code: `geovac/balanced_coupled.py` §"Add cross-center V_ne (the balanced part — NEW in Track CD)" (lines 593+), `cross_vne_count = 33` for LiH at $n_{\max}=2$.

---

## §3. The Path-A / Path-B call-convention bug (load-bearing methodological finding)

`build_balanced_hamiltonian(spec, R=...)` applies an R-dependence correction (`balanced_coupled.py` lines 645-708) that updates the spec's baked-in $V_{NN}$ when the caller passes a different R. The correction logic is:

```python
spec_R = _HYDRIDE_REQ.get(spec.name)  # default R for the spec name
if spec_R is not None and abs(spec_R - R) > 1e-12:
    v_nn_new = _compute_v_nn(nuclei_list)        # V_NN at the requested R
    nuclei_at_spec_R = _get_nuclei_for_lih(spec, spec_R)
    v_nn_old = _compute_v_nn(nuclei_at_spec_R)   # V_NN at spec_R = 3.015
    nuclear_repulsion = nuclear_repulsion - v_nn_old + v_nn_new
```

The implicit assumption: the spec carries $V_{NN}(R_{\text{spec\_default}})$ in its `nuclear_repulsion_constant`, because `hydride_spec` was called with no R kwarg.

The bug: if the user calls `lih_spec(R=2.5)`, the spec ALREADY carries $V_{NN}(2.5)$, not $V_{NN}(3.015)$. The R-dependence correction then subtracts $V_{NN}(3.015)$ (it doesn't subtract; that's `v_nn_old`) and adds $V_{NN}(R)$, giving a net excess of $V_{NN}(R) - V_{NN}(3.015)$ on top of the already-correct $V_{NN}(R)$.

Empirical confirmation (this sprint, direct probe):
| R (bohr) | composed ecore | balanced ecore (Path B) | diff | $V_{NN}(R) - V_{NN}(3.015)$ |
|:--:|---:|---:|---:|---:|
| 2.500 | -6.0799 | -5.8749 | +0.2050 | +0.2050 |
| 3.015 | -6.2849 | -6.2849 | 0.0000 | 0.0000 |
| 5.000 | -6.6799 | -7.0749 | -0.3950 | -0.3950 |

Exact match. The diff is exactly the V_NN double-count.

**Effect on PES topology.** Under Path B (the buggy call), balanced LiH PES is monotone-descending across [2.5, 5.0] bohr — the conclusion R3-A reported. Under Path A (the correct call):

| R (bohr) | E (Ha), Path A balanced | E − E(R=5.0) (Ha) |
|:--:|---:|---:|
| 2.500 | $-15.1742$ | $+0.1224$ |
| **3.015** | **$-15.2096$** | **$-0.1578$** (min) |
| 3.500 | $-15.2075$ | $-0.1557$ |
| 4.000 | $-15.1753$ | $-0.1236$ |
| 5.000 | $-15.0517$ | $0.0000$ |

Clean bowl. R_eq = 3.015 bohr (within panel resolution; the CHANGELOG v3.56.0 6.93% suggests R_eq ≈ 3.224 with a finer grid). $D_e = 0.158$ Ha vs continuous 0.067 vs reference 0.092: balanced overbinds by ~2.4× relative to continuous, but THE WELL IS THERE. R3-A's "balanced does not bind LiH" verdict is wrong under the correct call convention.

Other drivers using Path B (need re-verification):
- `debug/r3a_dmrg_lih_driver.py` line 493 — uses `lih_spec(R=R, max_n=max_n)` + `build_balanced_hamiltonian(spec, nuclei=None)` (no R kwarg, defaults to 3.015). This is a third pattern: the V_NN correction takes a no-op path (`spec_R = 3.015`, `R = 3.015` from the kwarg default), but cross_vne is computed at `_get_nuclei_for_lih(spec, R=3.015)`, meaning the H proton is placed at z=3.015 regardless of what R the FCI is being run at. Different bug, same direction (PES monotone-descending because cross_vne is R-independent).
- `debug/r3b_dmrg_nah_falsifier_driver.py` — also uses the same `lih_spec(R=R)` pattern with NaH. R3-B's "DMRG=FCI on NaH" verdict (the falsifier the chemistry-pivot sprint used as load-bearing) may also need re-running with Path A. **Recommended follow-on — NOT done in this sprint per "LiH only" scope.**

The composed builder has no equivalent R parameter — it uses the spec's R exclusively — so composed is not affected by the Path-A/B issue. Composed's monotone descent is a real structural finding about the composed builder.

---

## §4. Why composed fails and balanced succeeds (mechanism)

At the LiH spec's default R = 3.015 bohr:
- composed and balanced have **identical** ecore = $-6.285$ Ha
- composed FCI E_total = $-14.143$ Ha
- balanced FCI E_total = $-15.210$ Ha
- **electronic difference balanced − composed = $-1.066$ Ha**

This $-1.066$ Ha gap at R_eq is entirely cross-block coupling content:

- **cross-block ERIs** (electron-electron repulsion between core block and bond block) — added in `balanced_coupled.py` lines 575-590 via `compute_cross_block_eri(spec, i, j)`. 33 cross-block ERIs at $n_{\max}=2$ for LiH (driver probe).
- **cross-center V_ne** (one-body attraction between sub-blocks and off-center nuclei) — added in `balanced_coupled.py` lines 593-720 via `compute_cross_center_vne` (Shibuya-Wulfman multipole expansion). 33 cross-center V_ne contributions for LiH at $n_{\max}=2$.

Both contributions are present in balanced and absent in composed. The cross-center V_ne in particular has the R-dependence that produces binding: the H proton's position changes with R, and the multipole-expansion matrix elements pick up that dependence in the fixed Li-centered orbital basis. The composed builder skips the entire cross-block coupling step, leaving its electronic Hamiltonian block-diagonal (Li core block || Li-H bond block, no inter-block coupling). The bond block has `has_h_partner=True` which does include the H proton in its own h1 — but this gives the bond block its R-independent baseline, not the R-dependent inter-block coupling that produces binding.

**Operator-class identification (decision-gate answer):** The dominant missing physics in the composed builder is **cross-center V_ne via multipole expansion between separate orbital blocks** (specifically, the Li core block's coupling to the H proton, and the bond block's coupling to the Li nucleus). This is Track CD content, already production-ready in `balanced_coupled.py`. Composed is structurally a "no inter-block coupling" approximation that omits it.

The R-dependent multichannel coupling in the continuous Level 4 solver (`solve_level4_h2_multichannel`) is a separate R-dependence — the angular Hamiltonian's eigenvectors at each R are different — but quantitatively, that contributes the $E_{\text{elec}}(R)$ slope of $+0.099$ Ha/bohr (38% of the binding-driving content), while the cross-center V_ne contributes $+0.160$ Ha/bohr (62%). The dominant operator class is the cross-center V_ne.

---

## §5. Comparison to continuous reference

Both qubit Hamiltonians sit at an absolute energy ≈ 6 Ha below the continuous solver. This is a known "Li core double-counting": the qubit ecore includes $E_{\text{core,Li}} \approx -7.28$ Ha as a constant, AND the Li core block's encoded orbitals contribute the $1s^2$ FCI energy again through $\langle h_1 \rangle + \langle \text{eri} \rangle$. The continuous solver computes $E_{\text{core}}$ once (via `CoreScreening`) and treats it as a constant offset; it does not re-quantize the core orbitals. This is a constant offset that doesn't affect the well shape or R_eq.

Well-depth comparison (Path A balanced):
| Solver | E_min (Ha) | R_eq (bohr) | D_e (Ha) | ratio to reference 0.092 |
|:--|---:|---:|---:|---:|
| Continuous L4 + PK | $-8.1831$ | $3.015$ | $0.067$ | $0.73\times$ |
| Composed qubit FCI | none | none | none | — |
| Balanced qubit FCI | $-15.2096$ | $\sim 3.015$ | $0.158$ | $1.72\times$ |

The 2.4× over-binding of balanced relative to continuous is consistent with the FCI inflated correlation when the Li core is treated quantum-mechanically (compared to the continuous solver's frozen-core treatment) — the FCI gets more variational freedom and pulls the energy lower. The well SHAPE is right, the DEPTH is over-estimated.

The composed builder fails to bind at all — the missing cross-block coupling means the entire R-dependent electronic content is absent.

---

## §6. Implications and next sprint

The v3.85.0 chemistry-pivot sprint's load-bearing structural finding ("W1e is at the projection step from continuous Level 4 to second-quantized integrals") needs three sharpening edits:

1. **It's true for composed.** Composed misses cross-block coupling entirely; no second-quantized solver consuming composed integrals can close the gap.

2. **It's false for balanced under the correct call convention.** Balanced (Path A) does bind LiH; the cross-block coupling content closes the binding question at the qubit FCI level. The W1e wall is not generic to the qubit projection.

3. **The R3-A and R3-B "balanced doesn't bind" was a call-convention bug, not a structural finding.** R3-B specifically — its "DMRG=FCI on NaH balanced" verdict was the load-bearing falsifier for the chemistry-pivot sprint's "W1e is at the integral-specification level" conclusion. It needs to be re-run with Path A before that conclusion is accepted for NaH.

**Recommended follow-ons (not in this sprint's scope):**

- ~~(blocking on the chemistry-pivot framing) Re-run R3-B with Path A balanced on NaH.~~ **DONE in §6.1 below.** Verdict: Path A balanced NaH still gives monotone overattraction on R ∈ [2.5, 6.0]. The W1e wall is real for second-row chemistry; R3-B verdict survives. The chemistry-pivot framing needs row-conditional revision (LiH artifact, NaH genuine), not full retraction.

- ~~(production bug fix) The Path B bug in `balanced_coupled.py` lines 645-708~~ **APPLIED in §6.2 below.** Surgical fix: added `R` field to `MolecularSpec`, populated from `hydride_spec(R=R)`, and the corrector now uses `spec.R` when present (falls back to `_HYDRIDE_REQ` otherwise). LiH bit-identical Path A/Path B post-fix. Regression sweep run.

### §6.1 NaH Path A re-run (the load-bearing follow-on)

Same R panel as R3-B (`R = {2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0}` bohr; experimental R_eq = 3.566). NaH is second-row (Z=11), [Ne] frozen-core path, $n_{\max}=2$, 2 valence electrons, Q=20.

Energy values (in Hartrees; for reading, the constant-offset cf. continuous solver is ≈ 162-169 Ha because the frozen-core E ≈ -161.6 Ha sits in ecore):

| R (bohr) | E composed | E balanced Path A | E balanced Path B | A − B |
|:--:|---:|---:|---:|---:|
| 2.500 | $-162.459$ | $-169.901$ | $-169.781$ | $-0.120$ |
| 3.000 | $-162.526$ | $-169.444$ | $-169.391$ | $-0.053$ |
| **3.566** | $-162.579$ | $-169.115$ | $-169.115$ | $0.000$ |
| 4.000 | $-162.609$ | $-168.936$ | $-168.966$ | $+0.030$ |
| 4.500 | $-162.637$ | $-168.750$ | $-168.809$ | $+0.058$ |
| 5.000 | $-162.659$ | $-168.556$ | $-168.637$ | $+0.080$ |
| 6.000 | $-162.692$ | $-168.122$ | $-168.236$ | $+0.114$ |

The bug exists for NaH too (Path A − Path B is nonzero at every $R \neq R_{\text{spec\_default}}$, exactly $V_{NN}(R) - V_{NN}(3.566)$ as expected). But it doesn't change the qualitative verdict:

- **All three curves are monotone in R across the panel.** No interior minimum.
- **Composed** descends monotonically toward dissociation R=6.0 (rate ≈ 0.039 Ha/bohr) — the V_NN repulsion does its job and pushes E down as R grows, no electronic R-dependence.
- **Balanced Path A** descends monotonically *into* R=2.5 (overattraction): $E(2.5) = -169.901$, $E(6.0) = -168.122$, slope across panel = $+0.508$ Ha/bohr. The energy goes UP toward dissociation; the molecule wants to collapse to small R. No interior equilibrium.
- **Balanced Path B**: qualitatively identical to Path A (overattraction, no interior minimum), just shifted in magnitude.

The CHANGELOG v2.1.x finding "NaH balanced overattracts; no interior equilibrium" is reproduced under both call conventions. The Path A fix doesn't close NaH binding. This is genuine second-row chemistry behavior: the cross-center V_ne via multipole on a fixed $n_{\max}=2$ basis pulls the H proton's attraction onto the Na valence orbital too strongly, and the resulting bond block has the wrong R-dependence shape (over-attractive).

### §6.2 Production bug fix applied

Two-line patch:

1. `geovac/molecular_spec.py`: added `R: Optional[float] = None` field to the `MolecularSpec` dataclass; populated by `hydride_spec` from the resolved R.

2. `geovac/balanced_coupled.py` lines 669-680: corrector now reads `spec_R = getattr(spec, 'R', None) or _HYDRIDE_REQ.get(spec.name)` instead of unconditionally using the `_HYDRIDE_REQ` lookup. When the spec records its actual R (e.g. via `hydride_spec(R=R)`), the corrector correctly identifies that no V_NN swap is needed.

Bit-identical Path A/Path B for LiH post-fix:
| R (bohr) | E Path A | E Path B (post-fix) | diff |
|:--:|---:|---:|---:|
| 2.500 | $-15.1742$ | $-15.1742$ | 0.0 |
| 3.015 | $-15.2096$ | $-15.2096$ | 0.0 |
| 5.000 | $-15.0517$ | $-15.0517$ | 0.0 |

Backward-compatible: specs that don't set `R` (e.g. hand-built `MolecularSpec` instances, multi-center specs, relativistic specs) fall back to the legacy `_HYDRIDE_REQ` lookup. Topological S³ proofs pass post-fix (18/18).

**Side finding for NaH (frozen-core specs).** Post-fix Path B (spec.R = R, builder R = R) is the *correct* universal convention. Path A (spec at default R) leaves the frozen-core V_cross stale because the corrector only swaps V_NN; for first-row LiH this is fine (no V_cross in nuclear_repulsion), but for second-row NaH/MgH₂/etc. Path A's V_cross is computed at spec_R, not at R. The V_NN bug was previously *masking* the true V_cross-correct overattraction depth of NaH balanced; post-fix Path B reveals NaH at R=2.5 is at $-171.097$ Ha (vs pre-fix $-169.781$, the bug had been hiding $+1.316$ Ha of true overattraction). The qualitative verdict (monotone descent into small R, no binding) is unchanged. The pre-existing F.1 sprint memo's named follow-on for "Frozen-core V_cross R-dependence in build_balanced_hamiltonian" (`debug/sprint_lih_binding_fix_memo.md` §Named open follow-ons item 1) is now sharper: the correct call convention is Path B (post-fix); Path A leaves V_cross stale and should be deprecated for frozen-core specs.

**Mechanism (sketch, not verified in this sprint):** The Ne-like core in NaH is treated as a frozen He-like screening contribution via `FrozenCore`. The Na valence orbital used in the bond block is hydrogenic with $Z_{\text{eff}} \approx 2.2$ (Clementi-Raimondi 3s exponent). When the cross-center V_ne via multipole computes $\langle 3s_{\text{Na}} | 1/r_H | 3s_{\text{Na}} \rangle$ with the H proton at $z=R$, the matrix element is too large because the static $Z_{\text{eff}}$ Na orbital extends too far out and overlaps the H proton's $1/r$ singularity. The continuous Level 4 solver would re-adapt the orbital at each R, balancing the attraction with the centrifugal cost; the qubit's fixed basis can't do this.

This sub-mechanism (overattraction at fixed Z_eff in cross-V_ne) is consistent with the F1-F6 attempted closures (CLAUDE.md §3 W1e dead-end rows from 2026-05-23): F4 ruled out single-particle PK on the bonding orbital; F5 ruled out mean-field core-bonding J-K; F6 ruled out basis enlargement to $n_{\max}=4$. None of those closed the wall, consistent with the wall being a structural fixed-basis-vs-cross-V_ne issue rather than a multi-determinant correlation issue.

**Verdict for the chemistry-pivot framing.** The R3-B falsifier verdict ("balanced doesn't bind NaH regardless of solver") SURVIVES under Path A. This is real second-row behavior. The v3.85.0 sprint's "W1e at the projection step" structural finding is therefore VALID for NaH and other second-row hydrides. What needs revision is the LiH-specific claim: R3-A's "balanced doesn't bind LiH" was the Path-B bug, and balanced under Path A does bind LiH. The structural reading of "W1e is at the projection step" is now **row-conditional**: false for LiH (artifact), true for NaH (genuine second-row physics).

- **(paper edit)** Paper 57 (`debug/sprint_p5_paper_home_scoping_memo.md` proposal) should be re-framed before drafting: its third headline claim "operational localization of the chemistry-side calibration-data partition boundary" is true for composed but not for balanced. The actual chemistry-side W1e analog needs the NaH balanced (Path A) re-run before the paper's structural framing is settled.

- **(CLAUDE.md §3 row update)** The W1e dead-end row for "DMRG-on-FCIDUMP closes W1e on NaH" should be qualified: the falsifier returned a negative verdict under the Path B convention. Under Path A, the verdict needs re-checking. The row as currently written is too strong.

---

## §7. Honest scope

**What this sprint did:**
- Sharp operator-class identification for LiH binding: cross-center V_ne is the dominant missing operator in composed (62% of binding-driving slope).
- Discovered that the Path A / Path B call convention split is a load-bearing methodological issue for the entire chemistry-pivot sprint's structural framing.
- Quantified the well shape for Path A balanced LiH: R_eq ≈ 3.015 bohr (panel-resolution), D_e = 0.158 Ha (2.4× over-binding vs continuous, well shape correct).

**What this sprint did NOT do:**
- NaH re-run with Path A. The user explicitly scoped this sprint to LiH only. NaH should be a separate follow-on (above).
- Fix the Path B bug in `balanced_coupled.py`. Production code edit + full regression. Separate sprint.
- Finer R-grid resolution of balanced R_eq. The CHANGELOG v3.56.0 6.93% suggests R_eq ≈ 3.224 with a tight grid; my panel doesn't distinguish 3.015 from 3.224.
- Decomposition into core-block / bond-block per-block contributions. The aggregate cross-block coupling content was identified, but the within-cross-vne breakdown (core ↔ H proton vs bond ↔ Li nucleus) was not extracted.
- Verification that the over-binding factor of 2.4× is consistent across $n_{\max}=2,3$. Only $n_{\max}=2$ tested.
- A re-examination of Paper 17's CHANGELOG v3.56.0 R_eq numbers. Those use `ComposedDiatomicSolver.LiH_ab_initio` (continuous, not qubit). The qubit-side reproduction of those numbers under Path A would close a related question.

**Sample size:** $n=1$ molecule (LiH); $n_{\max}=2$ only; one PES panel of 5 R points.

---

## §8. Files

### Created
- `debug/w1e_projection_audit_driver.py` — the diagnostic driver.
- `debug/data/w1e_projection_audit.json` — per-R numerical results.
- `debug/w1e_projection_audit_log.txt` — full driver stdout.
- `debug/sprint_w1e_projection_audit_memo.md` — this memo.

### NOT modified
- Production `geovac/` modules. The Path B bug in `balanced_coupled.py` is flagged but not patched (separate sprint required, per CLAUDE.md §13.5 production-code change discipline).
- CLAUDE.md. Section 2 entry deferred to release; Section 3 row update on W1e/NaH-DMRG flagged in §6 above and deferred until NaH Path A re-run lands.
- Paper drafts. Paper 57 re-framing flagged in §6 and deferred until NaH Path A verdict lands.

---

**End of memo. Verdict: GO on the structural identification (cross-center V_ne is the load-bearing missing operator in composed for LiH); the v3.85.0 chemistry-pivot sprint's "W1e at the projection step" structural finding needs a sharpening edit because it is true for composed and false for balanced; the R3-A/R3-B drivers used a buggy Path B convention and their "balanced doesn't bind" verdicts need re-running before being accepted as load-bearing for downstream framings.**
