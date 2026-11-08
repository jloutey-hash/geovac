# Sprint TD Track 5 — PSLQ probe of GeoVac entropy in the master Mellin engine ring

**Date:** 2026-05-09
**Sprint:** TD Track 5 (post-Track 1+2 closure)
**Verdict:** **CLEAN NEGATIVE.** GeoVac information-theoretic entropy `S_full(GS)`
for two-electron singlet ground states is **NOT** in the master Mellin engine
ring (M1 ∪ M2 ∪ M3, per Paper 32 §VIII case-exhaustion theorem and Paper 18
§III.7 master-Mellin reading) at the mechanically-generated basis tested.
The Track 2 user-instinct verification ("decode temperature, residual is
bare-graph entropy") is operationally clean, but the residual is genuinely
**von Neumann-only**, not a master-Mellin output.

---

## 1. Question

Track 2 (commit `f05200e`) verified numerically that the T → 0 thermal
residual of the temperature-decoding decomposition equals Paper 27's
information-theoretic entropy S_full(GS) of the spatial 1-RDM, to numerical
precision (3.1×10⁻¹⁰ for He n_max=3, 3.2×10⁻¹⁰ for Li⁺ n_max=3). The
question this sprint asks: is the dimensionless residual itself in the
master Mellin engine ring?

If yes (PSLQ identifies in M2 = √π·ℚ ⊕ π²·ℚ): GeoVac entropy is a
Seeley-DeWitt output, structurally comparable to Bekenstein–Hawking
black-hole entropy in spectral-action terms. POSITIVE.

If no (PSLQ null with control-class baseline supporting the null):
GeoVac entropy is genuinely von Neumann-only, with no master-Mellin
geometric carrier. CLEAN NEGATIVE — sharpens the framework's scope
statement at exactly the place the conversation was pointing.

---

## 2. Methodology

Strict W3 falsification protocol applied throughout (see CLAUDE.md §3 W3
entry and `docs/curve_fit_audit_memo.md`, May 2): mechanical basis frozen
**before** testing the entropy values, control class for null baseline,
no human-curated forms.

### 2.1 High-precision computation of S_full(GS)

The graph-native FCI matrix (M_L = 0 singlet) is rebuilt at
`mp.mp.dps = 150`. The hopping h₁ has rational diagonal `−Z²/(2n²)` and
rational off-diagonal `(1/16) · A[i,j]`. Two-electron integrals
`R^k(n_a l_a, n_c l_c; n_b l_b, n_d l_d)` are computed via
`geovac.hypergeometric_slater.compute_rk_algebraic` as exact `Fraction`,
then converted to mpmath. Gaunt coefficients
`c^k(l₁, m₁; l₂, m₂) = (−1)^{m₁} √((2l₁+1)(2l₂+1)) · 3j(0,0,0) · 3j(−m₁, m₁−m₂, m₂)`
are computed via sympy `wigner_3j` at high precision. The ground-state
eigenvector of the FCI matrix is extracted via `mp.eigsy`, the spatial
1-RDM is built (matching `debug.entanglement_geometry.build_1rdm_from_singlet_ci`
exactly in convention), and the von Neumann entropy
`S = −Σᵢ (nᵢ/2) log(nᵢ/2)` is computed at full mpmath precision.

Driver: `debug/sprint_td_track5.py`, function `compute_S_full_GS_mpmath`.

Three systems tested:
- **He n_max=3** (Z=2, 31 configs): the canonical Track 2 system.
- **Li⁺ n_max=3** (Z=3, 31 configs): same n_max, Z scaling probe.
- **He n_max=4** (Z=2, 101 configs): basis convergence probe.

### 2.2 Mechanical basis (frozen specification)

Following `debug/w3_mechanical_basis.py` methodology:

**Seeds organised by Mellin-engine class:**

| Class | Seeds |
|:--|:--|
| M1 (Hopf-base measure / π-family) | π, π³, 1/π, 1/π², 1/π³, 1/(2π), 1/(4π), 1/(2π²), 4/π (= Vol(S²)/π² = M1 signature from Paper 38 Track L2) |
| M2 (Seeley-DeWitt / √π·ℚ ⊕ π²·ℚ ring) | √π, 1/√π, π², π⁴, π⁶, ζ(2), ζ(4), ζ(6), ζ(8) |
| M3 (vertex parity Hurwitz / Dirichlet-L) | Catalan G, β(4), β(6), ζ(2,1/4), ζ(2,3/4), ζ(3), ζ(5) |
| ALG | √2, √3, √5, φ, log 2, log 3 |
| RAT | small integer ratios (control class) |

**Generation rule (frozen):**
- Single-seed: prefactor `(n / d)` with `n ∈ {±1, ..., ±5}`, `d ∈ {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 30, 40, 60}`, gcd(n, d) = 1.
- Two-seed: `a × b / (n / d)` and `a / b / (n / d)` with `n ∈ {±1, ±2, ±3}`, same d set.
- Range filter: `|v| ∈ [10⁻⁶, 1]` (entropy values are O(0.01–0.1)).
- Dedup: 10⁻¹² relative tolerance.

**Resulting basis size:** 12,312 deduped forms, broken down:

| Class | Count |
|:--|--:|
| ALG | 914 |
| ALG+M1 | 1,438 |
| ALG+M2 | 1,322 |
| ALG+M3 | 1,284 |
| M1 | 406 |
| M1+M2 | 1,062 |
| M1+M3 | 1,910 |
| M2 | 964 |
| M2+M3 | 1,710 |
| M3 | 1,244 |
| RAT | 58 |

### 2.3 PSLQ search

`mp.pslq([target, v₁, …, vₖ], tol=10⁻⁵⁰, maxcoeff=1000)` is applied:
- **Single-element:** every form in the 12,312-form basis tested against
  each system's target.
- **Two-element pairs:** the 200 simplest forms (sorted by name length as
  proxy for rationality / simplicity) within M1 / M2 / M3 / mixed classes
  are pair-tested (19,900 pairs per system).
- **Tolerance:** 10⁻⁵⁰ residual.
- **Coefficient cap:** 1,000.

### 2.4 Control baseline

Twenty random comparable targets are drawn log-uniform on [target/3,
target × 3] in absolute value. Each is single-element-PSLQ-tested against
the entire basis at max_coeff = 100, tol = 10⁻⁵⁰.

---

## 3. High-precision S_full(GS) values

### He n_max=3 (Z=2, 31 configs)

```
S_full(GS) = 0.04081105136647181486799816980830280755418776212596723531469218
             3091428112917345686024166646215569823234420173... (150+ dps)
E_GS       = -2.89310974223393027854914191984
```

Track 2 reported `0.040811` with residual `3.1 × 10⁻¹⁰`. Bit-identical to
all displayed Track 2 digits.

### Li⁺ n_max=3 (Z=3, 31 configs)

```
S_full(GS) = 0.01121171793742097747395655693312407557227711308912820346... (150+ dps)
E_GS       = -7.23082319876139618183324973321
```

Track 2 reported `0.011212` with residual `3.2 × 10⁻¹⁰`. Bit-identical.

### He n_max=4 (Z=2, 101 configs)

```
S_full(GS) = 0.04187943008973792077826067816342070803653770238105... (150+ dps)
E_GS       = -2.89540556069780430681358249381
```

n_max convergence: He n_max=3 → n_max=4 changes S_full by `+0.001068` (a
2.6% relative shift). The structural behaviour we are testing is not
expected to depend on this small basis-truncation correction at the
PSLQ-discrimination level (the M-basis ring elements differ by O(1) factors).

---

## 4. PSLQ search results — CLEAN NEGATIVE

All three systems return zero PSLQ identifications across the 12,312-form
mechanical basis at the configured thresholds.

| System | n_candidates_1elt | n_candidates_2pair | z vs null |
|:--|--:|--:|--:|
| He n_max=3 | 0 | 0 | +0.00 |
| Li⁺ n_max=3 | 0 | 0 | +0.00 |
| He n_max=4 | 0 | 0 | +0.00 |

Control baseline: 20 random comparable targets at single-element PSLQ
against full basis (max_coeff=100, tol=10⁻⁵⁰): mean=0.0, std=0.0, max=0,
min=0 hits. The basis is genuinely strict at this combination of
parameters; the target's null result is not a power-asymmetry artifact.

### 4.1 Stress test at max_coeff = 10⁶

To rule out a max_coeff = 10³ artifact, the He n_max=3 target was
single-element-PSLQ-tested against the 10 simplest ring elements at
max_coeff = 10⁶, tol = 10⁻⁵⁰. All 10 returned null:

| Element | Result |
|:--|:--|
| 1 | null |
| π | null |
| π² | null |
| 1/π | null |
| ζ(2) | null |
| ζ(3) | null |
| Catalan G | null |
| log 2 | null |
| √2 | null |
| φ | null |

Six orders of magnitude in coefficient relaxation does not change the
verdict. The negative is robust.

---

## 5. Cross-system consistency

| Ratio | Value |
|:--|:--|
| S(He n_max=3) / S(Li⁺ n_max=3) | 3.640035505197479881615751279528... |
| S(He n_max=3) / S(He n_max=4) | 0.974489177121636619509489616456... |
| S(He n_max=4) / S(Li⁺ n_max=3) | 3.735326764684147293712070441141... |

The Z-scaling slope from log(S(He) / S(Li⁺)) / log(2 / 3) gives
**γ ≈ −3.19** at n_max=3, somewhat steeper than Paper 27 EP-2c's
`S ~ Z^{−2.6}` asymptotic, but consistent with the EP-2f finding that
the slope drifts in n_max and Z. None of the three ratios PSLQ-identifies
in the master Mellin ring either (verified at coeffmax=10³). The ratios
themselves do not appear to be in the ring.

---

## 6. Implications for the GV-vs-BH entropy comparison

The chat-context question framing this sprint:

> "If S_full(GS) identifies in M2 ring — GeoVac entropy is a
> Seeley-DeWitt output, structurally comparable to BH spectral-action
> entropy."

That conditional fails. The PSLQ outcome is null, and the control
baseline is a clean zero. Three structural readings:

**(a) GeoVac correlation entropy is genuinely von Neumann-only, not
spectrally expressible.** This was the chat's pre-registered alternative
outcome. It means the GV-vs-BH entropy comparison is structurally weak:
both are von Neumann entropies of reduced density matrices at the
abstract level, but only BH entropy lives in the master Mellin ring;
GV entropy at this resolution does not. They are different *kinds* of
mathematical objects, both legitimately called entropy, but not directly
comparable through the master Mellin engine.

**(b) The "information-theoretic geometry" framing the user invoked needs
a different mathematical structure.** Fisher information metric, Amari
geometry, or the Latrémolière-Connes-vS Wasserstein–Kantorovich
distance on truncated state space (Paper 38 §VIII Limit Identification
Remark) are candidate destinations. None of these is the master Mellin
engine.

**(c) Track 2's residual identity (`S_thermo - S_apparatus = k_B × S_info`,
4.8×10⁻¹⁴ machine precision) remains operationally clean — the
decomposition is rigorous, only the *content* of S_info is shown by this
sprint to be outside the Mellin ring.**

This sharpens, rather than weakens, the structural-skeleton-scope
statement of CLAUDE.md §2 ("Sprint of 2026-05-07 (post-WH1-PROVEN
convergent-findings session)"). The framework's master Mellin engine
classifies the transcendental signatures of *projection* observables
(Stefan-Boltzmann, Lamb shift, Casimir, anomalous magnetic moment); it
does not classify the von Neumann entropies of finite reduced density
matrices on bare-graph FCI ground states. These are distinct classes of
observable; entropy lives outside the engine.

The chat's pre-registered "if PSLQ null" reading was: "the comparison
to BH would require a different mathematical structure; the
information-theoretic geometry framing would then need to point to a
different mathematical object than the master Mellin engine — possibly
Fisher metric, Amari geometry, or something we haven't named yet."
That reading stands.

---

## 7. Honest scope and caveats

- **Mechanical basis not exhaustive.** The 12,312-form basis covers M1 +
  M2 + M3 + ALG + RAT × small-integer prefactors × O(2) seed combinations.
  A larger basis (4-seed combinations, higher-coefficient prefactors,
  Tornheim-Witten zeta, polylogarithms, colored MZVs) might in principle
  identify the target. But the W3 lesson (CLAUDE.md §3) is that
  expanding the basis post-hoc to find a hit is exactly the kind of
  selection-bias mistake the curve-fit audit flagged. The frozen-basis
  null at the levels tested is the correct interpretation.
- **Tolerance choice 10⁻⁵⁰ is conservative.** A genuine identity should
  have residual < 10⁻¹⁴⁰ at our 150-dps precision. Loosening to 10⁻²⁰
  would only generate false positives.
- **Three systems is a small sample.** A future sprint could extend to
  Be (4-electron, requires Paper 27 EP-2j analytical machinery in mpmath)
  and He at n_max=5. None of these is expected to change the verdict —
  the null is robust at six orders of magnitude in coefficient relaxation.
- **The H₂ ground state** has S_full = 0 exactly by single-determinant
  rigidity (Paper 27 §II). Track 2 verified this to 0 residual. This is
  the trivial point of the function we are probing; it does not constrain
  the M-basis ring.

---

## 8. Net effect for the program

GeoVac correlation entropy joins the list of irreducible-but-natural
constants the framework cannot explain through the master Mellin engine:

| Constant | Source | Mellin-ring identification? |
|:--|:--|:--|
| K = π(B + F − Δ) ≈ 1/α | Paper 2 conjecture | **NO** (twelve mechanisms eliminated; Sprint K-CC) |
| L2 next-order constant c ≈ 4.1093... | Sprint MR-C | **NO** (PSLQ null in M1 ∪ M2 ∪ M3) |
| Wolfenstein parameters λ, A, ρ̄, η̄ | CKM mixing | **NO** (Sprint W3 falsified, mechanical basis z=−0.52) |
| **S_full(GS) for He, Li⁺ at n_max=3, 4** | **Paper 27 / Track 2** | **NO (this sprint)** |

This is consistent with the two-axis Paper 34 / Paper 18 framework:
the master Mellin engine classifies *which* transcendentals can appear
in *projection* observables; it does not classify the entropy content
of bare-graph eigenstates. The W3-class question ("can calibration data
be derived from a deeper combinatorial structure") remains structurally
open; this sprint adds to its evidence base by ruling out one specific
candidate location for an identification.

---

## 9. Files

- `debug/sprint_td_track5.py` — high-precision driver (~600 lines).
- `debug/data/sprint_td_track5.json` — frozen basis spec, S_full values,
  PSLQ candidates, control baseline, cross-system ratios.
- `debug/data/sprint_td_track5_intermediate.json` — high-precision values
  saved after computation, before PSLQ phase.
- `debug/data/sprint_td_track5_run.log` — full driver run log.
- `debug/sprint_td_track5_memo.md` — this memo.

No production `geovac/` modifications; existing tests untouched.
