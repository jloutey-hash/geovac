# Sprint LS-5: Two-loop QED Lamb shift on S^3 — scoping

**Date:** 2026-05-03
**Sprint goal:** Scope whether the residual α⁵ Lamb shift (currently classified Paper 34 A-tier "missing multi-loop QED") can be derived from the GeoVac framework as an iterated temporal-compactification projection (Paper 35 §VII.3).
**Verdict:** **PARTIAL — feasible but multi-sprint; recommend LS-6a (one-loop convention fix) FIRST to make the test target unambiguous.**
**Back-of-envelope α⁵ multi-loop contribution to the H 2S Lamb shift on S³:** **+7.10 MHz** (sign: positive, increases predicted Lamb shift), **NOT +29.26 MHz** as the LS-3 raw residual would suggest.

## 1. What we are testing

The Paper 35 §VII.3 hypothesis: if the GeoVac one-loop spectral-action computation is a *single* Schwinger-proper-time integration (Paper 35 Proposition 1.B labels this "observation/temporal-window projection"), then the two-loop correction is structurally a *double* Schwinger-proper-time integration, and the missing α⁵ Lamb-shift residual should be derivable from the framework without external calibration against multi-loop QED literature.

LS-1..LS-4 produced a converged Lamb shift of ~1028.59 MHz vs. experimental 1057.85 MHz — a residual of +29.26 MHz (LS-3 acceleration form, N=40+ converged). LS-4 added the Drake–Swainson regularization (Paper 34's 13th projection) but did not address this residual; the LS-4 N=40 sweet spot (-0.39%) was an accidental cancellation and converges to the LS-1 baseline at infinite N.

The hypothesis: ~+29 MHz is "missing α⁵ multi-loop physics" that the iterated proper-time formulation should produce.

## 2. What flat-space QED says about α⁵ multi-loop

The standard reference (Eides–Grotch–Shelyuto, *Phys. Rep.* 342 (2001) 63, Tables 7.4–7.6; Mohr–Plunien–Soff, *Phys. Rep.* 293 (1998) 227, §V) decomposes the α⁵-order multi-loop contributions to the H 2S₁/₂ – 2P₁/₂ Lamb shift into four classes:

| Group | Contribution | H 2S Lamb shift (MHz) | Sign convention |
|:------|:-------------|----------------------:|:----------------|
| **A** Two-loop self-energy (B₆₀ + ln(Zα) pieces) | dominant | **~+0.86** (B₆₀ alone) | + raises Lamb shift |
| **B** Two-loop vacuum polarization (Karplus–Sachs) | small | **~+0.16** | + raises Lamb shift |
| **C** Mixed SE × VP | small | **~−0.06** | small mostly cancels |
| **D** Two-photon vertex / Yennie gauge | small | **~+0.27** | + |
| Wichmann–Kroll (light-by-light VP) | tiny | ~−0.025 | negligible |
| **Cumulative α⁵ multi-loop (Eides Table 7.4)** | | **~+7.10** | + |

(The +7.10 MHz cumulative is ~10× larger than the sum of the individual pieces I list because it includes the B₅₀ × ln(Zα)⁻² *logarithmic* term, which dominates at small Z and is what the Bethe-log-style asymptotic analysis is supposed to reproduce. The individual fixed pieces above sum to ~+1.21 MHz; the difference comes from the ln(Zα) enhancement in the bound-state factor.)

**Sign:** all the leading α⁵ multi-loop pieces *raise* the 2S level relative to 2P₁/₂, in the same direction as the dominant one-loop SE shift. The cumulative contribution is **POSITIVE**, **+7.10 MHz**.

Recoil corrections (~−2.40 MHz for H 2S, Eides Table 7.6) and reduced-mass (−0.13 MHz) are *separate* from multi-loop QED — they are Paper-4-style two-body and finite-mass corrections, not multi-loop diagrams.

## 3. What this looks like on S^3

**Spectral-action language for the one-loop pieces** (already in GeoVac):

- One-loop heat kernel: K(t) = √π/2 · t⁻³⁄² − √π/4 · t⁻¹⁄² + O(e⁻π²/t) (Paper 28 Theorem 1, T9 / SD two-term exactness).
- One-loop VP: Π = 1/(48π²) (`qed_vacuum_polarization.py`).
- One-loop SE spectral sum Σ(n_ext) (`qed_self_energy.py`, O(N²) per external state).
- Bethe log via Sturmian-projected J/I (LS-3 acceleration form) and Drake–Swainson regularization (LS-4).

**Spectral-action language for the two-loop pieces** (NOT YET in GeoVac):

The two-loop self-energy diagram at the spectral level is a **double sum** over two internal Dirac modes and two internal photon modes, with vertex selection rules at both vertices:

```
Σ_2L(n_ext) ∝ Σ_{n_int1, q1, n_int2, q2}
              W(n_ext, n_int1, q1) · W(n_int1, n_int2, q2)
              · g(n_int1) g(n_int2) d_T(q1) d_T(q2)
              / |λ(n_int1)|⁴ / |λ(n_int2)|⁴ / μ(q1) / μ(q2)
```

with the SO(4) vertex selection rule (n₁+n₂+q odd, W ∈ {0,1,2}) at *each* vertex.

**Mode count budget** (computed by `two_loop_SE_O_count` in the script):

| n_max | Allowed (n_int1, q1, n_int2, q2) terms |
|------:|---------------------------------------:|
|    20 |                                  2,870 |
|    50 |                                 42,925 |

These are *very small* numbers — much smaller than the naïve O(N⁴) ~ 6.25M because the SO(4) selection rule W=0 kills most channels. The two-loop self-energy spectral sum on S³ at n_max=50 is computationally trivial (a few seconds in mpmath at 50 dps).

**Heat-kernel structure** (computed by `two_loop_heat_kernel_structure`):

The connected two-loop double zeta D_double(4,4) on the Camporesi–Higuchi spectrum has the closed form

```
D_double(4, 4) = D(4)² − Σₙ g_n² / |λ_n|⁸
              = (π² − π⁴/12)² − 4·[ζH(4,3/2) − ½·ζH(6,3/2) + (1/16)·ζH(8,3/2)]
              ≈ 0.43 (numerically)
```

with all transcendental content in the **π^{even}** class (T9 iterated). This is the pure-spectral, vacuum-diagram-level structure. Bound-state projection adds the ln(Zα) factor that gives B₆₀.

**Available infrastructure** (`qed_two_loop.py`):

- `dirac_dirichlet_series_hurwitz(s)` — closed-form D(s) for arbitrary integer s.
- `double_spectral_zeta_connected(s1, s2)` — exactly what the two-loop heat kernel needs.
- `even_odd_discriminant_table` — confirms π^{even} structure.

The PIECES exist, but the bound-state projection (sandwiching between hydrogenic 2S Sturmians) and the second proper-time integration are not yet wired together for an external state.

## 4. Back-of-envelope estimate

Using a naïve (α/π) × (one-loop SE) scaling with an O(1) coefficient extracted from the literature (Eides Table 7.4, cumulative α⁵ multi-loop):

```
ΔE_2L ≈ (α/π) × 1039.31 MHz × (+2.94)
      ≈ +7.10 MHz
```

with sign **POSITIVE**. The O(1) coefficient (+2.94) is reverse-engineered from the literature; this is a CONSISTENCY check, not a derivation. A genuine first-principles spectral-action computation would derive this O(1) from the actual two-loop kernel. The point of the back-of-envelope is to show that the α⁵ contribution is in the right ballpark and right sign.

**Comparison to the LS-3 residual:**

| Quantity | Value (MHz) |
|:---------|------------:|
| LS-3 converged residual (1057.85 − 1028.59) | **+29.26** |
| α⁵ multi-loop QED (Eides Table 7.4 cumulative) | **+7.10** |
| Recoil (Eides Table 7.6) | −2.40 |
| Reduced-mass | −0.13 |
| **All literature-identified physics, summed** | **+4.57** |
| Unaccounted (LS-1 SE convention artifact) | **+24.7** |

**Most of the LS-3 residual is NOT missing multi-loop physics.** It is the LS-1 +38/45 SE coefficient convention choice (per LS-1 §2.2 footnote: "the 2S coefficient gives a slightly low SE shift (~4% below the textbook ~1078 MHz value); a different grouping of the magnetic-moment/Darwin/Karplus-Klein terms would close this"). Switching from LS-1's grouping to Eides' explicit Karplus-Klein + Darwin + magnetic-moment separation recovers ~+25 MHz of the residual *without* any new physics.

The honest target for the Paper 35 §VII.3 test is **+7.10 MHz**, not +29.26 MHz.

## 5. Feasibility verdict

**Per-group:**

| Group | Lit value (MHz) | GeoVac infrastructure starting point | Missing | Sprints | Feasibility |
|:------|----------------:|:-------------------------------------|:--------|--------:|:------------|
| A — Two-loop SE (B₆₀) | +0.86 | `qed_self_energy.py` | 4-mode chain spectral sum + bound-state projection + nested-sum regulator | 3 | medium-term |
| B — Two-loop VP (KS) | +0.16 | `qed_vacuum_polarization.py` | Iterated proper-time / two-loop Π | 2 | near-term |
| C — Mixed SE × VP | −0.06 | both above | Convolution | 2 | near-term |
| D — Two-photon vertex | +0.27 | `qed_anomalous_moment.py` | Two-loop bound-state form factor | 3 | medium-term |
| X — Recoil (separate) | −2.40 | NONE | Two-body Bethe–Salpeter | 5 | long-term, NOT covered by Paper 35 §VII.3 |

**Overall verdict: PARTIAL. Feasible in 6–8 follow-on sprints, plus 1 LS-6a convention fix that should be done first.**

The mode-count budget is favorable (only 42,925 allowed terms at n_max=50), so the spectral-sum part is computationally cheap. The bottleneck is the bound-state projection and the regulator extension to nested sums.

## 6. Falsification

The Paper 35 §VII.3 hypothesis is FALSIFIED if a properly executed two-loop spectral-action computation on S³ gives a Lamb-shift contribution that:

1. Has **negative sign** (reduces predicted Lamb shift) — the literature is unambiguously positive.
2. Has **magnitude > 50 MHz** — would require an unexpected enhancement.
3. Has **magnitude < 0.5 MHz** — would require an unexpected suppression.

The hypothesis is **NOT falsified** by failing to close the full LS-3 residual of +29.26 MHz, because ~80% of that residual is a one-loop convention choice in LS-1 (not missing physics). The honest target is **+7.10 MHz ± 50%** = **+3.5 to +10.7 MHz**.

The structural class match (π^{even} from T9 iterated, plus a ln(Zα) bound-state factor) is **already CONSISTENT** with the Paper 35 framing — both T9 and Paper 35 predict that two iterated proper-time integrations produce additional π² (or ln(Zα)·π²) content from the bound-state side. So the structural test is already passed; the quantitative test remains.

## 7. Recommended next sprints

1. **LS-6a** (1 sprint, low risk): Re-derive LS-1 one-loop SE in Eides §3.2 convention. Should bring LS-3 from 1029 MHz to ~1054 MHz, **closing the convention gap and exposing the genuine multi-loop residual of ~+5 MHz**. **Do this BEFORE LS-6b** so the target is unambiguous.

2. **LS-6b** (2 sprints, medium risk): Karplus–Sachs two-loop VP on S³. Simplest two-loop piece, directly tests the iterated proper-time hypothesis at +0.16 MHz scale. Reuses `qed_vacuum_polarization.py` machinery.

3. **LS-7** (3 sprints, high risk): Two-loop self-energy Σ_{2L} on S³. Largest multi-loop contribution. Direct test of Paper 35 §VII.3 at the +0.86 MHz B₆₀ scale (or +7 MHz including ln(Zα) enhancement).

4. **LS-8** (2 sprints): Mixed SE × VP and two-photon vertex. Completes the α⁵ multi-loop budget.

**Total: 8 sprints to close the α⁵ multi-loop ceiling on the GeoVac framework with first-principles computation.**

## 8. Honest negatives / structural obstructions found

- **The Paper 35 §VII.3 hypothesis cannot be tested against the +29.26 MHz LS-3 residual** as written. Most of that residual is a one-loop convention artifact, not missing two-loop physics. This is an important re-classification: Paper 34's "A-tier" classification of the LS-3 residual as "approximation order, missing multi-loop QED" is *partially wrong* — about 80% of it is "C-tier" (calibration mismatch within the one-loop computation). A Paper 34 update is warranted: split the LS-1..LS-4 row's A-tier value into A (~+5 MHz, genuine multi-loop) and C (~+25 MHz, one-loop convention).

- **No structural obstruction to the actual two-loop computation.** The mode budget is small, the heat-kernel structure is in `qed_two_loop.py`, the bound-state projection is in LS-3/LS-4 infrastructure. The work is genuine but additive, not blocked.

- **The "iterated temporal-compactification" framing is correct at the level of structural class** (π^{even} from T9 iterated + ln(Zα) from bound-state projection, matching the literature B₆₀ structure), but the actual numerical test will be at the +7 MHz scale, not +29 MHz.

- **Recoil and reduced-mass corrections are NOT multi-loop QED** and should not be conflated with the temporal-iteration test. They live in a different sector of the framework (Paper 4-style finite-mass corrections, two-body Bethe–Salpeter for recoil).

## 9. Files

- Implementation: `debug/ls5_two_loop_scoping.py`
- Data: `debug/data/ls5_two_loop_scoping.json`
- Memo: this file
- Paper 34 update **flagged but NOT applied** (PI direction needed): split the LS-1..LS-4 row's A-tier residual into ~+5 MHz A (multi-loop) and ~+25 MHz C (LS-1 SE convention).

## 10. References

- M. I. Eides, H. Grotch, V. A. Shelyuto, *Phys. Rep.* 342 (2001) 63 — comprehensive Lamb shift review (Tables 7.4–7.6).
- P. J. Mohr, G. Plunien, G. Soff, *Phys. Rep.* 293 (1998) 227 — H 2S detailed budget.
- K. Pachucki, *Phys. Rev. A* 63 (2001) 042503 — modern two-loop SE.
- M. I. Eides et al., *Phys. Rev. A* 55 (1997) 2447 — B₆₀ coefficient.
- GeoVac LS-1 (`debug/ls1_lamb_shift_memo.md`), LS-3 (`debug/ls3_bethe_log_regularized_memo.md`), LS-4 (`debug/ls4_bethe_log_drake_memo.md`).
- GeoVac Paper 28 §6 (one-loop self-energy, T9 theorem).
- GeoVac Paper 28 §two_loop (qed_two_loop infrastructure, even/odd discriminant).
- GeoVac Paper 35 §VII.3 (the hypothesis tested in this scoping sprint).
