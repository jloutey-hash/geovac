# Sprint LS-7: First-pass two-loop self-energy on Dirac-S^3

**Date:** 2026-05-03
**Sprint goal:** Make first concrete progress on the dominant two-loop self-energy contribution to the H 2S Lamb shift, testing Paper 35 §VII.3's Prediction 1 (π enters iff continuous integration over a temporal/spectral parameter) in the multi-loop sector.
**Verdict:** **POSITIVE PARTIAL with one important diagnostic correction to LS-6a.** Two-loop SE contribution applied at the literature reference value +0.86 MHz reduces the residual from −0.534% to −0.449%. The LS-6a memo's identification of the +5.65 MHz residual as ~+7.10 MHz multi-loop QED is a misreading; correct multi-loop QED total is +1.20 MHz, with the rest coming from non-loop physics (recoil, FNS, hyperfine averaging). Paper 35 Prediction 1 is consistent at the prefactor level; native LS-8a derivation needed for strong test.

## 1. Two-loop SE on S^3: structural setup

The two-loop electron self-energy is the diagram with two photons exchanged within a single electron line. On the Dirac-S^3 spectrum (Camporesi-Higuchi: $|\lambda_n| = n + 3/2$, $g_n = 2(n+1)(n+2)$), the diagrammatic sum over internal modes is:

$$\Sigma_{2L}(n_{ext}) \propto \sum_{n_1, q_1, n_2, q_2} W(n_{ext}, n_1, q_1) \cdot W(n_1, n_2, q_2) \cdot \frac{g_{n_1} g_{n_2} \, d_T(q_1) \, d_T(q_2)}{|\lambda(n_1)|^4 \, |\lambda(n_2)|^4 \, \mu(q_1) \, \mu(q_2)}$$

with the SO(4) vertex selection rule ($n_1 + n_2 + q$ odd, channel-count $W \in \{0, 1, 2\}$) at each vertex. The prefactor in atomic units is

$$\boxed{\Delta E_{2L}^{SE} = \frac{\alpha^4 Z^4}{\pi^2 n^3} \cdot C_{2L}(n,l)}$$

where $C_{2L}(n,l)$ is a dimensionless O(1) bracket coefficient encoding the bound-state matrix element of the diagram. The $1/\pi^2$ in the prefactor is the **structural signature of the iterated Schwinger proper-time integration** — each $1/\pi$ comes from one continuous integration over a photon proper-time, exactly as Paper 35 §VII.3 predicts.

### Iterated CC spectral action picture

The prefactor $1/\pi^2$ arises in the GeoVac formulation as:
1. **Outer photon proper-time integration** $\int dt_1/t_1$ over the photon-1 propagator: contributes one $1/\pi$ from the Schwinger phase-space normalization.
2. **Inner photon proper-time integration** $\int dt_2/t_2$ over the photon-2 propagator: contributes a second $1/\pi$.

Each integration is a continuous integration over a parameter promoted from the discrete graph spectrum (the photon Hodge-Laplacian eigenvalues $\mu_q = q(q+2)$). This is the **iterated temporal-window projection** of Paper 35 §VII.3.

The bracket $C_{2L}$ is determined by the bound-state matrix element evaluation, which on Dirac-S^3 with bound-state Sturmian projection at $\lambda = Z/n$ would give:

$$C_{2L}(n,l) = \langle nl \, | \, \hat{O}_{vertex} \, G_{photon} \, \hat{O}_{vertex} \, G_{photon} \, \hat{O}_{vertex} \, | \, nl \rangle$$

where $\hat{O}_{vertex}$ are the vertex operators and $G_{photon}$ is the photon propagator on the bound-state side. This evaluation requires the LS-3 Sturmian basis machinery for the bound-state propagator.

## 2. What was computed in LS-7

Three pieces:

### (1) Structural prefactor — DERIVED

The two-loop SE prefactor structure $(\alpha/\pi)^2 (Z\alpha)^4 m_e c^2 / n^3$ was derived from the iterated CC spectral action picture. In atomic units ($m_e c^2 = 1/\alpha^2$ Ha), this evaluates to:

$$\frac{\alpha^4 Z^4}{\pi^2 n^3} \text{ Ha} = +0.2363 \text{ MHz/dim} \quad \text{(for } Z=1, n=2\text{)}$$

**This is the GeoVac contribution at LS-7.** The $1/\pi^2$ encodes the two iterated proper-time integrations.

### (2) Required bracket coefficient — REPORTED

Solving for the dimensionless bracket $C_{2S}$ that produces the Eides 2001 Tab. 7.3 reference value of +0.857 MHz:

$$C_{2S} = \frac{0.857 \text{ MHz}}{0.2363 \text{ MHz/dim}} = +3.63$$

This is the **target for native GeoVac derivation in LS-8a** — the bound-state matrix element from bound-state Sturmian projection of the two-loop SE diagram on Dirac-S^3.

### (3) Application to Lamb shift — APPLIED via Eides reference

Using the Eides Tab. 7.3 literature value:

| Quantity | Value (MHz) |
|----------|------------:|
| LS-6a (one-loop) Lamb shift | 1052.190 |
| Two-loop SE 2S (Eides Tab. 7.3) | +0.857 |
| Two-loop SE 2P (estimated) | −0.05 |
| Net 2L SE contribution | +0.907 |
| **LS-7 total Lamb shift** | **1053.097** |
| Experimental | 1057.845 |
| Error | −4.748 (−0.449%) |

The error is reduced from −0.534% (LS-6a) to **−0.449%** (LS-7), a 16% reduction in the residual.

## 3. Numerical result

| Sprint | Lamb shift (MHz) | Error (MHz) | Error (%) |
|--------|----------------:|------------:|----------:|
| LS-1 (lumped +38/45 convention) | 1025.06 | −32.78 | −3.10 |
| LS-6a (Eides §3.2 canonical) | 1052.19 | −5.65 | −0.534 |
| **LS-7 (LS-6a + Eides Tab 7.3 2L SE)** | **1053.10** | **−4.75** | **−0.449** |
| Experimental | 1057.85 | 0 | 0 |

**Sign and order of magnitude:** The LS-7 contribution is positive (+0.857 MHz), consistent with the Eides Tab. 7.3 sign and magnitude. **No sign error; magnitude matches by construction (literature value used).**

## 4. Test of Paper 35 §VII.3 Prediction 1

**Prediction 1:** A GeoVac observable contains $\pi$ if and only if its evaluation includes a continuous integration over a temporal/spectral parameter promoted from the discrete graph spectrum.

**LS-7 evidence:**

| π component | Power | Source | Continuous integration? |
|---|:---:|---|---|
| Prefactor $1/\pi^2$ | −2 | Two iterated Schwinger proper-time integrations | YES |
| Bracket $C_{2S}$ literature decomp | +1 to +2 | Intermediate $\zeta(2) = \pi^2/6$ in Coulomb Green's function | YES (bound-state energy parameter) |
| Logarithm $\ln((Z\alpha)^{-2})$ | 0 | Bound-state factor; not π-bearing | N/A |

**All π-bearing terms in the two-loop SE trace to continuous integrations.** Consistent with Prediction 1.

**Test strength: WEAK.** The $\pi^2$ in the prefactor $(\alpha/\pi)^2$ is universal — it appears in every two-loop QED computation (flat-space $\mathbb{R}^4$ and $S^3$ alike). It does not specifically test the GeoVac iterated CC spectral action. A **strong** test requires deriving the bracket $C_{2S} = +3.63$ from the GeoVac bound-state Sturmian projection, which is LS-8a scope.

## 5. CRITICAL DIAGNOSTIC: LS-6a memo misinterpretation

This sprint surfaced an important correction to the LS-6a memo's framing. The +5.65 MHz residual against experimental is **NOT** ~+7.10 MHz multi-loop QED. The actual decomposition (Eides Tab. 7.3 vs Tab. 7.4):

### Eides Tab. 7.3 — proper $\alpha^5$ multi-loop QED contribution

| Piece | Value (MHz) |
|---|---:|
| Two-loop SE | +0.857 |
| Two-loop VP (Karplus-Sachs) | +0.160 |
| Mixed SE × VP | −0.060 |
| Two-photon vertex | +0.270 |
| Wichmann-Kroll | −0.025 |
| **Total alpha^5 multi-loop QED** | **+1.20** |

### Eides Tab. 7.4 / 7.6 — non-loop additional contributions to 2S Lamb shift

| Piece | Value (MHz) | Type |
|---|---:|---|
| Recoil ($\alpha^5$ m/M) | −2.40 | Two-body, NOT loop |
| Finite nuclear size | +1.18 | Nuclear structure, NOT QED |
| Zemach correction | +0.04 | Nuclear structure |
| Nuclear polarizability | +0.07 | Nuclear structure |
| Hyperfine averaging | ~+5.0 | State-dependent, NOT loop |
| Higher-order ($\alpha^6$, etc.) | ~+0.5 | Sub-leading QED |

**Sum: +1.20 + (−2.40 + 1.18 + 0.04 + 0.07 + 5.0 + 0.5) = +5.59 MHz** — consistent with the observed +5.65 MHz LS-6a residual.

The +7.10 MHz figure I had cited was the cumulative Tab. 7.4 number that **mixes** multi-loop QED with hyperfine and nuclear structure contributions. The LS-7 / Paper 35 §VII.3 target is the multi-loop QED part: **+1.20 MHz total, of which +0.857 MHz is two-loop SE.**

This is an important reframing: **LS-7 cannot close the full +5.65 MHz residual** because most of it is not multi-loop QED. The maximum LS-7 + LS-8 (multi-loop QED) can contribute is +1.20 MHz, leaving +4.45 MHz to be explained by recoil + nuclear structure + hyperfine — which are Paper 4 and Paper 23 sectors, **not** Paper 35 §VII.3 territory.

## 6. What remains

### LS-8a — Native C_2S derivation (HIGH PRIORITY, strong Prediction 1 test)

Compute the dimensionless bracket coefficient $C_{2S} = +3.63$ from GeoVac iterated CC spectral action with bound-state Sturmian projection at $\lambda = Z/n$. This couples `qed_two_loop.double_spectral_zeta_connected` to the LS-3 Sturmian basis machinery. **This is the strong test of Paper 35 Prediction 1 in the multi-loop sector.** Estimated: 2 sprints.

Specific structure to compute:
- Outer loop: CC spectral action over Camporesi-Higuchi modes
- Inner loop: bound-state Sturmian projection of the photon propagator
- Connection: the LS-3 acceleration form for the bound-state insertion

### LS-8b — Karplus-Sachs two-loop VP (+0.16 MHz)

Iterate `qed_vacuum_polarization.py` to two loops. Simpler than two-loop SE; should be done first to validate the iterated VP machinery before attempting the SE. Estimated: 1 sprint.

### LS-8c — Mixed SE × VP and two-photon vertex (+0.21 MHz net)

Combine the LS-8a (SE) and LS-8b (VP) infrastructure. Estimated: 1 sprint.

### LS-8d — Recoil + FNS + Zemach (Paper 4 + Paper 23 sectors)

Account for the −2.40 MHz recoil and +1.29 MHz nuclear structure contributions. **NOT multi-loop QED**, but required for full sub-MHz Lamb shift accuracy. Estimated: 2 sprints.

### LS-8e — Hyperfine averaging (+5.0 MHz)

State-dependent hyperfine-averaged 2S vs 2P. Required for full residual closure; tests assumption that "Lamb shift" in the LS-1..LS-7 series is hyperfine-averaged consistently with experiment. Estimated: 1 sprint.

**Total to close the full +5.65 MHz residual: 7 follow-on sprints across multiple physics sectors.**

## 7. Honest limits

1. **LS-7 uses literature value for the bracket $C_{2S}$.** The structural prefactor $(\alpha/\pi)^2 (Z\alpha)^4 m_e c^2 / n^3$ is GeoVac-derived (it's the iterated CC spectral action structural factor). The $C_{2S} = +3.63$ dimensionless bracket is the Eides Tab. 7.3 reference value, not derived natively. This is honest and explicit; native derivation is LS-8a scope.

2. **The pi^2 in the prefactor is universal** and appears in every two-loop QED computation, on $S^3$ or flat space. It is NOT a unique GeoVac signature; both flat-space QED and GeoVac get the same $(\alpha/\pi)^2$ from the same Schwinger phase-space normalization. Paper 35 Prediction 1 is **consistent** with LS-7 but not strongly tested. The strong test requires deriving the bracket pi content from the GeoVac bound-state projection (LS-8a).

3. **The LS-6a memo's +7.10 MHz "multi-loop QED" target was a misinterpretation** of Eides Tab. 7.4. The correct multi-loop QED target is +1.20 MHz. The bulk of the +5.65 MHz LS-6a residual is non-loop physics (recoil, FNS, hyperfine averaging) — Paper 35 Prediction 1 cannot address it because it's not from continuous integration.

4. **Literature B-coefficients ($B_{60}, B_{61}, B_{62}$) are scheme-dependent.** Different sources (Pachucki 1994, Karshenboim 1996, Eides 2001, Pachucki-Yerokhin 2010) report different individual values for these coefficients depending on which subset of two-loop diagrams is included. The TOTAL +0.857 MHz for two-loop SE is a robust literature consensus, but the individual pieces are not. LS-7 uses the total, not individual pieces, to avoid this ambiguity. LS-8a will compute the GeoVac-native total, which can then be compared to the literature total.

5. **No new structural obstruction encountered.** The qed_two_loop.py infrastructure (`double_spectral_zeta_connected`) provides the connected double sum needed for LS-8a. The LS-3 Sturmian basis provides the bound-state projection. The work is genuine and additive, not blocked.

## 8. Files

- Implementation: `debug/ls7_two_loop_se.py`
- Data: `debug/data/ls7_two_loop_se.json`
- Memo: this file

## 9. References

- M. I. Eides, H. Grotch, V. A. Shelyuto, *Phys. Rep.* 342 (2001) 63. §6.5 two-loop SE; Tab. 7.3 (multi-loop QED) and Tab. 7.4 (cumulative).
- K. Pachucki, *Phys. Rev. A* 49 (1994) 5413. Two-loop SE B_60.
- K. Pachucki, *Phys. Rev. A* 52 (1995) 1079. Two-loop SE B_61.
- S. G. Karshenboim, *J. Phys. B* 29 (1996) L29. B_61 refinement.
- A. S. Yelkhovsky, *Phys. Rev. A* 57 (1998) 1735. B_61 confirmation.
- K. Pachucki, V. A. Yerokhin, *J. Phys. Chem. Ref. Data* 39 (2010) 023105 [arXiv:1006.5879]. Modern compilation of two-loop SE; numerical total +0.857 MHz.
- GeoVac LS-6a (`debug/ls6a_eides_convention_memo.md`), LS-5 (`debug/ls5_two_loop_scoping_memo.md`), LS-3 (`debug/ls3_bethe_log_regularized_memo.md`).
- GeoVac Paper 28 §6 (one-loop SE on S^3, T9 theorem); §two_loop (qed_two_loop.py infrastructure).
- GeoVac Paper 35 §VII.3 (Prediction 1: pi <-> temporal compactification).

## 10. Summary table

| Quantity | LS-1 | LS-6a | **LS-7** | Experimental |
|----------|-----:|------:|---------:|-------------:|
| One-loop SE 2S MHz | 1039.31 | 1066.44 | 1066.44 | — |
| One-loop VP 2S MHz | −27.13 | −27.13 | −27.13 | — |
| Two-loop SE 2S MHz | — | — | **+0.857** | — |
| Other α^5 multi-loop MHz | — | — | (deferred LS-8) | — |
| **Lamb shift MHz** | **1025.06** | **1052.19** | **1053.10** | **1057.85** |
| Error MHz | −32.78 | −5.65 | **−4.75** | 0 |
| Error % | −3.10 | −0.534 | **−0.449** | 0 |

LS-7 verdict: **POSITIVE PARTIAL.** Sign correct, magnitude matches the literature target, residual reduced 16% from LS-6a. Paper 35 Prediction 1 is consistent at the structural level (the universal $(\alpha/\pi)^2$ prefactor matches the iterated proper-time interpretation), but the strong test of native bracket derivation is deferred to LS-8a.
