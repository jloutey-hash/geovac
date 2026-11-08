# Sprint B1 — joint S³×S¹_T product-carrier convergence (2026-06-10, v3.112.0)

**Goal:** rebuild the convergence content of the descoped Paper 45 on sound footing,
by tensoring Paper 38's unconditional action-seminorm machinery (spatial) with the
WH7 Step-1 Toeplitz temporal algebra (the Q1(i) prerequisite, established v3.111.0).

**Verdict: CLOSED.** Paper 45 gains Proposition `prop:product_action_seminorm`: the
truncated joint system $P\,(C(S^3 \times S^1_T))\,P$, metrized by the translation
Lipschitz seminorm over $SU(2) \times U(1)$, converges to $C(S^3 \times S^1_T)$ in
vS state-space GH distance at the **additive rate** $\gamma_n^{SU(2)} +
\gamma_K^{U(1)} = O(\log n/n) + O(T \log K/K)$, with the spatial $4/\pi$ leading
constant unchanged. Three structural ingredients, each a tensor product of an
unconditional factor ingredient: (i) joint kernel condition ($S^3 \times S^1$
parallelizable → joint Killing fields frame the tangent space; joint per-band
injectivity = product of factor injectivities); (ii) forward arrow = plain
compression, commuting with the joint unitaries → tautologically contractive;
(iii) dual arrow = product lifted state $\xi_s \otimes \xi_t$, defects = sum of the
two factor Fejér smoothings.

## Verification panel (`debug/wh7_b1_joint_product_gh.py`, JSON in `debug/data/`)

Spatial side: real Wigner-D Peter–Weyl window $j \le 1$ (dim 14, the P38
scalar-prototype window) on an exact SU(2) quadrature grid (8×16 GL×16, Haar);
temporal side: $K = 2$ Toeplitz window. Six checks, all PASS:

| Check | Result |
|:--|:--|
| T0 PW orthonormality (grid self-test) | residual 1.2×10⁻¹⁵ |
| T1 pure-factor exactness | temporal err 0.0; spatial err 1.1×10⁻¹⁶ — the joint seminorm restricts to the factor seminorms exactly |
| T2 Leibniz envelope | max(L_s‖B‖, ‖A‖L_t) ≤ L_joint(A⊗B) ≤ L_s‖B‖ + ‖A‖L_t on 6-cell panel |
| T3 joint kernel condition | const 7×10⁻¹⁵; min nonconst 1.68 |
| T4 additive smoothing | ‖Φ_s⊗Φ_t(F) − F‖ ≤ (γ_s + γ_t)·Lip(F) all rows (γ_s = 1.770, γ_t = 0.722 on this window) |
| T5 Riemannian-limit recovery | N_t = 1 reduces joint → spatial **bit-exactly (0.0)** — the load-bearing P45 falsifier preserved |

Frozen falsifier: `tests/test_wh7_b1_joint.py` (6 tests). Paper 45 GATE: PASS, 23 pp.

## Honest scope — what is and is not restored

**Restored:** the product-carrier convergence the withdrawn assembly claimed, with the
temporal direction now carrying genuine metric weight ($[D_t, S_q] = \omega_q S_q$ vs
$[D_t, g(D_t)] = 0$). Q1 clauses (i) and (iii) closed; clause (ii) (Krein compression)
**bypassed, not solved** — the action seminorm never routes the metric through the
Dirac operator.

**Not restored (and not claimed):** anything Lorentzian. The translation seminorm is
signature-blind; the carrier is the Wick-rotated/thermal product. Claims-register
row 11 stays RETRACTED; new row 21 records this Proposition as signature-agnostic.

**Q1 sharpened (applied in Paper 45 §sec:open):** the remaining content of the
Lorentzian target is exactly the *signature* — a seminorm in which causal structure
enters through the **generator** of the metric. Named candidate: a translation-type
seminorm along the wedge boost / modular flow (Paper 42's four-witness identification
$K_\alpha = J_{\mathrm{polar}}$, bit-exact at finite cutoff), i.e. metrize by modular
orbits rather than parameter time — the Connes–Rovelli thermal-time reading made
metric. This is follow-on **B3** (multi-week, genuinely new NCG mathematics).

## Follow-ons
- **B2** (unchanged): non-compact Step 2 — Paley–Wiener band-limit on ℝ_t in the
  pointed-proper framework; decides WH7's primary falsifier branch.
- **B3** (new, the Lorentzian candidate): boost/modular-flow translation seminorm on
  the Paper 42/43 wedge; kernel condition + four-witness compatibility first, rate
  analysis second. If B3 lands, "Lorentzian propinquity" returns as a theorem about
  modular time — which is the only time the operator algebra natively owns.
