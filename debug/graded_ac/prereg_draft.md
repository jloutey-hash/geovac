# Pre-registered prediction: the CCM Yukawa coefficient at finite cutoff on the S³ window (GA)

**DRAFT for PI review. Not registered. Not in docs/.** Drafted 2026-10-29, before any S³ window value of k_A, k_Φ or c is formed at n_max ≤ 7.
**Commit status:** to be committed on its own, before the scorer (`debug/graded_ac/ga_score.py`, not yet written) is run.
**Rests on:**
- `debug/lit_scan/ccm_spectral_action_yukawa_memo.md`: the relation Y₂(S) = 4g² (CCM07 eq. 4.17), and the reading that the outer factor cancels in the continuum.
- `debug/graded_ac/graded_product_memo.md`: the graded (Dabrowski–Dossena) product, the one-forms and the torus toy.
- `debug/graded_ac/continuum_ref.py`: the smooth-S³ reference.
- Paper 32 prop:odd(c): the true window is odd.

**Prior art:** the coefficient 4 is Chamseddine–Connes–Marcolli's (Adv. Theor. Math. Phys. 11, 991 (2007)). The product with the grading taken from the finite factor is Dabrowski–Dossena (IJGMMP 8, 1833 (2011), eq. 18).

## Why

In CCM the coefficient 4 in Y₂ = 4g² is a ratio of two same-order heat-kernel densities, the gauge F² term and the Higgs |dφ|² term. In the continuum the outer manifold cancels from that ratio. GeoVac's outer factor is a finite truncation: the true Camporesi–Higuchi window. The question is whether the truncation moves the coefficient. Equivalently: does the window's ratio of gauge to Higgs kinetic normalizations match the smooth sphere's?

## The derivation (done before registering)

Product: D = D_GV⊗γ_F + 1⊗D_F on the true window. D² = D_GV²⊗1 + 1⊗D_F² holds exactly (`controls.log`). The γ_GV⊗D_F form does not exist on the window, because no grading commutes with O^CH.

One multiplier pair supplies both fluctuations: a = P D^{1/2}_{−½,+½} P and b = P D^{1/2}_{+½,+½} P.
- Gauge: the one-form a[D_GV,b]⊗Tγ_F with T = P_matter. Since T commutes with D_F, the full-product ε² heat coefficient equals Tr_F(T²e^{−tD_F²}) times the outer coefficient, exactly (toy, to 3.6e-15).
- Higgs: the one-form (ab)⊗[D_F,β]. Its kinetic part is the Duhamel term of C = [D_GV, ab].

So c depends on the outer factor alone, and no D_F enters.

Outer functionals, with D_GV diagonal:
- k_A(t) = (1/2)∂²_ε Tr e^{−t(D_GV+εA)²}|₀, where A = ω + ω† and ω = a[D_GV,b];
- k_Φ(t) = (t²/2)∫₀¹ Tr(C e^{−stD²} C† e^{−(1−s)tD²}) ds.

The continuum reference k_X^∞(t) is the same functional built from the **uncompressed** functions a, b and the **true** product ab. It is summed over every Dirac level n ≤ 51, and a rigorous Gaussian bound covers the rest. This is route (R1).

Route (R2), the closed-form anchor, comes from the Gilkey a₄ density on the unit round S³. Its F-part is (s/6)|F|² with s = 2. The functions entering it are F = d(2 Im a db), with ∫|F|² = 8 vol = 16π², and ∫|dφ|² = (4/3) vol = 8π²/3. Together these give

    k_A^∞ ≈ (2/3)√π t^{1/2},   k_Φ^∞ ≈ (1/3)√π t^{1/2},   R(0⁺) = k_A^∞/k_Φ^∞ = 2.

Scored quantity, at t = τ/(n_max+½)²:

    c(n_max; τ) = 4 · [k_A(n_max,t)/k_A^∞(t)] / [k_Φ(n_max,t)/k_Φ^∞(t)].

The pure-gauge companion takes a = 1 and the same b: A_pg = [D,b] + h.c. Its continuum value is exactly 0, because D + ε i[D,χ] is unitarily equivalent to D. The companion ratio is ρ_pg = |k_A^pg(n_max,t) / k_A(n_max,t)|.

## What I know before comparing (disclosed)

- No S³ window value of k_A, k_Φ, c or ρ_pg has been formed at any n_max ≤ 7.
- One-form norms only, at n_max = 2, 3 (`controls.log`): gauge 4.797118 and 8.702251, Higgs 1.590869 and 2.964635. Both are zero in their controls (Y = 0, b = 1).
- Reference, from `continuum_ref.log`:
  - CG closed form against sympy: 1.1e-16.
  - The new builder reproduces the prototype's window one-forms: singular values agree to 2.2e-16 at n_max = 2 and 8.9e-16 at n_max = 3.
  - Pure gauge k^∞ is ≤ 3.6e-15 in absolute value. Its two terms are each up to 18.6 in size, and the tail bound is ≤ 1.6e-16.
  - (R1) approaches (R2): R = 2.352, 2.168, 2.082, 2.040 at t = 0.16, 0.08, 0.04, 0.02. Then (R−2)/t = 2.20, 2.10, 2.05, 2.02, so the O(t) correction is visible, with coefficient near 2. Over the same t, k_A^∞/(R2 leading) rises from 0.922 to 0.990 and k_Φ^∞/(R2 leading) from 0.784 to 0.970.
  - The denominators at the eight scored points have tail bounds no worse than 6.5e-30 relative (values in the log).
- S³ window at fixed t = 0.1, outside the scored range:

  | n_max | c | \|c/4−1\| | ρ_pg |
  |---|---|---|---|
  | 10 | 3.942427 | 1.4e-2 | 1.2e-2 |
  | 12 | 3.998697 | 3.3e-4 | 1.9e-4 |
  | 14 | 3.99998753 | 3.1e-6 | 1.4e-6 |

  So the S³ pipeline approaches 4 when the cutoff is far inside the window. Where it misses, ρ_pg is the same order as the miss.
- **Toy prior (flat 2-torus; `toy_prior.log`, `controls.log`).** At fixed τ the compressed gauge functional is not gauge invariant, and c is STOP-type. At the toy analogues of the scored points (R² = (n+½)², n = 2–5):

  | τ | c | \|c/4−1\| | ρ_pg |
  |---|---|---|---|
  | 1 | −55.3, −104.8, −161.3, −239.9 | 14.8, 27.2, 41.3, 61.0 | 1.40, 1.21, 1.13, 1.10 |
  | 4 | −2.86, −5.89, −6.86, −10.86 | 1.71, 2.47, 2.72, 3.72 | 3.23, 1.82, 1.67, 1.37 |

  At fixed t the same toy returns 4 to 1e-12.

## Cases (fresh; not computed before)

Eight scored points: τ ∈ {1, 4} × n_max ∈ {2, 3, 4, 5}, with t = τ/(n_max+½)². The construction is exactly as above, using `continuum_ref.Sphere` with window |λ| ≤ n_max + ½. The numerators are computed by the window functions of `continuum_ref.py` (`gauge_eps2`, `kinetic` with T = the whole window). The denominators are the logged (R1) values.

## Predictions

**GA-0 (controls; run FIRST).**
- (a) The scorer reproduces the fixed-t = 0.1 window values for n_max = 10, 12, 14 above, to 1e-12.
- (b) The tail bounds at the eight scored t are below 1e-12 relative.
- (c) The pure-gauge k^∞ is below 1e-12 in absolute value at the eight scored t.

If any of these fails, every item is void.

**GA1 (the relation fails at finite cutoff).** At every scored point, |c/4 − 1| > 1e-3.
- Falsifier: any scored point with |c/4 − 1| ≤ 1e-3.
- Bound: 1e-3 is about 1700 times below the smallest toy miss at a scored analogue (1.71). It is also nine orders above the precision the machinery reaches when finite equals continuum (toy at fixed t: 1e-12; reference tails ≤ 6.5e-30).
- Probability 0.85.

**GA2 (the mechanism is gauge breaking).** At every scored point with |c/4 − 1| > 1e-3, ρ_pg ≥ 0.1.
- Falsifier: any such point with ρ_pg < 0.1.
- Bound: the smallest toy ρ_pg at a scored analogue is 1.10, eleven times above 0.1. The S³ fixed-t points have ρ_pg/|c/4−1| = 0.81, 0.60 and 0.43, so the bound leaves a factor of about 4 even for that weaker coupling between ρ_pg and the miss.
- Probability 0.75.

**GA3 (trend).** At each τ, |c/4 − 1| increases strictly with n_max over 2, 3, 4, 5.
- Toy: strictly increasing at both τ, roughly ∝ (n+½)²/τ. The expected mechanism is that the gauge-breaking term is mass-like, of order t^{−1/2}, while the kinetic term is of order t^{1/2}, so the ratio scales as 1/t.
- Falsifier: any decrease at either τ.
- Probability 0.6. It is lower because the toy is flat and two-dimensional, and S³ has small-n_max level structure the torus lacks.

## Decision gate

- **GO (the registered picture holds):** GA-0 passes, and GA1 and GA2 hold. GA3 is reported, not gated.
- **BORDERLINE:** GA1 holds and GA2 fails at one or two points.
- **STOP:** GA-0 fails (void); or GA2 fails at three or more points (the miss is not gauge breaking, so its mechanism is unknown); or GA1 fails.
- **The interesting outcome, named:** GA1 failing means c is within 1e-3 of 4 at a finite τ where the toy misses by order 1 or more. That would mean something on S³ protects gauge invariance of the compressed functional, either the curvature or the odd window. It would be read as such and given a follow-up probe. It would not be absorbed into the gate as a plain miss.

## Caps

- No change to geovac/, tests/, papers/, CLAUDE.md or memory/.
- One multiplier pair and two τ values. The sizes of c are not predicted, only the threshold, the mechanism and the trend.
- **Scope:** the commutative-factor cutoff on the DD product only. It says nothing on convergence of the almost-commutative extension (WH1) and does not touch Paper 32's H1 non-selection theorem.
- **This is a 3D Euclidean cousin of CCM's 4D relation.** On S³ the shared heat order is t^{1/2}, not the f₀ moment.
- **The product is KO-1** (3 + 6, J = J_GV⊗J_Fγ_F), not CCM's KO-2 (4 + 6).
- The 4 is CCM's, not GeoVac's. Nothing here predicts a top or Higgs mass.
- The (R1) sums are float64. The per-term arithmetic is not exact: the CG square roots enter squared but are summed in floating point.

## Disposition

(Appended after comparing. The predictions above are not edited.)
