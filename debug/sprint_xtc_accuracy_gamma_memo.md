# Sprint: per-atom xTC accuracy + the γ-determination question (2026-08-24)

Canonical memo for the accuracy-axis follow-on to the monoatomic sparsity sweep. Question
chain (PI-driven): does xTC buy per-atom *accuracy* (not just sparsity)? → can the geminal
width γ be *determined* rather than fitted (Avery's "solve-and-tabulate")? → does a
*position-dependent* γ(r) (the literature's determined correlation length) beat the
fragility? Verdict: the determined correlation length is real and physics-only, but the
cheap non-Hermitian operator cannot carry it — two-walls / cost-conservation, confirmed on
the accuracy and position-dependence axes.

## 1. Per-atom xTC accuracy (does xTC buy accuracy?)
Largely already answered by the v5.0.2–5.0.8 arc (verified current-state, did NOT re-derive):
- **He:** plain s-only FCI ~17–25 mHa short; **R12-CI (accurate geminal) reaches 0.45–0.80 mHa** — beats chemical accuracy (Paper 59 §f12). Cheap native TC matches ~2.6 mHa only at a hand-picked γ.
- **Li:** cheap native TC ~3× (plain ~25 → physical-optimum ~8 mHa), not chemical.
- **NEW — the electron-count scaling** (`debug/xtc_accuracy_ns_converge.py`, ns-converged s-only):
  after separating basis-incompleteness (ns=3 was radial-incomplete: Li +227, Be +555 mHa)
  from the cusp via ns-convergence, at the oracle γ xTC recovers **94–98%** of the gap
  (He −1.4, Li +1.8, Be +2.0 mHa) — the cusp recovery **scales to Be (4e)**. BUT the
  γ-fragility **worsens with electron count**: E_TC crosses exact (non-variational), and the
  crossing sharpens |dE/dγ| = **40 → 43 → 253** mHa/γ across He→Li→Be (Be caveat: ns=6 not
  fully converged, so its absolute number is soft; the qualitative slope-explosion is robust).

## 2. Can γ be DETERMINED (not fitted)? — global γ (`debug/xtc_gamma_determined_probe.py`)
Framing: "determine-and-tabulate γ" = F12 (fixed transferable geminal). Tested physics-only
global γ candidates (no exact-energy peek): γ=1.0 (F12 std), γ=k, γ=k/2.
- **NO fixed/determined global γ works for the cheap operator.** γ=1.0 → 5/2/**59** mHa
  (He/Li/Be — Be 33× worse than Li); γ=k → 17/12/86; γ=k/2 → 0.6/−21/31. All scatter, Be
  always bad. The per-atom accurate γ is a non-transferable razor's-edge crossing (0.83/0.97/0.60).
- **Why F12 gets away with a fixed γ and we don't:** F12 uses the geminal VARIATIONALLY
  (R12-CI, CI coeffs adapt); the adaptive freedom IS the overlap matrix IS the conditioning
  wall (κ(S)≈200). The cheap non-variational dressing has no adaptive freedom → fixed γ fails.

## 3. Does a determined POSITION-DEPENDENT γ(r)∝n^{1/3} beat it? (`debug/xtc_posdep_gamma.py`)
Explorer literature report (`general-purpose` agent, verified): the physical correlation
length IS position-dependent and determinable — **local range-separation μ(r)=(√π/2)f(r)/n₂(r)**
(Giner–Toulouse–Savin, JCP 149 194301 2018) from the on-top pair density, no energy fitting;
and the **electron-avoidance radius** (Wagner–Gori-Giorgi, PRA 90 052512 2014), a one-electron
sphere (generalized Wigner–Seitz). Both → γ(r) ∝ n(r)^{1/3} in the uniform limit.

Baked that determined shape (bounded, density-weighted-mean-1 normalized) into the cheap
non-Hermitian TC operator, He/Li/Be:
- **Position-dependent γ makes it WORSE, not better.** best accuracy −9/−35/+59 mHa (vs global
  −0.2/+0.4/+2); slopes **80/134/348** (vs global 40/43/253). Three implementations (raw,
  bounded, mean-normalized) all agree.
- **Validated not-a-bug:** `build_posdep` with a flat shape reproduces the engine's global-γ
  ground **bit-for-bit** (diff 0.0 at γ₀=0.8, 1.0 on He).
- **Mechanism:** the cheap operator's non-Hermitian *convective* term feeds on the correlation
  factor's spatial GRADIENTS; position-dependence adds spatial structure → bigger gradients →
  bigger/more-sensitive non-Hermitian term → MORE fragile. Position-dependence amplifies the
  very thing that causes the fragility. The determined γ(r) works in the variational /
  range-separated-DFT framework, NOT as a cheap dressing.

## 4. Skeleton-native version (does equal-area GIVE γ(r)∝n^{1/3}?) — DEFERRED (honest)
Paper 0's equal-area premise is confirmed ("area-proportional state counting", "fundamental
area per state" σ₀=πd₀²/2). BUT its "shells" are the angular-momentum l-subshells (shell k has
2k−1=2l+1 states, l=k−1; grouping n gives the n² principal count), so its packing radius
r_k=k·d₀ is NOT real space. Connecting the equal-area cells to the physical density n(r) needs
the full Fock momentum→position map — a genuine derivation, not a session scaling calc. NOT
faked; flagged as a real open direction for the skeleton program.

## 5. Honest scope / verdict
- **Exact/validated:** the fragility + its electron-count scaling (§1); the fixed-global-γ
  failure (§2); the position-dependent-γ-makes-it-worse result (§3, machinery validated).
- **MEASURED-soft:** Be's absolute accuracy (ns=6 not fully converged).
- **Two-walls confirmed on new axes:** determined γ (global OR position-dependent) cannot rescue
  the cheap non-Hermitian operator; the determined correlation lives on the variational/
  conditioning-limited side. [[two_kinds_of_sparsity]], the CHEM-ACCURACY `/walls` cluster.
- **Open:** the skeleton-native γ(r) derivation (§4, Fock map); the determined γ(r) in a proper
  variational GeoVac framework (would re-incur conditioning).
- Positive kept: Avery's "bake in the formula, not a table" is right in spirit; the formula
  γ(r)∝n^{1/3} exists and is determined-not-fitted (Giner μ(r) / avoidance radius).

Drivers: `debug/xtc_accuracy_per_atom.py`, `xtc_accuracy_ns_converge.py`,
`xtc_gamma_determined_probe.py`, `xtc_posdep_gamma.py`; data `debug/data/xtc_*.json`.
