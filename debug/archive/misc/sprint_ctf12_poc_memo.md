# Sprint memo: CT-F12 PoC on a Coulomb-Sturmian basis (He)

**Date:** 2026-08-23  **Scope:** debug/ only; no production/paper edits.
**Question:** does a Kato-cusp-enforcing Slater geminal (f12 = e^{-γr12}) let a *materially
smaller* basis reach a target correlation-energy accuracy than plain FCI, on GeoVac's own
Slater/Coulomb-Sturmian basis (not Gaussian)? Ref: Motta et al., PCCP 22, 24270 (2020).

## Method
Coulomb-Sturmian s-basis S_{n0}(r)=hydrogenic_radial(r,n,0,Z=n·k) (shared decay e^{-kr};
size = #s-functions n_s; k = single scale). Reused the validated SturmianCI angular
(Gaunt c^k) + Löwdin + Slater–Condon FCI machinery; only the 2e radial kernel is swapped.
Two realizations of the geminal were tested:

- **Route A — Hermitian two-body *effective-Hamiltonian dressing* (the literal ask).**
  Canonical/Hermitian transcorrelation drops the non-Hermitian convective term K from
  H̃ = e^{-τ}He^{τ} = H + D + K (τ=f(r12)).  For one pair (He) the multiplicative part gives
  the effective interaction  w(r) = (1-e^{-γr})/r + (γ/2)e^{-γr} − (¼)e^{-2γr}  (finite at 0).
  Diagonalize H with 1/r12 → w(r12).  (A2 = symmetrized variant w=1/r − ¼e^{-2γr}.)
- **Route B — variational explicitly-correlated CI (R12-CI), the correct Hermitian route.**
  Augment the orbital 2e space with geminal pair functions G_ref = e^{-γr12}R_ref(r1)R_ref(r2)
  and Rayleigh–Ritz.  Variational ⇒ bounded below by E_exact.  s-only orbitals are isotropic,
  so the geminal alone carries the θ12 (angular) correlation.  Kinetic via the gradient form
  ½∫(∇₁A·∇₁B+∇₂A·∇₂B) (first derivatives only); all 2e integrals reduce to (r1,r2,x) quadrature.

## Geminal-integral validation (every geminal integral checked vs direct r12 quadrature)
- geminal multipole kernels g_k^w(r1,r2): Gauss–Legendre vs scipy.quad → **|Δ|<5e-6**.
- Coulomb multipole vs analytic r<^k/r>^{k+1} → **<1e-13**.
- end-to-end (1s1s|w|1s1s): grid assembly vs independent (r1,r2,x) Laguerre quad → **3.8e-6**.
- ⟨1s|T|1s⟩ grid = k²/2 **exact**; single-1s He grid = −2.847672 vs analytic −2.847656 (**16 µHa**).
- **R12-CI variational check: every reported energy lies ABOVE exact** (min margin 0.45 mHa) — genuine.

## Results — accuracy vs basis size (He, exact −2.903724; s-limit −2.879029)
| method | basis | q=2·n_orb | E (Ha) | err vs exact |
|---|---|---:|---:|---:|
| plain **s-only** FCI | 7 s-Sturmians | 14 | −2.879009 | **24.7 mHa (s-limit plateau, any size)** |
| plain FCI s,p,d | 14 orb (max_n 3) | 28 | −2.885277 | 18.4 mHa |
| plain FCI s,p,d,f,g | 55 orb (max_n 5) | 110 | −2.886637 | **17.1 mHa (still stuck)** |
| **R12-CI** | **3 s-Sturm + 1 geminal** | **6** | **−2.902926** | **0.80 mHa** |
| **R12-CI** | 4 s-Sturm + 2 geminal | 8 | −2.903275 | 0.45 mHa |

Route A (naive Hermitian dressing): **collapses** — H+D gives −2.97 … −3.24 Ha (25–330 mHa
*below* exact) and *worsens* with basis size (non-isospectral; dropping K breaks the spectrum).
Route A2 (symmetrized): bounded but **non-variational** — improves then *overshoots* exact
(−2.916 at max_n=3, γ=1), so untrustworthy as a target-accuracy method.

## Verdict — **SPLIT: GO on the geminal physics; STOP on the naive two-body dressing**
- **The cusp-enforcing geminal buys a large basis reduction (GO).** One Slater geminal added to
  a 3-s-orbital basis reaches **0.80 mHa** of exact (qubit-proxy 6). Plain FCI in the same
  Sturmian family is **stuck at ~17 mHa even at q=110** (55 orbitals) and s-only FCI plateaus at
  24.7 mHa forever. So plain FCI never reaches the ~1 mHa target at any tested size, while R12-CI
  reaches it at q≈6. **Basis-size / qubit-proxy ratio at matched (~1 mHa) accuracy: ≥18×
  (110/6) and effectively unbounded** (plain FCI plateaus above target). In raw variational-
  function counts: 7 functions (R12-CI) beat 1540 orbital pairs (plain, q=110). The geminal
  replaces the entire high-l partial-wave tower — the classic F12 win, reproduced on GeoVac's
  own Slater/Sturmian basis (Paper 59 §f12: geminal integrals are native here).
- **The naive Hermitian two-body *effective-Hamiltonian* dressing does NOT work (STOP).** Drop-K
  (H+D) collapses; the symmetrized variant overshoots. Neither is isospectral to H, so neither is
  a trustworthy Hermitian potential-swap — consistent with CLAUDE.md §3 (prior TC bolt-ons
  plateaued/failed). The working Hermitian realization keeps the geminal as an *explicit
  correlated basis function* (R12-CI), which is variational but does **not** fold into a simple
  2-body qubit term. True qubit-count reduction à la Motta therefore needs the non-Hermitian TC
  (deliberately avoided here) or the projector-based [2]_R12 with a CABS — not a naive dressing.

**Honest caveats.** (1) He-only, s-only orbitals + geminal (single pair ⇒ no 3-body terms; the
clean case). H2 (two pairs) would introduce 3-body TC terms and a genuine CABS/projector — not
attempted. (2) The plain-FCI "17 mHa floor" partly reflects the shared-k Sturmian basis being
radially under-complete at high l (compounding slow partial-wave convergence); a per-l-optimized
STO basis would do better, but R12-CI uses the *same* s-Sturmians, so the head-to-head is fair and
the geminal advantage is basis-family-internal. (3) The qubit-proxy=2·n_orb undercounts the
geminal (an extra correlated function, not an orbital); the ≥18× figure is the qubit-proxy ratio
under the R12-CI packaging, not a folded-Hamiltonian qubit count.

## Files
- debug/ctf12_validate.py — operator/cusp + geminal-multipole validation vs direct quad.
- debug/ctf12_poc_he.py — Route-A effective-Hamiltonian dressing (collapse) + grid ERI engine.
- debug/ctf12_r12ci_he.py — Route-B variational R12-CI (the working method).
- debug/data/ctf12_r12ci_he.json, debug/data/ctf12_summary.json — all numbers.
