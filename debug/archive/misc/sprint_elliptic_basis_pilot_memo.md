# Elliptic-basis pilot — do the bond's CM singular moduli give a better state-space basis? (2026-08-23)

**Question.** The two-center two-electron ERI is a genus-1 period on the Legendre/Γ(2)
family (Paper 59); its modulus is `m = 1 − c_min/c_max` (fock_f12_genus_probe.py), with
CM fibers at singular moduli **disc-4 m=1/2**, **disc-8 m=3−2√2**. That geometry lives on
the *interaction*. Pilot: is it also a useful *state-space basis* principle — if we place
the two 1s exponents so the dominant ERI's modulus is a CM value, is that exponent ratio
variationally distinguished (near the CI minimum), or generic (inert)?

**Bearish prior** (Explorer STOP + our own v4.98 F12 sprint): convergence is cusp-driven;
the period is orthogonal to the cusp. Expected: CM ratio generic.

## Structural finding (while building)
The **single-ζ H₂ basis sits at c₁=c₂ ⇒ m=0 ⇒ genus 0** — the elliptic structure is
*absent* from the minimal basis. It only switches on with ≥2 distinct scales, which is
exactly what the certified single-exponent-per-center closed form (`geovac/qfd_core`,
`two_center_eri`) cannot represent (`hybrid needs one Z on the shared center`). So the
pilot runs on the mixed-exponent **numeric** engine `geovac.sturmian_integrals`
(~1e-4; fast Lmax14/nr1600/nth120 vs high Lmax24/nr3000/nth200 consistent, ordering
preserved), on n=1 1s Slaters (exact per-orbital kinetic `T = a⟨i|1/r_c|j⟩ − ½a²S`).

## Method
H₂, Z=1 each, R=1.4 bohr. Basis: 2 1s-Slaters/center, exponents {z₁,z₂} shared by
symmetry. E_tot = FCI(E_elec) + 1/R. Sweep the ratio β=z_max/z_min (the axis the CM
condition constrains) and 2D-optimize (z₁,z₂). Both scale-map candidates marked
(c∝1/ζ² → disc-4 β=√2, disc-8 β=1.099; c∝1/ζ → disc-4 β=2, disc-8 β=1.207).

## Result — NEGATIVE (as expected), map-robust
| quantity | value |
|:--|:--|
| single-ζ (e=1.0) | E_tot = −1.10662 Ha |
| **true 2D optimum** | z*=(1.109, 1.418), **β_opt = 1.279**, E* = −1.15303 (exact H₂ −1.1744) |
| fixed-scale optimum | β* ≈ 1.35 (smooth single max) |
| CM markers | disc-4: √2=1.414 / 2.0 · disc-8: 1.099 / 1.207 |
| optimal modulus | m(β_opt) = 1−1/1.279² ≈ **0.389 (generic, non-singular)** |
| landscape flatness | E span over β∈[1.25,1.45] ≈ 2×10⁻⁴ Ha |
| E(√2) − E(β*=1.35), high grid | +1.5×10⁻⁴ Ha (√2 on the flat shoulder) |

**The variationally-optimal H₂ basis places the dominant ERI at a GENERIC modulus
(m≈0.39), not at a CM singular modulus.** β_opt=1.279 misses every CM ratio (nearest
disc-8/p1 1.207 off 6%, disc-4/p2 √2 off 10% — outside the ~1e-4 grid noise). The disc-4
ratio √2 is a near-miss sitting ~1.5e-4 Ha above the optimum on a landscape flat to
~2e-4, so its "nearness" confers no advantage; disc-8 and disc-4/p1 are clear misses.

**Verdict.** The interaction's period/CM (complex-multiplication) structure is
**variationally inert as a state-space basis-selection principle** for the two-electron
bond. Directly confirms the Explorer STOP + F12 prior with a measurement, and closes the
loop that sprint left open: convergence is set by ordinary variational optimization (→
ultimately the cusp), not by where the bond's elliptic curve has CM.

## Scope / what is NOT ruled out
2-ζ, 1s-only, homonuclear H₂, 2e — minimal. The **CM-scale-placement** lever is inert.
Two other "elliptic-adapted basis" levers remain UNTESTED: (i) uniformizing-coordinate
basis functions natural on the curve; (ii) graded/resummed max_n expansion along the
modular direction (BD Γ(2) Lambert route). This pilot does not touch either.

Drivers: `debug/elliptic_basis_pilot.py` (engine + FCI), `debug/elliptic_basis_sweep.py`,
`debug/elliptic_basis_focus.py`, `debug/elliptic_basis_focus2.py`.
Data: `debug/data/elliptic_basis_{sweep,focus2}.json`.

## Lever (ii) — modular resummation of the radial convergence series: NEGATIVE (wrong channel)
Diagnostic (2026-08-23): two-center H2 correlation vs basis size, even-tempered + nested 1s
ladders (numeric engine).
- **Grid-artifact flag:** the nested ladder at K=5 forced exponent 5.82 (scale ~0.17 bohr)
  past the radial grid resolution, giving E=-1.1867 **below** exact -1.1744 (variational
  violation) with a blown-up increment ratio ~10. This is numerical breakdown, NOT modular
  structure; discarded. (Caught by the variational-bound check.)
- **Trustworthy part (moderate exponents, K<=4):** the RADIAL s-only series converges fast
  to an s-limit E_inf(s) ~ -1.155 Ha; remaining radial tail ~1-2 mHa.
- **Channel split (the decider):** exact H2(R=1.4) = -1.17447; the ~18-19 mHa gap
  E_inf(s) - E_exact is **angular** correlation (p,d / Kato cusp). The ERI's angular part
  already closed to a **genus-0 Bessel j0**, so the elliptic/modular period lives ENTIRELY
  in the radial channel. A radial resummation can at most buy the ~1-2 mHa radial tail; the
  ~18 mHa that limits accuracy is in the angular/cusp channel the period cannot reach.

**Verdict (ii): NEGATIVE — wrong channel.** Consistent with F12 ("cusp is the lever,
orthogonal to the period"). Drivers: `debug/elliptic_convergence_{diag,nested}.py`.

## Thread verdict — "elliptic bond geometry -> chemical accuracy" CLOSED (measured)
Two independent, measured negatives: (i) CM-scale placement is variationally inert; (ii)
modular radial resummation is in the wrong channel (accuracy = angular/cusp; period =
radial). The bond's elliptic geometry is real and is a decidability/structure object; it is
**inert for chemical accuracy** in every way testable on the two-center engine. The
accuracy-limiting physics is the angular/cusp channel (genus-0, standard F12 territory),
structurally orthogonal to the radial elliptic period. Not a dead end -- a measured closure
of the "geometry doubles as accuracy engine for bonds" question.

## l>0 two-center engine + the PROPER (measured) test (2026-08-23, follow-on)
"Do it properly": built a validated mixed-exponent, l>0, cross-center two-center integral
engine so the angular question is MEASURED, not inferred.

**Engine** `debug/two_center_grid_lm.py` (`TwoCenterLM`): node-less STOs (zeta,l,m,center),
3D (r,u,phi) grid about A, graded radial grid, Coulomb-multipole ERI; one- and two-electron;
real or complex harmonics. **Validated to ~1e-6 against TWO independent background oracles**
(`debug/one_electron_lm_oracle.py` prolate-spheroidal + spherical; `debug/eri_crosscenter_oracle.py`
independent multipole) AND exact closed forms: H-atom T/⟨1/r⟩/S exact to 1e-6; p-density ERI vs
exact closed form 3e-7; cross-center l>0 ERI vs oracle 2.5e-6; m-selection zero at 1e-17.
Two real bugs caught by validation: an outer-radial-integral divergence (l>0 overflow; fixed via
reverse-cumulative ∫_r^∞) and a graded-grid near-origin resolution issue. CI assembly
`debug/two_center_ci_lm.py` reproduces the s-only pilot. Memoization = N^4 → N^2.

**Stage 1 — the gap IS angular (MEASURED).** H2 R=1.4, single-zeta:
| basis | E_tot | gap to exact |
|:--|:--|:--|
| s only | -1.148068 | +26.4 mHa |
| s + p_z (sigma) | -1.152964 | +21.5 mHa (-4.9) |
| s + full p (sigma+pi) | -1.162225 | +12.2 mHa (-9.3) |
Angular functions HALVE the gap to exact; the pi (p_x,p_y, left-right) drop (-9.3) is ~2x the
p_z drop (-4.9). Turns the earlier inference into a measurement.

**Stage 3 — elliptic/CM structure INERT even with angular correlation (MEASURED + coincidence
audited).** Swept beta=zeta_p/zeta_s in s+full-p; the coarse grid-min landed at beta=sqrt2 (the
disc-4 CM ratio, m=0.5) -- but the AUDIT (fine sweep, sqrt2 deliberately off-grid, parabola fit)
gives **beta*=1.475, m=0.541 (GENERIC)**; sqrt2 sits 8e-5 Ha off the optimum = grid noise, NOT
distinguished. The apparent CM hit was a grid-placement + flat-landscape artifact (caught by the
audit discipline). Robust to the c(zeta) map (beta*=1.48 -> generic modulus under either map).

## THREAD VERDICT (now fully MEASURED, not inferred)
The bond's elliptic geometry is INERT for chemical accuracy, established three independent ways:
(i) CM-scale placement variationally generic (s-only AND s+p); (ii) the accuracy gap is angular
(measured: p functions halve it, pi dominant); (iii) the elliptic period is radial, so it does
not reach the angular/cusp channel where the accuracy lives. The skeleton/geometry maps the
structure; chemical accuracy is the angular/cusp (Layer-2, standard F12) channel, which no
discrete geometry shortcuts. Drivers: stage1_angular.py, stage3_{elliptic_angular,audit,final}.py;
xcheck_{1e,eri}_oracle.py. Engine + oracles are reusable validated infrastructure.
