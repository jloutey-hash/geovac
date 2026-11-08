# Sprint memo — the cusp exploration: Fock expansion → F12 is native to the Fock momentum representation

Date: 2026-08-19 | Branch: work/sparsity-boundary (uncommitted) | Owning paper:
papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex (new sec:f12)

## One-line
PI-directed "let's try the Fock expansion" for the electron-electron cusp resolved into a
genuine, GeoVac-native result: (a) the cusp lever is explicit *linear* r₁₂ (Kato), not the
distinctive Fock *logarithms* (µHa); (b) explicit-r₁₂ (F12) correlation is **native to the
Fock momentum representation** — the same machinery as the Slater ERI with a regular kernel,
no Gaussian geminal expansion, no RI; (c) BUT the geminal does not lower the genus —
multi-center F12 integrals inherit the SAME elliptic wall as the multi-center ERI, because
**two-electron transcendence is set by density topology, not the operator**. Written up as
Paper 59 sec:f12.

## Context reframe (why the old TC negatives don't block this)
GeoVac's ~5 prior transcorrelated (TC/r₁₂) negatives (CLAUDE.md §3) were **sparsity-motivated**
— explicit r₁₂ densifies the qubit Hamiltonian (2.66× Pauli), killing the sparsity advantage.
In the PI's Avery/continuous pivot (sparsity no longer the constraint), those negatives do not
transfer, so explicit r₁₂ is back on the table. Diagnostic-before-engineering applied (≥2
negatives on the wall) → probes, not an engine.

## Track 1 — Fock-expansion diagnostic on He (driver debug/fock_expansion_he_diagnostic.py)
Minimal Hylleraas variational He calc in (r1,r2,r12), perimetric quadrature (validated exact:
single 1s² at ζ=1.6875 → −2.847656, err 2.5e-7). Isolated the linear-r₁₂ lever from the log lever:

| basis | content | E (Ha) | err |
|---|---|---|---|
| A | analytic, NO r₁₂ (orbital-product floor) | −2.87811 | 25.6 mHa |
| B | A + **linear r₁₂ (Kato)** | −2.90319 | 0.54 mHa |
| C2 | B + r₁₂ ln(r₁₂) | −2.90358 | 0.15 mHa |

**Levers:** linear r₁₂ (Kato) = **25.1 mHa** (the entire correlation); Fock logs = r₁₂ln(r₁₂)
**392 µHa**, R²lnR **18 µHa**, r₁₂ln(r₁+r₂) **6 µHa** (converged nlag 48↔64).
**Verdict:** the cusp/chemistry lever is explicit *linear* r₁₂; the distinctive Fock *logs* are a
µHa precision-physics refinement (Pekeris/Frankowski territory), invisible at chemical accuracy
(1.6 mHa). Vindicates the finite-representation intuition: ONE r₁₂ term = 25 mHa, vs orbital
products stuck at 25.6 mHa needing ℓ→∞.

## Track 2 — F12 is native to the Fock momentum representation
Drivers debug/fock_f12_momentum_probe.py (kernel-swap + pair), fock_f12_genus_probe.py (genus).

- **Kernel-swap exact [MEASURED 1e-14].** A correlated integral ⟨φ_aφ_b|f(r₁₂)|φ_cφ_d⟩ =
  (1/(2π)³)∫ ρ̃_ac* f̃ ρ̃_bd d³k — the IDENTICAL momentum object as the ERI with 4π/k² → f̃(k);
  the j₀ angular closure is kernel-independent. Coulomb same-center reproduces the exact 5a/8;
  the Slater geminal e^{−γr₁₂} has the EXACT regular FT 8πγ/(k²+γ²)² (no Gaussian expansion — the
  first F12 workaround is unnecessary), momentum = direct to 1e-14.
- **2-center pair is genus-0 elementary [SYMBOLIC].** Single-center densities ⇒ rational × j₀ ⇒
  ∫[k²/(k²+γ²)²]j₀(kR)dk = (π/4γ)e^{−γR} (sympy).
- **Geminal does NOT lower the genus [SYMBOLIC + MEASURED].** The two-scale curve period
  ∫dk/√((c₁k²+1)(c₂k²+1)) = (1/√c_max)K(m) exactly (elliptic, genus-1 signature; matches mpmath
  1e-16), set by the TWO two-center densities NOT the kernel. A rational kernel (Coulomb 4π/k² OR
  geminal 8πγ/(k²+γ²)²) is a rational function ON the curve ⇒ f̃·dk/y is still an elliptic integral,
  genus unchanged (geminal trades elliptic-1st-kind for -3rd-kind). Off-diagonal geminal moment
  deviates from any single-scale elementary form (2.8%, 19%); diagonal c₁=c₂ is elementary.

**The map (transcendence = density topology, not operator):**

| two-center transition densities | example | genus | class |
|---|---|---|---|
| 0 | ⟨AA\|f\|BB⟩ | 0 | elementary |
| 1 | ⟨AA\|f\|BC⟩ | 0 | {E₁, ln} |
| 2 | ⟨AB\|f\|CD⟩ | 1 | elliptic (the Paper-59 wall) |

## The unification (the actual result)
In the Fock momentum representation the correlation factor and the Coulomb operator are on
identical footing; the difficulty of any two-electron integral lives in the densities, never in
the operator. So F12 is native and costs the same as the ERI, and **Paper 59's elliptic frontier
gates explicitly-correlated accuracy too** — the same closed form that evaluates the 3-center ERI
evaluates the multi-center F12 correction. What correlation buys is not a cheaper integral but a
compact basis (few correlated integrals vs many repulsion integrals, cusp captured). The F12
thread and the Paper-59 elliptic thread are the SAME problem.

## Deliverables
- Paper 59 new **sec:f12** ("Explicitly-correlated integrals are native to this representation");
  compiles clean (11 pp). Compounds the Phase-4 re-review OWED.
- tests/test_routeC_momentum.py +3: `test_f12_kernel_swap_and_native_geminal`,
  `test_f12_geminal_pair_is_elementary`, `test_f12_multicenter_shares_eri_elliptic_genus` (green).
- drivers: fock_expansion_he_diagnostic.py, fock_f12_momentum_probe.py, fock_f12_genus_probe.py.

## Honest scope
Two-electron integrals only. NOT addressed: the many-electron / resolution-of-identity machinery
of a full F12 method; a resource comparison vs tuned Gaussian F12; whether the elliptic
multi-center integral is tractable at production scale (that is Paper 59's open frontier). Nothing
committed (PI controls commit/tag/Release/QA-recert).
