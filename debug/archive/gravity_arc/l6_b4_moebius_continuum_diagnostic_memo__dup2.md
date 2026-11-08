# B4 — Möbius α>1 is a finite-a artifact; the continuum is Sommerfeld–Cheeger

**Date:** 2026-05-29 (L6 bookkeeping task list, item B4)
**Type:** Diagnostic (continuum-extrapolation). **Verdict: clean NEGATIVE on
Möbius — the open mechanism question is RESOLVED.**

## Background
The α>1 conical slope was validated (v3.19.0 Track 5; task #25; Paper 51 §G4-5
follow-on) as the **Möbius closed form**
$$
\mathrm{slope}(\alpha)=-\tfrac{1}{12}\cdot\frac{\alpha}{2\alpha-1},\qquad
\mathrm{recovery}=\frac{\mathrm{slope}}{-1/12}=\frac{\alpha}{2\alpha-1},
$$
where $\mathrm{slope}=\Delta_K/(1/\alpha-\alpha)$,
$\Delta_K=K_{\rm wedge}(\alpha,t)-\alpha\,K_{\rm disk}(t)$. The form matched
sub-2% across α∈{1.5,2,3,4,5,10} but the **mechanism was OPEN**: the v3.19.0
Fursaev–Solodukhin attribution was falsified (Preitschopf, fabricated arXiv ID;
Fursaev–Miele 1996 has no Möbius modification), and Routes A/B/C plus the
harmonic-conjugate derivation all came up empty. GD-2 demoted the soft-IR
"mechanism" to a t≈1 coincidence; GD-3/GD-6 found the FORM robust across t-shape,
azimuthal discretization, and R — **all at fixed lattice spacing a≈0.05**.

## The diagnostic
The prior validation used the **FD azimuthal substrate at fixed a=0.05, single
t=1.0**. This sprint's B1 finding (the FD lattice carries an O(a) Dirichlet
boundary error; the α=1 tip converges cleanly to 1/6 only after a→0) motivated
re-measuring the α>1 recovery with the clean machinery: **SPECTRAL azimuthal
(exact $m_{\rm eff}=(k+\tfrac12)/\alpha$) + radial continuum extrapolation
(a→0) + multiple t**. `debug/l6_b4_moebius_continuum_diagnostic.py`.

## Result — recovery → 1 (SC), not α/(2α−1) (Möbius)
At α=2, t=1, deep radial sweep (a = R/N_ρ, N_ρ = 200…3200):

| N_ρ | a | recovery | (1−recovery) | ratio |
|----:|---:|---------:|-------------:|------:|
| 200 | 0.0500 | 0.6747 | 0.3253 | — |
| 400 | 0.0250 | 0.7690 | 0.2310 | 1.409 |
| 800 | 0.0125 | 0.8367 | 0.1633 | 1.414 |
| 1600 | 0.00625 | 0.8847 | 0.1153 | 1.416 |
| 3200 | 0.00313 | 0.9186 | 0.0814 | 1.417 |

- **Recovery climbs monotonically from ≈Möbius (0.675 at a=0.05) toward SC=1**,
  with no sign of stopping at Möbius (0.667).
- **(1−recovery) ~ √a exactly**: the ratio is bit-stable at 1.416 = √2 (a halves
  → (1−rec) drops by √2), order $a^{0.500}$. So recovery $= 1 - C\sqrt{a} \to 1$.
- The √a convergence is the signature of the **genuine cone singularity at
  α≠1** (the apex is a true deficit/excess-angle cone point), vs the smooth-apex
  O(a) at α=1 (B1). The FD heat-trace tip resolves the cone only as $a\to0$.
- Broader panel (α∈{2,3,5}, t∈{0.5,1,2}): every cell shows the same climb away
  from Möbius toward SC=1, t-stable (t-spread ~0.04); larger α is slower
  (further under-converged at N_ρ=800), consistent with the √a rate.

## Reading
**The Möbius α/(2α−1) form is a finite-a substrate artifact.** At the substrate
scale a≈0.05 the FD-regularized cone gives recovery ≈ α/(2α−1); the continuum
(a→0) limit is the standard **Sommerfeld–Cheeger scalar continuation**
(recovery → 1, slope → −1/12 constant in α), the exact cone heat-kernel tip
$\tfrac{1}{12}(1/\alpha-\alpha)$. The prior multi-sprint "robustness" held only at
fixed a — none of GD-2/3/6 took the radial continuum limit, which is exactly
where Möbius breaks. **The open mechanism question is answered: there is no
exotic Möbius mechanism. That is why FS, Routes A/B/C, and the harmonic-conjugate
attempt all found nothing — there was nothing to find.** The α>1 thread, the one
genuinely-open item in the gravity arc beyond bookkeeping, is **RETIRED** as a
finite-a artifact.

## Honest scope
- recovery→1 with (1−rec)~√a is confirmed cleanly at α=2 (√2 ratio stable over 4
  halvings to N_ρ=3200). At α=3,5 the climb-toward-SC is the same but slower
  (under-converged at the N_ρ used); I claim the **limit is SC** at α=2 and the
  **climb-away-from-Möbius is robust across the panel**, not that every α is
  separately converged to 1.
- This does not touch L6 (the α=1 replica derivative, B1/B2 closed at 1/6) — α=1
  is the smooth-apex case with clean O(a). B4 is the α≠1 cone case.
- Consistency: the SC continuation $\tfrac{1}{12}(1/\alpha-\alpha)$ at α=1
  derivative gives the same 1/6 that B2 confirmed — SC and the α=1 replica agree
  at α=1, as they must.

## Consequence for Paper 51
The §G4-5 follow-on Möbius paragraphs (the α>1 closed form, the α=10 "lock")
should be reframed: the Möbius form is the fixed-a appearance of an object whose
radial continuum limit is plain Sommerfeld–Cheeger. The "mechanism OPEN" caveat
is replaced by "RESOLVED: finite-a artifact, continuum = SC."

## Files
`debug/l6_b4_moebius_continuum_diagnostic.py`,
`debug/data/l6_b4_moebius_continuum_diagnostic.json`.
