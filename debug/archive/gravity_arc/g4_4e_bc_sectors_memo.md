# Sprint G4-4e — BC sectors comparison

**Date:** 2026-05-29
**Verdict:** **POSITIVE-G4-4e-DIAGNOSTIC-CONFIRMED.** Side-by-side comparison on the SAME wedge substrate confirms: spinor (anti-periodic BC) extracts SC at 99.4% mean recovery; scalar (periodic BC) extracts SC at 66.3% mean recovery. **Anti-periodic + half-integer angular momentum is structurally essential for clean conical-defect SC extraction on a discrete substrate.**

## 1. Side-by-side at t = 1.0

Same substrate ($N_\rho = 200$, $a = 0.05$, $N_0 = 120$, proper wedge lattice). Compute $\Delta_K^{\rm sector}/(1/\alpha - \alpha)$ for both BC sectors:

| $\alpha$ | spinor slope | spinor recovery | scalar slope | scalar recovery |
|---|---|---|---|---|
| 2/5 | $-0.08333$ | **99.99%** | $+0.0603$ | 72.4% |
| 1/2 | $-0.0832$ | 99.81% | $+0.0565$ | 67.8% |
| 3/5 | $-0.0827$ | 99.28% | $+0.0531$ | 63.7% |
| 2/3 | $-0.0822$ | 98.65% | $+0.0511$ | 61.3% |

**Mean spinor recovery: 99.4%. Mean scalar recovery: 66.3%. Ratio: 1.50×.**

## 2. F6 sanity at α = 1

- Wedge-Dirac(α=1) vs disk-Dirac: rel_err = 0.00e+00 (bit-exact)
- Wedge-Scalar(α=1) vs disk-Scalar: rel_err = 0.00e+00 (bit-exact)

Both BC sectors correctly reduce to the standard disk at the smooth limit.

## 3. Structural reading

The continuum Sommerfeld/Cheeger formulas give the same magnitude $1/12$ for the leading tip term in BOTH BC sectors (with opposite signs):
- Scalar (continuum): $+\frac{1}{12}(1/\alpha - \alpha)$
- Spinor (continuum): $-\frac{1}{12}(1/\alpha - \alpha)$

On the discrete substrate, the spinor extracts the continuum coefficient to 5 significant figures (G4-4c week 3 best 1.00002), while the scalar misses it by ~34%. The structural reason:

1. **Scalar has a zero-mode** ($m = 0$) with no centrifugal barrier (effective centrifugal $-1/4$ in $u$-representation, attractive). This zero-mode contribution is asymmetric in $\alpha \leftrightarrow 1/\alpha$ because the integer-$m$ discrete lattice samples differ between $\alpha < 1$ ($N_\phi = \alpha N_0$ sites) and $\alpha > 1$ ($N_\phi = \alpha N_0$ more sites).

2. **Spinor has NO zero-mode** (anti-periodic shifts to $m \geq 1/2$). The half-integer lattice samples $m_k = (k+1/2)/\alpha$ for $k = 0, 1, \ldots$ with no privileged origin. This makes the discrete eigenvalue structure SYMMETRIC under $\alpha \leftrightarrow 1/\alpha$ at sprint-scale finite $N_0$, allowing clean SC extraction.

**The half-integer angular momentum structure (forced by fermionic BC) is the structural reason the spinor extracts SC cleanly where the scalar fails.**

## 4. Implication for G4-5 / G4-6

The discrete replica method for $S_{\rm BH}$ (G4-5, multi-month) and the full discrete-substrate $S_{\rm BH}$ derivation (G4-6) operate on the **spinor** sector. G4-4e confirms the spinor SC extraction has clean sprint-scale precision; the load-bearing tip-coefficient identification works at the substrate level. The multi-month G4-5/G4-6 program has a structurally validated foundation.

The SCALAR case (G4-3c-proper / T1) fails to extract SC cleanly at sprint scale because of the BC asymmetry; this is consistent with the scalar case requiring more sophisticated approaches (replica method directly, continuum extrapolation) for $S_{\rm BH}$-class calculations.

## 5. Files
- `debug/g4_4e_bc_sectors.py` + JSON + this memo
