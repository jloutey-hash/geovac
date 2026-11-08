# Sprint GD-3 — Reading-B test: the Möbius form is substrate-class-universal

**Date:** 2026-05-29. **Type:** diagnostic (Reading-B follow-on to GD-2 + Route A). **Verdict:** Reading B CONFIRMED, strongly. The Möbius *form* reproduces across two independent azimuthal discretizations — FD-lattice and exact-spectral — to ~4 significant figures. The form is a robust feature of the substrate class (anti-periodic spinor wedge + discrete radial + α-scaling), not a discretization artifact.

## 1. The test

Route A found no continuum Möbius and named the decisive follow-up: does an *alternative* substrate discretization reproduce the Möbius? GD-2 showed the form is t-robust on the FD substrate. GD-3 compares two azimuthal discretizations of the spinor wedge, differing only in the azimuthal eigenvalue:
- **FD-azimuthal** (`DiscreteWedgeDirac`): m_eff = (2/h_φ) sin(π(k+½)/N_φ) — lattice dispersion.
- **spectral-azimuthal** (`DiscreteWedgeDiracSpectral`): m_eff = (k+½)/α — exact continuum half-integer.
Both share the Hermitian polar radial Laplacian (u = √ρ f) and the anti-periodic spinor BC. Slope = [K_wedge − α·K_disk]/(1/α−α). N_ρ=150, a=0.05, N_0=80.

## 2. Result (the two substrates nearly coincide)

| α | t | FD slope | spectral slope | Möbius | SC |
|:-:|:-:|:--------:|:--------------:|:------:|:--:|
| 2 | 1.0 | −0.05623 | −0.05622 | −0.05556 | +0.125 |
| 3 | 1.0 | −0.04883 | −0.04883 | −0.05000 | +0.222 |
| 5 | 2.0 | −0.04905 | −0.04904 | −0.04630 | +0.400 |

Both substrates are Möbius-shaped (closer to Möbius than SC) at **every** tested (α, t), and agree with each other to ~4 significant figures. SC is wrong-sign and wildly off throughout.

## 3. Reading

- **Reading B CONFIRMED:** the Möbius form is substrate-class-universal — independent of the azimuthal discretization scheme. It is set by the anti-periodic spinor BC + radial structure + α-scaling, which both discretizations share, not by any FD/lattice artifact.
- **Combined with Route A** (no continuum Möbius in the published spin-1/2 conical literature): the Möbius form is a *robust* feature of this discrete-substrate class with no published continuum analog.
- **Combined with GD-2** (form t-robust, coefficient + soft-IR t≈1-tuned): the honest picture is now sharp — **Möbius FORM = robust (in t AND across discretization); exact coefficient = t≈1 sweet-spot; soft-IR "mechanism" = t≈1 coincidence.**

## 4. Flagged for the mechanism question (open)

The slope uses the convention Δ_K = K_wedge − α·K_disk (bulk-normalized subtraction). Whether the Möbius-vs-SC difference is a *physical* excess-angle effect or an artifact of this specific subtraction convention is a candidate next probe (compare to an SC-convention subtraction). The mechanism — why this substrate class produces α/(2α−1) — remains open; GD-3 establishes it is real and discretization-robust, not why.

## 5. Documentation
- Paper 51 Möbius caveat: Reading-B confirmation sentence added (form substrate-universal).

## Files
- `debug/sprint_gd3_reading_b_substrate.py`
