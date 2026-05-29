# GeoVac's correspondence-physics position

**Date:** 2026-05-29
**Path:** Continuation thread, Task 10. Synthesis memo articulating where GeoVac actually sits in the correspondence-physics landscape.
**Verdict:** **The framework's actual position is a third structurally distinct category — "discrete spectral skeleton of boundary CFT₃-on-S³" — separated from AdS/CFT by absence of any bulk gravity and from dS/CFT by absence of a Hartle-Hawking wavefunctional mechanism. The propinquity convergence (Paper 38) IS the framework's analog of holographic correspondence, but operates between two abstract C*-algebraic objects rather than between two physical geometries.**

## 1. Three categories of "X/CFT-like" relations

Across mathematical physics, "X/CFT-like" relations come in three structural varieties:

### Category I — Holographic (bulk gravity ↔ boundary CFT)

**Examples:** AdS/CFT (Maldacena 1997), dS/CFT (Strominger 2001, Maldacena 2002), Kerr/CFT, flat-space/CFT (Bondi-Metzner-Sachs analog).

**Structural ingredients:**
- A bulk gravitational theory in $d+1$ dimensions.
- A boundary theory in $d$ dimensions (CFT, often).
- A dictionary mapping bulk fields to boundary operators.
- A partition-function or wavefunction identification: $Z_\text{bulk}[\phi_0] = Z_\text{boundary}[\phi_0]$ or $\Psi_\text{bulk}[\phi_0] = Z_\text{boundary}[\phi_0]$.

**Operationally:** boundary CFT data is RECONSTRUCTED from bulk gravity by holographic computation. Equivalently, bulk gravity is computed from boundary CFT via reconstruction (RT/HRT minimum surfaces, etc.).

### Category II — Topological (discrete bulk ↔ boundary state)

**Examples:** Chern-Simons / WZW on a boundary (Witten 1989), state-sum models / boundary TQFT, tensor networks / MERA.

**Structural ingredients:**
- A bulk in $d+1$ dimensions, typically with a discrete or topological structure.
- A boundary theory on a $d$-manifold, derived from bulk topological data.
- An identification of the boundary partition function with a bulk invariant.

**Operationally:** boundary observables are produced as TOPOLOGICAL invariants of the bulk discrete structure. The bulk is finite or discrete; the boundary CFT (or its analog) emerges in a continuum limit.

### Category III — Spectral-triple discretization (substrate ↔ continuum CFT)

**Examples — GeoVac is the headline case.** Connes-van Suijlekom 2021 spectral truncation of the circle $S^1$ converging to round-$S^1$ CFT; Paper 38's WH1 PROVEN proves the same for $S^3$ via Latrémolière propinquity; Paper 50's bit-exact KPS matches reproduce CFT_3 partition function data.

**Structural ingredients:**
- A discrete spectral triple $(A_{n_{\max}}, H_{n_{\max}}, D_{n_{\max}})$ at finite cutoff.
- A continuum spectral triple $(A_\infty, H_\infty, D_\infty)$ in the limit.
- A propinquity-style metric $\Lambda$ on the space of spectral triples.
- Convergence: $\Lambda(T_{n_{\max}}, T_\infty) \to 0$ at proven rate (Paper 38: $C_3 \cdot \gamma_{n_{\max}}$ with $\gamma_n \sim 4/\pi \cdot \log n / n$).
- Continuum CFT_d is the SPECTRAL CONTENT of the limit triple, extracted via zeta-function derivatives (KPS-style).

**Operationally:** boundary CFT data is COMPUTED DIRECTLY from spectral asymptotics of the discrete substrate, not reconstructed from any bulk. There is no "bulk" — the discrete substrate IS the framework's intrinsic object, and the continuum CFT is what the substrate looks like at the propinquity limit.

## 2. Why GeoVac is Category III, not Category I

Three structural reasons:

**Reason 1: No bulk geometry exists.** Paper 50 §6 explicitly states no AdS_4 / H^4 / S^4 machinery in `geovac/`. The framework's only geometric structure is the discrete substrate on $S^3$ (the Fock-projection image of the hydrogen spectrum). There is no candidate bulk to be holographic FROM.

**Reason 2: No holographic dictionary.** A Category-I theory maps bulk fields $\phi_\text{bulk}(x, z)$ to boundary operators $\mathcal{O}_\text{boundary}(x)$ via an asymptotic relation $\phi_\text{bulk}(x, z \to 0) \sim z^\Delta \mathcal{O}_\text{boundary}(x)$. GeoVac has no $z$ direction. The spectral triple has labels $(n, l, m, s)$ for the discrete substrate and the boundary CFT_3 emerges in the continuum limit of these spectral labels, NOT via an extra dimension being absorbed.

**Reason 3: No Witten / Hartle-Hawking functional integral.** Category-I theories produce boundary CFT data by INTEGRATING OVER BULK FIELDS (or geometries) with appropriate boundary conditions. GeoVac produces CFT_3 data by computing SPECTRAL TRACES on the discrete substrate: $F = -\frac{1}{2}\zeta'_\Delta(0)$ is a zeta-function-derivative computation, not a functional integral.

## 3. Why GeoVac is also distinct from Category II

GeoVac is not a topological / state-sum / tensor-network construction either:

**Distinction T1:** Category-II theories are typically TOPOLOGICAL on the bulk — observables depend on topology, not metric. GeoVac's discrete spectral triple has SPECTRAL content (Dirac eigenvalues $|\lambda_n| = n + 3/2$ on $S^3$) — explicitly metric-dependent. The propinquity bound depends on the metric of the round-$S^3$ continuum limit.

**Distinction T2:** State-sum models and tensor networks construct the boundary as an EMERGENT object from bulk discrete data. GeoVac's boundary CFT_3 is NOT emergent — it IS the substrate (the $S^3$ is the framework's natural geometric arena, not a boundary of something).

**Distinction T3:** Category-II theories have a "bulk = discrete, boundary = continuum" pattern. GeoVac has "substrate = discrete spectral triple, continuum limit = continuum spectral triple." Both are spectral-triple objects; one is a discretization of the other. No bulk/boundary asymmetry.

## 4. What IS GeoVac's "correspondence"?

The framework's analog of the holographic correspondence is:

$$
\boxed{\text{Propinquity convergence:}\quad \Lambda(T_{n_{\max}}, T_\infty) \to 0 \text{ at rate } C_3 \gamma_{n_{\max}}.}
$$

This is Paper 38 WH1 PROVEN. It is a correspondence between:
- $T_{n_{\max}}$: the discrete spectral triple at finite cutoff $n_{\max}$.
- $T_\infty$: the continuum spectral triple on round $S^3$.

Both objects are **operator-algebraic**, not geometric in the differential-geometry sense. The convergence is in the **Latrémolière propinquity metric**, which measures how close two spectral triples are in their algebra-Lipschitz-Dirac data.

The boundary CFT_3 partition function values (Paper 50 §3) are RECOVERED in the continuum limit because the zeta function $\zeta_\Delta(s)$ of the discrete spectral triple converges to the continuum zeta in the propinquity sense, and $F = -\frac{1}{2}\zeta'(0)$ is a derivative at zero that is preserved by the limit.

This is a different kind of object than any holographic correspondence. It is not bulk-boundary; it is discrete-continuum. The "correspondence" is between two C*-algebraic discretization levels, not between two physical theories.

## 5. The Wick-rotation thermal structure is the boundary side of any holographic theory

Papers 42, 43, 47, 49 establish the wedge KMS / Tomita-Takesaki / Bisognano-Wichmann / Hawking / Unruh / Sewell structure on the GeoVac wedge of $S^3$. This is the **modular-flow structure of the wedge KMS state**.

Any holographic correspondence involving thermal states (e.g., AdS-Schwarzschild ↔ thermal CFT, dS-Schwarzschild ↔ static patch CFT, BTZ ↔ thermal $S^1 \times \mathbb{R}$ CFT) has the SAME modular-flow structure on its boundary CFT side. The four-witness theorem (Paper 42) explicitly identifies this — Hawking, Sewell, Bisognano-Wichmann, and Unruh are four readings of one modular-flow structure.

So GeoVac shares the BOUNDARY-MODULAR-FLOW content with all holographic-correspondence theories. The difference is in what the BULK is:
- AdS/CFT: bulk is anti-de Sitter
- dS/CFT: bulk is de Sitter
- GeoVac: no bulk; the substrate is a discrete spectral triple

The Connes-Rovelli thermal-time identification (Paper 49) makes this concrete: thermal time IS the modular flow, regardless of whether you came at it from a bulk-gravity calculation or a direct boundary-CFT analysis. GeoVac is the latter.

## 6. The framework's correspondence position in one sentence

> **GeoVac is the discrete spectral-triple realization of the boundary CFT_3-on-S^3 partition function data and modular structure, distinguished from AdS/CFT and dS/CFT by absence of any bulk gravitational theory and presence of a propinquity-convergence "correspondence" between operator-algebraic discretization levels.**

## 7. Cross-correspondence catalogue

For Paper 50 §8 (or a future positioning paper), it is worth tabulating where GeoVac sits across the landscape:

| Theory | Bulk | Boundary CFT | Correspondence mechanism | GeoVac analog? |
|---|---|---|---|---|
| AdS_5/CFT_4 (Maldacena 1997) | AdS_5 IIB SUGRA | $\mathcal{N}=4$ SYM on $S^3 \times \mathbb{R}$ | Holographic dictionary | NO bulk; CFT_4 not in scope |
| AdS_4/CFT_3 (Aharony-Bergman-Jafferis-Maldacena 2008) | AdS_4 M-theory | ABJM on $S^3 \times \mathbb{R}$ | Holographic dictionary | NO bulk; partial CFT_3 match (Paper 50 KPS) but different CFT |
| dS_4/CFT_3 (Strominger 2001, Maldacena 2002) | dS_4 | Boundary CFT_3 on $I^+ \cong S^3$ | Hartle-Hawking | NO bulk; boundary CFT_3 ADJACENT (Paper 50 KPS reproduces standard data) |
| Kerr/CFT (Guica-Hartman-Song-Strominger 2009) | Kerr | Chiral CFT_2 | Holographic, asymptotic-symmetry | NO Kerr; no CFT_2 |
| Chern-Simons/WZW (Witten 1989) | 3D CS bulk | 2D WZW boundary | Topological / state-sum | NO topological CS; not in scope |
| Connes-vS truncations (2021) | Discrete spectral triple on $S^1$ | Round-$S^1$ spectral triple | Propinquity convergence | YES — GeoVac extends to $S^3$ |
| **GeoVac (Papers 38, 50)** | **Discrete spectral triple on $S^3$** | **Round-$S^3$ spectral triple + CFT_3 partition data** | **Latrémolière propinquity convergence at rate $C_3 \gamma_{n_{\max}}$** | **Self** |

The table places GeoVac in Category III as a natural extension of the Connes-vS Toeplitz $S^1$ truncation, with concrete CFT_3 partition function data as the witness that the continuum limit is recognizably a CFT (rather than just an abstract spectral triple).

## 8. What's publishable

Three options, in increasing scope:

**Option α — Paper 50 §8 extension (~1-2 days, low cost).**
Add §8.X "Structural distinction from holographic correspondences (AdS/CFT, dS/CFT)" with the §6 table + 1-2 paragraphs of structural reading. The headline goes in: GeoVac is Category III, not Category I or II.

**Option β — Standalone short paper (~1-2 weeks, math.OA-adjacent).**
Working title: *"Discrete spectral-triple realization of CFT_3-on-S^3: a non-holographic alternative to dS/CFT for boundary CFT data."* 12-15 pages, math.OA or math-ph venue. Builds on Paper 38 (propinquity convergence) and Paper 50 (bit-exact CFT_3 data); positions the framework explicitly against AdS/CFT and dS/CFT.

**Option γ — Full positioning paper (~1-2 months, hep-th / math-ph audience).**
Working title: *"Spectral-triple discretization as a correspondence mechanism: GeoVac, Connes-van Suijlekom, and the structure of non-holographic CFT realizations."* 25-30 pages. Surveys Category I/II/III; positions GeoVac as a paradigm of Category III; engages with broader landscape including Mondino-Sämann synthetic Lorentzian framework, etc.

## 9. Recommendation

**Option α (Paper 50 §8 extension) for the immediate move.** Cost is small; captures the headline structural distinction; locks in the priority statement.

Option β is the natural follow-on IF the PI wants to invest more in this framework-positioning direction. Option γ is multi-month and should be parked unless calibration arc lands a major new result.

## 10. Honest scope

This memo:
- **Synthesizes** the dS/CFT scoping (Task 9) plus Papers 38, 50, 42, 43, 47, 49 into a unified position statement.
- **Articulates** GeoVac as Category III, distinct from Category I (holographic) and Category II (topological).
- **Names** three publishable options ranging from §8 extension to full positioning paper.

Does NOT:
- Execute any of the options.
- Modify Paper 50 or any other paper.
- Claim GeoVac IS a new correspondence in the holographic sense.

## 11. Cross-references

- `debug/sprint_ds_cft_scoping_memo.md` — Task 9, structural-overlap finding
- Paper 50 §3 (Klebanov-Pufu-Safdi match), §4 (wedge KMS), §6 (bulk-side blocked), §8 (catalogue extension venue)
- Paper 38 — WH1 PROVEN: propinquity convergence (the "correspondence" mechanism)
- Paper 42 — four-witness Wick-rotation theorem at finite cutoff
- Paper 43 — Lorentzian extension
- Paper 47 — two-rate hybrid (non-compact extension)
- Paper 49 — Connes-Rovelli thermal-time stack
- Connes-van Suijlekom 2021 — Toeplitz $S^1$ truncation (Category III precedent)
- Strominger 2001 (arXiv:hep-th/0106113) — dS/CFT proposal
- Maldacena 2002 (arXiv:astro-ph/0210603) — Euclidean dS/CFT
- `memory/geovac_structural_skeleton_scope_pattern.md` — structural-skeleton scope

## 12. Files

- `debug/geovac_correspondence_position_memo.md` (this)
