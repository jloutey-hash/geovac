# Phase 2: Spectral Determinants & Second Selection Principle
# Date: 2026-03-22

## Executive Summary

Three main results:

1. **The spectral determinant near-miss is confirmed at 0.102%.** The gap between
   det'(S1)*det'(S3)/pi = 41.957 and B = 42 is 0.04276. No simple correction
   factor from the fibration geometry closes it exactly, BUT two additive
   corrections come within 0.003%: `pi/72` and `1/24`.

2. **CORRECTION: d_max(n_max=3) = 3, not 4.** The Phase 1 claim that d_max = 4
   is only true for n_max >= 4. At n_max = 3, all 14 nodes are boundary nodes
   and the maximum degree is 3. The "topological d_max = 4" is the degree of a
   fully interior node, which first appears at n_max = 4.

3. **K/pi is well-approximated by near_miss + F** (error 0.04%), suggesting
   the combination rule K = pi(B + F - Delta) may decompose as
   K ~ pi * (det'(S1)*det'(S3)/pi + zeta(2) - 1/40).


---

## Part 1: Spectral Determinant Values (50-digit precision)

| Quantity | Formula | Value |
|:---------|:--------|:------|
| det'(Delta_S1) | 4*pi^2 | 39.47841760435743... |
| det'(Delta_S2) | exp(1/2 - 4*zeta'(-1)) | 3.195311486059186... |
| det'(Delta_S3) | pi * exp(zeta(3)/(2*pi^2)) | 3.338851214151638... |
| zeta(3) | Apery's constant | 1.202056903159594... |
| zeta'(-1) | derivative at s=-1 | -0.165421143700451... |
| Glaisher A | exp(1/12 - zeta'(-1)) | 1.282427129100623... |

### The Near-Miss

```
det'(S1) * det'(S3) / pi = 4*pi^2 * exp(zeta(3)/(2*pi^2)) = 41.95724178323259...
B = 42
Gap = 0.042758216767411
Relative gap = 0.1018% (confirmed: matches Paper 2)
```

### Key Constants for Reference

```
F = zeta(2) = pi^2/6 = 1.644934066848226...
Delta = 1/40 = 0.025
K/pi = 42 + F - Delta = 43.619934066848226...
K = 137.036064414481541...
```


---

## Part 2: Bridge Search Results

### Best additive corrections to near_miss -> 42

| Correction | Value | Relative error |
|:-----------|:------|:---------------|
| **+ pi/72** | 42.0009 | **2.08e-5 (0.002%)** |
| **+ 1/24** | 41.9989 | **2.60e-5 (0.003%)** |
| + zeta(3)/(2*pi^2) | 42.0181 | 4.32e-4 (0.04%) |
| + 1/40 = Delta | 41.9822 | 4.23e-4 (0.04%) |
| + 1/14 = 1/N(3) | 42.0287 | 6.83e-4 (0.07%) |

**pi/72 is the closest.** Note: pi/72 = pi/(8 * 9) = pi/(8 * n_max^2) at n_max = 3.
This could be a finite-size spectral correction involving the maximum S3 eigenvalue
(|lambda_3| = 8) and the shell count (n_max^2 = 9). However, we have NOT derived
this from first principles.

**1/24 is also very close.** Note: 1/24 appears in the Euler-Maclaurin formula
and in heat kernel coefficients. The B2 Bernoulli number / 2 = 1/12, and 1/24
arises as a regularization correction in zeta function calculations on manifolds.

### Best multiplicative corrections

None come close. The required multiplier is 1.001019..., which doesn't match
any simple expression from the fibration data.

### Pairwise/triple arithmetic search

Only det'(S1)*det'(S3)/pi itself comes within 1% of B = 42.
For K/pi, the closest pairwise hit is det'(S1) + 4 = 43.478 (0.3% off).

### Can K/pi be expressed via spectral determinants?

Best candidates:

| Expression | Value | Error to K/pi |
|:-----------|:------|:--------------|
| near_miss * (1 + F/B) | 43.6005 | -0.04% |
| near_miss + F | 43.6022 | -0.04% |
| near_miss + F - Delta | 43.5772 | -0.098% |

**The pattern near_miss + F is suggestive:** if det'(S1)*det'(S3)/pi were exactly 42,
then K/pi = det'(S1)*det'(S3)/pi + F - Delta would be exact. The 0.04% residual
is the propagated gap from the 0.10% near-miss.

**The multiplicative form near_miss * (1 + F/B) is equally good** and has a nicer
structure: it says the spectral determinant product gets a fractional correction
equal to the fiber invariant divided by the base invariant.

### Analysis of the gap itself

```
gap = 42 - near_miss = 0.042758216767411
gap / F = 0.02599 (not a clean ratio)
gap / Delta = 1.7103 (not clean)
gap * B = 1.7958 (not clean)
1/gap = 23.387 (not recognizable)
gap / pi = 0.01361 (not clean)
```

**Conclusion: The gap does not encode any simple function of F, Delta, or other
project constants.** The near-miss remains unexplained.


---

## Part 3: d_max Analysis — CORRECTION TO PHASE 1

### Critical Finding: d_max(n_max=3) = 3, NOT 4

The Phase 1 summary stated d_max = 4, which is the degree of a **fully interior
node** (one with T+, T-, L+, L- all available). But at n_max = 3, NO NODE is
fully interior:

```
n_max=3 (14 states):
  Degree 1: 4 states  [(1,0,0), (3,0,0), (3,2,-2), (3,2,+2)]
  Degree 2: 8 states  [(2,0,0), (2,1,±1), (3,1,±1), (3,2,-1), (3,2,0), (3,2,+1)]
  Degree 3: 2 states  [(2,1,0), (3,1,0)]
  Degree 4: 0 states  ← NO FULLY INTERIOR NODES
```

The first degree-4 node appears at n_max = 4: state (2,1,0) gains both T± and L±.

### d_max vs n_max (verified 1-15)

| n_max | N_states | d_max | # degree-4 |
|:-----:|:--------:|:-----:|:----------:|
| 1 | 1 | 0 | 0 |
| 2 | 5 | 2 | 0 |
| 3 | 14 | 3 | 0 |
| 4 | 30 | 4 | 1 |
| 5 | 55 | 4 | 5 |
| 6 | 91 | 4 | 14 |
| 7 | 140 | 4 | 30 |
| 10 | 385 | 4 | 140 |
| 15 | 1240 | 4 | 650 |

**d_max = 4 for all n_max >= 4.** The fraction of degree-4 nodes grows monotonically
toward ~52% by n_max = 15.

### Revised Selection Principle

The Phase 1 equation `2 * d_max = n_max^2 - 1` still selects n_max = 3, BUT
the argument must be stated carefully:

- d_max = 4 is the **topological maximum degree** of the lattice type (2 transition
  types x 2 directions), not the realized maximum at n_max = 3.
- At n_max = 3, boundary effects reduce all degrees below 4.
- The equation should be read as: "the S3 eigenvalue at n_max equals the
  bipartite spectral radius bound for this lattice class."

That is: |lambda_{n_max}(S3)| = 2 * d_max(topology)
where d_max(topology) = 4 regardless of n_max.

**This is actually a stronger statement:** the topology of the lattice (2 transition
types giving d_max = 4) determines the asymptotic spectral radius = 8, and n_max = 3
is the unique shell at which the S3 Casimir eigenvalue matches this bound.


---

## Part 4: Transition Operators & Hopf Geometry

### The Two Transition Types

| Transition | Changes | Preserves | Geometric interpretation |
|:-----------|:--------|:----------|:------------------------|
| T± (radial) | n -> n±1 | l, m | Shell-changing = radial motion on S3 |
| L± (angular) | m -> m±1 | n, l | Within-shell = motion along Hopf fiber S1 |

### Hopf Fibration Correspondence

- **L± corresponds to the Hopf fiber S1:** The m quantum number is the U(1) charge
  (eigenvalue of the angular momentum projection J_z). Changing m is rotating
  around the Hopf fiber circle.

- **T± corresponds to radial motion on S3:** Changing n moves between eigenspaces
  of Delta_{S3}. In the Fock stereographic picture, this is changing the "latitude"
  (distance from the north pole). This is NOT motion on the base S2.

- **l is a selection rule, not a transition:** The angular momentum quantum number
  labels the SO(3) representation within each SO(4) shell. It is preserved by
  both transitions. Motion that changes l would correspond to motion on the base
  S2 of the Hopf fibration, but this does NOT appear as a lattice transition.

### Why Only 2 Transition Types (not 3)?

S3 is 3-dimensional, yet the lattice has only 2 independent transition types.
This is because the Hopf fibration collapses the 3 dimensions into:
1. Fiber direction (S1, 1-dim): L±
2. "Radial" direction (perpendicular to fiber, 1 effective dim): T±

The base S2 (2-dim) contributes through the state-labeling structure (which l, m
values exist within each shell) but not through direct transitions. The Fock
stereographic projection maps the 2D base directions into a single radial
variable (the principal quantum number n). This is why d_max = 4, not 6.

### Formula Chain

```
Hopf fibration has 2 independent motion types
  -> d_max = 2 * 2 = 4 (each type has ± direction)
  -> Bipartite bound: lambda_max(graph) -> 2 * d_max = 8
  -> S3 eigenvalue match: n^2 - 1 = 8 => n_max = 3
  -> kappa = E_1 / lambda_max = -0.5 / 8 = -1/16
```


---

## Open Questions for Phase 3

1. **Can pi/72 be derived from spectral geometry?** It equals pi/(|lambda_3| * n_max^2).
   Is there a heat kernel or zeta function correction on S3 that produces this?

2. **The multiplicative approximation near_miss * (1 + F/B) ~ K/pi** has nice
   structure (fiber/base ratio). Can this be made exact?

3. **The d_max = 3 at n_max = 3 issue:** The selection principle works with the
   topological d_max = 4, not the realized d_max. Is there a way to state the
   principle that makes this distinction rigorous?

4. **Why no l-transitions?** The absence of l-changing transitions is a design
   choice in the GeoVac lattice (preserving [H, L^2] = 0). What would happen
   if l-transitions were included? Would d_max = 6 and the selection give n_max = sqrt(13)?


---

## Files

- `phase2_spectral_determinants.py` — All computations
- `phase2_summary.md` — This file
- `degree_analysis.png` — d_max and interior fraction plots
- `bridge_search.png` — Bridge search best hits
