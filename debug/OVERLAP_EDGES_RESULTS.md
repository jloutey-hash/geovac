# Overlap-Induced Edges: Results

**Date:** March 12, 2026
**Version:** v0.9.38+
**Scripts:** `debug/validate_overlap_edges.py`
**Data:** `debug/data/overlap_edges_h2.txt`

---

## 1. Hypothesis

Adding cross-atom edges to the molecular adjacency matrix proportional to
orbital overlap would make the graph topology R-dependent. As atoms approach,
overlap grows → degree increases → kinetic energy rises → repulsive wall.

If the same edge formula works for both H2 and LiH, the repulsive wall would
be a universal structural property of the graph, not a fitted parameter.

---

## 2. Implementation

New parameters on `MolecularLatticeIndex`:
- `overlap_edges`: edge formula ('s2', 'fock', 'abs', 'kappa')
- `overlap_edge_scale`: multiplicative factor

Adds symmetric edges A[a,b] = f(S_ab) for all s-orbital pairs between atoms.
The degree matrix D[a,a] = Σ_b A[a,b] automatically increases.

---

## 3. Results: H2 (nmax=3, 1540 SDs)

| Config | R_min (bohr) | Error vs 1.401 | k_fit | D_e (Ha) | Has min? |
|:-------|:-------------|:----------------|:------|:---------|:---------|
| **baseline** | 1.147 | -18.2% | 1.15 | 0.39 | Yes |
| s2(1) | 1.147 | -18.2% | 1.15 | 0.39 | Yes |
| s2(10) | 1.050 | **-25.1%** | 1.58 | 0.78 | Yes |
| s2(50) | — | — | — | 4.41 | **No** |
| fock(1) | 1.147 | -18.2% | 1.15 | 0.39 | Yes |
| fock(10) | 1.005 | **-28.3%** | 1.95 | 0.71 | Yes |
| fock(50) | — | — | — | 2.69 | **No** |

**Observation:** Overlap edges push R_min to SHORTER distances (not longer).
At moderate scale (10): well deepens by 2× and R_min shifts 7-10% shorter.
At large scale (50): minimum disappears entirely — system becomes unbound.

## 4. Results: LiH Quick Test (5 R-points)

| R (bohr) | E_baseline | E_s2(10) | ΔE |
|:----------|:-----------|:---------|:---|
| 2.000 | -8.287 | -8.514 | -0.227 (more bound) |
| 2.500 | -8.178 | -8.333 | -0.155 |
| 3.015 | -8.118 | -8.203 | -0.086 |
| 3.500 | -8.092 | -8.129 | -0.037 |
| 4.000 | -8.079 | -8.090 | -0.010 |

**Same behavior:** overlap edges add MORE bonding at short R. The PES gradient
is steeper (more overbinding), not flatter.

---

## 5. Root Cause: Edges = Bonding, Not Repulsion

The approach fails because in graph theory, **edges create delocalization
(bonding), not localization (repulsion)**.

The graph Laplacian is L = D - A. Adding an edge with weight w:
- Increases A[a,b] by w → adds off-diagonal hopping → lowers energy (bonding)
- Increases D[a,a] by w → adds diagonal kinetic penalty → raises energy

For the bonding ground state of a diatomic, the hopping contribution ALWAYS
dominates because the ground state wavefunction has large amplitude on both
atoms. The net effect is **more bonding**, not less.

At very large scales, the overlap-edge hopping overwhelms all other terms,
destroying the PES minimum entirely.

### Why this is fundamental

In real quantum mechanics, Pauli repulsion arises from **orthogonalization**:
when same-spin orbitals overlap, they must be orthogonalized, which raises
kinetic energy. This is a constraint on the wavefunction, not an interaction.

In the graph framework:
- The graph encodes connectivity (who can hop where)
- More connectivity → more delocalization → lower kinetic energy
- There is no graph mechanism for "forced localization due to overlap"
- The Pauli principle is enforced by FCI antisymmetrization, but the
  one-electron kinetic penalty from orthogonalization is missing

### Comparison with t_corr_lambda (diagonal-only approach)

The earlier `t_corr_lambda` approach added corrections ONLY to the H1
diagonal: h1_diag[a] += λ Σ_b S²(a,b). This is equivalent to adding the
degree contribution WITHOUT the off-diagonal hopping. That approach DID
create repulsion (positive curvature in LiH), but was not universal
(opposite sign for H2 vs LiH).

| Approach | Degree (kinetic ↑) | Hopping (kinetic ↓) | Net effect |
|:---------|:-------------------|:--------------------|:-----------|
| t_corr_lambda | Yes | No | Repulsion (not universal) |
| overlap_edges | Yes | **Yes** | **Bonding** (wrong direction) |

---

## 6. Conclusions

1. **Overlap edges in the adjacency matrix do NOT create a repulsive wall.**
   They add bonding (hopping) that overwhelms the kinetic penalty (degree).

2. **The approach fails for both H2 and LiH** in the same way: more
   overbinding at short R, not less.

3. **The graph Laplacian has no natural mechanism for Pauli repulsion.**
   Edges = connectivity = delocalization = bonding. The kinetic penalty
   from orthogonalization is a one-body diagonal effect, not a graph edge.

4. **The diagonal-only approach (t_corr_lambda) was the right direction**
   for repulsion but was not universal. The overlap-edge approach is the
   wrong direction entirely.

5. **The R_eq problem remains open.** Neither approach produces universal
   equilibrium geometry without fitted parameters.

---

## 7. What's Left to Try

1. **Löwdin orthogonalization** of the two-center basis (tried in v0.9.12,
   catastrophically wrong — but the failure mode may have been an
   implementation issue, not a conceptual one).

2. **Explicit exchange repulsion** via off-diagonal V_ee matrix elements
   (currently only diagonal Coulomb J and limited exchange K are included).

3. **Counterpoise-corrected PES** — remove BSSE at every R point and
   check if the corrected PES has the right curvature. This is expensive
   but would isolate whether the curvature problem is from BSSE.

4. **Larger basis (nmax=4,5)** — check if R_min converges toward R_eq
   with increasing basis size. If so, the problem is basis truncation,
   not a missing mechanism.
