# Nuclear Lattice Derivation Notes

**Date:** March 12, 2026
**Status:** Complete draft
**Paper:** Paper 10 — The Nuclear Lattice
**Code:** `geovac/nuclear_lattice.py`
**Tests:** `tests/test_nuclear_lattice.py` (52/52 passing)

---

## 1. Position-Momentum Duality

The electron lattice and nuclear lattice are dual structures:

| Property | Electrons | Nuclei |
|:---------|:----------|:-------|
| **Native space** | Momentum **p** | Position **R** |
| **Dual space** | Position **r** (via Fock projection) | Momentum **P** (conjugate) |
| **S³ origin** | Fock projection of p-space | TBD — projection of R-space? |
| **Angular** | (l, m) on S², SO(3) | (J, M_J) on S², SO(3) |
| **Radial** | n = 1, 2, ..., ∞ | v = 0, 1, ..., v_max |
| **Algebra** | SU(2) ⊗ SU(1,1) | SU(2) ⊗ SU(2) |
| **Extent** | Infinite (ionization) | Finite (dissociation) |

**Key insight:** The electron S³ comes from Fock's stereographic projection of momentum space. The nuclear lattice should come from some analogous projection of position space, but this has not been derived yet. For now, we construct the nuclear graph directly from the algebraic structure of the Morse oscillator and rigid rotor.

---

## 2. The Rotational Paraboloid

### Construction
States |J, M_J⟩ with J = 0, 1, ..., J_max and -J ≤ M_J ≤ J.

**Degeneracy:** 2J + 1 per shell (identical to electron l-shell: 2l + 1)

**Cumulative count:** Σ(2J'+1) from J'=0 to J = (J+1)²

**Adjacency:** L±|J, M_J⟩ = √[J(J+1) - M_J(M_J ± 1)] |J, M_J ± 1⟩

This is structurally identical to the electron (l, m) lattice. Both are SO(3) angular momentum representations. The J=0 state is isolated (no M_J neighbors), just like the l=0 s-orbital.

### Selection rules from graph structure
- ΔM_J = 0, ±1 (from L± adjacency)
- ΔJ = ±1 for electric dipole (from inter-shell coupling)

---

## 3. The Vibrational Chain: Morse SU(2)

### The algebraic identification

The Morse oscillator has an exact SU(2) algebraic structure (Frank & Lemus 1992, Alhassid et al. 1986). The key map:

```
j = ω_e / (2 ω_e x_e) - 1/2        (SU(2) representation label)
m = j - v                             (angular momentum projection)
v = 0, 1, ..., v_max                  (vibrational quantum number)
```

### Ladder operators

The SU(2) lowering operator in the |j, m⟩ basis:

```
J₋|j, m⟩ = √[(j+m)(j-m+1)] |j, m-1⟩
```

Substituting m = j - v:

```
⟨v+1|J₋|v⟩ = √[(2j - v)(v + 1)]
```

**IMPORTANT:** The original implementation used √[(j-v)(j+v+1)] which is WRONG. This formula gives NaN for v > j because (j-v) becomes negative. The correct formula √[(2j-v)(v+1)] is always non-negative for 0 ≤ v ≤ 2j.

### Derivation of the correct formula

Starting from J₋|j, m⟩ = √[(j+m)(j-m+1)] |j, m-1⟩:

With m = j - v:
- j + m = j + (j - v) = 2j - v
- j - m + 1 = j - (j - v) + 1 = v + 1

Therefore: ⟨v+1|J₋|v⟩ = √[(2j - v)(v + 1)]

### Properties of the ladder elements

w(v) = √[(2j - v)(v + 1)] has the following properties:
- w(0) = √(2j) (first edge)
- Maximum at v ≈ j - 1/2 (middle of chain)
- w(v_max - 1) approaches zero near dissociation
- All elements are real and positive for 0 ≤ v < 2j

### v_max determination

The maximum vibrational quantum number must satisfy two conditions:
1. **SU(2) bound:** v ≤ floor(j) (so that m = j - v ≥ 0 in the physical representation)
2. **Energy bound:** E_v ≤ D_e (state must be below dissociation)

We use v_max = max{v : v ≤ floor(j) AND E_v ≤ D_e}.

For ideal Morse (where ω_e x_e = ω_e²/(4D_e) exactly), these two conditions are equivalent. For real molecules (where ω_e x_e is measured independently), condition 2 may be more restrictive.

### Resulting bound state counts

| Molecule | j | v_max | n_states | Expt ~states |
|:---------|:---|:------|:---------|:-------------|
| H₂ | 17.6 | 13 | 14 | ~14 |
| HCl | 27.8 | 17 | 18 | ~23 |
| CO | 81.1 | 81 | 82 | ~82 |
| LiH | 29.8 | 22 | 23 | ~23 |

---

## 4. The Nuclear Graph

### Product construction

G_nuc = G_vib × G_rot (Cartesian product graph)

**Adjacency:** A_nuc = A_vib ⊗ I_rot + I_vib ⊗ A_rot

**Laplacian:** L_nuc = L_vib ⊗ I_rot + I_vib ⊗ L_rot

**Dimension:** N_nuc = (v_max + 1) × (J_max + 1)²

### Vibration-rotation coupling

B_v = B_e - α_e(v + 1/2) modifies the rotational energy at each v:

E_{v,J} = ω_e(v+½) - ω_e x_e(v+½)² + B_v J(J+1) - D_v J²(J+1)²

This is encoded in the diagonal of the nuclear Hamiltonian, not in the graph adjacency. The adjacency encodes the topology (which states are connected); the Hamiltonian adds the energetics.

---

## 5. Electron-Nuclear Coupling

### Born-Oppenheimer fibration

The electron paraboloid sits above each nuclear vertex, with p₀ varying:

p₀(v) = Z_eff(v) / n

At v=0: p₀ = Z/n (standard Fock constraint)
At higher v: p₀ decreases (weaker binding as bond stretches)

### Franck-Condon factors

For vertical transitions between electronic surfaces:
- Same v: FC = 1 (by orthonormality)
- Different v: FC < 1, following a Poisson-like distribution

These factors become edge weights connecting nuclear vertices on different electronic fibers.

---

## 6. Connection to Missing Kinetic Repulsion

The LCAO FCI paper found:
- Graph Laplacian kinetic energy is R-independent
- Virial ratio η ≈ 30 (should be 1 at equilibrium)
- No equilibrium geometry in balanced Hamiltonian

The nuclear lattice provides context:
- R is no longer a classical parameter — it's determined by v
- R_eq corresponds to v = 0 ground state of nuclear graph
- The Fock-weighted correction Δh_aa = λ Σ_b S²_ab (Z_A/n_a)² (Z_B/n_b)²
  captures the cost of overlapping information structures
- Full resolution requires coupled electron-nuclear graph dynamics

---

## 7. ℏ as Discreteness Marker

Every ℏ in the rovibrational spectrum signals discrete lattice structure:

| Quantity | ℏ appearance | Lattice origin |
|:---------|:-------------|:---------------|
| E_v = ℏω(v+½) | ℏω = chain spacing | Vibrational chain edges |
| E_J = BJ(J+1) | B ∝ ℏ² | Rotational paraboloid L² |
| ZPE = ℏω/2 | ℏω = ground state | Chain Laplacian lowest eigenvalue |
| ωx = ℏω²/(4D_e) | ℏ = anharmonicity | Finite SU(2) representation |

---

## 8. Key Equations Summary

### Morse SU(2) identification
```
j = ω_e/(2 ω_e x_e) - 1/2
m = j - v
v_max = max{v : v ≤ floor(j) AND E_v ≤ D_e}
```

### Ladder operators
```
⟨v+1|J₋|v⟩ = √[(2j - v)(v + 1)]     (CORRECT)
⟨v-1|J₊|v⟩ = √[(2j - v + 1)(v)]      (raising)
```

### Vibrational adjacency
```
A_vib[v, v±1] = √[(2j - min(v,v±1))(min(v,v±1) + 1)]  (symmetric)
```

### Morse spectrum (two equivalent forms)
```
E_v = ω_e(v+½) - ω_e x_e(v+½)²                          (standard)
E_v = ω_e(v+½)[1 - (v+½)/(2j+1)]                        (SU(2) form)
```

### Rovibrational energy
```
E_{v,J} = ω_e(v+½) - ω_e x_e(v+½)² + B_v J(J+1) - D_v J²(J+1)²
B_v = B_e - α_e(v+½)
```

### Nuclear graph
```
A_nuc = A_vib ⊗ I_rot + I_vib ⊗ A_rot
L_nuc = L_vib ⊗ I_rot + I_vib ⊗ L_rot
N_nuc = (v_max + 1) × (J_max + 1)²
```

---

## 9. Open Questions

1. **Position-space S³?** Is there an analog of Fock's projection for nuclear position space that would give the nuclear graph the same rigorous S³ foundation as the electron graph?

2. **Graph Laplacian eigenvalues:** The current nuclear Hamiltonian is diagonal (energies are placed on the diagonal). Can the eigenvalues of L_nuc itself encode the rovibrational spectrum through some analog of κ = -1/16?

3. **Polyatomic generalization:** Normal modes → multiple SU(2) chains. Fermi resonance ↔ inter-chain coupling in the product graph.

4. **Non-adiabatic dynamics:** Coriolis coupling, vibronic coupling, and conical intersections as off-diagonal connections between different electronic fibers.

5. **Nuclear kinetic correction:** Can the nuclear graph's v-dependence of p₀ provide the missing R-dependent kinetic energy needed for equilibrium geometries?
