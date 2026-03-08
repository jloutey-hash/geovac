# Architecture: Lattice vs Hamiltonian

**Question:** Is the Hamiltonian "baked into" the lattice?

**Answer:** NO - but there are two different modes that use the relationship differently.

---

## Core Concept

```
GeometricLattice = BASIS STATES (quantum numbers + graph structure)
                        ↓
        MoleculeHamiltonian = OPERATOR (builds matrices using lattice)
```

The **lattice** defines the quantum numbers (n, l, m) and graph connectivity.
The **Hamiltonian** builds matrices (kinetic, potential) using that structure.

---

## Mode 1: Legacy "Spectral Delocalization" (Original H₂ Method)

### **Used in:** `tests/benchmark_suite.py`, original H₂ calculations

### **How it works:**

```python
# User manually creates lattices
lattice_A = GeometricLattice(max_n=5)  # Left atom: (1,0,0), (2,0,0), ...
lattice_B = GeometricLattice(max_n=5)  # Right atom: (1,0,0), (2,0,0), ...

# Connect them with "bridges" (graph edges)
mol = MoleculeHamiltonian(
    lattices=[lattice_A, lattice_B],
    connectivity=[(0, 1, 16)],  # 16 edges connect lattice A to B
    kinetic_scale=-1/16          # Universal constant
)
```

### **Inside MoleculeHamiltonian:**

```python
# Step 1: Build combined graph
#   Lattice A: states 0-54
#   Lattice B: states 55-109
#   Bridges: 16 edges connecting boundary states

# Step 2: Build Hamiltonian from graph topology
D = degree_matrix(graph)  # How many edges per node
A = adjacency_matrix(graph)  # Which nodes connect
T = kinetic_scale * (D - A)  # Graph Laplacian = kinetic energy

# Step 3: That's it!
H = T  # NO explicit nuclear potentials
```

### **Key insight:**
- The **graph topology** (nodes + edges) IS the physics
- "Bridges" allow electron delocalization between atoms
- NO spatial coordinates needed
- Potential energy is **implicitly** in the graph structure
- Elegant but **abstract** - can't easily add arbitrary nuclei

---

## Mode 2: Multi-Center Chemistry (He/H⁻/H₃⁺)

### **Used in:** `demo/chemistry_lab.py`, new chemistry benchmarks

### **How it works:**

```python
# User specifies PHYSICS (nuclear positions + charges)
mol = MoleculeHamiltonian(
    nuclei=[(0.0, 0.0, 0.0), (1.4, 0.0, 0.0)],  # Explicit 3D coords!
    nuclear_charges=[1, 1],                      # Z values
    max_n=5                                      # Auto-create lattices
)
```

### **Inside MoleculeHamiltonian (multi-step process):**

#### **Step 1: Auto-create lattices** (`__init__` line 1215)
```python
# Create one lattice per nucleus
self.lattices = [GeometricLattice(max_n) for _ in range(n_atoms)]
# Each has states: (1,0,0), (2,0,0), (2,1,-1), (2,1,0), (2,1,+1), ...
```

#### **Step 2: Compute spatial coordinates with Z-scaling** (`_compute_molecular_coordinates` line 1718-1719)
```python
for (n, l, m) in lattice.states:
    Z = nuclear_charges[atom_idx]

    # UNIVERSAL Z-SCALING (Bohr model)
    r_base = n**2         # Hydrogen Bohr radius
    r = r_base / Z        # Contracts for high-Z atoms!

    # Convert (r, l, m) → (x, y, z)
    x = r * sin(θ) * cos(φ)
    y = r * sin(θ) * sin(φ)
    z = r * cos(θ)

    # Offset by nuclear position
    coords[i] = nucleus_pos + (x, y, z)
```

**Example (Helium, Z=2):**
```
State (n=1, l=0, m=0):  r = 1²/2 = 0.5 Bohr  (compact!)
State (n=2, l=0, m=0):  r = 4/2 = 2.0 Bohr
State (n=3, l=0, m=0):  r = 9/2 = 4.5 Bohr
```

**Example (Hydrogen, Z=1):**
```
State (n=1, l=0, m=0):  r = 1²/1 = 1.0 Bohr  (standard)
State (n=2, l=0, m=0):  r = 4/1 = 4.0 Bohr
State (n=3, l=0, m=0):  r = 9/1 = 9.0 Bohr
```

#### **Step 3: Build Hamiltonian with EXPLICIT nuclear potential** (`_build_molecular_hamiltonian` line 1345-1352)
```python
# Graph Laplacian (kinetic energy)
T = kinetic_scale * (D - A)

# EXPLICIT nuclear attraction
V_nuc = np.zeros(n_states)
for i, (n, l, m) in enumerate(all_states):
    Z = nuclear_charges[atom_idx]
    V_nuc[i] = -Z / (n**2)  # Coulomb potential

# Total single-electron Hamiltonian
H_1 = T + V_nuc  # <-- EXPLICIT potential, not baked in!
```

#### **Step 4: Full CI for 2 electrons** (`compute_ground_state` with `method='full_ci'`)
```python
# Tensor product for 2-electron system
H_total = H_1 ⊗ I + I ⊗ H_1 + V_ee

# V_ee computed from spatial coordinates
for i, j in state_pairs:
    r_ij = ||coords[i] - coords[j]||  # Use computed coords!
    V_ee[i,j] = 1 / r_ij  # Electron-electron repulsion
```

---

## Comparison Table

| Feature | **Legacy (Spectral)** | **Multi-Center (Chemistry)** |
|---------|----------------------|------------------------------|
| **Lattice creation** | Manual by user | Auto-created from `max_n` |
| **Spatial coordinates** | NO (abstract graph) | YES (computed with Z-scaling) |
| **Nuclear potential** | Implicit in graph | Explicit V = -Z/r |
| **Electron repulsion** | Graph-based | Coordinate-based (1/r_ij) |
| **Physics input** | Graph topology | Nuclear positions + charges |
| **Best for** | H₂ molecule (known topology) | Arbitrary atoms/ions (unknown topology) |
| **Example** | `lattices=[A,B]` + bridges | `nuclei=[(0,0,0)]` + `nuclear_charges=[2]` |

---

## Visual Diagram

### **Legacy Mode (Abstract Graph)**
```
     Lattice A           Lattice B
    (1,0,0)             (1,0,0)
    (2,0,0)  ←bridges→  (2,0,0)
    (2,1,-1)            (2,1,-1)
      ...                 ...
        ↓
   Graph Laplacian = Hamiltonian
   (No coordinates needed!)
```

### **Multi-Center Mode (Coordinate-Based)**
```
User Input:
  nuclei = [(0,0,0), (1.4,0,0)]
  nuclear_charges = [1, 1]

        ↓

Auto-create lattices:
  Lattice A @ (0,0,0):   states (1,0,0), (2,0,0), ...
  Lattice B @ (1.4,0,0): states (1,0,0), (2,0,0), ...

        ↓

Compute spatial coords with Z-scaling:
  State A[0]: r=1²/1=1.0, coords=(0,0,1.0)
  State A[1]: r=4/1=4.0, coords=(0,0,4.0)
  State B[0]: r=1²/1=1.0, coords=(1.4,0,1.0)

        ↓

Build Hamiltonian matrices:
  T = kinetic_scale * (D - A)   [graph Laplacian]
  V_nuc[i] = -Z/n²              [explicit Coulomb]
  V_ee[i,j] = 1/||r_i - r_j||   [from coords]

        ↓

  H = T + V_nuc  (single-electron)
  H_total = H⊗I + I⊗H + V_ee  (Full CI)
```

---

## How Chemistry Benchmarks Use This

### **Helium (He):**
```python
mol = MoleculeHamiltonian(
    nuclei=[(0, 0, 0)],      # Single nucleus at origin
    nuclear_charges=[2],      # Z=2
    max_n=5
)

# Internally:
# - Creates 1 lattice with 55 states
# - Computes coords: r = n²/2 (contracted by 2x!)
# - Builds H_1 = T + V_nuc with V = -2/n²
# - Full CI: H_total for 2 electrons
```

### **H⁻ (Hydride anion):**
```python
mol = MoleculeHamiltonian(
    nuclei=[(0, 0, 0)],
    nuclear_charges=[1],      # Z=1
    max_n=5,
    kinetic_scale=2.789474    # SPECIAL: positive for anions!
)

# Internally:
# - Creates 1 lattice with 55 states
# - Computes coords: r = n²/1 (standard H scale)
# - Builds H_1 = T + V_nuc with V = -1/n²
# - CRITICAL: positive kinetic_scale prevents overbinding!
```

### **H₃⁺ (Trihydrogen cation):**
```python
nuclei = [
    (0.0, 0.953, 0.0),        # Triangle vertex 1
    (-0.825, -0.476, 0.0),    # Triangle vertex 2
    (0.825, -0.476, 0.0)      # Triangle vertex 3
]

mol = MoleculeHamiltonian(
    nuclei=nuclei,
    nuclear_charges=[1, 1, 1],  # Three protons
    max_n=5
)

# Internally:
# - Creates 3 lattices (one per nucleus)
# - Computes coords for each lattice relative to its nucleus
# - Connects lattices with bridges
# - Builds H_1 with V_nuc from ALL 3 nuclei
# - Full CI: 2 electrons in 165 states (3x55)
# - CRITICAL: Adds V_NN (nuclear-nuclear repulsion)
```

---

## Key Takeaways

1. **Lattice ≠ Hamiltonian**
   - Lattice = basis states (quantum numbers + graph)
   - Hamiltonian = operators built FROM lattice

2. **Two modes, two philosophies:**
   - **Legacy:** Physics in graph topology (abstract, elegant)
   - **Chemistry:** Physics in explicit coordinates (concrete, general)

3. **Why chemistry mode needed:**
   - Can't design graph topology for arbitrary atoms/ions
   - Z-scaling requires coordinates (can't be abstract)
   - Explicit potentials needed for different nuclear charges

4. **Both modes coexist:**
   - Legacy still used for H₂ benchmarks
   - Chemistry used for He/H⁻/H₃⁺/future atoms

5. **The "baking" happens at different stages:**
   - Legacy: Physics baked into **lattice connectivity**
   - Chemistry: Lattice is generic, physics baked into **Hamiltonian matrices**

---

## Summary

You were right to think about "baking" - but it's more nuanced:

- **Legacy mode:** Hamiltonian IS the graph structure (baked into lattice connectivity)
- **Chemistry mode:** Hamiltonian is COMPUTED from coordinates (NOT baked into lattice)

The chemistry mode is what enabled us to solve He/H⁻/H₃⁺ with universal Z-scaling and system-specific kinetic calibration. The lattices are now just a generic basis - the real physics comes from the explicit nuclear positions, charges, and computed potentials.
