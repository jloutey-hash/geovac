# Paper 15 Figures

Figures for "The Level 4 Natural Geometry" paper.

## Figures to Generate

### Figure 1: Level 4 Coordinate System Schematic
- **Type:** Diagram (TikZ or hand-drawn + scanned)
- **Content:** Two nuclei A, B on z-axis separated by R. Two electrons at
  positions r_1, r_2 from the molecular midpoint. Show:
  - R_e = sqrt(r_1^2 + r_2^2) as the electronic hyperradius
  - alpha = arctan(r_2/r_1) as the correlation angle
  - theta_1, theta_2 as body-frame polar angles from z-axis
  - Label rho = R/(2R_e) as the controlling parameter
- **Generate with:** TikZ in LaTeX or matplotlib schematic

### Figure 2: Convergence Plot — D_e % vs l_max
- **Type:** Line plot with data points
- **Data source:** `geovac/level4_multichannel.py` convergence study
- **Content:**
  - x-axis: l_max (0, 1, 2, 3, 4)
  - y-axis: D_e / D_e_exact * 100 (%)
  - Horizontal dashed line at 92.4% (Paper 12 Neumann V_ee)
  - Horizontal dashed line at 100% (exact)
  - Data points: (0, 30.8), (1, 37.3), (2, 87.8), (3, 88.5), (4, 95.5)
- **Generate with:**
  ```python
  from geovac.level4_multichannel import solve_level4_h2_multichannel
  # Run at R=1.4 for l_max = 0, 1, 2, 3, 4
  ```

### Figure 3: Adiabatic Potential U(R_e) at R = 1.4 bohr
- **Type:** Line plot
- **Content:**
  - x-axis: R_e (bohr), range [0.3, 10]
  - y-axis: U(R_e) (Ha)
  - Multiple curves for l_max = 0, 2, 4 showing how channels deepen the well
  - Mark the minimum for each curve
- **Generate with:**
  ```python
  from geovac.level4_multichannel import compute_adiabatic_curve_mc
  R_e_grid = np.linspace(0.5, 10, 100)
  for lm in [0, 2, 4]:
      U = compute_adiabatic_curve_mc(1.4, R_e_grid, l_max=lm)
  ```

### Figure 4: Fiber Bundle Diagram (Optional)
- **Type:** Conceptual diagram
- **Content:** The double-adiabatic bundle structure:
  Base B1 (R_nuc) -> Fiber F1 (electronic) -> Base B2 (R_e) -> Fiber F2 (angular)
  Show the connection (Berry phase / non-adiabatic coupling)
- **Generate with:** TikZ or hand-drawn
