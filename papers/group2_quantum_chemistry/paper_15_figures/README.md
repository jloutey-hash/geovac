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

### Figure 2: Convergence — D_e % vs l_max, both azimuthal sectors
- **Type:** Grouped bar chart with reference lines
- **Data source:** `tab:extended_convergence` in the paper itself. The figure
  adds no data of its own; every number is copied from that table, so the two
  cannot drift apart.
- **Content:**
  - x-axis: l_max (1, 2, 3, 4, 5, 6)
  - y-axis: D_e / D_e_exact * 100 (%)
  - sigma only (m_max = 0): 37.2, 79.5, 80.2, 87.0, 87.1, 88.8
  - sigma+pi (m_max = 1): 43.7, 86.2, 86.8, 93.6, 93.7, 96.0
  - Reference line at 92.4% — Paper 12, **sigma only**
  - Reference line at 99.1% — Paper 12, **|m| <= 1**
  - Reference line at 100% — exact
  - Every line is labelled with the sector it belongs to, in the figure's own
    legend. A reader must not be able to take an ordering between the two
    coordinate systems from the image; no such ordering is claimed anywhere.
- **Generate with:**
  ```
  python debug/regen_paper15_convergence_fig.py
  ```

> **Regenerated 2026-09-14.** The previous PNG (dated 2026-07-05) plotted
> 30.8 / 37.3 / 87.8 / 88.5 / 95.5 at l_max = 0..4, and **four of those five
> numbers appear nowhere in Paper 15** — its own table gives 37.2 / 79.5 /
> 80.2 / 87.0 for the same sector, so the figure showed 95.5% where the table
> says 87.0%. It also drew a single unqualified line labelled "Paper 12
> (92.4%)" that the bars crossed, which made the coordinate-system comparison
> the paper had withdrawn. The old generator
> (`debug/archive/misc/generate_paper15_figures.py`) wrote to
> `papers/core/paper_15_figures`, a path that stopped existing at the
> 2026-05-22 reorganisation, so the figure had become unregenerable in place —
> which is why it went stale silently and why nothing caught it.

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
