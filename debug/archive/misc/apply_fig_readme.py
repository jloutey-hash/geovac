r"""Bring paper_15_figures/README.md in line with the regenerated figure.

The README still described the obsolete figure -- the same five data points
(30.8, 37.3, 87.8, 88.5, 95.5) that appear in no table of Paper 15, and a
"generate with" recipe pointing at a solver call rather than the script that
actually draws it.  This is what let the figure drift: the documented source
was a sketch, not a runnable generator, and the only real generator wrote to
papers/core/, a path dead since the 2026-05-22 reorganisation.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

R = "papers/group2_quantum_chemistry/paper_15_figures/README.md"

OLD = """### Figure 2: Convergence Plot — D_e % vs l_max
- **Type:** Line plot with data points
- **Data source:** `geovac/level4_multichannel.py` convergence study
- **Content:**
  - x-axis: l_max (0, 1, 2, 3, 4)
  - y-axis: D_e / D_e_exact * 100 (%)
  - Horizontal dashed line at 92.4% (Paper 12 Neumann V_ee, SIGMA-ONLY --
    not a like-for-like reference line for a sigma+pi curve; Paper 12's own
    basis with |m| <= 1 reaches 99.1%)
  - Horizontal dashed line at 100% (exact)
  - Data points: (0, 30.8), (1, 37.3), (2, 87.8), (3, 88.5), (4, 95.5)
- **Generate with:**
  ```python
  from geovac.level4_multichannel import solve_level4_h2_multichannel
  # Run at R=1.4 for l_max = 0, 1, 2, 3, 4
  ```"""

NEW = """### Figure 2: Convergence — D_e % vs l_max, both azimuthal sectors
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
> which is why it went stale silently and why nothing caught it."""

with io.open(R, encoding="utf-8") as fh:
    t = fh.read()

if OLD not in t:
    print("FAILED to match")
    sys.exit(1)

with io.open(R, "w", encoding="utf-8") as fh:
    fh.write(t.replace(OLD, NEW, 1))

print("  + figure README describes the regenerated figure and its provenance")
