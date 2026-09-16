"""Regenerate paper_15_figures/convergence_lmax.png from the paper's own table.

Why this exists (round-3 DELTA finding B2).  The shipped PNG dates from
2026-07-05 and carries an obsolete solver vintage: it plots 30.8 / 37.3 / 87.8 /
88.5 / 95.5 at l_max = 0..4, and four of those five numbers appear NOWHERE in
Paper 15.  The live table (tab:extended_convergence) gives sigma-only
37.2 / 79.5 / 80.2 / 87.0 / 87.1 / 88.8 at l_max = 1..6.  The figure's 95.5% at
l_max = 4 is 8.5 points above the table's own sigma value for the same
truncation.

It also drew a single red reference line labelled just "Paper 12 (92.4%)", which
the bars cross -- so the image asserted exactly the coordinate-system ordering
that rounds 1-3 withdrew from the prose.  The caption carried the scope; the
image did not, and a figure is read before its caption.

This is a redraw from the paper's authoritative table, not a new measurement:
every number below is copied from tab:extended_convergence, and both Paper 12
reference lines are labelled with the sector they belong to so the figure cannot
be read as an ordering claim on its own.

The superseded generator (debug/archive/misc/generate_paper15_figures.py) also
wrote to papers/core/paper_15_figures, a path that stopped existing at the
2026-05-22 reorganisation -- so the figure had become unregenerable in place.
That is why it went stale silently.
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "papers", "group2_quantum_chemistry",
                   "paper_15_figures")

BLUE = "#0077BB"
ORANGE = "#EE7733"
RED = "#CC3311"
GRAY = "#BBBBBB"
DARK = "#333333"

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 7.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "serif",
    "mathtext.fontset": "cm",
})

# --- tab:extended_convergence, verbatim (2D solver + Schwartz cusp correction)
LMAX = [1, 2, 3, 4, 5, 6]
SIGMA = [37.2, 79.5, 80.2, 87.0, 87.1, 88.8]
SIGMA_PI = [43.7, 86.2, 86.8, 93.6, 93.7, 96.0]

# --- Paper 12 reference values, each tagged with its azimuthal sector
P12_SIGMA = 92.4
P12_M1 = 99.1


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.9))

    x = np.arange(len(LMAX))
    w = 0.38
    ax.bar(x - w / 2, SIGMA, width=w, color=BLUE, edgecolor="white",
           linewidth=0.5, zorder=3, label=r"$\sigma$ only ($m_{\max}=0$)")
    ax.bar(x + w / 2, SIGMA_PI, width=w, color=ORANGE, edgecolor="white",
           linewidth=0.5, zorder=3, label=r"$\sigma{+}\pi$ ($m_{\max}=1$)")

    # Reference lines. Every line and series is named in ONE legend placed
    # outside the axes, so the figure carries its own scope and nothing
    # overlaps -- the previous inline labels collided with each other.
    ax.axhline(y=100, color=GRAY, ls="--", lw=1.0, zorder=2, label="Exact")
    ax.axhline(y=P12_M1, color=RED, ls=":", lw=1.2, zorder=2,
               label=r"Paper 12, $|m|\leq1$ (99.1%)")
    ax.axhline(y=P12_SIGMA, color=RED, ls="--", lw=1.2, zorder=2,
               label=r"Paper 12, $\sigma$ only (92.4%)")

    # The quadrupole onset is the table's own largest step: +42.3 at l_max = 2.
    ax.annotate("quadrupole\nopens", xy=(1 - w / 2, SIGMA[1]), xytext=(0.30, 58),
                fontsize=7, color=DARK, ha="center",
                arrowprops=dict(arrowstyle="->", color=DARK, lw=0.8))

    ax.set_xlabel(r"$l_{\max}$")
    ax.set_ylabel(r"$D_e / D_e^{\mathrm{exact}}$ (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([str(v) for v in LMAX])
    ax.set_ylim(0, 106)
    ax.set_xlim(-0.75, len(LMAX) - 0.25)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2,
              frameon=False, handlelength=1.8, columnspacing=1.2)

    ax.yaxis.grid(True, alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    path = os.path.join(OUT, "convergence_lmax.png")
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print("  [OK] %s" % os.path.normpath(path))
    print("  sigma    :", SIGMA)
    print("  sigma+pi :", SIGMA_PI)
    print("  both columns copied from tab:extended_convergence")


if __name__ == "__main__":
    main()
