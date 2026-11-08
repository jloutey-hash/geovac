"""Summary figure for the BeH2 off-axis (bending) conical-intersection study (Josh 2026-08-25).
Left: log10 min-gap map over (d, theta) with the three CI branches overlaid.
Right: CI location d*(theta) per branch + the sigma/pi coupling driving the 1.7 promotion.
Reads beh2_bending_ci.json (gap map), _symmetry.json (census), _berry.json (robust Berry)."""
from __future__ import annotations
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = json.load(open("debug/data/beh2_bending_ci.json"))
CEN = json.load(open("debug/data/beh2_bending_ci_symmetry.json"))["census"]
BER = json.load(open("debug/data/beh2_bending_ci_berry.json"))

gm = M["gap_map"]
dg = np.array(gm["dgrid"]); th = np.array(gm["theta"]); G = np.array(gm["gap"], float)

fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))

im = ax[0].pcolormesh(dg, th, np.log10(np.clip(G, 1e-7, None)), shading="auto", cmap="magma")
plt.colorbar(im, ax=ax[0], label="log10  min F-gap  (dark valleys = degeneracies)")
# overlay CI branches (from census: dot per protected even<->odd CI, color = Berry phase)
for c in CEN:
    for e in c["cis"]:
        pi = e["berry"].get("0.1", e["berry"].get("0.06", "0"))
        col = "cyan" if pi == "PI" else "orange"
        ax[0].plot(e["dstar"], c["theta"], "o", ms=5, mfc=col, mec="k", mew=0.4)
ax[0].plot([], [], "o", mfc="cyan", mec="k", label="protected CI (Berry pi)")
ax[0].plot([], [], "o", mfc="orange", mec="k", label="even<->odd touch (Berry 0)")
ax[0].axhline(180, color="w", ls=":", alpha=0.5)
ax[0].set_xlabel("Be-H bond length d (bohr)")
ax[0].set_ylabel("H-Be-H angle theta (deg)")
ax[0].set_title("min F-gap over (d, theta) + CI census")
ax[0].legend(loc="lower left", fontsize=8)

# right: d*(theta) branches
sb = BER["sigma_branch"]; b2 = BER["branch_2p67"]; pr = BER["promotion_1p7"]
ax[1].plot([r["theta"] for r in sb], [r["dstar"] for r in sb], "o-", color="crimson",
           label="sigma-CI (<- linear 2.445)")
ax[1].plot([r["theta"] for r in b2], [r["dstar"] for r in b2], "s-", color="teal",
           label="2nd CI branch (~2.67)")
ax[1].plot([r["theta"] for r in pr], [r["dstar"] for r in pr], "^-", color="purple",
           label="promoted CI (~1.7, bend-only)")
ax[1].axvline(180, color="k", ls=":", alpha=0.4)
ax[1].annotate("LINEAR\n(equilib.)", (180, 2.445), textcoords="offset points",
               xytext=(-46, -2), fontsize=8, color="crimson")
ax[1].set_xlabel("H-Be-H angle theta (deg)  [180 = linear]")
ax[1].set_ylabel("CI location d* (bohr)")
ax[1].set_title("CI branches vs bend  (all Berry = pi)")
ax[1].legend(loc="upper left", fontsize=8)
ax[1].invert_xaxis()

# inset: sigma/pi coupling promotes the 1.7 touch into a CI
axi = ax[1].inset_axes([0.58, 0.08, 0.38, 0.34])
axi.plot([r["theta"] for r in pr], [r["sigpi_coupling"] for r in pr], ".-", color="purple")
axi.set_title("sig/pi coupling", fontsize=7)
axi.set_xlabel("theta", fontsize=6); axi.tick_params(labelsize=6)
axi.axhline(0, color="grey", lw=0.5)

fig.suptitle("BeH2 three-center conical intersection OFF the linear axis  "
             f"(engine validated vs mpmath linear driver: |dEig|=1.7e-13)", fontsize=11)
fig.tight_layout()
fig.savefig("debug/plots/beh2_bending_ci.png", dpi=130)
print("wrote debug/plots/beh2_bending_ci.png")
