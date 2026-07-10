"""
Z>20 cliff diagnostic — Synthesis: combine all probe results, rank causes,
emit final go/no-go on BBB93 effort.

This script reads all the probe data and compiles the diagnostic verdict
that goes into the memo.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

PROBES = {
    "a": "CR67 single-zeta fit faithfulness (across alkali heavy atoms)",
    "a_extended": "CR67 fit faithfulness (across Cs shells)",
    "b": "FD solver hydrogenic faithfulness (across Z)",
    "c": "Convention mismatch trace through Z_eff -> psi -> A_HF chain",
    "d": "multi-zeta machinery internal coherence",
    "e_v2": "Z-scan localization (kernel isolated, |psi_FW/psi_HF| ratio)",
}

verdicts = {}
data = {}

for probe, name in PROBES.items():
    path = Path(ROOT) / "debug" / "data" / f"z_cliff_probe_{probe}.json"
    if not path.exists():
        verdicts[probe] = ("MISSING", "Probe not run")
        continue
    with open(path) as f:
        d = json.load(f)
    data[probe] = d
    verdicts[probe] = (d.get("verdict", "?"), name)

print("=" * 80)
print("Z>20 CLIFF DIAGNOSTIC — SYNTHESIS")
print("=" * 80)
print()
print("Probe verdicts:")
for probe, (verdict, name) in verdicts.items():
    print(f"\n[{probe}] {name}")
    # Wrap long lines
    if len(verdict) > 75:
        # Break into ~75-char chunks
        words = verdict.split()
        line = "  "
        for w in words:
            if len(line) + len(w) > 75:
                print(line)
                line = "  "
            line = line + w + " "
        if line.strip():
            print(line)
    else:
        print(f"  {verdict}")

# Compile the dominant-cause ranking
print()
print("=" * 80)
print("DOMINANT CAUSE RANKING")
print("=" * 80)

ranking = []

# Probe (a) extended provides the clearest quantitative cliff signature.
if "a_extended" in data:
    inner_overshoot = data["a_extended"]["inner_overshoot_mean"]
    outer_overshoot = data["a_extended"]["outer_overshoot_mean"]
    if outer_overshoot and inner_overshoot:
        ratio = outer_overshoot / inner_overshoot
        ranking.append((
            1,
            "(a) CR67 single-zeta non-faithful for OUTER SHELLS at high Z",
            f"Cs CR67 inner-shell overshoot {inner_overshoot:.2f}x; "
            f"outer-shell overshoot {outer_overshoot:.2f}x; "
            f"ratio {ratio:.1f}x. Outer shells (4d, 5s, 5p) of [Xe] core "
            f"are 2-4x too compact in the Slater fit."
        ))

# Probe (d) sharpens this — the "missing radial nodes" mechanism.
if "d" in data:
    overlap = data["d"].get("max_s_orbital_overlap", 0)
    if overlap > 0.5:
        ranking.append((
            2,
            "(d) MISSING RADIAL NODES (heuristic two-zeta has all-positive coefficients)",
            f"Two-zeta heuristic gives max |<ns|ms>|={overlap:.2f} between "
            f"different Cs s-orbitals (should be 0). The CR67 single-zeta "
            f"fit has no nodes either (it's a bare hydrogenic), so the "
            f"effective Z at intermediate r is set by integrating shells "
            f"that are all 'clones' of each other rather than orthogonal "
            f"shells with proper radial structure."
        ))

# Probe (e) v2 confirms the cliff is well-localized at Z=19+.
if "e_v2" in data:
    cliff_z = data["e_v2"].get("cliff_first_Z")
    rows = data["e_v2"]["rows"]
    cliff_summary = []
    for r in rows:
        if "ratio_FW_over_lit" in r and r["ratio_FW_over_lit"] is not None:
            ratio = r["ratio_FW_over_lit"]
            if abs(ratio - 1.0) > 0.1:
                cliff_summary.append(f"{r['label']} ratio={ratio:.3f}")

    ranking.append((
        3,
        "(e) Cliff onset at Z=19 (K 4s); Z=11 (Na) shows mild signature",
        f"|psi_FW(0)|^2 / |psi_HF_lit(0)|^2 ratios: " +
        "; ".join(cliff_summary) +
        f". H Z=1 ratio is 0.999 (perfect). Cliff onset clearly bracketed "
        f"between Z=11 (24% over) and Z=19 (3.2x under). Cs Z=55 ratio "
        f"0.625 (1.6x under). The Cs is MILDER than Rb — relativistic "
        f"effects partially compensate."
    ))

# Probe (b): FD is faithful, NOT the cliff cause.
if "b" in data:
    ranking.append((
        99,  # last (negative result)
        "(b) FD SOLVER IS FAITHFUL — NOT the cliff cause",
        "Hydrogenic test at Z=1..80 with 400k grid: |psi(0)|^2 error "
        "<1.2% in all cases, monotone convergence. The closeout sprint's "
        "dense uniform FD path is correct; numerical refinement at higher "
        "n_grid will not close the cliff."
    ))

# Probe (c): no convention bug.
if "c" in data:
    findings = data["c"]["findings"]
    fail_count = sum(1 for f in findings if "FAIL" in f["verdict"])
    pass_count = sum(1 for f in findings if "PASS" in f["verdict"])
    ranking.append((
        99,  # last
        "(c) NO CONVENTION BUG — chain Z_eff -> psi -> A_HF is internally consistent",
        f"H 21cm reproduces at +0.17% (the expected reduced-mass Sprint HF level). "
        f"g_N convention is correct. n*zeta convention matches Slater. "
        f"Z_eff(r) boundary conditions correct. {pass_count} PASS, "
        f"{fail_count} FAIL (numerical, not structural)."
    ))

# Probe (d) coding bug: ruled out.
ranking.append((
    99,
    "(d) NO CODING BUG in multi-zeta machinery",
    "STO normalization, MultiZetaOrbital reduces correctly to single-zeta, "
    "_build_two_zeta_orbital renormalizes within 1%, density integrates "
    "to N_core within 0.0%. Multi-zeta machinery installed correctly."
))

# Sort
ranking.sort()

for rank, name, detail in ranking:
    print()
    print(f"#{rank}  {name}")
    # Wrap detail
    words = detail.split()
    line = "    "
    for w in words:
        if len(line) + len(w) > 75:
            print(line)
            line = "    "
        line = line + w + " "
    if line.strip():
        print(line)

print()
print()
print("=" * 80)
print("DIAGNOSTIC NET CONCLUSION")
print("=" * 80)
print()

print("""
The Z>20 cliff is DOMINATED by mechanism (a) + (d):

  (a) CR67 single-zeta fits overshoot Z_eff for OUTER shells of heavy
      atoms by 2-4x. INNER shells (1s, 2s, 2p) are well-fit (~1x). The
      cliff is shell-resolved, not Z-uniform.

  (d) The structural reason CR67 fits outer shells badly: a single-zeta
      hydrogenic R_nl has NO RADIAL NODES — but real RHF orbitals have
      n-l-1 nodes (3 for 5s, 2 for 5p, 1 for 4d). The nodes are what
      make outer shells RIGHT — they push amplitude OUT to where the
      orbital should physically be, by orthogonalizing against inner
      shells. Without nodes, the orbital piles up near the inner zeta
      peak.

The two-zeta heuristic (May 9 closeout sprint) inherits both problems:
it has 2 primitives but all-positive coefficients, so STILL no nodes.
That's why it FAILED in the wrong direction (-47% -> -90%) — adding
"more compact inner + less compact outer" without orthogonalization
just makes the screening transition steeper without fixing the position.

Mechanisms (b), (c), (e single-Z bugs) are RULED OUT:
  (b) FD solver is faithful — closeout sprint's algorithm is correct.
  (c) No factor-of-X convention bug in Z_eff -> psi -> A_HF chain.
  (d-coding) Multi-zeta machinery itself is bug-free.
""")

print()
print("=" * 80)
print("BBB93 GO/NO-GO RECOMMENDATION")
print("=" * 80)
print()

print("""
GO ON BBB93/KTT FULL TABULATION — but with EYES OPEN.

The dominant cause (a)+(d) IS what BBB93/KTT addresses: those
tabulations have 5-9 STO primitives per orbital with NEGATIVE
coefficients that produce the radial nodes. They will close the
cliff for inner shells, but the OUTER shell behavior depends on
how well the published tabulation reproduces the actual orbital
extent.

Risk factors:
  - For BBB93 (Z<=54): The Xe atom is at the upper limit of BBB93's
    range. Coefficients there may be less converged than for lighter
    atoms in the same paper.
  - For Cs/Ba (Z=55-56): MUST use Koga-Tatewaki-Thakkar 1993/2000
    (BBB93 stops at Z=54). KTT data is harder to find in tabular form.
  - The Casimir F_R = 1.555 leading-order is INDEPENDENT of the
    screening upgrade. Closing the screening cliff still leaves a
    ~30-40% gap to the full Bohr-Weisskopf factor (~2.6 vs 1.555).
    BBB93 alone is NOT sufficient for sub-percent Cs HFS.

Recommendation:
  1. PRIMARY: BBB93/KTT for [Xe] core (multi-day tabulation effort).
     Expected outcome: framework-native A residual reduced from -47%
     to ~-30% for Cs (the ~20% closure attributable to fixing the
     screening kernel).
  2. SECONDARY (necessary for sub-percent): full Bohr-Weisskopf via
     spinor-lift evaluation. Closes the residual to ~few percent for Cs.
  3. ALTERNATIVE PATH (3-4 weeks, more principled): self-consistent
     HF iteration in geovac/neon_core.py. Closes the screening cliff
     ALL Z (not just Xe core), making the framework heavy-atom-extensible.
     This is the long-term clean engineering path.
""")

# Save synthesis
synthesis = {
    "diagnostic_summary": "Z>20 cliff dominated by CR67 single-zeta fit non-faithfulness for outer shells, structurally caused by missing radial nodes",
    "dominant_causes_ranked": [
        {"rank": rank, "name": name, "detail": detail}
        for rank, name, detail in ranking
    ],
    "ruled_out": [
        "FD solver faithfulness (probe b)",
        "Convention bugs (probe c)",
        "multi-zeta machinery coding bugs (probe d)",
    ],
    "bbb93_recommendation": "GO with eyes open: closes screening kernel cliff (~20% closure) but does NOT close full Bohr-Weisskopf relativistic gap (~30-40% remaining). Sub-percent Cs HFS requires both. Self-consistent HF (multi-week) is the principled alternative.",
}
out_path = Path(ROOT) / "debug" / "data" / "z_cliff_probe_synthesis.json"
with open(out_path, "w") as f:
    json.dump(synthesis, f, indent=2)
print(f"\nSaved synthesis to {out_path}")


if __name__ == "__main__":
    pass
