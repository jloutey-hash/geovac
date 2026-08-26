"""
Track 4a / Probe 1 -- defect <-> continuum-coupling (polarizability) correlation.

Hypothesis under test (from WH7 "discreteness is compactness"): the max_n
geometry defect (v4.73.0 A/B/C: 100% orbital-basis, Layer-2/continuum-radial,
zero discrete-angular content) is a de-compactification cost, so |R_eq defect|
should track each system's continuum coupling -- proxied here by the static
dipole polarizability -- rather than bond topology / electron count.

Data source: debug/data/chem_error_atlas.md (15-row signed atlas, read in full).
Only rows with a REPORTED signed R_eq error qualify (atoms have no R_eq axis
per the atlas's own pattern-summary point 5). Polarizability values are AGENT
KNOWLEDGE, not looked up this sprint (no live literature search run) -- every
value is confidence-tagged and rounded to avoid false precision.

Method:
  1. Build the join table (system, method variant, |R_eq defect| %,
     N_electrons, N_centers, polarizability [a.u.], confidence tag).
  2. Spearman rank correlation of |defect| vs polarizability, N_electrons,
     N_centers -- on (a) the full row-level table (flagging pseudoreplication:
     one molecule's polarizability/N_e/N_c is identical across its method
     variants, so ties dominate) and (b) the system-level collapse (one row
     per distinct molecule, canonical variant only).
  3. Honest power statement: what |rho| would even be distinguishable from
     noise at n=4 (system-level) and n=9-11 (row-level)?
"""
from __future__ import annotations
import numpy as np
from scipy import stats

# --------------------------------------------------------------------- data
# Row-level join table. Source column = atlas row #. Defect = |R_eq error %|.
# polarizability_au: agent knowledge (isotropic static dipole polarizability,
# a.u. = bohr^3), ROUNDED, confidence-tagged. NOT a literature lookup this
# sprint -- see confidence column and the findings memo for caveats.
rows = [
    # system,    variant,                    |defect|%, N_e, N_c, pol_au, pol_confidence
    ("H2+",  "spectral Laguerre (row 1)",        0.25,   1,  2,   3.0,  "LOW (rough recall; lit. range ~2-4 a.u., strongly R/axis dependent, cation so << H2's 5.4)"),
    ("LiH",  "composed canonical l=2 (row 9)",   5.3,    4,  2,  26.0,  "MEDIUM (recalled range ~24-28 a.u. across sources)"),
    ("LiH",  "composed ab initio PK (row 9a)",   6.4,    4,  2,  26.0,  "MEDIUM (same molecule as above)"),
    ("LiH",  "composed fitted PK (row 9b)",      1.5,    4,  2,  26.0,  "MEDIUM (same molecule as above)"),
    ("LiH",  "balanced n_max=2 (memo)",          6.9,    4,  2,  26.0,  "MEDIUM (same molecule as above)"),
    ("LiH",  "balanced n_max=3 (memo)",          8.8,    4,  2,  26.0,  "MEDIUM (same molecule as above)"),
    ("LiH",  "composed l_max=3 (memo, drift)",  15.7,    4,  2,  26.0,  "MEDIUM (same molecule as above)"),
    ("LiH",  "composed l_max=4 (memo, drift)",  25.5,    4,  2,  26.0,  "MEDIUM (same molecule as above)"),
    ("LiH-4N","full 4e SO(12), l_max=2 (row 12)",63.5,   4,  2,  26.0,  "MEDIUM (same molecule; unbound D_e)"),
    ("BeH2", "composed l_max=2 (row 10)",       11.7,    6,  3,  15.0,  "LOW (no confident recall; order-of-magnitude estimate for a small linear 6e hydride, not a recalled literature figure)"),
    ("H2O",  "composed 5-block (row 11)",       19.4,   10,  3,   9.9,  "MEDIUM-HIGH (well-known experimental value ~1.47 A^3 = 9.9 a.u.; matches the value suggested in the task prompt)"),
]

print("=" * 100)
print("JOIN TABLE (row-level, n=%d)" % len(rows))
print("=" * 100)
hdr = f"{'system':8s} {'variant':34s} {'|defect|%':>9s} {'N_e':>4s} {'N_c':>4s} {'pol(au)':>8s}  confidence"
print(hdr)
for r in rows:
    print(f"{r[0]:8s} {r[1]:34s} {r[2]:9.2f} {r[3]:4d} {r[4]:4d} {r[5]:8.1f}  {r[6]}")

defect = np.array([r[2] for r in rows])
ne = np.array([r[3] for r in rows])
nc = np.array([r[4] for r in rows])
pol = np.array([r[5] for r in rows])

print("\n" + "=" * 100)
print("SPEARMAN CORRELATIONS -- ROW LEVEL (n=%d, WARNING: heavy ties -- 8/11 rows share" % len(rows))
print("LiH's single (N_e, N_c, pol) triple, so this is dominated by within-LiH method")
print("variance, not a genuine 11-point cross-system comparison)")
print("=" * 100)
for name, x in (("polarizability", pol), ("N_electrons", ne), ("N_centers", nc)):
    rho, p = stats.spearmanr(defect, x)
    print(f"  |defect| vs {name:15s}: rho={rho:+.3f}  p={p:.3f}  (n={len(rows)})")

# ---------------------------------------------------- system-level collapse
print("\n" + "=" * 100)
print("SYSTEM-LEVEL COLLAPSE (one row per distinct molecule, canonical variant only)")
print("=" * 100)
canonical = {
    "H2+": (0.25, 1, 2, 3.0),
    "LiH": (5.3, 4, 2, 26.0),        # canonical composed l=2 headline
    "BeH2": (11.7, 6, 3, 15.0),
    "H2O": (19.4, 10, 3, 9.9),
}
sysnames = list(canonical.keys())
sdefect = np.array([canonical[s][0] for s in sysnames])
sne = np.array([canonical[s][1] for s in sysnames])
snc = np.array([canonical[s][2] for s in sysnames])
spol = np.array([canonical[s][3] for s in sysnames])
print(f"{'system':8s} {'|defect|%':>9s} {'N_e':>4s} {'N_c':>4s} {'pol(au)':>8s}")
for s in sysnames:
    d, n1, n2, p = canonical[s]
    print(f"{s:8s} {d:9.2f} {n1:4d} {n2:4d} {p:8.1f}")
for name, x in (("polarizability", spol), ("N_electrons", sne), ("N_centers", snc)):
    rho, p = stats.spearmanr(sdefect, x)
    print(f"  |defect| vs {name:15s}: rho={rho:+.3f}  p={p:.3f}  (n=4)")

# also: including LiH-4N as a 5th "system" (same molecule, unbound extreme) to
# see how sensitive n=4/5 collapse is to that one inclusion decision
print("\n[sensitivity] adding LiH-4N (63.5%, same N_e/N_c/pol as LiH) as a 5th point:")
sdefect5 = np.append(sdefect, 63.5)
sne5 = np.append(sne, 4)
snc5 = np.append(snc, 2)
spol5 = np.append(spol, 26.0)
for name, x in (("polarizability", spol5), ("N_electrons", sne5), ("N_centers", snc5)):
    rho, p = stats.spearmanr(sdefect5, x)
    print(f"  |defect| vs {name:15s}: rho={rho:+.3f}  p={p:.3f}  (n=5)")

# --------------------------------------------------------------- power stmt
print("\n" + "=" * 100)
print("HONEST POWER ANALYSIS")
print("=" * 100)
print("Normal-approximation critical |rho| for two-tailed alpha=0.05: "
      "rho_crit ~= 1.96/sqrt(n-1)")
for n in (4, 5, 9, 11):
    if n - 1 > 0:
        rc = 1.96 / np.sqrt(n - 1)
        feasible = "IMPOSSIBLE (exceeds 1.0 -- no rho, even +-1, reaches significance)" if rc > 1.0 else f"rho_crit ~= {rc:.3f}"
        print(f"  n={n:2d}: {feasible}")

print("\nStandard Spearman critical-value tables (two-tailed alpha=0.05) confirm:")
print("  n=9  -> rho_crit ~ 0.683 ;  n=11 -> rho_crit ~ 0.618 ;  n=4 -> undefined/unreachable.")
