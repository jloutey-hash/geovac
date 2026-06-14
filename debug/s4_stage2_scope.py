"""Stage-2 scoping: classify the S^(4) multiple-t content by weight & depth,
and determine the assembly/identification track structure.

Reads the census + tables from debug/data/s4_decomp_tables.json and the
o-space relation coefficients, and reports:
  - weight x depth histogram of all distinct t-values
  - per-weight count of t-values that need identification
  - which weights have depth>=3 content (the depth-verdict frontier)
  - the maximum weight (ceiling) and the level-1 vs level-2 split heuristic
"""
import json
from collections import defaultdict
from pathlib import Path

DATA = Path(__file__).parent / "data"
d = json.loads((DATA / "s4_decomp_tables.json").read_text())

# all distinct t-value atoms across the cores
atoms = set()
for nm, tab in d["tables"].items():
    for ks in tab:
        s = eval(ks)  # noqa: S307
        if s[0] == "t":
            atoms.add(s[1])

# weight = sum of args ; depth = len
by_wd = defaultdict(int)
by_w = defaultdict(list)
for a in atoms:
    w, dep = sum(a), len(a)
    by_wd[(w, dep)] += 1
    by_w[w].append(a)

print("=== weight x depth histogram (distinct multiple-t atoms) ===")
weights = sorted(by_w)
depths = sorted({len(a) for a in atoms})
hdr = "  w\\d |" + "".join("  d=%d" % dd for dd in depths) + "   total"
print(hdr)
for w in weights:
    row = "  %3d |" % w
    tot = 0
    for dd in depths:
        c = by_wd.get((w, dd), 0)
        tot += c
        row += "  %3d" % c if c else "    ."
    row += "    %3d" % tot
    print(row)
print("  total atoms:", len(atoms))

print("\n=== depth>=3 content by weight (the depth-verdict frontier) ===")
for w in weights:
    d3 = [a for a in by_w[w] if len(a) >= 3]
    if d3:
        print("  w=%2d : %2d depth>=3 atoms; sample %s"
              % (w, len(d3), sorted(d3)[:4]))

print("\n=== ceiling & realized-depth question ===")
wmax = max(weights)
d4 = [a for a in atoms if len(a) == 4]
print("  weight ceiling:", wmax, " (odd-weight tower predicted 5..13)")
print("  depth-4 atoms:", len(d4), " max weight among them:",
      max(sum(a) for a in d4))
print("  -> realized-depth-<=3 test = whether all depth-4 atoms reduce to")
print("     depth<=3 in the assembly (k=3 pattern: depth<=k-1).")

# rough level-1 reducibility heuristic: a t-value is "level-1 trivially
# reducible" if it has <=1 non-trailing slot above the trailing block, etc.
# Here just report the trailing-1 census (the hard-to-sum, easy-to-identify set)
tr = [a for a in atoms if a[-1] == 1]
print("\n=== trailing-1 atoms (numerically hard, symbolically the target) ===")
print("  count:", len(tr))
for dep in depths:
    c = [a for a in tr if len(a) == dep]
    if c:
        print("  depth %d : %d  (weights %s)"
              % (dep, len(c), sorted({sum(a) for a in c})))
