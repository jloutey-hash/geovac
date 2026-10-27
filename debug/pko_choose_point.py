"""PK record, phase 1: choose the blind second parameter point (a, b, g1, g2).
Requirement: every candidate singular value of the kite or its system is distinct and nonzero, so that each
clause of the cut rule is tested separately.  Candidates (all as -z):
  opposite cuts  (a+-b)^2, (g1+-g2)^2, (a+-g2)^2, (b+-g1)^2       (predicted genuine)
  adjacent cuts  (a+-g1)^2, (b+-g2)^2                             (predicted ordinary)
  sums of squares of every pair, and every single m^2             (seen only in the system, or nowhere)
Also no mass equal to another, and none equal to a sum or difference of two others.
Among admissible points with a = 1 and the others in (1/2)Z or (1/3)Z up to 4, print those with the smallest
largest candidate.  Pure arithmetic; scores nothing."""
from fractions import Fraction as F
from itertools import product, combinations


def cands(a, b, g1, g2):
    m = {"a": a, "b": b, "g1": g1, "g2": g2}
    out = {}
    for (x, y), tag in [(("a", "b"), "opp"), (("g1", "g2"), "opp"), (("a", "g2"), "opp"), (("b", "g1"), "opp"),
                        (("a", "g1"), "adj"), (("b", "g2"), "adj")]:
        out[f"{tag} ({x}+{y})^2"] = (m[x] + m[y]) ** 2
        out[f"{tag} ({x}-{y})^2"] = (m[x] - m[y]) ** 2
    for x, y in combinations(m, 2):
        out[f"sq {x}^2+{y}^2"] = m[x] ** 2 + m[y] ** 2
    for x in m:
        out[f"one {x}^2"] = m[x] ** 2
    return out


def admissible(a, b, g1, g2):
    ms = [a, b, g1, g2]
    if len(set(ms)) < 4:
        return False
    for x, y, w in product(ms, repeat=3):
        if len({x, y, w}) == 3 and (x + y == w or abs(x - y) == w):
            return False
    v = list(cands(a, b, g1, g2).values())
    return 0 not in v and len(set(v)) == len(v)


grid = sorted({F(i, 2) for i in range(1, 9)} | {F(i, 3) for i in range(1, 13)})
best = []
for b, g1, g2 in product(grid, repeat=3):
    a = F(1)
    if admissible(a, b, g1, g2):
        best.append((max(cands(a, b, g1, g2).values()), max(x.denominator for x in (b, g1, g2)), (a, b, g1, g2)))
best.sort()
for mx, den, pt in best[:12]:
    print(mx, den, [str(x) for x in pt])
