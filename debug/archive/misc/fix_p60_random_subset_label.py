"""A sign confusion of mine, caught by the guard rather than by me.

`debug/p60_avery_102_probe.py` computed `best = max(best, locked(...))` over
ENERGIES.  Maximising E selects the LEAST-bound subset, so the reported figure
is the WORST of the 200 random 102-subsets, not the best.  Measured under the
driver's own seed:  best 13.862 mHa, worst 903.404 mHa (58% of random subsets
exclude the dominant 1s^2 configuration).

The argument is unaffected in direction -- random selection is still far worse
than weighted selection (13.9 vs 7.25 mHa) and both are far above Avery's 1.22
-- but the number was mislabelled in the paper and the registry.  Fixed in all
three places, and the driver now reports both ends.
"""
import io

DRV = "debug/p60_avery_102_probe.py"
PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
REG = "debug/qa/numeric_registry.py"

# ------------------------------------------------------------------- driver
d = io.open(DRV, encoding="utf-8").read()
D_OLD = """rng = np.random.default_rng(0)
best = -1e9
for _ in range(200):
    s = np.sort(rng.choice(K, size=nsel, replace=False))
    best = max(best, locked(M[np.ix_(s, s)]))
print("200 random %d-subsets: best locked E=%.7f  gap=%7.3f mHa" % (nsel, best, (best - EXACT) * 1000))"""
D_NEW = """rng = np.random.default_rng(0)
# NOTE the sign:  E is negative and deeper binding is LOWER, so the BEST random
# subset is the one with the MINIMUM energy.  An earlier version of this driver
# took max() and therefore reported the worst subset as the best.
energies = [locked(M[np.ix_(np.sort(rng.choice(K, size=nsel, replace=False)),
                            np.sort(rng.choice(K, size=nsel, replace=False)))])
            for _ in range(0)]  # placeholder, replaced below
energies = []
for _ in range(200):
    s = np.sort(rng.choice(K, size=nsel, replace=False))
    energies.append(locked(M[np.ix_(s, s)]))
best, worst = min(energies), max(energies)
print("200 random %d-subsets: BEST gap=%7.3f mHa   WORST gap=%7.3f mHa"
      % (nsel, (best - EXACT) * 1000, (worst - EXACT) * 1000))"""
assert D_OLD in d, "driver locus not found"
d = d.replace(D_OLD, D_NEW, 1)
d = d.replace('best_random=best,', 'best_random=best, worst_random=worst,', 1)
io.open(DRV, "w", encoding="utf-8").write(d)
print("driver: min/max corrected, both ends reported")

# -------------------------------------------------------------------- paper
t = io.open(PAP, encoding="utf-8").read()
P_OLD = r"""while $200$ random $102$-subsets reach only $903$~mHa."""
P_NEW = r"""while of $200$ random $102$-subsets the best reaches only
$13.9$~mHa and the worst $903$~mHa --- $58\%$ of them omit the dominant $1s^{2}$
configuration entirely."""
assert P_OLD in t, "paper locus not found"
t = t.replace(P_OLD, P_NEW, 1)
io.open(PAP, "w", encoding="utf-8").write(t)
print("paper: random-subset figure relabelled (best vs worst)")

# ----------------------------------------------------------------- registry
r = io.open(REG, encoding="utf-8").read()
R_OLD = '''"random 102-subsets reach only 903 mHa. The cited Avery figure "'''
R_OLD2 = '''        aliases={7.057: "the full K=244 pool", 903.4: "best of 200 random 102-subsets"}),'''
R_NEW2 = ('''        aliases={7.057: "the full K=244 pool",\n'''
          '''                 13.862: "BEST of 200 random 102-subsets",\n'''
          '''                 903.4: "WORST of 200 random 102-subsets -- the driver's max() over\\n'''
          '''                         energies selected the least-bound subset;  mislabelled as\\n'''
          '''                         'best' until 2026-09-08"}),''')
assert R_OLD2 in r, "registry alias locus not found"
r = r.replace(R_OLD2, R_NEW2, 1)
if R_OLD in r:
    r = r.replace(R_OLD, '''"random 102-subsets reach 13.9 mHa at best. The cited Avery figure "''', 1)
io.open(REG, "w", encoding="utf-8").write(r)
print("registry: aliases relabelled with the sign note")
