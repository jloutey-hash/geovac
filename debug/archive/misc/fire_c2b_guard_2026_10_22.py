"""Fire test for tests/test_paper59_c2b_elliptic.py (EC1-EC3): each planted wrong answer must make
the test aimed at it fail."""
import copy, json, sys
sys.path.insert(0, "tests")
import test_paper59_c2b_elliptic as G

def fires(fn):
    try:
        fn()
    except AssertionError:
        return True
    return False

orig = copy.deepcopy(G.DATA)
ladder = json.load(open("tests/data/paper59_ladder_2026_10_22.json"))["g32"]
plants = [
    ("wrong P (the (2,1,1,2) numerator)", "P", "16*x**2 - 88*x + 25",
     [G.test_family_and_j_invariant, G.test_c2b_is_the_image_of_the_elliptic_period_operator, G.test_source_solution_is_a_period_of_the_family]),
    ("wrong R (616 -> 612)", "R", "(4*x - 5)*(16*x**2 - 612*x + 25)",
     [G.test_family_and_j_invariant, G.test_source_solution_is_a_period_of_the_family]),
    ("wrong target (the chain's R2 in place of C2b)", "C2b_theta", ladder["chain_R2"],
     [G.test_c2b_is_the_image_of_the_elliptic_period_operator]),
    ("map numerators swapped", None, None, [G.test_c2b_is_the_image_of_the_elliptic_period_operator]),
]
for name, key, val, tests in plants:
    G.DATA.clear(); G.DATA.update(copy.deepcopy(orig))
    if key is None:
        G.DATA["N0"], G.DATA["N1"] = orig["N1"], orig["N0"]
    else:
        G.DATA[key] = val
    print(name, "->", [(t.__name__, fires(t)) for t in tests], flush=True)
G.DATA.clear(); G.DATA.update(copy.deepcopy(orig))
# a wrong period constant: compare with pi/sqrt(2) instead of pi/sqrt(3)
import mpmath as mp
real_sqrt = mp.sqrt
G.mp.sqrt = lambda v: real_sqrt(2) if v == 3 else real_sqrt(v)
print("wrong period constant ->", fires(G.test_source_solution_is_a_period_of_the_family), flush=True)
G.mp.sqrt = real_sqrt
print("clean ->", [(t.__name__, not fires(t)) for t in (G.test_family_and_j_invariant, G.test_c2b_is_the_image_of_the_elliptic_period_operator, G.test_source_solution_is_a_period_of_the_family)])
