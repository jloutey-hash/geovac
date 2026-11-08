import copy, sys, time
sys.path.insert(0, "tests")
import test_paper59_four_density_path as P

def fires(fn):
    try:
        fn()
    except AssertionError:
        return True
    return False

orig = copy.deepcopy(P.DATA)
G = P._grids()
P._grids = lambda *a, **k: G
P.DATA["Y4_Q"] = ["u**2"]          # a wrong top piece: theta^2 in place of the F_4-equivalent Y4
print("numeric fires on a wrong top piece (theta^2 for Y4):", fires(P.test_w12_annihilates_path_and_needs_f4), flush=True)
P.DATA.clear(); P.DATA.update(copy.deepcopy(orig))
P._grids = lambda *a, **k: (G[1], G[0])
print("numeric fires on the star's integrand:", fires(P.test_w12_annihilates_path_and_needs_f4), flush=True)
