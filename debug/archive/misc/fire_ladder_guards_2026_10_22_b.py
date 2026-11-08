"""Path guard: fixed middle test passes, and each slow test FAILS on the wrong answer it names."""
import copy, sys, time
sys.path.insert(0, "tests")
import sympy as sp
import test_paper59_four_density_path as P

def fires(fn):
    try:
        fn()
    except AssertionError:
        return True
    return False

orig = copy.deepcopy(P.DATA)
def reset():
    P.DATA.clear(); P.DATA.update(copy.deepcopy(orig))
t0 = time.time()
print("middle passes (fixed lookup):", not fires(P.test_path_middle_pieces_are_legendre), "%.0fs" % (time.time() - t0), flush=True)
# middle, wrong answer 1: Z2 and Z2' swapped (Z2 claimed as the right factor)
P.DATA["Z2p_Q"], P.DATA["Z2_Q"] = orig["Z2_Q"], orig["Z2p_Q"]
print("middle fires on Z2/Z2' swapped:", fires(P.test_path_middle_pieces_are_legendre), flush=True); reset()
# middle, wrong answer 2: Z2 map replaced by the Z2' map ('both pieces are the same Legendre image')
P.DATA["maps"]["Legendre -> Z2"] = orig["maps"]["Legendre -> Z2'"]
print("middle fires on a wrong Z2 map:", fires(P.test_path_middle_pieces_are_legendre), flush=True); reset()
# f4: wrong stored map (one coefficient off)
k = [k for k in orig["maps"] if k.startswith("F_4")][0]
m = list(orig["maps"][k]); m[2] = str(sp.sympify(m[2], locals={"t": P.t}) * sp.Rational(1001, 1000))
P.DATA["maps"][k] = m
t0 = time.time()
print("f4 fires on a wrong map:", fires(P.test_path_f4_piece_is_exact), "%.0fs" % (time.time() - t0), flush=True); reset()
# numerics: compute the grid once, reuse it
t0 = time.time()
G = P._grids()
print("grid %.0fs" % (time.time() - t0), flush=True)
P._grids = lambda *a, **k: G
print("numeric passes on the cached grid:", not fires(P.test_w12_annihilates_path_and_needs_f4), flush=True)
# wrong answer: 'the path sits below the F_4 piece' -> claim Y4 is trivial: W12 := R8 (Y4 replaced by identity)
P.DATA["Y4_Q"] = ["1"]
print("numeric fires on 'R8 already kills P':", fires(P.test_w12_annihilates_path_and_needs_f4), flush=True); reset()
# wrong integrand: the star's weight in place of the path's
P._grids = lambda *a, **k: (G[1], G[0])
print("numeric fires on the star's integrand:", fires(P.test_w12_annihilates_path_and_needs_f4), flush=True)
