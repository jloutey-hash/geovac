"""Fire tests: each guard must FAIL on the wrong answer it names."""
import copy, sys
sys.path.insert(0, "tests")
import sympy as sp
import test_paper59_correlated_triangle_g32 as G
import test_paper59_four_density_path as P

def fires(fn):
    try:
        fn()
    except AssertionError:
        return True
    return False

orig = copy.deepcopy(G.DATA)
# 1. 'Y sits below the K3 piece': pretend the right factor below F is L9 itself (it kills Y)
G._theta_to_Q = lambda cs: G._L9(G.DR.u)
print("GT Y-test fires on 'K6 Y = 0':", fires(G.test_g32_y_needs_the_k3_piece))
# 2. wrong order-three piece: F replaced by K6's top three theta coefficients is not F_3-equivalent
G.DATA["factor_F"] = orig["chain_V4"]
print("GT map-test fires on a wrong F:", fires(G.test_g32_k3_piece_is_exact))
# 3. 'a piece above order three': claim L9 = F o V5 (order-five quotient)
G.DATA.clear(); G.DATA.update(copy.deepcopy(orig)); G.DATA["chain_K6"] = orig["chain_V5"]
print("GT series-test fires on a wrong chain:", fires(G.test_g32_composition_series))
# 4. path bottom: R2 replaced by Legendre twisted wrongly (theta-form of R3's top two not R2)
P.DATA["R2_Q"] = P.DATA["Z2p_Q"]
print("path bottom-test fires on a wrong R2:", fires(P.test_path_bottom_pieces))
