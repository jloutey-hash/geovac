"""HG2: rederive Delta(u1) without LM Eqs. 32-33, as the D=3 Baikov polynomial (Gram determinant of the three
loop momenta) of the tetrahedron G(1,1,1,1,1,1) on its maximal cut, then compare with LM's transcribed form.

Propagators (LM Eq. 13): k1^2+u1^2, k2^2+u2^2, k3^2+u3^2, k32^2+w1^2, k13^2+w2^2, k21^2+w3^2 (Euclidean),
k_ij = k_i - k_j. On the maximal cut k_i.k_i = -u_i^2 and (k_i-k_j)^2 = -w_ij^2.
"""
import sympy as sp

u1, u2, u3, w1, w2, w3 = sp.symbols('u1 u2 u3 w1 w2 w3')
sq = {1: -u1**2, 2: -u2**2, 3: -u3**2}
wij = {(3, 2): w1, (1, 3): w2, (2, 1): w3}


def kk(i, j):
    if i == j:
        return sq[i]
    m = wij.get((i, j)) or wij.get((j, i))
    return (sq[i] + sq[j] + m**2)/2      # (k_i - k_j)^2 = -m^2


G = sp.Matrix(3, 3, lambda i, j: kk(i + 1, j + 1))
gram = sp.expand(G.det())

s1 = (-u2**2*w1**2 - u2**2*w2**2 + u2**2*w3**2 - u3**2*w1**2 + u3**2*w2**2 - u3**2*w3**2
      + w1**4 - w1**2*w2**2 - w1**2*w3**2)
s2 = (u2**4*w2**2 + u2**2*u3**2*w1**2 - u2**2*u3**2*w2**2 - u2**2*u3**2*w3**2 - u2**2*w1**2*w2**2
      + u2**2*w2**4 - u2**2*w2**2*w3**2 + u3**4*w3**2 - u3**2*w1**2*w3**2 - u3**2*w2**2*w3**2
      + u3**2*w3**4 + w1**2*w2**2*w3**2)
Delta = sp.expand(w1**2*u1**4 + s1*u1**2 + s2)

ratio = sp.simplify(gram/Delta)
print("Gram det (maximal cut) =", sp.factor(gram))
print("Gram / Delta_LM =", ratio)
print("HG2 match up to a constant:", ratio.free_symbols == set())
# which variable pairs with u1 in LM's labelling? also report the u1-polynomial structure of the Gram det
print("Gram det as polynomial in u1:", sp.Poly(gram, u1).all_coeffs()[:1], "... degree", sp.Poly(gram, u1).degree())
