"""
Numerical verification of the Casimir triangle inequality on SU(3).

Statement: For irreps π, π' with highest weights λ, λ', and any σ in
the decomposition of π ⊗ π'*, does

    |C(π) - C(π')| <= C(σ)

hold?

For SU(3) with highest weight (p,q):
    C(p,q) = (1/3)(p^2 + q^2 + p*q + 3p + 3q)
    dim(p,q) = (p+1)(q+1)(p+q+2)/2
    Dual: (p,q)^* = (q,p)

Tensor product (p1,q1) ⊗ (p2,q2) decomposition can be computed by direct
character computation (Weyl character formula) and product expansion.

We use sympy to enumerate weight vectors and compute character products
for verification on small examples.
"""

import sympy as sp
from sympy import Rational, symbols, expand, Poly, S, simplify
from itertools import product
from typing import List, Tuple, Dict


# =====================================================================
# SU(3) Casimir
# =====================================================================

def casimir_su3(p: int, q: int) -> Rational:
    """
    Quadratic Casimir of SU(3) irrep with Dynkin label (p, q).

    Formula: C(p,q) = (1/3)(p^2 + q^2 + p*q + 3p + 3q)

    This is the standard normalization with C(adjoint) = N_c = 3.
    Verify: (1,1) gives C = (1/3)(1+1+1+3+3) = 9/3 = 3. ✓
    """
    return Rational(p * p + q * q + p * q + 3 * p + 3 * q, 3)


def dim_su3(p: int, q: int) -> int:
    """Dimension of SU(3) irrep (p,q): (p+1)(q+1)(p+q+2)/2"""
    return (p + 1) * (q + 1) * (p + q + 2) // 2


# =====================================================================
# SU(3) tensor product decomposition via character computation
# =====================================================================
#
# We use the direct method: Weyl character formula gives
#   chi_lambda * chi_mu = sum_nu N^nu_{lambda,mu} chi_nu
#
# For SU(3), it's easier to use the explicit weight diagram method:
# expand chi_lambda in terms of characters chi_nu by enumerating
# weights with multiplicities.
#
# We use a known algorithmic approach: Brauer-Klimyk (Racah-Speiser)
# algorithm, which is finite and exact.
#
# For our verification needs, we'll use the simpler approach of using
# sympy's representation theory utilities or implement Brauer-Klimyk
# directly for small (p,q).
# =====================================================================


def weights_of_irrep_su3(p: int, q: int) -> List[Tuple[int, int]]:
    """
    Return all weights (a, b) of the SU(3) irrep (p,q) WITH multiplicity.

    We work in the basis of fundamental weights omega_1, omega_2.
    Weight = a*omega_1 + b*omega_2, encoded as (a, b).

    The highest weight is (p, q).

    For SU(3) the multiplicity is given by Freudenthal-like formulas;
    for small reps we can enumerate by descent through the weight lattice.

    Algorithm: starting from (p, q), apply lowering operators (subtract
    simple roots a_1 = 2*omega_1 - omega_2 and a_2 = -omega_1 + 2*omega_2),
    keep all weights in the convex hull of the Weyl orbit of (p,q).

    Multiplicity by Freudenthal: works but is recursive. For small cases
    we can use known formulas:
       - Each layer of the irrep is hexagonal (or triangular at corners)
       - On each shell, multiplicity equals shell index from outermost

    For our purposes we use the explicit Freudenthal recursion.
    """
    # We compute via Freudenthal's formula
    # Convert to alpha-basis (simple root basis)
    # alpha_1 = (2, -1), alpha_2 = (-1, 2) in omega-basis
    # So weight w in omega-basis can be lowered by alpha_1: (a-2, b+1)
    # and by alpha_2: (a+1, b-2)

    # Highest weight
    hw = (p, q)

    # Generate all weights by descending (BFS)
    # We keep going as long as we stay in the polytope (convex hull of W-orbit of hw)
    # Equivalently: a weight w is in the irrep iff w lies in W*hw + Q-, i.e., is
    # below hw in the dominance order on the W-orbit.

    # Use the simple criterion: weight (a, b) is in the irrep (p,q) iff it's in
    # the convex hull of the Weyl orbit. For SU(3), Weyl group is S_3 acting on
    # weights; orbit of (p,q) has up to 6 elements. The convex hull is a hexagon
    # (or triangle if p=0 or q=0).

    # Easier approach: rho = (1,1), compute weight orbit, then enumerate lattice
    # points in convex hull.

    # Convert weight (a,b) to standard 3D form (l1, l2, l3) with l1+l2+l3=0:
    # omega_1 corresponds to (2/3, -1/3, -1/3), omega_2 to (1/3, 1/3, -2/3)
    # So weight = a*omega_1 + b*omega_2 = ((2a+b)/3, (-a+b)/3, (-a-2b)/3)
    # which has sum 0. ✓

    def to_3d(a, b):
        return (Rational(2*a + b, 3), Rational(-a + b, 3), Rational(-a - 2*b, 3))

    def from_3d(v):
        # v = (l1, l2, l3); a = l1 - l2, b = l2 - l3
        return (v[0] - v[1], v[1] - v[2])

    # Weyl group S_3 orbit of hw in 3D
    hw_3d = to_3d(p, q)
    from itertools import permutations
    orbit_3d = set()
    for perm in permutations(range(3)):
        wp = tuple(hw_3d[perm[i]] for i in range(3))
        orbit_3d.add(wp)

    # Find lattice points (a, b) in omega-basis such that to_3d(a,b) is in convex hull
    # Bounds: |a|, |b| <= p + q
    bound = p + q

    weights_with_mult = {}

    # The set of weights of (p,q) consists of all weights mu with mu <= (p,q) in the
    # dominance order AND mu in the convex hull of the Weyl orbit.
    # Equivalent: weights of (p,q) = (W-orbit hull) ∩ (hw + root lattice).
    # Root lattice for SU(3): generated by alpha_1, alpha_2.
    # alpha_1 = (2, -1) in omega-basis, alpha_2 = (-1, 2) in omega-basis.

    # Both alpha_1 and alpha_2 have integer coords in omega basis. So all weights
    # mu with mu = hw - n1*alpha_1 - n2*alpha_2 for n1, n2 >= 0.

    # We use Freudenthal recursion to compute multiplicities.
    # Multiplicity of mu in V_lambda = mult_freudenthal(mu, lambda)

    # Freudenthal:
    #   ((|lambda + rho|^2 - |mu + rho|^2)) * mult(mu) =
    #     2 * sum_{alpha > 0} sum_{k >= 1} (mu + k*alpha, alpha) * mult(mu + k*alpha)

    # For SU(3): positive roots are alpha_1 = (2,-1), alpha_2 = (-1,2), alpha_1+alpha_2 = (1,1) in omega

    # Inner product on omega-basis: requires the Gram matrix (Cartan inverse)
    # For SU(3), the symmetrized Cartan matrix is
    #   <alpha_i, alpha_j>: alpha_1.alpha_1 = 2, alpha_1.alpha_2 = -1, alpha_2.alpha_2 = 2
    # In omega basis, <omega_i, omega_j> = (Cartan^-1)_{ij}
    # For SU(3): Cartan = [[2,-1],[-1,2]], Cartan^-1 = (1/3)[[2,1],[1,2]]

    G = sp.Matrix([[Rational(2,3), Rational(1,3)],
                   [Rational(1,3), Rational(2,3)]])

    def inner(w1, w2):
        a1, b1 = w1
        a2, b2 = w2
        v1 = sp.Matrix([a1, b1])
        v2 = sp.Matrix([a2, b2])
        return (v1.T @ G @ v2)[0, 0]

    # Positive roots in omega basis
    alpha_1 = (2, -1)
    alpha_2 = (-1, 2)
    alpha_3 = (1, 1)  # alpha_1 + alpha_2
    pos_roots = [alpha_1, alpha_2, alpha_3]

    rho = (1, 1)  # half sum of pos roots in omega basis: (a1+a2+a3)/2 = (1,1)

    lam = (p, q)

    # First, identify the set of weights via Weyl orbit hull
    def weight_le(mu, nu):
        """Check mu <= nu in the standard order (nu - mu in N*alpha_1 + N*alpha_2)."""
        d = (nu[0] - mu[0], nu[1] - mu[1])
        # Express d in alpha basis: alpha_1 = (2,-1), alpha_2 = (-1,2)
        # d = n1 * (2,-1) + n2 * (-1, 2)
        # 2 n1 - n2 = d_a
        # -n1 + 2 n2 = d_b
        # det = 4 - 1 = 3
        # n1 = (2*d_a + d_b)/3
        # n2 = (d_a + 2*d_b)/3
        n1 = Rational(2 * d[0] + d[1], 3)
        n2 = Rational(d[0] + 2 * d[1], 3)
        return n1.is_integer and n2.is_integer and n1 >= 0 and n2 >= 0

    # All candidate weights mu with mu <= lambda in dominance partial order
    # Bounds: weights live in {hw - n1*alpha_1 - n2*alpha_2 : n1, n2 >= 0}
    # n1, n2 are bounded by some function of p,q.
    # alpha_1 = (2,-1), alpha_2 = (-1,2). So a = p - 2*n1 + n2, b = q + n1 - 2*n2.
    # For mu in irrep, n1 and n2 are bounded; rough bound: n1, n2 <= p + q + 1 each
    # (very loose but safe).
    candidates = []
    n_bound = 2 * (p + q) + 2  # generous bound on lowering depth
    for n1 in range(n_bound + 1):
        for n2 in range(n_bound + 1):
            a_test = p - 2 * n1 + n2
            b_test = q + n1 - 2 * n2
            mu = (a_test, b_test)
            if mu not in candidates:
                candidates.append(mu)

    # Filter: also must lie in convex hull of Weyl orbit of lam.
    # In SU(3), this is automatic for mu <= lam in the orbit (Bourbaki / standard).
    # Actually NOT automatic: mu must be in convex hull AND in lam - Q+ AND mu must
    # have mult > 0. Freudenthal handles all this.

    # Compute multiplicities by Freudenthal
    mults: Dict[Tuple[int, int], Rational] = {}

    # Sort candidates by depth from lambda. Depth = n1 + n2 where lam - mu = n1*a1 + n2*a2.
    # n1 = (2*(p-a) + (q-b))/3, n2 = ((p-a) + 2*(q-b))/3. Depth = n1 + n2 = (p-a + q-b).
    # Process by ascending depth (highest weight first).
    def depth_from_lam(m):
        d = (lam[0] - m[0], lam[1] - m[1])
        n1 = Rational(2 * d[0] + d[1], 3)
        n2 = Rational(d[0] + 2 * d[1], 3)
        return int(n1 + n2)

    # Filter candidates that are reachable (n1, n2 non-negative integers)
    valid_candidates = []
    for m in candidates:
        d = (lam[0] - m[0], lam[1] - m[1])
        n1 = Rational(2 * d[0] + d[1], 3)
        n2 = Rational(d[0] + 2 * d[1], 3)
        if n1.is_integer and n2.is_integer and n1 >= 0 and n2 >= 0:
            valid_candidates.append(m)

    candidates_sorted = sorted(valid_candidates, key=depth_from_lam)

    mults[lam] = Rational(1)

    for mu in candidates_sorted:
        if mu == lam:
            continue
        # Freudenthal recursion:
        # ((<lam+rho, lam+rho> - <mu+rho, mu+rho>)) * mult(mu) =
        #   2 * sum_{alpha > 0} sum_{k>=1} <mu + k*alpha, alpha> * mult(mu + k*alpha)

        lhs_factor = inner((lam[0] + rho[0], lam[1] + rho[1]),
                           (lam[0] + rho[0], lam[1] + rho[1])) - \
                     inner((mu[0] + rho[0], mu[1] + rho[1]),
                           (mu[0] + rho[0], mu[1] + rho[1]))

        if lhs_factor == 0:
            continue  # Should only happen if mu is in W-orbit of lam, mult = 1

        rhs = Rational(0)
        for alpha in pos_roots:
            k = 1
            while True:
                mu_plus = (mu[0] + k * alpha[0], mu[1] + k * alpha[1])
                if mu_plus not in mults:
                    break
                rhs += inner(mu_plus, alpha) * mults[mu_plus]
                k += 1
        rhs *= 2

        mult_mu = rhs / lhs_factor
        if mult_mu > 0:
            mults[mu] = mult_mu

    # Convert to list with multiplicities
    result = []
    for mu, m in mults.items():
        if m > 0:
            for _ in range(int(m)):
                result.append(mu)

    return result


def tensor_product_su3(lam: Tuple[int, int], mu: Tuple[int, int]) -> Dict[Tuple[int, int], int]:
    """
    Decompose V_lambda ⊗ V_mu into irreducibles, return {nu: multiplicity}.

    Use Brauer-Klimyk (Racah-Speiser): for each weight mu_i of V_mu (with mult m_i),
    consider lambda + mu_i + rho. If this lies in the closure of the dominant Weyl
    chamber (boundary or interior), count contributions:
      - If on the boundary (some <lambda+mu_i+rho, alpha_simple> = 0): zero contribution
      - If interior: reflect into dominant chamber via shortest Weyl element w,
        contributes sign(w) * m_i to the irrep with HW (w(lambda+mu_i+rho) - rho)
    """
    p, q = lam
    weights_mu = weights_of_irrep_su3(*mu)

    rho = (1, 1)

    result: Dict[Tuple[int, int], int] = {}

    # Group weights of V_mu by value (dict with multiplicity)
    weight_mults: Dict[Tuple[int, int], int] = {}
    for w in weights_mu:
        weight_mults[w] = weight_mults.get(w, 0) + 1

    for w, m in weight_mults.items():
        # nu_shifted = lambda + w + rho
        a = lam[0] + w[0] + rho[0]
        b = lam[1] + w[1] + rho[1]

        # Check if (a, b) is in interior of dominant chamber: a > 0 and b > 0
        # If on boundary (a = 0 or b = 0): contributes 0
        # If a < 0 or b < 0: reflect via Weyl group

        # SU(3) Weyl group: S_3 with simple reflections
        #   s_1: (a, b) -> (-a, a+b)  [reflection across alpha_1 hyperplane]
        #   s_2: (a, b) -> (a+b, -b)  [reflection across alpha_2 hyperplane]
        # We need to find w in W and sign(w) such that w*(a,b) is in dominant chamber
        # OR detect that (a,b) is on a wall.

        # Use iterated reflections: keep reflecting until in dominant chamber or on wall.
        sign = 1
        on_wall = False
        cur_a, cur_b = a, b
        max_iter = 20
        for _ in range(max_iter):
            if cur_a == 0 or cur_b == 0:
                on_wall = True
                break
            if cur_a > 0 and cur_b > 0:
                break
            if cur_a < 0:
                cur_a, cur_b = -cur_a, cur_a + cur_b
                sign *= -1
                continue
            if cur_b < 0:
                cur_a, cur_b = cur_a + cur_b, -cur_b
                sign *= -1
                continue

        if on_wall:
            continue

        # Now (cur_a, cur_b) is in interior; HW is (cur_a - 1, cur_b - 1)
        nu = (cur_a - 1, cur_b - 1)
        if nu[0] < 0 or nu[1] < 0:
            continue
        result[nu] = result.get(nu, 0) + sign * m

    # Filter zero entries
    return {k: v for k, v in result.items() if v != 0}


# =====================================================================
# Verify Casimir triangle inequality
# =====================================================================

def verify_inequality(lam: Tuple[int, int], lam_prime: Tuple[int, int],
                      verbose: bool = True) -> Tuple[bool, List]:
    """
    Verify |C(lam) - C(lam')| <= C(sigma) for every sigma in lam ⊗ (lam')^*.

    Dual of (p, q) is (q, p) for SU(3).
    """
    p, q = lam
    p_prime, q_prime = lam_prime
    lam_prime_dual = (q_prime, p_prime)

    C_lam = casimir_su3(*lam)
    C_lam_prime = casimir_su3(*lam_prime)
    diff = abs(C_lam - C_lam_prime)

    decomp = tensor_product_su3(lam, lam_prime_dual)

    all_pass = True
    violations = []

    if verbose:
        print(f"  lambda = {lam}, C(lambda) = {C_lam}")
        print(f"  lambda' = {lam_prime}, C(lambda') = {C_lam_prime}")
        print(f"  (lambda')* = {lam_prime_dual}")
        print(f"  |C(lambda) - C(lambda')| = {diff}")
        print(f"  Decomposition lambda (x) (lambda')* = {decomp}")

    for sigma, mult in decomp.items():
        if mult > 0:  # only check actual representations (skip if any were anti-symmetric)
            C_sigma = casimir_su3(*sigma)
            holds = diff <= C_sigma
            status = "PASS" if holds else "FAIL"
            if verbose:
                print(f"    sigma = {sigma}: C(sigma) = {C_sigma}, |dC| = {diff}, {status}")
            if not holds:
                all_pass = False
                violations.append((sigma, C_sigma, diff))

    return all_pass, violations


# =====================================================================
# Run tests
# =====================================================================

def run_tests():
    print("=" * 70)
    print("SU(3) Casimir triangle inequality verification")
    print("=" * 70)

    # Sanity checks first
    print("\n--- Sanity: known Casimir values ---")
    print(f"C(0,0) = {casimir_su3(0, 0)} (expect 0, singlet)")
    print(f"C(1,0) = {casimir_su3(1, 0)} (expect 4/3, fundamental)")
    print(f"C(0,1) = {casimir_su3(0, 1)} (expect 4/3, antifund)")
    print(f"C(1,1) = {casimir_su3(1, 1)} (expect 3, adjoint)")
    print(f"C(2,0) = {casimir_su3(2, 0)} (expect 10/3, sextet)")
    print(f"C(3,0) = {casimir_su3(3, 0)} (expect 6, decuplet)")
    print(f"C(2,1) = {casimir_su3(2, 1)} (expect 16/3)")
    print(f"C(2,2) = {casimir_su3(2, 2)} (expect 8, [27])")

    # Sanity: known tensor product 3 (x) 3-bar = 1 (+) 8
    print("\n--- Sanity: 3 (x) 3-bar decomposition ---")
    decomp = tensor_product_su3((1, 0), (0, 1))
    print(f"  (1,0) (x) (0,1) = {decomp}")
    print(f"  Expected: {{(0,0): 1, (1,1): 1}}")
    expected = {(0, 0): 1, (1, 1): 1}
    if decomp == expected:
        print("  PASS correct")
    else:
        print("  FAIL MISMATCH")

    # Sanity: 3 (x) 3 = 6 (+) 3-bar
    print("\n--- Sanity: 3 (x) 3 decomposition ---")
    decomp = tensor_product_su3((1, 0), (1, 0))
    print(f"  (1,0) (x) (1,0) = {decomp}")
    print(f"  Expected: {{(2,0): 1, (0,1): 1}}")
    expected = {(2, 0): 1, (0, 1): 1}
    if decomp == expected:
        print("  PASS correct")
    else:
        print("  FAIL MISMATCH")

    # Sanity: 8 (x) 8 = 1 (+) 8 (+) 8 (+) 10 (+) 10-bar (+) 27
    print("\n--- Sanity: 8 (x) 8 decomposition ---")
    decomp = tensor_product_su3((1, 1), (1, 1))
    print(f"  (1,1) (x) (1,1) = {decomp}")
    print(f"  Expected: 1 + 8 + 8 + 10 + 10-bar + 27")
    expected = {(0, 0): 1, (1, 1): 2, (3, 0): 1, (0, 3): 1, (2, 2): 1}
    if decomp == expected:
        print("  PASS correct")
    else:
        print(f"  FAIL MISMATCH")
        print(f"     got: {decomp}")
        print(f"     expected: {expected}")

    # Now systematic verification
    print("\n--- Systematic: all (p,q), (p',q') with p+q <= 3, p'+q' <= 3 ---")

    test_cases = []
    for p in range(4):
        for q in range(4):
            if p + q > 3:
                continue
            for p2 in range(4):
                for q2 in range(4):
                    if p2 + q2 > 3:
                        continue
                    test_cases.append(((p, q), (p2, q2)))

    fail_count = 0
    pass_count = 0
    all_violations = []

    for lam, lam_prime in test_cases:
        passes, violations = verify_inequality(lam, lam_prime, verbose=False)
        if not passes:
            fail_count += 1
            all_violations.extend([(lam, lam_prime, v) for v in violations])
            print(f"  FAIL: lambda={lam}, lambda'={lam_prime}")
            for sigma, C_sigma, diff in violations:
                print(f"    sigma={sigma}, C(sigma)={C_sigma}, |dC|={diff}, ratio={float(diff)/float(C_sigma) if C_sigma > 0 else 'inf'}")
        else:
            pass_count += 1

    print(f"\n--- Summary ---")
    print(f"  Passed: {pass_count} / {len(test_cases)}")
    print(f"  Failed: {fail_count} / {len(test_cases)}")

    if fail_count > 0:
        print("\n  ALL VIOLATIONS:")
        for (lam, lam_prime, (sigma, C_sigma, diff)) in all_violations:
            ratio = float(diff) / float(C_sigma) if C_sigma > 0 else float('inf')
            print(f"    lambda={lam}, lambda'={lam_prime}, sigma={sigma}, C(sigma)={C_sigma}, |dC|={diff}, ratio={ratio:.4f}")

    # Detailed view of the most interesting cases
    print("\n--- Detailed view: 8 (x) 8-bar = (1,1) (x) (1,1) ---")
    verify_inequality((1, 1), (1, 1), verbose=True)

    print("\n--- Detailed view: extreme case (3,0) and (0,3) (decuplets) ---")
    verify_inequality((3, 0), (0, 3), verbose=True)

    print("\n--- Detailed view: (2,1) and (1,2) (mismatched) ---")
    verify_inequality((2, 1), (1, 2), verbose=True)

    print("\n--- Detailed view: (3,0) and (1,0) ---")
    verify_inequality((3, 0), (1, 0), verbose=True)


if __name__ == "__main__":
    run_tests()
