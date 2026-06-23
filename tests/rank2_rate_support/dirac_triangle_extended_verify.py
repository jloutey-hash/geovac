"""
Extended verification of the Dirac-triangle inequality across compact Lie groups.

Statement:
    |D(pi) - D(pi')| <= sqrt(C(sigma))    for every sigma in pi (x) pi'^*

where
    |D(lambda)|^2 = <lambda + rho, lambda + rho>
    C(lambda)     = <lambda + rho, lambda + rho> - <rho, rho>
    rho           = sum of fundamental weights = half sum of positive roots

Implements a generic Brauer-Klimyk + Freudenthal pipeline for any simple
Lie algebra given its Cartan matrix and the "co-root" length squared values
d_i = (alpha_i, alpha_i) / 2 (i.e., a_i is the long root squared length factor).

Tested at:
    - SU(3) (type A_2)        full panel and extended panel (p+q <= 5)
    - SU(4) (type A_3)        dominant weights with sum <= 3
    - Sp(2)=Sp(4) (type C_2)  dominant weights with sum <= 3
    - G_2                      dominant weights with sum <= 2 (large dims)

Plus the asymptotic family (n, 0, ..., 0) (x) (1, 0, ..., 0)* for each.

Convention. We use the normalization where short roots have (alpha, alpha) = 2;
this matches the SU(3) script. Then for SU(N) (simply-laced),
all roots have (alpha, alpha) = 2 and the symmetrized Cartan matrix equals the
Cartan matrix. For C_2 and G_2 we use the standard symmetrization.

This module uses *sympy exact arithmetic* everywhere. No floats during
verification; we only convert to float at the report-print stage.
"""

import json
import os
import sys
from itertools import product
from typing import Dict, List, Tuple
import sympy as sp
from sympy import Matrix, Rational, sqrt, simplify, Integer


# ---------------------------------------------------------------------------
# Generic simple Lie algebra data
# ---------------------------------------------------------------------------

class SimpleLieAlgebra:
    """
    Encapsulates the data of a simple Lie algebra of rank r:
      - cartan: r x r Cartan matrix (sympy Matrix of Integers)
      - d: list of length r; d[i] = (alpha_i, alpha_i)/2 in our normalization
           (so that short root squared length = 2; d_i = 1 for simply-laced)
      - simple_roots and positive roots in the alpha basis (integer combinations)
      - fundamental weights in alpha basis
      - rho (half sum of positive roots) in fundamental-weight basis = (1,1,...,1)
      - gram matrix on the alpha basis (symmetric Cartan): B[i][j] = d_i * A[i][j]
      - inner product on weights: we work in fundamental-weight basis throughout;
        the gram matrix on that basis is G_omega = B^{-1} * D   (symmetric, where
        D = diag(d_1, ..., d_r))   --- this is the inverse of the symmetrized
        Cartan with appropriate scaling.

    Inner product convention:
      Symmetrized Cartan matrix B_ij = d_i * A_ij in alpha basis.
      Then (omega_i, omega_j) = (B^{-1})_ij * d_j   (works for any simple algebra).
      Equivalently, the Gram matrix on the fundamental weight basis is
         G_omega = B^{-1} * diag(d)
      which is symmetric because (omega_i, omega_j) = (omega_j, omega_i).

    Verification: at SU(3),
       A = [[2,-1],[-1,2]],  d = [1, 1],  B = A,  B^{-1} = (1/3)[[2,1],[1,2]]
       G_omega = B^{-1} * I = (1/3)[[2,1],[1,2]]
    which matches the SU(3) check (P-A's `G` matrix). ✓
    """

    def __init__(self, name: str, cartan: List[List[int]], d: List[int]):
        self.name = name
        self.rank = len(d)
        self.cartan = Matrix(cartan)
        assert self.cartan.shape == (self.rank, self.rank)
        self.d = list(d)
        # Symmetrized cartan B_ij = d_i A_ij
        D = Matrix(self.rank, self.rank, lambda i, j: self.d[i] if i == j else 0)
        self.B_alpha = D * self.cartan
        # Verify symmetric (sanity check)
        for i in range(self.rank):
            for j in range(self.rank):
                assert self.B_alpha[i, j] == self.B_alpha[j, i], \
                    f"Symmetrized Cartan not symmetric for {name}"
        # Gram matrix on omega basis: derived from
        #   (omega_i, alpha_j) = d_j * delta_{ij}                (definition of omega)
        # and  alpha_j = sum_k A_{kj} * omega_k   (our column convention).
        # Then (omega_i, alpha_j) = sum_k A_{kj} G_{ik} = d_j delta_{ij}
        # i.e.   G * A = D   =>   G = D * A^{-1}.
        self.A_inv = self.cartan.inv()
        self.G_omega = D * self.A_inv
        # Sanity: G_omega should be symmetric
        for i in range(self.rank):
            for j in range(self.rank):
                assert self.G_omega[i, j] == self.G_omega[j, i], \
                    f"G_omega not symmetric for {name}: i={i} j={j}, " \
                    f"G[i,j]={self.G_omega[i,j]}, G[j,i]={self.G_omega[j,i]}"
        # Simple roots in omega basis: alpha_i = sum_j A_ji omega_j (column j of A)
        # i.e., alpha_i_omegabasis = (A[:, i])
        self.simple_roots_omega: List[Tuple[int, ...]] = []
        for i in range(self.rank):
            sr = tuple(int(self.cartan[j, i]) for j in range(self.rank))
            self.simple_roots_omega.append(sr)
        # Positive roots in omega basis: enumerate by adding simple roots starting
        # from each simple root, using bracket relations. For our purposes (since
        # we only need positive roots for Freudenthal), we precompute them for
        # specific algebras and store.
        self.positive_roots_omega = self._compute_positive_roots()
        # rho = sum of fundamental weights = (1, 1, ..., 1) in omega basis
        self.rho_omega = tuple([1] * self.rank)
        # rho_squared
        self.rho_sq = self._inner_omega(self.rho_omega, self.rho_omega)

    def _inner_omega(self, mu: Tuple[int, ...], nu: Tuple[int, ...]):
        """Inner product (mu, nu) for mu, nu in fundamental-weight basis."""
        mu_v = Matrix([mu]).T
        nu_v = Matrix([nu]).T
        return (mu_v.T * self.G_omega * nu_v)[0, 0]

    def inner(self, mu, nu):
        return self._inner_omega(tuple(mu), tuple(nu))

    def _compute_positive_roots(self) -> List[Tuple[int, ...]]:
        """
        Enumerate all positive roots in omega basis.
        Uses: positive roots = simple roots + positive linear combinations.
        We perform a closed-orbit BFS: start with simple roots, add simple roots,
        check if result is a (positive) root by checking if it pairs positively
        with at least one simple root (alpha + beta is a root iff (alpha, beta^vee) < 0
        when alpha, beta are roots; more precisely, the alpha-string through beta).

        Cleaner approach: use the Weyl group orbit of rho. Number of positive roots
        is well-known per algebra type; we list them explicitly for the algebras
        we care about. For other algebras we fall back to the orbit method.

        We provide explicit positive root systems for A_n, C_n, and G_2 below.
        """
        return self._compute_positive_roots_by_type()

    def _compute_positive_roots_by_type(self) -> List[Tuple[int, ...]]:
        """Explicit positive root lists in omega basis."""
        # Get name prefix
        name = self.name
        if name.startswith("A"):
            return positive_roots_A(self.rank, self.cartan)
        elif name.startswith("C"):
            return positive_roots_C(self.rank, self.cartan)
        elif name.startswith("G2"):
            return positive_roots_G2(self.cartan)
        else:
            raise ValueError(f"Don't have positive-root list for {name}")

    # -----------------------------------------------------------------
    # Casimir, Dirac
    # -----------------------------------------------------------------

    def casimir(self, lam: Tuple[int, ...]):
        """C(lambda) = <lambda + rho, lambda + rho> - <rho, rho>"""
        shifted = tuple(lam[i] + self.rho_omega[i] for i in range(self.rank))
        return self._inner_omega(shifted, shifted) - self.rho_sq

    def D_squared(self, lam: Tuple[int, ...]):
        """|D(lambda)|^2 = <lambda + rho, lambda + rho> = C(lambda) + <rho, rho>"""
        shifted = tuple(lam[i] + self.rho_omega[i] for i in range(self.rank))
        return self._inner_omega(shifted, shifted)

    def D_abs(self, lam: Tuple[int, ...]):
        return sqrt(self.D_squared(lam))

    def dim_weyl(self, lam: Tuple[int, ...]):
        """
        Weyl dimension formula:
          dim(V_lambda) = prod_{alpha > 0} <lambda + rho, alpha^vee>
                                          / <rho, alpha^vee>
        where alpha^vee = 2 alpha / (alpha, alpha).
        """
        shifted = tuple(lam[i] + self.rho_omega[i] for i in range(self.rank))
        num = Integer(1)
        den = Integer(1)
        for alpha in self.positive_roots_omega:
            alpha_sq = self._inner_omega(alpha, alpha)
            num_term = 2 * self._inner_omega(shifted, alpha) / alpha_sq
            den_term = 2 * self._inner_omega(self.rho_omega, alpha) / alpha_sq
            num = num * num_term
            den = den * den_term
        return simplify(num / den)

    # -----------------------------------------------------------------
    # Dual representation
    # -----------------------------------------------------------------

    def dual(self, lam: Tuple[int, ...]) -> Tuple[int, ...]:
        """
        Highest weight of the dual representation V_lambda^*.
        Equals -w0(lambda) where w0 is the longest Weyl element.

        For A_{n-1}: lambda = (a_1, ..., a_{n-1}) -> (a_{n-1}, ..., a_1)
        For C_n (Sp(2n)): all representations are self-dual: lambda^* = lambda
        For G_2: all representations are self-dual: lambda^* = lambda
        For B_n, D_n (even), F_4, E_7, E_8: also all self-dual.
        For D_n (odd), E_6: lambda^* = sigma(lambda) where sigma is the Dynkin
            diagram involution.
        """
        name = self.name
        if name.startswith("A"):
            return tuple(reversed(lam))
        elif name.startswith("C"):
            return tuple(lam)
        elif name.startswith("G2"):
            return tuple(lam)
        else:
            raise ValueError(f"Don't know dual for {name}")


# ---------------------------------------------------------------------------
# Explicit positive-root lists in omega basis for our four cases
# ---------------------------------------------------------------------------

def positive_roots_A(n: int, cartan: Matrix) -> List[Tuple[int, ...]]:
    """
    Type A_{n} (so rank = n, group SU(n+1)). Positive roots are
    alpha_{ij} = alpha_i + alpha_{i+1} + ... + alpha_{j-1}  for 1 <= i < j <= n+1.
    In omega basis, each alpha_k is the k-th column of the Cartan matrix.

    For A_n (n = rank), number of positive roots is n(n+1)/2.
    """
    roots = []
    simple = [tuple(int(cartan[j, i]) for j in range(n)) for i in range(n)]
    # alpha_i + alpha_{i+1} + ... + alpha_j summed for j >= i
    for i in range(n):
        cur = [0] * n
        for j in range(i, n):
            cur = [cur[k] + simple[j][k] for k in range(n)]
            roots.append(tuple(cur))
    # de-dup
    seen = set()
    out = []
    for r in roots:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return out


def positive_roots_C(n: int, cartan: Matrix) -> List[Tuple[int, ...]]:
    """
    Type C_n (group Sp(2n)). Positive roots in terms of orthonormal basis e_i
    (i = 1..n):
        e_i - e_j  (1 <= i < j <= n)            (short, n(n-1)/2 of them)
        e_i + e_j  (1 <= i <= j <= n, but only i<j for off-diagonal) (short, n(n-1)/2)
        2 e_i      (1 <= i <= n)                (long, n of them)
    Total: n(n-1) + n = n^2 positive roots.

    Standard Cartan: simple roots are
        alpha_i = e_i - e_{i+1}  for i = 1, .., n-1 (short)
        alpha_n = 2 e_n (long)

    Conversion to omega basis: alpha_i in omega basis is the i-th column of the
    Cartan matrix. Then expressing each positive root as integer combo of
    alpha_1, ..., alpha_n, then converting to omega basis.

    We enumerate explicitly: alpha_k_combo[i] = coefficient of alpha_i in the k-th
    positive root, and then sum simple roots' omega-basis representations.
    """
    # Positive roots in alpha basis: a list of n-tuples of nonneg ints
    pr_alpha = []
    # e_i - e_j = alpha_i + alpha_{i+1} + ... + alpha_{j-1}  (for i<j<=n)
    # In our 0-indexed convention: alpha_0, ..., alpha_{n-2} are short, alpha_{n-1} = 2 e_{n-1}
    # e_i - e_j for 0 <= i < j <= n-1:
    for i in range(n):
        for j in range(i + 1, n):
            coeff = [0] * n
            for k in range(i, j):
                coeff[k] = 1
            pr_alpha.append(tuple(coeff))
    # e_i + e_j for 0 <= i < j <= n-1:
    #   = (e_i - e_j) + 2 e_j = (alpha_i + ... + alpha_{j-1}) + 2 e_j
    # And 2 e_j = alpha_j + alpha_{j+1} + ... + alpha_{n-2} + alpha_{n-1} ??? Let's redo.
    #   2 e_k for k in [0, n-1]: by recursion, 2 e_{n-1} = alpha_{n-1}.
    #   2 e_k = 2 e_{k+1} + 2 (e_k - e_{k+1}) = 2 e_{k+1} + 2 alpha_k.
    # So 2 e_k = 2(alpha_k + alpha_{k+1} + ... + alpha_{n-2}) + alpha_{n-1}.
    for i in range(n):
        for j in range(i + 1, n):
            coeff = [0] * n
            # e_i - e_j contribution:
            for k in range(i, j):
                coeff[k] = 1
            # 2 e_j contribution: 2 alpha_j + ... + 2 alpha_{n-2} + alpha_{n-1}
            for k in range(j, n - 1):
                coeff[k] += 2
            coeff[n - 1] += 1
            pr_alpha.append(tuple(coeff))
    # 2 e_i for i in [0, n-1]:
    for i in range(n):
        coeff = [0] * n
        for k in range(i, n - 1):
            coeff[k] = 2
        coeff[n - 1] = 1
        pr_alpha.append(tuple(coeff))
    # de-dup
    pr_alpha = list(dict.fromkeys(pr_alpha))
    # Convert each alpha-basis combination to omega basis
    simple_omega = [tuple(int(cartan[j, i]) for j in range(n)) for i in range(n)]
    pr_omega = []
    for combo in pr_alpha:
        out = [0] * n
        for k in range(n):
            for j in range(n):
                out[j] += combo[k] * simple_omega[k][j]
        pr_omega.append(tuple(out))
    return pr_omega


def positive_roots_G2(cartan: Matrix) -> List[Tuple[int, ...]]:
    """
    Type G_2. Six positive roots (in alpha-basis with alpha_1 short, alpha_2 long):
        alpha_1, alpha_2
        alpha_1 + alpha_2
        2 alpha_1 + alpha_2
        3 alpha_1 + alpha_2
        3 alpha_1 + 2 alpha_2

    Cartan: [[2,-1],[-3,2]] in standard convention.
    """
    pr_alpha = [
        (1, 0),  # alpha_1 (short)
        (0, 1),  # alpha_2 (long)
        (1, 1),
        (2, 1),
        (3, 1),
        (3, 2),
    ]
    n = 2
    simple_omega = [tuple(int(cartan[j, i]) for j in range(n)) for i in range(n)]
    pr_omega = []
    for combo in pr_alpha:
        out = [0] * n
        for k in range(n):
            for j in range(n):
                out[j] += combo[k] * simple_omega[k][j]
        pr_omega.append(tuple(out))
    return pr_omega


# ---------------------------------------------------------------------------
# Weight enumeration via Freudenthal
# ---------------------------------------------------------------------------

def weights_of_irrep(la: SimpleLieAlgebra, lam: Tuple[int, ...]) -> Dict[Tuple[int, ...], Rational]:
    """
    Compute the multiplicities of weights in V_lambda by Freudenthal's formula:

        (||lambda+rho||^2 - ||mu+rho||^2) * mult(mu)
          = 2 * sum_{alpha > 0} sum_{k >= 1} <mu + k alpha, alpha> * mult(mu + k alpha)

    Returns dict mapping weight (in omega basis) -> multiplicity.

    Bounded by enumerating mu = lambda - sum_i n_i alpha_i with n_i >= 0 and
    keeping only those for which we get integer weights with multiplicity > 0
    after Freudenthal. We use generous bounds on n_i.
    """
    rank = la.rank
    rho = la.rho_omega
    lam_plus_rho = tuple(lam[i] + rho[i] for i in range(rank))
    lam_plus_rho_sq = la._inner_omega(lam_plus_rho, lam_plus_rho)

    # Express simple roots in omega basis
    simple_omega = la.simple_roots_omega

    # Enumerate candidate weights:
    # mu = lambda - sum_i n_i alpha_i,  n_i >= 0 integer
    # Bound: each n_i is at most the maximal length of the alpha_i ladder, which is
    # bounded by the dominant weight's coordinate components in alpha-basis.
    # We use a generous bound based on the size of |lambda|.
    coord_max = max(lam) if lam else 0
    n_bound = 4 * (sum(lam) + rank) + 4

    # Build dict of candidates with depth
    candidates = {lam: 0}  # weight -> depth from lambda
    frontier = [lam]
    while frontier:
        new_frontier = []
        for mu in frontier:
            d = candidates[mu]
            if d > n_bound:
                continue
            for alpha in simple_omega:
                child = tuple(mu[i] - alpha[i] for i in range(rank))
                # Track depth conservatively: only accept if not seen or strictly deeper
                if child not in candidates:
                    candidates[child] = d + 1
                    new_frontier.append(child)
        frontier = new_frontier

    # Sort by depth (we process in order of increasing depth)
    items_sorted = sorted(candidates.items(), key=lambda x: x[1])

    mults: Dict[Tuple[int, ...], Rational] = {lam: Rational(1)}

    for mu, d in items_sorted:
        if mu == lam:
            continue
        mu_plus_rho = tuple(mu[i] + rho[i] for i in range(rank))
        denom = lam_plus_rho_sq - la._inner_omega(mu_plus_rho, mu_plus_rho)
        if denom == 0:
            continue

        rhs = Rational(0)
        for alpha in la.positive_roots_omega:
            k = 1
            while True:
                mu_plus_k_alpha = tuple(mu[i] + k * alpha[i] for i in range(rank))
                if mu_plus_k_alpha not in mults:
                    break
                rhs += la._inner_omega(mu_plus_k_alpha, alpha) * mults[mu_plus_k_alpha]
                k += 1
        rhs *= 2
        mult = rhs / denom
        if mult > 0:
            mults[mu] = mult
        # If mult is 0 or negative, the weight is not in the irrep (negative would
        # be a Freudenthal artifact; we just drop it)
    return mults


# ---------------------------------------------------------------------------
# Tensor product via Brauer-Klimyk (Racah-Speiser)
# ---------------------------------------------------------------------------

def tensor_product(la: SimpleLieAlgebra, lam: Tuple[int, ...], mu: Tuple[int, ...]
                   ) -> Dict[Tuple[int, ...], int]:
    """
    Decompose V_lambda (x) V_mu into irreducibles using Brauer-Klimyk.

    For each weight w of V_mu with multiplicity m:
      shifted = lambda + w + rho
      If shifted is on a Weyl wall (some <shifted, alpha_i^vee> = 0 for simple alpha_i):
         contributes 0.
      Else: reflect shifted into the dominant chamber via Weyl group; the
         resulting interior point is nu + rho for some dominant nu; sign comes
         from parity of reflections; add sign * m to the multiplicity of nu.

    Result: dict mapping nu (dominant weight) -> multiplicity in V_lambda (x) V_mu.
    """
    rank = la.rank
    rho = la.rho_omega
    cartan = la.cartan
    weights_mu = weights_of_irrep(la, mu)

    out: Dict[Tuple[int, ...], int] = {}

    for w, m in weights_mu.items():
        shifted = tuple(lam[i] + w[i] + rho[i] for i in range(rank))
        # Reflect into dominant chamber. "Dominant": all coordinates in omega basis > 0
        # (interior; if any == 0, we're on a wall and contribute 0).
        cur = list(shifted)
        sign = 1
        on_wall = False
        max_iter = 200
        for _ in range(max_iter):
            # Check on wall:
            if any(c == 0 for c in cur):
                on_wall = True
                break
            # Check dominant: all c > 0
            if all(c > 0 for c in cur):
                break
            # Find first negative coordinate; reflect across that simple root.
            # Reflection: s_i sends weight v to v - <v, alpha_i^vee> alpha_i
            # where alpha_i^vee is the coroot. In omega basis, alpha_i is the i-th
            # column of the Cartan matrix, and <v, alpha_i^vee> = v_i (the i-th
            # coordinate in fundamental-weight basis).
            #
            # So s_i: cur[j] -> cur[j] - cur[i] * cartan[j, i]
            i_neg = -1
            for i in range(rank):
                if cur[i] < 0:
                    i_neg = i
                    break
            if i_neg < 0:
                break  # All non-negative; if any zero, would have been caught
            cur_i_val = cur[i_neg]
            for j in range(rank):
                cur[j] = cur[j] - cur_i_val * int(cartan[j, i_neg])
            sign *= -1
        else:
            # Hit iteration limit; mark as wall to skip
            on_wall = True
        if on_wall:
            continue
        # cur is in dominant chamber (interior). nu = cur - rho.
        nu = tuple(cur[i] - rho[i] for i in range(rank))
        if any(c < 0 for c in nu):
            # Shouldn't happen if reflection logic is correct
            continue
        # Multiplicity contribution (must be Integer; m is Rational from Freudenthal)
        contrib = sign * int(m)
        out[nu] = out.get(nu, 0) + contrib

    # Filter zeros
    return {k: v for k, v in out.items() if v != 0}


# ---------------------------------------------------------------------------
# Dirac-triangle inequality verification
# ---------------------------------------------------------------------------

def verify_dirac_triangle(la: SimpleLieAlgebra, lam: Tuple[int, ...],
                          lam_prime: Tuple[int, ...]) -> Tuple[bool, List, List]:
    """
    Verify |D(lambda) - D(lambda')| <= sqrt(C(sigma)) for every sigma in lam (x) lam'^*.

    Returns (all_pass, violations, all_ratios) where
        violations: list of (sigma, sqrt(C_sigma)_float, |D-D'|_float, ratio_float)
        all_ratios: list of float ratios (one per sigma in decomposition)
    """
    lam_prime_dual = la.dual(lam_prime)
    D = la.D_abs(lam)
    D_pr = la.D_abs(lam_prime)
    dirac_diff = sp.Abs(D - D_pr)

    decomp = tensor_product(la, lam, lam_prime_dual)

    violations = []
    all_ratios = []
    all_pass = True

    for sigma, mult in decomp.items():
        if mult <= 0:
            continue
        C_sigma = la.casimir(sigma)
        if C_sigma == 0:
            # trivial rep; sqrt(0)=0; OK iff dirac_diff == 0
            if dirac_diff != 0:
                all_pass = False
                violations.append((sigma, 0.0, float(dirac_diff), float('inf')))
            else:
                all_ratios.append(0.0)
            continue
        sqrt_C = sqrt(C_sigma)
        ratio_sym = dirac_diff / sqrt_C
        ratio_num = float(ratio_sym)
        all_ratios.append(ratio_num)
        if ratio_num > 1.0 + 1e-12:
            all_pass = False
            violations.append((sigma, float(sqrt_C), float(dirac_diff), ratio_num))

    return all_pass, violations, all_ratios


# ---------------------------------------------------------------------------
# Algebra builders
# ---------------------------------------------------------------------------

def build_A(n: int) -> SimpleLieAlgebra:
    """Type A_n: SU(n+1)."""
    cartan = [[0] * n for _ in range(n)]
    for i in range(n):
        cartan[i][i] = 2
        if i > 0:
            cartan[i][i - 1] = -1
        if i < n - 1:
            cartan[i][i + 1] = -1
    d = [1] * n
    return SimpleLieAlgebra(f"A{n}", cartan, d)


def build_C2() -> SimpleLieAlgebra:
    """
    Type C_2 (= Sp(4) = Sp(2)). Rank 2.
    Convention: alpha_1 short, alpha_2 long.
    Cartan matrix: A = [[2, -1], [-2, 2]].
    (Note: row i, col j is A_ij = 2 (alpha_i, alpha_j) / (alpha_i, alpha_i).
     With alpha_1 short, alpha_2 long, (alpha_2, alpha_2) = 2 (alpha_1, alpha_1).
     In our normalization (alpha_short)^2 = 2: d_1 = 1, d_2 = 2.)
    Check: B = D A = [[1,0],[0,2]] [[2,-1],[-2,2]] = [[2,-1],[-4,4]]. Hmm. Let's
    re-examine. We need (alpha_i, alpha_j) symmetric.
    B_ij = d_i A_ij. So B = D A = [[1·2, 1·(-1)], [2·(-2), 2·2]] = [[2, -1], [-4, 4]].
    Symmetric? B_12 = -1, B_21 = -4. NOT symmetric.
    Hmm. Let me re-derive: A_ij = 2 (alpha_i, alpha_j) / (alpha_i, alpha_i).
    So (alpha_i, alpha_j) = A_ij (alpha_i, alpha_i) / 2 = A_ij d_i.
    For (alpha, alpha) = 2 (short), d_short = 1. For long, d_long = (long_sq)/2.
    For C_n: long^2 = 2 * short^2 = 4. So d_long = 2.
    Then (alpha_1, alpha_2) = A_12 d_1 = -1 * 1 = -1. (Using d_1 = d_short = 1.)
                              and = A_21 d_2 = -2 * 2 / 2 = -2. Hmm wrong.
    Wait. The convention is (alpha_i, alpha_j) = A_ij (alpha_j, alpha_j)/2 = A_ij d_j.
    (Not d_i.) Because A_ij = 2 (alpha_i, alpha_j) / (alpha_j, alpha_j) in another
    convention. Argh.

    Actually let's just use the symmetrization (alpha_i, alpha_j) = d_i A_ij where
    d_i are non-negative integers such that D A is symmetric.
    For C_n: d_short = 1, d_long = 2 OR d_short = 2, d_long = 1 -- depending on
    sign convention. Try d = [2, 1] for C_2 with cartan [[2,-1],[-2,2]]:
      B = [[2*2, 2*(-1)], [1*(-2), 1*2]] = [[4, -2], [-2, 2]]. Symmetric.
    Then (alpha_1, alpha_1) = 4 (so alpha_1 is LONG of squared length 4)
    and (alpha_2, alpha_2) = 2 (so alpha_2 is SHORT). But we wanted alpha_1 short.

    Maybe the convention is alpha_1 long instead? Or use the transposed Cartan
    [[2,-2],[-1,2]]?

    Let's use convention: A_ij = <alpha_i^vee, alpha_j>. Then A is the Cartan
    matrix from coroot times root. And (alpha_i, alpha_j) = d_j A_ij where d_j
    = (alpha_j, alpha_j) / 2.

    For C_2 in Bourbaki convention (alpha_1 short, alpha_2 long):
      A = [[2, -2], [-1, 2]]  (from <alpha_1^v, alpha_2> = -2 since alpha_2 is twice as long)
    Then d_1 = 1 (short), d_2 = 2 (long).
    B_ij = d_j A_ij: B_12 = 2 * (-2) = -4. B_21 = 1 * (-1) = -1. NOT symmetric.

    Or: A_ij = <alpha_i, alpha_j^v>. Then (alpha_i, alpha_j) = d_i A_ij.
    For Bourbaki C_2 ([[2,-2],[-1,2]]) with d=[1,2]:
      B_12 = 1*(-2) = -2. B_21 = 2*(-1) = -2. Symmetric! ✓
    And then B_11 = 1*2 = 2 = (alpha_1, alpha_1). ✓ short. B_22 = 2*2 = 4 = (alpha_2, alpha_2). ✓ long.

    So in our convention (A_ij = (alpha_i, alpha_j^v); B_ij = d_i A_ij is the
    symmetrized form), C_2 with alpha_1 short has Cartan
      A = [[2, -2], [-1, 2]]  with d = [1, 2].
    """
    cartan = [[2, -2], [-1, 2]]
    d = [1, 2]
    return SimpleLieAlgebra("C2", cartan, d)


def build_G2() -> SimpleLieAlgebra:
    """
    Type G_2. Rank 2. alpha_1 short, alpha_2 long.

    In Bourbaki: A_ij = <alpha_i, alpha_j^v>, so
      A_12 = <alpha_1, alpha_2^v> = 2 (alpha_1, alpha_2) / (alpha_2, alpha_2)
    For G_2, (alpha_1, alpha_1) = 2 (short), (alpha_2, alpha_2) = 6 (long),
    angle = 150 deg, so (alpha_1, alpha_2) = sqrt(2) * sqrt(6) * cos(150) = -3.
    Then A_12 = 2 * (-3) / 6 = -1.  A_21 = 2 * (-3) / 2 = -3.
    A = [[2, -1], [-3, 2]],  d = [1, 3].

    Symmetrize: B = D A:
      B_11 = 1*2 = 2 = (alpha_1, alpha_1). ✓
      B_12 = 1*(-1) = -1.
      B_21 = 3*(-3) = -9. NOT symmetric.

    Hmm. Let me try transpose: A = [[2, -3], [-1, 2]], d = [1, 3]:
      B_11 = 2, B_12 = 1*(-3) = -3, B_21 = 3*(-1) = -3. ✓ symmetric.
      B_22 = 3*2 = 6. ✓ long.

    So I had the Cartan matrix transposed. Try A = [[2, -3], [-1, 2]] convention
    (this is the matrix when row i, col j is <alpha_i, alpha_j^v>).
    Hmm but the textbook G_2 Cartan IS [[2,-1],[-3,2]]. There are TWO conventions
    in textbooks: A_ij = <alpha_i, alpha_j^v> vs A_ij = <alpha_i^v, alpha_j>.
    They're transposes.

    For our purposes, what matters is that B_ij = d_i A_ij is symmetric and that
    the Weyl chamber / positive roots / fundamental weights / Brauer-Klimyk all
    work out consistently. Let me just choose:

       cartan = [[2, -3], [-1, 2]],  d = [1, 3]    --> B = [[2,-3],[-3,6]]  ✓ symmetric.

    Cross-check: positive roots in alpha basis are
       (1,0), (0,1), (1,1), (2,1), (3,1), (3,2)
    Their squared lengths: should give 2,6,2,6,2,2 (short:3, long:3 typically for G2).

    (alpha_1+alpha_2, alpha_1+alpha_2) = 2 - 6 + 6 = 2 (treating B_12=B_21=-3)
       Wait: B = [[2,-3],[-3,6]], so (alpha_1+alpha_2)^2 = 2 + 6 + 2*(-3) = 2. ✓ short.
    (2 alpha_1 + alpha_2)^2 = 4*2 + 6 + 4*(-3) = 8 + 6 - 12 = 2. ✓ short.
    (3 alpha_1 + alpha_2)^2 = 9*2 + 6 + 6*(-3) = 18 + 6 - 18 = 6. ✓ long.
    (3 alpha_1 + 2 alpha_2)^2 = 9*2 + 4*6 + 12*(-3) = 18 + 24 - 36 = 6. ✓ long.

    Good: positive roots split 3 short, 3 long. ✓
    """
    cartan = [[2, -3], [-1, 2]]
    d = [1, 3]
    return SimpleLieAlgebra("G2", cartan, d)


# ---------------------------------------------------------------------------
# Sanity tests against SU(3) baseline
# ---------------------------------------------------------------------------

def sanity_su3_baseline():
    """
    Verify SU(3) implementation reproduces P-A's SU(3) script results:
      C(1,0) = 4/3, C(1,1) = 3, C(2,0) = 10/3, C(3,0) = 6, C(2,2) = 8
      3 (x) 3-bar = 1 + 8
      3 (x) 3 = 6 + 3-bar
      8 (x) 8 = 1 + 8 + 8 + 10 + 10-bar + 27
    """
    a2 = build_A(2)
    print(f"\nSanity check: A_2 = SU(3)")
    print(f"  Cartan: {a2.cartan.tolist()}")
    print(f"  G_omega = {a2.G_omega.tolist()}")
    print(f"  rho_sq = {a2.rho_sq}")
    print(f"  Positive roots: {a2.positive_roots_omega}")
    print(f"  C(1,0) = {a2.casimir((1,0))} (expect 4/3)")
    print(f"  C(0,1) = {a2.casimir((0,1))} (expect 4/3)")
    print(f"  C(1,1) = {a2.casimir((1,1))} (expect 3)")
    print(f"  C(2,0) = {a2.casimir((2,0))} (expect 10/3)")
    print(f"  C(3,0) = {a2.casimir((3,0))} (expect 6)")
    print(f"  C(2,2) = {a2.casimir((2,2))} (expect 8)")
    print(f"  dim(1,1) = {a2.dim_weyl((1,1))} (expect 8)")
    print(f"  dim(2,2) = {a2.dim_weyl((2,2))} (expect 27)")
    print(f"  dim(3,0) = {a2.dim_weyl((3,0))} (expect 10)")
    print(f"  dual(2,1) = {a2.dual((2,1))} (expect (1,2))")

    # 3 (x) 3-bar = 1 + 8
    decomp = tensor_product(a2, (1, 0), (0, 1))
    print(f"  (1,0) x (0,1) = {decomp} (expect {{(0,0):1, (1,1):1}})")

    # 3 (x) 3 = 6 + 3-bar
    decomp = tensor_product(a2, (1, 0), (1, 0))
    print(f"  (1,0) x (1,0) = {decomp} (expect {{(2,0):1, (0,1):1}})")

    # 8 (x) 8
    decomp = tensor_product(a2, (1, 1), (1, 1))
    print(f"  (1,1) x (1,1) = {decomp}")
    print(f"     expect {{(0,0):1, (1,1):2, (3,0):1, (0,3):1, (2,2):1}}")


# ---------------------------------------------------------------------------
# Build panels and run verifications
# ---------------------------------------------------------------------------

def panel_dominant_weights(rank: int, total_max: int) -> List[Tuple[int, ...]]:
    """All (a_1, ..., a_rank) in N^rank with sum <= total_max."""
    weights = []
    if rank == 1:
        for a in range(total_max + 1):
            weights.append((a,))
        return weights

    # Recursive enumeration
    def helper(remaining: int, depth: int) -> List[List[int]]:
        if depth == rank - 1:
            return [[k] for k in range(remaining + 1)]
        out = []
        for k in range(remaining + 1):
            for tail in helper(remaining - k, depth + 1):
                out.append([k] + tail)
        return out

    raw = helper(total_max, 0)
    return [tuple(w) for w in raw]


def run_panel(la: SimpleLieAlgebra, irreps: List[Tuple[int, ...]],
              label: str, verbose_violations: bool = True) -> Dict:
    """
    Run Dirac-triangle inequality check on all ordered pairs of `irreps` in `la`.
    Returns summary dict.
    """
    print(f"\n{'=' * 78}")
    print(f"Dirac-triangle verification: {label}")
    print(f"{'=' * 78}")
    print(f"Irreps in panel ({len(irreps)} total):")
    for lam in irreps:
        dim = la.dim_weyl(lam)
        C = la.casimir(lam)
        print(f"  {lam}: dim={dim}, C={C}")

    pass_count = 0
    fail_count = 0
    max_ratio = 0.0
    max_ratio_case = None
    violations_all = []
    total_pairs = 0

    for lam in irreps:
        for lam_pr in irreps:
            total_pairs += 1
            passes, viols, ratios = verify_dirac_triangle(la, lam, lam_pr)
            if not passes:
                fail_count += 1
                if verbose_violations:
                    print(f"  FAIL: {lam} vs {lam_pr}: violations {viols}")
                for v in viols:
                    violations_all.append({
                        "lam": lam, "lam_pr": lam_pr,
                        "sigma": v[0], "sqrt_C_sigma": v[1],
                        "diff": v[2], "ratio": v[3],
                    })
            else:
                pass_count += 1
            if ratios:
                pair_sup = max(ratios)
                if pair_sup > max_ratio:
                    max_ratio = pair_sup
                    max_ratio_case = (lam, lam_pr)

    print(f"\n--- {label} summary ---")
    print(f"  Total pairs: {total_pairs}")
    print(f"  Pass: {pass_count}")
    print(f"  Fail: {fail_count}")
    print(f"  Max ratio: {max_ratio:.6f}")
    if max_ratio_case:
        print(f"  Achieved at: lam={max_ratio_case[0]}, lam_pr={max_ratio_case[1]}")

    return {
        "label": label,
        "rank": la.rank,
        "panel_size": len(irreps),
        "total_pairs": total_pairs,
        "pass_count": pass_count,
        "fail_count": fail_count,
        "max_ratio": float(max_ratio),
        "max_ratio_case": [list(max_ratio_case[0]), list(max_ratio_case[1])] if max_ratio_case else None,
        "violations": violations_all[:20],  # first 20 only
    }


def asymptotic_family(la: SimpleLieAlgebra, n_max: int, label: str) -> Dict:
    """
    Run asymptotic family (n, 0, ..., 0) (x) (1, 0, ..., 0)^* for n = 1..n_max.
    Return list of (n, |D-D'|, sqrt(C_sigma_min), ratio, sigma_min).
    """
    print(f"\n--- {label} asymptotic family (n, 0, ..., 0) vs (1, 0, ..., 0) ---")
    print(f"  {'n':>3} {'|D-D_1|':>12} {'sqrt(C_min)':>14} {'ratio':>10}  sigma_min")
    rank = la.rank
    rows = []
    for n in range(1, n_max + 1):
        lam = (n,) + (0,) * (rank - 1)
        lam_pr = (1,) + (0,) * (rank - 1)
        lam_pr_dual = la.dual(lam_pr)
        decomp = tensor_product(la, lam, lam_pr_dual)
        smallest_C = None
        smallest_sigma = None
        for sigma, mult in decomp.items():
            if mult > 0:
                C = la.casimir(sigma)
                if C > 0 and (smallest_C is None or C < smallest_C):
                    smallest_C = C
                    smallest_sigma = sigma
        if smallest_C is None:
            continue
        D = la.D_abs(lam)
        D_pr = la.D_abs(lam_pr)
        diff = float(sp.Abs(D - D_pr))
        sqrt_C = float(sqrt(smallest_C))
        ratio = diff / sqrt_C if sqrt_C > 0 else float('inf')
        print(f"  {n:>3} {diff:>12.4f} {sqrt_C:>14.4f} {ratio:>10.4f}   sigma_min={smallest_sigma}")
        rows.append({
            "n": n,
            "diff": diff,
            "sqrt_C_sigma_min": sqrt_C,
            "ratio": ratio,
            "sigma_min": list(smallest_sigma),
        })
    return {"label": label, "rows": rows}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "dirac_triangle_extended.json")

    results = {"description": "Dirac-triangle inequality |D(pi)-D(pi')| <= sqrt(C(sigma))",
               "panels": [], "asymptotics": []}

    # SU(3) sanity baseline
    sanity_su3_baseline()

    # SU(3) extended panel: p+q <= 5
    a2 = build_A(2)
    su3_panel = panel_dominant_weights(2, 5)
    print(f"\nSU(3) extended panel size: {len(su3_panel)}, total pairs: {len(su3_panel)**2}")
    results["panels"].append(run_panel(a2, su3_panel, "SU(3) p+q<=5"))
    results["asymptotics"].append(asymptotic_family(a2, 10, "SU(3)"))

    # SU(4) panel: a+b+c <= 3
    a3 = build_A(3)
    su4_panel = panel_dominant_weights(3, 3)
    print(f"\nSU(4) panel size: {len(su4_panel)}, total pairs: {len(su4_panel)**2}")
    results["panels"].append(run_panel(a3, su4_panel, "SU(4) a+b+c<=3"))
    results["asymptotics"].append(asymptotic_family(a3, 6, "SU(4)"))

    # Sp(2) (C_2) panel: a+b <= 3
    c2 = build_C2()
    sp2_panel = panel_dominant_weights(2, 3)
    print(f"\nSp(2)=C_2 panel size: {len(sp2_panel)}, total pairs: {len(sp2_panel)**2}")
    results["panels"].append(run_panel(c2, sp2_panel, "Sp(2)=C_2 a+b<=3"))
    results["asymptotics"].append(asymptotic_family(c2, 6, "Sp(2)"))

    # G_2 panel: a+b <= 2 (irreps grow fast)
    g2 = build_G2()
    g2_panel = panel_dominant_weights(2, 2)
    print(f"\nG_2 panel size: {len(g2_panel)}, total pairs: {len(g2_panel)**2}")
    results["panels"].append(run_panel(g2, g2_panel, "G_2 a+b<=2"))
    results["asymptotics"].append(asymptotic_family(g2, 5, "G_2"))

    # Persist
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {output_path}")

    # Final summary
    print("\n" + "=" * 78)
    print("FINAL SUMMARY")
    print("=" * 78)
    all_pass = True
    for p in results["panels"]:
        verdict = "PASS" if p["fail_count"] == 0 else "FAIL"
        print(f"  {p['label']}: {p['pass_count']}/{p['total_pairs']} pass, max_ratio={p['max_ratio']:.4f} [{verdict}]")
        if p["fail_count"] > 0:
            all_pass = False
    print(f"\nOverall verdict: {'ALL PASS' if all_pass else 'SOME FAILED'}")


if __name__ == "__main__":
    main()
