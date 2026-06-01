"""
Tensor-product Dirac operator for two particles on S³ ⊗ S³
============================================================
Step 1 of the head-on program: build D_total = D₁ ⊗ I₂ + γ₁ ⊗ D₂
at finite n_max and compute inner fluctuations.

The key question: do the cross-modulation terms
    [D_total, a₁ ⊗ a₂] = [D₁, a₁] ⊗ a₂ + γ₁·a₁ ⊗ [D₂, a₂]
generate an inter-particle potential when evaluated via the spectral action?

We use the SCALAR operator system (not spinor) for the initial diagnostic,
since the scalar Fock-projected S³ is the framework's core object and
the algebra structure (3-Y integrals, Gaunt selection rules) is cleanest there.

For the scalar case:
  - D₁ = Camporesi-Higuchi eigenvalues ±(n + 1/2) on the scalar basis
  - γ₁ = chirality grading (σ_x in the Weyl/anti-Weyl doubling)
  - The algebra O_{n_max} = P_{n_max} C(S³) P_{n_max} acts by 3-Y multipliers

The tensor product D_total acts on H₁ ⊗ H₂ = (2·dim₁) ⊗ (2·dim₂)
where the factor of 2 is the chirality doubling.
"""
import sys, json, time
import numpy as np
from itertools import product as iterproduct

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.full_dirac_operator_system import (
    full_dirac_basis, FullDiracLabel,
    camporesi_higuchi_full_dirac_matrix,
)
from geovac.operator_system import TruncatedOperatorSystem


def build_single_particle_scalar(n_max):
    """Build the single-particle system in the SCALAR sector.

    Uses the scalar Fock-projected basis |n, l, m> with a chirality doubling
    (dim = 2 * N_scalar). The Dirac is diagonal with eigenvalues ±(n + 1/2).
    The chirality γ swaps the two copies.

    This is simpler than the full spinor Dirac but captures the essential
    tensor-product structure.
    """
    from geovac.operator_system import HyperLabel

    # Scalar basis
    scalar_basis = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                scalar_basis.append(HyperLabel(n=n, l=l, m=m))

    N_scalar = len(scalar_basis)
    dim = 2 * N_scalar  # chirality doubling

    # Dirac: diagonal, +(n+1/2) on first half, -(n+1/2) on second half
    D = np.zeros((dim, dim))
    for i, b in enumerate(scalar_basis):
        eigenval = b.n + 0.5
        D[i, i] = +eigenval          # positive chirality
        D[N_scalar + i, N_scalar + i] = -eigenval  # negative chirality

    # Chirality γ: swaps the two sectors
    gamma = np.zeros((dim, dim))
    gamma[:N_scalar, N_scalar:] = np.eye(N_scalar)
    gamma[N_scalar:, :N_scalar] = np.eye(N_scalar)

    return {
        'scalar_basis': scalar_basis,
        'N_scalar': N_scalar,
        'dim': dim,
        'D': D,
        'gamma': gamma,
        'n_max': n_max,
    }


def build_tensor_product_dirac(sp1, sp2):
    """Build D_total = D₁ ⊗ I₂ + γ₁ ⊗ D₂ on H₁ ⊗ H₂.

    This is the standard tensor-product Dirac for a product of spectral triples.
    """
    d1, d2 = sp1['dim'], sp2['dim']
    dim_total = d1 * d2

    D1 = sp1['D']
    D2 = sp2['D']
    gamma1 = sp1['gamma']
    I1 = np.eye(d1)
    I2 = np.eye(d2)

    # D_total = D₁ ⊗ I₂ + γ₁ ⊗ D₂
    D_total = np.kron(D1, I2) + np.kron(gamma1, D2)

    return D_total, dim_total


def build_algebra_generators(n_max):
    """Build the algebra generators for the truncated operator system.

    Each generator M_{NLM} is the 3-Y multiplier matrix:
    (M_{NLM})_{(n,l,m),(n',l',m')} = <n,l,m|Y^{(3)}_{NLM}|n',l',m'>

    These are the SCALAR operator system generators (no chirality doubling).
    For the full Dirac, we embed them as a₁ = M ⊗ I_chirality (acting on
    both chirality sectors identically).
    """
    op_sys = TruncatedOperatorSystem(n_max)

    # The generators are already in the operator system
    # Get the basis matrices
    generators = op_sys._generators  # list of (N, L, M, matrix) tuples

    return op_sys, generators


def compute_inner_fluctuations(sp1, sp2, D_total, n_max):
    """Compute the inner fluctuation [D_total, a₁ ⊗ a₂] for algebra generators.

    For a = a₁ ⊗ a₂:
    [D_total, a₁ ⊗ a₂] = [D₁, a₁] ⊗ a₂ + γ₁·a₁ ⊗ [D₂, a₂]

    The first term is a "gauge field on particle 1 modulated by particle 2's state."
    The second term is the reverse with a chirality twist.

    We want to know:
    1. Are the cross terms nonzero? (Do the two particles interact?)
    2. What is the structure of the resulting operator? (Angular decomposition)
    3. Does summing a_i [D, b_i] over generators produce a recognizable potential?
    """
    d1, d2 = sp1['dim'], sp2['dim']
    D1, D2 = sp1['D'], sp2['D']
    gamma1 = sp1['gamma']
    I1, I2 = np.eye(d1), np.eye(d2)
    N_scalar = sp1['N_scalar']
    half = N_scalar

    print(f"  Single-particle: n_max={n_max}, dim_scalar={N_scalar}, dim_chirality_doubled={d1}")
    print(f"  Tensor product: dim_total = {d1*d2}")

    print(f"\n  --- Inner fluctuation diagnostic ---")

    # Build a 3-Y multiplier matrix (the algebra generator)
    from geovac.operator_system import build_multiplier_matrix, HyperLabel

    scalar_basis = sp1['scalar_basis']

    N_Y, L_Y, M_Y = 2, 1, 0
    M_scalar = build_multiplier_matrix(N_Y, L_Y, M_Y, scalar_basis).real

    # Embed into chirality-doubled space: [[M, 0], [0, M]]
    M_full = np.zeros((d1, d1))
    M_full[:half, :half] = M_scalar
    M_full[half:, half:] = M_scalar

    comm_D1_M = D1 @ M_full - M_full @ D1
    print(f"\n  3-Y multiplier (N=2,L=1,M=0):")
    print(f"  ||M_scalar|| = {np.linalg.norm(M_scalar):.6f}")
    print(f"  ||[D₁, M]|| = {np.linalg.norm(comm_D1_M):.6f}")
    print(f"  [D₁, M] ≠ 0: {np.linalg.norm(comm_D1_M) > 1e-10}")

    # Now compute the TENSOR PRODUCT inner fluctuation
    # [D_total, M₁ ⊗ I₂] = [D₁, M₁] ⊗ I₂ + γ₁·M₁ ⊗ [D₂, I₂]
    #                      = [D₁, M₁] ⊗ I₂    (since [D₂, I₂] = 0)
    # This is NOT interesting — it's just a gauge field on particle 1, no coupling.

    # The INTERESTING case is:
    # [D_total, M₁ ⊗ M₂] = [D₁, M₁] ⊗ M₂ + γ₁·M₁ ⊗ [D₂, M₂]
    # Both terms involve BOTH particles. The second term has γ₁ (chirality twist).

    M2_full = np.zeros((d2, d2))
    M2_full[:half, :half] = M_scalar
    M2_full[half:, half:] = M_scalar

    # Full tensor product commutator
    term1 = np.kron(comm_D1_M, M2_full)           # [D₁, M₁] ⊗ M₂
    term2 = np.kron(gamma1 @ M_full, D2 @ M2_full - M2_full @ D2)  # γ₁·M₁ ⊗ [D₂, M₂]

    full_comm = term1 + term2
    print(f"\n  Tensor product [D_total, M₁ ⊗ M₂]:")
    print(f"  ||term1|| = ||[D₁,M₁] ⊗ M₂|| = {np.linalg.norm(term1):.6f}")
    print(f"  ||term2|| = ||γ₁·M₁ ⊗ [D₂,M₂]|| = {np.linalg.norm(term2):.6f}")
    print(f"  ||full_comm|| = {np.linalg.norm(full_comm):.6f}")

    # The GAUGE FIELD is A = Σ a_i [D, b_i]. For the tensor product,
    # A_total = Σ (a₁_i ⊗ a₂_i) [D_total, (b₁_i ⊗ b₂_i)]
    # The spectral action Tr(f(D_total + A_total)²/Λ²) then gives the
    # effective action including the inter-particle interaction.

    # Key diagnostic: the SPECTRUM of D_total
    evals_D = np.sort(np.linalg.eigvalsh(D_total))

    print(f"\n  Spectrum of D_total (first 20 eigenvalues):")
    for i in range(min(20, len(evals_D))):
        print(f"    λ_{i} = {evals_D[i]:.4f}")

    # Compare with D₁ ⊗ I₂ + I₁ ⊗ D₂ (NO chirality twist)
    D_no_gamma = np.kron(D1, I2) + np.kron(I1, D2)
    evals_no_gamma = np.sort(np.linalg.eigvalsh(D_no_gamma))

    print(f"\n  Spectrum of D₁⊗I + I⊗D₂ (no chirality, first 20):")
    for i in range(min(20, len(evals_no_gamma))):
        print(f"    λ_{i} = {evals_no_gamma[i]:.4f}")

    # The DIFFERENCE between D_total and D_no_gamma is the chirality twist.
    # D_total = D₁ ⊗ I + γ₁ ⊗ D₂  vs  D_no_gamma = D₁ ⊗ I + I ⊗ D₂
    # Difference = (γ₁ - I₁) ⊗ D₂
    # This is the INTERACTION KERNEL in some sense.

    delta_D = D_total - D_no_gamma
    print(f"\n  ||D_total - D_no_gamma|| = ||(γ₁ - I₁) ⊗ D₂|| = {np.linalg.norm(delta_D):.4f}")

    # Heat trace comparison: Tr(exp(-t*D_total²)) vs Tr(exp(-t*D₁²)) * Tr(exp(-t*D₂²))
    # If they differ, the tensor product structure generates an interaction
    print(f"\n  --- Heat trace diagnostic ---")
    D_total_sq = D_total @ D_total
    D1_sq = D1 @ D1
    D2_sq = D2 @ D2

    for t in [0.1, 0.5, 1.0, 2.0]:
        ht_total = np.trace(np.diag(np.exp(-t * np.diag(D_total_sq))))
        # Wait, D_total² is not diagonal. Need eigenvalues.
        evals_sq = np.linalg.eigvalsh(D_total_sq)
        ht_total = np.sum(np.exp(-t * evals_sq))

        ht_1 = np.sum(np.exp(-t * np.diag(D1_sq)))
        ht_2 = np.sum(np.exp(-t * np.diag(D2_sq)))
        ht_product = ht_1 * ht_2

        ratio = ht_total / ht_product if ht_product > 0 else float('nan')
        print(f"  t={t:.1f}: Tr(e^{{-tD²_total}}) = {ht_total:.4f}, "
              f"Tr(e^{{-tD²₁}}) × Tr(e^{{-tD²₂}}) = {ht_product:.4f}, "
              f"ratio = {ratio:.6f}")

    return {
        'comm_norm': np.linalg.norm(full_comm),
        'term1_norm': np.linalg.norm(term1),
        'term2_norm': np.linalg.norm(term2),
        'delta_D_norm': np.linalg.norm(delta_D),
    }


def main():
    print("=" * 72)
    print("TENSOR-PRODUCT DIRAC FOR TWO PARTICLES ON S³ ⊗ S³")
    print("Step 1: Inner fluctuation diagnostic")
    print("=" * 72)

    for n_max in [2, 3]:
        print(f"\n{'='*60}")
        print(f"n_max = {n_max}")
        print(f"{'='*60}")

        t0 = time.time()
        sp1 = build_single_particle_scalar(n_max)
        sp2 = build_single_particle_scalar(n_max)

        D_total, dim_total = build_tensor_product_dirac(sp1, sp2)
        print(f"  Build time: {time.time()-t0:.2f} s")

        results = compute_inner_fluctuations(sp1, sp2, D_total, n_max)

        print(f"\n  VERDICT for n_max={n_max}:")
        if results['comm_norm'] > 1e-10:
            print(f"  ✓ Inner fluctuation [D_total, a₁⊗a₂] is NONZERO")
            print(f"    → The tensor product generates inter-particle coupling")
        else:
            print(f"  ✗ Inner fluctuation is zero → no coupling")

    print(f"\n{'='*72}")
    print("INTERPRETATION")
    print(f"{'='*72}")
    print("""
  The heat trace ratio Tr(e^{-tD²_total}) / [Tr(e^{-tD²₁}) × Tr(e^{-tD²₂})]
  measures whether the tensor-product spectral action factorizes.

  Ratio = 1: the two particles are independent (no interaction)
  Ratio ≠ 1: the spectral action of D_total contains cross terms
             that represent inter-particle coupling

  The chirality twist γ₁ in D_total = D₁⊗I + γ₁⊗D₂ is what generates
  the non-factorization. Without it (D = D₁⊗I + I⊗D₂), the heat trace
  factorizes exactly and there is no interaction.

  The cross terms in the spectral action are the framework-native
  inter-particle potential. Their angular structure should match the
  Gaunt selection rules, and their radial structure is the NEW prediction
  that we compare against 1/r₁₂.
""")


if __name__ == '__main__':
    main()
