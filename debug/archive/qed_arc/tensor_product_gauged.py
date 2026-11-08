"""
Gauged tensor-product spectral action: extracting the inter-particle potential
==============================================================================
Step 2: The free D_total² factorizes (Step 1 confirmed). The interaction
enters through gauging: D → D + A where A = Σ aᵢ[D, bᵢ].

The gauged operator squared:
  (D + A)² = D² + DA + AD + A²
           = D² + {D, A} + A²

D² factorizes. The interaction potential is in {D, A} + A².

For the tensor product with A built from algebra elements a₁ ⊗ a₂:
  A = Σ (a₁ᵢ ⊗ a₂ᵢ) [D_total, (b₁ᵢ ⊗ b₂ᵢ)]

The spectral action Tr(f((D+A)²/Λ²)) then contains:
  - Kinetic terms (from D²): single-particle energies
  - Gauge kinetic terms (from A²): propagator of the mediating field
  - INTERACTION terms (from {D, A}): the two-body potential

We extract the interaction by computing:
  V_interaction = Tr(f((D+A)²/Λ²)) - Tr(f(D²/Λ²))
at leading order in A.
"""
import sys, json, time
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.operator_system import build_multiplier_matrix, HyperLabel, allowed_multiplier_labels
from debug.tensor_product_dirac import build_single_particle_scalar, build_tensor_product_dirac


def build_gauge_field(sp, n_max):
    """Build the gauge field A = Σ aᵢ[D, bᵢ] from the algebra generators.

    For a single spectral triple, the gauge field (one-form) is:
      A = Σ aᵢ [D, bᵢ]
    where aᵢ, bᵢ ∈ A (the algebra).

    For the TENSOR PRODUCT, the gauge field couples the two factors:
      A_total = Σ (a₁ᵢ ⊗ a₂ᵢ) [D_total, (b₁ᵢ ⊗ b₂ᵢ)]

    The simplest nontrivial gauge field: use a single pair (a, b)
    where a and b are 3-Y multiplier matrices.

    We build A for the full tensor product and decompose it into
    single-particle and cross-particle parts.
    """
    D = sp['D']
    dim = sp['dim']
    N_scalar = sp['N_scalar']
    scalar_basis = sp['scalar_basis']

    # Get all allowed multiplier labels for the operator system
    labels = allowed_multiplier_labels(n_max)

    # Build all multiplier matrices (in chirality-doubled space)
    multipliers = []
    for (N, L, M) in labels:
        M_scalar = build_multiplier_matrix(N, L, M, scalar_basis)
        # Take real part (imaginary part cancels in Hermitian combinations)
        M_real = M_scalar.real
        if np.linalg.norm(M_real) < 1e-15:
            continue
        # Embed in chirality-doubled space
        M_full = np.zeros((dim, dim))
        M_full[:N_scalar, :N_scalar] = M_real
        M_full[N_scalar:, N_scalar:] = M_real
        multipliers.append({
            'label': (N, L, M),
            'M_scalar': M_real,
            'M_full': M_full,
        })

    # Build the single-particle gauge field: A₁ = Σ Mᵢ [D, Mⱼ]
    # This is the "generic" gauge field using all multiplier pairs
    # For the diagnostic, use A = M₁ [D, M₂] for a specific pair
    # that maximizes the commutator norm

    commutators = []
    for mult in multipliers:
        comm = D @ mult['M_full'] - mult['M_full'] @ D
        norm = np.linalg.norm(comm)
        if norm > 1e-10:
            commutators.append({
                **mult,
                'comm': comm,
                'comm_norm': norm,
            })

    commutators.sort(key=lambda x: -x['comm_norm'])

    return multipliers, commutators


def compute_gauged_spectral_action(sp1, sp2, n_max):
    """Compute the gauged spectral action and extract the interaction.

    Strategy:
    1. Build A_total from the tensor-product inner fluctuations
    2. Compute (D + εA)² for small ε
    3. Extract the ε-linear term {D, A} (first-order interaction)
    4. Extract the ε²-quadratic term A² (gauge kinetic / second-order)
    5. Compute the heat trace Tr(e^{-t(D+εA)²}) and compare with factorized
    """
    d1, d2 = sp1['dim'], sp2['dim']
    dim_total = d1 * d2
    D1, D2 = sp1['D'], sp2['D']
    gamma1 = sp1['gamma']
    I1, I2 = np.eye(d1), np.eye(d2)

    # Build D_total
    D_total = np.kron(D1, I2) + np.kron(gamma1, D2)

    # Build gauge fields for each particle
    print(f"  Building gauge fields...")
    mults1, comms1 = build_gauge_field(sp1, n_max)
    mults2, comms2 = build_gauge_field(sp2, n_max)

    print(f"  {len(mults1)} multiplier matrices, {len(comms1)} with nonzero [D, M]")
    if comms1:
        print(f"  Largest ||[D, M]|| = {comms1[0]['comm_norm']:.4f} "
              f"at label {comms1[0]['label']}")

    if not comms1 or not comms2:
        print("  No nonzero commutators — cannot build gauge field")
        return None

    # Build a representative tensor-product gauge field
    # A_total = (a₁ ⊗ a₂) [D_total, (b₁ ⊗ b₂)]
    # Use the strongest commutator pair for each particle
    a1 = comms1[0]['M_full']
    comm1 = comms1[0]['comm']  # [D₁, b₁]
    a2 = comms2[0]['M_full']
    comm2 = comms2[0]['comm']  # [D₂, b₂]

    # The tensor-product inner fluctuation:
    # [D_total, b₁ ⊗ b₂] = [D₁, b₁] ⊗ b₂ + γ₁·b₁ ⊗ [D₂, b₂]
    b1_full = comms1[0]['M_full']
    b2_full = comms2[0]['M_full']

    fluct_term1 = np.kron(comm1, b2_full)
    fluct_term2 = np.kron(gamma1 @ b1_full, comm2)
    fluct = fluct_term1 + fluct_term2

    # A_total = (a₁ ⊗ a₂) · fluctuation
    a_tensor = np.kron(a1, a2)
    A_total = a_tensor @ fluct

    # Make Hermitian: A_phys = (A + A†)/2
    A_total = (A_total + A_total.T) / 2

    print(f"\n  ||A_total|| = {np.linalg.norm(A_total):.4f}")
    print(f"  ||D_total|| = {np.linalg.norm(D_total):.4f}")
    print(f"  ||A|| / ||D|| = {np.linalg.norm(A_total)/np.linalg.norm(D_total):.4f}")

    # Compute (D + εA)² = D² + ε{D,A} + ε²A²
    D_sq = D_total @ D_total
    DA_anticomm = D_total @ A_total + A_total @ D_total  # {D, A}
    A_sq = A_total @ A_total

    print(f"\n  ||D²|| = {np.linalg.norm(D_sq):.4f}")
    print(f"  ||{{D, A}}|| = {np.linalg.norm(DA_anticomm):.4f}")
    print(f"  ||A²|| = {np.linalg.norm(A_sq):.4f}")

    # THE KEY TEST: does {D, A} factorize?
    # If {D, A} = X₁ ⊗ I₂ + I₁ ⊗ X₂ for some X₁, X₂, there's no interaction.
    # If {D, A} has a genuine cross-term (not factorizable), there IS interaction.

    # Partial trace test: if {D, A} factorizes as X₁⊗I + I⊗X₂, then
    # Tr₂({D,A}) = d₂ · X₁ + Tr₂(I₂) · X₂ = d₂ · X₁ + d₂ · X₂
    # and the "connected" part {D,A} - (1/d₂)Tr₂({D,A})⊗I - I⊗(1/d₁)Tr₁({D,A})
    # + (1/(d₁d₂))Tr({D,A})·I⊗I should be zero for factorizable.

    # Reshape {D,A} as (d1, d2, d1, d2) tensor
    DA_tensor = DA_anticomm.reshape(d1, d2, d1, d2)

    # Partial traces
    Tr2_DA = np.trace(DA_tensor, axis1=1, axis2=3)  # (d1, d1) matrix
    Tr1_DA = np.trace(DA_tensor, axis1=0, axis2=2)  # (d2, d2) matrix
    full_trace = np.trace(DA_anticomm)

    # Factorizable part: (1/d₂)Tr₂ ⊗ I₂ + I₁ ⊗ (1/d₁)Tr₁ - (1/(d₁d₂))Tr · I⊗I
    factorizable = (np.kron(Tr2_DA / d2, I2)
                    + np.kron(I1, Tr1_DA / d1)
                    - (full_trace / (d1 * d2)) * np.eye(dim_total))

    # Connected (non-factorizable) part = the interaction!
    connected = DA_anticomm - factorizable

    print(f"\n  --- Factorization test of {{D, A}} ---")
    print(f"  ||{{D,A}}|| = {np.linalg.norm(DA_anticomm):.6f}")
    print(f"  ||factorizable part|| = {np.linalg.norm(factorizable):.6f}")
    print(f"  ||connected part|| = {np.linalg.norm(connected):.6f}")
    ratio = np.linalg.norm(connected) / np.linalg.norm(DA_anticomm) if np.linalg.norm(DA_anticomm) > 0 else 0
    print(f"  Connected fraction = {ratio:.4f}")

    if ratio > 0.01:
        print(f"  ✓ CONNECTED PART IS NONZERO — interaction exists!")
    else:
        print(f"  ✗ Connected part negligible — no interaction at this order")

    # Same test for A²
    A_sq_tensor = A_sq.reshape(d1, d2, d1, d2)
    Tr2_Asq = np.trace(A_sq_tensor, axis1=1, axis2=3)
    Tr1_Asq = np.trace(A_sq_tensor, axis1=0, axis2=2)
    full_trace_Asq = np.trace(A_sq)
    factorizable_Asq = (np.kron(Tr2_Asq / d2, I2)
                         + np.kron(I1, Tr1_Asq / d1)
                         - (full_trace_Asq / (d1 * d2)) * np.eye(dim_total))
    connected_Asq = A_sq - factorizable_Asq
    ratio_Asq = np.linalg.norm(connected_Asq) / np.linalg.norm(A_sq) if np.linalg.norm(A_sq) > 0 else 0

    print(f"\n  --- Factorization test of A² ---")
    print(f"  ||A²|| = {np.linalg.norm(A_sq):.6f}")
    print(f"  ||connected part|| = {np.linalg.norm(connected_Asq):.6f}")
    print(f"  Connected fraction = {ratio_Asq:.4f}")

    # Heat trace at several ε values
    print(f"\n  --- Gauged heat trace ---")
    print(f"  {'eps':>6} {'Tr(e^-t(D+eA)^2)':>20} {'factorized':>12} {'ratio':>10}")

    t_test = 0.5
    evals_D_sq = np.linalg.eigvalsh(D_sq)
    ht_free = np.sum(np.exp(-t_test * evals_D_sq))
    ht_factorized_1 = np.sum(np.exp(-t_test * np.diag(D1 @ D1)))
    ht_factorized_2 = np.sum(np.exp(-t_test * np.diag(D2 @ D2)))
    ht_factorized = ht_factorized_1 * ht_factorized_2

    for eps in [0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]:
        D_gauged = D_total + eps * A_total
        D_gauged_sq = D_gauged @ D_gauged
        evals_gauged = np.linalg.eigvalsh(D_gauged_sq)
        ht_gauged = np.sum(np.exp(-t_test * evals_gauged))
        r = ht_gauged / ht_factorized if ht_factorized > 0 else float('nan')
        print(f"  {eps:>6.2f} {ht_gauged:>20.6f} {ht_factorized:>12.6f} {r:>10.6f}")

    return {
        'DA_anticomm_norm': float(np.linalg.norm(DA_anticomm)),
        'connected_norm': float(np.linalg.norm(connected)),
        'connected_fraction': float(ratio),
        'A_sq_connected_fraction': float(ratio_Asq),
        'A_norm': float(np.linalg.norm(A_total)),
    }


def main():
    print("=" * 72)
    print("GAUGED TENSOR-PRODUCT SPECTRAL ACTION")
    print("Extracting the inter-particle potential")
    print("=" * 72)

    for n_max in [2, 3]:
        print(f"\n{'='*60}")
        print(f"n_max = {n_max}")
        print(f"{'='*60}")

        t0 = time.time()
        sp1 = build_single_particle_scalar(n_max)
        sp2 = build_single_particle_scalar(n_max)

        result = compute_gauged_spectral_action(sp1, sp2, n_max)
        print(f"  Total time: {time.time()-t0:.1f} s")

    print(f"\n{'='*72}")
    print("INTERPRETATION")
    print(f"{'='*72}")
    print("""
  The factorization test decomposes {D, A} and A² into:
    - Factorizable part: X₁⊗I + I⊗X₂ (single-particle, no interaction)
    - Connected part: genuinely two-body (the interaction)

  If the connected fraction is nonzero, the gauged spectral action
  generates an inter-particle potential. The angular structure of
  the connected part should match the Gaunt selection rules from
  Paper 22 (angular sparsity theorem). The radial structure is
  the NEW prediction that we compare against 1/r₁₂.

  The heat trace at ε > 0 shows how the spectral action departs
  from the factorized (free) value as the gauge field strength grows.
""")


if __name__ == '__main__':
    main()
