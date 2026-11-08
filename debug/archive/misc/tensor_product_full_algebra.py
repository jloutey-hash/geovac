"""
Full-algebra gauge field on the tensor product T_{S³} ⊗ T_{S³}
================================================================
The previous sprint used a SINGLE multiplier pair for the gauge field,
giving 65% Gaunt-compatible at n_max=3. This script sums over ALL
multiplier pairs to build the rotationally invariant gauge field.

In NCG, the most general gauge field (one-form) is:
    A = Σᵢ aᵢ [D, bᵢ]
where aᵢ, bᵢ range over the full algebra. For the truncated operator
system O_{n_max}, the algebra generators are the 3-Y multiplier matrices
M_{NLM}. The general one-form is:
    A = Σ_{(N₁L₁M₁), (N₂L₂M₂)} c_{12} · M_{N₁L₁M₁} [D, M_{N₂L₂M₂}]

For the diagnostic, we use c_{12} = 1 for all pairs (the "democratic"
gauge field). This sums over ALL angular channels equally, producing a
rotationally invariant gauge field.

For the TENSOR PRODUCT:
    A_total = Σ_{i,j,k,l} (M_i ⊗ M_j) [D_total, M_k ⊗ M_l]

This is expensive (n_mult⁴ terms), so we use the factored form:
    A_total = A₁ ⊗ I₂ + γ₁ ⊗ A₂  (single-particle gauge fields)
plus CROSS terms:
    A_cross = Σ_{i,j} (M_i ⊗ M_j) [D_total, M_i ⊗ M_j]
which generate the inter-particle coupling.
"""
import sys, json, time
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.operator_system import (
    build_multiplier_matrix, HyperLabel, allowed_multiplier_labels,
)
from debug.tensor_product_dirac import build_single_particle_scalar


def build_all_multipliers(n_max, scalar_basis, N_scalar, dim):
    """Build all multiplier matrices for the truncated operator system."""
    labels = allowed_multiplier_labels(n_max)
    multipliers = []
    for (N, L, M) in labels:
        M_sc = build_multiplier_matrix(N, L, M, scalar_basis)
        # Use the FULL complex matrix (both real and imaginary parts matter
        # for the general one-form; Hermiticity is imposed on A, not on
        # individual terms)
        if np.linalg.norm(M_sc) < 1e-15:
            continue
        M_full = np.zeros((dim, dim), dtype=complex)
        M_full[:N_scalar, :N_scalar] = M_sc
        M_full[N_scalar:, N_scalar:] = M_sc
        multipliers.append({
            'label': (N, L, M),
            'M_full': M_full,
        })
    return multipliers


def build_single_particle_oneform(D, multipliers):
    """Build the single-particle one-form A = Σᵢ Mᵢ [D, Mᵢ].

    Using a = b = Mᵢ (self-paired) for simplicity. The general case
    A = Σᵢⱼ Mᵢ [D, Mⱼ] would give the same span but with more freedom
    in the coefficients.
    """
    dim = D.shape[0]
    A = np.zeros((dim, dim), dtype=complex)
    for mult in multipliers:
        M = mult['M_full']
        comm = D @ M - M @ D
        A += M @ comm  # a_i [D, b_i] with a_i = b_i = M_i
    # Hermitianize
    A = (A + A.conj().T) / 2
    return A


def build_cross_gauge_field(sp1, sp2, multipliers, D_total, gamma1, max_pairs=None):
    """Build the rotationally invariant cross gauge field.

    For a rotationally invariant one-form, use the conjugate pairing:
        A = Σ_{NLM} M†_{NLM} · [D, M_{NLM}]
    where M† = M.conj().T (Hermitian conjugate of the multiplier).

    This ensures:
    - Hermiticity of A (automatic: (M†·[D,M])† = [D,M]†·M = -[D,M†]·M)
    - m-conservation (M†_{NLM} has quantum number -M, [D, M_{NLM}] has +M,
      so the product has net M=0 → conserves total m)

    For the tensor product:
        A_cross = Σ_{i,j} (M†ᵢ ⊗ M†ⱼ) · [D_total, (Mᵢ ⊗ Mⱼ)]
    """
    d1, d2 = sp1['dim'], sp2['dim']
    D1, D2 = sp1['D'], sp2['D']

    n_mult = len(multipliers)
    if max_pairs is not None:
        n_use = min(n_mult, max_pairs)
    else:
        n_use = n_mult

    A_cross = np.zeros((d1 * d2, d1 * d2), dtype=complex)
    n_terms = 0

    for i in range(n_use):
        Mi = multipliers[i]['M_full']
        Mi_dag = Mi.conj().T
        comm_i = D1 @ Mi - Mi @ D1
        if np.linalg.norm(comm_i) < 1e-14:
            continue

        for j in range(n_use):
            Mj = multipliers[j]['M_full']
            Mj_dag = Mj.conj().T
            comm_j = D2 @ Mj - Mj @ D2
            if np.linalg.norm(comm_j) < 1e-14:
                continue

            # [D_total, Mᵢ ⊗ Mⱼ] = [D₁, Mᵢ] ⊗ Mⱼ + γ₁·Mᵢ ⊗ [D₂, Mⱼ]
            fluct = np.kron(comm_i, Mj) + np.kron(gamma1 @ Mi, comm_j)

            # a†·[D, b] = (M†ᵢ ⊗ M†ⱼ) · fluct
            a_dag = np.kron(Mi_dag, Mj_dag)
            term = a_dag @ fluct
            A_cross += term
            n_terms += 1

    # Hermitianize (should already be close to Hermitian by construction)
    A_cross = (A_cross + A_cross.conj().T) / 2

    return A_cross.real, n_terms


def factorization_test(matrix, d1, d2, label=""):
    """Test whether a matrix on H₁⊗H₂ factorizes."""
    dim_total = d1 * d2
    I1, I2 = np.eye(d1), np.eye(d2)

    M_tensor = matrix.reshape(d1, d2, d1, d2)
    Tr2 = np.trace(M_tensor, axis1=1, axis2=3)
    Tr1 = np.trace(M_tensor, axis1=0, axis2=2)
    full_tr = np.trace(matrix)

    factorizable = (np.kron(Tr2 / d2, I2)
                    + np.kron(I1, Tr1 / d1)
                    - (full_tr / (d1 * d2)) * np.eye(dim_total))
    connected = matrix - factorizable

    norm_total = np.linalg.norm(matrix)
    norm_connected = np.linalg.norm(connected)
    frac = norm_connected / norm_total if norm_total > 0 else 0

    print(f"  {label}: ||total||={norm_total:.6f}, ||connected||={norm_connected:.6f}, "
          f"fraction={frac:.4f}")
    return connected, frac


def angular_decomposition(connected, scalar_basis, d_full):
    """Decompose the connected part into angular channels."""
    N_scalar = len(scalar_basis)
    dim_full = d_full

    # Extract (+,+) chirality sector
    V_pp = np.zeros((N_scalar, N_scalar, N_scalar, N_scalar))
    for a in range(N_scalar):
        for b in range(N_scalar):
            for c in range(N_scalar):
                for d in range(N_scalar):
                    row = a * dim_full + b
                    col = c * dim_full + d
                    V_pp[a, b, c, d] = connected[row, col]

    # Check m-conservation
    m_conserving = 0.0
    m_violating = 0.0
    for a, ba in enumerate(scalar_basis):
        for b, bb in enumerate(scalar_basis):
            for c, bc in enumerate(scalar_basis):
                for d, bd in enumerate(scalar_basis):
                    val = V_pp[a, b, c, d]
                    if abs(val) < 1e-14:
                        continue
                    if ba.m + bb.m == bc.m + bd.m:
                        m_conserving += val**2
                    else:
                        m_violating += val**2

    total_m = m_conserving + m_violating
    m_cons_pct = m_conserving / total_m * 100 if total_m > 0 else 0

    # Gaunt compatibility
    gaunt_compat = 0.0
    gaunt_incompat = 0.0
    for a, ba in enumerate(scalar_basis):
        for b, bb in enumerate(scalar_basis):
            for c, bc in enumerate(scalar_basis):
                for d, bd in enumerate(scalar_basis):
                    val = V_pp[a, b, c, d]
                    if abs(val) < 1e-14:
                        continue
                    dl1 = abs(ba.l - bc.l)
                    dl2 = abs(bb.l - bd.l)
                    if dl1 == dl2:
                        gaunt_compat += val**2
                    else:
                        gaunt_incompat += val**2

    total_g = gaunt_compat + gaunt_incompat
    gaunt_pct = gaunt_compat / total_g * 100 if total_g > 0 else 0

    # Multipole decomposition
    by_k = {}
    for a, ba in enumerate(scalar_basis):
        for b, bb in enumerate(scalar_basis):
            for c, bc in enumerate(scalar_basis):
                for d, bd in enumerate(scalar_basis):
                    val = V_pp[a, b, c, d]
                    if abs(val) < 1e-14:
                        continue
                    dl1 = abs(ba.l - bc.l)
                    dl2 = abs(bb.l - bd.l)
                    if dl1 == dl2:
                        by_k[dl1] = by_k.get(dl1, 0) + val**2

    total_k = sum(by_k.values())

    print(f"  m-conservation: {m_cons_pct:.1f}%")
    print(f"  Gaunt-compatible: {gaunt_pct:.1f}%")
    if total_k > 0:
        print(f"  Multipole decomposition:")
        for k in sorted(by_k.keys()):
            print(f"    k={k}: {by_k[k]/total_k*100:.1f}%")

    return m_cons_pct, gaunt_pct, by_k


def main():
    print("=" * 72)
    print("FULL-ALGEBRA GAUGE FIELD ON T_{S³} ⊗ T_{S³}")
    print("=" * 72)

    for n_max in [2, 3]:
        print(f"\n{'='*60}")
        print(f"n_max = {n_max}")
        print(f"{'='*60}")

        t0 = time.time()
        sp1 = build_single_particle_scalar(n_max)
        sp2 = build_single_particle_scalar(n_max)
        d1, d2 = sp1['dim'], sp2['dim']
        N_scalar = sp1['N_scalar']
        scalar_basis = sp1['scalar_basis']

        D1, D2 = sp1['D'], sp2['D']
        gamma1 = sp1['gamma']
        I1, I2 = np.eye(d1), np.eye(d2)
        D_total = np.kron(D1, I2) + np.kron(gamma1, D2)

        multipliers = build_all_multipliers(n_max, scalar_basis, N_scalar, d1)
        print(f"  {len(multipliers)} multiplier matrices in O_{{n_max={n_max}}}")

        # Count multipliers with nonzero [D, M]
        n_active = sum(1 for m in multipliers
                       if np.linalg.norm(D1 @ m['M_full'] - m['M_full'] @ D1) > 1e-14)
        print(f"  {n_active} have nonzero [D, M] (cross-shell couplings)")

        # Build the full cross gauge field
        print(f"\n  Building full-algebra cross gauge field...")
        A_cross, n_terms = build_cross_gauge_field(
            sp1, sp2, multipliers, D_total, gamma1
        )
        print(f"  {n_terms} nonzero (i,j) pairs summed")
        print(f"  ||A_cross|| = {np.linalg.norm(A_cross):.6f}")
        print(f"  Build time: {time.time()-t0:.1f} s")

        # {D, A} for the full-algebra gauge field
        DA = D_total @ A_cross + A_cross @ D_total
        print(f"  ||{{D, A_cross}}|| = {np.linalg.norm(DA):.6f}")

        # Factorization test
        print(f"\n  --- Factorization test ---")
        connected_DA, frac_DA = factorization_test(DA, d1, d2, "{D, A}")
        connected_Asq, frac_Asq = factorization_test(
            A_cross @ A_cross, d1, d2, "A²")

        # Angular decomposition of the connected part
        if np.linalg.norm(connected_DA) > 1e-14:
            print(f"\n  --- Angular decomposition of connected {{D, A}} ---")
            m_pct, g_pct, by_k = angular_decomposition(
                connected_DA, scalar_basis, d1)

        # Compare with single-multiplier result
        print(f"\n  --- Comparison: single-multiplier vs full-algebra ---")
        # Single-multiplier: use strongest commutator
        best_mult = max(multipliers,
                        key=lambda m: np.linalg.norm(D1 @ m['M_full'] - m['M_full'] @ D1))
        M_best = best_mult['M_full']
        comm_best = D1 @ M_best - M_best @ D1

        fluct_single = np.kron(comm_best, M_best) + np.kron(gamma1 @ M_best, comm_best)
        a_single = np.kron(M_best, M_best)
        A_single = a_single @ fluct_single
        A_single = (A_single + A_single.conj().T).real / 2

        DA_single = D_total @ A_single + A_single @ D_total
        conn_single, frac_single = factorization_test(DA_single, d1, d2, "{D,A} single")

        if np.linalg.norm(conn_single) > 1e-14:
            print(f"  Single-multiplier angular:")
            angular_decomposition(conn_single, scalar_basis, d1)

        # Heat trace comparison
        print(f"\n  --- Heat trace (full algebra) ---")
        t_test = 0.5
        ht_1 = np.sum(np.exp(-t_test * np.diag(D1 @ D1)))
        ht_2 = np.sum(np.exp(-t_test * np.diag(D2 @ D2)))
        ht_fact = ht_1 * ht_2

        for eps in [0.0, 0.1, 0.5, 1.0]:
            D_g = D_total + eps * A_cross
            evals_g = np.linalg.eigvalsh(D_g @ D_g)
            ht_g = np.sum(np.exp(-t_test * evals_g))
            r = ht_g / ht_fact
            print(f"    eps={eps:.1f}: ratio = {r:.8f}")

    print(f"\n{'='*72}")
    print("VERDICT")
    print(f"{'='*72}")


if __name__ == '__main__':
    main()
