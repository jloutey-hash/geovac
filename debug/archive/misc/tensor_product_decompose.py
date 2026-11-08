"""
Decompose the connected two-body interaction into angular channels
=================================================================
Step 3: The connected part of {D, A} from Step 2 is a genuine two-body
operator on H₁ ⊗ H₂. Decompose it into the (l₁, l₂) → (l₁', l₂')
angular channel structure and compare with:

1. The Gaunt selection rules (Paper 22): which channels are nonzero?
2. The Neumann expansion of 1/r₁₂ (Paper 12): angular × radial form

The Neumann expansion:
  1/r₁₂ = Σ_k (r_<^k / r_>^{k+1}) * (4π/(2k+1)) * Σ_q Y*_kq(Ω₁) Y_kq(Ω₂)

In the S³ angular basis, this becomes:
  V_ee ~ Σ_k C_k * Σ_q <n₁l₁m₁|Y_kq|n₁'l₁'m₁'> <n₂l₂m₂|Y*_kq|n₂'l₂'m₂'>

The ANGULAR part is exactly the 3-Y integrals (Gaunt coefficients).
The RADIAL part C_k involves the Legendre Q functions (imported physics in Paper 12).

If the connected part of {D, A} has the SAME angular channel structure,
the framework derives the selection rules natively. If the radial
weights also match, the framework derives 1/r₁₂ itself.
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


def build_connected_interaction(sp1, sp2, n_max):
    """Build the connected (non-factorizable) part of {D, A}."""
    d1, d2 = sp1['dim'], sp2['dim']
    dim_total = d1 * d2
    N_scalar = sp1['N_scalar']
    D1, D2 = sp1['D'], sp2['D']
    gamma1 = sp1['gamma']
    I1, I2 = np.eye(d1), np.eye(d2)
    scalar_basis = sp1['scalar_basis']

    D_total = np.kron(D1, I2) + np.kron(gamma1, D2)

    # Build gauge field from strongest commutator
    labels = allowed_multiplier_labels(n_max)
    best_comm = None
    best_norm = 0
    best_M = None

    for (N, L, M) in labels:
        M_sc = build_multiplier_matrix(N, L, M, scalar_basis).real
        if np.linalg.norm(M_sc) < 1e-15:
            continue
        M_full = np.zeros((d1, d1))
        M_full[:N_scalar, :N_scalar] = M_sc
        M_full[N_scalar:, N_scalar:] = M_sc
        comm = D1 @ M_full - M_full @ D1
        n = np.linalg.norm(comm)
        if n > best_norm:
            best_norm = n
            best_comm = comm
            best_M = M_full

    # Build A_total
    fluct1 = np.kron(best_comm, best_M)
    fluct2 = np.kron(gamma1 @ best_M, best_comm)
    fluct = fluct1 + fluct2
    a_tensor = np.kron(best_M, best_M)
    A_total = a_tensor @ fluct
    A_total = (A_total + A_total.T) / 2

    # {D, A}
    DA = D_total @ A_total + A_total @ D_total

    # Connected part
    DA_tensor = DA.reshape(d1, d2, d1, d2)
    Tr2 = np.trace(DA_tensor, axis1=1, axis2=3)
    Tr1 = np.trace(DA_tensor, axis1=0, axis2=2)
    full_tr = np.trace(DA)

    factorizable = (np.kron(Tr2 / d2, I2)
                    + np.kron(I1, Tr1 / d1)
                    - (full_tr / (d1 * d2)) * np.eye(dim_total))
    connected = DA - factorizable

    return connected, DA, scalar_basis


def decompose_angular_channels(connected, scalar_basis, n_max):
    """Decompose the connected interaction into angular channels.

    The connected operator lives on (chirality-doubled H₁) ⊗ (chirality-doubled H₂).
    We extract the SCALAR SECTOR (positive chirality on both particles)
    and decompose into (n₁,l₁,m₁; n₂,l₂,m₂) → (n₁',l₁',m₁'; n₂',l₂',m₂') blocks.

    The key question: what (Δl₁, Δl₂) channels have nonzero matrix elements?
    For 1/r₁₂, the Gaunt selection rules give: Δl₁ = k, Δl₂ = k for each
    multipole order k, with the triangle inequality on the 3j symbols.
    """
    N_scalar = len(scalar_basis)
    dim_full = 2 * N_scalar  # chirality-doubled

    # Extract the (+,+) chirality sector: both particles in positive chirality
    # In our basis ordering: first N_scalar states are + chirality
    # Tensor product index: (i₁ * d2 + i₂) where i₁, i₂ ∈ [0, dim_full)
    # The (+,+) sector has i₁ ∈ [0, N_scalar) and i₂ ∈ [0, N_scalar)

    d = dim_full
    V_pp = np.zeros((N_scalar, N_scalar, N_scalar, N_scalar))
    for a in range(N_scalar):
        for b in range(N_scalar):
            for c in range(N_scalar):
                for d_idx in range(N_scalar):
                    row = a * dim_full + b     # (+,+) sector
                    col = c * dim_full + d_idx
                    V_pp[a, b, c, d_idx] = connected[row, col]

    # Now decompose V_pp into angular channels
    # V_pp[a, b, c, d] = <n₁l₁m₁, n₂l₂m₂ | V | n₁'l₁'m₁', n₂'l₂'m₂'>
    # where a = (n₁, l₁, m₁) and b = (n₂, l₂, m₂), etc.

    # Group by (l₁, l₂) → (l₁', l₂') channel
    channels = {}
    for a, ba in enumerate(scalar_basis):
        for b, bb in enumerate(scalar_basis):
            for c, bc in enumerate(scalar_basis):
                for d_idx, bd in enumerate(scalar_basis):
                    val = V_pp[a, b, c, d_idx]
                    if abs(val) < 1e-14:
                        continue
                    key = (ba.l, bb.l, bc.l, bd.l)
                    channels[key] = channels.get(key, 0) + val**2

    # Sort by magnitude
    sorted_channels = sorted(channels.items(), key=lambda x: -x[1])

    print(f"\n  Angular channel decomposition of connected {'{'}D, A{'}'} (++ chirality sector):")
    print(f"  {'(l1,l2)→(l1p,l2p)':>25} {'||V||²':>12} {'fraction':>10}")
    print(f"  " + "-" * 50)

    total_norm_sq = sum(v for _, v in sorted_channels)
    for (l1, l2, l1p, l2p), norm_sq in sorted_channels[:15]:
        frac = norm_sq / total_norm_sq if total_norm_sq > 0 else 0
        print(f"  ({l1},{l2})→({l1p},{l2p}){' ':>15} {norm_sq:>12.6e} {frac:>9.4f}")

    # Check Gaunt selection rules
    print(f"\n  Gaunt selection rule check:")
    print(f"  For 1/r₁₂ Neumann expansion, allowed channels have:")
    print(f"    |l₁ - l₁'| = |l₂ - l₂'| = k for some multipole k")
    print(f"    Triangle: |l₁ - k| ≤ l₁' ≤ l₁ + k, same for l₂")

    gaunt_compatible = 0
    gaunt_incompatible = 0
    for (l1, l2, l1p, l2p), norm_sq in sorted_channels:
        if norm_sq < 1e-20:
            continue
        dl1 = abs(l1 - l1p)
        dl2 = abs(l2 - l2p)
        if dl1 == dl2:
            gaunt_compatible += norm_sq
        else:
            gaunt_incompatible += norm_sq

    total = gaunt_compatible + gaunt_incompatible
    if total > 0:
        print(f"  Gaunt-compatible: {gaunt_compatible/total*100:.1f}%")
        print(f"  Gaunt-incompatible: {gaunt_incompatible/total*100:.1f}%")
    else:
        print(f"  No nonzero channels found")

    # m-selection rule check
    print(f"\n  m-conservation check (m₁ + m₂ = m₁' + m₂'):")
    m_conserving = 0
    m_violating = 0
    for a, ba in enumerate(scalar_basis):
        for b, bb in enumerate(scalar_basis):
            for c, bc in enumerate(scalar_basis):
                for d_idx, bd in enumerate(scalar_basis):
                    val = V_pp[a, b, c, d_idx]
                    if abs(val) < 1e-14:
                        continue
                    if ba.m + bb.m == bc.m + bd.m:
                        m_conserving += val**2
                    else:
                        m_violating += val**2
    total_m = m_conserving + m_violating
    if total_m > 0:
        print(f"  m-conserving: {m_conserving/total_m*100:.1f}%")
        print(f"  m-violating: {m_violating/total_m*100:.1f}%")

    return V_pp, channels


def compare_with_coulomb(V_pp, scalar_basis, n_max):
    """Compare the extracted interaction matrix elements with 1/r₁₂.

    At n_max=2, the scalar basis has states:
    |1,0,0>, |2,0,0>, |2,1,-1>, |2,1,0>, |2,1,1>

    The dominant Coulomb matrix elements at the lowest level are:
    <1s,1s|1/r₁₂|1s,1s> (direct Coulomb integral)

    For S³ at unit radius, the "Coulomb" integral involves the
    chordal distance, which in the Fock projection gives the
    standard Slater integrals F^k and exchange integrals G^k.

    At n_max=2, the nonzero angular channels for 1/r₁₂ are:
    k=0: (0,0)→(0,0) [monopole]
    k=1: (0,1)→(1,0) and permutations [dipole]
    """
    N = len(scalar_basis)

    # The "Coulomb-like" structure has specific symmetries:
    # 1. Hermitian: V[a,b,c,d] = V[c,d,a,b]*
    # 2. Exchange: V[a,b,c,d] = V[b,a,d,c] (for identical particles)
    # 3. m-conservation: m_a + m_b = m_c + m_d

    # Check symmetries
    hermitian_err = np.linalg.norm(V_pp - V_pp.transpose(2, 3, 0, 1))
    exchange_err = np.linalg.norm(V_pp - V_pp.transpose(1, 0, 3, 2))

    print(f"\n  Symmetry check of the connected interaction:")
    print(f"  Hermiticity error: {hermitian_err:.2e}")
    print(f"  Exchange symmetry error: {exchange_err:.2e}")

    # Extract the dominant matrix elements and their quantum numbers
    print(f"\n  Largest matrix elements (top 10):")
    print(f"  {'(n1,l1,m1)':>12} {'(n2,l2,m2)':>12} {'(n1p,l1p,m1p)':>14} {'(n2p,l2p,m2p)':>14} {'V':>12}")
    print(f"  " + "-" * 70)

    elements = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    val = V_pp[a, b, c, d]
                    if abs(val) > 1e-14:
                        elements.append((a, b, c, d, val))

    elements.sort(key=lambda x: -abs(x[4]))

    for a, b, c, d, val in elements[:10]:
        ba, bb, bc, bd = scalar_basis[a], scalar_basis[b], scalar_basis[c], scalar_basis[d]
        print(f"  ({ba.n},{ba.l},{ba.m:+d}) ({bb.n},{bb.l},{bb.m:+d}) "
              f"({bc.n},{bc.l},{bc.m:+d}) ({bd.n},{bd.l},{bd.m:+d})  {val:>+11.6e}")

    # Compare structure with known Slater integrals
    # The k=0 (monopole) channel: (l₁, l₂) → (l₁, l₂), Δm = 0
    # Corresponds to the F⁰ Slater integral (direct Coulomb)
    # The k=1 (dipole) channel: |Δl₁| = |Δl₂| = 1
    # Corresponds to the F¹ / G¹ Slater integrals

    # Count elements by multipole order k
    by_k = {}
    for a, b, c, d, val in elements:
        ba, bb, bc, bd = scalar_basis[a], scalar_basis[b], scalar_basis[c], scalar_basis[d]
        dl1 = abs(ba.l - bc.l)
        dl2 = abs(bb.l - bd.l)
        if dl1 == dl2:
            k = dl1
            by_k[k] = by_k.get(k, 0) + val**2

    print(f"\n  Multipole decomposition (assuming Δl₁ = Δl₂ = k):")
    total_by_k = sum(by_k.values())
    for k in sorted(by_k.keys()):
        frac = by_k[k] / total_by_k * 100 if total_by_k > 0 else 0
        print(f"    k={k}: {frac:.1f}% of ||V||²")

    return elements


def main():
    print("=" * 72)
    print("ANGULAR DECOMPOSITION OF THE CONNECTED INTERACTION")
    print("=" * 72)

    for n_max in [2, 3]:
        print(f"\n{'='*60}")
        print(f"n_max = {n_max}")
        print(f"{'='*60}")

        sp1 = build_single_particle_scalar(n_max)
        sp2 = build_single_particle_scalar(n_max)

        connected, DA, scalar_basis = build_connected_interaction(sp1, sp2, n_max)

        print(f"  ||connected|| = {np.linalg.norm(connected):.6f}")

        V_pp, channels = decompose_angular_channels(connected, scalar_basis, n_max)
        elements = compare_with_coulomb(V_pp, scalar_basis, n_max)

    print(f"\n{'='*72}")
    print("VERDICT")
    print(f"{'='*72}")
    print("""
  If the angular channels match the Gaunt selection rules (Δl₁ = Δl₂ = k),
  the framework's tensor-product gauged spectral action reproduces the
  ANGULAR structure of the two-body Coulomb interaction natively.

  The RADIAL weights (relative strength of k=0 vs k=1 vs k=2 multipoles)
  would then need to be compared with the Slater F^k / G^k integrals from
  Paper 12's Neumann expansion. If they match, the framework derives V_ee
  from spectral geometry. If only the selection rules match (but not the
  radial weights), the framework derives WHICH channels couple but not
  HOW STRONGLY — the radial coupling remains calibration data.
""")


if __name__ == '__main__':
    main()
