"""
Compare the tensor-product interaction matrix elements with Coulomb 1/r₁₂
==========================================================================
The full-algebra gauge field gives a pure k=0 monopole, 100% Gaunt-compatible,
100% m-conserving interaction. Now compare the MATRIX ELEMENT VALUES with
the actual Coulomb two-electron integrals from Paper 12's algebraic engine.

The connected part V_connected has elements:
    V[a,b,c,d] = <n₁l₁m₁, n₂l₂m₂ | V_connected | n₁'l₁'m₁', n₂'l₂'m₂'>

The Coulomb integrals are:
    C[a,b,c,d] = <n₁l₁m₁, n₂l₂m₂ | 1/r₁₂ | n₁'l₁'m₁', n₂'l₂'m₂'>

If V = λ·C (proportional), the framework derives the Coulomb interaction
with a single overall scale factor λ (which IS the coupling constant).

If V ≠ λ·C but has the same angular structure, the framework derives
selection rules but not the radial kernel.
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
from geovac.casimir_ci import two_electron_integral
from debug.tensor_product_dirac import build_single_particle_scalar


def build_coulomb_matrix(scalar_basis):
    """Build the exact Coulomb 1/r₁₂ matrix in the Fock-projected basis.

    C[a,b,c,d] = <a,b|1/r₁₂|c,d> in the Dirac (physics) convention.
    """
    N = len(scalar_basis)
    C = np.zeros((N, N, N, N))

    for a, ba in enumerate(scalar_basis):
        for b, bb in enumerate(scalar_basis):
            for c, bc in enumerate(scalar_basis):
                for d, bd in enumerate(scalar_basis):
                    val = two_electron_integral(
                        ba.n, ba.l, ba.m,
                        bb.n, bb.l, bb.m,
                        bc.n, bc.l, bc.m,
                        bd.n, bd.l, bd.m,
                    )
                    if abs(val) > 1e-15:
                        C[a, b, c, d] = val
    return C


def build_connected_interaction_full_algebra(n_max):
    """Build the connected interaction from the full rotationally invariant gauge field."""
    sp1 = build_single_particle_scalar(n_max)
    sp2 = build_single_particle_scalar(n_max)
    d1, d2 = sp1['dim'], sp2['dim']
    N_scalar = sp1['N_scalar']
    scalar_basis = sp1['scalar_basis']

    D1, D2 = sp1['D'], sp2['D']
    gamma1 = sp1['gamma']
    I1, I2 = np.eye(d1), np.eye(d2)
    D_total = np.kron(D1, I2) + np.kron(gamma1, D2)

    # Build multipliers
    labels = allowed_multiplier_labels(n_max)
    multipliers = []
    for (N, L, M) in labels:
        M_sc = build_multiplier_matrix(N, L, M, scalar_basis)
        if np.linalg.norm(M_sc) < 1e-15:
            continue
        M_full = np.zeros((d1, d1), dtype=complex)
        M_full[:N_scalar, :N_scalar] = M_sc
        M_full[N_scalar:, N_scalar:] = M_sc
        multipliers.append(M_full)

    # Build rotationally invariant A_cross = Σ M†ᵢ⊗M†ⱼ · [D_total, Mᵢ⊗Mⱼ]
    A_cross = np.zeros((d1 * d2, d1 * d2), dtype=complex)
    for Mi in multipliers:
        Mi_dag = Mi.conj().T
        comm_i = D1 @ Mi - Mi @ D1
        if np.linalg.norm(comm_i) < 1e-14:
            continue
        for Mj in multipliers:
            Mj_dag = Mj.conj().T
            comm_j = D2 @ Mj - Mj @ D2
            if np.linalg.norm(comm_j) < 1e-14:
                continue
            fluct = np.kron(comm_i, Mj) + np.kron(gamma1 @ Mi, comm_j)
            a_dag = np.kron(Mi_dag, Mj_dag)
            A_cross += a_dag @ fluct

    A_cross = ((A_cross + A_cross.conj().T) / 2).real

    # {D, A}
    DA = D_total @ A_cross + A_cross @ D_total

    # Connected part
    dim_total = d1 * d2
    DA_tensor = DA.reshape(d1, d2, d1, d2)
    Tr2 = np.trace(DA_tensor, axis1=1, axis2=3)
    Tr1 = np.trace(DA_tensor, axis1=0, axis2=2)
    full_tr = np.trace(DA)
    factorizable = (np.kron(Tr2 / d2, I2)
                    + np.kron(I1, Tr1 / d1)
                    - (full_tr / (d1 * d2)) * np.eye(dim_total))
    connected = DA - factorizable

    # Extract the (+,+) chirality sector
    V_pp = np.zeros((N_scalar, N_scalar, N_scalar, N_scalar))
    for a in range(N_scalar):
        for b in range(N_scalar):
            for c in range(N_scalar):
                for d in range(N_scalar):
                    row = a * d1 + b
                    col = c * d1 + d
                    V_pp[a, b, c, d] = connected[row, col]

    return V_pp, scalar_basis


def compare_interactions(V_spectral, V_coulomb, scalar_basis, label=""):
    """Compare the spectral-action interaction with Coulomb."""
    N = len(scalar_basis)

    # Collect nonzero elements from both
    spectral_vals = []
    coulomb_vals = []
    labels_list = []

    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    vs = V_spectral[a, b, c, d]
                    vc = V_coulomb[a, b, c, d]
                    if abs(vs) > 1e-14 or abs(vc) > 1e-14:
                        spectral_vals.append(vs)
                        coulomb_vals.append(vc)
                        ba = scalar_basis[a]
                        bb = scalar_basis[b]
                        bc = scalar_basis[c]
                        bd = scalar_basis[d]
                        labels_list.append(
                            f"({ba.n}{ba.l}{ba.m:+d},{bb.n}{bb.l}{bb.m:+d}|"
                            f"{bc.n}{bc.l}{bc.m:+d},{bd.n}{bd.l}{bd.m:+d})"
                        )

    spectral_vals = np.array(spectral_vals)
    coulomb_vals = np.array(coulomb_vals)

    # Overall correlation
    mask_both = (np.abs(spectral_vals) > 1e-14) & (np.abs(coulomb_vals) > 1e-14)
    n_both = np.sum(mask_both)
    n_spectral_only = np.sum((np.abs(spectral_vals) > 1e-14) & (np.abs(coulomb_vals) <= 1e-14))
    n_coulomb_only = np.sum((np.abs(spectral_vals) <= 1e-14) & (np.abs(coulomb_vals) > 1e-14))

    print(f"\n  {label} — element comparison:")
    print(f"    Nonzero in both: {n_both}")
    print(f"    Spectral only: {n_spectral_only}")
    print(f"    Coulomb only: {n_coulomb_only}")

    if n_both == 0:
        print(f"    No overlapping nonzero elements — cannot compare radial weights")
        return None

    # For the overlapping elements, compute the ratio V_spectral / V_coulomb
    ratios = spectral_vals[mask_both] / coulomb_vals[mask_both]
    mean_ratio = np.mean(ratios)
    std_ratio = np.std(ratios)
    cv = std_ratio / abs(mean_ratio) if abs(mean_ratio) > 1e-14 else float('inf')

    print(f"    Mean ratio V_spectral/V_coulomb = {mean_ratio:.6e}")
    print(f"    Std of ratio = {std_ratio:.6e}")
    print(f"    CV (std/mean) = {cv:.4f}")

    if cv < 0.01:
        print(f"    ✓ PROPORTIONAL — V_spectral = {mean_ratio:.6e} × V_coulomb (CV < 1%)")
        print(f"      The framework derives the Coulomb interaction with a single scale factor!")
    elif cv < 0.1:
        print(f"    ~ APPROXIMATELY proportional (CV < 10%) — same structure, small deviations")
    else:
        print(f"    ✗ NOT proportional (CV = {cv:.1%}) — radial weights differ")

    # Show the largest matrix elements side by side
    combined = list(zip(labels_list, spectral_vals, coulomb_vals))
    combined_both = [(lab, vs, vc) for lab, vs, vc in combined
                     if abs(vs) > 1e-14 and abs(vc) > 1e-14]
    combined_both.sort(key=lambda x: -abs(x[2]))

    print(f"\n    Top matrix elements (sorted by |V_coulomb|):")
    print(f"    {'element':>40} {'V_spectral':>12} {'V_coulomb':>12} {'ratio':>10}")
    print(f"    " + "-" * 78)

    for lab, vs, vc in combined_both[:15]:
        r = vs / vc if abs(vc) > 1e-14 else float('nan')
        print(f"    {lab:>40} {vs:>+11.6e} {vc:>+11.6e} {r:>+9.4f}")

    # Pearson correlation
    if n_both >= 3:
        corr = np.corrcoef(spectral_vals[mask_both], coulomb_vals[mask_both])[0, 1]
        print(f"\n    Pearson correlation: {corr:.6f}")
        if abs(corr) > 0.99:
            print(f"    ✓ Very strong linear correlation (|r| > 0.99)")
        elif abs(corr) > 0.9:
            print(f"    ~ Strong linear correlation (|r| > 0.9)")

    # Sign agreement
    sign_agree = np.sum(np.sign(spectral_vals[mask_both]) == np.sign(coulomb_vals[mask_both]))
    print(f"    Sign agreement: {sign_agree}/{n_both} ({sign_agree/n_both*100:.0f}%)")

    return {
        'n_both': int(n_both),
        'n_spectral_only': int(n_spectral_only),
        'n_coulomb_only': int(n_coulomb_only),
        'mean_ratio': float(mean_ratio),
        'std_ratio': float(std_ratio),
        'cv': float(cv),
        'pearson': float(corr) if n_both >= 3 else None,
    }


def main():
    print("=" * 72)
    print("RADIAL COMPARISON: SPECTRAL-ACTION V vs COULOMB 1/r₁₂")
    print("=" * 72)

    for n_max in [2, 3]:
        print(f"\n{'='*60}")
        print(f"n_max = {n_max}")
        print(f"{'='*60}")

        t0 = time.time()

        # Build the Coulomb matrix
        sp = build_single_particle_scalar(n_max)
        scalar_basis = sp['scalar_basis']
        print(f"  Building Coulomb matrix ({len(scalar_basis)}⁴ = "
              f"{len(scalar_basis)**4} elements)...")
        V_coulomb = build_coulomb_matrix(scalar_basis)
        n_nonzero_C = np.sum(np.abs(V_coulomb) > 1e-14)
        print(f"  {n_nonzero_C} nonzero Coulomb elements")

        # Build the connected spectral-action interaction
        print(f"  Building spectral-action interaction...")
        V_spectral, _ = build_connected_interaction_full_algebra(n_max)
        n_nonzero_S = np.sum(np.abs(V_spectral) > 1e-14)
        print(f"  {n_nonzero_S} nonzero spectral elements")
        print(f"  Build time: {time.time()-t0:.1f} s")

        # Compare
        result = compare_interactions(V_spectral, V_coulomb, scalar_basis,
                                       f"n_max={n_max}")

    print(f"\n{'='*72}")
    print("INTERPRETATION")
    print(f"{'='*72}")
    print("""
  If V_spectral = λ × V_coulomb (proportional with single scale factor):
    → The framework DERIVES the Coulomb potential from spectral geometry.
    → λ is the coupling constant (related to α or e²).
    → This would be the single biggest structural advance since WH1.

  If V_spectral has same angular structure but different radial weights:
    → The framework derives WHICH channels couple (selection rules)
    → but NOT how strongly (radial kernel is calibration data).
    → Still significant: Gaunt rules emerge from the algebra, not imported.

  If V_spectral has different structure entirely:
    → The gauge field generates a DIFFERENT interaction than Coulomb.
    → The tensor-product spectral action is genuine physics, just not V_ee.
""")


if __name__ == '__main__':
    main()
