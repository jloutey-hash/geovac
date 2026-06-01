"""
Resolvent-based two-body interaction from the spectral triple
=============================================================
Paper 54 showed: the spectral action (heat-kernel weighting) gives the
right angular selection rules but wrong radial coupling strengths.

The open question: does the RESOLVENT (D² - z)⁻¹ do better?

The Coulomb potential IS the Green's function of the Laplacian:
  ∇²G = -4πδ  =>  G ~ 1/r₁₂

On S³, the Laplacian has eigenvalues n²-1 (n=1,2,...), so:
  G_{S³}(Ω₁, Ω₂) = Σ_{N≥2} Σ_{LM} Y*_{NLM}(Ω₁) Y_{NLM}(Ω₂) / (N²-1)

The two-body matrix elements are:
  V^{res}_{abcd} = Σ_{NLM} M_{NLM}[a,c] × M_{NLM}[b,d] × w(N)

where M_{NLM} are the 3-Y multiplier matrices and w(N) is the weight.

Three candidate weightings:
  1. Laplacian resolvent:  w(N) = 1/(N²-1)     [Green's function of Δ on S³]
  2. Dirac resolvent:      w(N) = 1/(N+1/2)²   [Green's function of D²]
  3. Spectral action:      w(N) from {D, A}     [Paper 54's construction]

Compare all three to exact Coulomb Slater integrals.
"""
import sys
import json
import time
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.operator_system import (
    build_multiplier_matrix, HyperLabel, allowed_multiplier_labels,
)
from geovac.hypergeometric_slater import get_rk_float
from sympy.physics.wigner import wigner_3j
import sympy


def build_scalar_basis(n_max):
    """Build the scalar Fock-projected basis |n, l, m> up to n_max."""
    basis = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                basis.append(HyperLabel(n=n, l=l, m=m))
    return basis


def compute_resolvent_interaction(n_max, weight_fn, label=""):
    """Compute the two-body interaction V_{abcd} from the resolvent construction.

    The Green's function on S³ expands as:
      G(Ω₁, Ω₂) = Σ_{NLM} Y_{NLM}(Ω₁) Y*_{NLM}(Ω₂) / λ_N

    The two-body matrix element is:
      V_{abcd} = ∫∫ Y*_a(1) Y*_b(2) G(1,2) Y_c(1) Y_d(2) dΩ₁ dΩ₂
               = Σ_{NLM} M_{NLM}[a,c] × conj(M_{NLM}[d,b]) × w(N)

    where M_{NLM}[i,j] = ∫ Y*_i Y_{NLM} Y_j dΩ (the 3-Y integral).
    The conjugation pattern conj(M[d,b]) = M†[b,d] enforces m-conservation.
    """
    basis = build_scalar_basis(n_max)
    N_dim = len(basis)
    labels = allowed_multiplier_labels(n_max)

    V = np.zeros((N_dim, N_dim, N_dim, N_dim), dtype=np.complex128)

    n_multipliers = 0
    for (N, L, M) in labels:
        w = weight_fn(N)
        if abs(w) < 1e-30:
            continue

        M_mat = build_multiplier_matrix(N, L, M, basis)

        if np.linalg.norm(M_mat) < 1e-15:
            continue

        n_multipliers += 1

        # V_{abcd} = M[a,c] × conj(M[d,b]) × w
        # conj(M[d,b]) = M†[b,d] = (M.conj().T)[b,d]
        M_dag = M_mat.conj().T

        for a in range(N_dim):
            for c in range(N_dim):
                mac = M_mat[a, c]
                if abs(mac) < 1e-15:
                    continue
                for b in range(N_dim):
                    for d in range(N_dim):
                        m_dag_bd = M_dag[b, d]
                        if abs(m_dag_bd) < 1e-15:
                            continue
                        V[a, b, c, d] += mac * m_dag_bd * w

    V_real = V.real
    imag_norm = np.linalg.norm(V.imag)
    real_norm = np.linalg.norm(V_real)
    print(f"  [{label}] n_max={n_max}: {n_multipliers} nonzero multipliers, "
          f"dim={N_dim}, V nonzero = {np.count_nonzero(np.abs(V_real) > 1e-14)}, "
          f"imag/real = {imag_norm/real_norm:.2e}")
    return V_real, basis


def compute_coulomb_slater(n_max):
    """Compute exact Coulomb two-electron integrals in the Fock basis.

    The two-electron integral is:
      <ab|1/r₁₂|cd> = Σ_k X_k(a,c) × X_k(b,d) × R^k(a,c;b,d)

    where X_k are angular (Gaunt) coefficients and R^k are radial Slater integrals.

    For the S³ Fock-projected basis, the angular part is:
      X_k(a,c) = (-1)^{m_a} C(l_a, k, l_c; -m_a, m_a-m_c, m_c) × <l_a || C^k || l_c>
    where the reduced matrix element involves a 3j symbol.
    """
    basis = build_scalar_basis(n_max)
    N_dim = len(basis)
    V = np.zeros((N_dim, N_dim, N_dim, N_dim))

    for a, ba in enumerate(basis):
        for b, bb in enumerate(basis):
            for c, bc in enumerate(basis):
                for d, bd in enumerate(basis):
                    # m-conservation
                    if ba.m + bb.m != bc.m + bd.m:
                        continue

                    val = _two_electron_integral(ba, bb, bc, bd)
                    if abs(val) > 1e-15:
                        V[a, b, c, d] = val

    return V, basis


def _two_electron_integral(a, b, c, d):
    """Compute <ab|1/r₁₂|cd> using Slater-Condon rules.

    In the angular momentum basis, the two-electron integral factorizes:
      <ab|V|cd> = Σ_k c_k(a,c) × c_k(b,d) × R^k(n_a l_a, n_c l_c; n_b l_b, n_d l_d)

    where c_k(a,c) = (-1)^{m_a} × √((2l_a+1)(2l_c+1)) × 3j(l_a,k,l_c;0,0,0) × 3j(l_a,k,l_c;-m_a,m_a-m_c,m_c)
    """
    val = 0.0

    # Determine allowed k range from Gaunt selection rules
    k_min_ac = abs(a.l - c.l)
    k_max_ac = a.l + c.l
    k_min_bd = abs(b.l - d.l)
    k_max_bd = b.l + d.l

    k_min = max(k_min_ac, k_min_bd)
    k_max = min(k_max_ac, k_max_bd)

    for k in range(k_min, k_max + 1):
        # Parity selection: l_a + k + l_c must be even
        if (a.l + k + c.l) % 2 != 0:
            continue
        if (b.l + k + d.l) % 2 != 0:
            continue

        # Angular coefficient for (a,c) pair
        ck_ac = _gaunt_ck(a.l, a.m, c.l, c.m, k)
        if abs(ck_ac) < 1e-15:
            continue

        # Angular coefficient for (b,d) pair
        ck_bd = _gaunt_ck(b.l, b.m, d.l, d.m, k)
        if abs(ck_bd) < 1e-15:
            continue

        # Radial Slater integral R^k
        try:
            rk = get_rk_float(a.n, a.l, c.n, c.l, b.n, b.l, d.n, d.l, k)
        except (ValueError, ZeroDivisionError):
            continue

        val += ck_ac * ck_bd * rk

    return val


def _gaunt_ck(la, ma, lc, mc, k):
    """Gaunt angular coefficient c_k(a,c).

    c_k = (-1)^{m_a} × √((2l_a+1)(2l_c+1)) × 3j(l_a,k,l_c;0,0,0) × 3j(l_a,k,l_c;-m_a,q,m_c)
    where q = m_a - m_c.
    """
    q = ma - mc
    if abs(q) > k:
        return 0.0

    # 3j symbols
    j1 = float(wigner_3j(la, k, lc, 0, 0, 0))
    if abs(j1) < 1e-15:
        return 0.0

    j2 = float(wigner_3j(la, k, lc, -ma, q, mc))
    if abs(j2) < 1e-15:
        return 0.0

    prefactor = (-1)**ma * np.sqrt((2*la + 1) * (2*lc + 1))
    return prefactor * j1 * j2


def compare_interactions(V_test, V_coulomb, basis, label):
    """Compare a candidate interaction to exact Coulomb."""
    N = len(basis)

    # Collect nonzero pairs
    test_vals = []
    coul_vals = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    vt = V_test[a, b, c, d]
                    vc = V_coulomb[a, b, c, d]
                    if abs(vc) > 1e-14:
                        test_vals.append(vt)
                        coul_vals.append(vc)

    test_vals = np.array(test_vals)
    coul_vals = np.array(coul_vals)

    if len(coul_vals) == 0:
        print(f"  [{label}] No nonzero Coulomb elements found!")
        return {}

    # Pearson correlation
    if np.std(test_vals) > 1e-15 and np.std(coul_vals) > 1e-15:
        pearson = np.corrcoef(test_vals, coul_vals)[0, 1]
    else:
        pearson = 0.0

    # Proportionality check: if V_test ∝ V_coulomb, then V_test/V_coulomb = const
    mask = np.abs(coul_vals) > 1e-10
    if np.sum(mask) > 1:
        ratios = test_vals[mask] / coul_vals[mask]
        cv_ratio = np.std(ratios) / np.abs(np.mean(ratios)) if np.abs(np.mean(ratios)) > 1e-15 else float('inf')
        mean_ratio = np.mean(ratios)
    else:
        cv_ratio = float('inf')
        mean_ratio = 0.0

    # Sign agreement
    sign_agree = np.mean(np.sign(test_vals[mask]) == np.sign(coul_vals[mask])) if np.sum(mask) > 0 else 0

    # How many Coulomb-nonzero elements does the candidate also have nonzero?
    coverage = np.mean(np.abs(test_vals) > 1e-14)

    # m-conservation check
    m_conserving = 0
    m_total = 0
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    if abs(V_test[a, b, c, d]) > 1e-14:
                        m_total += 1
                        if basis[a].m + basis[b].m == basis[c].m + basis[d].m:
                            m_conserving += 1

    result = {
        'label': label,
        'n_coulomb_nonzero': int(np.sum(mask)),
        'n_test_nonzero': int(np.sum(np.abs(test_vals) > 1e-14)),
        'pearson': float(pearson),
        'cv_ratio': float(cv_ratio),
        'mean_ratio': float(mean_ratio),
        'sign_agreement': float(sign_agree),
        'coverage': float(coverage),
        'm_conservation': float(m_conserving / m_total) if m_total > 0 else 0,
    }

    print(f"\n  === {label} ===")
    print(f"  Coulomb nonzero elements: {result['n_coulomb_nonzero']}")
    print(f"  Test nonzero elements:    {result['n_test_nonzero']}")
    print(f"  Pearson correlation:      {result['pearson']:.6f}")
    print(f"  CV of ratio V_test/V_coul:{result['cv_ratio']:.4f}")
    print(f"  Mean ratio:               {result['mean_ratio']:.6f}")
    print(f"  Sign agreement:           {result['sign_agreement']:.1%}")
    print(f"  Coverage:                 {result['coverage']:.1%}")
    print(f"  m-conservation:           {result['m_conservation']:.1%}")

    return result


def main():
    print("=" * 72)
    print("RESOLVENT TWO-BODY DIAGNOSTIC")
    print("Paper 54 open question: does (D²)⁻¹ recover 1/r₁₂?")
    print("=" * 72)

    results = {}

    for n_max in [2, 3]:
        print(f"\n{'=' * 72}")
        print(f"n_max = {n_max}")
        print(f"{'=' * 72}")

        # 1. Exact Coulomb
        print(f"\n  Computing exact Coulomb Slater integrals...")
        t0 = time.time()
        V_coulomb, basis = compute_coulomb_slater(n_max)
        n_coul = np.count_nonzero(np.abs(V_coulomb) > 1e-14)
        print(f"  Done in {time.time()-t0:.1f}s. {n_coul} nonzero elements.")

        # 2. Laplacian resolvent: w(N) = 1/(N² - 1), skip N=1
        print(f"\n  Computing Laplacian resolvent interaction...")
        t0 = time.time()
        V_lap, _ = compute_resolvent_interaction(
            n_max,
            weight_fn=lambda N: 1.0 / (N**2 - 1) if N >= 2 else 0.0,
            label="Laplacian"
        )
        print(f"  Done in {time.time()-t0:.1f}s.")

        # 3. Dirac resolvent: w(N) = 1/(N + 1/2)²
        print(f"\n  Computing Dirac resolvent interaction...")
        t0 = time.time()
        V_dirac, _ = compute_resolvent_interaction(
            n_max,
            weight_fn=lambda N: 1.0 / (N + 0.5)**2,
            label="Dirac"
        )
        print(f"  Done in {time.time()-t0:.1f}s.")

        # 4. Uniform weight (no radial weighting): w(N) = 1
        print(f"\n  Computing uniform-weight interaction (angular-only)...")
        t0 = time.time()
        V_uniform, _ = compute_resolvent_interaction(
            n_max,
            weight_fn=lambda N: 1.0,
            label="Uniform"
        )
        print(f"  Done in {time.time()-t0:.1f}s.")

        # 5. Inverse-linear: w(N) = 1/N  (between Laplacian and Dirac)
        print(f"\n  Computing 1/N weight interaction...")
        t0 = time.time()
        V_invN, _ = compute_resolvent_interaction(
            n_max,
            weight_fn=lambda N: 1.0 / N,
            label="1/N"
        )
        print(f"  Done in {time.time()-t0:.1f}s.")

        # Compare all to Coulomb
        print(f"\n  --- Comparison with exact Coulomb ---")
        r = {}
        r['laplacian'] = compare_interactions(V_lap, V_coulomb, basis, f"Laplacian 1/(N²-1), n_max={n_max}")
        r['dirac'] = compare_interactions(V_dirac, V_coulomb, basis, f"Dirac 1/(N+½)², n_max={n_max}")
        r['uniform'] = compare_interactions(V_uniform, V_coulomb, basis, f"Uniform w=1, n_max={n_max}")
        r['invN'] = compare_interactions(V_invN, V_coulomb, basis, f"1/N weight, n_max={n_max}")

        results[n_max] = r

        # Print the top-10 Coulomb elements and their resolvent counterparts
        print(f"\n  --- Top-10 Coulomb matrix elements vs resolvent ---")
        elements = []
        N_dim = len(basis)
        for a in range(N_dim):
            for b in range(N_dim):
                for c in range(N_dim):
                    for d in range(N_dim):
                        vc = V_coulomb[a, b, c, d]
                        if abs(vc) > 1e-14:
                            elements.append((a, b, c, d, vc))

        elements.sort(key=lambda x: -abs(x[4]))

        hdr = f"  {'(a,b,c,d)':>20} {'Coulomb':>10} {'Lap-res':>10} {'Dirac-res':>10} {'Uniform':>10}"
        print(hdr)
        print("  " + "-" * len(hdr))
        for a, b, c, d, vc in elements[:10]:
            ba, bb, bc, bd = basis[a], basis[b], basis[c], basis[d]
            lbl = f"({ba.n}{ba.l},{bb.n}{bb.l},{bc.n}{bc.l},{bd.n}{bd.l})"
            vl = V_lap[a, b, c, d]
            vd = V_dirac[a, b, c, d]
            vu = V_uniform[a, b, c, d]
            print(f"  {lbl:>20} {vc:>10.6f} {vl:>10.6f} {vd:>10.6f} {vu:>10.6f}")

    # Summary
    print(f"\n{'=' * 72}")
    print("SUMMARY")
    print(f"{'=' * 72}")
    print(f"\n  {'Method':>25} {'n_max=2 Pearson':>16} {'n_max=3 Pearson':>16} {'n_max=2 CV':>12} {'n_max=3 CV':>12}")
    print("  " + "-" * 85)
    for method in ['laplacian', 'dirac', 'uniform', 'invN']:
        r2 = results[2].get(method, {})
        r3 = results[3].get(method, {})
        p2 = r2.get('pearson', 0)
        p3 = r3.get('pearson', 0)
        c2 = r2.get('cv_ratio', float('inf'))
        c3 = r3.get('cv_ratio', float('inf'))
        print(f"  {method:>25} {p2:>16.6f} {p3:>16.6f} {c2:>12.4f} {c3:>12.4f}")

    # Save results
    out = {str(k): {mk: mv for mk, mv in v.items()} for k, v in results.items()}
    with open('debug/data/resolvent_two_body_diagnostic.json', 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n  Results saved to debug/data/resolvent_two_body_diagnostic.json")


if __name__ == '__main__':
    main()
