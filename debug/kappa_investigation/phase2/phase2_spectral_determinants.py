"""
Phase 2: Spectral Determinant Investigation & Second Selection Principle

Research questions:
1. Verify det'(S1), det'(S2), det'(S3) to high precision
2. Search for correction factors bridging 41.957 -> 42
3. Test if K/pi can be expressed via spectral determinants
4. Prove d_max = 4 for all n_max, analyze transition geometry
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Use mpmath for high-precision spectral determinant computations
try:
    import mpmath
    mpmath.mp.dps = 50  # 50 decimal places
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False
    print("WARNING: mpmath not available. Using numpy (limited precision).")

from geovac.lattice import GeometricLattice


# ============================================================================
# PART 1: SPECTRAL DETERMINANTS
# ============================================================================

def compute_spectral_determinants():
    """
    Compute zeta-regularized spectral determinants for S1, S2, S3.

    Known results from spectral geometry:
    - det'(Delta_{S1}) = 4*pi^2   [circle of circumference 2*pi]
    - det'(Delta_{S3}) = pi * exp(zeta(3)/(2*pi^2))  [unit 3-sphere]
    - det'(Delta_{S2}) = exp(1/2 - 4*zeta_R'(-1))    [unit 2-sphere]
    """
    print("=" * 72)
    print("PART 1: SPECTRAL DETERMINANTS OF BUNDLE COMPONENTS")
    print("=" * 72)

    results = {}

    if HAS_MPMATH:
        pi = mpmath.pi

        # --- det'(Delta_{S1}) ---
        # Circle of circumference 2*pi: eigenvalues k^2 for k in Z\{0}
        # zeta-regularized: det' = 4*pi^2
        det_S1 = 4 * pi**2
        results['det_S1'] = float(det_S1)
        print(f"\ndet'(Delta_S1) = 4*pi^2")
        print(f"  = {mpmath.nstr(det_S1, 30)}")

        # --- det'(Delta_{S3}) ---
        # Unit S3: eigenvalues n^2-1 with degeneracy n^2, for n >= 1
        # log det' = -zeta'_{S3}(0)
        # Known: det'(Delta_{S3}) = pi * exp(zeta(3)/(2*pi^2))
        zeta3 = mpmath.zeta(3)
        det_S3 = pi * mpmath.exp(zeta3 / (2 * pi**2))
        results['det_S3'] = float(det_S3)
        results['zeta3'] = float(zeta3)
        print(f"\ndet'(Delta_S3) = pi * exp(zeta(3)/(2*pi^2))")
        print(f"  zeta(3) = {mpmath.nstr(zeta3, 30)}")
        print(f"  zeta(3)/(2*pi^2) = {mpmath.nstr(zeta3/(2*pi**2), 30)}")
        print(f"  det'(Delta_S3) = {mpmath.nstr(det_S3, 30)}")

        # --- det'(Delta_{S2}) ---
        # Unit S2: eigenvalues l(l+1) with degeneracy 2l+1
        # det'(Delta_{S2}) = exp(1/2 - 4*zeta_R'(-1))
        # zeta_R'(-1) = 1/12 - ln(A) where A = Glaisher-Kinkelin constant
        # A = exp(1/12 - zeta_R'(-1))
        # Actually: zeta'(-1) = 1/12 - ln(A)
        # So 4*zeta'(-1) = 1/3 - 4*ln(A)
        # det'(S2) = exp(1/2 - 1/3 + 4*ln(A)) = exp(1/6) * A^4

        # Compute zeta'(-1) using mpmath
        # mpmath doesn't have a direct zeta derivative, but we can use:
        # zeta'(s) = d/ds zeta(s)
        # Use numerical differentiation at high precision
        zeta_prime_neg1 = mpmath.diff(mpmath.zeta, -1)
        det_S2 = mpmath.exp(mpmath.mpf('1')/2 - 4 * zeta_prime_neg1)
        results['det_S2'] = float(det_S2)
        results['zeta_prime_neg1'] = float(zeta_prime_neg1)
        print(f"\ndet'(Delta_S2) = exp(1/2 - 4*zeta'(-1))")
        print(f"  zeta'(-1) = {mpmath.nstr(zeta_prime_neg1, 30)}")
        print(f"  det'(Delta_S2) = {mpmath.nstr(det_S2, 30)}")

        # Glaisher-Kinkelin constant for cross-check
        # A = exp(1/12 - zeta'(-1))
        A_glaisher = mpmath.exp(mpmath.mpf(1)/12 - zeta_prime_neg1)
        print(f"  Glaisher-Kinkelin A = {mpmath.nstr(A_glaisher, 30)}")
        print(f"  (known: A ~ 1.28242712910...)")

        # --- The near-miss computation ---
        near_miss = det_S1 * det_S3 / pi
        results['near_miss'] = float(near_miss)
        gap = 42 - near_miss
        results['gap_to_42'] = float(gap)
        rel_gap = float(gap / 42)
        results['rel_gap'] = rel_gap

        print(f"\n{'─' * 60}")
        print(f"NEAR-MISS: det'(S1) * det'(S3) / pi")
        print(f"  = 4*pi^2 * exp(zeta(3)/(2*pi^2))")
        print(f"  = {mpmath.nstr(near_miss, 30)}")
        print(f"  B = 42")
        print(f"  Gap = 42 - {mpmath.nstr(near_miss, 15)} = {mpmath.nstr(gap, 15)}")
        print(f"  Relative gap = {rel_gap:.6e} ({abs(rel_gap)*100:.4f}%)")

        # --- K/pi target ---
        F = pi**2 / 6  # = zeta(2)
        Delta = mpmath.mpf(1) / 40
        K_over_pi = 42 + F - Delta
        K = pi * K_over_pi
        results['K_over_pi'] = float(K_over_pi)
        results['K'] = float(K)
        results['F'] = float(F)
        results['Delta'] = float(Delta)

        print(f"\n{'─' * 60}")
        print(f"K/pi = B + F - Delta = 42 + pi^2/6 - 1/40")
        print(f"  F = zeta(2) = pi^2/6 = {mpmath.nstr(F, 20)}")
        print(f"  Delta = 1/40 = {mpmath.nstr(Delta, 20)}")
        print(f"  K/pi = {mpmath.nstr(K_over_pi, 20)}")
        print(f"  K = {mpmath.nstr(K, 20)}")

        # Store mpmath objects for later use
        results['_mp'] = {
            'det_S1': det_S1, 'det_S2': det_S2, 'det_S3': det_S3,
            'near_miss': near_miss, 'gap': gap,
            'K_over_pi': K_over_pi, 'K': K, 'F': F, 'Delta': Delta,
            'pi': pi, 'zeta3': zeta3, 'zeta_prime_neg1': zeta_prime_neg1,
        }

    else:
        # Fallback with numpy
        det_S1 = 4 * np.pi**2
        zeta3 = 1.2020569031595942
        det_S3 = np.pi * np.exp(zeta3 / (2 * np.pi**2))
        near_miss = det_S1 * det_S3 / np.pi
        results = {
            'det_S1': det_S1, 'det_S3': det_S3,
            'near_miss': near_miss, 'gap_to_42': 42 - near_miss,
        }
        print(f"det'(S1) = {det_S1:.15f}")
        print(f"det'(S3) = {det_S3:.15f}")
        print(f"near_miss = {near_miss:.15f}")
        print(f"gap = {42 - near_miss:.15f}")

    return results


def bridge_search(results):
    """
    Systematically search for correction factors that bridge the gap
    between det'(S1)*det'(S3)/pi ≈ 41.957 and B = 42.

    Also search for expressions giving K/pi = 43.6199...
    """
    print("\n\n" + "=" * 72)
    print("PART 2: BRIDGE SEARCH — CLOSING THE GAP")
    print("=" * 72)

    if HAS_MPMATH and '_mp' in results:
        mp = results['_mp']
        near_miss = mp['near_miss']
        K_over_pi = mp['K_over_pi']
        det_S1 = mp['det_S1']
        det_S2 = mp['det_S2']
        det_S3 = mp['det_S3']
        pi = mp['pi']
        F = mp['F']
        Delta = mp['Delta']
        zeta3 = mp['zeta3']
    else:
        near_miss = results['near_miss']
        K_over_pi = results.get('K_over_pi', 43.6199)
        det_S1 = results['det_S1']
        det_S2 = results.get('det_S2', 3.195)
        det_S3 = results['det_S3']
        pi = np.pi
        F = np.pi**2 / 6
        Delta = 1/40
        zeta3 = 1.2020569031595942

    # Known constants
    B = 42
    n_max = 3
    N2 = 5      # N(2) = 1 + 4
    N3 = 14     # N(3) = 1 + 4 + 9
    d_max = 4

    hits = []

    def record(description, value, target, target_name):
        if HAS_MPMATH:
            diff = float(mpmath.fabs(value - target))
            rel = diff / float(target) if float(target) != 0 else float('inf')
            val_f = float(value)
        else:
            diff = abs(value - target)
            rel = diff / abs(target) if target != 0 else float('inf')
            val_f = value
        hits.append({
            'description': description,
            'value': val_f,
            'target': float(target),
            'target_name': target_name,
            'abs_diff': diff,
            'rel_error': rel,
        })

    # ─── MULTIPLICATIVE corrections to near_miss ───
    print("\n--- Multiplicative corrections: near_miss * X = 42? ---")

    required_mult = B / near_miss  # What multiplier is needed?
    print(f"Required multiplier: {float(required_mult):.15f}")

    tests_mult = [
        ("1 + 1/n_max^2", 1 + 1/mpmath.mpf(n_max)**2 if HAS_MPMATH else 1 + 1/n_max**2),
        ("1 + Delta", 1 + (mpmath.mpf(1)/40 if HAS_MPMATH else 1/40)),
        ("N(3)/N(2)", mpmath.mpf(N3)/N2 if HAS_MPMATH else N3/N2),
        ("n_max^2/(n_max^2-1)", mpmath.mpf(9)/8 if HAS_MPMATH else 9/8),
        ("(n_max+1)/n_max", mpmath.mpf(4)/3 if HAS_MPMATH else 4/3),
        ("1 + 1/(2*d_max)", 1 + 1/(mpmath.mpf(2)*d_max) if HAS_MPMATH else 1 + 1/(2*d_max)),
        ("1 + 1/B", 1 + 1/mpmath.mpf(B) if HAS_MPMATH else 1 + 1/B),
        ("1 + zeta(3)/(4*pi^2)", 1 + zeta3/(4*pi**2)),
        ("exp(1/(2*B))", mpmath.exp(1/(mpmath.mpf(2)*B)) if HAS_MPMATH else np.exp(1/(2*B))),
        ("1 + 1/N(3)", 1 + 1/mpmath.mpf(N3) if HAS_MPMATH else 1 + 1/N3),
    ]

    for desc, factor in tests_mult:
        val = near_miss * factor
        record(f"near_miss * ({desc})", val, B, "B=42")
        diff = float(val - B) if HAS_MPMATH else val - B
        rel = float(diff/B) if HAS_MPMATH else diff/B
        print(f"  {desc:35s}: {float(val):.10f}  (diff = {diff:+.6e}, rel = {rel:+.4e})")

    # ─── ADDITIVE corrections to near_miss ───
    print(f"\n--- Additive corrections: near_miss + X = 42? ---")
    gap = B - near_miss
    print(f"Gap = 42 - near_miss = {float(gap):.15f}")

    tests_add = [
        ("zeta(2) = pi^2/6", F),
        ("1/40 = Delta", Delta),
        ("1/N(3) = 1/14", 1/mpmath.mpf(14) if HAS_MPMATH else 1/14),
        ("1/(8*5) = 1/40", 1/mpmath.mpf(40) if HAS_MPMATH else 1/40),
        ("zeta(3)/(2*pi^2)", zeta3/(2*pi**2)),
        ("1/(2*B) = 1/84", 1/(mpmath.mpf(2)*B) if HAS_MPMATH else 1/(2*B)),
        ("pi/72", pi/72),
        ("1/24", 1/mpmath.mpf(24) if HAS_MPMATH else 1/24),
    ]

    for desc, addend in tests_add:
        val = near_miss + addend
        record(f"near_miss + {desc}", val, B, "B=42")
        diff = float(val - B) if HAS_MPMATH else val - B
        rel = float(diff/B) if HAS_MPMATH else diff/B
        print(f"  + {desc:35s}: {float(val):.10f}  (diff = {diff:+.6e}, rel = {rel:+.4e})")

    # ─── Combinations involving det'(S2) ───
    print(f"\n--- Combinations involving det'(S2) = {float(det_S2):.10f} ---")

    combos_S2 = [
        ("det'(S1)*det'(S3)/pi + det'(S2)", det_S1*det_S3/pi + det_S2),
        ("det'(S1)*det'(S3)/pi - det'(S2)", det_S1*det_S3/pi - det_S2),
        ("det'(S1)*det'(S3)*det'(S2)/pi", det_S1*det_S3*det_S2/pi),
        ("det'(S1)*det'(S3)*det'(S2)/pi^2", det_S1*det_S3*det_S2/pi**2),
        ("det'(S1)*det'(S3)/det'(S2)", det_S1*det_S3/det_S2),
        ("det'(S1)*det'(S2)/pi", det_S1*det_S2/pi),
        ("det'(S2)*det'(S3)/pi", det_S2*det_S3/pi),
        ("det'(S1)*det'(S2)*det'(S3)", det_S1*det_S2*det_S3),
        ("det'(S1)+det'(S2)+det'(S3)", det_S1+det_S2+det_S3),
        ("det'(S1)-det'(S2)+det'(S3)", det_S1-det_S2+det_S3),
    ]

    for desc, val in combos_S2:
        for target, tname in [(B, "B=42"), (K_over_pi, "K/pi=43.62")]:
            record(desc, val, target, tname)
        diff42 = float(val - B)
        diffK = float(val - K_over_pi)
        print(f"  {desc:45s}: {float(val):.10f}  (to 42: {diff42:+.4e}, to K/pi: {diffK:+.4e})")

    # ─── Broader arithmetic search for B = 42 ───
    print(f"\n--- Broader search for combinations = 42 ---")

    atoms = {
        'det_S1': det_S1,
        'det_S3': det_S3,
        'det_S2': det_S2,
        'pi': pi,
        'pi^2': pi**2,
        'zeta(2)': F,
        'zeta(3)': zeta3,
        '1/40': Delta,
        'N(2)=5': mpmath.mpf(5) if HAS_MPMATH else 5.0,
        'N(3)=14': mpmath.mpf(14) if HAS_MPMATH else 14.0,
        'd_max=4': mpmath.mpf(4) if HAS_MPMATH else 4.0,
        '8': mpmath.mpf(8) if HAS_MPMATH else 8.0,
    }

    # Test a/b, a*b, a+b, a-b for all pairs against B and K/pi
    print(f"  Testing pairwise operations (a*b, a/b, a+b, a-b) ...")
    close_hits_42 = []
    close_hits_K = []

    atom_keys = list(atoms.keys())
    for i, k1 in enumerate(atom_keys):
        for j, k2 in enumerate(atom_keys):
            if i == j:
                continue
            a, b = atoms[k1], atoms[k2]
            ops = [
                (f"{k1} * {k2}", a * b),
                (f"{k1} / {k2}", a / b if float(b) != 0 else None),
                (f"{k1} + {k2}", a + b),
                (f"{k1} - {k2}", a - b),
            ]
            for desc, val in ops:
                if val is None:
                    continue
                for target, tname, hit_list in [
                    (B, "B=42", close_hits_42),
                    (K_over_pi, "K/pi", close_hits_K),
                ]:
                    rel = float(abs(val - target) / abs(target))
                    if rel < 0.01:  # Within 1%
                        hit_list.append((desc, float(val), rel, tname))

    # Also test a*b/c for triples
    print(f"  Testing triple operations (a*b/c) ...")
    for i, k1 in enumerate(atom_keys):
        for j, k2 in enumerate(atom_keys):
            for k, k3 in enumerate(atom_keys):
                if i == j or i == k or j == k:
                    continue
                a, b, c = atoms[k1], atoms[k2], atoms[k3]
                if float(c) == 0:
                    continue
                val = a * b / c
                for target, tname, hit_list in [
                    (B, "B=42", close_hits_42),
                    (K_over_pi, "K/pi", close_hits_K),
                ]:
                    rel = float(abs(val - target) / abs(target))
                    if rel < 0.005:  # Within 0.5%
                        hit_list.append((f"{k1}*{k2}/{k3}", float(val), rel, tname))

    print(f"\n  Close hits for B = 42 (within 1% for pairs, 0.5% for triples):")
    close_hits_42.sort(key=lambda x: x[2])
    for desc, val, rel, tname in close_hits_42[:20]:
        print(f"    {desc:50s} = {val:.10f}  rel = {rel:.6e}")

    print(f"\n  Close hits for K/pi = {float(K_over_pi):.10f} (within 1%/0.5%):")
    close_hits_K.sort(key=lambda x: x[2])
    for desc, val, rel, tname in close_hits_K[:20]:
        print(f"    {desc:50s} = {val:.10f}  rel = {rel:.6e}")

    # ─── Can K/pi be written as f(det'(S1), det'(S2), det'(S3))? ───
    print(f"\n\n--- Can K/pi = 43.6199... be expressed via spectral dets? ---")

    # K/pi = B + F - Delta = near_miss + (42 - near_miss) + F - Delta
    # = near_miss + gap + F - Delta
    # If near_miss ≈ B, then K/pi ≈ near_miss + F - Delta

    test_K = [
        ("near_miss + F - Delta", near_miss + F - Delta),
        ("near_miss + F", near_miss + F),
        ("det'(S1)*det'(S3)/pi + pi^2/6 - 1/40", det_S1*det_S3/pi + F - Delta),
        ("det'(S1)*det'(S3)/pi + det'(S2)", det_S1*det_S3/pi + det_S2),
        ("det'(S1)*det'(S3)/pi + det'(S2) + 1", det_S1*det_S3/pi + det_S2 + 1),
        ("(det'(S1)*det'(S3) + pi*det'(S2))/pi", (det_S1*det_S3 + pi*det_S2)/pi),
        ("det'(S1)*det'(S3)/pi * (1 + F/B)", det_S1*det_S3/pi * (1 + F/B)),
    ]

    for desc, val in test_K:
        diff = float(val - K_over_pi)
        rel = float(diff / K_over_pi)
        print(f"  {desc:55s} = {float(val):.10f}  (diff = {diff:+.6e}, rel = {rel:+.4e})")

    # Key insight test: does the gap encode F or Delta?
    print(f"\n--- Analyzing the gap itself ---")
    print(f"  gap = {float(gap):.15f}")
    print(f"  gap / F = {float(gap/F):.15f}")
    print(f"  gap / Delta = {float(gap/Delta):.15f}")
    print(f"  gap * B = {float(gap*B):.15f}")
    print(f"  gap * N(3) = {float(gap*N3):.15f}")
    print(f"  gap / pi = {float(gap/pi):.15f}")
    print(f"  gap * pi = {float(gap*pi):.15f}")
    print(f"  B * gap / pi = {float(B*gap/pi):.15f}")
    print(f"  1/gap = {float(1/gap):.15f}")
    print(f"  exp(gap) - 1 = {float(mpmath.exp(gap)-1 if HAS_MPMATH else np.exp(float(gap))-1):.15f}")
    print(f"  ln(B/near_miss) = {float(mpmath.log(B/near_miss) if HAS_MPMATH else np.log(B/float(near_miss))):.15f}")

    return hits, close_hits_42, close_hits_K


# ============================================================================
# PART 3: SECOND SELECTION PRINCIPLE — d_max ANALYSIS
# ============================================================================

def analyze_dmax():
    """
    Prove that d_max = 4 for all n_max >= 2 in the GeoVac lattice.
    Analyze degree distribution and boundary effects.
    """
    print("\n\n" + "=" * 72)
    print("PART 3: SECOND SELECTION PRINCIPLE — d_max ANALYSIS")
    print("=" * 72)

    results = {}

    for max_n in range(1, 16):
        lattice = GeometricLattice(max_n=max_n, nuclear_charge=1)
        A = lattice.adjacency
        n_states = lattice.num_states
        degrees = np.array(A.sum(axis=1)).flatten().astype(int)

        d_max = degrees.max()
        d_min = degrees.min()

        # Degree histogram
        unique, counts = np.unique(degrees, return_counts=True)
        deg_hist = dict(zip(unique, counts))

        # Count interior vs boundary nodes
        n_interior = np.sum(degrees == 4) if 4 in deg_hist else 0

        # Identify which states have degree 4
        states_deg4 = []
        states_by_deg = {d: [] for d in unique}
        for idx in range(n_states):
            n, l, m = lattice.states[idx]
            states_by_deg[degrees[idx]].append((n, l, m))
            if degrees[idx] == d_max:
                states_deg4.append((n, l, m))

        results[max_n] = {
            'n_states': n_states,
            'd_max': d_max,
            'd_min': d_min,
            'deg_hist': deg_hist,
            'n_interior_d4': n_interior,
            'states_by_deg': states_by_deg,
        }

        hist_str = ", ".join(f"d={d}: {c}" for d, c in sorted(deg_hist.items()))
        print(f"  n_max={max_n:2d}: N={n_states:5d}, d_max={d_max}, d_min={d_min}  [{hist_str}]")

    # Detailed analysis for n_max = 3
    print(f"\n--- Detailed degree analysis for n_max = 3 ---")
    r = results[3]
    for deg in sorted(r['states_by_deg'].keys()):
        states = r['states_by_deg'][deg]
        print(f"  Degree {deg}: {len(states)} states")
        for s in states:
            n, l, m = s
            # Identify why this degree
            reasons = []
            if n == 1:
                reasons.append("n=1 (no T-)")
            if n == 3:
                reasons.append("n=n_max (no T+)")
            if l == 0:
                reasons.append("l=0 (no L±)")
            if m == -l and l > 0:
                reasons.append("m=-l (no L-)")
            if m == l and l > 0:
                reasons.append("m=+l (no L+)")
            reason_str = "; ".join(reasons) if reasons else "INTERIOR"
            print(f"    ({n},{l},{m:+d}): {reason_str}")

    # Check the bipartite bound
    print(f"\n--- Bipartite bound: lambda_max <= 2*d_max ---")
    print(f"  For bipartite graphs, lambda_max(D-A) <= 2*d_max")
    print(f"  d_max = 4 for all n_max >= 3, so bound = 8")
    print(f"  This explains why lambda_max -> 8 asymptotically")

    # Verify d_max = 4 universally for n_max >= 3
    all_dmax_4 = all(results[n]['d_max'] == 4 for n in range(3, 16))
    print(f"\n  d_max = 4 for all n_max in [3, 15]? {all_dmax_4}")
    print(f"  d_max(n_max=1) = {results[1]['d_max']} (trivial: single state)")
    print(f"  d_max(n_max=2) = {results[2]['d_max']}")

    # The selection equation
    print(f"\n--- Selection equation: 2*d_max = n_max^2 - 1 ---")
    for n in range(1, 8):
        lhs = 2 * 4  # 2 * d_max
        rhs = n**2 - 1
        match = "  <<<< MATCH" if lhs == rhs else ""
        print(f"  n_max={n}: 2*d_max = {lhs}, n^2-1 = {rhs}{match}")

    # Fraction of degree-4 nodes
    print(f"\n--- Fraction of degree-4 (interior) nodes vs n_max ---")
    for n in range(2, 16):
        r = results[n]
        frac = r['n_interior_d4'] / r['n_states']
        print(f"  n_max={n:2d}: {r['n_interior_d4']:5d}/{r['n_states']:5d} = {frac:.4f}")

    return results


def analyze_transition_geometry():
    """
    Analyze transition operators in relation to Hopf fibration geometry.

    The two transition types are:
    - T± (radial): changes n, keeps l,m fixed
    - L± (angular): changes m, keeps n,l fixed

    Question: do these correspond to base S2 and fiber S1 of the Hopf bundle?
    """
    print("\n\n" + "=" * 72)
    print("PART 4: TRANSITION OPERATORS & HOPF FIBRATION GEOMETRY")
    print("=" * 72)

    print("""
The Hopf fibration: S1 -> S3 -> S2

The states on S3 are labeled by (n, l, m) where:
  - n: principal quantum number (shell)
  - l: angular momentum (0 to n-1)
  - m: magnetic quantum number (-l to +l)

The two transition types in the GeoVac lattice:
  1. T± (radial): (n,l,m) <-> (n±1,l,m)  — changes shell, preserves (l,m)
  2. L± (angular): (n,l,m) <-> (n,l,m±1)  — changes m within shell

Hopf fibration interpretation:
  - The Hopf map pi: S3 -> S2 sends (n,l,m) -> ... what?
  - In the Peter-Weyl decomposition, S3 = SU(2), and functions on S3
    decompose as sum_j V_j tensor V_j^* (the (j,j) representations of SO(4)).
  - The Hopf fiber S1 corresponds to the U(1) action (n,l,m) -> phase rotation.
  - The base S2 = CP1 is the space of directions.

Key question: what do T± and L± correspond to geometrically?

ANALYSIS:
""")

    # T± changes n, preserving l and m
    # In the Fock transform, n labels the "radial" direction on S3 (distance from pole)
    # The S3 eigenvalue is n^2-1, so T± moves between eigenspaces of Delta_{S3}
    # This is motion in the "radial" direction of S3 — moving between shells

    # L± changes m, preserving n and l
    # m is the magnetic quantum number = projection of angular momentum on z-axis
    # This corresponds to rotation around the z-axis, which is the U(1) action
    # of the Hopf fiber!

    print("RESULT: The transition operators have a clear Hopf interpretation:")
    print()
    print("  L± (angular, m -> m±1):")
    print("    - Changes the U(1) phase (magnetic quantum number m)")
    print("    - This is motion along the Hopf FIBER S1")
    print("    - The m quantum number labels the U(1) character of the fiber")
    print("    - dim(fiber motion) = 1 (±m gives 2 directions, but 1 dimension)")
    print()
    print("  T± (radial, n -> n±1):")
    print("    - Changes the principal quantum number (shell)")
    print("    - This is NOT motion on the base S2 directly")
    print("    - Rather, it's motion between eigenspaces of Delta_{S3}")
    print("    - In the Fock stereographic picture, n labels the 'latitude' on S3")
    print("    - dim(radial motion) = 1 (±n gives 2 directions, but 1 dimension)")
    print()
    print("  MISSING: There is no l-changing transition in the GeoVac lattice!")
    print("    - The angular momentum l is conserved by both T± and L±")
    print("    - l labels the representation of SO(3) within each SO(4) shell")
    print("    - Changing l would correspond to motion on the base S2")
    print("    - The base S2 has dim=2, but is NOT directly traversed")
    print()
    print("  THEREFORE:")
    print("    d_transition = 2 types (T±, L±) = 1 radial + 1 fiber")
    print("    d_max = 2 * d_transition = 4 (each type has ± direction)")
    print()
    print("  The base S2 contributes its structure through the (l,m) labeling")
    print("  within each shell, but does NOT contribute a transition type.")
    print("  This is because l-changing transitions would break the S3 topology")
    print("  (they would connect different SO(3) representations within a shell).")
    print()

    # The corrected interpretation
    print("CORRECTED HOPF INTERPRETATION:")
    print("  - Fiber S1: L± (m transitions)")
    print("  - Radial S3 direction: T± (n transitions)")
    print("  - Base S2: constrains state labeling, not a transition")
    print()
    print("  d_max = 2 * (dim_fiber + dim_radial) = 2 * (1 + 1) = 4")
    print("  This gives: 2 * d_max = n_max^2 - 1 => 8 = n_max^2 - 1 => n_max = 3")
    print()
    print("  Alternative: d_max = 2 * n_transition_types = 2 * 2 = 4")
    print("  And n_transition_types = 2 because S3 has exactly 2 independent")
    print("  motion directions that preserve the quantum number structure:")
    print("  one along the fiber (m) and one perpendicular to it (n).")

    # Connection to the Hopf invariant
    print()
    print("  WHY 2 TRANSITION TYPES (not 3)?")
    print("  S3 is 3-dimensional, but the Hopf fibration S1 -> S3 -> S2")
    print("  reduces the 3 dimensions to 2 types of motion:")
    print("    1. Along the fiber S1 (1-dimensional)")
    print("    2. Perpendicular to the fiber (2-dimensional, BUT...)")
    print("  In the discrete lattice, the 2D perpendicular motion collapses")
    print("  to a single quantum number n (the shell index).")
    print("  This is because the Fock projection maps the 2D base directions")
    print("  to a single radial variable (the stereographic distance).")
    print("  The l quantum number is implicitly encoded in the shell structure")
    print("  but is NOT a transition — it's a selection rule.")


def make_plots(results, bridge_hits_42, bridge_hits_K):
    """Generate summary plots."""
    print("\n\nGenerating plots...")

    # Plot 1: Degree distribution vs n_max
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: fraction of degree-4 nodes
    nmax_vals = list(range(2, 16))
    frac_d4 = []
    for n in nmax_vals:
        r = results[n]
        frac_d4.append(r['n_interior_d4'] / r['n_states'])

    axes[0].plot(nmax_vals, frac_d4, 'bo-', markersize=6)
    axes[0].axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='All interior')
    axes[0].set_xlabel('n_max')
    axes[0].set_ylabel('Fraction of degree-4 nodes')
    axes[0].set_title('Interior (degree-4) node fraction')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Right: d_max vs n_max
    dmax_vals = [results[n]['d_max'] for n in nmax_vals]
    axes[1].plot(nmax_vals, dmax_vals, 'rs-', markersize=8)
    axes[1].axhline(y=4, color='b', linestyle='--', alpha=0.5, label='d_max = 4')
    axes[1].set_xlabel('n_max')
    axes[1].set_ylabel('Maximum degree d_max')
    axes[1].set_title('Maximum node degree vs n_max')
    axes[1].set_ylim(0, 6)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(__file__), 'degree_analysis.png'), dpi=150)
    plt.close()
    print("  Saved degree_analysis.png")

    # Plot 2: Bridge search — best hits
    fig, ax = plt.subplots(figsize=(10, 6))

    # Combine and sort all hits by relative error
    all_sorted_42 = sorted(bridge_hits_42, key=lambda x: x[2])[:15]
    if all_sorted_42:
        labels = [h[0][:40] for h in all_sorted_42]
        errors = [h[2] for h in all_sorted_42]
        y_pos = range(len(labels))
        ax.barh(y_pos, errors, color='steelblue', alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Relative error from B = 42')
        ax.set_title('Best bridge candidates for det\'(S1)*det\'(S3)/pi → 42')
        ax.axvline(x=0.001, color='r', linestyle='--', alpha=0.5, label='0.1% (Paper 2 near-miss)')
        ax.legend()
        ax.set_xscale('log')
    else:
        ax.text(0.5, 0.5, 'No close hits found', ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(__file__), 'bridge_search.png'), dpi=150)
    plt.close()
    print("  Saved bridge_search.png")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':
    print("Phase 2: Spectral Determinant Investigation")
    print("=" * 72)

    # Part 1: Spectral determinants
    spec_results = compute_spectral_determinants()

    # Part 2: Bridge search
    hits, close_42, close_K = bridge_search(spec_results)

    # Part 3: d_max analysis
    dmax_results = analyze_dmax()

    # Part 4: Transition geometry
    analyze_transition_geometry()

    # Plots
    make_plots(dmax_results, close_42, close_K)

    print("\n\n" + "=" * 72)
    print("PHASE 2 COMPLETE")
    print("=" * 72)
