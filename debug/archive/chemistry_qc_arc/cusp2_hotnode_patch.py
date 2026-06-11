"""
Track CUSP-2: (1s,1s) hot-node patch.

Idea: the cusp's contribution to the FCI energy is missing from any
finite n_max basis. The Schwartz partial-wave extrapolation says this
missing energy scales as A/(l_max+2)^4 where A = (10/pi)|psi(0)|^2.
EP-2c/2g identified the (1s,1s) pair-state as the unique cusp hot-node.

So: ADD the Schwartz-extrapolated tail directly to the diagonal of
the (1s,1s) singlet config in the CI matrix.

The minimal "patch":
  H_patched[I_00, I_00] += DeltaE_Schwartz(l_max)

where I_00 is the index of the (1s,1s) singlet config, and
  DeltaE_Schwartz(l_max) = -A / (l_max + 2)^4.

We test:
1. He at multiple n_max ∈ {3, 4, 5, 6, 7} with the patch.
2. Compare to the empirical extrapolation: E(l_max) - E_exact ~ A/(l_max+2)^4.
3. Does the patch break the 0.20% floor at n_max=7?

Output: debug/data/cusp2_hotnode_patch.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.casimir_ci import (
    build_graph_native_fci, _build_orbital_basis,
)


E_EXACT_HE = -2.903724377


def _find_1s1s_index(n_max, m_total=0):
    """Locate the index of the (1s,1s) singlet config."""
    orbitals = _build_orbital_basis(n_max)
    configs = []
    n_spatial = len(orbitals)
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                configs.append((i, j))
    s1s_idx = None
    for I, (i, j) in enumerate(configs):
        if i == j and orbitals[i] == (1, 0, 0):
            s1s_idx = I
            break
    return s1s_idx, len(configs)


def _l_max_at_n_max(n_max):
    """Largest l in the basis at given n_max."""
    return n_max - 1


def main():
    Z = 2
    print(f"CUSP-2: hot-node patch on He (Z={Z})")
    print("=" * 60)

    # Reference: unpatched at each n_max
    rows = []
    print("\nUnpatched graph-native CI:")
    for n_max in (2, 3, 4, 5, 6, 7):
        t0 = time.time()
        H = build_graph_native_fci(Z=Z, n_max=n_max)
        E0 = float(np.linalg.eigvalsh(H)[0])
        err = E0 - E_EXACT_HE
        err_pct = 100 * abs(err) / abs(E_EXACT_HE)
        n_cfg = H.shape[0]
        s1s_idx, _ = _find_1s1s_index(n_max)
        l_max = _l_max_at_n_max(n_max)
        wall = time.time() - t0
        rows.append({
            'n_max': n_max, 'l_max': l_max, 'n_cfg': n_cfg,
            'E_unpatched': E0, 'err_pct': err_pct, 'err_abs': err,
            's1s_idx': s1s_idx, 'wall_s': wall, 'H': H,
        })
        print(f"  n_max={n_max} l_max={l_max} n_cfg={n_cfg:4d}  "
              f"E={E0:.6f}  err={err_pct:.4f}% ({wall:.1f}s)")

    # Calibrate Schwartz coefficient A from the data:
    # E(l_max) - E_exact ~ -A/(l_max+2)^4
    # Use n_max=4,5,6,7 (l_max=3,4,5,6) — high enough for asymptotic
    fit_rows = [r for r in rows if r['n_max'] >= 4]
    x = np.array([1.0 / (r['l_max'] + 2) ** 4 for r in fit_rows])
    y = np.array([r['err_abs'] for r in fit_rows])
    # E_exact = E_unpatched + DeltaE_cusp; DeltaE_cusp = E_exact - E_unpatched.
    # Negative because graph-native over-binds slightly at He (Z=2).
    # err_abs = E_unpatched - E_exact; so positive if over-binding (E_un < E_exact)
    # Schwartz: E_unpatched - E_exact = -A/(l_max+2)^4 if cusp is the dominant error
    # (which would mean E_unpatched > E_exact); but we have over-binding (E_unpatched < E_exact)
    # => err_abs is NEGATIVE. Let me check.
    print(f"\nFit calibration data:")
    for r in fit_rows:
        print(f"  l_max={r['l_max']}: err_abs = {r['err_abs']:.6e}  "
              f"(1/(l+2)^4 = {1/(r['l_max']+2)**4:.4e})")

    # Linear fit: err_abs = m * (1/(l+2)^4) + c
    if len(fit_rows) >= 2:
        m_fit, c_fit = np.polyfit(x, y, 1)
        print(f"  Linear fit slope (= -A): {m_fit:.6e}, intercept: {c_fit:.6e}")
        A_emp = -m_fit  # err_abs = -A/(l+2)^4 + c, so A = -slope
    else:
        m_fit = c_fit = A_emp = 0.0

    # Patch: subtract |err_abs at THIS l_max| from H[s1s_idx, s1s_idx]
    # Wait — the IDEA is that adding cusp correction to (1s,1s) diagonal
    # SHIFTS the GS energy by some amount that we can fit.
    # If we add Delta to H[s1s,s1s], the GS shift to first order is
    # |c_s1s_GS|^2 * Delta where c_s1s_GS is the GS amplitude on (1s,1s).
    # For He GS, |c_s1s|^2 ~ 0.95+. So we need
    # Delta = err_abs / |c_s1s_GS|^2 to shift GS by exactly err_abs.
    print("\nPatched results (target: shift E_unpatched by -err_abs):")
    for r in rows:
        H = r['H'].copy()
        # First find current c_s1s amplitude in unpatched GS
        eigs_u, vecs_u = np.linalg.eigh(H)
        c_s1s = float(vecs_u[r['s1s_idx'], 0])
        c2 = c_s1s ** 2
        # Schwartz tail at this l_max (using fit)
        if A_emp != 0 and r['l_max'] >= 2:
            DeltaE_Schwartz = -A_emp / (r['l_max'] + 2) ** 4
        else:
            DeltaE_Schwartz = 0.0
        # Apply patch
        H[r['s1s_idx'], r['s1s_idx']] += DeltaE_Schwartz / c2 if c2 > 0.1 else 0.0
        eigs_p = np.linalg.eigvalsh(H)
        E_p = float(eigs_p[0])
        err_p = 100 * abs(E_p - E_EXACT_HE) / abs(E_EXACT_HE)
        improvement = r['err_pct'] - err_p
        print(f"  n_max={r['n_max']} l_max={r['l_max']}  "
              f"|c_s1s|^2={c2:.4f}  Delta_S={DeltaE_Schwartz:.3e} Ha  "
              f"E_p={E_p:.6f}  err_p={err_p:.4f}%  (delta {improvement:+.4f}pp)")
        r['E_patched'] = E_p
        r['err_p_pct'] = err_p
        r['c_s1s_GS_squared'] = c2
        r['DeltaE_Schwartz'] = DeltaE_Schwartz
        del r['H']  # don't serialize

    out = {
        'Z': Z,
        'A_empirical': float(A_emp),
        'rows': rows,
    }
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'cusp2_hotnode_patch.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")


if __name__ == '__main__':
    main()
