"""Precision catalogue: Helium 2^1S - 2^3S exchange splitting.

The first MULTI-ELECTRON CORRELATION precision test in the catalogue. Tests the
graph-native CI machinery on a singlet-triplet splitting in the (1s)(2s)
configuration, which is determined ENTIRELY by the V_ee exchange integral
2 K(1s,2s) — common-mode kinetic and one-body terms cancel between S=0 and S=1.

This complements He 2^3P (Sprint 2026-05-08), which tested fine-structure
splittings driven by spin-orbit + spin-spin + spin-other-orbit. Here, no
relativistic operators enter — the splitting is pure non-relativistic
electron-electron correlation.

Architecture
------------
- Same orbital basis as He ground state graph-native CI (Track DI Sprint 3C):
  hydrogenic (n,l,m) at exponent k=Z=2, hybrid one-body Hamiltonian (exact
  -Z^2/(2n^2) diagonal + kappa = -1/16 graph adjacency off-diagonal),
  analytical Slater V_ee from `geovac.casimir_ci.two_electron_integral`
  (Wigner-3j angular + hypergeometric_slater radial, machine-precision).
- Singlet sector (S=0, symmetric spatial): pairs (i, j) with i <= j;
  bra/ket spatial wavefunctions [|ij> + |ji>]/sqrt(2) for i != j, |ii> for i==j.
- Triplet sector (S=1, antisymmetric spatial): pairs (i, j) with i < j;
  bra/ket spatial wavefunctions [|ij> - |ji>]/sqrt(2). Both M_S=±1 and M_S=0
  triplet sublevels are degenerate by spin rotational symmetry; we use the
  antisymmetric spatial sector directly.
- State labeling: filter to the (l1=0, l2=0) "ss" sub-block of M_L=0. This
  isolates the pure S-state sector (1^1S, 2^1S, 3^1S, ...; 2^3S, 3^3S, ...)
  cleanly. Diagnostic shows the ss-block result is bit-stable vs the full
  M_L=0 sector (which mixes ^1S, ^1P_{M_L=0}, ^1D_{M_L=0}). The naive
  "sort the M_L=0 eigenvalues" labeling used by compute_he_spectrum() is
  WRONG for n_max>=3: at n_max=3, the second M_L=0 singlet eigenvalue is
  dominantly (1s,3p) ^1P, NOT (1s,2s) ^1S. The ss-block is the correct
  identification.

Key finding: graph-native CI does NOT preserve a localized (1s,2s)
eigenstate as n_max grows. The graph kappa adjacency connects 1s<->2s<->3s<->...,
so the (n=2 s-character) is fragmented across many eigenstates. The lowest
ss-singlet excited eigenvalue is the framework's best approximation to
2^1S, but its (1s,2s) projection drops from 85% (n_max=2) to 30% (n_max=6).

Reference values
----------------
NIST Atomic Spectra Database (Drake high-precision NR limit):
  E(1^1S) = 0 cm^-1 (definition of zero)
  E(2^3S_1) at 159855.9743 cm^-1 above 1^1S
  E(2^1S_0) at 166277.4406 cm^-1 above 1^1S
  Splitting: nu(2^1S - 2^3S) = 6421.4663 cm^-1 = 192.510 THz = 0.029260 Ha

(Source: NIST ASD, accessed 2026; consistent with Drake & Yan 1992 high-precision
non-relativistic infinite-mass values listed in HE_NR_REFERENCE in casimir_ci.py:
  E(1_1S) = -2.903724377  Ha
  E(2_3S) = -2.175229378  Ha
  E(2_1S) = -2.145974046  Ha
  -> dE_st = 0.029255 Ha = 6420.4 cm^-1 (matches NIST to 1 cm^-1)

We use NIST 6421.4663 cm^-1 as the reference.

What the framework reproduces
-----------------------------
At leading order, the singlet-triplet splitting is dE_st = 2 K(1s,2s) where
K(a,b) = G^0(a,b) is the exchange (transition-density) Slater integral. From
casimir_ci.py: G^0(1s, 2s) at k=Z=2 = 2 * (16/729) Ha = 32/729 Ha = 0.04389 Ha
~ 9633 cm^-1 (a 50% overestimate vs NIST). This is the bare HF-level number;
graph-native CI adds correlation through:
  (i) the graph off-diagonal h1 (kappa * (-A) connecting (n,l,m) <-> (n',l',m))
      which mixes 1s with 2s, 3s, etc. and softens the orbital
  (ii) configuration mixing in the FCI sector (1s 2s) <-> (1s 3s) <-> (2s 2s)
       <-> (2p 2p, M_L=0) etc.
Both effects systematically REDUCE the exchange splitting from the bare
HF value, bringing it toward 6421 cm^-1.

What the framework does NOT include
-----------------------------------
Pure Coulomb non-relativistic limit only. No relativistic, QED, or recoil
corrections (which are <1 cm^-1 for the 2^1S-2^3S splitting and well below
framework precision). The 0.20% small-Z graph-validity-boundary artifact
(Z_c~1.84, He at Z=2 is just above) affects ABSOLUTE energies; the
splitting is a difference between two states with similar orbital character,
so common-mode errors largely cancel and the splitting may be MORE accurate
than the 0.20% absolute floor would suggest.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.casimir_ci import (
    HE_NR_REFERENCE,
    build_graph_native_fci,
    compute_he_spectrum,
)

# --- Physical constants ---
HA_TO_CM_INVERSE: float = 219474.6313632
HA_TO_HZ: float = 6.579683920502e15
HA_TO_THZ: float = HA_TO_HZ * 1.0e-12

# --- Reference values (NIST Atomic Spectra Database) ---
NIST_SPLITTING_CM_INVERSE: float = 6421.4663  # 2^1S_0 - 2^3S_1 (J-averaged)
DRAKE_E_1S_HA: float = -2.903724377  # 1^1S NR infinite-mass
DRAKE_E_2_3S_HA: float = -2.175229378
DRAKE_E_2_1S_HA: float = -2.145974046
DRAKE_SPLITTING_HA: float = DRAKE_E_2_1S_HA - DRAKE_E_2_3S_HA  # ~0.029255 Ha
DRAKE_SPLITTING_CM_INVERSE: float = DRAKE_SPLITTING_HA * HA_TO_CM_INVERSE


def _build_ss_only_singlet(n_max: int) -> np.ndarray:
    """Build the (l_1 = l_2 = 0) singlet sub-block FCI matrix only.

    This avoids constructing the full M_L=0 sector; the ss-block has
    dimension n_s_orbs * (n_s_orbs+1) / 2 where n_s_orbs = n_max
    (one s-orbital per shell n=1..n_max). At n_max=7: 28x28 matrix.
    """
    from geovac.casimir_ci import (_build_orbital_basis,
                                    _build_graph_h1, two_electron_integral)
    h1_mat, orbitals = _build_graph_h1(Z=2, n_max=n_max)

    # Indices of the s-orbitals (l=0, m=0 only) — these have one orbital per n
    s_indices = [i for i, (n, l, m) in enumerate(orbitals) if l == 0]

    # Singlet (l_1=l_2=0) configurations: pairs (i, j) with i <= j (M_L=0 trivially)
    configs = [(i, j) for i in s_indices for j in s_indices if i <= j]
    n_configs = len(configs)
    H = np.zeros((n_configs, n_configs))

    Z = 2
    k_orb = float(Z)

    def g_int(a, b, c, d):
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        return two_electron_integral(na, la, ma, nb, lb, mb,
                                     nc, lc, mc, nd, ld, md, k_orb)

    for I in range(n_configs):
        i, j = configs[I]
        for J in range(I, n_configs):
            p, q = configs[J]
            bra_perms = [(i, j, 1.0)]
            if i != j:
                bra_perms.append((j, i, 1.0))
            ket_perms = [(p, q, 1.0)]
            if p != q:
                ket_perms.append((q, p, 1.0))
            N_IJ = np.sqrt(float(len(bra_perms)))
            N_PQ = np.sqrt(float(len(ket_perms)))
            me = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    if b == d:
                        me += sign * h1_mat[a, c]
                    if a == c:
                        me += sign * h1_mat[b, d]
                    me += sign * g_int(a, b, c, d)
            me /= (N_IJ * N_PQ)
            H[I, J] = me
            H[J, I] = me
    return H


def _build_ss_only_triplet(n_max: int) -> np.ndarray:
    """Build the (l_1 = l_2 = 0) triplet sub-block FCI matrix only.

    Triplet requires i < j (antisymmetric spatial); diagonal pair (i,i) is
    Pauli-forbidden. At n_max=7: 21x21 matrix.
    """
    from geovac.casimir_ci import (_build_orbital_basis,
                                    _build_graph_h1, two_electron_integral)
    h1_mat, orbitals = _build_graph_h1(Z=2, n_max=n_max)
    s_indices = [i for i, (n, l, m) in enumerate(orbitals) if l == 0]
    configs = [(i, j) for k1, i in enumerate(s_indices)
                       for j in s_indices[k1+1:]]
    n_configs = len(configs)
    if n_configs == 0:
        return np.zeros((0, 0))
    H = np.zeros((n_configs, n_configs))

    Z = 2
    k_orb = float(Z)
    parity = -1.0  # triplet sign

    def g_int(a, b, c, d):
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        return two_electron_integral(na, la, ma, nb, lb, mb,
                                     nc, lc, mc, nd, ld, md, k_orb)

    for I in range(n_configs):
        i, j = configs[I]
        for J in range(I, n_configs):
            p, q = configs[J]
            bra_perms = [(i, j, 1.0), (j, i, parity)]
            ket_perms = [(p, q, 1.0), (q, p, parity)]
            N_IJ = np.sqrt(2.0)
            N_PQ = np.sqrt(2.0)
            me = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    if b == d:
                        me += sign * h1_mat[a, c]
                    if a == c:
                        me += sign * h1_mat[b, d]
                    me += sign * g_int(a, b, c, d)
            me /= (N_IJ * N_PQ)
            H[I, J] = me
            H[J, I] = me
    return H


def compute_st_at_nmax(n_max: int) -> Dict[str, Any]:
    """Compute E(1^1S), E(2^3S), E(2^1S), and the singlet-triplet splitting.

    Uses the (l_1 = l_2 = 0) ss-only sub-block to extract pure S-state
    eigenvalues cleanly. This avoids the L-mixing that occurs in the
    full M_L=0 sector and gives a STABLE labeling of 1^1S, 2^1S, 2^3S
    that does not depend on accidental energy ordering.

    Returns
    -------
    dict with keys:
      n_max, dim_singlet_ss, dim_triplet_ss, t_build, t_solve,
      E_1_1S_Ha, E_2_3S_Ha, E_2_1S_Ha,
      dE_st_Ha, dE_st_cm_inverse,
      err_vs_NIST_pct, err_vs_NIST_cm_inverse,
      err_vs_Drake_pct, err_vs_Drake_cm_inverse,
      err_E_1_1S_pct, err_E_2_3S_pct, err_E_2_1S_pct,
      variational_*_ok, ordering_ok, ...
    """
    # --- Singlet ss-block: (l1=l2=0, M_L=0, S=0) ---
    t0 = time.time()
    H_s = _build_ss_only_singlet(n_max)
    t_build_s = time.time() - t0
    t0 = time.time()
    evals_s = np.sort(np.linalg.eigvalsh(H_s))
    t_solve_s = time.time() - t0
    dim_singlet = H_s.shape[0]

    # --- Triplet ss-block: (l1=l2=0, M_L=0, S=1) ---
    t0 = time.time()
    H_t = _build_ss_only_triplet(n_max)
    t_build_t = time.time() - t0
    t0 = time.time()
    if H_t.shape[0] > 0:
        evals_t = np.sort(np.linalg.eigvalsh(H_t))
    else:
        evals_t = np.array([])
    t_solve_t = time.time() - t0
    dim_triplet = H_t.shape[0]

    # State extraction in the pure ss block:
    # Singlet ss: index 0 -> 1^1S, index 1 -> 2^1S, ...
    # Triplet ss: index 0 -> 2^3S, index 1 -> 3^3S, ...
    # (No 1^3S exists: would require both electrons in 1s with parallel spin,
    # forbidden by Pauli.)
    E_1_1S = float(evals_s[0])
    E_2_1S = float(evals_s[1]) if dim_singlet >= 2 else None
    E_2_3S = float(evals_t[0]) if dim_triplet >= 1 else None

    # Splitting
    dE_st_Ha = (E_2_1S - E_2_3S) if (E_2_1S is not None and E_2_3S is not None) else None
    dE_st_cm = dE_st_Ha * HA_TO_CM_INVERSE if dE_st_Ha is not None else None
    dE_st_THz = dE_st_Ha * HA_TO_THZ if dE_st_Ha is not None else None

    # Errors
    if dE_st_cm is not None:
        err_NIST_pct = (dE_st_cm - NIST_SPLITTING_CM_INVERSE) / NIST_SPLITTING_CM_INVERSE * 100.0
        err_NIST_cm = dE_st_cm - NIST_SPLITTING_CM_INVERSE
        err_Drake_pct = (dE_st_cm - DRAKE_SPLITTING_CM_INVERSE) / DRAKE_SPLITTING_CM_INVERSE * 100.0
        err_Drake_cm = dE_st_cm - DRAKE_SPLITTING_CM_INVERSE
    else:
        err_NIST_pct = err_NIST_cm = err_Drake_pct = err_Drake_cm = None

    # Per-state errors vs Drake NR reference
    err_E_1_1S_pct = abs((E_1_1S - DRAKE_E_1S_HA) / DRAKE_E_1S_HA) * 100.0 if E_1_1S is not None else None
    err_E_2_3S_pct = (
        abs((E_2_3S - DRAKE_E_2_3S_HA) / DRAKE_E_2_3S_HA) * 100.0 if E_2_3S is not None else None
    )
    err_E_2_1S_pct = (
        abs((E_2_1S - DRAKE_E_2_1S_HA) / DRAKE_E_2_1S_HA) * 100.0 if E_2_1S is not None else None
    )

    # Variational and ordering checks (these are the diagnostic asserts)
    variational_1_1S = (E_1_1S is None) or (E_1_1S > DRAKE_E_1S_HA)
    variational_2_3S = (E_2_3S is None) or (E_2_3S > DRAKE_E_2_3S_HA)
    variational_2_1S = (E_2_1S is None) or (E_2_1S > DRAKE_E_2_1S_HA)
    # Hund: E(2^1S) > E(2^3S) — singlet HIGHER than triplet
    ordering_ok = (E_2_1S is not None and E_2_3S is not None and E_2_1S > E_2_3S)
    # Both excited states above ground:
    excited_above_ground_ok = (
        E_2_3S is not None and E_2_1S is not None
        and E_2_3S > E_1_1S and E_2_1S > E_1_1S
    )
    # Ionization threshold: He+ ground state = -2.0 Ha. Both 2^3S and 2^1S
    # must be below this (both are bound He neutral states).
    below_ionization_ok = (
        E_2_3S is not None and E_2_1S is not None
        and E_2_3S < -2.0 and E_2_1S < -2.0
    )

    return {
        'n_max': n_max,
        'dim_singlet': dim_singlet,
        'dim_triplet': dim_triplet,
        'dim_singlet_ss': dim_singlet,
        'dim_triplet_ss': dim_triplet,
        't_build_singlet_s': t_build_s,
        't_solve_singlet_s': t_solve_s,
        't_build_triplet_s': t_build_t,
        't_solve_triplet_s': t_solve_t,
        'E_1_1S_Ha': E_1_1S,
        'E_2_3S_Ha': E_2_3S,
        'E_2_1S_Ha': E_2_1S,
        'dE_st_Ha': dE_st_Ha,
        'dE_st_cm_inverse': dE_st_cm,
        'dE_st_THz': dE_st_THz,
        'err_vs_NIST_pct': err_NIST_pct,
        'err_vs_NIST_cm_inverse': err_NIST_cm,
        'err_vs_Drake_pct': err_Drake_pct,
        'err_vs_Drake_cm_inverse': err_Drake_cm,
        'err_E_1_1S_pct': err_E_1_1S_pct,
        'err_E_2_3S_pct': err_E_2_3S_pct,
        'err_E_2_1S_pct': err_E_2_1S_pct,
        'variational_1_1S_ok': variational_1_1S,
        'variational_2_3S_ok': variational_2_3S,
        'variational_2_1S_ok': variational_2_1S,
        'ordering_ok': ordering_ok,
        'excited_above_ground_ok': excited_above_ground_ok,
        'below_ionization_ok': below_ionization_ok,
    }


def estimate_eta_for_nmax(n_max: int, prev_t: float, prev_dim: int, dim: int) -> float:
    """Heuristic ETA: graph-native FCI build is dominated by O(dim^2) integral
    work in the configuration loop, where dim ~ n_orb^2. So total time scales
    roughly as dim^2 / prev_dim^2 * prev_t.
    """
    if prev_dim == 0:
        return 0.0
    return prev_t * (dim / prev_dim) ** 2


def main() -> Dict[str, Any]:
    print("=" * 78)
    print("Precision catalogue: Helium 2^1S - 2^3S exchange splitting")
    print("=" * 78)
    print()
    print("System:        He neutral, configurations 1s 2s")
    print("Architecture:  graph-native CI (Track DI Sprint 3C)")
    print("Spin sectors:  singlet (S=0, symmetric spatial), triplet (S=1, antisym)")
    print(f"Reference:     NIST  6421.4663 cm^-1 (192.510 THz, 0.029260 Ha)")
    print(f"               Drake {DRAKE_SPLITTING_CM_INVERSE:.4f} cm^-1 ({DRAKE_SPLITTING_HA:.6f} Ha)")
    print()
    print(f"{'n_max':>6} {'dim_S':>7} {'dim_T':>7} {'t_build (s)':>12}  "
          f"{'E(2^3S) Ha':>13} {'E(2^1S) Ha':>13} "
          f"{'dE_st (cm-1)':>14} {'err NIST (%)':>14}")
    print("-" * 100)

    # Determine which n_max values to attempt. The ss-only block scales as
    # O(n_max^4) (not O(n_orbs^4 ~ n_max^8) of the full M_L=0 sector), so we
    # can go to much larger n_max cheaply. n_max=8 ss block has ~36 configs.
    #
    # 2026-05-08 cleanup-sprint REFRESH: the original cap at n_max=11 was due
    # to a now-fixed bug in `hypergeometric_slater.compute_rk_float` (the
    # float path silently degraded for n >= 6 and produced wrong-sign values
    # at n >= 12 due to catastrophic cancellation between O(1e25) Laguerre-
    # product terms summing to O(1e11)).  The fix in place dispatches at
    # ``_FLOAT_PATH_MAX_N = 4``: n <= 4 uses the fast pure-float path bit-
    # identically; n >= 5 delegates to ``compute_rk_algebraic`` (exact
    # ``Fraction``) and casts to float at the end.  Both paths are wrapped
    # in the module-level cache so each unique quartet is paid for once per
    # CI build.  We extend through n_max=13 with a per-step timeout so the
    # run terminates if R^k evaluation grows beyond 10 minutes per step.
    NMAX_VALUES = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    MAX_N_MAX_HARD_CAP = 13
    PER_STEP_TIMEOUT_S = 1200.0  # 20 minutes (n=12,13 R^k integrals are slower)

    rows: List[Dict[str, Any]] = []
    prev_t_total = 0.0
    prev_dim = 0

    for n_max in NMAX_VALUES:
        # ETA: ss-block dim is n_max*(n_max+1)/2; cost ~ dim^2 * (n_max^2)
        # for the integral evaluation. Roughly O(n_max^6) in total work.
        est_dim = n_max * (n_max + 1) // 2
        eta = estimate_eta_for_nmax(n_max, prev_t_total, prev_dim, est_dim) if prev_dim else 0.0

        if eta > PER_STEP_TIMEOUT_S:
            print(f"\n[STOP] Estimated time at n_max={n_max} is {eta:.0f}s > {PER_STEP_TIMEOUT_S:.0f}s budget.")
            print(f"       Stopping at n_max={n_max - 1}.")
            break

        try:
            row = compute_st_at_nmax(n_max)
        except MemoryError:
            print(f"\n[STOP] MemoryError at n_max={n_max}. Stopping.")
            break

        prev_t_total = (row['t_build_singlet_s'] + row['t_build_triplet_s']
                        + row['t_solve_singlet_s'] + row['t_solve_triplet_s'])
        prev_dim = row['dim_singlet']

        # Print row
        E_2_3S = row['E_2_3S_Ha']
        E_2_1S = row['E_2_1S_Ha']
        dE_cm = row['dE_st_cm_inverse']
        err_NIST = row['err_vs_NIST_pct']

        e_2_3s_str = f"{E_2_3S:>+13.6f}" if E_2_3S is not None else "         N/A "
        e_2_1s_str = f"{E_2_1S:>+13.6f}" if E_2_1S is not None else "         N/A "
        de_str = f"{dE_cm:>+14.2f}" if dE_cm is not None else "          N/A "
        err_str = f"{err_NIST:>+13.3f}%" if err_NIST is not None else "          N/A "

        print(f"{n_max:>6d} {row['dim_singlet']:>7d} {row['dim_triplet']:>7d} "
              f"{prev_t_total:>12.2f}  {e_2_3s_str} {e_2_1s_str} "
              f"{de_str} {err_str}")

        rows.append(row)

        # Early stop if convergence has clearly plateaued (relative change < 0.05%
        # over two consecutive steps — tighter than the 0.2% per-step threshold
        # so we record at least 3 plateau values).
        if len(rows) >= 4:
            last = rows[-1]['dE_st_cm_inverse']
            prev = rows[-2]['dE_st_cm_inverse']
            prev2 = rows[-3]['dE_st_cm_inverse']
            if all(x is not None for x in (last, prev, prev2)) and abs(prev) > 1e-10:
                rel_change_1 = abs(last - prev) / abs(prev) * 100.0
                rel_change_2 = abs(prev - prev2) / abs(prev2) * 100.0
                if rel_change_1 < 0.05 and rel_change_2 < 0.05:
                    print(f"\n[CONVERGED] Relative change at n_max={n_max} is {rel_change_1:.4f}% (and prev step {rel_change_2:.4f}%); stopping.")
                    break

    # --- Diagnostics ---
    print()
    print("--- Validation checks at largest n_max ---")
    final = rows[-1]
    print(f"  variational(1^1S):    E={final['E_1_1S_Ha']:+.6f} > exact={DRAKE_E_1S_HA:+.6f}: "
          f"{'OK' if final['variational_1_1S_ok'] else 'VIOLATED'}")
    print(f"  variational(2^3S):    E={final['E_2_3S_Ha']:+.6f} > exact={DRAKE_E_2_3S_HA:+.6f}: "
          f"{'OK' if final['variational_2_3S_ok'] else 'VIOLATED'}")
    print(f"  variational(2^1S):    E={final['E_2_1S_Ha']:+.6f} > exact={DRAKE_E_2_1S_HA:+.6f}: "
          f"{'OK' if final['variational_2_1S_ok'] else 'VIOLATED'}")
    print(f"  Hund's rule (E(2^1S) > E(2^3S)):  {'OK' if final['ordering_ok'] else 'VIOLATED'}")
    print(f"  Excited above ground state:       {'OK' if final['excited_above_ground_ok'] else 'VIOLATED'}")
    print(f"  Below He+ ionization (E < -2 Ha): {'OK' if final['below_ionization_ok'] else 'VIOLATED'}")

    # --- Compare splitting accuracy with absolute energy accuracy ---
    print()
    print("--- Splitting vs absolute-energy accuracy ---")
    for row in rows:
        nm = row['n_max']
        if row['err_vs_NIST_pct'] is None:
            continue
        err_split = abs(row['err_vs_NIST_pct'])
        err_gs = row['err_E_1_1S_pct']
        err_2_3S = row['err_E_2_3S_pct']
        err_2_1S = row['err_E_2_1S_pct']
        print(f"  n_max={nm}: |err_dE_st| = {err_split:>6.3f}%, "
              f"|err_E(1^1S)| = {err_gs:>6.3f}%, "
              f"|err_E(2^3S)| = {err_2_3S:>6.3f}%, "
              f"|err_E(2^1S)| = {err_2_1S:>6.3f}%")

    # Common-mode cancellation analysis: ratio of splitting error to typical
    # absolute error
    if len(rows) >= 2:
        last = rows[-1]
        if (last['err_vs_NIST_pct'] is not None and last['err_E_1_1S_pct'] is not None
                and last['err_E_1_1S_pct'] > 0):
            ratio = abs(last['err_vs_NIST_pct']) / last['err_E_1_1S_pct']
            print()
            print(f"  Splitting-error / abs-energy-error ratio at n_max={last['n_max']}: "
                  f"{ratio:.2f}")
            if ratio > 1.0:
                print("  Splitting is LESS accurate than absolute (no common-mode cancellation).")
            else:
                print("  Splitting is MORE accurate than absolute (common-mode cancellation observed).")

    # --- Save JSON ---
    out_path = PROJECT_ROOT / "debug" / "data" / "precision_catalogue_he_2s_singlet_triplet.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    save: Dict[str, Any] = {
        'system': 'He 2^1S_0 - 2^3S_1 exchange splitting (J-averaged)',
        'architecture': 'graph-native CI (Track DI Sprint 3C); hybrid h1 (exact diagonal + kappa graph adjacency); analytical Slater V_ee at k=Z=2',
        'method': 'build_graph_native_fci with explicit spin sectors',
        'reference_NIST_cm_inverse': NIST_SPLITTING_CM_INVERSE,
        'reference_NIST_THz': NIST_SPLITTING_CM_INVERSE * 29.9792458 * 1e-3,
        'reference_NIST_Ha': NIST_SPLITTING_CM_INVERSE / HA_TO_CM_INVERSE,
        'reference_Drake_NR_Ha': DRAKE_SPLITTING_HA,
        'reference_Drake_NR_cm_inverse': DRAKE_SPLITTING_CM_INVERSE,
        'reference_Drake_E_1_1S_Ha': DRAKE_E_1S_HA,
        'reference_Drake_E_2_3S_Ha': DRAKE_E_2_3S_HA,
        'reference_Drake_E_2_1S_Ha': DRAKE_E_2_1S_HA,
        'note': (
            'Hund first rule: E(2^1S) > E(2^3S); singlet higher than triplet '
            'because the symmetric spatial wavefunction (singlet) places electrons '
            'closer together, costing more Coulomb repulsion than the antisymmetric '
            'triplet. Splitting magnitude = 2 * K(1s, 2s) at leading order, where '
            'K = G^0(1s, 2s) is the exchange Slater integral.'
        ),
        'rows': rows,
    }
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(save, f, indent=2)
    print(f"\nSaved JSON to {out_path}")

    # --- Verdict summary ---
    print()
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    final = rows[-1]
    nm = final['n_max']
    de_cm = final['dE_st_cm_inverse']
    err_pct = final['err_vs_NIST_pct']
    print()
    print(f"At largest tractable n_max = {nm}:")
    print(f"  dE(2^1S - 2^3S)  = {de_cm:.2f} cm^-1")
    print(f"  vs NIST          = {NIST_SPLITTING_CM_INVERSE:.4f} cm^-1")
    print(f"  Residual         = {err_pct:+.3f}% ({final['err_vs_NIST_cm_inverse']:+.2f} cm^-1)")
    print()
    print(f"Variational bound respected for all three states: "
          f"{all([final['variational_1_1S_ok'], final['variational_2_3S_ok'], final['variational_2_1S_ok']])}")
    print(f"Hund's rule (E(2^1S) > E(2^3S)) verified:           {final['ordering_ok']}")
    print(f"Both states below He+ ionization threshold:         {final['below_ionization_ok']}")
    print()

    return save


if __name__ == "__main__":
    main()
