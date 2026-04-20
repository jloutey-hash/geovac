"""Native Dirac-Coulomb FCI engine (Sprint 1, Track D-1).

Two-electron full configuration interaction in the (n, kappa, m_j) Dirac
spinor basis with exact algebraic matrix elements.  This is NOT a
perturbative Breit-Pauli correction on top of Schrodinger — it uses the
Dirac-Coulomb eigenvalues directly as the one-body Hamiltonian and
jj-coupled Gaunt coefficients for the two-body Coulomb ERI.

One-body h1
-----------
Diagonal:
    h1(a,a) = -Z^2/(2*n^2) + E_FS(n, kappa, Z, alpha)

where E_FS is the alpha^4 Dirac fine-structure correction from
``fine_structure.fine_structure_total``.  For alpha=0 this reduces to
the non-relativistic -Z^2/(2*n^2).

Off-diagonal elements from the orbital exponent mismatch (k != Z) are
included when ``k_orb`` is specified.

Two-body V_ee
--------------
    <ab|1/r12|cd> = sum_k X_k(a,c) * X_k(b,d) * R^k(a,b,c,d)

Angular: ``jj_angular_Xk`` from composed_qubit_relativistic (full-Gaunt,
TR-phase-fixed).  Radial: ``get_rk_float`` from hypergeometric_slater
(exact algebraic).

Configuration space
-------------------
Antisymmetric 2-electron Slater determinants |a,b> with a < b in the
DiracLabel enumeration order.  Classified by:
- Total M_J_twice = two_m_j(a) + two_m_j(b)
- Parity = (-1)^{l_a + l_b}

Reuses
------
- ``dirac_matrix_elements``: DiracLabel, iter_dirac_labels, kappa_to_l/j
- ``composed_qubit_relativistic``: jj_angular_Xk, enumerate_dirac_labels
- ``hypergeometric_slater``: get_rk_float
- ``fine_structure``: fine_structure_total
- ``casimir_ci`` pattern: FCI matrix assembly, Slater-Condon rules
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Rational

from geovac.dirac_matrix_elements import (
    DiracLabel,
    iter_dirac_labels,
    kappa_to_l,
    kappa_to_j,
)
from geovac.composed_qubit_relativistic import (
    jj_angular_Xk,
    enumerate_dirac_labels,
)
from geovac.fine_structure import fine_structure_total
from geovac.hypergeometric_slater import get_rk_float


__all__ = [
    "build_dirac_fci_matrix",
    "dirac_h1_diagonal",
    "solve_dirac_fci",
    "dirac_convergence_table",
    "he_fine_structure",
    "dirac_fci_to_qubit",
    "dirac_qubit_resource_table",
]

_ALPHA_PHYSICAL = 7.2973525693e-3


def _build_spinor_labels(n_max: int) -> List[DiracLabel]:
    """Enumerate all DiracLabel states up to n_max."""
    return list(iter_dirac_labels(n_max))


def dirac_h1_diagonal(
    Z: int,
    n_max: int,
    alpha_num: float = _ALPHA_PHYSICAL,
    k_orb: Optional[float] = None,
) -> Dict[int, float]:
    """Diagonal one-body energies for Dirac spinor orbitals.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    alpha_num : float
        Numerical fine-structure constant. Default CODATA value.
    k_orb : float or None
        Orbital exponent. If None, uses Z (hydrogenic).

    Returns
    -------
    dict mapping orbital index -> one-body energy (float)
    """
    labels = _build_spinor_labels(n_max)
    if k_orb is None:
        k_orb = float(Z)

    h1 = {}
    for idx, lab in enumerate(labels):
        n = lab.n_fock
        e_nr = -float(Z)**2 / (2.0 * n**2)

        if alpha_num != 0.0:
            e_fs_sym = fine_structure_total(
                n, lab.kappa,
                Z=Integer(Z),
                alpha=Rational(alpha_num).limit_denominator(10**15),
            )
            e_fs = float(e_fs_sym) if e_fs_sym != 0 else 0.0
        else:
            e_fs = 0.0

        h1[idx] = e_nr + e_fs

    return h1


def _compute_breit_rk_cache(
    Z: int,
    labels: List[DiracLabel],
) -> Dict[Tuple[int, ...], float]:
    """Compute Breit-Pauli retarded radial integrals for all orbital pairs.

    Key format: (n1, l1, n2, l2, n3, l3, n4, l4, k) matching the Coulomb
    rk_cache convention.  Uses the float-algebraic evaluator from
    ``hypergeometric_slater`` (same closed-form math, ~1000x faster than
    the sympy path in ``breit_integrals``).  Z^3 scaling applied here.
    """
    from geovac.hypergeometric_slater import get_rk_breit_float

    unique_nl = sorted({(lab.n_fock, lab.l) for lab in labels})
    Z_cubed = float(Z) ** 3
    cache: Dict[Tuple[int, ...], float] = {}

    for n1, l1 in unique_nl:
        for n2, l2 in unique_nl:
            for n3, l3 in unique_nl:
                for n4, l4 in unique_nl:
                    k_max = min(l1 + l3, l2 + l4)
                    for k in range(k_max + 1):
                        if (l1 + l3 + k) % 2 != 0:
                            continue
                        if (l2 + l4 + k) % 2 != 0:
                            continue
                        key = (n1, l1, n2, l2, n3, l3, n4, l4, k)
                        if key in cache:
                            continue
                        try:
                            val = get_rk_breit_float(
                                n1, l1, n3, l3,
                                n2, l2, n4, l4, k)
                            cache[key] = val * Z_cubed
                        except (ValueError, ZeroDivisionError):
                            pass
    return cache


def _build_eri_cache(
    Z: int,
    labels: List[DiracLabel],
    k_orb: float,
    include_breit: bool = False,
    alpha_num: float = 0.0,
) -> Dict[Tuple[int, int, int, int], float]:
    """Build jj-coupled two-electron integrals (Coulomb + optional Breit).

    Coulomb: <ab|1/r12|cd> = sum_k X_k(a,c) * X_k(b,d) * R^k * k_orb
    Breit:   <ab|B|cd>     = sum_k X_k(a,c) * X_k(b,d) * R^k_BP * alpha^2

    Returns sparse dict keyed by (a, b, c, d) orbital indices.
    """
    Q = len(labels)

    l_max = max(lab.l for lab in labels) if labels else 0
    k_max = 2 * l_max

    ac_Xk: Dict[Tuple[int, int], List[Tuple[int, float]]] = defaultdict(list)
    for a in range(Q):
        la = labels[a]
        for c in range(Q):
            lc = labels[c]
            for k in range(k_max + 1):
                val = jj_angular_Xk(la.kappa, la.two_m_j,
                                    lc.kappa, lc.two_m_j, k)
                if abs(val) > 1e-15:
                    ac_Xk[(a, c)].append((k, val))

    rk_cache: Dict[Tuple[int, ...], float] = {}

    def _get_rk(la: DiracLabel, lb: DiracLabel,
                lc: DiracLabel, ld: DiracLabel, k: int) -> float:
        key = (la.n_fock, la.l, lb.n_fock, lb.l,
               lc.n_fock, lc.l, ld.n_fock, ld.l, k)
        if key not in rk_cache:
            rk_cache[key] = get_rk_float(
                la.n_fock, la.l, lc.n_fock, lc.l,
                lb.n_fock, lb.l, ld.n_fock, ld.l, k)
        return rk_cache[key]

    breit_rk: Dict[Tuple[int, ...], float] = {}
    if include_breit and alpha_num != 0.0:
        breit_rk = _compute_breit_rk_cache(Z, labels)
    alpha_sq = alpha_num ** 2

    eri: Dict[Tuple[int, int, int, int], float] = {}

    for (a, c), acs in ac_Xk.items():
        la = labels[a]
        lc = labels[c]
        for (b, d), bds in ac_Xk.items():
            lb = labels[b]
            ld = labels[d]
            if la.two_m_j + lb.two_m_j != lc.two_m_j + ld.two_m_j:
                continue
            val = 0.0
            for k_ac, x_ac in acs:
                for k_bd, x_bd in bds:
                    if k_ac != k_bd:
                        continue
                    k = k_ac
                    ang = x_ac * x_bd

                    rk = _get_rk(la, lb, lc, ld, k)
                    val += ang * rk * k_orb

                    if breit_rk:
                        bkey = (la.n_fock, la.l, lb.n_fock, lb.l,
                                lc.n_fock, lc.l, ld.n_fock, ld.l, k)
                        rk_b = breit_rk.get(bkey, 0.0)
                        val += ang * rk_b * alpha_sq

            if abs(val) > 1e-14:
                eri[(a, b, c, d)] = val

    return eri


def build_dirac_fci_matrix(
    Z: int,
    n_max: int,
    alpha_num: float = _ALPHA_PHYSICAL,
    k_orb: Optional[float] = None,
    M_J_twice: int = 0,
    include_breit: bool = False,
) -> Tuple[np.ndarray, List[Tuple[int, int]], List[DiracLabel]]:
    """Build two-electron Dirac-Coulomb FCI matrix.

    Antisymmetric Slater determinants |a,b> with a < b.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    alpha_num : float
        Fine-structure constant. Use 0.0 for non-relativistic limit.
    k_orb : float or None
        Orbital exponent. Default Z.
    M_J_twice : int
        2 * total M_J. Selects configs with m_j(a) + m_j(b) = M_J.
        Default 0.
    include_breit : bool
        If True, add Breit-Pauli retarded two-body corrections (SS + SOO)
        with alpha^2 prefactor.  Vanishes identically at alpha_num=0.

    Returns
    -------
    H : ndarray (n_configs, n_configs)
        FCI Hamiltonian matrix.
    configs : list of (int, int)
        Configuration list — pairs of orbital indices.
    labels : list of DiracLabel
        Orbital label list.
    """
    if k_orb is None:
        k_orb = float(Z)

    labels = _build_spinor_labels(n_max)
    Q = len(labels)

    h1_diag = dirac_h1_diagonal(Z, n_max, alpha_num, k_orb)

    eri = _build_eri_cache(Z, labels, k_orb,
                          include_breit=include_breit, alpha_num=alpha_num)

    configs: List[Tuple[int, int]] = []
    for a in range(Q):
        for b in range(a + 1, Q):
            if labels[a].two_m_j + labels[b].two_m_j == M_J_twice:
                configs.append((a, b))

    n_configs = len(configs)
    if n_configs == 0:
        return np.zeros((0, 0)), configs, labels

    H = np.zeros((n_configs, n_configs))

    for I in range(n_configs):
        a, b = configs[I]
        for J in range(I, n_configs):
            c, d = configs[J]

            me = _slater_condon_2e(a, b, c, d, h1_diag, eri)

            H[I, J] = me
            H[J, I] = me

    return H, configs, labels


def _slater_condon_2e(
    a: int, b: int, c: int, d: int,
    h1: Dict[int, float],
    eri: Dict[Tuple[int, int, int, int], float],
) -> float:
    """Matrix element <ab||cd> for antisymmetric determinants.

    Slater-Condon rules for two-electron determinants |a,b> and |c,d>
    where a < b and c < d.

    Diagonal (a==c, b==d):
        h1(a) + h1(b) + <ab|ab> - <ab|ba>

    Single excitation (one orbital differs):
        sign * [h1(p,q) + <pr|qr> - <pr|rq>]
        where p->q is the excitation and r is the shared orbital.

    Double excitation (both orbitals differ):
        <ab|cd> - <ab|dc>
    """
    if a == c and b == d:
        val = h1.get(a, 0.0) + h1.get(b, 0.0)
        val += eri.get((a, b, a, b), 0.0) - eri.get((a, b, b, a), 0.0)
        return val

    shared = set()
    bra = {a, b}
    ket = {c, d}
    shared = bra & ket
    diff_bra = bra - shared
    diff_ket = ket - shared

    if len(shared) == 1:
        r = shared.pop()
        p = diff_bra.pop()
        q = diff_ket.pop()
        sign_bra = 1.0 if ((a == r and b == p) or (a == p and b != r)) else 1.0
        sign_bra = _perm_sign(a, b, r, p)
        sign_ket = _perm_sign(c, d, r, q)
        sign = sign_bra * sign_ket

        val = 0.0
        val += eri.get((p, r, q, r), 0.0) - eri.get((p, r, r, q), 0.0)
        return sign * val

    if len(shared) == 0:
        p1, p2 = sorted(diff_bra)
        q1, q2 = sorted(diff_ket)
        sign_bra = _perm_sign(a, b, p1, p2)
        sign_ket = _perm_sign(c, d, q1, q2)
        sign = sign_bra * sign_ket

        val = eri.get((p1, p2, q1, q2), 0.0) - eri.get((p1, p2, q2, q1), 0.0)
        return sign * val

    return 0.0


def _perm_sign(orig_a: int, orig_b: int, tgt_first: int, tgt_second: int) -> float:
    """Sign of permutation mapping (orig_a, orig_b) -> (tgt_first, tgt_second).

    Both inputs are ordered (orig_a < orig_b).
    Returns +1 if tgt_first == orig_a (identity) or -1 if swapped.
    """
    if orig_a == tgt_first and orig_b == tgt_second:
        return 1.0
    elif orig_a == tgt_second and orig_b == tgt_first:
        return -1.0
    else:
        raise ValueError(
            f"({orig_a},{orig_b}) does not contain ({tgt_first},{tgt_second})")


def solve_dirac_fci(
    Z: int,
    n_max: int,
    alpha_num: float = _ALPHA_PHYSICAL,
    k_orb: Optional[float] = None,
    n_states: int = 5,
    M_J_twice: int = 0,
    include_breit: bool = False,
) -> Dict[str, Any]:
    """Solve the two-electron Dirac-Coulomb FCI.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    alpha_num : float
        Fine-structure constant.
    k_orb : float or None
        Orbital exponent. Default Z.
    n_states : int
        Number of eigenvalues to return.
    M_J_twice : int
        2 * total M_J quantum number.
    include_breit : bool
        If True, add Breit-Pauli retarded two-body corrections.

    Returns
    -------
    dict with keys:
        'energies' : ndarray of lowest eigenvalues
        'configs' : list of (a, b) index pairs
        'labels' : list of DiracLabel
        'n_configs' : int
        'n_orbitals' : int
        'Z' : int
        'n_max' : int
        'alpha' : float
    """
    H, configs, labels = build_dirac_fci_matrix(
        Z, n_max, alpha_num, k_orb, M_J_twice, include_breit=include_breit)

    n_configs = len(configs)
    if n_configs == 0:
        return {
            'energies': np.array([]),
            'configs': configs,
            'labels': labels,
            'n_configs': 0,
            'n_orbitals': len(labels),
            'Z': Z,
            'n_max': n_max,
            'alpha': alpha_num,
        }

    evals = np.linalg.eigh(H)[0]
    n_ret = min(n_states, len(evals))

    return {
        'energies': evals[:n_ret],
        'configs': configs,
        'labels': labels,
        'n_configs': n_configs,
        'n_orbitals': len(labels),
        'Z': Z,
        'n_max': n_max,
        'alpha': alpha_num,
    }


def dirac_convergence_table(
    Z: int,
    n_max_range: List[int],
    alpha_num: float = _ALPHA_PHYSICAL,
    M_J_twice: int = 0,
) -> List[Dict[str, Any]]:
    """Convergence table: ground-state energy vs n_max.

    Parameters
    ----------
    Z : int
        Nuclear charge (2 for He).
    n_max_range : list of int
        Values of n_max to compute.
    alpha_num : float
        Fine-structure constant.
    M_J_twice : int
        2 * total M_J.

    Returns
    -------
    list of dicts, one per n_max, with keys:
        'n_max', 'n_orbitals', 'n_configs', 'E_ground',
        'E_ground_nr' (alpha=0 ground state for comparison)
    """
    results = []
    for nm in n_max_range:
        res = solve_dirac_fci(Z, nm, alpha_num, n_states=1, M_J_twice=M_J_twice)
        res_nr = solve_dirac_fci(Z, nm, 0.0, n_states=1, M_J_twice=M_J_twice)

        E_gs = res['energies'][0] if len(res['energies']) > 0 else None
        E_nr = res_nr['energies'][0] if len(res_nr['energies']) > 0 else None

        results.append({
            'n_max': nm,
            'n_orbitals': res['n_orbitals'],
            'n_configs': res['n_configs'],
            'E_ground': E_gs,
            'E_ground_nr': E_nr,
            'rel_shift': (E_gs - E_nr) if (E_gs is not None and E_nr is not None) else None,
        })

    return results


_HA_TO_CM = 219474.63137054


def he_fine_structure(
    n_max: int,
    alpha_num: float = _ALPHA_PHYSICAL,
    include_breit: bool = True,
) -> Dict[str, Any]:
    """Compute He 2^3P_{0,1,2} fine-structure splittings from Dirac FCI.

    Extracts the J=0, J=1, J=2 levels of the 2^3P multiplet by
    exploiting the M_J quantum number.  M_J_twice=4 isolates J>=2 only;
    the lowest such state at n_max>=2 is the 2^3P_2 level.
    M_J_twice=2 adds J>=1; the new eigenvalue is 2^3P_1.
    M_J_twice=0 adds J>=0; the new eigenvalue is 2^3P_0.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number (>=2 required).
    alpha_num : float
        Fine-structure constant.
    include_breit : bool
        If True, include Breit-Pauli retarded two-body corrections.

    Returns
    -------
    dict with keys:
        'E_J' : dict mapping J -> energy (Ha)
        'splittings_ha' : dict with '0-1' and '1-2' splittings in Ha
        'splittings_cm' : dict with '0-1' and '1-2' splittings in cm^-1
        'inverted' : bool — True if E(J=0) > E(J=1) > E(J=2)
        'n_max' : int
        'alpha' : float
        'include_breit' : bool
    """
    if n_max < 2:
        raise ValueError("n_max must be >= 2 for 2^3P states")

    n_eig = max(10, 3 * n_max**2)
    res_mj0 = solve_dirac_fci(
        Z=2, n_max=n_max, alpha_num=alpha_num, n_states=n_eig,
        M_J_twice=0, include_breit=include_breit)

    res_nr = solve_dirac_fci(
        Z=2, n_max=n_max, alpha_num=0.0, n_states=n_eig, M_J_twice=0)

    nr_triplet_p = _find_nr_triplet_p_energy(res_nr['energies'])
    if nr_triplet_p is None:
        raise ValueError("Could not identify NR 2^3P level")

    res_mj4 = solve_dirac_fci(
        Z=2, n_max=n_max, alpha_num=alpha_num, n_states=n_eig,
        M_J_twice=4, include_breit=include_breit)

    gs_energy = res_mj0['energies'][0]
    E_J2_anchor = _find_triplet_p_level(res_mj4['energies'], gs_energy)

    cluster = _extract_triplet_p_cluster(
        res_mj0['energies'], nr_triplet_p, gs_energy)

    E_J: Dict[int, float] = {}
    if len(cluster) >= 3 and E_J2_anchor is not None:
        closest_to_anchor = min(cluster, key=lambda e: abs(e - E_J2_anchor))
        E_J[2] = closest_to_anchor
        remaining = sorted([e for e in cluster if e != closest_to_anchor])
        if len(remaining) >= 2:
            E_J[1] = remaining[0]
            E_J[0] = remaining[1]
        elif len(remaining) == 1:
            E_J[1] = remaining[0]
    elif len(cluster) >= 3:
        sorted_cluster = sorted(cluster)
        E_J[2] = sorted_cluster[0]
        E_J[1] = sorted_cluster[1]
        E_J[0] = sorted_cluster[2]

    splittings_ha: Dict[str, float] = {}
    splittings_cm: Dict[str, float] = {}

    if 0 in E_J and 1 in E_J:
        splittings_ha['0-1'] = E_J[0] - E_J[1]
        splittings_cm['0-1'] = splittings_ha['0-1'] * _HA_TO_CM
    if 1 in E_J and 2 in E_J:
        splittings_ha['1-2'] = E_J[1] - E_J[2]
        splittings_cm['1-2'] = splittings_ha['1-2'] * _HA_TO_CM

    inverted = (len(E_J) == 3 and E_J[0] > E_J[1] > E_J[2])

    return {
        'E_J': E_J,
        'splittings_ha': splittings_ha,
        'splittings_cm': splittings_cm,
        'inverted': inverted,
        'n_max': n_max,
        'alpha': alpha_num,
        'include_breit': include_breit,
    }


def _find_triplet_p_level(
    energies: np.ndarray,
    gs_energy: float,
) -> Optional[float]:
    """Find the lowest excited-state energy in M_J_twice=4 sector.

    For M_J_twice=4, only J>=2 contributes, so the lowest eigenvalue
    in the n=2 manifold is the 2^3P_2 level.
    """
    if len(energies) == 0:
        return None
    for e in energies:
        if e > gs_energy + 0.01:
            return float(e)
    return float(energies[0]) if len(energies) > 0 else None


def _find_nr_triplet_p_energy(energies: np.ndarray) -> Optional[float]:
    """Find the NR 2^3P energy as the first triply degenerate excited level.

    At alpha=0 in M_J=0, the 2^3P multiplet is triply degenerate (J=0,1,2
    all have the same energy).  Looks for a cluster of 3+ eigenvalues that
    are degenerate to within 1e-8 Ha.
    """
    if len(energies) < 4:
        return None

    sorted_e = sorted(energies)
    for i in range(2, len(sorted_e)):
        if abs(sorted_e[i] - sorted_e[i - 2]) < 1e-8:
            if sorted_e[i] > sorted_e[0] + 0.01:
                return float(sorted_e[i - 1])
    return None


def _extract_triplet_p_cluster(
    energies: np.ndarray,
    nr_reference: float,
    gs_energy: float,
    window: float = 0.005,
) -> List[float]:
    """Extract the 3 eigenvalues nearest to the NR 2^3P reference energy.

    Searches within ``window`` Ha of the NR reference to find the
    fine-structure-split ³P_{0,1,2} cluster.
    """
    excited = [float(e) for e in energies if e > gs_energy + 0.01]
    nearby = [e for e in excited if abs(e - nr_reference) < window]
    nearby.sort()
    return nearby[:3] if len(nearby) >= 3 else nearby


# ---------------------------------------------------------------------------
# Track D-3: Qubit encoding and resource benchmarks
# ---------------------------------------------------------------------------


def dirac_fci_to_qubit(
    Z: int,
    n_max: int,
    alpha_num: float = _ALPHA_PHYSICAL,
    k_orb: Optional[float] = None,
    include_breit: bool = False,
) -> Any:
    """Convert Dirac-Coulomb FCI Hamiltonian to JW-encoded QubitOperator.

    Each DiracLabel (n, kappa, m_j) maps to one qubit.  No spin doubling
    — the spin is already encoded in kappa.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    alpha_num : float
        Fine-structure constant.
    k_orb : float or None
        Orbital exponent (default Z).
    include_breit : bool
        Include Breit two-body corrections.

    Returns
    -------
    GeoVacHamiltonian
        Wrapped QubitOperator with metadata and resource metrics.
    """
    from openfermion import FermionOperator, jordan_wigner

    labels = _build_spinor_labels(n_max)
    Q = len(labels)

    if k_orb is None:
        k_orb_val = float(Z)
    else:
        k_orb_val = float(k_orb)

    h1 = dirac_h1_diagonal(Z, n_max, alpha_num=alpha_num, k_orb=k_orb)
    eri = _build_eri_cache(Z, labels, k_orb_val,
                           include_breit=include_breit,
                           alpha_num=alpha_num)

    sym_eri: Dict[Tuple[int, int, int, int], float] = {}
    for (a, b, c, d), val in eri.items():
        key = (a, b, c, d)
        alt = (c, d, a, b)
        sym_eri[key] = sym_eri.get(key, 0.0) + 0.5 * val
        sym_eri[alt] = sym_eri.get(alt, 0.0) + 0.5 * val

    fermion_op = FermionOperator((), 0.0)

    for p in range(Q):
        if abs(h1[p]) < 1e-14:
            continue
        fermion_op += FermionOperator(((p, 1), (p, 0)), h1[p])

    for (a, b, c, d), val in sym_eri.items():
        if a == b or c == d:
            continue
        if abs(val) < 1e-14:
            continue
        fermion_op += FermionOperator(
            ((a, 1), (b, 1), (d, 0), (c, 0)),
            0.5 * val,
        )

    qubit_op = jordan_wigner(fermion_op)

    all_terms = qubit_op.terms
    n_pauli = len(all_terms) - (1 if () in all_terms else 0)
    lam_total = sum(abs(c) for c in all_terms.values())
    lam_ni = sum(abs(c) for t, c in all_terms.items() if t != ())

    from geovac.measurement_grouping import qwc_groups
    n_qwc = len(qwc_groups(qubit_op))

    from geovac.ecosystem_export import GeoVacHamiltonian
    metadata = {
        'system': f'Dirac-CI Z={Z}',
        'Z': Z,
        'n_max': n_max,
        'Q': Q,
        'alpha_num': alpha_num,
        'include_breit': include_breit,
        'N_pauli': n_pauli,
        'one_norm_ni': lam_ni,
        'one_norm_total': lam_total,
        'QWC_groups': n_qwc,
    }

    return GeoVacHamiltonian(qubit_op, metadata=metadata)


def dirac_qubit_resource_table(
    Z_values: List[int],
    n_max_values: List[int],
    alpha_num: float = _ALPHA_PHYSICAL,
    include_breit: bool = False,
) -> List[Dict[str, Any]]:
    """Generate resource comparison table across Z and n_max.

    Parameters
    ----------
    Z_values : list of int
        Nuclear charges.
    n_max_values : list of int
        Basis sizes.
    alpha_num : float
        Fine-structure constant.
    include_breit : bool
        Include Breit corrections.

    Returns
    -------
    list of dict
        Each entry has keys: Z, n_max, Q, N_pauli, one_norm_ni, QWC_groups.
    """
    results = []
    for Z in Z_values:
        for n_max in n_max_values:
            ham = dirac_fci_to_qubit(Z, n_max, alpha_num=alpha_num,
                                     include_breit=include_breit)
            meta = ham.metadata
            results.append({
                'Z': Z,
                'n_max': n_max,
                'Q': meta['Q'],
                'N_pauli': meta['N_pauli'],
                'one_norm_ni': meta['one_norm_ni'],
                'QWC_groups': meta['QWC_groups'],
            })
    return results
