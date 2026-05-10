"""Calc Track He-Oscillator: Helium 2^1P -> 1^1S oscillator strength.

Paper 34 §V.C.5 NEW. Multi-electron extension of Sprint Calc-L (hydrogen
Lyman alpha at +0.055%) and Sprint Calc-P (hydrogen polarizability EXACT
at N_basis=2 in Sturmian basis).

Structural question
-------------------
Sprint Calc-P verified Sturmian's CONTINUUM-CLOSING property at exact
N_basis=2 for one-electron polarizability. The Dalgarno-Lewis function
(r^2 + 2r) e^(-r) lives in span{r e^(-r), r^2 e^(-r)} = span{S_2,1, S_3,1}
exactly, and the bound-state Sturmian basis absorbs the H continuum
structurally without continuum integration.

Does this property extend to multi-electron transitions? If yes, the
Sturmian basis becomes the natural basis for multi-electron transition
computations. If no, characterize where the closure breaks.

Architecture
------------
We compute the He 2^1P -> 1^1S oscillator strength via the LENGTH FORM:

    f = (2/3) omega |<1^1S | r_1 + r_2 | 2^1P>|^2

where omega = E(2^1P) - E(1^1S), and the matrix element is between
multi-electron antisymmetrized wavefunctions in the graph-native CI
sector at fixed M_L.

Two test paths
--------------

PATH A: Hydrogenic basis at fixed Z=Z_eff(1s)=2 (the "naive" baseline,
        same as Track P bound-state-only). Tests how the framework's
        STANDARD construction does on this transition.

PATH B: Sturmian-style basis at varying k_orb. Specifically, compare:
        - k_orb = Z = 2 (Coulomb Sturmian at Z, same exponent for 1s and 2p)
        - k_orb = 27/16 = 1.6875 (Slater 1s effective)
        - k_orb = 1 (the 2p-effective)
        and watch convergence with n_max.

The KEY DIAGNOSTIC: if the Sturmian closure property extends, we expect
rapid convergence (sub-percent by N=3) at the natural k_orb. If it fails,
the convergence will be either:
  (a) slow power-law (cusp problem in radial integrand of dipole)
  (b) wrong limit (multi-electron correlation not captured)

State labeling
--------------
Use subblock=(0, 0) for 1^1S and subblock=(0, 1) for 2^1P.

For 2^1P transition matrix elements in length form:
- 1^1S sits in (l1=0, l2=0, M_L=0) singlet subblock; lowest eigenvalue.
- 2^1P sits in (l1=0, l2=1, M_L=0) singlet subblock; lowest eigenvalue.
  At M_L=0, only the m=0 sublevel of the 2p shell contributes.

The dipole operator z = z_1 + z_2 connects these because z is rank-1
(angular Y_1^0). Selection rule: Delta L = +/-1 (Wigner 3j requires
|L_i - L_f| <= 1 <= L_i + L_f); here 0 -> 1, allowed. Spin: singlet
-> singlet (z is spin-independent).

Transition matrix element
-------------------------
For two-electron singlet-to-singlet:

  <Psi(1^1S) | z_1 + z_2 | Psi(2^1P)>

where each Psi is a CI expansion over Slater determinants in the
spatial-orbital basis. We expand:

  Psi(1^1S) = sum_{ij} c_ij^{1S} |ij,1S>
  Psi(2^1P) = sum_{ij} c_ij^{1P} |ij,1P>

where |ij,1S> and |ij,1P> are spin-symmetric singlet spatial functions.
The dipole matrix element between determinants reduces to:

  <ij,1S | z_1 + z_2 | pq,1P> =
      (1/N_IJ N_PQ) * sum over (a,b,c,d) sign-products of
        <a|z|c> delta_{b,d} + <b|z|d> delta_{a,c}

i.e. a sum of one-particle dipole matrix elements, each tagged by an
orbital-overlap delta on the spectator electron.

Reference values
----------------
Drake handbook (Drake & Yan 1992 high-precision NR):
  f(2^1P -> 1^1S, He) = 0.27616  (length form)

Lewis, Hessel et al. 1989 experimental: f ~ 0.276 (consistent)

Earlier Theodosiou 1987 / Schiff & Pekeris: 0.27616 +/- 0.00003

These are all infinite-mass NR limits; relativistic and recoil
corrections are sub-percent at He.

Alternative Z_eff conventions to test
-------------------------------------
- Z_eff(1s) = 27/16 = 1.6875 (variational, He ground-state minimum)
- Z_eff(1s) = 1.69 (Slater rules)
- Z_eff(1s) = 2 (bare nuclear charge, no screening)
- Z_eff(2p) = 1 (full screen by 1s for excited p-electron)

Date: 2026-05-09
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.casimir_ci import (
    HE_NR_REFERENCE,
    _build_graph_h1,
    _wigner3j,
    _gaunt_ck,
    two_electron_integral,
)
from geovac.dirac_matrix_elements import radial_matrix_element
from geovac.lattice import GeometricLattice
import sympy as sp


# Reference values
F_DRAKE_2P_1S: float = 0.27616  # Drake & Yan 1992 / Theodosiou 1987 / Schiff-Pekeris
F_DRAKE_LOWER: float = 0.27600
F_DRAKE_UPPER: float = 0.27632

HA_TO_EV: float = 27.21138386


# =====================================================================
# Single-particle dipole matrix elements <a|z|c>
# =====================================================================

def single_particle_dipole_z(
    n_a: int, l_a: int, m_a: int,
    n_c: int, l_c: int, m_c: int,
    Z_orb: float,
) -> float:
    """<n_a, l_a, m_a | z | n_c, l_c, m_c> in length form.

    z = sqrt(4 pi / 3) r Y_1^0 in spherical-tensor form.

    The matrix element factorizes as:
        <a|z|c> = c_1(l_a, m_a, l_c, m_c) * <R_{n_a, l_a} | r | R_{n_c, l_c}>

    where c_1(...) is a Gaunt coefficient (Wigner 3j combination).

    For real spherical harmonics with z direction, the angular factor is:
        (-1)^m_a sqrt((2l_a+1)(2l_c+1)) (l_a 1 l_c; 0 0 0)(l_a 1 l_c; -m_a 0 m_c)

    Selection rules:
        l_c = l_a +/- 1 (parity, Wigner 3j (l_a 1 l_c; 0 0 0) nonzero)
        m_c = m_a (z is m=0 photon)

    Parameters
    ----------
    n_a, l_a, m_a : int
        Initial orbital quantum numbers.
    n_c, l_c, m_c : int
        Final orbital quantum numbers.
    Z_orb : float
        Orbital exponent for hydrogenic radial functions (= Z for Coulomb
        Sturmian convention; = Z_eff for screened hydrogenic convention).

    Returns
    -------
    Matrix element value in atomic units (a_0 if Z_orb=1).
    """
    # Selection rules
    if m_a != m_c:
        return 0.0
    if abs(l_c - l_a) != 1:
        return 0.0

    # Angular factor (Gaunt c_k for k=1)
    c1 = _gaunt_ck(l_a, m_a, l_c, m_c, 1)
    if abs(c1) < 1e-15:
        return 0.0

    # Radial factor: <R_{n_a, l_a} | r | R_{n_c, l_c}> at orbital exponent Z_orb
    # radial_matrix_element returns symbolic; at Z=Z_orb gives the exact
    # rational/algebraic result.
    R_sym = radial_matrix_element(n_a, l_a, n_c, l_c, "r", Z=sp.Rational(1))
    # The function returns it for general Z; the Z scaling for <r> is 1/Z.
    # radial_matrix_element with Z=1 gives (...) * 1; for general Z it scales as 1/Z.
    R_at_Z = float(R_sym) / Z_orb

    return c1 * R_at_Z


# =====================================================================
# Multi-electron CI matrix in (l_1=0, l_2=L) sub-block
# =====================================================================

def build_singlet_LM_subblock(
    n_max: int,
    L_target: int,
    M_L_target: int,
    Z: int = 2,
    k_orb: float = None,
) -> Tuple[np.ndarray, List[Tuple[int, int]], List[Tuple[int, int, int]]]:
    """Build singlet sub-block FCI matrix at (l_1=0, l_2=L_target), M_L=M_L_target.

    Returns
    -------
    H : (n_configs, n_configs) FCI matrix
    configs : list of (i, j) orbital index pairs
    orbitals : list of (n, l, m) for index lookup
    """
    if k_orb is None:
        k_orb = float(Z)

    h1_mat, orbitals = _build_graph_h1(Z=Z, n_max=n_max)
    n_spatial = len(orbitals)

    # Indices for s-orbitals (l=0, m=0) and target-l-orbitals
    s_indices = [i for i, (n, l, m) in enumerate(orbitals)
                 if l == 0 and m == 0]
    target_indices = [i for i, (n, l, m) in enumerate(orbitals)
                      if l == L_target and m == M_L_target]

    if L_target == 0:
        # ss singlet: i <= j
        configs = [(i, j) for i in s_indices for j in s_indices if i <= j]
    else:
        # sl singlet (l != 0): i in s, j in target
        configs = [(i, j) for i in s_indices for j in target_indices]

    n_configs = len(configs)
    if n_configs == 0:
        return np.zeros((0, 0)), [], orbitals

    H = np.zeros((n_configs, n_configs))

    def g_int(a, b, c, d):
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        return two_electron_integral(na, la, ma, nb, lb, mb,
                                     nc, lc, mc, nd, ld, md, k_orb)

    parity = 1.0  # singlet

    for I in range(n_configs):
        i, j = configs[I]
        for J in range(I, n_configs):
            p, q = configs[J]
            bra_perms = [(i, j, 1.0)]
            if i != j:
                bra_perms.append((j, i, parity))
            ket_perms = [(p, q, 1.0)]
            if p != q:
                ket_perms.append((q, p, parity))
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
    return H, configs, orbitals


# =====================================================================
# Multi-electron transition matrix element
# =====================================================================

def transition_dipole(
    psi_initial: np.ndarray,
    configs_initial: List[Tuple[int, int]],
    psi_final: np.ndarray,
    configs_final: List[Tuple[int, int]],
    orbitals: List[Tuple[int, int, int]],
    Z_orb: float,
    parity_initial: float = 1.0,
    parity_final: float = 1.0,
) -> float:
    """<Psi_initial | z_1 + z_2 | Psi_final> for two-electron singlet states.

    Each Psi is a CI expansion sum_I c_I |I,singlet> where
    |I,singlet> = (|i_I j_I> + |j_I i_I>) / sqrt(2) for i != j and
    |I,singlet> = |i_I i_I> for i == j.

    The dipole operator is one-body. By Slater-Condon rules:

        <I,1S | z_1 + z_2 | J,1P> =
          (1/N_I N_J) sum_{(a,b)<-I, (c,d)<-J}
            sign * [<a|z|c> delta_{b,d} + <b|z|d> delta_{a,c}]

    where the sum is over the (1 or 2) permutations on each side.

    Returns
    -------
    Matrix element <Psi_init | r_1+r_2 (z) | Psi_final>.

    Note: r in length form means r-cosθ which is the z-component.
    """
    # Cache single-particle dipoles
    dipole_cache: Dict[Tuple[int, int, int, int, int, int], float] = {}

    def z_op(a: int, c: int) -> float:
        """<orbital a | z | orbital c>"""
        na, la, ma = orbitals[a]
        nc, lc, mc = orbitals[c]
        key = (na, la, ma, nc, lc, mc)
        if key not in dipole_cache:
            dipole_cache[key] = single_particle_dipole_z(
                na, la, ma, nc, lc, mc, Z_orb=Z_orb,
            )
        return dipole_cache[key]

    matrix_element = 0.0

    for I, c_I in enumerate(psi_initial):
        if abs(c_I) < 1e-15:
            continue
        i, j = configs_initial[I]
        bra_perms = [(i, j, 1.0)]
        if i != j:
            bra_perms.append((j, i, parity_initial))
        N_I = np.sqrt(float(len(bra_perms)))

        for J, c_J in enumerate(psi_final):
            if abs(c_J) < 1e-15:
                continue
            p, q = configs_final[J]
            ket_perms = [(p, q, 1.0)]
            if p != q:
                ket_perms.append((q, p, parity_final))
            N_J = np.sqrt(float(len(ket_perms)))

            me_IJ = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    # <a|z|c> delta_{b,d}
                    if b == d:
                        me_IJ += sign * z_op(a, c)
                    # <b|z|d> delta_{a,c}
                    if a == c:
                        me_IJ += sign * z_op(b, d)

            me_IJ /= (N_I * N_J)
            matrix_element += c_I * c_J * me_IJ

    return matrix_element


# =====================================================================
# Oscillator strength assembly
# =====================================================================

def compute_oscillator_strength(
    n_max: int,
    Z: int = 2,
    k_orb: float = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Compute He 2^1P -> 1^1S oscillator strength.

    f = (2/3) * omega * |<1^1S | r_1 + r_2 (z) | 2^1P>|^2

    where the matrix element is the z-component of the dipole; the (2/3)
    factor accounts for averaging over the three Cartesian components
    in length form for a transition between an L=0 initial state and
    L=1 final state. Specifically, in length form the standard formula is

        f = (2/3) * (E_f - E_i) * |<i | r | f>|^2  (atomic units)

    where |<i|r|f>|^2 = sum over all 2^1P sublevels (m=-1, 0, +1) of
    |<1^1S | r | 2^1P_m>|^2. By rotational invariance the three
    sublevels contribute equally, so

        |<i|r|f>|^2 = 3 * |<1^1S | z | 2^1P_{m=0}>|^2

    Hence
        f = 2 * (E_f - E_i) * |<1^1S | z | 2^1P_{m=0}>|^2

    Note: the conventional factor is (2/3) for a transition from a non-
    degenerate L=0 initial state. If we sum over the upper-state sublevels
    we get factor 2.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number (basis size parameter).
    Z : int
        Nuclear charge (default 2 for He).
    k_orb : float, optional
        Orbital exponent. If None, uses Z (Coulomb Sturmian convention).
    verbose : bool

    Returns
    -------
    Dict with: f_length, omega_Ha, dipole_z, E_1S_Ha, E_2P_Ha, n_max,
               n_configs_1S, n_configs_2P, k_orb, Z
    """
    if k_orb is None:
        k_orb = float(Z)

    if verbose:
        print(f"  n_max={n_max}, Z={Z}, k_orb={k_orb}")

    # Build 1^1S sub-block (l_1=0, l_2=0, M_L=0)
    H_1S, configs_1S, orbitals = build_singlet_LM_subblock(
        n_max=n_max, L_target=0, M_L_target=0, Z=Z, k_orb=k_orb,
    )
    if verbose:
        print(f"  1S configs: {len(configs_1S)}")
    eigvals_1S, eigvecs_1S = np.linalg.eigh(H_1S)
    psi_1S = eigvecs_1S[:, 0]
    E_1S = eigvals_1S[0]

    # Build 2^1P sub-block (l_1=0, l_2=1, M_L=0 for the m=0 sublevel)
    H_2P, configs_2P, _ = build_singlet_LM_subblock(
        n_max=n_max, L_target=1, M_L_target=0, Z=Z, k_orb=k_orb,
    )
    if verbose:
        print(f"  2P configs: {len(configs_2P)}")
    eigvals_2P, eigvecs_2P = np.linalg.eigh(H_2P)
    psi_2P = eigvecs_2P[:, 0]
    E_2P = eigvals_2P[0]

    omega = E_2P - E_1S

    # Compute transition dipole
    z_matrix_element = transition_dipole(
        psi_1S, configs_1S, psi_2P, configs_2P, orbitals,
        Z_orb=k_orb,
        parity_initial=1.0, parity_final=1.0,
    )

    # Oscillator strength: f = 2 * omega * |z_me|^2 (with sum over upper-state sublevels)
    f_length = 2.0 * omega * z_matrix_element ** 2

    return {
        "n_max": n_max,
        "Z": Z,
        "k_orb": k_orb,
        "n_configs_1S": len(configs_1S),
        "n_configs_2P": len(configs_2P),
        "E_1S_Ha": E_1S,
        "E_2P_Ha": E_2P,
        "omega_Ha": omega,
        "omega_eV": omega * HA_TO_EV,
        "dipole_z_au": z_matrix_element,
        "dipole_squared_au": z_matrix_element ** 2,
        "f_length": f_length,
    }


# =====================================================================
# Main convergence study
# =====================================================================

def main() -> Dict[str, Any]:
    print("=" * 70)
    print("He 2^1P -> 1^1S oscillator strength")
    print("Sturmian closure test for multi-electron transitions")
    print("=" * 70)

    results: Dict[str, Any] = {
        "reference": {
            "f_drake_yan_1992": F_DRAKE_2P_1S,
            "f_lower_uncertainty": F_DRAKE_LOWER,
            "f_upper_uncertainty": F_DRAKE_UPPER,
            "source": "Drake & Yan 1992 / Theodosiou 1987 / Schiff-Pekeris",
        },
        "convergence_study": {},
        "z_eff_study": {},
    }

    # Path A: standard Z_orb = Z = 2 (Coulomb Sturmian convention)
    print("\n" + "=" * 70)
    print("PATH A: Coulomb Sturmian at k_orb = Z = 2")
    print("=" * 70)
    print(f"  Drake & Yan reference: f = {F_DRAKE_2P_1S:.5f}")
    path_a_results = []
    for n_max in [2, 3, 4, 5]:
        t0 = time.time()
        try:
            r = compute_oscillator_strength(
                n_max=n_max, Z=2, k_orb=2.0, verbose=True,
            )
            elapsed = time.time() - t0
            r["wall_seconds"] = elapsed
            err_pct = (r["f_length"] - F_DRAKE_2P_1S) / F_DRAKE_2P_1S * 100.0
            r["error_pct"] = err_pct
            path_a_results.append(r)
            print(f"  n_max={n_max}: f = {r['f_length']:.6f}  "
                  f"(err = {err_pct:+.3f}%, omega = {r['omega_Ha']:.4f} Ha, "
                  f"|z|^2 = {r['dipole_squared_au']:.4f}, t = {elapsed:.1f}s)")
        except Exception as e:
            print(f"  n_max={n_max}: FAILED ({type(e).__name__}: {e})")
            path_a_results.append({"n_max": n_max, "error": str(e)})
    results["convergence_study"]["k_orb_Z"] = path_a_results

    # Path B: Z_eff conventions at fixed n_max=4
    print("\n" + "=" * 70)
    print("PATH B: Z_eff conventions at n_max = 4")
    print("=" * 70)
    z_eff_options = [
        (2.0, "bare Z"),
        (27.0 / 16.0, "variational 1s (27/16)"),
        (1.69, "Slater rules 1s"),
        (1.0, "screened 2p"),
    ]
    z_eff_results = []
    for k_orb, label in z_eff_options:
        t0 = time.time()
        try:
            r = compute_oscillator_strength(
                n_max=4, Z=2, k_orb=k_orb, verbose=False,
            )
            elapsed = time.time() - t0
            r["label"] = label
            r["wall_seconds"] = elapsed
            err_pct = (r["f_length"] - F_DRAKE_2P_1S) / F_DRAKE_2P_1S * 100.0
            r["error_pct"] = err_pct
            z_eff_results.append(r)
            print(f"  k_orb = {k_orb:.4f} ({label}):")
            print(f"    f = {r['f_length']:.6f}  (err = {err_pct:+.3f}%)")
            print(f"    omega = {r['omega_Ha']:.4f} Ha,  |z|^2 = {r['dipole_squared_au']:.4f}")
            print(f"    E_1S = {r['E_1S_Ha']:.6f}  E_2P = {r['E_2P_Ha']:.6f}")
        except Exception as e:
            print(f"  k_orb = {k_orb}: FAILED ({type(e).__name__}: {e})")
            z_eff_results.append({"k_orb": k_orb, "label": label, "error": str(e)})
    results["z_eff_study"] = z_eff_results

    # Diagnostic summary
    print("\n" + "=" * 70)
    print("STURMIAN CLOSURE DIAGNOSTIC")
    print("=" * 70)
    f_values = [r["f_length"] for r in path_a_results if "f_length" in r]
    if len(f_values) >= 2:
        print(f"\nPath A convergence (k_orb = Z = 2):")
        print(f"  n_max=2: f = {f_values[0]:.6f}  (err = {path_a_results[0].get('error_pct',0):+.3f}%)")
        for i in range(1, len(f_values)):
            df = f_values[i] - f_values[i-1]
            print(f"  n_max={i+2}: f = {f_values[i]:.6f}  delta = {df:+.6f}  "
                  f"(err = {path_a_results[i].get('error_pct',0):+.3f}%)")
        if len(f_values) >= 4:
            r1 = abs((f_values[-1] - f_values[-2]) / max(abs(f_values[-1]), 1e-12))
            r2 = abs((f_values[-2] - f_values[-3]) / max(abs(f_values[-2]), 1e-12))
            convergence = "rapid (sub-percent)" if r1 < 0.01 else (
                "moderate" if r1 < 0.05 else "slow / drift")
            print(f"\n  Successive ratios: {r2:.4f}, {r1:.4f}")
            print(f"  Convergence verdict: {convergence}")

    print("\n" + "=" * 70)
    print(f"Reference:  f(2^1P -> 1^1S) = {F_DRAKE_2P_1S:.5f} (Drake & Yan 1992)")
    print("=" * 70)

    # Save
    out_file = PROJECT_ROOT / "debug" / "data" / "he_oscillator_v1.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as fp:
        json.dump(results, fp, indent=2, default=float)
    print(f"\nSaved: {out_file}")

    return results


if __name__ == "__main__":
    main()
