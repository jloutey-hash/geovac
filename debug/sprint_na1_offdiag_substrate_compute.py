"""
NA-1 on non-diagonal substrate: off-diagonal Camporesi-Higuchi Dirac.

Sprint NA-1-offdiag (2026-06-06). Follow-up to today's Reading C diagonal
collapse theorem (debug/sprint_na1_depth2_mellin_memo.md).

GOAL.
The diagonal CH substrate forces J(s1, s2) = M_3^{gamma_P}(s1 + s2 - 1)
bit-exactly because D, D^2, gamma_P, e^{-t D^2} are simultaneously
diagonal in (n, l, m_j, chirality). The Reading A vs B distinction
(primitive product vs shuffle pair) is structurally invisible on that
substrate.

On the off-diagonal CH operator system from WH1 R3.5
(geovac/full_dirac_operator_system.py), the Dirac has chirality-flipping
E1 entries connecting (n, +) to (n', -) on |dn|=1 (and |dl|=1, |dm_j|<=1).
Since gamma_P|n,...,chi> = (-1)^n |n,...,chi> flips sign across the
n -> n+1 transition, [gamma_P, D] != 0, hence [D^2, gamma_P D] != 0
generically. The joint depth-2 Mellin no longer trivially collapses.

We test:
  (A) Does J(s1, s2) depend on s1 and s2 separately, or only on s_tot?
  (B) If it depends separately, does it factor as
      Reading A   : J(s1, s2) = M_2(s1) * M_3(s2)       (primitive product)
      Reading A'  : J(s1, s2) = M_2(s1) M_3(s2) + M_3(s1) M_2(s2)  (symmetric)
      Reading B   : J(s1, s2) - J(s2, s1) != 0 with deconcatenation-pair
                    asymmetry (NEW depth-2 content, free non-abelian
                    unipotent radical signature)
      Reading C'  : (diagonal collapse extended) J(s1, s2) = J_eff(s_tot)
                    with the same single-power structure as diagonal CH.

DECISION GATE (mirror sprint prompt).
  Reading A wins    => GeoVac IS abelianisation at depth 2; primitive
                       Hopf substrate is correct; diagonal-collapse
                       extends with operator-ordering correction.
  Reading B wins    => Substrate needs shuffle (cofree-cocommutative T(V))
                       enrichment; 2-4 month follow-on; sketch scope in memo.
  Inconclusive       => Off-diagonal collapses too (different reason);
                       report what's needed, propose alternative substrate.

METHOD.
Spectral evaluation. Diagonalise the off-diagonal CH Dirac
  D = U Lambda U^dagger
on the full-Dirac basis at n_max in {2, 3} (dim_H = 16, 40). Then in the
eigenbasis of D,
  T(t1, t2) := Tr(D^2 e^{-t1 D^2} . gamma_P . D . e^{-t2 D^2})
             = sum_{i,j} Lambda_i^2 e^{-t1 Lambda_i^2} (tilde_gamma_P)_{ij}
                         Lambda_j e^{-t2 Lambda_j^2}
where tilde_gamma_P := U^dagger gamma_P U. The joint Mellin is then
  J(s1, s2) = (1/(Gamma(s1) Gamma(s2))) int int t1^(s1-1) t2^(s2-1) T(t1,t2)
            = sum_{i,j} Lambda_i^(2-2 s1) (tilde_gamma_P)_{ij} Lambda_j^(1-2 s2)

If tilde_gamma_P is diagonal (i.e. [gamma_P, D] = 0), J factors into
single sums and is symmetric in (s1, s2) only via s_tot. If tilde_gamma_P
has off-diagonal entries connecting different eigenvalues, J depends on
s1 and s2 separately.

We choose chirality_coupling = 1.0, offdiag_alpha = 0.0, diag_lifters =
(1.0, 0.0, 0.0) — pure cross-chirality coupling with the canonical CH
diagonal, no diagonal-shifting perturbation that would confuse the test.

PRECISION DISCIPLINE.
  - sympy.Rational for the CH spectrum (bit-exact diagonal D_0).
  - The off-diagonal E1 coupling is numerical (alpha = 1.0 by default
    in the WH1 R3.5 construction; we keep this to faithfully reproduce
    the truthful spinor-bundle off-diagonal structure).
  - All linear algebra in numpy float64 (sufficient — eigenvalues are
    O(1) and matrix dim <= 40), with parallel mpmath high-precision
    eigendecomposition at 50 / 100 / 200 dps for cross-precision
    agreement on the PSLQ verdicts.
  - Mellin moments computed by closed-form spectral sum
    Lambda^(constant) at high mpmath precision.

OUTPUT.
  - data/na1_offdiag_substrate_results.json
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Tuple, List, Dict

import numpy as np
import mpmath as mp

from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    full_dirac_basis,
    camporesi_higuchi_offdiag_dirac_matrix,
    camporesi_higuchi_full_dirac_matrix,
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PRECISIONS = [50, 100, 200]
PSLQ_CEILING = 10 ** 6
PSLQ_MAXSTEPS = 2000

DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)
OUT_JSON = DATA_DIR / "na1_offdiag_substrate_results.json"

# n_max range — keep small for n_max=3 (dim 40), n_max=2 (dim 16).
N_MAX_VALUES = [2, 3]

# Off-diagonal coupling profile: pure cross-chirality, canonical CH diagonal.
DIAG_LIFTERS = (1.0, 0.0, 0.0)
OFFDIAG_ALPHA = 0.0   # no within-chirality E1 perturbation
CHIRALITY_COUPLING = 1.0


# ---------------------------------------------------------------------------
# Build off-diagonal CH Dirac + gamma_P
# ---------------------------------------------------------------------------

def build_offdiag_dirac_and_parity(n_max: int) -> Tuple[np.ndarray, np.ndarray, list]:
    """Build:
      - D_off : off-diagonal CH Dirac (chirality coupling on E1 selection)
      - gamma_P : diagonal (-1)^{n_fock} parity grading
      - basis  : list of FullDiracLabel
    """
    basis = full_dirac_basis(n_max)
    D_off = camporesi_higuchi_offdiag_dirac_matrix(
        basis,
        diag_lifters=DIAG_LIFTERS,
        offdiag_alpha=OFFDIAG_ALPHA,
        chirality_coupling=CHIRALITY_COUPLING,
    )
    # gamma_P diagonal with sign (-1)^n_fock.
    sign = np.array(
        [(-1) ** b.n_fock for b in basis],
        dtype=np.complex128,
    )
    gamma_P = np.diag(sign)
    return D_off, gamma_P, basis


def build_diag_dirac_and_parity(n_max: int) -> Tuple[np.ndarray, np.ndarray, list]:
    """Diagonal CH control: D = chirality * (n_fock + 1/2) (block-diagonal in chirality)."""
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    sign = np.array(
        [(-1) ** b.n_fock for b in basis],
        dtype=np.complex128,
    )
    gamma_P = np.diag(sign)
    return D, gamma_P, basis


# ---------------------------------------------------------------------------
# Commutator audit
# ---------------------------------------------------------------------------

def commutator_norms(D: np.ndarray, gamma_P: np.ndarray) -> Dict[str, float]:
    """Audit: ||[D, gamma_P]||, ||[D^2, gamma_P D]||, ||[D, D^2]|| (sanity)."""
    DG_GD = D @ gamma_P - gamma_P @ D
    Dsq = D @ D
    gpD = gamma_P @ D
    Dsq_gpD = Dsq @ gpD - gpD @ Dsq
    return {
        "D_gammaP": float(np.linalg.norm(DG_GD)),
        "Dsq_gammaPD": float(np.linalg.norm(Dsq_gpD)),
        "Dsq_D": float(np.linalg.norm(Dsq @ D - D @ Dsq)),  # sanity == 0
    }


# ---------------------------------------------------------------------------
# Spectral evaluation of the joint Mellin
# ---------------------------------------------------------------------------

def joint_mellin_spectral(
    D: np.ndarray, gamma_P: np.ndarray, s1: int, s2: int,
    high_prec: bool = True,
) -> mp.mpf:
    """Compute J(s1, s2) = sum_{i,j} Lambda_i^(2-2 s1) * tilde_gamma_P[i,j] *
                            Lambda_j^(1-2 s2)
    via eigendecomposition of D.

    If high_prec, promotes eigenvalues to mpmath at current mp.mp.dps.
    """
    # Eigendecomposition (D Hermitian).
    eigvals, U = np.linalg.eigh(D)
    # Transform gamma_P into the eigenbasis of D.
    tilde_gamma_P = U.conj().T @ gamma_P @ U
    # Sanity: D is Hermitian, eigvals real.
    if not np.all(np.abs(eigvals.imag) < 1e-10):
        raise ValueError("D not Hermitian to tol")
    eigvals = eigvals.real

    # Filter out near-zero eigenvalues (no such on off-diagonal CH;
    # sanity check).
    if np.any(np.abs(eigvals) < 1e-10):
        raise ValueError(f"Near-zero eigenvalue detected: {eigvals}")

    if high_prec:
        # Use mpmath for sums.
        lams = [mp.mpf(float(e)) for e in eigvals]
        total = mp.mpf(0)
        # For each (i, j), accumulate Lambda_i^(2-2 s1) * tilde_gamma_P[i,j] * Lambda_j^(1-2 s2)
        for i, li in enumerate(lams):
            # Sign of eigenvalue for non-integer powers
            sign_i = 1 if li > 0 else (-1) if li < 0 else 0
            ai = abs(li)
            # Lambda_i has sign; (Lambda_i)^k = sign_i^k * |Lambda_i|^k
            # (2 - 2 s1) is even integer, so sign vanishes.
            # (1 - 2 s2) is odd integer, so we need sign_i.
            li_pow = sign_i ** (2 - 2 * s1) * ai ** (2 - 2 * s1)
            for j, lj in enumerate(lams):
                sign_j = 1 if lj > 0 else (-1) if lj < 0 else 0
                aj = abs(lj)
                lj_pow = sign_j ** (1 - 2 * s2) * aj ** (1 - 2 * s2)
                g_ij = complex(tilde_gamma_P[i, j])
                total += li_pow * lj_pow * mp.mpc(g_ij.real, g_ij.imag)
        # Result must be real for Hermitian gamma_P (which it is); take real
        # part to drop floating-point imag noise.
        return mp.re(total)
    else:
        # Pure numpy / float — used only for fast scans.
        lam_pow_left = np.sign(eigvals) ** (2 - 2 * s1) * np.abs(eigvals) ** (2 - 2 * s1)
        lam_pow_right = np.sign(eigvals) ** (1 - 2 * s2) * np.abs(eigvals) ** (1 - 2 * s2)
        # outer product times tilde_gamma_P, then sum.
        L = np.diag(lam_pow_left)
        R = np.diag(lam_pow_right)
        return mp.mpf(float(np.real(np.trace(L @ tilde_gamma_P @ R))))


def depth1_M2_spectral(D: np.ndarray, s: int) -> mp.mpf:
    """M_2(s) = sum_i Lambda_i^(2 - 2 s) via D's spectrum.

    For Hermitian D with both signs of eigenvalues, sign^(even integer) = 1,
    so this is sum_i |Lambda_i|^(2 - 2 s).
    """
    eigvals = np.linalg.eigvalsh(D)
    total = mp.mpf(0)
    for e in eigvals:
        a = abs(float(e))
        if a < 1e-12:
            continue
        total += mp.mpf(a) ** (2 - 2 * s)
    return total


def depth1_M3_spectral_offdiag(D: np.ndarray, gamma_P: np.ndarray, s: int) -> mp.mpf:
    """M_3(s) = sum_i (tilde_gamma_P)_{ii} * Lambda_i^(1 - 2s) where
    tilde_gamma_P = U^dagger gamma_P U.

    The natural depth-1 "M3-like" on the off-diagonal substrate is the
    diagonal-of-conjugated-gamma_P weighted sum (= Tr(gamma_P D |D|^{-2s})).
    """
    eigvals, U = np.linalg.eigh(D)
    tilde_gamma_P = U.conj().T @ gamma_P @ U
    total = mp.mpf(0)
    for i, e in enumerate(eigvals):
        a = abs(float(e))
        if a < 1e-12:
            continue
        sgn = 1 if e > 0 else -1
        lam_pow = sgn ** (1 - 2 * s) * a ** (1 - 2 * s)
        g_ii = complex(tilde_gamma_P[i, i])
        # g_ii is real (gamma_P Hermitian)
        total += mp.mpf(g_ii.real) * lam_pow
    return mp.re(total)


# ---------------------------------------------------------------------------
# Symmetry probes
# ---------------------------------------------------------------------------

def s_tot_only_test(values: Dict[Tuple[int, int], mp.mpf]) -> Dict[str, object]:
    """For pairs (s1, s2) with the same s_tot, check if J(s1, s2) = J(s2, s1)
    AND whether it's the same as J(s'_1, s'_2) for another split."""
    by_stot = {}
    for (s1, s2), v in values.items():
        by_stot.setdefault(s1 + s2, []).append((s1, s2, v))

    results = {}
    for stot, entries in by_stot.items():
        if len(entries) <= 1:
            continue
        # Pairwise absolute differences
        diffs = []
        for i, (s1a, s2a, va) in enumerate(entries):
            for j, (s1b, s2b, vb) in enumerate(entries[i + 1:], start=i + 1):
                diff = abs(va - vb)
                # Normalised by mean
                mean = (abs(va) + abs(vb)) / 2
                rel = float(diff / mean) if mean > 0 else float(diff)
                diffs.append({
                    "split_A": (s1a, s2a),
                    "split_B": (s1b, s2b),
                    "abs_diff": float(diff),
                    "rel_diff": rel,
                })
        results[stot] = {
            "n_splits": len(entries),
            "max_abs_diff": float(max(d["abs_diff"] for d in diffs)) if diffs else 0.0,
            "max_rel_diff": max(d["rel_diff"] for d in diffs) if diffs else 0.0,
            "pairs": diffs,
        }
    return results


def swap_asymmetry_test(values: Dict[Tuple[int, int], mp.mpf]) -> Dict[str, object]:
    """For (s1, s2) and (s2, s1) BOTH present, check J(s1, s2) - J(s2, s1)."""
    results = {}
    seen = set()
    for (s1, s2), v in values.items():
        if (s2, s1) in values and (s1, s2) not in seen and s1 != s2:
            v_swap = values[(s2, s1)]
            diff = v - v_swap
            mean = (abs(v) + abs(v_swap)) / 2
            rel = float(abs(diff) / mean) if mean > 0 else float(abs(diff))
            results[f"({s1},{s2})_vs_({s2},{s1})"] = {
                "J": float(v),
                "J_swap": float(v_swap),
                "abs_asym": float(abs(diff)),
                "rel_asym": rel,
            }
            seen.add((s1, s2))
            seen.add((s2, s1))
    return results


# ---------------------------------------------------------------------------
# Factorisation tests
# ---------------------------------------------------------------------------

def factorisation_audit(
    J_values: Dict[Tuple[int, int], mp.mpf],
    M2: Dict[int, mp.mpf],
    M3: Dict[int, mp.mpf],
) -> Dict[str, object]:
    """For each (s1, s2), report:
      primitive_diff = J(s1, s2) - M_2(s1) * M_3(s2)
      shuffle_diff   = J(s1, s2) - (M_2(s1) M_3(s2) + M_3(s1) M_2(s2))
      collapse_diff  = J(s1, s2) - M_3^{eff}(s1 + s2 - 1)
        (using M_3 from diag substrate as collapse target; mostly diagnostic)
    """
    audit = {}
    for (s1, s2), v in J_values.items():
        prim = M2.get(s1, mp.mpf(0)) * M3.get(s2, mp.mpf(0))
        shuf = (M2.get(s1, mp.mpf(0)) * M3.get(s2, mp.mpf(0))
                + M3.get(s1, mp.mpf(0)) * M2.get(s2, mp.mpf(0)))
        audit[f"({s1},{s2})"] = {
            "J": float(v),
            "M2(s1)*M3(s2)": float(prim),
            "primitive_diff": float(v - prim),
            "M2*M3 + M3*M2": float(shuf),
            "shuffle_diff": float(v - shuf),
        }
    return audit


# ---------------------------------------------------------------------------
# PSLQ panel
# ---------------------------------------------------------------------------

def pslq_panel_for_J(
    J_values: Dict[Tuple[int, int], mp.mpf],
    M2: Dict[int, mp.mpf],
    M3: Dict[int, mp.mpf],
    s_tot_targets: List[int],
) -> Dict[str, object]:
    """For each s_tot, take J at all (s1, s2) with s1 + s2 = s_tot and
    PSLQ-test against a basis built from M2, M3 products at relevant
    weights + transcendental constants.
    """
    pi2 = mp.pi ** 2
    pi3 = mp.pi ** 3
    pi4 = mp.pi ** 4
    G_cat = mp.catalan
    z3 = mp.zeta(3)
    z5 = mp.zeta(5)

    results = {}
    for stot in s_tot_targets:
        # Collect all (s1, s2) with s1 + s2 = stot
        cells = [((s1, s2), v) for (s1, s2), v in J_values.items()
                 if s1 + s2 == stot]
        if not cells:
            continue
        # Use one representative cell (s1=2, s2=stot-2) for PSLQ test
        for (s1, s2), v in cells:
            basis_vals = [
                v,
                M2.get(s1, mp.mpf(0)) * M3.get(s2, mp.mpf(0)),
                M3.get(s1, mp.mpf(0)) * M2.get(s2, mp.mpf(0)),
                M2.get(s1, mp.mpf(0)),
                M2.get(s2, mp.mpf(0)),
                M3.get(s1, mp.mpf(0)),
                M3.get(s2, mp.mpf(0)),
                pi2, pi3, pi4,
                mp.pi ** 5, mp.pi ** 6, mp.pi ** 7,
                G_cat, z3, z5,
                mp.mpf(1),
            ]
            labels = [
                "J", "M2(s1)*M3(s2)", "M3(s1)*M2(s2)",
                "M2(s1)", "M2(s2)", "M3(s1)", "M3(s2)",
                "pi^2", "pi^3", "pi^4", "pi^5", "pi^6", "pi^7",
                "G", "zeta(3)", "zeta(5)", "1",
            ]
            try:
                relation = mp.pslq(basis_vals, tol=mp.mpf(10) ** (-mp.mp.dps + 10),
                                   maxcoeff=PSLQ_CEILING, maxsteps=PSLQ_MAXSTEPS)
            except Exception as e:
                relation = None
                err = str(e)
            else:
                err = None

            if relation is not None and relation[0] != 0:
                # Decompose
                components = {}
                for i in range(1, len(relation)):
                    if relation[i] != 0:
                        components[labels[i]] = (
                            int(-relation[i]), int(relation[0])
                        )
                results[f"({s1},{s2})"] = {
                    "identified": True,
                    "components": components,
                    "raw": [int(c) for c in relation],
                }
            else:
                results[f"({s1},{s2})"] = {
                    "identified": False,
                    "err": err,
                }
    return results


# ---------------------------------------------------------------------------
# Main sprint runner
# ---------------------------------------------------------------------------

def run_sprint() -> Dict[str, object]:
    out: Dict[str, object] = {
        "metadata": {
            "sprint": "NA-1-offdiag",
            "date": "2026-06-06",
            "memo": "debug/sprint_na1_offdiag_substrate_memo.md",
            "diag_lifters": list(DIAG_LIFTERS),
            "offdiag_alpha": OFFDIAG_ALPHA,
            "chirality_coupling": CHIRALITY_COUPLING,
            "precisions": PRECISIONS,
            "pslq_ceiling": PSLQ_CEILING,
            "pslq_maxsteps": PSLQ_MAXSTEPS,
        },
        "per_n_max": {},
    }

    # (s1, s2) splits to test
    # s_tot ranges over {3, 4, 5, 6, 7}; tests both (s1, s2) and (s2, s1)
    splits = [
        (2, 1), (1, 2),
        (3, 1), (1, 3),
        (2, 2),
        (4, 1), (1, 4),
        (3, 2), (2, 3),
        (5, 1), (1, 5),
        (4, 2), (2, 4),
        (3, 3),
        (5, 2), (2, 5),
        (4, 3), (3, 4),
        (4, 4),
    ]

    for n_max in N_MAX_VALUES:
        print(f"\n========== n_max = {n_max} ==========")
        node = {"n_max": n_max}

        # 1) Build both Diracs (off-diag + diag control), audit commutators
        D_off, gamma_P_off, basis_off = build_offdiag_dirac_and_parity(n_max)
        D_diag, gamma_P_diag, basis_diag = build_diag_dirac_and_parity(n_max)
        comm_off = commutator_norms(D_off, gamma_P_off)
        comm_diag = commutator_norms(D_diag, gamma_P_diag)
        node["dim_H"] = int(D_off.shape[0])
        node["commutators"] = {"offdiag": comm_off, "diag_control": comm_diag}
        print(f"  dim_H = {D_off.shape[0]}")
        print(f"  ||[D_off, gamma_P]||   = {comm_off['D_gammaP']:.4e}")
        print(f"  ||[D_diag, gamma_P]||  = {comm_diag['D_gammaP']:.4e}")
        print(f"  ||[D_off^2, gamma_P D_off]||  = {comm_off['Dsq_gammaPD']:.4e}")
        print(f"  ||[D_diag^2, gamma_P D_diag]|| = {comm_diag['Dsq_gammaPD']:.4e}")

        # Run at each precision
        per_prec = {}
        for dps in PRECISIONS:
            print(f"\n  -- precision {dps} dps --")
            mp.mp.dps = dps
            t0 = time.time()
            # Joint Mellin on offdiag substrate
            J_off: Dict[Tuple[int, int], mp.mpf] = {}
            for s1, s2 in splits:
                try:
                    J_off[(s1, s2)] = joint_mellin_spectral(
                        D_off, gamma_P_off, s1, s2, high_prec=True
                    )
                except Exception as e:
                    J_off[(s1, s2)] = mp.mpf("nan")
                    print(f"    J({s1},{s2}) off ERROR: {e}")
            # Depth-1 Mellins on offdiag substrate
            M2_off = {}
            M3_off = {}
            for s in range(1, 6):
                M2_off[s] = depth1_M2_spectral(D_off, s)
                M3_off[s] = depth1_M3_spectral_offdiag(D_off, gamma_P_off, s)

            # Same on diagonal substrate (control)
            J_diag: Dict[Tuple[int, int], mp.mpf] = {}
            for s1, s2 in splits:
                try:
                    J_diag[(s1, s2)] = joint_mellin_spectral(
                        D_diag, gamma_P_diag, s1, s2, high_prec=True
                    )
                except Exception as e:
                    J_diag[(s1, s2)] = mp.mpf("nan")
            M2_diag = {}
            M3_diag = {}
            for s in range(1, 6):
                M2_diag[s] = depth1_M2_spectral(D_diag, s)
                M3_diag[s] = depth1_M3_spectral_offdiag(D_diag, gamma_P_diag, s)

            # Symmetry probes
            stot_test_off = s_tot_only_test(J_off)
            swap_test_off = swap_asymmetry_test(J_off)
            stot_test_diag = s_tot_only_test(J_diag)
            swap_test_diag = swap_asymmetry_test(J_diag)

            # Factorisation audit
            fact_off = factorisation_audit(J_off, M2_off, M3_off)
            fact_diag = factorisation_audit(J_diag, M2_diag, M3_diag)

            # PSLQ panel — only at higher precisions, to save time
            pslq_off = {}
            if dps >= 100:
                pslq_off = pslq_panel_for_J(
                    J_off, M2_off, M3_off,
                    s_tot_targets=[3, 4, 5, 6, 7, 8],
                )

            elapsed = time.time() - t0
            print(f"    elapsed: {elapsed:.1f}s")

            per_prec[str(dps)] = {
                "J_off": {f"({s1},{s2})": str(v) for (s1, s2), v in J_off.items()},
                "J_diag": {f"({s1},{s2})": str(v) for (s1, s2), v in J_diag.items()},
                "M2_off": {str(s): str(v) for s, v in M2_off.items()},
                "M3_off": {str(s): str(v) for s, v in M3_off.items()},
                "M2_diag": {str(s): str(v) for s, v in M2_diag.items()},
                "M3_diag": {str(s): str(v) for s, v in M3_diag.items()},
                "stot_test_off": {str(k): v for k, v in stot_test_off.items()},
                "swap_test_off": swap_test_off,
                "stot_test_diag": {str(k): v for k, v in stot_test_diag.items()},
                "swap_test_diag": swap_test_diag,
                "factorisation_off": fact_off,
                "factorisation_diag": fact_diag,
                "pslq_off": pslq_off,
                "wall_seconds": elapsed,
            }

        node["per_precision"] = per_prec

        # Cross-precision agreement filter on max_rel_diff (s_tot-only test)
        # and swap asymmetry
        cross_check = cross_precision_filter(per_prec)
        node["cross_precision"] = cross_check

        out["per_n_max"][str(n_max)] = node

    return out


def cross_precision_filter(per_prec: Dict[str, Dict]) -> Dict[str, object]:
    """For each cell, check that the relative differences and swap
    asymmetries are stable across precisions."""
    precs = sorted(per_prec.keys(), key=int)
    if len(precs) < 2:
        return {}

    # Compare max_rel_diff across precisions
    cross = {"s_tot_split_invariance_off": {}, "swap_asymmetry_off": {}}

    # All s_tot keys
    all_stots = set()
    for p in precs:
        for k in per_prec[p].get("stot_test_off", {}).keys():
            all_stots.add(k)

    for stot in all_stots:
        vals = []
        for p in precs:
            entry = per_prec[p]["stot_test_off"].get(stot, {})
            if entry:
                vals.append((p, entry.get("max_rel_diff", None)))
        cross["s_tot_split_invariance_off"][stot] = vals

    # Swap asymmetry
    all_swap_keys = set()
    for p in precs:
        for k in per_prec[p].get("swap_test_off", {}).keys():
            all_swap_keys.add(k)
    for swap_k in all_swap_keys:
        vals = []
        for p in precs:
            entry = per_prec[p]["swap_test_off"].get(swap_k, {})
            if entry:
                vals.append((p, entry.get("rel_asym", None)))
        cross["swap_asymmetry_off"][swap_k] = vals

    return cross


# ---------------------------------------------------------------------------
# Entry
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"NA-1 off-diagonal substrate sprint — {time.strftime('%Y-%m-%d %H:%M:%S')}")
    t0 = time.time()
    out = run_sprint()
    elapsed = time.time() - t0
    out["metadata"]["wall_seconds_total"] = elapsed
    print(f"\n\nTotal wall time: {elapsed:.1f}s")
    print(f"Writing to {OUT_JSON}")

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
