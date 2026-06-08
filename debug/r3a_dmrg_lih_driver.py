"""
Sprint R3-A: DMRG-on-FCIDUMP for LiH composed at n_max=2.

Round 3 of the hybrid-pipeline arc. Goal: confirm that the P1 FCIDUMP
exporter unblocks chemical-accuracy DMRG-equivalent solving on the
production LiH system, and that the resulting PES reproduces the
CHANGELOG v3.56.0 R_eq ~2.82% match against CCCBDB experimental
~3.014 bohr.

----------------------------------------------------------------------
NAMED BUG (uncovered during this sprint, see memo §3):
  ``ecosystem_export.hamiltonian('LiH', ...)`` calls
  ``build_composed_hamiltonian(spec, pk_in_hamiltonian=False)``
  for ``core_method='pk'`` (the default), which surfaces a ``h1``
  that EXCLUDES the Phillips-Kleinman pseudopotential. The
  ``GeoVacHamiltonian.h1`` property therefore returns the bare
  one-electron matrix only; the ``to_fcidump`` exporter writes
  ``h1`` straight through; and downstream FCI on the FCIDUMP
  recovers ``E = -14.475`` Ha (PK-excluded, the wrong reference)
  instead of the Sprint P2 headline ``E = -14.143`` Ha
  (PK-included, the publication-grade number).

  WORKAROUND: bypass ``ecosystem_export.hamiltonian`` and build
  the FCIDUMP directly from
  ``build_composed_hamiltonian(spec, pk_in_hamiltonian=True)``.
  This driver does exactly that. With the workaround applied,
  the FCIDUMP-FCI matches the native FCI to machine zero on every
  R in the panel (max diff < 1e-13 Ha).

  RECOMMENDED FIX (not applied per "DO NOT MODIFY production code"):
  ``_build_hydride`` should expose the PK-included ``h1`` in the
  ``GeoVacHamiltonian.h1`` attribute, or ``to_fcidump`` should grow
  an ``include_pk: bool = True`` switch that lifts ``h1_pk`` into
  ``h1`` before writing.
----------------------------------------------------------------------

Pipeline at each R (PK-included workaround):
  1. ``build_composed_hamiltonian(spec, pk_in_hamiltonian=True)`` ->
     (h1, eri, ecore, n_electrons).
  2. Wrap by hand in a ``GeoVacHamiltonian`` and call ``to_fcidump``.
  3. Round-trip via ``read_fcidump``; verify bit-exact.
  4. ``coupled_fci_energy`` on the round-tripped integrals at
     n_electrons=4. This IS the DMRG-converged limit: the DMRG
     diagnostic sprint (2026-06-07) established that LiH composed
     is block-decoupled with state-side bond rank <=4 across every
     sub-block, so DMRG at chi_max=4 reaches sector FCI to machine
     zero.
  5. SVD-based Schmidt-rank profile on the FCI ground state at the
     production R, reporting the alpha-vs-beta cut bond dim
     (the coarsest single-cut characterization of DMRG-side bond
     dim that DMRG would have to navigate). The per-sub-block cuts
     of the DMRG diagnostic are the finer-grained version.

A real iterative DMRG sweep is omitted because no DMRG library
(block2 / tenpy / quimb / pyscf) is installed in this environment
AND the DMRG diagnostic sprint established that chi_max=4 reaches
FCI bit-exactly on this Hamiltonian's per-block ground state. A 30-
qubit Hamiltonian whose ground state has bond rank <=4 is decisively
non-blocking for any modern DMRG sweep; the load-bearing finding is
the bond-dim measurement, not the sweep convergence.

Output:
  debug/data/r3a_dmrg_lih.json -- full numerical results
  debug/data/lih_r*.fcidump   -- per-R FCIDUMP exports (round-tripped)

Author: PM agent
Date: 2026-06-07
"""
from __future__ import annotations

import json
import os
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Production headline from Paper 20 / Sprint P2 (Q=30, n_max=2, R=3.015):
#   E_FCI_composed = -14.143000 Ha (bit-stable per debug/data/p2_vqe_benchmark.json)
# CHANGELOG v3.56.0 production R_eq comparison: 2.82% vs CCCBDB ~3.014 bohr.
PAPER_20_E_LIH_FCI = -14.143000  # Ha at R=3.015, n_max=2, composed
CCCBDB_LIH_R_EQ = 3.014           # bohr, experimental

# R-grid: focused on the production minimum + flanks for PES topology
R_GRID = [2.5, 3.015, 3.5, 4.0, 5.0]

# Schmidt-rank cutoffs for the bond-dim sweep (state-side, alpha-vs-beta cut)
CHI_GRID = [1, 2, 4, 8, 16]

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Build + FCIDUMP export (PK-included workaround)
# ---------------------------------------------------------------------------

def build_and_export_fcidump(R: float, max_n: int = 2) -> Dict[str, Any]:
    """Build composed LiH with PK in h1, hand-wrap, export FCIDUMP, read back.

    Bypasses ``ecosystem_export.hamiltonian('LiH', ...)`` because it
    surfaces PK-EXCLUDED integrals (see file docstring's NAMED BUG).
    """
    from geovac.molecular_spec import lih_spec
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.ecosystem_export import GeoVacHamiltonian, read_fcidump

    spec = lih_spec(R=R, max_n=max_n)
    t0 = time.perf_counter()
    res = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
    build_wall = time.perf_counter() - t0

    n_electrons = sum(b.n_electrons for b in spec.blocks)

    H_wrap = GeoVacHamiltonian(
        res["qubit_op"],
        metadata={"system": "LiH", "R_bohr": R, "max_n": max_n,
                  "pk_in_hamiltonian": True,
                  "note": "PK-included via hand wrap (Sprint R3-A workaround)"},
        h1=res["h1"],
        eri=res["eri"],
        ecore=res["nuclear_repulsion"],
        n_electrons=n_electrons,
    )

    fcidump_path = str(DATA_DIR / f"lih_r{R:.3f}_nmax{max_n}.fcidump")
    t0 = time.perf_counter()
    write_meta = H_wrap.to_fcidump(fcidump_path)
    write_wall = time.perf_counter() - t0

    t0 = time.perf_counter()
    parsed = read_fcidump(fcidump_path)
    read_wall = time.perf_counter() - t0

    # Bit-exact round-trip diagnostics
    max_h1_diff = float(np.max(np.abs(res["h1"] - parsed["h1"])))
    max_eri_diff = float(np.max(np.abs(res["eri"] - parsed["eri"])))
    ecore_diff = float(abs(res["nuclear_repulsion"] - parsed["ecore"]))

    return {
        "R": R,
        "fcidump_path": fcidump_path,
        "M": int(res["M"]),
        "Q": int(res["Q"]),
        "n_pauli": int(res["N_pauli"]),
        "n_electrons": int(n_electrons),
        "ecore": float(res["nuclear_repulsion"]),
        "write_meta": {
            "n_one_body_terms": int(write_meta["n_one_body_terms"]),
            "n_two_body_terms": int(write_meta["n_two_body_terms"]),
        },
        "roundtrip_max_h1_diff": max_h1_diff,
        "roundtrip_max_eri_diff": max_eri_diff,
        "roundtrip_ecore_diff": ecore_diff,
        "build_wall_s": round(build_wall, 3),
        "write_wall_s": round(write_wall, 3),
        "read_wall_s": round(read_wall, 3),
        # Pass through the parsed integrals for downstream FCI
        "_h1_parsed": parsed["h1"],
        "_eri_parsed": parsed["eri"],
        "_ecore_parsed": parsed["ecore"],
        # And the native integrals for the consistency check
        "_h1_native": res["h1"],
        "_eri_native": res["eri"],
        "_ecore_native": res["nuclear_repulsion"],
    }


# ---------------------------------------------------------------------------
# Step 2: Sector-restricted FCI on the round-tripped integrals
# ---------------------------------------------------------------------------

def fci_from_roundtripped(rt: Dict[str, Any]) -> Dict[str, Any]:
    """Sector FCI on the FCIDUMP-round-tripped integrals."""
    from geovac.coupled_composition import coupled_fci_energy

    res = {
        "M": rt["M"],
        "h1": rt["_h1_parsed"],
        "eri": rt["_eri_parsed"],
        "nuclear_repulsion": rt["_ecore_parsed"],
    }
    t0 = time.perf_counter()
    out = coupled_fci_energy(res, n_electrons=rt["n_electrons"], verbose=False)
    wall = time.perf_counter() - t0
    return {
        "E_fci_from_fcidump": float(out["E_coupled"]),
        "sector_dim": int(out["n_det"]),
        "fci_wall_s": round(wall, 3),
    }


def fci_native(rt: Dict[str, Any]) -> Dict[str, Any]:
    """Native FCI on the unmodified integrals (consistency check)."""
    from geovac.coupled_composition import coupled_fci_energy

    res = {
        "M": rt["M"],
        "h1": rt["_h1_native"],
        "eri": rt["_eri_native"],
        "nuclear_repulsion": rt["_ecore_native"],
    }
    t0 = time.perf_counter()
    out = coupled_fci_energy(res, n_electrons=rt["n_electrons"], verbose=False)
    wall = time.perf_counter() - t0
    return {
        "E_fci_native": float(out["E_coupled"]),
        "fci_native_wall_s": round(wall, 3),
    }


# ---------------------------------------------------------------------------
# Schmidt-rank profile on the FCI ground state via the alpha-vs-beta cut
# ---------------------------------------------------------------------------

def fci_with_eigvec_via_diag(rt: Dict[str, Any]) -> Tuple[float, np.ndarray, int, int]:
    """Sector FCI returning the ground-state eigenvector.

    Uses the same matrix elements as ``coupled_fci_energy`` but exposes
    the eigenvector. Cross-checks: the energy returned here is bit-exact
    equal to ``coupled_fci_energy``'s, modulo numerical noise.
    """
    from itertools import combinations
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import eigsh
    from geovac.coupled_composition import _excitation_phase  # type: ignore

    M = rt["M"]
    h1 = rt["_h1_parsed"]
    eri = rt["_eri_parsed"]
    nuclear_repulsion = rt["_ecore_parsed"]
    n_electrons = rt["n_electrons"]
    n_up = n_electrons // 2
    n_down = n_electrons // 2

    alpha_strings = list(combinations(range(M), n_up))
    beta_strings = list(combinations(range(M), n_down))
    n_alpha = len(alpha_strings)
    n_beta = len(beta_strings)
    n_det = n_alpha * n_beta

    alpha_idx_map = {s: i for i, s in enumerate(alpha_strings)}
    beta_idx_map = {s: i for i, s in enumerate(beta_strings)}

    H_fci = lil_matrix((n_det, n_det))

    def det_index(a_idx: int, b_idx: int) -> int:
        return a_idx * n_beta + b_idx

    # Diagonal: one-body + Coulomb - Exchange
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            I = det_index(ai, bi)
            E_diag = nuclear_repulsion
            for p in alpha:
                E_diag += h1[p, p]
            for p in beta:
                E_diag += h1[p, p]
            for i_idx in range(n_up):
                for j_idx in range(i_idx + 1, n_up):
                    p, q = alpha[i_idx], alpha[j_idx]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for i_idx in range(n_down):
                for j_idx in range(i_idx + 1, n_down):
                    p, q = beta[i_idx], beta[j_idx]
                    E_diag += eri[p, p, q, q] - eri[p, q, q, p]
            for p in alpha:
                for q in beta:
                    E_diag += eri[p, p, q, q]
            H_fci[I, I] = E_diag

    # Singles within alpha
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for p in alpha:
            for r in range(M):
                if r in alpha_set:
                    continue
                new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                if new_alpha not in alpha_idx_map:
                    continue
                ai_new = alpha_idx_map[new_alpha]
                phase = _excitation_phase(alpha, p, r)
                for bi, beta in enumerate(beta_strings):
                    I = det_index(ai, bi)
                    J = det_index(ai_new, bi)
                    val = phase * h1[r, p]
                    for o in alpha_set - {p}:
                        val += phase * (eri[r, p, o, o] - eri[r, o, o, p])
                    for o in beta:
                        val += phase * eri[r, p, o, o]
                    H_fci[I, J] += val

    # Singles within beta
    for bi, beta in enumerate(beta_strings):
        beta_set = set(beta)
        for p in beta:
            for r in range(M):
                if r in beta_set:
                    continue
                new_beta = tuple(sorted((beta_set - {p}) | {r}))
                if new_beta not in beta_idx_map:
                    continue
                bi_new = beta_idx_map[new_beta]
                phase = _excitation_phase(beta, p, r)
                for ai, alpha in enumerate(alpha_strings):
                    I = det_index(ai, bi)
                    J = det_index(ai, bi_new)
                    val = phase * h1[r, p]
                    for o in beta_set - {p}:
                        val += phase * (eri[r, p, o, o] - eri[r, o, o, p])
                    for o in alpha:
                        val += phase * eri[r, p, o, o]
                    H_fci[I, J] += val

    # Doubles alpha-alpha
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        alpha_list = list(alpha)
        n_a = len(alpha_list)
        for i in range(n_a):
            for j in range(i + 1, n_a):
                p, q = alpha_list[i], alpha_list[j]
                for r in range(M):
                    if r in alpha_set:
                        continue
                    for s in range(r + 1, M):
                        if s in alpha_set:
                            continue
                        new_alpha = tuple(sorted((alpha_set - {p, q}) | {r, s}))
                        if new_alpha not in alpha_idx_map:
                            continue
                        ai_new = alpha_idx_map[new_alpha]
                        mid_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                        ph1 = _excitation_phase(alpha, p, r)
                        ph2 = _excitation_phase(mid_alpha, q, s)
                        phase = ph1 * ph2
                        val = phase * (eri[r, p, s, q] - eri[r, q, s, p])
                        for bi in range(n_beta):
                            I = det_index(ai, bi)
                            J = det_index(ai_new, bi)
                            H_fci[I, J] += val

    # Doubles beta-beta
    for bi, beta in enumerate(beta_strings):
        beta_set = set(beta)
        beta_list = list(beta)
        n_b = len(beta_list)
        for i in range(n_b):
            for j in range(i + 1, n_b):
                p, q = beta_list[i], beta_list[j]
                for r in range(M):
                    if r in beta_set:
                        continue
                    for s in range(r + 1, M):
                        if s in beta_set:
                            continue
                        new_beta = tuple(sorted((beta_set - {p, q}) | {r, s}))
                        if new_beta not in beta_idx_map:
                            continue
                        bi_new = beta_idx_map[new_beta]
                        mid_beta = tuple(sorted((beta_set - {p}) | {r}))
                        ph1 = _excitation_phase(beta, p, r)
                        ph2 = _excitation_phase(mid_beta, q, s)
                        phase = ph1 * ph2
                        val = phase * (eri[r, p, s, q] - eri[r, q, s, p])
                        for ai in range(n_alpha):
                            I = det_index(ai, bi)
                            J = det_index(ai, bi_new)
                            H_fci[I, J] += val

    # Doubles alpha-beta
    for ai, alpha in enumerate(alpha_strings):
        alpha_set = set(alpha)
        for p in alpha:
            for r in range(M):
                if r in alpha_set:
                    continue
                new_alpha = tuple(sorted((alpha_set - {p}) | {r}))
                if new_alpha not in alpha_idx_map:
                    continue
                ai_new = alpha_idx_map[new_alpha]
                ph_a = _excitation_phase(alpha, p, r)
                for bi, beta in enumerate(beta_strings):
                    beta_set = set(beta)
                    for q in beta:
                        for s in range(M):
                            if s in beta_set:
                                continue
                            new_beta = tuple(sorted((beta_set - {q}) | {s}))
                            if new_beta not in beta_idx_map:
                                continue
                            bi_new = beta_idx_map[new_beta]
                            ph_b = _excitation_phase(beta, q, s)
                            I = det_index(ai, bi)
                            J = det_index(ai_new, bi_new)
                            val = ph_a * ph_b * eri[r, p, s, q]
                            H_fci[I, J] += val

    H_csr = H_fci.tocsr()
    H_csr = (H_csr + H_csr.T) / 2.0

    eigvals, eigvecs = eigsh(H_csr, k=1, which="SA")
    E0 = float(eigvals[0])
    psi = np.asarray(eigvecs[:, 0])
    return E0, psi, n_alpha, n_beta


def schmidt_profile(
    psi: np.ndarray,
    n_alpha: int,
    n_beta: int,
    chi_grid: List[int],
) -> Dict[str, Any]:
    """Alpha-vs-beta SVD-based bond-dim characterization of the FCI state.

    This is the coarsest single-cut characterization of MPS bond rank
    that any DMRG sweep on this Hamiltonian would have to support at
    the alpha-vs-beta partition. The DMRG diagnostic sprint did finer
    per-orbital cuts on the per-block decomposition; the alpha-vs-beta
    cut here is a sprint-scale proxy on the full-system FCI vector.
    """
    psi_mat = psi.reshape(n_alpha, n_beta)
    U, S, Vt = np.linalg.svd(psi_mat, full_matrices=False)
    S2 = S ** 2
    total = float(S2.sum())
    truncation = []
    for chi in chi_grid:
        if chi >= len(S):
            retained = total
            truncated = 0.0
        else:
            retained = float(S2[:chi].sum())
            truncated = float(S2[chi:].sum())
        truncation.append({
            "chi": int(chi),
            "retained_weight": retained,
            "truncated_weight": truncated,
            "fidelity_loss_1mF": 1.0 - retained / total,
        })
    return {
        "singular_values": [float(s) for s in S],
        "num_significant_1em6": int(np.sum(S2 > 1e-6)),
        "num_significant_1em10": int(np.sum(S2 > 1e-10)),
        "num_significant_1em14": int(np.sum(S2 > 1e-14)),
        "truncation_panel": truncation,
    }


# ---------------------------------------------------------------------------
# PES analysis
# ---------------------------------------------------------------------------

def find_r_eq(R_arr: List[float], E_arr: List[float]) -> Tuple[float, float]:
    """Quadratic fit around the minimum for R_eq and E_min."""
    R_arr = np.asarray(R_arr)
    E_arr = np.asarray(E_arr)
    i_min = int(np.argmin(E_arr))
    if 0 < i_min < len(R_arr) - 1:
        Rs = R_arr[i_min - 1:i_min + 2]
        Es = E_arr[i_min - 1:i_min + 2]
        coef = np.polyfit(Rs, Es, 2)
        a, b, c = coef
        if abs(a) > 1e-12:
            R_eq = -b / (2.0 * a)
            E_eq = c - b ** 2 / (4.0 * a)
            return float(R_eq), float(E_eq)
    return float(R_arr[i_min]), float(E_arr[i_min])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_and_export_fcidump_balanced(R: float, max_n: int = 2) -> Dict[str, Any]:
    """Same as ``build_and_export_fcidump`` but using the balanced builder.

    The balanced builder includes cross-center V_ne and is the binding-
    accuracy target per CHANGELOG v3.56.0 ("balanced 6.93% R_eq"). The
    composed builder is the qubit-resource benchmark target. Both are
    exported as FCIDUMP files for completeness of the hybrid pipeline.
    """
    from geovac.molecular_spec import lih_spec
    from geovac.balanced_coupled import build_balanced_hamiltonian
    from geovac.ecosystem_export import GeoVacHamiltonian, read_fcidump

    spec = lih_spec(R=R, max_n=max_n)
    t0 = time.perf_counter()
    res = build_balanced_hamiltonian(spec, nuclei=None, verbose=False)
    build_wall = time.perf_counter() - t0

    n_electrons = sum(b.n_electrons for b in spec.blocks)

    H_wrap = GeoVacHamiltonian(
        res["qubit_op"],
        metadata={"system": "LiH", "R_bohr": R, "max_n": max_n,
                  "builder": "balanced_coupled"},
        h1=res["h1"], eri=res["eri"],
        ecore=res["nuclear_repulsion"],
        n_electrons=n_electrons,
    )
    fcidump_path = str(DATA_DIR / f"lih_balanced_r{R:.3f}.fcidump")
    write_meta = H_wrap.to_fcidump(fcidump_path)
    parsed = read_fcidump(fcidump_path)

    return {
        "R": R,
        "fcidump_path": fcidump_path,
        "M": int(res["M"]),
        "n_electrons": int(n_electrons),
        "ecore": float(res["nuclear_repulsion"]),
        "roundtrip_max_h1_diff": float(np.max(np.abs(res["h1"] - parsed["h1"]))),
        "roundtrip_max_eri_diff": float(np.max(np.abs(res["eri"] - parsed["eri"]))),
        "_h1_parsed": parsed["h1"],
        "_eri_parsed": parsed["eri"],
        "_ecore_parsed": parsed["ecore"],
        "build_wall_s": round(build_wall, 3),
    }


def fci_from_roundtripped_at(rt: Dict[str, Any]) -> float:
    from geovac.coupled_composition import coupled_fci_energy
    res = {"M": rt["M"], "h1": rt["_h1_parsed"], "eri": rt["_eri_parsed"],
           "nuclear_repulsion": rt["_ecore_parsed"]}
    out = coupled_fci_energy(res, n_electrons=rt["n_electrons"], verbose=False)
    return float(out["E_coupled"])


def main() -> Dict[str, Any]:
    print("Sprint R3-A: DMRG-on-FCIDUMP for LiH composed at n_max=2")
    print("=" * 72)
    print("Workaround: bypass ecosystem_export.hamiltonian wrapper (PK bug); ")
    print("use build_composed_hamiltonian(spec, pk_in_hamiltonian=True) directly.")
    print("=" * 72)

    per_R_results: List[Dict[str, Any]] = []
    schmidt_R = 3.015

    for R in R_GRID:
        print(f"\n[R = {R:.3f} bohr]")
        rt = build_and_export_fcidump(R)
        print(
            f"  FCIDUMP round-trip: max h1 diff = {rt['roundtrip_max_h1_diff']:.2e}, "
            f"max eri diff = {rt['roundtrip_max_eri_diff']:.2e}, "
            f"ecore diff = {rt['roundtrip_ecore_diff']:.2e}"
        )
        fci_rt = fci_from_roundtripped(rt)
        fci_nat = fci_native(rt)
        E_diff = abs(fci_rt["E_fci_from_fcidump"] - fci_nat["E_fci_native"])
        print(
            f"  E(FCIDUMP) = {fci_rt['E_fci_from_fcidump']:+.6f} Ha   "
            f"E(native)  = {fci_nat['E_fci_native']:+.6f} Ha   "
            f"|diff| = {E_diff:.2e} Ha"
        )

        record = {
            "R": R,
            **{k: v for k, v in rt.items() if not k.startswith("_")},
            **fci_rt,
            **fci_nat,
            "abs_diff_fcidump_vs_native": float(E_diff),
        }

        if abs(R - schmidt_R) < 1e-6:
            print(f"  [Schmidt-rank profile at R = {R:.3f}]")
            t0 = time.perf_counter()
            E0, psi, n_alpha, n_beta = fci_with_eigvec_via_diag(rt)
            sweep_wall = time.perf_counter() - t0
            print(
                f"    manual-FCI E0 = {E0:+.6f} Ha   "
                f"|diff vs coupled_fci| = {abs(E0 - fci_rt['E_fci_from_fcidump']):.2e}"
            )
            sp = schmidt_profile(psi, n_alpha, n_beta, CHI_GRID)
            print(
                f"    Schmidt rank @1e-6: {sp['num_significant_1em6']} / "
                f"@1e-10: {sp['num_significant_1em10']} / "
                f"@1e-14: {sp['num_significant_1em14']}"
            )
            for trow in sp["truncation_panel"]:
                print(
                    f"      chi = {trow['chi']:2d}  retained {trow['retained_weight']:.10f}   "
                    f"truncated {trow['truncated_weight']:.2e}   "
                    f"1-F = {trow['fidelity_loss_1mF']:.2e}"
                )
            sp["manual_fci_E0"] = float(E0)
            sp["manual_fci_wall_s"] = round(sweep_wall, 3)
            sp["n_alpha"] = int(n_alpha)
            sp["n_beta"] = int(n_beta)
            record["schmidt_profile"] = sp

        per_R_results.append(record)

    # ------------------------------------------------------------------
    # Companion: balanced builder over the same R-grid, FCIDUMP path
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("Companion: balanced builder PES (FCIDUMP path)")
    print("=" * 72)
    balanced_per_R = []
    for R in R_GRID + [3.224]:  # add the CHANGELOG balanced R_eq for reference
        rt_b = build_and_export_fcidump_balanced(R)
        E_b = fci_from_roundtripped_at(rt_b)
        balanced_per_R.append({
            "R": R,
            "M": rt_b["M"],
            "n_electrons": rt_b["n_electrons"],
            "ecore": rt_b["ecore"],
            "E_fci_from_fcidump": float(E_b),
            "roundtrip_max_h1_diff": rt_b["roundtrip_max_h1_diff"],
            "roundtrip_max_eri_diff": rt_b["roundtrip_max_eri_diff"],
            "fcidump_path": rt_b["fcidump_path"],
        })
        print(f"  R={R:.3f}  E_balanced_FCIDUMP = {E_b:+.6f} Ha")
    balanced_per_R.sort(key=lambda r: r["R"])

    # PES topology (composed)
    R_arr = [r["R"] for r in per_R_results]
    E_arr = [r["E_fci_from_fcidump"] for r in per_R_results]
    R_eq, E_eq = find_r_eq(R_arr, E_arr)
    err_pct = 100.0 * (R_eq - CCCBDB_LIH_R_EQ) / CCCBDB_LIH_R_EQ

    print("\n" + "=" * 72)
    print("PES summary (composed FCIDUMP -> FCI):")
    for r, e in zip(R_arr, E_arr):
        print(f"  R = {r:.3f} bohr   E = {e:+.6f} Ha")
    print(f"  R_eq (parabola fit) = {R_eq:.4f} bohr   E_eq = {E_eq:+.6f} Ha")
    print(f"  CCCBDB R_eq         = {CCCBDB_LIH_R_EQ:.4f} bohr")
    print(f"  error                = {err_pct:+.2f} %")
    print(f"  Paper 20 / Sprint P2 ref E@R=3.015 = {PAPER_20_E_LIH_FCI:+.6f} Ha")

    prod = next(r for r in per_R_results if abs(r["R"] - 3.015) < 1e-6)
    paper_20_match_abs = float(abs(prod["E_fci_from_fcidump"] - PAPER_20_E_LIH_FCI))
    print(f"  |E_FCIDUMP(R=3.015) - Paper20_FCI| = {paper_20_match_abs:.2e} Ha")
    print(f"  GO gate: < 1 mHa @ R=3.015 -> {'PASS' if paper_20_match_abs < 1e-3 else 'FAIL'}")
    err_abs = abs(err_pct)
    print(f"  GO gate: R_eq within 5% of CCCBDB -> {'PASS' if err_abs < 5.0 else 'FAIL'} ({err_abs:.2f}%)")

    # Balanced PES topology
    R_arr_b = [r["R"] for r in balanced_per_R]
    E_arr_b = [r["E_fci_from_fcidump"] for r in balanced_per_R]
    R_eq_b, E_eq_b = find_r_eq(R_arr_b, E_arr_b)
    err_pct_b = 100.0 * (R_eq_b - CCCBDB_LIH_R_EQ) / CCCBDB_LIH_R_EQ
    print("\nPES summary (balanced FCIDUMP -> FCI):")
    for r, e in zip(R_arr_b, E_arr_b):
        print(f"  R = {r:.3f} bohr   E = {e:+.6f} Ha")
    print(f"  R_eq (parabola fit) = {R_eq_b:.4f} bohr   error = {err_pct_b:+.2f}%")
    print(f"  CHANGELOG v3.56.0 balanced R_eq = 3.224 bohr (6.93% err)")

    payload = {
        "system": "LiH",
        "max_n": 2,
        "Q": 30,
        "n_pauli": 333,
        "R_grid": R_arr,
        "per_R": per_R_results,
        "balanced_per_R": balanced_per_R,
        "balanced_R_eq_fit": R_eq_b,
        "balanced_E_eq_fit": E_eq_b,
        "balanced_R_eq_err_pct": err_pct_b,
        "R_eq_fit": R_eq,
        "E_eq_fit": E_eq,
        "cccbdb_R_eq": CCCBDB_LIH_R_EQ,
        "R_eq_err_pct": err_pct,
        "paper_20_E_ref": PAPER_20_E_LIH_FCI,
        "paper_20_match_abs_Ha": paper_20_match_abs,
        "decision_gate": {
            "chemical_accuracy_at_prod_R_Ha": float(paper_20_match_abs),
            "chemical_accuracy_pass": bool(paper_20_match_abs < 1e-3),
            "R_eq_err_pct_abs": float(err_abs),
            "R_eq_pass": bool(err_abs < 5.0),
        },
        "named_bug": {
            "where": "ecosystem_export.hamiltonian('LiH', ...) for core_method='pk'",
            "what": "exposes h1 EXCLUDING the Phillips-Kleinman pseudopotential",
            "impact": "FCIDUMP file is PK-excluded -> downstream FCI gets wrong E by ~+0.33 Ha at R=3.015",
            "workaround": "bypass ecosystem wrapper; call build_composed_hamiltonian(spec, pk_in_hamiltonian=True) directly and hand-wrap into GeoVacHamiltonian for to_fcidump",
            "recommended_fix": "either expose pk-included h1 in _build_hydride, or grow to_fcidump(include_pk=True) switch",
        },
    }
    for r in payload["per_R"]:
        r.pop("_h1_parsed", None)
        r.pop("_eri_parsed", None)
        r.pop("_ecore_parsed", None)
        r.pop("_h1_native", None)
        r.pop("_eri_native", None)
        r.pop("_ecore_native", None)
    for r in payload["balanced_per_R"]:
        r.pop("_h1_parsed", None)
        r.pop("_eri_parsed", None)
        r.pop("_ecore_parsed", None)

    out_path = DATA_DIR / "r3a_dmrg_lih.json"
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"\nResults -> {out_path}")
    return payload


if __name__ == "__main__":
    main()
