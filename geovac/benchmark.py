"""
Benchmark Suite for 3-Electron Sparse Dynamics Engine
======================================================

Deliverable 2 — Two benchmark tests:

Test A — Unitarity (run_unitarity_test)
  Propagate a delta-kick initial state on the Lithium single-particle lattice
  for 10,000 Crank-Nicolson steps. Log |<psi|psi> - 1| at every step.
  Success criterion: norm deviation < 1e-12 throughout.

Test B — Scaling (run_scaling_benchmark)
  Sweep max_n from 3 to NMAX_MAX for the Lithium system. At each point,
  assemble the full N-electron Hamiltonian (via LatticeIndex), record
  assembly time and number of non-zero entries. Fit a power law
      assembly_time ~ NNZ^alpha
  Success criterion: alpha < 2 (sub-quadratic). Generate publishable
  log-log plot of assembly_time vs NNZ.

Both tests write figures to the output_dir provided (default: ./figures/).

Author: GeoVac Development Team
Date: February 2026
"""

import os
import time
import tracemalloc
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh as _eigsh

try:
    from .lattice import GeometricLattice
    from .dynamics import TimePropagator
    from .lattice_index import LatticeIndex
    from .atomic_solver import AtomicSolver
except ImportError:
    from lattice import GeometricLattice
    from dynamics import TimePropagator
    from lattice_index import LatticeIndex
    from atomic_solver import AtomicSolver


# Maximum nmax to attempt in the scaling benchmark.
# Increase only if memory > 8 GB and time budget > 10 min.
NMAX_MAX: int = 8

# Time budget per nmax in seconds. Assembly attempts beyond this are skipped.
ASSEMBLY_TIME_BUDGET_S: float = 120.0


# -----------------------------------------------------------------------
# Test A — Unitarity
# -----------------------------------------------------------------------

def run_unitarity_test(
    max_n: int = 5,
    n_steps: int = 10_000,
    dt: float = 0.01,
    output_dir: str = "./figures",
) -> Dict:
    """
    Test A: Unitarity of Crank-Nicolson propagation on Lithium lattice.

    Uses the Lithium single-particle Hamiltonian H1 (Z=3, max_n=max_n).
    Initial state: delta function on the 1s state (index 0, lowest energy).
    This represents a "delta-kick" preparation — an instantaneous excitation
    into the bare lattice state before time evolution.

    Propagates for n_steps at timestep dt (default 0.01 a.u. ≈ 0.24 as).

    Success criterion: max(|<psi|psi> - 1|) < 1e-12 over all steps.

    Parameters
    ----------
    max_n : int
        Lattice size for Li (default 5 gives 55 spatial states)
    n_steps : int
        Number of Crank-Nicolson steps (default 10,000)
    dt : float
        Time step in atomic units (default 0.01)
    output_dir : str
        Directory to write norm_vs_step.png

    Returns
    -------
    dict with keys: norm_deviations, max_deviation, passed, elapsed_s
    """
    print("\n" + "=" * 65)
    print("TEST A — Unitarity: Crank-Nicolson on Lithium Lattice")
    print("=" * 65)

    # --- Build Li single-particle Hamiltonian ---
    print(f"\n  Building Li lattice: max_n={max_n}, Z=3")
    solver = AtomicSolver(max_n=max_n, Z=3)
    H = solver.H
    n_states = H.shape[0]
    print(f"  States: {n_states},  H.nnz: {H.nnz}")

    # --- Delta-kick initial state: all weight on index 0 (1s) ---
    psi0 = np.zeros(n_states, dtype=complex)
    psi0[0] = 1.0 + 0j

    # --- Build time propagator ---
    prop = TimePropagator(H, dt=dt)

    # --- Run propagation, log norm at every step ---
    norm_deviations = np.empty(n_steps)
    psi = psi0.copy()

    t0 = time.perf_counter()
    for step in range(n_steps):
        psi = prop.step(psi)
        norm_dev = abs(np.dot(psi.conj(), psi).real - 1.0)
        norm_deviations[step] = norm_dev
    elapsed = time.perf_counter() - t0

    max_dev = float(np.max(norm_deviations))
    passed = max_dev < 1e-12

    print(f"\n  Steps: {n_steps:,}  |  dt = {dt} a.u.  |  elapsed = {elapsed:.2f}s")
    print(f"  Max |<psi|psi> - 1| = {max_dev:.3e}")
    print(f"  {'PASSED' if passed else 'FAILED'}: "
          f"criterion < 1e-12  (achieved {max_dev:.3e})")

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        os.makedirs(output_dir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 4))
        steps = np.arange(1, n_steps + 1)
        ax.semilogy(steps, norm_deviations, linewidth=0.6, color="steelblue",
                    label=r"$|\langle\psi|\psi\rangle - 1|$")
        ax.axhline(1e-12, color="red", linestyle="--", linewidth=1.2,
                   label="Success criterion $10^{-12}$")
        ax.set_xlabel("Step")
        ax.set_ylabel(r"$|\langle\psi|\psi\rangle - 1|$")
        ax.set_title(
            f"Crank-Nicolson Unitarity — Li lattice (max_n={max_n}, "
            f"{n_states} states, {n_steps:,} steps)"
        )
        ax.legend(fontsize=9)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()

        outpath = os.path.join(output_dir, "unitarity_norm_vs_step.png")
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"\n  Plot saved: {outpath}")
    except ImportError:
        print("\n  [WARNING] matplotlib not available — plot skipped.")

    return {
        "norm_deviations": norm_deviations,
        "max_deviation": max_dev,
        "passed": passed,
        "elapsed_s": elapsed,
        "n_states": n_states,
        "n_steps": n_steps,
        "dt": dt,
    }


# -----------------------------------------------------------------------
# Test B — Scaling
# -----------------------------------------------------------------------

def run_scaling_benchmark(
    nmax_range: Optional[List[int]] = None,
    output_dir: str = "./figures",
    time_budget_s: float = ASSEMBLY_TIME_BUDGET_S,
) -> Dict:
    """
    Test B: Assembly-time scaling of the Lithium 3-electron Hamiltonian.

    Sweeps max_n from 3 up to NMAX_MAX, assembles the full N=3 electron
    Hamiltonian at each point via LatticeIndex.assemble_hamiltonian(),
    records:
      - n_sd: number of Slater determinants
      - nnz: non-zero entries in assembled H
      - assembly_time_s: wall-clock seconds for assembly
      - peak_memory_mb: peak RSS change during assembly (tracemalloc)

    Fits a power law  assembly_time ~ NNZ^alpha  via log-log linear regression.
    Success criterion: alpha < 2 (sub-quadratic scaling).

    The log-log plot (assembly_time vs NNZ) is the primary publishable figure.

    Parameters
    ----------
    nmax_range : list of int, optional
        Which max_n values to benchmark. Default: range(3, NMAX_MAX+1).
    output_dir : str
        Directory for the scaling_log_log.png figure.
    time_budget_s : float
        Skip assembly if estimated time > budget (default 120s).

    Returns
    -------
    dict with keys: results (list of per-nmax dicts), alpha, r_squared, passed
    """
    if nmax_range is None:
        nmax_range = list(range(3, NMAX_MAX + 1))

    print("\n" + "=" * 65)
    print("TEST B — Scaling: Li 3-electron Hamiltonian (n=3 FCI)")
    print("=" * 65)
    print(f"\n  Sweeping max_n in {nmax_range}")
    print(f"  Time budget per point: {time_budget_s}s")

    results = []

    for max_n in nmax_range:
        n_sp = 2 * sum(n ** 2 for n in range(1, max_n + 1))
        # Quick combinatorics estimate
        n_sd_est = _comb3(n_sp)
        print(f"\n  max_n={max_n}: ~{n_sp} spin-orbitals, "
              f"~{n_sd_est:,} Slater determinants")

        # Suppress Mulliken warning during benchmark
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            try:
                # Enumerate SD basis (cheap even for large n_sp)
                t_enum0 = time.perf_counter()
                idx = LatticeIndex(
                    n_electrons=3, max_n=max_n, nuclear_charge=3
                )
                t_enum = time.perf_counter() - t_enum0

                # Estimate time using a tiny probe assembly
                # (skip if n_sd is huge and time would overflow budget)
                if n_sd_est > 5_000_000 and time_budget_s < 300:
                    print(f"    Skipping assembly: n_sd={n_sd_est:,} > 5M "
                          f"(exceeds time budget estimate)")
                    results.append({
                        "max_n": max_n,
                        "n_sd": idx.n_sd,
                        "nnz": None,
                        "assembly_time_s": None,
                        "peak_memory_mb": None,
                        "skipped": True,
                        "reason": "n_sd > 5M, time budget",
                    })
                    continue

                # Full assembly with memory tracking
                tracemalloc.start()
                t_asm0 = time.perf_counter()
                H = idx.assemble_hamiltonian()
                t_asm = time.perf_counter() - t_asm0
                _, peak_mem = tracemalloc.get_traced_memory()
                tracemalloc.stop()

                nnz = int(H.nnz)
                peak_mb = peak_mem / 1e6

                print(f"    n_sd={idx.n_sd:,}  nnz={nnz:,}  "
                      f"time={t_asm:.3f}s  peak={peak_mb:.1f} MB")

                results.append({
                    "max_n": max_n,
                    "n_sd": idx.n_sd,
                    "nnz": nnz,
                    "assembly_time_s": t_asm,
                    "peak_memory_mb": peak_mb,
                    "skipped": False,
                })

            except MemoryError:
                print(f"    MemoryError at max_n={max_n} — skipping.")
                results.append({
                    "max_n": max_n,
                    "n_sd": n_sd_est,
                    "nnz": None,
                    "assembly_time_s": None,
                    "peak_memory_mb": None,
                    "skipped": True,
                    "reason": "MemoryError",
                })
            except Exception as exc:
                print(f"    Exception at max_n={max_n}: {exc} — skipping.")
                results.append({
                    "max_n": max_n,
                    "n_sd": n_sd_est,
                    "nnz": None,
                    "assembly_time_s": None,
                    "peak_memory_mb": None,
                    "skipped": True,
                    "reason": str(exc),
                })

    # --- Power-law fit: log(time) = alpha * log(NNZ) + const ---
    valid = [r for r in results if not r["skipped"] and r["nnz"] and r["assembly_time_s"]]
    alpha = None
    r_squared = None
    passed = False

    if len(valid) >= 2:
        log_nnz = np.log10([r["nnz"] for r in valid])
        log_t = np.log10([r["assembly_time_s"] for r in valid])
        coeffs = np.polyfit(log_nnz, log_t, 1)
        alpha = float(coeffs[0])
        # R^2
        predicted = np.polyval(coeffs, log_nnz)
        ss_res = np.sum((log_t - predicted) ** 2)
        ss_tot = np.sum((log_t - log_t.mean()) ** 2)
        r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0
        passed = alpha < 2.0
        print(f"\n  Power-law fit: assembly_time ~ NNZ^{alpha:.3f}  "
              f"(R²={r_squared:.4f})")
        print(f"  {'PASSED' if passed else 'FAILED'}: "
              f"criterion alpha < 2.0  (achieved {alpha:.3f})")
    else:
        print("\n  WARNING: fewer than 2 valid data points — fit skipped.")

    # --- NNZ growth exponent: NNZ ~ max_n^beta ---
    nnz_beta = None
    nnz_beta_r2 = None
    if len(valid) >= 2:
        log_mn = np.log10([r["max_n"] for r in valid])
        log_nnz = np.log10([r["nnz"] for r in valid])
        beta_coeffs = np.polyfit(log_mn, log_nnz, 1)
        nnz_beta = float(beta_coeffs[0])
        pred_b = np.polyval(beta_coeffs, log_mn)
        ss_res_b = np.sum((log_nnz - pred_b) ** 2)
        ss_tot_b = np.sum((log_nnz - log_nnz.mean()) ** 2)
        nnz_beta_r2 = float(1.0 - ss_res_b / ss_tot_b) if ss_tot_b > 0 else 1.0
        print(f"  NNZ growth fit:  NNZ ~ max_n^{nnz_beta:.3f}  "
              f"(R²={nnz_beta_r2:.4f})")
        if nnz_beta > 2.0:
            print(f"  WARNING: super-quadratic NNZ growth (beta={nnz_beta:.2f} > 2). "
                  f"Sparsity mask bounds per-SD connections but cannot contain "
                  f"combinatorial N_SD growth for 3-electron FCI.")

    # --- Plot ---
    _plot_scaling(valid, alpha, r_squared, nnz_beta, nnz_beta_r2, output_dir)

    return {
        "results": results,
        "valid_results": valid,
        "alpha": alpha,
        "r_squared": r_squared,
        "passed": passed,
        "nnz_beta": nnz_beta,
    }


def _plot_scaling(
    valid: List[Dict],
    alpha: Optional[float],
    r_squared: Optional[float],
    nnz_beta: Optional[float],
    nnz_beta_r2: Optional[float],
    output_dir: str,
) -> None:
    """
    Two-panel publishable figure: (left) assembly_time vs NNZ, (right) NNZ vs max_n.

    Left panel shows the O(NNZ) assembly claim (alpha~1 = sub-quadratic).
    Right panel shows total NNZ growth with max_n — exposes whether the
    sparsity mask contains state-space growth or not. If NNZ ~ max_n^beta
    with beta > 2.0, a visible WARNING annotation is added.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if not valid:
            print("  [WARNING] No valid data points — plot skipped.")
            return

        nnz_arr = np.array([r["nnz"] for r in valid])
        time_arr = np.array([r["assembly_time_s"] for r in valid])
        nmax_arr = np.array([r["max_n"] for r in valid], dtype=float)

        os.makedirs(output_dir, exist_ok=True)
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(
            "Li 3-Electron FCI — Sparse Assembly Scaling (Graph-topology JOIN)",
            fontsize=12, fontweight="bold",
        )

        # ---- Left panel: assembly_time vs NNZ ----
        ax = axes[0]
        ax.loglog(nnz_arr, time_arr, "o-", color="navy", markersize=7,
                  linewidth=1.5, label="Measured")
        for nnz, t, mn in zip(nnz_arr, time_arr, nmax_arr):
            ax.annotate(f"$n_{{\\rm max}}={int(mn)}$", xy=(nnz, t),
                        xytext=(4, 4), textcoords="offset points", fontsize=8)

        if alpha is not None and len(nnz_arr) >= 2:
            nnz_fit = np.logspace(np.log10(nnz_arr.min()) - 0.2,
                                  np.log10(nnz_arr.max()) + 0.2, 200)
            log_nnz = np.log10(nnz_arr)
            C = 10 ** (np.mean(np.log10(time_arr) - alpha * log_nnz))
            label = (f"Fit: $t\\propto{{\\rm NNZ}}^{{{alpha:.2f}}}$"
                     + (f"  ($R^2={r_squared:.4f}$)" if r_squared else ""))
            ax.loglog(nnz_fit, C * nnz_fit ** alpha, "--", color="firebrick",
                      linewidth=1.5, label=label)

        # Reference slopes
        nnz_ref = np.logspace(np.log10(nnz_arr.min()), np.log10(nnz_arr.max()), 50)
        t0 = time_arr[0]
        ax.loglog(nnz_ref, t0 * (nnz_ref / nnz_arr[0]), ":", color="green",
                  linewidth=1.0, alpha=0.7, label="$\\alpha=1$ (linear)")
        ax.loglog(nnz_ref, t0 * (nnz_ref / nnz_arr[0]) ** 2, ":", color="orange",
                  linewidth=1.0, alpha=0.7, label="$\\alpha=2$ (quadratic)")

        ax.set_xlabel("NNZ in $H_{FCI}$", fontsize=11)
        ax.set_ylabel("Assembly time (s)", fontsize=11)
        ax.set_title("Assembly Time vs NNZ\n(sub-quadratic claim: $\\alpha<2$)", fontsize=10)
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(True, which="both", alpha=0.3)

        # ---- Right panel: NNZ vs max_n ----
        ax2 = axes[1]
        ax2.loglog(nmax_arr, nnz_arr, "s-", color="darkgreen", markersize=7,
                   linewidth=1.5, label="Measured NNZ")

        if nnz_beta is not None and len(nnz_arr) >= 2:
            mn_fit = np.logspace(np.log10(nmax_arr.min()) - 0.05,
                                 np.log10(nmax_arr.max()) + 0.05, 200)
            log_mn = np.log10(nmax_arr)
            C2 = 10 ** (np.mean(np.log10(nnz_arr) - nnz_beta * log_mn))
            label2 = (f"Fit: ${{\\rm NNZ}}\\propto n_{{\\rm max}}^{{{nnz_beta:.2f}}}$"
                      + (f"  ($R^2={nnz_beta_r2:.4f}$)" if nnz_beta_r2 else ""))
            ax2.loglog(mn_fit, C2 * mn_fit ** nnz_beta, "--", color="purple",
                       linewidth=1.5, label=label2)

            if nnz_beta > 2.0:
                ax2.text(
                    0.05, 0.07,
                    f"WARNING: super-quadratic NNZ growth\n"
                    f"($\\beta={nnz_beta:.2f}>2$) — sparsity mask\n"
                    f"bounds per-SD connections only;\n"
                    f"combinatorial $N_{{SD}}$ growth is unavoidable\n"
                    f"for full 3-electron FCI.",
                    transform=ax2.transAxes, fontsize=8,
                    color="red", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                              edgecolor="red", alpha=0.85),
                )

        # Reference: quadratic in max_n
        mn_ref = np.logspace(np.log10(nmax_arr.min()), np.log10(nmax_arr.max()), 50)
        n0_ref = nnz_arr[0]
        ax2.loglog(mn_ref, n0_ref * (mn_ref / nmax_arr[0]) ** 2, ":", color="orange",
                   linewidth=1.0, alpha=0.7, label="$\\beta=2$ reference")

        ax2.set_xlabel("$n_{\\rm max}$", fontsize=11)
        ax2.set_ylabel("NNZ in $H_{FCI}$", fontsize=11)
        ax2.set_title("NNZ Growth vs Basis Size\n(sparsity mask audit)", fontsize=10)
        ax2.legend(fontsize=8, loc="upper left")
        ax2.grid(True, which="both", alpha=0.3)

        fig.tight_layout()
        outpath = os.path.join(output_dir, "scaling_log_log.png")
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"\n  Primary publishable figure (two-panel) saved: {outpath}")

    except ImportError:
        print("\n  [WARNING] matplotlib not available — scaling plot skipped.")


def _comb3(n: int) -> int:
    """C(n, 3) = n*(n-1)*(n-2)//6."""
    if n < 3:
        return 0
    return n * (n - 1) * (n - 2) // 6


# -----------------------------------------------------------------------
# Gap 2 — Li ground state energy audit
# -----------------------------------------------------------------------

#: Exact non-relativistic Li ground state (Pekeris 1959, NIST-CODATA)
LI_EXACT_ENERGY_HA: float = -7.4781

#: Wall-time budget per max_n before issuing ResourceWarning and skipping (s)
ENERGY_AUDIT_TIME_BUDGET_S: float = 300.0

#: Memory budget before issuing ResourceWarning (bytes, 4 GB)
ENERGY_AUDIT_MEMORY_BUDGET_B: int = 4 * 1024 ** 3


def _run_audit_one_method(
    max_n_list: List[int],
    vee_method: str,
) -> List[Dict]:
    """
    Inner helper: run Li energy audit for a single vee_method.
    Returns list of row dicts (one per max_n).
    """
    rows: List[Dict] = []

    header = (f"{'max_n':>6} {'n_sd':>10} {'E_GeoVac (Ha)':>15} "
              f"{'E_exact (Ha)':>14} {'|dE| (Ha)':>12} {'err %':>8}")
    sep = "-" * len(header)
    print(f"\n  V_ee method: {vee_method.upper()}")
    print("  " + header)
    print("  " + sep)

    for max_n in max_n_list:
        n_sp_est = 2 * sum(n ** 2 for n in range(1, max_n + 1))
        n_sd_est = _comb3(n_sp_est)

        warnings.warn(
            f"[max_n={max_n}, vee={vee_method}] LatticeIndex V_ee approximation active. "
            f"Chordal: V_ee=κ/d²_chord (κ=1 Ha); Mulliken: ~12% overestimation. "
            f"Neither reproduces exact Hartree integrals.",
            UserWarning,
            stacklevel=2,
        )

        try:
            tracemalloc.start()
            t0 = time.perf_counter()

            idx = LatticeIndex(
                n_electrons=3, max_n=max_n, nuclear_charge=3,
                vee_method=vee_method,
            )
            H = idx.assemble_hamiltonian()

            t_asm = time.perf_counter() - t0
            _, peak_mem = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            if t_asm > ENERGY_AUDIT_TIME_BUDGET_S:
                msg = (f"max_n={max_n} vee={vee_method} assembly took {t_asm:.1f}s "
                       f"> budget {ENERGY_AUDIT_TIME_BUDGET_S:.0f}s.")
                warnings.warn(msg, ResourceWarning, stacklevel=1)
                print(f"  [ResourceWarning] {msg}")
                rows.append({
                    "max_n": max_n, "n_sd": idx.n_sd, "vee_method": vee_method,
                    "E_geovac": None, "abs_err": None, "pct_err": None,
                    "peak_mb": peak_mem / 1e6, "skipped": True,
                    "reason": "time budget exceeded",
                })
                continue

            if peak_mem > ENERGY_AUDIT_MEMORY_BUDGET_B:
                msg = (f"max_n={max_n} vee={vee_method} peak memory "
                       f"{peak_mem / 1e9:.2f} GB > budget.")
                warnings.warn(msg, ResourceWarning, stacklevel=1)
                print(f"  [ResourceWarning] {msg}")
                rows.append({
                    "max_n": max_n, "n_sd": idx.n_sd, "vee_method": vee_method,
                    "E_geovac": None, "abs_err": None, "pct_err": None,
                    "peak_mb": peak_mem / 1e6, "skipped": True,
                    "reason": "memory budget exceeded",
                })
                continue

            # Reuse already-assembled H (avoid double-assembly)
            k = min(1, idx.n_sd - 2)
            _eigvals, _ = _eigsh(H, k=k, which="SA")
            eigvals = _eigvals[np.argsort(_eigvals)]
            E_gv = float(eigvals[0])
            abs_err = abs(E_gv - LI_EXACT_ENERGY_HA)
            pct_err = abs_err / abs(LI_EXACT_ENERGY_HA) * 100.0

            print(f"  {max_n:>6} {idx.n_sd:>10,} {E_gv:>15.6f} "
                  f"{LI_EXACT_ENERGY_HA:>14.4f} {abs_err:>12.4f} {pct_err:>7.2f}%")

            rows.append({
                "max_n": max_n, "n_sd": idx.n_sd, "vee_method": vee_method,
                "E_geovac": E_gv, "abs_err": abs_err, "pct_err": pct_err,
                "peak_mb": peak_mem / 1e6, "skipped": False,
            })

        except MemoryError:
            tracemalloc.stop()
            msg = f"max_n={max_n} vee={vee_method} raised MemoryError."
            warnings.warn(msg, ResourceWarning, stacklevel=1)
            print(f"  [ResourceWarning] {msg}")
            rows.append({
                "max_n": max_n, "n_sd": n_sd_est, "vee_method": vee_method,
                "E_geovac": None, "abs_err": None, "pct_err": None,
                "peak_mb": None, "skipped": True, "reason": "MemoryError",
            })
        except Exception as exc:
            tracemalloc.stop()
            print(f"  [ERROR] max_n={max_n} vee={vee_method}: {exc}")
            rows.append({
                "max_n": max_n, "n_sd": n_sd_est, "vee_method": vee_method,
                "E_geovac": None, "abs_err": None, "pct_err": None,
                "peak_mb": None, "skipped": True, "reason": str(exc),
            })

    print("  " + sep)
    return rows


def run_li_energy_audit(
    max_n_list: Optional[List[int]] = None,
    output_dir: str = "./figures",
    vee_methods: Optional[List[str]] = None,
) -> List[Dict]:
    """
    Gap 2: Li ground state energy convergence vs exact non-relativistic value.

    Runs for each vee_method in vee_methods (default: ['chordal', 'mulliken']),
    prints a side-by-side comparison table, and writes the results to
    figures/li_energy_convergence.txt.

    Parameters
    ----------
    max_n_list : list of int, optional
        Basis sizes to audit. Default: [4, 6].
    output_dir : str
        Directory for the txt output.
    vee_methods : list of str, optional
        V_ee methods to compare. Default: ['chordal', 'mulliken'].

    Returns
    -------
    list of dicts for the PRIMARY method (chordal), one per max_n.
    """
    if max_n_list is None:
        max_n_list = [4, 6]
    if vee_methods is None:
        vee_methods = ['chordal', 'mulliken']

    print("\n" + "=" * 65)
    print("GAP 2 — Li Ground State Energy Audit (V_ee Method Comparison)")
    print("=" * 65)
    print(f"\n  Z=3, n_electrons=3, max_n in {max_n_list}")
    print(f"  Reference: E_exact = {LI_EXACT_ENERGY_HA} Ha (Pekeris 1959)")
    print(f"  Methods compared: {vee_methods}")
    print(f"  Memory budget: {ENERGY_AUDIT_MEMORY_BUDGET_B / 1e9:.1f} GB  |  "
          f"Time budget: {ENERGY_AUDIT_TIME_BUDGET_S:.0f}s")

    all_rows: Dict[str, List[Dict]] = {}
    for method in vee_methods:
        all_rows[method] = _run_audit_one_method(max_n_list, method)

    # --- Comparison summary ---
    primary_method = vee_methods[0]
    if len(vee_methods) >= 2:
        print(f"\n  {'max_n':>6}  {'method':>10}  {'E_GeoVac (Ha)':>15}  "
              f"{'err %':>8}  {'vs exact'}")
        print("  " + "-" * 55)
        for max_n in max_n_list:
            for method in vee_methods:
                r = next((x for x in all_rows[method]
                          if x["max_n"] == max_n and not x.get("skipped")), None)
                if r:
                    tag = ""
                    if method == primary_method:
                        other = vee_methods[1]
                        r2 = next((x for x in all_rows[other]
                                   if x["max_n"] == max_n and not x.get("skipped")), None)
                        if r2:
                            delta = r["pct_err"] - r2["pct_err"]
                            tag = (f"  [{delta:+.1f}pp vs {other}]")
                    print(f"  {max_n:>6}  {method:>10}  {r['E_geovac']:>15.6f}  "
                          f"{r['pct_err']:>7.2f}%{tag}")
                else:
                    print(f"  {max_n:>6}  {method:>10}  [SKIPPED]")

    # --- Write txt ---
    os.makedirs(output_dir, exist_ok=True)
    txt_path = os.path.join(output_dir, "li_energy_convergence.txt")
    with open(txt_path, "w") as f:
        f.write("Li 3-Electron Ground State Energy — V_ee Method Comparison\n")
        f.write("=" * 65 + "\n")
        f.write(f"Reference: E_exact = {LI_EXACT_ENERGY_HA} Ha (Pekeris 1959)\n")
        f.write(f"Methods: {vee_methods}\n\n")
        for method in vee_methods:
            f.write(f"\n--- {method.upper()} ---\n")
            header = (f"{'max_n':>6} {'n_sd':>10} {'E_GeoVac':>15} "
                      f"{'E_exact':>12} {'|dE|':>10} {'err%':>7}")
            f.write(header + "\n" + "-" * len(header) + "\n")
            for r in all_rows[method]:
                if r["skipped"]:
                    f.write(f"  {r['max_n']:>4}  [SKIPPED: {r.get('reason', '?')}]\n")
                else:
                    f.write(f"  {r['max_n']:>4} {r['n_sd']:>10,} "
                            f"{r['E_geovac']:>15.6f} "
                            f"{LI_EXACT_ENERGY_HA:>12.4f} "
                            f"{r['abs_err']:>10.4f} "
                            f"{r['pct_err']:>6.2f}%\n")
        f.write("\nNOTE: Errors do NOT converge to zero with max_n for either method.\n")
        f.write("Chordal: graph-native (Paper 7, Eq. 14), kappa_ee=1 Ha, d=2 antipodal.\n")
        f.write("Mulliken: ~12% overestimation; ~10x error for same-orbital pairs.\n")
    print(f"\n  Comparison table written: {txt_path}")

    # Summary line for the primary method
    primary_rows = all_rows[primary_method]
    valid_rows = [r for r in primary_rows if not r.get("skipped", True)]
    if valid_rows:
        best = min(valid_rows, key=lambda r: r["pct_err"])
        print(f"\n[LI ENERGY]  best ({primary_method}): max_n={best['max_n']}, "
              f"E={best['E_geovac']:.4f} Ha, "
              f"error={best['pct_err']:.2f}% vs exact")

    return primary_rows
