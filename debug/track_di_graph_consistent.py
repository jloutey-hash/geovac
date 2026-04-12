"""Sprint 3D: Graph-eigenbasis consistent CI — close the 0.075% gap.

Verifies FCI basis invariance: transforming V_ee to the graph eigenbasis
should give IDENTICAL energies to the graph-native approach (Sprint 3C).

If confirmed, the ~0.075% floor is NOT from basis mismatch — it's from
cusp/truncation (embedding exchange constant content).

Also characterizes the convergence floor via power-law fitting.
"""

import numpy as np
import json
import time

from geovac.casimir_ci import (
    build_graph_native_fci,
    build_graph_consistent_fci,
)

E_EXACT = -2.903724377


def run_comparison(n_max_range: range) -> dict:
    """Compare graph-native and graph-consistent FCI at each n_max."""
    results = {
        "metadata": {
            "track": "DI Sprint 3D",
            "description": "FCI basis invariance test: graph-native vs graph-consistent",
            "E_exact_He": E_EXACT,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        },
        "comparison": [],
    }

    for n_max in n_max_range:
        print(f"\n{'='*60}")
        print(f"n_max = {n_max}")
        print(f"{'='*60}")

        # Graph-native (Sprint 3C)
        t0 = time.time()
        H_native = build_graph_native_fci(Z=2, n_max=n_max)
        E_native = np.linalg.eigvalsh(H_native)[0]
        t_native = time.time() - t0

        # Graph-consistent (Sprint 3D)
        t0 = time.time()
        H_consistent = build_graph_consistent_fci(Z=2, n_max=n_max)
        E_consistent = np.linalg.eigvalsh(H_consistent)[0]
        t_consistent = time.time() - t0

        diff = abs(E_native - E_consistent)
        err_native = abs((E_native - E_EXACT) / E_EXACT) * 100
        err_consistent = abs((E_consistent - E_EXACT) / E_EXACT) * 100

        print(f"  n_configs     = {H_native.shape[0]}")
        print(f"  E_native      = {E_native:.12f}  ({err_native:.6f}%)")
        print(f"  E_consistent  = {E_consistent:.12f}  ({err_consistent:.6f}%)")
        print(f"  |difference|  = {diff:.2e} Ha")
        print(f"  Time: native {t_native:.2f}s, consistent {t_consistent:.2f}s")

        entry = {
            "n_max": n_max,
            "n_configs": int(H_native.shape[0]),
            "E_native": float(E_native),
            "E_consistent": float(E_consistent),
            "difference_Ha": float(diff),
            "error_native_pct": float(err_native),
            "error_consistent_pct": float(err_consistent),
            "time_native_s": float(t_native),
            "time_consistent_s": float(t_consistent),
        }
        results["comparison"].append(entry)

    return results


def characterize_floor(results: dict) -> dict:
    """Characterize the convergence floor from the graph-native sequence.

    Uses the Sprint 3C data (n_max=1-7) plus any new data points.
    Fits E(n) = E_inf + C / n^p (3-parameter) to characterize the floor.
    """
    # Load Sprint 3C data for the full n_max=1-7 sequence
    try:
        with open("debug/data/track_di_graph_native.json") as f:
            sprint3c = json.load(f)
        native_seq = sprint3c["graph_native_sequence"]
    except FileNotFoundError:
        print("Sprint 3C data not found; using comparison data only")
        native_seq = [
            {"n_max": e["n_max"], "E0": e["E_native"]}
            for e in results["comparison"]
        ]

    n_vals = np.array([e["n_max"] for e in native_seq])
    E_vals = np.array([e["E0"] for e in native_seq])
    errors = np.abs((E_vals - E_EXACT) / E_EXACT) * 100

    print(f"\n{'='*60}")
    print("Convergence characterization (graph-native sequence)")
    print(f"{'='*60}")
    print(f"{'n_max':>5} {'E0':>14} {'error%':>10} {'delta%':>10}")
    for i, (n, E, err) in enumerate(zip(n_vals, E_vals, errors)):
        delta = errors[i] - errors[i - 1] if i > 0 else float('nan')
        print(f"{int(n):5d} {E:14.8f} {err:10.6f} {delta:10.6f}")

    # 3-parameter fit: E(n) = E_inf + C / n^p
    # Use n_max >= 3 for the fit (avoid n=1,2 transient)
    from scipy.optimize import curve_fit

    def model(n, E_inf, C, p):
        return E_inf + C / n**p

    mask = n_vals >= 3
    try:
        popt, pcov = curve_fit(
            model, n_vals[mask], E_vals[mask],
            p0=[E_EXACT, 0.02, 1.0],
            maxfev=10000,
        )
        E_inf, C_fit, p_fit = popt
        floor_pct = abs((E_inf - E_EXACT) / E_EXACT) * 100

        print(f"\n3-parameter fit (n >= 3): E(n) = E_inf + C/n^p")
        print(f"  E_inf = {E_inf:.8f} Ha")
        print(f"  C     = {C_fit:.6f}")
        print(f"  p     = {p_fit:.4f}")
        print(f"  Floor = {floor_pct:.4f}%")
        print(f"  E_inf error = {(E_inf - E_EXACT)*1000:.4f} mHa")

        floor_data = {
            "E_inf": float(E_inf),
            "C": float(C_fit),
            "p": float(p_fit),
            "floor_pct": float(floor_pct),
            "E_inf_minus_exact_mHa": float((E_inf - E_EXACT) * 1000),
        }
    except Exception as e:
        print(f"\nFit failed: {e}")
        floor_data = {"status": "fit_failed", "error": str(e)}

    results["floor_characterization"] = floor_data
    return results


def main():
    print("Sprint 3D: Graph-Eigenbasis Consistent CI")
    print("=" * 60)
    print("Hypothesis: the 0.075% floor is from basis mismatch between")
    print("graph h1 (off-diagonal) and hydrogenic V_ee (exact radial).")
    print()
    print("Test: transform V_ee to graph eigenbasis and compare FCI.")
    print("If energies are identical -> floor is NOT from mismatch.")
    print()

    # Run comparison at n_max=1-4 (n_max=5 takes ~30s for integral tensor)
    results = run_comparison(range(1, 5))

    # Check basis invariance
    all_identical = all(
        e["difference_Ha"] < 1e-8 for e in results["comparison"]
    )
    results["basis_invariance_verified"] = all_identical

    print(f"\n{'='*60}")
    print("RESULT: FCI BASIS INVARIANCE", "VERIFIED" if all_identical else "FAILED")
    print(f"{'='*60}")

    if all_identical:
        print()
        print("The graph-consistent and graph-native FCI give IDENTICAL")
        print("energies (to machine precision).  The ~0.075% floor is")
        print("NOT from basis mismatch between h1 and V_ee operators.")
        print()
        print("The floor is from cusp + finite-n_max basis truncation")
        print("(embedding exchange constant content in the exchange")
        print("constant taxonomy of Paper 18).")

        results["conclusion"] = (
            "FCI is invariant under orbital rotation. The 0.075% floor "
            "is NOT from basis mismatch. It is from cusp + finite-n_max "
            "truncation (embedding exchange constant content)."
        )
    else:
        results["conclusion"] = "UNEXPECTED: basis invariance not verified"

    # Characterize the floor
    results = characterize_floor(results)

    # Save
    outpath = "debug/data/track_di_graph_consistent.json"
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
