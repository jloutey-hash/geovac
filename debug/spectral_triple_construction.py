"""
Driver script: build the Fock-projected S^3 spectral triple at n_max=2,3.

Prints full diagnostics including:
  - Hilbert space dimensions and sector structure
  - D spectrum (exact algebraic + numerical)
  - [D, pi(f)] non-triviality check
  - Connes distance lower bounds
  - Full axiom verification
  - Pi-free certificate
  - Supertrace values
  - Heat-kernel spectral action

Usage:
    python debug/spectral_triple_construction.py
"""

import json
import sys
import time

import numpy as np
from sympy import Integer, N, Rational

from geovac.spectral_triple import FockSpectralTriple


def run_diagnostics(n_max: int, skip_slow: bool = False) -> dict:
    """Run full spectral triple diagnostics at given n_max."""
    print(f"\n{'='*70}")
    print(f"  Fock Spectral Triple at n_max = {n_max}")
    print(f"{'='*70}")

    t0 = time.time()
    st = FockSpectralTriple(n_max=n_max, j_type="permutation")
    t_build = time.time() - t0

    print(f"\nConstruction time: {t_build:.2f}s")
    print(f"repr: {st}")
    print(f"dim_H = {st.dim_H}")
    print(f"n_sectors = {st.n_sectors}")
    print(f"sectors = {st.sectors}")

    hs = st.hilbert_space_summary()
    print(f"chi+ states = {hs['chi_plus']}")
    print(f"chi- states = {hs['chi_minus']}")

    # Spectrum
    print(f"\n--- D eigenvalues (numerical) ---")
    D = st.dirac_operator
    arr = np.array(D.tolist(), dtype=np.float64)
    evals = sorted(np.linalg.eigvalsh(arr))
    for ev in evals:
        print(f"  {ev:+.10f}")

    # Spectrum pairing check
    paired = True
    for a, b in zip(evals, sorted(-np.array(evals))):
        if abs(a - b) > 1e-8:
            paired = False
            break
    print(f"Spectrum is paired (lambda, -lambda): {paired}")

    # Commutator non-triviality
    print(f"\n--- Commutator [D, pi(f)] ---")
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    f_generic = [Integer(primes[k]) for k in range(st.n_sectors)]
    t0 = time.time()
    comm = st.commutator(f_generic)
    t_comm = time.time() - t0
    n_nz = sum(1 for i in range(st.dim_H) for j in range(st.dim_H) if comm[i, j] != 0)
    print(f"f = {f_generic[:st.n_sectors]}")
    print(f"[D, pi(f)] nonzero entries: {n_nz} / {st.dim_H**2}")
    print(f"Commutator time: {t_comm:.2f}s")

    # Connes distances
    print(f"\n--- Connes distance lower bounds ---")
    distances = {}
    n_pairs = st.n_sectors * (st.n_sectors - 1) // 2
    if n_pairs <= 36:  # only compute for small systems
        for i in range(st.n_sectors):
            for j in range(i + 1, st.n_sectors):
                t0 = time.time()
                d = st.connes_distance(i, j)
                t_d = time.time() - t0
                key = f"{st.sectors[i]}-{st.sectors[j]}"
                val = float(d)
                distances[key] = val
                print(f"  d({st.sectors[i]}, {st.sectors[j]}) >= {val:.6f}  ({t_d:.2f}s)")
    else:
        print(f"  Skipping ({n_pairs} pairs, too many)")

    # Axiom checks
    print(f"\n--- Axiom Checks ---")
    results = {}

    t0 = time.time()
    results["D_selfadjoint"] = st.check_selfadjoint()
    print(f"D = D^T:                  {results['D_selfadjoint']}")

    results["gamma_squared_I"] = st.check_grading_square()
    print(f"gamma^2 = I:              {results['gamma_squared_I']}")

    results["gamma_selfadjoint"] = st.check_grading_selfadjoint()
    print(f"gamma = gamma^T:          {results['gamma_selfadjoint']}")

    dg_result, dg_msg = st.check_D_gamma_anticommute()
    results["D_gamma_anticommute"] = dg_result
    print(f"{{D, gamma}} = 0:           {dg_result}")
    print(f"  -> {dg_msg}")

    results["Lambda_gamma_commute"] = st.check_D_gamma_commute_diagonal()
    print(f"[Lambda, gamma] = 0:      {results['Lambda_gamma_commute']}")

    j2_ok, j2_eps = st.check_J_squared()
    results["J_squared"] = j2_ok
    results["J_squared_epsilon"] = j2_eps
    print(f"J^2 = epsilon*I:          {j2_ok} (epsilon = {j2_eps})")

    jd_ok, jd_msg = st.check_J_D_relation()
    results["J_D_relation"] = jd_ok
    print(f"JD = DJ:                  {jd_ok}")
    print(f"  -> {jd_msg}")

    results["J_gamma_commute"] = st.check_J_gamma_commute()
    print(f"J*gamma = gamma*J:        {results['J_gamma_commute']}")

    t0 = time.time()
    results["order_zero"] = st.check_order_zero()
    t_o0 = time.time() - t0
    print(f"Order-zero:               {results['order_zero']}  ({t_o0:.2f}s)")

    if not skip_slow:
        t0 = time.time()
        o1_ok, o1_msg = st.check_order_one()
        t_o1 = time.time() - t0
        results["order_one"] = o1_ok
        print(f"Order-one:                {o1_ok}  ({t_o1:.2f}s)")
        print(f"  -> {o1_msg}")
    else:
        print("Order-one:                SKIPPED (slow)")
        results["order_one"] = "SKIPPED"

    results["pi_free"] = st.verify_pi_free()
    print(f"Pi-free certificate:      {results['pi_free']}")

    # Supertraces
    print(f"\n--- Supertraces ---")
    str_D = st.supertrace(st.dirac_operator)
    str_I = st.supertrace(st.grading * st.grading)
    str_D2 = st.supertrace(st.dirac_operator * st.dirac_operator)
    print(f"Str(D)   = {str_D}")
    print(f"Str(I)   = {str_I}")
    print(f"Str(D^2) = {str_D2}")
    print(f"Tr(gamma)= {st.grading.trace()}")

    # Heat-kernel spectral action
    print(f"\n--- Heat-kernel spectral action ---")
    for t_val in [Rational(1, 100), Rational(1, 10), Rational(1, 2), Integer(1)]:
        sa = st.spectral_action_heat(t_val)
        print(f"Tr(exp(-{t_val}*D^2)) = {float(N(sa, 10)):.8f}")

    # Kramers J comparison
    print(f"\n--- Kramers J comparison ---")
    st_k = FockSpectralTriple(n_max=n_max, j_type="kramers")
    j2_ok_k, j2_eps_k = st_k.check_J_squared()
    print(f"Kramers J^2 = epsilon*I:  {j2_ok_k} (epsilon = {j2_eps_k})")
    jd_ok_k, jd_msg_k = st_k.check_J_D_relation()
    print(f"Kramers JD = -DJ:         {jd_ok_k}")
    print(f"  -> {jd_msg_k}")
    print(f"Kramers J*gamma=gamma*J:  {st_k.check_J_gamma_commute()}")

    # Collect data for JSON
    data = {
        "n_max": n_max,
        "dim_H": st.dim_H,
        "n_sectors": st.n_sectors,
        "sectors": [list(s) for s in st.sectors],
        "chi_plus": hs["chi_plus"],
        "chi_minus": hs["chi_minus"],
        "eigenvalues_numerical": [float(e) for e in evals],
        "spectrum_paired": paired,
        "commutator_nonzero_entries": n_nz,
        "connes_distances": distances,
        "axiom_checks": {k: (v if isinstance(v, (bool, int, str)) else str(v))
                         for k, v in results.items()},
        "supertraces": {
            "Str_D": str(str_D),
            "Str_I": str(str_I),
            "Str_D2": str(str_D2),
            "Tr_gamma": str(st.grading.trace()),
        },
        "kramers_J_squared_epsilon": j2_eps_k,
    }
    return data


if __name__ == "__main__":
    all_data = {}

    # n_max = 2: full diagnostics
    all_data["n_max_2"] = run_diagnostics(2)

    # n_max = 3: full diagnostics (order-one is slow but feasible)
    all_data["n_max_3"] = run_diagnostics(3, skip_slow=("--skip-slow" in sys.argv))

    # Save data
    out_path = "debug/data/spectral_triple_construction.json"
    with open(out_path, "w") as f:
        json.dump(all_data, f, indent=2, default=str)
    print(f"\nData saved to {out_path}")
