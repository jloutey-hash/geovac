"""Sprint Q5'-OffDiag-Closure-nmax4 - L4 multi-year extension to n_max=4
for asymptotic fill-fraction verification of the OffDiag algebra closure.

L4 panel at v3.63.0 (sprint_q5p_offdiag_closure_memo):
  n_max=2: dim A_OD = 84, fill = 84/256 = 21/64 = 0.328
  n_max=3: dim A_OD = 592, fill = 592/1600 = 37/100 = 0.370

Predicted growth law (tentative, 2 data points): ~ dim(H)^2.13.
  At n_max=2, dim H = 16 -> H^2.13 = 547 (off by ~15% vs measured 84;
  ratio dim H^2.13 / measured ~6.5 — the 2.13 fit is to the LOG of the
  ratio over a small range, not to absolute scale).
  At n_max=3, dim H = 40 -> H^2.13 = 4136 (off by ~7x vs measured 592).

The growth-law candidate from the v3.63.0 memo is actually log-log slope
fit between n=2 and n=3:
  log(592/84) / log(40/16) = log(7.048) / log(2.5) = 2.123
This sprint adds n_max=4 (dim H = 80) to extend the log-log fit.

n_max=4 cutoff data:
  N(4) = 4*(4+3)/2 = 14 sectors
  dim H = 2*4*5*6/3 = 80
  dim H^2 = 6400
  N^2 = 196 (max possible block count)

Skip the Lie subalgebra test (cubic in n_commutators on a dim H^2 vector
space; at n_max=3 already showed Lie does NOT close in transition span,
no new structural information from third data point).

Discipline: bit-exact sympy.Rational throughout; no PSLQ; no floats.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Rational, Integer, Matrix, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple


# ----------------------------------------------------------------------
# Helpers (mirrored from v3.63.0 L4 driver)
# ----------------------------------------------------------------------

def serialize_rational(x) -> str:
    if isinstance(x, sp.Rational):
        return f"{x.p}/{x.q}" if x.q != 1 else str(x.p)
    if isinstance(x, sp.Integer):
        return str(int(x))
    return str(x)


def sector_states(triple: FockSpectralTriple) -> Dict[Tuple[int, int], List[int]]:
    out: Dict[Tuple[int, int], List[int]] = {s: [] for s in triple.sectors}
    for i, k in enumerate(triple._state_to_sector):
        out[triple.sectors[k]].append(i)
    return out


def basic_transitions(triple: FockSpectralTriple) -> List[Tuple[Tuple[int, int], Tuple[int, int]]]:
    Lambda = triple.diagonal_part
    A_part = triple.dirac_operator - Lambda  # = kappa * A_graph
    sec_states = sector_states(triple)
    out: List[Tuple[Tuple[int, int], Tuple[int, int]]] = []
    for s_from in triple.sectors:
        for s_to in triple.sectors:
            if s_to == s_from:
                continue
            nz = False
            for i in sec_states[s_to]:
                for j in sec_states[s_from]:
                    if A_part[i, j] != 0:
                        nz = True
                        break
                if nz:
                    break
            if nz:
                out.append((s_from, s_to))
    return out


def compute_closure_per_block(triple: FockSpectralTriple) -> Dict:
    """Per-block closure dimension via rank of (e_t A^k e_s)_{k=0..N-1} slabs."""
    N = triple.dim_H
    Lambda = triple.diagonal_part
    A_part = triple.dirac_operator - Lambda
    sec_states = sector_states(triple)

    # Pre-compute A^k for k = 0, ..., N-1
    print(f"  Building A^k for k = 0, ..., {N-1} (dim H = {N})", flush=True)
    A_powers: List[Matrix] = [sp.eye(N)]
    t_power_start = time.time()
    for k in range(1, N):
        if k % 10 == 0:
            print(f"    A^{k} after {time.time()-t_power_start:.1f}s", flush=True)
        A_powers.append(A_powers[-1] * A_part)
    print(f"  A^k powers complete in {time.time()-t_power_start:.1f}s", flush=True)

    per_block: Dict[str, Dict] = {}
    total = 0
    n_blocks_realized = 0
    t_block_start = time.time()
    block_idx = 0
    n_blocks_total = len(triple.sectors) ** 2

    for s_to in triple.sectors:
        for s_from in triple.sectors:
            block_idx += 1
            d_to = len(sec_states[s_to])
            d_from = len(sec_states[s_from])
            if d_to == 0 or d_from == 0:
                continue
            max_dim = d_to * d_from
            # Build (k = 0, ..., N-1) x (d_to * d_from) matrix
            rows: List[List[sp.Rational]] = []
            for k in range(N):
                Ak = A_powers[k]
                row = []
                for i in sec_states[s_to]:
                    for j in sec_states[s_from]:
                        row.append(Ak[i, j])
                rows.append(row)
            mat = Matrix(rows)
            rk = mat.rank()
            if rk > 0:
                n_blocks_realized += 1
                total += rk
                per_block[f"{s_to}<-{s_from}"] = dict(
                    from_=str(s_from),
                    to=str(s_to),
                    d_to=d_to,
                    d_from=d_from,
                    max_dim=max_dim,
                    realized_rank=rk,
                    fill_fraction=float(rk) / max_dim,
                )
            if block_idx % 20 == 0:
                print(f"    block {block_idx}/{n_blocks_total} after "
                      f"{time.time()-t_block_start:.1f}s "
                      f"(running total = {total})", flush=True)
    print(f"  Block scan complete in {time.time()-t_block_start:.1f}s", flush=True)

    return dict(
        total_closure_dim=total,
        n_blocks_realized=n_blocks_realized,
        per_block=per_block,
    )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def run_closure_for_nmax(n_max: int) -> Dict:
    print("=" * 70)
    print(f"Sprint Q5'-OffDiag-Closure-nmax4 -- n_max = {n_max}")
    print("=" * 70, flush=True)

    t0 = time.time()
    triple = FockSpectralTriple(n_max=n_max)
    print(f"dim H = {triple.dim_H}, n_sectors = {triple.n_sectors}", flush=True)
    print(f"sectors = {triple.sectors}", flush=True)

    trans_pairs = basic_transitions(triple)
    print(f"\n# idempotents: {triple.n_sectors}", flush=True)
    print(f"# single-step transitions: {len(trans_pairs)}", flush=True)

    print("\n--- Per-block closure dimension ---", flush=True)
    t_closure = time.time()
    closure = compute_closure_per_block(triple)
    print(f"\n  Closure computation time: {time.time()-t_closure:.1f}s", flush=True)
    closure_dim = closure["total_closure_dim"]
    print(f"  TOTAL CLOSURE DIM at n_max = {n_max}: {closure_dim}", flush=True)
    print(f"  Number of nonzero blocks: {closure['n_blocks_realized']} of "
          f"{triple.n_sectors ** 2} possible", flush=True)

    sec_states = sector_states(triple)
    sec_dims = {s: len(sec_states[s]) for s in triple.sectors}
    dim_H_sq = triple.dim_H * triple.dim_H
    upper_bound_full = sum(d_to * d_from for d_to in sec_dims.values()
                          for d_from in sec_dims.values())

    print(f"\n  Path-algebra upper bound: {upper_bound_full}", flush=True)
    print(f"  dim(H)^2: {dim_H_sq}", flush=True)
    print(f"  Fill fraction: {closure_dim}/{upper_bound_full} = "
          f"{Rational(closure_dim, upper_bound_full)} "
          f"~ {float(Rational(closure_dim, upper_bound_full)):.4f}", flush=True)

    elapsed = time.time() - t0
    print(f"\nWall time: {elapsed:.1f}s", flush=True)

    return dict(
        n_max=n_max,
        dim_H=triple.dim_H,
        n_sectors=triple.n_sectors,
        sectors=[str(s) for s in triple.sectors],
        n_idempotents=triple.n_sectors,
        n_single_step_transitions=len(trans_pairs),
        transitions=[f"{k[0]} -> {k[1]}" for k in sorted(trans_pairs, key=str)],
        closure_dim=closure_dim,
        n_blocks_realized=closure["n_blocks_realized"],
        per_block=closure["per_block"],
        path_algebra_upper_bound_full=upper_bound_full,
        dim_H_squared=dim_H_sq,
        fill_fraction_str=str(Rational(closure_dim, upper_bound_full)),
        wall_time_s=elapsed,
    )


def main():
    results = {
        "sprint": "Q5p-OffDiag-Closure-nmax4 (T3)",
        "date": "2026-06-06 (continuation of v3.63.0 L4 sub-sprint)",
        "discipline": "bit-exact sympy.Rational throughout; no Lie test",
    }

    res_4 = run_closure_for_nmax(n_max=4)
    results["n_max_4"] = res_4

    # Log-log fit with the v3.63.0 data points
    print("\n" + "=" * 70)
    print("3-point log-log fit of dim A_OD vs dim H")
    print("=" * 70, flush=True)
    panel = [
        ("n_max=2", 16, 84),
        ("n_max=3", 40, 592),
        ("n_max=4", res_4["dim_H"], res_4["closure_dim"]),
    ]
    print(f"  {'cutoff':>10} {'dim H':>6} {'dim A_OD':>10} "
          f"{'ratio':>10} {'log-log':>10}", flush=True)
    prev = None
    for label, dH, dA in panel:
        if prev:
            ratio_H = dH / prev[1]
            ratio_A = dA / prev[2]
            loglog = sp.log(Rational(dA, prev[2])) / sp.log(Rational(dH, prev[1]))
            print(f"  {label:>10} {dH:>6} {dA:>10} {ratio_A/ratio_H:>10.4f} "
                  f"{float(loglog):>10.4f}", flush=True)
        else:
            print(f"  {label:>10} {dH:>6} {dA:>10}", flush=True)
        prev = (label, dH, dA)

    results["panel_log_log"] = [
        dict(label=p[0], dim_H=p[1], dim_A_OD=p[2]) for p in panel
    ]

    out_path = Path(__file__).parent / "data" / "sprint_q5p_offdiag_closure_nmax4.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
