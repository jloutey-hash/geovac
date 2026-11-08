"""Sprint Q5'-T-Path-Generating (T2) — closed-form generating function for
the chirality-weighted two-step path count N^(2)_{s' -> s} on the
OffDiag substrate.

L3 (Bridge-Id, v3.63.0) identified the constructive ingredient
T_path(n_max=3, palindrome (e_2, e_3, e_3, e_2), interior) = 24 and
T_path(n_max=2, boundary) = 8 of the bridge identity
drift_{n_max>=3} = -kappa^4. The multi-year follow-on flagged in the L3
memo §6 and the v3.63.0 umbrella memo §5 is a CLOSED-FORM GENERATING
FUNCTION for these path counts as a function of (n_max, s, s').

Definition (from Bridge-Id memo §3.3):
    eta(T_{s' -> s}) = Tr(gamma * D * T_{s' -> s})
                     = kappa^2 * Tr(gamma * A * e_s * A * e_{s'})
                     =: kappa^2 * N^(2)_{s' -> s}

The chirality-weighted two-step path count N^(2)_{s' -> s} is an integer
counting chirality-weighted two-step E1 paths from sector s' back to s'
through intermediate sector s.

Bit-exact panel at n_max = 1, 2, 3, 4, 5:
    - Build FockSpectralTriple at each cutoff
    - For each adjacent sector pair (s', s) with eta(T_{s' -> s}) != 0,
      compute N^(2)_{s' -> s} = Tr(gamma * A * e_s * A * e_{s'})
    - Look for closed-form pattern as a function of (n', l', n, l, n_max)
    - Identify any "boundary" vs "interior" structural splits

Discipline: bit-exact sympy.Rational throughout; no floats; no PSLQ.
"""
from __future__ import annotations

import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Integer, Rational, Matrix, Symbol, sympify, simplify, factor, expand

from geovac.spectral_triple import FockSpectralTriple


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def sector_states(triple: FockSpectralTriple) -> Dict[Tuple[int, int], List[int]]:
    out: Dict[Tuple[int, int], List[int]] = {s: [] for s in triple.sectors}
    for i, k in enumerate(triple._state_to_sector):
        out[triple.sectors[k]].append(i)
    return out


def compute_n2_for_pair(
    triple: FockSpectralTriple,
    s_from: Tuple[int, int],
    s_to: Tuple[int, int],
) -> int:
    """N^(2)_{s' -> s} = Tr(gamma * A * e_s * A * e_{s'}) where A = D - Lambda
    is the off-diagonal Dirac (i.e. kappa * A_graph).

    Bit-exactly returned as a Python int (since the chirality-weighted
    two-step count is integer-valued by construction).
    """
    N = triple.dim_H
    Lambda = triple.diagonal_part
    A_part = triple.dirac_operator - Lambda  # = kappa * A_graph
    gamma = triple.grading

    sec_states = sector_states(triple)

    # Build e_s, e_{s'} as diagonal idempotents
    from sympy import zeros as sp_zeros
    e_s = sp_zeros(N, N)
    for i in sec_states[s_to]:
        e_s[i, i] = Integer(1)
    e_sp = sp_zeros(N, N)
    for i in sec_states[s_from]:
        e_sp[i, i] = Integer(1)

    # Tr(gamma * A * e_s * A * e_{s'})
    # The chirality-weighted two-step count is the trace WITHOUT
    # the kappa^2 prefactor on A. Since A_part = kappa * A_graph,
    # Tr(gamma * A_part * e_s * A_part * e_{s'})
    #     = kappa^2 * Tr(gamma * A_graph * e_s * A_graph * e_{s'})
    M = gamma * A_part * e_s * A_part * e_sp
    trace = sum(M[i, i] for i in range(N))
    # Strip kappa^2 factor
    kappa = Rational(-1, 16)
    kappa2 = kappa ** 2
    if trace == 0:
        return 0
    n2 = trace / kappa2
    # Should be integer
    if not n2.is_integer:
        raise RuntimeError(f"Non-integer N^(2): {trace}/{kappa2} = {n2}")
    return int(n2)


# ----------------------------------------------------------------------
# Main panel
# ----------------------------------------------------------------------


def panel_for_nmax(n_max: int) -> Dict:
    triple = FockSpectralTriple(n_max=n_max)
    print(f"\nn_max = {n_max}: dim H = {triple.dim_H}, "
          f"n_sectors = {triple.n_sectors}", flush=True)
    print(f"  sectors: {[str(s) for s in triple.sectors]}", flush=True)

    # Find adjacent sector pairs (s', s) with N^(2) != 0
    sec_states_d = sector_states(triple)
    pairs = []
    for s_from in triple.sectors:
        for s_to in triple.sectors:
            n2 = compute_n2_for_pair(triple, s_from, s_to)
            if n2 != 0:
                pairs.append((s_from, s_to, n2))

    # Sort by source sector lexicographically
    pairs.sort(key=lambda x: (x[0], x[1]))

    print(f"  # non-zero N^(2)_{{s'->s}} pairs: {len(pairs)}", flush=True)
    print(f"  {'(n,l) from':>12} -> {'(n,l) to':>12}  N^(2)", flush=True)
    for s_from, s_to, n2 in pairs:
        print(f"  {str(s_from):>12} -> {str(s_to):>12}  {n2:+d}", flush=True)

    return dict(
        n_max=n_max,
        dim_H=triple.dim_H,
        n_sectors=triple.n_sectors,
        sectors=[list(s) for s in triple.sectors],
        n_pairs_realized=len(pairs),
        pairs=[
            dict(
                s_from=list(s_from),
                s_to=list(s_to),
                N2=n2,
                n_from=s_from[0],
                l_from=s_from[1],
                n_to=s_to[0],
                l_to=s_to[1],
            )
            for s_from, s_to, n2 in pairs
        ],
        total_N2=sum(abs(p[2]) for p in pairs),
        total_N2_signed=sum(p[2] for p in pairs),
    )


# ----------------------------------------------------------------------
# Closed-form pattern detection
# ----------------------------------------------------------------------


def detect_closed_form(panels: List[Dict]) -> Dict:
    """Look for a closed-form polynomial N^(2)_{s' -> s}(n_max, n', l', n, l).

    Hypothesis 1 (cutoff-independent): N^(2)_{s' -> s} depends only on
    (n', l', n, l), not on n_max. Test by checking values at common
    pairs across cutoffs.

    Hypothesis 2 (selection rule): only |n - n'| = 1 and |l - l'| <= 1
    transitions appear (E1 dipole selection rule from sphere harmonics).

    Hypothesis 3 (closed form): N^(2)_{(n', l') -> (n, l)} is a polynomial
    in (n, l) for fixed transition type.
    """
    # Collect all distinct pair values across cutoffs
    pair_values: Dict[Tuple[int, int, int, int], List[Tuple[int, int]]] = defaultdict(list)
    for panel in panels:
        for entry in panel["pairs"]:
            key = (entry["n_from"], entry["l_from"], entry["n_to"], entry["l_to"])
            pair_values[key].append((panel["n_max"], entry["N2"]))

    # Check hypothesis 1: cutoff-independence
    cutoff_independent = True
    cutoff_dep_pairs = []
    for key, vals in pair_values.items():
        n2_vals = set(v[1] for v in vals)
        if len(n2_vals) > 1:
            cutoff_independent = False
            cutoff_dep_pairs.append(dict(
                pair=key,
                cutoff_values=vals,
            ))

    # Hypothesis 1 result
    h1 = dict(
        statement="N^(2)_{s' -> s} depends only on (s', s), not on n_max",
        passed=cutoff_independent,
        cutoff_dependent_pairs=cutoff_dep_pairs,
    )

    # Hypothesis 2: selection rule
    selection_violations = []
    for key in pair_values:
        n_from, l_from, n_to, l_to = key
        if abs(n_to - n_from) > 1:
            selection_violations.append(dict(pair=key, type="|Delta n| > 1"))
        if abs(l_to - l_from) > 1:
            selection_violations.append(dict(pair=key, type="|Delta l| > 1"))
    h2 = dict(
        statement="Only |Delta n| <= 1 and |Delta l| <= 1 transitions appear",
        passed=len(selection_violations) == 0,
        violations=selection_violations,
    )

    # Group pairs by transition type
    # Type A: (n, l) -> (n, l+1)   (same shell, l up)
    # Type B: (n, l) -> (n, l-1)   (same shell, l down)
    # Type C: (n, l) -> (n+1, l)   (shell up, same l)
    # Type D: (n, l) -> (n-1, l)   (shell down, same l)
    # Type E: (n, l) -> (n+1, l+1) and (n+1, l-1) (shell up, l +/- 1)
    # Type F: (n, l) -> (n-1, l-1) etc (shell down, l +/- 1)
    groups: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for key, vals in pair_values.items():
        n_from, l_from, n_to, l_to = key
        dn = n_to - n_from
        dl = l_to - l_from
        # Canonical N^(2) value
        n2 = vals[0][1] if cutoff_independent else None
        if n2 is None:
            continue
        groups[f"dn={dn},dl={dl}"].append((n_from, l_from, n2))

    # Try to fit a closed-form polynomial in (n, l) per group
    h3_results = {}
    for group_label, entries in groups.items():
        if len(entries) < 2:
            h3_results[group_label] = dict(
                n_entries=len(entries),
                entries=[dict(n=e[0], l=e[1], N2=e[2]) for e in entries],
                closed_form="(insufficient data)",
            )
            continue
        # Try ansatz: N^(2) = polynomial in (n, l). Sample N^(2) values
        # and fit symbolic polynomial degree (1, 2).
        # Start with: try N^(2) = a*n + b*l + c
        n_sym, l_sym = Symbol("n"), Symbol("l")
        # Build system
        if len(entries) >= 3:
            # Try degree-1 ansatz: a*n + b*l + c
            a, b, c = Symbol("a"), Symbol("b"), Symbol("c")
            eqs = []
            for n, l, n2 in entries[:3]:
                eqs.append(sp.Eq(a * n + b * l + c, n2))
            sol = sp.solve(eqs, [a, b, c])
            if sol:
                # Verify on remaining entries
                fit = sol[a] * n_sym + sol[b] * l_sym + sol[c]
                all_match = True
                for n, l, n2 in entries:
                    if fit.subs({n_sym: n, l_sym: l}) != n2:
                        all_match = False
                        break
                if all_match:
                    h3_results[group_label] = dict(
                        n_entries=len(entries),
                        entries=[dict(n=e[0], l=e[1], N2=e[2]) for e in entries],
                        closed_form_degree1=str(fit),
                        closed_form=str(fit),
                        verified_all_entries=True,
                    )
                    continue
            # Try degree-2 ansatz: a*n^2 + b*n*l + c*l^2 + d*n + e*l + f
            if len(entries) >= 6:
                aa, bb, cc, dd, ee, ff = sp.symbols("aa bb cc dd ee ff")
                eqs = []
                for n, l, n2 in entries[:6]:
                    eqs.append(sp.Eq(aa * n * n + bb * n * l + cc * l * l
                                     + dd * n + ee * l + ff, n2))
                sol = sp.solve(eqs, [aa, bb, cc, dd, ee, ff])
                if sol:
                    fit = (sol[aa] * n_sym ** 2 + sol[bb] * n_sym * l_sym
                           + sol[cc] * l_sym ** 2 + sol[dd] * n_sym
                           + sol[ee] * l_sym + sol[ff])
                    all_match = True
                    for n, l, n2 in entries:
                        if fit.subs({n_sym: n, l_sym: l}) != n2:
                            all_match = False
                            break
                    if all_match:
                        h3_results[group_label] = dict(
                            n_entries=len(entries),
                            entries=[dict(n=e[0], l=e[1], N2=e[2]) for e in entries],
                            closed_form_degree2=str(sp.simplify(fit)),
                            closed_form=str(sp.simplify(fit)),
                            verified_all_entries=True,
                        )
                        continue
        h3_results[group_label] = dict(
            n_entries=len(entries),
            entries=[dict(n=e[0], l=e[1], N2=e[2]) for e in entries],
            closed_form="(no closed form at degree <= 2)",
            verified_all_entries=False,
        )

    return dict(
        h1=h1,
        h2=h2,
        h3=h3_results,
    )


# ----------------------------------------------------------------------
# Generating function for total N^(2) at each cutoff
# ----------------------------------------------------------------------


def generating_function(panels: List[Dict]) -> Dict:
    """Compute the cumulative chirality-weighted two-step path count
    T_path(n_max) = sum_{pairs at cutoff n_max} N^(2)_{s' -> s}.

    Look for a closed-form polynomial in n_max.
    """
    samples = []
    for panel in panels:
        n_max = panel["n_max"]
        T_path_total = sum(p["N2"] for p in panel["pairs"])
        T_path_abs = sum(abs(p["N2"]) for p in panel["pairs"])
        samples.append((n_max, T_path_total, T_path_abs))
        print(f"  n_max = {n_max}: T_path_signed = {T_path_total}, "
              f"T_path_abs = {T_path_abs}", flush=True)

    # Try to fit T_path_signed = polynomial in n_max
    n_sym = Symbol("n")
    closed_signed = None
    if len(samples) >= 4:
        # Try degree-3 ansatz
        a, b, c, d = sp.symbols("a b c d")
        eqs = []
        for n_max, T, _ in samples[:4]:
            eqs.append(sp.Eq(a * n_max ** 3 + b * n_max ** 2 + c * n_max + d, T))
        sol = sp.solve(eqs, [a, b, c, d])
        if sol:
            fit = (sol[a] * n_sym ** 3 + sol[b] * n_sym ** 2
                   + sol[c] * n_sym + sol[d])
            ok = all(
                fit.subs({n_sym: n_max}) == T for n_max, T, _ in samples
            )
            if ok:
                closed_signed = str(sp.simplify(fit))

    closed_abs = None
    if len(samples) >= 4:
        a, b, c, d = sp.symbols("a b c d")
        eqs = []
        for n_max, _, T in samples[:4]:
            eqs.append(sp.Eq(a * n_max ** 3 + b * n_max ** 2 + c * n_max + d, T))
        sol = sp.solve(eqs, [a, b, c, d])
        if sol:
            fit = (sol[a] * n_sym ** 3 + sol[b] * n_sym ** 2
                   + sol[c] * n_sym + sol[d])
            ok = all(
                fit.subs({n_sym: n_max}) == T for n_max, _, T in samples
            )
            if ok:
                closed_abs = str(sp.simplify(fit))

    return dict(
        samples=[dict(n_max=s[0], T_path_signed=s[1], T_path_abs=s[2])
                 for s in samples],
        closed_form_signed=closed_signed or "(no degree-3 polynomial fit)",
        closed_form_abs=closed_abs or "(no degree-3 polynomial fit)",
    )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------


def main() -> None:
    t0 = time.time()

    print("=" * 70)
    print("Sprint Q5'-T-Path-Generating (T2)")
    print("=" * 70)

    panels = []
    for n_max in (1, 2, 3, 4, 5):
        panels.append(panel_for_nmax(n_max))

    print("\n" + "=" * 70)
    print("Closed-form pattern detection")
    print("=" * 70, flush=True)
    closed_form = detect_closed_form(panels)
    print(f"\nHypothesis 1 (cutoff-independence): {'PASS' if closed_form['h1']['passed'] else 'FAIL'}",
          flush=True)
    if not closed_form["h1"]["passed"]:
        for entry in closed_form["h1"]["cutoff_dependent_pairs"]:
            print(f"  cutoff-dependent pair: {entry}", flush=True)
    print(f"\nHypothesis 2 (E1 selection |dn|<=1, |dl|<=1): "
          f"{'PASS' if closed_form['h2']['passed'] else 'FAIL'}", flush=True)
    print(f"\nHypothesis 3 (closed form per (dn, dl) group):", flush=True)
    for group_label, info in closed_form["h3"].items():
        print(f"  {group_label}:", flush=True)
        if "closed_form_degree1" in info:
            print(f"    deg-1 closed form: N^(2) = {info['closed_form_degree1']}",
                  flush=True)
        elif "closed_form_degree2" in info:
            print(f"    deg-2 closed form: N^(2) = {info['closed_form_degree2']}",
                  flush=True)
        else:
            print(f"    {info['closed_form']}", flush=True)
        for e in info["entries"][:5]:
            print(f"      ({e['n']}, {e['l']}): N^(2) = {e['N2']}", flush=True)

    print("\n" + "=" * 70)
    print("Total generating function T_path(n_max)")
    print("=" * 70, flush=True)
    gen = generating_function(panels)
    print(f"\nClosed form T_path_signed(n_max) = {gen['closed_form_signed']}",
          flush=True)
    print(f"Closed form T_path_abs(n_max) = {gen['closed_form_abs']}", flush=True)

    wall = time.time() - t0
    print(f"\nWall: {wall:.2f} s", flush=True)

    out = dict(
        sprint="Q5p-T-Path-Generating (T2)",
        date="2026-06-06 (continuation of v3.63.0 L3 sub-sprint)",
        purpose="Closed-form generating function for chirality-weighted two-step path count N^(2)_{s' -> s} on the OffDiag substrate, extending v3.63.0 L3 Bridge-Id panel.",
        panels=panels,
        closed_form_per_transition=closed_form,
        generating_function=gen,
        wall_seconds=wall,
    )

    out_path = Path(__file__).parent / "data" / "sprint_q5p_t_path_generating.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWritten: {out_path}", flush=True)


if __name__ == "__main__":
    main()
