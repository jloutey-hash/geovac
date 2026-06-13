"""S^(4) scoping diagnostic — T1 collapse, T2 expansion inventory, T3 convergence.

Diagnostic-before-engineering pass for the k=4 rung of the S-tower.
Predictions pre-committed in debug/sprint_s4_scoping_memo.md BEFORE this ran
(W3 protocol).

    S^(4) = sum_{n1..n4>=1} I(n1,n2) I(n2,n3) I(n3,n4) phi(n1..n4),
    I(a,b) = 2 min(a,b) - 1 - [a=b],
    phi(n) = 2/(n+3/2)^2 - (1/2)/(n+3/2)^4.

T1: direct CG (production so4_channel_count) vs collapsed kernel, Fractions, N=8.
T2: mechanical 27-term expansion of prod_e (2m_e - 1 - delta_e), shape
    classification by union-find + connected components, bit-exact identity
    check, named-object inventory (the k=4 analog of the k=3 nine-term identity).
T3: O(N) factorized partial sums (numpy float64), monotone lower bounds,
    convergence-class check vs (ln N)^3/N tail, empirical limit estimate,
    Fraction cross-check at N=12.

Output: debug/data/s4_scoping_diag.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geovac.qed_vertex import so4_channel_count  # noqa: E402

OUT: dict = {}
t0 = time.time()

N_EXACT = 8


def phi(n: int) -> Fraction:
    x = Fraction(2 * n + 3, 2)
    return 2 / x ** 2 - Fraction(1, 2) / x ** 4


def I_closed(a: int, b: int) -> int:
    if min(a, b) == 0:
        return 0
    return 2 * min(a, b) - 1 - (1 if a == b else 0)


# ------------------------------------------------------------------ T1
print("T1: k=4 chain collapse against production so4_channel_count ...")
N = N_EXACT
I_direct = {}
for a in range(0, N + 1):
    for b in range(0, N + 1):
        I_direct[(a, b)] = sum(so4_channel_count(a, b, q)
                               for q in range(1, a + b + 1))
pair_exact = all(I_direct[(a, b)] == I_closed(a, b)
                 for a in range(N + 1) for b in range(N + 1))
PH = {n: phi(n) for n in range(0, N + 1)}
S4_direct = sum(I_direct[(a, b)] * I_direct[(b, c)] * I_direct[(c, d)]
                * PH[a] * PH[b] * PH[c] * PH[d]
                for a in range(N + 1) for b in range(N + 1)
                for c in range(N + 1) for d in range(N + 1))
S4_collapsed = sum(I_closed(a, b) * I_closed(b, c) * I_closed(c, d)
                   * PH[a] * PH[b] * PH[c] * PH[d]
                   for a in range(1, N + 1) for b in range(1, N + 1)
                   for c in range(1, N + 1) for d in range(1, N + 1))
t1_pass = pair_exact and (S4_direct == S4_collapsed)
print("  pairwise I closed form exact: %s" % pair_exact)
print("  S4 direct-CG == collapsed (N=%d): %s" % (N, S4_direct == S4_collapsed))
print("  S4(N=%d) = %.10f" % (N, float(S4_collapsed)))
OUT["T1_pair_exact"] = pair_exact
OUT["T1_collapse_bitexact"] = (S4_direct == S4_collapsed)
OUT["S4_partial_N8"] = str(S4_collapsed)

# ------------------------------------------------------------------ T2
print("T2: mechanical expansion of prod_e (2m_e - 1 - delta_e) ...")
EDGES = [(0, 1), (1, 2), (2, 3)]


def component_value(mults: tuple, N: int) -> Fraction:
    """Value of a path-shaped component: classes with phi-multiplicities
    `mults` along the path, min-kernel between consecutive classes."""
    r = len(mults)
    total = Fraction(0)
    if r == 1:
        return sum(PH[n] ** mults[0] for n in range(1, N + 1))
    idx = [1] * r
    rng = range(1, N + 1)
    if r == 2:
        for v0 in rng:
            for v1 in rng:
                total += min(v0, v1) * PH[v0] ** mults[0] * PH[v1] ** mults[1]
    elif r == 3:
        for v0 in rng:
            for v1 in rng:
                m01 = min(v0, v1) * PH[v0] ** mults[0] * PH[v1] ** mults[1]
                for v2 in rng:
                    total += m01 * min(v1, v2) * PH[v2] ** mults[2]
    elif r == 4:
        for v0 in rng:
            for v1 in rng:
                m01 = min(v0, v1) * PH[v0] ** mults[0] * PH[v1] ** mults[1]
                for v2 in rng:
                    m12 = m01 * min(v1, v2) * PH[v2] ** mults[2]
                    for v3 in rng:
                        total += m12 * min(v2, v3) * PH[v3] ** mults[3]
    else:
        raise ValueError("unexpected component size %d" % r)
    return total


NAMES = {(1,): "P", (2,): "Q", (3,): "R3", (4,): "R4",
         (1, 1): "S_min", (1, 2): "M2", (1, 3): "T31", (2, 2): "T22",
         (1, 1, 1): "M3", (1, 1, 2): "M3e", (1, 2, 1): "M3m",
         (1, 1, 1, 1): "M4"}


def classify(choice: tuple) -> tuple:
    """Return sorted tuple of canonical component names for one of the 27
    expansion terms, plus the coefficient."""
    coeff = 1
    parent = list(range(4))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    min_edges = []
    for (e, kind) in zip(EDGES, choice):
        if kind == "m":
            coeff *= 2
            min_edges.append(e)
        elif kind == "one":
            coeff *= -1
        else:  # delta
            coeff *= -1
            ra, rb = find(e[0]), find(e[1])
            parent[max(ra, rb)] = min(ra, rb)
    roots = {}
    for v in range(4):
        roots.setdefault(find(v), []).append(v)
    # adjacency among classes via min-edges
    class_of = {v: find(v) for v in range(4)}
    adj = {}
    for (u, v) in min_edges:
        cu, cv = class_of[u], class_of[v]
        if cu == cv:
            raise RuntimeError("min-edge inside a glued class (cycle?)")
        adj.setdefault(cu, set()).add(cv)
        adj.setdefault(cv, set()).add(cu)
        if len(adj[cu]) > 2 or len(adj[cv]) > 2:
            raise RuntimeError("non-path component encountered")
    # connected components over classes
    comps = []
    seen = set()
    for c in roots:
        if c in seen:
            continue
        stack, comp = [c], []
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            stack.extend(adj.get(x, ()))
        comps.append(comp)
    names = []
    for comp in comps:
        if len(comp) == 1:
            mults = (len(roots[comp[0]]),)
        else:
            ends = [c for c in comp if len(adj.get(c, ())) == 1]
            if len(ends) != 2:
                raise RuntimeError("non-path component (no 2 ends)")
            path = [ends[0]]
            while len(path) < len(comp):
                nxt = [x for x in adj[path[-1]] if x not in path]
                path.append(nxt[0])
            mults = tuple(len(roots[c]) for c in path)
            mults = min(mults, mults[::-1])
        names.append(NAMES[mults])
    return tuple(sorted(names)), coeff


inventory: dict = {}
for choice in product(["m", "one", "delta"], repeat=3):
    key, coeff = classify(choice)
    inventory[key] = inventory.get(key, 0) + coeff
inventory = {k: v for k, v in inventory.items() if v != 0}

VALUES = {}
for key in inventory:
    for name in key:
        if name not in VALUES:
            mults = [m for m, nm in NAMES.items() if nm == name][0]
            VALUES[name] = component_value(mults, N_EXACT)
assembled = Fraction(0)
for key, coeff in inventory.items():
    term = Fraction(coeff)
    for name in key:
        term *= VALUES[name]
    assembled += term
t2_pass = (assembled == S4_collapsed)
print("  expansion identity bit-exact (N=%d): %s" % (N_EXACT, t2_pass))
print("  identity: S4 = " + " + ".join(
    "%+d*%s" % (c, "*".join(k)) for k, c in sorted(inventory.items())))
new_objects = sorted({nm for key in inventory for nm in key}
                     - {"P", "Q", "R3", "S_min", "M2", "M3"})
print("  objects new at k=4: %s" % new_objects)
OUT["T2_identity_bitexact"] = t2_pass
OUT["T2_inventory"] = {"*".join(k): c for k, c in sorted(inventory.items())}
OUT["T2_new_objects"] = new_objects

# ------------------------------------------------------------------ T3
print("T3: factorized O(N) partial sums (float64) ...")


def s4_partial(N: int) -> float:
    n = np.arange(1, N + 1, dtype=np.float64)
    x = n + 1.5
    ph = 2.0 / x ** 2 - 0.5 / x ** 4
    P = ph.sum()
    c_nphi = np.cumsum(n * ph)
    c_phi = np.cumsum(ph)
    A = c_nphi + n * (P - c_phi)           # A[m-1] = sum_n min(n,m) phi(n)
    L = 2.0 * A - P - ph                   # L[m-1] = sum_n I(m,n) phi(n)
    F = ph * L
    SF = F.sum()
    c_nF = np.cumsum(n * F)
    c_F = np.cumsum(F)
    B = c_nF + n * (SF - c_F)              # B[m-1] = sum_n min(m,n) F(n)
    G = 2.0 * B - SF - F                   # G[m-1] = sum_n I(m,n) F(n)
    return float((F * G).sum())


# cross-check the prefix-sum pipeline against the exact Fraction value
N12 = 12
PH12 = {n: phi(n) for n in range(0, N12 + 1)}
S4_exact_12 = sum(I_closed(a, b) * I_closed(b, c) * I_closed(c, d)
                  * PH12[a] * PH12[b] * PH12[c] * PH12[d]
                  for a in range(1, N12 + 1) for b in range(1, N12 + 1)
                  for c in range(1, N12 + 1) for d in range(1, N12 + 1))
xerr = abs(s4_partial(12) - float(S4_exact_12))
print("  pipeline vs exact Fraction at N=12: |diff| = %.3e" % xerr)
OUT["T3_pipeline_crosscheck_N12"] = xerr

GRID = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 4 * 10 ** 6]
vals = {}
for Ng in GRID:
    tA = time.time()
    vals[Ng] = s4_partial(Ng)
    print("  S4(N=%-8d) = %.8f   (%.1f s)" % (Ng, vals[Ng], time.time() - tA))
OUT["T3_partial_sums"] = {str(k): v for k, v in vals.items()}

# monotonicity + convergence class: tail(N) ~ c (ln N)^3 / N predicts
# decade ratios r = tail(10N)/tail(N) = (ln 10N / ln N)^3 / 10.
mono = all(vals[GRID[i + 1]] > vals[GRID[i]] for i in range(len(GRID) - 1))
print("  monotone increasing: %s" % mono)
# fit S_inf and c on the last three points with tail model c*(ln N)^3/N
import itertools  # noqa: E402


def fit_tail(pts, j):
    # least squares for S_inf, c in  S(N) = S_inf - c (ln N)^j / N
    A_ = np.array([[1.0, -(np.log(Ng) ** j) / Ng] for Ng, _ in pts])
    b_ = np.array([v for _, v in pts])
    sol, *_ = np.linalg.lstsq(A_, b_, rcond=None)
    resid = float(np.max(np.abs(A_ @ sol - b_)))
    return float(sol[0]), float(sol[1]), resid


pts = [(Ng, vals[Ng]) for Ng in GRID[-3:]]
fits = {}
for j in (2, 3, 4):
    S_inf, c, resid = fit_tail(pts, j)
    fits[j] = (S_inf, c, resid)
    print("  tail model (ln N)^%d/N: S_inf = %.6f, c = %.4f, max-resid = %.2e"
          % (j, S_inf, c, resid))
best_j = min(fits, key=lambda j: fits[j][2])
S_inf_est, c_est, _ = fits[best_j]
tail_4e6 = c_est * np.log(4e6) ** best_j / 4e6
print("  best model j=%d; empirical S^(4) ~= %.4f; est. tail at N=4e6 ~= %.4f"
      % (best_j, S_inf_est, tail_4e6))
OUT["T3_monotone"] = mono
OUT["T3_fits"] = {str(j): {"S_inf": f[0], "c": f[1], "resid": f[2]}
                  for j, f in fits.items()}
OUT["T3_best_j"] = best_j
OUT["T3_S4_empirical"] = S_inf_est
OUT["T3_tail_est_4e6"] = float(tail_4e6)

OUT["runtime_s"] = time.time() - t0
out_path = Path(__file__).resolve().parent / "data" / "s4_scoping_diag.json"
out_path.write_text(json.dumps(OUT, indent=2))
print("wrote %s  (%.1f s total)" % (out_path, OUT["runtime_s"]))
